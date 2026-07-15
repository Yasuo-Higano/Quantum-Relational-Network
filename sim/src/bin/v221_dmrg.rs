//! v22.1 計算層: 重電荷 Schwinger の 2 サイト DMRG — Rust 移植 (PROMPT/4: 重い計算は Rust)
//!
//! explore/dmrg_schwinger.py + dmrg_heavy.py の完全な置き換え。numpy 版は einsum の
//! 一時割当てが律速で N=64 の 1 対に 20 時間超 — 本版は割当てゼロの内側ループで
//! 同一物理を計算する (期待 ~20 倍)。
//!
//! 模型 (python 版と同一規約): H = −x Σ (c†c + h.c.) + Σ_pairs w(k,l) q_k q_l,
//! w(k,l) = N−1−max(k,l) + λ_Q, q_k = qf·(n_k − [k odd]), λ_Q = 10 (全電荷罰則 —
//! run1 の教訓: 開鎖端の電荷は自由で荷電端状態が E₁ を偽装する)。
//! MPO 走和構成 D=5: W = [[I, q, cl·q²], [0, I, 2cl·q], [0, 0, I]] ⊕ ホップ,
//! cl = N−1−l + λ_Q。励起状態 = 直交罰則 λ_o|ψ₀⟩⟨ψ₀| (python 版と同一)。
//!
//! マルチスレッド: 6 ジョブ (アンカー N=10 × 2 + 本測定 N=64 × 4) を独立スレッドで
//! 実行。各ジョブは内部完全直列・固定シード — スレッド分割で結果は変わらない (決定的)。
//!
//! 自己ゲート (PASS/FAIL 内蔵):
//!   [G1] N=8, x=2, q=1: DMRG E₀ = 全空間 ED ± 1e-8 (dim 256 の jacobi)
//!   [G2] N=8 罰則法 E₁ = ED E₁ ± 1e-6
//!   [G3/G4] アンカー N=10 (q=1,2): DMRG E₀/E₁ = 全空間 ED (dim 1024) ± 1e-7/1e-5
//! 出力: results/v221_dmrg_gaps.json (v221_dmrgex が監査・判定する一次データ)

use std::io::Write as _;
use uft_sim::*;

const LAMQ: f64 = 10.0;
const LAM_O: f64 = 50.0;

// ---------------- MPO (D=5, 走和構成) ----------------
// W[l]: (dl, dr, 2, 2) — 端は 1×5 / 5×1。添字 [p][r][s_ket][s_bra]
struct Mpo {
    n: usize,
    w: Vec<Vec<f64>>, // per site: dl*dr*4, index ((p*dr + r)*2 + s)*2 + t
    dl: Vec<usize>,
    dr: Vec<usize>,
}

fn build_mpo(n: usize, x: f64, qf: f64, qpen: f64) -> Mpo {
    let mut w_all = Vec::new();
    let mut dls = Vec::new();
    let mut drs = Vec::new();
    for l in 0..n {
        let bg = if l % 2 == 1 { 1.0 } else { 0.0 };
        // q = qf*(n_op − bg·I): 対角 (s==t): q[0] = −qf·bg, q[1] = qf·(1−bg)
        let qd = [-qf * bg, qf * (1.0 - bg)];
        let cl = (n as f64 - 1.0 - l as f64).max(0.0) + qpen;
        let (dl, dr) = (if l == 0 { 1 } else { 5 }, if l == n - 1 { 1 } else { 5 });
        let mut w = vec![0.0f64; dl * dr * 4];
        // フル 5×5 を組んでから端をスライス
        let mut full = vec![0.0f64; 5 * 5 * 4];
        {
            let mut set = |p: usize, r: usize, s: usize, t: usize, v: f64| {
                full[((p * 5 + r) * 2 + s) * 2 + t] += v;
            };
            for s in 0..2 {
                set(0, 0, s, s, 1.0); // I
                set(4, 4, s, s, 1.0); // I
                set(0, 4, s, s, cl * qd[s] * qd[s]); // cl q²
                set(0, 1, s, s, qd[s]); // q
                set(1, 1, s, s, 1.0); // I
                set(1, 4, s, s, 2.0 * cl * qd[s]); // 2 cl q
            }
            // ホップ: W[0,2] = −x·c† (|1⟩⟨0| → s_ket=1, s_bra=0), W[2,4] = c
            set(0, 2, 1, 0, -x);
            set(2, 4, 0, 1, 1.0);
            set(0, 3, 0, 1, -x); // −x·c
            set(3, 4, 1, 0, 1.0); // c†
        }
        for p in 0..dl {
            let pp = if l == 0 { 0 } else { p };
            for r in 0..dr {
                let rr = if l == n - 1 { 4 } else { r };
                for s in 0..2 {
                    for t in 0..2 {
                        w[((p * dr + r) * 2 + s) * 2 + t] = full[((pp * 5 + rr) * 2 + s) * 2 + t];
                    }
                }
            }
        }
        w_all.push(w);
        dls.push(dl);
        drs.push(dr);
    }
    Mpo {
        n,
        w: w_all,
        dl: dls,
        dr: drs,
    }
}

// ---------------- MPS ----------------
#[derive(Clone)]
struct Mps {
    a: Vec<Vec<f64>>, // per site: dl*2*dr, index (al*2 + s)*dr + ar
    bond: Vec<usize>, // n+1
}

fn init_mps(n: usize, chi: usize, seed: u64, rand_init: bool) -> Mps {
    // Néel バイアス行 + 全 χ ランダム (右正準) — run5 の教訓: 小ボンド開始 (≤4) は
    // 積状態のランク障壁で χ 成長が阻害され x=16 系が −578.5 に停滞する
    // (python の全 χ ランダムは −582.3 到達)。ランク自由 + 物理方向バイアスの併用。
    let bond: Vec<usize> = (0..=n)
        .map(|l| {
            let a = 1usize << l.min(24);
            let b = 1usize << (n - l).min(24);
            chi.min(a).min(b)
        })
        .collect();
    let mut rng = Rng::new(seed);
    let mut a_all = Vec::new();
    for l in 0..n {
        let (dl, dr) = (bond[l], bond[l + 1]);
        let neel = if l % 2 == 1 { 1usize } else { 0 };
        // 行直交ランダム (右正準): dl 行 × (2dr) 列を Gram-Schmidt
        let cols = 2 * dr;
        let mut m = vec![0.0f64; dl * cols];
        for v in m.iter_mut() {
            *v = 0.05 * rng.gauss();
        }
        // 行 0 に積状態成分を注入 (GS が最初に正規化するので主成分として残る)。
        // rand_init = true (励起状態用) では注入せず純乱数 — ψ₀ との Krylov 崩壊を避ける
        if !rand_init {
            m[neel * dr] += 1.0;
        }
        for i in 0..dl {
            for j in 0..i {
                let mut ov = 0.0;
                for k in 0..cols {
                    ov += m[i * cols + k] * m[j * cols + k];
                }
                for k in 0..cols {
                    m[i * cols + k] -= ov * m[j * cols + k];
                }
            }
            let nr: f64 = (0..cols)
                .map(|k| m[i * cols + k] * m[i * cols + k])
                .sum::<f64>()
                .sqrt();
            if nr > 1e-12 {
                for k in 0..cols {
                    m[i * cols + k] /= nr;
                }
            }
        }
        // (dl, 2dr) 行優先 → a[(al*2+s)*dr + ar] = m[al, s*dr + ar]
        let mut a = vec![0.0f64; dl * 2 * dr];
        for al in 0..dl {
            for s in 0..2 {
                for ar in 0..dr {
                    a[(al * 2 + s) * dr + ar] = m[al * cols + s * dr + ar];
                }
            }
        }
        a_all.push(a);
    }
    Mps { a: a_all, bond }
}

// ---------------- 環境 ----------------
// Lc[l]: (a_ket, p, c_bra) 添字 (a*dp + p)*dc + c — サイト 0..l−1 を縮約済み
// Rc[l]: (a_ket, p, c_bra) — サイト l..n−1 を縮約済み (左端で評価)
struct Env {
    da: usize,
    dp: usize,
    dc: usize,
    v: Vec<f64>,
}

fn env_unit() -> Env {
    Env {
        da: 1,
        dp: 1,
        dc: 1,
        v: vec![1.0],
    }
}

// upd_L: Lc[l+1][b,q,d] = Σ Lc[l][a,p,c]·A[a,s,b]·W[p,q,s,t]·A[c,t,d]
fn upd_l(lc: &Env, a: &[f64], w: &[f64], dl: usize, dr: usize, wdl: usize, wdr: usize) -> Env {
    let (da, dp, dc) = (lc.da, lc.dp, lc.dc);
    debug_assert_eq!(da, dl);
    debug_assert_eq!(dc, dl);
    debug_assert_eq!(dp, wdl);
    // T1[p,c,s,b] = Σ_a Lc[a,p,c]·A[(a,s),b]
    let mut t1 = vec![0.0f64; dp * dc * 2 * dr];
    for aa in 0..da {
        for p in 0..dp {
            for c in 0..dc {
                let lv = lc.v[(aa * dp + p) * dc + c];
                if lv == 0.0 {
                    continue;
                }
                for s in 0..2 {
                    let arow = (aa * 2 + s) * dr;
                    let trow = ((p * dc + c) * 2 + s) * dr;
                    for b in 0..dr {
                        t1[trow + b] += lv * a[arow + b];
                    }
                }
            }
        }
    }
    // T2[q,t,c,b] = Σ_{p,s} T1[p,c,s,b]·W[(p,q),(s,t)]
    let mut t2 = vec![0.0f64; wdr * 2 * dc * dr];
    for p in 0..wdl {
        for q in 0..wdr {
            for s in 0..2 {
                for t in 0..2 {
                    let wv = w[((p * wdr + q) * 2 + s) * 2 + t];
                    if wv == 0.0 {
                        continue;
                    }
                    for c in 0..dc {
                        let trow = ((p * dc + c) * 2 + s) * dr;
                        let orow = ((q * 2 + t) * dc + c) * dr;
                        for b in 0..dr {
                            t2[orow + b] += wv * t1[trow + b];
                        }
                    }
                }
            }
        }
    }
    // Lc'[b,q,d] = Σ_{c,t} T2[q,t,c,b]·A[(c,t),d]
    let mut out = vec![0.0f64; dr * wdr * dr];
    for q in 0..wdr {
        for t in 0..2 {
            for c in 0..dc {
                let trow = ((q * 2 + t) * dc + c) * dr;
                let arow = (c * 2 + t) * dr;
                for b in 0..dr {
                    let tv = t2[trow + b];
                    if tv == 0.0 {
                        continue;
                    }
                    let orow = (b * wdr + q) * dr;
                    for d in 0..dr {
                        out[orow + d] += tv * a[arow + d];
                    }
                }
            }
        }
    }
    Env {
        da: dr,
        dp: wdr,
        dc: dr,
        v: out,
    }
}

// upd_R: Rc[l][a,p,c] = Σ Rc[l+1][b,q,d]·A[(a,s),b]·W[(p,q),(s,t)]·A[(c,t),d]
fn upd_r(rc: &Env, a: &[f64], w: &[f64], dl: usize, dr: usize, wdl: usize, wdr: usize) -> Env {
    let (db, dq, dd) = (rc.da, rc.dp, rc.dc);
    debug_assert_eq!(db, dr);
    debug_assert_eq!(dq, wdr);
    // T1[q,d,s,a] = Σ_b Rc[b,q,d]·A[(a,s),b]
    let mut t1 = vec![0.0f64; dq * dd * 2 * dl];
    for b in 0..db {
        for q in 0..dq {
            for d in 0..dd {
                let rv = rc.v[(b * dq + q) * dd + d];
                if rv == 0.0 {
                    continue;
                }
                for aa in 0..dl {
                    for s in 0..2 {
                        t1[((q * dd + d) * 2 + s) * dl + aa] += rv * a[(aa * 2 + s) * dr + b];
                    }
                }
            }
        }
    }
    // T2[p,t,d,a] = Σ_{q,s} T1[q,d,s,a]·W[(p,q),(s,t)]
    let mut t2 = vec![0.0f64; wdl * 2 * dd * dl];
    for p in 0..wdl {
        for q in 0..wdr {
            for s in 0..2 {
                for t in 0..2 {
                    let wv = w[((p * wdr + q) * 2 + s) * 2 + t];
                    if wv == 0.0 {
                        continue;
                    }
                    for d in 0..dd {
                        let trow = ((q * dd + d) * 2 + s) * dl;
                        let orow = ((p * 2 + t) * dd + d) * dl;
                        for aa in 0..dl {
                            t2[orow + aa] += wv * t1[trow + aa];
                        }
                    }
                }
            }
        }
    }
    // Rc'[a,p,c] = Σ_{d,t} T2[p,t,d,a]·A[(c,t),d]
    let mut out = vec![0.0f64; dl * wdl * dl];
    for p in 0..wdl {
        for t in 0..2 {
            for d in 0..dd {
                let trow = ((p * 2 + t) * dd + d) * dl;
                for c in 0..dl {
                    let av = a[(c * 2 + t) * dr + d];
                    if av == 0.0 {
                        continue;
                    }
                    for aa in 0..dl {
                        out[(aa * wdl + p) * dl + c] += t2[trow + aa] * av;
                    }
                }
            }
        }
    }
    Env {
        da: dl,
        dp: wdl,
        dc: dl,
        v: out,
    }
}

// ---------------- 2 サイト DMRG ----------------
struct OrthoEnv {
    // oL[l]: (a_cur, c_phi), oR[l]: (a_cur, c_phi) — 左端評価
    ol: Vec<Vec<f64>>,
    or_: Vec<Vec<f64>>,
    dims_ol: Vec<(usize, usize)>,
    dims_or: Vec<(usize, usize)>,
}

fn ov_l(
    prev: &[f64],
    pa: usize,
    pc: usize,
    a: &[f64],
    dl: usize,
    dr: usize,
    phi: &[f64],
    pdl: usize,
    pdr: usize,
) -> (Vec<f64>, usize, usize) {
    // out[b,d] = Σ_{a,c,s} prev[a,c]·A[(a,s),b]·φ[(c,s),d]
    debug_assert_eq!(pa, dl);
    debug_assert_eq!(pc, pdl);
    // T[c,s,b] = Σ_a prev[a,c]·A[(a,s),b]
    let mut t1 = vec![0.0f64; pc * 2 * dr];
    for aa in 0..pa {
        for c in 0..pc {
            let pv = prev[aa * pc + c];
            if pv == 0.0 {
                continue;
            }
            for s in 0..2 {
                let arow = (aa * 2 + s) * dr;
                let trow = (c * 2 + s) * dr;
                for b in 0..dr {
                    t1[trow + b] += pv * a[arow + b];
                }
            }
        }
    }
    // out[b,d] = Σ_{c,s} T[c,s,b]·φ[(c,s),d]
    let mut out = vec![0.0f64; dr * pdr];
    for c in 0..pc {
        for s in 0..2 {
            let trow = (c * 2 + s) * dr;
            let prow = (c * 2 + s) * pdr;
            for b in 0..dr {
                let tv = t1[trow + b];
                if tv == 0.0 {
                    continue;
                }
                for d in 0..pdr {
                    out[b * pdr + d] += tv * phi[prow + d];
                }
            }
        }
    }
    (out, dr, pdr)
}

fn ov_r(
    prev: &[f64],
    pb: usize,
    pd: usize,
    a: &[f64],
    dl: usize,
    dr: usize,
    phi: &[f64],
    pdl: usize,
    pdr: usize,
) -> (Vec<f64>, usize, usize) {
    // out[a,c] = Σ_{b,d,s} prev[b,d]·A[(a,s),b]·φ[(c,s),d]
    debug_assert_eq!(pb, dr);
    debug_assert_eq!(pd, pdr);
    // T[a,s,d] = Σ_b A[(a,s),b]·prev[b,d]
    let mut t1 = vec![0.0f64; dl * 2 * pd];
    for aa in 0..dl {
        for s in 0..2 {
            let arow = (aa * 2 + s) * dr;
            let trow = (aa * 2 + s) * pd;
            for b in 0..dr {
                let av = a[arow + b];
                if av == 0.0 {
                    continue;
                }
                for d in 0..pd {
                    t1[trow + d] += av * prev[b * pd + d];
                }
            }
        }
    }
    let mut out = vec![0.0f64; dl * pdl];
    for aa in 0..dl {
        for s in 0..2 {
            let trow = (aa * 2 + s) * pd;
            for c in 0..pdl {
                let prow = (c * 2 + s) * pdr;
                let mut acc = 0.0;
                for d in 0..pd {
                    acc += t1[trow + d] * phi[prow + d];
                }
                out[aa * pdl + c] += acc;
            }
        }
    }
    (out, dl, pdl)
}

// Hestenes 一側 Jacobi SVD: A (rows×cols) → (U rows×k, σ, Vᵀ k×cols), σ 降順。
// U も V も回転で直交構成 (σ 除算は保持列のみ) — Gram 法の尾部破壊 (σ < 1e-8 が
// 固有値 2 乗の床で壊れる) を除去する (run6 の教訓: 尾部の質が x=16 の収束を支配)。
fn hestenes_svd(
    a_in: &[f64],
    rows: usize,
    cols: usize,
    kmax: usize,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let mut b = a_in.to_vec(); // 列ベクトル集合として扱う (row-major: b[i*cols + j])
    let mut v = vec![0.0f64; cols * cols];
    for j in 0..cols {
        v[j * cols + j] = 1.0;
    }
    let tol = 1e-14;
    for _sw in 0..40 {
        let mut rotated = false;
        for pcol in 0..cols {
            for qcol in pcol + 1..cols {
                let mut app = 0.0;
                let mut aqq = 0.0;
                let mut apq = 0.0;
                for i in 0..rows {
                    let bp = b[i * cols + pcol];
                    let bq = b[i * cols + qcol];
                    app += bp * bp;
                    aqq += bq * bq;
                    apq += bp * bq;
                }
                if apq.abs() <= tol * (app * aqq).sqrt() {
                    continue;
                }
                rotated = true;
                let theta = (aqq - app) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let c = 1.0 / (t * t + 1.0).sqrt();
                let sn = t * c;
                for i in 0..rows {
                    let bp = b[i * cols + pcol];
                    let bq = b[i * cols + qcol];
                    b[i * cols + pcol] = c * bp - sn * bq;
                    b[i * cols + qcol] = sn * bp + c * bq;
                }
                for i in 0..cols {
                    let vp = v[i * cols + pcol];
                    let vq = v[i * cols + qcol];
                    v[i * cols + pcol] = c * vp - sn * vq;
                    v[i * cols + qcol] = sn * vp + c * vq;
                }
            }
        }
        if !rotated {
            break;
        }
    }
    // σ = 列ノルム、降順ソート
    let mut sig: Vec<(f64, usize)> = (0..cols)
        .map(|j| {
            let mut n2 = 0.0;
            for i in 0..rows {
                n2 += b[i * cols + j] * b[i * cols + j];
            }
            (n2.sqrt(), j)
        })
        .collect();
    sig.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
    let k = kmax.min(cols).min(rows);
    let mut u = vec![0.0f64; rows * k];
    let mut sv = vec![0.0f64; k];
    let mut vt = vec![0.0f64; k * cols];
    for (kk, &(sg, j)) in sig.iter().take(k).enumerate() {
        sv[kk] = sg;
        if sg > 1e-300 {
            let inv = 1.0 / sg;
            for i in 0..rows {
                u[i * k + kk] = b[i * cols + j] * inv;
            }
        }
        for i in 0..cols {
            vt[kk * cols + i] = v[i * cols + j];
        }
    }
    (u, sv, vt)
}

// 大域 Rayleigh 商 ⟨ψ|H|ψ⟩/⟨ψ|ψ⟩ — 局所 Ritz 値と違い環境破損の影響を受けず、
// 非変分値が原理的に出ない (run4 の教訓: q=1 の E₁ が純ホップ下界を割る値を報告していた)
fn mps_expect_h(mps: &Mps, mpo: &Mpo) -> f64 {
    let n = mpo.n;
    let mut env = env_unit();
    for l in 0..n {
        env = upd_l(
            &env,
            &mps.a[l],
            &mpo.w[l],
            mps.bond[l],
            mps.bond[l + 1],
            mpo.dl[l],
            mpo.dr[l],
        );
    }
    let e = env.v[0];
    // ノルム (単位 MPO なし転送)
    let mut t = vec![1.0f64];
    let mut dcur = 1usize;
    for l in 0..n {
        let (dl, dr) = (mps.bond[l], mps.bond[l + 1]);
        debug_assert_eq!(dcur, dl);
        let mut nt = vec![0.0f64; dr * dr];
        for a in 0..dl {
            for c in 0..dl {
                let tv = t[a * dl + c];
                if tv == 0.0 {
                    continue;
                }
                for st in 0..2 {
                    let arow = (a * 2 + st) * dr;
                    let crow = (c * 2 + st) * dr;
                    for b in 0..dr {
                        let av = mps.a[l][arow + b];
                        if av == 0.0 {
                            continue;
                        }
                        for dcol in 0..dr {
                            nt[b * dr + dcol] += tv * av * mps.a[l][crow + dcol];
                        }
                    }
                }
            }
        }
        t = nt;
        dcur = dr;
    }
    e / t[0]
}

fn mps_overlap(x: &Mps, y: &Mps) -> f64 {
    let n = x.a.len();
    let mut t = vec![1.0f64];
    let (mut dx, mut dy) = (1usize, 1usize);
    for l in 0..n {
        let (xl, xr) = (x.bond[l], x.bond[l + 1]);
        let (yl, yr) = (y.bond[l], y.bond[l + 1]);
        debug_assert_eq!(dx, xl);
        debug_assert_eq!(dy, yl);
        let mut nt = vec![0.0f64; xr * yr];
        for a in 0..xl {
            for c in 0..yl {
                let tv = t[a * yl + c];
                if tv == 0.0 {
                    continue;
                }
                for st in 0..2 {
                    let arow = (a * 2 + st) * xr;
                    let crow = (c * 2 + st) * yr;
                    for b in 0..xr {
                        let av = x.a[l][arow + b];
                        if av == 0.0 {
                            continue;
                        }
                        for d in 0..yr {
                            nt[b * yr + d] += tv * av * y.a[l][crow + d];
                        }
                    }
                }
            }
        }
        t = nt;
        dx = xr;
        dy = yr;
    }
    // 両者とも正規化済み前提 (DMRG 出力)
    t[0]
}

// 接続電荷相関 G(r) = ⟨q_i q_j⟩ − ⟨q_i⟩⟨q_j⟩ (q = n − bg) と有効質量
// m_eff(r) = −ln(G(r+2)/G(r))/2 — E₁ 不要の質量推定器 (run8 の教訓: q=1 の
// 罰則法 E₁ は密な多ボソン梯子で降下せず、両エンジンとも未収束)
fn charge_corr(mps: &Mps, i: usize, j: usize) -> f64 {
    let n = mps.a.len();
    let bg = |k: usize| if k % 2 == 1 { 1.0 } else { 0.0 };
    // ⟨O_i O_j⟩ / ⟨1⟩: O = n − bg (対角)
    let contract = |ops: &[(usize, bool)]| -> f64 {
        // ops: (site, insert_q)
        let mut t = vec![1.0f64];
        let mut dcur = 1usize;
        for l in 0..n {
            let (dl, dr) = (mps.bond[l], mps.bond[l + 1]);
            debug_assert_eq!(dcur, dl);
            let ins = ops.iter().any(|&(k, on)| on && k == l);
            let mut nt = vec![0.0f64; dr * dr];
            for a in 0..dl {
                for c in 0..dl {
                    let tv = t[a * dl + c];
                    if tv == 0.0 {
                        continue;
                    }
                    for st in 0..2 {
                        let w = if ins { st as f64 - bg(l) } else { 1.0 };
                        if w == 0.0 {
                            continue;
                        }
                        let arow = (a * 2 + st) * dr;
                        let crow = (c * 2 + st) * dr;
                        for b in 0..dr {
                            let av = mps.a[l][arow + b];
                            if av == 0.0 {
                                continue;
                            }
                            for d2 in 0..dr {
                                nt[b * dr + d2] += tv * w * av * mps.a[l][crow + d2];
                            }
                        }
                    }
                }
            }
            t = nt;
            dcur = dr;
        }
        t[0]
    };
    let nrm = contract(&[]);
    let qq = contract(&[(i, true), (j, true)]) / nrm;
    let qi = contract(&[(i, true)]) / nrm;
    let qj = contract(&[(j, true)]) / nrm;
    qq - qi * qj
}

#[allow(clippy::too_many_arguments)]
fn dmrg(
    mpo: &Mpo,
    chi: usize,
    sweeps: usize,
    seed: u64,
    ortho: Option<&Mps>,
    lam_o: f64,
    warm: Option<&Mps>,
    rand_init: bool,
) -> (f64, Mps) {
    let n = mpo.n;
    let mut mps = match warm {
        Some(m) => m.clone(),
        None => init_mps(n, chi, seed, rand_init),
    };
    // 環境
    let mut lc: Vec<Env> = (0..=n).map(|_| env_unit()).collect();
    let mut rc: Vec<Env> = (0..=n).map(|_| env_unit()).collect();
    for l in (1..n).rev() {
        rc[l] = upd_r(
            &rc[l + 1],
            &mps.a[l],
            &mpo.w[l],
            mps.bond[l],
            mps.bond[l + 1],
            mpo.dl[l],
            mpo.dr[l],
        );
    }
    // 直交罰則環境
    let mut oe = OrthoEnv {
        ol: vec![vec![1.0]; n + 1],
        or_: vec![vec![1.0]; n + 1],
        dims_ol: vec![(1, 1); n + 1],
        dims_or: vec![(1, 1); n + 1],
    };
    if let Some(phi) = ortho {
        for l in (1..n).rev() {
            let (v, da, dc) = ov_r(
                &oe.or_[l + 1],
                oe.dims_or[l + 1].0,
                oe.dims_or[l + 1].1,
                &mps.a[l],
                mps.bond[l],
                mps.bond[l + 1],
                &phi.a[l],
                phi.bond[l],
                phi.bond[l + 1],
            );
            oe.or_[l] = v;
            oe.dims_or[l] = (da, dc);
        }
    }
    let mut energy = 0.0f64;
    for _sw in 0..sweeps {
        // λ_o ランプ: 前半 sweep は弱め (脱出容易)、後半で全強度 (直交強制)
        let lam_eff = if _sw * 2 < sweeps { lam_o * 0.2 } else { lam_o };
        for dir in 0..2 {
            let bonds: Vec<usize> = if dir == 0 {
                (0..n - 1).collect()
            } else {
                (0..n - 1).rev().collect()
            };
            for &l in &bonds {
                let (d0, d1, d2) = (mps.bond[l], mps.bond[l + 1], mps.bond[l + 2]);
                let dim = d0 * 2 * 2 * d2;
                // θ0[a,s,t,b] = Σ_m A[l][(a,s),m]·A[l+1][(m,t),b]
                let mut th = vec![0.0f64; dim];
                for aa in 0..d0 {
                    for s in 0..2 {
                        let arow = (aa * 2 + s) * d1;
                        for m in 0..d1 {
                            let av = mps.a[l][arow + m];
                            if av == 0.0 {
                                continue;
                            }
                            for t in 0..2 {
                                let brow = (m * 2 + t) * d2;
                                let orow = ((aa * 2 + s) * 2 + t) * d2;
                                for b in 0..d2 {
                                    th[orow + b] += av * mps.a[l + 1][brow + b];
                                }
                            }
                        }
                    }
                }
                // 直交射影ベクトル
                let mut pvecs: Vec<Vec<f64>> = Vec::new();
                if let Some(phi) = ortho {
                    let (pa, pcp) = oe.dims_ol[l];
                    let (pb, pep) = oe.dims_or[l + 2];
                    let (p0, p1, p2) = (phi.bond[l], phi.bond[l + 1], phi.bond[l + 2]);
                    debug_assert_eq!(pcp, p0);
                    debug_assert_eq!(pep, p2);
                    // pv[a,s,t,b] = Σ oL[a,c]·φ[l][(c,s),d]·φ[l+1][(d,t),e]·oR[b,e]
                    // T1[a,s,d] = Σ_c oL[a,c]·φ[(c,s),d]
                    let mut t1 = vec![0.0f64; pa * 2 * p1];
                    for aa in 0..pa {
                        for c in 0..p0 {
                            let ov = oe.ol[l][aa * p0 + c];
                            if ov == 0.0 {
                                continue;
                            }
                            for s in 0..2 {
                                let prow = (c * 2 + s) * p1;
                                let trow = (aa * 2 + s) * p1;
                                for d in 0..p1 {
                                    t1[trow + d] += ov * phi.a[l][prow + d];
                                }
                            }
                        }
                    }
                    // T2[a,s,t,e] = Σ_d T1[a,s,d]·φ[l+1][(d,t),e]
                    let mut t2 = vec![0.0f64; pa * 2 * 2 * p2];
                    for aa in 0..pa {
                        for s in 0..2 {
                            let trow = (aa * 2 + s) * p1;
                            for d in 0..p1 {
                                let tv = t1[trow + d];
                                if tv == 0.0 {
                                    continue;
                                }
                                for t in 0..2 {
                                    let prow = (d * 2 + t) * p2;
                                    let orow = ((aa * 2 + s) * 2 + t) * p2;
                                    for e in 0..p2 {
                                        t2[orow + e] += tv * phi.a[l + 1][prow + e];
                                    }
                                }
                            }
                        }
                    }
                    // pv[a,s,t,b] = Σ_e T2[a,s,t,e]·oR[b,e]
                    let mut pv = vec![0.0f64; dim];
                    for aa in 0..pa {
                        for s in 0..2 {
                            for t in 0..2 {
                                let trow = ((aa * 2 + s) * 2 + t) * p2;
                                let orow = ((aa * 2 + s) * 2 + t) * d2;
                                for b in 0..pb {
                                    let mut acc = 0.0;
                                    for e in 0..p2 {
                                        acc += t2[trow + e] * oe.or_[l + 2][b * p2 + e];
                                    }
                                    pv[orow + b] += acc;
                                }
                            }
                        }
                    }
                    let nr: f64 = pv.iter().map(|z| z * z).sum::<f64>().sqrt();
                    if nr > 1e-12 {
                        for z in pv.iter_mut() {
                            *z /= nr;
                        }
                        pvecs.push(pv);
                    }
                }
                // hmul クロージャ
                let le = &lc[l];
                let re = &rc[l + 2];
                let w1 = &mpo.w[l];
                let w2 = &mpo.w[l + 1];
                let (wd0, wd1) = (mpo.dl[l], mpo.dr[l]);
                let (wd1b, wd2) = (mpo.dl[l + 1], mpo.dr[l + 1]);
                debug_assert_eq!(wd1, wd1b);
                let hmul = |v: &[f64], out: &mut Vec<f64>| {
                    out.clear();
                    out.resize(dim, 0.0);
                    // t1[p,c,u,v,b] = Σ_a Le[a,p,c]·θ[(a,u,v),b]
                    let mut t1 = vec![0.0f64; wd0 * d0 * 4 * d2];
                    for aa in 0..d0 {
                        for p in 0..wd0 {
                            for c in 0..d0 {
                                let lv = le.v[(aa * le.dp + p) * le.dc + c];
                                if lv == 0.0 {
                                    continue;
                                }
                                for uv in 0..4 {
                                    let vrow = (aa * 4 + uv) * d2;
                                    let trow = ((p * d0 + c) * 4 + uv) * d2;
                                    for b in 0..d2 {
                                        t1[trow + b] += lv * v[vrow + b];
                                    }
                                }
                            }
                        }
                    }
                    // t2[q,c,s,v,b] = Σ_{p,u} t1[p,c,(u,v),b]·W1[(p,q),(u,s)]
                    let mut t2 = vec![0.0f64; wd1 * d0 * 4 * d2];
                    for p in 0..wd0 {
                        for q in 0..wd1 {
                            for u in 0..2 {
                                for s in 0..2 {
                                    let wv = w1[((p * wd1 + q) * 2 + u) * 2 + s];
                                    if wv == 0.0 {
                                        continue;
                                    }
                                    for c in 0..d0 {
                                        for vv in 0..2 {
                                            let trow = ((p * d0 + c) * 4 + u * 2 + vv) * d2;
                                            let orow = ((q * d0 + c) * 4 + s * 2 + vv) * d2;
                                            for b in 0..d2 {
                                                t2[orow + b] += wv * t1[trow + b];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // t3[r,c,s,t,b] = Σ_{q,v} t2[q,c,(s,v),b]·W2[(q,r),(v,t)]
                    let mut t3 = vec![0.0f64; wd2 * d0 * 4 * d2];
                    for q in 0..wd1 {
                        for r in 0..wd2 {
                            for vv in 0..2 {
                                for t in 0..2 {
                                    let wv = w2[((q * wd2 + r) * 2 + vv) * 2 + t];
                                    if wv == 0.0 {
                                        continue;
                                    }
                                    for c in 0..d0 {
                                        for s in 0..2 {
                                            let trow = ((q * d0 + c) * 4 + s * 2 + vv) * d2;
                                            let orow = ((r * d0 + c) * 4 + s * 2 + t) * d2;
                                            for b in 0..d2 {
                                                t3[orow + b] += wv * t2[trow + b];
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                    // out[c,s,t,d] = Σ_{r,b} t3[r,c,(s,t),b]·Re[(b,r),d]
                    for r in 0..wd2 {
                        for c in 0..d0 {
                            for st in 0..4 {
                                let trow = ((r * d0 + c) * 4 + st) * d2;
                                let orow = (c * 4 + st) * d2;
                                for b in 0..d2 {
                                    let tv = t3[trow + b];
                                    if tv == 0.0 {
                                        continue;
                                    }
                                    let rrow = (b * re.dp + r) * re.dc;
                                    for d in 0..d2 {
                                        out[orow + d] += tv * re.v[rrow + d];
                                    }
                                }
                            }
                        }
                    }
                    // 直交罰則
                    for pv in &pvecs {
                        let mut ip = 0.0;
                        for i in 0..dim {
                            ip += pv[i] * v[i];
                        }
                        let lp = lam_eff * ip;
                        for i in 0..dim {
                            out[i] += lp * pv[i];
                        }
                    }
                };
                // Lanczos (m=25, 全再直交)
                let mut v0 = th.clone();
                let nr: f64 = v0.iter().map(|z| z * z).sum::<f64>().sqrt();
                for z in v0.iter_mut() {
                    *z /= nr;
                }
                let m_kry = 25usize;
                let mut basis: Vec<Vec<f64>> = vec![v0];
                let mut al: Vec<f64> = Vec::new();
                let mut be: Vec<f64> = Vec::new();
                let mut wbuf: Vec<f64> = Vec::new();
                for j in 0..m_kry {
                    hmul(&basis[j], &mut wbuf);
                    let a: f64 = basis[j].iter().zip(wbuf.iter()).map(|(x, y)| x * y).sum();
                    al.push(a);
                    // w −= a v_j + b v_{j−1}, 全再直交
                    for (bi, bv) in basis.iter().enumerate() {
                        let ip: f64 = bv.iter().zip(wbuf.iter()).map(|(x, y)| x * y).sum();
                        if ip != 0.0 {
                            for i in 0..dim {
                                wbuf[i] -= ip * bv[i];
                            }
                        }
                        let _ = bi;
                    }
                    let b: f64 = wbuf.iter().map(|z| z * z).sum::<f64>().sqrt();
                    if b < 1e-11 || j == m_kry - 1 {
                        break;
                    }
                    be.push(b);
                    basis.push(wbuf.iter().map(|z| z / b).collect());
                }
                let k = al.len();
                let mut tmat = vec![0.0f64; k * k];
                for i in 0..k {
                    tmat[i * k + i] = al[i];
                    if i + 1 < k {
                        tmat[i * k + (i + 1)] = be[i];
                        tmat[(i + 1) * k + i] = be[i];
                    }
                }
                let (evs, evec) = jacobi_eigh(&tmat, k);
                energy = evs[0];
                let mut gv = vec![0.0f64; dim];
                for aidx in 0..k {
                    let cc = evec[aidx]; // 最低固有ベクトルの成分 a
                    if cc == 0.0 {
                        continue;
                    }
                    for i in 0..dim {
                        gv[i] += cc * basis[aidx][i];
                    }
                }
                let nr: f64 = gv.iter().map(|z| z * z).sum::<f64>().sqrt();
                for z in gv.iter_mut() {
                    *z /= nr;
                }
                // 分解: Hestenes SVD (U も Vᵀ も回転直交 — 尾部の質を保つ)
                let rows = d0 * 2;
                let colsn = 2 * d2;
                // ランク下限 1・σ > 1e-10 で切る (python 実装と同一基準)
                let (umat, sv, vt) = hestenes_svd(&gv, rows, colsn, chi);
                let mut kkeep = 0usize;
                for kk2 in 0..sv.len() {
                    if sv[kk2] > 1e-10 {
                        kkeep += 1;
                    } else {
                        break;
                    }
                }
                let kkeep = kkeep.max(1);
                let snorm: f64 = sv[..kkeep].iter().map(|z| z * z).sum::<f64>().sqrt();
                if dir == 0 {
                    // A[l] = U (左正準), A[l+1] = diag(S/‖S‖)Vᵀ
                    let kfull = sv.len();
                    let mut a0 = vec![0.0f64; d0 * 2 * kkeep];
                    for i in 0..rows {
                        for kk2 in 0..kkeep {
                            a0[i * kkeep + kk2] = umat[i * kfull + kk2];
                        }
                    }
                    let mut a1 = vec![0.0f64; kkeep * 2 * d2];
                    for kk2 in 0..kkeep {
                        let w = sv[kk2] / snorm;
                        for j in 0..colsn {
                            a1[kk2 * colsn + j] = w * vt[kk2 * colsn + j];
                        }
                    }
                    mps.a[l] = a0;
                    mps.a[l + 1] = a1;
                    mps.bond[l + 1] = kkeep;
                    lc[l + 1] = upd_l(
                        &lc[l],
                        &mps.a[l],
                        &mpo.w[l],
                        mps.bond[l],
                        mps.bond[l + 1],
                        mpo.dl[l],
                        mpo.dr[l],
                    );
                    if let Some(phi) = ortho {
                        let (v2, da, dc) = ov_l(
                            &oe.ol[l],
                            oe.dims_ol[l].0,
                            oe.dims_ol[l].1,
                            &mps.a[l],
                            mps.bond[l],
                            mps.bond[l + 1],
                            &phi.a[l],
                            phi.bond[l],
                            phi.bond[l + 1],
                        );
                        oe.ol[l + 1] = v2;
                        oe.dims_ol[l + 1] = (da, dc);
                    }
                } else {
                    // A[l+1] = Vᵀ (右正準), A[l] = U diag(S/‖S‖)
                    let kfull = sv.len();
                    let mut a1 = vec![0.0f64; kkeep * 2 * d2];
                    for kk2 in 0..kkeep {
                        for j in 0..colsn {
                            a1[kk2 * colsn + j] = vt[kk2 * colsn + j];
                        }
                    }
                    let mut a0 = vec![0.0f64; d0 * 2 * kkeep];
                    for i in 0..rows {
                        for kk2 in 0..kkeep {
                            a0[i * kkeep + kk2] = umat[i * kfull + kk2] * sv[kk2] / snorm;
                        }
                    }
                    mps.a[l] = a0;
                    mps.a[l + 1] = a1;
                    mps.bond[l + 1] = kkeep;
                    rc[l + 1] = upd_r(
                        &rc[l + 2],
                        &mps.a[l + 1],
                        &mpo.w[l + 1],
                        mps.bond[l + 1],
                        mps.bond[l + 2],
                        mpo.dl[l + 1],
                        mpo.dr[l + 1],
                    );
                    if let Some(phi) = ortho {
                        let (v, da, dc) = ov_r(
                            &oe.or_[l + 2],
                            oe.dims_or[l + 2].0,
                            oe.dims_or[l + 2].1,
                            &mps.a[l + 1],
                            mps.bond[l + 1],
                            mps.bond[l + 2],
                            &phi.a[l + 1],
                            phi.bond[l + 1],
                            phi.bond[l + 2],
                        );
                        oe.or_[l + 1] = v;
                        oe.dims_or[l + 1] = (da, dc);
                    }
                }
            }
        }
    }
    let _ = energy;
    (mps_expect_h(&mps, mpo), mps)
}

fn gap_pair(n: usize, x: f64, qf: f64, chi: usize, sweeps: usize) -> (f64, f64, f64, Mps) {
    let mpo = build_mpo(n, x, qf, LAMQ);
    // x アニーリング (n > 10): 小 x では Néel がほぼ厳密 — 断熱継続で強結合の
    // 捕獲状態を回避 (run3 の教訓: x9_q1 は初期化非依存の吸着点 −345.685 を持つ)
    let (e0, mps0) = if n > 10 {
        let mut mps: Option<Mps> = None;
        let mut e = 0.0;
        for &xs in &[x / 9.0, x / 3.0, x] {
            let mpo_s = build_mpo(n, xs, qf, LAMQ);
            let sw = if (xs - x).abs() < 1e-12 { sweeps } else { 4 };
            let (ee, m2) = dmrg(&mpo_s, chi, sw, 7, None, 0.0, mps.as_ref(), false);
            e = ee;
            mps = Some(m2);
        }
        (e, mps.unwrap())
    } else {
        dmrg(&mpo, chi, sweeps, 7, None, 0.0, None, false)
    };
    // E₁: Néel 初期化 + λ ランプ (run7 の教訓: 乱数初期化は q=1 の N=64 で
    // 遠方盆地に停滞する [gap 55-102 の変分上界どまり]。Néel は ψ₀ 方向に近く
    // 罰則が早期から直交低状態へ押し出す。ψ₀ 厳密一致でないので Krylov 崩壊もない)
    let (e1, mps1) = dmrg(&mpo, chi, sweeps, 11, Some(&mps0), LAM_O, None, false);
    let ov = mps_overlap(&mps1, &mps0).abs();
    (e0, e1, ov, mps0)
}

// 全空間 ED (自己ゲート用, N ≤ 10)
fn ed_full(n: usize, x: f64, qf: f64, qpen: f64) -> (f64, f64) {
    let dim = 1usize << n;
    let mut h = vec![0.0f64; dim * dim];
    for s in 0..dim {
        let mut e = 0.0;
        let mut acc = 0.0;
        for k in 0..n {
            let occ = ((s >> k) & 1) as f64;
            let bg = if k % 2 == 1 { 1.0 } else { 0.0 };
            acc += qf * (occ - bg);
            if k < n - 1 {
                e += acc * acc;
            }
        }
        e += qpen * acc * acc;
        h[s * dim + s] += e;
        for k in 0..n - 1 {
            let b0 = (s >> k) & 1;
            let b1 = (s >> (k + 1)) & 1;
            if b0 != b1 {
                let t = s ^ (1 << k) ^ (1 << (k + 1));
                h[s * dim + t] += -x;
            }
        }
    }
    let (ev, _) = jacobi_eigh(&h, dim);
    (ev[0], ev[1])
}

fn main() {
    self_test();
    println!("=== v22.1 計算層: 重電荷 Schwinger DMRG (Rust, 決定的 6 スレッド) ===\n");
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!(
            "  [{}] {}  {}",
            if ok { "PASS" } else { "FAIL" },
            name,
            detail
        );
        if !ok {
            nfail += 1;
        }
    };
    let t0 = std::time::Instant::now();

    // ---- 自己ゲート [G1/G2]: N=8 で DMRG = ED ----
    {
        let (ed0, ed1) = ed_full(8, 2.0, 1.0, LAMQ);
        let (e0, e1, _ov8, _m8) = gap_pair(8, 2.0, 1.0, 32, 8);
        check(
            "[G1] N=8 DMRG E₀ = ED ± 1e-8",
            (e0 - ed0).abs() < 1e-8,
            format!(
                "DMRG {:.10} vs ED {:.10} (Δ {:.1e})",
                e0,
                ed0,
                (e0 - ed0).abs()
            ),
        );
        check(
            "[G2] N=8 罰則法 E₁ = ED ± 1e-6",
            (e1 - ed1).abs() < 1e-6,
            format!(
                "DMRG {:.10} vs ED {:.10} (Δ {:.1e}, {} s)",
                e1,
                ed1,
                (e1 - ed1).abs(),
                t0.elapsed().as_secs()
            ),
        );
    }

    // ---- 6 ジョブを並列実行 (各ジョブは内部直列・固定シード = 決定的) ----
    let jobs: Vec<(&str, usize, f64, f64, usize, usize)> = vec![
        ("anchor_q1", 10, 2.0, 1.0, 64, 8),
        ("anchor_q2", 10, 2.0, 2.0, 64, 8),
        ("n64_x9_q1", 64, 9.0, 1.0, 128, 12),
        ("n64_x9_q2", 64, 9.0, 2.0, 128, 12),
        ("n64_x16_q1", 64, 16.0, 1.0, 128, 12),
        ("n64_x16_q2", 64, 16.0, 2.0, 128, 12),
    ];
    let results: Vec<(String, f64, f64, f64, (f64, f64, f64))> = std::thread::scope(|sc| {
        let handles: Vec<_> = jobs
            .iter()
            .map(|&(name, n, x, qf, chi, sw)| {
                sc.spawn(move || {
                    let (e0, e1, ov, mps0) = gap_pair(n, x, qf, chi, sw);
                    // 相関質量: 主推定器 = 深プラトー m_eff @ r=20 (run9 の教訓:
                    // r=8..16 は q=1 が下から・q=2 が上から収束中で漸近が浅く、
                    // 比が +8% 偏る。基準 2±7% は不変 — 漸近域の精細化)。
                    // 記録: r=12 (run9 の主) と r=16。
                    let meff = if n >= 64 {
                        let g = |r: usize| charge_corr(&mps0, n / 2 - r / 2, n / 2 + r / 2).abs();
                        let m20 = -((g(22) / g(18)).ln()) / 4.0;
                        let m12 = -((g(14) / g(10)).ln()) / 4.0;
                        let m16 = -((g(18) / g(14)).ln()) / 4.0;
                        (m20, m12, m16)
                    } else {
                        (0.0, 0.0, 0.0)
                    };
                    (name.to_string(), e0, e1, ov, meff)
                })
            })
            .collect();
        handles.into_iter().map(|h| h.join().unwrap()).collect()
    });
    for (name, e0, e1, ov, meff) in &results {
        println!(
            "    [{}] E₀ = {:.8}, E₁ = {:.8}, gap = {:.6}, |⟨ψ₁|ψ₀⟩| = {:.1e}, m_eff = {:.4} [{:.4}, {:.4}] ({} s)",
            name,
            e0,
            e1,
            e1 - e0,
            ov,
            meff.0,
            meff.1,
            meff.2,
            t0.elapsed().as_secs()
        );
    }

    // ---- 自己ゲート [G3/G4]: アンカー N=10 = ED (dim 1024) ----
    for (qi, qf) in [(0usize, 1.0f64), (1, 2.0)] {
        let (ed0, ed1) = ed_full(10, 2.0, qf, LAMQ);
        let (_, e0, e1, _ov, _me) = &results[qi];
        check(
            &format!("[G3] アンカー q={} E₀ = ED ± 1e-7", qf as i64),
            (e0 - ed0).abs() < 1e-7,
            format!("Δ {:.1e}", (e0 - ed0).abs()),
        );
        check(
            &format!("[G4] アンカー q={} E₁ = ED ± 1e-5 (罰則法)", qf as i64),
            (e1 - ed1).abs() < 1e-5,
            format!("Δ {:.1e} (ED gap {:.6})", (e1 - ed1).abs(), ed1 - ed0),
        );
    }

    // ---- [G5] 交差エンジン錨: python 実装 (explore/dmrg_heavy.py, run2 完走分) の
    // E₀ と一致するか — 独立実装・独立初期化の変分アンカー
    {
        let py = [
            ("n64_x9_q2", -305.45320375377077f64),
            ("n64_x16_q2", -582.2989782722372),
        ];
        for (name, pye0) in py {
            let e0 = results.iter().find(|r| r.0 == name).unwrap().1;
            check(
                &format!("[G5] {} E₀ = python ± 0.02 (交差エンジン)", name),
                (e0 - pye0).abs() < 0.02,
                format!(
                    "rust {:.6} vs python {:.6} (Δ {:.1e})",
                    e0,
                    pye0,
                    (e0 - pye0).abs()
                ),
            );
        }
    }
    // ---- [G6] 対整合: E₁ ≥ E₀ − 1e-6 (罰則法が下を見つけたら E₀ が stuck の自己検出)
    for (name, e0, e1, ov, _me) in &results {
        if name.starts_with("n64") {
            check(
                &format!("[G6] {} 対整合 E₁ ≥ E₀", name),
                e1 >= &(e0 - 1e-6),
                format!("gap = {:.6}", e1 - e0),
            );
            check(
                &format!("[G7] {} 直交度 |⟨ψ₁|ψ₀⟩| < 0.1", name),
                *ov < 0.1,
                format!("{:.1e}", ov),
            );
        }
    }

    // ---- [G8] 相関質量 vs gap 質量 (q=2 — 両推定器が有効な場での内部整合) ----
    for (name, x) in [("n64_x9_q2", 9.0f64), ("n64_x16_q2", 16.0)] {
        let r = results.iter().find(|r| r.0 == name).unwrap();
        let m_gap = (r.2 - r.1) / (2.0 * x.sqrt()); // M/g (gap 法)
        let m_corr = (r.4).0 * x.sqrt(); // m_eff は格子単位 → M/g = m_eff·√x
        check(
            &format!("[G8] {} 相関質量 = gap 質量 ± 20%", name),
            (m_corr / m_gap - 1.0).abs() < 0.20,
            format!(
                "corr {:.4} vs gap {:.4} (比 {:.3})",
                m_corr,
                m_gap,
                m_corr / m_gap
            ),
        );
    }

    // ---- JSON 出力 (v221_dmrgex の一次データ) ----
    let mut kv: Vec<(String, f64)> = Vec::new();
    for (name, e0, e1, _ov, _me) in &results {
        if name.starts_with("anchor") {
            let q = &name[7..];
            kv.push((format!("anchor_n10_{}_e0", q), *e0));
            kv.push((format!("anchor_n10_{}_e1", q), *e1));
        } else {
            kv.push((format!("{}_e0", name), *e0));
            kv.push((format!("{}_e1", name), *e1));
            kv.push((format!("{}_gap", name), e1 - e0));
            kv.push((
                format!("{}_meff", name),
                _ov * 0.0 + results.iter().find(|r| &r.0 == name).unwrap().4 .0,
            ));
        }
    }
    let mut js = String::from("{\n");
    for (i, (k, v)) in kv.iter().enumerate() {
        js.push_str(&format!(
            " \"{}\": {:.12}{}\n",
            k,
            v,
            if i + 1 < kv.len() { "," } else { "" }
        ));
    }
    js.push('}');
    let mut f = std::fs::File::create("results/v221_dmrg_gaps.json").expect("write json");
    f.write_all(js.as_bytes()).expect("write");
    println!("\n[artifact] results/v221_dmrg_gaps.json");
    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
