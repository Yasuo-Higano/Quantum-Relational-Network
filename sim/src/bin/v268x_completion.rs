//! v26.8-X v268x_completion — X チャネル 2 比の連続極限 (PRED-016 の 4 比完成)
//!
//! 事前登録: spec §12.2-4/§12.4-5 (bc644d4)。D チャネル 2 比 (staggered 2×/Wilson 1×,
//! v26.8-B/C) に X (cross 偏極) の 2 比を加え、PRED-016 の 4 比
//!   A_D^stag/2A_or, A_X^stag/(Z²·2A_or), A_D^Wil/A_or, A_X^Wil/(A_or/2)
//! を完成させる。**正直な線引き**: PRED-016 の登録バーは「1% 以内・系統 0.5%」—
//! 現器械の外挿 spread は ~2% なので、本版は「4/4 比が 1 ± 2% で成立」までを
//! 確立し、登録バー (1%) への到達は残作業として interim 記録する (インフレ禁止)。
//!
//! X 頂点:
//!  staggered: v268z で tree-level 認証済みの 4 隅 point-split (両片 w = −σρ/16,
//!    **全 k で厳密 taste-singlet** — Clifford 像の元) をそのまま BZ 積分へ。
//!    node 極限 = 2·T_xz (凍結構成の tree 正規化 Z = 2) ⇒ raw の的 = Z²·A[T_xz]
//!    = 4·A[T_xz]。連続 2-taste の A[T_xz] = 2A_oracle/2 = A_oracle (σ_Txz = σ_D/2,
//!    v26.8-A で σ_X̂ = σ_D を機械証明済み) ⇒ **的 = 4·A_oracle(1 Dirac)**。
//!  Wilson: 混合ボンド転写 (η 不要 — taste なし): z ホップに −iαx/4 + x ホップに
//!    −iαz/4 (+h.c.) ⇒ V_X^Wil(k) = ½(αx sin kz + αz sin kx) — node で
//!    ½(αx pz + αz px) = T_xz^cont (Z = 1)。h_xz への Wilson-β 片の結合は 0 に取る
//!    (O(a) の scheme 自由度 — 文書化)。**的 = A[T_xz] 1 flavor = A_oracle/2**。
//!
//! 検査 (凍結):
//!  [S0] staggered X 頂点の taste-singlet 性の BZ 抜き取り: 遠ノード点でも
//!       V ∈ Clifford 像 (8 成分基底の (−1)^{sx}/(−1)^{sz} 対角×flip 構造 —
//!       v268z の構成的証明の回帰。数値: 行列が Γ 構造の張る形に一致 1e-12)
//!  [S1] 求積自己整合 (両エンジン, a_min, GL up): < 1e-3
//!  [S2] **staggered X**: A_raw(a) ladder の凍結 2 モデル外挿 (a²/a+a²) で
//!       A₀/(4A_oracle) = 1 ± 0.02
//!  [S3] **Wilson X**: 凍結 2 モデル (a+a²/a²+a⁴) で A₀/(A_oracle/2) = 1 ± 0.02
//!  [S4] PRED-016 総括: 4 比全て 1 ± 2% (D 2 比は v26.8-B/C の凍結値を引用) —
//!       登録バー 1% 未達の interim 判定を明記
//!  [S5] 変異: staggered X piece2 の η parity 破り (ε 3→2, v268z S8 と同一) →
//!       A が > 10% 変化
//!
//! 事前登録分岐: (a) S0–S3 PASS → 4/4 比の universality が 2% で成立 (PRED-016 は
//!   interim — 1% バーは残作業) / (b) S2/S3 miss → X 転写の one-loop 破綻 (税 =
//!   tree 認証だけでは不足だった事実を公表) / (c) S0/S1 FAIL → 器械。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

fn a_oracle_1d() -> f64 {
    -1.0 / (160.0 * PI * PI)
}

// ================= staggered 8 成分エンジン (v268b の写経) =================

fn sbit(s: usize, ax: usize) -> usize {
    (s >> ax) & 1
}

fn h8(k: [f64; 3], m: f64) -> Vec<f64> {
    let mut h = vec![0.0f64; 64];
    for s in 0..8usize {
        let cx = if sbit(s, 0) == 0 { 1.0 } else { -1.0 } * k[0].cos();
        h[s + s * 8] += cx;
        let s2 = s ^ 1;
        let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * k[1].cos();
        h[s2 + s * 8] += cy;
        let s3 = s ^ 3;
        let cz = if sbit(s, 2) == 0 { 1.0 } else { -1.0 } * k[2].cos();
        h[s3 + s * 8] += cz;
        let s4 = s ^ 7;
        h[s4 + s * 8] += m;
    }
    h
}

struct Term {
    eps: usize,
    d: [i32; 3],
    w: f64,
}

fn vertex8(terms: &[Term], k: [f64; 3], q: f64) -> Vec<f64> {
    let mut v = vec![0.0f64; 64];
    for t in terms {
        for s in 0..8usize {
            let s2 = s ^ t.eps;
            let mut ph1 = 0.0f64;
            let mut ph2 = 0.0f64;
            for ax in 0..3 {
                let ka = k[ax] + PI * sbit(s, ax) as f64;
                let ka2 = k[ax] + if ax == 1 { q } else { 0.0 } + PI * sbit(s2, ax) as f64;
                ph1 += ka * t.d[ax] as f64;
                ph2 += ka2 * t.d[ax] as f64;
            }
            let mid = 0.5 * q * t.d[1] as f64;
            v[s + s2 * 8] += t.w * ((ph1 + mid).cos() + (-ph2 + mid).cos());
        }
    }
    v
}

/// v268z で認証済みの point-split X (両片 w = −σρ/16)。mutate: piece2 ε 3→2
fn t_x_split(mutate: bool) -> Vec<Term> {
    let mut v = Vec::new();
    for sg in [1i32, -1] {
        for rh in [1i32, -1] {
            let c = (sg * rh) as f64 / 16.0;
            v.push(Term { eps: 0, d: [sg, 0, 2 * rh], w: -c });
            v.push(Term {
                eps: if mutate { 2 } else { 3 },
                d: [2 * rh, 0, sg],
                w: -c,
            });
        }
    }
    v
}

fn chi_x_integrand(k: [f64; 3], qs_lat: &[f64], w_null: &[f64], m: f64, mutate: bool) -> f64 {
    let hk = h8(k, m);
    let (wk, vk) = jacobi_eigh(&hk, 8);
    let mut acc = 0.0f64;
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let hq = h8(kq, m);
        let (wq, vq) = jacobi_eigh(&hq, 8);
        let vx = vertex8(&t_x_split(mutate), k, q);
        let mut chi = 0.0f64;
        for mu in 4..8 {
            for nu in 0..4 {
                let mut mel = 0.0f64;
                for r in 0..8 {
                    let mut s = 0.0f64;
                    for c in 0..8 {
                        s += vx[c + r * 8] * vk[c + nu * 8];
                    }
                    mel += vq[r + mu * 8] * s;
                }
                chi += 2.0 * mel * mel / (wq[mu] - wk[nu]);
            }
        }
        acc += w_null[qi] * chi;
    }
    acc
}

// ================= Wilson 4 成分エンジン (v268c の写経) =================

type M4 = [(f64, f64); 16];

fn mzero() -> M4 {
    [(0.0, 0.0); 16]
}

fn mmul(a: &M4, b: &M4) -> M4 {
    let mut o = mzero();
    for r in 0..4 {
        for k in 0..4 {
            let av = a[k + r * 4];
            if av.0 == 0.0 && av.1 == 0.0 {
                continue;
            }
            for c in 0..4 {
                let bv = b[c + k * 4];
                o[c + r * 4].0 += av.0 * bv.0 - av.1 * bv.1;
                o[c + r * 4].1 += av.0 * bv.1 + av.1 * bv.0;
            }
        }
    }
    o
}

fn mtrace_re(a: &M4) -> f64 {
    (0..4).map(|i| a[i + i * 4].0).sum()
}

fn alphas() -> [M4; 4] {
    let mut ax = mzero();
    let mut ay = mzero();
    let mut az = mzero();
    let mut b = mzero();
    let put = |m: &mut M4, r: usize, c: usize, re: f64, im: f64| {
        m[c + r * 4] = (re, im);
    };
    for blk in 0..2 {
        let (ro, co) = if blk == 0 { (0, 2) } else { (2, 0) };
        put(&mut ax, ro, co + 1, 1.0, 0.0);
        put(&mut ax, ro + 1, co, 1.0, 0.0);
        put(&mut ay, ro, co + 1, 0.0, -1.0);
        put(&mut ay, ro + 1, co, 0.0, 1.0);
        put(&mut az, ro, co, 1.0, 0.0);
        put(&mut az, ro + 1, co + 1, -1.0, 0.0);
    }
    for i in 0..4 {
        put(&mut b, i, i, if i < 2 { 1.0 } else { -1.0 }, 0.0);
    }
    [ax, ay, az, b]
}

fn mwil(k: [f64; 3], m: f64, r: f64) -> f64 {
    m + r * (3.0 - k[0].cos() - k[1].cos() - k[2].cos())
}

fn ewil(k: [f64; 3], m: f64, r: f64) -> f64 {
    let s2 = k[0].sin().powi(2) + k[1].sin().powi(2) + k[2].sin().powi(2);
    let mk = mwil(k, m, r);
    (s2 + mk * mk).sqrt()
}

fn projw(al: &[M4; 4], k: [f64; 3], m: f64, r: f64, sign: f64) -> M4 {
    let mut h = mzero();
    for ax in 0..3 {
        let s = k[ax].sin();
        for i in 0..16 {
            h[i].0 += s * al[ax][i].0;
            h[i].1 += s * al[ax][i].1;
        }
    }
    let mk = mwil(k, m, r);
    for i in 0..16 {
        h[i].0 += mk * al[3][i].0;
        h[i].1 += mk * al[3][i].1;
    }
    let e = ewil(k, m, r);
    let mut p = mzero();
    for i in 0..16 {
        p[i].0 = sign * h[i].0 / (2.0 * e);
        p[i].1 = sign * h[i].1 / (2.0 * e);
    }
    for i in 0..4 {
        p[i + i * 4].0 += 0.5;
    }
    p
}

/// Wilson X 頂点: V = ½(αx sin kz + αz sin kx) (h_xz への Wilson-β 結合は 0 — scheme)
fn vertex_x_wil(al: &[M4; 4], k: [f64; 3]) -> M4 {
    let mut v = mzero();
    let (sz, sx) = (k[2].sin(), k[0].sin());
    for i in 0..16 {
        v[i].0 = 0.5 * (sz * al[0][i].0 + sx * al[2][i].0);
        v[i].1 = 0.5 * (sz * al[0][i].1 + sx * al[2][i].1);
    }
    v
}

fn wil_x_integrand(al: &[M4; 4], k: [f64; 3], qs_lat: &[f64], w_null: &[f64], r: f64) -> f64 {
    let pm = projw(al, k, 0.0, r, -1.0);
    let e2 = ewil(k, 0.0, r);
    let vd = vertex_x_wil(al, k);
    let pmv = mmul(&pm, &vd);
    let mut acc = 0.0f64;
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let pp = projw(al, kq, 0.0, r, 1.0);
        let e1 = ewil(kq, 0.0, r);
        let tr = mtrace_re(&mmul(&mmul(&pp, &vd), &pmv));
        acc += w_null[qi] * tr * 2.0 / (e1 + e2);
    }
    acc
}

// ================= 共通: GL / BZ / 重み / 外挿 =================

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut xs = vec![0.0f64; n];
    let mut ws = vec![0.0f64; n];
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0f64, x);
            for kk in 2..=n {
                let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        xs[i] = x;
        let (mut p0, mut p1) = (1.0f64, x);
        for kk in 2..=n {
            let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
            p0 = p1;
            p1 = p2;
        }
        let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
        ws[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (xs, ws)
}

/// nodes1d: セル [lo, hi] を edges で分割した GL 節点列
fn make_nodes(edges: &[f64], gl: &(Vec<f64>, Vec<f64>)) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for w2 in edges.windows(2) {
        let (a, b) = (w2[0], w2[1]);
        let (cc, hh) = (0.5 * (a + b), 0.5 * (b - a));
        for (x, wgt) in gl.0.iter().zip(&gl.1) {
            out.push((cc + hh * x, wgt * hh));
        }
    }
    out
}

fn bz_sum<F: Fn([f64; 3]) -> f64 + Sync>(nodes: &[(f64, f64)], nthreads: usize, f: F) -> f64 {
    let n1 = nodes.len();
    let mut rows: Vec<Option<f64>> = Vec::new();
    rows.resize_with(n1, || None);
    let chunk = n1.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            let f = &f;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let ix = t * chunk + i;
                    let (kx, wx) = nodes[ix];
                    let mut acc = 0.0f64;
                    for &(ky, wy) in nodes.iter() {
                        for &(kz, wz) in nodes.iter() {
                            acc += wy * wz * f([kx, ky, kz]);
                        }
                    }
                    *slot = Some(acc * wx);
                }
            });
        }
    });
    rows.into_iter().map(|o| o.unwrap()).sum::<f64>() / (2.0 * PI).powi(3)
}

fn null_weights(qs: &[f64; 4]) -> [f64; 4] {
    let mut a = [[0.0f64; 4]; 4];
    let mut rhs = [0.0f64; 4];
    for (i, &q) in qs.iter().enumerate() {
        a[0][i] = 1.0;
        a[1][i] = q * q;
        a[2][i] = q.powi(4);
        a[3][i] = q.powi(4) * (q * q).ln();
    }
    rhs[3] = 1.0;
    let mut m = a;
    for col in 0..4 {
        let piv = (col..4)
            .max_by(|&r1, &r2| m[r1][col].abs().partial_cmp(&m[r2][col].abs()).unwrap())
            .unwrap();
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for row in col + 1..4 {
            let f = m[row][col] / d;
            for c in col..4 {
                m[row][c] -= f * m[col][c];
            }
            rhs[row] -= f * rhs[col];
        }
    }
    let mut w = [0.0f64; 4];
    for col in (0..4).rev() {
        let mut s = rhs[col];
        for c in col + 1..4 {
            s -= m[col][c] * w[c];
        }
        w[col] = s / m[col][col];
    }
    w
}

fn fit_a0(avals: &[(f64, f64)], basis: &dyn Fn(f64) -> Vec<f64>) -> f64 {
    let p = basis(avals[0].0).len();
    let mut ata = vec![0.0f64; p * p];
    let mut atb = vec![0.0f64; p];
    for &(a, y) in avals {
        let bs = basis(a);
        for i in 0..p {
            for j in 0..p {
                ata[j + i * p] += bs[i] * bs[j];
            }
            atb[i] += bs[i] * y;
        }
    }
    let mut m = ata;
    let mut rr = atb;
    for col in 0..p {
        let piv = (col..p)
            .max_by(|&r1, &r2| m[col + r1 * p].abs().partial_cmp(&m[col + r2 * p].abs()).unwrap())
            .unwrap();
        for c in 0..p {
            m.swap(c + col * p, c + piv * p);
        }
        rr.swap(col, piv);
        let d = m[col + col * p];
        for row in col + 1..p {
            let f = m[col + row * p] / d;
            for c in col..p {
                m[c + row * p] -= f * m[c + col * p];
            }
            rr[row] -= f * rr[col];
        }
    }
    let mut sol = vec![0.0f64; p];
    for col in (0..p).rev() {
        let mut s = rr[col];
        for c in col + 1..p {
            s -= m[c + col * p] * sol[c];
        }
        sol[col] = s / m[col + col * p];
    }
    sol[0]
}

fn main() {
    self_test();
    println!("=== v26.8-X v268x_completion — X チャネル 2 比 (PRED-016 の 4 比完成) ===\n");
    println!("事前登録: spec §12.2-4/§12.4-5 (bc644d4)。的: staggered X = 4·A_oracle (Z=2 の");
    println!("凍結 tree 正規化) / Wilson X = A_oracle/2。登録バー (1%) 未達なら interim (正直判定)。\n");
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
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
    let al = alphas();
    let qset = [0.3f64, 0.6, 0.9, 1.2];
    let w_null = null_weights(&qset);
    let a_or = a_oracle_1d();
    let gl14 = gauss_legendre(14);

    // ---- [S0] staggered X 頂点の構造回帰 (遠ノード点でも所定の対角/flip 構造) ----
    {
        // 構成的定理 (v268z): piece1 は s 対角 × (−1)^{sx}·sin kx·sin 2kz,
        // piece2 は (sx,sy) flip × (−1)^{sz}·(符号)·sin kz·sin 2kx。数値で全 64 成分照合。
        let k = [0.9f64, 0.4, -1.3];
        let q = 0.17f64;
        let v = vertex8(&t_x_split(false), k, q);
        let mut worst = 0.0f64;
        for s in 0..8usize {
            for s2 in 0..8usize {
                let got = v[s + s2 * 8];
                let sgx = if sbit(s, 0) == 0 { 1.0 } else { -1.0 };
                let sgz = if sbit(s, 2) == 0 { 1.0 } else { -1.0 };
                // 開発記録 (run1 → run2): piece2 の期待符号を −0.5 と誤記 (検査コード側の
                // 転記ミス — 頂点構成は v268z 認証のまま)。両片とも +(1/2):
                //   V = +(1/2)(−1)^{sx} sin kx sin 2kz [対角] + (1/2)(−1)^{sz} sin kz sin 2kx [flip]
                let want = if s2 == s {
                    0.5 * sgx * k[0].sin() * (2.0 * k[2]).sin()
                } else if s2 == (s ^ 3) {
                    0.5 * sgz * k[2].sin() * (2.0 * k[0]).sin()
                } else {
                    0.0
                };
                worst = worst.max((got - want).abs());
            }
        }
        check(
            "[S0] staggered X 頂点の閉形式構造 (対角 Γx·sinkx·sin2kz + flip Γz·sinkz·sin2kx)",
            worst < 1e-12,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- ladder 共通 ----
    let ladder = [0.35f64, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045];
    let stag_edges = |a: f64| -> Vec<f64> {
        let c = PI / 2.0;
        let r1 = (10.0 * a * qset[3]).min(1.2).max(0.4);
        let r2 = (r1 / 3.0).max(2.0 * (2.0 * a * qset[0]).min(r1 / 3.0));
        let r3 = (2.0 * a * qset[0]).min(r2 / 2.0);
        vec![0.0, c - r1, c - r2, c - r3, c, c + r3, c + r2, c + r1, PI]
    };
    let wil_edges = |a: f64| -> Vec<f64> {
        let r1 = (10.0 * a * qset[3]).min(1.2).max(0.4);
        let r2 = (r1 / 3.0).max(2.0 * (2.0 * a * qset[0]).min(r1 / 3.0));
        let r3 = (2.0 * a * qset[0]).min(r2 / 2.0);
        vec![-PI, -r1, -r2, -r3, 0.0, r3, r2, r1, PI]
    };

    // ---- staggered X ladder ----
    let mut avals_sx: Vec<(f64, f64)> = Vec::new();
    println!("    [A(a) 表 — staggered X (的 4A_or)] a | A/4A_or");
    for &a in &ladder {
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let nodes = make_nodes(&stag_edges(a), &gl14);
        let aval = bz_sum(&nodes, nthreads, |k| {
            chi_x_integrand(k, &qs_lat, &w_null, 0.0, false)
        }) / a.powi(4);
        println!("      a = {:.3}: {:.4} ({} s)", a, aval / (4.0 * a_or), t0.elapsed().as_secs());
        avals_sx.push((a, aval));
    }

    // ---- Wilson X ladder ----
    let mut avals_wx: Vec<(f64, f64)> = Vec::new();
    println!("    [A(a) 表 — Wilson X (的 A_or/2)] a | A/(A_or/2)");
    for &a in &ladder {
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let nodes = make_nodes(&wil_edges(a), &gl14);
        let aval = bz_sum(&nodes, nthreads, |k| {
            wil_x_integrand(&al, k, &qs_lat, &w_null, 1.0)
        }) / a.powi(4);
        println!("      a = {:.3}: {:.4} ({} s)", a, aval / (0.5 * a_or), t0.elapsed().as_secs());
        avals_wx.push((a, aval));
    }

    // ---- [S1] 求積自己整合 ----
    {
        let a = *ladder.last().unwrap();
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let gl18 = gauss_legendre(18);
        let f_s = bz_sum(&make_nodes(&stag_edges(a), &gl18), nthreads, |k| {
            chi_x_integrand(k, &qs_lat, &w_null, 0.0, false)
        }) / a.powi(4);
        let f_w = bz_sum(&make_nodes(&wil_edges(a), &gl18), nthreads, |k| {
            wil_x_integrand(&al, k, &qs_lat, &w_null, 1.0)
        }) / a.powi(4);
        let rel = (f_s / avals_sx.last().unwrap().1 - 1.0)
            .abs()
            .max((f_w / avals_wx.last().unwrap().1 - 1.0).abs());
        check(
            "[S1] 求積自己整合: a_min で GL 14 → 18 の変化 < 1e-3 (両エンジン)",
            rel < 1e-3,
            format!("max 相対変化 = {:.1e}", rel),
        );
    }

    // ---- [S2] staggered X 外挿 (凍結モデル a² / a+a²) ----
    let (sx_q, sx_l);
    {
        sx_q = fit_a0(&avals_sx, &|a: f64| vec![1.0, a * a]);
        sx_l = fit_a0(&avals_sx, &|a: f64| vec![1.0, a, a * a]);
        let (rq, rl) = (sx_q / (4.0 * a_or), sx_l / (4.0 * a_or));
        let dev = (rq - 1.0).abs().max((rl - 1.0).abs()).max((rq - rl).abs());
        println!(
            "    [S2 外挿 stag-X] 比 = {:.4} (a²) / {:.4} (a+a²), spread {:.3}",
            rq, rl, (rq - rl).abs()
        );
        check(
            "[S2] staggered X: A₀/(4A_oracle) = 1 ± 0.02 (凍結 2 モデル)",
            dev < 0.02,
            format!("偏差 max = {:.3}", dev),
        );
    }

    // ---- [S3] Wilson X 外挿 (凍結モデル a+a² / a²+a⁴) ----
    let (wx_l, wx_q);
    {
        wx_l = fit_a0(&avals_wx, &|a: f64| vec![1.0, a, a * a]);
        wx_q = fit_a0(&avals_wx, &|a: f64| vec![1.0, a * a, a.powi(4)]);
        let (rl, rq) = (wx_l / (0.5 * a_or), wx_q / (0.5 * a_or));
        let dev = (rl - 1.0).abs().max((rq - 1.0).abs()).max((rl - rq).abs());
        println!(
            "    [S3 外挿 Wil-X] 比 = {:.4} (a+a²) / {:.4} (a²+a⁴), spread {:.3}",
            rl, rq, (rl - rq).abs()
        );
        check(
            "[S3] Wilson X: A₀/(A_oracle/2) = 1 ± 0.02 (凍結 2 モデル)",
            dev < 0.02,
            format!("偏差 max = {:.3}", dev),
        );
    }

    // ---- [S4] PRED-016 総括 (4 比) ----
    {
        // D の 2 比は v26.8-B/C の公表値 (凍結)
        let d_stag = [1.0102f64, 0.9908]; // v268b
        let d_wil = [0.9872f64, 1.0058]; // v268c
        let x_stag = [sx_q / (4.0 * a_or), sx_l / (4.0 * a_or)];
        let x_wil = [wx_l / (0.5 * a_or), wx_q / (0.5 * a_or)];
        println!("    [S4 PRED-016 総括 (比, 2 モデル)]");
        println!("      A_D^stag/2A_or  = {:.4}/{:.4} (v26.8-B)", d_stag[0], d_stag[1]);
        println!("      A_D^Wil/A_or    = {:.4}/{:.4} (v26.8-C)", d_wil[0], d_wil[1]);
        println!("      A_X^stag/4A_or  = {:.4}/{:.4} (本版 — Z=2 の tree 正規化込み)", x_stag[0], x_stag[1]);
        println!("      A_X^Wil/(A_or/2)= {:.4}/{:.4} (本版)", x_wil[0], x_wil[1]);
        let all: Vec<f64> = [d_stag, d_wil, x_stag, x_wil].concat();
        let worst = all.iter().map(|r| (r - 1.0).abs()).fold(0.0f64, f64::max);
        println!(
            "      ⇒ 4/4 比が 1 ± {:.1}% — **登録バー (1%/系統 0.5%) は未達 → PRED-016 は interim** (正直判定)",
            100.0 * worst
        );
        check(
            "[S4] PRED-016 総括: 4/4 比が 1 ± 2% (登録バー 1% 未達の interim を明記)",
            worst < 0.02,
            format!("max 偏差 = {:.3}", worst),
        );
    }

    // ---- [S5] 変異 ----
    {
        let a = 0.125f64;
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let nodes = make_nodes(&stag_edges(a), &gl14);
        let bad = bz_sum(&nodes, nthreads, |k| {
            chi_x_integrand(k, &qs_lat, &w_null, 0.0, true)
        }) / a.powi(4);
        let good = avals_sx[3].1;
        check(
            "[S5] 変異: staggered X piece2 の η parity 破り (ε 3→2) → A が > 10% 変化",
            (bad / good - 1.0).abs() > 0.10,
            format!("相対変化 = {:.3}", (bad / good - 1.0).abs()),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.8-X".into())),
        ("kind".into(), Json::Str("x_channel_completion".into())),
        ("pred".into(), Json::Str("PRED-016 (4 比 — interim)".into())),
        (
            "x_stag_ladder".into(),
            Json::Arr(
                avals_sx
                    .iter()
                    .map(|&(a, v)| {
                        Json::Obj(vec![
                            ("a".into(), Json::Num(a)),
                            ("ratio_4Aor".into(), Json::Num(v / (4.0 * a_or))),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "x_wil_ladder".into(),
            Json::Arr(
                avals_wx
                    .iter()
                    .map(|&(a, v)| {
                        Json::Obj(vec![
                            ("a".into(), Json::Num(a)),
                            ("ratio_halfAor".into(), Json::Num(v / (0.5 * a_or))),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("x_stag_a0_ratio".into(), Json::Num(sx_q / (4.0 * a_or))),
        ("x_wil_a0_ratio".into(), Json::Num(wx_l / (0.5 * a_or))),
    ]);
    let p = write_artifact("results/v268x_completion.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **4/4 比の二離散化 universality が 2% で成立 — PRED-016 は interim (登録バー 1% への精緻化が残作業)**"
        } else {
            "FAIL — 分岐 (b) X 転写の one-loop 破綻 / (c) 器械。欄が一次ソース"
        }
    );
    println!(
        "\n総合判定: {} ({} s)",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" },
        t0.elapsed().as_secs()
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
