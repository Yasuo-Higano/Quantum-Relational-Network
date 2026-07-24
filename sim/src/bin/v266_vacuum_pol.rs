//! v26.6 重力真空偏極の監査 I — 接触項・停留背景・縦 Ward・projector 分解 (第二十七期, 経路 B)
//!
//! 事前登録: paper/grav-vacuum-polarization-spec.md (コミット 157ca53 で実装前に凍結)。
//! 数値ゲート・判定分岐は同 §4 が一次ソース。projector 代数と ŷ チャネル辞書は
//! proofs/Projector.lean (定理 13 本) — D/X = spin-2, S = P0s, **L = yy = 純ゲージ**。
//!
//! 結合則 scheme BOND-A: 方向 i のボンド振幅 t → t·(1 + h_ii(中点))^{−1/2}。
//!   一次頂点 V_A = −(1/2) T_A(q)、二次接触項 S_AB = (3/4) δ_dir Σ c_A c_B t (c†c+h.c.)、
//!   質量項は h 非結合 (正準規格化)、counterterm Γ_ren = E − Λ Σ_x √g (Λ = −e_y)。
//! 完全静的核 (モード規格化 w_q = 1/2 [q≠0] / 1 [q=0] で q=0 連続):
//!   k̂_AB(q) = (3/4) e_A δ_AB − (1/4) χ_AB(q) − (Λ/4)(J − 2I)_AB
//!
//! 検査 (凍結ゲート — spec §4):
//!  [S0] dense FD Hessian 完全性: exact ボンド因子の H[ε] を再対角化した Richardson
//!       中心差分 = ⟨S⟩ − χ_V (N=8, m∈{0,0.5}, 対角 5 モード + 交差 3 モード +
//!       cos/sin 等値)。abs Δ ≤ max(1e-5·|FD|, 1e-6)
//!  [S1] block 3×3 χ 行列の dense 照合 (abs 1e-9) + [S1b] cos モード写像
//!       2 χ^TT_cos = χ^complex (abs 1e-9, 交差込み)
//!  [S2] tadpole: block = 解析 k 和 (abs 1e-10, N∈{8,32}) / e_y = e_z (相対 1e-13)
//!  [S3] 背景停留 (要件 0): Λ = −e_y で |dΓ_ren/dε| ≤ 1e-6 (uni-yy/zz, 3 段 Richardson)
//!       / Δ_x(N) = |e_x/e_y − 1| が N∈{8,16,32,64} で単調減少 (m=0.5)
//!  [S4] 接触項の q 非依存の器械認証: |⟨T_i(2qŷ)⟩|/V ≤ 1e-12 (dense, j∈{1,2,3})
//!  [S5] χ_D の完全核成分としての再現: c₁^χ[D] = v26.5 公表値 ±0.0004 (N=64) +
//!       c₁^χ[xx](ŷ) = v26.4 横値 ±0.0004
//!  [S6] 縦 Ward 汚染 (主結果): R(m) = |c₁^χ[L]| / |c₁^χ[D]| — branch (a) R∈[1/3,3]
//!       (b) R<1/3 (c) R>3。分解能 |c₁^χ[L]| > 3×窓系統
//!  [S7] スカラー混合 K_SL (新測定): c₁^χ[SL] と質量走行の記録 (分解能条件)
//!  [S8] uniform 連続性: block 一様変形 E(ε) の FD = V[(3/4)e − ¼χ(0)] (相対 1e-6,
//!       N=32, xx/yy/交差) / m=0.5 の |k̂(q₁) − k̂(0)| ≤ 5% max|k̂(0)|
//!  [S9] 変異検出: (i) 接触項落とし → S0 (ii) Λ 落とし → S3 (iii) T_zz 折返し
//!       スワップ落とし → S1 (逸脱 > 1e-4)
//!
//! 判定分岐 (凍結): (a) S0–S5, S8 PASS → 完全核の器械認証、S6 の R が主結果 /
//!   (b) S0/S1 FAIL → 接触項導出・交差器械の誤り / (c) S5/S8 FAIL → 規格化・窓の再設計。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// v26.5 公表 c₁^(2) (graviton plus, N=64) と v26.4 横値 c₁[xx](ŷ) — S5 回帰の的
const REF265_C1D: [(f64, f64); 2] = [(0.0, -0.01935), (0.5, -0.01444)];
const REF264_C1XX_QY: [(f64, f64); 2] = [(0.0, -0.01556), (0.5, -0.01264)];

// ================= dense 機構 (N=8 の認証層) =================

/// 変形モード: dir 方向の h_dir(x) = eps·c(中点)。kind: 0 = 一様, 1 = cos, 2 = sin
#[derive(Clone, Copy)]
struct Mode {
    dir: usize,
    kind: u8,
    j: usize,
    eps: f64,
}

fn mode_weight(m: &Mode, n: usize, mid2: [usize; 3]) -> f64 {
    let q = 2.0 * PI * m.j as f64 / n as f64;
    let ph = 0.5 * q * mid2[1] as f64; // q ∥ ŷ — 位相は中点 y 座標
    match m.kind {
        0 => 1.0,
        1 => ph.cos(),
        _ => ph.sin(),
    }
}

struct Bond {
    i: usize,
    j: usize,
    t: f64,
    dir: usize,
    mid2: [usize; 3],
}

fn bonds(n: usize) -> Vec<Bond> {
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut out = Vec::new();
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                let twx = if x == n - 1 { -1.0 } else { 1.0 };
                out.push(Bond {
                    i,
                    j: idx((x + 1) % n, y, z),
                    t: 0.5 * twx,
                    dir: 0,
                    mid2: [2 * x + 1, 2 * y, 2 * z],
                });
                let ey = if x % 2 == 0 { 1.0 } else { -1.0 };
                out.push(Bond {
                    i,
                    j: idx(x, (y + 1) % n, z),
                    t: 0.5 * ey,
                    dir: 1,
                    mid2: [2 * x, 2 * y + 1, 2 * z],
                });
                let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                out.push(Bond {
                    i,
                    j: idx(x, y, (z + 1) % n),
                    t: 0.5 * ez,
                    dir: 2,
                    mid2: [2 * x, 2 * y, 2 * z + 1],
                });
            }
        }
    }
    out
}

/// H[h] dense: ボンド因子 (1 + Σ_modes h)^{−1/2} (scheme BOND-A の exact 形)
fn build_h_deformed(n: usize, m_stag: f64, modes: &[Mode]) -> Vec<f64> {
    let ns = n * n * n;
    let mut h = vec![0.0f64; ns * ns];
    for b in bonds(n) {
        let mut hb = 0.0;
        for md in modes {
            if md.dir == b.dir {
                hb += md.eps * mode_weight(md, n, b.mid2);
            }
        }
        assert!(1.0 + hb > 0.0, "ボンド因子が特異 (h = {})", hb);
        let t = b.t / (1.0 + hb).sqrt();
        h[b.j + b.i * ns] += t;
        h[b.i + b.j * ns] += t;
    }
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                let i = idx(x, y, z);
                h[i + i * ns] += m_stag * sgn; // 質量項は h 非結合 (spec §2)
            }
        }
    }
    h
}

/// 真空エネルギー (下半分の総和 — jacobi_eigh は昇順)
fn e0_of(w: &[f64], ns: usize) -> f64 {
    w[..ns / 2].iter().sum()
}

fn e0_dense(n: usize, m_stag: f64, modes: &[Mode]) -> f64 {
    let ns = n * n * n;
    let h = build_h_deformed(n, m_stag, modes);
    let (w, _) = jacobi_eigh(&h, ns);
    e0_of(&w, ns)
}

/// dense E₀ の並列バッチ評価 (独立ジョブの決定的スレッド分割 — 結果は並列度に依らない)
fn e0_dense_batch(n: usize, cfgs: &[(f64, Vec<Mode>)], nthreads: usize) -> Vec<f64> {
    let mut out: Vec<Option<f64>> = Vec::new();
    out.resize_with(cfgs.len(), || None);
    let chunk = cfgs.len().div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in out.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let (m, modes) = &cfgs[t * chunk + i];
                    *slot = Some(e0_dense(n, *m, modes));
                }
            });
        }
    });
    out.into_iter().map(|o| o.unwrap()).collect()
}

/// 真空相関 G_ij = Σ_{ν occ} v_iν v_jν
fn vac_corr(v: &[f64], ns: usize) -> Vec<f64> {
    let nocc = ns / 2;
    let mut g = vec![0.0f64; ns * ns];
    for nu in 0..nocc {
        for i in 0..ns {
            let vi = v[i + nu * ns];
            if vi == 0.0 {
                continue;
            }
            for j in 0..ns {
                g[j + i * ns] += vi * v[j + nu * ns];
            }
        }
    }
    g
}

/// 実モード頂点 T_A = Σ_{b∈dir} c(中点) t_b (c†c + h.c.)
fn mode_vertex(n: usize, md: &Mode) -> Vec<f64> {
    let ns = n * n * n;
    let mut o = vec![0.0f64; ns * ns];
    for b in bonds(n) {
        if b.dir != md.dir {
            continue;
        }
        let w = mode_weight(md, n, b.mid2);
        o[b.j + b.i * ns] += b.t * w;
        o[b.i + b.j * ns] += b.t * w;
    }
    o
}

/// 実演算子対の静的交差 Lehmann: Σ_{occ→unocc} 2 M_A M_B / ΔE (示量)
fn chi_cross_real(w: &[f64], v: &[f64], ns: usize, oa: &[f64], ob: &[f64]) -> f64 {
    let nocc = ns / 2;
    let tv = |o: &[f64]| -> Vec<f64> {
        let mut t = vec![0.0f64; ns * nocc];
        for c in 0..nocc {
            for r in 0..ns {
                let mut s = 0.0;
                for k in 0..ns {
                    s += o[k + r * ns] * v[k + c * ns];
                }
                t[r + c * ns] = s;
            }
        }
        t
    };
    let ta = tv(oa);
    let tb = tv(ob);
    let mut chi = 0.0f64;
    for mu in nocc..ns {
        for nu in 0..nocc {
            let (mut ma, mut mb) = (0.0f64, 0.0f64);
            for k in 0..ns {
                let vm = v[k + mu * ns];
                ma += vm * ta[k + nu * ns];
                mb += vm * tb[k + nu * ns];
            }
            chi += 2.0 * ma * mb / (w[mu] - w[nu]);
        }
    }
    chi
}

/// dense の複素頂点 T_A(qŷ) (which: 1 = x ボンド片, 2 = y ボンド片, 3 = z ボンド片;
/// 中点位相規約 — v26.4/26.5 と同一)
fn dense_vertex(n: usize, qy: f64, which: usize) -> (Vec<f64>, Vec<f64>) {
    let ns = n * n * n;
    let mut re = vec![0.0f64; ns * ns];
    let mut im = vec![0.0f64; ns * ns];
    for b in bonds(n) {
        if b.dir + 1 != which {
            continue;
        }
        let ph = 0.5 * qy * b.mid2[1] as f64;
        let (c, s) = (ph.cos(), ph.sin());
        re[b.j + b.i * ns] += b.t * c;
        re[b.i + b.j * ns] += b.t * c;
        im[b.j + b.i * ns] += b.t * s;
        im[b.i + b.j * ns] += b.t * s;
    }
    (re, im)
}

/// dense の複素 3 頂点交差 Lehmann → 3×3 対称行列 (示量)
fn chi_dense_matrix(
    w: &[f64],
    v: &[f64],
    ns: usize,
    ops: &[(Vec<f64>, Vec<f64>); 3],
) -> [[f64; 3]; 3] {
    let nocc = ns / 2;
    let tv = |o: &(Vec<f64>, Vec<f64>)| -> (Vec<f64>, Vec<f64>) {
        let mut tre = vec![0.0f64; ns * nocc];
        let mut tim = vec![0.0f64; ns * nocc];
        for c in 0..nocc {
            for r in 0..ns {
                let (mut sr, mut si) = (0.0, 0.0);
                for k in 0..ns {
                    let vv = v[k + c * ns];
                    sr += o.0[k + r * ns] * vv;
                    si += o.1[k + r * ns] * vv;
                }
                tre[r + c * ns] = sr;
                tim[r + c * ns] = si;
            }
        }
        (tre, tim)
    };
    let ts: Vec<(Vec<f64>, Vec<f64>)> = ops.iter().map(tv).collect();
    let mut chi = [[0.0f64; 3]; 3];
    for mu in nocc..ns {
        for nu in 0..nocc {
            let mut mm = [(0.0f64, 0.0f64); 3];
            for (a, t) in ts.iter().enumerate() {
                let (mut mr, mut mi) = (0.0f64, 0.0f64);
                for k in 0..ns {
                    let vm = v[k + mu * ns];
                    mr += vm * t.0[k + nu * ns];
                    mi += vm * t.1[k + nu * ns];
                }
                mm[a] = (mr, mi);
            }
            let de = w[mu] - w[nu];
            for a in 0..3 {
                for b in a..3 {
                    let x = 2.0 * (mm[a].0 * mm[b].0 + mm[a].1 * mm[b].1) / de;
                    chi[a][b] += x;
                    if a != b {
                        chi[b][a] += x;
                    }
                }
            }
        }
    }
    chi
}

// ================= block 機構 (q ∥ ŷ — v26.4/26.5 の認証済み実装を写経) =================

fn block_h_scaled(n: usize, m: f64, cky: f64, ckz: f64, s: [f64; 3]) -> Vec<f64> {
    let dim = 4 * n;
    let mut h = vec![0.0f64; dim * dim];
    let id = |x: usize, c: usize| x + n * c;
    let add = |h: &mut Vec<f64>, a: usize, b: usize, t: f64| {
        h[b + a * dim] += t;
        h[a + b * dim] += t;
    };
    let ysgn = [1.0, -1.0, 1.0, -1.0];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            add(&mut h, id(x, c), id((x + 1) % n, c), 0.5 * tw * s[0]);
            h[id(x, c) + id(x, c) * dim] += px * ysgn[c] * cky * s[1];
            if c == 0 || c == 2 {
                add(&mut h, id(x, c), id(x, c + 1), px * zsgn[c] * ckz * s[2]);
            }
        }
        add(&mut h, id(x, 0), id(x, 3), px * m);
        add(&mut h, id(x, 1), id(x, 2), px * m);
    }
    h
}

/// q∥ŷ の頂点 (ブロック対, 実行列, 格納 v[標的 + 元·dim])。which: 1 = T_xx, 2 = T_yy,
/// 3 = T_zz。sw = 折返しスワップ (v26.4 開発記録 — 縮約は「元」側)。
fn vertex_qy(n: usize, ky: f64, ckz: f64, q: f64, sw: bool, which: usize) -> Vec<f64> {
    let dim = 4 * n;
    let id = |x: usize, c: usize| x + n * c;
    let tc = |c: usize| -> usize {
        if sw {
            [1usize, 0, 3, 2][c]
        } else {
            c
        }
    };
    let mut v = vec![0.0f64; dim * dim];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let kyc = ky + if c == 1 || c == 3 { PI } else { 0.0 };
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            if which == 1 {
                v[id((x + 1) % n, tc(c)) + id(x, c) * dim] += 0.5 * tw;
                v[id(x, tc(c)) + id((x + 1) % n, c) * dim] += 0.5 * tw;
            }
            if which == 2 {
                v[id(x, tc(c)) + id(x, c) * dim] += px * (kyc + q / 2.0).cos();
            }
            if which == 3 && (c == 0 || c == 2) {
                let coef = px * zsgn[c] * ckz;
                v[id(x, tc(c + 1)) + id(x, c) * dim] += coef;
                v[id(x, tc(c)) + id(x, c + 1) * dim] += coef;
            }
        }
    }
    v
}

/// ブロック対の 3 頂点交差 Lehmann → 3×3 (縮約は「元」側 o[r + k·dim])
fn chi_pair_matrix(
    w1: &[f64],
    v1: &[f64],
    w2: &[f64],
    v2: &[f64],
    dim: usize,
    ops: &[Vec<f64>; 3],
) -> [[f64; 3]; 3] {
    let nocc = dim / 2;
    let tv = |o: &[f64]| -> Vec<f64> {
        let mut t = vec![0.0f64; dim * nocc];
        for ccol in 0..nocc {
            for r in 0..dim {
                let mut s = 0.0;
                for k in 0..dim {
                    s += o[r + k * dim] * v1[k + ccol * dim];
                }
                t[r + ccol * dim] = s;
            }
        }
        t
    };
    let ts: Vec<Vec<f64>> = ops.iter().map(|o| tv(o)).collect();
    let mut chi = [[0.0f64; 3]; 3];
    for mu in nocc..dim {
        for nu in 0..nocc {
            let mut mm = [0.0f64; 3];
            for (a, t) in ts.iter().enumerate() {
                let mut s = 0.0;
                for k in 0..dim {
                    s += v2[k + mu * dim] * t[k + nu * dim];
                }
                mm[a] = s;
            }
            let de = w2[mu] - w1[nu];
            for a in 0..3 {
                for b in a..3 {
                    let x = 2.0 * mm[a] * mm[b] / de;
                    chi[a][b] += x;
                    if a != b {
                        chi[b][a] += x;
                    }
                }
            }
        }
    }
    chi
}

/// 3×3 χ 行列の q 走査 (kz 行ごとに決定的スレッド分割)。
/// mutate = true: T_zz 頂点の折返しスワップを落とす (S9-iii 用)。
fn chi_scan_matrix(
    n: usize,
    m: f64,
    js: &[usize],
    nthreads: usize,
    mutate: bool,
) -> Vec<[[f64; 3]; 3]> {
    let nb = n / 2;
    let mut rows: Vec<Option<Vec<[[f64; 3]; 3]>>> = Vec::new();
    rows.resize_with(nb, || None);
    let chunk = nb.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let jz = t * chunk + i;
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    let mut eigs: Vec<(Vec<f64>, Vec<f64>)> = Vec::with_capacity(nb);
                    for jy in 0..nb {
                        let cky = (2.0 * PI * jy as f64 / n as f64).cos();
                        let h = block_h_scaled(n, m, cky, ckz, [1.0, 1.0, 1.0]);
                        eigs.push(jacobi_eigh(&h, 4 * n));
                    }
                    let dim = 4 * n;
                    let mut acc = vec![[[0.0f64; 3]; 3]; js.len()];
                    for (ji, &j) in js.iter().enumerate() {
                        let q = 2.0 * PI * j as f64 / n as f64;
                        for jy in 0..nb {
                            let ky = 2.0 * PI * jy as f64 / n as f64;
                            let mut jt = jy + j;
                            let mut sw = false;
                            while jt >= nb {
                                jt -= nb;
                                sw = !sw;
                            }
                            let ox = vertex_qy(n, ky, ckz, q, sw, 1);
                            let oy = vertex_qy(n, ky, ckz, q, sw, 2);
                            let oz = vertex_qy(n, ky, ckz, q, if mutate { false } else { sw }, 3);
                            let (w1, v1) = &eigs[jy];
                            let (w2, v2) = &eigs[jt];
                            let c = chi_pair_matrix(w1, v1, w2, v2, dim, &[ox, oy, oz]);
                            for a in 0..3 {
                                for b in 0..3 {
                                    acc[ji][a][b] += c[a][b];
                                }
                            }
                        }
                    }
                    *slot = Some(acc);
                }
            });
        }
    });
    let vol = (n * n * n) as f64;
    let mut out = vec![[[0.0f64; 3]; 3]; js.len()];
    for r in rows {
        for (ji, c) in r.unwrap().iter().enumerate() {
            for a in 0..3 {
                for b in 0..3 {
                    out[ji][a][b] += c[a][b] / vol;
                }
            }
        }
    }
    out
}

/// block tadpole: e_i = ⟨T_i(0)⟩/V (占有状態の対角和; 決定的スレッド分割)
fn tadpole_block(n: usize, m: f64, nthreads: usize) -> [f64; 3] {
    let nb = n / 2;
    let mut rows: Vec<Option<[f64; 3]>> = Vec::new();
    rows.resize_with(nb, || None);
    let chunk = nb.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let jz = t * chunk + i;
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    let dim = 4 * n;
                    let mut acc = [0.0f64; 3];
                    for jy in 0..nb {
                        let ky = 2.0 * PI * jy as f64 / n as f64;
                        let cky = ky.cos();
                        let h = block_h_scaled(n, m, cky, ckz, [1.0, 1.0, 1.0]);
                        let (_w, v) = jacobi_eigh(&h, dim);
                        for (wi, which) in [1usize, 2, 3].iter().enumerate() {
                            let o = vertex_qy(n, ky, ckz, 0.0, false, *which);
                            let mut s = 0.0;
                            for nu in 0..dim / 2 {
                                for r in 0..dim {
                                    let mut a = 0.0;
                                    for k in 0..dim {
                                        a += o[r + k * dim] * v[k + nu * dim];
                                    }
                                    s += v[r + nu * dim] * a;
                                }
                            }
                            acc[wi] += s;
                        }
                    }
                    *slot = Some(acc);
                }
            });
        }
    });
    let vol = (n * n * n) as f64;
    let mut out = [0.0f64; 3];
    for r in rows {
        let r = r.unwrap();
        for a in 0..3 {
            out[a] += r[a] / vol;
        }
    }
    out
}

/// 解析 tadpole: e_i = −(1/2V) Σ_k cos²k_i / √(Σcos² + m²) (twist-x 格子)。
/// 因子 1/2: N³ 個の k 全てに負エネルギー枝を数えると band pairing (k ↔ k+π 系)で
/// 占有状態 (N³/2 個) を厳密に 2 重計上する — 実装時にこの因子を落としており、
/// 走行前の独立 python 照合 (dense e_y = −0.198 vs 解析 −0.396) が検出した (開発記録)。
fn tadpole_analytic(n: usize, m: f64) -> [f64; 3] {
    let mut acc = [0.0f64; 3];
    for jx in 0..n {
        let kx = (2.0 * jx as f64 + 1.0) * PI / n as f64;
        for jy in 0..n {
            let ky = 2.0 * PI * jy as f64 / n as f64;
            for jz in 0..n {
                let kz = 2.0 * PI * jz as f64 / n as f64;
                let (cx, cy, cz) = (kx.cos(), ky.cos(), kz.cos());
                let e = (cx * cx + cy * cy + cz * cz + m * m).sqrt();
                acc[0] -= cx * cx / e;
                acc[1] -= cy * cy / e;
                acc[2] -= cz * cz / e;
            }
        }
    }
    let vol2 = 2.0 * (n * n * n) as f64;
    [acc[0] / vol2, acc[1] / vol2, acc[2] / vol2]
}

/// block 一様変形の真空エネルギー E(s) (s = ボンド因子 3 方向; 決定的スレッド分割)
fn e0_block_scaled(n: usize, m: f64, s: [f64; 3], nthreads: usize) -> f64 {
    let nb = n / 2;
    let mut rows: Vec<Option<f64>> = Vec::new();
    rows.resize_with(nb, || None);
    let chunk = nb.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let jz = t * chunk + i;
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    let dim = 4 * n;
                    let mut acc = 0.0;
                    for jy in 0..nb {
                        let cky = (2.0 * PI * jy as f64 / n as f64).cos();
                        let h = block_h_scaled(n, m, cky, ckz, s);
                        let (w, _) = jacobi_eigh(&h, dim);
                        acc += w[..dim / 2].iter().sum::<f64>();
                    }
                    *slot = Some(acc);
                }
            });
        }
    });
    // ブロック層は (ky, kz) の n/2 × n/2 のみ — 全 BZ は ±ky, ±kz の 4 重に折返し済みの
    // 4 成分基底で被覆される (v25.1/v26.3 の分解)。x 鎖 n × 4 成分 × (n/2)² = n³ ✓
    rows.into_iter().map(|r| r.unwrap()).sum()
}

/// 最小二乗 (v26.3–26.5 写経)
fn lstsq(xs: &[Vec<f64>], y: &[f64]) -> Vec<f64> {
    let p = xs.len();
    let mut a = vec![0.0f64; p * p];
    let mut b = vec![0.0f64; p];
    for i in 0..p {
        for jj in 0..p {
            a[jj + i * p] = xs[i].iter().zip(&xs[jj]).map(|(u, v)| u * v).sum();
        }
        b[i] = xs[i].iter().zip(y).map(|(u, v)| u * v).sum();
    }
    let mut idx: Vec<usize> = (0..p).collect();
    for col in 0..p {
        let piv = (col..p)
            .max_by(|&r1, &r2| {
                a[col + idx[r1] * p]
                    .abs()
                    .partial_cmp(&a[col + idx[r2] * p].abs())
                    .unwrap()
            })
            .unwrap();
        idx.swap(col, piv);
        let d = a[col + idx[col] * p];
        for r in col + 1..p {
            let f = a[col + idx[r] * p] / d;
            for cc in col..p {
                a[cc + idx[r] * p] -= f * a[cc + idx[col] * p];
            }
            b[idx[r]] -= f * b[idx[col]];
        }
    }
    let mut out = vec![0.0f64; p];
    for col in (0..p).rev() {
        let mut s = b[idx[col]];
        for cc in col + 1..p {
            s -= a[cc + idx[col] * p] * out[cc];
        }
        out[col] = s / a[col + idx[col] * p];
    }
    out
}

/// c₁ フィット (v26.3 プロトコル: 窓 (1,4)/(2,5) × モデル {1,q²} (+q²ln(1/q) [m=0]))
fn fit_c1(n: usize, m: f64, chi_of_j: &dyn Fn(usize) -> f64) -> (f64, f64) {
    let q = |j: usize| 2.0 * PI * j as f64 / n as f64;
    let mut c1s = Vec::new();
    for (lo, hi) in [(1usize, 4usize), (2, 5)] {
        let jr: Vec<usize> = (lo..=hi).collect();
        let y: Vec<f64> = jr
            .iter()
            .map(|&j| (chi_of_j(j) - chi_of_j(0)) / (q(j) * q(j)))
            .collect();
        let ones: Vec<f64> = jr.iter().map(|_| 1.0).collect();
        let q2: Vec<f64> = jr.iter().map(|&j| q(j) * q(j)).collect();
        let c = lstsq(&[ones.clone(), q2.clone()], &y);
        c1s.push(c[0]);
        if m == 0.0 {
            let ql: Vec<f64> = jr
                .iter()
                .map(|&j| q(j) * q(j) * (1.0 / q(j)).ln())
                .collect();
            let c = lstsq(&[ones, q2, ql], &y);
            c1s.push(c[0]);
        }
    }
    let mut sorted = c1s.clone();
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = sorted[sorted.len() / 2];
    let spread = sorted.last().unwrap() - sorted.first().unwrap();
    (med, spread)
}

/// Richardson 付き中心 2 階差分 (ε, ε/2)
fn d2_richardson(f: &dyn Fn(f64) -> f64, eps: f64) -> f64 {
    let f0 = f(0.0);
    let d = |e: f64| (f(e) + f(-e) - 2.0 * f0) / (e * e);
    (4.0 * d(eps / 2.0) - d(eps)) / 3.0
}

/// 交差 2 階差分 (4 点) + Richardson
fn d2_cross_richardson(f: &dyn Fn(f64, f64) -> f64, eps: f64) -> f64 {
    let d = |e: f64| (f(e, e) - f(e, -e) - f(-e, e) + f(-e, -e)) / (4.0 * e * e);
    (4.0 * d(eps / 2.0) - d(eps)) / 3.0
}

fn main() {
    self_test();
    println!(
        "=== v26.6 重力真空偏極の監査 I — 接触項・停留背景・縦 Ward・projector 分解 (第二十七期) ===\n"
    );
    println!("事前登録: paper/grav-vacuum-polarization-spec.md §4 (コミット 157ca53 で凍結)");
    println!("          (a) S0–S5, S8 PASS → 完全核の器械認証、S6 の R が主結果 /");
    println!("          (b) S0/S1 FAIL → 接触項・交差器械の誤り / (c) S5/S8 FAIL → 規格化の再設計");
    println!("チャネル辞書は proofs/Projector.lean (D/X = spin-2, S = P0s, L = yy = 純ゲージ)\n");
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
    let ms = [0.0f64, 0.5];

    // ================= dense 基底系 (N=8) =================
    let nd = 8usize;
    let nsd = nd * nd * nd;
    let vold = nsd as f64;
    let mut dense_eig: Vec<Option<(Vec<f64>, Vec<f64>)>> = Vec::new();
    dense_eig.resize_with(ms.len(), || None);
    std::thread::scope(|sc| {
        for (slot, &m) in dense_eig.iter_mut().zip(&ms) {
            sc.spawn(move || {
                let h = build_h_deformed(nd, m, &[]);
                *slot = Some(jacobi_eigh(&h, nsd));
            });
        }
    });
    let dense_eig: Vec<(Vec<f64>, Vec<f64>)> = dense_eig.into_iter().map(|o| o.unwrap()).collect();
    // 真空相関と tadpole (dense)
    let mut e_dense = [[0.0f64; 3]; 2]; // [mi][dir]
    let mut gs: Vec<Vec<f64>> = Vec::new();
    for (mi, _m) in ms.iter().enumerate() {
        let (_, v) = &dense_eig[mi];
        let g = vac_corr(v, nsd);
        for b in bonds(nd) {
            e_dense[mi][b.dir] += 2.0 * b.t * g[b.j + b.i * nsd];
        }
        for d in 0..3 {
            e_dense[mi][d] /= vold;
        }
        gs.push(g);
    }
    println!(
        "    [dense 基底] N=8 対角化 + 真空相関 完了 ({} s) — e = ({:.6}, {:.6}, {:.6}) [m=0]",
        t0.elapsed().as_secs(),
        e_dense[0][0],
        e_dense[0][1],
        e_dense[0][2]
    );

    // ---- dense FD 点のバッチ構築 (S0/S0c/S3a/S9 で共有 — 並列評価) ----
    let eps = 0.02;
    let diag_defs: [(usize, u8, usize, &str); 5] = [
        (0, 0, 0, "uni-xx"),
        (1, 0, 0, "uni-yy"),
        (0, 1, 1, "cos-xx j1"),
        (1, 1, 2, "cos-yy j2"),
        (2, 1, 2, "cos-zz j2"),
    ];
    let cross_defs: [(usize, u8, usize, usize, u8, usize, &str); 3] = [
        (0, 0, 0, 1, 0, 0, "uni xx-yy"),
        (0, 0, 0, 2, 0, 0, "uni xx-zz"),
        (0, 1, 2, 1, 1, 2, "cos xx-yy j2"),
    ];
    let mut cfgs: Vec<(f64, Vec<Mode>)> = Vec::new();
    let mut i_f0 = [0usize; 2];
    let mut i_diag = [[[0usize; 4]; 5]; 2]; // [mi][mode][pt: +ε −ε +ε/2 −ε/2]
    let mut i_cross = [[[0usize; 8]; 3]; 2]; // [mi][mode][pt: (±ε,±ε)×4, (±ε/2,±ε/2)×4]
    let mut i_stat = [[[0usize; 6]; 2]; 2]; // [mi][dir y/z][pt: ±ε ±ε/2 ±ε/4]
    let mut i_sin = [0usize; 4]; // m=0 のみ (sin-yy j2)
    for (mi, &m) in ms.iter().enumerate() {
        i_f0[mi] = cfgs.len();
        cfgs.push((m, vec![]));
        for (k, &(dir, kind, j, _)) in diag_defs.iter().enumerate() {
            for (p, s) in [1.0f64, -1.0, 0.5, -0.5].iter().enumerate() {
                i_diag[mi][k][p] = cfgs.len();
                cfgs.push((m, vec![Mode { dir, kind, j, eps: s * eps }]));
            }
        }
        for (k, &(d1, k1, j1, d2, k2, j2, _)) in cross_defs.iter().enumerate() {
            let mut p = 0usize;
            for lev in [1.0f64, 0.5] {
                for (sa, sb) in [(1.0f64, 1.0f64), (1.0, -1.0), (-1.0, 1.0), (-1.0, -1.0)] {
                    i_cross[mi][k][p] = cfgs.len();
                    cfgs.push((
                        m,
                        vec![
                            Mode { dir: d1, kind: k1, j: j1, eps: sa * lev * eps },
                            Mode { dir: d2, kind: k2, j: j2, eps: sb * lev * eps },
                        ],
                    ));
                    p += 1;
                }
            }
        }
        for (di, &dir) in [1usize, 2].iter().enumerate() {
            for (p, s) in [1.0f64, -1.0, 0.5, -0.5, 0.25, -0.25].iter().enumerate() {
                i_stat[mi][di][p] = cfgs.len();
                cfgs.push((m, vec![Mode { dir, kind: 0, j: 0, eps: s * eps }]));
            }
        }
    }
    for (p, s) in [1.0f64, -1.0, 0.5, -0.5].iter().enumerate() {
        i_sin[p] = cfgs.len();
        cfgs.push((0.0, vec![Mode { dir: 1, kind: 2, j: 2, eps: s * eps }]));
    }
    let evals = e0_dense_batch(nd, &cfgs, nthreads);
    println!(
        "    [dense FD] {} 点の並列バッチ評価 完了 ({} s, {} threads)",
        cfgs.len(),
        t0.elapsed().as_secs(),
        nthreads
    );
    // Richardson 組立 (点は評価済み — ここは純算術)
    let d2_from = |f0: f64, e: [f64; 4]| -> f64 {
        let de = (e[0] + e[1] - 2.0 * f0) / (eps * eps);
        let dh = (e[2] + e[3] - 2.0 * f0) / (0.25 * eps * eps);
        (4.0 * dh - de) / 3.0
    };
    let d2c_from = |e: [f64; 8]| -> f64 {
        let de = (e[0] - e[1] - e[2] + e[3]) / (4.0 * eps * eps);
        let dh = (e[4] - e[5] - e[6] + e[7]) / (eps * eps);
        (4.0 * dh - de) / 3.0
    };
    let d1_from = |e: [f64; 6]| -> f64 {
        let d1 = (e[0] - e[1]) / (2.0 * eps);
        let d2v = (e[2] - e[3]) / eps;
        let d3 = (e[4] - e[5]) / (0.5 * eps);
        let r1 = (4.0 * d2v - d1) / 3.0;
        let r2 = (4.0 * d3 - d2v) / 3.0;
        (16.0 * r2 - r1) / 15.0
    };
    let pick4 = |ix: [usize; 4]| -> [f64; 4] { [evals[ix[0]], evals[ix[1]], evals[ix[2]], evals[ix[3]]] };
    let pick6 = |ix: [usize; 6]| -> [f64; 6] {
        [evals[ix[0]], evals[ix[1]], evals[ix[2]], evals[ix[3]], evals[ix[4]], evals[ix[5]]]
    };
    let pick8 = |ix: [usize; 8]| -> [f64; 8] {
        [
            evals[ix[0]], evals[ix[1]], evals[ix[2]], evals[ix[3]], evals[ix[4]], evals[ix[5]],
            evals[ix[6]], evals[ix[7]],
        ]
    };

    // ---- [S0] dense FD Hessian 完全性 ----
    let mut fd_uni_yy_m0 = 0.0f64; // S9-i で再利用
    {
        let mut worst_rep = (0.0f64, String::new());
        for (mi, &m) in ms.iter().enumerate() {
            let (w, v) = &dense_eig[mi];
            for (k, &(dir, kind, j, name)) in diag_defs.iter().enumerate() {
                let md = Mode { dir, kind, j, eps: 0.0 };
                let fd = d2_from(evals[i_f0[mi]], pick4(i_diag[mi][k]));
                if mi == 0 && k == 1 {
                    fd_uni_yy_m0 = fd;
                }
                // 組立: ⟨S⟩ − χ_V = (3/4)Σ c² t 2G − (1/4)χ^TT
                let ov = mode_vertex(nd, &md);
                let mut s_exp = 0.0;
                for b in bonds(nd) {
                    if b.dir == dir {
                        let c = mode_weight(&md, nd, b.mid2);
                        s_exp += 0.75 * c * c * 2.0 * b.t * gs[mi][b.j + b.i * nsd];
                    }
                }
                let chi_tt = chi_cross_real(w, v, nsd, &ov, &ov);
                let asm = s_exp - 0.25 * chi_tt;
                let dev = (fd - asm).abs();
                let tol = (1e-5 * fd.abs()).max(1e-6);
                if dev / tol > worst_rep.0 {
                    worst_rep = (
                        dev / tol,
                        format!("{} m={:.1}: FD = {:.8}, 組立 = {:.8}", name, m, fd, asm),
                    );
                }
            }
            for (k, &(d1, k1, j1, d2, k2, j2, name)) in cross_defs.iter().enumerate() {
                let ma = Mode { dir: d1, kind: k1, j: j1, eps: 0.0 };
                let mb = Mode { dir: d2, kind: k2, j: j2, eps: 0.0 };
                let fd = d2c_from(pick8(i_cross[mi][k]));
                // 交差は接触項ゼロ (方向が異なる) — 組立 = −(1/4)χ^TT_AB
                let oa = mode_vertex(nd, &ma);
                let ob = mode_vertex(nd, &mb);
                let asm = -0.25 * chi_cross_real(w, v, nsd, &oa, &ob);
                let dev = (fd - asm).abs();
                let tol = (1e-5 * fd.abs()).max(1e-6);
                if dev / tol > worst_rep.0 {
                    worst_rep = (
                        dev / tol,
                        format!("{} m={:.1}: FD = {:.8}, 組立 = {:.8}", name, m, fd, asm),
                    );
                }
            }
        }
        check(
            "[S0] dense FD Hessian 完全性: d²E/dh² = ⟨S⟩ − χ_V (8 モード × m 2 種)",
            worst_rep.0 < 1.0,
            format!("worst = {:.3} × 許容 ({})", worst_rep.0, worst_rep.1),
        );
        // cos/sin 等値 (FD レベル)
        let fdc = d2_from(evals[i_f0[0]], pick4(i_diag[0][3]));
        let fds = d2_from(evals[i_f0[0]], pick4(i_sin));
        check(
            "[S0c] cos/sin モードの FD Hessian 等値 (y 並進不変性の端から端)",
            (fdc - fds).abs() < (1e-5 * fdc.abs()).max(1e-6),
            format!("cos = {:.8}, sin = {:.8}", fdc, fds),
        );
    }

    // ---- [S1] block 3×3 χ vs dense + [S1b] cos モード写像 ----
    {
        let mut worst = 0.0f64;
        let mut worst_map = 0.0f64;
        for (mi, &m) in ms.iter().enumerate() {
            let (w, v) = &dense_eig[mi];
            let blocky = chi_scan_matrix(nd, m, &[1, 2], nthreads, false);
            for (ji, &j) in [1usize, 2].iter().enumerate() {
                let qy = 2.0 * PI * j as f64 / nd as f64;
                let ops = [
                    dense_vertex(nd, qy, 1),
                    dense_vertex(nd, qy, 2),
                    dense_vertex(nd, qy, 3),
                ];
                let cd = chi_dense_matrix(w, v, nsd, &ops);
                for a in 0..3 {
                    for b in 0..3 {
                        worst = worst.max((blocky[ji][a][b] - cd[a][b] / vold).abs());
                    }
                }
                // [S1b] 実 cos モード頂点の χ^TT_cos ×2 = χ^complex (対角 3 + 交差 xx-yy)
                for (da, db) in [(0usize, 0usize), (1, 1), (2, 2), (0, 1)] {
                    let ma = Mode {
                        dir: da,
                        kind: 1,
                        j,
                        eps: 0.0,
                    };
                    let mb = Mode {
                        dir: db,
                        kind: 1,
                        j,
                        eps: 0.0,
                    };
                    let oa = mode_vertex(nd, &ma);
                    let ob = mode_vertex(nd, &mb);
                    let cc = chi_cross_real(w, v, nsd, &oa, &ob);
                    worst_map = worst_map.max((2.0 * cc - cd[da][db]).abs() / vold);
                }
            }
        }
        check(
            "[S1] block 3×3 χ 行列 = dense (N=8, j∈{1,2}, m 2 種, 全 9 成分)",
            worst < 1e-9,
            format!("max|Δ| = {:.1e} ({} s)", worst, t0.elapsed().as_secs()),
        );
        check(
            "[S1b] cos モード写像 2χ^TT_cos = χ^complex (対角 + 交差)",
            worst_map < 1e-9,
            format!("max|Δ|/V = {:.1e}", worst_map),
        );
    }

    // ---- [S2] tadpole の解析照合 + W 対称 ----
    let mut e_block = std::collections::BTreeMap::new(); // (n, mi) -> [f64;3]
    {
        let mut worst = 0.0f64;
        for &n in &[8usize, 32] {
            for (mi, &m) in ms.iter().enumerate() {
                let eb = tadpole_block(n, m, nthreads);
                let ea = tadpole_analytic(n, m);
                for d in 0..3 {
                    worst = worst.max((eb[d] - ea[d]).abs());
                }
                e_block.insert((n, mi), eb);
            }
        }
        // dense 側も解析と照合 (N=8)
        for (mi, _) in ms.iter().enumerate() {
            let ea = tadpole_analytic(nd, ms[mi]);
            for d in 0..3 {
                worst = worst.max((e_dense[mi][d] - ea[d]).abs());
            }
        }
        check(
            "[S2a] tadpole e_i: block/dense = 解析 k 和 (N∈{8,32}, m 2 種)",
            worst < 1e-10,
            format!("max|Δ| = {:.1e}", worst),
        );
        let mut worst_wz = 0.0f64;
        for &n in &[8usize, 32] {
            for mi in 0..ms.len() {
                let e = e_block[&(n, mi)];
                worst_wz = worst_wz.max((e[1] / e[2] - 1.0).abs());
            }
        }
        check(
            "[S2b] e_y = e_z (W 対称性)",
            worst_wz < 1e-13,
            format!("max|e_y/e_z − 1| = {:.1e}", worst_wz),
        );
    }

    // ---- [S3] 背景停留 (要件 0) ----
    {
        // [S3a] dense FD: Γ_ren = E − Λ Σ√g, Λ = −e_y^{dense} (バッチ点から純算術)
        let mut worst = 0.0f64;
        let mut detail = String::new();
        let stat_signs = [1.0f64, -1.0, 0.5, -0.5, 0.25, -0.25];
        for (mi, &m) in ms.iter().enumerate() {
            let lam = -e_dense[mi][1];
            for (di, _dir) in [1usize, 2].iter().enumerate() {
                let mut g6 = [0.0f64; 6];
                for (p, s) in stat_signs.iter().enumerate() {
                    let ev = evals[i_stat[mi][di][p]];
                    g6[p] = ev - lam * vold * (1.0 + s * eps).sqrt();
                }
                let dg = d1_from(g6);
                worst = worst.max(dg.abs());
                if di == 0 {
                    detail = format!("|dΓ/dε|(uni-yy, m={:.1}) = {:.1e}", m, dg.abs());
                }
            }
        }
        check(
            "[S3a] 背景停留: Λ = −e_y で dΓ_ren/dε = 0 (uni-yy/zz, dense FD 3 段 Richardson)",
            worst < 1e-6,
            format!("max = {:.1e} ({})", worst, detail),
        );
        // [S3b] x 残差の N 減少 (twist の有限サイズ圧力異方性)
        for &n in &[16usize, 64] {
            for (mi, &m) in ms.iter().enumerate() {
                let eb = tadpole_block(n, m, nthreads);
                e_block.insert((n, mi), eb);
            }
        }
        let dx = |n: usize, mi: usize| -> f64 {
            let e = e_block[&(n, mi)];
            (e[0] / e[1] - 1.0).abs()
        };
        println!(
            "    [S3b 表] Δ_x(N) = |e_x/e_y − 1|: m=0:   {:.2e} → {:.2e} → {:.2e} → {:.2e} (N=8,16,32,64)",
            dx(8, 0),
            dx(16, 0),
            dx(32, 0),
            dx(64, 0)
        );
        println!(
            "                                     m=0.5: {:.2e} → {:.2e} → {:.2e} → {:.2e}",
            dx(8, 1),
            dx(16, 1),
            dx(32, 1),
            dx(64, 1)
        );
        check(
            "[S3b] 停留残差 Δ_x(N) が単調減少 (m=0.5; twist 圧力異方性は有限サイズ)",
            dx(8, 1) > dx(16, 1) && dx(16, 1) > dx(32, 1) && dx(32, 1) > dx(64, 1),
            format!("{:.2e} → {:.2e} → {:.2e} → {:.2e}", dx(8, 1), dx(16, 1), dx(32, 1), dx(64, 1)),
        );
    }

    // ---- [S4] 接触項の q 非依存の器械認証: ⟨T_i(2qŷ)⟩ = 0 ----
    {
        let mut worst = 0.0f64;
        for (mi, _m) in ms.iter().enumerate() {
            for j in [1usize, 2, 3] {
                let q2y = 2.0 * PI * (2 * j) as f64 / nd as f64;
                for which in [1usize, 2, 3] {
                    // ⟨T(2q)⟩ = Σ_b 2·(t_b e^{i 2q·u}) G_ij — 頂点行列要素との縮約
                    let (re, im) = dense_vertex(nd, q2y, which);
                    let (mut tr_re, mut tr_im) = (0.0f64, 0.0f64);
                    for b in bonds(nd) {
                        if b.dir + 1 != which {
                            continue;
                        }
                        let g = gs[mi][b.j + b.i * nsd];
                        tr_re += 2.0 * re[b.j + b.i * nsd] * g;
                        tr_im += 2.0 * im[b.j + b.i * nsd] * g;
                    }
                    worst = worst.max((tr_re / vold).abs()).max((tr_im / vold).abs());
                }
            }
        }
        check(
            "[S4] y 並進不変性: |⟨T_i(2qŷ)⟩|/V = 0 (接触項 ⟨S(q)⟩ = (3/8)e_i V の厳密性)",
            worst < 1e-12,
            format!("max = {:.1e}", worst),
        );
    }

    // ================= block 走査 (N ∈ {32, 64}) と完全核 =================
    let ns_list = [32usize, 64];
    let jmax = 6usize;
    let js_all: Vec<usize> = (0..=jmax).collect();
    let mut tab: Vec<Vec<Vec<[[f64; 3]; 3]>>> = Vec::new(); // [ni][mi][j]
    for &n in &ns_list {
        let mut per_m = Vec::new();
        for &m in &ms {
            let got = chi_scan_matrix(n, m, &js_all, nthreads, false);
            println!(
                "    [走査] N={} m={:.1} 完了 ({} s) — χ_xx(0) = {:.6}, χ_yy(0) = {:.6}, χ_xy(0) = {:.6}",
                n,
                m,
                t0.elapsed().as_secs(),
                got[0][0][0],
                got[0][1][1],
                got[0][0][1]
            );
            per_m.push(got);
        }
        tab.push(per_m);
    }

    // チャネル射影 (proofs/Projector.lean の ŷ 辞書): D = (xx−zz)/√2, S = (xx+zz)/√2, L = yy
    let chan_d = |c: &[[f64; 3]; 3]| 0.5 * (c[0][0] + c[2][2]) - c[0][2];
    let chan_s = |c: &[[f64; 3]; 3]| 0.5 * (c[0][0] + c[2][2]) + c[0][2];
    let chan_l = |c: &[[f64; 3]; 3]| c[1][1];
    let chan_sl = |c: &[[f64; 3]; 3]| (c[0][1] + c[2][1]) / (2.0f64).sqrt();
    let chan_ds = |c: &[[f64; 3]; 3]| 0.5 * (c[0][0] - c[2][2]);
    let chan_dl = |c: &[[f64; 3]; 3]| (c[0][1] - c[2][1]) / (2.0f64).sqrt();

    // 完全核 k̂_AB(q) = (3/4)e_A δ_AB − ¼χ_AB(q) − (Λ/4)(J−2I)_AB (Λ = −e_y)
    let khat = |ni: usize, mi: usize, j: usize| -> [[f64; 3]; 3] {
        let e = e_block[&(ns_list[ni], mi)];
        let lam = -e[1];
        let chi = &tab[ni][mi][j];
        let mut k = [[0.0f64; 3]; 3];
        for a in 0..3 {
            for b in 0..3 {
                let jm2i = if a == b { -1.0 } else { 1.0 };
                k[a][b] = -0.25 * chi[a][b] - 0.25 * lam * jm2i;
                if a == b {
                    k[a][b] += 0.75 * e[a];
                }
            }
        }
        k
    };

    // 核の表 (N=64)
    println!("\n    [完全核 k̂ チャネル表 (N=64)] j | m | K_DD | K_SS | K_LL | K_SL | K_DS | K_DL");
    for mi in 0..ms.len() {
        for j in [0usize, 1, 2, 4, 6] {
            let k = khat(1, mi, j);
            println!(
                "      j={} m={:.1}:  {:+.6}  {:+.6}  {:+.6}  {:+.6}  {:+.2e}  {:+.2e}",
                j,
                ms[mi],
                chan_d(&k),
                chan_s(&k),
                chan_l(&k),
                chan_sl(&k),
                chan_ds(&k),
                chan_dl(&k)
            );
        }
    }
    {
        let k0 = khat(1, 0, 0);
        let e = e_block[&(64, 0)];
        println!(
            "    [縦核の解剖 m=0] k̂_LL(0) = ½e_y − ¼χ_yy(0) = {:+.6} (e_y = {:+.6}, χ_yy(0) = {:+.6})",
            chan_l(&k0),
            e[1],
            tab[1][0][0][1][1]
        );
        println!("      — 連続極限の diffeo 不変性はこの量の恒等消滅を要求する (q⁰ 繰り込み条件)");
    }

    // ---- c₁ フィット (χ 単位; k̂ の q² 係数は ×(−1/4)) ----
    // c1[ni][mi][ch]: ch 0 = D, 1 = S, 2 = L(yy), 3 = SL, 4 = xx(ŷ), 5 = zz(ŷ)
    let chans: [(&str, &dyn Fn(&[[f64; 3]; 3]) -> f64); 6] = [
        ("D (spin-2 plus)", &chan_d),
        ("S (P0s 横トレース)", &chan_s),
        ("L (yy 縦 = 純ゲージ)", &chan_l),
        ("SL (P0s×P0w 混合)", &chan_sl),
        ("xx(ŷ)", &|c| c[0][0]),
        ("zz(ŷ)", &|c| c[2][2]),
    ];
    let mut c1 = vec![vec![[(0.0f64, 0.0f64); 6]; ms.len()]; ns_list.len()];
    println!("\n    [c₁ 表 (χ 単位)] N | m | D | S | L(yy) | SL");
    for ni in 0..ns_list.len() {
        for mi in 0..ms.len() {
            for (ci, (_, f)) in chans.iter().enumerate() {
                c1[ni][mi][ci] = fit_c1(ns_list[ni], ms[mi], &|j: usize| f(&tab[ni][mi][j]));
            }
            println!(
                "      N={} m={:.1}:  {:+.5}(±{:.5})  {:+.5}(±{:.5})  {:+.5}(±{:.5})  {:+.5}(±{:.5})",
                ns_list[ni],
                ms[mi],
                c1[ni][mi][0].0,
                c1[ni][mi][0].1,
                c1[ni][mi][1].0,
                c1[ni][mi][1].1,
                c1[ni][mi][2].0,
                c1[ni][mi][2].1,
                c1[ni][mi][3].0,
                c1[ni][mi][3].1
            );
        }
    }

    // ---- [S5] 回帰: c₁^χ[D] = v26.5 / c₁^χ[xx](ŷ) = v26.4 ----
    {
        let mut worst = 0.0f64;
        for (mi, &(_m, r)) in REF265_C1D.iter().enumerate() {
            worst = worst.max((c1[1][mi][0].0 - r).abs());
        }
        for (mi, &(_m, r)) in REF264_C1XX_QY.iter().enumerate() {
            worst = worst
                .max((c1[1][mi][4].0 - r).abs())
                .max((c1[1][mi][5].0 - r).abs());
        }
        check(
            "[S5] 完全核成分の回帰: c₁^χ[D] = v26.5 公表値 / c₁^χ[xx,zz](ŷ) = v26.4 横値 (±0.0004)",
            worst < 4e-4,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [S6] 縦 Ward 汚染 (主結果) ----
    {
        let mut rs = [0.0f64; 2];
        let mut resolved = true;
        for mi in 0..ms.len() {
            let (cl, sl) = c1[1][mi][2];
            let (cd, _sd) = c1[1][mi][0];
            rs[mi] = (cl / cd).abs();
            if cl.abs() <= 3.0 * sl.abs() {
                resolved = false;
            }
        }
        let in_band = |r: f64| (1.0 / 3.0..=3.0).contains(&r);
        let branch = if in_band(rs[0]) && in_band(rs[1]) {
            "a: 縦汚染は spin-2 と同桁 — bare c₁ は非共変 regulator 汚染に支配される"
        } else if rs[0] < 1.0 / 3.0 && rs[1] < 1.0 / 3.0 {
            "b: 近似的共変性 (縦 ≪ spin-2)"
        } else {
            "c: 縦汚染が支配的 (R > 3)"
        };
        println!(
            "    [S6 branch] R(m) = |c₁^χ[L]|/|c₁^χ[D]| = {:.3} (m=0) / {:.3} (m=0.5) ⇒ branch {}",
            rs[0], rs[1], branch
        );
        println!(
            "      縦チャネル (yy@qŷ) は純ゲージ (Projector.lean: yhat_longitudinal_is_gauge) —"
        );
        println!(
            "      連続極限では k̂ の縦列は恒等的に消えるべき量。その q² 係数 c₁^χ[L] = {:+.5} (m=0.5)",
            c1[1][1][2].0
        );
        println!("      は繰り込み条件を課すべき純汚染であり、spin-2 の c₁ と同じ桁で走る。");
        check(
            "[S6] 縦 Ward 汚染比 R の分解能 (|c₁^χ[L]| > 3×窓系統) と branch 記録",
            resolved,
            format!("R = {:.3}/{:.3}, branch 判定は上の欄", rs[0], rs[1]),
        );
    }

    // ---- [S7] スカラー混合 K_SL (新測定) ----
    {
        let (c0, s0) = c1[1][0][3];
        let (c5, s5) = c1[1][1][3];
        let (cd, _) = c1[1][1][0];
        println!(
            "    [S7] c₁^χ[SL] = {:+.5}(±{:.5}) [m=0] → {:+.5}(±{:.5}) [m=0.5] — 質量走行 {:.1}%",
            c0,
            s0,
            c5,
            s5,
            100.0 * (c5 / c0 - 1.0).abs()
        );
        println!(
            "      χ_SL(0) = {:+.6} [m=0] / {:+.6} [m=0.5] (P0s×P0w 遷移核の静的値)",
            chan_sl(&tab[1][0][0]),
            chan_sl(&tab[1][1][0])
        );
        let resolved = c5.abs() > 3.0 * s5.abs() || 3.0 * s5.abs() < 0.3 * cd.abs();
        check(
            "[S7] スカラー混合 c₁^χ[SL] の測定が分解能を持つ (非零解決 or 上界 30% |c₁[D]|)",
            resolved,
            format!("|c₁[SL]| = {:.5}, 3σ = {:.5}", c5.abs(), 3.0 * s5.abs()),
        );
    }

    // ---- [S8] uniform 連続性 ----
    {
        let nb = 32usize;
        let volb = (nb * nb * nb) as f64;
        let mut worst_rel = 0.0f64;
        for (mi, &m) in ms.iter().enumerate() {
            let e = e_block[&(nb, mi)];
            let chi0 = &tab[0][mi][0];
            // 対角 xx / yy
            for dir in [0usize, 1] {
                let fd = d2_richardson(
                    &|eps: f64| {
                        let mut s = [1.0f64; 3];
                        s[dir] = 1.0 / (1.0 + eps).sqrt();
                        e0_block_scaled(nb, m, s, nthreads)
                    },
                    0.02,
                );
                let asm = volb * (0.75 * e[dir] - 0.25 * chi0[dir][dir]);
                worst_rel = worst_rel.max((fd / asm - 1.0).abs());
            }
            // 交差 xx-yy
            let fd = d2_cross_richardson(
                &|ea: f64, eb: f64| {
                    let s = [
                        1.0 / (1.0 + ea).sqrt(),
                        1.0 / (1.0 + eb).sqrt(),
                        1.0,
                    ];
                    e0_block_scaled(nb, m, s, nthreads)
                },
                0.02,
            );
            let asm = volb * (-0.25 * chi0[0][1]);
            worst_rel = worst_rel.max((fd / asm - 1.0).abs());
        }
        check(
            "[S8a] uniform 連続性: block 一様変形 FD = V[(3/4)e δ − ¼χ(0)] (N=32, xx/yy/交差)",
            worst_rel < 1e-6,
            format!("max 相対 = {:.1e} ({} s)", worst_rel, t0.elapsed().as_secs()),
        );
        // q₁ → 0 連続性 (m=0.5 ゲート / m=0 報告)
        let scale = {
            let k0 = khat(0, 1, 0);
            let mut s = 0.0f64;
            for a in 0..3 {
                for b in 0..3 {
                    s = s.max(k0[a][b].abs());
                }
            }
            s
        };
        let mut worst_gap = [0.0f64; 2];
        for mi in 0..ms.len() {
            let k0 = khat(0, mi, 0);
            let k1 = khat(0, mi, 1);
            for a in 0..3 {
                for b in 0..3 {
                    worst_gap[mi] = worst_gap[mi].max((k1[a][b] - k0[a][b]).abs());
                }
            }
        }
        println!(
            "    [S8b 表] max|k̂(q₁) − k̂(0)| (N=32): m=0: {:.2e} / m=0.5: {:.2e} (scale = {:.3})",
            worst_gap[0], worst_gap[1], scale
        );
        check(
            "[S8b] k̂ の q → 0 連続性 (m=0.5: ≤ 5% scale; m=0 は報告)",
            worst_gap[1] < 0.05 * scale,
            format!("{:.2e} < {:.2e}", worst_gap[1], 0.05 * scale),
        );
    }

    // ---- [S9] 変異検出 (破壊層) ----
    {
        // (i) 接触項落とし → S0 の uni-yy が崩れる (S0 の FD を再利用)
        let (mi, m) = (0usize, 0.0f64);
        let (w, v) = &dense_eig[mi];
        let md = Mode {
            dir: 1,
            kind: 0,
            j: 0,
            eps: 0.0,
        };
        let ov = mode_vertex(nd, &md);
        let chi_tt = chi_cross_real(w, v, nsd, &ov, &ov);
        let asm_mut = -0.25 * chi_tt; // ⟨S⟩ を落とした変異組立
        check(
            "[S9-i] 変異: 接触項 (3/4)⟨S⟩ 落とし → S0 が検出",
            (fd_uni_yy_m0 - asm_mut).abs() > 1e-3,
            format!("逸脱 {:.3} > 1e-3", (fd_uni_yy_m0 - asm_mut).abs()),
        );
        // (ii) Λ 落とし → 停留が破れる (S3a の点から Λ なしの dE/dε)
        let dg = d1_from(pick6(i_stat[0][0]));
        check(
            "[S9-ii] 変異: Λ counterterm 落とし → S3a が検出",
            dg.abs() > 1e-3,
            format!("|dE/dε| = {:.3} > 1e-3 (S3a の減算後値と対照)", dg.abs()),
        );
        // (iii) T_zz 折返しスワップ落とし → S1 (dense 照合) が検出
        let qy = 2.0 * PI * 2.0 / nd as f64;
        let ops = [
            dense_vertex(nd, qy, 1),
            dense_vertex(nd, qy, 2),
            dense_vertex(nd, qy, 3),
        ];
        let cd = chi_dense_matrix(w, v, nsd, &ops);
        let bad = chi_scan_matrix(nd, m, &[2], nthreads, true);
        let dev = (bad[0][0][2] - cd[0][2] / vold).abs();
        check(
            "[S9-iii] 変異: T_zz 折返しスワップ落とし → S1 が検出",
            dev > 1e-4,
            format!("逸脱 {:.2e} > 1e-4", dev),
        );
    }

    // ---- artifact ----
    let e64 = [e_block[&(64, 0)], e_block[&(64, 1)]];
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.6".into())),
        ("kind".into(), Json::Str("grav_vacuum_pol_kernel".into())),
        (
            "spec".into(),
            Json::Str("paper/grav-vacuum-polarization-spec.md (157ca53)".into()),
        ),
        (
            "tadpole_n64".into(),
            Json::Arr(
                (0..2)
                    .map(|mi| {
                        Json::Obj(vec![
                            ("m".into(), Json::Num(ms[mi])),
                            ("e_x".into(), Json::Num(e64[mi][0])),
                            ("e_y".into(), Json::Num(e64[mi][1])),
                            ("e_z".into(), Json::Num(e64[mi][2])),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "c1_chi_units_n64".into(),
            Json::Arr(
                (0..2)
                    .map(|mi| {
                        Json::Obj(vec![
                            ("m".into(), Json::Num(ms[mi])),
                            ("D_spin2".into(), Json::Num(c1[1][mi][0].0)),
                            ("S_p0s".into(), Json::Num(c1[1][mi][1].0)),
                            ("L_longitudinal".into(), Json::Num(c1[1][mi][2].0)),
                            ("SL_mixing".into(), Json::Num(c1[1][mi][3].0)),
                            (
                                "R_long_over_spin2".into(),
                                Json::Num((c1[1][mi][2].0 / c1[1][mi][0].0).abs()),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "khat_ll_q0_n64".into(),
            Json::Arr(
                (0..2)
                    .map(|mi| Json::Num(chan_l(&khat(1, mi, 0))))
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v266_vacuum_pol.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **完全核の器械認証が成立、縦 Ward 汚染比 R が主結果** — branch は [S6] の欄が一次ソース。解釈は docs/uft-v26.6.md へ"
        } else {
            "FAIL — 分岐 (b)/(c) は各検査の欄を一次ソースとする"
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
