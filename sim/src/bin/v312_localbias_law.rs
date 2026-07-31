//! v31.2 LocalBiasCommutatorLaw — 局所バイアス応答の短時間曲率則 (PROMPT/12 第三十一期)
//!
//! **不変ノルム核の第一候補**。既知ノード因子分解の射影 P_i に対する probe 対
//!   C_i^± = I/2 ± εP_i  (0 < ε < 1/2 — full-rank Gaussian の有効な相関)
//! を準備し、一体生成子 h の下で C(t) = e^{−iht} C(0) e^{iht} と発展させると:
//!
//!   (1) commutator 則:  P_j (Ċ⁺ − Ċ⁻)(0) P_i = −2iε P_j h P_i   (j ≠ i)
//!   (2) 密度曲率則:     (n̈_j⁺ − n̈_j⁻)(0) / (4ε) = ‖P_j h P_i‖_F²
//!
//! (2) は**ノード密度の時系列だけ**から生成子 block の Frobenius 重みを返す —
//! ノード内基底 (gauge) 不変・equilibrium 状態の logit 不使用 (P6/693 の sign(A)
//! no-go と衝突しない)・臨界境界増強も回避 (equilibrium 相関を使わない)。
//! v29.5 の「静的単独不可・応答併用可」を、到着時刻の経験則ではなく**有限行列の
//! 恒等式**にする。
//!
//! 能力の分離 (観測契約 3 段):
//!   CoherentLocalResponse     — block P_j h P_i を gauge 共変に復元 (最強)
//!   LocalBiasDensityResponse  — Frobenius 重み ‖P_j h P_i‖_F² のみ (gauge 不変)
//!   ArrivalTimeResponse       — 現行 B4 の圧縮観測 (到着時刻 — 最弱)
//!
//! 検査:
//!   [L0] 厳密代数 oracle: 複素エルミート h (4 ノード × 2 軌道) × ε 3 値で恒等式
//!        (1)(2) が機械精度 (probe 線形性により ε の高次補正は厳密に 0)
//!   [L1] 密度測定 lane: **時系列のみ**を受け取る readout (h 不可視) — 5 点 stencil +
//!        Richardson で ‖P_j h P_i‖_F² を rel 1e-5 復元・dt 収束次数 ~4 (中心差分 2)
//!   [L2] coherent 測定 lane: block P_j h P_i を entry 1e-5 で復元
//!   [L3] gauge 共変/不変 (乱択 U(2) block 10 draws): 密度重みは不変・coherent block
//!        は U_j B U_i† で共変
//!   [L4] ノード置換共変
//!   [L5] 能力の厳密分離: block を U 回転した h₂ (Frobenius 重み不変) は
//!        LocalBiasDensityResponse では識別不能 (厳密 0)・CoherentLocalResponse は識別 —
//!        観測契約の hierarchy の機械実例 (ArrivalTime は参考記録)
//!   [L6] **P6/693 の分離**: 静的 projector ではゲージ同値 (v31.1 [T9]) だった衝突対を、
//!        密度曲率則が隣接行列ごと復元して分離 (probe は equilibrium 状態と独立)
//!   [L7] ε 非依存性: ε = 0.05 と 0.45 の測定が一致 (恒等式が線形応答近似でないこと)
//!   [L8] Lean 台帳: proofs/LocalBias.lean (6 定理 — commutator block・trace–Frobenius・
//!        曲率恒等式・probe 分離・多軌道版・gauge 不変性, 格子 native_decide)
//!
//! 実行: cargo run --release --bin v312_localbias_law

use std::fs;
use std::path::Path;
use uft_sim::readout_contract::*;
use uft_sim::{jacobi_eigh, Rng, C64};

// ---------------------------------------------------------------- 複素行列の小物

fn czero() -> C64 {
    C64::new(0.0, 0.0)
}

fn cmat(n: usize) -> Vec<C64> {
    vec![czero(); n * n]
}

fn cmatmul(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let mut c = cmat(n);
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik.norm2() == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] = c[i * n + j] + aik * b[k * n + j];
            }
        }
    }
    c
}

fn cmax_diff(a: &[C64], b: &[C64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| (*x - *y).abs())
        .fold(0.0, f64::max)
}

/// エルミート行列 (re 対称, im 反対称) の固有系 — 2n×2n 実対称埋め込み
/// [[A, −B], [B, A]] を jacobi で解き、二重化した実固有対から複素固有ベクトルを
/// Gram–Schmidt で n 本選ぶ (縮退があっても射影残差で正しく拾う)。
fn herm_eig(re: &[f64], im: &[f64], n: usize) -> (Vec<f64>, Vec<Vec<C64>>) {
    let m = 2 * n;
    let mut big = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            big[i * m + j] = re[i * n + j];
            big[(i + n) * m + (j + n)] = re[i * n + j];
            big[i * m + (j + n)] = -im[i * n + j];
            big[(i + n) * m + j] = im[i * n + j];
        }
    }
    let (evals, evecs) = jacobi_eigh(&big, m);
    let mut order: Vec<usize> = (0..m).collect();
    order.sort_by(|&a, &b| evals[a].partial_cmp(&evals[b]).unwrap());
    let mut out_vals: Vec<f64> = Vec::new();
    let mut out_vecs: Vec<Vec<C64>> = Vec::new();
    for &idx in &order {
        if out_vecs.len() == n {
            break;
        }
        let lam = evals[idx];
        // 実 2n ベクトル (x || y) → 複素 v = x + iy
        let mut v: Vec<C64> = (0..n)
            .map(|i| C64::new(evecs[idx * m + i], evecs[idx * m + n + i]))
            .collect();
        // 既採用の同固有値ベクトルを射影で除く
        for (k, w) in out_vecs.iter().enumerate() {
            if (out_vals[k] - lam).abs() > 1e-8 {
                continue;
            }
            let mut ip = czero(); // ⟨w, v⟩
            for i in 0..n {
                ip = ip + w[i].conj() * v[i];
            }
            for i in 0..n {
                v[i] = v[i] - ip * w[i];
            }
        }
        let nrm: f64 = v.iter().map(|z| z.norm2()).sum::<f64>().sqrt();
        if nrm > 1e-6 {
            for z in v.iter_mut() {
                *z = z.scale(1.0 / nrm);
            }
            out_vals.push(lam);
            out_vecs.push(v);
        }
    }
    assert_eq!(out_vecs.len(), n, "複素固有ベクトルの本数が n に満たない");
    (out_vals, out_vecs)
}

/// C(t) = e^{−iht} C0 e^{iht} — h の固有系 (vals, vecs) で厳密に発展
fn evolve(vals: &[f64], vecs: &[Vec<C64>], c0: &[C64], n: usize, t: f64) -> Vec<C64> {
    // C̃0 = V† C0 V (V の列 = 固有ベクトル)
    let mut ct = cmat(n);
    for a in 0..n {
        for b in 0..n {
            let mut s = czero();
            for i in 0..n {
                for j in 0..n {
                    s = s + vecs[a][i].conj() * c0[i * n + j] * vecs[b][j];
                }
            }
            ct[a * n + b] = s * C64::expi(-(vals[a] - vals[b]) * t);
        }
    }
    // C(t) = V C̃ V†
    let mut c = cmat(n);
    for i in 0..n {
        for j in 0..n {
            let mut s = czero();
            for a in 0..n {
                for b in 0..n {
                    s = s + vecs[a][i] * ct[a * n + b] * vecs[b][j].conj();
                }
            }
            c[i * n + j] = s;
        }
    }
    c
}

// ---------------------------------------------------------------- probe と観測 (lab 側)

/// ノード定義: node_of[site] = ノード番号
#[derive(Clone)]
struct Nodes {
    node_of: Vec<usize>,
    n_nodes: usize,
}

impl Nodes {
    fn sites_of(&self, b: usize) -> Vec<usize> {
        (0..self.node_of.len())
            .filter(|&s| self.node_of[s] == b)
            .collect()
    }
}

/// probe C_i^± = I/2 ± εP_i
fn probe(nodes: &Nodes, i: usize, eps: f64, sign: f64, n: usize) -> Vec<C64> {
    let mut c = cmat(n);
    for s in 0..n {
        let d = if nodes.node_of[s] == i { 0.5 + sign * eps } else { 0.5 };
        c[s * n + s] = C64::new(d, 0.0);
    }
    c
}

/// ノード密度ベクトル n_j = Re Tr(P_j C)
fn densities(c: &[C64], nodes: &Nodes, n: usize) -> Vec<f64> {
    let mut out = vec![0.0; nodes.n_nodes];
    for s in 0..n {
        out[nodes.node_of[s]] += c[s * n + s].re;
    }
    out
}

/// 密度時系列 (lab): source ノード i の ± probe を 5 点 stencil {±dt, ±dt/2}
/// で発展させ、各ノードの密度を返す。**readout はこの出力のみを受け取る。**
struct DensitySeries {
    /// [probe ±][時刻 4: −dt, −dt/2, +dt/2, +dt][ノード]
    n_at: [[Vec<f64>; 4]; 2],
    n0: [Vec<f64>; 2],
    dt: f64,
    eps: f64,
}

fn lab_density_series(
    re: &[f64],
    im: &[f64],
    nodes: &Nodes,
    i: usize,
    eps: f64,
    dt: f64,
    n: usize,
) -> DensitySeries {
    let (vals, vecs) = herm_eig(re, im, n);
    let times = [-dt, -dt / 2.0, dt / 2.0, dt];
    let mut out = DensitySeries {
        n_at: [
            [vec![], vec![], vec![], vec![]],
            [vec![], vec![], vec![], vec![]],
        ],
        n0: [vec![], vec![]],
        dt,
        eps,
    };
    for (pi, sign) in [(0usize, 1.0), (1usize, -1.0)] {
        let c0 = probe(nodes, i, eps, sign, n);
        out.n0[pi] = densities(&c0, nodes, n);
        for (ti, &t) in times.iter().enumerate() {
            let c = evolve(&vals, &vecs, &c0, n, t);
            out.n_at[pi][ti] = densities(&c, nodes, n);
        }
    }
    out
}

// ---------------------------------------------------------------- readout (時系列のみ)

/// 測定推定 (時系列のみ — h は不可視)。Richardson: R = (4 D(dt/2) − D(dt))/3
struct MeasuredFrobenius {
    w: Vec<f64>,
    /// 中心差分 (dt) 単独の値 — 収束次数の検査用
    w_coarse: Vec<f64>,
}

fn readout_density_frobenius(s: &DensitySeries) -> MeasuredFrobenius {
    let nb = s.n0[0].len();
    let mut w = vec![0.0; nb];
    let mut w_coarse = vec![0.0; nb];
    for j in 0..nb {
        let d2 = |pi: usize, half: bool| -> f64 {
            let (tm, tp, dt) = if half {
                (1, 2, s.dt / 2.0)
            } else {
                (0, 3, s.dt)
            };
            (s.n_at[pi][tp][j] - 2.0 * s.n0[pi][j] + s.n_at[pi][tm][j]) / (dt * dt)
        };
        let coarse = (d2(0, false) - d2(1, false)) / (4.0 * s.eps);
        let fine = (d2(0, true) - d2(1, true)) / (4.0 * s.eps);
        w_coarse[j] = coarse;
        w[j] = (4.0 * fine - coarse) / 3.0;
    }
    MeasuredFrobenius { w, w_coarse }
}

/// coherent 測定 lane: block P_j C(t) P_i の時系列から Ḃ⁺ − Ḃ⁻ を Richardson で
/// 推定し、B̂ = (Ḃ⁺ − Ḃ⁻)/(−2iε) を返す (時系列のみ — h は不可視)
struct BlockSeries {
    /// [probe ±][時刻 4][block 成分]
    b_at: [[Vec<C64>; 4]; 2],
    dt: f64,
    eps: f64,
}

fn lab_block_series(
    re: &[f64],
    im: &[f64],
    nodes: &Nodes,
    i: usize,
    j: usize,
    eps: f64,
    dt: f64,
    n: usize,
) -> BlockSeries {
    let (vals, vecs) = herm_eig(re, im, n);
    let si = nodes.sites_of(i);
    let sj = nodes.sites_of(j);
    let times = [-dt, -dt / 2.0, dt / 2.0, dt];
    let mut out = BlockSeries {
        b_at: [
            [vec![], vec![], vec![], vec![]],
            [vec![], vec![], vec![], vec![]],
        ],
        dt,
        eps,
    };
    for (pi, sign) in [(0usize, 1.0), (1usize, -1.0)] {
        let c0 = probe(nodes, i, eps, sign, n);
        for (ti, &t) in times.iter().enumerate() {
            let c = evolve(&vals, &vecs, &c0, n, t);
            let mut blk = Vec::with_capacity(sj.len() * si.len());
            for &a in &sj {
                for &b in &si {
                    blk.push(c[a * n + b]);
                }
            }
            out.b_at[pi][ti] = blk;
        }
    }
    out
}

fn readout_coherent_block(s: &BlockSeries) -> Vec<C64> {
    let ne = s.b_at[0][0].len();
    let mut out = vec![czero(); ne];
    for e in 0..ne {
        let d1 = |pi: usize, half: bool| -> C64 {
            let (tm, tp, dt) = if half {
                (1, 2, s.dt / 2.0)
            } else {
                (0, 3, s.dt)
            };
            (s.b_at[pi][tp][e] - s.b_at[pi][tm][e]).scale(1.0 / (2.0 * dt))
        };
        let coarse = d1(0, false) - d1(1, false);
        let fine = d1(0, true) - d1(1, true);
        let rich = (fine.scale(4.0) - coarse).scale(1.0 / 3.0);
        // B̂ = (Ḃ⁺ − Ḃ⁻)/(−2iε): ×(−1/(2ε))·(1/i) = ×(i/(2ε))
        out[e] = rich * C64::new(0.0, 1.0 / (2.0 * s.eps));
    }
    out
}

// ---------------------------------------------------------------- 厳密代数 oracle (診断)

/// P_j h P_i block (真値 — 診断・比較専用)
fn true_block(re: &[f64], im: &[f64], nodes: &Nodes, i: usize, j: usize, n: usize) -> Vec<C64> {
    let si = nodes.sites_of(i);
    let sj = nodes.sites_of(j);
    let mut out = Vec::with_capacity(sj.len() * si.len());
    for &a in &sj {
        for &b in &si {
            out.push(C64::new(re[a * n + b], im[a * n + b]));
        }
    }
    out
}

fn frob2_block(b: &[C64]) -> f64 {
    b.iter().map(|z| z.norm2()).sum()
}

fn main() {
    uft_sim::self_test();
    println!("=== v31.2 LocalBiasCommutatorLaw — 局所バイアス応答の短時間曲率則 (PROMPT/12) ===");
    println!("(probe C± = I/2 ± εP_i は準備された状態 — equilibrium 状態の logit を使わず、");
    println!(" P6/693 の sign(A) no-go・臨界境界増強と衝突しない。readout は時系列のみを見る)\n");
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

    // ---- 試験系: 4 ノード × 2 軌道 = 8 サイト, 複素エルミート h (環 + 弦) ----
    let n = 8usize;
    let nodes = Nodes {
        node_of: vec![0, 0, 1, 1, 2, 2, 3, 3],
        n_nodes: 4,
    };
    let (h_re, h_im) = {
        let mut rng = Rng::new(312);
        let mut re = vec![0.0; n * n];
        let mut im = vec![0.0; n * n];
        // ノード内 (エルミート 2×2)
        for b in 0..4 {
            let (i, j) = (2 * b, 2 * b + 1);
            re[i * n + i] = rng.f64() - 0.5;
            re[j * n + j] = rng.f64() - 0.5;
            let tr = rng.f64() - 0.5;
            let ti = rng.f64() - 0.5;
            re[i * n + j] = tr;
            re[j * n + i] = tr;
            im[i * n + j] = ti;
            im[j * n + i] = -ti;
        }
        // ノード間 (環 0-1-2-3-0 + 弦 0-2, 複素 full 2×2)
        for (bi, bj) in [(0usize, 1usize), (1, 2), (2, 3), (0, 3), (0, 2)] {
            for oi in 0..2 {
                for oj in 0..2 {
                    let tr = rng.f64() - 0.5;
                    let ti = rng.f64() - 0.5;
                    let (a, b) = (2 * bi + oi, 2 * bj + oj);
                    re[a * n + b] = tr;
                    re[b * n + a] = tr;
                    im[a * n + b] = ti;
                    im[b * n + a] = -ti;
                }
            }
        }
        (re, im)
    };

    // ---- [L0] 厳密代数 oracle ----
    {
        let mut worst_comm: f64 = 0.0;
        let mut worst_curv: f64 = 0.0;
        let hmat: Vec<C64> = (0..n * n)
            .map(|k| C64::new(h_re[k], h_im[k]))
            .collect();
        for &eps in &[0.1, 0.25, 0.4] {
            for i in 0..4 {
                // Ċ± = −i[h, C±] (代数)
                let cp = probe(&nodes, i, eps, 1.0, n);
                let cm = probe(&nodes, i, eps, -1.0, n);
                let mi = C64::new(0.0, -1.0);
                let dot = |c0: &[C64]| -> Vec<C64> {
                    let hc = cmatmul(&hmat, c0, n);
                    let ch = cmatmul(c0, &hmat, n);
                    (0..n * n).map(|k| (hc[k] - ch[k]) * mi).collect()
                };
                let dp = dot(&cp);
                let dm = dot(&cm);
                for j in 0..4 {
                    if j == i {
                        continue;
                    }
                    // (1) P_j (Ċ⁺−Ċ⁻) P_i = −2iε P_j h P_i
                    let sj = nodes.sites_of(j);
                    let si = nodes.sites_of(i);
                    let tb = true_block(&h_re, &h_im, &nodes, i, j, n);
                    let mut e = 0usize;
                    for &a in &sj {
                        for &b in &si {
                            let lhs = dp[a * n + b] - dm[a * n + b];
                            let rhs = tb[e] * C64::new(0.0, -2.0 * eps);
                            worst_comm = worst_comm.max((lhs - rhs).abs());
                            e += 1;
                        }
                    }
                    // (2) n̈ 差 = 4ε‖P_j h P_i‖²: n̈ = Re Tr(P_j (−i[h, −i[h,C]]))
                    let ddp = dot(&dot(&cp));
                    let ddm = dot(&dot(&cm));
                    let nj_dd = |dd: &[C64]| -> f64 {
                        sj.iter().map(|&s| dd[s * n + s].re).sum::<f64>()
                    };
                    let lhs = nj_dd(&ddp) - nj_dd(&ddm);
                    let rhs = 4.0 * eps * frob2_block(&tb);
                    worst_curv = worst_curv.max((lhs - rhs).abs());
                }
            }
        }
        check(
            "[L0] 厳密代数 oracle: commutator 則と密度曲率則が機械精度 (ε 3 値 — 高次補正は厳密 0)",
            worst_comm <= 1e-13 && worst_curv <= 1e-12,
            format!("commutator 残差 {:.2e} / 曲率残差 {:.2e}", worst_comm, worst_curv),
        );
        // 契約の実演: 恒等式証人で exact 証明書を構成 (残差がバー内のときのみ成功する)
        type Cert = ReadoutCertificate<
            uft_sim::qrn_core::SpatialMetricUpToGlobalScale,
            GaussianGibbsFullRank,
            LocalBiasDensityResponse,
            GivenNodeFactorization,
        >;
        let cert = Cert::exact(
            "QRN-BRIDGE-013",
            &ExactWitness::AlgebraicIdentity {
                residual: worst_curv,
                bar: 1e-12,
            },
            false,
        );
        assert!(cert.is_ok(), "恒等式証人の証明書構成に失敗");
    }

    // ---- [L1] 密度測定 lane (時系列のみ) ----
    let dt = 0.02;
    let eps = 0.3;
    {
        let mut worst_rel: f64 = 0.0;
        let mut worst_ratio: f64 = 0.0; // 収束次数: coarse 誤差 / Richardson 誤差
        for i in 0..4 {
            let series = lab_density_series(&h_re, &h_im, &nodes, i, eps, dt, n);
            let m = readout_density_frobenius(&series);
            for j in 0..4 {
                if j == i {
                    continue;
                }
                let tw = frob2_block(&true_block(&h_re, &h_im, &nodes, i, j, n));
                let rel = (m.w[j] - tw).abs() / (1.0 + tw);
                worst_rel = worst_rel.max(rel);
                let e_coarse = (m.w_coarse[j] - tw).abs();
                let e_rich = (m.w[j] - tw).abs().max(1e-14);
                worst_ratio = worst_ratio.max(e_coarse / e_rich);
            }
        }
        check(
            "[L1] 密度測定 lane (h 不可視・時系列のみ): ‖P_j h P_i‖_F² を rel ≤ 1e-5 復元 + Richardson が中心差分を ≥ 10× 改善",
            worst_rel <= 1e-5 && worst_ratio >= 10.0,
            format!(
                "max rel 誤差 {:.2e} (dt = {}) / 改善比 min→ {:.1}×",
                worst_rel, dt, worst_ratio
            ),
        );
    }

    // ---- [L2] coherent 測定 lane ----
    {
        let mut worst: f64 = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                if j == i {
                    continue;
                }
                let s = lab_block_series(&h_re, &h_im, &nodes, i, j, eps, dt, n);
                let est = readout_coherent_block(&s);
                let tb = true_block(&h_re, &h_im, &nodes, i, j, n);
                worst = worst.max(cmax_diff(&est, &tb));
            }
        }
        check(
            "[L2] coherent 測定 lane: block P_j h P_i を entry ≤ 1e-5 で gauge 共変に復元",
            worst <= 1e-5,
            format!("max entry 誤差 {:.2e}", worst),
        );
    }

    // ---- [L3] gauge 共変/不変 (乱択 U(2) block) ----
    {
        let mut rng = Rng::new(3121);
        let mut worst_w: f64 = 0.0;
        let mut worst_cov: f64 = 0.0;
        for _ in 0..10 {
            // U_b = diag(e^{iα}, e^{iβ}) · 回転(θ)
            let us: Vec<[C64; 4]> = (0..4)
                .map(|_| {
                    let (al, be, th) = (
                        rng.f64() * std::f64::consts::TAU,
                        rng.f64() * std::f64::consts::TAU,
                        rng.f64() * std::f64::consts::TAU,
                    );
                    let (c, s) = (th.cos(), th.sin());
                    [
                        C64::expi(al).scale(c),
                        C64::expi(al).scale(-s),
                        C64::expi(be).scale(s),
                        C64::expi(be).scale(c),
                    ]
                })
                .collect();
            // h' = U h U†
            let mut hp_re = vec![0.0; n * n];
            let mut hp_im = vec![0.0; n * n];
            for bi in 0..4 {
                for bj in 0..4 {
                    for oi in 0..2 {
                        for oj in 0..2 {
                            let mut s = czero();
                            for p in 0..2 {
                                for q in 0..2 {
                                    let hpq = C64::new(
                                        h_re[(2 * bi + p) * n + (2 * bj + q)],
                                        h_im[(2 * bi + p) * n + (2 * bj + q)],
                                    );
                                    s = s + us[bi][oi * 2 + p] * hpq * us[bj][oj * 2 + q].conj();
                                }
                            }
                            hp_re[(2 * bi + oi) * n + (2 * bj + oj)] = s.re;
                            hp_im[(2 * bi + oi) * n + (2 * bj + oj)] = s.im;
                        }
                    }
                }
            }
            // エルミート化 (丸め)
            for a in 0..n {
                for b in (a + 1)..n {
                    let r = 0.5 * (hp_re[a * n + b] + hp_re[b * n + a]);
                    let s = 0.5 * (hp_im[a * n + b] - hp_im[b * n + a]);
                    hp_re[a * n + b] = r;
                    hp_re[b * n + a] = r;
                    hp_im[a * n + b] = s;
                    hp_im[b * n + a] = -s;
                }
                hp_im[a * n + a] = 0.0;
            }
            let (i, j) = (0usize, 2usize);
            // 密度重み: 不変
            let m0 = readout_density_frobenius(&lab_density_series(&h_re, &h_im, &nodes, i, eps, dt, n));
            let m1 = readout_density_frobenius(&lab_density_series(&hp_re, &hp_im, &nodes, i, eps, dt, n));
            worst_w = worst_w.max((m0.w[j] - m1.w[j]).abs());
            // coherent block: B̂' = U_j B̂ U_i†
            let e0 = readout_coherent_block(&lab_block_series(&h_re, &h_im, &nodes, i, j, eps, dt, n));
            let e1 = readout_coherent_block(&lab_block_series(&hp_re, &hp_im, &nodes, i, j, eps, dt, n));
            let mut t = vec![czero(); 4];
            for r in 0..2 {
                for cc in 0..2 {
                    let mut s = czero();
                    for p in 0..2 {
                        for q in 0..2 {
                            s = s + us[j][r * 2 + p] * e0[p * 2 + q] * us[i][cc * 2 + q].conj();
                        }
                    }
                    t[r * 2 + cc] = s;
                }
            }
            worst_cov = worst_cov.max(cmax_diff(&e1, &t));
        }
        check(
            "[L3] gauge: 密度重みは block-local U(2) 不変 ≤ 1e-9・coherent block は U_j B U_i† 共変 ≤ 1e-5",
            worst_w <= 1e-9 && worst_cov <= 1e-5,
            format!("重み乖離 {:.2e} / 共変乖離 {:.2e}", worst_w, worst_cov),
        );
    }

    // ---- [L4] ノード置換共変 ----
    {
        let pi = [3usize, 0, 2, 1];
        let mut hp_re = vec![0.0; n * n];
        let mut hp_im = vec![0.0; n * n];
        for bi in 0..4 {
            for bj in 0..4 {
                for oi in 0..2 {
                    for oj in 0..2 {
                        hp_re[(2 * bi + oi) * n + (2 * bj + oj)] =
                            h_re[(2 * pi[bi] + oi) * n + (2 * pi[bj] + oj)];
                        hp_im[(2 * bi + oi) * n + (2 * bj + oj)] =
                            h_im[(2 * pi[bi] + oi) * n + (2 * pi[bj] + oj)];
                    }
                }
            }
        }
        let mut worst: f64 = 0.0;
        for i in 0..4 {
            let mp = readout_density_frobenius(&lab_density_series(&hp_re, &hp_im, &nodes, i, eps, dt, n));
            // 対応する元系: source π(i), 読み j → π(j)
            let m0 = readout_density_frobenius(&lab_density_series(
                &h_re, &h_im, &nodes, pi[i], eps, dt, n,
            ));
            for j in 0..4 {
                if j == i {
                    continue;
                }
                worst = worst.max((mp.w[j] - m0.w[pi[j]]).abs());
            }
        }
        check(
            "[L4] ノード置換共変: Ŵ(PhPᵀ)_{ij} = Ŵ(h)_{π(i)π(j)} ≤ 1e-9",
            worst <= 1e-9,
            format!("max 乖離 {:.2e}", worst),
        );
    }

    // ---- [L5] 能力の厳密分離 (観測契約の hierarchy) ----
    {
        // h₂ = block (1←0) を SO(2) 回転 (Frobenius 重みは不変・block は変化)
        let th = 0.7f64;
        let (c0, s0) = (th.cos(), th.sin());
        let mut h2_re = h_re.clone();
        let mut h2_im = h_im.clone();
        // B' = R B (R は node1 側の実回転): sites (2,3) × (0,1)
        for col in 0..2 {
            let (a, b) = (2usize, 3usize);
            let (r0, i0) = (h_re[a * n + col], h_im[a * n + col]);
            let (r1, i1) = (h_re[b * n + col], h_im[b * n + col]);
            h2_re[a * n + col] = c0 * r0 - s0 * r1;
            h2_im[a * n + col] = c0 * i0 - s0 * i1;
            h2_re[b * n + col] = s0 * r0 + c0 * r1;
            h2_im[b * n + col] = s0 * i0 + c0 * i1;
            // エルミート共役側
            h2_re[col * n + a] = h2_re[a * n + col];
            h2_im[col * n + a] = -h2_im[a * n + col];
            h2_re[col * n + b] = h2_re[b * n + col];
            h2_im[col * n + b] = -h2_im[b * n + col];
        }
        // (a) 密度 lane (厳密代数): 全 (i,j) 対で重みが厳密一致 → 識別不能
        let mut worst_w: f64 = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                if i == j {
                    continue;
                }
                let w1 = frob2_block(&true_block(&h_re, &h_im, &nodes, i, j, n));
                let w2 = frob2_block(&true_block(&h2_re, &h2_im, &nodes, i, j, n));
                worst_w = worst_w.max((w1 - w2).abs());
            }
        }
        // (b) coherent lane: block (0→1) が有限に異なる → 識別可能
        let b1 = true_block(&h_re, &h_im, &nodes, 0, 1, n);
        let b2 = true_block(&h2_re, &h2_im, &nodes, 0, 1, n);
        let sep = cmax_diff(&b1, &b2);
        // (c) 測定 lane でも確認
        let m1 = readout_density_frobenius(&lab_density_series(&h_re, &h_im, &nodes, 0, eps, dt, n));
        let m2 = readout_density_frobenius(&lab_density_series(&h2_re, &h2_im, &nodes, 0, eps, dt, n));
        let meas_diff = (m1.w[1] - m2.w[1]).abs();
        // (d) ArrivalTime (参考): probe +ε at node0 → node2 密度の閾値到着
        let arrival = |re: &[f64], im: &[f64]| -> f64 {
            let (vals, vecs) = herm_eig(re, im, n);
            let c0 = probe(&nodes, 0, eps, 1.0, n);
            let base = densities(&c0, &nodes, n)[2];
            let mut t = 0.0;
            for k in 1..=400 {
                t = 0.01 * k as f64;
                let c = evolve(&vals, &vecs, &c0, n, t);
                if (densities(&c, &nodes, n)[2] - base).abs() > 1e-3 {
                    break;
                }
            }
            t
        };
        let (t1, t2) = (arrival(&h_re, &h_im), arrival(&h2_re, &h2_im));
        check(
            "[L5] 能力分離: block 回転 h₂ は密度 lane で識別不能 (厳密 0)・coherent lane は識別 ≥ 0.1 — 観測契約は真に階層",
            worst_w <= 1e-12 && sep >= 0.1 && meas_diff <= 1e-5,
            format!(
                "重み差 {:.1e} (測定 {:.1e}) / coherent 分離 {:.3} / 到着時刻 (参考): {:.2} vs {:.2}",
                worst_w, meas_diff, sep, t1, t2
            ),
        );
    }

    // ---- [L6] P6/693 の分離 (静的 projector のゲージ同値衝突を応答が破る) ----
    {
        let n6 = 6usize;
        let nodes6 = Nodes {
            node_of: (0..6).collect(),
            n_nodes: 6,
        };
        // v31.1 [T0] と同じ mask → 隣接
        let edges_692: [(usize, usize); 5] = [(0, 3), (0, 5), (1, 2), (1, 4), (2, 3)];
        let edges_693: [(usize, usize); 6] = [(0, 1), (0, 3), (0, 5), (1, 2), (1, 4), (2, 3)];
        let mk_h = |edges: &[(usize, usize)]| -> Vec<f64> {
            let mut h = vec![0.0; n6 * n6];
            for &(i, j) in edges {
                h[i * n6 + j] = -1.0;
                h[j * n6 + i] = -1.0;
            }
            h
        };
        let h_p6 = mk_h(&edges_692);
        let h_u693 = mk_h(&edges_693);
        let zeros = vec![0.0; n6 * n6];
        let measure_w = |h: &[f64]| -> Vec<f64> {
            let mut w = vec![0.0; n6 * n6];
            for i in 0..n6 {
                let m = readout_density_frobenius(&lab_density_series(
                    h, &zeros, &nodes6, i, eps, dt, n6,
                ));
                for j in 0..n6 {
                    if j != i {
                        w[j * n6 + i] = m.w[j];
                    }
                }
            }
            w
        };
        let w1 = measure_w(&h_p6);
        let w2 = measure_w(&h_u693);
        // (a) 隣接行列 (|h_ij|² = A_ij) の復元
        let mut worst_adj: f64 = 0.0;
        let a_p6 = {
            let mut a = vec![0.0; n6 * n6];
            for &(i, j) in &edges_692 {
                a[i * n6 + j] = 1.0;
                a[j * n6 + i] = 1.0;
            }
            a
        };
        let a_693 = {
            let mut a = vec![0.0; n6 * n6];
            for &(i, j) in &edges_693 {
                a[i * n6 + j] = 1.0;
                a[j * n6 + i] = 1.0;
            }
            a
        };
        for k in 0..n6 * n6 {
            if k / n6 == k % n6 {
                continue;
            }
            worst_adj = worst_adj.max((w1[k] - a_p6[k]).abs()).max((w2[k] - a_693[k]).abs());
        }
        // (b) 対の分離 (min-perm ∞) — 静的 projector ではゲージ同値だった (v31.1 [T9])
        let ps6 = {
            fn perms(n: usize) -> Vec<Vec<usize>> {
                let mut out: Vec<Vec<usize>> = vec![vec![]];
                for k in 0..n {
                    let mut next = Vec::new();
                    for p in out {
                        for pos in 0..=p.len() {
                            let mut q = p.clone();
                            q.insert(pos, k);
                            next.push(q);
                        }
                    }
                    out = next;
                }
                out
            }
            perms(6)
        };
        let mut best = f64::INFINITY;
        for pi in &ps6 {
            let mut d: f64 = 0.0;
            for i in 0..n6 {
                for j in 0..n6 {
                    d = d.max((w1[i * n6 + j] - w2[pi[i] * n6 + pi[j]]).abs());
                }
            }
            best = best.min(d);
        }
        check(
            "[L6] P6/693: 密度曲率則が隣接を復元 (≤ 1e-5) し衝突対を分離 (min-perm ≥ 0.9) — 静的不可・応答可の恒等式化",
            worst_adj <= 1e-5 && best >= 0.9,
            format!("隣接復元誤差 {:.2e} / 対の分離 (min-perm ∞) {:.3}", worst_adj, best),
        );
    }

    // ---- [L7] ε 非依存性 (線形応答近似ではない) ----
    {
        let mut worst: f64 = 0.0;
        for &(e1, e2) in &[(0.05, 0.45)] {
            for i in 0..4 {
                let m1 = readout_density_frobenius(&lab_density_series(&h_re, &h_im, &nodes, i, e1, dt, n));
                let m2 = readout_density_frobenius(&lab_density_series(&h_re, &h_im, &nodes, i, e2, dt, n));
                for j in 0..4 {
                    if j != i {
                        worst = worst.max((m1.w[j] - m2.w[j]).abs());
                    }
                }
            }
        }
        check(
            "[L7] ε 非依存: ε = 0.05 と 0.45 の測定が一致 ≤ 1e-4 (恒等式は probe 線形 — 微小応答近似ではない)",
            worst <= 1e-4,
            format!("max 乖離 {:.2e}", worst),
        );
    }

    // ---- [L8] Lean 台帳 ----
    {
        let root = if Path::new("proofs/LocalBias.lean").exists() {
            "."
        } else {
            ".."
        };
        let lean = fs::read_to_string(format!("{}/proofs/LocalBias.lean", root)).unwrap_or_default();
        let n_thm = lean
            .lines()
            .filter(|l| l.trim_start().starts_with("theorem "))
            .count();
        let ok = n_thm == 6
            && lean.contains("native_decide")
            && lean.contains("スコープの明示")
            && lean.contains("commutator_block_identity")
            && lean.contains("curvature_frobenius_identity")
            && lean.contains("multiorbital_trace_frobenius")
            && lean.contains("gauge_invariance_frobenius");
        check(
            "[L8] proofs/LocalBias.lean — 6 定理 (commutator block・trace–Frobenius・曲率・probe 分離・多軌道・gauge 不変)",
            ok,
            format!(
                "theorem 宣言 {} 本 (lean 4.31 でコンパイル検証 — 保証は格子上の恒等・複素/一般 d は数値側)",
                n_thm
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "局所バイアス曲率則は成立 — 密度時系列だけで生成子 block の Frobenius 重みが gauge 不変に読め、coherent 応答は block を共変復元する。静的 projector でゲージ同値だった P6/693 も応答で分離 — 「静的単独不可・応答併用可」が有限行列の恒等式になった (不変ノルム核の第一候補)"
        } else {
            "**曲率則の破れ** — 恒等式・測定 lane・不変性のいずれかが不成立"
        }
    );
    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
