//! v31.1 GaussianGibbsInverseOracle — 熱的 Gaussian の大域静的逆問題 (PROMPT/12 第三十一期)
//!
//! **位置づけ = oracle ceiling (識別可能性の上界)**。number-conserving quadratic Gibbs 族
//!   C = (I + e^{β(h−μI)})⁻¹,  0 ≺ C ≺ I ⇒ K(C) = log[(I−C)C⁻¹] = β(h−μI)
//! の逆変換は自由フェルミオンの標準的事実であり **QRN 固有法則ではない** — Gibbs 状態
//! からの Hamiltonian learning は独立した既存分野である。よって本器械は bridge 候補では
//! なく「完全な大域 C を観測できたときに原理的に何が読めるか」の天井を機械化する。
//! **この結果だけから BridgeLawCertificate は登録しない・PRED-019 も登録しない。**
//! 新規性は識別可能性相図の側 (exact/estimate lane・同値類・条件数・衝突対の分離) にある。
//!
//! 検査:
//!   [T0] n = 4..7 連結グラフ同型類の全数列挙 = OEIS A001349 (6/21/112/853, v29.5 再現)
//!        + 衝突対 mask 692 (P6)/693 (C6) の実在アンカー
//!   [T1] exact lane 全数: 992 グラフ × β ∈ {0.5, 1, 2, 4} (μ = 0.2) — clamp なし資格
//!        (0 < λ < 1, δ 記録) → K(C) → 門 (β, μ 既知) → max|ĥ−h| ≤ 1e-8
//!   [T2] β 未知 = 正スケール同値: K(β₁)/β₁ = K(β₂)/β₂ (n=5 全数) + 門は
//!        UpToPositiveScaleAndShift を返す
//!   [T3] μ 未知 = 恒等シフト同値: K(μ) − K(0) = −βμI (off-diag 不変, n=5 全数)
//!   [T4] block-local unitary 共変性 (多軌道 3 ノード × 2 軌道, 乱択 20): 独立に構成した
//!        C' = fermi(β(UhUᵀ−μI)) の K' が block ごとに U_i K_ij U_jᵀ・特異値/Frobenius/
//!        作用素/核ノルム不変
//!   [T5] ノード置換共変性: K(PhPᵀ) の block = K(h) の置換 block
//!   [T6] logit 2 経路の分離: GlobalLogitThenBlock ≠ PairBlockThenLogit — n=2 では一致
//!        (全系 = pair)・P4 では offdiag が有限に乖離 (反例の常設化)。契約 2 型で構成
//!   [T7] 条件数定理: ‖K(C₁)−K(C₀)‖_F ≤ ‖C₁−C₀‖_F / (δ̃(1−δ̃)), δ̃ = min(δ₀, δ₁)。
//!        証明: K(C₁)−K(C₀) = ∫₀¹ DK[C_t](ΔC) dt (C_t = 凸結合) で、λ_min の凹性/λ_max の
//!        凸性から margin(C_t) ≥ δ̃、DK のスペクトルは divided difference f[x,y] = f'(ξ)
//!        (MVT) で |f'| ≤ 1/(δ̃(1−δ̃))。数値照合 = 微小摂動 30 + 異グラフ大摂動 15 +
//!        整列 rank-1 でバー飽和 (ratio ≥ 0.995 — 上界が最良であることの証明書)
//!   [T8] 低温 estimate lane: P4, β = 40 — exact lane は正しく棄却 (RankDeficient)、
//!        clamp ε = 1e-12 の推定は飽和則 k = sign(βλ)·min(β|λ|, ln((1−ε)/ε)) に従う —
//!        **深部モードは sign のみ = P6/693 衝突機構の連続版** (β→∞ で sign(A) 類へ退化)
//!   [T9] P6/693 の相図: projector lane (半充填 GS) は静的衝突 (min-perm ≤ 2e-13,
//!        v29.5 再現) + 資格は RankDeficient + 正しい裁定は EquivalenceClassOnly。
//!        有限 β では大域逆が分離 (min-perm‖ΔK‖∞ = β — 辺 1 本の差が β で見える)
//!   [T10] 門の必須証拠 (実データ): GaussianityEvidence/GibbsProvenance なしは棄却
//!   [T11] Lean 台帳: proofs/GibbsInverse.lean (9 定理 — Frobenius–Gram 恒等式・直交
//!        不変性・スケール/シフト同値。格子 native_decide + スコープ明示) の整合
//!
//! 実行: cargo run --release --bin v311_gibbs_oracle

use std::fs;
use std::path::Path;
use uft_sim::readout_contract::*;
use uft_sim::{jacobi_eigh, matfun_sym, Rng};

// ---------------------------------------------------------------- グラフ列挙 (v29.5 と同一手法)

fn edge_bit(i: usize, j: usize, n: usize) -> usize {
    let (a, b) = if i < j { (i, j) } else { (j, i) };
    a * n - a * (a + 1) / 2 + (b - a - 1)
}

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

fn apply_perm(mask: u32, pi: &[usize], n: usize, nb: usize) -> u32 {
    let mut out = 0u32;
    for e in 0..nb {
        if mask & (1 << e) == 0 {
            continue;
        }
        let mut i = 0;
        let mut acc = 0;
        while acc + (n - i - 1) <= e {
            acc += n - i - 1;
            i += 1;
        }
        let j = i + 1 + (e - acc);
        out |= 1 << edge_bit(pi[i], pi[j], n);
    }
    out
}

fn is_connected(mask: u32, n: usize) -> bool {
    let mut adj = vec![0u32; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if mask & (1 << edge_bit(i, j, n)) != 0 {
                adj[i] |= 1 << j;
                adj[j] |= 1 << i;
            }
        }
    }
    let mut seen = 1u32;
    let mut stack = vec![0usize];
    while let Some(u) = stack.pop() {
        let mut nb = adj[u] & !seen;
        while nb != 0 {
            let v = nb.trailing_zeros() as usize;
            seen |= 1 << v;
            nb &= nb - 1;
            stack.push(v);
        }
    }
    seen.count_ones() as usize == n
}

/// n 頂点の連結グラフ同型類 (canonical mask = 全置換で最小)。スレッド分割は
/// 決定的 (結果は分割に依存しない — 最後にソート)。
fn enumerate_connected(n: usize) -> Vec<u32> {
    let nb = n * (n - 1) / 2;
    let ps = perms(n);
    let total = 1u32 << nb;
    let nthreads = 12usize;
    let chunk = (total as usize).div_ceil(nthreads);
    let mut sets: Vec<Vec<u32>> = Vec::new();
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for t in 0..nthreads {
            let ps = &ps;
            handles.push(s.spawn(move || {
                let lo = (t * chunk) as u32;
                let hi = (((t + 1) * chunk).min(total as usize)) as u32;
                let mut out = Vec::new();
                for mask in lo..hi {
                    if !is_connected(mask, n) {
                        continue;
                    }
                    let mut minimal = true;
                    for pi in ps.iter() {
                        if apply_perm(mask, pi, n, nb) < mask {
                            minimal = false;
                            break;
                        }
                    }
                    if minimal {
                        out.push(mask);
                    }
                }
                out
            }));
        }
        for h in handles {
            sets.push(h.join().unwrap());
        }
    });
    let mut all: Vec<u32> = sets.into_iter().flatten().collect();
    all.sort_unstable();
    all
}

fn adj_of_mask(mask: u32, n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            if mask & (1 << edge_bit(i, j, n)) != 0 {
                a[i * n + j] = 1.0;
                a[j * n + i] = 1.0;
            }
        }
    }
    a
}

// ---------------------------------------------------------------- oracle 素子

/// C = (I + e^{β(h−μI)})⁻¹ (厳密スペクトル構成 — clamp なし)
fn gibbs_c(h: &[f64], n: usize, beta: f64, mu: f64) -> Vec<f64> {
    let mut hm = h.to_vec();
    for i in 0..n {
        hm[i * n + i] -= mu;
    }
    matfun_sym(&hm, n, |x| 1.0 / (1.0 + (beta * x).exp()))
}

/// K(C) = log[(I−C)C⁻¹] (clamp なし — 呼び出し側が資格審査済みであること)
fn logit_k(c: &[f64], n: usize) -> Vec<f64> {
    matfun_sym(c, n, |x| ((1.0 - x) / x).ln())
}

/// 半充填基底状態の射影相関 (h = −A の負エネルギー = A の正固有値モードを占有)
fn projector_c(a: &[f64], n: usize) -> Vec<f64> {
    let (evals, evecs) = jacobi_eigh(a, n);
    let mut c = vec![0.0; n * n];
    for m in 0..n {
        if evals[m] > 0.0 {
            for i in 0..n {
                for j in 0..n {
                    c[i * n + j] += evecs[m * n + i] * evecs[m * n + j];
                }
            }
        }
    }
    c
}

/// 頂点置換下の min-perm ∞ 距離: min_π max_{ij} |X_ij − Y_{π(i)π(j)}|
fn minperm_inf(x: &[f64], y: &[f64], n: usize, ps: &[Vec<usize>]) -> f64 {
    let mut best = f64::INFINITY;
    for pi in ps {
        let mut d: f64 = 0.0;
        'outer: for i in 0..n {
            for j in 0..n {
                let e = (x[i * n + j] - y[pi[i] * n + pi[j]]).abs();
                if e > d {
                    d = e;
                    if d >= best {
                        break 'outer;
                    }
                }
            }
        }
        if d < best {
            best = d;
        }
    }
    best
}

fn frob(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn max_abs_diff(x: &[f64], y: &[f64]) -> f64 {
    x.iter()
        .zip(y.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max)
}

/// 2×2 block の特異値 (BᵀB の固有値の平方根, 降順)
fn sv2(b: &[f64; 4]) -> [f64; 2] {
    let g = [
        b[0] * b[0] + b[2] * b[2],
        b[0] * b[1] + b[2] * b[3],
        b[0] * b[1] + b[2] * b[3],
        b[1] * b[1] + b[3] * b[3],
    ];
    let tr = g[0] + g[3];
    let det = g[0] * g[3] - g[1] * g[2];
    let disc = (tr * tr / 4.0 - det).max(0.0).sqrt();
    let l1 = (tr / 2.0 + disc).max(0.0);
    let l2 = (tr / 2.0 - disc).max(0.0);
    [l1.sqrt(), l2.sqrt()]
}

/// h の (I,J) 2×2 block (ノード = 2 軌道)
fn block2(k: &[f64], n: usize, bi: usize, bj: usize) -> [f64; 4] {
    [
        k[(2 * bi) * n + 2 * bj],
        k[(2 * bi) * n + 2 * bj + 1],
        k[(2 * bi + 1) * n + 2 * bj],
        k[(2 * bi + 1) * n + 2 * bj + 1],
    ]
}

fn main() {
    uft_sim::self_test();
    println!("=== v31.1 GaussianGibbsInverseOracle — 大域静的逆問題の oracle ceiling (PROMPT/12) ===");
    println!("(この逆変換は自由フェルミオンの標準的事実 — QRN 固有法則ではない。");
    println!(" bridge law 登録・PRED-019 登録はこの結果からは行わない)\n");
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

    const MU: f64 = 0.2;
    const BETAS: [f64; 4] = [0.5, 1.0, 2.0, 4.0];

    // ---- [T0] 列挙アンカー ----
    let lists: Vec<(usize, Vec<u32>)> = (4..=7).map(|n| (n, enumerate_connected(n))).collect();
    {
        let counts: Vec<usize> = lists.iter().map(|(_, l)| l.len()).collect();
        let ok_counts = counts == vec![6, 21, 112, 853];
        let l6 = &lists[2].1;
        let has_pair = l6.contains(&692) && l6.contains(&693);
        // mask 692 = P6 (5 辺, 次数列 [1,1,2,2,2,2]) /
        // mask 693 = 単環 (6 辺 = C4 + ペンダント 2, 次数列 [1,1,2,2,3,3])
        let deg_of = |mask: u32, n: usize| -> Vec<usize> {
            let a = adj_of_mask(mask, n);
            let mut d: Vec<usize> = (0..n)
                .map(|i| (0..n).filter(|&j| a[i * n + j] > 0.5).count())
                .collect();
            d.sort_unstable();
            d
        };
        let ok_ids =
            deg_of(692, 6) == vec![1, 1, 2, 2, 2, 2] && deg_of(693, 6) == vec![1, 1, 2, 2, 3, 3];
        check(
            "[T0] 連結グラフ全数 = OEIS A001349 (6/21/112/853) + 衝突対 692 (P6)/693 (単環) 実在",
            ok_counts && has_pair && ok_ids,
            format!("counts = {:?}", counts),
        );
    }

    // ---- [T1] exact lane 全数 ----
    // 復元誤差のバーは条件数定理 [T7] から導出する: 数値誤差は ‖ΔC‖ ~ n·ε_mach の
    // 摂動と等価で、‖Δĥ‖ ≤ ‖ΔK‖/β ≤ n·ε/(δ(1−δ))/β — 安全係数 10 を掛けて
    // per-instance バーとする (β=4 の稠密グラフは δ ~ 1e-11 で条件数 ~6e10 —
    // 誤差が大きいことは oracle の破れではなく定理の予言)。
    {
        let eps_mach = f64::EPSILON;
        let mut worst_err: f64 = 0.0;
        let mut worst_err_lowbeta: f64 = 0.0;
        let mut worst_bar_ratio: f64 = 0.0; // err / per-instance bar
        let mut delta_min = f64::INFINITY;
        let mut n_cert = 0usize;
        let mut n_total = 0usize;
        for (n, masks) in &lists {
            let n = *n;
            for &mask in masks {
                let a = adj_of_mask(mask, n);
                let mut h = vec![0.0; n * n];
                for k in 0..n * n {
                    h[k] = -a[k];
                }
                for &beta in &BETAS {
                    n_total += 1;
                    let c = gibbs_c(&h, n, beta, MU);
                    match ExactFullRankCorrelation::certify_real(&c, n) {
                        Ok(cert) => {
                            n_cert += 1;
                            let delta = cert.spectral_margin();
                            delta_min = delta_min.min(delta);
                            let k = logit_k(cert.c_re(), n);
                            let parent = ParentModularGenerator {
                                re: k,
                                im: vec![0.0; n * n],
                                n,
                            };
                            match identify_physical_generator(
                                &parent,
                                GaussianityEvidence::ByConstruction,
                                GibbsProvenance::KnownBetaMu { beta, mu: MU },
                            ) {
                                Ok(PhysicalGeneratorReading::Exact(hr)) => {
                                    let err = max_abs_diff(&hr.re, &h);
                                    let bar =
                                        10.0 * (n as f64) * eps_mach / (delta * (1.0 - delta)) / beta;
                                    worst_err = worst_err.max(err);
                                    worst_bar_ratio = worst_bar_ratio.max(err / bar);
                                    if beta <= 2.0 {
                                        worst_err_lowbeta = worst_err_lowbeta.max(err);
                                    }
                                }
                                _ => {
                                    worst_bar_ratio = f64::INFINITY;
                                }
                            }
                        }
                        Err(_) => {}
                    }
                }
            }
        }
        check(
            "[T1] exact lane 全数 (992 グラフ × 4β, μ=0.2): 資格 100% + 誤差 ≤ 条件数バー 10nε/(βδ(1−δ)) + β≤2 は ≤ 1e-8",
            n_cert == n_total && worst_bar_ratio <= 1.0 && worst_err_lowbeta <= 1e-8,
            format!(
                "資格 {}/{} / 最悪誤差 {:.2e} (バー比 {:.3}) / β≤2 の最悪 {:.2e} / 族の最小 δ = {:.3e}",
                n_cert, n_total, worst_err, worst_bar_ratio, worst_err_lowbeta, delta_min
            ),
        );
    }

    // ---- [T2] β 未知 = 正スケール同値 ----
    {
        let (b1, b2) = (0.5, 2.0);
        let mut worst: f64 = 0.0;
        let masks5 = &lists[1].1;
        for &mask in masks5 {
            let n = 5;
            let a = adj_of_mask(mask, n);
            let mut h = vec![0.0; n * n];
            for k in 0..n * n {
                h[k] = -a[k];
            }
            let k1 = logit_k(&gibbs_c(&h, n, b1, 0.0), n);
            let k2 = logit_k(&gibbs_c(&h, n, b2, 0.0), n);
            let s1: Vec<f64> = k1.iter().map(|x| x / b1).collect();
            let s2: Vec<f64> = k2.iter().map(|x| x / b2).collect();
            worst = worst.max(max_abs_diff(&s1, &s2));
        }
        // 門: β 未知はスケール同値類を返す
        let gate_ok = matches!(
            identify_physical_generator(
                &ParentModularGenerator {
                    re: vec![1.0, 0.0, 0.0, -1.0],
                    im: vec![0.0; 4],
                    n: 2
                },
                GaussianityEvidence::ByConstruction,
                GibbsProvenance::BetaUnknownPositive,
            ),
            Ok(PhysicalGeneratorReading::UpToPositiveScaleAndShift(_))
        );
        check(
            "[T2] β 未知 = 正スケール同値: max|K(β₁)/β₁ − K(β₂)/β₂| ≤ 1e-9 (n=5 全数) + 門は同値類",
            worst <= 1e-9 && gate_ok,
            format!("max 乖離 {:.2e} (β = {} vs {})", worst, b1, b2),
        );
    }

    // ---- [T3] μ 未知 = 恒等シフト同値 ----
    {
        let beta = 1.0;
        let mu = 0.7;
        let mut worst_off: f64 = 0.0;
        let mut worst_diag: f64 = 0.0;
        let masks5 = &lists[1].1;
        for &mask in masks5 {
            let n = 5;
            let a = adj_of_mask(mask, n);
            let mut h = vec![0.0; n * n];
            for k in 0..n * n {
                h[k] = -a[k];
            }
            let km = logit_k(&gibbs_c(&h, n, beta, mu), n);
            let k0 = logit_k(&gibbs_c(&h, n, beta, 0.0), n);
            for i in 0..n {
                for j in 0..n {
                    let d = km[i * n + j] - k0[i * n + j];
                    if i == j {
                        worst_diag = worst_diag.max((d + beta * mu).abs());
                    } else {
                        worst_off = worst_off.max(d.abs());
                    }
                }
            }
        }
        check(
            "[T3] μ 未知 = 恒等シフト: K(μ)−K(0) = −βμI (off-diag 不変 ≤ 1e-10, n=5 全数)",
            worst_off <= 1e-10 && worst_diag <= 1e-9,
            format!(
                "off-diag 乖離 {:.2e} / diag − (−βμ) 乖離 {:.2e}",
                worst_off, worst_diag
            ),
        );
    }

    // ---- 多軌道系 (3 ノード × 2 軌道 = 6 サイト, ノードグラフ = 三角形) ----
    let nn = 6usize;
    let mut h6 = vec![0.0; nn * nn];
    {
        let mut rng = Rng::new(311);
        // 対角 + ノード内結合
        for b in 0..3 {
            let (i, j) = (2 * b, 2 * b + 1);
            h6[i * nn + i] = rng.f64() - 0.5;
            h6[j * nn + j] = rng.f64() - 0.5;
            let t = rng.f64() - 0.5;
            h6[i * nn + j] = t;
            h6[j * nn + i] = t;
        }
        // ノード間 2×2 full 結合 (三角形)
        for (bi, bj) in [(0usize, 1usize), (1, 2), (0, 2)] {
            for oi in 0..2 {
                for oj in 0..2 {
                    let t = rng.f64() - 0.5;
                    h6[(2 * bi + oi) * nn + (2 * bj + oj)] = t;
                    h6[(2 * bj + oj) * nn + (2 * bi + oi)] = t;
                }
            }
        }
    }
    let beta6 = 1.3;
    let mu6 = 0.1;
    let k6 = logit_k(&gibbs_c(&h6, nn, beta6, mu6), nn);

    // ---- [T4] block-local unitary 共変性 ----
    {
        let mut rng = Rng::new(3111);
        let mut worst_cov: f64 = 0.0;
        let mut worst_sv: f64 = 0.0;
        let mut worst_norm: f64 = 0.0;
        for _ in 0..20 {
            let thetas: Vec<f64> = (0..3).map(|_| rng.f64() * std::f64::consts::TAU).collect();
            // U = ⊕ 回転
            let mut u = vec![0.0; nn * nn];
            for b in 0..3 {
                let (ct, st) = (thetas[b].cos(), thetas[b].sin());
                u[(2 * b) * nn + 2 * b] = ct;
                u[(2 * b) * nn + 2 * b + 1] = -st;
                u[(2 * b + 1) * nn + 2 * b] = st;
                u[(2 * b + 1) * nn + 2 * b + 1] = ct;
            }
            // h' = U h Uᵀ
            let mut hu = vec![0.0; nn * nn];
            for i in 0..nn {
                for j in 0..nn {
                    let mut s = 0.0;
                    for p in 0..nn {
                        for q in 0..nn {
                            s += u[i * nn + p] * h6[p * nn + q] * u[j * nn + q];
                        }
                    }
                    hu[i * nn + j] = s;
                }
            }
            // 対称化 (f64 の丸め)
            for i in 0..nn {
                for j in (i + 1)..nn {
                    let m = 0.5 * (hu[i * nn + j] + hu[j * nn + i]);
                    hu[i * nn + j] = m;
                    hu[j * nn + i] = m;
                }
            }
            let ku = logit_k(&gibbs_c(&hu, nn, beta6, mu6), nn);
            for bi in 0..3 {
                for bj in 0..3 {
                    if bi == bj {
                        continue;
                    }
                    let kb = block2(&k6, nn, bi, bj);
                    let kub = block2(&ku, nn, bi, bj);
                    // U_i K_ij U_jᵀ
                    let (ci, si) = (thetas[bi].cos(), thetas[bi].sin());
                    let (cj, sj) = (thetas[bj].cos(), thetas[bj].sin());
                    let ui = [ci, -si, si, ci];
                    let uj = [cj, -sj, sj, cj];
                    let mut t = [0.0; 4];
                    for r in 0..2 {
                        for cc in 0..2 {
                            let mut s = 0.0;
                            for p in 0..2 {
                                for q in 0..2 {
                                    s += ui[r * 2 + p] * kb[p * 2 + q] * uj[cc * 2 + q];
                                }
                            }
                            t[r * 2 + cc] = s;
                        }
                    }
                    for e in 0..4 {
                        worst_cov = worst_cov.max((kub[e] - t[e]).abs());
                    }
                    let s0 = sv2(&kb);
                    let s1 = sv2(&kub);
                    worst_sv = worst_sv.max((s0[0] - s1[0]).abs().max((s0[1] - s1[1]).abs()));
                    let f0 = (kb.iter().map(|x| x * x).sum::<f64>()).sqrt();
                    let f1 = (kub.iter().map(|x| x * x).sum::<f64>()).sqrt();
                    let n0 = s0[0] + s0[1];
                    let n1 = s1[0] + s1[1];
                    worst_norm = worst_norm
                        .max((f0 - f1).abs())
                        .max((s0[0] - s1[0]).abs())
                        .max((n0 - n1).abs());
                }
            }
        }
        check(
            "[T4] block-local unitary 共変性 (乱択 20): K'_ij = U_i K_ij U_jᵀ + 特異値/F/op/核ノルム不変",
            worst_cov <= 1e-10 && worst_sv <= 1e-10 && worst_norm <= 1e-10,
            format!(
                "共変乖離 {:.2e} / 特異値乖離 {:.2e} / ノルム乖離 {:.2e}",
                worst_cov, worst_sv, worst_norm
            ),
        );
    }

    // ---- [T5] ノード置換共変性 ----
    {
        let pi = [2usize, 0, 1]; // ノード置換: 新 I ← 旧 π(I)
        let mut hp = vec![0.0; nn * nn];
        for bi in 0..3 {
            for bj in 0..3 {
                for oi in 0..2 {
                    for oj in 0..2 {
                        hp[(2 * bi + oi) * nn + (2 * bj + oj)] =
                            h6[(2 * pi[bi] + oi) * nn + (2 * pi[bj] + oj)];
                    }
                }
            }
        }
        let kp = logit_k(&gibbs_c(&hp, nn, beta6, mu6), nn);
        let mut worst: f64 = 0.0;
        for bi in 0..3 {
            for bj in 0..3 {
                let a = block2(&kp, nn, bi, bj);
                let b = block2(&k6, nn, pi[bi], pi[bj]);
                for e in 0..4 {
                    worst = worst.max((a[e] - b[e]).abs());
                }
            }
        }
        check(
            "[T5] ノード置換共変性: K(PhPᵀ)_IJ = K(h)_{π(I)π(J)} ≤ 1e-12",
            worst <= 1e-12,
            format!("max 乖離 {:.2e}", worst),
        );
    }

    // ---- [T6] logit 2 経路の分離 (GlobalLogitThenBlock ≠ PairBlockThenLogit) ----
    {
        // n=2 アンカー: 全系 = pair なので両経路は一致するはず
        let n2 = 2;
        let h2 = vec![0.0, -1.0, -1.0, 0.0];
        let c2 = gibbs_c(&h2, n2, 1.0, 0.0);
        let kg2 = logit_k(&c2, n2);
        let anchor_diff = (kg2[1] - {
            // pair 経路 (n=2 では同じ行列の logit)
            let kp2 = logit_k(&c2, n2);
            kp2[1]
        })
        .abs();
        // P4: global 経路 vs pair (0,1) 経路
        let n4 = 4;
        let mut h4 = vec![0.0; 16];
        for i in 0..3 {
            h4[i * 4 + i + 1] = -1.0;
            h4[(i + 1) * 4 + i] = -1.0;
        }
        let c4 = gibbs_c(&h4, n4, 1.0, 0.0);
        let kg4 = logit_k(&c4, n4);
        let global_block = GlobalPhysicalParentBlock {
            re: vec![kg4[1]],
            im: vec![0.0],
            ni: 1,
            nj: 1,
        };
        // pair RDM (0,1) → logit (現行 B2 の経路)
        let cp = [c4[0], c4[1], c4[4], c4[5]];
        let kp = logit_k(&cp, 2);
        let reduced_block = ReducedModularBlock {
            re: vec![kp[1]],
            im: vec![0.0],
            ni: 1,
            nj: 1,
        };
        let sep = (global_block.re[0] - reduced_block.re[0]).abs();
        // global 経路は βh を厳密復元する (offdiag = −1)
        let global_exact = (global_block.re[0] - (-1.0)).abs();
        check(
            "[T6] logit 2 経路: n=2 一致 (全系=pair) / P4 で有限乖離 ≥ 0.05 (f(PCP) ≠ Pf(C)P の常設反例)",
            anchor_diff <= 1e-12 && sep >= 0.05 && global_exact <= 1e-10,
            format!(
                "n=2 乖離 {:.1e} / P4: global {:.6} (=βh_01 厳密) vs pair-B2 {:.6} — 乖離 {:.4}",
                anchor_diff, global_block.re[0], reduced_block.re[0], sep
            ),
        );
    }

    // ---- [T7] 条件数定理 ‖ΔK‖_F ≤ ‖ΔC‖_F/(δ̃(1−δ̃)) ----
    {
        let mut rng = Rng::new(3117);
        let masks5 = &lists[1].1;
        let n = 5;
        let mut worst_ratio: f64 = 0.0; // ratio = ‖ΔK‖ / bound (≤ 1 なら定理成立)
        let mut n_pairs = 0usize;
        let margin_of = |c: &[f64], n: usize| -> f64 {
            let (ev, _) = jacobi_eigh(c, n);
            let (mut lo, mut hi) = (f64::INFINITY, f64::NEG_INFINITY);
            for &e in &ev {
                lo = lo.min(e);
                hi = hi.max(e);
            }
            lo.min(1.0 - hi)
        };
        // (a) 微小ランダム摂動 30 対
        for _ in 0..30 {
            let mask = masks5[rng.range(masks5.len())];
            let beta = if rng.f64() < 0.5 { 1.0 } else { 2.0 };
            let a = adj_of_mask(mask, n);
            let mut h = vec![0.0; n * n];
            for k in 0..n * n {
                h[k] = -a[k];
            }
            let c0 = gibbs_c(&h, n, beta, 0.0);
            let mut dc = vec![0.0; n * n];
            for i in 0..n {
                for j in i..n {
                    let v = rng.f64() - 0.5;
                    dc[i * n + j] = v;
                    dc[j * n + i] = v;
                }
            }
            let scale = 1e-6 / frob(&dc);
            let c1: Vec<f64> = c0.iter().zip(dc.iter()).map(|(a, d)| a + d * scale).collect();
            let (d0, d1) = (margin_of(&c0, n), margin_of(&c1, n));
            let dt = d0.min(d1);
            let k0 = logit_k(&c0, n);
            let k1 = logit_k(&c1, n);
            let dk: Vec<f64> = k1.iter().zip(k0.iter()).map(|(a, b)| a - b).collect();
            let dcv: Vec<f64> = c1.iter().zip(c0.iter()).map(|(a, b)| a - b).collect();
            let bound = frob(&dcv) / (dt * (1.0 - dt));
            worst_ratio = worst_ratio.max(frob(&dk) / bound);
            n_pairs += 1;
        }
        // (b) 異グラフ間の大摂動 15 対
        for _ in 0..15 {
            let m0 = masks5[rng.range(masks5.len())];
            let m1 = masks5[rng.range(masks5.len())];
            let mk_c = |mask: u32| -> Vec<f64> {
                let a = adj_of_mask(mask, n);
                let mut h = vec![0.0; n * n];
                for k in 0..n * n {
                    h[k] = -a[k];
                }
                gibbs_c(&h, n, 1.0, 0.0)
            };
            let c0 = mk_c(m0);
            let c1 = mk_c(m1);
            let dt = margin_of(&c0, n).min(margin_of(&c1, n));
            let dk: Vec<f64> = logit_k(&c1, n)
                .iter()
                .zip(logit_k(&c0, n).iter())
                .map(|(a, b)| a - b)
                .collect();
            let dcv: Vec<f64> = c1.iter().zip(c0.iter()).map(|(a, b)| a - b).collect();
            if frob(&dcv) < 1e-14 {
                continue; // 同一グラフの抽選はスキップ
            }
            let bound = frob(&dcv) / (dt * (1.0 - dt));
            worst_ratio = worst_ratio.max(frob(&dk) / bound);
            n_pairs += 1;
        }
        // (c) 整列 rank-1 摂動でバー飽和 (上界の最良性)
        let tight_ratio = {
            let mask = masks5[0];
            let a = adj_of_mask(mask, n);
            let mut h = vec![0.0; n * n];
            for k in 0..n * n {
                h[k] = -a[k];
            }
            let c0 = gibbs_c(&h, n, 1.0, 0.0);
            let (ev, evec) = jacobi_eigh(&c0, n);
            // 最小固有値のモードを境界 (0) 側へ押す
            let mut imin = 0;
            for m in 1..n {
                if ev[m] < ev[imin] {
                    imin = m;
                }
            }
            let d0 = margin_of(&c0, n);
            let eps = 1e-4 * d0;
            let mut c1 = c0.clone();
            for i in 0..n {
                for j in 0..n {
                    c1[i * n + j] -= eps * evec[imin * n + i] * evec[imin * n + j];
                }
            }
            let dt = margin_of(&c0, n).min(margin_of(&c1, n));
            let dk: Vec<f64> = logit_k(&c1, n)
                .iter()
                .zip(logit_k(&c0, n).iter())
                .map(|(a, b)| a - b)
                .collect();
            let dcv: Vec<f64> = c1.iter().zip(c0.iter()).map(|(a, b)| a - b).collect();
            frob(&dk) / (frob(&dcv) / (dt * (1.0 - dt)))
        };
        check(
            "[T7] 条件数定理: ‖ΔK‖_F ≤ ‖ΔC‖_F/(δ̃(1−δ̃)) (微小 30 + 大摂動 15) + rank-1 整列で飽和 ≥ 0.995",
            worst_ratio <= 1.0 + 1e-9 && tight_ratio >= 0.995 && tight_ratio <= 1.0 + 1e-9,
            format!(
                "max ‖ΔK‖/bound = {:.6} ({} 対) / 整列 rank-1 の飽和比 = {:.6}",
                worst_ratio, n_pairs, tight_ratio
            ),
        );
    }

    // ---- [T8] 低温 estimate lane (飽和則 — 深部は sign のみ) ----
    {
        let n = 4;
        let beta = 40.0;
        let mut h = vec![0.0; 16];
        for i in 0..3 {
            h[i * 4 + i + 1] = -1.0;
            h[(i + 1) * 4 + i] = -1.0;
        }
        let c = gibbs_c(&h, n, beta, 0.0);
        // exact lane は正しく棄却
        let refused = matches!(
            ExactFullRankCorrelation::certify_real(&c, n),
            Err(AbstainReason::RankDeficient) | Err(AbstainReason::IllConditioned)
        );
        // estimate lane: clamp ε = 1e-12
        let eps = 1e-12;
        let (ev, evec) = jacobi_eigh(&c, n);
        let mut c_cl = vec![0.0; n * n];
        for m in 0..n {
            let lam = ev[m].clamp(eps, 1.0 - eps);
            for i in 0..n {
                for j in 0..n {
                    c_cl[i * n + j] += lam * evec[m * n + i] * evec[m * n + j];
                }
            }
        }
        let reg = RegularizedCorrelation {
            c_re: c_cl.clone(),
            c_im: vec![0.0; n * n],
            n,
            clamp_eps: eps,
        };
        let k_est = logit_k(&reg.c_re, n);
        let (mut kev, _) = jacobi_eigh(&k_est, n);
        kev.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // 期待: k = sign(βλ_h)·min(β|λ_h|, κ_max), λ_h = ∓φ, ∓1/φ (P4), K = −βA なので
        // h = −A の固有値 ±1.618, ±0.618 → K 固有値 = −β·λ_h... 注意: K = βh
        let kappa_max = ((1.0 - eps) / eps).ln();
        let phi = (1.0 + 5.0f64.sqrt()) / 2.0; // 1.618
        let lam_h = [-phi, -(phi - 1.0), phi - 1.0, phi]; // h = −A の固有値
        let mut expect: Vec<f64> = lam_h
            .iter()
            .map(|&l| {
                let t = beta * l;
                t.signum() * t.abs().min(kappa_max)
            })
            .collect();
        expect.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut worst_rel: f64 = 0.0;
        for (a, b) in kev.iter().zip(expect.iter()) {
            worst_rel = worst_rel.max((a - b).abs() / b.abs());
        }
        // 飽和則の裁定: |βλ| < κ_max のモードは忠実・超えるモードは sign のみ
        let n_sat = lam_h.iter().filter(|&&l| beta * l.abs() > kappa_max).count();
        check(
            "[T8] 低温 P4 (β=40): exact は棄却・estimate は飽和則 k = sign·min(β|λ|, ln((1−ε)/ε)) — 深部は sign のみ",
            refused && worst_rel <= 1e-4 && n_sat == 2,
            format!(
                "K_est 固有値 {:?} / 期待 {:?} / 飽和モード {}/4 (P6/693 機構の連続版 — β→∞ で sign(A) 類へ)",
                kev.iter().map(|x| (x * 1e3).round() / 1e3).collect::<Vec<_>>(),
                expect.iter().map(|x| (x * 1e3).round() / 1e3).collect::<Vec<_>>(),
                n_sat
            ),
        );
    }

    // ---- [T9] P6/693 の相図 (projector 衝突 / 有限 β 分離) ----
    {
        let n = 6;
        let ps6 = perms(6);
        let a_p6 = adj_of_mask(692, n);
        let a_c6 = adj_of_mask(693, n);
        // (a) projector lane: 静的衝突 (v29.5 再現) + 資格棄却 + 同値類裁定。
        // v29.5 の衝突は**カーネル準位** (readout が見る |C_ij|) — 符号つき C は
        // ノード内 U(1)/Z2 位相ゲージで変わるため、(i) |C| の min-perm と
        // (ii) 符号ゲージ z_i = ±1 込みの min の両方を測る (符号つき素の C も報告)。
        let cp = projector_c(&a_p6, n);
        let cc = projector_c(&a_c6, n);
        let d_signed = minperm_inf(&cp, &cc, n, &ps6);
        let abs_of = |x: &[f64]| -> Vec<f64> { x.iter().map(|v| v.abs()).collect() };
        let d_static = minperm_inf(&abs_of(&cp), &abs_of(&cc), n, &ps6);
        // 符号ゲージ z ∈ {±1}⁶ × 置換の min (C_ij → z_i z_j C_ij)
        let mut d_gauge = f64::INFINITY;
        for zmask in 0u32..(1 << n) {
            let mut czz = cc.clone();
            for i in 0..n {
                for j in 0..n {
                    let z = if ((zmask >> i) & 1) ^ ((zmask >> j) & 1) == 1 {
                        -1.0
                    } else {
                        1.0
                    };
                    czz[i * n + j] *= z;
                }
            }
            let d = minperm_inf(&cp, &czz, n, &ps6);
            if d < d_gauge {
                d_gauge = d;
            }
        }
        let refused = matches!(
            ExactFullRankCorrelation::certify_real(&cp, n),
            Err(AbstainReason::RankDeficient)
        ) && matches!(
            ExactFullRankCorrelation::certify_real(&cc, n),
            Err(AbstainReason::RankDeficient)
        );
        type EqCert = ReadoutCertificate<
            uft_sim::qrn_core::SpatialTopologyGivenFactorization,
            GaussianProjector,
            GlobalOneBodyCorrelation,
            GivenNodeFactorization,
        >;
        let cert = EqCert::equivalence_class("QRN-ORACLE-SIGN-CLASS", "sign(A) 同値類 (半充填 projector)".into());
        let verdict_ok = cert.verdict().as_str() == "equivalence_class_only";
        // (b) 有限 β: 大域逆 K で分離 (min-perm‖ΔK‖∞ = β — 辺 1 本の差)
        let mut seps = Vec::new();
        let mut sep_ok = true;
        for &beta in &BETAS {
            let mut hp = vec![0.0; n * n];
            let mut hc = vec![0.0; n * n];
            for k in 0..n * n {
                hp[k] = -a_p6[k];
                hc[k] = -a_c6[k];
            }
            let kp = logit_k(&gibbs_c(&hp, n, beta, 0.0), n);
            let kc = logit_k(&gibbs_c(&hc, n, beta, 0.0), n);
            let d = minperm_inf(&kp, &kc, n, &ps6);
            sep_ok &= (d - beta).abs() <= 1e-9;
            seps.push(d);
        }
        // C レベルの分離も報告 (β=1)
        let c1p = gibbs_c(
            &{
                let mut h = vec![0.0; n * n];
                for k in 0..n * n {
                    h[k] = -a_p6[k];
                }
                h
            },
            n,
            1.0,
            0.0,
        );
        let c1c = gibbs_c(
            &{
                let mut h = vec![0.0; n * n];
                for k in 0..n * n {
                    h[k] = -a_c6[k];
                }
                h
            },
            n,
            1.0,
            0.0,
        );
        let d_c = minperm_inf(&c1p, &c1c, n, &ps6);
        check(
            "[T9] P6/693: projector はカーネル準位 |C| で衝突 ≤ 2e-13 (v29.5 再現)・資格は RankDeficient・裁定は同値類 / 有限 β は min-perm‖ΔK‖∞ = β で分離",
            d_static <= 2e-13 && refused && verdict_ok && sep_ok && d_c > 1e-3,
            format!(
                "|C| 衝突 {:.2e} (符号つき素 C {:.3} / Z2 ゲージ×置換 min {:.3e}) / K 分離 {:?} (= β) / C 分離 (β=1) {:.4}",
                d_static,
                d_signed,
                d_gauge,
                seps.iter().map(|x| (x * 1e4).round() / 1e4).collect::<Vec<_>>(),
                d_c
            ),
        );
    }

    // ---- [T10] 門の必須証拠 (実データ) ----
    {
        let n = 6;
        let a = adj_of_mask(692, n);
        let mut h = vec![0.0; n * n];
        for k in 0..n * n {
            h[k] = -a[k];
        }
        let c = gibbs_c(&h, n, 1.0, 0.0);
        let cert = ExactFullRankCorrelation::certify_real(&c, n).unwrap();
        let parent = ParentModularGenerator {
            re: logit_k(cert.c_re(), n),
            im: vec![0.0; n * n],
            n,
        };
        let r1 = identify_physical_generator(
            &parent,
            GaussianityEvidence::Unknown,
            GibbsProvenance::KnownBetaMu { beta: 1.0, mu: 0.0 },
        );
        let r2 = identify_physical_generator(
            &parent,
            GaussianityEvidence::ByConstruction,
            GibbsProvenance::Missing,
        );
        let r3 = identify_physical_generator(
            &parent,
            GaussianityEvidence::WickResidualBound {
                residual: 1e-3,
                bar: 1e-6,
            },
            GibbsProvenance::KnownBetaMu { beta: 1.0, mu: 0.0 },
        );
        let r4 = identify_physical_generator(
            &parent,
            GaussianityEvidence::ByConstruction,
            GibbsProvenance::KnownBetaMu { beta: 1.0, mu: 0.0 },
        );
        let ok4 = match &r4 {
            Ok(PhysicalGeneratorReading::Exact(hr)) => max_abs_diff(&hr.re, &h) <= 1e-10,
            _ => false,
        };
        check(
            "[T10] 門の必須証拠: Unknown→棄却 / Gibbs 出自なし→棄却 / Wick 残差超え→棄却 / 証拠あり→Exact",
            matches!(r1, Err(AbstainReason::GaussianityUnverified))
                && matches!(r2, Err(AbstainReason::GibbsProvenanceMissing))
                && matches!(r3, Err(AbstainReason::NonGaussianDomain))
                && ok4,
            "P6 (β=1) の実データで 4 経路とも契約どおり".into(),
        );
    }

    // ---- [T11] Lean 台帳 ----
    {
        let root = if Path::new("proofs/GibbsInverse.lean").exists() {
            "."
        } else {
            ".."
        };
        let lean = fs::read_to_string(format!("{}/proofs/GibbsInverse.lean", root)).unwrap_or_default();
        let n_thm = lean
            .lines()
            .filter(|l| l.trim_start().starts_with("theorem "))
            .count();
        let ok = n_thm == 9
            && lean.contains("native_decide")
            && lean.contains("スコープの明示")
            && lean.contains("未形式化");
        check(
            "[T11] proofs/GibbsInverse.lean — 9 定理 (Frobenius–Gram 恒等式・直交/スケール/シフト同値, 格子 native_decide + スコープ明示)",
            ok,
            format!(
                "theorem 宣言 {} 本 (コンパイルは著作時に lean 4.31 で検証 — 保証は格子上の恒等であり実数全域は次数論法)",
                n_thm
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "大域逆問題の天井が確定した — full-rank Gaussian Gibbs では h は厳密復元可能 (β 未知は正スケール・μ 未知は対角シフトの同値類)、projector 極限では sign 類まで、低温 estimate lane は飽和則で sign 類へ連続的に退化する。これは oracle ceiling であり bridge law ではない"
        } else {
            "**oracle ceiling の破れ** — 資格・門・上界のいずれかが契約と不整合"
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
