//! v34.3 有限データ昇格不能定理 (第四の no-go) と Robust Promotion (PROMPT/15 §4)
//!
//! 背景: 第三十三期の識別境界は「正しい入力 (exact 証明書)」の下で閉じた。しかし
//! addressability・cross-talk・glue overlap 等は実験では有限ショットからの推定で
//! あり、証明書は申告値だった。本版は「有限データから何が言えて何が言えないか」を
//! 定理と型で閉じる:
//!
//!   第四の no-go (Le Cam 二点下限の器械化 — Lean: proofs/FiniteDataNoGo.lean):
//!     異なる裁定を要する θ0, θ1 の N ショット分布の全変動距離を TV とすると、
//!     どんな強制回答器も平均誤り ≥ (1 − TV)/2。P0 = P1 なら ≥ 1/2。
//!     境界近傍では abstention (Straddled) か追加観測だけが正答である。
//!   Robust Promotion Theorem (正側):
//!     同時信頼集合 C_α(D) (Clopper–Pearson — 厳密 tail 反転) を exact reader に
//!     通し、像が単一クラスのときだけ昇格 → P(wrong promotion) ≤ α (無条件)。
//!     selective risk は ≤ α/P(answer) — 被覆単独からは出ない (反例を厳密計算)。
//!   禁止変換 22–29 (sim/src/finite_data.rs の型 — 各反例を厳密計算):
//!     22 PointEstimate ↛ Certificate / 23 Marginal ↛ Joint / 24 GoF ↛ Model /
//!     25 t0 ↛ t1 / 26 Mean ↛ Uniform / 27 Chart ↛ Glue / 28 ZeroHoldout ↛
//!     ZeroRisk / 29 ModelConditional ↛ ModelFree
//!
//! 方法: 全て解析的な厳密和 (二項・beta-二項・Markov 転送行列) と決定規則の全数
//! 列挙 — Monte Carlo は使わない (乱数なし・決定的)。
//!
//! 検証: [F1] 下限 = 全数列挙の最小 (規則 2^{N+1} 本) = Σmin/2 の三重一致 +
//!       Lean 実例 (N=4, 整数重み 147/512) との整数一致 / [F2] P0=P1 → 1/2 /
//!       [F3] 誤昇格 ≤ α (θ 走査の厳密和) + selective risk 反例 / [F4] 周辺直積の
//!       被覆破れ vs Bonferroni / [F5] 裁定 5 値の完備 (OutOfDomain・Insufficient・
//!       境界での Straddled ≥ 1−α・符号 orbit の EquivalenceClassOnly) /
//!       [F6] 禁止変換 22–29 の反例。

use uft_sim::finite_data::{
    binom_pmf, cp_interval, finite_data_self_test, zero_error_upper_bound, RobustEdgeReader,
    RobustVerdict, SignOrbitReader,
};
use uft_sim::{ln_gamma, self_test};

/// ln B(a,b) = lnΓa + lnΓb − lnΓ(a+b)
fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}

/// beta-二項 pmf (平均 a/(a+b) — 過分散: 同一平均の別ノイズモデル)
fn beta_binom_pmf(n: usize, k: usize, a: f64, b: f64) -> f64 {
    let ln_c = ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0);
    (ln_c + ln_beta(k as f64 + a, (n - k) as f64 + b) - ln_beta(a, b)).exp()
}

/// 2 状態 Markov 鎖 (定常周辺 θ・持続 ρ) の N ショット成功数分布 — 厳密 DP
fn markov_count_dist(n: usize, theta: f64, rho: f64) -> Vec<f64> {
    // 遷移: P(1→1) = θ + ρ(1−θ), P(0→1) = θ(1−ρ) — 定常周辺 = θ (ρ = 0 で iid)
    let p11 = theta + rho * (1.0 - theta);
    let p01 = theta * (1.0 - rho);
    // dp[k][s] = k 成功で現状態 s の確率
    let mut dp = vec![[0.0f64; 2]; n + 1];
    dp[1][1] = theta; // 初期ショットは定常分布から
    dp[0][0] = 1.0 - theta;
    for _ in 1..n {
        let mut nx = vec![[0.0f64; 2]; n + 1];
        for k in 0..=n {
            for s in 0..2 {
                let pr = dp[k][s];
                if pr == 0.0 {
                    continue;
                }
                let p_to1 = if s == 1 { p11 } else { p01 };
                if k + 1 <= n {
                    nx[k + 1][1] += pr * p_to1;
                }
                nx[k][0] += pr * (1.0 - p_to1);
            }
        }
        dp = nx;
    }
    dp.iter().map(|x| x[0] + x[1]).collect()
}

fn main() {
    self_test();
    finite_data_self_test().expect("finite_data self test");
    println!("=== v34.3 有限データ昇格不能定理と Robust Promotion (PROMPT/15 §4) ===\n");
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

    // ---------------- [F1] 第四の no-go: 下限の三重一致 ----------------
    println!("[F1] Le Cam 二点下限 — 全数列挙 = Σmin/2 = (1−TV)/2 の三重一致");
    {
        let (n, th0, th1) = (16usize, 0.25f64, 0.35f64);
        let w0: Vec<f64> = (0..=n).map(|k| binom_pmf(n, k, th0)).collect();
        let w1: Vec<f64> = (0..=n).map(|k| binom_pmf(n, k, th1)).collect();
        // 全 2^{n+1} 決定規則の平均誤りの最小
        let mut min_err = f64::INFINITY;
        let nmask = 1u64 << (n + 1);
        for mask in 0..nmask {
            let mut e = 0.0;
            for k in 0..=n {
                e += if (mask >> k) & 1 == 1 { w0[k] } else { w1[k] };
            }
            let e = e / 2.0;
            if e < min_err {
                min_err = e;
            }
        }
        let sum_min: f64 = (0..=n).map(|k| w0[k].min(w1[k])).sum::<f64>() / 2.0;
        let tv: f64 = (0..=n).map(|k| (w0[k] - w1[k]).abs()).sum::<f64>() / 2.0;
        let bound = (1.0 - tv) / 2.0;
        check(
            "[F1a] N=16: min_{全 131072 規則} 平均誤り = Σmin/2 = (1−TV)/2",
            (min_err - sum_min).abs() < 1e-14 && (sum_min - bound).abs() < 1e-14,
            format!("min = {:.12}, TV = {:.6}", min_err, tv),
        );
        // Lean 実例との整数一致 (N = 4, θ = 1/4 vs 1/2, スケール 256)
        let w0i: [i64; 5] = [81, 108, 54, 12, 1];
        let w1i: [i64; 5] = [16, 64, 96, 64, 16];
        let sum_min_i: i64 = (0..5).map(|k| w0i[k].min(w1i[k])).sum();
        let abs_diff_i: i64 = (0..5).map(|k| (w0i[k] - w1i[k]).abs()).sum();
        let mut min_err_i = i64::MAX;
        for mask in 0..32u32 {
            let e: i64 = (0..5)
                .map(|k| {
                    if (mask >> k) & 1 == 1 {
                        w0i[k]
                    } else {
                        w1i[k]
                    }
                })
                .sum();
            min_err_i = min_err_i.min(e);
        }
        check(
            "[F1b] Lean 実例 (proofs/FiniteDataNoGo.lean instance_n4): 整数一致",
            sum_min_i == 147 && abs_diff_i == 218 && min_err_i == 147 && 2 * 147 == 512 - 218,
            format!(
                "Σmin = {}, Σ|diff| = {}, 全 32 規則の最小 = {} (2·147 = 512 − 218)",
                sum_min_i, abs_diff_i, min_err_i
            ),
        );
        // 二項の f64 でも同じ実例が一致 (256 倍)
        let w0f: Vec<f64> = (0..=4).map(|k| binom_pmf(4, k, 0.25) * 256.0).collect();
        let dev: f64 = (0..5)
            .map(|k| (w0f[k] - w0i[k] as f64).abs())
            .fold(0.0, f64::max);
        check(
            "[F1c] f64 二項と整数重みの一致 (Lean ↔ Rust の橋)",
            dev < 1e-9,
            format!("max|Δ| = {:.2e}", dev),
        );
    }

    // ---------------- [F2] 識別不能の極限 ----------------
    println!("\n[F2] 観測契約が区別しなければ当てられない (TV = 0 → 誤り ≥ 1/2)");
    {
        let n = 12usize;
        let w: Vec<f64> = (0..=n).map(|k| binom_pmf(n, k, 0.4)).collect();
        let mut min_err = f64::INFINITY;
        for mask in 0..(1u64 << (n + 1)) {
            let mut e = 0.0;
            for k in 0..=n {
                e += w[k]; // どちらの仮説でも同じ重み — 規則に依らず総和は一定
                let _ = mask;
            }
            min_err = min_err.min(e / 2.0);
        }
        check(
            "[F2a] P0 = P1: 全規則の平均誤り = 1/2 (規則に依らない)",
            (min_err - 0.5).abs() < 1e-14,
            format!("min = {:.12}", min_err),
        );
        // 分解能はショット数の関数: 固定 (θ0, θ1) で TV は N に単調増加
        let (th0, th1) = (0.28f64, 0.32f64);
        let tv_of = |n: usize| -> f64 {
            (0..=n)
                .map(|k| (binom_pmf(n, k, th0) - binom_pmf(n, k, th1)).abs())
                .sum::<f64>()
                / 2.0
        };
        let (t4, t16, t64, t256) = (tv_of(4), tv_of(16), tv_of(64), tv_of(256));
        check(
            "[F2b] TV(N) は単調増加 — 境界の分解能は観測量の関数 (追加実験が下限を動かす)",
            t4 < t16 && t16 < t64 && t64 < t256 && t256 < 1.0,
            format!(
                "TV = {:.4} (N=4) → {:.4} (16) → {:.4} (64) → {:.4} (256); 強制誤り下限 {:.4} → {:.4}",
                t4,
                t16,
                t64,
                t256,
                (1.0 - t4) / 2.0,
                (1.0 - t256) / 2.0
            ),
        );
    }

    // ---------------- [F3] Robust Promotion: 誤昇格 ≤ α ----------------
    println!("\n[F3] Robust Promotion — 誤昇格 ≤ α (厳密和) と selective risk の区別");
    {
        let reader = RobustEdgeReader {
            tau: 0.3,
            alpha: 0.05,
            n_min: 10,
        };
        let n = 60usize;
        let mut max_wrong = 0.0f64;
        let mut worst_theta = 0.0;
        let grid: Vec<f64> = (1..=19)
            .map(|i| i as f64 * 0.05)
            .chain([0.29, 0.30, 0.31])
            .collect();
        for &theta in &grid {
            let mut p_wrong = 0.0;
            for k in 0..=n {
                let v = reader.read_counts(k, n);
                let wrong = match &v {
                    RobustVerdict::RobustExact { reading } => {
                        (reading == "edge" && theta <= reader.tau)
                            || (reading == "no_edge" && theta > reader.tau)
                    }
                    _ => false,
                };
                if wrong {
                    p_wrong += binom_pmf(n, k, theta);
                }
            }
            if p_wrong > max_wrong {
                max_wrong = p_wrong;
                worst_theta = theta;
            }
        }
        check(
            "[F3a] 全 θ 走査で P(wrong promotion) ≤ α (実測は片側 α/2 以下)",
            max_wrong <= 0.05 + 1e-12,
            format!("max = {:.5} at θ = {} (≤ α/2 = 0.025)", max_wrong, worst_theta),
        );
        // selective risk の反例: 被覆から出るのは無条件確率だけ
        let theta = 0.305;
        let (mut p_wrong, mut p_ans) = (0.0, 0.0);
        for k in 0..=n {
            let v = reader.read_counts(k, n);
            if let RobustVerdict::RobustExact { reading } = &v {
                let p = binom_pmf(n, k, theta);
                p_ans += p;
                if reading == "no_edge" {
                    p_wrong += p;
                }
            }
        }
        let cond = p_wrong / p_ans;
        check(
            "[F3b] selective risk 反例: P(wrong|answer) > α なのに P(wrong) ≤ α — 被覆単独から selective risk は出ない",
            p_wrong <= 0.05 && cond > 0.05,
            format!(
                "θ = 0.305: P(wrong) = {:.4} ≤ α, P(answer) = {:.4}, P(wrong|answer) = {:.3} > α (上限は α/P(answer) = {:.3})",
                p_wrong,
                p_ans,
                cond,
                0.05 / p_ans
            ),
        );
    }

    // ---------------- [F4] 禁止変換 23/27: 周辺直積 ≠ joint ----------------
    println!("\n[F4] MarginalIntervals ↛ JointConfidenceRegion / LocalChartCoverage ↛ GlobalGlueCoverage");
    {
        let (n, theta, m) = (40usize, 0.25f64, 6u32);
        let cov_at = |alpha: f64| -> f64 {
            (0..=n)
                .filter(|&k| {
                    let (lo, hi) = cp_interval(k, n, alpha);
                    lo <= theta && theta <= hi
                })
                .map(|k| binom_pmf(n, k, theta))
                .sum()
        };
        let c_marginal = cov_at(0.05);
        let joint_naive = c_marginal.powi(m as i32);
        let c_bonf = cov_at(0.05 / m as f64);
        let joint_bonf = c_bonf.powi(m as i32);
        check(
            "[F4a] 周辺 95% 区間の直積 (m = 6 独立パラメータ): joint 被覆 < 95% (破れ)",
            joint_naive < 0.95 && c_marginal >= 0.95,
            format!(
                "周辺被覆 {:.4} → 直積 {:.4} (miscoverage {:.3} > α — chart 局所被覆も同型で glue に昇格しない)",
                c_marginal,
                joint_naive,
                1.0 - joint_naive
            ),
        );
        check(
            "[F4b] Bonferroni (α/m): joint 被覆 ≥ 95% (正しい同時領域の最小構成)",
            joint_bonf >= 0.95,
            format!("α/6 周辺 {:.5} → 直積 {:.4}", c_bonf, joint_bonf),
        );
    }

    // ---------------- [F5] 裁定 5 値の完備 ----------------
    println!("\n[F5] 裁定 5 値 — 各値が到達可能で、境界では Straddled が保証つきの正答");
    {
        let reader = RobustEdgeReader {
            tau: 0.3,
            alpha: 0.05,
            n_min: 10,
        };
        let ood = reader.read(&[0u8, 1, 2, 0]);
        let insuf = reader.read(&[1u8, 0, 1]);
        check(
            "[F5a] OutOfDomain (支持 {0,1} 違反) と InsufficientObservation (n < n_min)",
            ood == RobustVerdict::OutOfDomain && insuf == RobustVerdict::InsufficientObservation,
            format!("{} / {}", ood.as_str(), insuf.as_str()),
        );
        // 境界 θ = τ: 被覆 ≥ 1−α ⇒ P(Straddled) ≥ 1−α (τ ∈ CI なら回答しない)
        let n = 60usize;
        let p_straddle: f64 = (0..=n)
            .filter(|&k| reader.read_counts(k, n) == RobustVerdict::Straddled)
            .map(|k| binom_pmf(n, k, reader.tau))
            .sum();
        check(
            "[F5b] θ = τ (裁定境界上): P(Straddled) ≥ 1 − α — 棄権が保証つきの正答",
            p_straddle >= 0.95,
            format!("P(Straddled) = {:.4}", p_straddle),
        );
        // 符号 orbit: TV(P₊, P₋) = 0 (契約は |θ| しか見ない) → 強制符号回答の誤り ≥ 1/2
        // EquivalenceClassOnly が唯一の回答経路 (クラス誤り ≤ α)
        let sr = SignOrbitReader {
            alpha: 0.05,
            n_min: 10,
        };
        let v = sr.read_counts(21, 60);
        let is_class = matches!(v, RobustVerdict::EquivalenceClassOnly { .. });
        // クラス誤り = |θ| が CI の外 ≤ α (被覆) — θ = ±0.35 の両方で厳密和
        let abs_theta = 0.35f64;
        let p_class_wrong: f64 = (0..=60)
            .filter(|&k| {
                let (lo, hi) = cp_interval(k, 60, 0.05);
                !(lo <= abs_theta && abs_theta <= hi)
            })
            .map(|k| binom_pmf(60, k, abs_theta))
            .sum();
        check(
            "[F5c] 符号 orbit (TV = 0): 回答は EquivalenceClassOnly のみ・クラス誤り ≤ α",
            is_class && p_class_wrong <= 0.05,
            format!(
                "verdict = {}, P(class wrong) = {:.4} (符号の強制回答は [F2a] より誤り ≥ 1/2)",
                v.as_str(),
                p_class_wrong
            ),
        );
    }

    // ---------------- [F6] 禁止変換 22/24/25/26/28/29 の反例 ----------------
    println!("\n[F6] 禁止変換の反例 (全て厳密和)");
    {
        // 22: PointEstimate ↛ Certificate
        let (n, tau, theta) = (50usize, 0.3f64, 0.28f64);
        let kmin = (tau * n as f64).floor() as usize + 1; // k/N > τ ⟺ k ≥ 16
        let p_point_wrong: f64 = (kmin..=n).map(|k| binom_pmf(n, k, theta)).sum();
        let reader = RobustEdgeReader {
            tau,
            alpha: 0.05,
            n_min: 10,
        };
        let p_robust_wrong: f64 = (0..=n)
            .filter(|&k| {
                matches!(reader.read_counts(k, n),
                    RobustVerdict::RobustExact { ref reading } if reading == "edge")
            })
            .map(|k| binom_pmf(n, k, theta))
            .sum();
        check(
            "[F6-22] PointEstimate ↛ Certificate: 点推定昇格の誤り ≫ α・robust は ≤ α",
            p_point_wrong > 0.2 && p_robust_wrong <= 0.05,
            format!(
                "θ = 0.28 (真は no_edge): 点推定 (k/N > τ で edge) の誤り = {:.3}, robust = {:.4}",
                p_point_wrong, p_robust_wrong
            ),
        );

        // 24: GoodnessOfFit ↛ NoiseModelValidity (同一平均の過分散 beta-二項)
        let (n, a, b) = (40usize, 0.5f64, 1.5f64); // 平均 a/(a+b) = 0.25
        let theta0 = 0.25f64;
        let bb_mean: f64 = (0..=n)
            .map(|k| k as f64 * beta_binom_pmf(n, k, a, b))
            .sum::<f64>()
            / n as f64;
        let cov_bb: f64 = (0..=n)
            .filter(|&k| {
                let (lo, hi) = cp_interval(k, n, 0.05);
                lo <= theta0 && theta0 <= hi
            })
            .map(|k| beta_binom_pmf(n, k, a, b))
            .sum();
        check(
            "[F6-24] GoodnessOfFit ↛ NoiseModelValidity: 一次積率は一致・被覆は崩壊",
            (bb_mean - theta0).abs() < 1e-12 && cov_bb < 0.60,
            format!(
                "beta-二項 (平均 {:.4} = θ0 厳密) の下で iid 95% 区間の実被覆 = {:.3}",
                bb_mean, cov_bb
            ),
        );

        // 25: CalibrationAt(t0) ↛ ValidAt(t1) (drift)
        let (n, th_t0, th_t1) = (100usize, 0.25f64, 0.40f64);
        let cov_drift: f64 = (0..=n)
            .filter(|&k| {
                let (lo, hi) = cp_interval(k, n, 0.05);
                lo <= th_t1 && th_t1 <= hi
            })
            .map(|k| binom_pmf(n, k, th_t0))
            .sum();
        check(
            "[F6-25] CalibrationAt(t0) ↛ ValidAt(t1): drift 下で t0 較正区間の t1 被覆が崩壊",
            cov_drift < 0.15,
            format!("θ(t0) = 0.25 → θ(t1) = 0.40: 被覆 = {:.4} (≥ 0.95 のはずが)", cov_drift),
        );

        // 26: MeanCrosstalk ↛ UniformCrosstalkBound
        let offdiag = [0.30f64, 0.01, 0.01, 0.01, 0.01, 0.01];
        let mean = offdiag.iter().sum::<f64>() / offdiag.len() as f64;
        let maxod = offdiag.iter().cloned().fold(0.0, f64::max);
        let bar = 0.1;
        check(
            "[F6-26] MeanCrosstalk ↛ UniformBound: 平均はバー内・単一対の漏れ 0.30 が独立性を破る",
            mean <= bar && maxod > bar,
            format!(
                "mean = {:.3} ≤ {} (平均証明書は通過) / max = {:.2} > {} (C2 の一様バーは拒否 — 漏れは操作のたびに 30% 乗る)",
                mean, bar, maxod, bar
            ),
        );

        // 28: ZeroHoldoutErrors ↛ ZeroPopulationRisk (PROMPT/15 §6 の数の機械化)
        let u9 = zero_error_upper_bound(9, 0.95);
        let u77 = zero_error_upper_bound(77, 0.95);
        let mut n_needed = 0usize;
        for n in 1..1000 {
            if zero_error_upper_bound(n, 0.95) <= 0.01 {
                n_needed = n;
                break;
            }
        }
        check(
            "[F6-28] ZeroHoldoutErrors ↛ ZeroPopulationRisk: 0/9 → 28.3%・0/77 → 3.8%・1% には 299 件",
            (u9 - 0.2831).abs() < 5e-4 && (u77 - 0.0382).abs() < 5e-4 && n_needed == 299,
            format!(
                "片側 95% 上限: 0/9 = {:.4}, 0/77 = {:.4}; ≤ 1% に必要な誤り 0 回答セル = {}",
                u9, u77, n_needed
            ),
        );

        // 29: ModelConditionalCertificate ↛ ModelFreeCertificate (相関鎖)
        let (n, theta0, rho) = (40usize, 0.25f64, 0.8f64);
        let dist = markov_count_dist(n, theta0, rho);
        let norm: f64 = dist.iter().sum();
        let mean_mk: f64 = dist.iter().enumerate().map(|(k, p)| k as f64 * p).sum::<f64>() / n as f64;
        let cov_mk: f64 = (0..=n)
            .filter(|&k| {
                let (lo, hi) = cp_interval(k, n, 0.05);
                lo <= theta0 && theta0 <= hi
            })
            .map(|k| dist[k])
            .sum();
        check(
            "[F6-29] ModelConditional ↛ ModelFree: 同一周辺 (定常 θ0) の相関鎖で iid 区間の被覆が崩壊",
            (norm - 1.0).abs() < 1e-12 && (mean_mk - theta0).abs() < 1e-12 && cov_mk < 0.75,
            format!(
                "Markov (ρ = 0.8, 周辺平均 {:.4} = θ0 厳密) の下で被覆 = {:.3} — 契約は iid を『登録』し、破れは OutOfDomain で拒否するのが正答",
                mean_mk, cov_mk
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "有限データの境界が閉じた — 境界近傍の強制回答は原理的に不可能 (第四の no-go)、\n       同時信頼集合の昇格は誤り ≤ α (Robust Promotion)、その保証は無条件確率であって\n       selective risk ではない。点推定・周辺区間・GoF・単一時点較正・平均バー・\n       holdout 0 件・モデル条件付き証明書は、いずれも証明書に昇格しない (22–29)。"
        } else {
            "**有限データ意味論の破れ** — no-go・被覆・反例計算を修復せよ"
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
