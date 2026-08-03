// finite_data — 有限データ証明書の型契約と同時信頼集合 (v34.3, PROMPT/15 §4)
//
// 第三十四期テーゼ「可アクセス性証明書は実験者の宣言ではなく、登録済み観測契約の
// 下で有限データから得られる同時信頼集合である」の型実装 (第一層):
//
//   1. **第四の no-go (有限データ昇格不能定理)**: 異なる裁定を要する二つの
//      interface θ0, θ1 の有限データ分布が近い (全変動距離 TV 小) とき、どんな
//      強制回答器の平均誤りも (1 − TV)/2 を下回れない
//      (Lean: proofs/FiniteDataNoGo.lean — le_cam_two_point・bayes_achieves で
//      下限は最良)。観測契約が二つを区別しない (P0 = P1) なら平均誤り ≥ 1/2。
//      境界近傍で常に正しい強制回答は原理的に存在せず、abstention (Straddled)
//      か追加観測が唯一の正答である。
//   2. **Robust Promotion Theorem (正側)**: 真値 θ を確率 ≥ 1 − α で含む同時
//      信頼集合 C_α(D) を構成し、exact reader の像 {q(θ) : θ ∈ C_α(D)} が単一
//      クラスのときだけ昇格する。誤昇格 ⊆ 被覆失敗 なので
//      P_θ(wrong promotion) ≤ α (Lean: promotion_exclusion + robust_promotion)。
//      **これは無条件の誤昇格確率であり selective risk ではない**:
//      P(wrong | answer) ≤ α / P(answer) — 回答条件付きの保証には回答率の下限が
//      別途要る (coverage 単独から selective risk ≤ α を導くのは禁止)。
//   3. **禁止変換 22–29** — 以下の型から証明書への変換は存在しない
//      (marker 型のみ定義し、反例は v343_finite_data_nogo が厳密計算で与える):
//        22 PointEstimate            ↛ AccessibilityCertificate
//        23 MarginalIntervals        ↛ JointConfidenceRegion
//        24 GoodnessOfFitPass        ↛ NoiseModelValidity
//        25 CalibrationAtT0          ↛ ValidAtT1
//        26 MeanCrosstalk            ↛ UniformCrosstalkBound
//        27 LocalChartCoverage       ↛ GlobalGlueCoverage
//        28 ZeroHoldoutErrors        ↛ ZeroPopulationRisk
//        29 ModelConditionalCertificate ↛ ModelFreeCertificate
//
// 有限データ意味論は OCS-1.0 (paper/operational-core-spec.md §14) のスコープ外と
// 宣言済み — 本モジュールが v34.3–v34.5 で確立する定理群が OCS-2.0 系の土台になる。

use crate::ln_gamma;

// ---------------------------------------------------------------- 裁定 5 値

/// 同時信頼集合上の裁定 (PROMPT/15 §4 の 5 値 — 凍結語彙)。
/// RobustExact は「C_α(D) 全体で exact reader の読みが単一 admissible クラス」の
/// ときのみ。集合が裁定境界を跨げば Straddled — 失敗ではなく有限データからの正答。
#[derive(Clone, Debug, PartialEq)]
pub enum RobustVerdict {
    /// 信頼集合全体で単一の読み — 昇格 (誤昇格確率 ≤ α)
    RobustExact { reading: String },
    /// 複数だが既知の同値関係の単一クラス内 (例: 符号 orbit {+θ, −θ})
    EquivalenceClassOnly { class_desc: String },
    /// 信頼集合が裁定境界を跨ぐ — 強制回答の禁止
    Straddled,
    /// 観測が契約の最低量に満たない (回答資格なし)
    InsufficientObservation,
    /// データが登録観測契約の外 (支持違反・棄却されたノイズモデル)
    OutOfDomain,
}

impl RobustVerdict {
    pub fn as_str(&self) -> &'static str {
        match self {
            RobustVerdict::RobustExact { .. } => "robust_exact",
            RobustVerdict::EquivalenceClassOnly { .. } => "equivalence_class_only",
            RobustVerdict::Straddled => "straddled",
            RobustVerdict::InsufficientObservation => "insufficient_observation",
            RobustVerdict::OutOfDomain => "out_of_domain",
        }
    }
}

// ---------------------------------------------------------------- 禁止変換 22–29 の marker 型

/// 禁止変換 22: 点推定 (k/N など) — 誤差領域を持たず、それ単体に操作的資格はない。
/// AccessibilityCertificate への変換は存在しない (v343 [F6-22]: 境界近傍で
/// 誤り確率が α を大きく超える)。
pub struct PointEstimate(pub f64);

/// 禁止変換 23/27: 周辺区間の束 — 各パラメータ (1−α) の直積は joint (1−α) 領域
/// ではない (被覆はおよそ (1−α)^m)。JointConfidenceRegion / GlobalGlueCoverage
/// への変換は存在しない (v343 [F4])。
pub struct MarginalIntervals(pub Vec<(f64, f64)>);

/// 禁止変換 24: 適合度検査の通過 — 一致した積率はモデルの妥当性を意味しない
/// (同一平均の過分散モデルで被覆が壊れる — v343 [F6-24])。
pub struct GoodnessOfFitPass(pub f64);

/// 禁止変換 25: 時刻 t0 の較正 — drift 下で t1 に持ち越せない (v343 [F6-25])。
pub struct CalibrationAtT0 {
    pub interval: (f64, f64),
    pub t0: f64,
}

/// 禁止変換 26: 平均 cross-talk — worst-case 上界ではない。UniformCrosstalkBound
/// への変換は存在しない (v343 [F6-26]: 平均がバー内でも単一対の漏れが独立
/// addressability を破る)。
pub struct MeanCrosstalk(pub f64);

/// 禁止変換 28: holdout 誤り 0 — 母集団リスク 0 ではない。片側 (1−conf) 上限は
/// 1 − (1 − conf)^{1/n} (v343 [F6-28]: 0/9 → 28.3%・0/77 → 3.8%・1% には 299 件)。
pub struct ZeroHoldoutErrors(pub usize);

/// 禁止変換 29: モデル条件付き証明書 — 登録ノイズモデル (iid 等) の下でのみ有効。
/// ModelFreeCertificate への変換は存在しない (v343 [F6-29]: 同一周辺分布の相関
/// 鎖で被覆が壊れる)。
pub struct ModelConditionalCertificate {
    pub interval: (f64, f64),
    pub model: &'static str,
}

// ---------------------------------------------------------------- 厳密二項計算

/// 二項 pmf (対数経由 — n ≤ 数千で安定)
pub fn binom_pmf(n: usize, k: usize, p: f64) -> f64 {
    if p <= 0.0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    if p >= 1.0 {
        return if k == n { 1.0 } else { 0.0 };
    }
    let ln_c = ln_gamma(n as f64 + 1.0) - ln_gamma(k as f64 + 1.0) - ln_gamma((n - k) as f64 + 1.0);
    (ln_c + (k as f64) * p.ln() + ((n - k) as f64) * (1.0 - p).ln()).exp()
}

/// P(X ≤ k)
pub fn binom_cdf_le(n: usize, k: usize, p: f64) -> f64 {
    (0..=k.min(n)).map(|j| binom_pmf(n, j, p)).sum()
}

/// P(X ≥ k)
pub fn binom_sf_ge(n: usize, k: usize, p: f64) -> f64 {
    (k..=n).map(|j| binom_pmf(n, j, p)).sum()
}

/// Clopper–Pearson 厳密区間 (両側 1−α): 二項 tail の反転を二分法で解く。
/// 片側性質: P_θ(θ < lo(K)) ≤ α/2, P_θ(θ > hi(K)) ≤ α/2 — 被覆 ≥ 1−α。
pub fn cp_interval(k: usize, n: usize, alpha: f64) -> (f64, f64) {
    let half = alpha / 2.0;
    let lo = if k == 0 {
        0.0
    } else {
        // P_θ(X ≥ k) は θ に単調増加 — = half を解く
        let (mut a, mut b) = (0.0f64, 1.0f64);
        for _ in 0..200 {
            let m = 0.5 * (a + b);
            if binom_sf_ge(n, k, m) < half {
                a = m;
            } else {
                b = m;
            }
        }
        0.5 * (a + b)
    };
    let hi = if k == n {
        1.0
    } else {
        // P_θ(X ≤ k) は θ に単調減少 — = half を解く
        let (mut a, mut b) = (0.0f64, 1.0f64);
        for _ in 0..200 {
            let m = 0.5 * (a + b);
            if binom_cdf_le(n, k, m) > half {
                a = m;
            } else {
                b = m;
            }
        }
        0.5 * (a + b)
    };
    (lo, hi)
}

/// 誤り 0 件 n 回の片側母集団リスク上限: 1 − (1 − conf)^{1/n}
/// (「ZeroHoldoutErrors ↛ ZeroPopulationRisk」の定量形)
pub fn zero_error_upper_bound(n: usize, conf: f64) -> f64 {
    1.0 - (1.0 - conf).powf(1.0 / n as f64)
}

// ---------------------------------------------------------------- Robust reader

/// 二値読み (edge / no_edge) の robust reader — 登録観測契約:
/// iid Bernoulli(θ)・支持 {0,1}・最低ショット数 n_min。
/// 裁定: CP 区間が τ の上に完全 → edge / 下に完全 → no_edge / 跨ぎ → Straddled。
pub struct RobustEdgeReader {
    pub tau: f64,
    pub alpha: f64,
    pub n_min: usize,
}

impl RobustEdgeReader {
    /// 生データからの裁定 (契約検査込み)
    pub fn read(&self, shots: &[u8]) -> RobustVerdict {
        if shots.iter().any(|&s| s > 1) {
            return RobustVerdict::OutOfDomain;
        }
        self.read_counts(shots.iter().filter(|&&s| s == 1).count(), shots.len())
    }
    /// 十分統計量 (成功数 k, ショット数 n) からの裁定
    pub fn read_counts(&self, k: usize, n: usize) -> RobustVerdict {
        if n < self.n_min {
            return RobustVerdict::InsufficientObservation;
        }
        let (lo, hi) = cp_interval(k, n, self.alpha);
        if hi <= self.tau {
            RobustVerdict::RobustExact {
                reading: "no_edge".into(),
            }
        } else if lo > self.tau {
            RobustVerdict::RobustExact {
                reading: "edge".into(),
            }
        } else {
            RobustVerdict::Straddled
        }
    }
}

/// 符号 orbit reader — 観測契約が |θ| しか見ない (Bernoulli(|θ|)) とき、
/// 符号は原理的に識別不能 (TV(P₊, P₋) = 0)。回答は常に同値類
/// {+|θ|, −|θ|} — RobustExact(符号) を返す経路は存在しない。
pub struct SignOrbitReader {
    pub alpha: f64,
    pub n_min: usize,
}

impl SignOrbitReader {
    pub fn read_counts(&self, k: usize, n: usize) -> RobustVerdict {
        if n < self.n_min {
            return RobustVerdict::InsufficientObservation;
        }
        let (lo, hi) = cp_interval(k, n, self.alpha);
        RobustVerdict::EquivalenceClassOnly {
            class_desc: format!("{{+m, -m}} for |theta| in [{:.4}, {:.4}]", lo, hi),
        }
    }
}

// ---------------------------------------------------------------- 自己検査

pub fn finite_data_self_test() -> Result<(), String> {
    // 二項の正規化
    let s: f64 = (0..=20).map(|k| binom_pmf(20, k, 0.3)).sum();
    if (s - 1.0).abs() > 1e-12 {
        return Err(format!("binom_pmf 正規化 {}", s));
    }
    // CP の実被覆 ≥ 1−α (θ 走査の厳密和)
    let (n, alpha) = (25usize, 0.10f64);
    for &theta in &[0.1, 0.3, 0.5, 0.7] {
        let cov: f64 = (0..=n)
            .filter(|&k| {
                let (lo, hi) = cp_interval(k, n, alpha);
                lo <= theta && theta <= hi
            })
            .map(|k| binom_pmf(n, k, theta))
            .sum();
        if cov < 1.0 - alpha - 1e-12 {
            return Err(format!("CP 被覆 {} < 1−α at θ={}", cov, theta));
        }
    }
    // 誤り 0 上限の既知値 (n = 3 の rule of three: ≈ 63.2%? 1−0.05^{1/3} = 0.6316)
    let u3 = zero_error_upper_bound(3, 0.95);
    if (u3 - (1.0 - 0.05f64.powf(1.0 / 3.0))).abs() > 1e-15 {
        return Err("zero_error_upper_bound".into());
    }
    Ok(())
}
