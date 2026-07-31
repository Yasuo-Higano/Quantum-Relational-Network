// readout_contract — 幾何読み出しの識別可能性契約 (v31.0, PROMPT/12)
//
// 第三十一期の主題「何が状態に符号化され、何が許された観測から読め、何が原理的に
// 読めないか」を型で分離する。qrn_core (v27.2/v30.0) が「QRN 本体 vs 測定器」の
// 境界を守るのに対し、本契約は読み出しそのものの 4 軸を守る:
//
//   能力 (BridgeCapability — qrn_core の 8 タグを再利用)
//   × 状態領域 (StateDomain — Gaussian Gibbs full-rank / projector / interacting / unknown)
//   × 観測契約 (ObservationContract — 大域 C / pair RDM / 静的局所 / 局所バイアス密度応答
//                / coherent 応答 / retarded 応答)
//   × 因子分解状態 (FactorizationStatus — given / operationally-inferred / unknown)
//
// 設計の核 (qrn_core と同型の封鎖 3 段構え):
//   1. **別型 = 暗黙変換禁止**。ParentModularGenerator (状態と整合する Gaussian parent
//      K(C) = log[(I−C)C⁻¹]) と PhysicalGenerator (模型の実生成子 h) は別型 —
//      変換の唯一の門は identify_physical_generator で、GaussianityEvidence と
//      GibbsProvenance の証拠を要求する (禁止変換 8: 証拠なしの同一視は abstain)。
//   2. **RegularizedCorrelation (clamp/正則化済み) から exact 資格は構成不能** —
//      ExactWitness に RegularizedCorrelation の variant が存在せず、
//      ExactFullRankCorrelation の唯一のコンストラクタは clamp なしの
//      スペクトル証明書 (0 < λ < 1, margin δ ≥ DELTA_EXACT_FLOOR) を内部計算で要求する。
//   3. **sealed** — 状態領域・観測契約・因子分解状態のタグは外部から追加できない。
//
// 裁定 (IdentifiabilityVerdict) は「非識別セルの正しい棄却」を一級市民にする:
// Abstain / EquivalenceClassOnly は失敗ではなく正しい読み出し結果である (HOLD-7 は
// 非識別セルで無理に回答したら FAIL と採点する — PROMPT/12 §7)。
//
// **ReadoutCertificate は昇格の門ではない** — bridge law の門は qrn_core の
// BridgeLawCertificate (全能力で登録簿空) のみ。本契約の証明書は個々の読み出し実験の
// 裁定記録であり、layer: bridge の主張を core に昇格させる効力を持たない。
//
// 一次ソース: docs/uft-v31.0.md (意味論凍結) / core.schema.yml (概念登録 + 禁止変換
// 8–10)。整合は v310_readout_semantics が機械検査する (ALWAYS_RUN)。

use crate::jacobi_eigh;
use std::marker::PhantomData;

pub use crate::qrn_core::BridgeCapability;

mod sealed {
    pub trait Sealed {}
}

// ---------------------------------------------------------------- 状態領域 (4 タグ)

/// 読み出しが資格を持つ状態族。**full-rank Gaussian の結果を pure / interacting へ
/// 暗黙拡張しない** (PROMPT/12 絶対禁止) — 拡張は当該 domain での再資格のみ。
pub trait StateDomain: sealed::Sealed {
    const NAME: &'static str;
}

/// 熱的 Gaussian (0 ≺ C ≺ I) — 生成子のスペクトル情報を保持する領域
pub enum GaussianGibbsFullRank {}
/// 純粋 Gaussian (C² = C) — 静的核は sign 情報まで (v29.5 P6/693 の領域)
pub enum GaussianProjector {}
/// 相互作用フェルミオン — one-body C は状態を決定しない
pub enum InteractingFermion {}
/// 出自不明 — Gaussian 性・Gibbs 出自の証拠なし
pub enum UnknownStateDomain {}

macro_rules! impl_tag {
    ($tr:ident, $t:ty, $name:literal) => {
        impl sealed::Sealed for $t {}
        impl $tr for $t {
            const NAME: &'static str = $name;
        }
    };
}
impl_tag!(StateDomain, GaussianGibbsFullRank, "gaussian_gibbs_full_rank");
impl_tag!(StateDomain, GaussianProjector, "gaussian_projector");
impl_tag!(StateDomain, InteractingFermion, "interacting_fermion");
impl_tag!(StateDomain, UnknownStateDomain, "unknown_state_domain");

// ---------------------------------------------------------------- 観測契約 (6 タグ)

/// 読み出しに許された観測資源。強さの半順序はここでは定義しない —
/// v31.3 (観測予算 hierarchy) が同一 hidden generator に対して実測する。
pub trait ObservationContract: sealed::Sealed {
    const NAME: &'static str;
}

/// 大域 one-body 相関行列 C 全体 (oracle ceiling の観測)
pub enum GlobalOneBodyCorrelation {}
/// 二ノード reduced 状態 (現行 B2 の観測 — 環境で renormalize される)
pub enum PairReducedStates {}
/// 静的局所観測量 (密度・局所相関のみ)
pub enum StaticLocalObservables {}
/// 局所バイアス probe 後の密度時系列 (v31.2 の密度曲率則)
pub enum LocalBiasDensityResponse {}
/// 局所バイアス probe 後の full coherence 応答 (v31.2 の commutator 則)
pub enum CoherentLocalResponse {}
/// retarded 応答核 (B6 系)
pub enum RetardedResponse {}

impl_tag!(ObservationContract, GlobalOneBodyCorrelation, "global_one_body_correlation");
impl_tag!(ObservationContract, PairReducedStates, "pair_reduced_states");
impl_tag!(ObservationContract, StaticLocalObservables, "static_local_observables");
impl_tag!(ObservationContract, LocalBiasDensityResponse, "local_bias_density_response");
impl_tag!(ObservationContract, CoherentLocalResponse, "coherent_local_response");
impl_tag!(ObservationContract, RetardedResponse, "retarded_response");

// ---------------------------------------------------------------- 因子分解状態 (3 タグ)

/// ノード因子分解の出自。現行の全 bridge 成果は GivenNodeFactorization (入力) —
/// 状態・操作からの選定 (RelationalDecompositionGoal) は未構成 (v29.5 [C5])。
pub trait FactorizationStatus: sealed::Sealed {
    const NAME: &'static str;
}

/// 因子分解を入力として与えられている (観測資源であり読み出しではない)
pub enum GivenNodeFactorization {}
/// 観測・操作だけから推定した因子分解 (OperationalAlgebra が必要 — v31.4)
pub enum OperationallyInferredFactorization {}
/// 因子分解不明 (mode 基底の任意性の下では静的一意選定は不可能 — v31.4 no-go)
pub enum UnknownFactorization {}

impl_tag!(FactorizationStatus, GivenNodeFactorization, "given_node_factorization");
impl_tag!(
    FactorizationStatus,
    OperationallyInferredFactorization,
    "operationally_inferred_factorization"
);
impl_tag!(FactorizationStatus, UnknownFactorization, "unknown_factorization");

// ---------------------------------------------------------------- 裁定と棄却理由

/// 棄却理由 — 8 種で凍結 (PROMPT/12 §v31.0)。棄却は失敗ではなく正しい読み出し結果。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AbstainReason {
    /// 固有値が {0, 1} に触れる (projector 系 — logit が発散)
    RankDeficient,
    /// spectral margin δ が f64 の exact lane 床を下回る
    IllConditioned,
    /// Gaussian 性の証拠がない (Unknown のまま logit を物理生成子と読まない)
    GaussianityUnverified,
    /// Gibbs 出自 (β, μ の意味) の証拠がない
    GibbsProvenanceMissing,
    /// 因子分解が不明 (mode 基底の任意性で静的一意選定は不可)
    UnknownFactorization,
    /// 観測契約が要求精度に足りない
    InsufficientObservation,
    /// 非 Gaussian 領域 (Wick 残差がバー超え)
    NonGaussianDomain,
    /// 局所 Gram 行列が rank 欠損 (無条件擬似逆は禁止 — support 制限か棄却)
    RankDeficientLocalGram,
}

impl AbstainReason {
    pub const ALL: [AbstainReason; 8] = [
        AbstainReason::RankDeficient,
        AbstainReason::IllConditioned,
        AbstainReason::GaussianityUnverified,
        AbstainReason::GibbsProvenanceMissing,
        AbstainReason::UnknownFactorization,
        AbstainReason::InsufficientObservation,
        AbstainReason::NonGaussianDomain,
        AbstainReason::RankDeficientLocalGram,
    ];
    pub fn as_str(self) -> &'static str {
        match self {
            AbstainReason::RankDeficient => "rank_deficient",
            AbstainReason::IllConditioned => "ill_conditioned",
            AbstainReason::GaussianityUnverified => "gaussianity_unverified",
            AbstainReason::GibbsProvenanceMissing => "gibbs_provenance_missing",
            AbstainReason::UnknownFactorization => "unknown_factorization",
            AbstainReason::InsufficientObservation => "insufficient_observation",
            AbstainReason::NonGaussianDomain => "non_gaussian_domain",
            AbstainReason::RankDeficientLocalGram => "rank_deficient_local_gram",
        }
    }
}

/// 識別可能性の裁定 — 5 値。「読めた」だけでなく「どこまで読めたか」を型で言わせる。
#[derive(Clone, Debug, PartialEq)]
pub enum IdentifiabilityVerdict {
    /// ノード内基底 (block-local unitary) の gauge を除いて厳密復元
    ExactUpToGauge,
    /// さらに正の大域スケール (β 未知) を除いて厳密復元
    ExactUpToGlobalScale,
    /// 正則化つき安定推定 (条件数上界を明示)
    StableEstimate { condition_bound: f64 },
    /// 同値類までしか識別できない (例: 半充填 projector の sign(A) 類 — v29.5)
    EquivalenceClassOnly { class_desc: String },
    /// 棄却 (理由つき) — 非識別セルでの正しい読み出し結果
    Abstain(AbstainReason),
}

impl IdentifiabilityVerdict {
    pub fn as_str(&self) -> &'static str {
        match self {
            IdentifiabilityVerdict::ExactUpToGauge => "exact_up_to_gauge",
            IdentifiabilityVerdict::ExactUpToGlobalScale => "exact_up_to_global_scale",
            IdentifiabilityVerdict::StableEstimate { .. } => "stable_estimate",
            IdentifiabilityVerdict::EquivalenceClassOnly { .. } => "equivalence_class_only",
            IdentifiabilityVerdict::Abstain(_) => "abstain",
        }
    }
}

// ---------------------------------------------------------------- 生成子 2 型 (禁止変換 8)

/// 完全な相関行列 C と整合する Gaussian parent K(C) = log[(I−C)C⁻¹]。
/// **物理生成子ではない** — 状態が (i) Gaussian で (ii) 物理生成子の Gibbs 状態で
/// (iii) β, μ の意味が与えられて初めて h と同定できる (門は identify_physical_generator)。
/// エルミート行列を (re 対称, im 反対称) の実部虚部で保持する。
#[derive(Clone, Debug, PartialEq)]
pub struct ParentModularGenerator {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub n: usize,
}

/// 模型の実時間発展生成子 h (一体エルミート)。読み出し側では同定対象であり、
/// ParentModularGenerator からの impl From は存在しない (禁止変換 8)。
#[derive(Clone, Debug, PartialEq)]
pub struct PhysicalGenerator {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub n: usize,
}

/// Gaussian 性の証拠 — 証拠なしの logit 読み出しは abstain (v31.5 で witness 実装)
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GaussianityEvidence {
    /// 構成により Gaussian (自由フェルミオン solver の出力)
    ByConstruction,
    /// Wick 残差の実測上界 (bar は走行前に凍結すること)
    WickResidualBound { residual: f64, bar: f64 },
    /// 証拠なし
    Unknown,
}

/// Gibbs 出自の証拠 — C = (I + e^{β(h−μI)})⁻¹ の β, μ の意味
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GibbsProvenance {
    /// β, μ が物理的に与えられている
    KnownBetaMu { beta: f64, mu: f64 },
    /// β > 0 だが値は未知 — 正の大域スケール同値類まで
    BetaUnknownPositive,
    /// 出自の証拠なし
    Missing,
}

/// 物理生成子の読み出し結果 — どこまで確定したかを型で区別する
#[derive(Clone, Debug, PartialEq)]
pub enum PhysicalGeneratorReading {
    /// β, μ 既知: h = K/β + μI を一意復元
    Exact(PhysicalGenerator),
    /// β 未知: 正の大域スケールを除いて確定 (代表元 = K そのもの)。
    /// μ の不定性は一様対角シフトのみで空間隣接に影響しない。
    UpToPositiveScaleAndShift(PhysicalGenerator),
}

/// ParentModularGenerator → PhysicalGenerator の**唯一の門** (禁止変換 8 の解錠条件)。
/// Gaussian 性と Gibbs 出自の証拠がなければ棄却を返す — 無条件変換は書けない。
pub fn identify_physical_generator(
    parent: &ParentModularGenerator,
    gaussianity: GaussianityEvidence,
    provenance: GibbsProvenance,
) -> Result<PhysicalGeneratorReading, AbstainReason> {
    match gaussianity {
        GaussianityEvidence::Unknown => return Err(AbstainReason::GaussianityUnverified),
        GaussianityEvidence::WickResidualBound { residual, bar } => {
            if !(residual.is_finite() && bar.is_finite()) || residual > bar {
                return Err(AbstainReason::NonGaussianDomain);
            }
        }
        GaussianityEvidence::ByConstruction => {}
    }
    match provenance {
        GibbsProvenance::Missing => Err(AbstainReason::GibbsProvenanceMissing),
        GibbsProvenance::KnownBetaMu { beta, mu } => {
            assert!(beta > 0.0 && beta.is_finite(), "β は正の有限値であること");
            let n = parent.n;
            let mut re = vec![0.0; n * n];
            let mut im = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    re[i * n + j] = parent.re[i * n + j] / beta + if i == j { mu } else { 0.0 };
                    im[i * n + j] = parent.im[i * n + j] / beta;
                }
            }
            Ok(PhysicalGeneratorReading::Exact(PhysicalGenerator { re, im, n }))
        }
        GibbsProvenance::BetaUnknownPositive => {
            Ok(PhysicalGeneratorReading::UpToPositiveScaleAndShift(PhysicalGenerator {
                re: parent.re.clone(),
                im: parent.im.clone(),
                n: parent.n,
            }))
        }
    }
}

// ---------------------------------------------------------------- 相関 2 型 (禁止変換 9)

/// exact lane の f64 床: spectral margin δ = min(λ_min, 1−λ_max) がこれ未満なら
/// exact を名乗れない (モジュラー核 κ = ln((1−c)/c) の f64 分解能床 ~1e-14 の系譜 —
/// CLAUDE.md 既知の落とし穴)。これは物理ではなく器械の契約。
pub const DELTA_EXACT_FLOOR: f64 = 1e-13;

/// 厳密 full-rank 相関 — **唯一のコンストラクタが clamp なしのスペクトル証明書を
/// 内部計算で要求する** (0 < λ < 1 かつ δ ≥ DELTA_EXACT_FLOOR)。
/// RegularizedCorrelation からの変換は存在しない (禁止変換 9)。
#[derive(Clone, Debug)]
pub struct ExactFullRankCorrelation {
    c_re: Vec<f64>,
    c_im: Vec<f64>,
    n: usize,
    delta: f64,
}

impl ExactFullRankCorrelation {
    /// 実対称 C の資格審査 (clamp なし)。失敗は棄却理由を返す:
    /// λ ∉ (0,1) → RankDeficient / 0 < δ < 床 → IllConditioned。
    pub fn certify_real(c: &[f64], n: usize) -> Result<Self, AbstainReason> {
        assert_eq!(c.len(), n * n, "C の形が n×n でない");
        for i in 0..n {
            for j in (i + 1)..n {
                assert!(
                    (c[i * n + j] - c[j * n + i]).abs() <= 1e-12,
                    "C が対称でない (呼び出し側のバグ)"
                );
            }
        }
        let (evals, _) = jacobi_eigh(c, n);
        Self::from_margin(c.to_vec(), vec![0.0; n * n], n, &evals)
    }

    /// エルミート C = A + iB の資格審査 — 2n×2n 実対称埋め込み [[A, −B], [B, A]]
    /// (固有値は各 2 重化されるが margin は不変)。
    pub fn certify_herm(re: &[f64], im: &[f64], n: usize) -> Result<Self, AbstainReason> {
        assert_eq!(re.len(), n * n);
        assert_eq!(im.len(), n * n);
        for i in 0..n {
            assert!(im[i * n + i].abs() <= 1e-12, "エルミートなら対角虚部は 0");
            for j in 0..n {
                assert!((re[i * n + j] - re[j * n + i]).abs() <= 1e-12, "re は対称");
                assert!((im[i * n + j] + im[j * n + i]).abs() <= 1e-12, "im は反対称");
            }
        }
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
        let (evals, _) = jacobi_eigh(&big, m);
        Self::from_margin(re.to_vec(), im.to_vec(), n, &evals)
    }

    fn from_margin(
        c_re: Vec<f64>,
        c_im: Vec<f64>,
        n: usize,
        evals: &[f64],
    ) -> Result<Self, AbstainReason> {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for &e in evals {
            lo = lo.min(e);
            hi = hi.max(e);
        }
        let delta = lo.min(1.0 - hi);
        if !(lo > 0.0 && hi < 1.0) {
            return Err(AbstainReason::RankDeficient);
        }
        if delta < DELTA_EXACT_FLOOR {
            return Err(AbstainReason::IllConditioned);
        }
        Ok(ExactFullRankCorrelation { c_re, c_im, n, delta })
    }

    pub fn c_re(&self) -> &[f64] {
        &self.c_re
    }
    pub fn c_im(&self) -> &[f64] {
        &self.c_im
    }
    pub fn n(&self) -> usize {
        self.n
    }
    /// スペクトル余裕 δ = min(λ_min, 1 − λ_max) — 条件数上界 1/(δ(1−δ)) の証人
    pub fn spectral_margin(&self) -> f64 {
        self.delta
    }
}

/// clamp / 正則化済み相関 — estimate lane 専用。数値推定としては正当だが、
/// **exact 識別可能性証明書の根拠にはできない** (禁止変換 9)。
#[derive(Clone, Debug)]
pub struct RegularizedCorrelation {
    pub c_re: Vec<f64>,
    pub c_im: Vec<f64>,
    pub n: usize,
    /// 施した clamp の ε (固有値を [ε, 1−ε] に切った等)
    pub clamp_eps: f64,
}

// ---------------------------------------------------------------- logit 2 経路の型分離

/// 全系 C に logit をかけてから block (i,j) を抽出したもの —
/// GlobalOneBodyCorrelation 観測の親生成子 block。一般に ReducedModularBlock と
/// 一致しない (f(P C P) ≠ P f(C) P) — v31.1 が反例を常設検査する。
#[derive(Clone, Debug, PartialEq)]
pub struct GlobalPhysicalParentBlock {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub ni: usize,
    pub nj: usize,
}

/// 二ノード RDM に logit をかけたもの (現行 B2) — 環境によって renormalize された
/// reduced modular coupling。観測予算の両端の一方として保持する。
#[derive(Clone, Debug, PartialEq)]
pub struct ReducedModularBlock {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub ni: usize,
    pub nj: usize,
}

// ---------------------------------------------------------------- patch 2 型 (循環防止)

/// 診断専用: **真の隣接半径**で選んだ patch。OperationalPatch への変換は存在しない —
/// 真の幾何で patch を選んで読み出しに使うと循環する (PROMPT/12 絶対禁止)。
#[derive(Clone, Debug)]
pub struct OraclePatch {
    pub center: usize,
    pub members: Vec<usize>,
}

/// 観測だけから構築した patch (B3 等の観測 graph からの拡張)。provenance に
/// 構築に使った観測契約を記録する。
#[derive(Clone, Debug)]
pub struct OperationalPatch {
    pub center: usize,
    pub members: Vec<usize>,
    pub provenance: &'static str,
}

// ---------------------------------------------------------------- 支持証明書・操作代数

/// 局所 Gram 行列が rank 欠損のときに返す支持証明書 — 無条件の擬似逆行列を禁止し、
/// support 制限つき読み出しか棄却 (RankDeficientLocalGram) を強制する (v31.4)。
#[derive(Clone, Debug)]
pub struct ObservableSupportCertificate {
    pub rank: usize,
    pub threshold: f64,
    pub nullspace_dim: usize,
}

/// FactorizationGivenObservables が成立するために必要な操作的構造 (v31.4 で実装)。
/// 静的状態だけでは因子分解は一意選定できない — 選定には準備・介入・測定・両立性の
/// 資源が要る (v29.5 [C5] の空隙を埋める側の定義)。
#[derive(Clone, Debug, Default)]
pub struct OperationalAlgebra {
    pub preparations: Vec<&'static str>,
    pub interventions: Vec<&'static str>,
    pub measurements: Vec<&'static str>,
    pub compatibility: Vec<&'static str>,
}

// ---------------------------------------------------------------- 読み出し証明書

/// exact 資格の証人 — **RegularizedCorrelation の variant は存在しない** (禁止変換 9)。
pub enum ExactWitness<'a> {
    /// 静的逆問題: clamp なしの full-rank スペクトル証明書
    FullRankCorrelation(&'a ExactFullRankCorrelation),
    /// 応答恒等式: 代数的恒等式の照合残差 (バーは走行前に凍結)
    AlgebraicIdentity { residual: f64, bar: f64 },
}

/// 読み出し実験 1 件の裁定記録 — 能力 × 状態領域 × 観測契約 × 因子分解状態を型に持つ。
/// **昇格の門ではない** (門は qrn_core::BridgeLawCertificate のみ)。verdict は私有 —
/// コンストラクタ 4 本 (exact / stable_estimate / equivalence_class / abstain) だけが
/// 裁定を書ける。
pub struct ReadoutCertificate<Cap, D, O, F>
where
    Cap: BridgeCapability,
    D: StateDomain,
    O: ObservationContract,
    F: FactorizationStatus,
{
    claim_id: &'static str,
    verdict: IdentifiabilityVerdict,
    _marker: PhantomData<(Cap, D, O, F)>,
}

impl<Cap, D, O, F> ReadoutCertificate<Cap, D, O, F>
where
    Cap: BridgeCapability,
    D: StateDomain,
    O: ObservationContract,
    F: FactorizationStatus,
{
    /// exact 裁定 — 証人が必要 (witness の検証に失敗したら Err)。
    /// up_to_scale = true なら ExactUpToGlobalScale (β 未知の同値類)。
    pub fn exact(
        claim_id: &'static str,
        witness: &ExactWitness,
        up_to_scale: bool,
    ) -> Result<Self, AbstainReason> {
        match witness {
            ExactWitness::FullRankCorrelation(cert) => {
                // 証明書の margin は構成時に検査済み — ここでは実在のみ要求
                debug_assert!(cert.spectral_margin() >= DELTA_EXACT_FLOOR);
            }
            ExactWitness::AlgebraicIdentity { residual, bar } => {
                if !(residual.is_finite() && *residual <= *bar) {
                    return Err(AbstainReason::InsufficientObservation);
                }
            }
        }
        Ok(ReadoutCertificate {
            claim_id,
            verdict: if up_to_scale {
                IdentifiabilityVerdict::ExactUpToGlobalScale
            } else {
                IdentifiabilityVerdict::ExactUpToGauge
            },
            _marker: PhantomData,
        })
    }

    /// 正則化つき安定推定 — RegularizedCorrelation を証人に取る (exact は名乗れない)
    pub fn stable_estimate(
        claim_id: &'static str,
        witness: &RegularizedCorrelation,
        condition_bound: f64,
    ) -> Self {
        let _ = witness.clamp_eps;
        ReadoutCertificate {
            claim_id,
            verdict: IdentifiabilityVerdict::StableEstimate { condition_bound },
            _marker: PhantomData,
        }
    }

    /// 同値類までの識別 (例: projector lane の sign(A) 類)
    pub fn equivalence_class(claim_id: &'static str, class_desc: String) -> Self {
        ReadoutCertificate {
            claim_id,
            verdict: IdentifiabilityVerdict::EquivalenceClassOnly { class_desc },
            _marker: PhantomData,
        }
    }

    /// 棄却 — 非識別セルでの正しい読み出し結果
    pub fn abstain(claim_id: &'static str, reason: AbstainReason) -> Self {
        ReadoutCertificate {
            claim_id,
            verdict: IdentifiabilityVerdict::Abstain(reason),
            _marker: PhantomData,
        }
    }

    pub fn claim_id(&self) -> &'static str {
        self.claim_id
    }
    pub fn verdict(&self) -> &IdentifiabilityVerdict {
        &self.verdict
    }
    pub fn capability(&self) -> &'static str {
        Cap::NAME
    }
    pub fn state_domain(&self) -> &'static str {
        D::NAME
    }
    pub fn observation(&self) -> &'static str {
        O::NAME
    }
    pub fn factorization(&self) -> &'static str {
        F::NAME
    }
}

// ---------------------------------------------------------------- 自己検査

/// readout_contract の不変条件 (v310_readout_semantics が呼ぶ):
/// タグの命名・封鎖経路・門の棄却挙動・exact lane の床。
pub fn readout_contract_self_test() -> Result<(), String> {
    // 1. タグ名の非空・重複なし (4 + 6 + 3 = 13)
    let names = [
        GaussianGibbsFullRank::NAME,
        GaussianProjector::NAME,
        InteractingFermion::NAME,
        UnknownStateDomain::NAME,
        GlobalOneBodyCorrelation::NAME,
        PairReducedStates::NAME,
        StaticLocalObservables::NAME,
        LocalBiasDensityResponse::NAME,
        CoherentLocalResponse::NAME,
        RetardedResponse::NAME,
        GivenNodeFactorization::NAME,
        OperationallyInferredFactorization::NAME,
        UnknownFactorization::NAME,
    ];
    for (i, a) in names.iter().enumerate() {
        if a.is_empty() {
            return Err("空のタグ名".into());
        }
        for b in names.iter().skip(i + 1) {
            if a == b {
                return Err(format!("タグ名の重複: {}", a));
            }
        }
    }
    // 2. タグは居住不能 (空 enum)
    if std::mem::size_of::<GaussianGibbsFullRank>() != 0
        || std::mem::size_of::<UnknownFactorization>() != 0
    {
        return Err("タグ型が居住可能になっている".into());
    }
    // 3. exact lane: 良条件 C は資格・projector は RankDeficient・薄い margin は IllConditioned
    let good = ExactFullRankCorrelation::certify_real(&[0.3, 0.0, 0.0, 0.7], 2);
    match &good {
        Ok(c) if (c.spectral_margin() - 0.3).abs() < 1e-12 => {}
        _ => return Err("良条件 C (λ = 0.3, 0.7) の資格に失敗".into()),
    }
    match ExactFullRankCorrelation::certify_real(&[1.0, 0.0, 0.0, 0.0], 2) {
        Err(AbstainReason::RankDeficient) => {}
        _ => return Err("projector C が RankDeficient にならない".into()),
    }
    match ExactFullRankCorrelation::certify_real(&[1e-15, 0.0, 0.0, 0.5], 2) {
        Err(AbstainReason::IllConditioned) => {}
        _ => return Err("margin 1e-15 が IllConditioned にならない".into()),
    }
    // 4. 門: 証拠なしの変換は棄却・証拠ありは β, μ を正しく解く
    let k = ParentModularGenerator {
        re: vec![2.0, -1.0, -1.0, 4.0],
        im: vec![0.0; 4],
        n: 2,
    };
    match identify_physical_generator(&k, GaussianityEvidence::Unknown, GibbsProvenance::Missing) {
        Err(AbstainReason::GaussianityUnverified) => {}
        _ => return Err("Gaussian 性未検証で棄却しない".into()),
    }
    match identify_physical_generator(
        &k,
        GaussianityEvidence::WickResidualBound { residual: 1e-2, bar: 1e-8 },
        GibbsProvenance::KnownBetaMu { beta: 1.0, mu: 0.0 },
    ) {
        Err(AbstainReason::NonGaussianDomain) => {}
        _ => return Err("Wick 残差バー超えで棄却しない".into()),
    }
    match identify_physical_generator(
        &k,
        GaussianityEvidence::ByConstruction,
        GibbsProvenance::Missing,
    ) {
        Err(AbstainReason::GibbsProvenanceMissing) => {}
        _ => return Err("Gibbs 出自なしで棄却しない".into()),
    }
    match identify_physical_generator(
        &k,
        GaussianityEvidence::ByConstruction,
        GibbsProvenance::KnownBetaMu { beta: 2.0, mu: 0.5 },
    ) {
        Ok(PhysicalGeneratorReading::Exact(h)) => {
            // K = β(h − μI) → h = K/β + μI: h = [[1.5, −0.5], [−0.5, 2.5]]
            let want = [1.5, -0.5, -0.5, 2.5];
            for (a, b) in h.re.iter().zip(want.iter()) {
                if (a - b).abs() > 1e-12 {
                    return Err("β, μ 既知の復元が誤り".into());
                }
            }
        }
        _ => return Err("証拠ありの門が開かない".into()),
    }
    match identify_physical_generator(
        &k,
        GaussianityEvidence::ByConstruction,
        GibbsProvenance::BetaUnknownPositive,
    ) {
        Ok(PhysicalGeneratorReading::UpToPositiveScaleAndShift(_)) => {}
        _ => return Err("β 未知がスケール同値類にならない".into()),
    }
    // 5. AbstainReason は 8 種 (名は一意)
    for (i, a) in AbstainReason::ALL.iter().enumerate() {
        for b in AbstainReason::ALL.iter().skip(i + 1) {
            if a.as_str() == b.as_str() {
                return Err("棄却理由の名の重複".into());
            }
        }
    }
    // 6. 証明書: 恒等式証人の残差バー超えは exact を拒否
    type RcTest = ReadoutCertificate<
        crate::qrn_core::SpatialMetricUpToGlobalScale,
        GaussianGibbsFullRank,
        GlobalOneBodyCorrelation,
        GivenNodeFactorization,
    >;
    match RcTest::exact(
        "TEST",
        &ExactWitness::AlgebraicIdentity { residual: 1e-3, bar: 1e-12 },
        false,
    ) {
        Err(AbstainReason::InsufficientObservation) => {}
        _ => return Err("残差バー超えの exact が拒否されない".into()),
    }
    let ok_cert = RcTest::exact(
        "TEST",
        &ExactWitness::FullRankCorrelation(good.as_ref().unwrap()),
        true,
    );
    match ok_cert {
        Ok(c) if *c.verdict() == IdentifiabilityVerdict::ExactUpToGlobalScale => {}
        _ => return Err("full-rank 証人の exact 構成に失敗".into()),
    }
    Ok(())
}
