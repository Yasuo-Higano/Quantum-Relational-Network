// qrn_core — 型付き Core Contract (v27.2, PROMPT/10 §4)
//
// 「何が QRN 本体で、何が既知 QFT の測定器か」を型で混同不能にする。
// 一次ソース: core.schema.yml (層分類) / docs/qrn-terminology.md (型境界) /
// docs/qrn-core-v1-spec.md (状態表示)。整合は v272_core_contract が機械検査し、
// 禁止 impl From の不在は v271_core_audit [S12] が常設監査する。
//
// 設計の核:
//   1. 別型 = 暗黙変換の禁止 (RegulatorSiteId → RelationalNodeId 等は From を書かない)。
//   2. 変換の唯一の門は BridgeLawCertificate — 登録簿 REGISTERED_BRIDGE_LAWS は
//      現在【空】(bridge law は一つも確立していない)。register() は常に None。
//   3. 禁止昇格の到達先 (QrnEvidence 等) は【居住不能型】(空 enum) — 値を構成する
//      コード自体が書けない。bridge law が確立したときに初めて variant を追加する。
//   4. 未定義部分は ContractStatus::{Conjectural, Undefined} を明示し、暗黙の
//      fallback や推測値を入れない。

// ---------------------------------------------------------------- 状態表示

/// 契約の状態 — 「Core 完成」の水増しを型で防ぐ (spec §7 の完了条件)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ContractStatus {
    /// 定義済み (実装が toy に限られる場合は note に明示)
    Defined,
    /// 固有原理は未定義 — 模型族だけがある
    ModelFamilyOnly,
    /// 仮説段階 (toy 実演のみ・bridge law 未確立)
    Conjectural,
    /// 支持する結果 0 件
    Unsupported,
    /// 未構成
    Undefined,
}

impl ContractStatus {
    pub fn as_str(self) -> &'static str {
        match self {
            ContractStatus::Defined => "defined",
            ContractStatus::ModelFamilyOnly => "model_family_only",
            ContractStatus::Conjectural => "conjectural",
            ContractStatus::Unsupported => "unsupported",
            ContractStatus::Undefined => "undefined",
        }
    }
}

/// QRN-Core v1 の構成要素 (PROMPT/10 §4 の 9 + GravityBridge)
pub struct ContractComponent {
    pub name: &'static str,
    pub layer: &'static str,
    pub status: ContractStatus,
    pub note: &'static str,
}

/// QRN-Core v1 の現在の状態 (docs/qrn-core-v1-spec.md §7 / core.schema.yml と
/// v272_core_contract が照合する)。**空の dynamics/bridge を Defined と書かないこと。**
pub const QRN_CORE_V1: [ContractComponent; 10] = [
    ContractComponent {
        name: "StateSpace",
        layer: "core",
        status: ContractStatus::Defined,
        note: "有限次元 Hilbert 空間の族と状態。実装は GaussianFermionState (toy, C3) のみ",
    },
    ContractComponent {
        name: "ObservableAlgebra",
        layer: "core",
        status: ContractStatus::Defined,
        note: "局所演算子代数。実装は相関行列の部分行列 (ガウス系限定)",
    },
    ContractComponent {
        name: "RelationalDecomposition",
        layer: "core",
        status: ContractStatus::Defined,
        note: "テンソル分解は入力でなく読み出し (v11.4)。型は RelationalNodeId",
    },
    ContractComponent {
        name: "EvolutionLaw",
        layer: "dynamics",
        status: ContractStatus::ModelFamilyOnly,
        note: "QRN 固有の発展原理は未定義 — GaussianToyModel / ConstrainedToyCoreV2 の模型族のみ",
    },
    ContractComponent {
        name: "ConstraintAlgebra",
        layer: "dynamics",
        status: ContractStatus::ModelFamilyOnly,
        note: "拘束 core 族 (v20–v22) の個別実装のみ — 一般代数は未定義",
    },
    ContractComponent {
        name: "GeometryBridge",
        layer: "bridge",
        status: ContractStatus::Conjectural,
        note: "B1–B4 候補の toy 実演のみ — bridge law (regulator 非依存の一意な読み出し) は未確立",
    },
    ContractComponent {
        name: "CausalBridge",
        layer: "bridge",
        status: ContractStatus::Conjectural,
        note: "光円錐の創発 (v1.1, toy)。Lorentzian 因果構造の再構成は未",
    },
    ContractComponent {
        name: "ClockBridge",
        layer: "bridge",
        status: ContractStatus::Undefined,
        note:
            "時計の読み出しは未構成 (ModularParameter/OperationalClockReading/ProperTime は型のみ)",
    },
    ContractComponent {
        name: "MatterBridge",
        layer: "bridge",
        status: ContractStatus::Conjectural,
        note: "k_F・密度の読み出し (v6.7, toy)",
    },
    ContractComponent {
        name: "GravityBridge",
        layer: "bridge",
        status: ContractStatus::Unsupported,
        note: "支持する結果 0 件 — v27.0 fork: graviton pole なし・c₂ は regulator 量",
    },
];

// ---------------------------------------------------------------- ノード 3 型

/// QRN 存在論上の関係要素 (layer: core — 実在の仮説)。
/// RegulatorSiteId からの変換は BridgeLawCertificate なしには存在しない。
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RelationalNodeId(pub u64);

/// 数値離散化上の格子点 (layer: adapter — regulator)。v27.0 fork で存在論から分離確定。
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct RegulatorSiteId(pub u64);

/// 有効理論上の点 (layer: adapter — 連続極限の座標)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ContinuumPoint(pub [f64; 4]);

// ---------------------------------------------------------------- 計量 2 型

/// 外部から入力する source 計量 h_μν (layer: adapter — v26.2–v27.0 の全応答測定の対象)
#[derive(Clone, Debug, PartialEq)]
pub struct ExternalMetricSource {
    /// source の摂動振幅 (BOND-A 規約での h)
    pub amplitude: f64,
    /// source の説明 (チャネル・運動量など — 器械台帳 instruments.yml の語彙)
    pub description: &'static str,
}

/// 状態から読み出される計量の候補 (layer: bridge — 現在 unsupported)。
/// ExternalMetricSource からの変換は BridgeLawCertificate なしには存在しない。
/// 構成には bridge law の成功条件 7 項 (spec §5) を満たす読み出しが要る。
#[derive(Clone, Debug, PartialEq)]
pub struct EmergentMetricCandidate {
    pub components: Vec<f64>,
    pub bridge_claim_id: &'static str,
}

// ---------------------------------------------------------------- 時間 4 型

/// シミュレーションの外部発展パラメータ t (layer: dynamics — 器械の入力)。
/// **A1 の創発時間ではない** — ProperTime への変換は BridgeLawCertificate が要る。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EvolutionParameter(pub f64);

/// モジュラー流のパラメータ s (layer: bridge — 熱時間仮説の変数)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModularParameter(pub f64);

/// 部分系時計の読み (layer: bridge — 操作的時間。未構成)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OperationalClockReading(pub f64);

/// 固有時間の候補 (layer: bridge — 未構成)
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ProperTime(pub f64);

// ---------------------------------------------------------------- 証拠の型 (居住不能な昇格先)

/// 証拠クラス — 禁止昇格 (qrn-terminology.md §2.4) の**到達先を居住不能型にする**。
/// 空 enum は値を構成するコード自体が書けない: 「QRN の証拠」「自然の観測量の的中」
/// 「独立外部再現」は現在、型レベルで 0 件である。bridge law・自然の的中・外部再現が
/// 成立したとき、対応する variant (certificate を運ぶ) を追加する改訂が唯一の門。
pub mod evidence {
    /// 較正・器械認証の証拠 (構成可能 — 現在の全 PASS はここに住む)
    #[derive(Clone, Copy, Debug)]
    pub struct CalibrationEvidence {
        pub source: &'static str,
    }
    /// 計算実験内の凍結予言の的中 (構成可能 — PRED-013〜018)
    #[derive(Clone, Copy, Debug)]
    pub struct InternalPrediction {
        pub pred_id: &'static str,
    }
    /// QRN (関係的量子ネットワーク存在論) の証拠 — **居住不能**
    pub enum QrnEvidence {}
    /// 自然の観測量の的中 (future_observation × scored-hit) — **居住不能**
    pub enum NaturalObservation {}
    /// 独立外部再現 (replications.yml の 6 条件を満たすもの) — **居住不能**
    pub enum IndependentReplication {}
    /// regulator 非依存が証明された普遍量 (bridge 側) — **居住不能**
    pub enum UniversalQuantityForBridge {}
}

// ---------------------------------------------------------------- Bridge law の門 (v30.0: 能力別)

/// 確立済み bridge law の登録簿の**統合ビュー** — 全能力の和集合で、**現在は空**
/// (qrn-core-v1-spec.md §5: 成功条件 7 項を満たす bridge law は一つもない)。
/// v30.0 で登録簿は能力ごと (BridgeCapability::REGISTERED) に分割された —
/// 空間位相の証拠と proper time の証拠を単一 claim id で混同できた旧設計の是正
/// (PROMPT/11)。追加は証拠 (claims.yml の id) を伴う改訂としてのみ行う。
/// v272_core_contract が全能力で空であることを検査する。
pub const REGISTERED_BRIDGE_LAWS: &[&str] = &[];

mod sealed_cap {
    pub trait Sealed {}
}

/// bridge 能力タグ (v30.0, PROMPT/11) — 証明書は能力ごとに別型で、ある能力の
/// 証拠を別の能力に流用する経路が型レベルで存在しない。実装は本モジュール内に
/// 封印 (sealed) — 外部 crate/モジュールから能力を追加できない。
pub trait BridgeCapability: sealed_cap::Sealed {
    const NAME: &'static str;
    /// この能力で確立済みの bridge law (claims.yml の id) — 現在は全能力で空
    const REGISTERED: &'static [&'static str];
}

// 8 能力タグ (空 enum — 値は存在せず、型パラメータとしてのみ使う)。
// 現行証拠で将来到達しうる上限は SpatialTopologyGivenFactorization /
// SpatialMetricUpToGlobalScale / CausalOrderGivenExternalClock の 3 つ (PROMPT/11)。
pub enum FactorizationGivenObservables {}
pub enum SpatialTopologyGivenFactorization {}
pub enum SpatialMetricUpToGlobalScale {}
pub enum CausalOrderGivenExternalClock {}
pub enum ConformalLorentzianStructure {}
pub enum VolumeMeasure {}
pub enum ClockCalibration {}
pub enum FullLorentzianMetric {}

macro_rules! impl_capability {
    ($t:ty, $name:literal) => {
        impl sealed_cap::Sealed for $t {}
        impl BridgeCapability for $t {
            const NAME: &'static str = $name;
            const REGISTERED: &'static [&'static str] = &[];
        }
    };
}
impl_capability!(FactorizationGivenObservables, "factorization_given_observables");
impl_capability!(SpatialTopologyGivenFactorization, "spatial_topology_given_factorization");
impl_capability!(SpatialMetricUpToGlobalScale, "spatial_metric_up_to_global_scale");
impl_capability!(CausalOrderGivenExternalClock, "causal_order_given_external_clock");
impl_capability!(ConformalLorentzianStructure, "conformal_lorentzian_structure");
impl_capability!(VolumeMeasure, "volume_measure");
// **ClockCalibration / FullLorentzianMetric には BridgeCapability を実装しない** —
// BridgeLawCertificate<ClockCalibration> はトレイト境界を満たせず**型レベルで
// 構成不能** (別の証拠が立つまで封鎖 — PROMPT/11「ClockCalibration、ProperTime、
// FullLorentzianMetric は別の証拠がない限り構成不能にすべき」)。sealed のみ
/// 実装し、外部からの後付け実装も封じる。
impl sealed_cap::Sealed for ClockCalibration {}
impl sealed_cap::Sealed for FullLorentzianMetric {}

/// 型境界を越える変換の唯一の門 (能力別)。能力 C の登録簿が空である限り構成不能。
pub struct BridgeLawCertificate<C: BridgeCapability> {
    claim_id: &'static str,
    _capability: std::marker::PhantomData<C>,
}

impl<C: BridgeCapability> BridgeLawCertificate<C> {
    /// 能力 C で登録済みの claim id に対してのみ証明書を発行する (現在は常に None)
    pub fn register(claim_id: &'static str) -> Option<Self> {
        if C::REGISTERED.contains(&claim_id) {
            Some(BridgeLawCertificate {
                claim_id,
                _capability: std::marker::PhantomData,
            })
        } else {
            None
        }
    }
    pub fn claim_id(&self) -> &'static str {
        self.claim_id
    }
    pub fn capability(&self) -> &'static str {
        C::NAME
    }
}

/// 格子点 → 関係要素 (禁止変換 1 の唯一の門 — factorization 能力の証明書必須。
/// v29.5 [C5] が機械記録したとおり readout は factorization を選定できない —
/// この能力の登録には FactorizationBridge の確立が要る)
pub fn promote_site_to_node(
    site: RegulatorSiteId,
    cert: &BridgeLawCertificate<FactorizationGivenObservables>,
) -> RelationalNodeId {
    let _ = cert.claim_id;
    RelationalNodeId(site.0)
}

/// 外部計量応答 → 創発計量候補 (禁止変換 3 の唯一の門 — 空間計量能力の証明書必須)
pub fn promote_external_to_emergent(
    source: &ExternalMetricSource,
    cert: &BridgeLawCertificate<SpatialMetricUpToGlobalScale>,
) -> EmergentMetricCandidate {
    EmergentMetricCandidate {
        components: vec![source.amplitude],
        bridge_claim_id: cert.claim_id(),
    }
}

// 旧 promote_evolution_to_proper_time (発展パラメータ → 固有時間) は v30.0 で
// **削除** — 対応する能力 ClockCalibration が BridgeCapability 未実装のため、
// この門は型として書けない (禁止変換 2 は関数の不在で強制される)。

// ---------------------------------------------------------------- 自己検査

/// qrn_core の不変条件 (v272_core_contract と lib.rs::self_test から呼ぶ):
/// 登録簿が空・証明書が発行不能・状態表示が水増しされていないこと。
pub fn qrn_core_self_test() -> Result<(), String> {
    if !REGISTERED_BRIDGE_LAWS.is_empty() {
        return Err(format!(
            "REGISTERED_BRIDGE_LAWS が空でない: {:?} — bridge law の確立は claims.yml と spec の改訂を伴うこと",
            REGISTERED_BRIDGE_LAWS
        ));
    }
    // 能力別登録簿が全て空・全能力で発行不能 (v30.0)
    fn cap_locked<C: BridgeCapability>() -> bool {
        C::REGISTERED.is_empty() && BridgeLawCertificate::<C>::register("QRN-META-030").is_none()
    }
    if !(cap_locked::<FactorizationGivenObservables>()
        && cap_locked::<SpatialTopologyGivenFactorization>()
        && cap_locked::<SpatialMetricUpToGlobalScale>()
        && cap_locked::<CausalOrderGivenExternalClock>()
        && cap_locked::<ConformalLorentzianStructure>()
        && cap_locked::<VolumeMeasure>())
    {
        return Err("能力別登録簿に登録がある / 未登録 id に証明書が発行された".into());
    }
    // 状態表示の水増し検査: dynamics/bridge に Defined が混入していないこと
    for c in &QRN_CORE_V1 {
        if (c.layer == "dynamics" || c.layer == "bridge") && c.status == ContractStatus::Defined {
            return Err(format!(
                "{} (layer: {}) が Defined — 固有原理・bridge law なしに Defined を名乗れない",
                c.name, c.layer
            ));
        }
    }
    // 居住不能型のサイズ 0 (空 enum は値を持たない)
    if std::mem::size_of::<evidence::QrnEvidence>() != 0 {
        return Err("QrnEvidence が居住可能になっている".into());
    }
    Ok(())
}

// ---------------------------------------------------------------- 減衰長の型 (v29.1)

/// 指数減衰率 κ > 0 (単位: 1/ノード間隔)。構成は TryFrom<LinearFit> のみ —
/// ln w vs d の**傾き**が負の有限値であることを型が強制する (v28.2/28.3 の
/// HOLD-3 採点器が切片を減衰率に使い裁定を反転させた事故 [PROMPT/11] の恒久対策)。
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct DecayRate(f64);

impl DecayRate {
    pub fn value(self) -> f64 {
        self.0
    }
}

impl TryFrom<crate::LinearFit> for DecayRate {
    type Error = &'static str;
    fn try_from(fit: crate::LinearFit) -> Result<Self, Self::Error> {
        if !fit.slope.is_finite() || fit.slope >= 0.0 {
            return Err("指数減衰を示す負の有限傾きではない");
        }
        Ok(DecayRate(-fit.slope))
    }
}

/// 核減衰長 ξ = 1/κ (単位: ノード間隔)。**物理的相関長と呼ばない** — カーネル冪
/// (|C| vs |C|²) や代数的前因子が未分離のため、これは「当該核の実効減衰長」である
/// (v29.1 の命名修正: 旧称「ξ (相関長)」を KernelDecayLength に限定)。
#[derive(Debug, Clone, Copy, PartialEq, PartialOrd)]
pub struct KernelDecayLength(f64);

impl KernelDecayLength {
    pub fn value(self) -> f64 {
        self.0
    }
}

impl TryFrom<crate::LinearFit> for KernelDecayLength {
    type Error = &'static str;
    fn try_from(fit: crate::LinearFit) -> Result<Self, Self::Error> {
        let rate = DecayRate::try_from(fit)?;
        let xi = 1.0 / rate.value();
        if !xi.is_finite() || xi <= 0.0 {
            return Err("正の有限な核減衰長を構成できない");
        }
        Ok(KernelDecayLength(xi))
    }
}
