/-
QrnPromotion.lean — 昇格禁止と QRN-Core v1 状態表示の形式化 (v27.2, PROMPT/10 §4)

検証内容 (全て有限・決定可能 — decide / rfl のみ、native_decide 不使用):
  [1] 証拠クラスの昇格グラフ: 許される辺は「弱化」(強い証拠で弱い主張を支持する)
      のみで、bridge law 登録簿 (certified) は空。推移閉包 (reach 8) の下で
        calibration          ↛ qrnEvidence
        internalPrediction   ↛ naturalObservation
        sameAuthorReplication↛ independentReplication
        regulatorQuantity    ↛ universalQuantity
      であり、いかなる経路でも証拠の強さ (strength) は増えない。
  [2] QRN-Core v1 の状態表示 (sim/src/qrn_core.rs の QRN_CORE_V1 の鏡像):
      kinematics 層のみ defined であり、dynamics/bridge 層に defined は存在しない
      (「空の dynamics/bridge を Core 完成と呼ぶ」ことの形式的禁止)。

注意 (PROMPT/10 §4 の明示要件): ここで証明されるのは昇格グラフと状態表示の
**形式的性質**であって、分類の物理的正しさではない (それは ASM-LAYER-SEMANTICS)。
Rust 側の対応物: sim/src/qrn_core.rs (REGISTERED_BRIDGE_LAWS = [] / 居住不能型)。
機械照合: v272_core_contract [T5] と v151_audit [8] (定理数 12)。
-/

namespace QrnPromotion

/-- 証拠クラス: 上段 4 つ = 現在構成可能な証拠 / 下段 4 つ = 禁止昇格の到達先
    (Rust 側では居住不能型)。 -/
inductive EvidenceClass where
  | calibration            -- 較正・器械認証
  | internalPrediction     -- 計算実験内の凍結予言の的中 (PRED-013〜018)
  | sameAuthorReplication  -- 同一作者の別実装照合 (algorithmic diversity)
  | regulatorQuantity      -- regulator 依存量 (bare c₁/c₂ 型)
  | qrnEvidence            -- QRN (関係的存在論) の証拠
  | naturalObservation     -- 自然の観測量の的中 (future_observation × hit)
  | independentReplication -- 独立外部再現 (replications.yml の 6 条件)
  | universalQuantity      -- regulator 非依存が証明された普遍量 (bridge 側)
  deriving DecidableEq, BEq, Repr

open EvidenceClass

/-- 全クラスの列挙 (有限探索用) -/
def all : List EvidenceClass :=
  [calibration, internalPrediction, sameAuthorReplication, regulatorQuantity,
   qrnEvidence, naturalObservation, independentReplication, universalQuantity]

/-- 証拠の強さ: 下段 (QRN 証拠・自然の的中・独立再現・普遍量) = 1、上段 = 0。 -/
def strength : EvidenceClass → Nat
  | calibration => 0
  | internalPrediction => 0
  | sameAuthorReplication => 0
  | regulatorQuantity => 0
  | qrnEvidence => 1
  | naturalObservation => 1
  | independentReplication => 1
  | universalQuantity => 1

/-- 弱化 (許される方向): 強い証拠は弱い主張を支持できる。逆は不可。 -/
def weaken : EvidenceClass → EvidenceClass → Bool
  | qrnEvidence, calibration => true
  | naturalObservation, internalPrediction => true
  | independentReplication, sameAuthorReplication => true
  | universalQuantity, regulatorQuantity => true
  | _, _ => false

/-- bridge law 登録簿 — **空** (sim/src/qrn_core.rs の REGISTERED_BRIDGE_LAWS = [] の鏡像)。
    昇格辺の追加は、成功条件 7 項 (docs/qrn-core-v1-spec.md §5) を満たす証拠を
    伴う本ファイルの改訂としてのみ行う。 -/
def certified : List (EvidenceClass × EvidenceClass) := []

/-- 1 ステップの許容遷移 = 弱化 ∪ 登録済み昇格 -/
def step (a b : EvidenceClass) : Bool := weaken a b || certified.contains (a, b)

/-- n ステップ以内の到達可能性 (8 クラスなので n = 8 で閉包) -/
def reach : Nat → EvidenceClass → EvidenceClass → Bool
  | 0, a, b => a == b
  | n + 1, a, b => reach n a b || all.any (fun c => reach n a c && step c b)

/-- [定理 1] 登録簿は空である (昇格の門は閉じている) -/
theorem empty_registry : certified = [] := rfl

/-- [定理 2] 弱化は強さを増やさない -/
theorem weaken_sound :
    (all.all fun a => all.all fun b =>
      !(weaken a b) || decide (strength b ≤ strength a)) = true := by decide

/-- [定理 3] 許容遷移 1 ステップは強さを増やさない (登録簿が空である限り) -/
theorem step_no_upgrade :
    (all.all fun a => all.all fun b =>
      !(step a b) || decide (strength b ≤ strength a)) = true := by decide

/-- [定理 4] いかなる経路 (推移閉包) でも強さは増えない -/
theorem reach_no_upgrade :
    (all.all fun a => all.all fun b =>
      !(reach 8 a b) || decide (strength b ≤ strength a)) = true := by decide

/-- [定理 5] 較正は QRN の証拠に到達しない -/
theorem calibration_never_qrn : reach 8 calibration qrnEvidence = false := by decide

/-- [定理 6] 計算実験内の的中は自然の観測量の的中に到達しない -/
theorem internal_never_natural :
    reach 8 internalPrediction naturalObservation = false := by decide

/-- [定理 7] 同一作者の別実装照合は独立外部再現に到達しない -/
theorem same_author_never_independent :
    reach 8 sameAuthorReplication independentReplication = false := by decide

/-- [定理 8] regulator 量は普遍量に到達しない -/
theorem regulator_never_universal :
    reach 8 regulatorQuantity universalQuantity = false := by decide

-- ---------------------------------------------------------------- 状態表示

/-- QRN-Core v1 の構成要素 (qrn_core.rs QRN_CORE_V1 の鏡像) -/
inductive Component where
  | stateSpace | observableAlgebra | relationalDecomposition
  | evolutionLaw | constraintAlgebra
  | geometryBridge | causalBridge | clockBridge | matterBridge | gravityBridge
  deriving DecidableEq, BEq

inductive Layer where
  | core | dynamics | bridge
  deriving DecidableEq, BEq

inductive CoreStatus where
  | defined | modelFamilyOnly | conjectural | unsupported | undefined
  deriving DecidableEq, BEq

open Component

def comps : List Component :=
  [stateSpace, observableAlgebra, relationalDecomposition,
   evolutionLaw, constraintAlgebra,
   geometryBridge, causalBridge, clockBridge, matterBridge, gravityBridge]

def layerOf : Component → Layer
  | stateSpace => .core
  | observableAlgebra => .core
  | relationalDecomposition => .core
  | evolutionLaw => .dynamics
  | constraintAlgebra => .dynamics
  | _ => .bridge

def statusOf : Component → CoreStatus
  | stateSpace => .defined
  | observableAlgebra => .defined
  | relationalDecomposition => .defined
  | evolutionLaw => .modelFamilyOnly
  | constraintAlgebra => .modelFamilyOnly
  | geometryBridge => .conjectural
  | causalBridge => .conjectural
  | clockBridge => .undefined
  | matterBridge => .conjectural
  | gravityBridge => .unsupported

/-- [定理 9] kinematics (core 層) は全成分 defined -/
theorem kinematics_defined :
    (comps.all fun c => !(layerOf c == .core) || (statusOf c == .defined)) = true := by
  decide

/-- [定理 10] dynamics / bridge 層に defined は存在しない —
    「空の dynamics/bridge を Core 完成と呼ぶ」ことの形式的禁止 -/
theorem no_silent_completion :
    (comps.all fun c =>
      !(layerOf c == .dynamics || layerOf c == .bridge) || !(statusOf c == .defined)) = true := by
  decide

/-- [定理 11] 重力 bridge は unsupported (v27.0 fork: graviton pole なし) -/
theorem gravity_bridge_unsupported : statusOf gravityBridge = .unsupported := rfl

/-- [定理 12] 時計 bridge は未構成 (ProperTime への変換の門は閉) -/
theorem clock_bridge_undefined : statusOf clockBridge = .undefined := rfl

end QrnPromotion
