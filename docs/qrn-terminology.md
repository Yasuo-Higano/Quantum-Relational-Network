# QRN 用語法と型境界 (v27.1 制定 — 第二十八期)

**Version**: v27.1 (2026-07-28)
**位置づけ**: PROMPT/10 §2 の成果物。本文書は「何が QRN 本体で、何が既知 QFT の
測定器か」を**コード上も文書上も混同不能にする**ための用語法・層分類・型境界を
定義する。機械可読の一次ソースは [core.schema.yml](../core.schema.yml)
(監査: `v271_core_audit` — スイートの常時実行層)。本文書と schema が食い違ったら
schema と監査バイナリが正。

---

## 1. 層 (layer) の定義 — 6 + 1 分類

すべての概念・主張は次のいずれか**一つ**の層に属する。「core」という語を
無修飾で使うことを廃止する (v27.0 §D.3 の「QRN-Core v1」は下記 4 層の複合体
だったことが監査で判明 — uft-v27.0.md §D.3 の監査注記を参照)。

| 層 | 内容 | 現在の状態 |
|---|---|---|
| `core` (QRN-Kinematics) | 状態空間・演算子代数・関係的分解・許容同値 — QRN の存在論的骨格 | **defined** (定義のみ。実装はガウス toy 1 族 = C3) |
| `dynamics` (QRN-Dynamics) | 状態を選び発展させる固有原理 (初期条件則・発展則・拘束代数) | **model_family_only** (固有原理は未定義 — RingChain 等の模型族だけがある) |
| `adapter` (Matter-on-Background) | 外部計量と物質の結合規約 — BOND-A・中点変調・Belinfante λ=−1/8・接触項・counterterm | **certified** (v26.2–v27.0 で認証。ただし QRN の公理ではなく格子離散化の契約) |
| `instrument` (Metrology) | 較正済み測定器 — null 結合・殻積分・Matsubara Ward・導出モデル外挿・監査台帳群 | **certified** (認証記録は v27.3 の instruments.yml に凍結) |
| `bridge` (QRN-Bridge) | 状態・相関から計量・因果・時計・物質を読み出す規則 | **conjectural** (toy 実演 C3 のみ。重力 bridge は unsupported) |
| `phenomenology` | 自然のデータへの適用 — SM 走査・湯川・flavor・宇宙論・BMV | C4 fit / C2 窓付き定理 (自然の観測量の的中 0) |
| `meta` | プログラム運営層 — 期統合・台帳・監査・凍結宣言 | 常設 |

**重要な帰結 (v27.1 の主判定)**: v27.0 §D.3 が「QRN-Core v1 の公理系」と呼んだ
BOND-A・中点変調・Belinfante λ=−1/8 は、QRN の存在論的公理ではなく
**`adapter` 層の格子契約**である。これらが定義するのは

> 外部計量上の Dirac 物質を、格子正則化で正しく測定するための認証済み計量応答カーネル

であり、QRN の中心命題「量子相関網から時空・重力・物質・因果が読み出される」
(= `bridge` 層) はまだ仮説段階のままである。

## 2. 型境界 — 同一視してはならないもの

以下の各組は**別の型**であり、無証明の変換 (暗黙のキャスト) を禁止する。
Rust 型としての実装は v27.2 (`sim/src/qrn_core.rs`)、機械検査は
`v271_core_audit` (禁止 `impl From` の不在検査)。

### 2.1 ノード 3 型 — 存在論と正則化の分離

```rust
struct RelationalNodeId(u64); // QRN 存在論上の関係要素 (実在の仮説, layer: core)
struct RegulatorSiteId(u64);  // 数値離散化上の格子点 (regulator, layer: adapter)
struct ContinuumPoint;        // 有効理論上の点 (連続極限の座標, layer: adapter)
```

v27.0 の fork 執行で「格子 = regulator」(存在論ではない) が確定した。一方
旧文書 (v1.0 A0, v6.7) はテンソル因子ノードを物理的実在として扱う。両者を
同じ型に潰すと「格子由来の構造の観測」が「自然が量子ネットワークである証拠」
に無音で昇格する。**禁止変換: `RegulatorSiteId -> RelationalNodeId`**
(許すには明示的な `BridgeLawCertificate` が要る — v27.2 で型化)。

### 2.2 計量 2 型 — 外部と創発の分離

```rust
struct ExternalMetricSource;     // 外部から入力する source h_μν (layer: adapter — v26.2–v27.0 の全て)
struct EmergentMetricCandidate;  // 状態から読み出される計量の候補 (layer: bridge — 現状 unsupported)
```

v27.0-C の確定事項: 自由場 matter loop に有限留数の graviton pole はなく、
c₂ は a⁻² 走行の regulator 量。経路 B の全成果は **ExternalMetricSource への
応答測定**であり、EmergentMetricCandidate の証拠ではない。
**禁止変換: `ExternalMetricSource -> EmergentMetricCandidate`** (external metric
response ≠ emergent metric)。

### 2.3 時間 4 型 — 発展パラメータと創発時間の分離

```rust
struct EvolutionParameter(f64);      // シミュレーションの発展パラメータ t (layer: dynamics — 器械の入力)
struct ModularParameter(f64);        // モジュラー流のパラメータ s (layer: bridge — 熱時間仮説の変数)
struct OperationalClockReading(f64); // 部分系時計の読み (layer: bridge — 操作的時間)
struct ProperTime(f64);              // 固有時間の候補 (layer: bridge — 未構成)
```

公理 A1 は「時間は外部パラメータではなくモジュラー流などから創発する」と
述べるが、旧 `QrnModel::evolve(&self, s, t: f64)` の `t` は**外部発展パラメータ**
であり、A1 の創発時間ではない。現状は両者が同じ `f64` に潰れている
(uft-v1.0.md A1 の監査注記を参照)。**禁止変換:
`EvolutionParameter -> ProperTime`** (bridge law なしの同一視は A1 の先取り)。
> v30.0 追記: 証明書は能力別 `BridgeLawCertificate<C>` に分割され、ProperTime へ
> 至る能力 (ClockCalibration) は BridgeCapability 未実装で**型レベル構成不能** —
> この変換の門は関数ごと削除された (docs/uft-v30.0.md §1)。

### 2.4 証拠 4 型 — 昇格の禁止

```text
CalibrationEvidence<T>  -> EvidenceFor<QRN>        禁止 (較正は理論の証拠でない)
InternalPrediction<T>   -> NaturalObservation<T>   禁止 (計算実験内の的中は自然の的中でない)
SameAuthorReplication   -> IndependentReplication  禁止 (同一作者・同一 AI の別言語実装は
                                                    algorithmic diversity どまり)
RegulatorQuantity<T>    -> UniversalQuantity<T>    禁止 (bare c₁/c₂ 型の scheme 量)
```

paper/grav-vacuum-polarization-spec.md §5 の禁止暗黙変換 6 種 (bond modulation ≠
vierbein / T00 保存 ≠ 保存 Tμν / plus ≠ spin-2 / 質量依存 ≠ regulator 依存 /
1/χ ≠ propagator / q² 係数 ≠ Newton 定数) はそのまま恒久有効で、本表はその
存在論・証拠論への拡張である。

## 3. 旧名の限定 (改名台帳)

| 旧名 | 正名 | 限定の内容 |
|---|---|---|
| `QrnState` | `GaussianFermionState` | ガウスフェルミオン相関行列の toy 状態 (C3)。QRN の一般状態空間ではない |
| `QrnModel` | `GaussianToyModel` | ガウス toy 模型族のインターフェース。QRN-Dynamics の固有原理ではない |
| 「QRN core」(v6.7–v27.0) | `GaussianFermionState` + 読み出し群 | kinematics の toy 実装 + bridge 読み出しの複合体。層を混同した旧称 |
| 「QRN-Core v1 の公理系」(v27.0 §D.3) | Matter-on-Background Adapter contract | adapter 層の格子契約 (§1 参照) |

コード上の改名は **v27.2 で実施済み** (`sim/src/qrn_core.rs` の導入と同一コミット。
旧名の残存ゼロは `v272_core_contract` [T2] が常設検査 — コメント内の改名台帳への
歴史的言及のみ許す)。加えて V2 系も同時に改名した: `QrnStateV2` →
`ConstrainedToyStateV2` / `QrnDynamicsV2` → `ConstrainedToyDynamicsV2`
(拘束模型族の toy — QRN の一般動力学ではない)。sim/src 共有部の変更に伴う
全スイート再検証の儀式は v28.0 で一括実施する。

## 4. 主張台帳の多軸 (v27.1 で claims.yml 全 214 件に付与)

C0–C5 (確立度) と直交する 6 軸。語彙は core.schema.yml が一次ソース、
書式は `v61_ledger` [7]、意味論は `v271_core_audit` が機械検査する。

| 軸 | 語彙 | 意味 |
|---|---|---|
| `layer` | core / dynamics / adapter / instrument / bridge / phenomenology / meta | §1 の層 |
| `evidence_kind` | theorem / reproduction / calibration / mechanism_demo / internal_holdout / interpretation / **natural_observation** / **external_replication** | 証拠の種類 (太字 2 種は現在 **0 件** — 機械監査) |
| `independence` | same_implementation / algorithmically_diverse / same_author_clean_room / **independent_author** | 証拠の独立性 (independent_author は現在 **0 件**) |
| `universality` | not_applicable / regulator_specific / scheme_dependent / continuum_universal / unknown | regulator 依存性 |
| `data_relation` | not_applicable / fitted / calibration_data_reused / preregistered_holdout / future_observation | データとの関係 |
| `physical_scope` | toy / effective_model / laboratory / natural | 主張が語る領域 |

### 昇格禁止規則 (v271_core_audit が全件検査)

```text
R1  evidence_kind = natural_observation の件数 = 0 (自然の観測量の的中 0 の機械化)
R2  evidence_kind = external_replication の件数 = 0 (独立外部再現 0 の機械化)
R3  independence = independent_author の件数 = 0 (同上)
R4  layer ∈ {core, dynamics, bridge} かつ C0 以外 ⇒ physical_scope ∈ {toy, effective_model}
    (toy mechanism → theory of nature の禁止)
R5  data_relation = future_observation を持つ主張は存在しない (未来のデータは主張を支えない)
R6  universality = continuum_universal ⇒ layer ∈ {adapter, instrument} または C0
    (連続普遍性は認証済み測定器の専有 — bridge/core への漏出禁止)
R7  evidence_kind = internal_holdout ⇒ data_relation = preregistered_holdout
```

R1 の運用: predictions.yml の scored-hit のうち「自然の観測量の的中」と数えるのは
physical_scope = natural **かつ** data_relation = future_observation のものだけ
(PRED ごとの分類表は core.schema.yml)。γ_UT (PRED-007) や θ13 (PRED-011) の hit は
**測定値が公知の holdout** (preregistered_holdout) であり、この定義で自然の的中に
数えない — README の「自然の観測量の的中 0」と整合する。

## 4b. 二層分離 (v29.2 追記): FactorizationBridge と GeometryBridge

bridge 層を二層に分離する (PROMPT/11 の意味論監査):

```text
global algebra/state → [FactorizationBridge] → local observable subalgebras
                     → [GeometryBridge]      → adjacency / metric / causal order
```

- **FactorizationBridge** (未定義): 大域代数・状態から局所部分代数 (テンソル因子)
  を抽出する規則。現行の全 bridge 成果はこれを**入力**に取る — NodeState はモードを
  ノードへ群化して受け取り、隠しているのは座標・隣接・ラベルであって分解ではない。
- **GeometryBridge** (現行成果 = C3): 与えられた部分代数間の相関・応答から幾何を
  抽出する。第二十九期〜v29.1 の確定主張の正確な scope は「**与えられたノード因子
  分解の下で**、ラベル・座標・隣接を入力せず、隠れた 1 次元隣接・位相・因果順を
  回復した」。
- 減衰長の正名 = **KernelDecayLength** (v29.1) — カーネル冪・前因子が未分離のため
  「物理的相関長」と呼ばない。

## 5. 検証と限界

- 本文書は分類の**規約**であり、分類自体の正しさは仮定 ASM-LAYER-SEMANTICS
  (assumptions.yml) に登録した。誤分類は v271_core_audit の網羅性を破る。
- 分類は「主張の役割」を一意に潰す簡約であり、複数層にまたがる主張
  (例: 期統合) は meta に逃がした。境界例は core.schema.yml の note に記録する。
