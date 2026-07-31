# QRN v31.0 — 第三十一期 開始: 幾何読み出しの識別可能性相図 — 意味論・型・protocol の凍結

**Version**: v31.0 (2026-07-31)
**Sim**: `sim/src/readout_contract.rs` (識別可能性契約) + `sim/src/bin/v310_readout_semantics.rs`
→ `results/v310_readout_semantics.txt` (8 検査 PASS)
**位置づけ**: PROMPT/12 ロードマップの第一版。第三十一期の主題は**高 genus・3D への
単純拡張ではなく、「幾何読み出しの識別可能性相図」の完成** — 何が状態に符号化され、
何が許された観測から読め、何が原理的に読めないかを、能力・状態領域・観測契約ごとに
確定する。

**基準**: main = e8f9daf (v30.0-B)。第三十期の全履歴・凍結節・SECRET・失敗記録は
改変しない。総計 PASS 1227 / FAIL 0 を baseline manifest とする。bridge law 登録簿は
全能力で空・external replication = 0 を維持。PRED-019 は未登録のまま開始する。

---

## 1. 中心命題 — 「読み出し」は一語ではない

> **状態に幾何が符号化されていること、許された観測から局所的に読み出せること、
> 物理的生成子として同定できること、因子分解自体を発見できることは、
> すべて別の能力である。**

v29.5 の P6/693 衝突 (半充填 projector の静的核は sign(A) まで) と v29.6/HOLD-6 の
熱的 Gaussian 成功は別現象ではない — 前者は純粋基底状態が生成子のスペクトル強度を
捨てる極限、後者は有限温度状態がそれを保持する領域である。臨界基底状態だけが
境界増強・Friedel 共鳴で崩れた事実 ([T2b]) もこの境界を補強する。今後は次の 5 段階を
混同しない:

| 段階 | 問題 | v30.0 時点の残高 |
|---|---|---|
| **E0**: Encoding | 完全な状態の中に生成子・幾何の情報が存在するか | 熱的 Gaussian で強い正例 |
| **E1**: Global inversion | 大域状態を完全に知れば復元できるか | logit 逆変換でほぼ自明 (v31.1 で oracle ceiling として機械化) |
| **E2**: Operational readout | pair RDM・密度・局所応答など制限観測で読めるか | HOLD-5/6 が強い正例 — ただし与えられた因子分解・指定状態族の下で |
| **E3**: Factorization discovery | ノード分解そのものを状態・操作から選べるか | 未定義・原理的非一意性 (v29.5 [C5] — v31.4 で no-go を定理化) |
| **E4**: Physical bridge | 自然系の時空・計量との対応か | 外部再現 0・自然観測の的中 0 |

v30.0 までの成果は主として **E2 (与えられた因子分解と指定状態族の下)** である。
E3・E4 との混同を型で禁止するのが本版の仕事。

## 2. 契約の 4 軸 (`sim/src/readout_contract.rs`)

読み出し証明書は能力だけでなく状態領域・観測契約・因子分解状態を型に持つ:

```rust
ReadoutCertificate<Capability, StateDomain, ObservationContract, FactorizationStatus>
```

- **Capability**: qrn_core の 8 能力タグを再利用 (ClockCalibration / FullLorentzianMetric
  は BridgeCapability 未実装のため読み出し証明書も構成不能)。
- **StateDomain** (4, sealed): `GaussianGibbsFullRank` / `GaussianProjector` /
  `InteractingFermion` / `UnknownStateDomain`。**full-rank Gaussian の結果を pure /
  interacting へ暗黙拡張しない** (絶対禁止)。
- **ObservationContract** (6, sealed): `GlobalOneBodyCorrelation` / `PairReducedStates` /
  `StaticLocalObservables` / `LocalBiasDensityResponse` / `CoherentLocalResponse` /
  `RetardedResponse`。強さの半順序はここでは定義しない — v31.3 が実測する。
- **FactorizationStatus** (3, sealed): `GivenNodeFactorization` /
  `OperationallyInferredFactorization` / `UnknownFactorization`。

**裁定 5 値** (`IdentifiabilityVerdict`): `ExactUpToGauge` / `ExactUpToGlobalScale` /
`StableEstimate {condition_bound}` / `EquivalenceClassOnly` / `Abstain(理由 8 種)`。
**非識別セルの正しい棄却は失敗ではなく一級の読み出し結果** — HOLD-7 (v32.0) は
非識別セルで無理に回答したら FAIL、正しい `EquivalenceClassOnly`/`Abstain` を PASS と
採点する。棄却理由 8 種: RankDeficient / IllConditioned / GaussianityUnverified /
GibbsProvenanceMissing / UnknownFactorization / InsufficientObservation /
NonGaussianDomain / RankDeficientLocalGram。

**ReadoutCertificate は昇格の門ではない** — bridge law の門は qrn_core の
`BridgeLawCertificate` (全能力で登録簿空) のみ。

## 3. 生成子 2 型と門 (禁止変換 8)

熱的 Gaussian 生成器 C = (I + e^{β(h−μI)})⁻¹ に対し 0 ≺ C ≺ I なら

```text
K(C) = log[(I−C)C⁻¹] = β(h−μI)
```

だが、**logit が返すのはまず「親モジュラー生成子」である**:

- `ParentModularGenerator` — 完全な C と整合する Gaussian parent K(C)。
- `PhysicalGenerator` — 模型の実時間発展生成子 h。
- 変換の唯一の門 = `identify_physical_generator`。**GaussianityEvidence**
  (ByConstruction / WickResidualBound / Unknown) と **GibbsProvenance**
  (KnownBetaMu / BetaUnknownPositive / Missing) の証拠を要求し、証拠がなければ
  棄却を返す。β 未知は正の大域スケール同値類 (`UpToPositiveScaleAndShift`)、
  μ 未知は一様対角シフトのみで空間隣接に影響しない。

この逆変換は**新しい QRN 法則ではない** — 自由フェルミオンの標準的事実であり、
Gibbs 状態からの Hamiltonian learning は独立した既存分野である。よって v31.1 の実装は
`GaussianGibbsInverseOracle` と命名し、**oracle ceiling (識別可能性の上界)** として
登録する。QRN 固有予言は不可・PRED-019 は登録不可・この式だけで
BridgeLawCertificate は登録不可。新規性が生じるのは**静的・動的・局所・大域・純粋・
熱的・Gaussian・非 Gaussian を横断した識別可能性相図、厳密な no-go 証明書、
観測資源別の読み出し限界**である。

## 4. 相関 2 型 (禁止変換 9) と logit 2 経路

- `ExactFullRankCorrelation` — 唯一のコンストラクタが **clamp なしのスペクトル証明書**
  (0 < λ < 1 かつ margin δ = min(λ_min, 1−λ_max) ≥ 1e-13 [f64 の器械床]) を内部計算で
  要求。エルミートは 2n×2n 実対称埋め込みで審査。
- `RegularizedCorrelation` — clamp/正則化済み。estimate lane 専用。**exact 識別可能性
  証明書の根拠にできない** (現行 B2 の eigenvalue clamp は数値推定として正当だが
  exact の根拠ではない)。

一般に f(P C P) ≠ P f(C) P なので、logit の 2 経路は別型:

- `GlobalPhysicalParentBlock` — 全系 C に logit をかけてから block 抽出。
- `ReducedModularBlock` — 二ノード RDM に logit (現行 B2) — 環境で renormalize された
  reduced modular coupling。

両者は観測予算の両端として v31.1 で比較し、一致しないことを反例で常設検査する。

## 5. patch 2 型 (禁止変換 10)・支持証明書・操作代数

- `OraclePatch` (診断専用: 真の隣接半径) と `OperationalPatch` (観測から構築) は別型 —
  真の幾何半径で patch を選ぶ読み出しは循環する。
- `ObservableSupportCertificate {rank, threshold, nullspace_dim}` — 局所 Gram の rank
  欠損時に無条件擬似逆を禁止し、support 制限か棄却を強制 (v31.4)。
- `OperationalAlgebra {preparations, interventions, measurements, compatibility}` —
  E3 (因子分解の選定) に必要な操作的構造の定義 (v31.4 で no-go 側を定理化)。

## 6. schema・台帳への波及

- **core.schema.yml**: 新概念 15 種を登録。**RelationalDecomposition の意味論差を是正**
  — 旧 note「テンソル分解は入力でなく読み出し (v11.4)」は、現行 bridge が因子分解を
  入力に取る事実 (v29.5 [C5]) と矛盾していた。新 note は「状態からの選定は設計目標
  (RelationalDecompositionGoal, 未構成) — 現行成果は GivenNodeFactorization を入力に
  取る」と分離する (v11.4 の MI マッチングは正例デモで選定則ではない)。禁止変換
  8–10 (親→物理生成子 / 正則化→exact / oracle→operational patch) を追加。
  注: `bridge_candidates.yml` の旧文言は v29 の実行前凍結ファイルのため**改変しない**
  (歴史的記録)。
- **replications.yml**: claim/capability scoped に拡張 (Unit D schema — D1 数値再現 /
  D2 end-to-end / D3 負定理 / D4 応答法則)。**gauge の Unit A/B/C の成功は geometry
  能力の blocker を解除しない — 解除は matching D2 のみ**。external_replications = 0
  維持。
- **qrn_core.rs**: RelationalDecomposition の note のみ是正 (追加以外の変更なし)。
  v30.0 の封鎖 (ProperTime 門の不在・ClockCalibration 構成不能・登録簿空) は
  v310 [R5] が不変を検査。

## 7. 絶対禁止 (期間中 — PROMPT/12)

1. ParentModularGenerator を無条件に PhysicalGenerator へ変換しない。
2. Full-rank Gaussian の結果を pure / interacting state へ暗黙拡張しない。
3. Regularized/clamped matrix から ExactIdentifiabilityCertificate を発行しない。
4. 真の adjacency・coordinate・graph distance・hidden h を readout に渡さない。
5. operational patch の選択に真の幾何半径を使わない。
6. P6/693 の static-projector collision を tuning で消そうとしない。
7. holdout 開封後に kernel・bar・generator・scorer を変更しない。
8. 外部再現を自己実装・同一 AI 実装で水増ししない。
9. PRED-019 を内部 toy theorem から登録しない。
10. open-boundary regulator mismatch の bar を緩めて成立扱いしない。

## 8. 儀式の扱い (suite 台帳)

共有部 (`lib.rs` への mod 宣言 1 行 + `qrn_core.rs` note 1 箇所 + 新規
`readout_contract.rs`) の変更は suite 台帳の全バイナリを無効化する。第三十一期は
v30.0 の前例に従い、**共有変更を本版 1 回に集約**し、各版では新規バイナリ + 監査層を
手動実行、**期末 (v32.0) の完全儀式で全バイナリの無波及を台帳比較で検証**する。
既存コードへの変更は追加のみ (qrn_core は note 文字列 1 箇所) — コンパイル互換は
`cargo build --release` 全数で確認済み。

## 9. ロードマップ (PROMPT/12)

| 版 | 主要課題 |
|---|---|
| v31.0 | 意味論・型・protocol の凍結 (本版) |
| v31.1 | GaussianGibbsInverseOracle — exact/estimate lane・条件数定理・n≤7 全数・P6/693 有限温度分離 |
| v31.2 | LocalBiasCommutatorLaw — 密度曲率 Frobenius 恒等式 (不変ノルム核の第一候補) |
| v31.3 | 観測予算 hierarchy — 7 lane の精度・coverage・abstention 曲線 |
| v31.4 | invariant operator response atlas / factorization no-go |
| v31.5 | 非 Gaussian transfer (spinless t-V + Z2) |
| v31.6 | 計量 VR persistence・高 genus・3D (β₃ には 4-simplex) |
| v31.7 | 開放境界 regulator 不一致の機構分類 (bar は変更しない) |
| Track X | 外部再現 Unit D (内部研究と並列) |
| v32.0 | HOLD-7 — 正解の的中と**非識別セルの正しい棄却**の両方を採点 |

## 10. 開発記録 (v31.0)

- 契約の封鎖は qrn_core と同型の 3 段構え (別型 / 唯一の門 / sealed) + 監査 4 重
  (v310 [R1] source 走査・[R2] schema・[R5] qrn_core 不変・自己検査)。
- [R7] は熱的 Gaussian round-trip の最小実演 (P4 鎖, β = 1.7, μ = 0.3): 資格審査 →
  K(C) → 門 → max|ĥ−h| < 1e-10 復元、projector は RankDeficient で正しく棄却。
  v31.1 が n ≤ 7 全数 × 複数 β へ拡張する。
- 「非識別の正しい棄却が一級市民」という採点原則を型 (Abstain が verdict の
  正規 variant) で先に凍結した — HOLD-7 の採点器を後から歪めないため。
