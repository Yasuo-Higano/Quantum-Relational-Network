# QRN v27.2 — 型付き Core Contract: 別型・居住不能型・bridge law の門

**Version**: v27.2
**Date**: 2026-07-29
**Sim**: `sim/src/qrn_core.rs` (共有モジュール) + `sim/src/bin/v272_core_contract.rs`
→ `results/v272_core_contract.txt` (6 検査 PASS)
**Lean**: `proofs/QrnPromotion.lean` (12 定理 — decide/rfl のみ、native_decide 不使用)
**位置づけ**: PROMPT/10 §4。v27.1 の層分類 (意味論) を**型 (コンパイラが拒否する境界)**
に降ろす。物理 run なし。

---

## 0. 設計の核 — 3 段の防御

1. **別型**: regulator/存在論 (`RegulatorSiteId` ≠ `RelationalNodeId`)・外部/創発計量
   (`ExternalMetricSource` ≠ `EmergentMetricCandidate`)・発展/固有時間
   (`EvolutionParameter` ≠ `ProperTime`) — 暗黙変換 (`impl From`) は存在せず、
   不在は `v271_core_audit` [S12] が常設走査。
2. **唯一の門**: `BridgeLawCertificate::register(claim_id)` — 登録簿
   `REGISTERED_BRIDGE_LAWS` は**空**なので常に `None`。bridge law の確立
   (成功条件 7 項, spec §5) を伴う改訂だけが門を開ける。
3. **居住不能型**: 禁止昇格の到達先 `QrnEvidence`・`NaturalObservation`・
   `IndependentReplication`・`UniversalQuantityForBridge` は**空 enum** — 値を
   構成するコード自体が書けない。「自然の観測量の的中 0・独立外部再現 0」が
   カウンタ (v27.1 R1–R3) に加えて**型レベルの事実**になった。

## 1. 改名 — toy の名から Qrn を外す

| 旧 (v6.7/v15.3) | 新 (v27.2) | 意味 |
|---|---|---|
| `QrnState` | `GaussianFermionState` | ガウスフェルミオン相関行列の toy 状態 (C3) |
| `QrnModel` | `GaussianToyModel` | ガウス toy 模型族 (RingChain/TfdPair/GrowingChain/PacketRing) |
| `QrnStateV2` | `ConstrainedToyStateV2` | 拘束を解いた基底上の多体 toy 状態 |
| `QrnDynamicsV2` | `ConstrainedToyDynamicsV2` | 拘束 toy 模型族の発展 |

`evolve(s, t)` / `step(s, dt)` の `t` は `EvolutionParameter` 型に変更 —
**f64 の裸渡しは型エラー**になり、「シミュレーションの発展パラメータ」と
「物理的に読み出された時間」の混同 (uft-v1.0.md A1 の監査注記) がコンパイル時に
遮断される。更新バイナリ: v67_core / v74_core2 / v94_qneccore / v114_tensorcore /
v153_corev2 / v154_continuum (旧名の残存ゼロは [T2] が常設検査)。

## 2. QRN_CORE_V1 — 状態表示の Rust const 化

`qrn_core::QRN_CORE_V1` (10 成分: StateSpace/ObservableAlgebra/
RelationalDecomposition/EvolutionLaw/ConstraintAlgebra/GeometryBridge/
CausalBridge/ClockBridge/MatterBridge/GravityBridge) が core.schema.yml の
状態表示と一致することを [T1] が照合。`qrn_core_self_test()` は lib.rs の
`self_test()` に組み込み (成功時無音) — **全バイナリが起動時に「登録簿が空・
dynamics/bridge に Defined の混入なし」を再検査する**。

## 3. Lean 形式化 (`QrnPromotion.lean`, 12 定理)

証拠クラス 8 頂点の昇格グラフ (許容辺 = 弱化のみ・登録簿 = 空) に対し:

- [1] 登録簿は空 / [2–4] 弱化・1 ステップ・**推移閉包 (reach 8) のいずれでも
  証拠の強さは単調非増加**
- [5–8] 4 禁止対の到達不能性: calibration ↛ qrnEvidence・internalPrediction ↛
  naturalObservation・sameAuthorReplication ↛ independentReplication・
  regulatorQuantity ↛ universalQuantity
- [9–12] 状態表示の鏡像: kinematics 全成分 defined・**dynamics/bridge 層に
  defined 不在 (空の Core を「完成」と呼ぶことの形式的禁止)**・gravity =
  unsupported・clock = undefined

**「Lean 証明済み」が意味するのは昇格グラフと状態表示の形式的性質であって、
分類の物理的正しさではない** (それは ASM-LAYER-SEMANTICS の規約)。
Lean 定理は 50 → 62 本 (8 ファイル)。

## 4. 改名バイナリの出力検証

| バイナリ | 判定 |
|---|---|
| v114_tensorcore / v154_continuum | **byte 同一** |
| v67_core / v74_core2 | 差分 = println 内の型名のみ (数値完全同一) |
| v94_qneccore | 差分 = タイミング行のみ (数値完全同一) |
| v153_corev2 | 数値 9 桁表示は全て同一。**発見**: 1e-15 級診断 |Δ| と JSON 下位桁が **run 間でドリフト** (下記) |

**開発記録 (発見): v153_corev2 の run 間非決定性** — 純粋性双対の診断
|Δ| = |S(A)−S(Aᶜ)| が同一バイナリの連続実行で 1.0e-15 → 1.9e-15 → 1.4e-15 と
揺れる。切り分け: (i) 旧ソース (v25.2 と同一) を今日ビルドしても committed 値と
不一致 (8.9e-16 vs 3.8e-15) ⇒ **v27.2 のパッチ起因ではない**。(ii) v67/v154 は
run 間 byte 同一 ⇒ lib 全般でなく v153 固有。(iii) 単スレッド・固定シード・
HashMap/時刻非依存で**機構は未特定**。判定閾値 1e-9 に対し 6 桁の余裕があり
PASS は不変だが、「並列化によらず結果が変わらない」という決定性規約の穴として
**v28.0 の残高に登録**する。

## 5. 検証 (2026-07-29)

- `v272_core_contract`: [T0] 登録簿空・水増しなし・居住不能型 / [T1] const ↔
  schema / [T2] 旧名残存 0 / [T3] evolve/step の型分離 / [T4] 証明書発行拒否 /
  [T5] Lean 12 定理 — **6 検査 PASS**。
- `lean QrnPromotion.lean` — 検証成功 (exit 0)。
- 監査層 (v61 [7 検査]・v151 [10 検査]・v271 [13 検査]・v217 ほか) 全 PASS。
- **スイート注記**: sim/src 共有部 (lib.rs + qrn_core.rs) の変更により全 164 本が
  「要再実行」となる。約 105 CPU 時間の完全再走の儀式 (数期に一度) は **v28.0 で
  一括実施**し、本版では変更 6 バイナリ + 新規 2 本 + 監査層の個別再検証で代える
  (開発記録: 素の cargo fmt 事故 [v27.1] の再発防止として、整形は対象ファイル
  指定に限定した)。

## 6. 残り (第二十八期)

- v27.3: instruments.yml (adapter/metrology の凍結)。
- v27.4: reproducer/ + replications.yml。
- v28.0: 完全再走の儀式 + Core v1 完了条件の判定 + 期統合
  (残高: v153 の run 間ドリフトの機構特定)。
