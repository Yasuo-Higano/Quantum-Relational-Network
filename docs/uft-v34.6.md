# QRN v34.6 — real-data lane: synthetic と実測の構造的分離

**Version**: v34.6 (2026-08-03)
**Sim**: `sim/src/bin/v346_real_data_lane.rs` → `results/v346_real_data_lane.txt` (7 検査 PASS)
**受け皿**: `reproducer/real_data/` (RECORDED_LANE.md / RECORD.schema.json /
REAL_DATA_LEDGER.yml / FIXTURES) + `reproducer/D2R_PACKET.md` (配布用)
**位置づけ**: PROMPT/15 §7。「synthetic shot noise を実測ノイズと呼ばない」を
型と受け皿で守り、v35.0 の科学的完成条件の判定土台を作る。

---

## 1. 二 lane の分離 (型)

`DataProvenance` (`sim/src/finite_data.rs`): `SyntheticCoverage` (真の interface
既知・coverage を直接採点・adversarial 生成可能 — HOLD 系と v34.5 robust atlas は
この lane) と `RecordedExperimental` (実装置の記録・latent 未知・model fit /
drift / 事前登録予測を採点) — **二つの間に変換は存在しない**。

## 2. recorded experimental lane の受け皿

- **RECORD.schema.json**: provenance const・チャネル shot 配列 {0,1}・
  **事前登録コミットメント** (未使用チャネル予測の sha256 — データ開示前の
  commit が検証可能であること)・**vendor topology commitment** (freeze 後開示)・
  失敗込み公開 = true。
- **drift gate**: split-half Clopper–Pearson 区間 (Bonferroni α/2m) が disjoint
  なら **OutOfDomain** — iid 契約の破れの正検出 (禁止変換 25/29 の運用形:
  破れたデータに iid 証明書を発行しない)。fixture: stationary (6/20 vs 6/20)
  通過・drifting (2/20 vs 14/20) 拒否。
- **事前登録機構**: sha256(開示予測) = commitment の機械照合 (HOLD の SECRET
  機構と同型)。
- **台帳**: recorded_runs = externally_operated_runs =
  preregistered_prediction_hits = **0** (数は実記録のみ・fixture は数えない)。

## 3. D2-R 配布パケット

PROMPT/15 §8:「campaign layer の追加実装はもう不要。次の仕事は schema の改良では
なく、**仕様パケットの配布と実施者の獲得**」— `reproducer/D2R_PACKET.md` に
実施者向けの単一案内を用意 (6 条件・読む順・提出方法・凍結された約束 [一件で
足りる・失敗も公開]・OCS-1.0 を clean-room の唯一の入力として参照)。
campaign layer (d2r-v1) は凍結のまま不変。**実配布と実施者の獲得は人間の作業**
であり、リポジトリは受け皿の完備までを担う。

## 4. v35.0 の科学的完成条件 (凍結)

次の**いずれか**が無い限り、期末表現は「instrumental closure」に留める
(PROMPT/15 §7):

1. externally operated D2-R report
2. 実データ上の事前登録された未使用応答チャネル予測の的中
3. 広い観測契約族を排除する新しい厳密 no-go

## 5. 残高

bridge law 空・PRED-019 未登録・自然の的中 0・external_replications = 0・
recorded_runs = 0 — 全て正直な 0 のまま。

## 6. 次 (v35.0-A/B)

HOLD-10 の二層: **HOLD-10S** (30 セル semantic adversarial — noise/provenance・
addressability/synthesis・factorization/context・graded/J・structured/resource の
5 群 × 6) と **HOLD-10C** (凍結生成分布から 600+ セル・300+ 回答セルの coverage
campaign — wrong_promotion_probability_upper_95 ≤ 0.01 ほか 8 指標)。
selective risk = 0.000 だけを完成条件にしない (v34.3 [F3b]/[F6-28] が根拠)。
