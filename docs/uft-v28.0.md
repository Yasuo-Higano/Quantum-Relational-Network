# QRN v28.0 — 第二十八期 完結: QRN-Core v1 の意味論的閉包・型付き境界・外部再現化

**Version**: v28.0 (期統合)
**Date**: 2026-07-29
**位置づけ**: 第二十八期 (PROMPT/10) の期統合。事前指令 = 「QRN-Core v1 の意味論的
閉包・型付き境界・外部再現化」— v27.0 §D.3 の成果物を無条件に QRN-Core v1 と
呼ばず、理論本体・測定器・橋渡し仮説を型レベルで分離した上で、Core v1 の完了条件
(spec §7) を判定する。**物理 run なしの期** — 成果は全て意味論・型・台帳・外部固定
であり、新しい物理主張は一つもない (それが目的である)。

---

## A. 総括 (v27.1 → v27.4)

**期のテーゼ: 「名を正すと、負けられる問いが立つ」**

- **v27.1 (核意味論監査)**: 全概念を 7 層に分類 (core.schema.yml, 67 概念・
  未分類 0)。**中心判定: §D.3 の「公理系」の正名は Matter-on-Background
  Adapter v1 + Metrology Suite v1** — QRN の公理ではない。全主張に 6 軸付与
  (**layer: core の主張は 0 件**が正直な姿)。昇格禁止 R1–R7 を機械化 —
  「自然の観測量の的中 0・独立外部再現 0」が文言からカウンタになった。
  **QRN-C0-001 を条件付き定理に限定・QRN-GRAV-001 (v0.7) を modular-response
  benchmark に降格** (Jacobson 1995/2015・Faulkner 2013 の前提と結論を区別)。
- **v27.2 (型付き Core Contract)**: 別型 (ノード 3・計量 2・時間 4)・bridge law
  の門 (登録簿 = **空**)・**居住不能型** (QrnEvidence 等 — 値を構成するコードが
  書けない)。旧 QrnState/QrnModel → GaussianFermionState/GaussianToyModel
  (toy の名から Qrn を外した)。Lean 12 定理 (推移閉包で 4 禁止対が到達不能・
  dynamics/bridge に defined 不在 — 水増しの形式的禁止)。
- **v27.3 (器械台帳の凍結)**: instruments.yml — 22 器械 × 16 フィールド
  (較正記録・陰性対照・故障モード = 教訓集・禁止解釈・**認証 SHA-256** —
  器械の無断変更が監査 FAIL になる体制)。常設回帰 5 本。
- **v27.4 (外部再現単位)**: reproducer/ — 3 単位 (アノマリー分類・λ 15 桁・
  真空偏極 universality) を共有コードなしで再実装可能に。replications.yml —
  独立外部再現の 6 条件を機械可読化 (**同一 AI は独立でない**と明記)。

## B. QRN-Core v1 完了条件の判定 (spec §7 — 全 10 項)

| 条件 | 判定 | 証拠 |
|---|---|---|
| ontology / regulator / background / instrument / bridge の分離 | 成立 | core.schema.yml (未分類 0) + v271 [S1] |
| 状態空間・観測量・発展則の状態の明示 | 成立 | QRN_CORE_V1 (10 成分) + spec §2–§3 + v272 [T1] |
| 外部時間と創発時間が別型 | 成立 | EvolutionParameter ≠ ProperTime (v272 [T3]) |
| 外部計量と創発計量が別型 | 成立 | ExternalMetricSource ≠ EmergentMetricCandidate (v271 [S12]) |
| 旧 Gaussian core の適用範囲の限定 | 成立 | 改名 + C3 toy 明記 (v272 [T2] 残存 0) |
| v0.7 / A3 の過大解釈の修正 | 成立 | QRN-C0-001 限定・GRAV-001 降格・監査注記 3 本 (v271 [S11]) |
| 自然観測 0・独立外部再現 0 の保持 | 成立 | v271 [S4/S5/S10] + v274 [P4] (常設) |
| 全監査 PASS | 成立 | 監査層 11 本 全 PASS (§E) |
| make suite PASS | 成立 | 完全再走の儀式 + 二相化再検証 (§E) |
| 負の結果と未定義部分が消されていない | 成立 | v26.9-D 凍結・undefined/unsupported 明示・Lean no_silent_completion |

**判定: QRN-Core v1 完了。** ただしこの「完了」が意味するのは「**QRN 固有の理論核
を構築できるかを、初めて失敗可能な形で問える状態**」であって統一理論ではない。
現時点で完成しているのは強い計算物理・理論監査プラットフォームである
(この一文を弱める改訂は spec §7 が禁止する)。

## C. 正直な残高 (不変の分類)

- 的中 = 計算実験内の機構予言 (PRED-013–018)。**自然の観測量の的中 0・
  独立外部再現 0** — v27.1 以降は宣言でなく機械監査 (R1–R3, v274 [P4])。
- **layer: core の主張 0 件・bridge law 0 本・empirical_prediction: none** —
  「QRN 本体」を支える結果はまだ一つもない。第二十八期はこの 0 を隠しうる構造
  (名称・暗黙変換・台帳の一次元性) をすべて撤去した。それが本期の成果である。

## D. 未解決問題の残高 (第二十九期への持ち越し)

1. **relational geometry bridge の比較** (次期の本丸): 外部計量を入力せず、
   状態と演算子代数だけから計量・因果・時計を読み出せるか。候補 B1–B4・
   訓練 2 系・holdout 4 系・失敗条件 (無調整で 2 微視的模型 × 2 regulator 一致
   しなければ棄却) を **bridge_candidates.yml に実行前凍結済み**。
   PRED-019 は解析的導出まで登録しない。
2. v153_corev2 の run 間ドリフトの機構特定 (1e-15 級診断 — 単スレッド・固定
   シードで機構未特定。決定性規約の穴として v27.2 で発見・登録)。
3. 相互作用系の q² 項生成 (PRED-013 の続き) — composite 路線は CG0–CG8
   (qrn-core-v1-spec.md §6) の封印のまま。
4. flavor 前方予言 τ = 1/12 + i/2 (待ち)・massive ρ₀ の 4D 分解・Wilson 温度チャネル。
5. 外部 (最優先): anomaly-search / modular-BW の投稿・**reproducer/ を添えた
   第三者再実装の公募**。

## E. スイート — 完全再走の儀式 (2026-07-29)

`make suite-full OUT=results/v280_full_suite.txt JOBS=12` — 直近の儀式 (v25.2,
136 本 905 PASS) 以来の完全再計算。sim/src 共有部の変更 (qrn_core.rs 導入・改名)
の波及と固定シード決定性・末桁ドリフトを一括検査。rustc 1.94.0。

- **全 166 本 実行: PASS 1161 / FAIL 2** — FAIL は物理ではなく**監査層の読み書き
  競合** (常時実行の v273 が、並列再走中の v268p の書き込み途中の結果を読み
  凍結判定文を見失った偽 FAIL)。**儀式が suite 設計の欠陥を発見した**。
- **修正**: 監査層 (ALWAYS_RUN) を全バイナリの後段に回す二相実行に変更
  (tools/suite.sh)。再検証: **実行 11 本 (監査層) + 引用 155 本 = 総計
  PASS 1163 / FAIL 0**。
- **末桁ドリフト検査: 既存 163 本の PASS/FAIL 数は儀式前後で完全一致**
  (差分は新規 3 本のみ)。v25.2 凍結値 (v252_manifest)・常設回帰 5 本の凍結
  判定文も無傷。既知の v153 run 間揺れ (1e-15 級診断) 以外の drift なし。
- 計算量: CPU 総和 539,181 秒 ≈ 150 CPU 時間 (最長 v192_boost3d 45,527 秒 —
  12 並列の overcommit で単独記録 30,508 秒より延伸)。

## 1. 残り

- なし (第二十八期 完結)。第二十九期 = bridge 比較 (bridge_candidates.yml の
  凍結プロトコル) から — 期の開始はユーザーの PROMPT を待つ。

## 2. 成果物

- v27.1: docs/qrn-core-v1-spec.md / docs/qrn-terminology.md / core.schema.yml /
  `v271_core_audit` (claims QRN-META-029/030)
- v27.2: sim/src/qrn_core.rs / proofs/QrnPromotion.lean (12 定理) /
  `v272_core_contract` (QRN-META-031)
- v27.3: instruments.yml / `v273_instruments` (QRN-META-032)
- v27.4: reproducer/ / replications.yml / `v274_reproducer` (QRN-META-033)
- v28.0: 本文書 / bridge_candidates.yml (実行前凍結) /
  results/v280_full_suite.txt / suite 二相化 (QRN-META-034)
