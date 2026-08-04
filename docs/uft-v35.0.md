# QRN v35.0 — HOLD-10: 二層 holdout (semantic 30 + coverage 640) と第三十四期統合

**Version**: v35.0-A (2026-08-03 — 凍結半) / v35.0-B (開封 — 本文書に追記)
**Sim**: `sim/src/bin/v350a_hold10_freeze.rs` → `results/v350a_hold10_freeze.txt`
(4 検査 PASS + train 満票 + 設計 12 シード満票)
**位置づけ**: PROMPT/15 §6。第三十四期の全器械 (有限データ意味論 v34.3・enumerator
v34.4・robust atlas v34.5・lane 分離 v34.6) を凍結し、**HOLD-9 型の「満票」を超えて
population-risk の上限を報告する**二層 holdout。

---

## 1. なぜ二層か (PROMPT/15 §6)

HOLD-9 の 9 answerable セルと FollowUp の 77 Answer は観測誤り 0 — しかし誤り 0/n
の片側 95% 上限は 0/9 → **28.3%**・0/77 → **3.82%** (v343 [F6-28] で機械化)。
母集団誤り率の上限を 1% 未満にするには誤り 0 の回答セルが **299 件以上**必要。
よって:

- **HOLD-10S** (30 セル・5 群 × 6) — 意味論 adversarial holdout。セル一覧は
  QRN-BRIDGE-040 と v350a 出力。回答 9 / 非識別 21 の設計。
- **HOLD-10C** (640 セル・凍結生成分布) — coverage campaign。回答 ≥ 300 を設計し
  **wrong_promotion の片側 95% 上限 ≤ 0.01** を含む 12 指標で採点。

`selective risk = 0.000` だけを完成条件にしない — これは第三十四期の定理
(v34.3 [F3b]: coverage は selective risk を含意しない) の採点への反映。

## 2. 凍結 (v35.0-A)

```text
sha256(SECRET) = c3a62f1b55d708e50915abd1c634f65d619c05ea52c2ab5283b8c1804d7851c4
holdout シード = SECRET から導出 (開封時に公開・機械照合)
lib pin (8 モジュール): operational_net 7898d244… / laboratory_interface a11f188f… /
  resource_profile f2fe4b96… / contextual_factorization e540eea6… /
  graded_recovery c1ebb02c… / structured_backend f57d9bdf… (以上 6 本は HOLD-9 pin
  と同値 — 第三十四期は v33 器械に触れていない) + finite_data 66352869… +
  factorization_enumerator e7a3e312…
凍結バー: α_cell = 0.0002 (促進/区間)・α_gate = 0.001 (drift/相関/ブロックの 3 契約
  ゲート)・τ = 0.3・N = 150 (misspec セルは 300–450)・answers ≥ 300・
  wrong_upper_95 ≤ 0.01・coverage_lower_95 ≥ 0.98・境界棄却 recall 1・
  misspec recall 1・回答率 ≥ 0.95・marginal→joint 0・窓再利用 0・lane drift 0
train (シード 35001, 可視): S 30 セル満票 + C 640 セル 12 指標満票
  (回答 469・W = 0 → upper 0.0064・coverage 554/554・境界棄却 10/10・
  misspec 32/32・回答率 0.9736)
設計走行: 12 シード (35002–35013) 全て満票
```

**相関 lane (PROMPT/15 §6 必須)**: iid 契約の破れを split-half だけでなく
**遷移数ゲート** (t ~ Binomial(n−1, 2p(1−p)) の CP 区間 vs p の CP 区間からの伝播
区間の disjoint 判定) で検出する — 同一周辺分布の持続的 Markov 鎖 (v343 [F6-29]
の反例) を契約側で拒否する器械。

## 3. 設計走行が発見した器械故障 (凍結前に修復 — 開発記録)

設計走行 (12 シード) が v34.4 enumerator の故障モードを 2 つ検出した:

1. **commutant の数値零空間の不安定** — 稠密 (無理数成分) の共役族で
   closure_center_basis が近零特異値の junk 方向を拾い、traceless 次元が 3 → 4 に
   膨れて orbit dedup が破れる (overlap = 3/4 = 0.750 ちょうど)。
   修復 = `commutant_basis` が返す基底の可換性を**機械再検証**して塵を落とす
   (証明書は再実行で検証する、の規律)。
2. **合成 commutant の冗長経路** — 被覆する明示メンバーが存在するのに singleton
   部分集合の補因子 (合成 commutant) が別候補として生き残る。
   修復 = **候補は族の極大可換 simple 部分集合のみ** (明示のメンバーがあるなら
   それを使う)。

修復後: probe 300 走で故障 0・既存 v344/v345 の出力 byte 一致 (裁定・数値無波及)。
「設計走行はシード頑健性の確認」(v34.0-A) が「設計走行は器械の故障モード検出」
としても機能した初の期。

## 4. 開封 (v35.0-B) — 本採点の確定表

SECRET を開示 (`HOLD10-51335c5…` — [H0] が sha256 = コミットメント c3a62f1b… と
FROZEN-HOLD10 区間の v350a との逐語一致・**lib pin 8 モジュール不変**を機械照合) し、
holdout シード 15384046040323535768 (SECRET 導出) で S 30 + C 640 セルを初生成・
本採点した (調整なし):

```text
HOLD-10S: 30/30 — selective risk 0・impossibility recall 21/21・
          answerable recall 9/9・強制回答 0
HOLD-10C: 回答 476 (設計 ≥ 300) — W (誤昇格) = 0
  wrong_promotion の片側 95% 上限 = 0.00627 ≤ 0.01   ← HOLD-9 に無かった報告
  coverage 554/554 (下限 0.99461 ≥ 0.98)
  boundary_abstention_recall 10/10・misspecification_recall 32/32
  (drift 16 + 相関 Markov 16 — 遷移数ゲート)
  answerable_recall 0.9904 (412/416)・強制符号回答 0・insufficient 22/22
  marginal_to_joint 0・窓再利用 0・structured/dense 裁定 drift 0
```

観測誤り 0 は「population リスク 0」ではない — 上限 0.63% が本 holdout が言える
上界であり、それを**指標として初めて要求・報告した**のが HOLD-10 の新しさ
(HOLD-9 の 9 回答の上限は 28.3% だった)。

## 5. 第三十四期の統合 — 期末判定は「instrumental closure」

**期テーゼ「可アクセス性証明書は実験者の宣言ではなく、登録済み観測契約の下で
有限データから得られる同時信頼集合である。集合が裁定境界を跨ぐ場合、Straddled /
EquivalenceClassOnly / Abstain は失敗ではなく有限データから導かれる唯一の正答で
ある」— 定理 (第四の no-go + Robust Promotion, Lean 10 定理)・型 (禁止変換 22–29・
DataProvenance)・器械 (robust atlas・enumerator・real-data lane)・holdout
(HOLD-10 二層満票 + population 上限) の全てで閉じた。**

| 版 | 成果 (確定) |
|---|---|
| v34.1 | FollowUp 受理 (REP-001 / Partially Replicated / external 0 維持)・Yukawa erratum (QRN-YUK-034 refuted_as_stated・族内定理・位相 ↛ 計量) |
| v34.2 | OCS-1.0 paper-closed spec (sha256 凍結・出力値なし)・closure manifest 20・probe 型分離の機械実証 (quench null / 1−Δ² / Busch) |
| v34.3 | 第四の no-go (Le Cam, Lean) + Robust Promotion (誤昇格 ≤ α)・selective risk ≠ coverage・禁止変換 22–29・0/n 上限の機械化 |
| v34.4 | sector-aware factorization enumerator (Wedderburn 証明書・出力 6 型) |
| v34.5 | robust atlas (σ_min 下界・worst-case xtalk・glue 区間・spectral-gap J・interval cost・lane 一致) |
| v34.6 | real-data lane (synthetic ↛ experimental)・D2-R 配布パケット・完成条件の凍結 |
| v35.0 | HOLD-10 凍結 → 開封 **S 30/30 + C 16 指標満票・population 上限 0.63%** (調整なし) |

**期末判定 (v34.6 で凍結した規約による)**: 外部 D2-R 報告 0・実データ上の事前登録
予測 0・完成条件に数える「新しい」厳密 no-go なし (第四の no-go は本期の計画内
成果) — よって期末表現は **instrumental closure** であり「QRN の物理的前進」とは
呼ばない。**正直な残高 (不変)**: bridge law 登録簿は全能力で空・PRED-019 未登録・
自然の観測量の的中 0・external_replications = 0・recorded_runs = 0。
blocker は v34.6 の受け皿に対する**外部の実施者と実データ**のみ。

**開発記録 (第三十四期)**:
- 反例の受理で期を開けた (v34.1 — 本プロジェクト初の外部起源 erratum)。
- 設計走行が器械の故障モードを 2 つ検出し凍結前に修復した (v35.0-A §3) —
  「設計走行はシード頑健性の確認」から「故障検出器」へ役割が広がった。
- 期を通じて v33 の 6 モジュールの lib pin が不変 — 第三十四期は第三十三期の
  器械に触れずにその上へ積んだ (加法的な期)。

**期末完全儀式 (追補で確定 — 2026-08-04)**: 全 211 本の完全再計算
(`make suite-full OUT=results/v350_full_suite.txt JOBS=12`, v35.0-B コミット後に起動):

```text
→ 実行 211 本: 総計 PASS 1547 / FAIL 0 (引用 0 — 全数再計算)
→ ドリフト検査 (台帳の期前後比較): 既存 203 本の PASS/FAIL は完全一致
  (ドリフト 0 件)・新規 8 本の +118 PASS (1429 + 118 = 1547 会計一致)・消失 0
→ 出力差分は壁時計の計時文字列のみ (数値は完全一致 — v336 の教訓どおり計時は
  合否に含まれない)。第三十四期の共有部変更 (finite_data /
  factorization_enumerator 新設・enumerator の凍結前修復) は既存物理に無波及
```

past 期は儀式を B コミットに同梱したが、本期は起動 (v35.0-B) と記録 (追補) を
分離して手続きを明示した — 儀式中の器械不変が git 履歴で機械検証できる。

## 6. 未解決 (第三十五期への課題 — PROMPT/15 の優先順位の残り)

1. **外部独立再現 (最優先・据え置き)** — D2R_PACKET の実配布と実施者獲得。
2. 実データ (recorded lane) の初回受理と事前登録予測。
3. profile の set-valued 関手性と安定性定理 (定義は v34.5 で器械化済み)。
4. 一般 GKLS 応答 (jump gauge・Kossakowski 同値類)・BCS 型 witness (本期は
   OutOfDomain 負制御まで)。
5. OCS-2.0 (有限データ意味論の spec 化 — v34.3–v34.5 の定理が材料)。
6. gravity・PRED-019 (凍結継続)。
