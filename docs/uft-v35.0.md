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

## 4. 開封 (v35.0-B — 開封時に追記)

*(凍結時点では空欄 — SECRET 開示・holdout 本採点・期末統合はここに追記される)*
