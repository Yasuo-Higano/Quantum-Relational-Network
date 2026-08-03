# QRN v34.1 — FollowUp 受理と Yukawa erratum: 反証の正式受理 (第三十四期 開幕)

**Version**: v34.1 (2026-08-03)
**Sim**: `sim/src/bin/v341_yukawa_erratum.rs` → `results/v341_yukawa_erratum.txt` (14 検査 PASS)
**位置づけ**: PROMPT/15 (第三十四期「申告証明書から観測証明書へ」) の必須第一課題。
新機能ではなく**反証の正式受理**から期を開ける — 第三十三期が「監査が版を開いた」
(v33.1) なら、第三十四期は「反例の受理が版を開いた」。

---

## 1. FollowUp 最終報告の受理 (replications.yml REP-001)

Quantum-Relational-Network-FollowUp (Julia, cross-model clean-room) は 2026-08-02 に
phase 7 の凍結 holdout 120 セルを**一度だけ開封**し、最終判定を確定した:

```text
総合判定: Partially Replicated
凍結: phase 6 (commit 5557eb9c…)・SECRET はコミットメント照合つき
120 セル: Answer 77 / EquivalenceClassOnly 12 / InsufficientObservation 11 /
          OutOfDomain 10 / Abstain 10
安全計量: selective risk 0・answerable recall 1・impossibility recall 1
最大数値誤差: ~2.78e-16
再現された能力: 有限次元応答・支持分解・代数診断 (中心/可換子)・homology
specification-limited: YUK-001..004 (格子行列・事前分布の仕様不足)・MSR/PREQ/
          CP-002,003/GRAV-002..007 など — paper だけでは独立実装に閉じない
```

**登録の裁定**: FollowUp 自身が「shared human operator による cross-model
clean-room であり、organizationally external human replication ではない」と
自己限定している。`replications.yml` の 6 条件 (different_author ほか) を
満たさないため **external_replications = 0 は維持** (v274 [P4]/v310 [R3] 監査下)。
ただし一次元の「external か否か」では情報を失うため、**多次元独立性 profile**
(replication_kind / human_operator / implementation_agent / repository /
numerical_kernel / generator / freeze_before_holdout / secret_holder /
organization / verdict) を v34.1 で導入し、REP-001 として登録した —
第三十三期テーゼ「属性は対象単体ではなく関係である」の再現独立性への適用。

**FollowUp が主リポジトリに返した宿題** (第三十四期の版計画に登録):
1. **YUK-005 反例** — 本版 (v34.1) で正式受理。下記 §2。
2. **operational core の spec 不在** — 多数の claim が specification-limited。
   paper だけで独立実装可能な core specification が存在しない → v34.2。
3. **probe 型の区別** — SignedInitialCovarianceProbe ≠ HamiltonianQuench
   (二準位反例あり)・NumberConservingResponse ≠ BCS/PairingResponse (応答則が
   1−Δ² に変わる) — 曖昧な「局所摂動」の語の再使用禁止 → v34.2 spec の型分離。

## 2. 反例の厳密受理 [E1] — 素数性は積を固定するが不等性を強制しない

v13.2 (QRN-YUK-015 発見 2) と論文 geometric-yukawa は

> 「Pf(F) = 3 (素数) が第 2 磁気固有面を平坦化し、2 タワーの抑制の掛け算を
> 構造的に禁じる」("index 3, being prime, forces one magnetic eigenplane to flatten")

という**普遍推論**を含んでいた。FollowUp の整数反対称行列

```text
F = ⎡  0   1   1   1 ⎤     Pf(F) = 1·1 − 1·(−1) + 1·1 = 3   (i64 厳密)
    ⎢ −1   0   1  −1 ⎥     FᵀF = 3·I₄                       (i64 厳密)
    ⎢ −1  −1   0   1 ⎥     det F = 9 = Pf²                  (i64 厳密)
    ⎣ −1   1  −1   0 ⎦  ⇒  f₊ = f₋ = √3 — 両固有面は等スケール
```

は 4×4 反対称の恒等式 f₊²+f₋² = Σ_{i<j}F²・f₊f₋ = |Pf| により両 skew singular
value がともに √3 — **素数性は積 f₊f₋ = 3 を固定するが不等性を強制しない**
(f± は整数である必要がない)。素朴推論の予測 (整数分解 {3,1} 型 ⇒ FᵀF 固有値
{9,9,1,1}) と厳密値 {3,3,3,3} の乖離 2.0 [E1e]。**普遍命題は 1 つの厳密証人で
反証される** — FollowUp が organizationally external でないことと反例の有効性は
別問題で、主リポジトリの整数再検算により数学的に決着した。

## 3. 位相 ↛ 計量 [E2] と族内定理 [E3] — なぜ誤り、なぜ族内観測は正しかったか

**[E2] Pfaffian は skew scale をどちらの向きにも決めない。** unimodular 剪断
S (det = 1) による同じ交代形式の別基底表示 F' = SᵀFS は Pf 不変 (= 3) のまま
Σ' = 8 → f₊ − f₋ = √2 (固有値 4±√7, jacobi 照合 < 1e-12)。同一の位相・代数
データ (Pf = 3) が等スケール (Σ = 6) とも不等スケール (Σ = 8) とも両立する —
**Pfaffian (基底不変の代数・位相量) から skew scale (計量依存の幾何量) への
昇格には metric/lattice compatibility bridge が必要**。

**[E3] 走査族の内側では不等が整数論的に強制される (窓なしの定理)。** v13.2 の
走査族 (座標 2-平面 + 傾き対, Pf = Q₁Q₂ + ts, f13 = f24 = 0) では:

```text
等スケール ⟺ Σ = a²+b²+c²+d² = 2|Pf| = 6 — この整数解は存在しない
  (Σ = 6 ⇒ |成分| ≤ 2 なので全数 96 点で完全列挙 — 窓仮定なし)
⇒ 族内 Pf = ±3 なら Σ ≥ 7 (AM-GM の 6 と整数性) ⇒ (f₊−f₋)² = Σ − 6 ≥ 1
最小ギャップ: (2,1,1,1) 型 16 点で f₊ − f₋ = 1 厳密 (f± = (√13±1)/2)
v13.2 の 7 走査点: gap² = {4, 4, 13, 49, 25, 25, 49} — 全て深い不等
反例 F: 6 成分全て非零 — 族の外 (f13 = 1, f24 = −1 ≠ 0)
```

つまり **族内観測 (QRN-YUK-015) は正しく、その普遍化だけが誤りだった**。
「素数だから」に見えた平坦化は、実際には走査族のパラメータ化 (2 成分を 0 に
固定) の整数論的帰結であり、素数性は無関係 (どんな固定積でも Σ ≥ 2|Pf| は
成り立ち、等号到達可能性だけが族に依存する)。

## 4. 台帳と論文の分割 [E4]

| 対象 | 処理 |
|---|---|
| QRN-YUK-015 | **走査族限定に分割** — 発見 2 を「admissible family の範囲内の有限結果 + 族内不等定理」に書き換え、原文 (普遍形) は正誤表 limitation に逐語保存。status の「最終決着」→「走査族内の決着」 |
| QRN-YUK-034 (新規, C2) | **refuted_as_stated の正式記録** — 反例の厳密検証 [E1]・分離の両方向 [E2]・族内定理 [E3] を 1 claim に登録 (evidence_kind: theorem・independence: algorithmically_diverse) |
| QRN-META-013 | Erratum limitation — 「傾き代替が素数 3 の構造論で落ちた」を「走査族の有限結果で落ちた」に訂正 (住所替えの結論は走査済み範囲で維持) |
| paper/geometric-yukawa-full.md | Abstract 差し替え (PROMPT/15 の推奨文)・§5c 再構成 (「why the prime 3 kills it」→ 有限走査 + erratum)・Erratum 節新設 (撤回文は blockquote 引用のみ — [E4b] が機械検査)・§6 限界の更新 |
| paper/geometric-yukawa.md | Abstract 案・§5c・claims 表の同期訂正 |
| replications.yml | REP-001 (cross_model_clean_room / partially_replicated / counts_as_external_replication: false)・多次元独立性 profile スキーマ |
| assumptions.yml | ASM-EXACT-INTEGER (i64 厳密性と列挙完全性の信頼) |

## 5. 一般教訓と残高

**一般教訓 (第三十四期の型設計への入力)**: 「位相・代数的不変量 → 計量依存の
幾何量」の無証明昇格は、QRN 全体で繰り返し警戒すべき族である (今回: Pf → skew
scale。既知の同族: 大域閉包 → marking [v32.2]・宣言 → 資格 [v33.2])。v34.2 の
core spec と v34.3 の禁止変換 22–29 に反映する。

**正直な残高 (不変)**: bridge law 登録簿は全能力で空・PRED-019 未登録・自然の
観測量の的中 0・**external_replications = 0** (REP-001 は cross-model clean-room
であり 6 条件を満たさない — FollowUp 自身の自己限定と一致)。

**開発記録**: 反例の受理で期を開けるのは本プロジェクト初 (これまでの正誤表
v9.2/v24.6/v29.1/v32.1 は全て内部発見)。「有効報告は一件で足りる」(Track X の
凍結された約束) の反例版が機能した — 報告の出自が cross-model でも、厳密証人は
整数再検算で主リポジトリの資産になる。

## 6. 次 (v34.2)

standalone operational core paper/spec — 第三十三期の certified interface /
context atlas / resource profile / graded recovery を、リポジトリの型名を知らずに
独立実装できる形で定義する論文。observation contract・裁定順序・全 no-go・
falsifier・publication closure manifest (paper_closed / repository_replay_only の
区別) を含み、結果値は未掲載とする (FollowUp の specification-limited 判定への
応答)。probe 型の分離 (SignedInitialCovarianceProbe ↛ HamiltonianQuench /
NumberConservingResponse ↛ BCSResponse) を明文化する。
