# QRN v34.3 — 有限データ昇格不能定理 (第四の no-go) と Robust Promotion Theorem

**Version**: v34.3 (2026-08-03)
**Sim**: `sim/src/bin/v343_finite_data_nogo.rs` → `results/v343_finite_data_nogo.txt`
(17 検査 PASS)・共有契約 `sim/src/finite_data.rs` 新設
**Lean**: `proofs/FiniteDataNoGo.lean` (10 定理 — 計 87 定理 / 11 ファイル)
**位置づけ**: PROMPT/15 §4。期テーゼ「可アクセス性証明書は宣言ではなく同時信頼
集合である」の定理層。exact な第三十三期の型と定理を**現実の裁定に使うために何が
観測されればよいか**を閉じ、閉じられない領域を第四の no-go として切り出す。

---

## 1. 第四の no-go: 有限データ昇格不能定理

登録観測モデルのパラメータ θ・有限データ分布 P_θ・exact 裁定 q(θ) について、
異なる裁定を要する θ₀, θ₁ の分布が近いとき:

```text
inf_δ [ P_θ₀(δ ≠ q(θ₀)) + P_θ₁(δ ≠ q(θ₁)) ] / 2  ≥  (1 − TV(P_θ₀, P_θ₁)) / 2
```

**Lean 機械証明** (`proofs/FiniteDataNoGo.lean` — 有限結果空間・整数重み):

| 定理 | 内容 |
|---|---|
| errMass_ge_minSum | どんな決定規則も点ごと最小の和を下回れない |
| two_minSum_eq | 2 Σmin = W₀ + W₁ − Σ\|w₀−w₁\| (TV との恒等式) |
| le_cam_two_point | 上 2 つの合成 = Le Cam 二点下限 |
| bayes_achieves | 尤度比規則が下限を**達成** — 下限は最良 (改良不能) |
| indistinguishable_half | P₀ = P₁ ⇒ 平均誤り ≥ 1/2 (観測契約が区別しなければ当てられない) |
| promotion_exclusion | 誤昇格 ⇒ 真値は信頼集合の外 (排除補題) |
| wrongMass_le_missMass / robust_promotion | 誤昇格質量 ≤ 被覆失敗質量 ≤ α |
| instance_n4_min / instance_n4_all_rules | N=4 実例 (整数 147/512)・全 32 規則の decide 全数 |

**Rust 側の三重一致** [F1]: N=16 の**全 131072 決定規則の全数列挙**の最小平均誤り
= Σmin/2 = (1−TV)/2 (1e-14)。Lean 実例とは整数 (147, 218) で一致 — 形式証明と
数値器械の橋。[F2]: P₀ = P₁ で誤り = 1/2 厳密・TV(N) は単調増加 (0.07 → 0.52,
N = 4 → 256) — **境界の分解能は観測量の関数であり、境界近傍では abstention か
追加実験だけが正答**。

## 2. Robust Promotion Theorem (正側)

データ D から真値を確率 ≥ 1−α で含む同時信頼集合 C_α(D) (Clopper–Pearson —
二項 tail の厳密反転) を構成し、**点推定ではなく集合全体**を exact reader に通す:

```text
Q_α(D) = { q(θ) : θ ∈ C_α(D) }
  単一クラス            → RobustExact      (誤昇格 ≤ α — Lean: robust_promotion)
  既知同値関係の単一クラス → EquivalenceClassOnly
  裁定境界を跨ぐ         → Straddled       (強制回答の禁止)
  観測が最低量未満        → InsufficientObservation
  契約外データ           → OutOfDomain
```

機械検証 (全て厳密和 — Monte Carlo なし):
- [F3a] 全 θ 走査で P(wrong promotion) ≤ α — 実測 max 0.0196 (片側 α/2 以下)。
- **[F3b] selective risk の区別 (反例)**: θ = 0.305 で P(wrong) = 0.011 ≤ α なのに
  **P(wrong | answer) = 0.318 > α** — 被覆単独から回答条件付き保証は出ない
  (P(wrong|answer) ≤ α/P(answer) — 分母の下限が別途要る)。HOLD-10 の採点設計
  (selective risk = 0.000 だけを完成条件にしない) の定理的根拠。
- [F5b] **裁定境界上では P(Straddled) ≥ 1−α** — 棄権が保証つきの正答 (被覆の系)。
- [F5c] 符号 orbit (契約が |θ| しか見ない — TV(P₊,P₋) = 0): 回答は
  EquivalenceClassOnly のみ (クラス誤り ≤ α)・強制符号回答は [F2a] より誤り ≥ 1/2。

## 3. 禁止変換 22–29 (型 `sim/src/finite_data.rs` + 厳密反例)

| # | 変換 (存在しない) | 反例 (厳密和) |
|---|---|---|
| 22 | PointEstimate ↛ AccessibilityCertificate | 境界近傍 θ=0.28 で点推定昇格の誤り 0.312 ≫ α (robust は 0.005) |
| 23 | MarginalIntervals ↛ JointConfidenceRegion | 周辺 95% × 6 の直積: joint 被覆 0.845 < 0.95 (Bonferroni α/m は 0.984) |
| 24 | GoodnessOfFit ↛ NoiseModelValidity | 同一平均の beta-二項 (一次積率は厳密一致) で被覆 0.329 |
| 25 | CalibrationAt(t₀) ↛ ValidAt(t₁) | drift θ 0.25→0.40 で t₁ 被覆 0.104 |
| 26 | MeanCrosstalk ↛ UniformCrosstalkBound | 平均 0.058 ≤ 0.1 でも単一対 0.30 が独立 addressability を破る |
| 27 | LocalChartCoverage ↛ GlobalGlueCoverage | 23 と同型 (glue は全 overlap の同時被覆) |
| 28 | ZeroHoldoutErrors ↛ ZeroPopulationRisk | 誤り 0 の片側 95% 上限: 0/9 → **28.3%**・0/77 → **3.8%**・≤1% には **299 セル** |
| 29 | ModelConditionalCertificate ↛ ModelFreeCertificate | 同一周辺の相関鎖 (ρ=0.8) で被覆 0.513 — iid は「登録」する仮定 |

28 は PROMPT/15 §6 の数の機械化 — HOLD-9 の 9 answerable / FollowUp の 77 Answer が
「観測誤り 0」であって母集団リスク 0 でないこと、HOLD-10C が 300+ 回答セルを
要する理由の一次ソース。

## 4. 接続

- **OCS-1.0 との関係**: spec §14 が「有限データ意味論はスコープ外・定理の後に
  spec 化」と宣言済み — 本版がその定理。v34.5 で実 reader (addressability σ_min・
  cross-talk・glue overlap・J) を信頼集合に持ち上げ、OCS-2.0 系の材料にする。
- **第三十三期との関係**: v33 の Straddled (区間がバーを跨ぐ) は「観測された区間」
  への棄権だった。本版はそれを「真値の信頼集合」への棄権に**根拠づけ**る —
  跨ぎ棄権は美徳ではなく、強制回答の誤り下限 (no-go) からの必然。
- **残高**: bridge law 空・PRED-019 未登録・自然の的中 0・external 0 — 不変。

## 5. 次 (v34.4)

sector-aware complete finite factorization enumerator — A, A′, Z(A) から central
projector 列挙 → Wedderburn block 証明 → multiplicity 分離 → 候補列挙 + 同値
witness。出力型 UniqueFactorization / FactorizationCandidateSet /
SectorwiseFactorization / IncompletePrimitiveSet / NontrivialCenterObstruction /
ScopeExceeded。
