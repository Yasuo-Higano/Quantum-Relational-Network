# QRN v33.4 — Contextual factorization: chart 局所復元と overlap glue

**Version**: v33.4 (2026-08-02)
**Sim**: `sim/src/bin/v334_contextual_factorization.rs` →
`results/v334_contextual_factorization.txt` (7 検査 PASS) / 共有契約 =
`sim/src/contextual_factorization.rs` (新設)
**位置づけ**: PROMPT/14 第四版。実際の実験室は**持ち場 (chart)** の集まりであり、
全系を一望する制御器はどこにもない。期テーゼの「role-typed context atlas 上で
**整合する**因子分解」の「整合」を機械化する: context ごとの局所復元 → overlap 上の
algebra matching → glue 裁定。

---

## 1. chart 局所復元 — 持ち場しか語らない

`ChartSpec` (実験者が宣言する primitive の部分集合) の局所復元 (凍結手順):

1. chart 部分 net の非可換グラフ成分 — **可換子証明書は大域 net から継承**する
   (証明書は測定データであり、部分 net で再計算・捏造しない)。
2. **chart 内証人ゲート** (v33.1 の局所版): 成分対ごとに chart 内共有文脈。
3. 各成分の閉包 = full matrix factor M_d (d ≥ 2)・成分中心自明。

出力 `ChartLocalFactorization` は chart 内の因子 (次元・traceless 部分代数・構成
primitive) のみ — 大域 fullness は chart の資格要件ではない。**chart の局所 Exact を
大域因子分解に昇格する変換は存在しない (禁止変換 19)**。

## 2. atlas glue の凍結手順と glue 定理 [A1][A2]

成分ノード (chart, factor) 対の**連結** = 共有 primitive ∨ 交差 NonCommuting
(certified)。連結対は**同一の因子部分代数** (traceless ONB の overlap ≥ 1 − 1e-9)
を指していなければならない (overlap 上の algebra matching):

- **全 matching 整合** → 因子クラスを束ね、次元一致・被覆 (Π d = n)・**大域証人**
  (クラス対ごとに共有文脈 — v33.1 の証人規律の atlas 版) を検査 → **GluedExact**。
- **glue 定理 (機械照合)**: chart A = {X₁,Z₁,X₂,Z₂}・B = {X₂,Z₂,X₃,Z₃} (qubit 2
  共有・bridge 文脈 {X₁,X₃}) の atlas glue [2,2,2] は、**直接大域復元 (v33.1 入口)
  と読み・gauge orbit が一致**する。
- **変成不変**: chart B の qubit-2 操作を局所 unitary u₂ で回した net でも glue は
  同じ [2,2,2] — 因子部分代数は集合として frame 非依存 (u M₂ u† = M₂)。

## 3. cocycle 不整合 → Abstain [A3]

chart B を entangler W = CZ₂₃ で捻ると、B の因子は W M₂(2) W† になり A の M₂(2)
との overlap は 1/3 — matching が破れる。A も B も**局所的には Exact** (因子 2 つ
ずつ) なのに、全被覆 (Π d = 8) を達成する整合 chart 群が存在しないため
**Abstain(GlueInconsistent)**。宣言された chart を黙って捨てて残りで「大域」を
名乗ることはしない。

## 4. 複数 glue → EquivalenceClassOnly [A4]

site charts (S1+S2) と DFT chart D の atlas: site 群と D 群は**それぞれ全被覆で
内部整合**・相互の matching は破れる (overlap ≈ 0.56) — 整合 atlas が 2 つあるとき
**EquivalenceClassOnly{2}** を返し、無制約 tie-break で 1 つを選ばない (v32.3 [F3]
site×DFT 裁定の atlas 版)。

## 5. witness 境界の両 lane 一致 [A5]

bridge 文脈 {X₁,X₃} を外すと: glue は Abstain(CompatibilityUnwitnessed)・直接大域
復元も Abstain(OperationalCompatibilityUnwitnessed) — **v33.1 の証人規律は atlas を
経由しても緩まない**。glue が資格の抜け穴にならないことの機械記録。

## 6. 型契約の登録 [A6]

- `sim/src/contextual_factorization.rs` — `ChartSpec` / `ChartLocalFactorization` /
  `ChartFailure` 4 種 / `AtlasReading` (GluedExact / EquivalenceClassOnly / Abstain
  4 理由)。`OperationalNet::commutator_certificate` (証明書の継承用アクセサ) を
  追加 (v32.2 契約は不変)。
- `core.schema.yml` に概念 3 種 + **禁止変換 19** を登録。

## 7. 正直な残高

- chart 復元は control lane の Exact 因子まで — chart 内の**中心非自明**
  (superselection) atlas・graded chart は未実装 (前者は sector つき matching、
  後者は v33.5 の主題)。
- matching バーは exact 証明書域 (overlap ≥ 1 − 1e-9) — ノイズ下の overlap 区間と
  Straddled 裁定は HOLD-9 変成セルの主題。
- 「複数 glue」の数え上げは chart 群の連結分割に基づく凍結規則 — 一般の極大整合
  部分 atlas の列挙 (指数的) はスコープ外 (toy の chart 数 ≤ 3 で厳密)。
- MeasurementContext / PreparationFamily / DriftRegime (v33.2) の atlas への統合は
  据え置き — 本版の chart は control 文脈の集まり。
- bridge law 登録簿は全能力で空・PRED-019 未登録・自然の的中 0・
  `external_replications = 0` のまま — Unit D2-R の公募が引き続き最優先。
