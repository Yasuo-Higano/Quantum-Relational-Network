# QRN v33.2 — Certified Laboratory Interface: 操作の出自証明と controller-free no-go

**Version**: v33.2 (2026-08-02)
**Sim**: `sim/src/bin/v332_certified_interface.rs` → `results/v332_certified_interface.txt`
(9 検査 PASS) / 共有契約 = `sim/src/laboratory_interface.rs` (新設)
**位置づけ**: PROMPT/14 第二版。v33.1 が固定した循環 (primitive 選別が答えを入力
する) の正側修復。目標は「物理的に可能な操作を状態や Hamiltonian だけから自動導出
する」ことでは**ない** — それは原理的に不可能であり、本版 [C5] が第三の no-go
として機械化する。目標は

> 「accessible operations を入力として与えた」から、「各 operation が、どの
> command・校正・合成列・誤差・資源によって accessible なのかを証明した」へ

進めることである。

---

## 1. 宣言 ≠ 資格 — 出自の 3 つの門 (禁止変換 14) [C1]

`DeclaredOperation` (行列 + 意図) は**資格なし** — `AccessibleOperation` への直接
変換は存在しない (**禁止変換 14**)。通過できる門は出自証明書 3 種のみ:

| 門 | 証明書 | 機械検証の内容 |
|---|---|---|
| 較正 | `IndependentAddressabilityCertificate` | command Jacobian の証明つき rank・σ_min 区間・cross-talk 区間・較正記録 sha256 |
| 合成 | `SynthesisCertificate` | bracket/線形結合の列の**機械実行** (レシピの主張ではなく実行結果が資格)・深さ・残差区間 |
| トモグラフィ | `TomographyCertificate` | 情報完全状態族 (design σ_min 資格) からの線形再構成・残差区間・effect 資格 |

- **sha256 結束**: 各証明書は対象行列の sha256 を記録し、別の行列への流用は
  `CertificateTargetMismatch` で構成時拒否 (X₁ の証明書を X₂ に付けられない)。
- **文字列 provenance の廃止**: 新 lane の出自は型 + sha256 で運ぶ。v32.2 凍結層の
  `PrimitiveOperation.provenance` 文字列は legacy フィールドとして定数で中立化。
- `ResourceBudget` (時間・振幅・帯域・深さ・誤差の 5 成分) は**成分半順序**のみ —
  Ord/PartialOrd を実装しない (恣意的な重み付き和での全順序化は新しい選択バイアス)。
  比較不能対 (1,2,·) vs (2,1,·) の実在を機械記録。

## 2. 数学的分解は独立 addressability を与えない (禁止変換 15) [C2][C3]

較正の資格審査 (`certify_addressability`): 標的族を HS 正規直交化し、正規化
command との重なり行列 M̂ の**全特異値 ≥ σ_bar** (証明つき rank) と cross-talk
(自標的以外への漏れ) の区間バーを検査する。

- 独立 site knobs (6 command 1:1, dim 8): rank 6・σ_min = 1.000000000・
  cross-talk ≤ 1e-12 → 資格。
- cross-coupling ε = 0.05 (G = X_k + εX_{k+1}): 記録つき資格 (cross-talk 0.049938
  ≤ バー 0.1)。ε = 0.3 → `CrosstalkExcess` 拒否。区間がバーを跨ぐ場合は
  `CrosstalkMarginStraddled` (強制判定の禁止 — HOLD-9 の cross-talk セルの器械)。
- **tied control no-go**: 装置が u(t)(X₁+X₂) しか持たないとき、{X₁, X₂} への
  数学的分解の較正申告は **rank 1 < 2 で構成時拒否** (**禁止変換 15**)。正直な
  interface (tied 1 本を 1 標的として較正) の net の読みは
  Abstain(InsufficientOperationalGenerators) — 「二つの独立 primitive」は立たない。

## 3. 可アクセス性は interface との関係 — 合成の門 [C4]

同じ作用素 X₁ が:

- **interface A = {X₁+X₂}**: Lie 閉包 (bracket + 線形結合) への相対残差 0.707107 —
  合成路なし (`NoSynthesisPath` 域)。
- **interface B = {X₁+X₂, Z₂}**: (1/i)-bracket 2 手 + 線形 1 手 (depth 3) で
  X₁ を合成 — 残差 ≤ 1e-12 の機械実行検証つきで `Synthesized` 資格。

**可アクセス性は作用素単体の属性ではない** (期テーゼの直接の機械実証)。tied だから
永久に不可なのではなく、合成路が実在すれば証明書つきで資格が立つ — 禁止されるのは
**証明書なしの数学的分解**である。

## 4. controller-free decomposition no-go — E3-A (第三の no-go) [C5]

> 同一の (H, H, ρ) に対して、異なる laboratory interface が非同値な OperationalNet
> と因子分解を生成し得る。したがって状態・Hamiltonian・大域作用素代数だけから
> 「物理的に accessible な操作」を一意に選ぶ写像は存在しない。

機械実証 (dim 8・drift H = 0・ρ = I/8 を固定し、certified interface だけを変える):

| interface | 読み |
|---|---|
| site knobs {X_i, Z_i} | Exact [2,2,2] — gauge orbit α |
| DFT knobs {VX_iV†, VZ_iV†} | Exact [2,2,2] — gauge orbit β (matching 不在, overlap 0.5625) |
| site + entangler X₁X₂ | Exact [2,4] |
| tied (X₁+X₂) | Abstain(InsufficientOperationalGenerators) |

三段の no-go が完成: **状態単独では選べない** (v31.4)・**global closure は
marking を消す** (v32.2 禁止変換 11)・**controller を消すと accessibility が消える**
(本版)。

## 5. context は役割ごとに別型 — 禁止変換 16 [C6]

「context = 可換集合」は役割意味論を混ぜる。v33.2 で 4 つの role-typed 文脈を
別型化し、いずれも代数的可換性 (CertifiedCommutator) から構成**できない**
(**禁止変換 16** — 禁止変換 12 の役割別展開):

| 文脈 | 要求する証明書 | 機械実証 |
|---|---|---|
| `ControlContext` | IndependentAddressabilityCertificate | 証明書が member を sha256 結束しない登録は拒否 |
| `MeasurementContext` | JointMeasurementCertificate | **joint measurability は可換性より広い**: 非可換 unsharp 対 (η = 0.6, ‖[E^X, E^Z]‖ = 0.2546 ≠ 0) が明示的 joint POVM (各元 PSD・総和 I・marginal 厳密一致) で資格 — η = 0.8 は正値性破れで拒否 (不偏 qubit 対の canonical 構成 G±± = ¼(I ± ηX ± ηZ) は Busch の iff を飽和する器械) |
| `PreparationFamily` | ConvexReachabilityCertificate | z = 0.2 は重み (0.375, 0.625) で機械検証・z = 0.9 (凸包外) は拒否 |
| `DriftRegime` | StabilityCertificate | 変動 0.02 ≤ バー 0.05 は資格・0.5 は拒否 |

測定側の出自 = トモグラフィの門: 情報完全 6 状態 (design σ_min 資格) から effect
を再構成 (残差 ≤ 1e-12)・偏りデータ (+0.02) は残差バーで拒否。

## 6. end-to-end と型契約の登録 [C7][C8]

- `AccessibleOperationalNet` — **DeclaredOperation を受け付ける口が存在しない** net。
  全 primitive が出自証明書つきで入り (admit のみ)、制御文脈は addressability
  証明書つきで v32.2 の可換子証明書要求に加えて登録され、v33.1 修復入口
  (contexts 必須・role 純度・被覆) で Exact [2,2,2] を復元する。
- `core.schema.yml` に概念 13 種 + **禁止変換 14/15/16** を登録。
- source 監査: `impl From<DeclaredOperation…>` / `impl Ord for ResourceBudget` /
  `impl PartialOrd for ResourceBudget` の不在。

## 7. 正直な残高

- 本版の較正・合成・トモグラフィ証明書は **exact toy** (誤差区間は ±1e-12 級の
  数値幅) — 実測定ノイズ下の区間較正・shot ノイズからの資格判定は未走査
  (HOLD-9 の変成セルで走査する)。
- ResourceBudget は AccessibleOperation に付くだけで、**budget を動かしたときの
  読みの変化 (resource-indexed factorization profile) は v33.3** の主題。
  合成の深さ・時間コストと budget の相互作用もそこで扱う。
- MeasurementContext / PreparationFamily / DriftRegime は型つき保持まで — 因子分解
  への統合 (context overlap の整合・glue) は v33.4。
- joint measurability の拒否は canonical 構成の正値性破れ (この構成では立たない)
  であり、一般の非存在証明は不偏 qubit 対の Busch の iff (教科書的事実) に依る。
- bridge law 登録簿は全能力で空・PRED-019 未登録・自然の的中 0・
  `external_replications = 0` のまま — Unit D2-R の公募が引き続き最優先。
