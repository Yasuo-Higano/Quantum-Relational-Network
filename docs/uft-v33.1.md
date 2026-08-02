# QRN v33.1 — 境界監査と型スコープ修復: contexts は定理入力・accessible primitive は未証明

**Version**: v33.1 (2026-08-02)
**Sim**: `sim/src/bin/v331_scope_repair.rs` → `results/v331_scope_repair.txt`
(7 検査 PASS) / 共有契約の増補 = `sim/src/operational_net.rs` (v33.1 節)
**位置づけ**: PROMPT/14 第一版 — 第三十三期の開幕。期テーゼ:

> **可アクセス性は作用素単体の属性ではなく、系・制御器・測定器・資源・誤差証明書の
> 関係である。局所性は、証明付き laboratory interface が生成する role-typed context
> atlas 上で整合する因子分解の、資源スケールにわたる安定な同値類としてのみ識別される。**

第一歩は新機能ではなく**境界監査**である: v32.3 の復元器が実際には何を入力に
取っていたかを反例で固定し、型スコープの空隙を修復する。旧署名は
`v323_factorization_recovery.rs` に凍結保存し (黙って直さない — v32.1 の
superseded 原則)、修復は `operational_net.rs` の新入口として並置する。

---

## 1. 監査 — v32.3 復元器の型スコープの空隙

### 1.1 contexts 盲目性 [B1]

`OperationalNet` は `contexts` (可換分解の区画) を保持し `add_context` は可換性を
検証するが、**v32.3 の復元手順は contexts を参照しない**。機械実証: 同一 primitive・
同一可換子証明書で contexts だけが異なる (∅ vs atlas {X₁,X₂,X₃}/{Z₁,Z₂,Z₃}) 2 つの
site net に対し、v32.3 参照手順 (旧署名の再現 — 歴史的原本は v323) は非可換グラフ
成分も最終読みも**完全一致** (両方 Exact [2,2,2])。contexts は入力でなかった。

これは**現行定理の誤りではない** — 定理の仮定「各ノードの局所生成子だけが選別
された primitive family」が型に存在せず、復元対象行列が net と**別引数 `gens`**
として並行に渡され、両者の整合が呼び出し側の慣行任せだった、という型スコープの
空隙である。

### 1.2 primitive 選別の循環 [B2]

site 6 primitive {X_i, Z_i} → Exact [2,2,2]。ここに **independently accessible な
entangler X₁X₂ を 1 本加えるだけで読みは [2,4] に併合**する (新旧 lane とも一致・
大域閉包は同一の M₈ [dim 64] のまま)。細粒度の [2,2,2] を保持すべきか [2,4] へ
粗視化すべきかは、操作コスト・文脈・合成経路なしには決まらない — つまり

> **primitive の選別自体が答えを入力してしまう循環**

が第三十三期の敵である。この循環は本版の型修復では解けない (どの操作が
physically accessible かの証明 — `DeclaredOperation → AccessibleOperation` の資格 —
は v33.2 Certified Laboratory Interface の主題)。

## 2. 修復 — MarkedRecoveryInput (gens 別渡しの廃止) [B3]

復元器の唯一の型付き入口を `MarkedRecoveryInput` に固定する:

- **gens 別渡しの廃止**: 生成子行列は net 自身の primitive からのみ取る
  (`generator_matrices` — 外部から行列を注入する経路は存在しない)。公開
  コンストラクタは `OperationalNet::recovery_input` のみ。
- **構成時資格審査** (拒否は型エラーであって Abstain ではない — 「不正な入力で
  走らない」と「走った上で棄却する」を混ぜない):
  - `RoleMixedRecovery` — Control 以外の役割 (測定 effect 等) の混入を拒否
    (**role-mixed recovery の禁止**)。測定・準備・drift の文脈意味論 (joint
    measurability 等) の型化は v33.2。
  - `NoDeclaredContexts` — 文脈 0 を拒否。**contexts が復元定理の入力になった**。
  - `ContextCoverageIncomplete` — 文脈 atlas が primitive family を被覆しない
    限り走らない。
- 資格を満たす net では v32.3 の凍結決定手順 (成分 → joint 閉包 → 中心 → 三値
  裁定・dust guard 込み) を再現する: site → Exact [2,2,2]・qutrit×qubit (C⁶) →
  Exact [2,3]。

## 3. 禁止変換 12 — 代数的可換 ↛ 操作的両立 [B4]

「context = 可換集合」は役割ごとに分解すべき概念であり、その第一歩として
**代数の層と操作の層を型で分離**する:

- `CertifiedCommutator` (bracket ノルムの区間証明書) は**代数的**事実。
- `JointContextWitness` (成分対の共同 addressability) は**操作的**事実 — 唯一の
  構成は宣言済み文脈の共有 (`joint_context_witness`)。一般 POVM の joint
  measurability は可換性で特徴づけられないため、前者から後者への暗黙変換は
  存在しない (**禁止変換 12**)。
- 機械実証: singleton 文脈のみの site net は**全対に definite な Commuting/
  NonCommuting 証明書があっても**、成分間の証人がなく
  `Abstain(OperationalCompatibilityUnwitnessed)` — 同じ素材を旧手順は Exact
  [2,2,2] と読む (修復が変えた点の対照)。

## 4. 禁止変換 13 — Liouvillian lane の型分離 [B5]

v32.4 の実装対象は L = −i[H,·] (Hamiltonian commutator lane) であり、一般の GKLS
生成子ではない。これを型で固定する:

- `HamiltonianCommutatorLiouvillian` — v32.4 応答恒等式 R⁽¹⁾ = −i Tr(B[H,A])・
  R⁽²⁾ = Tr([H,B][H,A]) の資格域。lane への資格審査は**導分 (Leibniz) 証明書**
  (`classify_generator`): L(AB) = L(A)B + A·L(B)・L(A†) = L(A)†・L(I) = 0 の全対
  検査。可換子 lane は全て**厳密 0** で通過し、Ĥ を中心を除いて一意に復元する
  (traceless gauge・復元残差 厳密 0)。
- `GklsLiouvillian` — GKLS 生成子 (γ_μ ≥ 0 を構成時資格審査)。**dissipator は
  Leibniz を破る**: 破れは γ に厳密比例 (γ = 0.15 → 0.3 で比 2.000000)・γ = 0 の縮退点は
  可換子 lane の資格を通る (連続性の対照)。
- **応答負制御**: R⁽¹⁾ 公式を GKLS 発展の測定 (RK4 + 4 次 stencil) に当てると
  γ = 0.3 で乖離 1.200 (測定自身は GKLS の線形応答 Tr(B L(A)) と ≤ 1.2e-9 で
  一致 — 乖離は数値誤差ではない)・γ = 0 で ≤ 1.3e-9 に回復。
- 帰結: **可換子 lane の応答証明書は GKLS lane に昇格しない** (**禁止変換 13**)。
  一般 GKLS の応答理論 (jump 表現の gauge — 推定対象は個々の L_μ でなく
  superoperator/Kossakowski 表現の同値類) は未構成で、次々期の独立テーマ。

## 5. 型契約の登録 [B6]

- `sim/src/operational_net.rs` (v33.1 節) — `MarkedRecoveryInput` /
  `RecoveryInputRejection` (3 種) / `JointContextWitness` /
  `FactorizationAbstainReason::OperationalCompatibilityUnwitnessed` /
  `HamiltonianCommutatorLiouvillian` / `GklsLiouvillian` / `classify_generator` /
  v32.3 kernel の lib 移植 (`closure_center_basis` / `closure_central_projectors` —
  dust guard 込み)。v32.2 契約 (禁止変換 11・odd 拒否ゲート) は不変。
- `core.schema.yml` に概念 4 種 + **禁止変換 12/13** を登録。
- 禁止 impl From の source 監査 (CertifiedCommutator / Liouvillian 両型 /
  MarkedRecoveryInput / JointContextWitness) — v331 [B6] が常設検査。
- v323 の旧署名 (`gens: &[Vec<C64>]`) は**凍結保存** — 歴史的原本として監査対象。

## 6. 正直な残高

- **primitive 選別の循環は未解決** (本版は境界を型で固定しただけ)。accessible
  operations の出自 — 校正・独立 addressability・合成経路・資源上限の証明書 —
  は v33.2 (Certified Laboratory Interface)、資源依存の factorization profile は
  v33.3 の主題。
- 本版の witness 意味論は最小形 (成分対ごとに共有文脈 1 つ) — context overlap の
  整合・cocycle・複数 glue の裁定は v33.4 (Contextual Factorization)。
- role-mixed recovery は拒否したが、`add_context` 自体の役割型付け
  (ControlContext / MeasurementContext / PreparationFamily / DriftRegime) は
  v33.2 で別型化する (v32.2 の add_context は可換性のみ検査 — 凍結のまま)。
- graded lane の recovery (odd CAR → Majorana frame / Dirac no-go) は v33.5。
- bridge law 登録簿は全能力で空・PRED-019 未登録・自然の的中 0・
  `external_replications = 0` のまま — **Unit D2-R の公募が引き続き最優先**
  (Track X: campaign layer は本期の並行トラック)。
