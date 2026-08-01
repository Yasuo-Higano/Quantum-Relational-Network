# QRN v32.2 — OperationalNet: 操作的文脈の型分離と global-algebra erasure no-go

**Version**: v32.2 (2026-08-01)
**Sim**: `sim/src/bin/v322_operational_net.rs` → `results/v322_operational_net.txt`
(8 検査 PASS) / 共有契約 = `sim/src/operational_net.rs`
**位置づけ**: PROMPT/13 第二版 — 第三十二期の主研究の開幕。v31.4 の E3 no-go
(状態単独では因子分解を選定できない — 因子分解は操作代数が運ぶ) の「操作代数」を、
文字列リストの雛形 (旧 OperationalAlgebra) から **marked family of subalgebras
(OperationalNet)** に精密化し、**第二の no-go** を機械化する。

---

## 1. 目標定理 A — global-algebra erasure no-go

> 異なる tensor factorization に対応する二つの operational generator family が、
> **同一の大域閉包代数 B(H) を生成し得る**。したがって大域代数の同型類だけから
> 因子分解は識別不能である。

**機械実証 [N1]** ((C²)⊗³, dim 8): site 生成族 {X_i, Z_i} と DFT₈ 共役族
{V X_i V†, V Z_i V†} の *-閉包 (HS Gram–Schmidt 成長 — 決定的) は**ともに複素次元
64 = M₈ 全体**。一方、su(2) 部分代数の対応は存在しない — subspace overlap は
対角自己 1.000000000000 に対し site×mode 最大 **0.617851** (< 0.9 バー)。
同一の閉包・異なる因子分解。

**帰結の型化 (禁止変換 11)**: `GlobalClosure` は閉包の同型類 (次元) しか運ばず、
OperationalNet・因子分解への変換は存在しない — 「閉包は marking を消す」が
コンパイル不能性として固定される (v322 [N6] が impl From の不在を source 監査)。

**[N2] marking は因子分解を運ぶ** (v32.3 recovery の最小 preview): 両 net の
非可換グラフ連結成分 = 3 × {X_i, Z_i} → 各成分閉包 = M₂ (dim 4)・local_dims
[2,2,2] — **同一の閉包から、net ごとに別の因子分解が復元される**。v31.4
(同一状態 → 3 幾何) の状態側 no-go と対になり、「何が因子分解を運ぶか」の答えが
**OperationalNet (marking)** に確定した。

## 2. 役割は同じ数学的型ではない [N3]

| 役割 | 型 | 閉包構造 | 機械証拠 |
|---|---|---|---|
| 準備 | `Preparation` (tr = 1・PSD を構成時資格審査) | **凸結合**で閉じる | mix は資格通過・行列積は tr = 1/8 で棄却 |
| 介入 | `ControlGenerator` (エルミート) | **Lie bracket** で閉じる | [X₁,Z₁]/i = −2Y₁ ∈ su(2)₁ (残差 0) |
| 測定 | `MeasurementEffect` (0 ≤ E ≤ I) | 作用素系 — **積閉包を要求しない** | span{I, n₁, n₂} = 3 → 積 n₁n₂ で 4 |
| drift | `DriftGenerator` (制御不能な発展) | — | Control と別型 (相互 From なし) |

## 3. 可換性は証明書 [N4]

`CertifiedCommutator` は (graded) bracket ノルムの**区間** [lo, hi] で持ち、閾値 τ
との 3 値裁定: hi < τ → Commuting / lo > τ → NonCommuting / **跨ぎ → Abstain**。
Abstain 対を含む文脈構成・非可換グラフの成分分解は**拒否**される (辺の強制禁止 —
HOLD-7 の棄却原則を可換構造の読み出しに継承)。実測: ‖[X₁,Z₁]‖_F = 2√8 =
5.656854 (明確な辺)・可換対 0・跨ぎ区間 [0.3τ, 1.5τ] → abstain。

## 4. grading は型 — Jordan–Wigner 弦の幾何誤読の遮断 [N5]

独立 3 フェルミオンモード (JW: γ₁ = XII, γ₃ = ZXI, γ₅ = ZZX) の site-local odd
演算子は、**ordinary 可換子では完全グラフ K₃ に見える** (‖[γ_a, γ_b]‖_F = 2√8 —
odd 対は反可換するから [γ,γ'] = 2γγ' ≠ 0)。真の構造は graded bracket (odd×odd は
反可換子) で読む: ‖{γ_a, γ_b}‖ = 0 (厳密) — **空グラフ = 真の独立**。

- `OperationalNet<OrdinaryCommutation>` は **odd primitive を構成時に拒否**する
  (型ゲート — 幾何の捏造が起こる前に遮断)。
- parity-even 双線形 (iγ₁γ₂ 等) は ordinary で厳密可換 (0.0) — **現行資格は
  parity-even lane** であり、graded lane は odd を扱う唯一の入口。

## 5. 型契約の登録 (v322 [N6][N7])

- `sim/src/operational_net.rs` — 役割 4 型 (相互 From なし)・OperatorParity・
  CommutationGrading 2 タグ (sealed)・CertifiedCommutator・OperationalNet
  (contexts は全対 Commuting 証明書を要求)・**FactorizationReading (v32.3 の出力型
  を先凍結**: Exact / SuperselectionSectors [中心非自明で tensor を強制しない] /
  EquivalenceClassOnly / Abstain 4 理由)・**InteractionHypergraph (v32.5 の出力型を
  先凍結**: w_S = ‖H_S‖_F²)・共用素子 (algebra_closure・commutant_dim)。
- core.schema.yml に概念 4 種 (OperationalNet / PrimitiveOperation /
  CertifiedCommutator / GlobalClosure) + **禁止変換 11** を登録。
- 旧 OperationalAlgebra (readout_contract) は互換のため凍結保存 + 後継注記。

## 6. 正直な残高

- 本版の系は toy (dim 8・可換子は exact ノルムの区間)。**ノイズ下の可換子測定・
  復元の gauge orbit 裁定・中心非自明系・graded recovery は v32.3** (目標定理 B)。
- 「操作可能な観測・相互作用が tensor product structure を規定する」という一般原理
  は既存研究にある (Zanardi らの operational tensor product / operational mereology
  の subalgebra–commutant 対)。**QRN の新規性は原理でなく、識別可能性の fiber・
  棄却型 (Abstain 一級市民)・証明書つき可換構造・holdout を含む構成的実装**に置く
  (PROMPT/13 §7)。
- bridge law 登録簿は全能力で空・external_replications = 0 のまま。
