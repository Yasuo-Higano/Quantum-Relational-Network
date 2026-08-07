# QRN v35.3 — set-valued resource profile の数学的地位: RPF 定理群

**位置づけ**: PROMPT/16 (第三十五期「外部化と観測商」) §7・優先度 3。v34.5 で
器械化した set-valued profile (budget ↦ 信頼集合上の読み) と v33.3 で凍結した
昇格規則 (禁止変換 17: stable ⟺ chain ≥ 2) の数学的地位を確定する。期テーゼの
「追加資源はその商を精密化し」の部分を、**精密化は自動ではなく構成で買う**と
正確化した版。

**一次ソース**: `results/v353_profile_refinement.txt` (8 検査 PASS)・
`results/v353_lean.txt` (`proofs/ResourceProfile.lean` 8 定理・終了コード 0)。

---

## 1. 問題の設定

budget poset B 上の識別問題。各 b ∈ B に:

- 同値関係 ~_b (その予算で識別できる粒度)。b ≤ b' で ~_{b'} ⊆ ~_b (細分) なら
  商写像の因子化 q_b = r_{b'b} ∘ q_{b'} がある。
- 有限データの信頼集合 C_b(D) ⊆ Θ (被覆 ≥ 1 − α)。
- 集合値の読み Q_b(D) = q_b(C_b(D))。

素朴な期待:「資源を増やせば読みは自動的に精密化する — profile は functor」。

## 2. RPF-1: 非関手性 no-go — 被覆から nesting は出ない

**Lean 最小反例** (`rpf1_coverage`/`rpf1_not_nested`/`rpf1_no_functoriality`):
Θ = {θ₁, θ₂}・データ 2 値・α = 1/2 の toy 水準で、cbLow(D) = {θ = D} と
cbHigh(D) = {θ ≠ D} はどちらも被覆 1/2 を満たすが、標本 D = true で
Q_{b'} = {θ₂} ⊄ {θ₁} = Q_b — 恒等の細分 (r = id) でも naturality が破れる。

**数値実測** ([R1], CP 区間・p = 0.42・n 100 → 300・4000 標本):

| 構成 | nesting 破れ率 |
|---|---|
| 独立標本 (予算ごとに別取得) | **0.437** |
| 累積標本 (restriction — 同じ master data の先頭) | **0.266** |
| intersection 構成 (α/2 同時配分 + 逐次交差) | **0** (被覆 0.970 ≥ 0.95) |

累積標本 (PROMPT/16 §3 の「common master data の restriction」) **単独では
不十分**であることが分かった — CP 区間は k の揺らぎで広がり方が変わるため、
n を増やした区間が前の区間からはみ出すことが 4 回に 1 回起こる。samplewise
nesting は **intersection (または anytime-valid 系) の構成で買う**しかない。

## 3. RPF-2: 条件付き refinement 定理 (正側)

`rpf2_refinement` (型多相の純論理 — 有限性・可測性の仮定なし):

> q_b = r ∘ q_{b'} (細分) かつ ∀θ (C_{b'}(D) θ → C_b(D) θ) (samplewise nesting)
> ならば r(Q_{b'}(D)) ⊆ Q_b(D)。

数値側 [R2]: intersection 構成 + verdict 語彙の商で違反 0/4000 (読みが実際に
狭まった標本 647 — 精密化は起きた上で単調)。**profile が lax refinement に
なるのは nesting を構成したときであり、かつそのときに限る** (RPF-1 が逆側)。

## 4. RPF-3/RPF-4: 安定性は margin 条件つき

- **RPF-3** (`rpf3_margin_stability`, 全整数値の一般定理): 区間 (lo ≤ hi) の
  端点摂動 d が margin(τ, lo, hi) 未満なら離散 verdict (edge/noEdge/straddled)
  は不変。margin は verdict 側ごとに lo−τ / τ−hi+1 / min(τ−lo+1, hi−τ)。
  数値 [R3]: margin 未満のランダム摂動 4000 件で flip 0・margin 超の flip 証人。
- **RPF-4** (`rpf4_boundary_instability`): boundary 上では最小摂動が verdict を
  変える — 摂動 1/n・verdict 距離 1 の対が全ての n に存在し ([R4]: 比 10 → 10⁵)、
  **離散 verdict への global Lipschitz 安定性は存在しない**。RPF-3 の margin
  条件が最良である (margin なしの安定性は買えない)。Straddled が一級市民で
  あることの定理側の根拠。

## 5. RPF-5: chain ≥ 2 昇格規則の導出 (v33.3 の宿題)

v33.3 は「stable ⟺ chain ≥ 2」(禁止変換 17) を**凍結された規則**として導入した。
本版でこれが refinement 構造の系になる:

- **健全性** (`rpf5_chain_soundness` = RPF-2 の系): nesting + 細分の下で、精密
  予算の singleton 読み v' の粗視化 r(v') は粗い予算の読み集合から排除されない —
  2-chain で一致した読みの昇格は矛盾を作らない。数値 [R5]: nested 2-chain の
  健全性違反 0。
- **単点は transient** (`rpf5_transient` + [R5] 実測): 被覆有効でも、独立標本の
  次予算で singleton 読みは消え得る — boundary 近傍 (p = 0.36, τ = 0.3) では
  singleton 371 件中 **273 件 (74%) が消失/変化**した。単点昇格の禁止は
  この実測率が物語る。

## 6. 検査一覧 (v353_profile_refinement — 8 PASS)

[R1a/b/c] nesting 破れ率 (独立/累積/intersection) / [R2] lax refinement 違反 0 /
[R3] margin 安定性 / [R4] Lipschitz 不在の発散表 / [R5] transient + 2-chain 健全性 /
[R6] Lean 反例の整数橋。

## 7. 限界と非主張

- set-valued 安定性・interleaving 自体には先行研究がある — 本版の主張は
  **観測商 + 有限データ信頼集合 + budget poset の同時扱い**に限る (安定性一般の
  新規性を主張しない — PROMPT/16 §7 の注意の遵守)。
- Lean の被覆は有限結果空間の整数重み (toy 水準 α = 1/2 の最小反例)。実 CP 区間
  の破れ率は数値実測 (シード固定・4000 標本)。
- 2 予算 chain の定理まで — budget poset 全体の多パラメータ persistence・
  比較不能対の一般論は v33.3 の器械の領分 (語彙の凍結は継続)。
- anytime-valid 信頼列・全予算同時 α 配分など intersection 以外の nesting 構成は
  設計選択肢として記録するに留める (OCS-2.0 の材料)。
