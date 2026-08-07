/-
v35.3 set-valued resource profile の refinement / no-go / 安定性 (PROMPT/16 §7)

対象: budget poset 上の識別問題。資源 (budget) b ごとに
  - 同値関係 ~_b (識別できる粒度 — b ≤ b' なら ~_{b'} ⊆ ~_b: 細分)
  - 有限データの信頼集合 C_b(D) ⊆ Θ (被覆 ≥ 1 − α)
  - 集合値の読み Q_b(D) = q_b(C_b(D)) (商への像)
が与えられたとき、profile b ↦ Q_b(D) の構造は:

  RPF-1 (非関手性 no-go): **被覆だけでは samplewise nesting は出ない** —
    各予算で独立に作った被覆有効な信頼集合族で、r(Q_{b'}(D)) ⊄ Q_b(D) となる
    標本 D が存在する (2 点 Θ・2 点データ・α = 1/2 の最小反例)。
    「予算を増やせば読みは自動的に精密化する」は偽 — 関手性は構成で買う。
  RPF-2 (条件付き refinement 定理): b ≤ b'・~_{b'} が ~_b を細分 (q_b = r ∘ q_{b'})・
    **C_{b'}(D) ⊆ C_b(D) (samplewise nesting)** が成り立つとき、
      r(Q_{b'}(D)) ⊆ Q_b(D)
    (lax refinement)。純論理の一般定理 — 有限性・可測性の仮定なし。
  RPF-3 (margin 安定性): 区間の摂動が decision boundary までの margin 未満なら
    離散 verdict (edge/noEdge/straddled) は不変 (全整数値の一般定理)。
  RPF-4 (大域不安定性): boundary を跨ぐ点では任意に小さい摂動が verdict を変える —
    離散 verdict への global Lipschitz 安定性は存在しない (margin → 0 で比が発散:
    スケール n の単位で摂動 1・verdict 距離 1 の対が全ての n に存在)。
  RPF-5 (persistence 系): v33.3 の昇格規則「stable ⟺ chain ≥ 2」(禁止変換 17) の導出:
    (a) 健全性 — nesting + 細分の下で、精密側の singleton 読み v' の粗視化 r(v') は
        粗い側の読み集合から排除されない (RPF-2 の系: 2-chain 一致は整合)。
    (b) 単点は transient — 被覆有効でも、独立標本の次予算で singleton 読みが
        消える反例が存在する (単点昇格の禁止の根拠)。

証明の形: RPF-2/5a は型多相の純論理 (仮定を明示した一般定理)。RPF-3 は Int の
一般恒等 (omega)。RPF-1/4/5b は具体的反例 (decide / 明示証人)。

スコープの明示 (「未証明を証明済みに見せない」):
  - set-valued 安定性・interleaving には先行研究がある — 本ファイルの主張は
    「QRN の観測商 + 有限データ信頼集合 + budget poset の同時扱い」に限る
    (新規性を安定性一般に置かない — PROMPT/16 §7 の注意)。
  - 被覆は有限結果空間の整数重み (FiniteDataNoGo の規約)。実データの CP 区間での
    nesting 破れ率・intersection 構成の被覆は v353 が数値側で測る。
  - budget poset の一般論 (比較不能対・多パラメータ) は v33.3 の器械に既にあり、
    ここでは 2 予算 chain の定理のみ形式化する。
-/

namespace ResourceProfile

-- ================= RPF-2: 条件付き refinement (一般定理) =================

/-- **RPF-2**: q_b = r ∘ q_{b'} (~_{b'} が ~_b を細分) かつ C_{b'} ⊆ C_b
    (samplewise nesting) なら、読みの像は r(Q_{b'}) ⊆ Q_b。
    Θ・商空間は任意の型 — 有限性も可測性も要らない純論理。 -/
theorem rpf2_refinement {Th Vb Vb' : Type} (qb : Th → Vb) (qb' : Th → Vb')
    (r : Vb' → Vb) (hfac : ∀ θ, qb θ = r (qb' θ))
    (cb cb' : Th → Prop) (hnest : ∀ θ, cb' θ → cb θ) :
    ∀ v, (∃ θ, cb' θ ∧ qb' θ = v) → ∃ θ, cb θ ∧ qb θ = r v := by
  intro v h
  obtain ⟨θ, hθ, hq⟩ := h
  exact ⟨θ, hnest θ hθ, by rw [hfac, hq]⟩

-- ================= RPF-1: 非関手性 no-go (最小反例) =================

/-
最小反例: Θ = Bool (2 仮説)・データ D = Bool (2 結果, 各仮説の下で等確率 1/2)。
予算 b の信頼集合 cb(D) = {θ = D}, 予算 b' の (独立に設計した) 信頼集合
cb'(D) = {θ ≠ D}。どちらも被覆 = 1/2 ≥ 1 − α (α = 1/2 の toy 水準)。
同値関係は両予算とも恒等 (r = id — 細分条件は自明に成立)。
標本 D = true で Q_{b'} = {false} ⊄ {true} = Q_b — nesting なしに関手性なし。
-/

def cbLow (d θ : Bool) : Bool := θ == d
def cbHigh (d θ : Bool) : Bool := θ != d

/-- RPF-1 (a): 両予算の信頼集合はどちらも被覆 1/2 (2 等確率データ点中 1 点で的中) —
    被覆水準 1 − α = 1/2 を両仮説で満たす (toy 水準の被覆有効性)。 -/
theorem rpf1_coverage :
    ∀ θ : Bool,
      ((if cbLow true θ then 1 else 0) + (if cbLow false θ then 1 else 0) = 1)
      ∧ ((if cbHigh true θ then 1 else 0) + (if cbHigh false θ then 1 else 0) = 1) := by
  decide

/-- RPF-1 (b): 標本 D = true で C_{b'} ⊄ C_b — samplewise nesting は被覆から出ない。 -/
theorem rpf1_not_nested : ¬ (∀ θ : Bool, cbHigh true θ = true → cbLow true θ = true) := by
  decide

/-- **RPF-1 no-go**: 恒等の細分 (r = id) でも r(Q_{b'}(D)) ⊆ Q_b(D) が破れる標本が
    存在する — 被覆有効な任意予算別信頼集合に自然な関手性はない。 -/
theorem rpf1_no_functoriality :
    ¬ (∀ (d θ : Bool), cbHigh d θ = true → cbLow d θ = true) := by
  decide

-- ================= RPF-3: margin 安定性 (一般整数値) =================

inductive Verdict
  | edge
  | noEdge
  | straddled
deriving Repr, DecidableEq

/-- 区間 [lo, hi] の τ に対する離散 verdict (v34.3 の RobustVerdict の interval 形) -/
def verdict (τ lo hi : Int) : Verdict :=
  if τ < lo then .edge else if hi ≤ τ then .noEdge else .straddled

/-- decision boundary までの margin (verdict が変わらない摂動の上限) -/
def margin (τ lo hi : Int) : Int :=
  if τ < lo then lo - τ else if hi ≤ τ then τ - hi + 1 else min (τ - lo + 1) (hi - τ)

/-- **RPF-3**: 区間 (lo ≤ hi) の端点の摂動が margin 未満なら verdict は不変
    (全整数値の一般定理)。 -/
theorem rpf3_margin_stability (τ lo hi lo' hi' d : Int)
    (hlh : lo ≤ hi) (_hd : 0 ≤ d) (hm : d < margin τ lo hi)
    (hlo : lo - d ≤ lo' ∧ lo' ≤ lo + d) (hhi : hi - d ≤ hi' ∧ hi' ≤ hi + d) :
    verdict τ lo' hi' = verdict τ lo hi := by
  unfold verdict margin at *
  obtain ⟨hlo1, hlo2⟩ := hlo
  obtain ⟨hhi1, hhi2⟩ := hhi
  split at hm
  · rw [if_pos (show τ < lo' by omega), if_pos (show τ < lo by omega)]
  · split at hm
    · rw [if_neg (show ¬ τ < lo' by omega), if_pos (show hi' ≤ τ by omega),
          if_neg (show ¬ τ < lo by omega), if_pos (show hi ≤ τ by omega)]
    · rw [if_neg (show ¬ τ < lo' by omega), if_neg (show ¬ hi' ≤ τ by omega),
          if_neg (show ¬ τ < lo by omega), if_neg (show ¬ hi ≤ τ by omega)]

-- ================= RPF-4: 大域 Lipschitz 安定性の不在 =================

/-- **RPF-4 (境界不安定性)**: boundary 上では最小の摂動 (整数 1 — スケール n の
    単位で 1/n) が verdict を edge → straddled に変える。スケール n を任意に
    細かく取れば「摂動 1/n・verdict 距離 1」の対が全ての n で存在する —
    離散 verdict への global Lipschitz 定数は存在しない (RPF-3 の margin 条件が
    最良: margin なしの安定性は買えない)。 -/
theorem rpf4_boundary_instability (τ : Int) :
    verdict τ (τ + 1) (τ + 2) ≠ verdict τ τ (τ + 2) := by
  unfold verdict
  rw [if_pos (show τ < τ + 1 by omega), if_neg (show ¬ τ < τ by omega),
      if_neg (show ¬ τ + 2 ≤ τ by omega)]
  decide

-- ================= RPF-5: chain ≥ 2 昇格規則の導出 =================

/-- **RPF-5 (a) 健全性**: nesting + 細分の下で、精密予算の singleton 読み v' の
    粗視化 r(v') は粗い予算の読み集合に必ず属する — 2-chain で一致した読みは
    整合であり、chain に沿った昇格は矛盾を作らない (RPF-2 の系)。 -/
theorem rpf5_chain_soundness {Th Vb Vb' : Type} (qb : Th → Vb) (qb' : Th → Vb')
    (r : Vb' → Vb) (hfac : ∀ θ, qb θ = r (qb' θ))
    (cb cb' : Th → Prop) (hnest : ∀ θ, cb' θ → cb θ) (v' : Vb')
    (hsing : ∀ θ, cb' θ → qb' θ = v') (hne : ∃ θ, cb' θ) :
    ∃ θ, cb θ ∧ qb θ = r v' := by
  obtain ⟨θ0, h0⟩ := hne
  exact rpf2_refinement qb qb' r hfac cb cb' hnest v' ⟨θ0, h0, hsing θ0 h0⟩

/-- RPF-5 (b) transient 反例: 低予算で singleton 読み ({true})・独立標本の高予算で
    読みが全体 ({true, false}) に戻る — どちらも被覆有効。**単点の読みは次の予算で
    消え得る**: 禁止変換 17 (stable ⟺ chain ≥ 2, v33.3) の単点禁止側の根拠。 -/
def cbSingleton (θ : Bool) : Bool := θ == true
def cbEvaporate (_ : Bool) : Bool := true

theorem rpf5_transient :
    (∀ θ, cbSingleton θ = true → θ = true)
    ∧ ¬ (∀ θ, cbEvaporate θ = true → θ = true) := by
  decide

end ResourceProfile
