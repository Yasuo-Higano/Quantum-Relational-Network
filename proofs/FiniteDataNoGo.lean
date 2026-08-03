/-
v34.3 有限データ昇格不能定理 (第四の no-go) と Robust Promotion Theorem (PROMPT/15 §4)

対象: 有限結果空間上の二仮説判別と、同時信頼集合による昇格の誤り上限。

  (1) **Le Cam 二点下限 (第四の no-go の有限版)** — 結果空間の各点 x に
      仮説 0/1 の確率重み (w0 x, w1 x) を置く (共通分母でスケール済みの整数重み —
      有理確率は常にこの形に書ける)。任意の決定規則 d : x ↦ {0,1} の誤り質量
        err(d) = Σ_{x : d x = 1} w0 x + Σ_{x : d x = 0} w1 x
      は Σ_x min(w0 x, w1 x) を下回れない。さらに
        2 Σ min = W0 + W1 − Σ|w0 − w1|
      なので、W0 = W1 = W (正規化) では平均誤り err/(2W) ≥ (1 − TV)/2、
      TV = Σ|w0−w1|/(2W)。観測契約が二つの interface を区別しない (P0 = P1)
      なら平均誤り ≥ 1/2 — **強制回答器は原理的に当てられない**。
      尤度比規則 (Bayes) が下限を達成する (下限は最良)。

  (2) **Robust Promotion Theorem (正側)** — 信頼集合 C(データ) が真値 θ を
      確率 ≥ 1 − α で含み (被覆)、昇格規則が「C 上で読みが一定のときだけ回答」
      なら、誤った昇格の確率 ≤ α。核は排除補題: 回答が r(θ) と異なるなら
      θ は C の外 (C 内で一定な値は θ ∈ C なら r(θ) に一致するから)。
      よって {誤昇格} ⊆ {被覆失敗} で質量単調性から上限が従う。
      **これは無条件の誤昇格確率であり、回答条件付き risk ≤ α ではない**
      (selective risk ≤ α/P(answer) — 分母は別途要る。数値側 v343 [F3b] が
      条件付き risk > α の実例を厳密計算で与える)。

証明の形: 整数重みのリスト帰納 + omega (min/natAbs は omega が扱う)。
実数係数への拡張は共通分母スケーリングの標準論法 (有理確率で完全一般)。
決定規則の全数実例 (N = 4 Bernoulli 対, 32 規則) は decide で機械確認。

スコープの明示 (PROMPT/12「未証明を証明済みに見せない」の継承):
  - 確率は有限結果空間上の有理数 (整数重み) — 連続結果空間・一般測度は
    論文側の一般版 (v34.3 文書) に委ね、ここでは有限版のみを機械保証する。
  - ランダム化決定規則は扱わない (決定的規則の凸結合なので下限は同じ —
    論法はコメントのみ、形式化は決定的規則)。
  - 被覆仮説そのもの (Clopper–Pearson 等の構成) は形式化しない — ここで
    証明するのは「被覆 ⇒ 誤昇格 ≤ α」の含意。構成の厳密性は v343 が
    二項分布の厳密和で機械検査する。
-/

namespace FiniteDataNoGo

/-- 結果空間の 1 点: (仮説 0 の重み, 仮説 1 の重み, 規則の回答)。
    回答 d = true は「仮説 1」、false は「仮説 0」。 -/
abbrev Row := Int × Int × Bool

/-- 誤り質量: 仮説 0 の下で true と答えた重み + 仮説 1 の下で false と答えた重み -/
def errMass : List Row → Int
  | [] => 0
  | (w0, w1, d) :: rest => (if d then w0 else w1) + errMass rest

/-- 点ごとの最小重みの和 (達成可能な最小誤り質量) -/
def minSum : List Row → Int
  | [] => 0
  | (w0, w1, _) :: rest => min w0 w1 + minSum rest

/-- 仮説 0 の総重み -/
def wSum0 : List Row → Int
  | [] => 0
  | (w0, _, _) :: rest => w0 + wSum0 rest

/-- 仮説 1 の総重み -/
def wSum1 : List Row → Int
  | [] => 0
  | (_, w1, _) :: rest => w1 + wSum1 rest

/-- 全変動距離の 2W 倍: Σ|w0 − w1| -/
def absDiffSum : List Row → Int
  | [] => 0
  | (w0, w1, _) :: rest => (w0 - w1).natAbs + absDiffSum rest

/-- [1] どんな決定規則も点ごと最小の和を下回れない (重みの符号仮定は不要) -/
theorem errMass_ge_minSum (rows : List Row) : minSum rows ≤ errMass rows := by
  induction rows with
  | nil => simp [errMass, minSum]
  | cons r rest ih =>
    obtain ⟨w0, w1, d⟩ := r
    simp only [errMass, minSum]
    cases d <;> simp <;> omega

/-- [2] 恒等式 2 Σ min = W0 + W1 − Σ|w0 − w1| -/
theorem two_minSum_eq (rows : List Row) :
    2 * minSum rows = wSum0 rows + wSum1 rows - absDiffSum rows := by
  induction rows with
  | nil => simp [minSum, wSum0, wSum1, absDiffSum]
  | cons r rest ih =>
    obtain ⟨w0, w1, d⟩ := r
    simp only [minSum, wSum0, wSum1, absDiffSum]
    omega

/-- [3] Le Cam 二点下限 (第四の no-go の有限版):
    2 err ≥ W0 + W1 − Σ|w0 − w1|。正規化 W0 = W1 = W で
    平均誤り err/(2W) ≥ (1 − TV)/2。 -/
theorem le_cam_two_point (rows : List Row) :
    wSum0 rows + wSum1 rows - absDiffSum rows ≤ 2 * errMass rows := by
  have h1 := errMass_ge_minSum rows
  have h2 := two_minSum_eq rows
  omega

/-- [4] 識別不能なら当てられない: 観測契約が二仮説を区別しない (w0 = w1 点ごと)
    とき、どんな規則の誤り質量も片側の総重み W を下回れない (平均誤り ≥ 1/2)。 -/
theorem indistinguishable_half (rows : List Row)
    (h : ∀ r ∈ rows, r.1 = r.2.1) : wSum0 rows ≤ errMass rows := by
  induction rows with
  | nil => simp [errMass, wSum0]
  | cons r rest ih =>
    obtain ⟨w0, w1, d⟩ := r
    have hr : w0 = w1 := h (w0, w1, d) List.mem_cons_self
    have hrest : ∀ r ∈ rest, r.1 = r.2.1 := fun r hm => h r (List.mem_cons_of_mem _ hm)
    have hmono := ih hrest
    simp only [errMass, wSum0]
    cases d <;> simp <;> omega

/-- 尤度比 (Bayes) 規則: w1 が真に大きい点でのみ「仮説 1」と答える -/
def bayesRows (pairs : List (Int × Int)) : List Row :=
  pairs.map fun p => (p.1, p.2, decide (p.1 < p.2))

/-- [5] Bayes 規則は下限を達成する — Le Cam 下限は最良 (改良不能) -/
theorem bayes_achieves (pairs : List (Int × Int)) :
    errMass (bayesRows pairs) = minSum (bayesRows pairs) := by
  induction pairs with
  | nil => rfl
  | cons p rest ih =>
    obtain ⟨w0, w1⟩ := p
    simp only [bayesRows, List.map_cons] at ih ⊢
    simp only [errMass, minSum]
    rcases Decidable.em (w0 < w1) with hlt | hge
    · rw [decide_eq_true hlt]
      simp
      omega
    · rw [decide_eq_false hge]
      simp
      omega

/-- [6] 排除補題 (Robust Promotion の核): 読み r が信頼集合 C 上で一定値 v を
    取り、真値 θ の読みが v と異なる (= 誤った昇格) なら、θ は C の外にある。 -/
theorem promotion_exclusion {Θ V : Type} (C : Θ → Prop) (r : Θ → V) (v : V)
    (hconst : ∀ θ', C θ' → r θ' = v) (θ : Θ) (hwrong : r θ ≠ v) : ¬C θ :=
  fun hC => hwrong (hconst θ hC)

/-- 被覆と誤昇格の質量: 行 = (重み, θ ∈ C(データ) か, 誤昇格か) -/
abbrev CovRow := Int × Bool × Bool

/-- 誤昇格の質量 -/
def wrongMass : List CovRow → Int
  | [] => 0
  | (w, _, wrong) :: rest => (if wrong then w else 0) + wrongMass rest

/-- 被覆失敗 (θ ∉ C) の質量 -/
def missMass : List CovRow → Int
  | [] => 0
  | (w, inC, _) :: rest => (if inC then 0 else w) + missMass rest

/-- [7] 質量単調性: 各データ点で「誤昇格 → 被覆失敗」(排除補題) が成り立ち
    重みが非負なら、誤昇格質量 ≤ 被覆失敗質量。 -/
theorem wrongMass_le_missMass (rows : List CovRow)
    (hpos : ∀ r ∈ rows, 0 ≤ r.1)
    (hexcl : ∀ r ∈ rows, r.2.2 = true → r.2.1 = false) :
    wrongMass rows ≤ missMass rows := by
  induction rows with
  | nil => simp [wrongMass, missMass]
  | cons r rest ih =>
    obtain ⟨w, inC, wrong⟩ := r
    have hw : 0 ≤ w := hpos (w, inC, wrong) List.mem_cons_self
    have hrest₁ : ∀ r ∈ rest, 0 ≤ r.1 := fun r hm => hpos r (List.mem_cons_of_mem _ hm)
    have hrest₂ : ∀ r ∈ rest, r.2.2 = true → r.2.1 = false := fun r hm =>
      hexcl r (List.mem_cons_of_mem _ hm)
    have hmono := ih hrest₁ hrest₂
    simp only [wrongMass, missMass]
    cases hwr : wrong
    · cases inC <;> simp <;> omega
    · have hin : inC = false := hexcl (w, inC, wrong) List.mem_cons_self hwr
      subst hin
      simp
      omega

/-- [8] Robust Promotion Theorem (有限版): 被覆失敗質量 ≤ A (被覆 ≥ 1 − α の
    スケール形) なら誤昇格質量 ≤ A。selective risk (回答条件付き) の上限では
    ないことに注意 — それは別途 P(answer) の下限を要する。 -/
theorem robust_promotion (rows : List CovRow) (A : Int)
    (hpos : ∀ r ∈ rows, 0 ≤ r.1)
    (hexcl : ∀ r ∈ rows, r.2.2 = true → r.2.1 = false)
    (hcov : missMass rows ≤ A) : wrongMass rows ≤ A :=
  Int.le_trans (wrongMass_le_missMass rows hpos hexcl) hcov

/-
具体実例 (N = 4 Bernoulli 対, θ0 = 1/4 vs θ1 = 1/2, 共通スケール 256):
  w0 = (3/4 + 1/4)^4 の二項重み × 256 = [81, 108, 54, 12, 1]   (k = 0..4)
  w1 = (1/2 + 1/2)^4 の二項重み × 256 = [16, 64, 96, 64, 16]
  Σ min = 16 + 64 + 54 + 12 + 1 = 147, Σ|diff| = 65+44+42+52+15 = 218,
  検算: 2·147 = 294 = 256 + 256 − 218 ✓ (最小平均誤り 147/512, TV = 218/512)
-/

/-- N = 4 実例の行を規則 (d0..d4) から組む -/
def rows4 (d0 d1 d2 d3 d4 : Bool) : List Row :=
  [(81, 16, d0), (108, 64, d1), (54, 96, d2), (12, 64, d3), (1, 16, d4)]

/-- [9] 実例の下限値: Σ min = 147 (2·147 = 512 − 218 の検算込み) -/
theorem instance_n4_min :
    minSum (rows4 false false true true true) = 147 ∧
    wSum0 (rows4 false false true true true) = 256 ∧
    wSum1 (rows4 false false true true true) = 256 ∧
    absDiffSum (rows4 false false true true true) = 218 := by decide

/-- [10] 全 32 決定規則の機械全数: どの規則も誤り質量 147 を下回れない
    (Le Cam 下限のこの実例での最良性 — Bayes 規則 (F,F,T,T,T) が 147 を達成) -/
theorem instance_n4_all_rules :
    ∀ d0 d1 d2 d3 d4 : Bool, 147 ≤ errMass (rows4 d0 d1 d2 d3 d4) := by decide

end FiniteDataNoGo
