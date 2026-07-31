/-
v31.1 GaussianGibbsInverseOracle の同値類構造の形式化 (PROMPT/12 第三十一期)

対象: 熱的 Gaussian の大域逆変換 K(C) = log[(I−C)C⁻¹] = β(h−μI) が定める同値類:
  (1) **block-local unitary ノルム不変性** — K の (i,j) block B に対し
      ‖U B Vᵀ‖_F = ‖B‖_F (U, V 直交)。核は Frobenius–Gram 恒等式
        ‖U B‖² = ‖B‖² + (G₁₁−1)·(1行寄与) + (G₂₂−1)·(2行寄与) + (G₁₂+G₂₁)·(交差項)
      (G = UᵀU) — 直交性 G = I が余剰項を消す。
  (2) **β 未知の正スケール同値** — 符号 (支持) とエントリ交差比の不変性
  (3) **μ 未知の恒等シフト無関係性** — off-diagonal は μ に依存しない

証明の形 (Projector.lean の規約): 恒等式は成分の整数多項式で per-variable 次数 ≤ 2。
よって 4 点/変数の格子 {−1,0,1,2} 上の全数検証 (native_decide) が多項式恒等を確定し
(3 点で十分のところ 4 点)、ℤ → ℝ への拡張は多項式恒等式の標準論法による。
この完了論法自体は未形式化 (mathlib 級) — **Lean が保証するのは格子上の恒等**。

スコープの明示 (PROMPT/12「未証明を証明済みに見せない」):
  - block は 2×2 (最小の多軌道ノード)。一般 d×d は同じ Gram 展開の帰納 — 未形式化
    (数値側 v311_gibbs_oracle [T4] が d=2 ノードの実乱択で検査する)。
  - matrix logarithm そのもの (K(C) の解析的構成・固有値分解の存在) は形式化しない —
    ここで証明するのは「K = β(h−μI) 型の行列」の同値類構造。
  - 直交行列の条件付き定理 (frob_orthogonal_*) も格子上 — 格子内の直交行列は符号付き
    置換 8 個で、実回転群への拡張は恒等式定理 + G = I 代入 (コメントの論法) による。
  - 正スケール同値の実数版 (β > 0) は初等順序論 — ここでは整数窓 β ∈ [1,8] の全数。
-/

namespace GibbsInverse

/-- 2×2 整数行列 (a b / c d) — 最小の多軌道ノード block -/
structure M2 where
  a : Int
  b : Int
  c : Int
  d : Int
deriving Repr, DecidableEq

def mul (x y : M2) : M2 :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d,
   x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩

/-- 転置 -/
def tp (x : M2) : M2 := ⟨x.a, x.c, x.b, x.d⟩

/-- Frobenius ノルムの 2 乗 -/
def frob2 (x : M2) : Int := x.a * x.a + x.b * x.b + x.c * x.c + x.d * x.d

/-- trace -/
def trM (x : M2) : Int := x.a + x.d

/-- Gram 行列 G = UᵀU -/
def gram (u : M2) : M2 := mul (tp u) u

def idm : M2 := ⟨1, 0, 0, 1⟩

/-- スカラー倍 (β·H) -/
def scaleM (k : Int) (x : M2) : M2 := ⟨k * x.a, k * x.b, k * x.c, k * x.d⟩

/-- 対角シフト (H + c·I) — 化学ポテンシャル μ の入り方 -/
def addDiag (x : M2) (c : Int) : M2 := ⟨x.a + c, x.b, x.c, x.d + c⟩

-- ---------------- 格子 ----------------

/-- 4 点/変数の格子 (per-variable 次数 ≤ 2 の恒等式には 3 点で十分 — 余裕 1 点) -/
def grid : List Int := [-1, 0, 1, 2]

def gridM2 : List M2 :=
  grid.flatMap fun a => grid.flatMap fun b =>
    grid.flatMap fun c => grid.map fun d => ⟨a, b, c, d⟩

def onPairs (p : M2 → M2 → Bool) : Bool :=
  gridM2.all fun u => gridM2.all fun b => p u b

def onTriples (p : M2 → M2 → M2 → Bool) : Bool :=
  gridM2.all fun u => gridM2.all fun b => gridM2.all fun v => p u b v

-- ---------------- (1) Frobenius–Gram 恒等式と直交不変性 ----------------

/-- 定理 1 (左乗恒等式): ‖UB‖² = ‖B‖² + Gram 偏差項。
    per-variable 次数 ≤ 2 の 8 変数多項式恒等式 — 格子 4⁸ = 65536 対で確定 -/
theorem frob_left_gram_identity :
    onPairs (fun u b =>
      frob2 (mul u b)
        == frob2 b
           + ((gram u).a - 1) * (b.a * b.a + b.b * b.b)
           + ((gram u).d - 1) * (b.c * b.c + b.d * b.d)
           + ((gram u).b + (gram u).c) * (b.a * b.c + b.b * b.d)) = true := by
  native_decide

/-- 定理 2 (右乗恒等式): ‖BVᵀ‖² = ‖B‖² + Gram 偏差項 (行側) -/
theorem frob_right_gram_identity :
    onPairs (fun v b =>
      frob2 (mul b (tp v))
        == frob2 b
           + ((gram v).a - 1) * (b.a * b.a + b.c * b.c)
           + ((gram v).d - 1) * (b.b * b.b + b.d * b.d)
           + ((gram v).b + (gram v).c) * (b.a * b.b + b.c * b.d)) = true := by
  native_decide

/-- 定理 3 (直交不変性, 左): G(U) = I ⇒ ‖UB‖² = ‖B‖² (格子上 — 実回転群へは
    定理 1 + G = I 代入の論法で拡張。格子内の直交行列は符号付き置換) -/
theorem frob_orthogonal_left :
    onPairs (fun u b =>
      !(gram u == idm) || (frob2 (mul u b) == frob2 b)) = true := by
  native_decide

/-- 定理 4 (直交不変性, 両側): G(U) = G(V) = I ⇒ ‖U B Vᵀ‖² = ‖B‖² —
    ノード内基底変換 K_ij ↦ U_i K_ij U_j† の Frobenius 不変性の格子版 -/
theorem frob_orthogonal_both :
    onTriples (fun u b v =>
      !(gram u == idm) || !(gram v == idm)
        || (frob2 (mul u (mul b (tp v))) == frob2 b)) = true := by
  native_decide

/-- 定理 5 (trace 相似不変性): tr(U B Uᵀ) = tr(B G(U)) — 直交なら tr(B)。
    Frobenius 以外の相似不変量 (固有値対称関数) の最低次の代表 -/
theorem trace_conj_gram_identity :
    onPairs (fun u b =>
      trM (mul u (mul b (tp u))) == trM (mul b (gram u))) = true := by
  native_decide

-- ---------------- (2) β 未知 — 正スケール同値 ----------------

/-- 定理 6 (符号窓): β ∈ [1,8] で sign(β·x) = sign(x) (x ∈ [−8,8]) —
    正スケールは K の支持 (隣接) と符号構造を変えない。実数 β > 0 版は初等順序論 -/
theorem scale_sign_window :
    ((List.range 8).all fun bm1 =>
      ((List.range 17).all fun xp8 =>
        let beta : Int := Int.ofNat bm1 + 1
        let x : Int := Int.ofNat xp8 - 8
        (beta * x).sign == x.sign)) = true := by
  native_decide

/-- 定理 7 (交差比不変): (βH) のエントリ交差積は H のそれと一致 —
    「正の大域スケールを除いて確定」の等式的内容 (a·b′ = b·a′ 型)。全対で成立 -/
theorem scale_cross_ratio (k : Int) (h : M2) :
    (scaleM k h).a * h.b = (scaleM k h).b * h.a
      ∧ (scaleM k h).c * h.d = (scaleM k h).d * h.c
      ∧ (scaleM k h).a * h.d = (scaleM k h).d * h.a := by
  refine ⟨?_, ?_, ?_⟩ <;> simp only [scaleM] <;> ac_rfl

/-- 定理 8 (零支持窓): β ∈ [1,8] で (β·x = 0 ⇔ x = 0) — スケール同値類の中で
    「結合が存在するか」は不変 (隣接の読み出しは β を知らずに可能) -/
theorem scale_support_zero_window :
    ((List.range 8).all fun bm1 =>
      ((List.range 17).all fun xp8 =>
        let beta : Int := Int.ofNat bm1 + 1
        let x : Int := Int.ofNat xp8 - 8
        ((beta * x) == 0) == (x == 0))) = true := by
  native_decide

-- ---------------- (3) μ 未知 — 恒等シフト無関係性 ----------------

/-- 定理 9 (μ の off-diagonal 無関係性): K = β(H − μI) の off-diagonal は μ に
    依存しない — 空間隣接の読み出しに μ の知識は不要 (定義的に成立する型の事実) -/
theorem mu_shift_offdiag (beta mu : Int) (h : M2) :
    (scaleM beta (addDiag h (-mu))).b = beta * h.b
      ∧ (scaleM beta (addDiag h (-mu))).c = beta * h.c := by
  exact ⟨rfl, rfl⟩

end GibbsInverse
