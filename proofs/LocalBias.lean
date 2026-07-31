/-
v31.2 LocalBiasCommutatorLaw の形式化 (PROMPT/12 第三十一期)

対象: 既知ノード因子分解の射影 P_i に対する probe 対 C_i^± = I/2 ± εP_i の応答則:
  (1) **commutator block 恒等式**: j ≠ i (P_j P_i = 0) で
        P_j (Ċ⁺ − Ċ⁻)(0) P_i = −2iε P_j h P_i   (Ċ = −i[h, C])
      核となるのは P_j (h P_i − P_i h) P_i = P_j h P_i (射影の直交性が第 2 項を消す)。
  (2) **密度曲率 Frobenius 恒等式**: n_j = Tr(P_j C) に対し
        (n̈_j⁺ − n̈_j⁻)(0) / (4ε) = ‖P_j h P_i‖_F²
      核となるのは −Tr(P_j [h,[h,P_i]]) = 2 Tr(P_j h P_i h) = 2‖P_j h P_i‖_F² (j ≠ i)。
  (3) **probe 差の分離**: (I/2 + εP) − (I/2 − εP) = 2εP — 状態非依存部が厳密に消える。
  (4) **局所ゲージ不変性**: ‖U B‖_F = ‖B‖_F (G(U) = I) — 読み出しがノード内基底に
      依存しないこと (GibbsInverse.lean の Frobenius–Gram 恒等式と同じ機構 —
      本ファイルでも自己完結に再掲する)。

証明の形 (Projector.lean / GibbsInverse.lean の規約): 全て成分の整数多項式恒等式
(per-variable 次数 ≤ 2)。3 点/変数で確定するところを格子 {−1, 0, 1, 2} (4 点) の全数
native_decide で検証。ℤ → ℝ/ℂ への拡張は多項式恒等式の標準論法 (未形式化 — Lean が
保証するのは格子上の恒等)。

スコープの明示 (「未証明を証明済みに見せない」):
  - (1)(2) は実対称 h・1 次元ノード (3 サイト: P₁ = diag(1,0,0), P₂ = diag(0,1,0))
    および 2 軌道ノード (4 サイト: P₁ = diag(1,1,0,0), P₂ = diag(0,0,1,1)) の格子恒等。
    複素エルミート h・一般 d 軌道は未形式化 — 数値側 v312_localbias_law が複素乱択で
    検査する。
  - 時間微分そのもの (Ċ = −i[h,C] の導出・Taylor 展開) は形式化しない — ここで証明
    するのは微分が生む**行列恒等式**の側。
-/

namespace LocalBias

-- ---------------- 3 サイト (1 次元ノード × 2 + 傍観サイト 1) ----------------

/-- 3×3 実対称行列 (h11 h12 h13 / h12 h22 h23 / h13 h23 h33) — 6 自由度 -/
structure S3 where
  a : Int -- h11
  b : Int -- h12
  c : Int -- h13
  d : Int -- h22
  e : Int -- h23
  f : Int -- h33
deriving Repr, DecidableEq

def ent (h : S3) : Nat → Nat → Int
  | 0, 0 => h.a
  | 0, 1 => h.b
  | 0, 2 => h.c
  | 1, 0 => h.b
  | 1, 1 => h.d
  | 1, 2 => h.e
  | 2, 0 => h.c
  | 2, 1 => h.e
  | _, _ => h.f

/-- ノード射影 (1 次元): P₁ = e₀e₀ᵀ, P₂ = e₁e₁ᵀ -/
def p1 (i j : Nat) : Int := if i == 0 && j == 0 then 1 else 0
def p2 (i j : Nat) : Int := if i == 1 && j == 1 then 1 else 0

def sum3 (g : Nat → Int) : Int := g 0 + g 1 + g 2

/-- (P₂ (h P₁ − P₁ h) P₁)_{ij} -/
def commBlock (h : S3) (i j : Nat) : Int :=
  sum3 fun x => sum3 fun y => sum3 fun z =>
    p2 i x * (ent h x y * p1 y z - p1 x y * ent h y z) * p1 z j

/-- (P₂ h P₁)_{ij} -/
def hBlock (h : S3) (i j : Nat) : Int :=
  sum3 fun x => sum3 fun y => p2 i x * ent h x y * p1 y j

def grid : List Int := [-1, 0, 1, 2]

def gridS3 : List S3 :=
  grid.flatMap fun a => grid.flatMap fun b => grid.flatMap fun c =>
    grid.flatMap fun d => grid.flatMap fun e => grid.map fun f => ⟨a, b, c, d, e, f⟩

def idx : List (Nat × Nat) :=
  [0, 1, 2].flatMap fun i => [0, 1, 2].map fun j => (i, j)

/-- 定理 1 (commutator block 恒等式): P₂ (hP₁ − P₁h) P₁ = P₂ h P₁ —
    P₂P₁ = 0 が第 2 項を消す。probe 差の commutator から h block が裸で出る機構 -/
theorem commutator_block_identity :
    (gridS3.all fun h => idx.all fun ij =>
      commBlock h ij.1 ij.2 == hBlock h ij.1 ij.2) = true := by native_decide

/-- Tr(P₂ h P₁ h) -/
def trPhPh (h : S3) : Int :=
  sum3 fun i => sum3 fun j => sum3 fun k => sum3 fun l =>
    p2 i j * ent h j k * p1 k l * ent h l i

/-- ‖P₂ h P₁‖_F² = Σ (P₂hP₁)_{ij}² -/
def frobBlock (h : S3) : Int :=
  sum3 fun i => sum3 fun j => hBlock h i j * hBlock h i j

/-- 定理 2 (trace–Frobenius 恒等式): Tr(P₂ h P₁ h) = ‖P₂ h P₁‖_F² (実対称 h) -/
theorem trace_frobenius_identity :
    (gridS3.all fun h => trPhPh h == frobBlock h) = true := by native_decide

/-- −Tr(P₂ [h,[h,P₁]]) の整数化: [h,P₁] = hP₁ − P₁h, [h,[h,P₁]] = h(hP₁−P₁h) − (hP₁−P₁h)h -/
def negTrDoubleComm (h : S3) : Int :=
  -(sum3 fun i => sum3 fun j => sum3 fun k => sum3 fun l =>
    p2 i j *
      (ent h j k * (sum3 fun m => ent h k m * p1 m l - p1 k m * ent h m l)
        - (sum3 fun m => ent h j m * p1 m k - p1 j m * ent h m k) * ent h k l)
      * (if l == i then 1 else 0))

/-- 定理 3 (密度曲率恒等式の核): −Tr(P₂[h,[h,P₁]]) = 2 Tr(P₂hP₁h) = 2‖P₂hP₁‖_F²
    (j ≠ i)。probe 差 ×ε と合わせ (n̈⁺−n̈⁻)/(4ε) = ‖P₂hP₁‖_F² を与える -/
theorem curvature_frobenius_identity :
    (gridS3.all fun h => negTrDoubleComm h == 2 * frobBlock h) = true := by native_decide

/-- 定理 4 (probe 差の分離): (I/2 + εP) − (I/2 − εP) = 2εP — 成分恒等式
    (I/2 は整数化のため 2 倍して I とする: (I + 2εP) − (I − 2εP) = 4εP) -/
theorem probe_difference_isolation :
    (grid.all fun eps => idx.all fun ij =>
      let plus := (if ij.1 == ij.2 then 1 else 0) + 2 * eps * p1 ij.1 ij.2
      let minus := (if ij.1 == ij.2 then 1 else 0) - 2 * eps * p1 ij.1 ij.2
      plus - minus == 4 * eps * p1 ij.1 ij.2) = true := by native_decide

-- ---------------- 4 サイト (2 軌道ノード × 2) ----------------

/-- 4×4 実対称 h = [[D₁, B], [Bᵀ, D₂]] — D₁ (3 自由度) + B (4) + D₂ (3) = 10 自由度 -/
structure S4 where
  d1a : Int
  d1b : Int
  d1d : Int
  b11 : Int
  b12 : Int
  b21 : Int
  b22 : Int
  d2a : Int
  d2b : Int
  d2d : Int
deriving Repr, DecidableEq

def ent4 (h : S4) : Nat → Nat → Int
  | 0, 0 => h.d1a
  | 0, 1 => h.d1b
  | 1, 0 => h.d1b
  | 1, 1 => h.d1d
  | 0, 2 => h.b11
  | 0, 3 => h.b12
  | 1, 2 => h.b21
  | 1, 3 => h.b22
  | 2, 0 => h.b11
  | 3, 0 => h.b12
  | 2, 1 => h.b21
  | 3, 1 => h.b22
  | 2, 2 => h.d2a
  | 2, 3 => h.d2b
  | 3, 2 => h.d2b
  | _, _ => h.d2d

/-- 2 軌道ノード射影: Q₁ = diag(1,1,0,0), Q₂ = diag(0,0,1,1) -/
def q1 (i j : Nat) : Int := if i == j && i < 2 then 1 else 0
def q2 (i j : Nat) : Int := if i == j && 2 ≤ i && i < 4 then 1 else 0

def sum4 (g : Nat → Int) : Int := g 0 + g 1 + g 2 + g 3

def trQhQh (h : S4) : Int :=
  sum4 fun i => sum4 fun j => sum4 fun k => sum4 fun l =>
    q2 i j * ent4 h j k * q1 k l * ent4 h l i

/-- ‖Q₂ h Q₁‖_F² = Σ B 成分² (= b11² + b12² + b21² + b22²) -/
def frobB (h : S4) : Int :=
  h.b11 * h.b11 + h.b12 * h.b12 + h.b21 * h.b21 + h.b22 * h.b22

def gridS4 : List S4 :=
  -- 10 変数格子は 4¹⁰ = 1M 点で重い — 対角 6 変数は {0,1}, block 4 変数は 4 点
  -- (検査対象の恒等式は block 変数に次数 2・対角変数に次数 ≤ 1 で、2 点/1 次数で十分)
  [0, 1].flatMap fun a => [0, 1].flatMap fun b => [0, 1].flatMap fun c =>
    grid.flatMap fun x11 => grid.flatMap fun x12 =>
      grid.flatMap fun x21 => grid.flatMap fun x22 =>
        [0, 1].flatMap fun d => [0, 1].flatMap fun e => [0, 1].map fun f =>
          (⟨a, b, c, x11, x12, x21, x22, d, e, f⟩ : S4)

/-- 定理 5 (多軌道 trace–Frobenius): Tr(Q₂ h Q₁ h) = ‖B‖_F² — 2 軌道ノードでも
    密度曲率が inter-node block の Frobenius 重みを裸で返す (対角 block は寄与しない)。
    対角変数は次数 ≤ 1 なので 2 点格子で確定 (コメントの次数勘定) -/
theorem multiorbital_trace_frobenius :
    (gridS4.all fun h => trQhQh h == frobB h) = true := by native_decide

-- ---------------- 局所ゲージ不変性 (自己完結の再掲) ----------------

structure M2 where
  a : Int
  b : Int
  c : Int
  d : Int
deriving Repr, DecidableEq

def mul2 (x y : M2) : M2 :=
  ⟨x.a * y.a + x.b * y.c, x.a * y.b + x.b * y.d,
   x.c * y.a + x.d * y.c, x.c * y.b + x.d * y.d⟩

def tp2 (x : M2) : M2 := ⟨x.a, x.c, x.b, x.d⟩

def frob2 (x : M2) : Int := x.a * x.a + x.b * x.b + x.c * x.c + x.d * x.d

def gram2 (u : M2) : M2 := mul2 (tp2 u) u

def gridM2 : List M2 :=
  grid.flatMap fun a => grid.flatMap fun b => grid.flatMap fun c => grid.map fun d =>
    (⟨a, b, c, d⟩ : M2)

/-- 定理 6 (読み出しの局所ゲージ不変性): G(U) = G(V) = I ⇒ ‖U B Vᵀ‖_F² = ‖B‖_F² —
    密度曲率読み出し ‖P_j h P_i‖_F² はノード内基底の取り方に依らない
    (GibbsInverse.frob_orthogonal_both と同じ機構 — 自己完結の格子版) -/
theorem gauge_invariance_frobenius :
    (gridM2.all fun u => gridM2.all fun b => gridM2.all fun v =>
      !(gram2 u == (⟨1, 0, 0, 1⟩ : M2)) || !(gram2 v == (⟨1, 0, 0, 1⟩ : M2))
        || (frob2 (mul2 u (mul2 b (tp2 v))) == frob2 b)) = true := by native_decide

end LocalBias
