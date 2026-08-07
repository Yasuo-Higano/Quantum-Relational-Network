/-
v35.2 開放系 signed-response 観測商の形式化 (PROMPT/16 第三十五期 §5)

対象: 有限モード・number-conserving・quasi-free (Markov) 開放系の normal covariance
  Ċ = X C + C X† + Y   (X = −ih − ½(ℓℓ† + gg†) は有効 drift, Y = gg† は affine 注入)
に対する signed probe 対 C±(0) = C₀ ± εP_i の応答が識別するもの:

  GQF-1 (affine 消去): 符号つき差分 D = C⁺ − C⁻ から Y が**恒等的に**消え、
    D は同次方程式 Ḋ = XD + DX† に従う。核は flow 右辺の双線形性。
  GQF-2 (有効 drift 曲率): (n̈_j⁺ − n̈_j⁻)(0)/(4ε) = ‖P_j X P_i‖_F²。核は
    (i) 第二 jet の第 1・第 3 項の trace 消去 (P_i P_j = 0),
    (ii) block Frobenius 恒等式 Tr(B B†) = ‖B‖²_F (B = P_j X P_i)。
    閉鎖系 (X = −ih) では既存の ‖P_j h P_i‖² (LocalBias.lean) に戻るが、
    開放系で曲率が読むのは **Hamiltonian ではなく有効 drift** である。
  GQF-3 (Hamiltonian 昇格 no-go): coherent hopping (h ≠ 0, ℓ = 0) と collective
    loss (h = 0, ℓℓ† = 2(1,−1)(1,−1)ᵀ) が**同一の曲率統計**を返す厳密反例対 —
    CurvatureOnlyOpenResponse 契約から Hamiltonian support への写像は存在しない。
  GQF-4 (正側の門): dissipator が node-block 対角 (cross-node 散逸 drift = 0) なら
    有効 drift の off-diagonal block は −i h の block に厳密一致 — 曲率は
    ‖P_j h P_i‖² へ還元し、Hamiltonian 昇格が解錠される。
  GQF-5 (charge attribution no-go): 電荷応答 (登録統計 = flow の総数 trace) が
    loss (散逸的電荷破れ・pairing なし) と BdG pairing (coherent 電荷破れ・
    散逸なし) で**同一値 −8** を取る厳密反例対 — ChargeNonconservingResponse
    単独から HamiltonianPairing は選べない。正側: H_B が block 対角 (Δ = 0) なら
    閉鎖 Nambu の電荷統計は恒等的に 0 — 散逸ゼロ証明書の下で ≠ 0 ⇒ Δ ≠ 0。

証明の形 (本ファイルの規約 — LocalBias.lean の格子 4 点法からの強化):
  行列を成分構造体 (2×2 Gauss 整数 / その 2×2 block = 4×4) で表し、恒等式を
  **全ての整数値で成立する一般定理**として simp (積の分配・可換の正規化) +
  omega (原子化線形算術) で証明する。全整数点で成立する多項式恒等式は
  (ℤⁿ の Zariski 稠密性により) 係数の恒等式であり、ℝ/ℂ への拡張は標準論法 —
  格子有限窓より強い形で同じ結論を与える。反例対 (GQF-3/5) は具体的整数
  インスタンスの decide。

スコープの明示 (「未証明を証明済みに見せない」):
  - ODE の解の存在・e^{Xt} の解析・Taylor 展開の収束は形式化しない — ここで
    証明するのは微分が生む**行列恒等式** (jet 代数) の側。数値側 v352 が
    dense Lindblad (2^N 次元の量子発展) との一致で covariance 閉包を検査する。
  - 行列は 2×2 (1 次元 node 対) と 2×2 block (多軌道 node 対 = 4×4)。一般
    n×n・一般 block 分割は同じ展開 — 数値側が複素乱択で検査。
  - GQF-5 の登録統計のスケール規約 (R2 = 2R, QN = 2Q) はコメントに固定 —
    統計が物理の dN/dt に一致することは v352 が dense lane で較正する。
    状態は物理的 (loss: C = diag(1,0), Λ = ℓℓ† ⪰ 0 / pairing: pure BCS)。
  - ℤ[i] → ℂ の拡張 (稠密性論法) 自体は未形式化 (mathlib 級)。
-/

namespace OpenQuotient

-- ================= Gauss 整数 ℤ[i] =================

structure GI where
  re : Int
  im : Int
deriving Repr, DecidableEq

def gadd (x y : GI) : GI := ⟨x.re + y.re, x.im + y.im⟩
def gsub (x y : GI) : GI := ⟨x.re - y.re, x.im - y.im⟩
def gmul (x y : GI) : GI := ⟨x.re * y.re - x.im * y.im, x.re * y.im + x.im * y.re⟩
def gconj (x : GI) : GI := ⟨x.re, -x.im⟩
def gzero : GI := ⟨0, 0⟩
def gone : GI := ⟨1, 0⟩
/-- |z|² (整数値) -/
def gnorm2 (x : GI) : Int := x.re * x.re + x.im * x.im

-- ================= 2×2 複素行列 (成分構造体) =================

/-- 2×2 Gauss 整数行列 (a b / c d) -/
structure M2 where
  a : GI
  b : GI
  c : GI
  d : GI
deriving Repr, DecidableEq

def madd (x y : M2) : M2 := ⟨gadd x.a y.a, gadd x.b y.b, gadd x.c y.c, gadd x.d y.d⟩
def msub (x y : M2) : M2 := ⟨gsub x.a y.a, gsub x.b y.b, gsub x.c y.c, gsub x.d y.d⟩
def mmul (x y : M2) : M2 :=
  ⟨gadd (gmul x.a y.a) (gmul x.b y.c), gadd (gmul x.a y.b) (gmul x.b y.d),
   gadd (gmul x.c y.a) (gmul x.d y.c), gadd (gmul x.c y.b) (gmul x.d y.d)⟩
/-- 共役転置 X† -/
def adjM (x : M2) : M2 := ⟨gconj x.a, gconj x.c, gconj x.b, gconj x.d⟩
def traceM (x : M2) : GI := gadd x.a x.d
/-- Frobenius ノルム² (整数値) -/
def frob2 (x : M2) : Int := gnorm2 x.a + gnorm2 x.b + gnorm2 x.c + gnorm2 x.d
/-- Gauss 整数スカラー倍 -/
def msmul (s : GI) (x : M2) : M2 := ⟨gmul s x.a, gmul s x.b, gmul s x.c, gmul s x.d⟩
def mzero : M2 := ⟨gzero, gzero, gzero, gzero⟩
def mid : M2 := ⟨gone, gzero, gzero, gone⟩

/-- 開放系 flow の右辺 (affine): XC + CX† + Y -/
def rhsAff (x c y : M2) : M2 := madd (madd (mmul x c) (mmul c (adjM x))) y
/-- 同次 flow の右辺: XD + DX† -/
def rhsHom (x d : M2) : M2 := madd (mmul x d) (mmul d (adjM x))

/-- 1 次元 node 対の射影 (2×2): P₁ = diag(1,0), P₂ = diag(0,1) -/
def p1 : M2 := ⟨gone, gzero, gzero, gzero⟩
def p2 : M2 := ⟨gzero, gzero, gzero, gone⟩

/-- 第二 jet の差分核: X²P + 2XPX† + P(X†)² (D(t) = 2ε e^{Xt} P e^{X†t} の t² 係数 ×2/(2ε)) -/
def jet2 (x p : M2) : M2 :=
  madd (madd (mmul (mmul x x) p)
             (madd (mmul (mmul x p) (adjM x)) (mmul (mmul x p) (adjM x))))
       (mmul p (mmul (adjM x) (adjM x)))

attribute [local simp] gadd gsub gmul gconj gzero gone gnorm2
  madd msub mmul adjM traceM frob2 msmul mzero mid rhsAff rhsHom p1 p2 jet2

-- ================= GQF-1: affine 消去 (一般整数値の恒等式) =================

/-- **GQF-1**: 符号つき差分から affine 項 Y が恒等的に消える —
    rhsAff(X, C⁺, Y) − rhsAff(X, C⁻, Y) = rhsHom(X, C⁺ − C⁻)。
    全ての X, C⁺, C⁻, Y ∈ M₂(ℤ[i]) で成立 (係数の多項式恒等式 → ℂ へ稠密性で拡張)。 -/
theorem gqf1_affine_cancellation (x cp cm y : M2) :
    msub (rhsAff x cp y) (rhsAff x cm y) = rhsHom x (msub cp cm) := by
  simp [M2.mk.injEq, GI.mk.injEq, Int.mul_sub, Int.sub_mul]
  omega

/-- GQF-1 系: probe 対 C± = C₀ ± εP の差分は εP + εP (= 2εP) — 状態依存部 C₀ の
    厳密消去 (ε は整数スケール — 有理/実への拡張は稠密性)。 -/
theorem gqf1_probe_difference (c0 p : M2) (eps : Int) :
    msub (madd c0 (msmul ⟨eps, 0⟩ p)) (msub c0 (msmul ⟨eps, 0⟩ p))
      = madd (msmul ⟨eps, 0⟩ p) (msmul ⟨eps, 0⟩ p) := by
  simp [M2.mk.injEq, GI.mk.injEq]
  omega

-- ================= GQF-2: 有効 drift 曲率 (一般整数値の恒等式) =================

/-- **GQF-2 核 — block Frobenius 恒等式**: Tr(B B†) = ‖B‖²_F (実数値・虚部 0)。
    曲率の中央項 Tr(P_j X P_i X†) は B = P_j X P_i でこの形に落ちる。 -/
theorem gqf2_block_frobenius (b : M2) :
    traceM (mmul b (adjM b)) = ⟨frob2 b, 0⟩ := by
  simp [GI.mk.injEq, Int.mul_comm, Int.mul_neg]
  omega

/-- GQF-2 一階消滅: Δṅ_j(0) ∝ Tr(P₂(XP₁ + P₁X†)) = 0 (j ≠ i) —
    signed 差分の一階応答は恒等的に 0 で、曲率 (二階) が主要項。 -/
theorem gqf2_first_derivative_zero (x : M2) :
    traceM (mmul p2 (madd (mmul x p1) (mmul p1 (adjM x)))) = gzero := by
  simp

/-- **GQF-2 (2×2・1 次元 node)**: Tr(P₂ · jet2(X, P₁)) = 2‖P₂XP₁‖²_F —
    非エルミート X (開放系有効 drift) を含む全整数値で成立。 -/
theorem gqf2_curvature_2x2 (x : M2) :
    traceM (mmul p2 (jet2 x p1)) = ⟨2 * frob2 (mmul p2 (mmul x p1)), 0⟩ := by
  simp [GI.mk.injEq, Int.mul_comm, Int.mul_neg, Int.neg_neg]
  omega

-- ================= 2×2 block (= 4×4, 多軌道 node) =================

/-- 4×4 行列を 2×2 block (各 M2) で表す: (A B / C D) -/
structure M4 where
  A : M2
  B : M2
  C : M2
  D : M2
deriving Repr, DecidableEq

def madd4 (x y : M4) : M4 := ⟨madd x.A y.A, madd x.B y.B, madd x.C y.C, madd x.D y.D⟩
def msub4 (x y : M4) : M4 := ⟨msub x.A y.A, msub x.B y.B, msub x.C y.C, msub x.D y.D⟩
def mmul4 (x y : M4) : M4 :=
  ⟨madd (mmul x.A y.A) (mmul x.B y.C), madd (mmul x.A y.B) (mmul x.B y.D),
   madd (mmul x.C y.A) (mmul x.D y.C), madd (mmul x.C y.B) (mmul x.D y.D)⟩
def adjM4 (x : M4) : M4 := ⟨adjM x.A, adjM x.C, adjM x.B, adjM x.D⟩
def traceM4 (x : M4) : GI := gadd (traceM x.A) (traceM x.D)
def msmul4 (s : GI) (x : M4) : M4 := ⟨msmul s x.A, msmul s x.B, msmul s x.C, msmul s x.D⟩

/-- 2 次元 node 対の射影 (block): P₁ = (I 0 / 0 0), P₂ = (0 0 / 0 I) -/
def q1 : M4 := ⟨mid, mzero, mzero, mzero⟩
def q2 : M4 := ⟨mzero, mzero, mzero, mid⟩

def jet2b (x p : M4) : M4 :=
  madd4 (madd4 (mmul4 (mmul4 x x) p)
               (madd4 (mmul4 (mmul4 x p) (adjM4 x)) (mmul4 (mmul4 x p) (adjM4 x))))
        (mmul4 p (mmul4 (adjM4 x) (adjM4 x)))

attribute [local simp] madd4 msub4 mmul4 adjM4 traceM4 msmul4 q1 q2 jet2b

/-- **GQF-2 (4×4・2 軌道 node block 版)**: Tr(P₂ · jet2(X, P₁)) = 2‖X₂₁ block‖²_F —
    多軌道 node でも曲率は有効 drift の cross block ノルムを厳密に読む。 -/
theorem gqf2_curvature_block (x : M4) :
    traceM4 (mmul4 q2 (jet2b x q1)) = ⟨2 * frob2 x.C, 0⟩ := by
  simp [GI.mk.injEq, Int.mul_comm, Int.mul_neg, Int.neg_neg]
  omega

-- ================= GQF-3: Hamiltonian 昇格 no-go (厳密反例対) =================

/-- coherent hopping: h = (0 1 / 1 0), 散逸なし → X_A = −ih = (0 −i / −i 0) -/
def XA : M2 := ⟨gzero, ⟨0, -1⟩, ⟨0, -1⟩, gzero⟩
/-- collective loss: h = 0, Λ = ℓℓ† = 2(1,−1)(1,−1)ᵀ → X_B = −½Λ = (−1 1 / 1 −1) -/
def XB : M2 := ⟨⟨-1, 0⟩, ⟨1, 0⟩, ⟨1, 0⟩, ⟨-1, 0⟩⟩

/-- 登録された曲率統計 (CurvatureOnlyOpenResponse 契約): w₂₁ = ‖P₂XP₁‖²_F -/
def curvStat (x : M2) : Int := frob2 (mmul p2 (mmul x p1))
/-- Hamiltonian 部の off-diagonal (× 2 スケール): 2h₂₁ = i(X₂₁ − conj X₁₂) -/
def hamOff2 (x : M2) : GI := gmul ⟨0, 1⟩ (gsub x.c (gconj x.b))

/-- GQF-3 (a): 両モデルの曲率統計は同一 (= 1) -/
theorem gqf3_same_curvature : curvStat XA = curvStat XB ∧ curvStat XA = 1 := by
  decide

/-- GQF-3 (b): Hamiltonian off-diagonal は異なる (coherent: 2h₂₁ = 2 / loss: 0) -/
theorem gqf3_different_hamiltonian : hamOff2 XA = ⟨2, 0⟩ ∧ hamOff2 XB = gzero := by
  decide

/-- **GQF-3 no-go**: 曲率統計から Hamiltonian off-diagonal への写像は存在しない —
    任意の f について、f が XA で正しければ XB で必ず誤る (同一入力・異なる正解)。 -/
theorem gqf3_no_map (f : Int → GI) :
    ¬(f (curvStat XA) = hamOff2 XA ∧ f (curvStat XB) = hamOff2 XB) := by
  intro ⟨hA, hB⟩
  have hc : curvStat XA = curvStat XB := gqf3_same_curvature.1
  rw [hc] at hA
  rw [hA] at hB
  have h2 : hamOff2 XA = ⟨2, 0⟩ := gqf3_different_hamiltonian.1
  have h0 : hamOff2 XB = gzero := gqf3_different_hamiltonian.2
  rw [h2, h0] at hB
  exact absurd hB (by decide)

-- ================= GQF-4: 正側の門 (block 対角 dissipator での還元) =================

/-- **GQF-4 (i)**: dissipator Γ が node-block 対角 (cross-node 散逸 drift = 0) なら、
    スケール有効 drift 2X = −2ih − Γ の cross block は −2i h₂₁ に厳密一致する —
    Γ は cross block に一切触れない (OffDiagonalDissipatorZero 証明書の中身)。 -/
theorem gqf4_block_diagonal_reduction (h : M4) (g1 g2 : M2) :
    (msub4 (msmul4 ⟨0, -2⟩ h) ⟨g1, mzero, mzero, g2⟩).C = msmul ⟨0, -2⟩ h.C := by
  simp

/-- **GQF-4 (ii)**: ‖−i·B‖²_F = ‖B‖²_F — 単位スカラー −i は Frobenius 等距
    (「−ih か h か」はノルムに影響しない)。 -/
theorem gqf4_neg_i_isometry (b : M2) :
    frob2 (msmul ⟨0, -1⟩ b) = frob2 b := by
  simp [Int.mul_neg, Int.neg_mul, Int.neg_neg]
  omega

/-- GQF-4 補題: ‖B + B‖²_F = 4‖B‖²_F (整数スケール 2 の帳尻 — 2B = B + B) -/
theorem gqf4_scale_two (b : M2) :
    frob2 (madd b b) = 4 * frob2 b := by
  simp [Int.add_mul, Int.mul_add]
  omega

/-- **GQF-4 (iii)**: ‖−2i·B‖²_F = 4‖B‖²_F — (i) の block 等式と合わせて、block
    対角 dissipator の下で曲率統計は ‖P_j h P_i‖² へ厳密還元する
    (OffDiagonalDissipatorZero 証明書が Hamiltonian 昇格を解錠する門)。 -/
theorem gqf4_scale_norm (b : M2) :
    frob2 (msmul ⟨0, -2⟩ b) = 4 * frob2 b := by
  have hcomp : msmul ⟨0, -2⟩ b = msmul ⟨0, -1⟩ (madd b b) := by
    simp [M2.mk.injEq, GI.mk.injEq]
    omega
  rw [hcomp, gqf4_neg_i_isometry, gqf4_scale_two]

-- ================= GQF-5: charge attribution no-go (厳密反例対) =================

/-
登録統計 (スケール規約を固定):
  loss lane (number-conserving covariance): S = Re Tr(XC + CX†) (Y = 0) —
    d(Tr C)/dt の flow 表現。モデル: Λ = ℓℓ† = diag(8,0) (ℓ = 2√2 e₁ — Λ は整数
    PSD), C = diag(1,0) (物理的: 0 ≤ C ≤ I)。X = −½Λ = diag(−4,0)。
  pairing lane (閉鎖 Nambu): S = Re Tr(QN · (−i)[H_B, R2]) — スケール規約
    R2 = 2R (R = Nambu covariance), QN = 2Q (Q = number 観測量) を宣言。
    モデル: h = 0, Δ = (0 1 / −1 0), 状態 = pure BCS (C = I/2, |κ| = ½ — 物理的)。
  この二つの登録統計が物理の dN/dt に一致することは v352 が dense lane で較正する。
  ここで証明するのは: 両モデルが**同一の統計値 −8** を返し、pairing の有無が
  異なること — 値から pairing への写像の不在。
-/

/-- loss モデルの電荷応答統計: Re Tr(XC + CX†), X = diag(−4,0), C = diag(1,0) -/
def chargeStatLoss : Int :=
  (traceM (rhsHom ⟨⟨-4, 0⟩, gzero, gzero, gzero⟩ p1)).re

/-- BdG Hamiltonian (h = 0, Δ = (0 1 / −1 0)): H_B = (0 Δ / Δ† 0) -/
def hBdG : M4 :=
  ⟨mzero, ⟨gzero, gone, ⟨-1, 0⟩, gzero⟩, ⟨gzero, ⟨-1, 0⟩, gone, gzero⟩, mzero⟩
/-- pure BCS 状態 (スケール R2 = 2R): 対角 block I, 異常 block 2κ = (0 i / −i 0) -/
def rBCS : M4 :=
  ⟨mid, ⟨gzero, ⟨0, 1⟩, ⟨0, -1⟩, gzero⟩, ⟨gzero, ⟨0, 1⟩, ⟨0, -1⟩, gzero⟩, mid⟩
/-- number 観測量 (スケール QN = 2Q): diag(1,1,−1,−1) -/
def qN : M4 := ⟨mid, mzero, mzero, msmul ⟨-1, 0⟩ mid⟩

attribute [local simp] hBdG rBCS qN

/-- pairing モデルの電荷応答統計: Re Tr(QN · (−i)[H_B, R2]) -/
def chargeStatPairing : Int :=
  (traceM4 (mmul4 qN (msmul4 ⟨0, -1⟩
    (msub4 (mmul4 hBdG rBCS) (mmul4 rBCS hBdG))))).re

/-- GQF-5 (a): 両モデルの電荷応答統計は同一 (= −8)。
    loss は pairing なし (Δ = 0)・pairing は散逸なし (Λ = 0) — 電荷非保存の
    出所が散逸か Hamiltonian pairing かは統計値から区別できない。 -/
theorem gqf5_same_charge_response :
    chargeStatLoss = -8 ∧ chargeStatPairing = -8 := by
  decide

/-- **GQF-5 no-go**: 電荷応答値から「pairing の有無」への写像は存在しない —
    loss モデル (pairing = false) と pairing モデル (pairing = true) が同一値。 -/
theorem gqf5_no_map (f : Int → Bool) :
    ¬(f chargeStatLoss = false ∧ f chargeStatPairing = true) := by
  intro ⟨hA, hB⟩
  have h1 : chargeStatLoss = -8 := gqf5_same_charge_response.1
  have h2 : chargeStatPairing = -8 := gqf5_same_charge_response.2
  rw [h1] at hA
  rw [h2] at hB
  rw [hA] at hB
  exact Bool.false_ne_true hB

/-- **GQF-5 正側 (解錠条件)**: H_B が block 対角 (Δ = 0 — pairing なしの閉鎖
    Nambu) なら電荷統計は**全ての状態 R で恒等的に 0** — よって散逸ゼロ証明書の
    下で統計 ≠ 0 は Δ ≠ 0 (pairing witness) を含意する。 -/
theorem gqf5_closed_charge_conservation (h1 h2 : M2) (r : M4) :
    (traceM4 (mmul4 qN (msmul4 ⟨0, -1⟩
      (msub4 (mmul4 ⟨h1, mzero, mzero, h2⟩ r) (mmul4 r ⟨h1, mzero, mzero, h2⟩))))).re
      = 0 := by
  simp [Int.mul_comm, Int.mul_neg, Int.neg_neg]
  omega

end OpenQuotient
