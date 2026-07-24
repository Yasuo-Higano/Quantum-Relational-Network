/-
v26.6 Barnes–Rivers 型 projector 代数の機械検証 (PROMPT/8 — 重力真空偏極監査)

対象: 3 次元静的計量摂動 h_ij(q) の対称 2 階テンソル空間 (6 次元) 上の
spin 射影演算子 (Barnes–Rivers の静的・空間版):
  θ_ij = δ_ij − q_i q_j / q²   (横射影),  ω_ij = q_i q_j / q²   (縦射影)
  P2   = ½(θ⊗θ + θ⊗θ) − ½ θθ   (spin-2, TT — 3D の横平面は 2 次元なので 1/(d−1)=1/2)
  P1   = ½(θ⊗ω 対称 4 項)        (spin-1, 縦横混合)
  P0s  = ½ θ_ij θ_kl              (spin-0, 横トレース)
  P0w  = ω_ij ω_kl                (spin-0, 縦縦)

**整数化**: 分母を避けるため θ' := q²δ − q q, ω' := q q (多項式) を使い、
全射影演算子を m := 2 q⁴ P で表す (m2, m1, m0s, m0w は Int 値の多項式)。
  P∘P = P      ⇔  m∘m = 2q⁴ m
  P_a∘P_b = 0  ⇔  m_a∘m_b = 0
  ΣP = 1_sym   ⇔  m2+m1+m0s+m0w = q⁴ (δ_ik δ_jl + δ_il δ_jk)
  tr P = rank  ⇔  tr m = 2q⁴ rank

**証明の形**: 各恒等式は q ∈ ℤ³ の多項式恒等式で、成分ごとの per-variable 次数は
≤ 8 (合成 = 4次×4次)。したがって格子 {0,…,8}³ (9 点/変数) 上の全数検証が
多項式の恒等的零を含意する (1 変数ずつ「次数 d の多項式が d+1 点で零なら零」を
適用する標準論法)。この完了論法自体は mathlib 級のため未形式化 — Lean が保証する
のは格子上の恒等 (native_decide)、実数全域への拡張は上記次数勘定による
(AnomalyUpstream の「窓 native_decide + 限界の明示」と同じ規約)。

**ŷ 方向インスタンス** (q = (0,1,0)): v26.6 の数値器械 (v266_vacuum_pol) が使う
チャネル辞書を定理化する:
  D = E_xx − E_zz (plus)      → spin-2 の固有ベクトル (固有値 1)
  X = E_xz + E_zx (cross)     → spin-2 の固有ベクトル
  S = E_xx + E_zz (横トレース) → P0s の固有ベクトル
  L = E_yy (縦)               → P0w の固有ベクトル = 純ゲージ (L ∝ q⊗ŷ + ŷ⊗q)
  E_xy+E_yx, E_yz+E_zy        → P1 の固有ベクトル
ゲージモード h = q⊗ξ + ξ⊗q は P2 と P0s に消され (gauge_annihilated)、
P1 + P0w が再現する (gauge_reproduced) — 「spin-2 と横トレースが gauge 不変
チャネル、縦 (yy@qŷ) が純ゲージチャネル」の代数的根拠。

実行: cd proofs && ~/.elan/bin/lean Projector.lean  (終了コード 0 = 全定理検証)
-/

-- ---------------- 添字と格子 ----------------

abbrev Q3 := Int × Int × Int

def qc (q : Q3) : Nat → Int
  | 0 => q.1
  | 1 => q.2.1
  | _ => q.2.2

def q2 (q : Q3) : Int := qc q 0 * qc q 0 + qc q 1 * qc q 1 + qc q 2 * qc q 2

def del (i j : Nat) : Int := if i == j then 1 else 0

/-- 横射影の整数化 θ'_ij = q²δ_ij − q_i q_j -/
def th (q : Q3) (i j : Nat) : Int := q2 q * del i j - qc q i * qc q j

/-- 縦射影の整数化 ω'_ij = q_i q_j -/
def om (q : Q3) (i j : Nat) : Int := qc q i * qc q j

-- ---------------- 射影演算子 (m = 2q⁴ P の整数化) ----------------

def m2 (q : Q3) (i j k l : Nat) : Int :=
  th q i k * th q j l + th q i l * th q j k - th q i j * th q k l

def m1 (q : Q3) (i j k l : Nat) : Int :=
  th q i k * om q j l + th q i l * om q j k + om q i k * th q j l + om q i l * th q j k

def m0s (q : Q3) (i j k l : Nat) : Int := th q i j * th q k l

def m0w (q : Q3) (i j k l : Nat) : Int := 2 * om q i j * om q k l

/-- 対称テンソル空間上の恒等演算子の 2 倍: i2 = δ_ik δ_jl + δ_il δ_jk -/
def i2 (i j k l : Nat) : Int := del i k * del j l + del i l * del j k

def sum3 (f : Nat → Int) : Int := f 0 + f 1 + f 2

/-- 演算子合成 (A∘B)_{ij,kl} = Σ_{ab} A_{ij,ab} B_{ab,kl} -/
def comp (a b : Q3 → Nat → Nat → Nat → Nat → Int) (q : Q3) (i j k l : Nat) : Int :=
  sum3 fun x => sum3 fun y => a q i j x y * b q x y k l

/-- trace: tr m = Σ_{ij} m_{ij,ij} -/
def trm (a : Q3 → Nat → Nat → Nat → Nat → Int) (q : Q3) : Int :=
  sum3 fun i => sum3 fun j => a q i j i j

def idx3 : List Nat := [0, 1, 2]

def idx4 : List (Nat × Nat × Nat × Nat) :=
  idx3.flatMap fun i => idx3.flatMap fun j => idx3.flatMap fun k => idx3.map fun l => (i, j, k, l)

/-- 格子 {0,…,8}³ — per-variable 次数 ≤ 8 の多項式恒等式を確定する 9 点格子 -/
def grid9 : List Q3 :=
  (List.range 9).flatMap fun a =>
    (List.range 9).flatMap fun b =>
      (List.range 9).map fun c => ((Int.ofNat a, Int.ofNat b, Int.ofNat c) : Q3)

def onGrid (p : Q3 → Bool) : Bool := grid9.all p

def allIdx (p : Nat → Nat → Nat → Nat → Bool) : Bool :=
  idx4.all fun t => p t.1 t.2.1 t.2.2.1 t.2.2.2

-- ---------------- 定理 1: 完全性 ΣP = 1_sym ----------------

theorem completeness :
    onGrid (fun q => allIdx fun i j k l =>
      m2 q i j k l + m1 q i j k l + m0s q i j k l + m0w q i j k l
        == q2 q * q2 q * i2 i j k l) = true := by native_decide

-- ---------------- 定理 2–5: 冪等性 P∘P = P (m∘m = 2q⁴ m) ----------------

theorem idem_m2 :
    onGrid (fun q => allIdx fun i j k l =>
      comp m2 m2 q i j k l == 2 * q2 q * q2 q * m2 q i j k l) = true := by native_decide

theorem idem_m1 :
    onGrid (fun q => allIdx fun i j k l =>
      comp m1 m1 q i j k l == 2 * q2 q * q2 q * m1 q i j k l) = true := by native_decide

theorem idem_m0s :
    onGrid (fun q => allIdx fun i j k l =>
      comp m0s m0s q i j k l == 2 * q2 q * q2 q * m0s q i j k l) = true := by native_decide

theorem idem_m0w :
    onGrid (fun q => allIdx fun i j k l =>
      comp m0w m0w q i j k l == 2 * q2 q * q2 q * m0w q i j k l) = true := by native_decide

-- ---------------- 定理 6: 直交性 P_a∘P_b = 0 (a ≠ b, 両順序) ----------------

theorem orthogonality :
    onGrid (fun q => allIdx fun i j k l =>
      comp m2 m1 q i j k l == 0 && comp m1 m2 q i j k l == 0 &&
      comp m2 m0s q i j k l == 0 && comp m0s m2 q i j k l == 0 &&
      comp m2 m0w q i j k l == 0 && comp m0w m2 q i j k l == 0 &&
      comp m1 m0s q i j k l == 0 && comp m0s m1 q i j k l == 0 &&
      comp m1 m0w q i j k l == 0 && comp m0w m1 q i j k l == 0 &&
      comp m0s m0w q i j k l == 0 && comp m0w m0s q i j k l == 0) = true := by native_decide

-- ---------------- 定理 7: 対称性 (自己共役性) m_{ij,kl} = m_{ji,kl} = m_{ij,lk} = m_{kl,ij} ----------------

theorem symmetry :
    onGrid (fun q => allIdx fun i j k l =>
      m2 q i j k l == m2 q j i k l && m2 q i j k l == m2 q i j l k &&
      m2 q i j k l == m2 q k l i j &&
      m1 q i j k l == m1 q j i k l && m1 q i j k l == m1 q i j l k &&
      m1 q i j k l == m1 q k l i j &&
      m0s q i j k l == m0s q j i k l && m0s q i j k l == m0s q i j l k &&
      m0s q i j k l == m0s q k l i j &&
      m0w q i j k l == m0w q j i k l && m0w q i j k l == m0w q i j l k &&
      m0w q i j k l == m0w q k l i j) = true := by native_decide

-- ---------------- 定理 8: trace = rank (tr m = 2q⁴ rank; rank 2,2,1,1 — 計 6) ----------------

theorem trace_ranks :
    onGrid (fun q =>
      trm m2 q == 2 * q2 q * q2 q * 2 && trm m1 q == 2 * q2 q * q2 q * 2 &&
      trm m0s q == 2 * q2 q * q2 q * 1 && trm m0w q == 2 * q2 q * q2 q * 1) = true := by
  native_decide

-- ---------------- 定理 9: Ward 収縮 q_i (P2)_{ij,kl} = 0, q_i (P0s)_{ij,kl} = 0 ----------------

theorem ward_transverse :
    onGrid (fun q => allIdx fun _ j k l =>
      sum3 (fun i => qc q i * m2 q i j k l) == 0 &&
      sum3 (fun i => qc q i * m0s q i j k l) == 0) = true := by native_decide

-- ---------------- ゲージモード h_kl = q_k ξ_l + q_l ξ_k ----------------

def hg (q xi : Q3) (k l : Nat) : Int := qc q k * qc xi l + qc q l * qc xi k

def apply4 (a : Q3 → Nat → Nat → Nat → Nat → Int) (q : Q3) (v : Nat → Nat → Int)
    (i j : Nat) : Int :=
  sum3 fun k => sum3 fun l => a q i j k l * v k l

def xis : List Q3 := [(0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (1, 0, 1), (0, 1, 1), (1, 1, 1)]

-- ---------------- 定理 10: ゲージ消去 P2 h_gauge = 0, P0s h_gauge = 0 ----------------

theorem gauge_annihilated :
    onGrid (fun q => xis.all fun xi =>
      idx3.all fun i => idx3.all fun j =>
        apply4 m2 q (hg q xi) i j == 0 && apply4 m0s q (hg q xi) i j == 0) = true := by
  native_decide

-- ---------------- 定理 11: ゲージ再現 (P1 + P0w) h_gauge = h_gauge ----------------

theorem gauge_reproduced :
    onGrid (fun q => xis.all fun xi =>
      idx3.all fun i => idx3.all fun j =>
        apply4 m1 q (hg q xi) i j + apply4 m0w q (hg q xi) i j
          == 2 * q2 q * q2 q * hg q xi i j) = true := by native_decide

-- ---------------- ŷ 方向インスタンス (q = (0,1,0), q² = 1, 2q⁴ = 2) ----------------
-- v266_vacuum_pol.rs のチャネル辞書はこの節の定理に従う。

def hy : Q3 := (0, 1, 0)

def eD (i j : Nat) : Int := del i 0 * del j 0 - del i 2 * del j 2  -- plus:  E_xx − E_zz
def eX (i j : Nat) : Int := del i 0 * del j 2 + del i 2 * del j 0  -- cross: E_xz + E_zx
def eS (i j : Nat) : Int := del i 0 * del j 0 + del i 2 * del j 2  -- 横トレース: E_xx + E_zz
def eL (i j : Nat) : Int := del i 1 * del j 1                      -- 縦: E_yy
def eG1 (i j : Nat) : Int := del i 0 * del j 1 + del i 1 * del j 0 -- E_xy + E_yx
def eG2 (i j : Nat) : Int := del i 1 * del j 2 + del i 2 * del j 1 -- E_yz + E_zy

/-- チャネル辞書: D, X は spin-2 / S は P0s / L は P0w / G1, G2 は P1 の固有ベクトル
    (固有値 1 ⇔ apply4 m = 2·v)、かつ他の 3 射影に消される。 -/
theorem yhat_dictionary :
    (idx3.all fun i => idx3.all fun j =>
      apply4 m2 hy eD i j == 2 * eD i j && apply4 m1 hy eD i j == 0 &&
        apply4 m0s hy eD i j == 0 && apply4 m0w hy eD i j == 0 &&
      apply4 m2 hy eX i j == 2 * eX i j && apply4 m1 hy eX i j == 0 &&
        apply4 m0s hy eX i j == 0 && apply4 m0w hy eX i j == 0 &&
      apply4 m0s hy eS i j == 2 * eS i j && apply4 m2 hy eS i j == 0 &&
        apply4 m1 hy eS i j == 0 && apply4 m0w hy eS i j == 0 &&
      apply4 m0w hy eL i j == 2 * eL i j && apply4 m2 hy eL i j == 0 &&
        apply4 m1 hy eL i j == 0 && apply4 m0s hy eL i j == 0 &&
      apply4 m1 hy eG1 i j == 2 * eG1 i j && apply4 m2 hy eG1 i j == 0 &&
      apply4 m1 hy eG2 i j == 2 * eG2 i j && apply4 m2 hy eG2 i j == 0) = true := by
  native_decide

/-- 縦チャネルは純ゲージ: 2·E_yy = q⊗ŷ + ŷ⊗q (ξ = ŷ)。
    v26.6 の縦 Ward 監査 (「連続極限では K の縦列が消えるべき」) の代数的根拠。 -/
theorem yhat_longitudinal_is_gauge :
    (idx3.all fun k => idx3.all fun l => hg hy (0, 1, 0) k l == 2 * eL k l) = true := by
  native_decide
