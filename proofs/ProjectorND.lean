/-
v26.8-0 抽象 Barnes–Rivers 代数 (PROMPT/9 §3.1 — spec §12.2-1)

Projector.lean (v26.6) は具体的な q ∈ {0..8}³ の格子 + 次数勘定による certificate
だった。本ファイルはこれを卒業し、**公理**

  θ² = θ,  ω² = ω,  θω = ωθ = 0,  θ + ω = I,  tr θ = τ (= d−1),  tr ω = 1
  (θ, ω 対称 — ω = q q^T/q² は rank 1。q² ≠ 0 は ω の存在が構造的に担う)

から Barnes–Rivers 射影の代数を**任意次元**で証明する。方法:

1. 対称 2 階テンソル空間上の演算子で BR 射影に現れる語は、rank-1 ω の下で
   6 個で閉じる (E7 = ω⊗ω-sym は E4 に退化):
     E1 = ½(θ⊗θ + θ⊗θ)   E2 = θ_{μν}θ_{ρσ}   E3 = ½(θω 対称 4 項)
     E4 = ω_{μν}ω_{ρσ}    E5 = θ_{μν}ω_{ρσ}   E6 = ω_{μν}θ_{ρσ}
2. 合成 (A∘B)_{μν,ρσ} = A_{μν,αβ}B_{αβ,ρσ} の乗積表は縮約計算 (各積に閉じた
   θ ループは高々 1 個 → 係数は τ の**高々 1 次**) で:
        |  E1    E2    E3   E4   E5    E6
     E1 |  E1    E2    0    0    E5    0
     E2 |  E2    τE2   0    0    τE5   0
     E3 |  0     0     E3   0    0     0
     E4 |  0     0     0    E4   0     E6
     E5 |  0     0     0    E5   0     E2
     E6 |  E6    τE6   0    0    τE4   0
3. この表を **d = 3 (τ = 2) と d = 4 (τ = 4−1 = 3) の具体行列インスタンスで
   全数認証**する (table_cert_d3 / table_cert_d4)。係数は τ に高々 1 次なので
   2 点 {2, 3} での一致が任意 τ での表を確定する (アフィン補間 — 上記 2 の
   縮約構造による次数上界が根拠)。
4. 以降の全定理は表の上の **ℤ[τ] 係数の形式代数の多項式恒等式** — 有限窓では
   なく任意次元 (任意 τ)、任意 q (q² ≠ 0) の定理になる。Projector.lean は
   「具体的な θ(q), ω(q) が公理を満たす」side の certificate として併存する。

射影 (P₂ = E1 − E2/τ, P₁ = E3, P₀s = E2/τ, P₀w = E4, 転送 P₀sw = E5/√τ,
P₀ws = E6/√τ) は分母回避のためスケール形で定理化する:
  m2 := τE1 − E2 (= τP₂): m2∘m2 = τ·m2
  m0s := E2: m0s∘m0s = τ·m0s / m1 := E3, m0w := E4: 冪等
  完全性: m2 + m0s + τE3 + τE4 = τ·(E1+E3+E4) = τ·I_sym
  rank (2·trace): tr2 = [τ²+τ, 2τ, 2τ, 2, 0, 0] →
  rank P₂ = (τ²+τ)/2 − 1, rank P₁ = τ, rank P₀s = rank P₀w = 1
  ⇒ d=4: (5, 3, 1, 1) / d=3: (2, 2, 1, 1) — 合計 d(d+1)/2 ✓

4D インスタンス (座標 (τ̂,x,y,z) = 添字 (0,1,2,3), q ∥ ŷ = 添字 2):
  S₃ := E_xx + E_zz は 4D では純 spin-0 ではない —
  **3·S₃ = 2·u0 + u20** (u0 = E_τ̂τ̂+E_xx+E_zz = θ [P₀s 固有], u20 = −2E_τ̂τ̂+E_xx+E_zz
  [P₂ 固有]) が正規化前の分解 (正規化後 S₃ = √(2/3)e₀ + e₂₀/√3)。
  D = E_xx − E_zz と X = E_xz + E_zx は 4D P₂ の固有ベクトル (真の TT)。

実行: cd proofs && ~/.elan/bin/lean ProjectorND.lean
-/

-- ================= τ の多項式 (ℤ 係数, 低次 → 高次) =================

def padd : List Int → List Int → List Int
  | [], b => b
  | a, [] => a
  | a :: as, b :: bs => (a + b) :: padd as bs

def pscale (c : Int) (p : List Int) : List Int := p.map (c * ·)

def pmul : List Int → List Int → List Int
  | [], _ => []
  | a :: as, b => padd (pscale a b) (0 :: pmul as b)

/-- 末尾の 0 を落とす正規形 -/
def ptrim (p : List Int) : List Int := (p.reverse.dropWhile (· == 0)).reverse

/-- τ での評価 -/
def peval (p : List Int) (t : Int) : Int := p.foldr (fun c acc => c + t * acc) 0

-- ================= 形式代数: V6 = Σ c_i(τ) E_i =================

abbrev V6 := List (List Int) -- 長さ 6

def zeroV : V6 := [[], [], [], [], [], []]

def basisE (i : Nat) : V6 := (List.range 6).map (fun k => if k == i then [1] else [])

def vadd (a b : V6) : V6 := List.zipWith padd a b

def vscale (p : List Int) (a : V6) : V6 := a.map (pmul p)

def vnorm (a : V6) : V6 := a.map ptrim

/-- 乗積表: (i, j) ↦ [(基底添字, 係数多項式)]。tau = [0,1]。 -/
def tbl (i j : Nat) : List (Nat × List Int) :=
  let one : List Int := [1]
  let tau : List Int := [0, 1]
  match i, j with
  | 0, 0 => [(0, one)]
  | 0, 1 => [(1, one)]
  | 0, 4 => [(4, one)]
  | 1, 0 => [(1, one)]
  | 1, 1 => [(1, tau)]
  | 1, 4 => [(4, tau)]
  | 2, 2 => [(2, one)]
  | 3, 3 => [(3, one)]
  | 3, 5 => [(5, one)]
  | 4, 3 => [(4, one)]
  | 4, 5 => [(1, one)]
  | 5, 0 => [(5, one)]
  | 5, 1 => [(5, tau)]
  | 5, 4 => [(3, tau)]
  | _, _ => []

def mulE (a b : V6) : V6 :=
  (List.range 6).foldl
    (fun acc i =>
      (List.range 6).foldl
        (fun acc2 j =>
          let cij := pmul (a.getD i []) (b.getD j [])
          (tbl i j).foldl
            (fun acc3 kc => vadd acc3 (vscale (pmul cij kc.2) (basisE kc.1)))
            acc2)
        acc)
    zeroV

def veq (a b : V6) : Bool := vnorm a == vnorm b

-- 単位元 (完全性の対象): I_sym = E1 + E3 + E4
def unitI : V6 := vadd (basisE 0) (vadd (basisE 2) (basisE 3))

-- スケール射影: m2 = τE1 − E2, m1 = E3, m0s = E2, m0w = E4
def m2f : V6 := vadd (vscale [0, 1] (basisE 0)) (vscale [-1] (basisE 1))
def m1f : V6 := basisE 2
def m0sf : V6 := basisE 1
def m0wf : V6 := basisE 3

-- ---------------- 定理 1: 乗積表の結合律 (代数の無矛盾性) ----------------

theorem assoc_table :
    ((List.range 6).all fun i => (List.range 6).all fun j => (List.range 6).all fun k =>
      veq (mulE (mulE (basisE i) (basisE j)) (basisE k))
        (mulE (basisE i) (mulE (basisE j) (basisE k)))) = true := by native_decide

-- ---------------- 定理 2: 単位元 I = E1 + E3 + E4 ----------------

theorem unit_law :
    ((List.range 6).all fun i =>
      veq (mulE unitI (basisE i)) (basisE i) && veq (mulE (basisE i) unitI) (basisE i)) = true := by
  native_decide

-- ---------------- 定理 3: 冪等性 (スケール形, 任意 τ) ----------------

theorem idempotency :
    (veq (mulE m2f m2f) (vscale [0, 1] m2f) &&
     veq (mulE m1f m1f) m1f &&
     veq (mulE m0sf m0sf) (vscale [0, 1] m0sf) &&
     veq (mulE m0wf m0wf) m0wf) = true := by native_decide

-- ---------------- 定理 4: 直交性 (12 対, 任意 τ) ----------------

theorem orthogonality_nd :
    (veq (mulE m2f m1f) zeroV && veq (mulE m1f m2f) zeroV &&
     veq (mulE m2f m0sf) zeroV && veq (mulE m0sf m2f) zeroV &&
     veq (mulE m2f m0wf) zeroV && veq (mulE m0wf m2f) zeroV &&
     veq (mulE m1f m0sf) zeroV && veq (mulE m0sf m1f) zeroV &&
     veq (mulE m1f m0wf) zeroV && veq (mulE m0wf m1f) zeroV &&
     veq (mulE m0sf m0wf) zeroV && veq (mulE m0wf m0sf) zeroV) = true := by native_decide

-- ---------------- 定理 5: 完全性 τ·(P₂+P₁+P₀s+P₀w) = τ·I ----------------

theorem completeness_nd :
    veq (vadd m2f (vadd m0sf (vadd (vscale [0, 1] m1f) (vscale [0, 1] m0wf))))
      (vscale [0, 1] unitI) = true := by native_decide

-- ---------------- 定理 6: scalar block の matrix units (転送演算子) ----------------
-- P₀sw = E5/√τ, P₀ws = E6/√τ: E5∘E6 = E2 (= τP₀s), E6∘E5 = τE4 (= τP₀w),
-- E5∘E5 = E6∘E6 = 0, 吸収則 m0s∘E5 = τE5, E5∘m0w = E5, m0w∘E6 = E6, E6∘m0s = τE6

theorem matrix_units :
    (veq (mulE (basisE 4) (basisE 5)) m0sf &&
     veq (mulE (basisE 5) (basisE 4)) (vscale [0, 1] m0wf) &&
     veq (mulE (basisE 4) (basisE 4)) zeroV &&
     veq (mulE (basisE 5) (basisE 5)) zeroV &&
     veq (mulE m0sf (basisE 4)) (vscale [0, 1] (basisE 4)) &&
     veq (mulE (basisE 4) m0wf) (basisE 4) &&
     veq (mulE m0wf (basisE 5)) (basisE 5) &&
     veq (mulE (basisE 5) m0sf) (vscale [0, 1] (basisE 5)) &&
     veq (mulE m2f (basisE 4)) zeroV && veq (mulE (basisE 4) m2f) zeroV &&
     veq (mulE m1f (basisE 4)) zeroV && veq (mulE (basisE 5) m1f) zeroV) = true := by
  native_decide

-- ---------------- 2·trace 汎関数と rank 公式 ----------------
-- tr2(E) = [τ²+τ, 2τ, 2τ, 2, 0, 0] (= 2 Tr)。縮約計算による (具体認証は下の
-- table_cert_d3/d4 が兼ねる)。rank P = Tr P:
--   2τ·rank P₂ = tr2(m2) = τ³+τ²−2τ, rank P₁ = tr2(E3)/2 = τ,
--   rank P₀s: tr2(m0s)/(2τ) = 1, rank P₀w = tr2(E4)/2 = 1

def tr2v : List (List Int) := [[0, 1, 1], [0, 2], [0, 2], [2], [], []]

def tr2 (a : V6) : List Int :=
  ptrim ((List.range 6).foldl (fun acc i => padd acc (pmul (a.getD i []) (tr2v.getD i []))) [])

theorem rank_formulas :
    (tr2 m2f == ptrim [0, -2, 1, 1] &&      -- τ³+τ²−2τ = 2τ·((τ²+τ)/2 − 1)
     tr2 m1f == [0, 2] &&                    -- 2τ ⇒ rank P₁ = τ
     tr2 m0sf == [0, 2] &&                   -- 2τ ⇒ rank P₀s = 1 (÷2τ)
     tr2 m0wf == [2] &&                      -- 2 ⇒ rank P₀w = 1
     -- 総和 = 2·dim Sym² = (τ+1)(τ+2): tr2(I) = τ²+3τ+2
     tr2 unitI == [2, 3, 1] &&
     -- インスタンス: d=4 (τ=3): rank P₂ = 5, P₁ = 3 / d=3 (τ=2): 2, 2
     peval (tr2 m2f) 3 == 2 * 3 * 5 && peval (tr2 m1f) 3 == 2 * 3 &&
     peval (tr2 m2f) 2 == 2 * 2 * 2 && peval (tr2 m1f) 2 == 2 * 2) = true := by
  native_decide

-- ================= 具体インスタンスによる乗積表の認証 =================
-- d 次元, 縦方向 = axis (q ∥ 当該軸)。θ = diag(1,…,0,…,1), ω = E_{axis,axis}。
-- 語 w_i = 2E_i (半整数回避のため 2 倍)。合成 (2A)∘(2B) = 4(A∘B) なので
--   comp(w_i, w_j) = 2·Σ_k c_k(τ_d)·w_k  を全成分で検査する。
-- 表の係数は τ の高々 1 次 (縮約で閉じる θ ループは高々 1 個) なので、
-- τ = 2 (d=3) と τ = 3 (d=4) の 2 点一致が任意 τ の表を確定する。

def delc (i j : Nat) : Int := if i == j then 1 else 0
def omc (axis i j : Nat) : Int := if i == axis && j == axis then 1 else 0
def thc (axis i j : Nat) : Int := delc i j - omc axis i j

def wrd (axis : Nat) (w : Nat) (i j k l : Nat) : Int :=
  match w with
  | 0 => thc axis i k * thc axis j l + thc axis i l * thc axis j k
  | 1 => 2 * thc axis i j * thc axis k l
  | 2 => thc axis i k * omc axis j l + thc axis i l * omc axis j k +
         omc axis i k * thc axis j l + omc axis i l * thc axis j k
  | 3 => 2 * omc axis i j * omc axis k l
  | 4 => 2 * thc axis i j * omc axis k l
  | _ => 2 * omc axis i j * thc axis k l

def sumd (d : Nat) (f : Nat → Int) : Int := (List.range d).foldl (fun a i => a + f i) 0

def compc (d axis wi wj i j k l : Nat) : Int :=
  sumd d fun a => sumd d fun b => wrd axis wi i j a b * wrd axis wj a b k l

/-- 表の RHS (2 倍スケール): 2·Σ_k c_k(τ)·w_k -/
def tblRhs (_d axis wi wj : Nat) (tau : Int) (i j k l : Nat) : Int :=
  2 * (tbl wi wj).foldl (fun acc kc => acc + peval kc.2 tau * wrd axis kc.1 i j k l) 0

def certTable (d axis : Nat) (tau : Int) : Bool :=
  (List.range 6).all fun wi => (List.range 6).all fun wj =>
    (List.range d).all fun i => (List.range d).all fun j =>
      (List.range d).all fun k => (List.range d).all fun l =>
        compc d axis wi wj i j k l == tblRhs d axis wi wj tau i j k l

def certTrace (d axis : Nat) (tau : Int) : Bool :=
  (List.range 6).all fun wi =>
    sumd d (fun i => sumd d fun j => wrd axis wi i j i j) == peval (tr2v.getD wi []) tau

theorem table_cert_d3 : (certTable 3 1 2 && certTrace 3 1 2) = true := by native_decide

theorem table_cert_d4 : (certTable 4 2 3 && certTrace 4 2 3) = true := by native_decide

-- ================= 4D インスタンスの辞書 (座標 (τ̂,x,y,z), q ∥ ŷ = 添字 2) =================

def applyW (d axis w : Nat) (v : Nat → Nat → Int) (i j : Nat) : Int :=
  sumd d fun k => sumd d fun l => wrd axis w i j k l * v k l

/-- スケール射影の 4D 具体形: 2m2 = τ·w0 − w1 (τ = 3) の適用 -/
def applyM2d4 (v : Nat → Nat → Int) (i j : Nat) : Int :=
  3 * applyW 4 2 0 v i j - applyW 4 2 1 v i j

def eS3 (i j : Nat) : Int := delc i 1 * delc j 1 + delc i 3 * delc j 3 -- E_xx + E_zz
def eU0 (i j : Nat) : Int := delc i 0 * delc j 0 + delc i 1 * delc j 1 + delc i 3 * delc j 3
def eU20 (i j : Nat) : Int := -2 * (delc i 0 * delc j 0) + delc i 1 * delc j 1 + delc i 3 * delc j 3
def eD4 (i j : Nat) : Int := delc i 1 * delc j 1 - delc i 3 * delc j 3   -- E_xx − E_zz
def eX4 (i j : Nat) : Int := delc i 1 * delc j 3 + delc i 3 * delc j 1   -- E_xz + E_zx

def idx4d : List Nat := [0, 1, 2, 3]

/-- 定理: S₃ の 4D 分解 3·S₃ = 2·u0 + u20 (正規化後 S₃ = √(2/3)e₀ + e₂₀/√3)。
    S₃ は 4D の純 spin-0 ではない — spin-0 (u0) と spin-2 (u20) の混合。 -/
theorem s3_4d_decomposition :
    (idx4d.all fun i => idx4d.all fun j =>
      3 * eS3 i j == 2 * eU0 i j + eU20 i j) = true := by native_decide

/-- 定理: 分解の成分は正しい 4D spin 固有ベクトル —
    u0 は P₀s 固有 (E2 u0 = τ u0, P₂ u0 = 0)、u20 は P₂ 固有 (m2 u20 = τ u20,
    E2/E3/E4 で消える)、u0 ⊥ u20 ⊥ (u0, u20 とも横空間内)。 -/
theorem s3_components_eigen :
    ((idx4d.all fun i => idx4d.all fun j =>
      applyW 4 2 1 eU0 i j == 2 * 3 * eU0 i j &&      -- w1 = 2E2: 2·E2 u0 = 2τ u0
      applyM2d4 eU0 i j == 0 &&
      applyM2d4 eU20 i j == 2 * 3 * eU20 i j &&        -- 2m2 u20 = 2τ u20
      applyW 4 2 1 eU20 i j == 0 &&
      applyW 4 2 2 eU20 i j == 0 &&
      applyW 4 2 3 eU20 i j == 0) &&
     -- 直交性 ⟨u0, u20⟩ = 0
     (sumd 4 (fun i => sumd 4 fun j => eU0 i j * eU20 i j) == 0)) = true := by
  native_decide

/-- 定理: D = E_xx − E_zz と X = E_xz + E_zx は 4D P₂ の固有ベクトル (真の TT —
    q₀ = 0, q ∥ ŷ でも 4D transverse-traceless)。他の 3 射影に消される。 -/
theorem dx_are_4d_tt :
    (idx4d.all fun i => idx4d.all fun j =>
      applyM2d4 eD4 i j == 2 * 3 * eD4 i j &&
      applyW 4 2 1 eD4 i j == 0 && applyW 4 2 2 eD4 i j == 0 &&
      applyW 4 2 3 eD4 i j == 0 &&
      applyM2d4 eX4 i j == 2 * 3 * eX4 i j &&
      applyW 4 2 1 eX4 i j == 0 && applyW 4 2 2 eX4 i j == 0 &&
      applyW 4 2 3 eX4 i j == 0) = true := by native_decide

/-- 定理: 4D Ward 収縮 — P₂ と P₀s の縦添字成分 (μ = 2 = q 方向) は恒等的に零。 -/
theorem ward_4d :
    (idx4d.all fun j => idx4d.all fun k => idx4d.all fun l =>
      applyM2d4 (fun a b => if a == k && b == l then 1 else 0) 2 j == 0 &&
      applyW 4 2 1 (fun a b => if a == k && b == l then 1 else 0) 2 j == 0) = true := by
  native_decide

/-- 定理: 4D ゲージモード h = q⊗ξ + ξ⊗q (q ∥ ŷ) は P₂ と P₀s に消される。 -/
theorem gauge_annihilated_4d :
    (idx4d.all fun xi =>
      let hg := fun a b => (if a == 2 then delc b xi else 0) + (if b == 2 then delc a xi else 0)
      idx4d.all fun i => idx4d.all fun j =>
        applyM2d4 hg i j == 0 && applyW 4 2 1 hg i j == 0) = true := by native_decide
