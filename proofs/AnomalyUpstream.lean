/-
v21.6 アノマリー探索の上流形式化 — 表の係数を定理にする (PROMPT/3 §6)

これまでの Lean 定理 (Anomaly*.lean) は「与えた係数表で探索が完全」を検証した。
本ファイルはその上流: 表にハードコードされていた SU(3) の係数
(Dynkin 指数 2T・立方アノマリー係数 A) を、既約表現 (p,q) の閉形式
  dim(p,q) = (p+1)(q+1)(p+q+2)/2
  C₂(p,q)  = (p² + q² + pq + 3p + 3q)/3
  T(p,q)   = dim·C₂/8          (dim adj = 8)
  A(p,q)   = (p−q)(p+2q+3)(q+2p+3)·dim/60   (規格化 A(1,0) = 1)
から導出される定理に昇格させる。検算アンカー: A(6)=(2,0) → 7, A(10)=(3,0) → 27,
A(8)=(1,1) → 0 (実表現), T(6) → 5/2 (いずれも教科書値)。

定理:
  - table_entries_from_formulas: 探索表の {1, 3, 3̄} の (2T, A) が閉形式の値に一致 (decide)
  - anchors: A(6)=7, A(10)=27, A(8)=0, 2T(6)=5, 2T(8)=6 (decide)
  - a_conj_antisym_window / t_conj_sym_window / dim_conj_sym_window:
      窓 p,q ≤ 12 の全対で A(p,q) = −A(q,p), T 対称, dim 対称 (native_decide)
      — 「共役で数える」同値類規約の上流正当化
  - a_real_zero_window: 自己共役 (p=q) で A = 0 (native_decide, 窓全域)
SU(2) 側: 2T(doublet) = 1 を weight 列 {±1} (2m 単位) の Σ(2m)²/2 から (decide)。

未形式化 (限界として台帳へ): アノマリー方程式そのものの群論的導出 (一般 N)・
J=0 構造零の代数証明 (実スペクトル定理 = mathlib 級)・index-3 数論 no-go・
B−L Plücker 分類。

実行: ~/.elan/bin/lean proofs/AnomalyUpstream.lean  (終了コード 0 = 全定理検証)
-/

-- ---------------- SU(3) 既約表現 (p,q) の閉形式 (60 倍整数化で Rat を回避) ----------------

def dim2 (p q : Int) : Int := (p + 1) * (q + 1) * (p + q + 2)  -- 2·dim

-- 整数化: T(p,q) = dim·C₂/8 = (dim2/2)·(poly/3)/8 = dim2·poly/48 — t48 := dim2·poly。
--   検算 (1,0): dim2 = 6, poly = 4 → t48 = 24 → T = 1/2 ✓
def t48 (p q : Int) : Int := dim2 p q * (p * p + q * q + p * q + 3 * p + 3 * q)

theorem t48_normalization : t48 1 0 = 24 := by decide  -- T = 24/48 = 1/2 (基本表現)

-- 120·A = (p−q)(p+2q+3)(q+2p+3)·dim2 (A = .../120): 検算 (1,0): 1·4·5·6 = 120 → A = 1 ✓
def a120 (p q : Int) : Int := (p - q) * (p + 2 * q + 3) * (q + 2 * p + 3) * dim2 p q

-- ---------------- 探索表 {1, 3, 3̄} の係数が公式から出る ----------------
-- 表 (Anomaly.lean の reps): 単一項 (色) の a3 ∈ {0, +1, −1}, 2T(色) ∈ {0, 1}
theorem table_singlet : a120 0 0 = 0 ∧ t48 0 0 = 0 := by decide
theorem table_fund : a120 1 0 = 120 ∧ t48 1 0 = 24 := by decide      -- A=+1, 2T=1
theorem table_antifund : a120 0 1 = -120 ∧ t48 0 1 = 24 := by decide -- A=−1, 2T=1

-- ---------------- 教科書アンカー (公式の外部検証) ----------------
theorem anchor_sextet : a120 2 0 = 7 * 120 ∧ 2 * t48 2 0 = 5 * 48 := by decide  -- A(6)=7, 2T=5
theorem anchor_decuplet : a120 3 0 = 27 * 120 := by decide                       -- A(10)=27
theorem anchor_adjoint : a120 1 1 = 0 ∧ 2 * t48 1 1 = 6 * 48 := by decide        -- A(8)=0, 2T=6

-- ---------------- SU(2): doublet の 2T = 1 が weight 列から出る ----------------
def su2Weights2m : List Int := [1, -1]  -- doublet の 2m
theorem su2_doublet_index : (su2Weights2m.map (fun m => m * m)).foldl (· + ·) 0 = 2 := by
  decide  -- Σ(2m)² = 2 → 2T = Σ(2m)²/2 = 1 (表の t2 = 1)

-- ---------------- 窓全域の構造定理 (native_decide) ----------------
def window : List (Nat × Nat) :=
  (List.range 13).flatMap (fun p => (List.range 13).map (fun q => (p, q)))

-- 共役 (p,q) ↔ (q,p) で A は反対称, T と dim は対称 — 「共役軌道で数える」規約の上流
theorem a_conj_antisym_window :
    window.all (fun pq => a120 (pq.1 : Int) (pq.2 : Int) = -a120 (pq.2 : Int) (pq.1 : Int)) = true := by native_decide

theorem t_conj_sym_window :
    window.all (fun pq => t48 (pq.1 : Int) (pq.2 : Int) = t48 (pq.2 : Int) (pq.1 : Int)) = true := by native_decide

theorem dim_conj_sym_window :
    window.all (fun pq => dim2 (pq.1 : Int) (pq.2 : Int) = dim2 (pq.2 : Int) (pq.1 : Int)) = true := by native_decide

-- 自己共役表現 (p=q) はアノマリー・フリー
theorem a_real_zero_window :
    (List.range 13).all (fun p => a120 (p : Int) (p : Int) = 0) = true := by native_decide

-- dim の正値性 (窓) — 公式の健全性
theorem dim_pos_window :
    window.all (fun pq => 0 < dim2 (pq.1 : Int) (pq.2 : Int)) = true := by native_decide
