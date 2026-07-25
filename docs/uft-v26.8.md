# 量子情報網理論 v26.8 — 連続極限 universality 監査 (0: source matching / A: 解析 oracle / B: staggered 連続極限)

**第二十七期第 8 版 (進行中)。判定 (a) — 0a: Lean 14 定理 / 0b: 8 検査 / A: 13 検査 / B: 6 検査 PASS。**
PROMPT/9 (spec §12) の中心命題の再定義に従う:

> 次の中心課題は q⁴ln q² の測定ではない。**BOND-A の格子 source が、連続極限で
> 2-taste Dirac の保存された Hilbert/Belinfante stress tensor に流れるかを証明する**。
> 失敗した場合、経路 B は「重力真空偏極」ではなく格子フェルミオンの弾性・strain
> response であり、gravitational interpretation を終了する (Gate 1)。

## 0a. ProjectorND.lean — 公理からの抽象 Barnes–Rivers 代数 (Lean 定理 14 本)

有限格子 {0..8}³ certificate (Projector.lean) を卒業し、公理 θ²=θ, ω²=ω, θω=ωθ=0,
θ+ω=I, trθ=τ, trω=1 から証明する。rank-1 ω の下で BR の語は 6 個 {E1..E6} で閉じ、
乗積表 (係数は τ の高々 1 次) の形式代数 (ℤ[τ]) として: 結合律 (216 triple)・単位元・
冪等・直交 12 対・完全性・**scalar block の matrix units (P₀sw/P₀ws)**・rank 公式
(d=4 で **(P₂,P₁,P₀s,P₀w) = (5,3,1,1)**、d=3 で (2,2,1,1) — 総和 d(d+1)/2) を機械検証。
表は d=3 (τ=2)/d=4 (τ=3) の具体行列で全数認証 (係数がアフィンなので 2 点で確定)。

**4D 辞書** (座標 (τ̂,x,y,z), q∥ŷ): **S₃ = (E_xx+E_zz)/√2 は 4D の純 spin-0 ではない** —
3·S₃ = 2·u0 + u20 (u0 = θ: P₀s 固有 / u20 = −2E_τ̂τ̂+E_xx+E_zz: P₂ 固有, 直交) を定理化
(正規化形 S₃ = √(2/3)e₀ + e₂₀/√3)。**D = E_xx−E_zz と X = E_xz+E_zx は 4D P₂ の
固有ベクトル (真の TT)** — 最小 TT universality test は D/X で実施可能。以後
「S を 4D spin-0 と呼ぶこと」「h_00, h_0i 実装前に full gravitational polarization と
呼ぶこと」は禁止 (spec §12.2-2)。完全共変核の取得後にのみ K_SS = (2/3)F₀ + (1/3)F₂。

## 0b. v268z_source_matching — Gate 1 の判定 (PRED-015 hit)

器械: 8 成分折込み基底で H(k) = Σcos k_i·Γ_i + m·Γ_m (Γ は厳密に反交換 — [S0])。
Clifford 単項 16 個 = spin⊗I_taste の全体、taste 可換子環は twirl 像 rank = 4 = **M₂
(2 tastes)** — [S1]。taste-singlet 射影 = Clifford 像への Frobenius 射影
(I_taste⊗tr_taste/2 と等価)。一般頂点公式 (parity ε, 変位 d, 重み w) を導出し、
N=8 周期格子の位置空間演算子の平面波 sandwich と **5.1e-15** で照合 — [S2]
(解析形と v26.x 器械系の接続)。matching ladder ε ∈ {0.4, 0.2, 0.1, 0.05}
(|p| = ε, |q| = 0.7ε, 4 方向):

| ε | r_D | m_D | Z_D | r_X naive | r_X split | m_X split | Z_X |
|---|---|---|---|---|---|---|---|
| 0.40 | 5e-16 | 2.5e-3 | 0.983 | **1.000** | 6e-16 | 4.8e-3 | 1.846 |
| 0.05 | 4e-15 | **3.9e-5** | **0.9997** | **1.000** | 5e-15 | **7.8e-5** | **1.9975** |

- **[S3] D は TreeLevelMatchedTTSource**: 対角頂点は Γ の多項式そのもの —
  **taste-singlet 性は構成的定理** (r_D ~ 4e-15)。shape 残差はきれいな O(ε²)
  (improvement 級の許容差)、Z_D → 1。**Gate 1 開通 — 経路 B は「strain response」
  ではなく tree-level で本物の TT source を持つ。**
- **[S4] 素朴 BOND-A off-diagonal 転写 (spec §2 の暫定則) は棄却**: r = **1.000**
  (全 ε) — η_x(≡1) 位相 z ホップは Λ_z = I⊗I⊗σ₃ 構造で **Clifford 像に直交
  (100% taste-nonsinglet)**。「認証前に物理を主張しない」とした §2 の警戒が的中。
- **[S5] 修正転写 point-split X の構成と認証**: 4 隅 point-split
  (d = ±x̂±2ẑ / ±ẑ±2x̂, σρ 符号交代 + η parity) は cos 項を消して sin·sin 積のみを
  残し、**厳密に taste-singlet** (r ~ 5e-15): Γx sin kx sin 2kz + Γz sin kz sin 2kx →
  ノードで αx pz + αz px ∝ T_xz。m = O(ε²) → 7.8e-5。
- [S6] **Z_D = 0.9997 ≠ Z_X = 1.9975** — D と X は別の正規化 (立方既約表現の分裂,
  spec §12.2-4 どおり bare 一致は要求しない — 個別繰り込み後の比較が universality)。
- [S8] 変異: η parity 破り (piece2 ε: (1,1,0)→(0,1,0)) → r = 0.707 で検出。

**PRED-015 hit** (D 通過 ∧ 認証済み X 転写の存在; naive 棄却は登録済みの分岐)。

## 0c. 開発記録

1. X_split run1 は 2 片の相対符号を誤実装 — 結果は反対称結合 −αx pz + αz px
   (回転生成子) で、**taste 残差 r は両片厳密 0 のため不感、shape ゲート m ≈ 2.0 /
   Z ≈ −0.25 が検出**した。symmetric stress には両片 w = −σρ/16。「r だけでは
   足りない — shape (m) と Z が matching certificate の別の目」。
2. 変異ゲート初版 (corner 符号反転) は r に不感と判明 — piece1 の全 Fourier 成分が
   Γx⊗I 構造 (1-link x 変位は常に (−1)^{sx} 重み) で singlet のまま。taste 構造を
   守るのは **η parity** — 変異は parity 破りに再設計 (r = 0.707 で検出)。

## A. v268a_oracle — 解析 one-loop oracle の三重経路一致 (Gate 2 開通)

規約を凍結 (Belinfante 頂点 Γ_ij = ¼[α_i(2p+q)_j + (i↔j)]、u†u = 1、
χ = ∫dE 2σ/E、null-combination Σw = Σwq² = Σwq⁴ = 0, Σwq⁴ln q² = 1 —
v268z の matching certificate と同一規約チェーン) した上で、1 Dirac flavor の
TT (D) チャネル oracle を**三重経路**で導出した (13 検査 PASS):

- **閉形式の導出と認証**: スピン和 tr[P₊(a·α)P₋(a·α)] =
  (1/E₁E₂)[(E₁E₂+k₁·p+m²)a² − 2(p·a)²] が明示 4×4 行列トレースと 8.9e-16。
  φ 平均で **D と X̂ = (T_xz+T_zx)/√2 の被積分が恒等に一致** (spin-2 の 2 偏極が
  同一 form factor — σ_X̂ = σ_D を独立実装で 6.7e-16)。
- **Lorentz 不変性の器械証明**: σ_D(E;q) が s = E²−q² のみに依存 (異 q 同 s で
  1.6e-15) ⇒ KL 表示 χ(q) = ∫ds ρ(s)/(s+q²) が厳密に従う。massless の
  **ρ_D(s) = s²/(160π²) (閉形式・厳密)**、massive も閉形式
  (pE_p/4π²)[(2/3)p² − (4/15)p⁴/E_p²] を導出・認証 (5.6e-16 / 1e-8 級)。
- **三重経路の一致 (Gate 2)**: Route I (直接ループの球座標求積) = Route II
  (吸収部→分散, 安定通分形 K(s) = n₀/Π(s+qᵢ²)) = Route III (解析):
  **A = −1/(160π²) = −6.33257e-4、すなわち 16π²A = −1/10 (厳密)** —
  I/II 相対差 4.6e-8、II/III 3.1e-9。**A(2 taste) = −1/(80π²) が v26.8-B の
  格子比較の分母 (PRED-016)**。
- **λ スケール不変 1.8e-15 (branch α)**: 非局所形は純 q⁴ln q² — v26.7.1 の
  M² 修正と同様、q³ 型の可能性は KL 表示 + ρ の解析性が排除。
- **massive decoupling (PRED-018 の oracle 側)**: A(m)/A(0) = 0.61 → 0.023
  (m/q̄ = 0.5 → 4)、大質量冪 **m⁻¹·⁹⁰ ≈ (q/m)²** — 教科書どおりの decoupling。
- **スカラー和則 (自前規約)**: ∫ρ_θ(s)/s³ ds = **1/(80π²)** (3.5e-12,
  m ∈ {0.5,1,2} で不変)。文献規約 σ_f = ρ_θ/(3s³) では 1/(240π²)、2 taste で
  1/(120π²) — PRED-017 の的の oracle 側。
- 開発記録 (run1→run3): (i) 素朴な Σwᵢσ(E;qᵢ) は大 E 域の **f64 桁落ち**
  (E⁴ 項の解析的相殺が数値では 1e-16×|w|×積分域で崩壊) で O(10⁷) 倍の誤差 —
  **通分形 K(s) = n₀/Π(s+qᵢ²) (n₃ = n₂ = n₁ = 0 が Σw 拘束で恒等)** に書き換えて
  桁落ちゼロ。(ii) Route I は円錐折れ目と遠方 Jacobian で発散級 — 球座標 +
  閾値パネル + 冪外挿 tail に再設計。(iii) スカラー照合の引数ミス (q=0 の対は
  k₁ = p であって −p でない)。「解析的に消える項は数値でも消してから積分せよ」。

## B. v268b_continuum — staggered TT の連続極限 (PRED-016 前半 hit)

**経路 B の最重要 falsifier の前半が通った。** 8 成分折込み基底の H(k) を無限体積
格子の厳密な 1 粒子 Hamiltonian として、χ_D の BZ 積分 (ノード中心セル [0,π)³,
3 段入れ子箱 + ノード点分割の GL 求積) を continuum trajectory
(q^lat = a·Q, m^lat = a·m_phys — lattice-unit m 固定の a→0 は禁止) で測った。
null 結合は k 点ごとに直接累積 (相殺深さ (aQ)⁴ ≫ f64 床 — v26.8-A の教訓の適用)。

- **[S0] 器械回帰が決定的**: BZ 積分の χ_D(2π/16) が v26.7-II の有限体積公表値と
  **4.1e-7** で一致 (0.154068 = 0.154068) — 新エンジンの規格化・頂点・占有の全てが
  認証済み器械系と接続された。
- **[S2] 主結果 — PRED-016 前半 hit**:

  | a | A(a)/2A_oracle |
  |---|---|
  | 0.50 | 1.1653 |
  | 0.25 | 1.0660 |
  | 0.125 | 1.0222 |
  | 0.0625 | 1.0065 |
  | **0.032** | **1.0014** |

  単調な O(a²) 級の収束。2 モデル外挿 (a² / a+a²) = 1.0102 / 0.9908 —
  **A₀/(2A_oracle) = 1 ± 0.02 (spread 1.9%)**。**staggered の繰り込み後 TT form
  factor は連続 2-taste Dirac に収束する** (taste 数 2 も同時に測定されたことになる
  — 1 taste なら比は 0.5 に収束したはず)。
- **[S3] q̂ = (2/a)sin(aq/2) 変種**: 外挿一致 1.1% (spec §12.1 の二重運動量規則)。
- **[S4] massive decoupling (PRED-018 の格子側)**: A₀(m_phys=1)/A₀(0) = 0.1793 vs
  oracle 比 0.1741 (閉形式 ρ から計算) — **3.0% 一致**。oracle 側の冪 (q/m)²
  (v26.8-A S8) と合わせ PRED-018 は両側から立った。
- [S1] 求積自己整合 4.1e-4 / [S5] 変異 (η parity 破り) 検出。
- 開発記録 (run1→run5): (i) ノード構造 (スケール aQ_min) を解像しない箱設計を
  S1 が 35% で検出 → 3 段入れ子 + ノード点分割へ。(ii) massive 対照の初版は
  m_phys = 2.0 で am = 0.7 の artifact 支配 → m_phys = 1.0 + 細 a 系列に変更。
  (iii) oracle 比較値を placeholder 定数で置いた不備 → 閉形式 ρ_D(s;m) から
  binary 内で計算する形に修正 (「出所不明の定数は禁じ手」)。

## 1. 登録済みの残り (spec §12.5 — 本版の C)

- **v26.8-C**: Wilson 独立離散化 (時間連続・空間格子 Hamiltonian — taste trap
  回避、1 flavor なので的は 1×A_oracle)・外挿 2 モデル事前登録 — **PRED-016 の
  完成** (A^stag_D/A_or = 1 は達成済み; A^Wil_D, X チャネル 2 種が残り)。
  spin-0 sum rule (PRED-017) は operator benchmark として並走。
- 三者一致の残り (lattice spectral density / dispersion 再構成の 2 経路 —
  本版 B は直接 null projection のみ)。

## 2. 成果物

0a: `proofs/ProjectorND.lean` (定理 14 本 — Lean 計 50 本/7 ファイル)。
0b: `sim/src/bin/v268z_source_matching.rs` / `results/v268z_source_matching.txt`
(8 検査 PASS) / `results/v268z_source_matching.json` / PRED-015 (**scored-hit**)。
claims: QRN-GRAV-043。
A: `sim/src/bin/v268a_oracle.rs` / `results/v268a_oracle.txt` (13 検査 PASS) /
`results/v268a_oracle.json` — **A_oracle(1 Dirac) = −1/(160π²) 凍結**。
claims: QRN-GRAV-044。
B: `sim/src/bin/v268b_continuum.rs` / `results/v268b_continuum.txt` (6 検査 PASS) /
`results/v268b_continuum.json` — **A₀/(2A_oracle) = 1 ± 2% (PRED-016 前半 hit)**。
claims: QRN-GRAV-045。
