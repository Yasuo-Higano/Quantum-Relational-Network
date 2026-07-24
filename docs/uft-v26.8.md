# 量子情報網理論 v26.8 — 連続極限 universality 監査 (0: source matching と定義閉包 — Gate 1 開通)

**第二十七期第 8 版 (進行中 — 本節は v26.8-0)。判定 (a) — 0a: Lean 14 定理 / 0b: 8 検査 PASS。**
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

## 1. 登録済みの残り (spec §12.3–12.5 — 本版の A/B/C)

- **v26.8-A**: 解析 one-loop oracle の二重導出 (直接 Feynman 積分 vs 曲率 form
  factor/Weyl anomaly) — Gate 2: 一致まで数値実装禁止。規約凍結リストは spec §12.3。
- **v26.8-B**: staggered TT continuum limit — null-combination 推定器 (q⁰/q²/q⁴
  counterterm を代数的に消す)・三者一致 (spectral/dispersion/直接)・continuum
  trajectory (am = a·m_phys)。massless/massive 二系列。
- **v26.8-C**: Wilson 独立離散化 (時間連続 Hamiltonian — taste trap 回避)・
  外挿 2 モデル事前登録 — PRED-016。spin-0 sum rule (PRED-017) は operator
  benchmark として並走。

## 2. 成果物

0a: `proofs/ProjectorND.lean` (定理 14 本 — Lean 計 50 本/7 ファイル)。
0b: `sim/src/bin/v268z_source_matching.rs` / `results/v268z_source_matching.txt`
(8 検査 PASS) / `results/v268z_source_matching.json` / PRED-015 (**scored-hit**)。
claims: QRN-GRAV-043。
