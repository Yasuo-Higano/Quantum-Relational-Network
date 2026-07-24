# v26.6–v26.8 事前登録: 重力真空偏極監査の仕様凍結

**登録日: 2026-07-24。この文書はコード実装より先にコミットされる (凍結)。**
PROMPT/8 の指示に従い、v26.6 以降の誘導重力経路 B を「pole search」ではなく
**full gravitational vacuum polarization の構成・Ward 認証・spectral 分解・
continuum universality audit** として仕様から凍結する。数値ゲート・判定分岐・
禁止暗黙変換は本文書が一次ソースであり、実装・実行後に変更しない
(ゲートの前提誤りが見つかった場合は v24.3/v26.2 の先例に従い開発記録として
両 run を公表する)。

## 0. 現状の正名 (v26.5 までの到達点の正確な表現)

> v26.5 は「誘導重力を発見した段階」ではなく、**自由フェルミオン格子に対する
> 外部計量応答の測定器を較正し、素朴な Newton 定数解釈を棄却した段階**である。

v26.2–26.5 が実際に検証したのは: 局所エネルギー密度 T_00 の分解と保存 (時間成分の
較正)、一本のボンド変調に対する Feynman–Hellmann 応答、静的感受率の等方性、
plus 偏極 1 成分の分解 — であり、完全な T_μν・完全な二点 Ward 恒等式・接触項では
ない。「要件 1–3 完遂」「spin-2 分解完了」という表現は過大であり、以後使わない。

## 1. 格子の存在論 (二者択一の凍結) — regulator を選ぶ

QRN はこれまで格子を「regulator」と「物理的微視構造」の二つの意味で使ってきた。
両者を混ぜると「非普遍だから失敗」と「microscopic model dependent だから予言」が
同時に使われてしまう。**v26.x 誘導重力経路では前者に固定する:**

- **格子は regulator である。** staggered 格子は試験 matter 理論 + UV 切断であり、
  QRN Core ではない。
- 帰結: (i) 連続極限 a→0 が必要。(ii) bare c₁ は scheme dependent であり、
  繰り込み条件なしに Newton 定数は予言できない。(iii) 比較・予言の対象は
  **共通の繰り込み条件を課した後の応答、非局所な q⁴ln q² の係数、spectral
  density、C_T/Weyl anomaly 対応係数**である。(iv) 独立離散化に bare 値の一致は
  要求しない (要件 10 の書き換え — §6)。
- 「格子/ネットワークが物理的微視構造」の路線は、**同じ Core が lattice
  spacing・graph・matter spectrum を決める QRN-Core v1 が定義されるまで主張資格を
  持たない**。その場合の検証様式 (bare 量の物理化、同一 RG basin 内の irrelevant
  変形への頑健性) は Core v1 仕様の側に登録する。

## 2. 結合の定義 (scheme BOND-A) — H[h] を二次まで

模型は v26.2–26.5 と同一の 3+1D staggered (η_x=1, η_y=(−1)^x, η_z=(−1)^{x+y},
staggered 質量 m(−1)^{x+y+z}, twist-x)。静的空間計量 g_ij(x) = δ_ij + h_ij(x) への
結合を次で**定義**する (正準規格化 ψ = g^{−1/4}ψ̃ 後の Hamiltonian 形式):

- **対角成分**: 方向 i のボンド b (中点 u_b) の振幅
  t_b → t_b · (1 + h_ii(u_b))^{−1/2}。
  展開 (1+h)^{−1/2} = 1 − h/2 + (3/8)h² + O(h³) により
  - 一次頂点: V_A = ∂H/∂ε_A = −(1/2) T_A(q) (T_A = v26.2–26.5 の stress 頂点、
    中点位相規約)
  - 二次接触項 (seagull): S_AB = ∂²H/∂ε_A∂ε_B
    = (3/4) δ_{dir(A),dir(B)} Σ_{b∈dir} c_A(u_b) c_B(u_b) t_b (c†c + h.c.)。
    **接触項は同一方向のボンドにのみ働く** (対角計量ではボンド因子が自方向の
    h_ii のみに依存するため)。
- **質量項は h に結合しない**: √g m ψ†βψ は正準規格化で m ψ̃†βψ̃ (√g が消える)。
  格子では on-site 項 (長さ要素を持たない) — scheme の定義として凍結。
- **off-diagonal 成分 (定義のみ — 実装は v26.7 登録)**: e = g^{−1/2} の行列則。
  h_xz = s のとき e = 1 − (s/2)(E_xz+E_zx) + (3/8)s²(E_xx+E_zz) + O(s³)。
  一次頂点は η-swap 混合ボンド (η_a 位相 × ĵ 方向ホップ、中点 sample) の対称化、
  二次接触項は x/z 対角ボンドに (3/8)s²。ブロック基底では π シフト対を結ぶ
  「奇運動量頂点」になる。**taste 汚染の有無 (連続極限で taste-singlet の
  ½{α,∂} に落ちるか) は v26.7 の認証課題** — 認証前に cross 偏極の物理は主張しない。
- **counterterm (tadpole 減算)**: Γ_ren[h] = E₀[h] − Λ Σ_x √(det g(x)),
  √g = Π_i (1+h_ii(x))^{1/2} (site sample)。**Λ := −e_⊥** (e_⊥ = y ボンド
  エネルギー密度。y↔z は W 対称性で厳密同値、x は twist の有限サイズ残差)。

**完全静的核** (モード規格化 w_q = 1/2 (q≠0 cos モード) / 1 (q=0 一様) で q=0 連続):

> k̂_AB(q) := (V w_q)^{−1} ∂²Γ_ren/∂ε_A∂ε_B
>          = (3/4) e_A δ_AB − (1/4) χ_AB(q) − (Λ/4)(J − 2I)_AB
>
> A, B ∈ {xx, yy, zz}, q ∥ ŷ。χ_AB は v26.3–26.5 規約の複素頂点 intensive 感受率
> (χ_AB = Σ 2 M_A M_B/ΔE / V)。J = 全 1 行列, I = 単位行列。

接触項と counterterm は**この scheme では q に依存しない** (y 並進不変性が
⟨T(2qŷ)⟩ = 0 を厳密に保証するため — [S4] で器械認証)。したがって q² 係数 c₁ は
接触完備化で不変であり、v26.3–26.5 の c₁ 結論は完全核でもそのまま成立する
(これ自体を [S5] で回帰認証する)。

## 3. Γ とは何か (このシリーズでの意味)

Γ[h] := E₀[h] (基底状態エネルギー汎関数 = T=0 静的 Euclidean 有効作用 / 単位時間)。
∂²Γ/∂h² は **ω=0 の静的核**である。ω≠0 (Lehmann 分母 2ΔE/(ΔE²−ω²)) と
有限サイズ spectral measure は v26.7 (§7)。「graviton propagator」の語は §6 の
昇格条件を満たすまで使わない。

## 4. v26.6 の検査と数値ゲート (凍結)

実装バイナリ: `sim/src/bin/v266_vacuum_pol.rs`。証明: `proofs/Projector.lean`
(Barnes–Rivers 型 projector 代数 — 冪等・直交・完全・trace(=rank)・Ward 収縮・
ŷ 方向インスタンスのチャネル辞書。**数値コードはこの証明済み規格に従う**)。

- **[S0] dense FD Hessian 完全性 (数学的恒等式)**: N=8, m∈{0,0.5}。
  モード集合 = {uni-xx, uni-yy, cos-xx(j=1), cos-yy(j=2), cos-zz(j=2),
  交差 uni(xx,yy), 交差 uni(xx,zz), 交差 cos(xx,yy)(j=2), cos-yy(j=2) vs
  sin-yy(j=2) の等値}。exact ボンド因子 (1+εc)^{−1/2} の H[ε] を再対角化した
  Richardson 中心差分 (ε = 0.02, 0.01) と、組立側 ⟨S⟩ − χ_V (dense Lehmann) の
  一致: **abs Δ ≤ max(1e-5·|FD|, 1e-6)**。
- **[S1] block 3×3 χ 行列の dense 照合 (器械)**: (χ_xx, χ_yy, χ_zz と全交差)(qŷ),
  N=8, j∈{1,2}, m∈{0,0.5}: **abs ≤ 1e-9**。[S1b] dense 側で χ_V(cos モード) =
  (1/8)χ_complex の写像恒等式: **abs ≤ 1e-10**。
- **[S2] tadpole の解析照合 (器械)**: e_i の block 読みと解析 k 和
  (−cos²k_i/E⁺ の BZ 和, twist-x 格子) の一致 N∈{8,32}: **abs ≤ 1e-10**;
  e_y = e_z: **相対 ≤ 1e-13**。
- **[S3] 背景停留 (要件 0)**: Λ = −e_y で dense FD |dΓ_ren/dε| (uni-yy, uni-zz,
  N=8): **abs ≤ 1e-6**; x 残差 Δ_x(N) = |e_x/e_y − 1| が N∈{8,16,32,64} で
  **単調減少 (m=0.5)**、m=0 は報告のみ。
- **[S4] 接触項の q 非依存の器械認証**: |⟨T_i(2qŷ)⟩|/V (i∈{x,y,z}, j∈{1,2,3},
  N=32, m∈{0,0.5}): **≤ 1e-12**。
- **[S5] χ_D の完全核成分としての再現 (回帰)**: k̂ の D=(xx−zz)/√2 チャネル
  q² 係数 ×(−4) が v26.5 公表 c₁^(2) (−0.01935 [m=0] / −0.01444 [m=0.5]) と
  **±0.0004** (N=64, v26.3 と同一の窓プロトコル)。
- **[S6] 縦チャネル汚染 (主結果)**: 連続極限の diffeo 不変性は縦 (longitudinal =
  yy@qŷ) チャネルの繰り込み後核の恒等消滅を要求する。lattice scheme の縦 q² 係数
  c₁^k[L] ≠ 0 は**純粋な非共変汚染**である。R(m) = |c₁^k[L]|/|c₁^k[D]| を測り、
  branch (a) R ∈ [1/3, 3] (両質量) — 汚染は spin-2 と同桁 = bare c₁ は非共変
  artifact に支配される / (b) R < 1/3 (近似的共変性) / (c) R > 3。
  分解能条件: |c₁^k[L]| > 3×窓系統。
  (v26.4 の実測 c₁[yy](ŷ) は既知なので R の値自体は新情報ではない — 新しいのは
  「接触完備核の中で縦チャネルが担う意味」の確定と、これを繰り込み条件の個数
  (少なくとも Λ・q⁰ 核・q² 縦の 3 つ) に翻訳すること。)
- **[S7] スカラー混合 K_SL (新測定)**: S=(xx+zz)/√2 と L=yy の混合核
  k̂_SL(q) の q² 係数と質量走行の記録 (branch 記録 — P0s×P0w transfer の実測)。
- **[S8] uniform 連続性 (端から端)**: block エンジンの一様変形 E(ε) (cky 等の
  rescale) の FD Hessian vs k̂(0) 組立 (N=32): **相対 ≤ 1e-6**; m=0.5 で
  |k̂_AB(q_min) − k̂_AB(0)| ≤ 0.05·max|k̂(0)|、m=0 は報告のみ。
- **[S9] 変異検出 (破壊層)**: (i) 接触項 (3/4) 落とし → S0 が **> 1e-3** で検出 /
  (ii) Λ counterterm 落とし → S3 が検出 / (iii) dense T_zz 折返しスワップ落とし
  (v26.5 と同型) → S1 が **> 1e-4** で検出。

**判定分岐 (凍結)**: (a) S0–S5, S8 PASS → 完全核の器械認証が成立、S6 の R が
主結果 (branch a/b/c は記録) / (b) S0 か S1 FAIL → 接触項導出または交差器械の誤り
(dense が真・物理値非公表で再設計) / (c) S5 か S8 FAIL → 規格化・窓の再設計。

## 5. 禁止暗黙変換 (主張生成の型制約 — 恒久)

以下の左辺から右辺への昇格は、対応する証明 (Lean) または器械 certificate
(dense 照合・FD 完全性・繰り込み条件の明示) を伴わない限り**禁止**する。
文書・claims.yml・論文のいずれでも適用する。

```text
bond modulation        ≠ full vierbein coupling
T00 conservation       ≠ conserved Tμν
plus component         ≠ complete spin-2
mass dependence        ≠ regulator dependence
inverse susceptibility ≠ graviton propagator
nonzero q² coefficient ≠ Newton constant
```

- 「mass dependence ≠ regulator dependence」の意味: v26.3–26.5 の m 走行は
  matter mass/threshold 依存の実証であり、**UV 離散化を変えた試験ではない**。
  regulator 依存の直接試験は v26.8 (独立離散化・a→0) まで主張しない。
  同様に N=32→64 は q_min を下げる有限体積極限であり、固定 (m_phys, q_phys,
  L_phys) の a→0 連続極限ではない。

## 6. 要件表の改訂 (凍結)

- **要件 0 (新設): 背景の停留性** — tadpole 減算後に平坦背景が
  δΓ/δg|_η = 0 を満たすこと。停留でない背景の周りの Hessian 反転は物理的
  伝播関数ではない。
- **要件 10 (書き換え)**: 「独立離散化に対し、**同じ繰り込み条件の下で非局所
  form factor と spectral density が同一の continuum limit に収束すること**。
  bare c₁ の一致は要求しない。」
- 「graviton propagator」の語の使用条件 (すべて満たすまで禁止):
  (1) h_μν を積分する動的変数として定義 (2) 背景停留 (3) cosmological/tadpole
  繰り込み (4) gauge redundancy と gauge fixing の定義 (5) TT residue 正・
  余分な scalar/ghost なし。
- composite-graviton 路線 (metric が microscopic state から創発) を採る場合は、
  **Weinberg–Witten no-go のどの仮定 (Lorentz covariance / 局所性 / covariant
  conserved stress tensor / 通常の粒子状態) を破るのかを計算開始前に明記**する。
  「ネットワークだから回避」は不可。

## 7. v26.7 の登録 (動的 spectral response)

- 主成果物は broadened plot ではなく**有限サイズ spectral measure**
  ρ_AB^(N)(s) = Σ_n Z_n,AB^(N) δ(s − s_n^(N)) (s = ΔE², 一粒子—正孔対の離散和)。
- **pole 判定 (凍結)**: massless pole の主張は s_n^(N) → 0 **かつ**
  Z_n^(N) → Z_* > 0 (residue の体積 scaling) と連続 spectral weight からの分離を
  要する。「分母が零に近い」は証拠にならない。
- 事前登録の自由場 benchmark: (i) massless spin-2 の δ(s=0) は**存在しない**
  (raw χ に pole なし — v26.5 の登録予想を踏襲) (ii) massive 系は二粒子
  threshold s ≥ (2m_gap)² (iii) 物理 TT チャネルの spectral density 非負
  (iv) Euclidean と retarded response の spectral 表現整合
  (v) plus/cross の連続極限縮退 (cross 器械認証後)。
- **副実験 (凍結・モデル交換禁止)**: v26.3 の登録予言「相互作用で χ_00 の q⁴ 保護
  が破れ q² に戻る」の検査。仕様: 1+1D staggered 鎖 (twist 境界, 半充填,
  many-body 厳密対角化/CG)、相互作用 V Σ_⟨ij⟩ n_i n_j、**V ∈ {0, 0.5, 1.0} 固定・
  符号 + (斥力) 固定**、サイズ列 N ∈ {10, 12, 14} (可能なら 16)、fit 規則 =
  最小 2 運動量点の log–log 傾き p。**判定: p(V=1.0, 最大 N) < 3.0 → 的中 /
  ≥ 3.0 → 外れ (外れてもモデル・V・fit 規則を交換しない)**。バイナリ名
  `v267b_q4break` も凍結。predictions.yml に PRED として登録する。

## 8. v26.8 の登録 (continuum universality — 経路 B の最重要 falsifier)

- 明示的格子間隔 a を導入し、固定 (m_phys, q_phys, L_phys) で a→0
  (N = L_phys/a → ∞)。staggered の taste 数 (3+1D で 2) を明示。
- 2-taste Dirac の解析的 1-loop 重力分極 form factor と照合: **q⁴ ln q² の係数**
  を抽出し、解析値と比較。独立離散化 (Wilson fermion) で再計算し、**同じ
  繰り込み条件を課した後の form factor を比較**する。
- **一致しなければ source・contact 項・taste 規格化・continuum scaling の
  いずれかが誤っている。その状態で 1/Π や graviton pole に進んではならない。**

## 9. v26.9 の登録 (動的 metric の分岐 — 実施前判断)

三つを混同しない: (1) 外部背景 metric への matter response (v26.6–26.8 の全て) /
(2) 補助 metric を導入し matter loop で kinetic term を誘導する EFT /
(3) metric 自体が microscopic state の複合励起として創発する理論。
v26.x で到達可能なのは (2) まで。(2) に進む場合も「重力場の kinetic term が
matter loop により誘導された」以上を主張しない ((3) には §6 の Weinberg–Witten
明記が要る)。massless spin-2 が本当に存在する場合、普遍結合 (要件 8) は美点では
なく整合性条件である (soft-graviton 定理)。

## 10. QRN-Core v1 との接続 (型制約)

Core v1 は新しい大規模シミュレーションではなく**型付き仕様**から始める:

```text
AtlasResult → ExternalMetricResponse → CovariantRenormalizedKernel
  → DynamicalMetric → MasslessTTMode → UniversalCoupling → BackreactingGravity
```

各矢印は proof certificate (Lean 証明オブジェクトまたは器械 certificate) なしには
生成不能とする。§5 の禁止暗黙変換はこの型系の公理である。Core v1 の最小成果物:
{H, H(θ), C_gauge, R_b, O_matter, O_geometry, B_backreaction} を単一の有限パラメタ
集合 θ から定義し、**同じ θ で少なくとも一つの未使用 cross-domain observable を
予測する**こと。

## 11. 論文・外部検証の優先順位 (登録)

1. v25.2-freeze の GitHub Release 化と Zenodo DOI (タグ・manifest は存在 —
   Release/DOI は手動手順)。
2. bounded anomaly-search と modular-BW は v26.x を待たず投稿する (独立査読単位)。
3. 第三者再実装を最優先 milestone とする (現状: 凍結新規予言の的中 0・独立外部
   再現 0)。
4. v26 誘導重力論文は **v26.8 の後** (それまでは「格子 stress-response
   instrumentation」の技術報告に留まる)。
