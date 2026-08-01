# D2-R — 応答 lane (protocol v32.1, end_to_end)

**geometry 能力の blocker 解除候補はこの lane のみ。**

曲率則 (n̈_j⁺(0) − n̈_j⁻(0))/(4ε) = ‖P_j h P_i‖_F² は任意の一体生成子に対する
**有限行列恒等式**であり (docs/uft-v31.2.md)、静的 B3 の長距離相関混入
([d2-static.md](d2-static.md) の反例) を受けない。D2-R はこれを単なる恒等式確認
(Unit D4) で終わらせず、**トポロジーまで接続した end-to-end** の再現単位にする。

## プロトコル (言語中立 — 本リポジトリの seed・コードを使わない)

1. 連結な一体生成子 h (エルミート; 重みつき可) を独立に生成する。
2. probe 対 C± = I/2 ± εP_i (0 < ε < 1/2) を準備し C(t) = e^{−iht} C± e^{iht} で
   発展させる。
3. **ノード密度時系列 n_j(t) = (C(t))_jj のみ**を readout に渡す (h は readout に不可視)。
4. 曲率 (n̈_j⁺(0) − n̈_j⁻(0))/(4ε) から ŵ_ij ≈ |h_ij|² を復元する
   (推奨: 5 点 stencil + Richardson, dt ≤ 0.02/‖h‖₁)。
5. 凍結 gap 則 ([d2-static.md](d2-static.md) と同一規則) で支持を得る。
6. clique complex (4-単体まで) → Z2 homology β₀..β₃ (**β₃ には ∂₄ が必須**) →
   vertex link 多様体性 (S²/D²/特異)。

## 必須セル (6 — schema の cell enum と同一)

| cell | 期待 |
|---|---|
| `kuhn_t3_L4plus` | Kuhn T³ (**L ≥ 4**): 支持 欠0余0・β = (1,3,3,1)・link closed |
| `cell16_s3` | 16-cell (S³): β = (1,0,0,1)・link closed |
| `solid_kuhn_ball` | 中実 Kuhn 3-ball: β = (1,0,0,0)・link boundary |
| `t3_L3_flag_break_negative` | T³ L=3: 周期軸線の巻き付き 3-クリークで clique complex ≠ 三角化 — **負制御** (位相一致を主張しないこと) |
| `independent_sparse_weighted` | 独立生成の疎な重みつきグラフ: 支持 欠0余0・ŵ_ij = \|h_ij\|² (rel ≤ 1e-4) |
| `high_noise_abstain` | 測定ノイズ大: **Abstain(InsufficientObservation)** — 強制回答は fail |

## ノイズ下の裁定 (v32.1 で凍結 — 2 段)

- **重み段** (凍結決定規則 4): 誤差見積り noise_error_bound = σ·17√6/(3·dt²·4ε) が
  0.1 を超えたら棄却。
- **支持段** (SupportNoiseCertificate — v321 [E5] で導出・凍結):

  ```text
  noise_error_bound(σ, ‖h‖₁) · √(2·ln(m/10⁻⁶)) ≤ 10⁻³ · max ŵ
  (m = n(n−1) = 閾値判定される順序対の数)
  ```

  不成立なら支持は Abstain(InsufficientObservation)。**重みバーの通過は支持の保証では
  ない** — F}oXO の σ = 1e-9 は重みバー内 (見積り 4.6e-4 ≤ 0.1) なのに、ノイズ最大値が
  gap 則の窓ガード (max·10⁻³) を跨いで余剰辺 1 が実測される (v321 [E5])。非辺の裁定は
  「点推定が閾値以下」ではなく、ノイズ最大値の信頼上界がガードを跨がないことで行う。

## 解除されるもの・されないもの

成功 (NO_SHARED_CODE.md の 6 条件 + matching claim scope) で解除されるのは
**spatial_topology_given_factorization の独立性 blocker のみ**。E3 (因子分解の選定)・
PRED-019・自然の観測量・bridge law 登録は別問題のまま残る。D2-R が対象とする能力は
現状 GivenNodeFactorization を入力に取る (readout は因子分解を選定できない — v31.4)。

## 報告

[unit-d-report.schema.json](unit-d-report.schema.json) の `units.D2R` に従う
(`claimed_capabilities` は `["spatial_topology_given_factorization"]` または `[]`
のみ適合 — 他の能力名は無効な能力昇格として不適合)。
失敗・不一致も同じ形式で提出する。
