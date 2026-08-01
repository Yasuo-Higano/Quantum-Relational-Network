# 単位 D — Geometry Identifiability (第三十一期, PROMPT/12 Track X)

gauge 系の単位 A/B/C は **geometry 能力の blocker を解除しない** (replications.yml)。
geometry 系能力 (spatial_topology_given_factorization 等) の独立外部再現は本単位で行う。
共有コード禁止・凍結プロトコル・失敗も公開 — の 6 条件は SPEC.md / NO_SHARED_CODE.md と
同一。**同一 AI による再実装は independence を満たさない (algorithmic diversity どまり)。**

## D1 — 凍結 C 行列からの kernel/readout 数値再現 (replication_level: numerical_only)

入力: `INPUTS/unit_d1_frozen_c.json` — u693 (単環 6 頂点, 辺表つき), β = 1,
h = −隣接行列, C = (I + e^{βh})⁻¹ の参照値 (1e-9 丸め)。

再現すべき量 (許容誤差 1e-6):
1. C 行列の再計算が参照値と一致する。
2. B3 核 |C_ij|² のスケールガード付き gap 支持 (max·10⁻³ 窓・窓内有意段差 ≥ ln 3、
   なければ窓 = 単一クラスタ) = 辺表そのもの (欠0余0)。
3. 大域 logit K(C) = log[(I−C)C⁻¹] = β·(−A) (off-diagonal 誤差 ≤ 1e-9)。

## D2 — 独立モデル生成・独立 state solver・独立 readout の end-to-end (replication_level: end_to_end)

**geometry 能力の blocker を解除できるのはこの水準のみ。**

**プロトコル改訂 (v32.1)**: 初版 D2 (v31.7,「任意の連結グラフ」) は外部走行 0 件の
時点でプロトコル反例が見つかり `superseded_before_external_run` として版分離された。
熱的相関 C = (I+e^{−A})⁻¹ は A の解析関数なので、連結グラフでは一般に非辺にも歩道
由来の非零相関が乗り、凍結 gap 則は任意グラフの支持を厳密再現しない — 反例 graph6
`F}oXO` (7 頂点 11 辺 → 13 辺の余剰報告)・n=7 連結同型類 853 中 22 で同型故障
(全て余剰のみ・欠落 0)。これは外部再現の失敗ではなく実行前に発見・登録された
プロトコル欠陥である (発見器械: `sim/src/bin/v321_d2_erratum.rs`・supersession 台帳:
[protocols/v32.1/protocol-index.yml](protocols/v32.1/protocol-index.yml)・旧文面の
逐語保存: [protocols/v31.7/d2-v1-superseded.md](protocols/v31.7/d2-v1-superseded.md))。

現行 D2 は 2 lane に分離される。報告契約 = [protocols/v32.1/unit-d-report.schema.json](protocols/v32.1/unit-d-report.schema.json)
(JSON Schema draft 2020-12 — 正直な失敗は適合・能力の水増しは不適合) /
許容誤差と凍結参照値 = [protocols/v32.1/unit-d-tolerances.yml](protocols/v32.1/unit-d-tolerances.yml)。

### D2-S — 静的 B3 lane ([protocols/v32.1/d2-static.md](protocols/v32.1/d2-static.md); StableEstimate scoped — blocker は解除しない)

事前証明された分離マージン (**B3SupportMarginCertificate** — 走行前に真値グラフから
計算でき、certificate 成立 ⟺ 凍結則 exact が n=4..7 全 992 グラフで例外 0) を持つ
グラフ族に対してのみ支持と下流位相の再現を主張する。必須負制御 = `F}oXO` (凍結則の
13 辺誤報告の再現 + certificate 不成立の確認)。certificate 不成立グラフの正答は
rejected_no_certificate — 強制回答は fail。

### D2-R — 応答 lane ([protocols/v32.1/d2-response.md](protocols/v32.1/d2-response.md); end_to_end — **blocker 解除の本命**)

独立に選んだ h → 独立 probe 準備 (C± = I/2 ± εP_i) → **密度時系列のみ** → 曲率則で
ŵ = |h_ij|² 復元 (有限行列恒等式 — 静的 lane の反例を受けない) → gap 支持 →
clique complex → Z2 homology (β₀..β₃ — **β₃ には ∂₄ が必須**) → vertex link。
必須セル 6: Kuhn T³ (**L ≥ 4**)・16-cell (S³)・中実 Kuhn 3-ball・T³ L=3 の flag
破れ負制御・独立生成 sparse weighted graph・高ノイズでの Abstain(InsufficientObservation)。
ノイズ下の支持は SupportNoiseCertificate (ガード比 — 重みバー通過は支持の保証では
ない) で裁定する。成功で解除されるのは **spatial_topology_given_factorization の
独立性 blocker のみ** — E3 因子分解・PRED-019・bridge law 登録は別問題のまま。

## D3 — 負定理の独立確認 (replication_level: negative_theorem)

1. **P6/693 の projector ゲージ同値**: P6 (path 6) と u693 の半充填基底状態射影 C は、
   Z2 ゲージ (z_i = ±1) × 頂点置換で min ‖C_P6 − z π C_693 π z‖∞ ≤ 2e-13 —
   full C の観測でも原理的に識別不能。有限 β では大域逆 K の min-perm‖ΔK‖∞ = β で分離。
2. **projector の資格棄却**: 半充填射影 C (固有値 {0,1}) は clamp なし full-rank 資格
   (0 < λ < 1) を通らない — 正しい読み出しは同値類/棄却であって強制回答ではない。
3. **factorization 非一意性**: 同一の熱的 C (重み場 ring12) に対し、site 基底 →
   ring 支持 / eigenmode 基底 → 支持なし (モード間相関 = 厳密 0) / pair 回転基底 →
   別の支持 — 静的 one-body 相関だけからの一意選定は不可能。疎性基準は mode 基底
   (幾何自明) を選ぶ。

## D4 — LocalBiasCurvatureLaw の独立実装 (replication_level: response_law)

probe 対 C± = I/2 ± εP_i (0 < ε < 1/2) を C(t) = e^{−iht}C±e^{iht} で発展させ、
ノード密度 n_j(t) = Tr(P_j C(t)) の**時系列のみ**から

```text
(n̈_j⁺(0) − n̈_j⁻(0)) / (4ε) = ‖P_j h P_i‖_F²
```

を検証する (j ≠ i)。期待:
1. 一様重み ring: Ŵ = 隣接行列 (rel ≤ 1e-4, 5 点 stencil + Richardson,
   dt ≤ 0.02/‖h‖₁ 推奨)。
2. ε ∈ {0.05, 0.45} で読みが一致する (恒等式 — 線形応答近似ではない)。
3. block-local unitary で Frobenius 重みが不変。
4. (拡張) spinless t-V (密度対角相互作用): 読みが V に依存しない (厳密転移)。

## 台帳への登録

報告は replications.yml の v31.0 拡張書式で登録する (claim_ids / capabilities /
replication_level / independence_scope / protocol_commit / generator_hash / input_hash)。
対象 claim: D1/D2-R → QRN-BRIDGE-013, QRN-BRIDGE-017 / D2-S → QRN-BRIDGE-021 /
D3 → QRN-BRIDGE-012, QRN-BRIDGE-015 / D4 → QRN-BRIDGE-013, QRN-BRIDGE-016。
**geometry 能力の解除効力を持つのは matching D2-R のみ** — D1/D2-S/D3/D4 は数値再現・
scoped 静的再現・負定理確認・応答法則確認として登録される (それ自体に価値があるが
blocker は残る)。
