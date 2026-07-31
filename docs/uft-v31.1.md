# QRN v31.1 — GaussianGibbsInverseOracle: 大域静的逆問題の oracle ceiling

**Version**: v31.1 (2026-07-31)
**Sim**: `sim/src/bin/v311_gibbs_oracle.rs` → `results/v311_gibbs_oracle.txt` (12 検査 PASS)
**Lean**: `proofs/GibbsInverse.lean` (9 定理 — 格子 native_decide + スコープ明示)
**位置づけ**: PROMPT/12 ロードマップ第二版。E1 (global inversion) の天井を機械化する。

---

## 1. 何を確定したか

number-conserving quadratic Gibbs 族 C = (I + e^{β(h−μI)})⁻¹ に対し 0 ≺ C ≺ I なら

```text
K(C) = log[(I−C)C⁻¹] = β(h−μI)
```

**この逆変換は自由フェルミオンの標準的事実であり QRN 固有法則ではない** (Gibbs 状態
からの Hamiltonian learning は独立した既存分野)。本器械は oracle ceiling — 「完全な
大域 C を観測できたとき原理的に何が読めるか」の上界 — として登録する。
**この結果だけから BridgeLawCertificate は登録しない・PRED-019 も登録しない。**

## 2. 確定表 (12 検査)

| 検査 | 結果 |
|---|---|
| [T0] 列挙 | n=4..7 連結グラフ = 6/21/112/853 (OEIS A001349, v29.5 再現)・衝突対 mask 692 (P6)/693 (単環 = C4+ペンダント 2) 実在 |
| [T1] exact lane 全数 | 992 グラフ × β ∈ {0.5,1,2,4} (μ=0.2): **資格 3968/3968・復元誤差は条件数バー 10nε/(βδ(1−δ)) 内 (比 0.463)**・β≤2 は max 1.3e-11。族の最小 δ = 1.7e-11 (β=4 の稠密極限 — 誤差増大は定理の予言どおり) |
| [T2] β 未知 | K(β₁)/β₁ = K(β₂)/β₂ (乖離 1.1e-13, n=5 全数) — **正の大域スケール同値類**。門は UpToPositiveScaleAndShift |
| [T3] μ 未知 | K(μ)−K(0) = −βμI (off-diag 乖離 1.8e-14) — **恒等シフトのみ・空間隣接に無影響** |
| [T4] gauge 共変性 | 多軌道 (3 ノード × 2 軌道) 乱択 20: K'_ij = U_i K_ij U_jᵀ (7e-15)・**特異値/Frobenius/作用素/核ノルムは block-local unitary 不変** (4e-15) |
| [T5] 置換共変性 | K(PhPᵀ)_IJ = K_{π(I)π(J)} (3.4e-15) |
| [T6] logit 2 経路 | **GlobalLogitThenBlock ≠ PairBlockThenLogit の常設反例**: n=2 は一致 (全系=pair)・P4 では global −1.000000 (=βh 厳密) vs pair-B2 −0.927576 — 乖離 0.0724。現行 B2 は環境 renormalize された reduced modular coupling |
| [T7] 条件数定理 | **‖ΔK‖_F ≤ ‖ΔC‖_F/(δ̃(1−δ̃))** (δ̃ = min margin)。微小 30 + 異グラフ大摂動 15 対で max 比 0.779・**整列 rank-1 摂動で飽和比 0.99996 — 上界は最良** |
| [T8] 低温 estimate lane | P4 β=40: exact は RankDeficient で正しく棄却。clamp ε=1e-12 の推定は**飽和則 k = sign(βλ)·min(β\|λ\|, ln((1−ε)/ε))** に従う (乖離 ≤ 1e-4)。**深部モードは sign のみ — β→∞ で sign(A) 類へ連続退化 (P6/693 機構の連続版)** |
| [T9] P6/693 相図 | 下記 §3 |
| [T10] 門 | Gaussian 性証拠なし / Gibbs 出自なし / Wick 残差超え → 全て棄却。証拠あり → Exact |
| [T11] Lean | GibbsInverse.lean 9 定理 (§4) |

## 3. 発見: projector 衝突は「ゲージ同値」である (v29.5 の強化)

v29.5 は P6/693 の**静的カーネル** (readout が見る |C_ij| 準位) の衝突を発見した。
本版はそれを 3 準位に分解して機械確定した:

| 準位 | min-perm 距離 | 意味 |
|---|---|---|
| 符号つき素の C | 0.194 | 行列としては異なる |
| \|C\| (カーネル準位) | 1.1e-15 | v29.5 の衝突の再現 |
| **Z2 ゲージ × 置換** | **1.1e-15** | **full C が厳密にゲージ同値** |

つまり静的衝突は「カーネルが偶然一致する」のではなく、**ノード内位相ゲージ
z_i = ±1 を許すと一体相関行列全体が厳密に一致する** — 与えられた因子分解の下でも、
GlobalOneBodyCorrelation (最強の静的観測) をもってしても、projector lane の P6/693 は
**原理的に識別不能** (ノード内基底は観測で固定できない自由度)。裁定
`EquivalenceClassOnly("sign(A) 同値類")` が正しい読み出しであることの機械的根拠。
有限 β では大域逆が min-perm‖ΔK‖∞ = β で分離する (辺 1 本の差が β で見える)。

**識別可能性相図の第一断面 (E1, GlobalOneBodyCorrelation, GivenNodeFactorization)**:

```text
GaussianGibbsFullRank (δ ≥ 床):  ExactUpToGauge (β,μ 既知) / ExactUpToGlobalScale (β 未知)
低温 (δ < 床):                    StableEstimate — 飽和則で深部は sign のみ
GaussianProjector:                EquivalenceClassOnly (sign(A) 類 — ゲージ同値)
証拠なし (Gaussianity/Gibbs):     Abstain
```

## 4. Lean 形式化 (proofs/GibbsInverse.lean, 9 定理)

同値類構造の等式的骨格を格子 native_decide で確定 (Projector.lean の規約 —
per-variable 次数 ≤ 2 の恒等式を 4 点/変数格子で検証・実数拡張は次数論法):

1–2. **Frobenius–Gram 恒等式** (左/右): ‖UB‖² = ‖B‖² + Gram 偏差項 — 直交性が
     余剰項を消す機構そのもの。
3–4. **直交不変性** (左/両側): G(U) = G(V) = I ⇒ ‖U B Vᵀ‖² = ‖B‖²。
5. **trace 相似不変性**: tr(U B Uᵀ) = tr(B·G(U))。
6. **符号窓**: β ∈ [1,8] で sign(βx) = sign(x)。
7. **交差比不変** (propositional): (βH) のエントリ交差積 = H のそれ。
8. **零支持窓**: βx = 0 ⇔ x = 0 — 隣接の有無は β 不要で読める。
9. **μ シフトの off-diag 無関係性** (propositional, 定義的)。

スコープの明示: block は 2×2・matrix log 自体は未形式化・格子上の保証
(「未証明を証明済みに見せない」— PROMPT/12)。

## 5. 正直な残高

- **これは E1 (大域逆) の天井であって E2 (操作的読み出し) ではない** — 完全な C の
  観測は実際の観測契約では与えられない。E2 の hierarchy は v31.3。
- [T6] が示すとおり、現行 B2 (pair RDM 経路) は大域親生成子と系統的に異なる —
  どちらが「良い」かではなく**別の観測契約の別の量**。
- 条件数定理の δ̃ は数値的に計算した margin — 区間演算での厳密化は未実施。
- Lean は格子上の恒等 + スコープ明示 — 実数全域は次数論法 (コメント)。

## 6. 開発記録

- [T0] の初期期待「693 = C6」は誤り — 単環 = unicyclic (C4 + ペンダント 2)。
  次数列 [1,1,2,2,3,3] で確定。
- [T1] の初期バー 1e-8 は β=4 で破れた — これは**条件数定理の予言どおり**
  (δ = 1.7e-11 → 増幅 ~6e10)。バーを定理由来の per-instance バーに是正し、
  誤差が定理に従うこと自体を検査に変えた (比 0.463)。
- [T9] の素の C 距離 0.194 の解剖が §3 のゲージ同値の発見につながった —
  「衝突はカーネルの偶然」ではなく「ゲージの必然」。
