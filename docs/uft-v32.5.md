# QRN v32.5 — Interaction hypergraph: H_S 直交分解・相関ホッピング・条件付き probe

**Version**: v32.5 (2026-08-01)
**Sim**: `sim/src/bin/v325_interaction_hypergraph.rs` →
`results/v325_interaction_hypergraph.txt` (7 検査 PASS)
**位置づけ**: PROMPT/13 第五版 — 候補 3 (相関ホッピング) と候補 5 (相互作用次数) を
一つの研究プログラムに統合する。v31.5 の未走査 ([H_int, n_j] ≠ 0 の相関 hopping) を
**「曲率則が壊れるか」ではなく「何が読めるか」**で採点する (v32.4 の応答階層の適用)。

---

## 1. 相互作用ハイパーグラフの正準定義 [G0][G1]

因子分解確定後 (v32.3)、H を**局所 Hilbert–Schmidt 条件期待値**で直交分解する:

```text
H = Σ_{S⊆V} H_S,   H_S = Σ_{T⊆S} (−1)^{|S\T|} E_T(H),   w_S = ‖H_S‖²_F
E_T(H) = (Tr_{T^c} H) ⊗ I_{T^c} / 2^{|T^c|}
```

- **[G0] 二重定義一致**: Möbius 条件期待値 = Pauli 支持射影 (5.6e-16)・直交性
  6.9e-18・完全性 3.9e-16 — E_T は基底自由の定義なので、w_S は**局所演算子基底の
  選び方に依存しない**。
- **[G1] 局所 unitary 不変** (max|Δw_S| = 1.8e-15)・非局所 unitary (DFT₈) は変える
  (max|Δw_S| = 2.83 — 負制御)。|S| = 1: on-site / 2: edge / 3: correlated-hopping 級
  hyperedge / ≥ 4: 高体。

## 2. 中心化分離 [G2] — 「三体項」の半分は二体に住む

V·n₃h₁₂ (相関ホッピング) は n₃ = 1/2 + (2n₃−1)/2 の中心化により
**二体 (V/2)h₁₂ ⊕ 真の三体 −(V/2)Z₃h₁₂** に正準分離される — w_{12} = w_{123} = V²
(等重み, 機械値 0.36 = V²)。二体部は「半充填平均で dressing された hopping」であり、
模型 H = t(h₁₂+h₂₃) + V n₃h₁₂ + μn₁ の正準重みは
w_{01} = 4(t+V/2)² / w_{12} = 4t² / w_{012} = V² / w_{0} = 2μ² (全て厳密一致)。

## 3. 三つの観測 lane (何が読めるか) [G3][G4][G5]

| lane | 読み | 機械値 |
|---|---|---|
| **1. density curvature** (非条件付き) | **遷移率和則**: K_uncond(j←i) = Σ_{S⊇{i,j}} w_S/4 — 対を含む hyperedge 重みの**和** (条件付き遷移率の Gram 核。「破れ」ではない) | K(2←1) = 1.300000 = (w₀₁+w₀₁₂)/4 厳密 |
| **2. conditional density probes** (Möbius/Boolean-Fourier) | 補助ノードの占有条件付けで次数を分離: K(v) = \|t + vV\|² 厳密・**K(1)−K(0) は {i,j,k} hyperedge 検出器** (V=0 で厳密 0 = 負制御)・混合恒等式 K_uncond = (K(0)+K(1))/2 | K(0) = 0.64 = t²・K(1) = 1.96 = (t+V)² |
| **3. coherent parity-even probes** | 一階応答で符号・位相を回復: t ↔ −t は密度曲率に**厳密不可視** — coherent 電流 R⁽¹⁾ は奇 (+0.8/−0.8)。密度 lane 単独の正答は**符号同値類** (強制回答しない) | 密度差 0.0・R⁽¹⁾ = ±0.8000 |

## 4. 位置づけ

- v31.5「曲率則の厳密転移 (密度対角 V)」+ 本版「相関 hopping では hyperedge 和を
  読む」で、**相互作用次数 × 観測契約の相図**が閉じた: 密度二階 = w の対和 /
  条件付き密度 = 次数分離 / coherent 一階 = 符号・位相 / (大域 oracle = Gaussian
  のみ — v31.5)。HOLD-8 の相互作用セル (quadratic/t-V/相関/pair/三体混合/H↔−H/
  位相対) はこの相図で採点する。
- w_S の定義は Zanardi 型の operational tensor product の上の標準的な直交分解 —
  新規性は**因子分解 (v32.3) と hypergraph の同時識別・観測契約別の読み・棄却型**
  の側に置く (PROMPT/13 §7)。

## 5. 正直な残高

- toy (3 モード JW, dim 8)。pair hopping・4 体・複数 hyperedge の交差項は HOLD-8
  セルで走査する (枠は本版で確定)。
- 遷移率和則の係数 1/4 は JW bond 演算子の規約 (‖h_bond‖²_F = 4) — 一般格子への
  規約非依存の正規化は d2-response 型の凍結時に固定する。
- 条件付き probe は補助ノードの占有固定 (射影準備) を仮定 — 準備誤差・部分条件付け
  のノイズ裁定は SupportNoiseCertificate (v32.1) の変成で扱う。
