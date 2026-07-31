# QRN v31.2 — LocalBiasCommutatorLaw: 局所バイアス応答の短時間曲率則

**Version**: v31.2 (2026-07-31)
**Sim**: `sim/src/bin/v312_localbias_law.rs` → `results/v312_localbias_law.txt` (9 検査 PASS)
**Lean**: `proofs/LocalBias.lean` (6 定理)
**位置づけ**: PROMPT/12 第三版 — **不変ノルム核の第一候補**。B3 (静的共分散) と B4
(到着時刻) を「operator block の不変ノルム」として統一する基礎量。

---

## 1. 法則

既知ノード因子分解の射影 P_i に対し probe 対 **C_i^± = I/2 ± εP_i** (0 < ε < 1/2)
を準備し、一体生成子 h の下で C(t) = e^{−iht}C(0)e^{iht} と発展させると、j ≠ i で:

```text
(1) commutator 則:  P_j (Ċ⁺ − Ċ⁻)(0) P_i = −2iε P_j h P_i
(2) 密度曲率則:     (n̈_j⁺(0) − n̈_j⁻(0)) / (4ε) = ‖P_j h P_i‖_F²
```

導出の核は P_jP_i = 0 (commutator の第 2 項が消える) と
−Tr(P_j[h,[h,P_i]]) = 2‖P_jhP_i‖_F² — **有限行列の恒等式**であり、線形応答近似では
ない (n̈ は C(0) に線形なので ± 差は εP_i 部分を厳密に単離する — [L7] で ε = 0.05 と
0.45 の測定一致 3.4e-11 を機械確認)。

## 2. なぜこれが第三十一期の主役か

1. **ノード内基底 (gauge) 不変** — h_ji ↦ U_j h_ji U_i† で Frobenius ノルム不変
   (Lean 定理 6 + [L3] 乱択 U(2) で 9.6e-12)。
2. **密度だけで読める** — full tomography 不要。5 点 stencil + Richardson で
   rel 4.0e-9 ([L1], dt = 0.02, 改善比 ~1e5)。readout は**時系列のみ**を受け取る
   (h は署名レベルで不可視)。
3. **純粋状態の rank 欠損に依存しない** — equilibrium 状態の logit を使わないため、
   P6/693 の sign(A) no-go と衝突しない。**[L6]: 静的 projector でゲージ同値だった
   P6/693 を、密度曲率則が隣接行列ごと復元して分離 (誤差 6.1e-8・min-perm 分離
   1.000)** — v29.5 の「静的単独不可・応答併用可」が到着時刻の経験則から**恒等式**に
   昇格した。
4. **臨界境界増強を回避** — 長距離 equilibrium 相関を使わず、生成子の短時間局所
   応答を直接読む (probe は準備状態 — 系の熱的/projector/臨界の別と独立)。

## 3. 能力の分離 (観測契約 hierarchy の機械実例, [L5])

block (0→1) を SO(2) 回転した h₂ (Frobenius 重み不変・block は変化) に対し:

| 観測契約 | h vs h₂ | 結論 |
|---|---|---|
| LocalBiasDensityResponse | 重み差 **厳密 0** (測定 6.3e-11) | 識別不能 — 契約準位の情報損失 |
| CoherentLocalResponse | block 差 0.300 | 識別可能 (gauge 共変 [L3]: U_jBU_i† 1.9e-14) |
| ArrivalTimeResponse (参考) | 0.09 vs 0.09 | 到着時刻も識別不能 (Frobenius 重みの下流) |

**LocalBiasDensityResponse < CoherentLocalResponse は真の階層** — v31.3 の観測予算
hierarchy の最初の確定点。

## 4. 検査一覧 (9 PASS)

| 検査 | 結果 |
|---|---|
| [L0] 厳密代数 oracle | 複素エルミート h (4 ノード × 2 軌道) × ε 3 値: commutator 残差 7.9e-17・曲率残差 2.2e-16 |
| [L1] 密度測定 lane | rel 4.0e-9 (時系列のみ)・Richardson 改善 1.0e5 倍 |
| [L2] coherent 測定 lane | block entry 1.4e-8 |
| [L3] gauge | 密度重み不変 9.6e-12・coherent 共変 1.9e-14 (U(2) 乱択 10) |
| [L4] ノード置換 | 9.7e-12 |
| [L5] 能力分離 | §3 |
| [L6] P6/693 | 隣接復元 6.1e-8・分離 1.000 |
| [L7] ε 非依存 | 3.4e-11 (0.05 vs 0.45) |
| [L8] Lean | 6 定理 |

## 5. Lean 形式化 (proofs/LocalBias.lean, 6 定理)

1. `commutator_block_identity` — P₂(hP₁−P₁h)P₁ = P₂hP₁ (射影直交性が第 2 項を消す)
2. `trace_frobenius_identity` — Tr(P₂hP₁h) = ‖P₂hP₁‖_F²
3. `curvature_frobenius_identity` — −Tr(P₂[h,[h,P₁]]) = 2‖P₂hP₁‖_F²
4. `probe_difference_isolation` — (I+2εP)−(I−2εP) = 4εP (状態非依存部の厳密相殺)
5. `multiorbital_trace_frobenius` — 2 軌道ノード版 (対角 block は寄与しない)
6. `gauge_invariance_frobenius` — G(U) = G(V) = I ⇒ ‖UBVᵀ‖² = ‖B‖²

格子 native_decide (Projector.lean 規約)。スコープ明示: 実対称・1/2 次元ノードの
格子恒等 — 複素エルミート・一般 d は数値側 [L0]–[L3] が乱択検査。

## 6. 正直な残高

- **観測資源は「局所 probe 準備 + 密度時系列」** — 状態だけからの読み出しではない。
  証明書型は SpatialMetric × GaussianGibbsFullRank (probe) × LocalBiasDensityResponse
  × GivenNodeFactorization。因子分解は依然入力 (E3 は未解決)。
- 読めるのは **h の block Frobenius 重み** (計量そのものではない) — 計量への昇格は
  B3 との整合 (v31.3) と VR persistence (v31.6) の後。
- probe C± = I/2 ± εP_i の準備は理想化 (瞬時準備・デコヒーレンスなし)。
- 到着時刻の h/h₂ 一致 ([L5] 参考) は「ArrivalTime が Frobenius 重みの下流」を示唆
  するがまだ 1 例 — B4 の圧縮性の系統評価は v31.3。
- HOLD-7 の holdout 検証はまだ (本版は資格・設計区画)。bridge_candidates.yml への
  凍結登録は v32.0 の HOLD-7 凍結時に行う。

## 7. 開発記録

- [L5] の到着時刻が h/h₂ で完全一致 (0.09 = 0.09) したのは想定外の副産物 — 到着
  時刻チャネル (B4 系) が「重みの下流」であることの最初の直接証拠。密度曲率則は
  その上流の基礎量を直接読む。
- herm_eig は 2n×2n 実対称埋め込み + 複素 Gram–Schmidt (縮退があっても射影残差で
  正しく拾う)。評価は h の固有系での厳密位相回転 — 時間積分誤差なし (誤差は
  有限差分のみ = Richardson で 1e5 倍改善が効く理由)。
