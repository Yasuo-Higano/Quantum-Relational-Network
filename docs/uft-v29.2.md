# QRN v29.2 — 意味論再基礎化: 候補の再登録・二層分離・HOLD-5 の凍結

**Version**: v29.2 (第三十期 — PROMPT/11 ロードマップ)
**Date**: 2026-07-30
**Sim**: なし (物理 run なし — 定義・凍結・文書のみ。lib にも触れない)
**位置づけ**: v29.1 (Integrity Erratum) の帰結を凍結記録に反映し、次段 (v29.3
合成 bridge / v29.4 定量計量) の**実行前凍結**を完了する。holdout を 3 回観測した
問題への応答 — 新しい検証力は本版で凍結した HOLD-5 から。

---

## 1. bridge_candidates.yml 追補 (append-only — 原文は as-run 不変)

1. **正名の記録**: B3 → B3-COV (density connected-covariance) / B4 →
   B4-DENSITY-FRONT (空間核×因果チャネルの**合成** readout — 単独候補ではない)。
2. **未試験候補の定義凍結** (実装・較正は v29.3 以降):
   - **B5-QFI**: w(i,j) = ‖ρ_ij − ρ_i⊗ρ_j‖_BKM — BKM 内積は固有分解の閉形式
     c(λ_a, λ_b) = (λ_a − λ_b)/(ln λ_a − ln λ_b) で厳密評価 (事前登録 B3 の本来形)。
   - **B6-COMMUTATOR**: t*(i,j) = 最初の t で |⟨[n_i(t), n_j]⟩| ≥ ε。Gaussian
     運動学変種 |{c_i(t), c†_j}| (状態非依存の運動学的円錐) を併記し区別
     (事前登録 B4 の本来形)。
3. **v29.3 の採点原則**: S (B1, B2, B3-COV, B5-QFI) × C (DENSITY-FRONT,
   COMMUTATOR) の**全組合せを独立セル**で採点。距離供給のフォールバック禁止・
   棄却候補の得点流用禁止。

## 2. HOLD-5 生成器の凍結 (走行前 — 生成分布・分割・seed 規則)

- **系族**: 不均一 1D Dirac — 滑らかな正値速度プロファイル
  v(x) = 1 + Σ A_k exp(−((x−c_k)/w_k)²)、3 regulator (単鎖 / 細分格子 a/2 /
  Wilson 2 成分)、m_phys·a ∈ {0, 0.05}、リング/開放端。
- **seed 分離**: train-0..3 (設計可) / val-4..7 (中間確認 1 回) / hold-8..15
  (**走行版まで生成すらしない**)。
- **コミットメント**: instance seed = sha256(SECRET + ":" + id) の先頭 8 バイト。
  **sha256(SECRET) = `cfeb1959f3ba63b17c8ea7d1261f0f24570d33f7f8a2b3956a3df761bd301ab2`**
  を公開し、SECRET は holdout 走行版の結果公開と同時に開示 — 第三者が
  コミットメントと全生成系列を検証できる。
- **v29.4 の採点法の予告**: 静的核の再構成距離が未使用の動的到着時刻
  τ(x,y) ≈ Σ 1/v を無調整予言。距離一致は Δ∞ = inf_α max |ln(d₁/(α d₂))|
  (スケールのみ自由)。バー値は走行版の凍結時に導出して事前登録。

## 3. 二層分離 (FactorizationBridge / GeometryBridge)

`global algebra/state → [FactorizationBridge] → local subalgebras →
[GeometryBridge] → adjacency/metric/causal order`。現行の全 bridge 成果は
**分解を入力に取る** (NodeState はモードをノードへ群化して受け取る — 隠したのは
座標・隣接・ラベル)。FactorizationBridge は未定義の別課題として schema に登録
(status: undefined)。spec §2 の「分解は読み出し」は kinematics の**設計目標**で
あり現行成果の記述ではない、と明確化した。

## 4. 変更一覧

- bridge_candidates.yml: 追補 §v29.2 (上記 1–2)。
- core.schema.yml: 概念 3 追加 (FactorizationBridge / BkmFisherKernel /
  CommutatorFrontKernel — 計 70 概念)。
- docs/qrn-terminology.md §4b・docs/qrn-core-v1-spec.md §2/§5: 二層分離。
- claims: QRN-META-036 (C5)。バイナリ追加なし・sim/src 変更なし (儀式不要)。
