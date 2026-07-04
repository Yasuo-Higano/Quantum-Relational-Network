# 論文骨子 2: 幾何からの湯川 — 乱雑係数なしの階層と余剰次元の数のベイズ選択

**状態**: 骨子 (v8.3)。v7.2/v8.1 の論文化。
**対象誌の候補**: JHEP / PRD (現象論) — ベイズ模型選択の方法論も売りになる。

---

## Title (案)

> **Yukawa Hierarchies from Magnetized Tori without Order-One Coefficients:
> Bayesian Selection of the Number of Extra Dimensions**

## Abstract (案, ~150 words)

> Froggatt–Nielsen fits of fermion mass hierarchies involve chosen integer charges
> and marginalized order-one coefficients. We replace both by geometry: generations
> are the exact zero modes of a lattice Dirac operator on a magnetized torus
> (degeneracy = flux = 3), sectors are distinguished by discrete Wilson lines, and
> Yukawa matrices are computed as overlap integrals with no random coefficients.
> A single torus is shown to be intrinsically too shallow (its hierarchy depth is
> set by the flux alone), failing the up-quark hierarchy by two orders of magnitude.
> The factorizable product T²×T² squares the suppression and reproduces five of six
> mass ratios and two of three out-of-sample CKM elements within a factor of five,
> with Bayesian evidence exceeding the anarchic bound by 15 nats. Treating the
> geometry itself as a hypothesis space, evidence from mass ratios alone selects
> three tori over two, while coarse Wilson-line lattices are always preferred.
> Including the CKM elements in the evidence — computed exactly by a triple sum
> over the shared Wilson line — resolves the apparent mass–mixing tension:
> three tori win on all nine observables, with the joint MAP placing eight of
> nine within a factor of five.

## 構成

1. **Introduction** — M0/M1/M2/M3 の階層 (改良方針 §6)。v6.5 のベイズ比較
   (M1: lnB≥23 だが電荷は選択) からの問題設定。
2. **Construction** — 磁束トーラスの厳密ゼロモード (縮退=Q、v2.3)、位置演算子による
   局在化 (一般位相での縮退回避)、Wilson 線によるセクターの住所 (中心が k サイト
   厳密シフト)、重なり積分。**開発記録**: 格子磁気並進が N|2Q でしか閉じない障害と
   Wilson 線への転換 (方法論として価値)。
3. **Single-torus no-go** — 階層深度は磁束のみで決まり格子サイズ不変。到達下限
   ~3×10⁻³、lnZ = −53.8 (アナーキー上界未満)。**原理的な陰性結果**。
4. **T²×T²** — 抑制の 2 乗。質量比 5/6 (m_u/m_t 比 2.2)、out-of-sample CKM 2/3
   (\|V_cb\| 比 1.03)。lnZ = −20.4 (アナーキー +15)。M1 との残差 −8.2 の分解。
5. **Geometry selection** — 模型空間 {トーラス数}×{Wilson 格子}: 質量の証拠は T³
   (+3.0)、格子細分化は常に Occam 罰負け。**質量と混合の緊張とその解消** (v8.1→v9.1):
   質量のみの MAP は CKM を失うが、CKM 込みの証拠 (共有 Wilson 線 K_Q の三重和で厳密
   計算、e セクターは因子化) でも T³ が勝つ (+1.6)。緊張は MAP 点評価と証拠 (Occam
   積分) の混同による人工物 — **方法論の節として独立の価値** (点評価で模型を捨てるな)。
6. **Limitations / Outlook (M3)** — 対角世代対の仮定、同一トーラス、Higgs 簡略化、
   Wilson 線・磁束・トーラス数の起源 (モジュライ安定化)。
7. **Reproducibility** — 全て決定論 (乱数なし)、外部依存なし、v6.5 と同一の尤度で
   の直接比較可能性。

## 図表計画

- 図 1: ゼロモードの局在と Wilson 線シフト (figures/fig_zeromode_wilson.svg — v9.3 作成済み)
- 図 2: 単一 T² vs T²×T² の到達可能な質量比集合 (figures/fig_attainable_ratios.svg — v9.3 作成済み。
  観測の up 対が雲の下に外れるのは m_c/m_t の既知の弱点 (比 14) の正直な可視化)
- 図 3: 幾何模型の lnZ 比較 (figures/fig_geometry_lnz.svg — v9.3 作成済み。両ラベル規約 (v9.2) と
  緊張の解消 (v9.1) を 3 パネルで)
- 表 1: MAP 幾何と 9 量の予測/実測
- 表 2: M0/M1/M2文献/M2geo/M2geo² の証拠とパラメータ数

## 主張と証拠の対応 (claims.yml)

| 主張 | id | 証拠 |
|---|---|---|
| ベイズ比較の枠組 (M0/M1) | QRN-YUK-002 | v65_bayes |
| M2geo (単一 no-go + T²×T²) | QRN-YUK-003 | v72_geomfn |
| 幾何の選択 (T³, 粗い格子) | QRN-YUK-004 | v81_geoselect |
| CKM 込み証拠 (緊張の解消) | QRN-YUK-005 | v91_ckmselect |
| 前提 (FN 機構, 磁束=世代) | QRN-C0-006, QRN-MATTER-001 | v23_matter |
