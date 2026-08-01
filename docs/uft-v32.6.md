# QRN v32.6 — VR exactness: 離散円環 bar 定理・端点規約・H2 persistence

**Version**: v32.6 (2026-08-01)
**Sim**: `sim/src/bin/v326_vr_exactness.rs` → `results/v326_vr_exactness.txt`
(6 検査 PASS)
**位置づけ**: PROMPT/13 第六版 — v32.0-B K3-holes の教訓「定量バーは観測量の応答則
から導出する」の履行。円の VR 複体が S¹, S³, S⁵… を経ること (Adamaszek–Adams) は
既知であり、本版の寄与は**発見ではなく、QRN の離散グラフ測地・有限標本・filtration
規約に対する厳密 bar 端点の導出と機械化**である。

---

## 1. 規約は型で運ぶ [V0]

```text
RipsConvention = DiameterLessThan | DiameterLessOrEqual
BarEndpoint    = Open | Closed        (本器械の bar は [birth, death))
```

整数 filtration では VR_<(r) ≡ VR_≤(r−1) — 規約差は bar 端点をちょうど +1 シフト
する (n = 5..20 全数機械照合)。**L/3 − s をそのまま全離散系の exact formula として
凍結してはいけない** — 規約と離散性が端点を動かす。

## 2. 離散円環 H1 bar 定理 [V1]

> C_n (測地距離・VR_≤・整数 filtration) の H1 bar は **[1, ⌈n/3⌉) のただ 1 本**。
> persistence = ⌈n/3⌉ − 1。

n = 4..30 全数で機械化 (persistence bar = per-step β₁ [GF2 rank] と全 r で一致):

- n ≡ 0 (mod 3): death/n = 1/3 **厳密** (9 件)
- n ≢ 0 (mod 3): death/n は 1/3 に**上から**近づく (超過 ≤ 2/(3n), n = 4 で飽和)

**連続極限の L/3 − s と離散 exact 値の分離**がこれで確定 — 連続式は floor/ceiling
を消した極限であり、有限標本の採点バーには離散式を使う。

## 3. H2 persistence — sparse reduction + column clearing [V2]

境界行列の Z2 sparse 削減 (列 = ソート済み行 index の対称差) + **column clearing**
(次元降順に削減し、対になった creator 列を下位次元でスキップ — Chen–Kerber twist):

- 8 面体 (S²): H2 bar = [1, 2) ちょうど 1 本・H1 なし。
- **wedge-S² 遷移**: n ≡ 0 (mod 3) の r = n/3 で β₂ = n/3 − 1 の短命 bar
  [n/3, n/3+1) — (n, β₂) = (9,2), (12,3), (15,4) 機械確認 (S¹ → S³ の狭間で
  2 球面の wedge を経る — 縮退点の正体)。
- clearing on/off で bar **完全一致**・列削減演算 25% 節約 (13744 → 10310) —
  v31.6 の「大規模 3D で brute force が blocker」への器械側の解。

## 4. アフィン則の導出と K3-holes の retrodiction [V3]

v32.0-B は「persistence はアフィン則 L/3 − birth」を事後分析として記録した。
本版の離散定理はこれを閉形式にする: **pers = ⌈n/3⌉ − 1** (birth = 格子間隔)。

2 穴の理論比: (16, 6) → (⌈16/3⌉−1)/(⌈6/3⌉−1) = 5/1 = **5.0 — HOLD-7 K3-holes の
実測 5.00 を厳密に retrodict** (凍結バー 2.67 は周長比例モデルの誤りだったことの
機構確認)。train 対 (14, 8) → 2.0 (実測 2.20 — 重み場・2D 埋め込みの系統込み)。
**バーと不成立記録は変更しない** (v31.0 絶対禁止 10 —「後から機構が分かったから
成立」にしない)。次の holdout の VR バーはこの離散式から導出する。

## 5. S³ 帯の直接測定 [V4]

VR_≤(C₁₁, 3) = S¹ (1,1,0,0) → VR_≤(C₁₁, 4) = S³ (1,0,0,1) (1/3 < 4/11 < 2/5,
β₃ は ∂₄ 必須 — v31.6 K5 anchor の教訓を再確認)。円の VR が 3 球面帯を経ることの
from-scratch 機械確認。

## 6. 正直な残高

- 対象は一様重み円環 + 8 面体の整数 filtration — **重みつき円環 (実数 filtration)
  の exact bar は未導出** (birth = max 辺重み・death = 最短 1/3 弦の実数化が予想 —
  次の holdout バー導出時に凍結)。
- レンズ空間・大型 3D・高次元は**意図的に後回し** (PROMPT/13 §6: Z2 Betti と link
  だけではレンズ空間は分類できない — torsion 不変量が必要。能力を
  Closed3ManifoldWitness / IntegralHomologySignature / TorsionLinkingSignature /
  LensSpaceEquivalenceClass / HomeomorphismClass に分けてから着手する)。
  現時点では例の拡張であり、能力の上限を変えない (優先順位 4)。
- H0 persistence は報告しない (連結性は既存 pipeline の別檢査)。
