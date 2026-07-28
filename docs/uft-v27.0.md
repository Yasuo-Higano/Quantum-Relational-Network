# QRN v27.0 — full 4D kernel と dynamic metric fork (spec §13)

**位置づけ**: 第二十七期の最終 arc。事前登録 = spec §13 (b46bce3, 改訂 3 —
実装前凍結): (13.1) lapse-shift 型 temporal 結合 (N·ε + Nⁱ·πᵢ_Bel,
λ = −1/8)、(13.2) full 4D kernel の厳密 Ward と連続 universality、
(13.3) **dynamic metric fork の判断プロトコル — 1/Π 解禁 3 条件・三分岐・
Weinberg–Witten 破れ仮定の事前明記義務**、(13.4) ユニット計画 (0/A/B/C/D —
D = 期統合)。前提: v26.9 arc 完了 (動的 Ward 144 恒等式・Belinfante 完成・
massless σ = ρ₂P₂ 完全崩壊)。

## A. v270a_kernel4d — full 4D kernel の厳密 Ward (5 検査 PASS)

**構成** (v26.9-B の恒等式を核の定義に転化):

  k̂^{0ν,B}(iq₀,q) := C_{T⁰ν_Bel, B}(iq₀)
  k̂^{yν,B}(iq₀,q) := C_{J_ν, B}(iq₀) − ⟨[T⁰ν_Bel(q), B(−q)]⟩ / q̂

(T⁰ν_Bel = λ = −1/8 改良込み Belinfante 密度・J_ν = 厳密流束・第 2 項 =
接触完備化)。恒等式 iq₀C − q̂C_J = −⟨[A,B]⟩ により **iq₀k̂^{0ν} − q̂k̂^{yν} = 0
がカットオフ有限のまま構成的に厳密** — v26.6 静的核 (17 検査) の 4D 拡張。

| 検査 | 結果 |
|---|---|
| A0 接触完備化 X_ν = −⟨[T⁰ν,B]⟩/q̂ の q → 0 正則性 | 成長比 1.000 (非発散) |
| **A1 厳密 4D Ward (4 行 × 10 列 × 2 周波数, 12³)** | **6.8e-14** |
| A2 Onsager 対称性 (反対称改良項の微小破れ込み) | 7.5e-9 |
| **A3 流束行 → Belinfante stress 行 (全 4 行)** | O(ε²): 例 ν=x 1.5e-2 → 9.5e-4 |
| A4 変異 (λ → 0) | ν=x 行 1.30 停留 (正版の 1370 倍) |

A3 が核の**対称テンソル整合性**: Ward-厳密な流束行が連続極限で Belinfante
stress 行に一致 → 核は連続極限で「対称・保存・Ward-厳密」の三条件を同時に
満たす。開発記録: (i) A0 run1 は「定数収束」を誤要求 — 交換子が q³ で消える
対 (X ∝ q² → 0) も正則 (正しい判定 = 非発散)。(ii) Onsager は反対称改良頂点
により厳密には成立しない (破れ ~1e-8 — 高次微小)。

**型名**: spec §13.2-A の Ward 要件は充足。ただし
`FullGravitationalVacuumPolarization` の発行は §13.3 の全条件 (v27.0-B の
連続 universality + fork 予言の事前登録) 通過まで保留。**1/Π 禁止維持**。

## 1. 残り (v27.0 arc)

- **v27.0-B**: k̂ の繰り込み後 form factor の連続極限 — 10×10 全チャネルの
  2 関数 (ρ₂/ρ₀) 崩壊と q⁴ln q² 係数の oracle 照合 (外挿モデルは §12.9 規律)。
- **v27.0-C**: fork 判定の執行 (1/k̂ の pole 構造 — 予言: 自由場 = pole なし
  = 分岐 (b) 既定)。
- **v27.0-D**: 期統合 (第二十七期統合文書・全スイート儀式・QRN-Core v1 着手・
  CLAUDE.md 到達点更新)。

## 2. 成果物

A: `sim/src/bin/v270a_kernel4d.rs` / `results/v270a_kernel4d.txt`
(5 検査 PASS) / `results/v270a_kernel4d.json`。claims: QRN-GRAV-057。
