# 論文骨子 1: アノマリー探索 — SM の最小性・一意性・孤立性の再現可能な全数探索

**状態**: 骨子 (v8.3)。改良方針 §16-2「最も強い一本」の論文化。
**対象誌の候補**: SciPost Physics / Computer Physics Communications / PRD。
再現性重視なら ReScience C 系も選択肢 (改良方針 §10)。

---

## Title (案)

> **The Standard Model as the Unique Minimal Anomaly-Free Chiral Spectrum:
> An Exhaustive, Independently Cross-Checked, and Machine-Verified Search**

## Abstract (案, ~150 words)

> We present an exhaustive search for anomaly-free chiral fermion spectra of
> SU(3)×SU(2)×U(1) gauge theories within explicitly bounded domains (representations
> up to the octet/triplet, hypercharges |Y| ≤ 3, up to 8 multiplets), subject to the
> five perturbative anomaly cancellation conditions, the Witten SU(2) anomaly,
> chirality, and charge under all factors. Within every domain scanned, the minimal
> solution is the 15-component Standard Model generation, and it is unique and
> isolated (no solutions with 15 < components < 16 beyond a small set of excluded
> neighbors). The search is implemented three times independently (depth-first,
> linear-elimination grid, meet-in-the-middle), certified by SHA-256 canonical-form
> certificates, and the two core domains are additionally verified as theorems in
> Lean 4 via native_decide with a poisoned-fuel design that makes completeness of
> the enumeration part of the theorem. Extending the gauge group, we prove within
> windows that the only rank-2 chiral extension is B−L (requiring ν_R) and that no
> rank-3 chiral extension exists. All results are reproducible from a
> dependency-free Rust repository with machine-checked claim ledger.

## 構成

1. **Introduction** — 「なぜこの物質内容か」の問い。既存研究 (アノマリー消去の教科書的
   結果、charge quantization 文献) との差分表 (docs/uft-v6.1.md §4 を転用)。
2. **Setup** — 領域と条件の機械可読な定義 (certificates/v62_domains.json)。
   探索領域の物理的正当化 = EFT 仮定 (脱結合 v5.3) であることの明示 (C2 の枠組)。
3. **Methods** — 3 独立実装 (DFS / 線形消去+グリッド / MITM) の一致、SHA-256 証明書、
   正準形 (gcd・共役・反転)。Lean 4 形式化: 毒値設計 (軌道外 +10⁶・燃料切れ +10⁹) に
   より「定理の成立が列挙の完全性を含意」する構成。信頼基盤の区別 (核 vs native_decide)。
4. **Results**
   - 最小性・一意性: {15: SM のみ}。スペクトル {15:1, 16:8, 17:1(弱三重項), 18:18, 22:2, 24:459}。
   - 頑健性の軸: 大表現 (6,6̄,8,三重項) / |Y|≤2,3 / ≤8 多重項 / ν_R (16 成分に SM+ν_R)。
   - **対照地図**: 最小性を担うのはカイラル性+全因子帯電、一意性を担うのは U(1)³。
     線形 3 アノマリーは領域内で冗長 (これ自体が発見)。
   - **U(1) の階段**: 2 本目は B−L が唯一 (rank-2 平面の Plücker 分類)、3 本目は
     存在しない (rank-3 全数 0 + 対照 355 で装置検証)。E6 の 27 のカイラル芯 = SM。
   - **例外群の完結** (v11.1): G₂/F₄/E₇/E₈ の代表表現は自己共役 (整数ウェイト系の
     機械検証; 対照 SU(3) 3 は複素) → カイラル芯は空。カイラル物質の入口は E₆ のみ。
   - GUT 降下の厳密整数検査 (SU(5)/SO(10)/Pati–Salam)。
5. **Machine verification** — Lean 定理 6 本の内容・実行時間 (991 s / 9734 s)・
   コスト予測法 (全再帰呼び出し数 × 0.6 µs)。
6. **Limitations** — 窓の有限性 (rank-3 の大電荷域、例外群、巨大表現)、16 成分隣人の
   排除が観測入力 (分数電荷ハドロン) に依ること、v4.3 域の Lean 未形式化。
7. **Reproducibility** — 外部依存なし Rust、固定シード、claims.yml 台帳、CI。

## 図表計画

- 図 1: 理論空間スペクトル (figures/v62_landscape.svg — 既存)
- 図 2: 対照地図 (figures/fig_controls_map.svg — v9.3 作成済み、数値は v62_atlas.json と照合)
- 図 3: U(1) の階段 (figures/fig_u1_staircase.svg — v9.3 作成済み、数値は v71/v82 JSON と照合)
- 表 1: 領域定義 (certificates から生成)
- 表 2: 3 実装の一致と SHA-256
- 表 3: Lean 定理と信頼基盤

## 主張と証拠の対応 (claims.yml)

| 主張 | id | 証拠 |
|---|---|---|
| 最小性・一意性 (v3.1 域) | QRN-GAUGE-003 | v31_gauge / v62_atlas / Anomaly.lean |
| 孤立スペクトル | QRN-GAUGE-006 | v43_landscape / v62_atlas |
| 大表現への頑健性 | QRN-GAUGE-007, -013 | v52_bigreps / v62_atlas / AnomalyBig.lean |
| 拡張軸・証明書・対照地図 | QRN-GAUGE-008, -010 | v62_atlas + certificates/ |
| 第 2 U(1) = B−L | QRN-GAUGE-009, -012 | v62_atlas / v71_twou1 |
| 第 3 U(1) なし | QRN-GAUGE-014 | v82_threeu1 |
| 機械検証 | QRN-GAUGE-011, -013 | proofs/*.lean |
