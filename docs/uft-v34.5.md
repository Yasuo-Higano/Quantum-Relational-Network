# QRN v34.5 — robust atlas: exact reader の同時信頼集合への持ち上げ

**Version**: v34.5 (2026-08-03)
**Sim**: `sim/src/bin/v345_robust_atlas.rs` → `results/v345_robust_atlas.txt` (13 検査 PASS)
**位置づけ**: PROMPT/15 §5 (v34.5)。v34.3 で確立した有限データ意味論
(区間・Bonferroni 同時割当・集合の全域で裁定一致のときだけ回答・跨ぎは Straddled)
を、第三十三期の各 exact reader に持ち上げる。**申告証明書 → 観測証明書**の器械。

---

## 1. 持ち上げの一覧

| reader | exact (v33) | robust (v34.5) |
|---|---|---|
| addressability rank/σ_min | 申告値 σ_min ≥ 0.5 | **同時下界** σ_lo = σ_min(M̂) − ‖ΔM‖_F ≥ 0.5 (Weyl) |
| cross-talk | 申告区間 ≤ 0.1 | **worst-case 上界** max+w ≤ 0.1 (平均は禁止 — 26) |
| glue orbit matching | overlap ≥ 0.9 (exact) | overlap **区間** vs 0.9 — 跨ぎは Straddled |
| charge witness → J | J² = −I ≤ 1e-9 (exact) | **spectral-gap 証明書**: 区間全域で σ_min(K) > 0 のときのみ J = K(−K²)^{−1/2} 構成 |
| resource cost | 申告 budget 半順序 | **interval cost** の 3 値採用 (確実採用/確実排除/跨ぎ) |
| profile | budget ↦ 読み | **set-valued profile** (跨ぎ点は class から除外して記録) |
| structured lane | dense と裁定一致 (exact) | **同じ区間意味論で一致** (graph 裁定・Straddled 込み) |

## 2. 主結果

- **[R0] 被覆の実例**: 合成 reader (2 統計量, Bonferroni α/2 ずつ) の全結果空間
  31² の厳密列挙 — 資格が偽の全境界セルで P(誤証明書) ≤ α (max 0.020)。
  v34.3 Robust Promotion の多統計量版の器械確認。
- **[R1] 下界の厳密性**: box 全 16 角で σ_min ≥ σ_lo (Weyl 束縛が有効)。tied
  (rank 1) は区間ごと確実拒否。平均 cross-talk 0.155 が通る box を worst-case
  上界 0.32 が拒否 — **禁止変換 26 の robust 版**。
- **[R2] glue**: site×CNOT 共役の overlap = **1/3 厳密** (v34.4 の値の再現) —
  区間裁定は matching (1.0) / 非同値 (1/3) / 境界 0.88±0.05 → Straddled。
- **[R3] J の構成条件**: 完全 charge は σ_min − ‖ΔK‖_F = 0.88 > 0 で
  J = K(−K²)^{−1/2} を構成 (J² = −I 厳密 0)。部分 charge (縮退 witness) は
  σ 区間が 0 を含み**構成拒否** — 「zero crossing がない場合のみ構成」の器械化。
- **[R4] interval cost**: grid 6 点 — 資源不足 → [2,2] (chain 2 = stable) →
  **跨ぎ点 Straddled** → 併合 [4] (chain 2 = stable)。**中点潰しの負制御**:
  点推定は跨ぎ点 b = 1.9 で確定読みを返す (禁止変換 18/22 の合流点)。
  set-valued profile の単調性を器械化 (関手性定理は成立まで語彙を凍結)。
- **[R5] dense/structured の同一区間意味論**: 統計量の lane 一致 (15 対)・
  **N = 8000 で site 族の graph 裁定確定 — N = 5000 では Bonferroni α/15 の
  CP 上限 1.3e-3 が τ = 1e-3 を跨ぎ確定不能** (分解能は観測量の関数, v343 [F2b]
  の実例)・不足ショット (N = 100) は両 lane とも Straddled が正答・対応原理
  (閉包/中心 = 2^{dim V}/2^{radical})。

## 3. 正直な境界

- **synthetic lane である**: 登録契約 = iid Bernoulli 統計量・決定的代表カウント。
  これを「実測ノイズ」とは呼ばない (PROMPT/15 §7 非交渉)。drift・相関・モデル
  不適合は禁止変換 24/25/29 が型で守り、real-data lane は v34.6。
- Bonferroni は保守的 — gauge orbit 上のタイトな同時領域 (GST 型) は未構成。
- set-valued profile は定義と単調性まで — 関手性・安定性定理は次期。

## 4. 次 (v34.6)

real-data lane の分離 (synthetic coverage lane / recorded experimental lane) と
clean-room freeze・D2-R 実募集の整備。
