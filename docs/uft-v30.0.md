# QRN v30.0 — 第三十期 統合: 能力別 certificate・HOLD-6・期の残高

**Version**: v30.0-A (2026-07-30 — 凍結半) / v30.0-B (holdout 開封 + 期統合 — 未実行)
**Sim**: `sim/src/qrn_core.rs` (能力別 certificate) → `results/v272_core_contract.txt`
(7 検査 PASS) / `sim/src/bin/v300a_hold6_freeze.rs` → `results/v300a_hold6_freeze.txt`
(4 検査 PASS)
**位置づけ**: PROMPT/11 ロードマップの最終版。qrn_core の共有部変更を伴うため、
本コミット後に**完全再走の儀式** (全 176 本) を実施し、その後 v30.0-B が HOLD-6 を
開封して期を閉じる。

---

## 1. 能力別 BridgeLawCertificate (v30.0-A — PROMPT/11 の型改訂)

旧 `BridgeLawCertificate` は単一の claim id しか持たず、**空間位相の証拠と
proper time の証拠を型で区別できなかった**。v30.0 で能力別に分割した:

```rust
pub trait BridgeCapability: sealed { const NAME; const REGISTERED: &[&str]; }
pub struct BridgeLawCertificate<C: BridgeCapability> { claim_id, PhantomData<C> }
```

- **8 能力タグ** (空 enum): FactorizationGivenObservables /
  SpatialTopologyGivenFactorization / SpatialMetricUpToGlobalScale /
  CausalOrderGivenExternalClock / ConformalLorentzianStructure / VolumeMeasure /
  ClockCalibration / FullLorentzianMetric。
- 現行証拠で将来到達しうる上限 = SpatialTopology / SpatialMetric / CausalOrder
  の 3 能力 (PROMPT/11)。**登録簿は全能力で空** — v29.4b の生存 (真値照合 24/24
  等) をもってしても登録しない (独立外部再現 0 のため — R2/R3 の機械化と整合)。
- **ClockCalibration / FullLorentzianMetric は BridgeCapability 未実装** —
  `BridgeLawCertificate<ClockCalibration>` はトレイト境界で**型レベル構成不能**。
  旧 `promote_evolution_to_proper_time` (発展パラメータ → 固有時間の門) は
  **関数ごと削除** — 禁止変換 2 は「門となる関数が書けない」ことで強制される。
- sealed パターンで外部からの能力追加も封鎖。検査: `qrn_core_self_test`
  (全能力の施錠) + `v272_core_contract` [T3]/[T4] (source 検査 + 6 能力 × 5 id の
  発行拒否)。schema に概念 CapabilityIndexedCertificate を登録。

## 2. HOLD-6 の凍結 (v30.0-A) — 2D topology pipeline の新鮮 holdout

v29.6 の資格は設計区画で検証力がない。HOLD-5 と同じ開封順序を 2D で敷いた:

```text
コミットメント公表 + 生成器/採点器/バー凍結 (v30.0-A = 本コミット)
  → SECRET 開示 + hold-0..7 初開封・本採点 (v30.0-B, 調整なし)
sha256(SECRET) = fe3c9cbd0c2d733f852422734ca4f212cd01e488fac763142ceef07d598d6b62
```

- **生成器** (FROZEN TOPO v32, SHA-256 = e4a1472f… — v30.0-B が照合): クラス ∈
  {torus, cylinder, disk, two-holes, sphere} (holdout は seed が決定)・サイズ乱択・
  **滑らかな速度場 v = 1 + Σ A_k exp(−r²/w²) の重みつきボンド** (v ∈ [0.65, 1.45]
  棄却保証 — スケールガード W_FRAC = 0.15 に対し最悪 nn 核比 0.20 の余裕。周期軸は
  min-image で場も周期化。球は角距離ガウス束)・ノード置換。状態 = 熱的 Gaussian
  (βt = 1, v29.6 の資格状態族)。採点 pipeline は v296 の逐語コピー。
- **バー (離散・事前登録)**: (Betti, 曲面性) がクラス期待と厳密一致 + 窓全会一致 +
  採用窓 ≥ 2。
- **[G1] train 5 系 (クラス別 1 系) 全て資格**: torus (1,2,1)/cylinder (1,1,0)/
  disk (1,0,0)/two-holes (1,2,0)/sphere (1,0,1) — 重み場・置換つきで満票。
- **[G1b] 生成器健全性 (新設ゲート)**: train 全系の**真の複体**が多様体かつ期待
  Betti。設計走行の発見: 貪欲成長の穴クラスタは凹形になり**真の複体にピンチ点
  (非多様体頂点) を作る** — pipeline の「not-surface」裁定が正しく、期待側が
  誤りだった。穴を凸ブロックに変更して根治 (凸集合の除去は各境界頂点の失う近傍が
  連続 → link は弧のまま)。「生成器の期待プロファイル自体を検査する」ゲートとして
  常設。

## 3. v30.0-B (本コミット時点で未実行)

儀式 (全 176 本再走 — qrn_core 変更の検証) の完走後に: SECRET 開示 → hold-0..7
初開封・本採点 (調整なし・[採点] 行) → 期統合 (第三十期の確定残高・教訓・
第三十一期への課題) を本文書に追記して期を閉じる。

## 4. 開発記録 (v30.0-A)

- ピンチ点の発見が [G1b] を生んだ: holdout の「期待プロファイル」も検査対象で
  ある — 採点器だけ検査して生成器を信じるのは非対称だった。
- 能力別 certificate の「構成不能」は 3 段構え: (i) 能力トレイト未実装 (型レベル)
  (ii) sealed で外部実装封鎖 (iii) 登録簿空 (値レベル) — Lean の昇格不能定理
  (QrnPromotion.lean) と合わせて 4 重。
