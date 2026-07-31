# QRN v30.0 — 第三十期 統合: 能力別 certificate・HOLD-6・期の残高

**Version**: v30.0-A (2026-07-30, コミット 8024f17 — 凍結半) / v30.0-B (2026-07-31 —
儀式 + holdout 開封 + 期統合)
**Sim**: `sim/src/qrn_core.rs` (能力別 certificate) → `results/v272_core_contract.txt`
(7 検査 PASS) / `sim/src/bin/v300a_hold6_freeze.rs` → `results/v300a_hold6_freeze.txt`
(4 検査 PASS) / `sim/src/bin/v300b_hold6_open.rs` → `results/v300b_hold6_open.txt`
(器械 2 検査 PASS + [採点])
**儀式**: qrn_core 共有部変更に伴う完全再走 — **全 176 本 実行, 総計 PASS 1224 /
FAIL 0** (`results/v300_full_suite.txt`, 壁時計 ~15.3 h, JOBS=12)。**末桁ドリフト:
既存 173 本の PASS/FAIL は儀式前後で台帳比較により完全一致** (差分は新規 3 本
[v295 11 + v296 6 + v300a 5 = +22] のみ) — 能力別改訂は既存物理に無波及。
**位置づけ**: PROMPT/11 ロードマップの最終版 = 第三十期の期統合。

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

## 3. v30.0-B — HOLD-6 の開封と本採点 (確定表)

SECRET を開示 (`HOLD6-b4b84a54679589ad2773200a3251c43d` — [H0] が sha256 =
コミットメント fe3c9cbd… と train seed の一致を機械照合) し、hold-0..7 を初生成・
本採点した (調整なし):

| instance | class (seed 決定) | n | β 実測 | 曲面裁定 | 裁定 |
|---|---|---|---|---|---|
| hold-0 | cylinder | 35 | (1,1,0) | boundary | バー内 |
| hold-1 | two-holes | 149 | (1,2,0) | boundary | バー内 |
| hold-2 | cylinder | 45 | (1,1,0) | boundary | バー内 |
| hold-3 | cylinder | 48 | (1,1,0) | boundary | バー内 |
| hold-4 | torus | 64 | **(1,2,1)** | closed | バー内 |
| hold-5 | two-holes | 116 | (1,2,0) | boundary | バー内 |
| hold-6 | sphere | 42 | **(1,0,1)** | closed | バー内 |
| hold-7 | sphere | 42 | (1,0,1) | closed | バー内 |

**全 8 系バー内 (満票)** — 全会一致・採用窓 4/4。dimension-agnostic pipeline が
**新鮮な重み場つき 2D holdout で生存**した。注記: disk クラスは holdout 抽選に
出なかった (seed の選択 — train でのみ検証済み。恣意ではないことは seed 系列の
第三者検証で確認可能)。

## 4. 第三十期の統合 — 確定残高

**期テーゼ: 「監査が裁定を正し、凍結が誠実を守り、読み出しが次元を跨いだ」**

| 版 | 成果 (確定) |
|---|---|
| v29.1 | Integrity Erratum — linfit 切片バグの再現・訂正。**B2 復活・B1 のみ棄却**。LinearFit 型 + 変成テスト常設 |
| v29.2 | 意味論再基礎化 — B5/B6 定義凍結・S×C 採点原則・二層分離・HOLD-5 コミットメント |
| v29.3 | S×C 合成 train 288 セル満票 (Wilson 実対称化・橋なし優先 pipeline) |
| v29.4a | val 1 回使用 — 資格健全 + 定量バー機械導出 + **SECRET 開示** |
| v29.4b | **HOLD-5 本採点**: 資格 576 満票・**v̂ 真値照合 24/24 (円環 Δ∞ ≤ 7%)・τ 予言 24/24・regulator 間 23/24** (1 対不成立を確定 — 開放鎖境界系統) |
| v29.5 | collision atlas — **静的核衝突対 (P6, 693) 発見** (静的核は sign(A) まで — KIN が 0.5 で分離 = 「静的単独不可・応答併用可」の最小実例)・Petersen 誤認の機構同定・factorization 選定不能の機械記録 |
| v29.6 | **dimension-agnostic pipeline** — torus (1,2,1)/cylinder/disk/sphere (1,0,1) を状態から end-to-end 同定。臨界 GS の境界増強 2D 版を発見 |
| v30.0 | 能力別 certificate (ProperTime への門は関数ごと不在)・**HOLD-6 満票**・儀式 PASS 1224/0 |

**正直な残高 (変わらないもの)**:
- **bridge law 登録簿は全能力で空のまま** — HOLD-5/6 の生存をもってしても登録
  しない。blocker は独立外部再現 0 (R2/R3 の機械化)。登録に将来最も近い能力は
  SpatialTopologyGivenFactorization (HOLD-6 満票 + v29.5 の同値類証明書が裏付け)。
- PRED-019 未登録 (QRN 固有の数値の解析的導出なし)・自然の観測量の的中 0。
- scope はすべて「**与えられたノード因子分解の下で**」— FactorizationBridge は
  未定義のまま (v29.5 [C5] がその空隙自体を機械記録)。

**未解決 (第三十一期への課題)**: (i) 不変ノルム核 (PROMPT/11 第二課題の残り —
B3-COV と B4 を同一応答 atlas の断面として統一する superoperator ノルム)。
(ii) 臨界状態に頑健な 2D 核 (境界増強・Friedel 共鳴を除去する読み出し)。
(iii) 計量 Vietoris–Rips persistence (穴スケールの寿命分離)・高 genus・3D。
(iv) regulator 間不一致 1 対 (開放鎖境界) の機構の定量化。
(v) **外部独立再現の公募 (最優先)** — reproducer/ は 3 単位を公開済み、
replications.yml は空のまま。

## 5. 開発記録 (v30.0)

- ピンチ点の発見が [G1b] を生んだ: holdout の「期待プロファイル」も検査対象で
  ある — 採点器だけ検査して生成器を信じるのは非対称だった。
- 能力別 certificate の「構成不能」は 3 段構え: (i) 能力トレイト未実装 (型レベル)
  (ii) sealed で外部実装封鎖 (iii) 登録簿空 (値レベル) — Lean の昇格不能定理
  (QrnPromotion.lean) と合わせて 4 重。
- 儀式は 4 回目 (v28.0, v29.1, 本版 ×1 + v25 期) — 二相実行 (監査層後段) 以降、
  偽 FAIL ゼロが続いている。台帳比較 (manifest tsv の git 前後 diff) による
  ドリフト検査を儀式の標準手順に昇格。
