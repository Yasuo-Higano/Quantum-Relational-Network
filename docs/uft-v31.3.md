# QRN v31.3 — 観測予算 hierarchy: 7 lane の識別可能性相図

**Version**: v31.3 (2026-07-31)
**Sim**: `sim/src/bin/v313_budget_hierarchy.rs` → `results/v313_budget_hierarchy.txt`
(10 検査 PASS)
**位置づけ**: PROMPT/12 第四版 — 「global に符号化されている (E0/E1)」と「その観測
契約で読める (E2)」の分離を、同一 hidden generator への 7 lane 一斉適用で機械化。

---

## 1. 設定

隠れ族 3 系 (重み場つき ring12 / 開鎖 chain12 / 単環 u693 — 決定的解析形) に対し、
観測契約の異なる 7 lane を同時適用。真値は採点でのみ参照 (readout には不可視)。

## 2. 確定した hierarchy (β = 1, 重み復元誤差)

| lane | 観測契約 | 重み誤差 | 位相 (支持) | 何が読めるか |
|---|---|---|---|---|
| 1 oracle | GlobalOneBodyCorrelation | 1.2e-14 | 欠0余0 | h 厳密 (E1 の天井) |
| 6 coherent | CoherentLocalResponse | 1.1e-7 | — | h block (Z2 ゲージ木伝播つき符号込み) |
| 5 密度応答 | LocalBiasDensityResponse | 3.7e-8 (重み²) | — | \|h_ij\|² (gauge 不変・状態非依存) |
| 2 patch | OperationalPatch (B3 支持の半径 2) | **1.7%** | — | 環境 renormalization が patch で減衰 |
| 3 pair-B2 | PairReducedStates | **20–30%** (系統) | 欠0余0 | 支持は正・重みは renormalized |
| 4 B3 | StaticLocalObservables | 単調 proxy (Spearman 0.973) | 欠0余0 | 重み順位のみ |
| 7 到着時刻 | ArrivalTimeResponse | 距離 proxy (Spearman 0.991) | — | 重みつき最短路の単調像 |

- **patch inversion の中間性が確定**: pair (2 サイト) → patch (半径 2) で系統偏差
  20.4% → 1.7%。観測を広げるほど環境 renormalization が減る (patchwise inversion
  仮説の最初の定量点)。OperationalPatch は B3 観測だけから構築し、OraclePatch
  (真の半径 — 診断専用) と型で分離 — 本走行では B3 支持が正しいため両者は一致
  (一致すること自体を検査)。
- pair-B2 の 30% 系統偏差は**失敗ではなく契約準位の性質** (v31.1 [T6] の
  f(PCP) ≠ Pf(C)P の定量化)。

## 3. 相図スライス (状態領域 × lane, ring12)

| 状態領域 | global exact | global estimate | pair-B2 | B3 | 応答 lane 5/6 |
|---|---|---|---|---|---|
| β=1 熱的 | Exact | — | 支持正・重み 20% | 支持正・単調 | Exact |
| β=25 深部 | **正しく棄却** (margin < 床) | gap 支持 欠0余0 | **gap 支持 破れ (余22)** — ranking は生存 | gap 支持 欠0余0 | **不変** (3.7e-8) |
| projector GS | RankDeficient 棄却 (sign 類 — v31.1 [T9]) | 同値類のみ | — | — | **不変** |

**発見 ([H8]): 「encoded but not operationally readable」の実例は pair 準位に出る** —
β=25 で pair rdm は縮約混合により full-rank のまま資格を通る (棄却が起きない) のに、
スケールガード付き gap 抽出器では支持が読めない (余22)。一方 ranking 準位
(precision@12) では全 lane 1.0 — **情報は符号化されているのに、その観測契約の
operational 抽出器が読めない**。失敗を消さず機械記録した。

## 4. ノイズの代価 ([H9])

σ = 1e-4 のノイズで、静的 lane (K 誤差 1.7e-3) に対し応答 lane は **5903×** の誤差
増幅 (9.8) — 二階差分 1/dt² の代価 (予言オーダと整合)。**応答 lane は状態非依存性を
ノイズ増幅で買っている** — trade-off の機械記録。HOLD-7 のノイズ軸の設計入力。

## 5. 正直な残高

- 隠れ族は 1D 系 + u693 の 3 系 — 2D・高次元の hierarchy は v31.6 以降。
- gap 抽出器はスケールガード付き最大対数ギャップ則 (v29.6 と同思想) — 抽出器を
  変えれば「readable」の境界は動く (ranking との差が示すとおり)。抽出器の凍結は
  HOLD-7 で行う。
- ノイズは観測量への iid Gaussian — 相関ノイズ・デコヒーレンスは未走査。
- 開発記録: 初版の gap 則 (ガードなし) は f64 ノイズ床の尾部 (対数比が発散) を
  カットし全系で偽の余剰辺を出した — v29.6 のスケールガードの必然性を独立に
  再発見した形。β=25 の pair 全対棄却という初期予想も誤り (縮約混合で pair rdm は
  常に full-rank) — 「棄却が起きないのに読めない」がこのセルの正しい記述。
