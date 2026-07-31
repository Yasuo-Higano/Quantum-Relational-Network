# QRN v32.0 — HOLD-7: identifiability 相図の holdout と第三十一期統合

**Version**: v32.0-A (2026-07-31 — 凍結半) / v32.0-B (開封 + 期統合)
**Sim**: `sim/src/bin/v320a_hold7_freeze.rs` → `results/v320a_hold7_freeze.txt`
(3 検査 PASS + train 17 セル満票)
**位置づけ**: PROMPT/12 の最終版 — 第三十一期の全器械 (v31.0–v31.7) を凍結し、
identifiability 相図**全体**を新鮮 holdout で検証する。

---

## 1. HOLD-7 が HOLD-5/6 と違う点

**非識別セルの正しい棄却を採点する** (v31.0 で型として先凍結した採点原則の執行):

- selective risk (回答セルの誤り) = 0 を要求
- coverage (回答可能セルの回答率) ≥ 0.9
- **impossibility recall (非識別セルの正しい棄却率) = 1.0 — 強制回答は FAIL**
- 変成対 (置換ゲージ) の読み出し一致・真値の readout 非流入 (採点器のみ真値参照)

## 2. 開封順序 (HOLD-5/6 と同一)

```text
v32.0-A (本コミット) = 生成器・採点器・バー・観測契約・棄却規則の凍結
                      + SECRET コミットメント公表 + train 採点 (可視シード 32001)
  → v32.0-B = SECRET 開示・holdout 初生成・本採点 (調整なし) + 期末完全儀式
sha256(SECRET) = c50c4c30993b1bb7113734609e87f116e8f280a30315571261f92377bd3ec9ea
```

凍結カーネル (FROZEN-HOLD7 区間, sha256 = e9d4079e…) は v31 期の器械の逐語コピー:
曲率則測定 lane (スペクトル半径スケール dt)・最終 gap 則 (窓内有意段差 ≥ ln 3)・
Z2 homology (∂₄ 込み)・link 分類・VR H1・GaussianGibbsInverseOracle・Wick witness・
t-V ED・凍結棄却規則 4 本 (因子分解なし → 棄却 / witness 超え → 棄却 / RankDeficient →
同値類 / ノイズ誤差見積り超え → 棄却)。

## 3. セルクラス (10 種 17 セル — 相図の全軸)

| クラス | 状態 × 観測 × 軸 | 要求 |
|---|---|---|
| K1 ×3 | 熱的 × 応答 × topo2d (torus/genus2/two-holes, 重み場+置換) | (β, 曲面種) 回答 |
| K2 ×2 | 熱的 × 応答 × topo3d (T³ L=4/16-cell/ball) | (β₀..β₃, link 種) 回答 |
| K3 ×2 | 熱的 × 応答 × 計量 (ring 1/3 法則・2 穴寿命比) | 定量回答 |
| K4 ×2 | Gibbs β=1 / **β=30** × oracle | 回答 / **棄却・同値類** |
| K5 ×1 | **projector (P6)** × oracle | **同値類** |
| K6 ×2 | **t-V 相互作用** × oracle / 応答 | **棄却** / 回答 (厳密転移) |
| K7 ×1 | **因子分解不明** (稠密回転基底) | **棄却** |
| K8 ×2 | ノイズ σ=1e-9 / **σ=1e-3** × 応答 | 回答 / **棄却** |
| K9 ×1 | 変成対 (置換) | 読み出し一致 |
| K10 ×1 | regulator (R-A×R-C 開放鎖) | 端バイアス a 減衰 ≥ 2× |

train (シード 32001, 可視) は **17 セル満票**: 回答 12/12 (coverage 1.00)・誤り 0・
非識別 5/5 正棄却 (recall 1.00)・強制回答 0。内訳: torus5 (1,2,1)・twoholes (1,2,0)・
genus2 (1,4,1)・S³ ×2 (1,0,0,1)・ring24 death/perim 0.3428・holes 比 2.00 (真 1.75)・
oracle 2.0e-14・cold/projector → EquivClass・t-V oracle → NonGaussian 棄却・
unknownfact/highnoise → 棄却・変成一致・regulator 減衰 6.57。

## 4. v32.0-B (開封) — 本採点の結果

(v32.0-B で追記)
