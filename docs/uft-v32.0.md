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

## 4. v32.0-B (開封) — 本採点の確定表

SECRET を開示 (`HOLD7-f6618f3e15fc1566479e431de0bf6c59` — [H0] が sha256 = コミット
メント c50c4c30… と FROZEN 区間の逐語一致 [SHA e9d4079e…] を機械照合) し、holdout
17 セル (シード = SECRET 導出) を初生成・本採点した (調整なし):

| セル | 種 | 裁定 | 結果 |
|---|---|---|---|
| K1-0 twoholes | 応答→位相 | β=(1,2,0) open | ✓ |
| K1-1 genus2 | 応答→位相 | β=(1,4,1) closed | ✓ |
| K1-2 torus4 | 応答→位相 | β=(1,2,1) closed | ✓ |
| K2-0 **T³** | 応答→3D | **β=(1,3,3,1) closed** | ✓ (3D が holdout 抽選に出て生存) |
| K2-1 S³ | 応答→3D | β=(1,0,0,1) closed | ✓ |
| K3-ring21 | 計量 VR | death/perim = 0.3488 | ✓ (バー [0.28, 0.37]) |
| **K3-holes16:6** | 計量 VR | 比 5.00 vs バー 2.67±35% | **✗ 不成立** |
| K4-warm | oracle | err 2.3e-14 | ✓ |
| K4-cold (非識別) | oracle | EquivClass | ✓ 正棄却 |
| K5-projector (非識別) | oracle | EquivClass | ✓ 正棄却 |
| K6-oracle V=1.17 (非識別) | oracle | Abstain(NonGaussian) | ✓ 正棄却 |
| K6-resp V=1.17 | 応答 (t-V) | 支持一致 (厳密転移) | ✓ |
| K7-unknownfact (非識別) | oracle | Abstain(UnknownFactorization) | ✓ 正棄却 |
| K8-lownoise | 応答 | 支持一致 | ✓ |
| K8-highnoise (非識別) | 応答 | Abstain(InsufficientObservation) | ✓ 正棄却 |
| K9-metamorphic | 応答 変成対 | 次数列 + 重み多重集合 一致 | ✓ |
| K10-regulator | 到着時刻 | a 減衰比 6.57 (バー ≥ 2) | ✓ |

**確定: 16/17 生存 — coverage 1.00・selective risk = 1/12 (K3-holes のみ)・
impossibility recall 5/5 = 1.0・強制回答 0。**

### K3-holes 不成立の機構 (事後分析 — 調整・再採点なし)

凍結バーは「persistence 比 = 周長比」を仮定した。実際の VR H1 persistence は
**アフィン則 (L/3 − birth)**: holdout の (16, 6) ではアフィン予測 4.33 (周長比
2.67 の 1.63 倍 — バー外)、train の (14, 8) は 2.20 (1.26 倍 — 偶然バー内)。
**採点器のバーモデル誤りを holdout が捕捉した** — これが holdout の機能であり、
不成立は確定・バーは変更しない。教訓 = 「定量バーは観測量の応答則から導出する」
(v29.4a の教訓の再確認 — 次期の VR バーはアフィン則から導出する)。

なお **非識別 5 セルの正しい棄却 (impossibility recall 1.0) は初の holdout 通過** —
v31.0 で型として先凍結した「棄却を一級市民にする」採点原則が新鮮データで機能した。

---

## 5. 期末完全儀式 (v32.0-B)

共有部 (`lib.rs` の mod 宣言・`qrn_core.rs` の note・新規 `readout_contract.rs`) の
変更は v31.0 に集約したため、期末に一度だけ全数再計算を行った:

```text
make suite-full OUT=results/v320_full_suite.txt JOBS=12
→ 実行 187 本 / 総計 PASS 1314 / FAIL 0 (壁時計 ~11.3 h)
```

**末桁ドリフト検査 (台帳の期前後比較)**: 旧台帳 (e8f9daf, 177 本) と新台帳 (187 本)
を突き合わせ、**既存 177 本の PASS/FAIL は完全一致 (ドリフト 0 件)**。差分は新規
10 本 (v310–v317, v320a/b) の +87 PASS のみ — **第三十一期の型改訂・器械最終化は
既存物理に無波及**。

## 6. 第三十一期の統合 — 確定残高

**期テーゼ: 「読めることと、読めないと言えることは、同じ資格の両面である」**

| 版 | 成果 (確定) |
|---|---|
| v31.0 | 意味論・型・protocol の凍結 — E0–E4 の型分離・`ReadoutCertificate<能力,状態領域,観測契約,因子分解>`・裁定 5 値 + 棄却理由 8 種・**禁止変換 8–10** (親→物理生成子 / 正則化→exact / oracle→operational patch)・schema の RelationalDecomposition 意味論差是正・replications.yml の capability scoped 化 |
| v31.1 | **GaussianGibbsInverseOracle** (E1 の天井, oracle ceiling): n≤7 全数 992 グラフ × 4β 資格満票・**条件数定理 ‖ΔK‖≤‖ΔC‖/(δ(1−δ))** (rank-1 整列で飽和 0.99996)・**発見: P6/693 の projector 衝突は Z2 ゲージ×置換の厳密同値** (full C でも識別不能)・低温飽和則 (深部は sign のみ) |
| v31.2 | **LocalBiasCommutatorLaw** (不変ノルム核): (n̈⁺−n̈⁻)/(4ε) = ‖P_jhP_i‖_F² の恒等式・密度時系列のみで rel 4e-9・gauge 不変・**P6/693 を応答が分離 (静的不可・応答可の恒等式化)**・能力階層の実証 (ArrivalTime = 重みの下流) |
| v31.3 | **観測予算 hierarchy** (7 lane): oracle → coherent → 密度応答 → **patch 1.7%** → pair-B2 20–30% → B3 単調 → 到着時刻。**「encoded but not operationally readable」の実例は pair 準位**・応答 lane の代価 = ノイズ 5903× 増幅 |
| v31.4 | **不変 atlas と E3 no-go**: whitened 不変量は演算子基底の任意可逆再結合に不変・**B3/B4 = atlas の成分/圧縮 (統一)**・**同一状態が 3 因子分解で 3 幾何**・疎性基準は幾何を自明化 (負制御)・**因子分解は OperationalAlgebra が運ぶ** |
| v31.5 | **非 Gaussian transfer**: **曲率則は密度対角相互作用に厳密転移** (V ∈ {0,2,4} で 1.000000000000)・oracle は Gaussian-only 確定 (witness 超え → 棄却)・強結合で静的 B3 が破れ応答は無傷・Z2 拘束系は位相のみ転移 |
| v31.6 | **3D + 計量**: β₃ に im ∂₄ 必須 (K5 anchor)・**from-state で T³/S³/3-ball を link 込み同定 (初の 3 次元)**・VR persistence の円環 1/3 法則・**2 穴の寿命分離 (v29.6 限界の解消)**・**genus-2 β₁=4**・安定性定理の bottleneck 照合・gap 抽出器の最終凍結 |
| v31.7 | **開放境界機構**: hold-12 の不一致 = **振幅消失型 reflection artifact** (壁位置差 Δδ=0.60 格子・振幅 ~O(a) 消失) — 境界普遍類の差ではない。**バーと不成立記録は不変** |
| Track X | **外部再現 Unit D** (D1 数値再現 / **D2 = geometry blocker 解除の唯一路** / D3 負定理 / D4 応答法則) を公開・`v310 [R8]` が常設監査 |
| v32.0 | **HOLD-7**: 凍結 (train 満票) → 開封 **16/17 生存**・**非識別 5/5 正棄却 (impossibility recall 1.0 の初通過)**・K3-holes 1 セル不成立確定 (バーモデル誤りを holdout が捕捉)・儀式 187 本 PASS 1314/0 (ドリフト 0) |

**正直な残高 (変わらないもの)**:
- **bridge law 登録簿は全能力で空のまま** — HOLD-7 の 16/17 生存をもってしても登録
  しない。blocker は独立外部再現 0 (R2/R3) — **Unit D2 の公募が最優先**。
- PRED-019 未登録・自然の観測量の的中 0・`external_replications = 0`。
- scope は「**与えられたノード因子分解の下で**」— ただし v31.4 で E3 の no-go が
  定理化され、FactorizationBridge には OperationalAlgebra が必要と確定した。
- **不成立の記録 2 件**: v29.4b hold-12 (regulator 1 対) と v32.0-B K3-holes
  (VR バーモデル) — どちらも調整せず確定し、機構を解剖して次期の設計入力にした。

**未解決 (第三十二期への課題)**:
1. **外部独立再現 (最優先・据え置き)** — Unit D2 の実施者募集。
2. VR バーのアフィン則導出 (K3-holes の教訓) と H2 persistence。
3. 相関ホッピング型相互作用での曲率則の破れ検証 ([H_int, n_j] ≠ 0 の系)。
4. 非自明 3 多様体 (レンズ空間等)・大型 3D 系・4-simplex 超の次元。
5. OperationalAlgebra の構成的定義 (E3 の正の側 — no-go の裏返し)。

## 7. 開発記録 (第三十一期)

- **器械の誤りを holdout が 2 度捕捉した**: v31.6 の gap 抽出器 (窓跨ぎ段差) は
  設計走行で、v32.0-B の VR バー (persistence の比例仮定) は holdout 本番で。
  前者は器械訂正として全バイナリに統一適用 (v31.5 [N5] の裁定が更新された —
  監査注記で追跡可能)、後者は不成立として確定 (バーは変更しない)。
- **統計量の参照系を誤る故障が 3 回** (v31.7): 局所 v̂ の離散化ノイズ・bulk 差の
  未除去・LDOS の bulk DOS 差。教訓 =「regulator 比較は自分の bulk を参照してから
  境界を見る」。
- **1D 調整値の暗黙拡張**を 2 度自ら踏んだ: 固定 dt (次数 14 の T³ で O(1) 誤差 →
  スペクトル半径スケールに昇格)・固定 gap 則。本期の禁止事項がそのまま自分に
  返ってきた形で、器械契約への昇格で根治した。
- **flag 性の破れ** (T³ L=3 の巻き付き 3-クリーク) は、幾何が小さすぎると
  clique complex が三角化と別物になるという系サイズ下限の発見 — v29.5 Petersen
  縮退と同族の教訓。
