# QRN v34.0 — HOLD-9: 出自 × 文脈整合 × 資源依存局所性の holdout と第三十三期統合

**Version**: v34.0-A (2026-08-02 — 凍結半) / v34.0-B (開封 + 期統合)
**Sim**: `sim/src/bin/v340a_hold9_freeze.rs` → `results/v340a_hold9_freeze.txt`
(3 検査 PASS + train 20 セル満票) / `sim/src/bin/v340b_hold9_open.rs` →
`results/v340b_hold9_open.txt` (開封 3 検査 PASS + holdout 本採点)
**位置づけ**: PROMPT/14 の最終版 — 第三十三期の全器械 (v33.1 修復入口 / v33.2
laboratory interface / v33.3 resource profile / v33.4 atlas glue / v33.5 graded
recovery / v33.6 structured backend) を凍結し、**操作の出自 × 文脈整合 ×
資源依存局所性**の識別可能性境界を新鮮 holdout で採点する。

---

## 1. HOLD-9 が採点するもの (HOLD-8 の後継)

従来の 3 計量:

```text
selective risk        = 0   (回答セルの誤り 0 — 強制回答も誤りに数える)
impossibility recall  = 1   (非識別セルの正しい棄却/同値類/構成時拒否)
answerable recall     = 1   (回答可能セルは全て答える)
```

に、第三十三期の 5 計量 (PROMPT/14) を加える:

```text
origin_certificate_coverage          = 1.0  (admit された全操作に出自証明書)
context_witness_coverage             = 1.0  (全 Exact 読みが証人ゲート通過)
raw_operation_promotions             = 0    (証明書なし/流用の昇格 0)
scope_violations                     = 0    (scope 外の裁定強行 0)
transient_factorization_promotions   = 0    (単点読みの局所性昇格 0)
```

## 2. 開封順序 (HOLD-5..8 と同一)

```text
v34.0-A (本コミット) = 生成器・採点器・バー・lib pin の凍結
                      + SECRET コミットメント公表 + train 採点 (可視シード 34001)
  → v34.0-B = SECRET 開示・holdout 初生成・本採点 (調整なし) + 期末完全儀式
sha256(SECRET) = ef6a8cd97b5d7693f4a4ffdb11ccdf42dbf5a971b0c79ba5f0f72f7b15739fcd
```

凍結カーネル (FROZEN-HOLD9 区間 — v340b と逐語一致を [H0] が照合) は第三十三期の
器械を**呼ぶ**駆動部であり、器械本体は lib 6 モジュール。**lib pin**: 凍結時の
sha256-16 (operational_net `7898d244…` / laboratory_interface `a11f188f…` /
resource_profile `f2fe4b96…` / contextual_factorization `e540eea6…` /
graded_recovery `c1ebb02c…` / structured_backend `f57d9bdf…`) を [A0]/[H0] が
照合する — 凍結から開封の間に器械が変わらないことの機械保証 (HOLD-8 の逐語コピー
方式の lib 版)。凍結バー: τ = 1e-3・σ_bar = 0.5・xtalk_bar = 0.1・合成/CAR/J バー
1e-9・orbit 非同値 < 0.9・transient 規則 chain ≥ 2 — 全て第三十三期の各版で凍結済み
の値の再固定。

## 3. セル 5 群 20 セル (隠しパラメータ = 変成 U×置換・対の抽選・γ・コスト・O(6) 回転・汚染振幅・qubit 数)

| 群 | セル | 要求 |
|---|---|---|
| 入力完全性 | IC1 raw + 流用証明書 | **CertificateTargetMismatch 拒否** (raw_operation_promotions = 0) |
| | IC2 net-owned gens | 復元行列 = net primitive **byte 同一** + Exact [2,2,2] |
| | IC3 role-mixed | **RoleMixedRecovery 構成時拒否** |
| | IC4 GKLS → 可換子 lane | **NonDerivation 拒否** (Leibniz 破れ) |
| accessibility | AC1 独立 knobs (変成) | Exact [2,2,2] (出自証明書つき) |
| | AC2 tied (隠し対) | 分解は rank 拒否 + 正直 net は **Abstain** |
| | AC3 同一閉包・別 interface | 両方 Exact・**orbit matching 不在** (非同値) |
| | AC4 合成証明書の流用 | **CertificateTargetMismatch 拒否** |
| context・resource | CR1 budget profile | 資源不足 → [2,2,2] → [2,4] (隠しコスト) |
| | CR2 一意 glue | GluedExact [2,2,2] (frame 回転不変) |
| | CR3 複数 glue | **EquivalenceClassOnly{2}** (tie-break なし) |
| | CR4 cross-talk 跨ぎ | **Straddled 棄却** (強制判定なし) |
| graded | GR1 odd CAR only | **MajoranaFrameOnly** (O(2N) orbit で止まる) |
| | GR2 charge witness | complex modes (**Σâ†â = Q ≤ 1e-9**) |
| | GR3 汚染 witness | **WitnessNotLinearOnFrame 棄却** |
| | GR4 ordinary odd | **構成時拒否** (JW 誤読の遮断) |
| scale・変成 | SC1 dense = Pauli | 隠しセル + 隠し置換で**裁定完全一致** |
| | SC2 大型 structured | Pauli Exact (隠し 32–60 qubit・行列なし)・**dense は ScopeExceeded** |
| | SC3 大域共役の共変 | dims 不変・orbit は W-共役で厳密対応 |
| | SC4 transient 非昇格 | stable = [2,2,2] のみ・単点 [2,4]/[8] は transient (**昇格 0**) |

**train (シード 34001, 可視) は 20 セル満票**: selective risk 0.000・impossibility
recall 1.000 (11/11)・answerable recall 1.000 (9/9)・強制回答 0・出自被覆 31/31・
証人被覆 7/7・昇格/違反 0。設計走行で 12 の独立シードでも満票 (生成器の
シード頑健性 — 隠しパラメータの縮退による器械故障がないことの設計時確認)。
*(訂正注記 [v34.0-B]: 本文書の A 版は区画を 12/8 と記していたが誤記 — 凍結
カーネルの answer_cell 宣言から機械的に定まる区画は 非識別 11 / 回答 9。機械出力は
区画数を印字せず、採点・裁定は無関係。)*

**注**: 外部再現 (Track X, D2-R) の成否は HOLD-9 のセルに**含まれない** — 内部
SECRET から作ったセルは外部独立性ではない (PROMPT/14 非交渉)。

## 4. v34.0-B (開封) — 本採点の確定表

SECRET を開示 (`HOLD9-44eed4c2f12b4019ed013ec44607b6b5` — [H0] が sha256 =
コミットメント ef6a8cd9… と FROZEN-HOLD9 区間の逐語一致 [kernel sha 317ccd3a…]・
**lib pin 6 モジュール不変**を機械照合) し、holdout 20 セル (シード
17251756188301620883 = SECRET 導出) を初生成・本採点した (調整なし):

| セル | 裁定 | 結果 |
|---|---|---|
| IC1 raw + 流用証明書 | CertificateTargetMismatch 拒否 | ✓ 正棄却 |
| IC2 net-owned gens | 復元行列 = primitive byte 同一 + Exact [2,2,2] | ✓ |
| IC3 role-mixed | RoleMixedRecovery 構成時拒否 | ✓ 正棄却 |
| IC4 GKLS → 可換子 lane | NonDerivation 拒否 | ✓ 正棄却 |
| AC1 独立 knobs (変成) | Exact [2,2,2] (出自つき) | ✓ |
| AC2 tied (隠し対) | 分解 rank 拒否 + Abstain(Insufficient) | ✓ 正棄却 |
| AC3 同一閉包・別 interface | 両 Exact・orbit matching 不在 (0.5625) | ✓ 正棄却 |
| AC4 合成証明書の流用 | CertificateTargetMismatch 拒否 | ✓ 正棄却 |
| CR1 budget profile | 資源不足 → [2,2,2] → [2,4] | ✓ |
| CR2 一意 glue | GluedExact [2,2,2] (frame 回転不変) | ✓ |
| CR3 複数 glue | EquivalenceClassOnly{2} | ✓ 正棄却 |
| CR4 cross-talk 跨ぎ | Straddled (強制判定なし) | ✓ 正棄却 |
| GR1 odd CAR only | MajoranaFrameOnly | ✓ 正棄却 |
| GR2 charge witness | 3 モード回復 (Σâ†â = Q) | ✓ |
| GR3 汚染 witness | WitnessNotLinearOnFrame | ✓ 正棄却 |
| GR4 ordinary odd | 構成時拒否 | ✓ 正棄却 |
| SC1 dense = Pauli | 隠しセル + 置換で裁定一致 | ✓ |
| SC2 大型 structured | Pauli Exact (行列なし)・dense は ScopeExceeded | ✓ |
| SC3 大域共役の共変 | dims 不変・orbit 厳密対応 | ✓ |
| SC4 transient 非昇格 | stable = [2,2,2] のみ | ✓ |

**確定: 20/20 満票 — selective risk 0.000・impossibility recall 1.000 (11/11)・
answerable recall 1.000 (9/9)・強制回答 0・origin_certificate_coverage 1.000
(31/31)・context_witness_coverage 1.000 (7/7)・raw_operation_promotions 0・
scope_violations 0・transient_factorization_promotions 0。**

**器械訂正 1 件 (開封走行で発見・凍結区間に触れない)**: v340b [H2] (開封手続きの
自己整合検査) の区画定数の誤記 8/12 → 9/11 — 凍結カーネルの answer_cell 宣言から
機械的に定まる値で、初版の [H2] が**正しく FAIL して検出**した (採点器・裁定・
FROZEN 区間・lib pin は不変。v31.6 の gap 抽出器訂正と同族の scaffolding 修正)。

## 5. 期末完全儀式 (v34.0-B)

共有部の変更 (lib 5 モジュール新設 + operational_net 増補) が第三十三期全体に
わたるため、期末に一度だけ全数再計算を行った:

```text
make suite-full OUT=results/v340_full_suite.txt JOBS=12   (壁時計 ~18.2 h)
→ 完全再計算 203 本: PASS 1427 / FAIL 2 — FAIL は v336 [S3] の壁時計バー
  (並列負荷で超過)。**壁時計を合否条件にしたのは「並列化で結果が変わらない」規約
  (PROMPT/4) に反する器械バグで、儀式が設計どおり検出した** → 計時を合否から除外
  (裁定・数値は不変)・儀式記録は保持
→ 末桁ドリフト検査 (台帳の期前後比較): **既存 183 非監査バイナリの PASS/FAIL は
  完全一致 (ドリフト 0 件)・監査層 12 本も一致** — 第三十三期の共有部変更
  (laboratory_interface / resource_profile / contextual_factorization /
  graded_recovery / structured_backend 新設・operational_net 増補) は既存物理に無波及
→ 台帳確定後の再検証 (v336 訂正の再走 + 監査層 + 引用集約):
  results/v340_full_suite_reverify.txt — 総計 PASS 1429 / FAIL 0
  (v33.0/v28.0 の「儀式 + 再検証」と同型の記録)
```

## 6. 第三十三期の統合 — 確定残高

**期テーゼ「可アクセス性は作用素単体の属性ではなく、系・制御器・測定器・資源・
誤差証明書の関係である。局所性は、証明付き laboratory interface が生成する
role-typed context atlas 上で整合する因子分解の、資源スケールにわたる安定な
同値類としてのみ識別される」— が、型 (禁止変換 12–21)・定理 (第三の no-go =
controller-free decomposition)・負制御 (entangler/tied/スカラー潰し/CZ 捻り/汚染
witness)・holdout (HOLD-9 20/20 満票 + 5 新計量) の全てで閉じた。**

| 版 | 成果 (確定) |
|---|---|
| v33.1 | **境界監査と型スコープ修復** — contexts 盲目性の機械実証・entangler 負制御 (選別の循環)・MarkedRecoveryInput (gens 別渡し廃止・構成時資格 3 種)・禁止変換 12 (代数的可換 ↛ 操作的両立)・13 (可換子 lane ↛ GKLS — Leibniz 導分証明書) |
| v33.2 | **Certified Laboratory Interface** — 宣言 ≠ 資格 (出自 3 証明書・sha256 結束・文字列 provenance 廃止)・tied no-go・**controller-free no-go (E3-A, 第三の no-go)**・role-typed 文脈 4 型 (joint measurability は可換性より広い)・禁止変換 14/15/16 |
| v33.3 | **Resource-Filtered OperationalNet** — profile: budget ↦ 読み・poset は barcode でない・昇格規則 chain ≥ 2 (transient 非昇格)・スカラー潰しの裁定反転・頂は経路を消す・禁止変換 17/18 |
| v33.4 | **Contextual factorization** — chart 局所復元 + overlap glue・**glue 定理 (= 直接大域復元と orbit 一致)**・cocycle 不整合 Abstain・複数 glue EquivClassOnly・witness 境界の両 lane 一致・禁止変換 19 |
| v33.5 | **Graded recovery の境界** — **Majorana locality ≠ Dirac locality**・O(2N) 不変・charge witness → 複素構造 J (J²=−I 4.4e-16) → U(N) を除くモード回復・縮退/汚染 witness の棄却・禁止変換 20 |
| v33.6 | **構造化スケーリング** — Pauli GF(2) backend (dense と全セル裁定一致・48 qubit 証明書を行列なしで)・Majorana quadratic backend (対応原理 2^{2m−1})・ScopeExceeded は正答・禁止変換 21 |
| Track X | **D2-R campaign layer (d2r-v1)** — 数値 kernel なしの公募受け皿 (事前登録 schema・validator fixture・AMBIGUITIES・OUTREACH/FUNNEL [数は実記録のみ]・「一件で足りる」の約束凍結) — v310 [R10] 常設監査 |
| v34.0 | **HOLD-9**: 凍結 (train 満票・lib pin) → 開封 **20/20 満票** — risk 0.000・recall 1.000/1.000・出自被覆 31/31・証人被覆 7/7・昇格/scope 違反/transient 昇格 0 (調整なし) |

**正直な残高 (変わらないもの)**:
- **bridge law 登録簿は全能力で空のまま** — HOLD-9 の満票をもってしても登録しない。
  blocker は独立外部再現 0 (R2/R3) — **Unit D2-R の公募が最優先** (campaign layer
  で受け皿を拡充済み・funnel は全段 0 の正直な記録)。
- PRED-019 未登録・自然の観測量の的中 0・`external_replications = 0`。
- E3 は A (controller-free no-go) と B (certified interface 下の復元) に分割して
  閉じたが、**「校正・独立 addressability・合成・資源上限」自体は依然として宣言
  入力**である (証明書の検証は機械・証明書の内容は実験者の申告 — 実測ノイズ下の
  区間較正は未走査)。
- toy スコープ: dense は dim ≤ 8 実演・structured lane は宣言 Pauli/quadratic 構造。

**未解決 (第三十四期への課題)**:
1. **外部独立再現 (最優先・据え置き)** — D2-R campaign の実施者募集 (funnel 稼働)。
2. 実測ノイズ下の証明書 (addressability/J/glue matching の区間資格と Straddled) —
   HOLD-9 は exact 証明書域だった。
3. structured lane の統合 — Pauli/quadratic backend 上の resource profile・atlas
   glue・graded witness (大型 complex mode 回復)。
4. 一般 GKLS 応答 (jump 表現 gauge・Kossakowski 同値類)・BCS 型 witness。
5. profile の関手性と安定性定理 (multiparameter persistence の語彙解禁の条件)。

## 7. 開発記録 (第三十三期)

- **監査が版を開いた**: 期の最初のコミットは新機能でなく境界監査 (v33.1) —
  「contexts が入力でなかった」ことの機械実証から型修復へ、という順序が
  PROMPT/14 の指示どおり機能した。
- **[H2] の区画誤記**: 開封 scaffolding の自己整合検査が自分の誤記 (8/12) を
  正しく FAIL して検出 — 検査が検査自身の入力ミスを捕まえる設計の実証。
- **壁時計バーの教訓**: v336 の計時合否は儀式の並列負荷 (JOBS=12) で初めて破れた —
  「並列化で結果が変わらない」規約 (PROMPT/4) は check 条件にも適用される。
  規模の主張は構成 (行列を生成しない) で担い、計時は合否にしない。
- **lib pin 方式**: HOLD-8 の「カーネル逐語コピー」から「lib モジュール sha256-16
  pin」へ — 器械が lib に住む期の凍結手法として機能した (開封時 6 pin 一致)。
