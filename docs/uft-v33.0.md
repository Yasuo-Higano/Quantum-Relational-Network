# QRN v33.0 — HOLD-8: 識別可能性境界の holdout と第三十二期統合

**Version**: v33.0-A (2026-08-01 — 凍結半) / v33.0-B (開封 + 期統合)
**Sim**: `sim/src/bin/v330a_hold8_freeze.rs` → `results/v330a_hold8_freeze.txt`
(3 検査 PASS + train 20 セル満票) / `sim/src/bin/v330b_hold8_open.rs` →
`results/v330b_hold8_open.txt` (開封 3 検査 PASS + holdout 本採点)
**位置づけ**: PROMPT/13 の最終版 — 第三十二期の全器械 (v32.1 SupportNoiseCertificate /
v32.2 OperationalNet / v32.3 recovery / v32.4 応答階層 / v32.5 hypergraph) を凍結し、
**factorization × interaction order × observation contract** の識別可能性境界を
新鮮 holdout で採点する。

---

## 1. HOLD-8 が採点するもの (HOLD-7 の後継)

「正解グラフを何本当てたか」ではなく**識別可能性の境界**:

```text
selective risk        = 0   (回答セルの誤り 0 — 強制回答も誤りに数える)
impossibility recall  = 1   (非識別セルの正しい棄却/同値類)
answerable recall     = 1   (回答可能セルは全て答える — coverage 単独は最大化しない)
```

## 2. 開封順序 (HOLD-5/6/7 と同一)

```text
v33.0-A (本コミット) = 生成器・採点器・バー・棄却規則の凍結
                      + SECRET コミットメント公表 + train 採点 (可視シード 33001)
  → v33.0-B = SECRET 開示・holdout 初生成・本採点 (調整なし) + 期末完全儀式
sha256(SECRET) = fb05a4a07dd7feef78e0bcecdaeeff933a656d36a38ceb26c71047002621c582
```

凍結カーネル (FROZEN-HOLD8 区間 — v330b と逐語一致を [H0] が照合) は第三十二期の
器械の逐語コピー: OperationalNet 構築 (可換子の区間証明書, τ = 1e-3)・recovery
決定手順 (成分 → 閉包 → 中心 → 三値裁定)・R⁽¹⁾/R⁽²⁾ 応答核・条件付き接ベクトル・
最終 gap 則・SupportNoiseCertificate (σ·z 対 窓ガード)・hyperedge 検出バー
(条件差 > 1e-6·max)・coherent 分離バー (0.05, 電流 + 実 hopping の 2 チャネル)。

**器械訂正 (設計走行で発見・v32.3 kernel に統一適用)**: 共役射影の数値塵
(‖PbP‖ ≈ 0 の候補が正規化されて基底に混入し restricted 次元を膨らませる) を
dust guard (閾値 1e-9) で除外 — v323 の既存 8 検査は全て不変 (exact 射影では
no-op)。v31.6 の gap 抽出器訂正と同じ「器械契約への昇格」型の修正。

## 3. セルクラス (8 種 20 セル — 隠しパラメータは重み・角度・置換・V・θ・σ)

| クラス | 内容 | 要求 |
|---|---|---|
| F1/F2 ×2 | site net / mode (DFT) net (隠し局所 U × 置換) | **Exact [2,2,2]** + 正しい gauge orbit |
| F3 ×1 | number operator のみ | **Abstain(InsufficientOperationalGenerators)** |
| F4 ×1 | 完全だが非互換な 2 net | **EquivalenceClassOnly** (強制一致 = FAIL) |
| F5 ×1 | 中心非自明 {X_a, Z_a, Z_b} | **SuperselectionSectors [(2,2),(2,2)]** (tensor 強制禁止) |
| F6 ×1 | odd 入力 + 偶双線形 (隠し係数) | ordinary **構成時拒否** + graded 化で [(4,1),(4,1)] |
| I1–I5 ×6 | quadratic (ring/chain 隠し) / t-V / 相関 hopping / pair hopping / 三体 ± null | 支持 欠0余0・**ŵ = t² (V 不可視)**・hyperedge 検出 (null 対で不発)・条件付き K(v) = \|t+vV\|² |
| I6/I7 ×4 | H↔−H 対・磁束 ±θ 対 | **密度 lane = 同値裁定 (分離主張 = FAIL)**・coherent lane = 分離 |
| M1/M3 ×2 | 変成対 (局所 U×置換)・基底の可逆再結合 | 同一 orbit (overlap 1)・recovery 不変 |
| M4/M5 ×2 | 可換子 margin 以下のノイズ・support margin 以下の弱辺 | **Abstain** (Straddled / InsufficientObservation) |

**train (シード 33001, 可視) は 20 セル満票**: selective risk 0.000・
impossibility recall 1.000・answerable recall 1.000。回答 14/14 (F1/F2/F5/F6・
I1–I5・I6c/I7c・M1/M3)・正棄却 6/6 (F3/F4・I6d/I7d・M4/M5)。

## 4. v33.0-B (開封) — 本採点の確定表

SECRET を開示 (`HOLD8-7cb8c1a325d20044e4f87f89c2bd7892` — [H0] が sha256 =
コミットメント fb05a4a0… と FROZEN-HOLD8 区間の逐語一致 [kernel sha 44d14839…] を
機械照合) し、holdout 20 セル (シード 18088044487616495343 = SECRET 導出) を
初生成・本採点した (調整なし):

| セル | 裁定 | 結果 |
|---|---|---|
| F1-site / F2-mode | Exact [2,2,2] + orbit overlap 1.000000000 | ✓✓ |
| F3-numberonly | Abstain(InsufficientOperationalGenerators) | ✓ 正棄却 |
| F4-incompatible | EquivalenceClassOnly (強制一致なし) | ✓ 正棄却 |
| F5-superselect | SuperselectionSectors [(2,2),(2,2)] | ✓ |
| F6-graded | ordinary 構成時拒否 + [(4,1),(4,1)] | ✓ |
| I1-quadratic-ring | 支持 欠0余0 (holdout 抽選は ring)・hyper なし | ✓ |
| I2-tv | 支持 欠0余0・ŵ = t² (V 不可視) | ✓ |
| I3-corrhop | hyper (0,1,2) 検出・K(v) = \|t+vV\|² | ✓ |
| I4-pairhop | 支持 {(0,1),(1,2)}・hyper なし | ✓ |
| I5a / I5b | hyper 検出 / null 対は不発 | ✓✓ |
| I6d / I7d (非識別) | 密度 lane = 同値裁定 (分離主張なし) | ✓✓ 正棄却 |
| I6c / I7c | coherent 2 チャネルが分離 | ✓✓ |
| M1 / M3 | 変成対 同一 orbit (1.000000000000)・再結合不変 | ✓✓ |
| M4 / M5 (非識別) | Abstain(Straddled / InsufficientObservation) | ✓✓ 正棄却 |

**確定: 20/20 満票 — selective risk 0.000・impossibility recall 1.000 (6/6)・
answerable recall 1.000 (14/14)・強制回答 0。** HOLD-7 で初通過した「棄却の採点」が、
因子分解 × 相互作用次数 × 観測契約の全軸で新鮮データを通過した。

## 5. 期末完全儀式 (v33.0-B)

共有部の変更 (lib.rs の `operational_net` 追加・readout_contract の後継注記) は
v32.2 に集約したため、期末に一度だけ全数再計算を行った:

```text
make suite-full OUT=results/v330_full_suite.txt JOBS=12   (壁時計 ~16.7 h)
→ 完全再計算 195 本: 物理・器械 194 本は全 PASS。監査層 v151_audit のみ 3 FAIL —
  儀式走行中の台帳が期末更新の中間状態 (v33.0-B 材料の登録前) にあったことを
  [3][9][10] が正しく検出した (監査が設計どおり機能した中間観測)
→ 台帳確定後の再検証 (監査層の再走 + 195 本引用集約):
  results/v330_full_suite_reverify.txt — 総計 PASS 1371 / FAIL 0
  (v28.0 の「儀式 + 再検証」と同型の記録)
```

**末桁ドリフト検査 (台帳の期前後比較)**: 旧台帳 (v32.0-B 時点) と新台帳 (195 本)
を突き合わせ、**既存 175 非監査バイナリの PASS/FAIL は完全一致 (ドリフト 0 件)** —
第三十二期の共有部変更 (operational_net 新設・readout_contract 後継注記) は既存物理に
無波及。差分は新規 8 本 (v321–v326, v330a/b) の +56 PASS と、監査層の設計上の変化
(v151 の claims 増加分) のみ。

## 6. 第三十二期の統合 — 確定残高

**期テーゼ: 「局所性は状態に宿るのではなく、操作的文脈の可換分解と、その上の
Liouvillian 応答の同値類としてのみ識別される」— が、型・定理・負制御・holdout の
全てで閉じた。**

| 版 | 成果 (確定) |
|---|---|
| v32.1 | **Unit D2 プロトコルの型修復と実行前反例** — `F}oXO` (凍結 gap 則が 11→13 辺の余剰・n=7 全数 22/853 故障 [欠落 0])・B3SupportMarginCertificate ⟺ exact (992 全数例外 0)・応答 lane が全反例修復・SupportNoiseCertificate (重みバー通過でも支持は落ちる)・報告契約を実 JSON Schema 化・superseded_before_external_run (黙って修正しない) |
| v32.2 | **OperationalNet + 第二の no-go** — 閉包は marking を消す (site/DFT 両族が M₈, su(2) overlap 0.618)・禁止変換 11・役割 4 型の分離・可換子の区間証明書 (跨ぎ Abstain)・JW 幾何誤読の型遮断 (odd は ordinary で K₃ 捏造) |
| v32.3 | **目標定理 B (marked recovery)** — 成分→閉包→中心の三値裁定・Exact [2,2,2]/[2,3]・**SuperselectionSectors (tensor 強制禁止・測定だけの軸 = 超選択ラベル・未 address = 多重度)**・gauge orbit 裁定 (site×DFT → EquivalenceClassOnly)・noise abstention・**パリティ超選択 [(4,1),(4,1)] の機械発見** |
| v32.4 | **Liouvillian 応答階層** — R⁽¹⁾/R⁽²⁾ 恒等式 (ε 非依存 9e-15)・一階×情報完全 = H を中心を除いて一意 (1.2e-16)・**R⁽²⁾ は H↔−H 不可視 (厳密)**・磁束は密度に不可視/coherent が分離・保存則和則・PSD・**v31.2 曲率則 = 本階層の対角特殊化** |
| v32.5 | **Interaction hypergraph** — H_S 直交分解 (Möbius = Pauli 支持 5.6e-16)・w_S 局所 unitary 不変・**中心化分離 (V·n₃h₁₂ = 二体⊕三体 等重み)**・**遷移率和則 K = Σ_{S⊇{i,j}} w_S/4 (「破れ」ではない)**・条件付き probe の次数分離・coherent 符号回復 |
| v32.6 | **VR exactness** — 離散円環 bar 定理 [1, ⌈n/3⌉) (n=4..30 全数)・規約の型 (VR_< は +1)・H2 persistence (sparse+clearing 25% 節約)・wedge-S² 遷移・**K3-holes 実測 5.00 の閉形式 retrodiction** |
| Track X | 外部公募の受け皿修復 — protocols/ 版分離・実 JSON Schema・v310 [R9] 常設 sha256 認証 (公募は v32.1 から継続中) |
| v33.0 | **HOLD-8**: 凍結 (train 20 セル満票) → 開封 **20/20 満票** — selective risk 0・impossibility recall 1.0 (6/6)・answerable recall 1.0 (14/14)・強制回答 0 (調整なし) |

**正直な残高 (変わらないもの)**:
- **bridge law 登録簿は全能力で空のまま** — HOLD-8 の満票をもってしても登録しない。
  blocker は独立外部再現 0 (R2/R3) — **Unit D2-R の公募が最優先** (v32.1 で受け皿は
  修復済み)。
- PRED-019 未登録・自然の観測量の的中 0・`external_replications = 0`。
- E3 は「操作的 fiber の機械化」まで — **どの operations が physically accessible
  かは依然入力** (完全な RelationalDecompositionGoal ではない)。
- toy スコープ: dim ≤ 8 (3 qubit/3 モード)・graded recovery は parity-even lane。

**未解決 (第三十三期への課題)**:
1. **外部独立再現 (最優先・据え置き)** — Unit D2-R の実施者募集。
2. 重みつき円環の実数 filtration exact bar (次の VR バー導出)・レンズ空間は
   torsion 不変量の能力分離後に。
3. graded lane の recovery (odd 演算子からのモード構造)・4 体以上の hypergraph 交差項。
4. OperationalNet の「物理的に accessible な操作」の出自 — E3 の残り半分。
5. 大型系 (dim > 8) への recovery のスケーリング (閉包成長の疎化)。

## 7. 開発記録 (第三十二期)

- **実行前反例がプロトコルを直した**: D2-v1 の「任意の連結グラフ」は外部走行 0 件の
  時点で反例が見つかり、版分離 + certificate scope に修復された。「主張域は器械が
  証明できる域に事前に切る」— v29.4a/K3-holes の教訓がプロトコル設計に到達した。
- **noise 裁定の 2 段化**: 重みバー (0.1) を通っても支持は落ちる (ガード跨ぎ) —
  SupportNoiseCertificate が観測契約の一部になった。
- **共役射影の数値塵** (HOLD-8 設計走行で発見): ‖PbP‖ ≈ 0 の候補が正規化されて
  基底に混入し restricted 次元を膨らませる — dust guard を v32.3 kernel に統一適用
  (既存検査は不変)。v31.6 の gap 抽出器訂正と同族の「器械契約への昇格」。
- **応答の符号規約**: R⁽¹⁾ の符号は発展の向き (e^{−iHt} vs e^{+iHt}) で反転する —
  v32.4 で Schrödinger 規約 (v31.2 と同一) に固定し、恒等式照合で機械確定した。
- **coherent 分離は 1 チャネルでは足りない**: 電流 J の応答は cosθ (偶) — 磁束対
  ±θ は J だけでは分離できず、実 hopping チャネル (sinθ 系) の併用で分離した
  (HOLD-8 設計走行の発見)。
