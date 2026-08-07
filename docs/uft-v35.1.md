# QRN v35.1 — 実験前仕様修復: RECORD v2・相関粒度型・status 一本化・outreach kit

**位置づけ**: PROMPT/16 (第三十五期「外部化と観測商」) §3–4。優先度 0 (器械保守) と
1-A (D2-R 外部化の配布準備)。**本版は instrument maintenance であり、物理的前進に
数えない** (PROMPT/16 §2「maintenance commit を physical advance と呼ばない」)。

**一次ソース**: `results/v351_record_v2.txt` (12 検査 PASS)。

---

## 1. なぜ「データ取得前」に直すのか

v34.6 の recorded lane には prose と schema の食い違いがあった:
`RECORDED_LANE.md` は「shot 配列**または counts** + 取得時刻」を許すが、
`RECORD.schema.json` (v1) は各 channel に `shots` を必須とし counts を表現できない。
また source/target/probe 符号/ε/発展時間/program hash/取得時刻 (チャネル別) が
型になっていなかった。

これを**実データ受領後に**直すと「結果を見た後の schema 変更」に見える —
事前登録の信頼はそこで失われる。よって v35.1 は、実データが 0 件である今
(REAL_DATA_LEDGER: recorded_runs = 0)、v1 を保存したまま v2 を別 hash で凍結する。
取得後の v2 変更は `RECORD_V2/MANIFEST.sha256` の drift として `v351_record_v2`
が機械検出する ([P2])。

## 2. RECORD v2 — 三段 schema (凍結)

単一の可変 JSON ではなく、役割の異なる三段に分離した
(`reproducer/real_data/RECORD_V2/`):

| 段 | ファイル | 内容 |
|---|---|---|
| A | `preregistration-v2.schema.json` | 観測契約 hash・解析 commit・calibration/held-out チャネル・予測コミットメント sha256・採点規則・misspecification gates (drift α・相関 α・二 ε の signed-linearity)・topology+mapping コミットメント |
| B | `acquisition-v2.schema.json` | provenance (`recorded_experimental` const)・**evolution_provenance** (native_analog / compiled_digital / emulator)・チャネルごとの source/target/probe_sign/ε/evolution_time/**evolution_family_id**/job_id/**program_hash**/compiler 版/calibration snapshot hash/**used_in_calibration (必須 boolean)**/raw_artifact_sha256・データ本体は相関粒度 3 型の **oneOf** |
| C | `disclosure-v2.schema.json` | 予測 plaintext・hidden topology・label mapping の開示・コミットメント照合 (boolean のまま記録)・採点 (hits+misses+straddled = total の会計)・失敗込み final_verdict |

設計原則:

- **全 object に `additionalProperties: false`** — 「採点後の追記」の型レベル禁止
  (負 fixture `acquisition_invalid_extra_field.json` が較正 — [P4])。
- **三段は sha256 で結束** — acquisition は prereg の実 hash を、disclosure は
  prereg + acquisition の実 hash を持つ。fixture は実 hash 連鎖で較正済み ([P4])。
- **± probe 対の evolution program 同一性** — 同一 (source, target, ε, t, family) で
  符号だけ異なる 2 チャネルは program_hash が一致しなければ不適合
  (PROMPT/16 §5.3 の器械化 — compiled lane で ± が別プログラムにコンパイルされる
  事故を構成時に拒否)。
- **evolution_provenance の 3 値** — `NativeAnalogEvolution` 以外から native
  hardware Hamiltonian topology を主張しない (§8 の型準備)。
- 凍結 hash: `RECORD_V2/MANIFEST.sha256`
  (prereg fe931d97… / acquisition 0b9bc397… / disclosure 262d81d0…)。
  v1 は不変 (sha256 = 780cd728… を [P1] がピン)。

## 3. 相関粒度型 (`sim/src/record_v2.rs`) — 禁止変換 30/31

HOLD-10 の教訓: 同一周辺分布の持続的 Markov 鎖は split-half を通過する
(v343 [F6-29]) ため、遷移数ゲートが必須になった。ところが**順序を失った
aggregate counts では遷移数ゲートが原理的に実行できない**。検査の可能性は
データの粒度の属性であり、後から昇格できない:

| 粒度 | 実行可能な検査 | 最良到達裁定 |
|---|---|---|
| `OrderedShots` | split-half (α/2) + 遷移数 (α/2) | `IidConsistent` (証明書) |
| `TimestampedBatches` | batch 間 drift (対ごと CP, α/2) + overdispersion (Pearson T vs χ²_{m−1}, α/2) | `CorrelationUnresolved` |
| `AggregateCounts` | なし | `CorrelationUnassessed` |

- **禁止変換 30**: `AggregateCounts` ↛ `IidCertificate`。
- **禁止変換 31**: `TimestampedBatches` ↛ `IidCertificate` — batch 間検査を全て
  通過しても「batch 内の serial correlation は観測不能」なので証明書は出ない。
- `IidCertificate` は private フィールドを持ち、`assess()` が OrderedShots で
  全ゲート通過したときのみ構成する (v33.2「門は較正」の規律)。
- OrderedShots の登録最小長 40 (split-half 各 20) — これ未満は
  `CorrelationUnresolved` (検査できないのに資格を出さない)。
- χ²_{m−1} は**近似であることを登録** (Pearson) — 較正 [P5d] が固定シードで
  偽検出率 ≤ 2α を機械検査する (実測 0.0070, α = 0.01)。

**較正結果** ([P5], 固定シード):

- 8 セル裁定 (iid/Markov/drift/短列/良性バッチ/batch drift/過分散/aggregate)
  全て期待どおり。aggregate は n = 10⁶ でも `CorrelationUnassessed`。
- ordered level: iid 偽検出 0/2000 (CP ゲートの保守性)。
- ordered power: 持続 Markov (stay 0.9, n = 800) の非 IidConsistent 率
  **200/200** (遷移数 175 + split-half 25 — 持続鎖は drift ゲートでも正しく
  棄却され得る。どちらのモードでも証明書は発行されない)。
- batches power: 交互 p = 0.38/0.62 (m = 12) の検出 200/200。

**発見 (較正が教えたこと)**: 初版バーは「遷移数ゲート単独で ≥ 0.95」としたが
168/200 で FAIL — 片側に張り付いた持続鎖は周辺分布ごと変わるため split-half が
先に検出し、また張り付き軌道 (p̂ ≈ 0.05) は q = 2p(1−p) ≈ 遷移率となり iid と
統計的に区別しにくい。**較正の主旨は「iid 証明書を誤発行しない」**なので、
検出 = 非 IidConsistent と登録し直した (n = 800 で誤資格 0/200)。

## 4. FollowUp status の一本化 (別リポジトリの統治修復)

FollowUp final report は Phase 7/8 完了 (2026-08-02) を記録する一方、
`docs/replication_plan.md` は「Phase 7 holdout has not been opened」と表示し
続けていた。**単一の機械可読 source of truth が無いことによる統治不全**であり、
QRN-FollowUp リポジトリに以下を実装した (commit `1435c1e` — 科学的内容・凍結
規則・事前登録された plan 本体には触れない governance 修復):

- `FINAL_STATUS.json` — 唯一の状態ソース (phase・verdict・holdout 消費済み・
  shared_human_operator: true / organizationally_external: false)。
- `tools/sync_status.jl` (std lib のみ) — README / replication_plan /
  final_report のマーカー区間を生成。README と final_report は**内容 byte 不変**
  でマーカー化のみ (plan の陳腐化 status だけが更新された)。
- `.github/workflows/status-check.yml` — drift で CI が落ちる。

主リポジトリ側は replications.yml の REP-001 (partially_replicated /
human_operator: shared / external_replications = 0) との整合を [P7] が常設監査。

## 5. D2-R outreach kit (`reproducer/campaigns/d2r-v1-outreach/`)

D2-R の blocker は validator でも定理でもなく**配布**である (funnel 全段 0)。
kit は配布実験の設計を**反応を見る前に**固定する:

- **3 cohort × 各 10 件 = 30 件** (Hamiltonian/Lindbladian learning・quantum
  network topology・TDA/reproducibility/formal methods)。
- 配布物 4 点凍結 (`MANIFEST.sha256`): OUTREACH_KIT / CHALLENGE_SUMMARY (一枚もの・
  英語・先行研究に対する位置づけの正直な記載つき) / QUAL_CHECKLIST / EMAIL_TEMPLATE。
  **数値 kernel・実装骨格・翻訳コードは渡さない** ([P6] がコード混入も監査)。
- **停止規則の事前登録**: 30 件で preregistered = 0 なら
  `protocol_externalizability: failed_at_current_burden` を記録 (物理の反証では
  ない)。その後にのみ D2R-v2 を設計し、v1 は変更しない。
- d2r-v1 の凍結物は逐語不変 — MANIFEST 自体の sha256 (af967b3e…) を [P6] がピン。
  OUTREACH_LEDGER への記録は**実接触のみ** (台帳行の正当な更新時は MANIFEST の
  該当行を同一コミットで更新)。
- **送信は人間の作業** — 本版が担うのは配布物の凍結と記録の器械化まで。
  funnel は全段 0 のまま (この行自体が正直な残高の記録)。

## 6. 検査一覧 (v351_record_v2 — 12 PASS)

| 検査 | 内容 | 結果 |
|---|---|---|
| [P0] | record_v2 自己検証 (χ² 分位の教科書値・粒度到達集合) | PASS |
| [P1] | v1 凍結ピン + lane 文書の v2 節アンカー | PASS |
| [P2] | v2 MANIFEST (三段 schema の hash 一致 — 取得後変更の drift 検出) | PASS |
| [P3] | schema 構造 (全 object additionalProperties:false・oneOf 3 粒度・必須フィールド) | PASS |
| [P4] | fixture 正 4/負 3・三段 sha256 結束・コミットメント照合・held-out 整合・± 対規則 | PASS |
| [P5a] | 裁定較正 8 セル (IidConsistent は OrderedShots のみ) | PASS |
| [P5b] | ordered level (iid 偽検出 0/2000 ≤ 1.5α) | PASS |
| [P5c] | ordered power (Markov stay 0.9: 非 IidConsistent 200/200) | PASS |
| [P5d] | batches level (偽検出 0.0070 ≤ 2α — χ² 近似の登録 slack) | PASS |
| [P5e] | batches power (交互 0.38/0.62: 200/200) | PASS |
| [P6] | outreach kit (d2r-v1 不変ピン・kit MANIFEST・停止規則・コードなし) | PASS |
| [P7] | FollowUp 整合 (REP-001 = partially_replicated・external 0 維持) | PASS |

## 7. 限界と非主張

- 本版は**器械保守**である。recorded_runs = 0・external_replications = 0・
  bridge law 空・PRED-019 未登録は全て不変 — どれも代用しない。
- 相関ゲートは iid 契約の**破れの検出器**であり、通過は iid の証明ではない
  (検出力は n と代替仮説の関数 — [P5c] の解剖を §3 に記録)。
- overdispersion の χ² は近似 (登録済み)。厳密化 (正確な条件付き分布・
  bounded martingale lane) は将来課題。
- outreach の実施 (送信) と応答は人間と第三者の領分 — 台帳が空である限り
  funnel は 0 のまま報告する。
