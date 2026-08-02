# D2-R 外部再現 campaign (d2r-v1) — protocol-only replication note

**これは外部再現の「受け皿」であって数値 kernel ではない。** 本ディレクトリは
コードを一切含まない (v310 [R10] が機械監査)。再現対象の凍結プロトコルは
[`reproducer/protocols/v32.1/d2-response.md`](../../protocols/v32.1/d2-response.md)
(sha256-16 = `ce27f04b56303110` — v310 [R9] が常設認証) であり、**逐語凍結のまま
変更しない**。

## なぜ D2-R か

`external_replications = 0` が本プロジェクトの最大の未解決残高である。D2-R
(応答 end-to-end) は **`spatial_topology_given_factorization` 能力の独立性 blocker
を解除する唯一の単位** — ただしそれだけであり、**bridge law・PRED-019・自然の
観測量は D2-R の成功でも解除されない** (別問題)。

## 数える条件 (6 つ全て必須 — replications.yml と同一)

1. 別作者 (同一作者・同一 AI の別言語再実装は数えない)
2. 別 repository
3. 数値 kernel 非共有 (本リポジトリのコードを読んでもよいが、移植・翻訳は不可)
4. 走行前凍結 (実装 commit → 予告 → 走行の順序が公開記録で確認できる)
5. commit 記録 (時系列が第三者検証可能)
6. **失敗を含む公開** (負の結果も報告する — 正直な失敗は schema 適合である)

## 2 段の参加区分

- **D2-R-QUAL** (練習): schema・hash・負制御の練習走行。**外部再現に数えない**。
  参加障壁を下げるための区分であり、資格・独立性の審査は行わない。
- **D2-R-FULL**: 凍結プロトコルの必須セルを全て満たす。**これだけを数える**。

## 凍結された約束 (事後変更の禁止)

- **有効な報告は一件で足りる** — 有効な一件が来た場合、対応する claim scope の
  独立性を更新する。**後から二件要求へ変更しない**。
- 判定規則は走行前凍結の schema
  ([`unit-d-report.schema.json`](../../protocols/v32.1/unit-d-report.schema.json))
  と本 campaign の validator 手順のみ — 報告受領後に規則を追加しない。
- **外部再現の成否を内部 holdout (HOLD-9 以降) のセルに入れない** — 内部 SECRET
  から作ったセルは外部独立性ではない。`external_replications` は第三者からの
  有効報告のみで変わる。
- 本 campaign は**プロトコルが外部化可能かを測るメタ実験**でもある —
  [`REPLICATION_FUNNEL.yml`](REPLICATION_FUNNEL.yml) が段階別の実数を記録する
  (物理証拠ではない)。

## 手順

1. **事前登録**: [`PREREGISTRATION.schema.json`](PREREGISTRATION.schema.json) に
   適合する JSON を自分の repository に commit し、issue 等で予告する。
2. **実装**: 凍結プロトコル文書だけから独立に実装する (質問は AMBIGUITIES へ —
   回答は凍結文の書き換えではなく追記の明確化として公開される)。
3. **走行 → 報告**: 結果 (成功・失敗とも) を report schema に適合する JSON で公開。
4. **検証**: [`REPORT_VALIDATOR/README.md`](REPORT_VALIDATOR/README.md) の手順で
   機械検証される — 能力の水増し (D2-S で geometry 解除を主張する等) は不適合。

疑義・曖昧さの発見はそれ自体が寄与である ([`AMBIGUITIES.yml`](AMBIGUITIES.yml))。
