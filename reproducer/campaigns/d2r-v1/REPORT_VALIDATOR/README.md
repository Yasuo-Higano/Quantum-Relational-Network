# REPORT_VALIDATOR — D2-R 報告の機械検証手順 (d2r-v1)

**このディレクトリは手順書であって数値 kernel ではない。** 検証は書式・独立性・
版結束の審査であり、物理値の再計算を含まない (それは報告者の実装の仕事)。

## 検証手順 (順に全て)

1. **事前登録の適合**: 事前登録 JSON が
   [`../PREREGISTRATION.schema.json`](../PREREGISTRATION.schema.json) (draft
   2020-12) に適合すること。`planned.protocol_sha256_16` は凍結値
   `ce27f04b56303110` (d2-response.md) と一致すること — 版の結束。
2. **報告の適合**: 報告 JSON が凍結 schema
   [`../../../protocols/v32.1/unit-d-report.schema.json`](../../../protocols/v32.1/unit-d-report.schema.json)
   (sha256-16 = `f578816c54db3d23`) に適合すること。
   **正直な失敗は適合である** (負の結果のフィールドは schema が受理する)。
   **能力の水増しは不適合である** (例: D2-S 報告で geometry blocker 解除を主張・
   語彙外の capability 文字列)。
3. **独立性 6 条件**: replications.yml の 6 条件それぞれに公開証拠 (URL) が
   対応すること。走行前凍結は「実装 commit → 予告 → 走行」の時系列が第三者
   検証可能であること。
4. **版結束**: 参照している凍結ファイルの sha256-16 が v310 [R9] の認証値と
   一致すること:
   - `d2-response.md` = `ce27f04b56303110`
   - `d2-static.md` = `f858fa4bdeaa3554`
   - `unit-d-report.schema.json` = `f578816c54db3d23`
   - `unit-d-tolerances.yml` = `0ebc7098c6961355`
5. **数える範囲**: D2-R-FULL の必須セルを全て満たす有効報告のみが
   `external_replications` と [`../REPLICATION_FUNNEL.yml`](../REPLICATION_FUNNEL.yml)
   の `full_D2R_valid` を進める。D2-R-QUAL は funnel の練習欄までで、独立性欄には
   進めない。
6. **記録**: 判定 (適合・不適合とも) は理由つきで公開し、不適合の場合も報告自体は
   funnel (`report_submitted`) に数える — 提出の事実を消さない。

## 負制御 (validator 自身の較正)

[`../INVALID_REPORT_FIXTURES/`](../INVALID_REPORT_FIXTURES/) の fixture が期待通りに
適合/不適合と判定されることを v310 [R10] が常設監査する — validator が緩む・
締まりすぎる drift の検出。
