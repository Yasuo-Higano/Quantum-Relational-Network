# D2-R 実施パケット (配布用) — Unit D2-R: 応答 end-to-end の独立外部再現

**これは何か**: 本リポジトリの geometry 読み出し主張の唯一の外部検証単位
Unit D2-R を、第三者がそのまま実施できる形にまとめた案内。**有効な報告は
一件で足りる** (凍結された約束 — 二件要求へ変更しない)。失敗・不一致の報告も
成功と同じ形式で受理し公開する。

## 実施者に求めるもの (6 条件 — replications.yml)

1. different_author (本リポジトリの作者と別の人間/組織)
2. independent_repository
3. no_shared_numerical_kernel (共有コードなし — 同一 AI による再実装は独立でない)
4. protocol_frozen_before_run (実行前に凍結プロトコルへコミット)
5. commit_hash_recorded
6. result_including_failures_public

## 読むもの (この順)

| ファイル | 内容 |
|---|---|
| `reproducer/SPEC.md` | 再現単位の総則 (失敗も同じ形式で提出) |
| `reproducer/protocols/v32.1/` | D2-R プロトコル本体 (凍結・sha256 認証) |
| `reproducer/campaigns/d2r-v1/PREREGISTRATION.schema.json` | 事前登録の機械可読形式 |
| `reproducer/campaigns/d2r-v1/AMBIGUITIES.yml` | 既知の曖昧点 (凍結文は追記明確化のみ) |
| `reproducer/campaigns/d2r-v1/INVALID_REPORT_FIXTURES/` | 適合/不適合の較正例 |
| `paper/operational-core-spec.md` (OCS-1.0) | **paper-closed core spec** — ソース・出力値なしで operational core を独立実装するための唯一の入力 (sha256 凍結) |
| `reproducer/real_data/RECORDED_LANE.md` | 実装置の記録を持つ実施者向けの recorded lane |

## 提出方法

GitHub issue または PR で、事前登録 (schema 適合 JSON) → 実行 → 報告。
報告は成功/失敗を問わず `replications.yml` に登録される (外部条件を満たす
報告のみが external_replications を動かす)。

## 約束 (凍結)

- 有効報告は**一件で足りる**。
- 外部再現の成否を内部 holdout の採点に入れない。
- 報告数・接触数の水増しはしない (数は実記録のみ — OUTREACH_LEDGER)。
- プロトコルの曖昧さが見つかった場合は、結果に合わせず版を上げる。
