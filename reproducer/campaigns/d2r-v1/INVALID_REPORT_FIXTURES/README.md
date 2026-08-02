# INVALID_REPORT_FIXTURES — validator の負制御 fixture (d2r-v1)

事前登録 schema に対する正例 1 + 負例 2。v310 [R10] が「正例は適合・負例は
それぞれ宣言された理由で不適合」を常設監査する (validator の drift 検出)。

| fixture | 期待 | 理由 |
|---|---|---|
| `prereg_valid_minimal.json` | 適合 | 必須フィールド完備・凍結 sha 結束・独立性宣言 |
| `prereg_invalid_missing_independence.json` | 不適合 | `independence` ブロック欠落 (6 条件の宣言なし) |
| `prereg_invalid_capability_inflation.json` | 不適合 | `unit` が語彙外 (`D2-S-geometry-unlock` — D2-S は blocker を解除しない) + 凍結 sha 不一致 |

fixture は**架空の記入例**であり実在の再現者ではない (author は placeholder)。
