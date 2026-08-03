# FIXTURES — scorer/schema 較正用 (実記録ではない)

ここにあるファイルは **fixture** であり、実験装置の実記録ではない。
`REAL_DATA_LEDGER.yml` の entries に数えることは禁止 (捏造の禁止)。
v346_real_data_lane が「台帳 0 件」と「fixture の較正結果」を独立に検査する。

- `example_stationary.json` — 適合 + drift gate 通過 (split-half 6/20 vs 6/20)
  + 事前登録コミットメントの開示一致 (sha256 照合)
- `example_drifting.json` — 適合だが drift gate で **OutOfDomain**
  (split-half 2/20 vs 14/20 — iid 契約の破れの正検出)
- `invalid_provenance_synthetic.json` — **不適合** (synthetic を本 lane に提出)
- `invalid_missing_prereg.json` — **不適合** (事前登録コミットメント欠落)
