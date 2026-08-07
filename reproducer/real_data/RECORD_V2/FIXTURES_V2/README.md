# FIXTURES_V2 — RECORD v2 validator 較正用 (実記録ではない)

ここにあるファイルは **fixture** であり、実験装置の実記録ではない。
`REAL_DATA_LEDGER.yml` の entries に数えることは禁止 (捏造の禁止 — v1 FIXTURES と同じ規約)。
hash 類 (program_hash / raw_artifact_sha256 等) は `fixture-...` 文字列の実 sha256 を
形式充足のために置いたダミーである。`v351_record_v2` が正負判定を較正する。

## 正例 (4)

- `prereg_valid.json` — A 段。prediction_commitment_sha256 は
  disclosure の `disclosed_prediction` 文字列の実 sha256 (三段結束の較正)。
- `acquisition_ordered_valid.json` — B 段・`ordered_shots` 粒度 (2 チャネル,
  ± probe 対が同一 program_hash / evolution_family_id)。
  preregistration_sha256 = prereg_valid.json の実 sha256。
- `acquisition_batches_valid.json` — B 段・`timestamped_batches` 粒度
  (4 batch, drift/overdispersion なし → 到達点は CorrelationUnresolved)。
- `disclosure_valid.json` — C 段。topology+mapping コミットメントの結合規則は
  sha256(disclosed_topology + "\n" + disclosed_label_mapping)。

## 負例 (3)

- `acquisition_invalid_extra_field.json` — トップレベルに `post_hoc_note`
  (additionalProperties:false の拒否対象 — 「採点後の追記」の型レベル禁止)。
- `acquisition_invalid_shots_and_counts.json` — data に shots と counts が同居
  (oneOf のどの粒度にも適合しない)。
- `prereg_invalid_no_heldout.json` — held_out_channels 欠落
  (held-out なしの事前登録は事前登録ではない)。
