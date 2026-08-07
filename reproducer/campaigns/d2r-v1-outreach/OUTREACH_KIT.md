# D2R-v1 OUTREACH KIT — 固定 campaign の設計 (v35.1, PROMPT/16 §4)

**位置づけ**: D2-R の blocker はもはや validator や定理ではなく、**実際に第三者へ
届けて参加者を得ること**である (funnel 全段 0 の正直な記録)。本 kit は配布実験の
設計文書であり、**d2r-v1 の protocol・schema・判定規則には一切触れない**
(凍結ピン: `sha256(campaigns/d2r-v1/MANIFEST.sha256) =
af967b3eb9a34511bd93785c05b4b6ca1899cd9459f0f9e2680219767e22de28` —
`v351_record_v2` が不変を機械監査する)。

## 固定 campaign (事前登録された設計 — 反応を見て変えない)

- **3 cohort × 各 10 件 = 合計 30 件**の適格な直接 outreach を記録する。
  - **Cohort H**: Hamiltonian / Lindbladian learning (量子系同定・トモグラフィ・
    短時間応答反転を扱う研究者/グループ)
  - **Cohort N**: quantum network topology (局所測定からの結線推定・量子ネット
    ワーク実験・デバイス較正)
  - **Cohort F**: TDA / reproducibility / formal methods (persistence・再現性
    工学・定理証明器による検証)
- **適格条件**: (i) 上記分野の公開研究実績がある個人/グループ、(ii) 連絡経路が
  公開されている (公開メール・issue・フォーラム)、(iii) 本プロジェクトの作者と
  組織的に独立 (different_author を満たし得る)。
- **配布物は以下に限定** (数値 kernel・実装骨格・翻訳コードは渡さない):
  1. `CHALLENGE_SUMMARY.md` (一枚もの)
  2. 凍結 protocol への参照 (`reproducer/D2R_PACKET.md` と
     `reproducer/campaigns/d2r-v1/` — リポジトリ URL + sha256 ピン)
  3. 事前登録 schema (`PREREGISTRATION.schema.json`) と較正 fixture
  4. QUAL 用チェックリスト (`QUAL_CHECKLIST.md`)
- **質問対応**: 全ての問い合わせは公開 ambiguity ledger
  (`d2r-v1/AMBIGUITIES.yml`) に集約する。凍結文の変更ではなく
  **non-normative clarification (追記明確化) のみ**。
- **記録**: 実際に行った接触のみを `d2r-v1/OUTREACH_LEDGER.yml` に記録する
  (id / date / channel / target_kind [cohort を含める] / response)。
  個人特定情報は書かない。**予定・シミュレーション・水増しの記載は禁止**。
  台帳が正当に成長したときは `MANIFEST.sha256` の該当行を同一コミットで更新する
  (台帳以外の行の変更は凍結違反 — v310 [R10] が drift を検出する)。

## 事前登録された停止規則 (結果を見てから決めない)

- **30 件の記録完了時点で preregistered = 0 の場合**:
  `REPLICATION_FUNNEL.yml` の dropout_reasons に
  `protocol_externalizability: failed_at_current_burden`
  を記録する。**これは物理や定理の反証ではない** — 現行の実施負担では
  プロトコルが外部化できなかった、というメタ実験の結果である。
- その後にのみ D2R-v2 (負担軽減版) を設計してよい。**v1 の条件は変更しない**
  (有効報告一件で足りる・判定規則の事後変更なし・成功で解除されるのは
  `spatial_topology_given_factorization` の外部独立性のみ)。
- 将来の v2 では独立性を「同じ AI か否か」の一軸ではなく、human/org・source
  access・shared code/kernel・shared prompt/workspace・model/revision・secret
  provenance の複数軸で記録し、`external verification` と `source-blind
  clean-room` を分離する (現行 v1 の判定は遡及変更しない)。

## 実施上の注意

- outreach の送信自体は人間 (リポジトリ管理者) の作業である。本リポジトリが
  担うのは配布物の凍結と記録の器械化まで。
- 位置づけの正直さ: 「局所測定から network topology を推論する」こと自体には
  先行研究がある。QRN の新規性は **signed curvature の厳密則・
  given-factorization scope・gauge-invariant block weight・no-go 群・
  fail-closed 有限データ証明書**に置く (最初の topology inference と主張しない)。
