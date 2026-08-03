# Recorded Experimental Lane (v34.6, PROMPT/15 §7)

**目的**: synthetic shot noise を「実測ノイズ」と呼ばないための構造的分離。
本 lane は**実際に物理装置から得た repeated finite-shot record** の受け皿であり、
synthetic coverage lane (HOLD 系・v34.5 robust atlas) とは別の型を持つ。

## 二 lane の分離 (凍結)

| | Synthetic coverage lane | Recorded experimental lane |
|---|---|---|
| 真の interface | **既知** (生成器が作る) | **未知** (latent) |
| 採点できるもの | 信頼集合の coverage を直接採点・adversarial boundary / model misspecification を生成可能 | model fit・drift 検査・予測 calibration・**未使用チャネル予測** |
| データの名前 | synthetic (これを実測と呼ぶことは禁止) | recorded experimental |
| 型 | `DataProvenance::SyntheticCoverage` | `DataProvenance::RecordedExperimental` |

二つの provenance の間に変換は存在しない (synthetic ↛ experimental)。
provenance を偽った提出はスキーマ段階で不適合 (v346 [L2] の負 fixture)。

## 提出物 (RECORD.schema.json 準拠)

1. **記録**: チャネルごとの shot 配列 (0/1) または counts + 取得時刻。
2. **事前登録**: 未使用応答チャネルの予測の sha256 コミットメント —
   **データ開示より前に** commit されたことが git 履歴等で検証可能であること。
3. **vendor topology commitment**: 装置の結線・トポロジ情報の sha256 —
   freeze 後に開示 (予測が topology を先に知らないことの証明)。
4. **operator 宣言**: 氏名/組織・独立性 (organizationally external か否か)。
5. **公開**: 失敗を含む全結果の公開に同意 (D2-R と同じ約束)。

## 採点 (v346 が器械化する規則 — 凍結)

- **契約検査**: 記録が登録ノイズモデル (iid) と整合するか —
  split-half の Clopper–Pearson 区間 (Bonferroni α/2m) が **disjoint なチャネルが
  あれば OutOfDomain** (drift の検出は失敗ではなく正答 — 禁止変換 25/29 の運用形。
  iid が破れたデータに iid 証明書を発行しない)。
- **予測採点**: 事前登録された未使用チャネル予測を、開示後の記録の
  信頼区間と照合 (RobustExact 意味論 — 区間が予測バーを跨げば Straddled)。
- **selective risk と coverage の区別** (v34.3): 観測誤り 0 は母集団リスク 0 では
  ない — 回答セル数 n に対し片側 95% 上限 1 − 0.05^{1/n} を常に併記する。

## 完成条件への接続 (v35.0 の科学的完成条件 — PROMPT/15 §7)

次の**いずれか**が無い限り、期末表現は「instrumental closure」に留める:
1. externally operated D2-R report (reproducer/D2R_PACKET.md)
2. 実データ上の事前登録された未使用応答チャネル予測の的中 (本 lane)
3. 広い観測契約族を排除する新しい厳密 no-go

## 正直な台帳

`REAL_DATA_LEDGER.yml` — recorded_runs / externally_operated_runs は
**実記録のみ** (捏造・シミュレーションでの水増しは禁止)。現在どちらも 0。
