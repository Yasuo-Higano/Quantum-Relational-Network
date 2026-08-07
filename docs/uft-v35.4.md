# QRN v35.4 — FollowUp adversary packet A/B と OCS-2.0 の gating

**位置づけ**: PROMPT/16 (第三十五期「外部化と観測商」) §9・優先度 4。FollowUp を
外部再現としてではなく **paper-closed theorem adversary** として使う体制の整備。
論文更新を実装 commit ごとではなく atomic packet 単位に固定する。

**一次ソース**: `results/v354_packets.txt` (6 検査 PASS)。

---

## 1. FollowUp の正しい使い方 (方針 §6 の器械化)

REP-001 (v34.1 受理) は cross-model clean-room だが **shared human operator**
であり、組織的外部再現ではない — これは FollowUp の final report 自身が明示して
いる。よって FollowUp の価値は「外部性」ではなく**敵対的検証**にある:
独立導出・反例探索・scope narrowing を要求し、`RefutedAsStated` と
`Inconclusive` を一級の成果として受理する。

## 2. Packet A — Open Response (`paper/packet-a-open-response.md`)

sha256 = e23e82db… (v354 [K1] が pin)。9 claims:

- **OQ-0** (較正課題): covariance ODE の X, Y の導出 — **規約 (どちらの散逸行列が
  転置で入るか) を意図的に渡さない** (v35.2 [G1] で我々自身が転置を誤った既知の
  失敗モードであり、独立導出の試金石)。
- **OQ-1..OQ-6**: GQF-1..6 の paper-closed 版 (status / domain / 量化子 /
  入出力型 / gauge / falsifier / forbidden interpretation を各 claim に付与。
  反例対は構成のみ与え、値は adversary が計算する)。
- **OQ-7** (open): gauge 族 {local phase, global frequency, 複素共役} は
  full-time 契約の観測同値の全てか — 完全分類の攻撃依頼。
- **OQ-8** (open): full-time 契約での Hamiltonian support 識別可能性 —
  **OQ-3 の scope discipline を最初に攻撃させる** (方針 §2.3 の指示)。

## 3. Packet B — Resource Profile (`paper/packet-b-resource-profile.md`)

sha256 = 76d7267c… (v354 [K1] が pin)。6 claims: RPF-1..5 の paper-closed 版 +
**RP-2q** (累積 master data の restriction だけでは nesting が出ないことの独立
検証依頼) + **RP-6** (anytime-valid 系など代替 nesting 構成の設計 open 問題 —
OCS-2.0 の材料)。

## 4. Packet C — OCS-2.0 は gating (凍結しない)

PROMPT/16 §9 のとおり、**OCS-2.0 は Packet A/B の反証と修正を受理した後にだけ**
作る。現時点で `paper/` に ocs-2.0 ファイルは存在しない — v354 [K5] が
「存在しないこと」を常設監査する (反証受理前の凍結は規律違反として FAIL する)。
OCS-1.0 は不変。盛り込む予定の内容 (PROMPT/16 §9 の凍結リスト): sample space・
iid/相関仮定・信頼集合構成・同時 α 配分・misspecification gates・
optional-stopping policy・exact verdict map・resource refinement maps
(v35.3 の intersection 構成が材料)・normative constructions・falsifiers・
forbidden interpretations。記録済み出力値は書かない。

## 5. FollowUp への配布 (staging)

packet A/B は QRN-FollowUp リポジトリの `papers/` に byte 同一で staging した
(FollowUp の AGENTS.md が許可する唯一の入力経路)。**clean-room rerun の実施は
人間 (ユーザー) の作業** — 別モデルの adversary agent が packet だけから独立
導出・反例探索を行い、報告は本リポジトリの受理判定 (RefutedAsStated 含む) で
記録する。原実装・結果値・Lean 証明は渡していない (paper-closed — v354 [K3])。

## 6. 検査一覧 (v354_packets — 6 PASS)

[K1] packet 実在 + sha256 pin / [K2] claim ID 全数 + 必須フィールド /
[K3] paper-closed (内部参照なし・expected_outputs_in_packet: false) /
[K4] adversary 規律 (RefutedAsStated/Inconclusive 受理・OQ-0 規約非開示) /
[K5] OCS-2.0 gating (ファイル不在 + 宣言) / [K6] scope discipline
(no-go の契約 scope 明示・full-time 問題は open)。

## 7. 限界と非主張

- packet の発行は**検証の依頼**であって検証ではない。GQF/RPF の no-go を
  完成条件 3 に数える資格は、FollowUp の反例探索を通過するまで発生しない。
- FollowUp が RefutedAsStated を返した場合、該当 claim は v34.1 の Yukawa
  erratum と同じ手続きで受理・訂正する (結果に合わせて packet を変えない —
  版を上げる)。
- 本版も instrument/paper 整備であり、実データ・外部再現・自然の的中は
  増えていない。
