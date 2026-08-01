# 再現の帰結 — 単位 ↔ 主張の対応表 (v27.4 凍結)

再現の成否は本表に従って claims.yml に反映する。**成功しても失敗しても登録する** —
失敗の隠蔽は本プログラムの反証可能性を破壊する。

## 成功 (6 条件を満たす独立外部再現) の場合

replications.yml に登録し、該当主張の `independence` を `independent_author` に、
関係する箇所の `evidence_kind` を `external_replication` に**意識的に**更新する。
このとき v271_core_audit の R2/R3 期待値 (現在 0 件) の更新を同一コミットで行う —
監査の期待値を書き換えずに数だけ増やすことはできない (それが設計である)。

| 単位 | 昇格対象の主張 |
|---|---|
| A | QRN-GAUGE-003 (v3.1 全数探索の最小・一意), QRN-GAUGE-008 (独立実装照合), QRN-GAUGE-011 (Lean 定理化) |
| B | QRN-GRAV-030 (1D 質量関数への厳密還元), QRN-GRAV-031 (閉形式同定), QRN-GRAV-032 (区間証明) |
| C | QRN-GRAV-044 (解析 oracle), QRN-GRAV-049 (PRED-016 universality), QRN-GRAV-048 (和則) |

これが成立した時点で「独立外部再現 0」の残高表示 (README・claims.yml・
core.schema.yml) を更新する。**それまでは 0 のまま維持する。**

## 失敗 (バー外の不一致) の場合

1. 不一致報告を replications.yml に **status: failed で登録** (削除しない)。
2. 仕様の曖昧さ・正準化規約の差など「実装等価性の問題」をまず切り分ける
   (単位 A はハッシュでなく多重項集合で判定してよい — SPEC.md)。
3. 切り分け後も不一致が残る場合、上表の主張を**降格** (C2 → 要再検証) し、
   README の到達点行に不一致を明記する。FAL-SUITE の発火に相当する扱い。

## 単位 D (Geometry Identifiability — v31.0/PROMPT-12 で追加・v32.1 で D2 を版分離, 詳細は UNIT_D.md)

| 単位 | 対象主張 (claims.yml) | 効力 |
|---|---|---|
| D1 | QRN-BRIDGE-013, QRN-BRIDGE-017 | 数値再現のみ (blocker は残る) |
| D2-S | QRN-BRIDGE-021 | 静的 lane の scoped 再現 (StableEstimate — blocker は解除しない) |
| D2-R | QRN-BRIDGE-013, QRN-BRIDGE-017 | **geometry 能力の blocker 解除候補 (これのみ)** |
| D3 | QRN-BRIDGE-012, QRN-BRIDGE-015 | 負定理の独立確認 |
| D4 | QRN-BRIDGE-013, QRN-BRIDGE-016 | 応答法則の独立確認 |

D2-v1 (v31.7 の「任意の連結グラフ」) は外部走行前の反例 (graph6 F}oXO) により
`superseded_before_external_run` — 版台帳は
[protocols/v32.1/protocol-index.yml](protocols/v32.1/protocol-index.yml)。
gauge 単位 A/B/C の成功は geometry 能力を解除しない (capability scoped —
replications.yml v31.0 拡張)。

## 注意

- 同一作者・同一 AI・同一数値カーネルの別言語実装は algorithmic diversity
  (claims.yml の independence: algorithmically_diverse / same_author_clean_room)
  であり、本表の昇格を発動しない。
- 部分的成功 (単位 A のみ等) は単位ごとに扱う。
