# QRN v27.4 — 外部再現単位の公開: reproducer/ と replications.yml

**Version**: v27.4
**Date**: 2026-07-29
**Sim**: `sim/src/bin/v274_reproducer.rs` → `results/v274_reproducer.txt` (7 検査 PASS)
**新規**: [reproducer/](../reproducer/SPEC.md) (SPEC・INPUTS×3・EXPECTED_SCHEMA・
TOLERANCES・CLAIMS・NO_SHARED_CODE) + [replications.yml](../replications.yml)
**位置づけ**: PROMPT/10 §6。最大の不足は内部計算量ではなく**外部に固定された
再現単位と独立の検証** (PROMPT/9 の結論) — その不足を埋める形を凍結する。
物理 run なし。

---

## 1. 3 つの再現単位 (共有コードなしで再実装可能)

| 単位 | 対象 | 期待値 (凍結一次ソースと機械照合) | 難度 |
|---|---|---|---|
| A | 有界領域アノマリー全数分類 (D₁ 窓) | SM 15 成分がただ一つ・SHA-256 = certificates/v62 | 整数厳密・最易 |
| B | 半空間モジュラー規格化 λ | 15 桁 + 証明付き区間 = v25.2 凍結値 (v252_bz_certificate) | 中 |
| C | Dirac 真空偏極 universality | oracle 16π²A = −1/10 (厳密) + 4 比 1 ± 1% (PRED-016 バー) | 高 |

各単位に**過去のバグからの警告** (弱二重項の二重計上・κ 床・正準化規約の差 等) を
正直に添付した — 再現者が同じ穴に落ちることは検証にならない。

## 2. 独立外部再現の 6 条件 (機械可読 — replications.yml)

different_author / independent_repository / no_shared_numerical_kernel /
protocol_frozen_before_run / commit_hash_recorded /
result_including_failures_public。**同一作者・同一 AI・同一数値カーネルの
別言語再実装は algorithmic diversity であり計数しない** (本リポジトリの実装は
AI 支援で書かれているため、同一 AI による「独立実装」は独立性の実質を欠く —
NO_SHARED_CODE.md に明記)。

## 3. 昇格・降格の手続き (reproducer/CLAIMS.md)

- **成功** → replications.yml 登録 + 該当主張 (A: QRN-GAUGE-003/008/011,
  B: QRN-GRAV-030/031/032, C: QRN-GRAV-044/048/049) の independence を昇格。
  **このとき v271_core_audit R2/R3 の期待値 (0 件) を同一コミットで更新する** —
  監査の期待値を書き換えずに数だけ増やすことはできない設計。
- **失敗** → status: failed で登録 (削除しない)・切り分け後も残る不一致は主張の
  降格 + README 到達点行への明記。

## 4. 検査 (v274_reproducer, 7 検査 PASS — 常時実行層)

[P0] 一式 9 ファイルの実在 / [P1] JSON 構文 (自前最小パーサ) /
[P2] **TOLERANCES ↔ 凍結一次ソースの一致** (単位 A の SHA-256 = certificates・
単位 B の区間 = v25.2 凍結証明書・oracle −1/10・PRED-016 バー — 外部への約束と
内部の凍結が同一の数であることの機械保証) / [P3] 単位 A の領域・期待解 ↔
certificates / [P4] 計数 0 = claims 台帳の external_replication 件数 /
[P5] CLAIMS.md の id 実在 / [P6] クリーンルーム条項の文書アンカー。

開発記録: (i) TOLERANCES の λ_⊥ 区間を 15 桁丸めで書いて P2 が検出 —
凍結値は shortest-repr の完全桁で転写する (「Json::Num は shortest-repr」の教訓の
再演)。(ii) 文書アンカーは markdown の改行折返しを跨ぐと不発 — 短い一意句を使う。

## 5. 残り (第二十八期)

- v28.0: 完全再走の儀式 (~105 CPU 時間) + Core v1 完了条件の判定 + 期統合。
- 外部作業 (期をまたぐ継続): anomaly-search / modular-BW の投稿・
  第三者再実装の公募 (reproducer/ を公募の添付物にする)。
