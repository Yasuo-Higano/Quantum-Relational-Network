# 統一場理論 v17.1 — 論文本文 3: CP は複素構造を要求する (英語完全原稿)

**Version**: v17.1
**Date**: 2026-07-09
**成果物**: `paper/cp-complex-structure-full.md` (英語完全原稿 Draft v1) + `paper/geometric-yukawa-full.md` への追記 (Addendum v17.1 — companion 相互参照)

---

## 0. 版の型

| 項目 | 内容 |
|---|---|
| **Question** | 第十七期の両アーク (幾何 v15.7→v16.10・原理 v16.11→v16.13) を、companion (geometric-yukawa) と同じ規律の英語完全原稿にできるか — 「追補 or 独立短報」の判断込み |
| **Claim to promote** | なし (執筆版 — v14.1/v14.2 の慣例どおり台帳追加なし。全主張は QRN-YUK-016〜025 / QRN-SEL-002〜004 が既に保持) |
| **Null models** | — (数値は一切新造しない: 本文の全数値は results/ の一次ソースからの転記) |
| **Falsifier** | 本文数値と results/ の不一致が見つかれば正誤表対象。書誌 [2]–[7] は**未照合** (照合は v14.3 の慣例どおり別版) |
| **Artifacts** | paper 2 ファイル + 本文書 |
| **Downgrade rule** | — |

## 1. 判断: 独立短報 (companion 方式)

追補ではなく**独立論文**とした。理由: (i) アークは自立した物語を持つ (構造零 → 失敗した修理 → 発見 → 反転 → 3 つの無罪 → 診断 → 縮退 → 生存原理 → 測度) — 追補にすると companion の後半を圧倒する。(ii) companion の主結果 (対の既約性) と本稿の主結果 (複素構造の必須性・測度問題) は台帳上も別の主張群。(iii) 相互参照は companion §6 への 1 段落 (Addendum v17.1) で足りる。

## 2. 原稿の構成 (paper/cp-complex-structure-full.md)

Abstract + 12 節: §2 構造零 (J=0 と \|V_td\| miss — PRED-002/003)、§3 傾き磁束の不成立と装置の罠 4 つ (development record)、§4 シアーが J を 11 桁立てる、§5 値段 1.27 nats と +306 nats 反転・\|V_td\| hit、§6 12/12 スコアカード、§7 三つの無罪 (Wilson/格子/τ_re — τ の量子化が N とともに細かくなる指摘含む)、§8 非対称地図 (Cabibbo = 対称の病気)、§9 証拠の分解能と Higgs 幅 1/18、§10 生存原理とその解剖 (測度補正クラス — 事後解釈の旗も明記)、§11 限界、§12 再現性 (事前登録分岐・回帰ゲート・キャッシュの非一次ソース宣言・holdout 簿記)。

## 3. 規律の記録

- 執筆中の自己検査で誤り 3 件を修正: 「eleven successes」(companion の rect は 8/9 — 曖昧表現に置換)、Wilson 細分の倍率 2⁵ (誤 — 数値を外した)、「ten years of model drift」(修辞の暴走 — window の記述に置換)。
- 書誌 [5] (Kobayashi 系列 — 複素構造/モジュラー CP) は**プレースホルダ**であることを本文に明記。照合版 (v14.3 方式の Web 照合) が済むまで投稿不可。
- companion への追記は §6 末尾の 1 段落のみ (本文の歴史的数値には触れない)。

## 4. 次の的

- 書誌照合版 (Web 照合 — [2][3][4][5][6][7]、特に [5] の確定)。
- LaTeX 整形・図版 (τ_re 地形・21 幾何ヒートマップ・生存曲線 Δ(β) は figures/ 化の価値が高い)。
- 数値転記の独立照合 (numpy 側での spot check — v9.2 の教訓「独立再実装は装置バグ検出器」)。
