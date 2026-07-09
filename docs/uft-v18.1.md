# 統一場理論 v18.1 — 出版整形: 4 論文の LaTeX 化とコンパイル検証

**Version**: v18.1
**Date**: 2026-07-10
**成果物**: `paper/tex/{anomaly-search, geometric-yukawa, cp-complex-structure, measure-dissolution}.tex` (+ 検証済み PDF 4 本) — 出版作業の残高「LaTeX 整形」の返済

---

## 0. 版の型

| 項目 | 内容 |
|---|---|
| **Question** | 4 論文 (.md 完全原稿) を投稿可能な LaTeX (revtex4-2, PRD preprint) に落とし、コンパイルまで検証できるか |
| **Claim to promote** | なし (出版整形版) |
| **Null models / 検証** | tectonic (自己完結 LaTeX エンジン, brew 導入) で 4 本ともコンパイル成功・紙面を PNG 検収 (abstract・数式・表・引用の描画確認) |
| **Falsifier** | .tex と .md の数値不一致は正誤表対象 (数値は全て .md 経由で results/ から) |
| **Artifacts** | tex 4 本 + 本文書 (.pdf は生成物 — リポジトリには .tex のみ) |

## 1. 内容

- **revtex4-2 (aps, prd, preprint)** で統一。著者・所属は**仮置き** (TODO コメント明示 —
  投稿前に要確認)。図は cp 論文の Figures 節がパスで参照 (投稿時に SVG → PDF 化)。
- 変換はコピーエディットを兼ねた。**発見・修正した原稿の欠陥**:
  - anomaly-search-full.md: 英語原稿への日本語混入 2 箇所 (§6「全称」→ universal /
    参考文献 [7] の括弧内)。
  - cp-complex-structure: §11 が v17.5 以降の展開より古い — **Addendum (v17.12)** を
    .md/.tex 両方に追加 (τ 谷・測度判定・orientation の解決を companion 参照で明示)。
  - measure-dissolution: Limitations の orientation 項に v17.13 の結果
    (γ = +66.8° が測定誤差内) を後日譚として追記。
- ツール: tectonic 0.16.9 + poppler (紙面検収用) を導入 — CLAUDE.md 追加規則
  「人間向け資料にフリーソフトウェア可」の範囲。

## 2. 開発記録 (小さな発見)

v18.0 スイート再実行で `results/v153_corev2.json` / `v154_continuum.json` の
末桁が揺れた (~1e-13 級)。原因はスレッド集約順序の非決定性 (Lanczos/並列和) —
[PASS] 判定・物理値は不変。**一次ソースは初出値を保持する規約**に従い revert し、
ここに記録する (v153/v154 の JSON 末桁はスレッド順序に依存する)。

## 3. 残作業 (出版)

- 著者・所属・謝辞の確定 (ユーザー判断)。
- 投稿先の決定 (草稿は SciPost/CPC/PRD [anomaly]・JHEP/PRD [他 3 本] を想定)。
- cp 論文の図 3 点の PDF 化と \includegraphics 組み込み・arXiv 用パッケージング。
