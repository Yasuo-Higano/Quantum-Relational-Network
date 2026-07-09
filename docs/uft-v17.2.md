# 統一場理論 v17.2 — 書誌の照合: 先行研究が §4 の隣に立った

**Version**: v17.2
**Date**: 2026-07-09
**成果物**: `paper/cp-complex-structure-full.md` の書誌確定 (Draft v2) + §4 への先行研究明示

---

## 0. 版の型

| 項目 | 内容 |
|---|---|
| **Question** | v17.1 の書誌プレースホルダ ([5] Kobayashi 系列・[4] PDG) を Web 照合で確定できるか。J=3.08×10⁻⁵ の出典は現行 PDG と整合するか |
| **Claim to promote** | なし (書誌版 — v14.3 の慣例) |
| **Falsifier** | 書誌の再照合で不一致が出れば正誤表対象 |
| **Artifacts** | paper 更新 + 本文書 (照合ソース一覧) |

## 1. 照合結果

| 参照 | 確定 | ソース |
|---|---|---|
| [3] Jarlskog | Phys. Rev. Lett. **55 (1985) 1039** — 表題・巻・頁を APS で確認 | link.aps.org/doi/10.1103/PhysRevLett.55.1039 |
| [4] PDG | **S. Navas et al., Phys. Rev. D 110, 030001 (2024)**。2024 CKM review の J ≈ 3.12×10⁻⁵。エンジン固定値 3.08×10⁻⁵ (事前登録時の当時値) との差 1.3% は σ=ln2 に対して無視可能 — この注記を参照文献に明記 | pdg.lbl.gov/2024/reviews/rpp2024-rev-ckm-matrix.pdf (+2016 版で 3.04, 2012 版で 2.96 の系列も確認) |
| [5] | **T. Kobayashi, K. Nishiwaki, Y. Tatsuta, "CP-violating phase on magnetized toroidal orbifolds," JHEP 04 (2017) 080, arXiv:1609.08608** — 著者・誌名・巻号を arXiv abs で確認 | arxiv.org/abs/1609.08608 |
| [2][6][7] | companion の v14.3 Web 照合を引き継ぐ (同一文献) | docs/uft-v14.3.md |

## 2. 科学的に重要な副産物 — 先行研究の明示

[5] の主張は「磁化トーラス orbifold の CP 位相には **複素構造モジュラス τ の実部の
非零が必須**」— 本計画 v16.2 の発見 (シアー = Re τ が J を作る) の**連続極限での
先行研究**である。論文 §4 に 1 文を追加した: 先行するのは方向であり、本計画の
寄与は (i) 格子厳密な実現 (全シアーで指数保護)、(ii) 仮定でなく**反証によって
強制された**発見経路、(iii) ベイズ的な値段の勘定 (+306 nats・holdout 修復) — と
分界を明示。発見の物語は変わらないが、「発見」の請求範囲は狭く正確になった。

## 3. 次の的

- LaTeX 整形・図版 (τ_re 地形・21 幾何ヒートマップ・生存曲線)。
- 数値転記の独立照合 (numpy spot check — v9.2 の教訓)。
- 投稿判断は出版作業の残高へ。
