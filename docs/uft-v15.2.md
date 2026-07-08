# 統一場理論 v15.2 — 依存グラフの Prolog 独立推論

**Version**: v15.2
**Date**: 2026-07-08
**パイプライン**: `dag/` (結果: `results/v152_dag.txt`)
**生成物**: `dag/facts.pl` / `dag/report.md` / `dag/report.json` / `dag/qrn_dag.dot` / `dag/qrn_dag.mmd`

---

## 0. 版の型

| 項目 | 内容 |
|---|---|
| **Question** | v15.1 の依存グラフ監査 (Rust) の導出量は、独立実装で再現できるか。グラフは推論・可視化の道具として使える形になっているか |
| **Claim to promote** | 監査の導出量 (閉包・深さ・単調性) を単一実装から**二重独立実装**へ — 新 C2 主張 QRN-TOOL-002 |
| **Null models** | なし (照合そのものが対照 — 1 件でも不一致なら FAIL) |
| **Falsifier** | FAL-SUITE (照合が崩れたら TOOL-002 撤回、原因を Rust/Prolog のどちらかのバグとして記録) |
| **Artifacts** | Python 変換器 + Prolog 規則/レポート + 4 形式の生成物 + result + 本文書 + CI ジョブ |
| **Downgrade rule** | 不一致発見時は一致していた範囲を明記して C3 (実装の相互不整合の記録) へ |

## 1. 動機 — CLAUDE.md「DAG/依存グラフ」規則

v15.1 で主張の依存構造は機械可読 (claims.graph.json) になったが、導出 (推移閉包・影響範囲)
は Rust の単一実装だった。本プロジェクトの規約 (CLAUDE.md 追記) に従い、

```text
JSON → Python で読み込み・正規化 → facts.pl → Prolog で推論 → DOT/Mermaid/Markdown/JSON
```

のパイプラインを導入する。狙いは二つ:

1. **独立実装の相互検証** — v6.2 が探索に課した三重照合の流儀を、監査自体に適用する。
   Prolog の表付き推移閉包 (`:- table depends_tc/2.`) は Rust の逆辺 BFS と
   アルゴリズムが異なり、両者の一致は実装バグへの強い防御になる。
2. **推論の道具化** — 「ASM-X を抜くと何が落ちるか」が Prolog の一行クエリ
   (`falls_by_asm('ASM-GAUSS', C)`) になり、可視化 (DOT/Mermaid) が付く。

## 2. 実装

| 段 | ファイル | 内容 |
|---|---|---|
| 正規化 | `dag/json_to_facts.py` | claims.graph.json → `facts.pl` (ID をアトム化、等級を c0..c5 に、Rust の導出値を `rust_*` 事実として併記) |
| 規則 | `dag/rules.pl` | `depends_tc/2` (表付き推移閉包)、`depth_of/2`、`blast_asm/2`、`blast_fal/2`、`mono_violation/1`、`cycle_node/1`、孤児検査 |
| レポート | `dag/report.pl` | 構造検査 + **Rust との全数照合** + 4 形式出力。不一致で exit 1 |
| 駆動 | `dag/run.sh` | `python3 dag/json_to_facts.py && swipl -q dag/report.pl` |

## 3. 結果 (results/v152_dag.txt)

7 検査すべて PASS:

- 構造 (Prolog 独立導出): 非循環・等級単調・孤児なし
- 照合 (Prolog = Rust): **深さ 全主張一致 / 被支持閉包 全主張一致 / 仮定の影響範囲
  全 37 件一致 / 反証条件の射程 全 15 件一致** (v15.2 時点 97 主張・112 辺)

生成されたレポートの仮定影響範囲 (上位) は v15.1 §4 と同一の順位:
ASM-LATTICE 51 / ASM-SEED 35 / ASM-LOWDIM 30 / ASM-PDG 27 / ASM-GAUSS 27。

## 4. 運用への組み込み

- **CI**: `dag` ジョブを追加 (swi-prolog を導入して `sh dag/run.sh` を実行し、
  生成物が committed 内容と一致することを `git diff --exit-code dag/` で検査)。
- **CONTRIBUTING**: 主張追加のたびに `v151_audit --write` → `sh dag/run.sh` で
  生成物を更新する手順を明記。

## 5. 限界

- 照合するのは**導出量**であり、辺の物理的正しさは v15.1 と同じく ASM-EDGE-SEMANTICS
  (人手判定) のまま。
- Mermaid 図は 97 節点・112 辺の全体図 — 論文図としては部分グラフの切り出しが要る
  (Prolog クエリで機械的に可能; 将来の図版作業)。
