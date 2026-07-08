#!/bin/sh
# QRN 依存グラフの Prolog パイプライン (CLAUDE.md「DAG/依存グラフ」規則):
#   claims.graph.json → (Python 正規化) → facts.pl → (Prolog 推論) →
#   DOT / Mermaid / Markdown / JSON レポート + Rust 監査との全数照合
# 前提: v151_audit --write 済みの claims.graph.json が最新であること。
set -e
cd "$(dirname "$0")/.."
python3 dag/json_to_facts.py
swipl -q dag/report.pl
