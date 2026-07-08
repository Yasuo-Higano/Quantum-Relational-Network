#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""claims.graph.json → dag/facts.pl (Prolog 事実) への変換・正規化。

CLAUDE.md「DAG/依存グラフ」規則のパイプライン第一段:
    JSON → Python で読み込み・正規化 → facts.pl → Prolog で推論 → レポート出力

正規化の内容:
  - ID は単一引用の Prolog アトムに (そのまま; ' はエスケープ)
  - 等級 C0..C5 は小文字アトム c0..c5 に
  - Rust 監査 (v151_audit) の導出値 (深さ・被支持閉包・影響範囲) を rust_* 事実として
    併記する — Prolog 側の独立推論と全数照合するため (独立実装の相互検証)
決定論: 入力 JSON の順序を保存し、タイムスタンプ等は書かない。
"""

import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


def atom(s: str) -> str:
    return "'" + s.replace("\\", "\\\\").replace("'", "\\'") + "'"


def main() -> None:
    src = ROOT / "claims.graph.json"
    dst = ROOT / "dag" / "facts.pl"
    data = json.loads(src.read_text(encoding="utf-8"))

    lines = []
    a = lines.append
    a("% ============================================================")
    a("% facts.pl — claims.graph.json から自動生成 (dag/json_to_facts.py)")
    a("% 手で編集しない。再生成: sh dag/run.sh")
    a("% ============================================================")
    for p in [
        "claim/3", "dep/2", "asm_of/2", "fal_of/2",
        "assumption/4", "falsifier/2",
        "rust_depth/2", "rust_closure/2", "rust_blast_asm/2", "rust_blast_fal/2",
    ]:
        a(f":- discontiguous {p}.")
    a("")

    a("% claim(Id, Version, Level).")
    for c in data["claims"]:
        a(f"claim({atom(c['id'])}, {atom(c['version'])}, {c['level'].lower()}).")
    a("")
    a("% dep(X, Y) — X は Y に依存する (Y が落ちれば X も落ちる)。")
    for c in data["claims"]:
        for d in c["deps"]:
            a(f"dep({atom(c['id'])}, {atom(d)}).")
    a("")
    a("% asm_of(Claim, Assumption) / fal_of(Claim, Falsifier)。")
    for c in data["claims"]:
        for s in c["asm"]:
            a(f"asm_of({atom(c['id'])}, {atom(s)}).")
        for f in c["fal"]:
            a(f"fal_of({atom(c['id'])}, {atom(f)}).")
    a("")
    a("% assumption(Id, Type, Scope, Status) / falsifier(Id, Status)。")
    for s in data["assumptions"]:
        a(
            f"assumption({atom(s['id'])}, {s['type']}, {s['scope']}, {s['status']})."
        )
    for f in data["falsifiers"]:
        a(f"falsifier({atom(f['id'])}, {f['status']}).")
    a("")
    a("% Rust 監査 (v151_audit) の導出値 — Prolog 独立推論との照合用。")
    for c in data["claims"]:
        a(f"rust_depth({atom(c['id'])}, {c['depth']}).")
        a(f"rust_closure({atom(c['id'])}, {c['supported_by_closure']}).")
    for s in data["assumptions"]:
        a(f"rust_blast_asm({atom(s['id'])}, {s['blast_radius']}).")
    for f in data["falsifiers"]:
        a(f"rust_blast_fal({atom(f['id'])}, {f['blast_radius']}).")
    a("")

    dst.write_text("\n".join(lines) + "\n", encoding="utf-8")
    n_dep = sum(len(c["deps"]) for c in data["claims"])
    print(
        f"[facts] claims {len(data['claims'])} / deps {n_dep} / "
        f"assumptions {len(data['assumptions'])} / falsifiers {len(data['falsifiers'])}"
        f" → {dst.relative_to(ROOT)}"
    )


if __name__ == "__main__":
    sys.exit(main())
