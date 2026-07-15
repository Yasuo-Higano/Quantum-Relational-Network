#!/usr/bin/env python3
"""v22.1 探索層の単一測定ランナー — dmrg_heavy.gap() を 1 点だけ計算して JSON に書く。

直列 4 測定 (~32h) を 4 並列 (~8h) に割るための道具。シード・χ・sweep は
dmrg_heavy.py と完全同一 (kind=anchor: N=10 χ=64 / kind=main: N=64 χ=128)。
各 JSON は explore/dmrg_merge.py が explore/dmrg_heavy_out.json に併合する。
"""
import json
import sys
sys.path.insert(0, "explore")
from dmrg_heavy import gap

def main():
    kind, x, qf, out = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4]
    if kind == "anchor":
        e0, e1 = gap(10, x, qf, chi=64, sweeps=8)
        d = {f"anchor_n10_q{int(qf)}_e0": e0, f"anchor_n10_q{int(qf)}_e1": e1}
    else:
        e0, e1 = gap(64, x, qf, chi=128, sweeps=8)
        d = {f"n64_x{int(x)}_q{int(qf)}_gap": e1 - e0,
             f"n64_x{int(x)}_q{int(qf)}_e0": e0, f"n64_x{int(x)}_q{int(qf)}_e1": e1}
    json.dump(d, open(out, "w"), indent=1)
    print(json.dumps(d), file=sys.stderr)

if __name__ == "__main__":
    main()
