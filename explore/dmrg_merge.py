#!/usr/bin/env python3
"""v22.1 探索層: dmrg_one.py の部品 JSON を一次ソース JSON に併合する。"""
import json
import sys

def main():
    out = {}
    for f in sys.argv[1:-1]:
        out.update(json.load(open(f)))
    json.dump(out, open(sys.argv[-1], "w"), indent=1)
    print(f"merged {len(sys.argv)-2} files -> {sys.argv[-1]} ({len(out)} keys)")

if __name__ == "__main__":
    main()
