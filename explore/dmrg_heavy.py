#!/usr/bin/env python3
"""v22.1 探索層: 重電荷ボソンの励起状態 DMRG — v21.1 の格子窓を開く。

開鎖には巻き付き ε が無い → v21.1 (環 ED) でボソンを隠したフラックス状態が
存在せず、E₁ − E₀ が中性セクターのボソン質量そのもの。
大 x (= 細かい格子) + 大 N (体積) で M·a ≪ 1 の軽域に入り、
同一 x での質量比 M(q=2)/M(q=1) を測る (比 = 残留格子効果の相殺 — v21.1 の教訓)。

出力 JSON は監査層 sim/src/bin/v221_dmrgex.rs が ED アンカーと照合して一次ソース化。
"""
import json
import sys
sys.path.insert(0, "explore")
from dmrg_schwinger import mpo, dmrg

def gap(N, x, qf, chi, sweeps, qpen=10.0):
    Ws = mpo(N, x, 0.0, {}, qf=qf, qpen=qpen)
    e0, _d0, A0 = dmrg(Ws, chi=chi, sweeps=sweeps, seed=7)
    e1, _d1, _A1 = dmrg(Ws, chi=chi, sweeps=sweeps, seed=11, ortho=[A0], lam_o=50.0)
    return e0, e1

def main():
    out = {}
    # ---- アンカー (監査層の ED と照合): N=10, x=2, q ∈ {1, 2} ----
    for qf in (1, 2):
        e0, e1 = gap(10, 2.0, float(qf), chi=64, sweeps=8)
        out[f"anchor_n10_q{qf}_e0"] = e0
        out[f"anchor_n10_q{qf}_e1"] = e1
        print(f"[anchor] N=10 q={qf}: E0={e0:.8f} E1={e1:.8f} gap={e1-e0:.6f}", file=sys.stderr)
    # ---- 本測定: N=64, x ∈ {9, 16}, q ∈ {1, 2} ----
    for x in (9.0, 16.0):
        for qf in (1, 2):
            e0, e1 = gap(64, x, float(qf), chi=128, sweeps=8)
            mg = (e1 - e0) / (2.0 * x ** 0.5)
            out[f"n64_x{int(x)}_q{qf}_gap"] = e1 - e0
            print(f"[main] N=64 x={x} q={qf}: gap={e1-e0:.5f} M/g={mg:.4f}", file=sys.stderr)
    path = sys.argv[1] if len(sys.argv) > 1 else "explore/dmrg_heavy_out.json"
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
