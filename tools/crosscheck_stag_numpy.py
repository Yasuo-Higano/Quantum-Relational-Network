#!/usr/bin/env python3
# tools/crosscheck_stag_numpy.py — 独立数値バックエンド照合 (PROMPT/6 PR9, v25.0)
#
# 自作 Jacobi (Rust, sim/src/dd.rs・stag.rs) と成熟した線形代数 (numpy/LAPACK) の
# 独立照合。3+1D staggered (x 開放・y/z 周期) N=16 の半空間について
#   - S (エンタングルメントエントロピー)      … f64 で頑健 (κe^{−κ} 抑制)
#   - λ_A(ξ=1,2) (符号付き NN 推定器)          … f64 信頼域 ξ ≤ 2 (v24.2 の較正)
# を results/v243_bwwindow.json (Rust 一次ソース) と突き合わせる。
# 重い計算ではない (単発 4096² eigh) — 検証層なので python を許可 (追加規則)。

import json
import sys

import numpy as np

N = 16
ns = N * N * N
half = N // 2


def idx(x, y, z):
    return x + N * (y + N * z)


# H 構築 (v234/stag.rs と同一規約)
H = np.zeros((ns, ns))
for x in range(N):
    for y in range(N):
        for z in range(N):
            i = idx(x, y, z)
            if x + 1 < N:
                j = idx(x + 1, y, z)
                H[i, j] += 0.5
                H[j, i] += 0.5
            ey = 0.5 if x % 2 == 0 else -0.5
            j = idx(x, (y + 1) % N, z)
            H[i, j] += ey
            H[j, i] += ey
            ez = 0.5 if (x + y) % 2 == 0 else -0.5
            j = idx(x, y, (z + 1) % N)
            H[i, j] += ez
            H[j, i] += ez

ev, vv = np.linalg.eigh(H)
nocc = ns // 2
assert ev[nocc] - ev[nocc - 1] > 1e-6, "閉殻ギャップ"

sel = np.array([idx(x, y, z) for z in range(N) for y in range(N) for x in range(half)])
V = vv[sel][:, :nocc]
C = V @ V.T

cw = np.linalg.eigvalsh(C)
cwc = np.clip(cw, 1e-14, 1 - 1e-14)
S = float(np.sum(-cwc * np.log(cwc) - (1 - cwc) * np.log(1 - cwc)))

# K = ln((1−C)/C) (f64, クランプ 1e-14) — ξ = 1, 2 の x-NN 平均 (符号付き)
w, u = np.linalg.eigh(C)
wc = np.clip(w, 1e-14, 1 - 1e-14)
kappa = np.log((1 - wc) / wc)
K = (u * kappa) @ u.T


def sel_a(x, y, z):
    return x + half * (y + N * z)


lam = {}
for xi in (1, 2):
    xb = half - 1 - xi
    vals = [K[sel_a(xb, y, z), sel_a(xb + 1, y, z)] for y in range(N) for z in range(N)]
    lam[xi] = float(np.pi * xi / np.mean(vals))

# Rust 一次ソースとの照合
with open("results/v243_bwwindow.json") as f:
    j = json.load(f)
row16 = next(r for r in j["rows"] if r["n"] == 16)
lam_rust = row16["lam_a"]
s_rust = row16["s_total"]

ok = True
print("=== numpy/LAPACK 独立照合 (N=16 半空間) ===")
d_s = abs(S - s_rust)
print(f"S:      numpy {S:.6f}  vs Rust {s_rust:.6f}  (Δ {d_s:.2e})")
ok &= d_s < 1e-4
for xi in (1, 2):
    d = abs(lam[xi] - lam_rust[xi - 1])
    print(f"λ(ξ={xi}): numpy {lam[xi]:.6f}  vs Rust {lam_rust[xi - 1]:.6f}  (Δ {d:.2e})")
    ok &= d < 1e-4  # f64 信頼域 (v24.2: ξ≤2 で relΔ < 1e-3; LAPACK↔Jacobi は同床)
print("総合判定:", "[PASS] 独立バックエンドは一次ソースを再現" if ok else "[FAIL]")
sys.exit(0 if ok else 1)
