#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v17.3 キャッシュ・モードの独立検証 — Rust の実埋め込みと numpy の複素経路の照合

sim/cache/ のモード表 (v16.11 導入・v16.11〜v16.13 の一次計算が使用) を、
**別実装・別の定式化**で検証する:

  Rust:  2N² 次元の実対称埋め込み [[R,−I],[I,R]] を循環ヤコビで対角化
  ここ:  N² 次元の複素エルミート H (H_ij = −e^{−iθ}) を numpy.linalg.eigh で対角化

照合は基底非依存の量のみで行う (縮退部分空間の基底は任意のため):
  (1) 縮退: numpy の最低 3 固有値の幅 < 1e-10、キャッシュの spread 記録と整合
  (2) ギャップ: E₄−E₃ がキャッシュの gap 記録と 1e-9 で一致
  (3) 部分空間: 射影子 ‖P_numpy − P_rust‖_F < 1e-7 (モードの張る空間が一致)
  (4) 正規直交: キャッシュ・モードのグラム行列が単位行列 (1e-10)

1 つでも不一致なら [FAIL] を出して exit(1) (make_figures.py と同じ規約)。
v9.2 の教訓「独立再実装は装置バグ検出器」の cache 版である。

再現手順: figenv/bin/python figures/verify_cache_modes.py
"""

import struct
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CACHE = ROOT / "sim" / "cache"
MAGIC = 0x51524E4D4F444531
Q = 3

nfail = 0


def check(name, ok, detail=""):
    global nfail
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        nfail += 1


def load_modes(tag, n, q, s, k):
    p = CACHE / f"modes_t{tag}_n{n}_q{q}_s{s}_k{k}.bin"
    buf = p.read_bytes()
    ns = n * n
    assert len(buf) == 8 * (6 + 2 + q * ns * 2), f"サイズ不一致: {p}"
    hdr = struct.unpack_from("<6Q", buf, 0)
    assert hdr == (MAGIC, tag, n, q, s, k), f"ヘッダ不一致: {p}"
    gap, spread = struct.unpack_from("<2d", buf, 48)
    arr = np.frombuffer(buf, dtype="<f8", offset=64).reshape(q, ns, 2)
    modes = arr[:, :, 0] + 1j * arr[:, :, 1]  # 実埋め込み (Re; Im) → 複素
    return modes, gap, spread


def hamiltonian(n, k_half, s):
    """flux_modes_shear_n の複素エルミート版 (独立の定式化)。"""
    ns = n * n
    phi = 2.0 * np.pi * Q / ns
    wl = phi * k_half / 2.0
    H = np.zeros((ns, ns), dtype=complex)

    def idx(x, y):
        return x + y * n

    def hop(i, j, th):
        H[i, j] += -np.exp(-1j * th)
        H[j, i] += -np.exp(1j * th)

    for x in range(n):
        for y in range(n):
            th_y = phi * x + wl
            if y == n - 1:
                hop(idx(x, y), idx((x + s) % n, 0), th_y)
            else:
                hop(idx(x, y), idx(x, y + 1), th_y)
            th_x = -phi * n * y if x == n - 1 else 0.0
            hop(idx(x, y), idx((x + 1) % n, y), th_x)
    return H


def main():
    n, tag = 36, 1
    samples = [(0, 0), (1, 0), (3, 7), (5, 11)]  # (シアー s, Wilson 半歩 k)
    print(f"=== キャッシュ・モードの独立検証: N={n}, {len(samples)} 標本 ===")
    for s, k in samples:
        modes, gap0, spread0 = load_modes(tag, n, Q, s, k)
        H = hamiltonian(n, k, s)
        w, v = np.linalg.eigh(H)
        spread = w[Q - 1] - w[0]
        gap = w[Q] - w[Q - 1]
        # (1) 縮退
        check(
            f"(s={s},k={k}) 縮退幅 (numpy)",
            spread < 1e-10 and abs(spread - spread0) < 1e-9,
            f"numpy {spread:.2e} / cache {spread0:.2e}",
        )
        # (2) ギャップ
        check(
            f"(s={s},k={k}) ギャップ一致",
            abs(gap - gap0) < 1e-9,
            f"numpy {gap:.6f} / cache {gap0:.6f}",
        )
        # (3) 射影子 (基底非依存)
        P_np = v[:, :Q] @ v[:, :Q].conj().T
        P_ru = modes.T @ modes.conj()  # Σ_i ψ_i ψ_i† : (ns,q)@(q,ns)
        dP = np.linalg.norm(P_np - P_ru)
        check(f"(s={s},k={k}) 射影子一致 ‖ΔP‖_F", dP < 1e-7, f"{dP:.2e}")
        # (4) 正規直交
        G = modes.conj() @ modes.T
        dG = np.linalg.norm(G - np.eye(Q))
        check(f"(s={s},k={k}) グラム = I", dG < 1e-10, f"{dG:.2e}")
    print(f"\n総合判定: {'[PASS] キャッシュは独立経路と一致' if nfail == 0 else '[FAIL]'}")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
