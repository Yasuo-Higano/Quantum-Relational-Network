#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v17.3 論文 3 (cp-complex-structure) の図 3 点

規約は make_figures.py と同一 (知的誠実性):
  - 図中の数値は results/ の機械可読 JSON (一次ソース) のみから読む。
  - 内蔵照合: 文書・論文に公表済みのアンカー値と突き合わせ、不一致なら exit(1)。

再現手順: figenv/bin/python figures/make_figures_cp.py
出力: figures/fig_cp_asym_map.svg    (§8: 21 幾何の lnZ₁₀ / |V_us| ヒートマップ)
      figures/fig_cp_taure.svg       (§7: τ_re 走査 — 谷 1/12 と峰 1/9)
      figures/fig_cp_depth_curve.svg (§10: Depth 生存曲線 Δ(β) と崩壊点)
図中ラベルは投稿用 (英語)。
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("SVG")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
nfail = 0


def check(name, ok, detail=""):
    global nfail
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}  {detail}")
    if not ok:
        nfail += 1


def jload(rel):
    return json.loads((ROOT / "results" / rel).read_text())


# ---------- 図 1: 非対称シアーの地図 (v16.9) ----------
def fig_asym_map():
    d = jload("v169_asymshear.json")
    rows = d["rows"]
    # アンカー照合 (docs/uft-v16.9.md §1 の表)
    by = {(r["s1"], r["s2"]): r for r in rows}
    check("v169 (3,3) lnZ₁₀ = −21.757", abs(by[(3, 3)]["lnz10"] - (-21.756581)) < 1e-4)
    check("v169 (0,0) rect 床 = −269.10", abs(by[(0, 0)]["lnz10"] - (-269.103744)) < 1e-4)
    check("v169 (2,5) |V_us| factor = 1.03", abs(by[(2, 5)]["vus_factor"] - 1.03) < 0.005)

    Z = np.full((6, 6), np.nan)
    V = np.full((6, 6), np.nan)
    for r in rows:
        Z[r["s1"], r["s2"]] = r["lnz10"]
        V[r["s1"], r["s2"]] = r["vus_factor"]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.4))
    # lnZ₁₀ (rect 床は色域を潰すのでクリップし注記)
    Zc = np.clip(Z, -30, None)
    im0 = axes[0].imshow(Zc, origin="lower", cmap="viridis", vmin=-30, vmax=-21)
    axes[0].set_title(r"ln $Z_{10}$ (10-observable evidence)")
    for (i, j), val in np.ndenumerate(Z):
        if not np.isnan(val):
            txt = f"{val:.1f}" if val > -100 else "floor"
            axes[0].text(j, i, txt, ha="center", va="center", fontsize=7,
                         color="white" if Zc[i, j] < -25 else "black")
    im1 = axes[1].imshow(V, origin="lower", cmap="RdYlGn_r", vmin=1.0, vmax=3.4)
    axes[1].set_title(r"$|V_{us}|$ MAP factor (1 = measured)")
    for (i, j), val in np.ndenumerate(V):
        if not np.isnan(val):
            axes[1].text(j, i, f"{val:.2f}", ha="center", va="center", fontsize=7)
    for ax, im in ((axes[0], im0), (axes[1], im1)):
        ax.set_xlabel(r"$s_2$")
        ax.set_ylabel(r"$s_1$")
        ax.set_xticks(range(6))
        ax.set_yticks(range(6))
        fig.colorbar(im, ax=ax, shrink=0.85)
    axes[0].annotate("evidence best (3,3)", xy=(3, 3), xytext=(0.2, 4.7),
                     arrowprops=dict(arrowstyle="->"), fontsize=8)
    axes[1].annotate("Cabibbo best (2,5)", xy=(5, 2), xytext=(0.2, 4.7),
                     arrowprops=dict(arrowstyle="->"), fontsize=8)
    fig.suptitle("Asymmetric shear map at N = 36: evidence prefers the diagonal, "
                 "Cabibbo prefers asymmetry", fontsize=11)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "fig_cp_asym_map.svg")
    plt.close(fig)


# ---------- 図 2: τ_re 走査 (v16.8) ----------
def fig_taure():
    d = jload("v168_taure.json")
    rows = d["rows"]
    check("v168 s=3 lnZ₁₀ = −21.757", abs(rows[2]["lnz10"] - (-21.757)) < 5e-4)
    check("v168 s=2 回帰 = −22.263", abs(rows[1]["lnz10"] - (-22.263)) < 5e-4)
    tau = [r["tau_re"] for r in rows]
    z = [r["lnz10"] for r in rows]
    vus = [r["vus_factor"] for r in rows]
    fig, ax1 = plt.subplots(figsize=(6.4, 4.2))
    ax1.plot(tau, z, "o-", color="tab:blue", label=r"ln $Z_{10}$")
    ax1.set_xlabel(r"$\tau_{\rm re} = s/N$  (N = 36)")
    ax1.set_ylabel(r"ln $Z_{10}$", color="tab:blue")
    ax1.tick_params(axis="y", labelcolor="tab:blue")
    ax2 = ax1.twinx()
    ax2.plot(tau, vus, "s--", color="tab:red", label=r"$|V_{us}|$ factor")
    ax2.axhline(1.8, color="tab:red", lw=0.8, ls=":", alpha=0.6)
    ax2.set_ylabel(r"$|V_{us}|$ MAP factor", color="tab:red")
    ax2.tick_params(axis="y", labelcolor="tab:red")
    ax1.annotate("valley 1/12\n(unreachable at N=18)", xy=(1 / 12, z[2]),
                 xytext=(0.095, -23.3), arrowprops=dict(arrowstyle="->"), fontsize=8)
    ax1.annotate("ridge 1/9", xy=(1 / 9, z[3]), xytext=(0.115, -24.6),
                 arrowprops=dict(arrowstyle="->"), fontsize=8)
    ax1.set_title(r"$\tau_{\rm re}$ scan: no Cabibbo valley on the symmetric section")
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "fig_cp_taure.svg")
    plt.close(fig)


# ---------- 図 3: Depth 生存曲線 (v16.12) ----------
def fig_depth_curve():
    d = jload("v1612_depthscan.json")
    uni = d["uniform"]
    betas = d["betas"]
    dz = [z - uni for z in d["lnz_beta"]]
    d_inf = d["lnz_inf"] - uni
    check("v1612 uniform = −23.1214", abs(uni - (-23.1214)) < 5e-4)
    check("v1612 峰 β=1 Δ=+3.338", abs(dz[2] - 3.338) < 5e-3)
    check("v1612 hard Δ = −8386", abs(d_inf - (-8385.912)) < 0.5)
    fig, ax = plt.subplots(figsize=(6.4, 4.2))
    ax.plot(betas, dz, "o-", color="tab:purple")
    ax.axhline(1.0, color="gray", lw=0.9, ls="--")
    ax.text(6.6, 1.25, "survival line +1", fontsize=8, color="gray")
    ax.axhline(0.0, color="gray", lw=0.6)
    ax.annotate(f"peak β=1.0 (Δ=+{dz[2]:.2f})", xy=(1.0, dz[2]), xytext=(2.2, 3.1),
                arrowprops=dict(arrowstyle="->"), fontsize=8)
    ax.annotate("collapse β*=2.0", xy=(2.0, dz[4]), xytext=(3.1, 0.4),
                arrowprops=dict(arrowstyle="->"), fontsize=8)
    ax.annotate(f"hard selection (β→∞): Δ = {d_inf:.0f}", xy=(8.0, dz[-1]),
                xytext=(3.6, -8.5), fontsize=8,
                arrowprops=dict(arrowstyle="->", ls=":"))
    ax.set_xlabel(r"depth tilt β  (prior ∝ $e^{\beta\,\cdot\,{\rm depth}}$)")
    ax.set_ylabel(r"Δ ln Z vs uniform prior (nats)")
    ax.set_title("Survival curve of the depth-tilted prior: "
                 "what survives is a tilt of order one")
    ax.set_ylim(-14, 4.6)
    fig.tight_layout()
    fig.savefig(ROOT / "figures" / "fig_cp_depth_curve.svg")
    plt.close(fig)


def main():
    print("=== 論文 3 の図生成 (一次ソース照合内蔵) ===")
    fig_asym_map()
    fig_taure()
    fig_depth_curve()
    print(f"\n総合判定: {'[PASS] 図 3 点を生成 (アンカー照合済み)' if nfail == 0 else '[FAIL]'}")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
