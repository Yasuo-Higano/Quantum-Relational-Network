#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""v9.2 論文図の生成 — 5 点 (paper/anomaly-search.md 図 2–3, paper/geometric-yukawa.md 図 1–3)

CLAUDE.md 追加規則「図表のような人間向け資料の作成には python などを利用して構わない」に基づく。
シミュレーション本体は従来どおり Rust (sim/) であり、本スクリプトは可視化のみを担う。

規約 (知的誠実性):
  - 図中の数値は results/ の機械可読 JSON (一次ソース) から読む。
  - results/ にない曲線 (波動関数の形・到達可能集合の点群) は numpy で再計算するが、
    その場合は必ず一次ソースに記録された数値 (縮退幅・ギャップ・中心・lnZ) と
    照合し、1 つでも不一致なら [FAIL] を出して exit(1) する。
    図が Rust の計算と食い違ったまま論文に載ることを防ぐための内蔵検証である。
  - 世代ラベルは v9.2 の安定規約 (中心を 0.5 サイト格子にスナップ後ソート)。本スクリプトの
    numpy 実装は生の昇順ソートだが、中心が厳密に格子上に乗る (偏差 ~1e-12) ため両者は
    同値 — その同値性こそ v9.2 で公表値の綱渡りを発見した照合である (docs/uft-v9.2.md §1)。

再現手順:
  python3 -m venv figenv && figenv/bin/pip install matplotlib numpy
  figenv/bin/python figures/make_figures.py

出力: figures/fig_controls_map.svg      (anomaly  図 2: 陰性対照の地図)
      figures/fig_u1_staircase.svg      (anomaly  図 3: U(1) の階段)
      figures/fig_zeromode_wilson.svg   (yukawa   図 1: ゼロモード局在と Wilson 線)
      figures/fig_attainable_ratios.svg (yukawa   図 2: 到達可能な質量比集合)
      figures/fig_geometry_lnz.svg      (yukawa   図 3: 幾何模型の証拠比較と緊張の解消)
図中のラベルは投稿用 (英語)。
"""

import json
import sys
from pathlib import Path

import numpy as np
import matplotlib

matplotlib.use("SVG")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
RES = ROOT / "results"
OUT = ROOT / "figures"
PNG_DIR = None  # 検収用 PNG の出力先 (コマンドラインで --png DIR)

plt.rcParams.update(
    {
        "font.size": 9,
        "axes.titlesize": 9.5,
        "axes.labelsize": 9,
        "svg.fonttype": "none",
        "figure.dpi": 110,
    }
)

FAILURES = []


def check(name, ok):
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}")
    if not ok:
        FAILURES.append(name)


def load(name):
    with open(RES / name, encoding="utf-8") as f:
        return json.load(f)


def savefig(fig, name):
    fig.savefig(OUT / name, bbox_inches="tight")
    if PNG_DIR is not None:
        fig.savefig(Path(PNG_DIR) / (name.replace(".svg", ".png")), bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"  -> figures/{name}")


# ======================================================================
# 物理の再計算 (図 1, 2 of yukawa 用) — sim/src/bin/v72_geomfn.rs と同一の模型
# ======================================================================

N = 18  # 格子
Q = 3  # 磁束 = 世代数
NS = N * N
PHI = 2.0 * np.pi * Q / NS
PHI0 = 0.83  # 局在化の一般位相 (縮退回避; v72 と同一)


def flux_modes(k_wilson):
    """磁束 Q・Wilson 線 k (Z₆, 1 サイト刻み) のゼロモード帯 (最低 Q 状態)。

    Rust 側 (v72_geomfn.rs / v91_ckmselect.rs の flux_modes) と同一の Landau ゲージ:
    y ホップに位相 φx + wl、x 境界に −φNy のツイスト。向き (転置) の規約は
    「Wilson 線 k で中心が +k サイト動く」という一次ソースの記録に合わせて固定し、
    下の照合検査がそれを検証する。戻り値: (modes (NS,Q), gap, spread)。"""
    wl = PHI * k_wilson
    H = np.zeros((NS, NS), dtype=complex)
    idx = lambda x, y: x + y * N
    for x in range(N):
        for y in range(N):
            i = idx(x, y)
            j = idx(x, (y + 1) % N)
            th = PHI * x + wl
            H[i, j] += -np.exp(+1j * th)
            H[j, i] += -np.exp(-1j * th)
            j = idx((x + 1) % N, y)
            th = -PHI * N * y if x == N - 1 else 0.0
            H[i, j] += -np.exp(+1j * th)
            H[j, i] += -np.exp(-1j * th)
    w, v = np.linalg.eigh(H)
    gap = w[Q] - w[Q - 1]
    spread = w[Q - 1] - w[0]
    return v[:, :Q], gap, spread


X_PHASE = np.exp(2j * np.pi * (np.arange(NS) % N) / N)  # 位置演算子 e^{2πix/N}


def localize(modes):
    """ゼロモード帯を位置演算子で局在化し、中心の昇順に並べる (v72 の localize と同一)。"""
    U = modes.conj().T @ (X_PHASE[:, None] * modes)  # 3×3
    M = np.exp(-1j * PHI0) * U
    H1 = (M + M.conj().T) / 2
    _, vec = np.linalg.eigh(H1)
    psi = modes @ vec  # (NS,3)
    z = np.einsum("si,s,si->i", psi.conj(), X_PHASE, psi)
    centers = (np.angle(z) / (2 * np.pi) * N) % N
    order = np.argsort(centers)
    return psi[:, order], centers[order]


def higgs_profile(sigma_h):
    """原点中心のガウス型 Higgs プロファイル (トーラス距離; v72/v91 と同一)。"""
    x = np.arange(NS) % N
    y = np.arange(NS) // N
    dx = np.minimum(x, N - x)
    dy = np.minimum(y, N - y)
    return np.exp(-(dx**2 + dy**2) / (2.0 * sigma_h**2))


def yukawa(psi_a, psi_b, phih):
    """重なり積分 Y_ij = Σ_s ψa_i(s)* ψb_j(s) φ_H(s)。"""
    return psi_a.conj().T @ (phih[:, None] * psi_b)


def ratios(m):
    """3×3 質量行列の特異値比 (σ1/σ3, σ2/σ3)。"""
    s = np.linalg.svd(m, compute_uv=False)
    return s[2] / s[0], s[1] / s[0]


# ======================================================================
# 図 A (anomaly 図 2): 陰性対照の地図 — どの条件が最小性・一意性を担うか
# ======================================================================

def fig_controls_map():
    print("[A] 対照地図 (data: results/v62_atlas.json)")
    atlas = load("v62_atlas.json")
    controls = atlas["controls"]
    labels = {
        "N1_nochiral": "$-$ chirality",
        "N2_nofactors": "$-$ charged under all factors",
        "N3_nowitten": "$-$ Witten SU(2)",
        "N4_nosu3cub": "$-$ SU(3)$^3$",
        "N5_nosu3sq": "$-$ SU(3)$^2$U(1)",
        "N6_nosu2sq": "$-$ SU(2)$^2$U(1)",
        "N7_nograv": "$-$ grav$^2$U(1)",
        "N8_nocubic": "$-$ U(1)$^3$",
    }
    check("対照 8 種が揃っている", [c["run"] for c in controls] == list(labels))
    rows = [("baseline (all conditions)", 15, 1, 468)] + [
        (labels[c["run"]], c["min_components"], c["solutions_at_15"], c["total_solutions"])
        for c in controls
    ]
    # 基線スペクトル {15:1,16:8,24:459} = 総解 468 (results/v62_atlas.txt R2 域)
    ypos = np.arange(len(rows))[::-1]
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.2, 3.1), sharey=True)

    for ax, col, base, title, xlab in (
        (ax1, 1, 15, "minimality: smallest chiral spectrum", "minimal number of components"),
        (ax2, 2, 1, "uniqueness: solutions with 15 components", "# solutions at 15 components"),
    ):
        vals = [r[col] for r in rows]
        broken = [v != base for v in vals]
        colors = ["#666666"] + ["#c0392b" if b else "#8aa8c6" for b in broken[1:]]
        ax.barh(ypos, vals, color=colors, height=0.62)
        ax.axvline(base, color="k", lw=0.8, ls="--")
        ax.set_title(title)
        ax.set_xlabel(xlab)
        for y, v in zip(ypos, vals):
            ax.text(v + 0.15 * base, y, str(v), va="center", fontsize=8)
        ax.set_xlim(0, max(vals) * 1.25)
    ax1.set_yticks(ypos)
    ax1.set_yticklabels([r[0] for r in rows])
    ax1.text(
        0.02, -0.32,
        "domain: small reps, $|Y|\\leq 3/2$, $\\leq 6$ multiplets, $\\leq 24$ components; "
        "red = removing the condition breaks the property (baseline dashed)",
        transform=ax1.transAxes, fontsize=7.5, color="#444444",
    )
    savefig(fig, "fig_controls_map.svg")
    # 一次ソースとの照合: 最小性を壊すのは {カイラル性, 全因子帯電, Witten, SU(3)³} のみ
    minb = {c["run"] for c in controls if c["min_components"] < 15}
    check(
        "最小性を壊す対照 = {N1,N2,N3,N4} (線形アノマリーは冗長)",
        minb == {"N1_nochiral", "N2_nofactors", "N3_nowitten", "N4_nosu3cub"},
    )
    uniqb = {c["run"] for c in controls if c["solutions_at_15"] != 1}
    check("一意性を壊す対照 ⊇ {U(1)³} (9 解)", "N8_nocubic" in uniqb)


# ======================================================================
# 図 B (anomaly 図 3): U(1) の階段 — Y → B−L → 存在しない
# ======================================================================

def fig_u1_staircase():
    print("[B] U(1) の階段 (data: results/v62_atlas.json, v71_twou1.json, v82_threeu1.json)")
    two = load("v71_twou1.json")
    three = load("v82_threeu1.json")
    min16 = two["min_rank2_components"][2]
    n3a, n3b, ctrl = (
        three["rank3_solutions_u3a"],
        three["rank3_solutions_u3b"],
        three["control_solutions"],
    )
    check("rank-2 の最小解は 16 成分 (SM+ν_R)", min16 == 16)
    check("rank-3 カイラル解は両窓で 0", n3a == 0 and n3b == 0)
    check("対照 (カイラル性オフ) は 355 解", ctrl == 355)

    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    # 階段: 段 = U(1) の本数。3 段目は「存在しない」ので破線で描く。
    for x0, y in ((0, 1), (1, 2)):
        ax.plot([x0, x0 + 1], [y, y], color="#2c3e50", lw=2.5)
        if y > 1:
            ax.plot([x0, x0], [y - 1, y], color="#2c3e50", lw=2.5)
    ax.plot([2, 2], [2, 3], color="#c0392b", lw=1.6, ls="--")
    ax.plot([2, 3], [3, 3], color="#c0392b", lw=1.6, ls="--")
    ax.text(2.5, 3.06, "✗", ha="center", color="#c0392b", fontsize=13)

    ax.text(0.5, 1.10, "U(1)$_Y$", ha="center", fontsize=11, weight="bold")
    ax.text(
        0.5, 0.82,
        "unique minimal chiral spectrum:\nSM generation (15 components)\n"
        "second charge on SM alone: all $\\propto Y$",
        ha="center", va="top", fontsize=7.6,
    )
    ax.text(1.5, 2.10, "U(1)$_{B-L}$", ha="center", fontsize=11, weight="bold")
    ax.text(
        1.5, 1.82,
        f"requires $\\nu_R$ ({min16} components);\nunique charge plane "
        "$\\{Y, B\\!-\\!L\\}$\n(Plücker classification of rank-2 planes)",
        ha="center", va="top", fontsize=7.6,
    )
    ax.text(2.5, 2.62, "no third U(1)", ha="center", fontsize=10, weight="bold", color="#c0392b")
    ax.text(
        2.5, 2.44,
        f"rank-3 chiral solutions in window: {n3a}\n"
        f"(control without chirality: {ctrl} —\nthe instrument is not blind)\n"
        "E$_6$: chiral core of $\\mathbf{27}$ = SM generation",
        ha="center", va="top", fontsize=7.6, color="#7f3030",
    )
    ax.set_xlim(-0.15, 3.15)
    ax.set_ylim(0.3, 3.45)
    ax.set_xticks([0.5, 1.5, 2.5])
    ax.set_xticklabels(["1st U(1)", "2nd U(1)", "3rd U(1)"])
    ax.set_yticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.set_title("the U(1) staircase ends after two steps: anomaly consistency counts the forces")
    savefig(fig, "fig_u1_staircase.svg")


# ======================================================================
# 図 C (yukawa 図 1): ゼロモードの局在と Wilson 線による平行移動
# ======================================================================

def fig_zeromode_wilson():
    print("[C] ゼロモード局在と Wilson 線 (再計算 + results/v72_geomfn.txt の数値と照合)")
    data = {}
    gaps = {}
    for k in range(6):
        modes, gap, spread = flux_modes(k)
        psi, centers = localize(modes)
        data[k] = (psi, centers)
        gaps[k] = gap
        if k == 0:
            # 一次ソース (results/v72_geomfn.txt [1]): 縮退幅 1.8e-13 ≪ ギャップ 0.115,
            # 中心 [6.0, 12.0, 18.0] (= 等間隔 6 サイト), 幅 ≈ 2.93
            check("k=0: 縮退幅 < 1e-10 (一次ソース 1.8e-13)", spread < 1e-10)
            check("k=0: ギャップ = 0.115 ±0.005 (一次ソース 0.115)", abs(gap - 0.115) < 0.005)
            check(
                "k=0: 中心が 6 サイト等間隔 (一次ソース [6,12,18])",
                np.allclose(np.diff(np.sort(centers)), 6.0, atol=0.05),
            )
            p0 = np.abs(psi[:, 0]) ** 2
            x = np.arange(NS) % N
            c0 = data[0][1][0]
            dx = (x - c0 + N / 2) % N - N / 2
            width = np.sqrt((p0 * dx**2).sum())
            check("k=0: モード幅 ≈ 2.93 ±0.1 (一次ソース 2.93)", abs(width - 2.93) < 0.1)
    # Wilson 線 1 目盛 = 1 サイトの平行移動 (縮退・ギャップ不変):
    # 中心集合 {c} が k=0 の集合を −k した {c₀−k mod N} と一致すること
    # (方向はゲージの向きの規約; 一次ソース results/v72_geomfn.txt は「1 サイト刻み」)
    shift_ok = all(
        np.allclose(np.sort(data[k][1]), np.sort((data[0][1] - k) % N), atol=0.02)
        for k in range(1, 6)
    )
    check("Wilson 線 k で中心が k サイト平行移動 (一次ソース: 1 サイト刻み)", shift_ok)
    check("ギャップは k に依らず不変 (±1e-6)", max(gaps.values()) - min(gaps.values()) < 1e-6)

    x = np.arange(N)
    fig, (ax1, ax2) = plt.subplots(
        1, 2, figsize=(7.2, 2.9), gridspec_kw={"width_ratios": [1.5, 1]}
    )
    colors = ["#1f6f8b", "#c98a00", "#7d4a94"]
    for g in range(3):
        prof0 = (np.abs(data[0][0][:, g]) ** 2).reshape(N, N).sum(axis=0)
        prof3 = (np.abs(data[3][0][:, g]) ** 2).reshape(N, N).sum(axis=0)
        ax1.plot(x, prof0, color=colors[g], lw=1.6, label=f"generation {g+1}")
        ax1.plot(x, prof3, color=colors[g], lw=1.3, ls="--")
    ax1.plot([], [], color="k", lw=1.6, label="$k=0$")
    ax1.plot([], [], color="k", lw=1.3, ls="--", label="$k=3$ (shifted 3 sites)")
    ax1.set_xlabel("lattice coordinate $x$")
    ax1.set_ylabel("$\\sum_y |\\psi_g(x,y)|^2$")
    ax1.set_title(f"zero modes of the flux-$Q{{=}}3$ torus (${N}\\times{N}$)")
    ax1.legend(fontsize=7.2, loc="upper right")
    ax1.set_xlim(0, N - 1)

    # 各世代の移動距離 |Δc|(k): 円周距離で測ると厳密に k サイト
    for g in range(3):
        c0 = data[0][1][g]
        dist = []
        for k in range(6):
            # k=0 の世代 g から出たモード (期待位置 c0−k) に最も近い中心を追跡
            exp_pos = (c0 - k) % N
            d = np.abs((data[k][1] - exp_pos + N / 2) % N - N / 2)
            j = int(np.argmin(d))
            dd = (c0 - data[k][1][j]) % N
            dist.append(min(dd, N - dd))
        ax2.plot(range(6), dist, "o", color=colors[g], ms=4.5, alpha=0.85)
    ax2.plot([0, 5], [0, 5], "-", color="#888888", lw=1.0, zorder=0, label="$|\\Delta c| = k$ (exact)")
    ax2.set_xlabel("Wilson line $k \\in \\mathbb{Z}_6$")
    ax2.set_ylabel("center displacement (sites)")
    ax2.set_title("rigid translation by $k$ sites;\ngap and degeneracy invariant")
    ax2.legend(fontsize=7.2, loc="upper left")
    savefig(fig, "fig_zeromode_wilson.svg")
    return data


# ======================================================================
# 図 D (yukawa 図 2): 到達可能な質量比集合 — 単一 T² の床と T²×T² の 2 乗
# ======================================================================

def t2_mass_lnz_stable(mode_data):
    """安定ラベルでの T²×T² 質量のみ証拠 (v9.2 の eval_mass と同じ模型) を再計算する。
    図の点群が v9.2 の模型空間そのものであることの端到端検査に使う。"""
    sig = np.log(2.0)
    norm2 = -np.log(2 * np.pi * sig * sig)
    tgt = {
        "u": (np.log(1.3e-5), np.log(3.7e-3)),
        "d": (np.log(1.1e-3), np.log(2.2e-2)),
        "e": (np.log(2.9e-4), np.log(5.9e-2)),
    }

    def lse(v):
        v = np.asarray(v).ravel()
        m = np.max(v)
        return m + np.log(np.sum(np.exp(v - m)))

    terms = []
    for sh in (1.0, 1.5, 2.0, 2.5):
        ph = higgs_profile(sh)
        yt = np.zeros((6, 6, 3, 3), complex)
        for a in range(6):
            for b in range(6):
                yt[a, b] = yukawa(mode_data[a][0], mode_data[b][0], ph)
        A = np.arange(36)
        M = yt[A[:, None] % 6, A[None, :] % 6] * yt[A[:, None] // 6, A[None, :] // 6]
        s = np.linalg.svd(M, compute_uv=False)
        r1 = np.log(s[..., 2] / s[..., 0])
        r2 = np.log(s[..., 1] / s[..., 0])
        ll = {
            q: -((r1 - tgt[q][0]) ** 2 + (r2 - tgt[q][1]) ** 2) / (2 * sig * sig) + norm2
            for q in "ude"
        }
        per_q = [lse(ll["u"][kq]) + lse(ll["d"][kq]) for kq in range(36)]
        terms.append(lse(per_q) + lse(ll["e"]))
    return lse(terms) - (10 * np.log(6) + np.log(4))


def fig_attainable_ratios(mode_data):
    print("[D] 到達可能な質量比集合 (再計算 + v7.2/v9.2 の一次ソースと照合)")
    phih = higgs_profile(1.0)  # v72 の MAP: σ_H = 1
    ytab = {}
    for ka in range(6):
        for kb in range(6):
            ytab[(ka, kb)] = yukawa(mode_data[ka][0], mode_data[kb][0], phih)
    singles = np.array([ratios(ytab[p]) for p in ytab])
    prods = np.array(
        [
            ratios(ytab[(a1, b1)] * ytab[(a2, b2)])
            for a1 in range(6)
            for b1 in range(6)
            for a2 in range(6)
            for b2 in range(6)
        ]
    )
    # 一次ソースとの照合:
    #  [i] 単一 T² MAP (k_Q=k_u=3; results/v72_geomfn.txt [3] — ラベル置換に不変):
    #      m_u/m_t = 2.98e-3, m_c/m_t = 5.48e-2
    r1s, r2s = ratios(ytab[(3, 3)])
    check("単一 T² MAP (3,3): r1 = 2.98e-3 ±1% (v7.2, ラベル不変)", abs(r1s / 2.98e-3 - 1) < 0.01)
    check("単一 T² MAP (3,3): r2 = 5.48e-2 ±1% (v7.2, ラベル不変)", abs(r2s / 5.48e-2 - 1) < 0.01)
    #  [ii] 単一トーラスの床 (到達下限) ~3e-3
    floor = singles[:, 0].min()
    check("単一 T² の床: min r1 ~ 3e-3 (到達下限, 一次ソース 2.98e-3)", 2.5e-3 < floor < 3.5e-3)
    #  [iii] 端到端: この点群から作った T²×T² 質量のみ証拠が v9.2 の安定ラベル値と一致
    lnz2 = t2_mass_lnz_stable(mode_data)
    ref = load("v92_labelstab.json")["lnZ_mass"]["T2xZ6"]["stable_label"]
    check(
        f"T²×T² 質量のみ lnZ = {ref:.2f} ±0.02 (v9.2 安定ラベルの一次ソース)",
        abs(lnz2 - ref) < 0.02,
    )
    # 積模型 MAP 点 (安定ラベル側の代表として散布図の注釈に使う)
    r1p, r2p = min(((ratios(ytab[(a1, b1)] * ytab[(a2, b2)]))
                    for a1 in range(6) for b1 in range(6)
                    for a2 in range(6) for b2 in range(6)),
                   key=lambda r: abs(np.log(r[0] / 1.3e-5)) + abs(np.log(r[1] / 3.7e-3)))

    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    ax.scatter(
        prods[:, 0], prods[:, 1], s=7, color="#8aa8c6", alpha=0.45, lw=0,
        label="$T^2{\\times}T^2$ (1296 Wilson configs)",
    )
    ax.scatter(
        singles[:, 0], singles[:, 1], s=26, color="#c0392b", alpha=0.9, lw=0,
        label="single $T^2$ (36 Wilson configs)", zorder=3,
    )
    obs = {  # 実測 (v6.5 以来共通のデータ; results/v72_geomfn.txt の実測列)
        "up ($m_u/m_t,\\ m_c/m_t$)": (1.3e-5, 3.7e-3, "*", 130),
        "down ($m_d/m_b,\\ m_s/m_b$)": (1.1e-3, 2.2e-2, "P", 60),
        "lepton ($m_e/m_\\tau,\\ m_\\mu/m_\\tau$)": (2.9e-4, 5.9e-2, "X", 55),
    }
    for lab, (r1, r2, mk, sz) in obs.items():
        ax.scatter([r1], [r2], marker=mk, s=sz, color="#1a7a3a", zorder=4, label="observed " + lab)
    ax.axvline(floor, color="#c0392b", lw=0.9, ls=":")
    ax.text(
        floor * 0.85, 1.05e-2,
        "single-torus floor $\\approx 3{\\times}10^{-3}$:\ndepth set by flux alone",
        fontsize=7.2, color="#c0392b", ha="right",
    )
    ax.annotate(
        "squared suppression\nreaches the up hierarchy",
        xy=(r1p, r2p), xytext=(0.03, 0.42), textcoords="axes fraction",
        fontsize=7.4, color="#33506b",
        arrowprops=dict(arrowstyle="->", color="#33506b", lw=0.9),
    )
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlim(6e-6, 3e-2)
    ax.set_ylim(2.2e-3, 1.6e-1)
    ax.set_xlabel("lightest / heaviest singular value  $\\sigma_1/\\sigma_3$")
    ax.set_ylabel("middle / heaviest  $\\sigma_2/\\sigma_3$")
    ax.set_title("attainable mass-ratio sets at $\\sigma_H{=}1$ (no random coefficients)")
    ax.legend(fontsize=6.8, loc="lower right", framealpha=0.9)
    savefig(fig, "fig_attainable_ratios.svg")


# ======================================================================
# 図 E (yukawa 図 3): 幾何模型の証拠比較 — 質量のみ vs 全 9 量、緊張の解消
# ======================================================================

def fig_geometry_lnz():
    print(
        "[E] 幾何模型の証拠比較 (data: results/v92_labelstab.json, v81_geoselect.json, "
        "v91_ckmselect.json, v72_geomfn.json)"
    )
    g92 = load("v92_labelstab.json")
    g81 = load("v81_geoselect.json")
    g91 = load("v91_ckmselect.json")
    g72 = load("v72_geomfn.json")
    order = ["T1xZ6", "T1xZ12", "T2xZ6", "T2xZ12", "T3xZ6", "T3xZ12"]
    disp = {
        "T1xZ6": "$T^1{\\times}Z_6$", "T1xZ12": "$T^1{\\times}Z_{12}$",
        "T2xZ6": "$T^2{\\times}Z_6$", "T2xZ12": "$T^2{\\times}Z_{12}$",
        "T3xZ6": "$T^3{\\times}Z_6$", "T3xZ12": "$T^3{\\times}Z_{12}$",
    }
    mass = g92["lnZ_mass"]
    nine = g92["lnZ_nine"]
    nine_order = [g for g in order if g in nine]
    check("v9.2: 勝者と順位は両ラベル規約で不変", g92["winner_invariant"] and g92["ranking_invariant"])
    for key, tbl in (("mass", mass), ("nine", nine)):
        win_stb = max(tbl, key=lambda g: tbl[g]["stable_label"])
        win_pub = max(tbl, key=lambda g: tbl[g]["published_label"])
        check(f"{key}: 勝者 = T³×Z₆ (安定側/公表側とも)", win_stb == win_pub == "T3xZ6")

    fig, (ax1, ax2, ax3) = plt.subplots(
        1, 3, figsize=(10.2, 3.2), gridspec_kw={"width_ratios": [1.5, 1.1, 0.9], "wspace": 0.34}
    )
    # (a) 質量のみ (6 幾何): 安定ラベルの棒 + 公表側の菱形マーカー
    xs = np.arange(len(order))
    stb = [mass[g]["stable_label"] for g in order]
    pub = [mass[g]["published_label"] for g in order]
    ax1.bar(xs, stb, 0.55, color="#8aa8c6", label="stable labelling (v9.2)")
    ax1.plot(xs, pub, "D", ms=4, color="#444444", mfc="none",
             label="published convention (v8.1)")
    ax1.axhline(g72["lnZ_m1_from_v65"], color="#888888", lw=0.9, ls=":",
                label="M1: free FN charges + 18 rand. coeffs")
    ax1.axhline(g72["lnZ_m0_upper_bound_from_v65"], color="#aaaaaa", lw=0.9, ls="--",
                label="M0 anarchy (upper bound)")
    star = order.index("T3xZ6")
    ax1.text(star, stb[star] + 1.2, "★", ha="center", fontsize=11, color="#b8860b")
    ax1.set_xticks(xs)
    ax1.set_xticklabels([disp[g] for g in order], fontsize=7.2)
    ax1.set_ylabel("ln evidence (6 mass ratios)")
    ax1.set_title("mass-only evidence:\n$T^3$ wins under both label conventions")
    ax1.legend(fontsize=6.6, loc="lower right")
    ax1.set_ylim(min(stb + pub) * 1.22, 0)

    # (b) 全 9 量 (4 幾何; T³×Z₁₂ は三重和 5×10⁹ のため対象外)
    xs9 = np.arange(len(nine_order))
    stb9 = [nine[g]["stable_label"] for g in nine_order]
    pub9 = [nine[g]["published_label"] for g in nine_order]
    ax2.bar(xs9, stb9, 0.5, color="#2c5f8a", label="stable labelling (v9.2)")
    ax2.plot(xs9, pub9, "D", ms=4, color="#444444", mfc="none",
             label="published convention (v9.1)")
    star9 = nine_order.index("T3xZ6")
    ax2.text(star9, stb9[star9] + 1.6, "★", ha="center", fontsize=11, color="#b8860b")
    ax2.set_xticks(xs9)
    ax2.set_xticklabels([disp[g] for g in nine_order], fontsize=7.2)
    ax2.set_ylabel("ln evidence (6 masses + 3 CKM)")
    ax2.set_title("full 9-observable evidence:\n$T^3$ keeps masses and CKM")
    ax2.legend(fontsize=6.6, loc="lower right")
    ax2.set_ylim(min(stb9 + pub9) * 1.25, 0)

    # (c) 緊張の解消 — 9 量の点評価: 「質量のみ MAP」→「全 9 量 MAP」(公表側規約の記録値)
    m81 = {m["label"]: m for m in g81["models"]}
    m91 = {m["label"]: m for m in g91["models"]}
    pairs = [("T² × Z₆", "$T^2{\\times}Z_6$", "#8a6d3b"), ("T³ × Z₆", "$T^3{\\times}Z_6$", "#2c5f8a")]
    for lab, dl, col in pairs:
        a = m81[lab]["map_lnL9_point"]
        b = m91[lab]["map_lnL"]
        ax3.plot([0, 1], [a, b], "o-", color=col, lw=1.6, ms=4.5, label=dl)
        ax3.text(1.07, b, f"{b:.1f}", fontsize=7, va="center", color=col)
        ax3.text(-0.07, a, f"{a:.1f}", fontsize=7, va="center", ha="right", color=col)
    check(
        "緊張の入れ替わり: 質量 MAP 点評価は T²>T³、全 9 量 MAP は T³>T²",
        m81["T² × Z₆"]["map_lnL9_point"] > m81["T³ × Z₆"]["map_lnL9_point"]
        and m91["T³ × Z₆"]["map_lnL"] > m91["T² × Z₆"]["map_lnL"],
    )
    ax3.set_xlim(-0.6, 1.6)
    ax3.set_xticks([0, 1])
    ax3.set_xticklabels(["mass-only MAP\n(v8.1 point eval)", "joint MAP\n(v9.1)"], fontsize=7.2)
    ax3.set_ylabel("$\\ln L$ on all 9 observables")
    ax3.set_title("the mass–mixing 'tension'\nwas a point-estimate artifact")
    ax3.legend(fontsize=7.0, loc="center right")
    savefig(fig, "fig_geometry_lnz.svg")


# ======================================================================

def main():
    global PNG_DIR
    if len(sys.argv) >= 3 and sys.argv[1] == "--png":
        PNG_DIR = sys.argv[2]
    print("=== v9.2 論文図の生成 (数値は results/ と照合) ===")
    fig_controls_map()
    fig_u1_staircase()
    mode_data = fig_zeromode_wilson()
    fig_attainable_ratios(mode_data)
    fig_geometry_lnz()
    if FAILURES:
        print(f"総合判定: [FAIL] ({len(FAILURES)} 件の照合が一次ソースと不一致)")
        sys.exit(1)
    print("総合判定: [PASS] (全図の数値が一次ソースと一致)")


if __name__ == "__main__":
    main()
