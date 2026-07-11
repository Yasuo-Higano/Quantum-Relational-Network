#!/usr/bin/env python3
"""v21.3 探索層: Schwinger 模型の 2 サイト DMRG (numpy のみ)。

二層構造 (PROMPT/3 §7) の初適用。探索層の出力 JSON は監査層
(sim/src/bin/v213_dmrgaudit.rs) が ED アンカーと照合して一次ソース化する。
自己検定を内蔵: N=8 で (i) MPO 縮約 = 独立構成の H (機械精度),
(ii) DMRG 基底エネルギー = 厳密対角化 (1e-8)。

模型 (占有基底, v20.1 の U(1) core と同一物理・開鎖):
  H = −x Σ (c†_n c_{n+1} + h.c.) + Σ_links L_n² + μ Σ (−1)^n n_n + const(質量規約)
  L_n = Σ_{k≤n} q_k,  q_k = n_k − [k odd] + probe_k
MPO (走和構成, 電場部 bond 3 + ホップ 2 = 全体 5):
  W = [[I, q, c q² + mass], [0, I, 2c q], [0, 0, I]] ⊕ ホップキャリア
"""
import json
import sys
import numpy as np

def site_ops():
    n_op = np.diag([0.0, 1.0])
    cdag = np.array([[0.0, 0.0], [1.0, 0.0]])  # |1⟩⟨0|
    c = cdag.T
    return n_op, cdag, c

def mpo(N, x, mu, probes):
    n_op, cdag, c = site_ops()
    I = np.eye(2)
    Ws = []
    for l in range(N):
        bg = 1.0 if l % 2 == 1 else 0.0
        q = n_op - (bg - probes.get(l, 0.0)) * I
        cl = float(max(N - 1 - l, 0))
        mass = mu * (1.0 if l % 2 == 0 else -1.0) * n_op
        D = 5
        W = np.zeros((D, D, 2, 2))
        W[0, 0] = I
        W[4, 4] = I
        W[0, 4] = cl * (q @ q) + mass
        W[0, 1] = q
        W[1, 1] = I
        W[1, 4] = 2.0 * cl * q
        W[0, 2] = -x * cdag
        W[2, 4] = c
        W[0, 3] = -x * c
        W[3, 4] = cdag
        Ws.append(W)
    Ws[0] = Ws[0][0:1]
    Ws[-1] = Ws[-1][:, 4:5]
    return Ws

def mpo_to_dense(Ws):
    N = len(Ws)
    H = None
    for l in range(N):
        W = Ws[l]  # (Dl, Dr, 2, 2)
        if H is None:
            H = W
        else:
            # H[a?, b, S, S'] ⊗ W: 縮約はボンドで
            H = np.einsum('abst,bcuv->acsutv', H, W)
            s1 = H.shape
            H = H.reshape(s1[0], s1[1], s1[2] * s1[3], s1[4] * s1[5])
    return H[0, 0]

def h_direct(N, x, mu, probes):
    """MPO と独立の直接構成 (自己検定用)"""
    n_op, cdag, c = site_ops()
    I = np.eye(2)
    dim = 2 ** N
    H = np.zeros((dim, dim))
    def kron_at(ops):
        M = np.array([[1.0]])
        for l in range(N):
            M = np.kron(M, ops.get(l, I))
        return M
    for l in range(N - 1):
        # フェルミオン: 開鎖最近接は JW 符号なし (占有基底の kron は site 順)
        H += -x * kron_at({l: cdag, l + 1: c})
        H += -x * kron_at({l: c, l + 1: cdag})
    for l in range(N):
        H += mu * (1.0 if l % 2 == 0 else -1.0) * kron_at({l: n_op})
    # 電場: Σ_links (Σ_{k≤n} q_k)²
    qs = []
    for l in range(N):
        bg = 1.0 if l % 2 == 1 else 0.0
        qs.append(kron_at({l: n_op - (bg - probes.get(l, 0.0)) * I}))
    for n in range(N - 1):
        Ln = sum(qs[k] for k in range(n + 1))
        H += Ln @ Ln
    return H

def lanczos_min(hmul, v0, m=25):
    v = v0 / np.linalg.norm(v0)
    Vs = [v]
    al, be = [], []
    for j in range(m):
        w = hmul(Vs[j])
        a = float(np.dot(Vs[j], w))
        al.append(a)
        w = w - a * Vs[j] - (be[-1] * Vs[j - 1] if j > 0 else 0.0)
        for vv in Vs:  # 再直交
            w -= np.dot(vv, w) * vv
        b = float(np.linalg.norm(w))
        if b < 1e-11 or j == m - 1:
            break
        be.append(b)
        Vs.append(w / b)
    T = np.diag(al)
    if be:
        T += np.diag(be, 1) + np.diag(be, -1)
    ev, evec = np.linalg.eigh(T)
    gv = sum(evec[j, 0] * Vs[j] for j in range(len(Vs)))
    return float(ev[0]), gv / np.linalg.norm(gv)

def dmrg(Ws, chi, sweeps, seed=7):
    N = len(Ws)
    rng = np.random.default_rng(seed)
    # 右正準ランダム MPS — ボンド b_l (site l の左) = min(χ, 2^l, 2^{N−l})
    bond = [int(min(chi, 2 ** min(l, 24), 2 ** min(N - l, 24))) for l in range(N + 1)]
    A = [None] * N
    for l in reversed(range(N)):
        dl, dr = bond[l], bond[l + 1]
        M = rng.normal(size=(dl, 2 * dr))
        Q = np.linalg.qr(M.T)[0].T  # (min(dl, 2dr), 2dr) 行直交
        if Q.shape[0] < dl:
            Q = np.vstack([Q, np.zeros((dl - Q.shape[0], 2 * dr))])
        A[l] = Q.reshape(dl, 2, dr)
    Lc = [None] * (N + 1)
    Rc = [None] * (N + 1)
    Lc[0] = np.ones((1, 1, 1))
    Rc[N] = np.ones((1, 1, 1))
    def upd_L(l):
        Lc[l + 1] = np.einsum('apc,asb,pqst,ctd->bqd', Lc[l], A[l], Ws[l], A[l], optimize=True)
    def upd_R(l):
        Rc[l] = np.einsum('bqd,asb,pqst,ctd->apc', Rc[l + 1], A[l], Ws[l], A[l], optimize=True)
    for l in reversed(range(1, N)):
        upd_R(l)
    energy = None
    for _sw in range(sweeps):
        for direction in (0, 1):
            rng_l = range(N - 1) if direction == 0 else reversed(range(N - 1))
            for l in rng_l:
                Le, Re = Lc[l], Rc[l + 2]
                th0 = np.einsum('asb,btc->astc', A[l], A[l + 1])
                sh = th0.shape
                W1, W2 = Ws[l], Ws[l + 1]
                def hmul(v):
                    # 規約: Lc[a_ket, p, a_bra] / W[p, q, s_ket, s_bra] / Rc[b_ket, r, b_bra]
                    th = v.reshape(sh)  # (a_ket, u, v, b_ket)
                    t = np.einsum('apc,auvb->pcuvb', Le, th, optimize=True)
                    t = np.einsum('pcuvb,pqus->qcsvb', t, W1, optimize=True)
                    t = np.einsum('qcsvb,qrvt->rcstb', t, W2, optimize=True)
                    out = np.einsum('rcstb,brd->cstd', t, Re, optimize=True)
                    return out.reshape(-1)
                e0, gv = lanczos_min(hmul, th0.reshape(-1))
                energy = e0
                th = gv.reshape(sh[0] * 2, 2 * sh[3])
                U, S, Vt = np.linalg.svd(th, full_matrices=False)
                k = int(min(chi, (S > 1e-10).sum() or 1))
                U, S, Vt = U[:, :k], S[:k], Vt[:k]
                S /= np.linalg.norm(S)
                if direction == 0:
                    A[l] = U.reshape(sh[0], 2, k)
                    A[l + 1] = (np.diag(S) @ Vt).reshape(k, 2, sh[3])
                    upd_L(l)
                else:
                    A[l + 1] = Vt.reshape(k, 2, sh[3])
                    A[l] = (U @ np.diag(S)).reshape(sh[0], 2, k)
                    upd_R(l + 1)
    # 密度 ⟨n_l⟩ (恒等環境, O(N²χ³))
    dens = []
    n_op = np.diag([0.0, 1.0])
    for l in range(N):
        Lb = np.ones((1, 1))
        for k2 in range(l):
            Lb = np.einsum('ac,asb,csd->bd', Lb, A[k2], A[k2], optimize=True)
        Rb = np.ones((1, 1))
        for k2 in reversed(range(l + 1, N)):
            Rb = np.einsum('bd,asb,csd->ac', Rb, A[k2], A[k2], optimize=True)
        val = np.einsum('ac,asb,st,ctd,bd->', Lb, A[l], n_op, A[l], Rb, optimize=True)
        nrm = np.einsum('ac,asb,csd,bd->', Lb, A[l], A[l], Rb, optimize=True)
        dens.append(float(val / nrm))
    return energy, dens

def e_field(dens, N, probes):
    e, p = [], 0.0
    for l in range(N):
        bg = 1.0 if l % 2 == 1 else 0.0
        p += dens[l] - bg + probes.get(l, 0.0)
        e.append(p)
    return e

def main():
    out = {}
    # ---- 自己検定 (N=8, x=2, μ=0.5, プローブ付き) ----
    N, x, mu = 8, 2.0, 0.5
    probes = {1: 0.5, 5: -0.5}
    Ws = mpo(N, x, mu, probes)
    Hm = mpo_to_dense(Ws)
    Hd = h_direct(N, x, mu, probes)
    dev = float(np.max(np.abs(Hm - Hd)))
    out["selftest_mpo_dev"] = dev
    ev = np.linalg.eigvalsh(Hd)
    e_dmrg, _ = dmrg(Ws, chi=64, sweeps=8)
    out["selftest_ed_e0"] = float(ev[0])
    out["selftest_dmrg_e0"] = e_dmrg
    print(f"[selftest] MPO dev = {dev:.2e}, ED = {ev[0]:.10f}, DMRG = {e_dmrg:.10f}", file=sys.stderr)
    # ---- アンカー: N=16, x=2, プローブ ±1/2 @ (0,6) の E_mid (v20.2 実測と照合) ----
    N, x = 16, 2.0
    probes = {0: 0.5, 6: -0.5}
    for mu in (0.0, 1.0):
        Wp = mpo(N, x, mu, probes)
        _e, densp = dmrg(Wp, chi=96, sweeps=8)
        W0 = mpo(N, x, mu, {})
        _e0, dens0 = dmrg(W0, chi=96, sweeps=8)
        emid = 0.5 * ((e_field(densp, N, probes)[2] - e_field(dens0, N, {})[2])
                      + (e_field(densp, N, probes)[3] - e_field(dens0, N, {})[3]))
        out[f"anchor_n16_emid_mu{int(mu)}"] = float(emid)
        print(f"[anchor] N=16 mu={mu}: E_mid = {emid:.4f}", file=sys.stderr)
    # ---- 本測定: N=48, プローブ ±1/2 @ (4, 28), 中点ボンド 15/16 ----
    N = 48
    probes = {4: 0.5, 28: -0.5}
    for mu in (0.0, 1.0):
        Wp = mpo(N, x, mu, probes)
        _e, densp = dmrg(Wp, chi=128, sweeps=10)
        W0 = mpo(N, x, mu, {})
        _e0, dens0 = dmrg(W0, chi=128, sweeps=10)
        efp, ef0 = e_field(densp, N, probes), e_field(dens0, N, {})
        emid = 0.5 * ((efp[15] - ef0[15]) + (efp[16] - ef0[16]))
        out[f"n48_emid_mu{int(mu)}"] = float(emid)
        print(f"[main] N=48 mu={mu}: E_mid = {emid:.4f}", file=sys.stderr)
    path = sys.argv[1] if len(sys.argv) > 1 else "explore/dmrg_out.json"
    json.dump(out, open(path, "w"), indent=1)
    print(json.dumps(out, indent=1))

if __name__ == "__main__":
    main()
