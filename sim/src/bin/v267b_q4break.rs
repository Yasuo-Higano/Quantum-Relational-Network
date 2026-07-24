//! v26.7 (I/II) v267b_q4break — 相互作用による χ_00 の q⁴ 保護の破れ (凍結副実験, PRED-013)
//!
//! 事前登録: paper/grav-vacuum-polarization-spec.md §7 (コミット 157ca53) +
//! predictions.yml PRED-013。プロトコルは凍結済み — **外れてもモデル・V・fit 規則を
//! 交換しない**。v26.3 §3 の機構予言:
//!   自由場では全エネルギー流 J^E_tot が保存量 → エネルギー双極子の行列要素が消え
//!   χ_00 ∝ q⁴。相互作用で J^E_tot の保存が壊れると χ_00 は q² に戻る。
//!
//! 凍結プロトコル (spec §7): 1+1D staggered 鎖 (twist 境界, 半充填, many-body 厳密),
//! 相互作用 V Σ_x n_x n_{x+1} (斥力固定), V ∈ {0, 0.5, 1.0}, N ∈ {10, 12, 14, 16},
//! fit 規則 p = ln[χ_00(q₂)/χ_00(q₁)] / ln(q₂/q₁) (q_j = 2πj/N — 最小 2 運動量点)。
//! **判定 [S4]: p(V=1.0, 最大 N) < 3.0 → PRED-013 的中 / ≥ 3.0 → 外れ (どちらも公表)**。
//!
//! 実行前凍結の補足宣言 (spec §7 が未指定だったパラメタ — 走査しない):
//!   staggered 質量 m = 0.5 を採点対象とする。理由: m = 0 の鎖 + 最近接密度相互作用は
//!   Jordan–Wigner で XXZ (+ 一様場) となり**可積分 — エネルギー流 J^E が厳密保存**
//!   されるため、「相互作用を入れれば保存が壊れる」という予言の前提 (非可積分な一般の
//!   相互作用) が成立しない。staggered 質量 (交替場) は可積分性を破り前提を満たす。
//!   m = 0 の可積分線は**機構の対照 [S5]** として記録する (予想: J^E 保存により
//!   相互作用があっても q⁴ が生き残る — 保存 ⟺ q⁴ の両面検査。採点対象外)。
//!
//! 模型 (JW 前のフェルミオン表示):
//!   H = Σ_x t_x (c†_x c_{x+1} + h.c.)/1 + m Σ_x (−1)^x n_x + V Σ_x n_x n_{x+1}
//!   t_x = 1/2 (x < N−1), t_{N−1} = −1/2 (反周期 twist — 1 粒子ノード k=π/2 を外す)
//!   T_00(x): ボンド項は中点位相の規約 (v26.2 と同一) —
//!   T_00(q) = Σ_bond e^{iq(x+1/2)} [t_x(c†c+h.c.) + V n_x n_{x+1}] + Σ_x e^{iqx} m(−1)^x n_x
//!   構成的に T_00(0) = H。χ_00(q) = 2 Σ_n |⟨n|T_00(q)|0⟩|²/(E_n−E_0) / N (示強)。
//!   解法: Lanczos 基底状態 + 共役勾配 (H−E₀)|x⟩ = P⊥ T|0⟩ (決定的)。
//!
//! 検査 (凍結):
//!  [S0a] V=0 の E₀: many-body Lanczos = 1 粒子 Dirac 海 (dense jacobi) — abs 1e-9
//!  [S0b] V=0 の χ_00(q_{1,2}): many-body CG = 1 粒子 Lehmann — 相対 1e-6
//!        (基底・JW 符号・CG の端から端の器械証明)
//!  [S1] T_00(0) = H: 決定的テストベクトル 3 本で ‖(ΣT−H)v‖ < 1e-12 (V 項込み)
//!  [S2] Lanczos 残差 < 1e-8 / CG 残差 < 1e-10 / 多体ギャップ > 0
//!  [S3] 自由対照: p(V=0, m=0.5, N=16 [採点サイズ]) > 3.3 — q⁴ の分解能。
//!       走行前較正 (自由場の厳密参照のみ・相互作用データ不使用): p_free =
//!       3.23 (N=10) / 3.83 (12) / 3.57 (14) / 3.80 (16) — N=10 の有限 q 補正が
//!       大きいため、ゲートは採点サイズに置く (他は報告)
//!  [S4] **判定 (採点)**: p(V=1.0, m=0.5, N=16) < 3.0 → PRED-013 hit / else miss
//!  [S5] 可積分対照 (記録): m=0 line。走行前較正の発見: 自由場 m=0 は本サイズで
//!       漸近域外 — N≡0 (mod 4) で χ(q₁) が**厳密零** (選択則)、N≡2 (mod 4) で
//!       p≈1.1。よって m=0 は p でなく χ 値を記録する (ゲートなし)
//!  [S6] 変異: T_00 から V 項を落とす → S1 が検出 (逸脱 > 0.01)
//!
//! 事前登録分岐: (a) S0–S2 PASS → S4 が PRED-013 を採点 (hit/miss どちらも判定 a) /
//!   (b) S0/S1/S2 FAIL → 器械の誤り (採点せず) / (c) S3 FAIL → 分解能不足 (採点せず)。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

// ---------------- 多体基底 (半充填の占有ビットマスク, 昇順) ----------------

fn basis_states(n: usize, np: usize) -> Vec<u32> {
    let mut out = Vec::new();
    for s in 0u32..(1u32 << n) {
        if s.count_ones() as usize == np {
            out.push(s);
        }
    }
    out
}

fn find(states: &[u32], s: u32) -> usize {
    states.binary_search(&s).expect("基底外の状態")
}

/// c†_a c_b + c†_b c_a (a < b) の JW 符号: a, b の狭義間の占有パリティ
fn hop_sign(s: u32, a: usize, b: usize) -> f64 {
    let mask = ((1u32 << b) - 1) & !((1u32 << (a + 1)) - 1);
    if (s & mask).count_ones() % 2 == 0 {
        1.0
    } else {
        -1.0
    }
}

struct Chain {
    n: usize,
    m_stag: f64,
    v_int: f64,
    states: Vec<u32>,
}

impl Chain {
    fn new(n: usize, m_stag: f64, v_int: f64) -> Self {
        Chain {
            n,
            m_stag,
            v_int,
            states: basis_states(n, n / 2),
        }
    }
    fn dim(&self) -> usize {
        self.states.len()
    }
    fn tbond(&self, x: usize) -> f64 {
        if x == self.n - 1 {
            -0.5
        } else {
            0.5
        }
    }
    /// H x (実 matvec)
    fn hx(&self, x: &[f64]) -> Vec<f64> {
        self.t00_apply(x, |_| 1.0, |_| 1.0)
    }
    /// 重みつき T_00 適用: ボンド重み wb(x_bond), サイト重み ws(x)。
    /// wb ≡ ws ≡ 1 で H。cos/sin モードで T_00(q) の実部/虚部。
    fn t00_apply(&self, xin: &[f64], wb: impl Fn(usize) -> f64, ws: impl Fn(usize) -> f64) -> Vec<f64> {
        let n = self.n;
        let mut out = vec![0.0f64; self.dim()];
        for (i, &s) in self.states.iter().enumerate() {
            let xi = xin[i];
            if xi == 0.0 {
                continue;
            }
            // 対角: 質量 + 相互作用 (ボンド重み)
            let mut d = 0.0;
            for x in 0..n {
                let nx = (s >> x) & 1;
                if nx == 1 {
                    d += self.m_stag * if x % 2 == 0 { 1.0 } else { -1.0 } * ws(x);
                }
                let x1 = (x + 1) % n;
                if nx == 1 && (s >> x1) & 1 == 1 {
                    d += self.v_int * wb(x);
                }
            }
            out[i] += d * xi;
            // ボンドホップ
            for x in 0..n {
                let x1 = (x + 1) % n;
                let (a, b) = if x < x1 { (x, x1) } else { (x1, x) };
                let na = (s >> a) & 1;
                let nb = (s >> b) & 1;
                if na == nb {
                    continue;
                }
                let s2 = s ^ ((1u32 << a) | (1u32 << b));
                let j = find(&self.states, s2);
                out[j] += self.tbond(x) * hop_sign(s, a, b) * wb(x) * xi;
            }
        }
        out
    }
    /// T_00(q) の (cos, sin) モードを |ψ⟩ に適用 (中点位相規約)
    fn t00q(&self, psi: &[f64], q: f64, sin_mode: bool) -> Vec<f64> {
        let f = move |u: f64| if sin_mode { (q * u).sin() } else { (q * u).cos() };
        self.t00_apply(psi, move |x| f(x as f64 + 0.5), move |x| f(x as f64))
    }
}

/// 共役勾配: (H − e0) x = r を |gs⟩ 直交補空間で解く。戻り値 (⟨r|x⟩, 残差)
fn cg_resolvent(ch: &Chain, e0: f64, gs: &[f64], r0: &[f64]) -> (f64, f64) {
    let dim = ch.dim();
    let proj = |v: &mut Vec<f64>| {
        let p: f64 = gs.iter().zip(v.iter()).map(|(a, b)| a * b).sum();
        for i in 0..dim {
            v[i] -= p * gs[i];
        }
    };
    let mut r = r0.to_vec();
    proj(&mut r);
    let rr0: f64 = r.iter().map(|x| x * x).sum();
    if rr0 < 1e-28 {
        return (0.0, 0.0);
    }
    let mut x = vec![0.0f64; dim];
    let mut p = r.clone();
    let mut rr = rr0;
    for _it in 0..5000 {
        let mut ap = ch.hx(&p);
        for i in 0..dim {
            ap[i] -= e0 * p[i];
        }
        proj(&mut ap);
        let pap: f64 = p.iter().zip(ap.iter()).map(|(a, b)| a * b).sum();
        let alpha = rr / pap;
        for i in 0..dim {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        proj(&mut r);
        let rr_new: f64 = r.iter().map(|x| x * x).sum();
        if rr_new < 1e-26 * rr0.max(1.0) {
            rr = rr_new;
            break;
        }
        let beta = rr_new / rr;
        rr = rr_new;
        for i in 0..dim {
            p[i] = r[i] + beta * p[i];
        }
    }
    let val: f64 = r0.iter().zip(x.iter()).map(|(a, b)| a * b).sum();
    (val, rr.sqrt())
}

/// 基底状態 (Lanczos, k=2 でギャップも) — 実 matvec を複素に包む
fn ground_state(ch: &Chain) -> (f64, f64, Vec<f64>, f64) {
    let dim = ch.dim();
    let mv = |v: &[(f64, f64)]| -> Vec<(f64, f64)> {
        let re: Vec<f64> = v.iter().map(|x| x.0).collect();
        let im: Vec<f64> = v.iter().map(|x| x.1).collect();
        let hre = ch.hx(&re);
        let him = ch.hx(&im);
        hre.into_iter().zip(him).collect()
    };
    let krylov = 240.min(dim);
    let (w, vecs, resid) = lanczos_lowest_herm(&mv, dim, 2, krylov, 12345);
    // 実ベクトル化 (実対称 H — 位相を実に回す)
    let v0 = &vecs[0];
    let (mut nr, mut ni) = (0.0f64, 0.0f64);
    for x in v0 {
        nr += x.0 * x.0;
        ni += x.1 * x.1;
    }
    let gs: Vec<f64> = if nr >= ni {
        v0.iter().map(|x| x.0).collect()
    } else {
        v0.iter().map(|x| x.1).collect()
    };
    let nrm: f64 = gs.iter().map(|x| x * x).sum::<f64>().sqrt();
    let gs: Vec<f64> = gs.iter().map(|x| x / nrm).collect();
    (w[0], w[1], gs, resid)
}

/// χ_00(q)/N: CG で 2[⟨r_c|(H−E0)⁻¹|r_c⟩ + ⟨r_s|(H−E0)⁻¹|r_s⟩]/N
fn chi00(ch: &Chain, e0: f64, gs: &[f64], q: f64) -> (f64, f64) {
    let mut chi = 0.0;
    let mut worst = 0.0f64;
    for sin_mode in [false, true] {
        let mut r = ch.t00q(gs, q, sin_mode);
        // ⟨0|T|0⟩ 成分は CG 内の射影で除かれる
        let (val, res) = cg_resolvent(ch, e0, gs, &r);
        chi += 2.0 * val;
        worst = worst.max(res);
        r.clear();
    }
    (chi / ch.n as f64, worst)
}

// ---------------- 1 粒子参照 (V = 0) ----------------

fn sp_hamiltonian(n: usize, m_stag: f64) -> Vec<f64> {
    let mut h = vec![0.0f64; n * n];
    for x in 0..n {
        let x1 = (x + 1) % n;
        let t = if x == n - 1 { -0.5 } else { 0.5 };
        h[x1 + x * n] += t;
        h[x + x1 * n] += t;
        h[x + x * n] += m_stag * if x % 2 == 0 { 1.0 } else { -1.0 };
    }
    h
}

/// 1 粒子 T_00(q) 頂点 (複素, 中点位相) と Dirac 海 Lehmann
fn sp_chi00(n: usize, m_stag: f64, q: f64) -> f64 {
    let h = sp_hamiltonian(n, m_stag);
    let (w, v) = jacobi_eigh(&h, n);
    let mut re = vec![0.0f64; n * n];
    let mut im = vec![0.0f64; n * n];
    for x in 0..n {
        let x1 = (x + 1) % n;
        let t = if x == n - 1 { -0.5 } else { 0.5 };
        let ph = q * (x as f64 + 0.5);
        re[x1 + x * n] += t * ph.cos();
        re[x + x1 * n] += t * ph.cos();
        im[x1 + x * n] += t * ph.sin();
        im[x + x1 * n] += t * ph.sin();
        let phs = q * x as f64;
        let ms = m_stag * if x % 2 == 0 { 1.0 } else { -1.0 };
        re[x + x * n] += ms * phs.cos();
        im[x + x * n] += ms * phs.sin();
    }
    let nocc = n / 2;
    let mut chi = 0.0;
    for mu in nocc..n {
        for nu in 0..nocc {
            let (mut mr, mut mi) = (0.0f64, 0.0f64);
            for k in 0..n {
                for l in 0..n {
                    mr += v[k + mu * n] * re[l + k * n] * v[l + nu * n];
                    mi += v[k + mu * n] * im[l + k * n] * v[l + nu * n];
                }
            }
            chi += 2.0 * (mr * mr + mi * mi) / (w[mu] - w[nu]);
        }
    }
    chi / n as f64
}

fn sp_e0(n: usize, m_stag: f64) -> f64 {
    let h = sp_hamiltonian(n, m_stag);
    let (w, _) = jacobi_eigh(&h, n);
    w[..n / 2].iter().sum()
}

fn main() {
    self_test();
    println!("=== v26.7 (I/II) v267b_q4break — 相互作用による χ_00 の q⁴ 保護の破れ (PRED-013) ===\n");
    println!("凍結: spec §7 (157ca53) + PRED-013。判定 [S4]: p(V=1.0, m=0.5, N=16) < 3.0 → hit。");
    println!("補足宣言 (実行前凍結): 採点対象は m = 0.5 (可積分性破壊)。m = 0 は XXZ 可積分線");
    println!("(J^E 厳密保存) の機構対照 [S5] — 予想は q⁴ 残存。走査・交換はしない。\n");
    let t0 = std::time::Instant::now();
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!(
            "  [{}] {}  {}",
            if ok { "PASS" } else { "FAIL" },
            name,
            detail
        );
        if !ok {
            nfail += 1;
        }
    };

    let ns_list = [10usize, 12, 14, 16];
    let vs = [0.0f64, 0.5, 1.0];
    let ms = [0.5f64, 0.0]; // [0] = 採点対象 (非可積分), [1] = 可積分対照

    // ---- [S1] T_00(0) = H (V 項込み) + [S6] 変異 ----
    {
        let ch = Chain::new(10, 0.5, 1.0);
        let dim = ch.dim();
        let mut worst = 0.0f64;
        let mut worst_mut = f64::INFINITY;
        for seed in [0usize, dim / 2, dim - 1] {
            let mut v = vec![0.0f64; dim];
            v[seed] = 1.0;
            let hv = ch.hx(&v);
            let tv = ch.t00q(&v, 0.0, false);
            let d: f64 = hv
                .iter()
                .zip(&tv)
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt();
            worst = worst.max(d);
            // 変異: T_00 から相互作用ボンド項を落とす (wb を質量のみに効かせる)
            let tv_bad = ch.t00_apply(&v, |_| 1.0, |_| 1.0); // 正
            let tv_bad2 = {
                // V 項を落とした T_00(0): ボンド重みで V だけ 0 にはできないので
                // V=0 の鎖の T_00(0) を V≠0 の H と比べる (エネルギー密度の取りこぼし)
                let ch0 = Chain::new(10, 0.5, 0.0);
                ch0.t00q(&v, 0.0, false)
            };
            let dmut: f64 = hv
                .iter()
                .zip(&tv_bad2)
                .map(|(a, b)| (a - b) * (a - b))
                .sum::<f64>()
                .sqrt();
            worst_mut = worst_mut.min(dmut);
            let _ = tv_bad;
        }
        check(
            "[S1] T_00(0) = H (V 項込み, 決定的ベクトル 3 本)",
            worst < 1e-12,
            format!("max‖Δv‖ = {:.1e}", worst),
        );
        check(
            "[S6] 変異: T_00 から V 項を落とす → S1 が検出",
            worst_mut > 0.01,
            format!("min 逸脱 {:.3} > 0.01", worst_mut),
        );
    }

    // ---- 走査 (全条件) ----
    // tab[mi][vi][ni] = (E0, gap, chi1, chi2, p)
    let mut tab = vec![vec![vec![(0.0f64, 0.0f64, 0.0f64, 0.0f64, 0.0f64); ns_list.len()]; vs.len()]; 2];
    let mut worst_lanczos = 0.0f64;
    let mut worst_cg = 0.0f64;
    let mut min_gap_scored = f64::INFINITY; // 採点線 (m=0.5) のみ — m=0 対照線は縮退あり
    let mut min_gap_ctrl = f64::INFINITY;
    println!("    [走査] m | V | N | dim | E₀ | gap | χ(q₁)/N | χ(q₂)/N | p");
    for (mi, &m) in ms.iter().enumerate() {
        for (vi, &v_int) in vs.iter().enumerate() {
            for (ni, &n) in ns_list.iter().enumerate() {
                let ch = Chain::new(n, m, v_int);
                let (e0, e1, gs, resid) = ground_state(&ch);
                worst_lanczos = worst_lanczos.max(resid);
                if mi == 0 {
                    min_gap_scored = min_gap_scored.min(e1 - e0);
                } else {
                    min_gap_ctrl = min_gap_ctrl.min(e1 - e0);
                }
                let q1 = 2.0 * PI / n as f64;
                let q2 = 4.0 * PI / n as f64;
                let (c1, r1) = chi00(&ch, e0, &gs, q1);
                let (c2, r2) = chi00(&ch, e0, &gs, q2);
                worst_cg = worst_cg.max(r1).max(r2);
                // 厳密零 (m=0 の選択則) で ln が発散しないようクランプ (m=0.5 系は無関係)
                let p = (c2.max(1e-300) / c1.max(1e-300)).ln() / (q2 / q1).ln();
                tab[mi][vi][ni] = (e0, e1 - e0, c1, c2, p);
                println!(
                    "      m={:.1} V={:.1} N={:2} dim={:5}: E₀={:+.6} gap={:.4} χ₁={:.6e} χ₂={:.6e} p={:.3}",
                    m,
                    v_int,
                    n,
                    ch.dim(),
                    e0,
                    e1 - e0,
                    c1,
                    c2,
                    p
                );
            }
        }
        println!();
    }
    println!("    ({} s)", t0.elapsed().as_secs());

    // ---- [S0] V=0 の器械照合 ----
    // 開発記録 (run1 → run2): S0b 初版は純相対比較 — m=0 の χ(q₁) 厳密零 (選択則) で
    // |0/1e-31 − 1| = 1.0 と誤発報した (ゲートの前提誤り)。混合許容
    // |Δ| < max(1e-6·|参照|, 1e-9) に較正。物理数値は run1/run2 で同一 (決定的)。
    {
        let mut worst_e = 0.0f64;
        let mut worst_chi = 0.0f64;
        for (mi, &m) in ms.iter().enumerate() {
            for (ni, &n) in ns_list.iter().enumerate() {
                let (e0, _, c1, c2, _) = tab[mi][0][ni];
                worst_e = worst_e.max((e0 - sp_e0(n, m)).abs());
                let s1 = sp_chi00(n, m, 2.0 * PI / n as f64);
                let s2 = sp_chi00(n, m, 4.0 * PI / n as f64);
                worst_chi = worst_chi
                    .max((c1 - s1).abs() / (1e-6 * s1.abs()).max(1e-9))
                    .max((c2 - s2).abs() / (1e-6 * s2.abs()).max(1e-9));
            }
        }
        check(
            "[S0a] V=0 の E₀: many-body Lanczos = 1 粒子 Dirac 海 (abs 1e-9)",
            worst_e < 1e-9,
            format!("max|Δ| = {:.1e}", worst_e),
        );
        check(
            "[S0b] V=0 の χ_00(q₁,q₂): many-body CG = 1 粒子 Lehmann (混合許容 max(1e-6 相対, 1e-9))",
            worst_chi < 1.0,
            format!("worst = {:.3} × 許容", worst_chi),
        );
    }

    // ---- [S2] 数値残差 ----
    // 開発記録 (run1 → run2): gap ゲート初版は m=0 対照線 (縮退基底あり — twist 下の
    // 中立線) を含んでいた。採点線 (m=0.5, gap ~ 1.0–1.8) に限定し、対照線は報告のみ。
    check(
        "[S2] Lanczos 残差 < 1e-8 / CG 残差 < 1e-10 / 採点線 (m=0.5) の多体ギャップ > 0.1",
        worst_lanczos < 1e-8 && worst_cg < 1e-10 && min_gap_scored > 0.1,
        format!(
            "Lanczos {:.1e} / CG {:.1e} / 採点線 min gap {:.4} (対照線 {:.1e} — 縮退あり・報告のみ)",
            worst_lanczos, worst_cg, min_gap_scored, min_gap_ctrl
        ),
    );

    // ---- [S3] 自由対照 (q⁴ が採点サイズで見えるか — 走行前較正: 厳密参照 3.80) ----
    {
        let p_free = tab[0][0][ns_list.len() - 1].4;
        println!(
            "    [S3 表] p(V=0, m=0.5) = {:.3} / {:.3} / {:.3} / {:.3} (N=10,12,14,16 — 厳密参照 3.23/3.83/3.57/3.80)",
            tab[0][0][0].4, tab[0][0][1].4, tab[0][0][2].4, tab[0][0][3].4
        );
        check(
            "[S3] 自由対照: p(V=0, m=0.5, N=16) > 3.3 (採点サイズの q⁴ 分解能)",
            p_free > 3.3,
            format!("p = {:.3}", p_free),
        );
    }

    // ---- [S4] 判定 (PRED-013 の採点) ----
    let p_scored = tab[0][2][ns_list.len() - 1].4;
    {
        let hit = p_scored < 3.0;
        println!(
            "    [S4 採点] p(V=1.0, m=0.5, N=16) = {:.3} — 凍結判定 p < 3.0 ⇒ **PRED-013 {}**",
            p_scored,
            if hit { "hit (的中)" } else { "miss (外れ)" }
        );
        println!(
            "      V 系列 (m=0.5, N=16): p = {:.3} (V=0) → {:.3} (V=0.5) → {:.3} (V=1.0)",
            tab[0][0][3].4, tab[0][1][3].4, tab[0][2][3].4
        );
        check(
            "[S4] PRED-013 の採点が確定 (hit/miss どちらも判定 a — モデル交換なし)",
            true,
            format!("p = {:.3} ⇒ {}", p_scored, if hit { "hit" } else { "miss" }),
        );
    }

    // ---- [S5] 可積分対照 (記録 — 走行前較正の発見により χ 値で記録) ----
    {
        println!("    [S5 可積分対照] m=0 (XXZ 線): 本サイズでは漸近域外 (自由場で χ(q₁) の厳密零");
        println!("      [N≡0 mod 4] / p≈1.1 [N≡2 mod 4] — 走行前の厳密参照で判明)。χ 値の記録:");
        for (vi, &v_int) in vs.iter().enumerate() {
            println!(
                "      V={:.1}: χ(q₁)/N = {:.3e} (N=14) / {:.3e} (N=16), χ(q₂)/N = {:.3e} (N=16)",
                v_int, tab[1][vi][2].2, tab[1][vi][3].2, tab[1][vi][3].3
            );
        }
        println!("      — 可積分線の q⁴ 判定は大 N の別器械が要る (記録のみ・採点外)");
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.7-I".into())),
        ("kind".into(), Json::Str("q4_protection_break".into())),
        ("pred".into(), Json::Str("PRED-013".into())),
        ("p_scored_v1_m05_n16".into(), Json::Num(p_scored)),
        (
            "p_table_m05".into(),
            Json::Arr(
                (0..vs.len())
                    .map(|vi| {
                        Json::Obj(vec![
                            ("V".into(), Json::Num(vs[vi])),
                            (
                                "p_by_n".into(),
                                Json::Arr(
                                    (0..ns_list.len()).map(|ni| Json::Num(tab[0][vi][ni].4)).collect(),
                                ),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "chi_table_m0_integrable".into(),
            Json::Arr(
                (0..vs.len())
                    .map(|vi| {
                        Json::Obj(vec![
                            ("V".into(), Json::Num(vs[vi])),
                            (
                                "chi1_by_n".into(),
                                Json::Arr(
                                    (0..ns_list.len()).map(|ni| Json::Num(tab[1][vi][ni].2)).collect(),
                                ),
                            ),
                            (
                                "chi2_by_n".into(),
                                Json::Arr(
                                    (0..ns_list.len()).map(|ni| Json::Num(tab[1][vi][ni].3)).collect(),
                                ),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v267b_q4break.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **S4 の採点が PRED-013 の status を確定** (hit/miss どちらも公表 — モデル交換なし)"
        } else {
            "FAIL — 分岐 (b) 器械 / (c) 分解能。採点せず"
        }
    );
    println!(
        "\n総合判定: {} ({} s)",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" },
        t0.elapsed().as_secs()
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
