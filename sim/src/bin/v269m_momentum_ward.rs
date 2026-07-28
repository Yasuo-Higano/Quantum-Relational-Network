//! v26.9-A v269m_momentum_ward — 運動量セクターの 4D Ward と Belinfante 対称性
//!
//! 事前登録: spec §12.6 (4D covariance closure) — v26.9-0 (energy 行) に続く
//! momentum 行。h₀ᵢ の source T₀ᵢ を構成し、(i) 運動量連続の式から**保存流由来の
//! stress** を導出して BOND-A (metric 変分由来) の stress と tree 照合、
//! (ii) **Belinfante 対称性 T₀y = T_y0 (= J_E)** を on-shell で検証、
//! (iii) 動的 Ward を momentum 行でも機械証明する。
//!
//! 構成 (taste 安全な 2 サイト分割 — 1 サイトシフトは taste 混合、2 サイトは
//! (η² = 1) で自明に singlet):
//!   T₀ᵢ(q) := Σ_x e^{iq(x_y+δ)}·(−i/4)[ψ†_{x+2î}ψ_x − ψ†_xψ_{x+2î}]
//!   ⇒ 折込み頂点 **V₀ᵢ(k;q) = −(1/2)sin(2kᵢ + qδ_{iy})·𝟙** (s 非依存 —
//!   自明に taste-singlet)。node 極限: −(1/2)sin(π+2pᵢ+…) → pᵢ + (q/2)δ_iy ✓
//!   (連続の T₀ᵢ 頂点 = 中点運動量)。q = 0 で [h(k), V₀y] = 0 厳密
//!   (スカラー×恒等) = 全運動量保存。
//! 連続の式 (q ∥ ŷ): ∂_t T⁰ʸ = −∂_y T^{yy} ⇒
//!   **V_S(k;q) := i·[h(k+qŷ)V₀y − V₀y h(k)]/(2 sin(q/2))** が局所 (Fourier 台
//!   有界) なら、V_S = 保存流由来の T_yy 頂点。その on-shell node 極限が
//!   BOND-A の T_yy (= V₀₀ の y 片, つまり h_yy ボンド変調) と Z = 1 で一致
//!   するか — **自由パラメータなしの正準規格同士の照合** (Gate 5 の成分:
//!   「BOND-A source = 保存 stress」の tree 判定)。
//! Belinfante: T₀y (運動量密度) と T_y0 = J_E^y (エネルギー流, v26.9-0 の
//!   V_J) は連続では on-shell で等しい (差は spin 流の発散)。格子では
//!   ‖P₊(k+q)[V₀y − V_J]P₋(k)‖_F / ‖P₊V₀yP₋‖_F → 0 (ε-ladder) を検査
//!   (v26.9-0 の教訓: 縮退があるので基底不変量 = 射影ブロック Frobenius)。
//!
//! 検査 (凍結):
//!  [M0] V₀y(k;0) = −(1/2)sin(2k_y)𝟙・転送ブロック整合・[h, V₀y(0)] = 0 (1e-15)
//!  [M1] 連続の式: V_S の Fourier 台 |n| ≤ 3 (残差 < 1e-12, q ∈ {0.1,0.5,1.0})
//!       かつ q → 0 極限存在
//!  [M2] 保存 stress = BOND-A stress (tree): node ε-ladder で
//!       ‖P₊(V_S − V_yy^A)P₋‖/‖P₊V_S P₋‖ → 0 — V_yy^A は BOND-A の y ボンド
//!       変調頂点 (= V₀₀ の y 片)。Z 補正なし
//!  [M3] Belinfante: ‖P₊(V₀y − V_J)P₋‖/‖P₊V₀yP₋‖ → 0 (ε-ladder, O(ε))
//!  [M4] 動的 Ward (momentum 行): q₀²χ_{0y,0y} + q̂²χ_SS = M₁^{0y} (BZ 積分 < 1e-8)
//!  [M5] 変異: V₀y の h.c. 片を ×1.02 (非 Hermite 化) → [h, V₀y(0)] ≠ 0 で
//!       可分性が破れ V_S(q→0) 発散 (比 > 10)
//!
//! 事前登録分岐: (a) 全 PASS → momentum 行開通 + BOND-A = 保存 stress (tree) +
//!   Belinfante 成立 — 残りは 10×10 kernel と contact 項込み full 4D Ward /
//!   (b) M2 で Z ≠ 1 → BOND-A の yy 正規化と保存 stress の不一致 (公表 — scheme
//!   の再定義が必要になる) / (c) M0/M1 FAIL → 器械。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

fn sbit(s: usize, ax: usize) -> usize {
    (s >> ax) & 1
}

fn h8(k: [f64; 3], m: f64) -> Vec<f64> {
    let mut h = vec![0.0f64; 64];
    for s in 0..8usize {
        let cx = if sbit(s, 0) == 0 { 1.0 } else { -1.0 } * k[0].cos();
        h[s + s * 8] += cx;
        let s2 = s ^ 1;
        let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * k[1].cos();
        h[s2 + s * 8] += cy;
        let s3 = s ^ 3;
        let cz = if sbit(s, 2) == 0 { 1.0 } else { -1.0 } * k[2].cos();
        h[s3 + s * 8] += cz;
        let s4 = s ^ 7;
        h[s4 + s * 8] += m;
    }
    h
}

/// T₀y 頂点: V₀y(k;q) = −(1/2)sin(2k_y + q)·𝟙。
/// mutate: s_y = 1 成分だけ ×1.02 (taste-nonsinglet 運動量密度) — s_y 依存
/// diag は h の Γz/Γm 片 (s_y flip) と非可換になり保存が破れる
/// (開発記録: run1 の「h.c. 片 ×1.02」変異は複素スカラー×𝟙 のままで
/// 依然 h と可換 — 変異として不発だった)。
fn v0y(k: [f64; 3], q: f64, mutate: bool) -> Vec<(f64, f64)> {
    let mut v = vec![(0.0f64, 0.0f64); 64];
    let val = -0.5 * (2.0 * k[1] + q).sin();
    for s in 0..8usize {
        let w = if mutate && sbit(s, 1) == 1 { 1.02 } else { 1.0 };
        v[s + s * 8] = (val * w, 0.0);
    }
    v
}

/// BOND-A の T_yy 頂点 (= V₀₀ の y 片): (−1)^{s_y}cos(k_y + q/2) at s → s^1
fn vyy_bond_a(k: [f64; 3], q: f64) -> Vec<f64> {
    let mut v = vec![0.0f64; 64];
    for s in 0..8usize {
        let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * (k[1] + 0.5 * q).cos();
        v[(s ^ 1) + s * 8] = cy;
    }
    v
}

type C8 = Vec<(f64, f64)>;

fn c8_from_real(a: &[f64]) -> C8 {
    a.iter().map(|&x| (x, 0.0)).collect()
}

fn c8_mul(a: &C8, b: &C8) -> C8 {
    let mut o = vec![(0.0f64, 0.0f64); 64];
    for r in 0..8 {
        for kk in 0..8 {
            let av = a[kk + r * 8];
            if av.0 == 0.0 && av.1 == 0.0 {
                continue;
            }
            for c in 0..8 {
                let bv = b[c + kk * 8];
                o[c + r * 8].0 += av.0 * bv.0 - av.1 * bv.1;
                o[c + r * 8].1 += av.0 * bv.1 + av.1 * bv.0;
            }
        }
    }
    o
}

/// 保存流由来 stress: 格子連続の式 ∂_tρ(q) = i[H,ρ(q)] = +2i sin(q/2)·J(q)
/// (J は ∂_tρ_x + div J = 0 の中点 Fourier 規約) ⇒
///   **V_S(k;q) = C(k;q)/(2 sin(q/2)), C = h(k+qŷ)V₀y − V₀y h(k)** — 実行列。
/// (開発記録: run1 は V_S := iC/q̂ とし、実の BOND-A 頂点との射影ブロック差が
/// rel = √2 で停留 — 位相直交の指紋。i は規約の取り違え。)
fn vs(k: [f64; 3], q: f64, m: f64, mutate: bool) -> C8 {
    let a = c8_from_real(&h8([k[0], k[1] + q, k[2]], m));
    let b = c8_from_real(&h8(k, m));
    let v = v0y(k, q, mutate);
    let am = c8_mul(&a, &v);
    let mb = c8_mul(&v, &b);
    let d = 2.0 * (0.5 * q).sin();
    let mut o = vec![(0.0f64, 0.0f64); 64];
    for i in 0..64 {
        o[i] = ((am[i].0 - mb[i].0) / d, (am[i].1 - mb[i].1) / d);
    }
    o
}

/// v26.9-0 のエネルギー流 (同じ実規約): V_J = [h(k+qŷ)V₀₀ − V₀₀h(k)]/(2sin(q/2))
fn vj_energy(k: [f64; 3], q: f64, m: f64) -> C8 {
    let a = c8_from_real(&h8([k[0], k[1] + q, k[2]], m));
    let b = c8_from_real(&h8(k, m));
    let v = c8_from_real(&h8([k[0], k[1] + 0.5 * q, k[2]], m));
    let am = c8_mul(&a, &v);
    let mb = c8_mul(&v, &b);
    let d = 2.0 * (0.5 * q).sin();
    let mut o = vec![(0.0f64, 0.0f64); 64];
    for i in 0..64 {
        o[i] = ((am[i].0 - mb[i].0) / d, (am[i].1 - mb[i].1) / d);
    }
    o
}

/// 射影ブロック P₊(k+qŷ)·V·P₋(k) の Frobenius ノルム (基底不変量)。
/// P± は固有分解から (縮退部分空間全体の射影なので一意)。
fn proj_block_norm(v: &C8, k: [f64; 3], q: f64, m: f64) -> f64 {
    let hk = h8(k, m);
    let (wk, vk) = jacobi_eigh(&hk, 8);
    let kq = [k[0], k[1] + q, k[2]];
    let hq = h8(kq, m);
    let (wq, vq) = jacobi_eigh(&hq, 8);
    let _ = (wk, wq);
    // ⟨μ (unocc, k+q) | V | ν (occ, k)⟩ を全対で
    let mut sum = 0.0f64;
    for mu in 4..8 {
        for nu in 0..4 {
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for r in 0..8 {
                let (mut sre, mut sim) = (0.0f64, 0.0f64);
                for c in 0..8 {
                    let vkc = vk[c + nu * 8];
                    sre += v[c + r * 8].0 * vkc;
                    sim += v[c + r * 8].1 * vkc;
                }
                let vqr = vq[r + mu * 8];
                re += vqr * sre;
                im += vqr * sim;
            }
            sum += re * re + im * im;
        }
    }
    sum.sqrt()
}

/// 差の射影ブロックノルム (v1 − v2)
fn proj_block_diff(v1: &C8, v2: &C8, k: [f64; 3], q: f64, m: f64) -> f64 {
    let d: C8 = (0..64)
        .map(|i| (v1[i].0 - v2[i].0, v1[i].1 - v2[i].1))
        .collect();
    proj_block_norm(&d, k, q, m)
}

fn main() {
    self_test();
    println!("=== v26.9-A v269m_momentum_ward — 運動量セクターの 4D Ward と Belinfante ===\n");
    println!("V₀y = −(1/2)sin(2k_y+q)𝟙 (2 サイト分割 — 自明 taste-singlet)。保存流由来");
    println!("stress V_S vs BOND-A の T_yy (Z なし正準照合)、Belinfante T₀y = J_E (on-shell)、");
    println!("動的 Ward q₀²χ₀₀ʸ + q̂²χ_SS = M₁ を機械検査。Gate 5 の momentum 行。\n");
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
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
    let m0 = 0.3f64;

    // ---- [M0] 基本恒等式 ----
    {
        let ks = [[0.4f64, -1.0, 2.1], [1.3, 0.6, -0.5]];
        let mut worst = 0.0f64;
        for &k in &ks {
            let v = v0y(k, 0.0, false);
            for s in 0..8usize {
                for s2 in 0..8usize {
                    let expect = if s == s2 { -0.5 * (2.0 * k[1]).sin() } else { 0.0 };
                    worst = worst.max((v[s2 + s * 8].0 - expect).abs().max(v[s2 + s * 8].1.abs()));
                }
            }
            // [h, V₀y(0)] = 0
            let hb = c8_from_real(&h8(k, m0));
            let vv = v0y(k, 0.0, false);
            let c1 = c8_mul(&hb, &vv);
            let c2 = c8_mul(&vv, &hb);
            for i in 0..64 {
                worst = worst.max((c1[i].0 - c2[i].0).abs().max((c1[i].1 - c2[i].1).abs()));
            }
            // 転送ブロック整合: V(k;q)ᵀ = V(k+qŷ;−q) (実対称スカラー)
            let q = 0.7f64;
            let va = v0y(k, q, false);
            let vb = v0y([k[0], k[1] + q, k[2]], -q, false);
            for r in 0..8 {
                for c in 0..8 {
                    worst = worst.max((va[c + r * 8].0 - vb[r + c * 8].0).abs());
                }
            }
        }
        check(
            "[M0] V₀y(k;0) = −(1/2)sin(2k_y)𝟙・[h, V₀y(0)] = 0・転送ブロック整合",
            worst < 1e-15,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [M1] 連続の式: V_S の局所性 ----
    {
        let mut worst = 0.0f64;
        for &q in &[0.1f64, 0.5, 1.0] {
            for ax in 0..3usize {
                let ngrid = 16usize;
                let base = [0.37f64, -0.81, 1.13];
                let vsl: Vec<C8> = (0..ngrid)
                    .map(|j| {
                        let mut k = base;
                        k[ax] = 2.0 * PI * j as f64 / ngrid as f64;
                        vs(k, q, m0, false)
                    })
                    .collect();
                for r in 0..8 {
                    for c in 0..8 {
                        let mut power_hi = 0.0f64;
                        let mut power_all = 0.0f64;
                        for nmode in 0..ngrid {
                            let nsig = if nmode <= ngrid / 2 {
                                nmode as i32
                            } else {
                                nmode as i32 - ngrid as i32
                            };
                            let (mut cre, mut cim) = (0.0f64, 0.0f64);
                            for (j, vm) in vsl.iter().enumerate() {
                                let th = -2.0 * PI * (nmode * j) as f64 / ngrid as f64;
                                let (vre, vim) = vm[c + r * 8];
                                cre += vre * th.cos() - vim * th.sin();
                                cim += vre * th.sin() + vim * th.cos();
                            }
                            let p = (cre * cre + cim * cim) / (ngrid * ngrid) as f64;
                            power_all += p;
                            if nsig.abs() > 3 {
                                power_hi += p;
                            }
                        }
                        if power_all > 1e-20 {
                            worst = worst.max((power_hi / power_all).sqrt());
                        }
                    }
                }
            }
        }
        let k = [0.9f64, 0.2, -1.4];
        let va = vs(k, 1e-4, m0, false);
        let vb = vs(k, 2e-4, m0, false);
        let mut dmax = 0.0f64;
        let mut vnorm = 0.0f64;
        for i in 0..64 {
            dmax = dmax.max((va[i].0 - vb[i].0).abs().max((va[i].1 - vb[i].1).abs()));
            vnorm = vnorm.max(va[i].0.abs().max(va[i].1.abs()));
        }
        check(
            "[M1] 連続の式: V_S = C/(2sin(q/2)) の Fourier 台 |n| ≤ 3 かつ q → 0 極限存在",
            worst < 1e-12 && dmax / vnorm < 1e-3,
            format!("台残差 = {:.1e}, 収束差/ノルム = {:.1e}", worst, dmax / vnorm),
        );
    }

    // ---- [M2] 保存 stress = BOND-A stress (tree, Z なし) ----
    {
        let nvec = [0.55f64, 0.2, -0.81];
        let nn = (nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]).sqrt();
        let n = [nvec[0] / nn, nvec[1] / nn, nvec[2] / nn];
        let c = PI / 2.0;
        let mut rels = Vec::new();
        for &eps in &[0.4f64, 0.2, 0.1, 0.05] {
            let mphys = 0.5 * eps;
            let k = [c + eps * n[0], c + eps * n[1], c + eps * n[2]];
            let q = 0.8 * eps;
            let v_s = vs(k, q, mphys, false);
            let v_a = c8_from_real(&vyy_bond_a(k, q));
            let dnum = proj_block_diff(&v_s, &v_a, k, q, mphys);
            let dden = proj_block_norm(&v_s, k, q, mphys);
            rels.push((eps, dnum / dden));
        }
        let mut msg = String::new();
        for &(e, r) in &rels {
            msg = format!("{} rel({}) = {:.2e}", msg, e, r);
        }
        let ok = rels[3].1 < 0.05 && rels[3].1 < rels[2].1 && rels[2].1 < rels[1].1;
        check(
            "[M2] 保存 stress = BOND-A T_yy (tree, Z なし): 射影ブロック差 → 0 (単調)",
            ok,
            format!("{}", msg),
        );
    }

    // ---- [M3] Belinfante: T₀y vs J_E (on-shell) ----
    {
        let nvec = [0.55f64, 0.2, -0.81];
        let nn = (nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]).sqrt();
        let n = [nvec[0] / nn, nvec[1] / nn, nvec[2] / nn];
        let c = PI / 2.0;
        let mut rels = Vec::new();
        for &eps in &[0.4f64, 0.2, 0.1, 0.05] {
            let mphys = 0.5 * eps;
            let k = [c + eps * n[0], c + eps * n[1], c + eps * n[2]];
            let q = 0.8 * eps;
            let v_p = v0y(k, q, false);
            let v_e = vj_energy(k, q, mphys);
            let dnum = proj_block_diff(&v_p, &v_e, k, q, mphys);
            let dden = proj_block_norm(&v_p, k, q, mphys);
            rels.push((eps, dnum / dden));
        }
        let mut msg = String::new();
        for &(e, r) in &rels {
            msg = format!("{} rel({}) = {:.2e}", msg, e, r);
        }
        let ok = rels[3].1 < 0.05 && rels[3].1 < rels[2].1 && rels[2].1 < rels[1].1;
        check(
            "[M3] Belinfante: on-shell で T₀y = J_E^y (射影ブロック差 → 0, 単調)",
            ok,
            format!("{}", msg),
        );
    }

    // ---- [M4] 動的 Ward (momentum 行, BZ 積分) ----
    {
        let q = 0.4f64;
        let qhat = 2.0 * (0.5 * q).sin();
        let ngrid = 24usize;
        let q0s = [0.3f64, 0.9];
        let mut sums = vec![0.0f64; 5];
        let chunk = ngrid.div_ceil(nthreads);
        let mut partials: Vec<Option<Vec<f64>>> = Vec::new();
        partials.resize_with(nthreads, || None);
        std::thread::scope(|sc| {
            for (t, slot) in partials.iter_mut().enumerate() {
                sc.spawn(move || {
                    let mut acc = vec![0.0f64; 5];
                    for jx in (t * chunk)..(((t + 1) * chunk).min(ngrid)) {
                        let kx = PI * (jx as f64 + 0.5) / ngrid as f64;
                        for jy in 0..ngrid {
                            let ky = PI * (jy as f64 + 0.5) / ngrid as f64;
                            for jz in 0..ngrid {
                                let kz = PI * (jz as f64 + 0.5) / ngrid as f64;
                                let k = [kx, ky, kz];
                                let hk = h8(k, m0);
                                let (wk, vk) = jacobi_eigh(&hk, 8);
                                let kq = [k[0], k[1] + q, k[2]];
                                let hq = h8(kq, m0);
                                let (wq, vq) = jacobi_eigh(&hq, 8);
                                let v0 = v0y(k, q, false);
                                let vsm = vs(k, q, m0, false);
                                for mu in 4..8 {
                                    for nu in 0..4 {
                                        let (mut p_re, mut p_im) = (0.0f64, 0.0f64);
                                        let (mut s_re, mut s_im) = (0.0f64, 0.0f64);
                                        for r in 0..8 {
                                            let (mut a_re, mut a_im, mut b_re, mut b_im) =
                                                (0.0f64, 0.0f64, 0.0f64, 0.0f64);
                                            for cc in 0..8 {
                                                let vkc = vk[cc + nu * 8];
                                                a_re += v0[cc + r * 8].0 * vkc;
                                                a_im += v0[cc + r * 8].1 * vkc;
                                                b_re += vsm[cc + r * 8].0 * vkc;
                                                b_im += vsm[cc + r * 8].1 * vkc;
                                            }
                                            let vqr = vq[r + mu * 8];
                                            p_re += vqr * a_re;
                                            p_im += vqr * a_im;
                                            s_re += vqr * b_re;
                                            s_im += vqr * b_im;
                                        }
                                        let de = wq[mu] - wk[nu];
                                        let a00 = p_re * p_re + p_im * p_im;
                                        let ass = s_re * s_re + s_im * s_im;
                                        acc[0] += 2.0 * a00 * de;
                                        for (ii, &q0) in q0s.iter().enumerate() {
                                            acc[1 + 2 * ii] +=
                                                2.0 * a00 * de / (de * de + q0 * q0);
                                            acc[2 + 2 * ii] +=
                                                2.0 * ass * de / (de * de + q0 * q0);
                                        }
                                    }
                                }
                            }
                        }
                    }
                    *slot = Some(acc);
                });
            }
        });
        for p in partials.into_iter() {
            let p = p.unwrap();
            for i in 0..5 {
                sums[i] += p[i];
            }
        }
        let mut worst = 0.0f64;
        let mut msg = String::new();
        for (ii, &q0) in q0s.iter().enumerate() {
            let lhs = q0 * q0 * sums[1 + 2 * ii] + qhat * qhat * sums[2 + 2 * ii];
            let rel = (lhs / sums[0] - 1.0).abs();
            worst = worst.max(rel);
            msg = format!("{} q₀={}: 残差 {:.1e}", msg, q0, rel);
        }
        check(
            "[M4] 動的 Ward (momentum 行): q₀²χ_{0y,0y} + q̂²χ_SS = M₁ (q = 0.4)",
            worst < 1e-8,
            format!("{} ({} s)", msg, t0.elapsed().as_secs()),
        );
    }

    // ---- [M5] 変異 ----
    {
        let k = [0.9f64, 0.2, -1.4];
        let norm = |v: &C8| -> f64 {
            v.iter()
                .map(|x| (x.0 * x.0 + x.1 * x.1).sqrt())
                .fold(0.0f64, f64::max)
        };
        let r_mut = norm(&vs(k, 1e-3, m0, true)) / norm(&vs(k, 0.1, m0, true));
        let r_good = norm(&vs(k, 1e-3, m0, false)) / norm(&vs(k, 0.1, m0, false));
        check(
            "[M5] 変異: s_y 依存重み ×1.02 (taste-nonsinglet) → 保存破れで V_S(q→0) 発散 (比 > 10; 正版 < 2)",
            r_mut > 10.0 && r_good < 2.0,
            format!("変異比 = {:.1}, 正版比 = {:.2}", r_mut, r_good),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.9-A".into())),
        ("kind".into(), Json::Str("momentum_sector_4d_ward".into())),
        (
            "v0y_construction".into(),
            Json::Str("V0y(k;q) = -(1/2)sin(2k_y+q)·1 (2 サイト分割, 自明 taste-singlet)".into()),
        ),
        (
            "belinfante".into(),
            Json::Str("on-shell P+ (V0y - V_J) P- → 0 (ε-ladder)".into()),
        ),
    ]);
    let p = write_artifact("results/v269m_momentum_ward.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **momentum 行開通 — 保存 stress = BOND-A stress (tree, Z なし)・Belinfante T₀y = J_E 成立・動的 Ward 厳密** (残り = 10×10 kernel と contact 込み full 4D Ward)"
        } else {
            "FAIL あり — 分岐 (b) BOND-A yy 正規化の不一致 (公表) / (c) 器械。欄が一次ソース"
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
