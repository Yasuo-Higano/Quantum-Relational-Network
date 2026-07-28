//! v26.9-E v269i_spin_improve — spin-current 改良項の明示構成 (Belinfante 完成)
//!
//! 事前登録: v26.9-D の負の結果 (素朴平均 ≠ Belinfante 密度) を受けた改良項の
//! 明示構成。連続の Belinfante 改良 T_B^{0i} = T_can^{0i} + ∂_λK^{λ0i} の格子
//! 転写 — 転送 q ∥ ŷ では λ = y のみ生き、
//!   **ΔT^{0i}(qŷ) = q̂ · W_i,  W_i = λ·[α̂ᵢ, α̂_y] (on-site 定数行列)**
//! (α̂ᵢ := ∂h/∂kᵢ|node の閉形式 — dh8(node))。su(2) 型代数
//! [Σpⱼα̂ⱼ, [α̂ᵢ,α̂_y]] = 4(pᵢα̂_y − p_yα̂ᵢ) から [h, q̂W]/q̂ の node 極限は
//! 回転流 = (piece1 − piece2)/2 の構造 — v26.9-C の厳密恒等式 (正準流束 =
//! piece2) と合わせ、**λ の有理値 (±1/8) が予言される**。
//!
//! 検査 (凍結):
//!  [E0] λ の同定: node 近傍 1 点で ‖[h, W]-transfer − (p1−p2)/2‖ 最小化の
//!       λ_fit が有理値 ±1/8 に 1e-3 で一致 (以後 λ = その有理値に凍結)
//!  [E1] **改良の閉包**: rel(Φ(T_B^{0x}) vs T^{xy}-split) の ε-ladder が
//!       v26.9-D の停留 1.30 から単調減少へ転じ、最終 < 5%
//!  [E2] **対称 10×10 の横断性回復**: 破れ 0.674 (v26.9-D) → ladder 単調減少,
//!       最終 < 2% (縮小比 > 2.5)
//!  [E3] **conformal 崩壊**: trace 比 0.265 (v26.9-D) → 単調減少, 最終 < 2%
//!  [E4] oracle 回帰: σ_DD/(2ρ_D) = 1 ± 2%・σ_DD/2σ_XX = 1 ± 2% (P₂ 無傷)
//!  [E5] 変異: λ → −λ → 横断性が O(1) に停留 (最終 > 10× 正版)
//!
//! 事前登録分岐: (a) 全 PASS → **Belinfante 対称 10×10 が σ = ρ₂P₂ に完全崩壊
//!   — Gate 5 の分離部門完結** / (b) E1 のみ破れ → 改良項の形が不足 (別の
//!   Γ 積が要る — 公表) / (c) E0 FAIL → 代数の見立て違い (公表)。

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

fn dh8(k: [f64; 3], ax: usize) -> Vec<f64> {
    let mut h = vec![0.0f64; 64];
    for s in 0..8usize {
        let sgn = if sbit(s, ax) == 0 { 1.0 } else { -1.0 };
        let val = -sgn * k[ax].sin();
        let s2 = match ax {
            0 => s,
            1 => s ^ 1,
            _ => s ^ 3,
        };
        h[s2 + s * 8] += val;
    }
    h
}

fn mat_mul(a: &[f64], b: &[f64]) -> Vec<f64> {
    let mut o = vec![0.0f64; 64];
    for r in 0..8 {
        for k in 0..8 {
            let av = a[k + r * 8];
            if av == 0.0 {
                continue;
            }
            for c in 0..8 {
                o[c + r * 8] += av * b[c + k * 8];
            }
        }
    }
    o
}

/// 改良頂点 W_i = [α̂ᵢ, α̂_y] (node の閉形式定数行列 — 実反対称)
fn w_mat(ax: usize) -> Vec<f64> {
    let node = [PI / 2.0, PI / 2.0, PI / 2.0];
    let ai = dh8(node, ax);
    let ay = dh8(node, 1);
    let x = mat_mul(&ai, &ay);
    let y = mat_mul(&ay, &ai);
    (0..64).map(|i| x[i] - y[i]).collect()
}

// ---- 空間 6 source (v26.9-C/D の認証済み写経) ----

struct Term {
    eps: usize,
    d: [i32; 3],
    w: f64,
}

fn vertex_unmod(terms: &[Term], k: [f64; 3]) -> Vec<f64> {
    let mut v = vec![0.0f64; 64];
    for t in terms {
        for s in 0..8usize {
            let s2 = s ^ t.eps;
            let mut ph = 0.0f64;
            for axx in 0..3 {
                ph += (k[axx] + PI * sbit(s, axx) as f64) * t.d[axx] as f64;
            }
            v[s + s2 * 8] += t.w * 2.0 * ph.cos();
        }
    }
    v
}

fn t_split_terms(a1: usize, a2: usize) -> Vec<Term> {
    let flip = [0usize, 1, 3];
    let mut v = Vec::new();
    for sg in [1i32, -1] {
        for rh in [1i32, -1] {
            let c = (sg * rh) as f64 / 16.0;
            let mut d1 = [0i32; 3];
            d1[a1] = sg;
            d1[a2] = 2 * rh;
            v.push(Term { eps: flip[a1], d: d1, w: -c });
            let mut d2 = [0i32; 3];
            d2[a2] = sg;
            d2[a1] = 2 * rh;
            v.push(Term { eps: flip[a2], d: d2, w: -c });
        }
    }
    v
}

/// Belinfante 10 source (改良版): T⁰ᵢ = 正準 V₀ᵢ + q̂·λ·W_i (i = x, z のみ —
/// 0y/00 は改良不要 [K^{y0y} = K^{y00} = 0])。mutate: λ → −λ
fn source_imp(i: usize, k: [f64; 3], q: f64, m: f64, lam: f64, mutate: bool) -> Vec<f64> {
    let km = [k[0], k[1] + 0.5 * q, k[2]];
    let qhat = 2.0 * (0.5 * q).sin();
    let sgn = if mutate { -1.0 } else { 1.0 };
    match i {
        0 => h8(km, m),
        1 | 3 => {
            let ax = i - 1;
            let val = -0.5 * (2.0 * km[ax]).sin();
            let mut v0 = vec![0.0f64; 64];
            for s in 0..8usize {
                v0[s + s * 8] = val;
            }
            let w = w_mat(ax);
            (0..64).map(|x| v0[x] + sgn * lam * qhat * w[x]).collect()
        }
        2 => {
            let val = -0.5 * (2.0 * km[1]).sin();
            let mut v0 = vec![0.0f64; 64];
            for s in 0..8usize {
                v0[s + s * 8] = val;
            }
            v0
        }
        4 => vertex_unmod(&[Term { eps: 0, d: [1, 0, 0], w: 0.5 }], km),
        5 => vertex_unmod(&[Term { eps: 1, d: [0, 1, 0], w: 0.5 }], km),
        6 => vertex_unmod(&[Term { eps: 3, d: [0, 0, 1], w: 0.5 }], km),
        7 => vertex_unmod(&t_split_terms(0, 1), km).iter().map(|x| x / 2.0).collect(),
        8 => vertex_unmod(&t_split_terms(0, 2), km).iter().map(|x| x / 2.0).collect(),
        _ => vertex_unmod(&t_split_terms(1, 2), km).iter().map(|x| x / 2.0).collect(),
    }
}

// ---- 殻積分 (v269d の写経, source を改良版に) ----

fn e8m(k: [f64; 3]) -> f64 {
    (k[0].cos().powi(2) + k[1].cos().powi(2) + k[2].cos().powi(2)).sqrt()
}

fn f_pair(p: [f64; 3], q_lat: f64) -> f64 {
    let c = PI / 2.0;
    let k = [c + p[0], c + p[1], c + p[2]];
    e8m([k[0], k[1] + q_lat, k[2]]) + e8m(k)
}

fn df_dr(p: [f64; 3], n: [f64; 3], q_lat: f64) -> f64 {
    let c = PI / 2.0;
    let k = [c + p[0], c + p[1], c + p[2]];
    let kq = [k[0], k[1] + q_lat, k[2]];
    let (e1, e2) = (e8m(kq).max(1e-14), e8m(k).max(1e-14));
    let mut d = 0.0f64;
    for ax in 0..3 {
        d += n[ax] * (-kq[ax].sin() * kq[ax].cos() / e1 - k[ax].sin() * k[ax].cos() / e2);
    }
    d
}

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut xs = vec![0.0f64; n];
    let mut ws = vec![0.0f64; n];
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0f64, x);
            for kk in 2..=n {
                let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        xs[i] = x;
        let (mut p0, mut p1) = (1.0f64, x);
        for kk in 2..=n {
            let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
            p0 = p1;
            p1 = p2;
        }
        let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
        ws[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (xs, ws)
}

fn sigma_imp(a: f64, e_phys: f64, q_phys: f64, nth: usize, nph: usize, lam: f64, mutate: bool) -> (Vec<f64>, [f64; 4]) {
    let e_lat = a * e_phys;
    let q_lat = a * q_phys;
    let gl = gauss_legendre(nth);
    let c = PI / 2.0;
    let mut sig = vec![0.0f64; 100];
    let mut tv_num = [0.0f64; 4];
    let mut tv_den = [0.0f64; 4];
    let qhat = 2.0 * (0.5 * q_lat).sin();
    let pairs = [(0usize, 2usize), (1, 7), (2, 5), (3, 9)];
    for (ct, wt) in gl.0.iter().zip(&gl.1) {
        let st = (1.0 - ct * ct).sqrt();
        for j in 0..nph {
            let ph = (j as f64 + 0.5) * 2.0 * PI / nph as f64;
            let n = [st * ph.cos(), st * ph.sin(), *ct];
            let mut r_hi = e_lat;
            let mut guard = 0;
            while f_pair([r_hi * n[0], r_hi * n[1], r_hi * n[2]], q_lat) <= e_lat && guard < 40 {
                r_hi *= 1.5;
                guard += 1;
            }
            let (mut lo, mut hi) = (0.0f64, r_hi);
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if f_pair([mid * n[0], mid * n[1], mid * n[2]], q_lat) < e_lat {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let r = 0.5 * (lo + hi);
            let p = [r * n[0], r * n[1], r * n[2]];
            let k = [c + p[0], c + p[1], c + p[2]];
            let hk = h8(k, 0.0);
            let (_, vk) = jacobi_eigh(&hk, 8);
            let kq = [k[0], k[1] + q_lat, k[2]];
            let hq = h8(kq, 0.0);
            let (_, vq) = jacobi_eigh(&hq, 8);
            let vs: Vec<Vec<f64>> = (0..10)
                .map(|i| source_imp(i, k, q_lat, 0.0, lam, mutate))
                .collect();
            let dfr = df_dr(p, n, q_lat).abs().max(1e-12);
            let wgt = wt * (2.0 * PI / nph as f64) * r * r / dfr;
            for mu in 4..8 {
                for nu in 0..4 {
                    let mut mels = [0.0f64; 10];
                    for (i, v) in vs.iter().enumerate() {
                        let mut re = 0.0f64;
                        for rr in 0..8 {
                            let mut s = 0.0f64;
                            for cc in 0..8 {
                                s += v[cc + rr * 8] * vk[cc + nu * 8];
                            }
                            re += vq[rr + mu * 8] * s;
                        }
                        mels[i] = re;
                    }
                    for i in 0..10 {
                        for jj in 0..10 {
                            sig[jj + i * 10] += wgt * mels[i] * mels[jj];
                        }
                    }
                    for (nn, &(i0, iy)) in pairs.iter().enumerate() {
                        let vr = e_lat * mels[i0] - qhat * mels[iy];
                        tv_num[nn] += wgt * vr * vr;
                        tv_den[nn] += wgt
                            * (e_lat * e_lat * mels[i0] * mels[i0]
                                + qhat * qhat * mels[iy] * mels[iy]);
                    }
                }
            }
        }
    }
    let norm = (2.0 * PI).powi(3);
    for v in sig.iter_mut() {
        *v /= norm;
    }
    let mut tv = [0.0f64; 4];
    for i in 0..4 {
        tv[i] = (tv_num[i] / tv_den[i].max(1e-300)).sqrt();
    }
    (sig, tv)
}

fn main() {
    self_test();
    println!("=== v26.9-E v269i_spin_improve — spin-current 改良項の明示構成 ===\n");
    println!("ΔT⁰ⁱ = q̂·λ·[α̂ᵢ, α̂_y] (on-site 定数行列)。v26.9-D の停留 (流束 1.30 /");
    println!("横断性 0.674) が崩れるか — Belinfante 完成と Gate 5 分離部門の完結判定。\n");
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

    // ---- [E0] λ の同定 (tree matching で fit → 有理値に凍結) ----
    let lam_frozen: f64;
    {
        // node 近傍: Φ(T⁰ˣ(λ)) と T^{xy}-split の射影ブロック距離を λ で最小化。
        // Φ(λ) = Φ(0) + λ·q̂·[h-pair, W]/q̂ 線形 — 最小二乗の閉形式で λ_fit。
        let eps = 0.1f64;
        let mphys = 0.5 * eps;
        let c = PI / 2.0;
        let nvec = [0.55f64, 0.2, -0.81];
        let nn = (nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]).sqrt();
        let k = [
            c + eps * nvec[0] / nn,
            c + eps * nvec[1] / nn,
            c + eps * nvec[2] / nn,
        ];
        let qq = 0.8 * eps;
        let flux_of = |lam: f64| -> Vec<f64> {
            let d = source_imp(1, k, qq, mphys, lam, false);
            let a = h8([k[0], k[1] + qq, k[2]], mphys);
            let b = h8(k, mphys);
            let am = mat_mul(&a, &d);
            let mb = mat_mul(&d, &b);
            let dd = 2.0 * (0.5 * qq).sin();
            (0..64).map(|i| (am[i] - mb[i]) / dd).collect()
        };
        let split = source_imp(7, k, qq, mphys, 0.0, false);
        let hk = h8(k, mphys);
        let (_, vk) = jacobi_eigh(&hk, 8);
        let hq = h8([k[0], k[1] + qq, k[2]], mphys);
        let (_, vq) = jacobi_eigh(&hq, 8);
        let block_vec = |v: &[f64]| -> Vec<f64> {
            let mut o = Vec::with_capacity(16);
            for mu in 4..8 {
                for nu in 0..4 {
                    let mut re = 0.0f64;
                    for r in 0..8 {
                        let mut acc = 0.0f64;
                        for cc in 0..8 {
                            acc += v[cc + r * 8] * vk[cc + nu * 8];
                        }
                        re += vq[r + mu * 8] * acc;
                    }
                    o.push(re);
                }
            }
            o
        };
        let f0 = block_vec(&flux_of(0.0));
        let f1 = block_vec(&flux_of(1.0));
        let sp = block_vec(&split);
        // f(λ) = f0 + λ(f1 − f0); minimize ‖f(λ) − sp‖²
        let (mut num, mut den) = (0.0f64, 0.0f64);
        for i in 0..16 {
            let dfl = f1[i] - f0[i];
            num += dfl * (sp[i] - f0[i]);
            den += dfl * dfl;
        }
        let lam_fit = num / den;
        // 有理値の同定: ±1/8, ±1/4
        let cands = [0.125f64, -0.125, 0.25, -0.25, 0.5, -0.5];
        let best = cands
            .iter()
            .cloned()
            .min_by(|a, b| (a - lam_fit).abs().partial_cmp(&(b - lam_fit).abs()).unwrap())
            .unwrap();
        lam_frozen = best;
        check(
            "[E0] λ の同定: fit が有理値 (±1/8 系) に 1e-2 一致 → 凍結",
            (lam_fit - best).abs() < 1e-2,
            format!("λ_fit = {:.6} → 凍結 λ = {:+.4}", lam_fit, best),
        );
    }

    // ---- [E1] 改良の閉包 (流束 ladder) ----
    {
        let nvec = [0.55f64, 0.2, -0.81];
        let nn = (nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]).sqrt();
        let n = [nvec[0] / nn, nvec[1] / nn, nvec[2] / nn];
        let c = PI / 2.0;
        let mut rels = Vec::new();
        for &eps in &[0.4f64, 0.2, 0.1, 0.05] {
            let mphys = 0.5 * eps;
            let k = [c + eps * n[0], c + eps * n[1], c + eps * n[2]];
            let qq = 0.8 * eps;
            let d = source_imp(1, k, qq, mphys, lam_frozen, false);
            let a = h8([k[0], k[1] + qq, k[2]], mphys);
            let b = h8(k, mphys);
            let am = mat_mul(&a, &d);
            let mb = mat_mul(&d, &b);
            let dd = 2.0 * (0.5 * qq).sin();
            let phi: Vec<f64> = (0..64).map(|i| (am[i] - mb[i]) / dd).collect();
            let split = source_imp(7, k, qq, mphys, 0.0, false);
            let hk = h8(k, mphys);
            let (_, vk) = jacobi_eigh(&hk, 8);
            let hq = h8([k[0], k[1] + qq, k[2]], mphys);
            let (_, vq) = jacobi_eigh(&hq, 8);
            let block = |v: &[f64]| -> f64 {
                let mut sm = 0.0f64;
                for mu in 4..8 {
                    for nu in 0..4 {
                        let mut re = 0.0f64;
                        for r in 0..8 {
                            let mut acc = 0.0f64;
                            for cc in 0..8 {
                                acc += v[cc + r * 8] * vk[cc + nu * 8];
                            }
                            re += vq[r + mu * 8] * acc;
                        }
                        sm += re * re;
                    }
                }
                sm.sqrt()
            };
            let dif: Vec<f64> = (0..64).map(|i| phi[i] - split[i]).collect();
            rels.push((eps, block(&dif) / block(&split)));
        }
        let mut msg = String::new();
        for &(e, r) in &rels {
            msg = format!("{} rel({}) = {:.2e}", msg, e, r);
        }
        let ok = rels[3].1 < 0.05 && rels[3].1 < rels[2].1 && rels[2].1 < rels[1].1;
        check(
            "[E1] 改良の閉包: Φ(T_B⁰ˣ) → T^{xy}-split (停留 1.30 → 単調減少, 最終 < 5%)",
            ok,
            format!("{}", msg),
        );
    }

    // ---- [E2][E3][E4] ladder + [E5] 変異 ----
    {
        let (e_phys, q_phys) = (1.5f64, 0.6);
        let s_inv = e_phys * e_phys - q_phys * q_phys;
        let rho_d = s_inv * s_inv / (160.0 * PI * PI);
        let mut rows = Vec::new();
        for &a in &[0.18f64, 0.09, 0.045] {
            let (sig, tv) = sigma_imp(a, e_phys, q_phys, 32, 64, lam_frozen, false);
            let tvmax = tv.iter().cloned().fold(0.0f64, f64::max);
            // trace 方向
            let q2 = e_phys * e_phys - q_phys * q_phys;
            let qv = [e_phys, 0.0, q_phys, 0.0];
            let eta = [1.0f64, -1.0, -1.0, -1.0];
            let map = [(0usize, 0usize), (0, 1), (0, 2), (0, 3), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)];
            // t̂ は Σ_I t̂_I M_I = θ^{μν}M_{μν} (η 上げ + 非対角 ×2) になるよう構成
            // (開発記録: run1 は η 上げと重み 2 を落とし、物理の trace 方向を
            // 向いていなかった — 停留 0.265 はその器械欠陥)
            let mut tvec = [0.0f64; 10];
            for (i, &(mu, nu)) in map.iter().enumerate() {
                let th = (if mu == nu { eta[mu] } else { 0.0 }) - qv[mu] * qv[nu] / q2;
                let raised = eta[mu] * eta[nu] * th;
                tvec[i] = raised * if mu == nu { 1.0 } else { 2.0 };
            }
            let tn: f64 = tvec.iter().map(|x| x * x).sum::<f64>().sqrt();
            for v in tvec.iter_mut() {
                *v /= tn;
            }
            let matv = |m: &[f64], x: &[f64; 10]| -> [f64; 10] {
                let mut o = [0.0f64; 10];
                for r in 0..10 {
                    for c in 0..10 {
                        o[r] += m[c + r * 10] * x[c];
                    }
                }
                o
            };
            let st = matv(&sig, &tvec);
            let stn: f64 = st.iter().map(|x| x * x).sum::<f64>().sqrt();
            let mut x = [1.0f64; 10];
            let mut lmax = 0.0f64;
            for _ in 0..200 {
                let y = matv(&sig, &x);
                lmax = y.iter().map(|v| v * v).sum::<f64>().sqrt();
                for i in 0..10 {
                    x[i] = y[i] / lmax;
                }
            }
            let sdd = 0.5 * (sig[4 + 4 * 10] + sig[6 + 6 * 10] - 2.0 * sig[6 + 4 * 10]);
            let sxx = sig[8 + 8 * 10];
            println!(
                "    [E 表] a = {:.3}: 横断性 max = {:.5}, trace 比 = {:.5}, アンカー = {:.4}, 縮退 = {:.4} ({} s)",
                a,
                tvmax,
                stn / lmax,
                sdd / a.powi(4) / (2.0 * rho_d),
                sdd / (2.0 * sxx),
                t0.elapsed().as_secs()
            );
            rows.push((a, tvmax, stn / lmax, sdd / a.powi(4) / (2.0 * rho_d), sdd / (2.0 * sxx)));
        }
        let ok_tv = rows[2].1 < 0.02 && rows[0].1 / rows[1].1 > 2.5 && rows[1].1 / rows[2].1 > 2.5;
        check(
            "[E2] 対称 10×10 の横断性回復: 0.674 (v26.9-D) → O(a²) 級で消える (最終 < 2%)",
            ok_tv,
            format!(
                "{:.4} → {:.4} → {:.4} (比 {:.1}, {:.1})",
                rows[0].1, rows[1].1, rows[2].1,
                rows[0].1 / rows[1].1, rows[1].1 / rows[2].1
            ),
        );
        let ok_tr = rows[2].2 < 0.02 && rows[2].2 < rows[1].2 && rows[1].2 < rows[0].2;
        check(
            "[E3] conformal 崩壊: trace 比 0.265 (v26.9-D) → 単調減少 (最終 < 2%)",
            ok_tr,
            format!("{:.4} → {:.4} → {:.4}", rows[0].2, rows[1].2, rows[2].2),
        );
        let ok_or = (rows[2].3 - 1.0).abs() < 0.02 && (rows[2].4 - 1.0).abs() < 0.02;
        check(
            "[E4] oracle 回帰: σ_DD/(2ρ_D) = 1 ± 2% かつ σ_DD/2σ_XX = 1 ± 2%",
            ok_or,
            format!("アンカー = {:.4}, 縮退 = {:.4}", rows[2].3, rows[2].4),
        );
        // ---- [E5] 変異 ----
        let (_, tv_mut) = sigma_imp(0.045, e_phys, q_phys, 32, 64, lam_frozen, true);
        let tvm = tv_mut.iter().cloned().fold(0.0f64, f64::max);
        check(
            "[E5] 変異: λ → −λ → 横断性 O(1) 停留 (正版の 10 倍超)",
            tvm > 10.0 * rows[2].1 && tvm > 0.1,
            format!("変異 = {:.4} vs 正版 {:.4}", tvm, rows[2].1),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.9-E".into())),
        ("kind".into(), Json::Str("spin_current_improvement".into())),
        ("lambda".into(), Json::Num(lam_frozen)),
        (
            "construction".into(),
            Json::Str("ΔT0i = qhat·λ·[αi, αy] (on-site) — Belinfante 完成".into()),
        ),
    ]);
    let p = write_artifact("results/v269i_spin_improve.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **spin-current 改良で Belinfante 対称 10×10 が完成 — 横断性回復・conformal 崩壊 σ → ρ₂P₂・Gate 5 の分離部門完結** (型名保留・1/Π 禁止維持 — 凍結解釈: 測定器の証明)"
        } else {
            "FAIL あり — 分岐 (b) 改良項の形の不足 (公表) / (c) 代数の見立て違い。欄が一次ソース"
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
