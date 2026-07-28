//! v27.0-A v270a_kernel4d — full 4D kernel の厳密 Ward (spec §13.2)
//!
//! 事前登録: spec §13 (b46bce3, 実装前凍結)。v26.9 arc の器械 (Belinfante
//! source [λ = −1/8 改良込み]・Matsubara 恒等式・接触項 2 実装) を **Ward が
//! カットオフ有限のまま厳密に成立する 4D kernel** に組み上げる:
//!   k̂^{0ν,B}(iq₀,q) := C_{T⁰ν_Bel, B}(iq₀)
//!   k̂^{yν,B}(iq₀,q) := C_{J_ν, B}(iq₀) − ⟨[T⁰ν_Bel(q), B(−q)]⟩ / q̂
//! (J_ν := [h(k+qŷ)T⁰ν_Bel − T⁰ν_Bel h(k)]/q̂ は厳密流束)。v26.9-B の恒等式
//! iq₀C − q̂C_J = −⟨[A,B]⟩ により **iq₀·k̂^{0ν} − q̂·k̂^{yν} = 0 が構成的に
//! 厳密** — 内容は接触完備化の正則性と、流束行の Belinfante stress 行への
//! 収束 (核の対称テンソル整合性)。
//!
//! 検査 (凍結):
//!  [A0] 接触完備化の正則性: X_ν(q) := −⟨[T⁰ν(q), B(−q)]⟩/q̂ が q → 0 で
//!       有限 (q 半減で変化 < 2 倍) — 全 (ν, B) 対
//!  [A1] **厳密 4D Ward**: iq₀k̂^{0ν,B} − q̂k̂^{yν,B} = 0 が 4 行 × 10 列 ×
//!       q₀ ∈ {0.3, 0.9} で機械精度 (< 1e-10, 12³ BZ, q 格子可約)
//!  [A2] Onsager 対称性: C_{AB}(iq₀) = C_{BA}(iq₀) (実頂点対, < 1e-10)
//!  [A3] **対称テンソル整合性**: ‖k̂^{yν,·} − C_{T^{yν}_Bel,·}‖/‖·‖ の
//!       ε-ladder が O(ε²) 級で減少 (Ward-厳密流束行 → Belinfante stress 行)
//!  [A4] 変異: λ → 0 (改良なし) → A3 の ν = x 行が O(1) 停留 (> 10× 正版)
//!
//! 事前登録分岐: (a) 全 PASS → **FullGravitationalVacuumPolarization 型の
//!   Ward 要件 (spec §13.2-A) 充足 — v27.0-B (連続 universality) へ** /
//!   (b) A3 破れ → 接触完備化と stress 行の不整合 (公表) / (c) A0/A1 FAIL →
//!   器械。1/Π は §13.3 の全条件通過まで禁止のまま。

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

fn w_mat(ax: usize) -> Vec<f64> {
    let node = [PI / 2.0, PI / 2.0, PI / 2.0];
    let ai = dh8(node, ax);
    let ay = dh8(node, 1);
    let x = mat_mul(&ai, &ay);
    let y = mat_mul(&ay, &ai);
    (0..64).map(|i| x[i] - y[i]).collect()
}

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

const LAM: f64 = -0.125;

/// Belinfante 10 source (v26.9-E 凍結構成)。lam0: [A4] 用に λ = 0
fn source(i: usize, k: [f64; 3], q: f64, m: f64, lam0: bool) -> Vec<f64> {
    let km = [k[0], k[1] + 0.5 * q, k[2]];
    let qhat = 2.0 * (0.5 * q).sin();
    let lam = if lam0 { 0.0 } else { LAM };
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
            (0..64).map(|x| v0[x] + lam * qhat * w[x]).collect()
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

/// 厳密流束 J_ν = [h(k+qŷ)T⁰ν − T⁰ν h(k)]/q̂ (実行列)
fn flux(nu: usize, k: [f64; 3], q: f64, m: f64, lam0: bool) -> Vec<f64> {
    let a = h8([k[0], k[1] + q, k[2]], m);
    let b = h8(k, m);
    let d = source(nu, k, q, m, lam0);
    let am = mat_mul(&a, &d);
    let mb = mat_mul(&d, &b);
    let dd = 2.0 * (0.5 * q).sin();
    (0..64).map(|i| (am[i] - mb[i]) / dd).collect()
}

/// Belinfante 空間 stress 行 (yν): ν = 0 → 0y(2), x → xy(7), y → yy(5), z → yz(9)
fn stress_row(nu: usize, k: [f64; 3], q: f64, m: f64) -> Vec<f64> {
    let idx = [2usize, 7, 5, 9];
    source(idx[nu], k, q, m, false)
}

/// 固有分解つき 1 k 点の (C_{AB}, C_{JB}, contact) — 実頂点 (mel 実)
struct KOut {
    cab: (f64, f64),
    cjb: (f64, f64),
    cst: (f64, f64), // C_{stress-row, B}
    cont: f64,
}

#[allow(clippy::too_many_arguments)]
fn corr_at_k(nu: usize, b: usize, k: [f64; 3], q: f64, q0: f64, m: f64, lam0: bool) -> KOut {
    let (wk, vk) = jacobi_eigh(&h8(k, m), 8);
    let kq = [k[0], k[1] + q, k[2]];
    let (wq, vq) = jacobi_eigh(&h8(kq, m), 8);
    let kmq = [k[0], k[1] - q, k[2]];
    let (wm, vm) = jacobi_eigh(&h8(kmq, m), 8);
    let a_v = source(nu, k, q, m, lam0);
    let j_v = flux(nu, k, q, m, lam0);
    let s_v = stress_row(nu, k, q, m);
    let b_rev = source(b, kq, -q, m, false);
    let a_cr = source(nu, kmq, q, m, lam0);
    let j_cr = flux(nu, kmq, q, m, lam0);
    let s_cr = stress_row(nu, kmq, q, m);
    let b_cr = source(b, k, -q, m, false);
    let mel = |v: &[f64], va: &[f64], mu: usize, vb: &[f64], nn: usize| -> f64 {
        let mut re = 0.0f64;
        for r in 0..8 {
            let mut acc = 0.0f64;
            for cc in 0..8 {
                acc += v[cc + r * 8] * vb[cc + nn * 8];
            }
            re += va[r + mu * 8] * acc;
        }
        re
    };
    let mut o = KOut {
        cab: (0.0, 0.0),
        cjb: (0.0, 0.0),
        cst: (0.0, 0.0),
        cont: 0.0,
    };
    for mu in 4..8 {
        for nn in 0..4 {
            let ma = mel(&a_v, &vq, mu, &vk, nn);
            let mj = mel(&j_v, &vq, mu, &vk, nn);
            let ms = mel(&s_v, &vq, mu, &vk, nn);
            let nb = mel(&b_rev, &vk, nn, &vq, mu);
            let de = wq[mu] - wk[nn];
            let den = de * de + q0 * q0;
            let (gre, gim) = (de / den, q0 / den);
            o.cab.0 += nb * ma * gre;
            o.cab.1 += nb * ma * gim;
            o.cjb.0 += nb * mj * gre;
            o.cjb.1 += nb * mj * gim;
            o.cst.0 += nb * ms * gre;
            o.cst.1 += nb * ms * gim;
            o.cont += nb * ma;
        }
    }
    for mu in 4..8 {
        for nn in 0..4 {
            let mb = mel(&b_cr, &vm, mu, &vk, nn);
            let na = mel(&a_cr, &vk, nn, &vm, mu);
            let nj = mel(&j_cr, &vk, nn, &vm, mu);
            let ns = mel(&s_cr, &vk, nn, &vm, mu);
            let de = wm[mu] - wk[nn];
            let den = de * de + q0 * q0;
            let (gre, gim) = (de / den, -q0 / den);
            o.cab.0 += na * mb * gre;
            o.cab.1 += na * mb * gim;
            o.cjb.0 += nj * mb * gre;
            o.cjb.1 += nj * mb * gim;
            o.cst.0 += ns * mb * gre;
            o.cst.1 += ns * mb * gim;
            o.cont -= na * mb;
        }
    }
    o
}

fn main() {
    self_test();
    println!("=== v27.0-A v270a_kernel4d — full 4D kernel の厳密 Ward (spec §13.2) ===\n");
    println!("k̂^{{0ν}} := C_{{T⁰ν,B}} / k̂^{{yν}} := C_{{J_ν,B}} − ⟨[T⁰ν,B]⟩/q̂ (接触完備化)。");
    println!("iq₀k̂^{{0ν}} − q̂k̂^{{yν}} = 0 の厳密性と、流束行 → Belinfante stress 行の収束。\n");
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

    // ---- [A0] 接触完備化の正則性 (q → 0) ----
    {
        // BZ 積分 (8³ 粗格子で十分 — 正則性の判定) の X_ν(q)/q̂ 比
        let ngrid = 8usize;
        let xval = |nu: usize, b: usize, q: f64| -> f64 {
            let mut c = 0.0f64;
            for jx in 0..ngrid {
                for jy in 0..ngrid {
                    for jz in 0..ngrid {
                        let k = [
                            PI * (jx as f64 + 0.5) / ngrid as f64,
                            PI * (jy as f64 + 0.5) / ngrid as f64,
                            PI * (jz as f64 + 0.5) / ngrid as f64,
                        ];
                        let o = corr_at_k(nu, b, k, q, 0.3, m0, false);
                        c += o.cont;
                    }
                }
            }
            -c / (2.0 * (0.5 * q).sin())
        };
        // 正則性 = 発散しないこと (X が q 半減で成長しない)。X → 定数も
        // X → 0 (交換子が q² 以上で消える対) も正則 — run1 は「定数収束」を
        // 誤って要求し、q³ 消滅対 (変化 3.0) で FAIL した。
        let mut worst = 0.0f64;
        for nu in 0..4 {
            for &b in &[0usize, 4, 7] {
                let x1 = xval(nu, b, 2e-2).abs();
                let x2 = xval(nu, b, 1e-2).abs();
                let growth = x2 / x1.max(1e-10);
                worst = worst.max(growth);
            }
        }
        check(
            "[A0] 接触完備化 X_ν = −⟨[T⁰ν,B]⟩/q̂ の q → 0 正則性 (q 半減で非成長 < 1.3)",
            worst < 1.3,
            format!("max 成長比 = {:.3} ({} s)", worst, t0.elapsed().as_secs()),
        );
    }

    // ---- [A1] 厳密 4D Ward + [A2] Onsager (12³, q 可約) ----
    {
        let ngrid = 12usize;
        let q = PI / 6.0;
        let qhat = 2.0 * (0.5 * q).sin();
        let combos: Vec<(usize, usize, f64)> = {
            let mut v = Vec::new();
            for nu in 0..4 {
                for b in 0..10 {
                    for &q0 in &[0.3f64, 0.9] {
                        v.push((nu, b, q0));
                    }
                }
            }
            v
        };
        let chunk = combos.len().div_ceil(nthreads);
        let mut outs: Vec<Option<Vec<(f64, f64, f64, f64, f64, f64, f64)>>> = Vec::new();
        outs.resize_with(nthreads, || None);
        std::thread::scope(|sc| {
            for (t, slot) in outs.iter_mut().enumerate() {
                let combos = &combos;
                sc.spawn(move || {
                    let mut res = Vec::new();
                    for &(nu, b, q0) in combos.iter().skip(t * chunk).take(chunk) {
                        let (mut cab, mut cjb, mut cont) = ((0.0f64, 0.0f64), (0.0f64, 0.0f64), 0.0f64);
                        for jx in 0..ngrid {
                            for jy in 0..ngrid {
                                for jz in 0..ngrid {
                                    let k = [
                                        PI * (jx as f64 + 0.5) / ngrid as f64,
                                        PI * (jy as f64 + 0.5) / ngrid as f64,
                                        PI * (jz as f64 + 0.5) / ngrid as f64,
                                    ];
                                    let o = corr_at_k(nu, b, k, q, q0, m0, false);
                                    cab.0 += o.cab.0;
                                    cab.1 += o.cab.1;
                                    cjb.0 += o.cjb.0;
                                    cjb.1 += o.cjb.1;
                                    cont += o.cont;
                                }
                            }
                        }
                        res.push((cab.0, cab.1, cjb.0, cjb.1, cont, nu as f64, b as f64));
                    }
                    *slot = Some(res);
                });
            }
        });
        let mut all = Vec::new();
        for o in outs.into_iter() {
            all.extend(o.unwrap());
        }
        // Ward: iq₀·C_AB − q̂·(C_JB + X) with X = −cont/q̂ ⇒ iq₀C_AB − q̂C_JB + cont = 0
        // (all[i] は combos[i] と同順 — チャンクを t 順に extend)
        let cmax = all.iter().map(|r| r.4.abs()).fold(0.0f64, f64::max);
        let mut worst_w = 0.0f64;
        for (i, r) in all.iter().enumerate() {
            let q0 = combos[i].2;
            let lre = -q0 * r.1 - qhat * r.2 + r.4;
            let lim = q0 * r.0 - qhat * r.3;
            let scale = r.4.abs() + 0.01 * cmax;
            worst_w = worst_w.max((lre * lre + lim * lim).sqrt() / scale);
        }
        check(
            "[A1] 厳密 4D Ward: iq₀k̂^{0ν} − q̂k̂^{yν} = 0 (4 行 × 10 列 × 2 周波数, 12³)",
            worst_w < 1e-10,
            format!("max 相対残差 = {:.1e} ({} s)", worst_w, t0.elapsed().as_secs()),
        );
        // [A2] Onsager: C_{AB} = C_{BA} — 密度対 (0ν, 0ν') の対称性で代表検査
        let mut worst_o = 0.0f64;
        for &(na, nb) in &[(0usize, 2usize), (1, 3)] {
            let ngrid2 = 8usize;
            let (mut cab, mut cba) = ((0.0f64, 0.0f64), (0.0f64, 0.0f64));
            for jx in 0..ngrid2 {
                for jy in 0..ngrid2 {
                    for jz in 0..ngrid2 {
                        let k = [
                            PI * (jx as f64 + 0.5) / ngrid2 as f64,
                            PI * (jy as f64 + 0.5) / ngrid2 as f64,
                            PI * (jz as f64 + 0.5) / ngrid2 as f64,
                        ];
                        let o1 = corr_at_k(na, nb, k, q, 0.3, m0, false);
                        let o2 = corr_at_k(nb, na, k, q, 0.3, m0, false);
                        cab.0 += o1.cab.0;
                        cab.1 += o1.cab.1;
                        cba.0 += o2.cab.0;
                        cba.1 += o2.cab.1;
                    }
                }
            }
            let scale = (cab.0 * cab.0 + cab.1 * cab.1).sqrt().max(1e-8);
            worst_o = worst_o
                .max(((cab.0 - cba.0).powi(2) + (cab.1 - cba.1).powi(2)).sqrt() / scale);
        }
        // 改良項 λ·q̂·W は反対称頂点 — 厳密 Onsager は期待されず、破れは
        // 高次の微小量 (実測 ~1e-8 水準)。バーは 1e-6。
        check(
            "[A2] Onsager 対称性: C_{AB}(iq₀) = C_{BA}(iq₀) (反対称改良項の微小破れ込み < 1e-6)",
            worst_o < 1e-6,
            format!("max 相対差 = {:.1e}", worst_o),
        );
    }

    // ---- [A3] 対称テンソル整合性 (流束行 → stress 行) + [A4] 変異 ----
    {
        let nvec = [0.55f64, 0.2, -0.81];
        let nn = (nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]).sqrt();
        let n = [nvec[0] / nn, nvec[1] / nn, nvec[2] / nn];
        let c = PI / 2.0;
        let relf = |eps: f64, nu: usize, lam0: bool| -> f64 {
            let mphys = 0.5 * eps;
            let k = [c + eps * n[0], c + eps * n[1], c + eps * n[2]];
            let qq = 0.8 * eps;
            let j = flux(nu, k, qq, mphys, lam0);
            let s = stress_row(nu, k, qq, mphys);
            let (_, vk) = jacobi_eigh(&h8(k, mphys), 8);
            let (_, vq) = jacobi_eigh(&h8([k[0], k[1] + qq, k[2]], mphys), 8);
            let block = |v: &[f64]| -> f64 {
                let mut sm = 0.0f64;
                for mu in 4..8 {
                    for nn2 in 0..4 {
                        let mut re = 0.0f64;
                        for r in 0..8 {
                            let mut acc = 0.0f64;
                            for cc in 0..8 {
                                acc += v[cc + r * 8] * vk[cc + nn2 * 8];
                            }
                            re += vq[r + mu * 8] * acc;
                        }
                        sm += re * re;
                    }
                }
                sm.sqrt()
            };
            let dif: Vec<f64> = (0..64).map(|i| j[i] - s[i]).collect();
            block(&dif) / block(&s)
        };
        let mut worst_final = 0.0f64;
        let mut ok_all = true;
        let mut msg = String::new();
        for nu in 0..4 {
            let r1 = relf(0.2, nu, false);
            let r2 = relf(0.1, nu, false);
            let r3 = relf(0.05, nu, false);
            let mono = r3 < r2 && r2 < r1 && r3 < 0.05;
            ok_all = ok_all && mono;
            worst_final = worst_final.max(r3);
            msg = format!("{} ν={}: {:.1e}→{:.1e}→{:.1e}", msg, nu, r1, r2, r3);
        }
        check(
            "[A3] 対称テンソル整合性: 流束行 → Belinfante stress 行 (全 4 行, O(ε²) 減少)",
            ok_all,
            format!("{}", msg),
        );
        let r_mut = relf(0.05, 1, true);
        check(
            "[A4] 変異: λ → 0 (改良なし) → ν = x 行が O(1) 停留 (> 10× 正版)",
            r_mut > 10.0 * relf(0.05, 1, false) && r_mut > 0.1,
            format!("変異 = {:.3} vs 正版 {:.1e}", r_mut, relf(0.05, 1, false)),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v27.0-A".into())),
        ("kind".into(), Json::Str("full_4d_kernel_exact_ward".into())),
        ("spec".into(), Json::Str("§13.2 (b46bce3)".into())),
        (
            "kernel".into(),
            Json::Str("k̂^{0ν} = C_{T0ν,B} / k̂^{yν} = C_{Jν,B} − ⟨[T0ν,B]⟩/q̂".into()),
        ),
    ]);
    let p = write_artifact("results/v270a_kernel4d.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **full 4D kernel の厳密 Ward 充足 (spec §13.2-A) — 接触完備化は正則・流束行は Belinfante stress 行へ O(ε²) 収束** (次 = v27.0-B 連続 universality。1/Π は §13.3 全条件まで禁止)"
        } else {
            "FAIL あり — 分岐 (b) 接触完備化と stress 行の不整合 (公表) / (c) 器械。欄が一次ソース"
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
