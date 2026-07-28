//! v26.9-B v269w_ward4d — full 4D Ward: 4 行 × 8 列 × 2 周波数の一括機械検査
//!
//! 事前登録: spec §12.6 (Gate 5 = full 4D Ward [h₀₀, h₀ᵢ 込み])。energy 行
//! (v26.9-0) と momentum 行 (v26.9-A) で確立した器械を全 Ward 行に一般化する。
//!
//! 器械 (一様構成):
//!  密度 4 種 D_ν: V₀₀ = h(k+q/2ŷ) / V₀ᵢ = −(1/2)sin(2kᵢ + qδ_iy)·𝟙 (i = x,y,z)
//!  保存フラックス 4 種 (連続の式から自動生成・v26.9-A の実規約):
//!    **V_Fν(k;q) := [h(k+qŷ)·D_ν − D_ν·h(k)] / (2 sin(q/2))**
//!  (F0 = J_E [v26.9-0]・Fy = 保存 T_yy [v26.9-A]・Fx/Fz = 保存 T_yx/T_yz —
//!   局所性は Fourier 台ゲートで個別に機械証明)
//!  列 8 種 B: 密度 4 + BOND-A stress 4 (V_xx, V_yy [中点変調], V_zz,
//!   V_xz [v268z 認証済み 4 隅 point-split — d_y = 0 なので一般頂点公式が有効])
//!
//! Matsubara Ward (本ユニットの核): C_{AB}(iq₀,q) := Σ_direct N_B M_A/(Δ−iq₀)
//!  + Σ_crossed N'_A M'_B/(Δ'+iq₀) (M_A = ⟨μ,k+q|Â|ν,k⟩, N_B = (B̂(k+q;−q))_{νμ})
//!  に対し、対ごとの厳密関係 M_{J_A} = Δ·M_A/q̂ (連続の式) と分数分解
//!  iq₀/(Δ∓iq₀) = ∓1 ± Δ/(Δ∓iq₀) から
//!    **iq₀·C_{AB}(iq₀) − q̂·C_{J_A B}(iq₀) = −⟨[A(q), B(−q)]⟩**
//!  が k 点ごとに厳密。接触項は独立 2 実装で照合:
//!    (i) 対和: Σ_d N_B M_A − Σ_c N'_A M'_B
//!    (ii) 占有トレース: tr[P_occ(k)·(Â(k−qŷ;q)B̂(k;−q) − B̂(k+qŷ;−q)Â(k;q))]
//!
//! Ward の物理の所在 (run1 の教訓 — 重要): J_A := [h(k+qŷ)D − Dh(k)]/q̂ と
//! **定義**した時点で上の Matsubara 恒等式は任意の双線形 A に対する厳密再配列
//! であり、A の保存性を要しない (変異 A でも恒等式は成立する)。**Ward の物理的
//! 内容は (i) フラックスの局所性 [W0] — 保存が破れると J_A は q → 0 で 1/q̂
//! 発散し有界局所頂点でなくなる — と (ii) 接触項の構造 [W1] に宿る** (連続でも
//! Ward = 保存則 + カレントの局所性)。もう一つの定理 (run1 で発見・導出):
//! 接触項の対和 (= −⟨[A,B]⟩) と占有トレース (= +⟨[A,B]⟩ — run2 で符号の
//! 指紋 |差|/|和| = 2.0 から確定) は **k 点ごとには一致しない** — 差は
//! occ-occ 片 tr[P_occ(k)ÂP_occ(k∓q)B̂] で、k → k+q の相対ラベル替えにより
//! **BZ 和でのみ相殺**する。相殺を格子上で厳密にするため、BZ 和の検査は
//! **q を格子と可約に取る** (q = π/6 = 一様 12³ 格子の 2 刻み — シフトが
//! 格子自己同型になり occ-occ 相殺が機械精度で成立)。
//!
//! 検査 (凍結):
//!  [W0] フラックス局所性: V_Fν (ν = 0,x,y,z) の Fourier 台 |n| ≤ 3 (< 1e-12)
//!  [W1] 接触項の 2 実装一致 (BZ 和): 対和 = 占有トレース (4×8 組, 12³, 混合
//!       正規化 < 1e-10)
//!  [W2] **full 4D Ward**: 4 行 × 8 列 × q₀ ∈ {0.3, 0.9} = 64 恒等式が
//!       k 点ごと (3 サンプル) と 12³ BZ 積分の両方で < 1e-10 (混合正規化 —
//!       分母は |接触項| + 0.01·max|接触項| [小分母増幅の防止])
//!  [W3] BOND-A ↔ 保存フラックスの off-shell 距離: ‖V_Fy − V_yy^A‖_F の
//!       ε-ladder が単調減少 (指数は測って公表 — on-shell O(ε²) は v26.9-A 済み)
//!  [W4] 変異: (i) V₀x に s_x 依存重み ×1.02 (保存破れ) → **フラックスの局所性
//!       破れ = V_F(q → 0) 発散 (比 > 10)** / (ii) 接触項の符号反転 → Ward 破れ
//!
//! 事前登録分岐: (a) W0–W2 PASS → **局所カレント + 計算可能な接触項で full 4D
//!   Ward が格子上で厳密に閉じる** (Gate 5 の Ward 機構は確立 — 残り = BOND-A
//!   差の連続極限と spin-0/2 分離 [v26.9-C] の後に Gate 5 総括) / (b) W0 の特定
//!   行のみ破れ → その密度の連続の式の破れ (公表) / (c) W1 FAIL → 器械 (索引)。

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

// ---------------- 密度 4 種 (mutate は V₀x の taste 破り) ----------------

fn density(nu: usize, k: [f64; 3], q: f64, m: f64, mutate: bool) -> C8 {
    match nu {
        0 => c8_from_real(&h8([k[0], k[1] + 0.5 * q, k[2]], m)),
        _ => {
            let ax = nu - 1; // 1,2,3 → x,y,z
            let val = -0.5 * (2.0 * k[ax] + if ax == 1 { q } else { 0.0 }).sin();
            let mut v = vec![(0.0f64, 0.0f64); 64];
            for s in 0..8usize {
                let w = if mutate && nu == 1 && sbit(s, 0) == 1 { 1.02 } else { 1.0 };
                v[s + s * 8] = (val * w, 0.0);
            }
            v
        }
    }
}

/// 保存フラックス: V_Fν = [h(k+qŷ)D_ν − D_ν h(k)]/(2 sin(q/2))
fn flux(nu: usize, k: [f64; 3], q: f64, m: f64, mutate: bool) -> C8 {
    let a = c8_from_real(&h8([k[0], k[1] + q, k[2]], m));
    let b = c8_from_real(&h8(k, m));
    let d = density(nu, k, q, m, mutate);
    let am = c8_mul(&a, &d);
    let mb = c8_mul(&d, &b);
    let dd = 2.0 * (0.5 * q).sin();
    (0..64)
        .map(|i| ((am[i].0 - mb[i].0) / dd, (am[i].1 - mb[i].1) / dd))
        .collect()
}

// ---------------- BOND-A stress 4 種 ----------------

struct Term {
    eps: usize,
    d: [i32; 3],
    w: f64,
}

/// 一般頂点公式 (v268z 認証 — d_y = 0 の項に有効)
fn vertex8(terms: &[Term], k: [f64; 3], q: f64) -> Vec<f64> {
    let mut v = vec![0.0f64; 64];
    for t in terms {
        for s in 0..8usize {
            let s2 = s ^ t.eps;
            let mut ph1 = 0.0f64;
            let mut ph2 = 0.0f64;
            for ax in 0..3 {
                let ka = k[ax] + PI * sbit(s, ax) as f64;
                let ka2 = k[ax] + if ax == 1 { q } else { 0.0 } + PI * sbit(s2, ax) as f64;
                ph1 += ka * t.d[ax] as f64;
                ph2 += ka2 * t.d[ax] as f64;
            }
            let mid = 0.5 * q * t.d[1] as f64;
            v[s + s2 * 8] += t.w * ((ph1 + mid).cos() + (-ph2 + mid).cos());
        }
    }
    v
}

/// BOND-A 列: 0 = xx, 1 = yy (中点変調 — v26.9-0 導出), 2 = zz, 3 = xz (point-split)
fn bond_a(col: usize, k: [f64; 3], q: f64) -> C8 {
    match col {
        0 => c8_from_real(&vertex8(&[Term { eps: 0, d: [1, 0, 0], w: 0.5 }], k, q)),
        1 => {
            let mut v = vec![0.0f64; 64];
            for s in 0..8usize {
                let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * (k[1] + 0.5 * q).cos();
                v[(s ^ 1) + s * 8] = cy;
            }
            c8_from_real(&v)
        }
        2 => c8_from_real(&vertex8(&[Term { eps: 3, d: [0, 0, 1], w: 0.5 }], k, q)),
        _ => {
            let mut ts = Vec::new();
            for sg in [1i32, -1] {
                for rh in [1i32, -1] {
                    let c = (sg * rh) as f64 / 16.0;
                    ts.push(Term { eps: 0, d: [sg, 0, 2 * rh], w: -c });
                    ts.push(Term { eps: 3, d: [2 * rh, 0, sg], w: -c });
                }
            }
            c8_from_real(&vertex8(&ts, k, q))
        }
    }
}

/// 列 B (0..8): 0-3 = 密度 D_ν, 4-7 = BOND-A stress
fn column(b: usize, k: [f64; 3], q: f64, m: f64) -> C8 {
    if b < 4 {
        density(b, k, q, m, false)
    } else {
        bond_a(b - 4, k, q)
    }
}

// ---------------- 対和と接触項 ----------------

/// 固有分解 (占有 = 0..4, 非占有 = 4..8; jacobi_eigh は昇順)
fn eig(k: [f64; 3], m: f64) -> (Vec<f64>, Vec<f64>) {
    jacobi_eigh(&h8(k, m), 8)
}

/// ⟨μ|V|ν⟩ (μ: 左固有ベクトル列 vq, ν: 右固有ベクトル列 vk) — 複素
fn mel(v: &C8, vq: &[f64], mu: usize, vk: &[f64], nu: usize) -> (f64, f64) {
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
    (re, im)
}

/// 1 k 点の Ward 三点セット: 行 nu, 列 b, 周波数 q₀ に対する
/// (iq₀·C_AB − q̂·C_JB, 接触項 [対和], 接触項 [トレース]) — 全て複素
#[allow(clippy::too_many_arguments)]
fn ward_at_k(
    nu: usize,
    b: usize,
    k: [f64; 3],
    q: f64,
    q0: f64,
    m: f64,
    mutate: bool,
) -> ((f64, f64), (f64, f64), (f64, f64)) {
    let qhat = 2.0 * (0.5 * q).sin();
    let (wk, vk) = eig(k, m);
    let kq = [k[0], k[1] + q, k[2]];
    let (wq, vq) = eig(kq, m);
    let kmq = [k[0], k[1] - q, k[2]];
    let (wm, vm) = eig(kmq, m);
    // 頂点
    let a_v = density(nu, k, q, m, mutate); // Â(k;q): k → k+q
    let j_v = flux(nu, k, q, m, mutate); // J_A(k;q)
    let b_rev = column(b, kq, -q, m); // B̂(k+q;−q): k+q → k
    let a_cr = density(nu, kmq, q, m, mutate); // Â(k−q;q): k−q → k
    let j_cr = flux(nu, kmq, q, m, mutate);
    let b_cr = column(b, k, -q, m); // B̂(k;−q): k → k−q
    let mut lhs = (0.0f64, 0.0f64);
    let mut cont_pair = (0.0f64, 0.0f64);
    // 直接セクター: |ν,k⟩ → |μ,k+q⟩
    for mu in 4..8 {
        for nuo in 0..4 {
            let ma = mel(&a_v, &vq, mu, &vk, nuo);
            let mj = mel(&j_v, &vq, mu, &vk, nuo);
            let nb = mel(&b_rev, &vk, nuo, &vq, mu); // (B̂)_{ν μ}
            let de = wq[mu] - wk[nuo];
            // 1/(Δ − iq₀) = (Δ + iq₀)/(Δ²+q₀²)
            let den = de * de + q0 * q0;
            let (gre, gim) = (de / den, q0 / den);
            // C_AB 直接片: N_B·M_A·g / C_JB: N_B·M_J·g
            let prod = |x: (f64, f64), y: (f64, f64)| (x.0 * y.0 - x.1 * y.1, x.0 * y.1 + x.1 * y.0);
            let am = prod(nb, ma);
            let jm = prod(nb, mj);
            let cab = (am.0 * gre - am.1 * gim, am.0 * gim + am.1 * gre);
            let cjb = (jm.0 * gre - jm.1 * gim, jm.0 * gim + jm.1 * gre);
            // iq₀·C − q̂·C_J
            lhs.0 += -q0 * cab.1 - qhat * cjb.0;
            lhs.1 += q0 * cab.0 - qhat * cjb.1;
            cont_pair.0 += am.0;
            cont_pair.1 += am.1;
        }
    }
    // 交差セクター: |ν,k⟩ → |μ,k−q⟩ (B が励起, A が脱励起)
    for mu in 4..8 {
        for nuo in 0..4 {
            let mb = mel(&b_cr, &vm, mu, &vk, nuo);
            let na = mel(&a_cr, &vk, nuo, &vm, mu);
            let nj = mel(&j_cr, &vk, nuo, &vm, mu);
            let de = wm[mu] - wk[nuo];
            // 1/(Δ + iq₀) = (Δ − iq₀)/(Δ²+q₀²)
            let den = de * de + q0 * q0;
            let (gre, gim) = (de / den, -q0 / den);
            let prod = |x: (f64, f64), y: (f64, f64)| (x.0 * y.0 - x.1 * y.1, x.0 * y.1 + x.1 * y.0);
            let am = prod(na, mb);
            let jm = prod(nj, mb);
            let cab = (am.0 * gre - am.1 * gim, am.0 * gim + am.1 * gre);
            let cjb = (jm.0 * gre - jm.1 * gim, jm.0 * gim + jm.1 * gre);
            lhs.0 += -q0 * cab.1 - qhat * cjb.0;
            lhs.1 += q0 * cab.0 - qhat * cjb.1;
            cont_pair.0 -= am.0;
            cont_pair.1 -= am.1;
        }
    }
    // 接触項 (占有トレース): tr[P_occ(k)(Â(k−q;q)B̂(k;−q) − B̂(k+q;−q)Â(k;q))]
    let mut cont_tr = (0.0f64, 0.0f64);
    {
        let t1 = c8_mul(&a_cr, &b_cr);
        let t2 = c8_mul(&b_rev, &a_v);
        for nuo in 0..4 {
            for r in 0..8 {
                for c in 0..8 {
                    let p = vk[r + nuo * 8] * vk[c + nuo * 8];
                    cont_tr.0 += p * (t1[c + r * 8].0 - t2[c + r * 8].0);
                    cont_tr.1 += p * (t1[c + r * 8].1 - t2[c + r * 8].1);
                }
            }
        }
    }
    (lhs, cont_pair, cont_tr)
}

fn main() {
    self_test();
    println!("=== v26.9-B v269w_ward4d — full 4D Ward: 4 行 × 8 列 × 2 周波数 ===\n");
    println!("iq₀·C_AB(iq₀) − q̂·C_{{J_A B}}(iq₀) = −⟨[A(q),B(−q)]⟩ を保存カレント系で");
    println!("一括機械検査。接触項は対和/占有トレースの独立 2 実装。Gate 5 の Ward 機構。\n");
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
    let q = 0.4f64;
    let names_r = ["ν=0", "ν=x", "ν=y", "ν=z"];
    let names_c = ["T00", "T0x", "T0y", "T0z", "Axx", "Ayy", "Azz", "Axz"];

    // ---- [W0] フラックス局所性 ----
    {
        let mut worst = 0.0f64;
        for nu in 0..4 {
            for &qv in &[0.1f64, 0.7] {
                for ax in 0..3usize {
                    let ngrid = 16usize;
                    let base = [0.37f64, -0.81, 1.13];
                    let vsl: Vec<C8> = (0..ngrid)
                        .map(|j| {
                            let mut k = base;
                            k[ax] = 2.0 * PI * j as f64 / ngrid as f64;
                            flux(nu, k, qv, m0, false)
                        })
                        .collect();
                    for r in 0..8 {
                        for c in 0..8 {
                            let mut hi = 0.0f64;
                            let mut all = 0.0f64;
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
                                all += p;
                                if nsig.abs() > 3 {
                                    hi += p;
                                }
                            }
                            if all > 1e-20 {
                                worst = worst.max((hi / all).sqrt());
                            }
                        }
                    }
                }
            }
        }
        check(
            "[W0] 保存フラックス V_Fν (ν = 0,x,y,z) の Fourier 台 |n| ≤ 3",
            worst < 1e-12,
            format!("max 台残差 = {:.1e}", worst),
        );
    }

    // ---- [W2a] Ward (サンプル k, 混合正規化) ----
    {
        let samples = [[0.5f64, 1.2, -0.8], [1.7, -0.4, 0.9], [0.23, 2.4, 1.55]];
        let mut worst_w = 0.0f64;
        for &k in &samples {
            // 2 パス: cmax を先に求め、小分母増幅を防ぐ
            let mut recs = Vec::new();
            let mut cmax = 0.0f64;
            for nu in 0..4 {
                for b in 0..8 {
                    for &q0 in &[0.3f64, 0.9] {
                        let (lhs, cp, _) = ward_at_k(nu, b, k, q, q0, m0, false);
                        cmax = cmax.max((cp.0 * cp.0 + cp.1 * cp.1).sqrt());
                        recs.push((lhs, cp));
                    }
                }
            }
            for (lhs, cp) in recs {
                let scale = (cp.0 * cp.0 + cp.1 * cp.1).sqrt() + 0.01 * cmax;
                worst_w = worst_w
                    .max(((lhs.0 + cp.0).powi(2) + (lhs.1 + cp.1).powi(2)).sqrt() / scale);
            }
        }
        check(
            "[W2a] full 4D Ward (k 点ごと): iq₀C − q̂C_J = −⟨[A,B]⟩ (64 恒等式 × 3 k)",
            worst_w < 1e-10,
            format!("max 相対残差 = {:.1e}", worst_w),
        );
    }

    // ---- [W1] 接触項 2 実装 (BZ 和) & [W2b] Ward (BZ 積分) ----
    {
        let ngrid = 12usize;
        let q = PI / 6.0; // 格子可約 (2 刻み) — occ-occ 相殺が厳密になる
        let combos: Vec<(usize, usize, f64)> = {
            let mut v = Vec::new();
            for nu in 0..4 {
                for b in 0..8 {
                    for &q0 in &[0.3f64, 0.9] {
                        v.push((nu, b, q0));
                    }
                }
            }
            v
        };
        let chunk = combos.len().div_ceil(nthreads);
        let mut outs: Vec<Option<Vec<(f64, f64, f64, f64, f64, f64)>>> = Vec::new();
        outs.resize_with(nthreads, || None);
        std::thread::scope(|sc| {
            for (t, slot) in outs.iter_mut().enumerate() {
                let combos = &combos;
                sc.spawn(move || {
                    let mut res = Vec::new();
                    for &(nu, b, q0) in combos.iter().skip(t * chunk).take(chunk) {
                        let (mut l, mut c, mut ctr) =
                            ((0.0f64, 0.0f64), (0.0f64, 0.0f64), (0.0f64, 0.0f64));
                        for jx in 0..ngrid {
                            for jy in 0..ngrid {
                                for jz in 0..ngrid {
                                    let k = [
                                        PI * (jx as f64 + 0.5) / ngrid as f64,
                                        PI * (jy as f64 + 0.5) / ngrid as f64,
                                        PI * (jz as f64 + 0.5) / ngrid as f64,
                                    ];
                                    let (lh, cp, ct) = ward_at_k(nu, b, k, q, q0, m0, false);
                                    l.0 += lh.0;
                                    l.1 += lh.1;
                                    c.0 += cp.0;
                                    c.1 += cp.1;
                                    ctr.0 += ct.0;
                                    ctr.1 += ct.1;
                                }
                            }
                        }
                        res.push((l.0, l.1, c.0, c.1, ctr.0, ctr.1));
                    }
                    *slot = Some(res);
                });
            }
        });
        let mut all: Vec<(f64, f64, f64, f64, f64, f64)> = Vec::new();
        for o in outs.into_iter() {
            all.extend(o.unwrap());
        }
        let cmax = all
            .iter()
            .map(|r| (r.2 * r.2 + r.3 * r.3).sqrt())
            .fold(0.0f64, f64::max);
        let mut worst_c = 0.0f64;
        let mut worst_w = 0.0f64;
        for r in &all {
            let scale = (r.2 * r.2 + r.3 * r.3).sqrt() + 0.01 * cmax;
            // 対和 = −⟨[A,B]⟩, トレース = +⟨[A,B]⟩ ⇒ 照合は cp + ct = 0
            worst_c = worst_c.max(((r.2 + r.4).powi(2) + (r.3 + r.5).powi(2)).sqrt() / scale);
            worst_w = worst_w.max(((r.0 + r.2).powi(2) + (r.1 + r.3).powi(2)).sqrt() / scale);
        }
        check(
            "[W1] 接触項 2 実装 (BZ 和のみで一致する定理): 対和 = −占有トレース (12³, q 可約)",
            worst_c < 1e-10,
            format!("max 相対差 = {:.1e}", worst_c),
        );
        check(
            "[W2b] full 4D Ward (12³ BZ 積分): 64 恒等式すべて",
            worst_w < 1e-10,
            format!("max 相対残差 = {:.1e} ({} s)", worst_w, t0.elapsed().as_secs()),
        );
    }

    // ---- [W3] BOND-A ↔ 保存フラックスの off-shell 距離 ----
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
            let vf = flux(2, k, qq, mphys, false);
            let va = bond_a(1, k, qq);
            let (mut d2, mut n2) = (0.0f64, 0.0f64);
            for i in 0..64 {
                d2 += (vf[i].0 - va[i].0).powi(2) + (vf[i].1 - va[i].1).powi(2);
                n2 += va[i].0.powi(2) + va[i].1.powi(2);
            }
            rels.push((eps, (d2 / n2).sqrt()));
        }
        let mut msg = String::new();
        for &(e, r) in &rels {
            msg = format!("{} rel({}) = {:.2e}", msg, e, r);
        }
        let expo = (rels[2].1 / rels[3].1).log2();
        let ok = rels[3].1 < rels[2].1 && rels[2].1 < rels[1].1 && rels[1].1 < rels[0].1;
        check(
            "[W3] ‖V_Fy − V_yy^A‖_F (off-shell 全行列) の ε-ladder 単調減少 (指数は公表)",
            ok,
            format!("{} — 実効指数 ~ε^{:.1}", msg, expo),
        );
    }

    // ---- [W4] 変異 ----
    {
        let k = [0.5f64, 1.2, -0.8];
        // (i) V₀x の taste 破り (s_x 依存重み) → 保存破れ = フラックス局所性破れ
        //     = V_F(q → 0) の 1/q̂ 発散 (恒等式自体は任意の A で成立 — 物理は局所性)
        let norm = |v: &C8| -> f64 {
            v.iter()
                .map(|x| (x.0 * x.0 + x.1 * x.1).sqrt())
                .fold(0.0f64, f64::max)
        };
        let r_mut = norm(&flux(1, k, 1e-3, m0, true)) / norm(&flux(1, k, 0.1, m0, true));
        let r_good = norm(&flux(1, k, 1e-3, m0, false)) / norm(&flux(1, k, 0.1, m0, false));
        // (ii) 接触項の符号反転 (正版でも lhs − cp ≠ 0)
        let (lhs2, cp2, _) = ward_at_k(0, 4, k, q, 0.3, m0, false);
        let r_sign = ((lhs2.0 - cp2.0).powi(2) + (lhs2.1 - cp2.1).powi(2)).sqrt()
            / (cp2.0 * cp2.0 + cp2.1 * cp2.1).sqrt().max(1e-6);
        check(
            "[W4] 変異: (i) V₀x taste 破り → V_Fx(q→0) 発散 (比 > 10; 正版 < 2) / (ii) 接触項符号反転 → 破れ",
            r_mut > 10.0 && r_good < 2.0 && r_sign > 0.1,
            format!(
                "発散比 = {:.1} (正版 {:.2}), 符号反転残差 = {:.2}",
                r_mut, r_good, r_sign
            ),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.9-B".into())),
        ("kind".into(), Json::Str("full_4d_ward_conserved".into())),
        (
            "identity".into(),
            Json::Str("iq0·C_AB − qhat·C_JB = −⟨[A(q),B(−q)]⟩ (4 行 × 8 列 × 2 周波数)".into()),
        ),
        ("rows".into(), Json::Arr(names_r.iter().map(|s| Json::Str(s.to_string())).collect())),
        ("cols".into(), Json::Arr(names_c.iter().map(|s| Json::Str(s.to_string())).collect())),
    ]);
    let p = write_artifact("results/v269w_ward4d.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **保存カレント系での full 4D Ward (4 行 × 8 列 × 2 周波数 = 64 恒等式) が格子上で厳密に閉じる — Gate 5 の Ward 機構確立** (残り = spin-0/2 の 4D 分離と BOND-A 差の連続極限 [v26.9-C] → Gate 5 総括)"
        } else {
            "FAIL あり — 分岐 (b) 特定行の連続の式破れ (公表) / (c) 器械 (索引)。欄が一次ソース"
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
