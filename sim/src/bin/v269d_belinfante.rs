//! v26.9-D v269d_belinfante — Belinfante 対称化の素朴平均は失敗する (負の結果)
//!
//! 事前登録: v26.9-C の Gate 5 総括で登録した「Belinfante 対称 10×10 の完全
//! 崩壊」の第一試行。**結果は負** — 知的誠実性の規約 (CLAUDE.md #5: 失敗した
//! 検証も削除せず記録) に従い、負の結果と原因の同定を本ユニットの成果として
//! 凍結する。
//!
//! 試行した構成: エネルギー流の閉形式 J^i_unmod = ½{h, ∂ᵢh} (x 方向連続の式の
//! p → 0 極限から導出 — これ自体は正しい器械 [D0]) と置換則で
//! **T⁰ᵢ_naive := (V₀ᵢ + J^i)/2** (運動量密度とエネルギー流の素朴平均)。
//!
//! 負の発見 3 件 (全て ε/a-ladder で停留 = 収束した非零定数):
//!  [D1] T⁰ˣ_naive の**厳密流束** Φ := [h(k+qŷ)T⁰ˣ − T⁰ˣh(k)]/q̂ は point-split
//!       T^{xy} (Belinfante stress) に収束しない — rel = 1.30 で停留 (正準
//!       V₀x 単独の流束と同じ距離)。**素朴平均の流束は正準 stress のまま**。
//!  [D2] 対称 10×10 (T⁰ᵢ_naive + BOND-A/split 空間 6) の横断性破れ 0.674 は
//!       a 非依存 — v26.9-C run1 の正準混合と同水準。
//!  [D3] **J 片は殻上の横断性破れの組合せに寄与しない**: tv(naive) と
//!       tv(正準のみ) の差 < 1% — (E₁−E₂) 型の J 行列要素が対称殻で相殺する
//!       ことの器械的証拠。
//! 結論: **格子 Belinfante 0i 密度は密度・流束の平均では作れない — spin-current
//! 改良項 (σ^{0i} 型 Γ 積の局所頂点, ∂_λK^{λ0i} の格子転写) の明示構成が必要**。
//! 次ユニットとして登録。P₂ (spin-2) 部門は無傷 (oracle 回帰 [D4] PASS —
//! v26.9-C の縮退・アンカーは D/X のみで閉じているため影響なし)。
//!
//! 検査 (凍結 — 負の結果の認証):
//!  [D0] J 構成の器械: J^y_subst − 厳密流 = (q/8)[∂²_yh, h] + O(q²) — O(q)
//!       交換子項 (q 0.2 → 0.1 で ~1/2) かつ J^x の Fourier 台 ≤ 2 (< 1e-12)
//!  [D1] Φ(T⁰ˣ_naive) vs T^{xy}-split: ε-ladder で定数 (spread < 5%) かつ > 0.5
//!  [D2] 対称 10×10 の横断性: a-ladder で定数 (spread < 5%) かつ > 0.1
//!  [D3] J 片の殻寄与消滅: |tv(naive) − tv(正準)| / tv < 1% (最細 a)
//!  [D4] oracle 回帰: σ_DD/(2ρ_D) = 1 ± 2% かつ σ_DD/2σ_XX = 1 ± 2% (P₂ 無傷)
//!
//! 事前登録分岐: 本ユニットは判定 (b) = 負の結果の公表として構成 —
//!   「素朴平均は Belinfante 密度でない」が確定すれば PASS (堂々巡り防止:
//!   spin-current 改良の設計は独立ユニットで行い、本ユニットでは試行しない)。

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

/// ∂h/∂kᵢ (閉形式): 該当 Γᵢ 片の cos → −sin
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

/// J^i(k; qŷ) = ½{h(mid), ∂ᵢh(mid)} (置換則) — 実行列
fn j_subst(k: [f64; 3], q: f64, m: f64, ax: usize) -> Vec<f64> {
    let mid = [k[0], k[1] + 0.5 * q, k[2]];
    let h = h8(mid, m);
    let d = dh8(mid, ax);
    let hd = mat_mul(&h, &d);
    let dh = mat_mul(&d, &h);
    (0..64).map(|i| 0.5 * (hd[i] + dh[i])).collect()
}

/// 厳密 y 流 (v26.9-0): [h(k+qŷ)h(mid) − h(mid)h(k)]/(2sin(q/2)) — 実行列
fn j_exact_y(k: [f64; 3], q: f64, m: f64) -> Vec<f64> {
    let a = h8([k[0], k[1] + q, k[2]], m);
    let b = h8(k, m);
    let mm = h8([k[0], k[1] + 0.5 * q, k[2]], m);
    let am = mat_mul(&a, &mm);
    let mb = mat_mul(&mm, &b);
    let dd = 2.0 * (0.5 * q).sin();
    (0..64).map(|i| (am[i] - mb[i]) / dd).collect()
}

// ---- 空間 6 source (v26.9-C の認証済み写経) ----

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
            for ax in 0..3 {
                ph += (k[ax] + PI * sbit(s, ax) as f64) * t.d[ax] as f64;
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

/// Belinfante 10 source (実行列)。順: 00, 0x, 0y, 0z, xx, yy, zz, xy, xz, yz。
/// mutate: T⁰ˣ を正準 V₀x に置換 (v26.9-C run1 の病理 = 負対照)
fn source_bel(i: usize, k: [f64; 3], q: f64, m: f64, mutate: bool) -> Vec<f64> {
    let km = [k[0], k[1] + 0.5 * q, k[2]];
    match i {
        0 => h8(km, m),
        1 | 3 => {
            let ax = i - 1;
            let val = -0.5 * (2.0 * km[ax]).sin();
            let mut v0 = vec![0.0f64; 64];
            for s in 0..8usize {
                v0[s + s * 8] = val;
            }
            if mutate && i == 1 {
                return v0; // 正準のみ (対称化なし)
            }
            let j = j_subst(k, q, m, ax);
            (0..64).map(|x| 0.5 * (v0[x] + j[x])).collect()
        }
        2 => {
            let val = -0.5 * (2.0 * km[1]).sin();
            let mut v0 = vec![0.0f64; 64];
            for s in 0..8usize {
                v0[s + s * 8] = val;
            }
            let j = j_exact_y(k, q, m);
            (0..64).map(|x| 0.5 * (v0[x] + j[x])).collect()
        }
        4 => vertex_unmod(&[Term { eps: 0, d: [1, 0, 0], w: 0.5 }], km),
        5 => vertex_unmod(&[Term { eps: 1, d: [0, 1, 0], w: 0.5 }], km),
        6 => vertex_unmod(&[Term { eps: 3, d: [0, 0, 1], w: 0.5 }], km),
        7 => vertex_unmod(&t_split_terms(0, 1), km).iter().map(|x| x / 2.0).collect(),
        8 => vertex_unmod(&t_split_terms(0, 2), km).iter().map(|x| x / 2.0).collect(),
        _ => vertex_unmod(&t_split_terms(1, 2), km).iter().map(|x| x / 2.0).collect(),
    }
}

// ---- 殻積分 ----

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

/// 殻積分: 10×10 σ 行列と 4 行の横断性破れ (Belinfante 対)
/// 返値: (σ[100], tviol[4] 相対, λ 用の対角情報は σ に含む)
fn sigma_bel(a: f64, e_phys: f64, q_phys: f64, nth: usize, nph: usize, mutate: bool) -> (Vec<f64>, [f64; 4]) {
    let e_lat = a * e_phys;
    let q_lat = a * q_phys;
    let gl = gauss_legendre(nth);
    let c = PI / 2.0;
    let mut sig = vec![0.0f64; 100];
    let mut tv_num = [0.0f64; 4];
    let mut tv_den = [0.0f64; 4];
    let qhat = 2.0 * (0.5 * q_lat).sin();
    // Belinfante 横断性の (0ν, yν) 対: ν=0: (00, 0y) / ν=x: (0x, xy) /
    // ν=y: (0y, yy) / ν=z: (0z, yz)
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
            // 行列要素
            let hk = h8(k, 0.0);
            let (_, vk) = jacobi_eigh(&hk, 8);
            let kq = [k[0], k[1] + q_lat, k[2]];
            let hq = h8(kq, 0.0);
            let (_, vq) = jacobi_eigh(&hq, 8);
            let vs: Vec<Vec<f64>> = (0..10).map(|i| source_bel(i, k, q_lat, 0.0, mutate)).collect();
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
    println!("=== v26.9-D v269d_belinfante — Belinfante 対称 10×10 の完全崩壊 ===\n");
    println!("T⁰ᵢ_Bel = (V₀ᵢ + J^i)/2, J^i = ½{{h(mid), ∂ᵢh(mid)}} (y は厳密流)。massless の");
    println!("conformal 崩壊 σ → ρ₂P₂ のみ (trace セクター消滅) を殻積分 ladder で判定。\n");
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

    // ---- [D0] J 構成 ----
    {
        // (i) J^y_subst vs 厳密流: 差が q² 級
        let k = [0.7f64, -0.4, 1.9];
        let m = 0.2f64;
        let diff = |q: f64| -> f64 {
            let a = j_subst(k, q, m, 1);
            let b = j_exact_y(k, q, m);
            let mut d = 0.0f64;
            let mut nn = 0.0f64;
            for i in 0..64 {
                d += (a[i] - b[i]).powi(2);
                nn += b[i].powi(2);
            }
            (d / nn).sqrt()
        };
        let (d1, d2) = (diff(0.2), diff(0.1));
        let ratio = d1 / d2;
        // (ii) J^x の Fourier 台
        let mut worst = 0.0f64;
        for ax in 0..3usize {
            let ngrid = 16usize;
            let base = [0.37f64, -0.81, 1.13];
            let vsl: Vec<Vec<f64>> = (0..ngrid)
                .map(|j| {
                    let mut kk = base;
                    kk[ax] = 2.0 * PI * j as f64 / ngrid as f64;
                    j_subst(kk, 0.5, m, 0)
                })
                .collect();
            for r in 0..8 {
                for c in 0..8 {
                    let (mut hi2, mut all) = (0.0f64, 0.0f64);
                    for nmode in 0..ngrid {
                        let nsig = if nmode <= ngrid / 2 {
                            nmode as i32
                        } else {
                            nmode as i32 - ngrid as i32
                        };
                        let (mut cre, mut cim) = (0.0f64, 0.0f64);
                        for (j, vm) in vsl.iter().enumerate() {
                            let th = -2.0 * PI * (nmode * j) as f64 / ngrid as f64;
                            cre += vm[c + r * 8] * th.cos();
                            cim += vm[c + r * 8] * th.sin();
                        }
                        let pw = (cre * cre + cim * cim) / (ngrid * ngrid) as f64;
                        all += pw;
                        if nsig.abs() > 2 {
                            hi2 += pw;
                        }
                    }
                    if all > 1e-20 {
                        worst = worst.max((hi2 / all).sqrt());
                    }
                }
            }
        }
        // 導出: J_exact = ½{∂_yh, h}(mid) + (q/8)[∂²_yh, h] + O(q²) — 差は
        // O(q) の交換子項 (run1 は O(q²) と誤想定: 半減比の正解は ~2)
        check(
            "[D0] J 構成: J^y_subst − 厳密流 = O(q) 交換子項 (半減比 ~2) かつ J^x の Fourier 台 ≤ 2",
            (1.8..2.7).contains(&ratio) && worst < 1e-12,
            format!("q 半減比 = {:.2} (差 {:.1e} → {:.1e}), 台残差 = {:.1e}", ratio, d1, d2, worst),
        );
    }

    // ---- [D0b] 診断: T⁰ˣ_Bel の厳密流束 Φ vs point-split T^{xy} ----
    {
        let nvec = [0.55f64, 0.2, -0.81];
        let nn = (nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]).sqrt();
        let n = [nvec[0] / nn, nvec[1] / nn, nvec[2] / nn];
        let c = PI / 2.0;
        let mut msg = String::new();
        let mut rels = Vec::new();
        for &eps in &[0.4f64, 0.2, 0.1, 0.05] {
            let mphys = 0.5 * eps;
            let k = [c + eps * n[0], c + eps * n[1], c + eps * n[2]];
            let qq = 0.8 * eps;
            let dsrc = source_bel(1, k, qq, mphys, false);
            let a = h8([k[0], k[1] + qq, k[2]], mphys);
            let b = h8(k, mphys);
            let am = mat_mul(&a, &dsrc);
            let mb = mat_mul(&dsrc, &b);
            let dd = 2.0 * (0.5 * qq).sin();
            let phi: Vec<f64> = (0..64).map(|i| (am[i] - mb[i]) / dd).collect();
            let split = source_bel(7, k, qq, mphys, false);
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
            let r = block(&dif) / block(&split);
            msg = format!("{} rel({}) = {:.2e}", msg, eps, r);
            rels.push(r);
        }
        let mx = rels.iter().cloned().fold(0.0f64, f64::max);
        let mn = rels.iter().cloned().fold(f64::INFINITY, f64::min);
        check(
            "[D1] 負の発見: Φ(T⁰ˣ_naive) は T^{xy}-split に収束しない (定数 1.30 で停留)",
            (mx - mn) / mx < 0.05 && mn > 0.5,
            format!("{} (spread {:.1e})", msg, (mx - mn) / mx),
        );
    }

    // ---- [D1][D2][D3] ladder ----
    {
        let (e_phys, q_phys) = (1.5f64, 0.6);
        let s_inv = e_phys * e_phys - q_phys * q_phys;
        let rho_d = s_inv * s_inv / (160.0 * PI * PI);
        let mut rows = Vec::new();
        for &a in &[0.18f64, 0.09, 0.045] {
            let (sig, tv) = sigma_bel(a, e_phys, q_phys, 32, 64, false);
            let tvmax = tv.iter().cloned().fold(0.0f64, f64::max);
            // trace 方向: t̂ ∝ θ_{μν} 成分 (E,Q 面): θ = η − qq/q²
            let q2 = e_phys * e_phys - q_phys * q_phys;
            let qv = [e_phys, 0.0, q_phys, 0.0];
            let mut eta = [1.0f64, -1.0, -1.0, -1.0];
            let mut tvec = [0.0f64; 10];
            let map = [(0usize, 0usize), (0, 1), (0, 2), (0, 3), (1, 1), (2, 2), (3, 3), (1, 2), (1, 3), (2, 3)];
            for (i, &(mu, nu)) in map.iter().enumerate() {
                let th = (if mu == nu { eta[mu] } else { 0.0 }) - qv[mu] * qv[nu] / q2;
                tvec[i] = th * if mu == nu { 1.0 } else { (2.0f64).sqrt() };
            }
            let tn: f64 = tvec.iter().map(|x| x * x).sum::<f64>().sqrt();
            for v in tvec.iter_mut() {
                *v /= tn;
            }
            // σ·t̂ と λ_max (べき乗法)
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
            let trace_frac = stn / lmax;
            // oracle: σ_DD = (σ44 + σ66 − 2σ46)/2, σ_XX = σ88
            let sdd = 0.5 * (sig[4 + 4 * 10] + sig[6 + 6 * 10] - 2.0 * sig[6 + 4 * 10]);
            let sxx = sig[8 + 8 * 10];
            let anchor = sdd / a.powi(4) / (2.0 * rho_d);
            let deg = sdd / (2.0 * sxx);
            println!(
                "    [D 表] a = {:.3}: 横断性 max = {:.5}, trace 比 = {:.5}, アンカー = {:.4}, σ_DD/2σ_XX = {:.4} ({} s)",
                a, tvmax, trace_frac, anchor, deg, t0.elapsed().as_secs()
            );
            let _ = &mut eta;
            rows.push((a, tvmax, trace_frac, anchor, deg));
        }
        let sp = (rows[0].1 - rows[2].1).abs() / rows[2].1;
        check(
            "[D2] 負の発見: 対称 10×10 の横断性破れ 0.674 は a 非依存 (素朴平均 ≠ Belinfante 密度)",
            sp < 0.05 && rows[2].1 > 0.1,
            format!(
                "{:.4} → {:.4} → {:.4} (spread {:.1e}), trace 比も停留 {:.3}",
                rows[0].1, rows[1].1, rows[2].1, sp, rows[2].2
            ),
        );

        // ---- [D3] J 片の殻寄与消滅 ----
        {
            let (_, tv_can) = sigma_bel(0.045, e_phys, q_phys, 32, 64, true);
            let rel = (tv_can[1] - rows[2].1).abs() / rows[2].1;
            check(
                "[D3] 負の発見: J 片は殻の破れに寄与しない — tv(naive) = tv(正準) < 1%",
                rel < 0.01,
                format!("相対差 = {:.1e} (naive {:.4} vs 正準 {:.4})", rel, rows[2].1, tv_can[1]),
            );
        }

        let ok_or = (rows[2].3 - 1.0).abs() < 0.02 && (rows[2].4 - 1.0).abs() < 0.02;
        check(
            "[D4] oracle 回帰: σ_DD/(2ρ_D) = 1 ± 2% かつ σ_DD/2σ_XX = 1 ± 2% (P₂ 無傷)",
            ok_or,
            format!("アンカー = {:.4}, 縮退 = {:.4}", rows[2].3, rows[2].4),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.9-D".into())),
        ("kind".into(), Json::Str("belinfante_symmetric_collapse".into())),
        (
            "construction".into(),
            Json::Str("T0i_Bel = (V0i + J^i)/2, J^i = (1/2){h(mid), ∂ᵢh(mid)} (y は厳密流)".into()),
        ),
    ]);
    let p = write_artifact("results/v269d_belinfante.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "判定 (b) 確定: **負の結果 — 素朴平均 T⁰ᵢ = (密度+流束)/2 は格子 Belinfante 密度ではない** (厳密流束は正準 stress のまま・J 片は殻の破れに寄与しない)。正しい構成 = spin-current 改良頂点の明示転写 — 次ユニットに登録 (P₂ 部門は無傷)"
        } else {
            "FAIL あり — 負の結果の認証自体が破れた (器械を疑う)。欄が一次ソース"
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
