//! v22.4 R₂ スケーリング崩壊の点予言検定 — G 窓の障害を法則として固定する (第二十三期)
//!
//! v22.2 の発見: ブースト比 R₂ が横質量 m_T でも体積 N でもなく無次元深さ
//! d·gap (= d/ξ) のスケーリング関数に乗る。カナリア有効 3 点:
//!   (d·gap, R₂) = (0.121, 0.79255), (0.362, 0.86452) [N=12], (0.277, 0.82075) [N=16]
//! N=16 点は N=12 曲線の線形補間から −2.2% — 崩壊精度 ~2%。
//!
//! 本版は崩壊仮説の**事前登録された点予言**を検定する:
//!   N=16, d=2.0 は d·gap = 2.0×0.1845 = 0.369 で N=12, d=1.5 (0.362) とほぼ同深さ
//!   → 崩壊が真なら R₂(16, 2.0) = 0.8667 ± 0.025 (N=12 曲線の補間 + 観測散らばり)。
//!   異なる (N, d) で同じ d·gap → 同じ R₂ — これが崩壊の定義そのもの。
//! 同時に N=12, d=2.0 (d·gap = 0.482, 外挿予測 ~0.90 — 記録のみ) と、カナリア死境界
//! (v22.2: d=2.5 で両 N とも死, R₁ ≈ 0.945 が N 非依存 = f64 κ 床の d-局所性) の
//! d=2.0 での前哨を測る。
//!
//! 事前登録: (a) R₁(16, 2.0) 有効 (|R₁−1| < 1%) ∧ |R₂(16, 2.0) − 0.8667| ≤ 0.025
//!   = 崩壊確認 — R₂ = f(d/ξ) が法則に昇格 (G 窓の障害 = 「f≳0.93 に必要な
//!   d/ξ ≳ 0.7 が f64 κ 床の向こう」が定量文になる) /
//!   (a′) 両点ともカナリア死 = 床が d=2.0 まで前進 (それ自体を記録) /
//!   (b) カナリア有効なのに点予言外れ = 崩壊反証 — 記録。
//! (N=12, d=2.0 は外挿域なので記録のみ — ゲートには使わない。)

use uft_sim::*;

// 3D staggered H (x 開放, y/z 周期) — v22.2 と同一
fn build_h3d_periodic(n: usize) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let mut add = |i: usize, j: usize, t: f64| {
        h[j + i * ns] += t;
        h[i + j * ns] += t;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                if x + 1 < n {
                    add(i, idx(x + 1, y, z), 0.5);
                }
                let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                add(i, idx(x, (y + 1) % n, z), ey);
                let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                add(i, idx(x, y, (z + 1) % n), ez);
            }
        }
    }
    h
}

fn entropy_sym_local(c: &[f64], m: usize) -> f64 {
    let (w, _) = jacobi_eigh(c, m);
    w.iter()
        .map(|&x| {
            let x = x.clamp(1e-14, 1.0 - 1e-14);
            -x * x.ln() - (1.0 - x) * (1.0 - x).ln()
        })
        .sum()
}

fn restrict_local(c: &[f64], ns: usize, sel: &[usize]) -> Vec<f64> {
    let m = sel.len();
    let mut o = vec![0.0f64; m * m];
    for a in 0..m {
        for b in 0..m {
            o[a + b * m] = c[sel[a] + sel[b] * ns];
        }
    }
    o
}

// rank-2 実回転 C' = OCOᵀ (v19.1 系, v22.2 と同一)
fn rotate_c_local(c: &[f64], ns: usize, u: &[f64], w: &[f64], alpha: f64) -> Vec<f64> {
    let (s, cs) = (alpha.sin(), alpha.cos());
    let cu: Vec<f64> = (0..ns)
        .map(|i| (0..ns).map(|j| c[j + i * ns] * u[j]).sum())
        .collect();
    let cw: Vec<f64> = (0..ns)
        .map(|i| (0..ns).map(|j| c[j + i * ns] * w[j]).sum())
        .collect();
    let uu: f64 = (0..ns).map(|i| u[i] * cu[i]).sum();
    let uw: f64 = (0..ns).map(|i| u[i] * cw[i]).sum();
    let ww: f64 = (0..ns).map(|i| w[i] * cw[i]).sum();
    let mut out = c.to_vec();
    let a = cs - 1.0;
    for i in 0..ns {
        let du = a * u[i] + s * w[i];
        let dw = a * w[i] - s * u[i];
        for j in 0..ns {
            out[j + i * ns] += du * cu[j] + dw * cw[j];
            out[i + j * ns] += du * cu[j] + dw * cw[j];
        }
    }
    for i in 0..ns {
        let dui = a * u[i] + s * w[i];
        let dwi = a * w[i] - s * u[i];
        for j in 0..ns {
            let duj = a * u[j] + s * w[j];
            let dwj = a * w[j] - s * u[j];
            out[j + i * ns] += dui * (uu * duj + uw * dwj) + dwi * (uw * duj + ww * dwj);
        }
    }
    out
}

fn main() {
    self_test();
    println!("=== v22.4 R₂ スケーリング崩壊の点予言検定 (第二十三期) ===\n");
    println!("事前登録: (a) R₁(16,2.0) 有効 ∧ R₂(16,2.0) = 0.8667 ± 0.025 = 崩壊確認 /");
    println!("          (a′) 両点カナリア死 = 床の前進を記録 / (b) 有効なのに外れ = 崩壊反証\n");
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
    let t0 = std::time::Instant::now();
    let d = 2.0f64;
    let mut rows: Vec<(usize, f64, f64, f64, bool)> = Vec::new(); // (n, dgap, r1, r2, valid)

    for &n in &[12usize, 16] {
        let ns = n * n * n;
        let h = build_h3d_periodic(n);
        let (ev, vv) = jacobi_eigh(&h, ns);
        let nocc = ns / 2;
        let gap = ev[nocc] - ev[nocc - 1];
        check(
            &format!("N={} 閉殻ギャップ (v22.2 の再現)", n),
            gap > 1e-6,
            format!(
                "gap = {:.4}, d·gap = {:.4} ({} s)",
                gap,
                d * gap,
                t0.elapsed().as_secs()
            ),
        );
        let mut c = vec![0.0f64; ns * ns];
        for k in 0..nocc {
            for i in 0..ns {
                let vi = vv[i + k * ns];
                if vi == 0.0 {
                    continue;
                }
                for j in 0..ns {
                    c[j + i * ns] += vi * vv[j + k * ns];
                }
            }
        }
        let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
        let mut sel = Vec::new();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n / 2 {
                    sel.push(idx(x, y, z));
                }
            }
        }
        let m = sel.len();
        let ca0 = restrict_local(&c, ns, &sel);
        let (cw, cvec) = jacobi_eigh(&ca0, m);
        let kappa: Vec<f64> = cw
            .iter()
            .map(|&ck| {
                let x2 = ck.clamp(1e-14, 1.0 - 1e-14);
                ((1.0 - x2) / x2).ln()
            })
            .collect();
        let dk_of = |ca: &[f64]| -> f64 {
            let mut dk = 0.0;
            for k in 0..m {
                let mut acc = 0.0;
                for i in 0..m {
                    let ui = cvec[i + k * m];
                    if ui == 0.0 {
                        continue;
                    }
                    for j in 0..m {
                        acc += ui * (ca[j + i * m] - ca0[j + i * m]) * cvec[j + k * m];
                    }
                }
                dk += kappa[k] * acc;
            }
            dk
        };
        let xplane = n as f64 / 2.0 - 0.5;
        let mut bonds: Vec<(usize, usize, f64, f64)> = Vec::new();
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let i = idx(x, y, z);
                    if x + 1 < n {
                        let xi = xplane - (x as f64 + 0.5);
                        if xi > 1e-9 {
                            bonds.push((i, idx(x + 1, y, z), 0.5, xi));
                        }
                    }
                    let xi_site = xplane - x as f64;
                    if xi_site > 1e-9 {
                        let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                        bonds.push((i, idx(x, (y + 1) % n, z), ey, xi_site));
                        let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                        bonds.push((i, idx(x, y, (z + 1) % n), ez, xi_site));
                    }
                }
            }
        }
        let two_pi = 2.0 * std::f64::consts::PI;
        let kboost_of = |cc: &[f64]| -> f64 {
            let mut acc = 0.0;
            for &(i, j, t, xi) in &bonds {
                acc += two_pi * xi * t * 2.0 * (cc[j + i * ns] - c[j + i * ns]);
            }
            acc
        };
        // Gauss 対 (v22.2 と同一規約, d = 2.0)
        let wid = 1.0f64;
        let mid = n as f64 / 2.0;
        let mk = |x0: f64| -> Vec<f64> {
            let mut p = vec![0.0f64; ns];
            for x in 0..n {
                for y in 0..n {
                    for z in 0..n {
                        let d2 = (x as f64 - x0).powi(2)
                            + (y as f64 - mid).powi(2)
                            + (z as f64 - mid).powi(2);
                        p[idx(x, y, z)] = (-d2 / (2.0 * wid * wid)).exp();
                    }
                }
            }
            let nr: f64 = p.iter().map(|a| a * a).sum::<f64>().sqrt();
            for a in p.iter_mut() {
                *a /= nr;
            }
            p
        };
        let u = mk(xplane - d);
        let mut w = mk(xplane - d - 1.0);
        let ov: f64 = u.iter().zip(w.iter()).map(|(a, b)| a * b).sum();
        for i in 0..ns {
            w[i] -= ov * u[i];
        }
        let nr: f64 = w.iter().map(|a| a * a).sum::<f64>().sqrt();
        for a in w.iter_mut() {
            *a /= nr;
        }
        let a1 = 0.02f64;
        let mut r1e = [0.0f64; 2];
        let mut r2e = [0.0f64; 2];
        for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
            let cp = rotate_c_local(&c, ns, &u, &w, al);
            let cm = rotate_c_local(&c, ns, &u, &w, -al);
            let cap = restrict_local(&cp, ns, &sel);
            let cam = restrict_local(&cm, ns, &sel);
            let ds = (entropy_sym_local(&cap, m) - entropy_sym_local(&cam, m)) / 2.0;
            let dka = (dk_of(&cap) - dk_of(&cam)) / 2.0;
            let dkb = (kboost_of(&cp) - kboost_of(&cm)) / 2.0;
            r1e[ii] = ds / dka;
            r2e[ii] = ds / dkb;
        }
        let r1 = (4.0 * r1e[0] - r1e[1]) / 3.0;
        let r2 = (4.0 * r2e[0] - r2e[1]) / 3.0;
        let valid = (r1 - 1.0).abs() < 0.01;
        println!(
            "    N={} d=2.0 (d·gap={:.4}): R₁={:.5} [{}], R₂={:.5} ({} s)",
            n,
            d * gap,
            r1,
            if valid { "有効" } else { "除外" },
            r2,
            t0.elapsed().as_secs()
        );
        rows.push((n, d * gap, r1, r2, valid));
    }

    // ---- 判定 ----
    let pred = 0.8667f64;
    let tol = 0.025f64;
    let p16 = rows.iter().find(|r| r.0 == 16).unwrap();
    let p12 = rows.iter().find(|r| r.0 == 12).unwrap();
    let branch_a = p16.4 && (p16.3 - pred).abs() <= tol;
    let branch_ap = !p16.4 && !p12.4;
    if p16.4 {
        println!(
            "    [採点] 点予言 R₂(16, 2.0) = 0.8667 ± 0.025: 実測 {:.5} (差 {:+.4}) → {}",
            p16.3,
            p16.3 - pred,
            if (p16.3 - pred).abs() <= tol {
                "的中"
            } else {
                "外れ"
            }
        );
    } else {
        println!(
            "    [記録] N=16 d=2.0 はカナリア死 (R₁ = {:.5}) — 点予言は採点不能",
            p16.2
        );
    }
    println!(
        "    [記録] N=12 d=2.0 (外挿域): R₁ = {:.5} [{}], R₂ = {:.5} (外挿予測 ~0.90)",
        p12.2,
        if p12.4 { "有効" } else { "除外" },
        p12.3
    );
    println!(
        "\n[判定] {}",
        if branch_a {
            "事前登録 (a): 崩壊確認 — R₂ = f(d/ξ) が法則に昇格。G 窓の障害 = 「f ≳ 0.93 に必要な d/ξ ≳ 0.7 が f64 κ 床の向こう」が定量文になった"
        } else if branch_ap {
            "事前登録 (a′): 両点カナリア死 — f64 床は d=2.0 まで前進 (境界の記録)"
        } else {
            "事前登録 (b): 崩壊反証または部分的カナリア死 — 記録"
        }
    );

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v22.4".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(n, dg, r1, r2, va)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("dgap".into(), Json::Num(dg)),
                            ("r1".into(), Json::Num(r1)),
                            ("r2".into(), Json::Num(r2)),
                            ("valid".into(), Json::Bool(va)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("pred".into(), Json::Num(pred)),
        ("branch_a".into(), Json::Bool(branch_a)),
    ]);
    let p = write_artifact("results/v224_collapse.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(a′)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
