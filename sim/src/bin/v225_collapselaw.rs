//! v22.5 R₂ スケーリング崩壊法則 — 生存窓内 6 点の合併単調性検定 (第二十三期)
//!
//! v22.2 の発見 (R₂ が無次元深さ d·gap = d/ξ のみの関数) は有効 3 点、
//! v22.4 の点予言は床の前進 (d=2.0 でカナリア死, R₁ 劣化が N 非依存) で採点不能だった。
//! 本版は f64 床の生存窓 (d ≤ 1.5) の内側で法則を検定する:
//!   N ∈ {12, 14} × d ∈ {0.5, 1.0, 1.5} の 6 点 — d·gap は
//!   N=12: {0.121, 0.241, 0.362} / N=14 (gap ≈ 0.207): {~0.103, ~0.207, ~0.310}
//!   と交互に噛み合う。崩壊が真なら N を混ぜて d·gap で並べた列が単調になり、
//!   N=14 点は N=12 折れ線の ±0.025 に乗る (v22.2 で観測した崩壊散らばり)。
//! (v22.2 の N=16 d=1.5 点 (0.277, 0.82075) は引用記録 — 本版のゲートには使わない。)
//!
//! 事前登録: (a) 有効点 ≥ 5/6 ∧ 合併列が単調 (許容 −0.01) ∧ N=14 有効点の
//!   N=12 折れ線からの残差 max ≤ 0.025 = 崩壊法則の成立 /
//!   (a′) 有効点 < 5 = 器械不足の記録 / (b) 単調破れ or 残差超過 = 崩壊反証。

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
    println!("=== v22.5 R₂ スケーリング崩壊法則 — 生存窓 6 点の合併単調性 (第二十三期) ===\n");
    println!("事前登録: (a) 有効 ≥5/6 ∧ 合併単調 ∧ N=14 残差 ≤ 0.025 = 崩壊法則成立 /");
    println!("          (a′) 有効 < 5 = 器械不足 / (b) 単調破れ・残差超過 = 崩壊反証\n");
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
    let dists = [0.5f64, 1.0, 1.5];
    let mut rows: Vec<(usize, f64, f64, f64, f64, bool)> = Vec::new(); // (n, d, dgap, r1, r2, valid)

    for &n in &[12usize, 14] {
        let ns = n * n * n;
        let h = build_h3d_periodic(n);
        let (ev, vv) = jacobi_eigh(&h, ns);
        let nocc = ns / 2;
        let gap = ev[nocc] - ev[nocc - 1];
        check(
            &format!("N={} 閉殻ギャップ (m_T=0 構成)", n),
            gap > 1e-6,
            format!("gap = {:.4} ({} s)", gap, t0.elapsed().as_secs()),
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
        // Gauss 対 (v22.2 と同一規約) — 生存窓 d ∈ {0.5, 1.0, 1.5}
        let wid = 1.0f64;
        let mid = n as f64 / 2.0;
        for &d in &dists {
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
                "    N={} d={:.1} (d·gap={:.4}): R₁={:.5} [{}], R₂={:.5} ({} s)",
                n,
                d,
                d * gap,
                r1,
                if valid { "有効" } else { "除外" },
                r2,
                t0.elapsed().as_secs()
            );
            rows.push((n, d, d * gap, r1, r2, valid));
        }
    }

    // ---- 判定 ----
    let valid_rows: Vec<&(usize, f64, f64, f64, f64, bool)> = rows.iter().filter(|r| r.5).collect();
    let nv = valid_rows.len();
    // 合併列 (d·gap 昇順) の単調性
    let mut sorted: Vec<&&(usize, f64, f64, f64, f64, bool)> = valid_rows.iter().collect();
    sorted.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap());
    let mut mono = true;
    for i in 1..sorted.len() {
        if sorted[i].4 < sorted[i - 1].4 - 0.01 {
            mono = false;
        }
    }
    // N=14 有効点の N=12 折れ線からの残差
    let base: Vec<(f64, f64)> = rows
        .iter()
        .filter(|r| r.0 == 12 && r.5)
        .map(|r| (r.2, r.4))
        .collect();
    let interp = |x: f64| -> Option<f64> {
        if base.len() < 2 {
            return None;
        }
        for w in base.windows(2) {
            if x >= w[0].0 - 0.02 && x <= w[1].0 + 0.02 {
                let t = (x - w[0].0) / (w[1].0 - w[0].0);
                return Some(w[0].1 + t * (w[1].1 - w[0].1));
            }
        }
        None
    };
    let mut max_res = 0.0f64;
    let mut n_scored = 0usize;
    for r in rows.iter().filter(|r| r.0 == 14 && r.5) {
        if let Some(f) = interp(r.2) {
            let res = (r.4 - f).abs();
            println!(
                "    [採点] N=14 d={:.1} (d·gap={:.4}): R₂={:.5} vs N=12 折れ線 {:.5} (残差 {:+.4})",
                r.1,
                r.2,
                r.4,
                f,
                r.4 - f
            );
            max_res = max_res.max(res);
            n_scored += 1;
        } else {
            println!(
                "    [記録] N=14 d={:.1} (d·gap={:.4}): R₂={:.5} — N=12 折れ線の外 (採点なし)",
                r.1, r.2, r.4
            );
        }
    }
    println!("    [引用記録] v22.2 の N=16 d=1.5: (d·gap, R₂) = (0.277, 0.82075)");
    let branch_a = nv >= 5 && mono && n_scored >= 1 && max_res <= 0.025;
    println!(
        "\n[判定] {}",
        if branch_a {
            "事前登録 (a): 崩壊法則の成立 — R₂ = f(d/ξ) が N を跨いで単調・±0.025 で一本の曲線"
        } else if nv < 5 {
            "事前登録 (a′): 有効点 < 5 — 器械不足の記録"
        } else {
            "事前登録 (b): 単調破れまたは残差超過 — 崩壊反証の記録"
        }
    );
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v22.5".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(n, d, dg, r1, r2, va)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("d".into(), Json::Num(d)),
                            ("dgap".into(), Json::Num(dg)),
                            ("r1".into(), Json::Num(r1)),
                            ("r2".into(), Json::Num(r2)),
                            ("valid".into(), Json::Bool(va)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("max_res".into(), Json::Num(max_res)),
        ("branch_a".into(), Json::Bool(branch_a)),
    ]);
    let p = write_artifact("results/v225_collapselaw.json", &j.render());
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
