//! v22.6 K_boost 離散化の変種比較 — R₂ 不足の残る説明系を検定する (第二十三期)
//!
//! v22.5 で R₂ 不足の幾何仮説は全滅し、残る説明系は「ブースト重み K_boost の格子
//! 離散化がシェル依存の O(1) 補正を持つ」だけになった。本版はそれを直接検定する:
//! 生存窓の 4 点 (N=12 × d ∈ {0.5, 1.0, 1.5} + N=14 × d=1.0) で、K_boost の
//! 離散化変種ごとに R₂ を測る。
//!   [A] 現行 (v19.2 凍結): 右楔のみ (ξ > 0)、x ボンドは中点 ξ・y/z ボンドはサイト ξ
//!   [B] 両側反対称: ξ の符号ごと全ボンドを含める (左楔は負の重み — Rindler の完全形)
//!   [C] セル平均: staggered 2 セル単位で ξ をセル中心にスナップ (倍加子の平均化)
//! 事前登録: (a) いずれかの変種が 4 点全てで |R₂ − 1| ≤ 0.05 = 離散化が原因と確定
//!   + 改良形の同定 / (a′) N=12 d=1.0 と N=14 d=1.0 の変種内差が [A] の半分以下に
//!   縮む = 部分説明 / (b) どの変種も改善しない = K_boost 重みは原因でない (これも
//!   残る説明系を一つ消す決定的情報)。
//! 器械: v22.5 と同一 (Gauss 対 rank-2 回転・R₁ カナリア |R₁−1| < 1%)。

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
    println!("=== v22.6 K_boost 離散化の変種比較 (第二十三期) ===\n");
    println!("事前登録: (a) ある変種が 4 点全て |R₂−1| ≤ 0.05 = 離散化が原因 /");
    println!("          (a′) 変種内の N 差が [A] の半分以下 = 部分説明 / (b) 改善なし\n");
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
    let mut rows: Vec<(usize, f64, usize, f64, f64, bool)> = Vec::new(); // (n, d, variant, r1, r2, valid)

    for &n in &[12usize, 14] {
        let dists: Vec<f64> = if n == 12 {
            vec![0.5, 1.0, 1.5]
        } else {
            vec![1.0]
        };
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
        // 3 変種のボンド集合: [A] 現行 (右楔・サイト ξ), [B] 両側反対称, [C] セル平均
        let mut variants: Vec<Vec<(usize, usize, f64, f64)>> =
            vec![Vec::new(), Vec::new(), Vec::new()];
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let i = idx(x, y, z);
                    if x + 1 < n {
                        let xi = xplane - (x as f64 + 0.5);
                        let j = idx(x + 1, y, z);
                        if xi > 1e-9 {
                            variants[0].push((i, j, 0.5, xi));
                            variants[2].push((i, j, 0.5, xi));
                        }
                        if xi.abs() > 1e-9 {
                            variants[1].push((i, j, 0.5, xi));
                        }
                    }
                    let xi_site = xplane - x as f64;
                    let xc = (2 * (x / 2)) as f64 + 0.5; // staggered 2 セルの中心
                    let xi_cell = xplane - xc;
                    let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                    let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                    let jy = idx(x, (y + 1) % n, z);
                    let jz = idx(x, y, (z + 1) % n);
                    if xi_site > 1e-9 {
                        variants[0].push((i, jy, ey, xi_site));
                        variants[0].push((i, jz, ez, xi_site));
                    }
                    if xi_site.abs() > 1e-9 {
                        variants[1].push((i, jy, ey, xi_site));
                        variants[1].push((i, jz, ez, xi_site));
                    }
                    if xi_cell > 1e-9 {
                        variants[2].push((i, jy, ey, xi_cell));
                        variants[2].push((i, jz, ez, xi_cell));
                    }
                }
            }
        }
        let two_pi = 2.0 * std::f64::consts::PI;
        let kboost_of = |cc: &[f64], vi: usize| -> f64 {
            let mut acc = 0.0;
            for &(i, j, t, xi) in &variants[vi] {
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
            let mut r2e = [[0.0f64; 2]; 3];
            for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
                let cp = rotate_c_local(&c, ns, &u, &w, al);
                let cm = rotate_c_local(&c, ns, &u, &w, -al);
                let cap = restrict_local(&cp, ns, &sel);
                let cam = restrict_local(&cm, ns, &sel);
                let ds = (entropy_sym_local(&cap, m) - entropy_sym_local(&cam, m)) / 2.0;
                let dka = (dk_of(&cap) - dk_of(&cam)) / 2.0;
                r1e[ii] = ds / dka;
                for vi in 0..3 {
                    let dkb = (kboost_of(&cp, vi) - kboost_of(&cm, vi)) / 2.0;
                    r2e[vi][ii] = ds / dkb;
                }
            }
            let r1 = (4.0 * r1e[0] - r1e[1]) / 3.0;
            let valid = (r1 - 1.0).abs() < 0.01;
            let mut line = format!(
                "    N={} d={:.1}: R₁={:.5} [{}]",
                n,
                d,
                r1,
                if valid { "有効" } else { "除外" }
            );
            for (vi, name) in ["A", "B", "C"].iter().enumerate() {
                let r2 = (4.0 * r2e[vi][0] - r2e[vi][1]) / 3.0;
                line.push_str(&format!(", R₂[{}]={:.5}", name, r2));
                rows.push((n, d, vi, r1, r2, valid));
            }
            println!("{} ({} s)", line, t0.elapsed().as_secs());
        }
    }

    // ---- 判定 ----
    let mut spread_a = 0.0f64;
    let mut branch_a = false;
    let mut branch_ap = false;
    for vi in 0..3 {
        let pts: Vec<&(usize, f64, usize, f64, f64, bool)> =
            rows.iter().filter(|r| r.2 == vi && r.5).collect();
        if pts.len() < 4 {
            continue;
        }
        let maxdev = pts.iter().map(|r| (r.4 - 1.0).abs()).fold(0.0f64, f64::max);
        let r12 = pts
            .iter()
            .find(|r| r.0 == 12 && (r.1 - 1.0).abs() < 0.01)
            .map(|r| r.4)
            .unwrap_or(0.0);
        let r14 = pts.iter().find(|r| r.0 == 14).map(|r| r.4).unwrap_or(0.0);
        let spread = (r12 - r14).abs();
        if vi == 0 {
            spread_a = spread;
        }
        println!(
            "    [変種 {}] max|R₂−1| = {:.4}, N 間差 (d=1.0) = {:.4}",
            ["A", "B", "C"][vi],
            maxdev,
            spread
        );
        if maxdev <= 0.05 {
            branch_a = true;
        }
        if vi > 0 && spread <= 0.5 * spread_a {
            branch_ap = true;
        }
    }
    println!(
        "\n[判定] {}",
        if branch_a {
            "事前登録 (a): ある変種が 4 点 |R₂−1| ≤ 0.05 — K_boost 離散化が原因と確定・改良形を同定"
        } else if branch_ap {
            "事前登録 (a′): 変種で N 間差が半減 — 部分説明の記録"
        } else {
            "事前登録 (b): どの変種も改善せず — K_boost 重みは原因でない (残る説明系が一つ消えた)"
        }
    );
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v22.6".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(n, d, vi, r1, r2, va)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("d".into(), Json::Num(d)),
                            ("variant".into(), Json::Int(vi as i64)),
                            ("r1".into(), Json::Num(r1)),
                            ("r2".into(), Json::Num(r2)),
                            ("valid".into(), Json::Bool(va)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(branch_a)),
    ]);
    let p = write_artifact("results/v226_kboost.json", &j.render());
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
