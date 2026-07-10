//! v19.4 ブースト同定の壁距離掃引 — 第二段の逸脱の機構分離
//!
//! v19.2 (N 走査) は核なしブースト比 R₂ = δS/δ⟨K_boost⟩ が 1 へ収束しない逸脱を
//! 示した (N=8: 0.80, N=12: 0.77 — 壁密着波束 wid=N/8)。本版はその機構を分離する:
//! 格子エンタングルメント・ハミルトニアンは壁から O(a) の距離で非局所項を持つ
//! (Eisler–Peschel / Arias–Casini–Huerta–Pontello の格子系の知見) ため、壁密着の
//! 摂動では局所ブースト形 2π Σ ξ T₀₀ が過大応答する — なら壁から離せば回復するはず。
//! N=12 固定・格子幅 wid=1.0 固定で、x 隣接波束対を壁距離 d ∈ {0.5,1.5,2.5,3.5} に
//! 置き、R₁ (第一法則 — 核床カナリア)・R₂・R_K を測る。
//! 事前登録は v19.2 の N=16 完走前に凍結 (盲検性): (a) 有効点で R₂(d) 単調増加かつ
//! 最遠有効点 > 0.90 = 壁近傍格子効果 / (b) プラトー・非単調 = 局所形不整合の記録。

use uft_sim::*;

/// 3+1D staggered ハミルトニアン (実対称, x 開放・y,z 反周期)
fn build_h3d(n: usize) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let mut add = |i: usize, j: usize, v: f64| {
        h[j + i * ns] += v;
        h[i + j * ns] += v;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                if x + 1 < n {
                    add(i, idx(x + 1, y, z), 0.5);
                }
                let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                if y + 1 < n {
                    add(i, idx(x, y + 1, z), ey);
                } else {
                    add(i, idx(x, 0, z), -ey);
                }
                let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                if z + 1 < n {
                    add(i, idx(x, y, z + 1), ez);
                } else {
                    add(i, idx(x, y, 0), -ez);
                }
            }
        }
    }
    h
}

/// フェルミオン・エントロピー (実対称 C の固有値から)
fn entropy_sym(c: &[f64], m: usize) -> f64 {
    let (w, _v) = jacobi_eigh(c, m);
    let mut s = 0.0;
    for &ck in &w {
        let x = ck.clamp(1e-14, 1.0 - 1e-14);
        s += -x * x.ln() - (1.0 - x) * (1.0 - x).ln();
    }
    s
}

/// 領域 A の C_A
fn restrict(c: &[f64], ns: usize, sel: &[usize]) -> Vec<f64> {
    let m = sel.len();
    let mut ca = vec![0.0f64; m * m];
    for (i, &si) in sel.iter().enumerate() {
        for (j, &sj) in sel.iter().enumerate() {
            ca[j + i * m] = c[sj + si * ns];
        }
    }
    ca
}

/// span{u, w} 平面の直交回転 O(α) による C(α) = O C Oᵀ (u ⊥ w, 正規化済み前提)。
/// ΔO = (cosα−1)(uuᵀ + wwᵀ) + sinα (uwᵀ − wuᵀ) — rank 2 なので外積で厳密に構成。
fn rotate_c(c: &[f64], ns: usize, u: &[f64], w: &[f64], alpha: f64) -> Vec<f64> {
    let (sa, cosa) = alpha.sin_cos();
    let (p, q) = (cosa - 1.0, sa);
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
    let add_outer = |a: &[f64], b: &[f64], s: f64, out: &mut Vec<f64>| {
        if s == 0.0 {
            return;
        }
        for i in 0..ns {
            let ai = a[i] * s;
            if ai == 0.0 {
                continue;
            }
            for j in 0..ns {
                out[j + i * ns] += ai * b[j];
            }
        }
    };
    // ΔO C + (ΔO C)ᵀ
    add_outer(u, &cu, p, &mut out);
    add_outer(&cu, u, p, &mut out);
    add_outer(w, &cw, p, &mut out);
    add_outer(&cw, w, p, &mut out);
    add_outer(u, &cw, q, &mut out);
    add_outer(&cw, u, q, &mut out);
    add_outer(w, &cu, -q, &mut out);
    add_outer(&cu, w, -q, &mut out);
    // ΔO C ΔOᵀ (span{u,w} 内のスカラーで閉じる)
    let r1u = p * uu + q * uw;
    let r1w = p * uw + q * ww;
    let r2u = p * uw - q * uu;
    let r2w = p * ww - q * uw;
    let s1u = p * r1u + q * r1w;
    let s1w = p * r1w - q * r1u;
    let s2u = p * r2u + q * r2w;
    let s2w = p * r2w - q * r2u;
    add_outer(u, u, s1u, &mut out);
    add_outer(u, w, s1w, &mut out);
    add_outer(w, u, s2u, &mut out);
    add_outer(w, w, s2w, &mut out);
    out
}

fn main() {
    self_test();
    println!(
        "=== v19.4 ブースト同定の壁距離掃引: 逸脱は壁近傍の格子効果か (第二段の機構分離) ===\n"
    );
    println!("事前登録 (v19.2 の N=16 完走前に凍結 — N=8,12 の逸脱 R₂≈0.77-0.80 のみ既知):");
    println!(
        "  R₁ = 1 ± 1% を点ごとのカナリアとする (外れた点は核床汚染として除外 — 除外も記録)。"
    );
    println!("  (a) 有効点で R₂(d) が単調増加 (許容 −0.01) かつ最遠有効点 R₂ > 0.90");
    println!("      = 逸脱は壁近傍の格子非局所性 — ブースト同定は壁から離れて回復 /");
    println!("  (b) それ以外 (プラトー・非単調) = 局所ブースト形自体の不整合を記録\n");
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

    let two_pi = 2.0 * std::f64::consts::PI;
    let n = 12usize;
    let ns = n * n * n;
    let t0 = std::time::Instant::now();
    let h = build_h3d(n);
    let (w_ev, v) = jacobi_eigh(&h, ns);
    let nocc = ns / 2;
    let gap = w_ev[nocc] - w_ev[nocc - 1];
    check(
        &format!("N={} の閉殻ギャップ", n),
        gap > 1e-6,
        format!("gap = {:.4} ({} s)", gap, t0.elapsed().as_secs()),
    );
    let mut c = vec![0.0f64; ns * ns];
    for k in 0..nocc {
        for i in 0..ns {
            let vi = v[i + k * ns];
            if vi == 0.0 {
                continue;
            }
            for j in 0..ns {
                c[j + i * ns] += vi * v[j + k * ns];
            }
        }
    }
    let mut dev: f64 = 0.0;
    for i in 0..ns {
        for j in 0..ns {
            let mut acc = 0.0;
            for k in 0..ns {
                acc += c[k + i * ns] * c[j + k * ns];
            }
            dev = dev.max((acc - c[j + i * ns]).abs());
        }
    }
    check(
        &format!("N={} の C²=C", n),
        dev < 1e-10,
        format!("最大偏差 {:.1e}", dev),
    );
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
    let ca0 = restrict(&c, ns, &sel);
    let (cw_ev, cv) = jacobi_eigh(&ca0, m);
    let kappa: Vec<f64> = cw_ev
        .iter()
        .map(|&ck| {
            let x = ck.clamp(1e-14, 1.0 - 1e-14);
            ((1.0 - x) / x).ln()
        })
        .collect();
    let dk_of = |ca: &[f64]| -> f64 {
        let mut dk = 0.0;
        for k in 0..m {
            let mut acc = 0.0;
            for i in 0..m {
                let ui = cv[i + k * m];
                if ui == 0.0 {
                    continue;
                }
                for j in 0..m {
                    acc += ui * (ca[j + i * m] - ca0[j + i * m]) * cv[j + k * m];
                }
            }
            dk += kappa[k] * acc;
        }
        dk
    };
    // ブースト重みつきボンド一覧 — 規約は v19.2 から凍結 (同一コード)
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
                    if y + 1 < n {
                        bonds.push((i, idx(x, y + 1, z), ey, xi_site));
                    } else {
                        bonds.push((i, idx(x, 0, z), -ey, xi_site));
                    }
                    let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                    if z + 1 < n {
                        bonds.push((i, idx(x, y, z + 1), ez, xi_site));
                    } else {
                        bonds.push((i, idx(x, y, 0), ez * -1.0, xi_site));
                    }
                }
            }
        }
    }
    let kboost_of = |cc: &[f64]| -> f64 {
        let mut acc = 0.0;
        for &(i, j, t, xi) in &bonds {
            acc += two_pi * xi * t * 2.0 * (cc[j + i * ns] - c[j + i * ns]);
        }
        acc
    };
    // 波束対: 格子幅 wid = 1.0 固定 (格子距離を掃く — 物理幅スケーリングではない)
    let wid = 1.0f64;
    let make_packet = |x0: f64, y0: f64, z0: f64| -> Vec<f64> {
        let mut p = vec![0.0f64; ns];
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let d2 =
                        (x as f64 - x0).powi(2) + (y as f64 - y0).powi(2) + (z as f64 - z0).powi(2);
                    p[idx(x, y, z)] = (-d2 / (2.0 * wid * wid)).exp();
                }
            }
        }
        let nrm: f64 = p.iter().map(|a| a * a).sum::<f64>().sqrt();
        for a in p.iter_mut() {
            *a /= nrm;
        }
        p
    };
    let orth = |u: &[f64], w: &mut Vec<f64>| {
        let ov: f64 = u.iter().zip(w.iter()).map(|(a, b)| a * b).sum();
        for i in 0..u.len() {
            w[i] -= ov * u[i];
        }
        let nrm: f64 = w.iter().map(|a| a * a).sum::<f64>().sqrt();
        for a in w.iter_mut() {
            *a /= nrm;
        }
    };
    let mid = n as f64 / 2.0;
    let dists = [0.5f64, 1.5, 2.5, 3.5];
    let mut rows: Vec<(f64, f64, f64, f64, bool)> = Vec::new(); // (d, R1, R2, RK, valid)
    for &d in &dists {
        let u = make_packet(xplane - d, mid, mid);
        let mut wv = make_packet(xplane - d - 1.0, mid, mid);
        orth(&u, &mut wv);
        let eval_at = |alpha: f64| -> (f64, f64, f64) {
            let cc = rotate_c(&c, ns, &u, &wv, alpha);
            let ca = restrict(&cc, ns, &sel);
            (entropy_sym(&ca, m), dk_of(&ca), kboost_of(&cc))
        };
        let a1 = 0.02f64;
        let mut r1e = [0.0f64; 2];
        let mut r2e = [0.0f64; 2];
        let mut rke = [0.0f64; 2];
        for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
            let (sp, kp, bp) = eval_at(al);
            let (sm, km, bm) = eval_at(-al);
            let ds = (sp - sm) / 2.0;
            let dka = (kp - km) / 2.0;
            let dkb = (bp - bm) / 2.0;
            r1e[ii] = ds / dka;
            r2e[ii] = ds / dkb;
            rke[ii] = dkb / dka;
        }
        let r1 = (4.0 * r1e[0] - r1e[1]) / 3.0;
        let r2 = (4.0 * r2e[0] - r2e[1]) / 3.0;
        let rk = (4.0 * rke[0] - rke[1]) / 3.0;
        let valid = (r1 - 1.0).abs() < 0.01;
        println!(
            "    d={:.1}: R₁={:.5} [{}], R₂={:.5}, R_K={:.5} ({} s)",
            d,
            r1,
            if valid {
                "有効"
            } else {
                "除外 (核床カナリア)"
            },
            r2,
            rk,
            t0.elapsed().as_secs()
        );
        rows.push((d, r1, r2, rk, valid));
    }
    check(
        "有効点が 2 点以上 (掃引が成立)",
        rows.iter().filter(|r| r.4).count() >= 2,
        format!("{}/4 点が有効", rows.iter().filter(|r| r.4).count()),
    );

    // ---- 判定 (記録) ----
    let valid_rows: Vec<&(f64, f64, f64, f64, bool)> = rows.iter().filter(|r| r.4).collect();
    let mut mono = true;
    for i in 1..valid_rows.len() {
        if valid_rows[i].2 < valid_rows[i - 1].2 - 0.01 {
            mono = false;
        }
    }
    let last_ok = valid_rows.last().map(|r| r.2 > 0.90).unwrap_or(false);
    println!(
        "\n[判定] R₂(d) の壁距離依存 (有効 {} 点):",
        valid_rows.len()
    );
    for r in &valid_rows {
        println!("    d={:.1}: R₂ = {:.4}", r.0, r.2);
    }
    println!(
        "    => {}",
        if mono && last_ok {
            "事前登録 (a): 逸脱は壁近傍の格子非局所性 — ブースト同定は壁から離れて回復 (第二段の機構が立った)"
        } else {
            "事前登録 (b): 局所ブースト形の不整合を記録 (プラトー/非単調)"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v19.4".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(d, r1, r2, rk, va)| {
                        Json::Obj(vec![
                            ("d".into(), Json::Num(d)),
                            ("r1".into(), Json::Num(r1)),
                            ("r2".into(), Json::Num(r2)),
                            ("rk".into(), Json::Num(rk)),
                            ("valid".into(), Json::Bool(va)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(mono && last_ok)),
    ]);
    let p = write_artifact("results/v194_boostdist.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
