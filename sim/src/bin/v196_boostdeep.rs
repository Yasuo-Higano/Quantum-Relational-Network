//! v19.6 粒子・正孔摂動によるブースト同定の深部掃引 — 第二段の器械 2 号
//!
//! v19.4 は Gauss 波束対の回転で壁距離 d ≥ 2.5 の δS 信号が消える (両波束の占有内容が
//! 縮退) ことを示した。本版は摂動を分枝分割に替える: 複素波束 φ = G e^{iπx/4} の
//! 非占有部 u = (1−C₀)φ / 占有部 w = C₀φ の回転は占有内容の差が最大で、深部でも
//! δS が消えない (v19.5 の QNEC 器械で較正済みの機構)。しかもコヒーレンス項の
//! ξ 重みつき応答 δ⟨K_boost⟩ は d と共に成長する。モジュラー核は一切使わない
//! (v19.4 の核カナリアの代替 = 事前登録した有効性判定: 信号量 + Richardson 線形性)。
//! 閉形式 ΔC(α) = sin²α (u*uᵀ − w*wᵀ) + sinα cosα (u*wᵀ + w*uᵀ) — rank-2 厳密・
//! 純粋性は構成的に保たれる (最浅点で S(A)=S(Aᶜ) ゲート)。
//! 事前登録: (a) 有効点で |R₂−1| 単調減少かつ最遠有効点 < 0.05 = 第二段成立 /
//! (a′) 単調接近だが未達 = 傾向の記録 / (b) プラトー・非単調 = 局所形不整合の記録。

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

fn main() {
    self_test();
    println!("=== v19.6 粒子・正孔摂動によるブースト同定の深部掃引 (第二段の器械 2 号) ===\n");
    println!(
        "事前登録: v19.4 の深部信号消失を分枝分割摂動 (u = 波束の非占有部, w = 占有部) で回避し、"
    );
    println!("  核なしブースト比 R₂ = δS/δ⟨K_boost⟩ を壁距離 d ∈ {{0.5,1.5,2.5,3.5,4.5}} (N=12) で測る。");
    println!("  点の有効性 (事前登録): |δS_odd(α)| > 1e-9 かつ Richardson 線形性 |R₂(α)−R₂(2α)| < 0.1|R₂|。");
    println!("  (a) 有効点で単調接近かつ最遠有効点 |R₂−1| < 0.05 = ブースト同定が壁から離れて完成 (第二段成立) /");
    println!(
        "  (a′) 単調接近だが 5% に未達 = 傾向の記録 / (b) プラトー・非単調 = 局所形不整合の記録\n"
    );
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
    let (eps, v) = jacobi_eigh(&h, ns);
    let nocc = ns / 2;
    let gap = eps[nocc] - eps[nocc - 1];
    check(
        &format!("N={} の閉殻ギャップ", n),
        gap > 1e-6,
        format!("gap = {:.4} ({} s)", gap, t0.elapsed().as_secs()),
    );
    let mut c0 = vec![0.0f64; ns * ns];
    for k in 0..nocc {
        for i in 0..ns {
            let vi = v[i + k * ns];
            if vi == 0.0 {
                continue;
            }
            for j in 0..ns {
                c0[j + i * ns] += vi * v[j + k * ns];
            }
        }
    }
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut sel = Vec::new();
    let mut sel_c = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n {
                if x < n / 2 {
                    sel.push(idx(x, y, z));
                } else {
                    sel_c.push(idx(x, y, z));
                }
            }
        }
    }
    let m = sel.len();
    // ブースト重みつきボンド (規約は v19.2 から凍結)
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
    // S₀ (α=0) は d に依らない真空値 — 1 回だけ実対称で評価
    {
        let ca = restrict(&c0, ns, &sel);
        println!(
            "    S₀(A) = {:.6} ({} s)",
            entropy_sym(&ca, m),
            t0.elapsed().as_secs()
        );
    }
    let sig = 1.0f64;
    let kx0 = std::f64::consts::PI / 4.0;
    let mid = n as f64 / 2.0;
    let dists = [0.5f64, 1.5, 2.5, 3.5, 4.5];
    // (d, r2, valid, 理由)
    let mut rows: Vec<(f64, f64, bool, String)> = Vec::new();
    for (di, &d) in dists.iter().enumerate() {
        // 複素波束 → 分枝分割
        let (xc, yc, zc) = (xplane - d, mid, mid);
        let mut phr = vec![0.0f64; ns];
        let mut phi = vec![0.0f64; ns];
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let d2 =
                        (x as f64 - xc).powi(2) + (y as f64 - yc).powi(2) + (z as f64 - zc).powi(2);
                    let g = (-d2 / (2.0 * sig * sig)).exp();
                    phr[idx(x, y, z)] = g * (kx0 * x as f64).cos();
                    phi[idx(x, y, z)] = g * (kx0 * x as f64).sin();
                }
            }
        }
        let mut ptr = vec![0.0f64; ns];
        let mut pti = vec![0.0f64; ns];
        for k in 0..ns {
            let mut ar = 0.0;
            let mut ai = 0.0;
            for i in 0..ns {
                let vv = v[i + k * ns];
                ar += vv * phr[i];
                ai += vv * phi[i];
            }
            ptr[k] = ar;
            pti[k] = ai;
        }
        let nrm_u: f64 = (nocc..ns).map(|k| ptr[k] * ptr[k] + pti[k] * pti[k]).sum();
        let nrm_w: f64 = (0..nocc).map(|k| ptr[k] * ptr[k] + pti[k] * pti[k]).sum();
        let ratio = nrm_u / (nrm_u + nrm_w);
        check(
            &format!("d={:.1} の占有/非占有分割が非退化", d),
            ratio > 0.15 && ratio < 0.85,
            format!("非占有比 = {:.3}", ratio),
        );
        // サイト基底の u, w (正規化)
        let (su, sw) = (nrm_u.sqrt(), nrm_w.sqrt());
        let mut ur = vec![0.0f64; ns];
        let mut ui = vec![0.0f64; ns];
        let mut wr = vec![0.0f64; ns];
        let mut wi = vec![0.0f64; ns];
        for k in 0..ns {
            let (a, b) = (ptr[k], pti[k]);
            if a == 0.0 && b == 0.0 {
                continue;
            }
            if k >= nocc {
                for i in 0..ns {
                    let vv = v[i + k * ns];
                    ur[i] += vv * a / su;
                    ui[i] += vv * b / su;
                }
            } else {
                for i in 0..ns {
                    let vv = v[i + k * ns];
                    wr[i] += vv * a / sw;
                    wi[i] += vv * b / sw;
                }
            }
        }
        // ΔC(α) = sin²α (uu† − ww†) + sinα cosα (uw† + wu†)
        let dc = |alpha: f64, i: usize, j: usize| -> (f64, f64) {
            let (s, c) = (alpha.sin(), alpha.cos());
            let s2 = s * s;
            let sc = s * c;
            // uu†: re = ur_i ur_j + ui_i ui_j, im = ui_i ur_j − ur_i ui_j? — C = x*xᵀ 規約:
            // (u*uᵀ)_{ij} = u*_i u_j: re = ur_i ur_j + ui_i ui_j, im = ur_i ui_j − ui_i ur_j
            let uu_re = ur[i] * ur[j] + ui[i] * ui[j];
            let uu_im = ur[i] * ui[j] - ui[i] * ur[j];
            let ww_re = wr[i] * wr[j] + wi[i] * wi[j];
            let ww_im = wr[i] * wi[j] - wi[i] * wr[j];
            let uw_re = ur[i] * wr[j] + ui[i] * wi[j] + wr[i] * ur[j] + wi[i] * ui[j];
            let uw_im = ur[i] * wi[j] - ui[i] * wr[j] + wr[i] * ui[j] - wi[i] * ur[j];
            (
                s2 * (uu_re - ww_re) + sc * uw_re,
                s2 * (uu_im - ww_im) + sc * uw_im,
            )
        };
        // S(α; sel) と δK_boost(α)
        let s_of = |alpha: f64, ss: &[usize]| -> f64 {
            let mm = ss.len();
            let mut cre = vec![0.0f64; mm * mm];
            let mut cim = vec![0.0f64; mm * mm];
            for a in 0..mm {
                for b in 0..mm {
                    let (re, im) = dc(alpha, ss[a], ss[b]);
                    cre[a + b * mm] = c0[ss[a] + ss[b] * ns] + re;
                    cim[a + b * mm] = im;
                }
            }
            entropy_corr_herm(&cre, &cim, mm)
        };
        let kb_of = |alpha: f64| -> f64 {
            let mut acc = 0.0;
            for &(i, j, t, xi) in &bonds {
                let (re, _) = dc(alpha, i, j);
                acc += two_pi * xi * t * 2.0 * re;
            }
            acc
        };
        // 純粋性+機構ゲート (最浅点のみ): S(A) = S(Aᶜ)
        if di == 0 {
            let sa = s_of(0.04, &sel);
            let sc_ = s_of(0.04, &sel_c);
            check(
                "d=0.5 の補集合一致 (純粋性+機構)",
                (sa - sc_).abs() < 1e-6,
                format!("|ΔS| = {:.2e}", (sa - sc_).abs()),
            );
        }
        let a1 = 0.02f64;
        let mut r2e = [0.0f64; 2];
        let mut ds_min = f64::INFINITY;
        for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
            let ds = (s_of(al, &sel) - s_of(-al, &sel)) / 2.0;
            let dkb = (kb_of(al) - kb_of(-al)) / 2.0;
            ds_min = ds_min.min(ds.abs());
            r2e[ii] = ds / dkb;
        }
        let r2 = (4.0 * r2e[0] - r2e[1]) / 3.0;
        let lin = (r2e[0] - r2e[1]).abs();
        let valid = ds_min > 1e-9 && lin < 0.1 * r2.abs();
        let reason = if valid {
            "有効".to_string()
        } else if ds_min <= 1e-9 {
            "除外 (信号不足)".to_string()
        } else {
            format!("除外 (非線形 {:.3})", lin)
        };
        println!(
            "    d={:.1}: R₂={:.5} [R₂(α)={:.5}, R₂(2α)={:.5}] |δS|min={:.2e} [{}] ({} s)",
            d,
            r2,
            r2e[0],
            r2e[1],
            ds_min,
            reason,
            t0.elapsed().as_secs()
        );
        rows.push((d, r2, valid, reason));
    }
    check(
        "有効点が 3 点以上 (掃引が成立)",
        rows.iter().filter(|r| r.2).count() >= 3,
        format!("{}/5 点が有効", rows.iter().filter(|r| r.2).count()),
    );

    // ---- 判定 (記録) ----
    let valid_rows: Vec<&(f64, f64, bool, String)> = rows.iter().filter(|r| r.2).collect();
    let mut mono = true;
    for i in 1..valid_rows.len() {
        if (valid_rows[i].1 - 1.0).abs() > (valid_rows[i - 1].1 - 1.0).abs() + 0.01 {
            mono = false;
        }
    }
    let last_dev = valid_rows.last().map(|r| (r.1 - 1.0).abs()).unwrap_or(9.9);
    println!(
        "\n[判定] R₂(d) (粒子・正孔摂動, 有効 {} 点):",
        valid_rows.len()
    );
    for r in &valid_rows {
        println!(
            "    d={:.1}: R₂ = {:.4} (|R₂−1| = {:.4})",
            r.0,
            r.1,
            (r.1 - 1.0).abs()
        );
    }
    println!(
        "    => {}",
        if mono && last_dev < 0.05 {
            "事前登録 (a): ブースト同定が壁から離れて完成 — K_A = 2π Σ ξ T₀₀ (第二段成立)"
        } else if mono {
            "事前登録 (a′): 単調接近だが 5% に未達 — 傾向の記録"
        } else {
            "事前登録 (b): プラトー・非単調 — 局所ブースト形の不整合を記録"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v19.6".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|(d, r2, va, _)| {
                        Json::Obj(vec![
                            ("d".into(), Json::Num(*d)),
                            ("r2".into(), Json::Num(*r2)),
                            ("valid".into(), Json::Bool(*va)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(mono && last_dev < 0.05)),
    ]);
    let p = write_artifact("results/v196_boostdeep.json", &j.render());
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
