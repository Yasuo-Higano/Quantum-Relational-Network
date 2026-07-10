//! v19.7 有限振幅の S_rel — 非線形不等式と二次近似の崩壊点 (Einstein 橋の第五段)
//!
//! 第三段 (v19.3) は α→0 の二次構造 (Fisher) を立てた。本版は有限振幅 α ∈
//! {0.2, 0.5, 1.0, π/2} で S_rel(α) = Δ⟨K_A⟩ − ΔS を厳密評価し、
//!   [1] 正値性 (厳密ゲート)・[2] 入れ子単調性 (α=1.0, 厳密ゲート) —
//! 非線形領域でも情報側の不等式が立つこと — と、
//!   [3] 比 r(α) = S_rel/((F/2)α²) の二次近似崩壊点 α* と四次構造の符号 (記録)
//! を測る。摂動族は 2 系: 族 A = 壁波束対の実直交回転 (v19.3 系)、族 B = 分枝分割
//! 粒子・正孔 (v19.5/19.6 系, 複素・閉形式 ΔC)。族 B の α=π/2 は 1 量子の完全励起 —
//! v0.7 の「1 量子は無限小でない」の 3+1D 定量化。
//! 事前登録: (a) ゲート全通過かつ四次符号が族間一致 = 第五段の到達可能形成立 /
//! (b) 族間不一致 = 非普遍の記録。「非線形」= 固定格子上の有限振幅状態空間の
//! 不等式であり、動的重力 (背反作用) は自由格子の外 (限界として明示)。

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
        "=== v19.7 有限振幅の S_rel: 非線形不等式と二次近似の崩壊点 (Einstein 橋の第五段) ===\n"
    );
    println!("事前登録: 正値性 (全 α)・入れ子単調性 (α=1.0) は厳密な数学的性質 = 装置ゲート。");
    println!("  記録: 比 r(α) = S_rel/((F/2)α²) の崩壊点 α* (|r−1| > 0.25 の最小 α) と");
    println!("  四次構造の符号 sign(r(0.2)−1)。分岐: (a) ゲート全通過かつ四次符号が 2 摂動族で");
    println!("  一致 = 第五段の到達可能形成立 (非線形構造の族間普遍性) / (b) 族間不一致 = 記録。");
    println!("  注: ここでの「非線形」は固定格子上の有限振幅状態空間の不等式 — 動的重力");
    println!("  (背反作用) は自由格子の外 (限界として明示)。族 B の α=π/2 は 1 量子の完全");
    println!("  励起 — v0.7 の「1 量子は無限小でない」が 3+1D の測定量になる\n");
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
    // 入れ子領域 A_k = {x < n/2 − k} と各領域のモジュラー核固有系
    struct Reg {
        sel: Vec<usize>,
        m: usize,
        ca0: Vec<f64>,
        s0: f64,
        cv: Vec<f64>,
        kappa: Vec<f64>,
    }
    let mut regs = Vec::new();
    for k in 0..3usize {
        let mut sel = Vec::new();
        for z in 0..n {
            for y in 0..n {
                for x in 0..n / 2 - k {
                    sel.push(idx(x, y, z));
                }
            }
        }
        let m = sel.len();
        let ca0 = restrict(&c0, ns, &sel);
        let s0 = entropy_sym(&ca0, m);
        let (cw_ev, cv) = jacobi_eigh(&ca0, m);
        let kappa: Vec<f64> = cw_ev
            .iter()
            .map(|&ck| {
                let x = ck.clamp(1e-14, 1.0 - 1e-14);
                ((1.0 - x) / x).ln()
            })
            .collect();
        regs.push(Reg {
            sel,
            m,
            ca0,
            s0,
            cv,
            kappa,
        });
    }
    let xplane = n as f64 / 2.0 - 0.5;
    let mid = n as f64 / 2.0;
    let wid = n as f64 / 8.0;
    // ---- 族 A: 壁波束対の実直交回転 (v19.3 と同系) ----
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
    let ua = make_packet(xplane - wid, mid - wid, mid);
    let mut wa = make_packet(xplane, mid - wid, mid);
    {
        let ov: f64 = ua.iter().zip(wa.iter()).map(|(a, b)| a * b).sum();
        for i in 0..ns {
            wa[i] -= ov * ua[i];
        }
        let nrm: f64 = wa.iter().map(|a| a * a).sum::<f64>().sqrt();
        for a in wa.iter_mut() {
            *a /= nrm;
        }
    }
    // ---- 族 B: 分枝分割 (壁の複素波束 e^{iπx/4}) ----
    let kx0 = std::f64::consts::PI / 4.0;
    let mut phr = vec![0.0f64; ns];
    let mut phi = vec![0.0f64; ns];
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let d2 = (x as f64 - xplane).powi(2)
                    + (y as f64 - mid).powi(2)
                    + (z as f64 - mid).powi(2);
                let g = (-d2 / (2.0 * wid * wid)).exp();
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
    check(
        "族 B の占有/非占有分割が非退化",
        nrm_u / (nrm_u + nrm_w) > 0.15 && nrm_u / (nrm_u + nrm_w) < 0.85,
        format!("非占有比 = {:.3}", nrm_u / (nrm_u + nrm_w)),
    );
    let (su, sw) = (nrm_u.sqrt(), nrm_w.sqrt());
    let mut ub_r = vec![0.0f64; ns];
    let mut ub_i = vec![0.0f64; ns];
    let mut wb_r = vec![0.0f64; ns];
    let mut wb_i = vec![0.0f64; ns];
    for k in 0..ns {
        let (a, b) = (ptr[k], pti[k]);
        if a == 0.0 && b == 0.0 {
            continue;
        }
        if k >= nocc {
            for i in 0..ns {
                let vv = v[i + k * ns];
                ub_r[i] += vv * a / su;
                ub_i[i] += vv * b / su;
            }
        } else {
            for i in 0..ns {
                let vv = v[i + k * ns];
                wb_r[i] += vv * a / sw;
                wb_i[i] += vv * b / sw;
            }
        }
    }
    // ΔC 閉形式。族 A (実): ΔC = (cosα−1)(uuᵀ+wwᵀ)C₀(...) → rotate_c と同値の
    // rank-2 形を直接使う: ΔC_A = OC₀Oᵀ − C₀。ここでは C(α) 全体を作らず、
    // 制限блок上で v19.3 同様に rotate_c (body) を用いる。
    // 族 B (複素): ΔC = sin²α (u*uᵀ − w*wᵀ) + sinα cosα (u*wᵀ + w*uᵀ)。
    let dcb = |alpha: f64, i: usize, j: usize| -> (f64, f64) {
        let (s, c) = (alpha.sin(), alpha.cos());
        let s2 = s * s;
        let sc = s * c;
        let uu_re = ub_r[i] * ub_r[j] + ub_i[i] * ub_i[j];
        let uu_im = ub_r[i] * ub_i[j] - ub_i[i] * ub_r[j];
        let ww_re = wb_r[i] * wb_r[j] + wb_i[i] * wb_i[j];
        let ww_im = wb_r[i] * wb_i[j] - wb_i[i] * wb_r[j];
        let uw_re = ub_r[i] * wb_r[j] + ub_i[i] * wb_i[j] + wb_r[i] * ub_r[j] + wb_i[i] * ub_i[j];
        let uw_im = ub_r[i] * wb_i[j] - ub_i[i] * wb_r[j] + wb_r[i] * ub_i[j] - wb_i[i] * ub_r[j];
        (
            s2 * (uu_re - ww_re) + sc * uw_re,
            s2 * (uu_im - ww_im) + sc * uw_im,
        )
    };
    // S_rel(α; 領域) — 族 A は実 (rotate_c + entropy_sym)、族 B は複素 (閉形式 + herm)
    let srel_a = |alpha: f64, rg: &Reg| -> f64 {
        let cc = rotate_c(&c0, ns, &ua, &wa, alpha);
        let ca = restrict(&cc, ns, &rg.sel);
        let s_a = entropy_sym(&ca, rg.m);
        let mut dk = 0.0;
        for k in 0..rg.m {
            let mut acc = 0.0;
            for i in 0..rg.m {
                let ui = rg.cv[i + k * rg.m];
                if ui == 0.0 {
                    continue;
                }
                for j in 0..rg.m {
                    acc += ui * (ca[j + i * rg.m] - rg.ca0[j + i * rg.m]) * rg.cv[j + k * rg.m];
                }
            }
            dk += rg.kappa[k] * acc;
        }
        dk - (s_a - rg.s0)
    };
    let srel_b = |alpha: f64, rg: &Reg| -> f64 {
        let mm = rg.m;
        let mut cre = vec![0.0f64; mm * mm];
        let mut cim = vec![0.0f64; mm * mm];
        for a in 0..mm {
            for b in 0..mm {
                let (re, im) = dcb(alpha, rg.sel[a], rg.sel[b]);
                cre[a + b * mm] = c0[rg.sel[a] + rg.sel[b] * ns] + re;
                cim[a + b * mm] = im;
            }
        }
        let s_a = entropy_corr_herm(&cre, &cim, mm);
        // Δ⟨K⟩ = Σ κ_k ⟨v_k|ΔC|v_k⟩ — v_k 実なので ΔC_im は消える
        let mut dk = 0.0;
        for k in 0..mm {
            let mut acc = 0.0;
            for i in 0..mm {
                let ui = rg.cv[i + k * mm];
                if ui == 0.0 {
                    continue;
                }
                for j in 0..mm {
                    acc += ui * (cre[i + j * mm] - rg.ca0[i + j * mm]) * rg.cv[j + k * mm];
                }
            }
            dk += rg.kappa[k] * acc;
        }
        dk - (s_a - rg.s0)
    };
    // ---- 測定 ----
    let alphas = [0.2f64, 0.5, 1.0, std::f64::consts::FRAC_PI_2];
    let mut fam_rows: Vec<(String, Vec<f64>, f64, f64, f64)> = Vec::new(); // (族, r(α)列, F, α*, 四次符号)
    for fam in 0..2usize {
        let name = if fam == 0 {
            "族A 壁波束対"
        } else {
            "族B 粒子・正孔"
        };
        let sr = |al: f64, rg: &Reg| -> f64 {
            if fam == 0 {
                srel_a(al, rg)
            } else {
                srel_b(al, rg)
            }
        };
        // 小 α アンカー: F = 2·S_even(0.02)/0.02²
        let a0 = 0.02f64;
        let se = (sr(a0, &regs[0]) + sr(-a0, &regs[0])) / 2.0;
        let f_est = 2.0 * se / (a0 * a0);
        println!(
            "    {}: F(アンカー α=0.02) = {:.5} ({} s)",
            name,
            f_est,
            t0.elapsed().as_secs()
        );
        let mut ratios = Vec::new();
        let mut all_pos = true;
        let mut min_sr = f64::INFINITY;
        for &al in &alphas {
            let s = sr(al, &regs[0]);
            min_sr = min_sr.min(s);
            if s <= 0.0 {
                all_pos = false;
            }
            let r = s / (0.5 * f_est * al * al);
            ratios.push(r);
            println!(
                "    {} α={:.4}: S_rel = {:+.6} r = S_rel/((F/2)α²) = {:.4} ({} s)",
                name,
                al,
                s,
                r,
                t0.elapsed().as_secs()
            );
        }
        check(
            &format!("{} の正値性 (厳密, 全 α)", name),
            all_pos && se > 0.0,
            format!("最小 S_rel = {:.3e}", min_sr.min(se)),
        );
        // 入れ子単調性 (α = 1.0)
        let s_k: Vec<f64> = (0..3).map(|k| sr(1.0, &regs[k])).collect();
        check(
            &format!("{} の入れ子単調性 (厳密, α=1.0)", name),
            s_k[0] >= s_k[1] - 1e-12 && s_k[1] >= s_k[2] - 1e-12,
            format!("{:.4} ≥ {:.4} ≥ {:.4}", s_k[0], s_k[1], s_k[2]),
        );
        let astar = alphas
            .iter()
            .zip(ratios.iter())
            .find(|(_, r)| (**r - 1.0).abs() > 0.25)
            .map(|(a, _)| *a)
            .unwrap_or(f64::INFINITY);
        let q_sign = (ratios[0] - 1.0).signum();
        println!(
            "    {}: α* = {} (|r−1|>0.25 の最小 α), 四次符号 = {:+.0}",
            name,
            if astar.is_finite() {
                format!("{:.4}", astar)
            } else {
                "> π/2".to_string()
            },
            q_sign
        );
        fam_rows.push((name.to_string(), ratios, f_est, astar, q_sign));
    }

    // ---- 判定 (記録) ----
    let sign_match = fam_rows[0].4 == fam_rows[1].4;
    println!(
        "\n[判定] {}",
        if nfail == 0 && sign_match {
            "事前登録 (a): 有限振幅の非線形不等式が全て成立し四次構造の符号が族間で一致 — 第五段の到達可能形成立"
        } else if nfail == 0 {
            "事前登録 (b): 不等式は成立・四次構造は族間で不一致 — 非普遍の記録"
        } else {
            "装置ゲート破れ — 判定無効"
        }
    );
    println!(
        "    族A: α* = {:.4}, 符号 {:+.0} / 族B: α* = {:.4}, 符号 {:+.0}",
        fam_rows[0].3, fam_rows[0].4, fam_rows[1].3, fam_rows[1].4
    );
    println!("    族B α=π/2 (1 量子の完全励起): S_rel = 有限 — v0.7 の「1 量子は無限小でない」の 3+1D 定量化");

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v19.7".into())),
        (
            "families".into(),
            Json::Arr(
                fam_rows
                    .iter()
                    .map(|(nm, rs, f, astar, sg)| {
                        Json::Obj(vec![
                            ("name".into(), Json::Str(nm.clone())),
                            (
                                "ratios".into(),
                                Json::Arr(rs.iter().map(|&x| Json::Num(x)).collect()),
                            ),
                            ("f".into(), Json::Num(*f)),
                            (
                                "astar".into(),
                                if astar.is_finite() {
                                    Json::Num(*astar)
                                } else {
                                    Json::Num(-1.0)
                                },
                            ),
                            ("qsign".into(), Json::Num(*sg)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(sign_match)),
    ]);
    let p = write_artifact("results/v197_nonlinear.json", &j.render());
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
