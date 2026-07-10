//! v19.3 相対エントロピーの三性質 — 正値・二次 (Fisher)・単調 (Einstein 橋の第三段)
//!
//! v2.1 (1+1D) の 3+1D 版: S_rel(ρ_α‖ρ₀) = Δ⟨K_A⟩ − ΔS の
//!   [1] 正値性 (厳密な数学的性質 — 破れたら実装誤り = 装置ゲート)
//!   [2] 二次性 S_rel = (F/2)α² (Fisher 計量 = 正準エネルギーの情報側) —
//!       偶部推定器 [S_rel(α)+S_rel(−α)]/2 が核床の奇数次汚染を厳密に消す
//!   [3] 入れ子単調性 A₂ ⊂ A₁ ⊂ A₀ ⟹ S_rel(A₂) ≤ S_rel(A₁) ≤ S_rel(A₀) (厳密)
//! を 3+1D 格子 Rindler (v19.1 の構成) で検証する。
//! 事前登録: [1][3] はゲート、[2] は局所傾き p(α→2α) の最小対 p₁ = p(0.01→0.02) が
//! 4 点全て 2 ± 0.05 で (a) Fisher 成立 — |p−2| の α 依存で四次尾の制御性を記録。
//! F を Richardson 抽出し N 依存 (物理幅固定波束) を記録 — 第四段 (QNEC)・
//! 第五段 (非線形) の規格化点になる。

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
        "=== v19.3 相対エントロピーの三性質: 正値・二次 (Fisher)・単調 (Einstein 橋の第三段) ===\n"
    );
    println!("事前登録: 正値性・入れ子単調性は厳密な数学的性質 = 装置ゲート (破れたら実装誤り)。");
    println!("          二次性は局所傾き p(α→2α) で測る — 最小対 p₁ = p(0.01→0.02) が 4 点全て");
    println!("          2 ± 0.05 = (a) Fisher 成立 (|p−2| の α 依存で四次尾の制御性を記録) / 外れ = (b) 記録。");
    println!("          F (Fisher 計量) を Richardson 抽出し N 依存を記録\n");
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

    // (n, pair, [p12, p23, p34], F)
    let mut fisher_rows: Vec<(usize, usize, [f64; 3], f64)> = Vec::new();
    for &n in &[8usize, 12] {
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
        let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
        // 入れ子領域 A_k = {x < n/2 − k}, k = 0, 1, 2
        let region = |k: usize| -> Vec<usize> {
            let mut sel = Vec::new();
            for z in 0..n {
                for y in 0..n {
                    for x in 0..n / 2 - k {
                        sel.push(idx(x, y, z));
                    }
                }
            }
            sel
        };
        // 波束対 (v19.1/19.2 と同系, 物理幅 wid = N/8)
        let wid = n as f64 / 8.0;
        let make_packet = |x0: f64, y0: f64, z0: f64| -> Vec<f64> {
            let mut p = vec![0.0f64; ns];
            for x in 0..n {
                for y in 0..n {
                    for z in 0..n {
                        let d2 = (x as f64 - x0).powi(2)
                            + (y as f64 - y0).powi(2)
                            + (z as f64 - z0).powi(2);
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
        let xplane = n as f64 / 2.0 - 0.5;
        let mid = n as f64 / 2.0;
        let centers = [
            (
                (xplane - wid / 2.0, mid, mid),
                (xplane + wid / 2.0, mid, mid),
            ),
            ((xplane - wid, mid - wid, mid), (xplane, mid - wid, mid)),
        ];
        // 領域ごとの (S₀, K 固有系) を用意
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
            let sel = region(k);
            let m = sel.len();
            let ca0 = restrict(&c, ns, &sel);
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
        for (pi, &((ax, ay, az), (bx, by, bz))) in centers.iter().enumerate() {
            let u = make_packet(ax, ay, az);
            let mut wv = make_packet(bx, by, bz);
            orth(&u, &mut wv);
            // S_rel(α; A_k)
            let srel = |alpha: f64, rg: &Reg| -> f64 {
                let cc = rotate_c(&c, ns, &u, &wv, alpha);
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
                            acc += ui
                                * (ca[j + i * rg.m] - rg.ca0[j + i * rg.m])
                                * rg.cv[j + k * rg.m];
                        }
                    }
                    dk += rg.kappa[k] * acc;
                }
                dk - (s_a - rg.s0)
            };
            // 偶部 S_even(α) = [S_rel(α)+S_rel(−α)]/2 を α ∈ {0.01,0.02,0.04,0.08} で測定。
            // 同じ評価から [1] 正値性 (8 点全て、厳密) も判定する。
            let als = [0.01f64, 0.02, 0.04, 0.08];
            let mut evens = Vec::new();
            let mut min_sr = f64::INFINITY;
            let mut sr08 = 0.0;
            for &al in &als {
                let sp = srel(al, &regs[0]);
                let sm = srel(-al, &regs[0]);
                min_sr = min_sr.min(sp).min(sm);
                if al == 0.08 {
                    sr08 = sp;
                }
                evens.push((sp + sm) / 2.0);
            }
            check(
                &format!("N={} 対{} の S_rel 正値性 (厳密)", n, pi + 1),
                min_sr > 0.0,
                format!("最小 S_rel = {:.3e} (α = ±0.01..±0.08 の 8 点)", min_sr),
            );
            // [2] 二次性: 局所傾き p(α→2α) = ln(S_even(2α)/S_even(α))/ln 2
            let ln2 = std::f64::consts::LN_2;
            let p_loc = [
                (evens[1] / evens[0]).ln() / ln2,
                (evens[2] / evens[1]).ln() / ln2,
                (evens[3] / evens[2]).ln() / ln2,
            ];
            // F の Richardson (α² 主導項): F(α) = 2 S_even/α² → (4F(0.01)−F(0.02))/3
            let f1 = 2.0 * evens[0] / (als[0] * als[0]);
            let f2 = 2.0 * evens[1] / (als[1] * als[1]);
            let f_rich = (4.0 * f1 - f2) / 3.0;
            println!(
                "    N={} 対{}: 局所傾き p = [{:.4}, {:.4}, {:.4}] (0.01→0.02→0.04→0.08), F = {:.5}",
                n,
                pi + 1,
                p_loc[0],
                p_loc[1],
                p_loc[2],
                f_rich
            );
            println!(
                "      |p−2| 列 = [{:.4}, {:.4}, {:.4}] (四次尾なら α² で縮む)",
                (p_loc[0] - 2.0).abs(),
                (p_loc[1] - 2.0).abs(),
                (p_loc[2] - 2.0).abs()
            );
            fisher_rows.push((n, pi, p_loc, f_rich));
            // [3] 入れ子単調性 (厳密): S_rel(A₀) ≥ S_rel(A₁) ≥ S_rel(A₂) at α = 0.08
            let sr_k = [sr08, srel(0.08, &regs[1]), srel(0.08, &regs[2])];
            let mono = sr_k[0] >= sr_k[1] - 1e-12 && sr_k[1] >= sr_k[2] - 1e-12;
            check(
                &format!("N={} 対{} の入れ子単調性 (厳密)", n, pi + 1),
                mono,
                format!("{:.3e} ≥ {:.3e} ≥ {:.3e}", sr_k[0], sr_k[1], sr_k[2]),
            );
        }
    }

    // ---- 判定 (記録) ----
    let all_quad = fisher_rows
        .iter()
        .all(|&(_, _, p, _)| (p[0] - 2.0).abs() < 0.05);
    println!(
        "\n[判定] {}",
        if all_quad {
            "事前登録 (a): 三性質 (正値・二次・単調) が 3+1D で成立 — Fisher 計量が立った (第三段成立)"
        } else {
            "事前登録 (b): 最小 α 対でも二次帯の外 — 逸脱の記録"
        }
    );
    println!("    F の N 依存 (物理幅固定):");
    for &(n, p, _, f) in &fisher_rows {
        println!("      N={} 対{}: F = {:.5}", n, p + 1, f);
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v19.3".into())),
        (
            "fisher".into(),
            Json::Arr(
                fisher_rows
                    .iter()
                    .map(|&(n, p, pl, f)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("pair".into(), Json::Int(p as i64)),
                            ("p12".into(), Json::Num(pl[0])),
                            ("p23".into(), Json::Num(pl[1])),
                            ("p34".into(), Json::Num(pl[2])),
                            ("f".into(), Json::Num(f)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(all_quad)),
    ]);
    let p = write_artifact("results/v193_fisher3d.json", &j.render());
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
