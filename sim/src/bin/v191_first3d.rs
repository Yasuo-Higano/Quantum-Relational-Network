//! v19.1 3+1D の第一法則 — 格子 Rindler の δS = δ⟨K⟩ (Einstein 橋の第一段)
//!
//! PROMPT/3 v20 の梯子 (QRN 状態 → モジュラー → QNEC → T_μν → Einstein → 3+1D) の
//! 第一段を 3+1D で立てる: エンタングルメント第一法則 δS = δ⟨K_A⟩。これは v0.7
//! (1+1D) の 3+1D 版であり、線形化 Einstein 方程式の情報理論的な芯である
//! (Jacobson / Faulkner ら — 半空間のモジュラー・ハミルトニアンはブースト)。
//!
//! 構成 (全て実対称): 3+1D staggered fermion (x 開放 = 単一の絡み面, y/z 反周期)、
//! 半充填閉殻の基底射影 C、領域 A = {x < N/2}、K_A は C_A の固有基底の
//! κ = ln((1−c)/c)。第一法則: R = δS_odd / δ⟨K_A⟩_odd → 1 (Richardson で α² 消去)。
//!
//! 事前登録: R = 1 ± 1% を「波束対 2 種 × 格子 N ∈ {8, 12}」の 4 点全てで
//! 満たせば (a) 第一段成立 / 系統的にずれれば (b) 破れの記録 — どちらも正当な
//! 結果であり、装置ゲート (閉殻・C²=C・面積則比) のみが PASS/FAIL。
//!
//! 装置改訂 (初走の教訓 — 開発記録): 拡がった Fermi 固有モード対の回転は、真の
//! κ ~ 100 の深部モードに重みを持ち、モジュラー核の f64 分解能床 (c の 1e-14
//! クランプ → κ ≤ ln(1e14) ≈ 32.2 に飽和) が α 非依存の定数バイアス b を作る
//! (初走: N=8 で b≈1.0%, N=12 で b≈5.4% — R = 1 + b + cα² 分解で同定)。本走は
//! 主摂動を「壁局在の波束対の直交回転」(任意の反対称生成子 J で C→OCOᵀ は厳密に
//! 射影子) に変更 — 深部モードへの重みが指数的に消え、核の床に触れない。
//! 固有モード対の旧摂動は診断として併記し b を明示する (モジュラー核の f64 床 =
//! 新しい数値の落とし穴として記録)。

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
    println!("=== v19.1 3+1D の第一法則: 格子 Rindler の δS = δ⟨K_A⟩ ===\n");
    println!("事前登録: R = 1 ± 1% (波束対 2 × N ∈ {{8,12}} の 4 点) = (a) 第一段成立 /");
    println!("          系統ずれ = (b) 破れの記録 — 分岐は記録であり装置ゲートのみ PASS/FAIL\n");
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

    let mut s_half = Vec::new();
    let mut ratios: Vec<(usize, usize, f64)> = Vec::new();
    let mut diag_bias: Vec<(usize, f64)> = Vec::new();
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
        let s0 = entropy_sym(&ca0, m);
        s_half.push((n, s0));
        println!("    N={}: S(A) = {:.4}", n, s0);
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
        // ---- 主摂動: 壁局在の波束対 ----
        let make_packet = |x0: f64, y0: f64, z0: f64, wid: f64| -> Vec<f64> {
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
        let xc = n as f64 / 2.0 - 0.5;
        let mid = n as f64 / 2.0;
        let centers = [
            ((xc - 0.5, mid, mid), (xc + 0.5, mid, mid)),
            ((xc - 1.0, mid - 1.0, mid), (xc, mid - 1.0, mid)),
        ];
        for (pi, &((ax, ay, az), (bx, by, bz))) in centers.iter().enumerate() {
            let u = make_packet(ax, ay, az, 1.2);
            let mut wv = make_packet(bx, by, bz, 1.2);
            orth(&u, &mut wv);
            let eval_at = |alpha: f64| -> (f64, f64) {
                let cc = rotate_c(&c, ns, &u, &wv, alpha);
                let ca = restrict(&cc, ns, &sel);
                (entropy_sym(&ca, m), dk_of(&ca))
            };
            let a1 = 0.02f64;
            let mut r_est = [0.0f64; 2];
            for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
                let (sp, kp) = eval_at(al);
                let (sm, km) = eval_at(-al);
                r_est[ii] = ((sp - sm) / 2.0) / ((kp - km) / 2.0);
            }
            let r = (4.0 * r_est[0] - r_est[1]) / 3.0;
            println!(
                "    N={} 波束対{}: R(α)={:.5}, R(2α)={:.5} → R = {:.5}",
                n,
                pi + 1,
                r_est[0],
                r_est[1],
                r
            );
            ratios.push((n, pi, r));
        }
        // ---- 診断: 固有モード対 (核床バイアス b の展示 — 主測定には不使用) ----
        {
            let wall_w = |k: usize| -> f64 {
                let mut wsum = 0.0;
                for z in 0..n {
                    for y in 0..n {
                        for x in (n / 2).saturating_sub(2)..(n / 2 + 2).min(n) {
                            let a = v[idx(x, y, z) + k * ns];
                            wsum += a * a;
                        }
                    }
                }
                wsum
            };
            let mut occs: Vec<usize> = (nocc.saturating_sub(20)..nocc).collect();
            let mut emps: Vec<usize> = (nocc..(nocc + 20).min(ns)).collect();
            occs.sort_by(|&a, &b| wall_w(b).partial_cmp(&wall_w(a)).unwrap());
            emps.sort_by(|&a, &b| wall_w(b).partial_cmp(&wall_w(a)).unwrap());
            let (ka, kb) = (occs[0], emps[0]);
            let ua: Vec<f64> = (0..ns).map(|i| v[i + ka * ns]).collect();
            let ub: Vec<f64> = (0..ns).map(|i| v[i + kb * ns]).collect();
            let eval_at = |alpha: f64| -> (f64, f64) {
                let cc = rotate_c(&c, ns, &ua, &ub, alpha);
                let ca = restrict(&cc, ns, &sel);
                (entropy_sym(&ca, m), dk_of(&ca))
            };
            let a1 = 0.02f64;
            let mut r_est = [0.0f64; 2];
            for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
                let (sp, kp) = eval_at(al);
                let (sm, km) = eval_at(-al);
                r_est[ii] = ((sp - sm) / 2.0) / ((kp - km) / 2.0);
            }
            let b = (4.0 * r_est[0] - r_est[1]) / 3.0 - 1.0;
            println!(
                "    N={} [診断] 固有モード対の核床バイアス b = {:+.4} (深部重みが f64 床に触れる)",
                n, b
            );
            diag_bias.push((n, b));
        }
    }
    let ratio_area = s_half[1].1 / s_half[0].1;
    check(
        "面積則スケール比 S(12)/S(8) ∈ [1.8, 2.8]",
        (1.8..=2.8).contains(&ratio_area),
        format!("{:.3} (期待 2.25)", ratio_area),
    );

    // ---- 判定 (記録 — 装置 FAIL ではない) ----
    let all_ok = ratios.iter().all(|&(_, _, r)| (r - 1.0).abs() < 0.01);
    let max_dev = ratios
        .iter()
        .map(|&(_, _, r)| (r - 1.0).abs())
        .fold(0.0f64, f64::max);
    println!(
        "\n[判定] 波束対 4 点の最大 |R−1| = {:.4} → {}",
        max_dev,
        if all_ok {
            "事前登録 (a): 3+1D の第一法則が立った — Einstein 橋の第一段成立"
        } else {
            "事前登録 (b): 系統ずれの記録 (次版で核床/摂動の追い込み)"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v19.1".into())),
        (
            "ratios".into(),
            Json::Arr(
                ratios
                    .iter()
                    .map(|&(n, p, r)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("pair".into(), Json::Int(p as i64)),
                            ("r".into(), Json::Num(r)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "diag_kernel_bias".into(),
            Json::Arr(
                diag_bias
                    .iter()
                    .map(|&(n, b)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("b".into(), Json::Num(b)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("s_half_8".into(), Json::Num(s_half[0].1)),
        ("s_half_12".into(), Json::Num(s_half[1].1)),
        ("area_ratio".into(), Json::Num(ratio_area)),
        ("branch_a".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v191_first3d.json", &j.render());
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
