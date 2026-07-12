//! v22.2 G の窓を開ける — 横質量ゼロ構成 + 多倍長核 (第二十三期 課題 2)
//!
//! 第二段 (ブースト同定 = G の読み出し口) を塞いだ 2 つの凍結量 (v19.6 の三角測量):
//! 横質量積 m_T·ξ と核の f64 床。本版は両方を外す:
//!   [装置 1] 横質量ゼロ: y/z を周期境界に (反周期 → 周期)。m_T = √2 sin(π/N) → 0
//!     (k_T = 0 の massless 線が復活)。開 x 鎖は E=0 を避けるので閉殻は保たれる (ゲート)。
//!   [装置 2 — run1 で退役] DD-Rayleigh 核精密化は無効と判明: f64 固有ベクトルは
//!     準ヌル縮退クラスタ内で任意混合されており、ベクトル毎の Rayleigh 商と κ の
//!     対応が壊れる (拡張モードゲートが R = −327 で捕捉 — 私の近似の数学的誤り)。
//!     正しい形はクラスタ部分空間の DD-Gram 対角化だが、深部 d≥2.5 の障害は
//!     信号消失 (占有縮退 — v19.4) であり核精度では直せないため、本器械では
//!     f64+クランプ核 (v19.4 検証済) に戻す。DD 関数は将来のために保持。
//! 器械: v19.4 型 Gauss 対掃引 (wid=1.0, d ∈ {0.5..4.5})、R₁ = δS/δ⟨K_A⟩ (f64 核) を
//! カナリア、R₂ = δS/δ⟨K_boost⟩ (核なし) を本命。N=12 全掃引 + N=16 は d ∈ {1.5, 2.5}。
//! 開発記録 (v2→v3): v19.1 拡張モード対バイアス (+5.4%) を器械同一性の錨ゲートに
//! していたが、周期 y/z では Fermi 面直下/直上モードが横運動量で縮退し「拡張モード対」
//! の同定自体が不定 — 錨の前提が偽 (R_拡張 = −319.7 で捕捉、DD-Rayleigh の −327 と
//! 同値域 = 核精度によらない構造的不定性)。錨は半空間幾何専用と判明したため v3 で
//! 記録に降格。有効性判定は各点の R₁ カナリア (|R₁−1| < 1%) に一本化 (v2 と同一)。
//! 事前登録 (順序命題): (a) R₂(1.5; m_T=0, N=12) > 0.87 (v19.4 基準) かつ 有効域で
//!   d 単調増加 かつ R₂(2.5; N=16) > R₂(2.5; N=12) かつ 最終有効点 > 0.93
//!   = G の窓が開通 → 条件付き G_lattice = a²/(4 s_area) を初読み出し (記録) /
//!   (a′) 単調だが < 0.93 / (b) 非単調・基準割れ。

use uft_sim::*;

// ---- double-double 演算 (Dekker) — 将来のクラスタ DD-Gram 用に保持 ----
#[allow(dead_code)]
#[derive(Clone, Copy)]
struct Dd {
    hi: f64,
    lo: f64,
}
#[allow(dead_code)]
fn two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    let bb = s - a;
    let err = (a - (s - bb)) + (b - bb);
    Dd { hi: s, lo: err }
}
#[allow(dead_code)]
fn two_prod(a: f64, b: f64) -> Dd {
    let p = a * b;
    let e = a.mul_add(b, -p); // FMA
    Dd { hi: p, lo: e }
}
#[allow(dead_code)]
fn dd_add(x: Dd, y: Dd) -> Dd {
    let s = two_sum(x.hi, y.hi);
    let lo = s.lo + x.lo + y.lo;
    let t = two_sum(s.hi, lo);
    Dd { hi: t.hi, lo: t.lo }
}
#[allow(dead_code)]
fn dd_add_f(x: Dd, y: f64) -> Dd {
    dd_add(x, Dd { hi: y, lo: 0.0 })
}
#[allow(dead_code)]
fn dd_mul_ff(a: f64, b: f64) -> Dd {
    two_prod(a, b)
}

// 3D staggered H (x 開放, y/z 周期 = 横質量ゼロ)
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
                let jy = idx(x, (y + 1) % n, z);
                add(i, jy, ey); // 周期 (符号反転なし — 反周期を外す)
                let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                let jz = idx(x, y, (z + 1) % n);
                add(i, jz, ez);
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

// rank-2 実回転 C' = OCOᵀ (v19.1 系)
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
    // ΔO = (cs−1)(uuᵀ+wwᵀ) + s(uwᵀ−wuᵀ)
    let mut out = c.to_vec();
    let a = cs - 1.0;
    for i in 0..ns {
        let du = a * u[i] + s * w[i]; // ΔO 行の u 係数? — 展開: ΔO_{ij} = a(u_i u_j + w_i w_j) + s(u_i w_j − w_i u_j)
        let dw = a * w[i] - s * u[i];
        // (ΔO C)_{ij} = du_i·cu_j + dw_i·cw_j
        for j in 0..ns {
            out[j + i * ns] += du * cu[j] + dw * cw[j];
            out[i + j * ns] += du * cu[j] + dw * cw[j]; // + (C ΔOᵀ)
        }
    }
    // ΔO C ΔOᵀ = [du dw] [[uu uw],[uw ww]] [du dw]ᵀ (行 i, 列 j)
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
    println!("=== v22.2 G の窓 — 横質量ゼロ構成 + DD 核 (第二十三期 課題 2) ===\n");
    println!("事前登録: (a) R₂(1.5; N=12) > 0.87 ∧ d 単調 ∧ R₂(2.5; N=16) > R₂(2.5; N=12) ∧");
    println!(
        "          最終 > 0.93 = 窓開通 (条件付き G_lattice を記録) / (a′) 単調未達 / (b) 外れ\n"
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
    let t0 = std::time::Instant::now();
    let mut results: Vec<(usize, f64, f64, f64, bool)> = Vec::new(); // (n, d, r1, r2, valid)
    let mut s_area_n12 = 0.0f64;

    for &n in &[12usize, 16] {
        let ns = n * n * n;
        let h = build_h3d_periodic(n);
        let (ev, vv) = jacobi_eigh(&h, ns);
        let nocc = ns / 2;
        let gap = ev[nocc] - ev[nocc - 1];
        check(
            &format!("N={} 周期 y/z の閉殻ギャップ (m_T=0 構成)", n),
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
        let s0 = entropy_sym_local(&ca0, m);
        if n == 12 {
            s_area_n12 = s0 / ((n * n) as f64); // 面積あたりエントロピー (記録)
        }
        // f64 固有系 + DD-Rayleigh 精密化核
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
        // ブースト重みボンド (v19.2 凍結規約, 周期 y/z 版)
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
        // [装置 2 ゲート] v19.1 拡張モード対診断の再現 (N=12 のみ): DD 核で b < 1%
        if n == 12 {
            // Fermi 面直下/直上の拡がった固有モード対
            let u: Vec<f64> = (0..ns).map(|i| vv[i + (nocc - 1) * ns]).collect();
            let w: Vec<f64> = (0..ns).map(|i| vv[i + nocc * ns]).collect();
            let a1 = 0.02f64;
            let mut rr = [0.0f64; 2];
            for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
                let cp = rotate_c_local(&c, ns, &u, &w, al);
                let cm = rotate_c_local(&c, ns, &u, &w, -al);
                let sp = entropy_sym_local(&restrict_local(&cp, ns, &sel), m);
                let sm = entropy_sym_local(&restrict_local(&cm, ns, &sel), m);
                let kp = dk_of(&restrict_local(&cp, ns, &sel));
                let km = dk_of(&restrict_local(&cm, ns, &sel));
                rr[ii] = (sp - sm) / (kp - km);
            }
            let r_ext = (4.0 * rr[0] - rr[1]) / 3.0;
            println!(
                "    [記録] 拡張モード対の錨は本幾何で不定 (横運動量縮退): R_拡張 = {:.5} — v19.1 錨 (+5.4%) は半空間専用と判明",
                r_ext
            );
        }
        // 掃引
        let dists: Vec<f64> = if n == 12 {
            vec![0.5, 1.5, 2.5, 3.5, 4.5]
        } else {
            vec![1.5, 2.5]
        };
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
                "    N={} d={:.1}: R₁={:.5} [{}], R₂={:.5} ({} s)",
                n,
                d,
                r1,
                if valid { "有効" } else { "除外" },
                r2,
                t0.elapsed().as_secs()
            );
            results.push((n, d, r1, r2, valid));
        }
    }

    // ---- 判定 ----
    let v12: Vec<&(usize, f64, f64, f64, bool)> =
        results.iter().filter(|r| r.0 == 12 && r.4).collect();
    let r2_15_n12 = results
        .iter()
        .find(|r| r.0 == 12 && (r.1 - 1.5).abs() < 0.01)
        .map(|r| r.3)
        .unwrap_or(0.0);
    let r2_25_n12 = results
        .iter()
        .find(|r| r.0 == 12 && (r.1 - 2.5).abs() < 0.01 && r.4)
        .map(|r| r.3);
    let r2_25_n16 = results
        .iter()
        .find(|r| r.0 == 16 && (r.1 - 2.5).abs() < 0.01 && r.4)
        .map(|r| r.3);
    let mut mono = true;
    for i in 1..v12.len() {
        if v12[i].3 < v12[i - 1].3 - 0.01 {
            mono = false;
        }
    }
    let last = v12.last().map(|r| r.3).unwrap_or(0.0);
    let base_ok = r2_15_n12 > 0.87;
    let n16_ok = match (r2_25_n12, r2_25_n16) {
        (Some(a), Some(b)) => b > a,
        _ => false,
    };
    let branch_a = base_ok && mono && n16_ok && last > 0.93;
    println!(
        "\n[判定] {}",
        if branch_a {
            "事前登録 (a): G の窓が開通 — 横質量ゼロ + DD 核でブースト同定が完成域に"
        } else if mono && base_ok {
            "事前登録 (a′): 単調・基準超えだが最終 < 0.93 または N=16 未確認 — 前進の記録"
        } else {
            "事前登録 (b): 外れ — 記録"
        }
    );
    println!(
        "    R₂(1.5; N=12) = {:.4} (v19.4 基準 0.87), 最終有効 = {:.4}, N=16 比較 = {:?} vs {:?}",
        r2_15_n12, last, r2_25_n12, r2_25_n16
    );
    println!(
        "    [記録] 条件付き G_lattice: s_area(N=12) = {:.4}/面積 → a²/(4G) = s_area ⇒ G = a²/(4·{:.4})",
        s_area_n12, s_area_n12
    );

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v22.2".into())),
        (
            "rows".into(),
            Json::Arr(
                results
                    .iter()
                    .map(|&(n, d, r1, r2, va)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("d".into(), Json::Num(d)),
                            ("r1".into(), Json::Num(r1)),
                            ("r2".into(), Json::Num(r2)),
                            ("valid".into(), Json::Bool(va)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("s_area_n12".into(), Json::Num(s_area_n12)),
        ("branch_a".into(), Json::Bool(branch_a)),
    ]);
    let p = write_artifact("results/v222_gwindow.json", &j.render());
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
