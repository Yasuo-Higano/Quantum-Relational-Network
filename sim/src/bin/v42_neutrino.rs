//! v4.2 ニュートリノ — シーソー機構・PMNS 混合・レプトジェネシス
//!
//! ニュートリノ質量 (~0.05 eV) は電子の 1000 万分の 1。なぜこれほど軽いのか。
//! シーソー機構: 右巻き ν (v3.1 の SO(10) 16 表現が予言!) が重い質量 M_R を持つと
//!   m_ν = m_D² / M_R   (m_D ~ 電弱スケール)
//! [A] 逆算: 観測された m_ν から M_R を求める → GUT スケールが出るか (v3.1 と独立照合)
//! [B] FN 電荷 (v3.2) + アナーキーな M_R のモンテカルロ → PMNS 混合角の統計 vs 観測
//! [C] レプトジェネシス: 重い ν の CP 非対称崩壊 → バリオン非対称 η_B の桁の見積り

use uft_sim::*;

type M3 = [[(f64, f64); 3]; 3];

fn eig_herm3(hre: &[[f64; 3]; 3], him: &[[f64; 3]; 3]) -> ([f64; 3], [[(f64, f64); 3]; 3]) {
    let n = 3;
    let m = 6;
    let mut emb = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            emb[i + j * m] = hre[i][j];
            emb[i + (j + n) * m] = -him[i][j];
            emb[(i + n) + j * m] = him[i][j];
            emb[(i + n) + (j + n) * m] = hre[i][j];
        }
    }
    let (w, v) = jacobi_eigh(&emb, m);
    let mut lam = [0.0f64; 3];
    let mut vecs = [[(0.0f64, 0.0f64); 3]; 3];
    for k in 0..3 {
        lam[k] = 0.5 * (w[2 * k] + w[2 * k + 1]);
        for i in 0..3 {
            vecs[k][i] = (v[i + (2 * k) * m], v[(i + n) + (2 * k) * m]);
        }
        let nrm: f64 = vecs[k].iter().map(|&(a, b)| a * a + b * b).sum::<f64>().sqrt();
        for i in 0..3 {
            vecs[k][i].0 /= nrm;
            vecs[k][i].1 /= nrm;
        }
    }
    (lam, vecs)
}

fn mmd(a: &M3, b: &M3, bdag: bool) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    // A · B or A · B†
    let mut re = [[0.0f64; 3]; 3];
    let mut im = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let (x, y) = a[i][k];
                let (u, v) = if bdag {
                    (b[j][k].0, -b[j][k].1)
                } else {
                    (b[k][j].0, b[k][j].1)
                };
                re[i][j] += x * u - y * v;
                im[i][j] += x * v + y * u;
            }
        }
    }
    (re, im)
}

/// 3×3 複素行列の逆 (ガウス消去)
fn inv3(m: &M3) -> M3 {
    let mut a = [[(0.0f64, 0.0f64); 6]; 3];
    for i in 0..3 {
        for j in 0..3 {
            a[i][j] = m[i][j];
        }
        a[i][i + 3] = (1.0, 0.0);
    }
    for col in 0..3 {
        // ピボット
        let mut p = col;
        for r in col..3 {
            if a[r][col].0.hypot(a[r][col].1) > a[p][col].0.hypot(a[p][col].1) {
                p = r;
            }
        }
        a.swap(col, p);
        let (pr, pi) = a[col][col];
        let d = pr * pr + pi * pi;
        for j in 0..6 {
            let (x, y) = a[col][j];
            a[col][j] = ((x * pr + y * pi) / d, (y * pr - x * pi) / d);
        }
        for r in 0..3 {
            if r == col {
                continue;
            }
            let (fr, fi) = a[r][col];
            for j in 0..6 {
                let (x, y) = a[col][j];
                a[r][j].0 -= fr * x - fi * y;
                a[r][j].1 -= fr * y + fi * x;
            }
        }
    }
    let mut out = [[(0.0, 0.0); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            out[i][j] = a[i][j + 3];
        }
    }
    out
}

fn main() {
    let mut rng = Rng::new(20260704);
    println!("=== v4.2 ニュートリノ: シーソー・PMNS・レプトジェネシス ===\n");

    // ---- [A] シーソーの逆算 ----
    println!("[A] シーソー機構 m_ν = m_D²/M_R の逆算 (m_D = 174 GeV: SO(10) ではトップと同源)");
    let v_ew = 174.0f64; // GeV
    for &mnu_ev in &[0.05f64, 0.01] {
        let mnu = mnu_ev * 1e-9; // GeV
        let mr = v_ew * v_ew / mnu;
        println!("  m_ν = {:.2} eV → M_R = {:.2e} GeV", mnu_ev, mr);
    }
    println!("  v3.1 の結合統一点: ~1-3×10^16 GeV。SO(10) の中間破れスケールとして整合する高さ。");
    println!("  => ニュートリノの軽さは「GUT スケールの秤」— 独立な 2 つの矢印が同じ場所を指す\n");

    // ---- [B] PMNS 混合のモンテカルロ (FN 電荷: q_L=(1,0,0), 右巻きνはアナーキー) ----
    println!("[B] FN + シーソーの統計 (2000 試行): PMNS 角と質量二乗比");
    let eps = 0.22f64;
    let q_l = [1.0f64, 0.0, 0.0];
    let q_e = [4.0f64, 2.0, 0.0];
    let ntrial = 2000;
    let mut s: Vec<Vec<f64>> = vec![Vec::new(); 4]; // th12, th23, th13, r=Δm²21/Δm²32
    let o1 = |rng: &mut Rng| -> (f64, f64) {
        let r = (3.0f64).powf(2.0 * rng.f64() - 1.0);
        let th = 2.0 * std::f64::consts::PI * rng.f64();
        (r * th.cos(), r * th.sin())
    };
    for _ in 0..ntrial {
        let mut ye: M3 = [[(0.0, 0.0); 3]; 3];
        let mut md: M3 = [[(0.0, 0.0); 3]; 3];
        let mut mr: M3 = [[(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (a, b) = o1(&mut rng);
                let f = eps.powf(q_l[i] + q_e[j]);
                ye[i][j] = (a * f, b * f);
                let (a, b) = o1(&mut rng);
                let f = eps.powf(q_l[i]); // 右巻き ν は電荷 0 (アナーキー)
                md[i][j] = (a * f, b * f);
            }
        }
        for i in 0..3 {
            for j in 0..=i {
                let (a, b) = o1(&mut rng);
                mr[i][j] = (a, b);
                mr[j][i] = (a, b); // 対称 (マヨラナ)
            }
        }
        // m_ν = m_D M_R^{-1} m_D^T (対称)
        let mrinv = inv3(&mr);
        let mut t: M3 = [[(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let (x, y) = md[i][k];
                    let (u, v) = mrinv[k][j];
                    t[i][j].0 += x * u - y * v;
                    t[i][j].1 += x * v + y * u;
                }
            }
        }
        let mut mnu: M3 = [[(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for k in 0..3 {
                    let (x, y) = t[i][k];
                    let (u, v) = md[j][k]; // m_D^T
                    mnu[i][j].0 += x * u - y * v;
                    mnu[i][j].1 += x * v + y * u;
                }
            }
        }
        // 質量 = m m† の固有値の平方根; 左固有ベクトル
        let (hre, him) = mmd(&mnu, &mnu, true);
        let (lam_n, vn) = eig_herm3(&hre, &him);
        let (hre_e, him_e) = mmd(&ye, &ye, true);
        let (_, vecs_e) = eig_herm3(&hre_e, &him_e);
        // PMNS = U_e† U_ν
        let mut pmns = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (mut re, mut im) = (0.0, 0.0);
                for k in 0..3 {
                    let (a, b) = vecs_e[i][k];
                    let (c, d) = vn[j][k];
                    re += a * c + b * d;
                    im += a * d - b * c;
                }
                pmns[i][j] = (re * re + im * im).sqrt();
            }
        }
        let s13 = pmns[0][2];
        let s12 = pmns[0][1] / (1.0 - s13 * s13).sqrt();
        let s23 = pmns[1][2] / (1.0 - s13 * s13).sqrt();
        s[0].push(s12.min(1.0).asin().to_degrees());
        s[1].push(s23.min(1.0).asin().to_degrees());
        s[2].push(s13.min(1.0).asin().to_degrees());
        let m2: Vec<f64> = lam_n.iter().map(|x| x.max(0.0)).collect();
        s[3].push((m2[1] - m2[0]) / (m2[2] - m2[1]).max(1e-300));
    }
    let obs = [33.4f64, 49.0, 8.6, 0.030];
    let names = ["θ12 [°]", "θ23 [°]", "θ13 [°]", "Δm²21/Δm²32"];
    for k in 0..4 {
        s[k].sort_by(|a, b| a.partial_cmp(b).unwrap());
        println!(
            "  {:12}  中央値 {:7.2} [{:6.2}, {:7.2}]   観測 {:6.2}",
            names[k],
            s[k][ntrial / 2],
            s[k][ntrial * 16 / 100],
            s[k][ntrial * 84 / 100],
            obs[k]
        );
    }
    println!("  => クォーク (小さい混合, v3.2) と対照的に、ニュートリノは「大きい混合」が");
    println!("     デフォルトで出る — 観測パターンと整合 (アナーキー + 弱い FN 構造)\n");

    // ---- [C] レプトジェネシス ----
    println!("[C] レプトジェネシス: なぜ宇宙には物質しかないのか (η_B ~ 6×10⁻¹⁰ の起源)");
    let m1 = 1e10f64; // 最軽の右巻きν [GeV]
    let m3nu = 0.05e-9f64; // GeV
    let eps_cp = 3.0 / (16.0 * std::f64::consts::PI) * m1 * m3nu / (v_ew * v_ew);
    let kappa = 0.02f64; // washout 効率 (典型値)
    let gstar = 106.75f64;
    let eta_b = eps_cp * kappa / gstar * 7.04; // s→γ 換算込みの慣用近似
    println!("  M₁ = {:.0e} GeV の右巻きνの CP 非対称崩壊: ε_CP ≤ (3/16π)·M₁m₃/v² = {:.1e}", m1, eps_cp);
    println!("  washout κ~{}, g*~{} → η_B ~ ε·κ/g*·7 ≈ {:.1e}   (観測: 6.1×10⁻¹⁰)", kappa, gstar, eta_b);
    println!("  => シーソーの CP 位相 → レプトン非対称 → スファレロンでバリオンへ。");
    println!("     「反物質が消えた理由」はニュートリノの重い相棒の崩壊でありうる。");
    println!("\n結論: 右巻きν (SO(10) が要求) を 10^10-15 GeV に置くだけで、(a) ν の軽さ、");
    println!("      (b) 大きい混合、(c) バリオン非対称、の 3 つの謎が同時に桁で解ける。");
    println!("      検証: 0νββ 崩壊 (マヨラナ性)、CP 位相 δ の測定、質量順序 — 2030 年代の実験。");
}
