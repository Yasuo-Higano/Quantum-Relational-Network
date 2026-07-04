//! v0.5 幾何=情報 — 調和格子のエンタングルメント・エントロピー (厳密計算)
//! H = ½Σp² + ½Σ(x_i - x_j)² + ½m²Σx²  (最近接結合の調和格子 = 自由スカラー場の格子化)
//!
//! 基底状態はガウス状態なので、部分領域 A のエントロピーは相関行列から厳密に求まる:
//!   X_ij = ⟨x_i x_j⟩, P_ij = ⟨p_i p_j⟩ を A に制限し、
//!   シンプレクティック固有値 ν_a = sqrt(eig(X_A P_A)) から
//!   S = Σ_a [(ν+½)ln(ν+½) - (ν-½)ln(ν-½)]
//!
//! A: 1D 臨界 (m→0):  S(ℓ) = (c/3)·ln[コード長] + 定数, c=1 → 対数則
//! B: 1D 質量あり:    S(ℓ) → 定数 (1次元の面積則)
//! C: 2D:            S ∝ 境界の長さ (面積則) — ブラックホールエントロピー S=A/4G の縮図

use uft_sim::*;

fn entropy_from_xp(xa: &[f64], pa: &[f64], n: usize) -> f64 {
    // sqrt(X_A) を作り M = √X P √X (対称) の固有値 = ν²
    let sx = matfun_sym(xa, n, |lam| lam.max(0.0).sqrt());
    let m1 = matmul(&sx, pa, n);
    let m2 = matmul(&m1, &sx, n);
    // 対称化 (数値誤差の掃除)
    let mut msym = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            msym[i + j * n] = 0.5 * (m2[i + j * n] + m2[j + i * n]);
        }
    }
    let (w, _) = jacobi_eigh(&msym, n);
    w.iter()
        .map(|&mu| {
            let nu = mu.max(0.25).sqrt();
            let (a, b) = (nu + 0.5, nu - 0.5);
            if b < 1e-14 {
                0.0
            } else {
                a * a.ln() - b * b.ln()
            }
        })
        .sum()
}

/// 1D 鎖 (PBC): 距離 d の相関 ⟨xx⟩, ⟨pp⟩ テーブル
fn corr_1d(n: usize, m: f64) -> (Vec<f64>, Vec<f64>) {
    let mut xc = vec![0.0; n];
    let mut pc = vec![0.0; n];
    for d in 0..n {
        let (mut sx, mut sp) = (0.0, 0.0);
        for k in 0..n {
            let q = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
            let w = (m * m + 2.0 - 2.0 * q.cos()).sqrt();
            let c = (q * d as f64).cos();
            sx += c / (2.0 * w);
            sp += c * w / 2.0;
        }
        xc[d] = sx / n as f64;
        pc[d] = sp / n as f64;
    }
    (xc, pc)
}

fn s_interval_1d(xc: &[f64], pc: &[f64], n: usize, l: usize) -> f64 {
    let mut xa = vec![0.0; l * l];
    let mut pa = vec![0.0; l * l];
    for i in 0..l {
        for j in 0..l {
            let d = (i as isize - j as isize).unsigned_abs() % n;
            xa[i + j * l] = xc[d];
            pa[i + j * l] = pc[d];
        }
    }
    entropy_from_xp(&xa, &pa, l)
}

fn main() {
    self_test();
    println!("=== v0.5 幾何=情報: エンタングルメント・エントロピーのスケーリング ===\n");

    // ---- A: 1D 臨界(質量 0, 開境界=Dirichlet でゼロモード回避) 対数則と中心電荷 ----
    let n = 1024;
    println!(
        "[A] 1D 臨界鎖 (N={}, m=0, 開境界): S(ℓ) = (c/6)·ln[(2N/π)sin(πℓ/N)] + k",
        n
    );
    // OBC の厳密モード: v_k(i)=√(2/(N+1))sin(πki/(N+1)), λ_k = 2-2cos(πk/(N+1))
    let lmax = 256usize;
    let nf = (n + 1) as f64;
    let mut xa_full = vec![0.0; lmax * lmax];
    let mut pa_full = vec![0.0; lmax * lmax];
    for k in 1..=n {
        let lam = 2.0 - 2.0 * (std::f64::consts::PI * k as f64 / nf).cos();
        let sq = lam.sqrt();
        let mut vk = vec![0.0; lmax];
        for i in 0..lmax {
            vk[i] =
                (2.0 / nf).sqrt() * (std::f64::consts::PI * k as f64 * (i + 1) as f64 / nf).sin();
        }
        for i in 0..lmax {
            for j in 0..=i {
                let t = vk[i] * vk[j];
                xa_full[i + j * lmax] += t / (2.0 * sq);
                pa_full[i + j * lmax] += t * sq / 2.0;
            }
        }
    }
    for i in 0..lmax {
        for j in (i + 1)..lmax {
            xa_full[i + j * lmax] = xa_full[j + i * lmax];
            pa_full[i + j * lmax] = pa_full[j + i * lmax];
        }
    }
    let ells = [4usize, 8, 16, 32, 64, 128, 256];
    let mut lnchord = Vec::new();
    let mut svals = Vec::new();
    println!("  ℓ     S(ℓ)      共形長");
    for &l in &ells {
        // 端に接する区間 [1..ℓ] (エンタングリング点は 1 つ → 係数 c/6)
        let mut xa = vec![0.0; l * l];
        let mut pa = vec![0.0; l * l];
        for i in 0..l {
            for j in 0..l {
                xa[i + j * l] = xa_full[i + j * lmax];
                pa[i + j * l] = pa_full[i + j * lmax];
            }
        }
        let s = entropy_from_xp(&xa, &pa, l);
        let chord =
            (2.0 * nf / std::f64::consts::PI) * (std::f64::consts::PI * l as f64 / nf).sin();
        println!("  {:4}  {:.5}   {:8.2}", l, s, chord);
        lnchord.push(chord.ln());
        svals.push(s);
    }
    let (_, slope) = linfit(&lnchord, &svals);
    let c_eff = 6.0 * slope;
    println!(
        "  => フィット勾配 = {:.4} → 中心電荷 c = 6×勾配 = {:.3}  (自由ボソン場の理論値 c=1)  {}",
        slope,
        c_eff,
        pass((c_eff - 1.0).abs() < 0.1)
    );
    println!("     臨界(質量ゼロ)ではもつれが全スケールに広がり、S は対数発散する。\n");

    // ---- B: 1D 質量あり → 飽和 (1D の「面積則」) ----
    println!("[B] 1D 質量あり (m=0.5): 相関長 1/m を超えると S が飽和する");
    let (xc2, pc2) = corr_1d(n, 0.5);
    println!("  ℓ     S(ℓ)");
    let mut last = 0.0;
    for &l in &[2usize, 4, 8, 16, 32, 64] {
        let s = s_interval_1d(&xc2, &pc2, n, l);
        println!("  {:4}  {:.6}", l, s);
        last = s;
    }
    let s8 = s_interval_1d(&xc2, &pc2, n, 8);
    println!(
        "  => S(8)={:.5} と S(64)={:.5} の差 {:.1e}: 飽和 [面積則]  {}",
        s8,
        last,
        (last - s8).abs(),
        pass((last - s8).abs() < 0.01)
    );
    println!("     ギャップ系ではもつれは境界近傍に局在する。\n");

    // ---- C: 2D 面積則 ----
    let l2 = 32usize;
    let m2d = 0.2;
    println!(
        "[C] 2D 格子 ({}×{}, m={}): ℓ×ℓ ブロックの S vs 境界長 4ℓ",
        l2, l2, m2d
    );
    // 相関テーブル (dx,dy)
    let mut xc2d = vec![0.0; l2 * l2];
    let mut pc2d = vec![0.0; l2 * l2];
    for dy in 0..l2 {
        for dx in 0..l2 {
            let (mut sx, mut sp) = (0.0, 0.0);
            for ky in 0..l2 {
                for kx in 0..l2 {
                    let qx = 2.0 * std::f64::consts::PI * kx as f64 / l2 as f64;
                    let qy = 2.0 * std::f64::consts::PI * ky as f64 / l2 as f64;
                    let w = (m2d * m2d + 4.0 - 2.0 * qx.cos() - 2.0 * qy.cos()).sqrt();
                    let c = (qx * dx as f64 + qy * dy as f64).cos();
                    sx += c / (2.0 * w);
                    sp += c * w / 2.0;
                }
            }
            xc2d[dx + dy * l2] = sx / (l2 * l2) as f64;
            pc2d[dx + dy * l2] = sp / (l2 * l2) as f64;
        }
    }
    let mut per = Vec::new();
    let mut s2d = Vec::new();
    println!("  ℓ×ℓ   境界長  S        S/境界長");
    for &l in &[2usize, 3, 4, 5, 6, 7, 8] {
        let nn = l * l;
        let mut xa = vec![0.0; nn * nn];
        let mut pa = vec![0.0; nn * nn];
        for a in 0..nn {
            for b in 0..nn {
                let (ax, ay) = (a % l, a / l);
                let (bx, by) = (b % l, b / l);
                let dx = (ax as isize - bx as isize).unsigned_abs();
                let dy = (ay as isize - by as isize).unsigned_abs();
                xa[a + b * nn] = xc2d[dx + dy * l2];
                pa[a + b * nn] = pc2d[dx + dy * l2];
            }
        }
        let s = entropy_from_xp(&xa, &pa, nn);
        println!(
            "  {}×{}   {:4}   {:.4}   {:.4}",
            l,
            l,
            4 * l,
            s,
            s / (4 * l) as f64
        );
        per.push((4 * l) as f64);
        s2d.push(s);
    }
    let (icpt, slope2) = linfit(&per, &s2d);
    // 決定係数
    let mean_s: f64 = s2d.iter().sum::<f64>() / s2d.len() as f64;
    let ss_tot: f64 = s2d.iter().map(|s| (s - mean_s).powi(2)).sum();
    let ss_res: f64 = per
        .iter()
        .zip(&s2d)
        .map(|(&p, &s)| (s - icpt - slope2 * p).powi(2))
        .sum();
    let r2 = 1.0 - ss_res / ss_tot;
    println!(
        "  => S = {:.4}·(境界長) + {:.3}, 決定係数 R² = {:.5}  {}",
        slope2,
        icpt,
        r2,
        pass(r2 > 0.999)
    );
    println!("     エントロピーは体積 ℓ² でなく境界 4ℓ に比例する = 面積則。");
    println!("\n結論: 基底状態のもつれは境界に宿る。「エントロピー ∝ 面積」はブラックホール");
    println!("      S = A/4G と同じスケーリング — 幾何(面積)と情報(もつれ)の等式の最初の実証。");
    println!("      時空の「面積」とは、そこを横切る量子相関の量なのではないか (→ v0.7)。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
