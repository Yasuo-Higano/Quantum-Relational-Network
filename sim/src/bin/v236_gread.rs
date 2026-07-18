//! v23.6 G_lattice の条件付き初読み出し — 第二十期以来の G 窓、開通 (第二十四期)
//!
//! 土台: v23.1 (K_A の NN 形 = 線形 BW ± 1%) + v23.5 (温度繰り込み λ∞ = 32/27,
//! 0.006%)。Jacobson 鎖 δS = δQ/T に T = λ∞/(2π) を通すと **G = a²/(4·λ∞·s_area)**。
//! massless フェルミオン半空間は Gioev–Klich の対数面積則 (S ~ A·ln L) を持ち得る
//! ので、s_area(N) = S(C_A)/N² の N 系列で「定数 vs 対数走行」を判別する:
//!   s(N) = c₀ + c₁·ln N (最小二乗) vs 定数 (平均) — SSR 比で裁く。
//! S は f64+クランプで十分 (深部 κ の切り捨ては S に κe^{−κ} でしか効かない —
//! DD 不要。λ∞ は v23.5 の確定値を定数として使う)。
//!
//! 装置ゲート: [G] 各 N の閉殻ギャップ > 1e-6 (v22.2 の再現)。
//! 事前登録: (a) 対数則が勝つ (SSR_log < 0.5·SSR_const ∧ c₁ > 0) = **G は cutoff
//!   走行結合として読める** — G(N) 表 + G(N=16) を条件付き初読み出しとして記録
//!   (誘導重力 [Sakharov] 描像と接続: 重力結合は cutoff の物性) /
//!   (a′) 定数則が勝つ = G∞ = 27a²/(128·s∞) の一点読み出し /
//!   (b) どちらも合わない (SSR 比 ∈ [0.5, 2]) = 特徴付け不足の記録。
//! いずれも「条件付き」= 自由フェルミオン toy の格子 G (QRN の重力読み出し口の
//! 最初の数値、実重力定数の主張ではない)。

use uft_sim::*;

// 3D staggered H (x 開放, y/z 周期) — v22.2 と同一 (成分は正確な半整数 → DD 誤差ゼロ)
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

fn main() {
    self_test();
    println!("=== v23.6 G_lattice の条件付き初読み出し (第二十四期) ===\n");
    println!("事前登録: (a) 対数面積則 → G は cutoff 走行として読む / (a′) 定数 → G∞ /");
    println!("          (b) どちらも合わない。G = a²/(4·λ∞·s_area), λ∞ = 32/27 (v23.5)\n");
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
    let lam_inf = 32.0 / 27.0; // v23.5 の確定値 (0.006% 一致の有理形)
    let mut rows: Vec<(usize, f64, f64)> = Vec::new(); // (N, s_area, G/a²)

    for &n in &[8usize, 10, 12, 14, 16] {
        let ns = n * n * n;
        let nocc = ns / 2;
        let h = build_h3d_periodic(n);
        let (ev, vv) = jacobi_eigh(&h, ns);
        let gap = ev[nocc] - ev[nocc - 1];
        check(
            &format!("[G] N={} 閉殻ギャップ", n),
            gap > 1e-6,
            format!("gap = {:.4} ({} s)", gap, t0.elapsed().as_secs()),
        );
        // C_A (f64) と S
        let idx3 = |x: usize, y: usize, z: usize| x + n * (y + n * z);
        let half = n / 2;
        let mut sel = Vec::new();
        for z in 0..n {
            for y in 0..n {
                for x in 0..half {
                    sel.push(idx3(x, y, z));
                }
            }
        }
        let m = sel.len();
        let mut ca = vec![0.0f64; m * m];
        for k in 0..nocc {
            for (a, &ia) in sel.iter().enumerate() {
                let va = vv[ia + k * ns];
                if va == 0.0 {
                    continue;
                }
                for (b, &ib) in sel.iter().enumerate() {
                    if b < a {
                        continue;
                    }
                    ca[a * m + b] += va * vv[ib + k * ns];
                }
            }
        }
        for a in 0..m {
            for b in 0..a {
                ca[a * m + b] = ca[b * m + a];
            }
        }
        let (cw, _) = jacobi_eigh(&ca, m);
        let s: f64 = cw
            .iter()
            .map(|&c| {
                let c = c.clamp(1e-14, 1.0 - 1e-14);
                -c * c.ln() - (1.0 - c) * (1.0 - c).ln()
            })
            .sum();
        let s_area = s / ((n * n) as f64);
        let g_over_a2 = 1.0 / (4.0 * lam_inf * s_area);
        println!(
            "    N={:2}: S = {:.4}, s_area = {:.5}, G/a² = {:.5} ({} s)",
            n,
            s,
            s_area,
            g_over_a2,
            t0.elapsed().as_secs()
        );
        rows.push((n, s_area, g_over_a2));
    }

    // ---- フィット: s(N) = c₀ + c₁ ln N vs 定数 ----
    let xs: Vec<f64> = rows.iter().map(|r| (r.0 as f64).ln()).collect();
    let ys: Vec<f64> = rows.iter().map(|r| r.1).collect();
    let np = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / np;
    let my = ys.iter().sum::<f64>() / np;
    let sxy: f64 = xs.iter().zip(&ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let c1 = sxy / sxx;
    let c0 = my - c1 * mx;
    let ssr_log: f64 = xs
        .iter()
        .zip(&ys)
        .map(|(x, y)| {
            let e = y - (c0 + c1 * x);
            e * e
        })
        .sum();
    let ssr_const: f64 = ys.iter().map(|y| (y - my) * (y - my)).sum();
    println!(
        "\n    フィット: s(N) = {:.5} + {:.5}·ln N (SSR {:.3e}) vs 定数 {:.5} (SSR {:.3e}) — 比 {:.3}",
        c0,
        c1,
        ssr_log,
        my,
        ssr_const,
        ssr_log / ssr_const
    );
    let log_wins = ssr_log < 0.5 * ssr_const && c1 > 0.0;
    let const_wins =
        ssr_log > 0.5 * ssr_const || c1.abs() * (16.0f64.ln() - 8.0f64.ln()) < 0.005 * my;
    let g16 = rows.iter().find(|r| r.0 == 16).unwrap().2;
    println!(
        "\n[判定] {}",
        if nfail > 0 {
            "装置ゲート故障 — 記録"
        } else if log_wins {
            "事前登録 (a): 対数面積則 (Gioev–Klich) — G は cutoff 走行結合として開通。G(N) 表が読み出し"
        } else if const_wins {
            "事前登録 (a′): 定数面積則 — G∞ の一点読み出し"
        } else {
            "事前登録 (b): 特徴付け不足 — 記録"
        }
    );
    println!(
        "    [読み出し] G(N=16)/a² = {:.5} (T = 2π/λ∞, λ∞ = 32/27) — 条件付き (自由フェルミオン toy)",
        g16
    );
    if log_wins {
        println!(
            "    [誘導重力接続] c₁/c₀ = {:.4} — G⁻¹ ∝ s_area は cutoff の物性 (Sakharov 描像)",
            c1 / c0
        );
    }

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v23.6".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(n, s, g)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("s_area".into(), Json::Num(s)),
                            ("g_over_a2".into(), Json::Num(g)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("c0".into(), Json::Num(c0)),
        ("c1".into(), Json::Num(c1)),
        ("ssr_ratio".into(), Json::Num(ssr_log / ssr_const)),
        ("g16_over_a2".into(), Json::Num(g16)),
        ("branch_a_log".into(), Json::Bool(nfail == 0 && log_wins)),
    ]);
    let p = write_artifact("results/v236_gread.json", &j.render());
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
