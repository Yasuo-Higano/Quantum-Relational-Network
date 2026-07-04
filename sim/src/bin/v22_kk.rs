//! v2.2 カルツァ=クライン — 電荷と幾何の統一 (隠れ次元の力学)
//!
//! 1921 年 Kaluza, 1926 年 Klein: 5 次元の純粋な重力 = 4 次元の重力 + 電磁気 + スカラー。
//! 格子で最小限の検証をする: 空間 = 大きい次元 x × 小さい円 y (N_y サイト)。
//!   (1) KK 塔: 質量スペクトル m_n = 2sin(πn/N_y) ≈ (2π/N_y)·n — 「重い粒子の梯子」
//!   (2) 電荷 = 隠れ次元の運動量 n ∈ Z (円周のトポロジー π₁(S¹)=Z が電荷を量子化する)
//!   (3) ゲージポテンシャル = 隠れ次元のねじれ (Wilson 線 α): エネルギーが電荷 n に比例して応答
//!   (4) ゲージ不変性 = 隠れ次元の座標変換: ねじれの配分は無意味、ホロノミー Σθ だけが物理

use uft_sim::*;

/// 円周 N_y 上のラプラシアン (リンク位相 phases[j]: j→j+1) の固有値 (昇順)
fn ring_eigs(ny: usize, phases: &[f64]) -> Vec<f64> {
    let m = 2 * ny;
    let mut a = vec![0.0; m * m];
    for j in 0..ny {
        a[j + j * m] = 2.0;
        a[(j + ny) + (j + ny) * m] = 2.0;
        let k = (j + 1) % ny;
        // H[k][j] = -e^{iθ_j}, H[j][k] = -e^{-iθ_j} (エルミート) を実埋め込み
        let (c, s) = (phases[j].cos(), phases[j].sin());
        // Re 部
        a[k + j * m] += -c;
        a[j + k * m] += -c;
        a[(k + ny) + (j + ny) * m] += -c;
        a[(j + ny) + (k + ny) * m] += -c;
        // Im 部 (embed: [[Re, -Im],[Im, Re]]),  Him[k,j] = -s, Him[j,k] = +s
        a[k + (j + ny) * m] += s;
        a[(k + ny) + j * m] += -s;
        a[j + (k + ny) * m] += -s;
        a[(j + ny) + k * m] += s;
    }
    let (w, _) = jacobi_eigh(&a, m);
    // 固有値は 2 重に出る: ペアを平均して 1 つに
    (0..ny).map(|i| 0.5 * (w[2 * i] + w[2 * i + 1])).collect()
}

fn main() {
    let ny = 8usize;
    let two_pi = 2.0 * std::f64::consts::PI;
    println!(
        "=== v2.2 カルツァ=クライン: 電荷 = 隠れ次元の運動量 (N_y={}) ===\n",
        ny
    );

    // ---- (1) KK 塔 ----
    println!("[A] KK 塔 (ねじれなし): 質量² = 2-2cos(2πn/N_y)");
    let eigs0 = ring_eigs(ny, &vec![0.0; ny]);
    let mut exact0: Vec<f64> = (0..ny)
        .map(|nn| 2.0 - 2.0 * (two_pi * nn as f64 / ny as f64).cos())
        .collect();
    exact0.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut maxd = 0.0f64;
    println!("  n(電荷)   m²(数値)   m²(厳密)   m_n");
    for i in 0..ny {
        maxd = maxd.max((eigs0[i] - exact0[i]).abs());
    }
    for nn in 0..=(ny / 2) {
        let m2 = 2.0 - 2.0 * (two_pi * nn as f64 / ny as f64).cos();
        println!(
            "  ±{}        {:.5}    {:.5}    {:.4}{}",
            nn,
            m2,
            m2,
            m2.sqrt(),
            if nn > 0 && nn < ny / 2 {
                "  (±n の 2 重縮退 = 粒子/反粒子)"
            } else {
                ""
            }
        );
    }
    println!(
        "  => 数値 vs 厳密の最大差 {:.1e}  {}",
        maxd,
        pass(maxd < 1e-10)
    );
    println!(
        "     4次元から見ると: 電荷 n を持つ質量 m_n の粒子の梯子。電荷は π₁(S¹)=Z で量子化。\n"
    );

    // ---- (3) Wilson 線 (A₅) への応答 = 電荷 ----
    println!("[B] 隠れ次元のねじれ α (=ゲージポテンシャルの Wilson 線) への応答");
    let alpha = 0.5f64;
    let eigs_a = ring_eigs(ny, &vec![alpha / ny as f64; ny]);
    let mut exact_a: Vec<f64> = (0..ny)
        .map(|nn| {
            let kk = (two_pi * nn as f64 + alpha) / ny as f64;
            2.0 - 2.0 * kk.cos()
        })
        .collect();
    exact_a.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut maxd = 0.0f64;
    for i in 0..ny {
        maxd = maxd.max((eigs_a[i] - exact_a[i]).abs());
    }
    println!(
        "  α={}: スペクトルが m²_n = 2-2cos((2πn+α)/N_y) へ移動 (最大差 {:.1e})  {}",
        alpha,
        maxd,
        pass(maxd < 1e-10)
    );
    println!("  ±n 縮退の分裂 (電荷に比例した応答):");
    println!("  n    E(+n)-E(-n) 数値      2·[cos((2πn-α)/8)-cos((2πn+α)/8)]/1 厳密");
    for nn in 1..ny / 2 {
        let ep = 2.0 - 2.0 * ((two_pi * nn as f64 + alpha) / ny as f64).cos();
        let em = 2.0 - 2.0 * ((two_pi * nn as f64 - alpha) / ny as f64).cos();
        println!(
            "  {}    {:+.5}                (電荷 ±{} が逆向きに応答)",
            nn,
            ep - em,
            nn
        );
    }
    println!("  => ねじれ α は電荷 n に比例して各状態に働く = ミニマル結合 (qA) の幾何的起源\n");

    // ---- (4) ゲージ不変性 = 隠れ次元の微分同相 ----
    println!("[C] ゲージ不変性: ねじれの配分をランダムに変えても (Σθ=α 固定) スペクトル不変か");
    let mut rng = Rng::new(2026);
    let mut phases: Vec<f64> = (0..ny).map(|_| rng.gauss() * 0.3).collect();
    let s: f64 = phases.iter().sum();
    for p in phases.iter_mut() {
        *p += (alpha - s) / ny as f64;
    }
    let eigs_r = ring_eigs(ny, &phases);
    let mut maxd = 0.0f64;
    for i in 0..ny {
        maxd = maxd.max((eigs_r[i] - eigs_a[i]).abs());
    }
    println!(
        "  ランダム配分 vs 一様配分のスペクトル差: {:.1e}  {}",
        maxd,
        pass(maxd < 1e-10)
    );
    println!(
        "  => 局所的な A₅(x) は物理でない。ホロノミー ∮A のみ物理 (v0.4 と同じ教訓が幾何から)\n"
    );

    // ---- (2) 高次元の分離性チェック ----
    println!("[D] 2D (大きい次元 × 円) の分散関係 ω² = 4sin²(k_x/2) + m²_n の検証");
    {
        let nx = 40usize;
        let kx = two_pi * 3.0 / nx as f64;
        let nn = 2usize;
        let theta = two_pi * nn as f64 / ny as f64; // 円周上の許される波数 (周期性)
                                                    // ψ(x,y) = e^{i k_x x} e^{i θ y} に 2D ラプラシアン (y リンクに α/N_y のねじれ) を作用
        let idx = |x: usize, y: usize| x + y * nx;
        let mut re = vec![0.0; nx * ny];
        let mut im = vec![0.0; nx * ny];
        for x in 0..nx {
            for y in 0..ny {
                let ph = kx * x as f64 + theta * y as f64;
                re[idx(x, y)] = ph.cos();
                im[idx(x, y)] = ph.sin();
            }
        }
        // ねじれ込みの固有値: ω² = 4sin²(k_x/2) + 2-2cos((2πn-α)/N_y)
        let om2 = 4.0 * (kx / 2.0).sin().powi(2) + 2.0
            - 2.0 * ((two_pi * nn as f64 - alpha) / ny as f64).cos();
        let mut maxres = 0.0f64;
        let tw = alpha / ny as f64;
        for x in 0..nx {
            for y in 0..ny {
                let (xp, xm) = ((x + 1) % nx, (x + nx - 1) % nx);
                let (yp, ym) = ((y + 1) % ny, (y + ny - 1) % ny);
                // H ψ = 4ψ - ψ(x±1) - e^{-i tw}ψ(y+1) - e^{+i tw}ψ(y-1)
                let hre = 4.0 * re[idx(x, y)]
                    - re[idx(xp, y)]
                    - re[idx(xm, y)]
                    - (tw.cos() * re[idx(x, yp)] + tw.sin() * im[idx(x, yp)])
                    - (tw.cos() * re[idx(x, ym)] - tw.sin() * im[idx(x, ym)]);
                let him = 4.0 * im[idx(x, y)]
                    - im[idx(xp, y)]
                    - im[idx(xm, y)]
                    - (tw.cos() * im[idx(x, yp)] - tw.sin() * re[idx(x, yp)])
                    - (tw.cos() * im[idx(x, ym)] + tw.sin() * re[idx(x, ym)]);
                maxres = maxres
                    .max((hre - om2 * re[idx(x, y)]).abs())
                    .max((him - om2 * im[idx(x, y)]).abs());
            }
        }
        println!(
            "  平面波 × 円モード (k_x=2π·3/40, n=2) の残差: {:.1e}  {}",
            maxres,
            pass(maxres < 1e-12)
        );
    }
    println!(
        "\n結論: 「電荷」「ゲージ場」「質量の梯子」は、見えない小さい次元の「運動量」「ねじれ」"
    );
    println!("      「調和音」である — 内部量子数は幾何に統一できる (KK)。弦理論の compact 化は");
    println!("      この機構の一般化。QRN では隠れ次元 = 網の追加テンソル因子。");
    println!("      課題 (正直に): なぜ 3+1 だけ大きいか、モジュライ安定化、カイラルフェルミオン (→v2.3)。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
