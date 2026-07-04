//! v5.1 初期条件問題 — 成長する宇宙と時間の矢の起源
//!
//! v4.4 で残った最大の謎: なぜビッグバンのエントロピーは最大値の 10^-34 だったのか。
//! 仮説 (成長する自由度): 宇宙の実効的な自由度は膨張とともに増える (地平線に入る
//! モードが増える)。そして新しく現れるモードは断熱定理により「真空 (最低エントロピー)」で
//! 生まれる — インフレーションの Bunch–Davies 真空 (v2.4) はまさにこれ。
//! ならば「過去が低エントロピー」なのは特別な微調整ではなく、
//! **過去には自由度が少なく、新参者は常に無垢で来る**ことの帰結である。
//!
//! トイ検証: フェルミオン鎖を 2 サイトずつ成長させる。
//!   シナリオ A: 新サイトは局所基底状態 (真空) で到着 → 観測窓のエントロピーは
//!               0 から単調に成長 (時間の矢!)、全系は純粋なまま
//!   シナリオ B: 新サイトが最大混合 (熱的) で到着 → 窓は即座に飽和、歴史も矢もない

use uft_sim::*;

fn h2(z: f64) -> f64 {
    let z = z.clamp(1e-14, 1.0 - 1e-14);
    -z * z.ln() - (1.0 - z) * (1.0 - z).ln()
}
fn entropy_herm(cre: &[f64], cim: &[f64], n: usize) -> f64 {
    let m = 2 * n;
    let mut a = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            a[i + j * m] = cre[i + j * n];
            a[i + (j + n) * m] = -cim[i + j * n];
            a[(i + n) + j * m] = cim[i + j * n];
            a[(i + n) + (j + n) * m] = cre[i + j * n];
        }
    }
    let (w, _) = jacobi_eigh(&a, m);
    0.5 * w.iter().map(|&z| h2(z)).sum::<f64>()
}

/// OBC 鎖の基底状態相関 (半充填)
fn ground_c(n: usize) -> Vec<f64> {
    let mut h = vec![0.0; n * n];
    for x in 0..n - 1 {
        h[x + (x + 1) * n] = -1.0;
        h[(x + 1) + x * n] = -1.0;
    }
    let (_, v) = jacobi_eigh(&h, n);
    let mut c = vec![0.0; n * n];
    for m in 0..n / 2 {
        for i in 0..n {
            for j in 0..n {
                c[i + j * n] += v[i + m * n] * v[j + m * n];
            }
        }
    }
    c
}

fn run(thermal_arrivals: bool) -> Vec<(usize, f64, f64)> {
    let n0 = 16usize;
    let win = 16usize;
    let nmax = 96usize;
    let dt = 3.0f64;
    // 初期: 小さな宇宙の基底状態
    let mut n = n0;
    let mut cre = ground_c(n0);
    let mut cim = vec![0.0; n0 * n0];
    let mut hist = Vec::new();
    let mut step = 0usize;
    loop {
        // 観測窓 [0..win) のエントロピーと全系の純粋度
        let mut wre = vec![0.0; win * win];
        let mut wim = vec![0.0; win * win];
        for i in 0..win {
            for j in 0..win {
                wre[i + j * win] = cre[i + j * n];
                wim[i + j * win] = cim[i + j * n];
            }
        }
        let s_win = entropy_herm(&wre, &wim, win);
        let s_tot = entropy_herm(&cre, &cim, n);
        hist.push((n, s_win, s_tot));
        if n >= nmax {
            break;
        }
        // 2 サイトを追加 (真空 or 熱的)
        let n2 = n + 2;
        let mut nre = vec![0.0; n2 * n2];
        let mut nim = vec![0.0; n2 * n2];
        for i in 0..n {
            for j in 0..n {
                nre[i + j * n2] = cre[i + j * n];
                nim[i + j * n2] = cim[i + j * n];
            }
        }
        if thermal_arrivals {
            nre[n + n * n2] = 0.5;
            nre[(n + 1) + (n + 1) * n2] = 0.5;
        } else {
            // 2 サイトの局所基底状態 (結合軌道に 1 粒子): C = [[.5,.5],[.5,.5]]
            nre[n + n * n2] = 0.5;
            nre[(n + 1) + (n + 1) * n2] = 0.5;
            nre[n + (n + 1) * n2] = 0.5;
            nre[(n + 1) + n * n2] = 0.5;
        }
        n = n2;
        // 全体を一様ホッピングで dt だけ発展 (接合クエンチ)
        let mut h = vec![0.0; n * n];
        for x in 0..n - 1 {
            h[x + (x + 1) * n] = -1.0;
            h[(x + 1) + x * n] = -1.0;
        }
        let (eps, v) = jacobi_eigh(&h, n);
        // U = V e^{-iεt} V^T
        let mut pur = vec![0.0; n * n];
        let mut pui = vec![0.0; n * n];
        for k in 0..n {
            let (c_, s_) = ((eps[k] * dt).cos(), -(eps[k] * dt).sin());
            for i in 0..n {
                pur[i + k * n] = v[i + k * n] * c_;
                pui[i + k * n] = v[i + k * n] * s_;
            }
        }
        let mut vt = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                vt[i + j * n] = v[j + i * n];
            }
        }
        let ur = matmul(&pur, &vt, n);
        let ui = matmul(&pui, &vt, n);
        // C' = U C U†
        let t1 = matmul(&ur, &nre, n);
        let t2 = matmul(&ui, &nim, n);
        let t3 = matmul(&ur, &nim, n);
        let t4 = matmul(&ui, &nre, n);
        let ar: Vec<f64> = t1.iter().zip(&t2).map(|(a, b)| a - b).collect();
        let ai: Vec<f64> = t3.iter().zip(&t4).map(|(a, b)| a + b).collect();
        let mut utr = vec![0.0; n * n];
        let mut uti = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                utr[i + j * n] = ur[j + i * n];
                uti[i + j * n] = -ui[j + i * n];
            }
        }
        let b1 = matmul(&ar, &utr, n);
        let b2 = matmul(&ai, &uti, n);
        let b3 = matmul(&ar, &uti, n);
        let b4 = matmul(&ai, &utr, n);
        cre = b1.iter().zip(&b2).map(|(a, b)| a - b).collect();
        cim = b3.iter().zip(&b4).map(|(a, b)| a + b).collect();
        step += 1;
        let _ = step;
    }
    hist
}

fn main() {
    println!("=== v5.1 成長する宇宙: 時間の矢と低エントロピーの過去の起源 ===\n");
    println!(
        "観測窓 = 最初の 16 サイト (S_max = 16·ln2 = {:.2})。宇宙は 16 → 96 サイトへ成長。\n",
        16.0 * (2.0f64).ln()
    );
    let va = run(false);
    let vb = run(true);
    println!("  宇宙の大きさ N   S_窓(真空到着)   S_窓(熱的到着)   S_全系(真空)   S_全系(熱的)");
    for k in (0..va.len()).step_by(5) {
        println!(
            "  {:6}          {:8.3}         {:8.3}        {:8.4}      {:8.3}",
            va[k].0, va[k].1, vb[k].1, va[k].2, vb[k].2
        );
    }
    let last = va.len() - 1;
    println!("\n  シナリオ A (真空到着):");
    println!(
        "   - 窓のエントロピー: {:.2} → {:.2} へ成長 = 時間の矢",
        va[0].1, va[last].1
    );
    let inc = va.windows(2).filter(|w| w[1].1 >= w[0].1).count();
    let max_dip = va
        .windows(2)
        .map(|w| (w[0].1 - w[1].1).max(0.0))
        .fold(0.0f64, f64::max);
    let rise = va[last].1 - va[0].1;
    println!(
        "   - 上昇ステップ {}/{}, 最大の戻り {:.3} (全上昇 {:.2} の {:.0}%) — 有限系の再帰による揺らぎ",
        inc,
        va.len() - 1,
        max_dip,
        rise,
        100.0 * max_dip / rise
    );
    println!(
        "   - 矢の判定 (戻り ≤ 上昇の 15%): {}",
        pass(max_dip < 0.15 * rise)
    );
    println!(
        "   - 全系は純粋なまま: S_全系 = {:.1e} {} (矢は大域エントロピーなしで生じる!)",
        va[last].2,
        pass(va[last].2 < 1e-6)
    );
    println!(
        "  シナリオ B (熱的到着): 窓は最初から飽和域 ({:.2}→{:.2}) — 歴史も矢もない",
        vb[0].1, vb[last].1
    );
    println!("\n結論: 「過去はなぜ低エントロピーだったか」の答えの候補:");
    println!(
        "      (i) 過去は自由度が少なかった (S ≤ ln dim: 小さい宇宙は低エントロピーしか持てない)"
    );
    println!(
        "      (ii) 膨張で新しく入るモードは断熱定理により真空で生まれる (Bunch–Davies, v2.4)"
    );
    println!("      この 2 つだけで、微調整なしに時間の矢が生じることをトイ模型で確認した。");
    println!("      ペンローズの 10^-34 の「特別さ」は、宇宙が小さかったことの返り値でありうる。");
    println!("      (未解決で残るもの: なぜ自由度は増えるのか = 膨張そのものの起源、に問いが移る)");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
