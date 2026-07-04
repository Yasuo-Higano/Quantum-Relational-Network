//! v13.1 Lanczos 疎ソルバの検証 — 大格子への道具 (残高 11)
//!
//! 稠密ヤコビ法は実埋め込み ~2600 次元 (N=6 の T⁴) が実用上限で、N=18 の T⁴
//! (10 万サイト) には届かない (v12.3 の道具の限界)。lib.rs に自作した
//! `lanczos_lowest_herm` (複素エルミート・完全再直交化・残差検証つき) を 3 段で検証する:
//!
//!  [1] 小さな稠密ランダムエルミート行列 (n=40): 最低 6 固有値がヤコビ法と厳密一致
//!  [2] 2D 磁束トーラス (N=18, Q=3): 縮退 3 の最低バンド — 幅とギャップが一次ソース
//!      (results/v72_geomfn.txt: ギャップ 0.115) と一致。縮退クラスタの分解能の検証
//!  [3] N=6 傾き T⁴ (2,2,1,−1): v12.1 (ヤコビ経路) のギャップ 0.3195 への回帰 —
//!      同じ物理を疎経路で再現し、実行時間の比 (~10 分 → 秒) を記録
//!  [4] 決定論と初期値独立性: 同シードで再現、異シードで固有値一致
//!
//! これが通れば v13.2 (N=18 傾き T⁴ の深い磁束族) が開く。

use uft_sim::*;

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

// ---------- [2] 用: 2D 磁束トーラス (v72 と同じ物理; 複素で直接) ----------
const N2: usize = 18;
const NS2: usize = N2 * N2;
const Q2D: usize = 3;

fn hops_2d() -> Vec<(usize, usize, f64)> {
    let phi = 2.0 * std::f64::consts::PI * Q2D as f64 / NS2 as f64;
    let idx = |x: usize, y: usize| x + y * N2;
    let mut hops = Vec::new();
    for x in 0..N2 {
        for y in 0..N2 {
            hops.push((idx(x, y), idx(x, (y + 1) % N2), phi * x as f64));
            let th = if x == N2 - 1 {
                -phi * (N2 as f64) * y as f64
            } else {
                0.0
            };
            hops.push((idx(x, y), idx((x + 1) % N2, y), th));
        }
    }
    hops
}

// ---------- [3] 用: N=6 傾き T⁴ (v12.1 と同一) ----------
const N4: usize = 6;
const NS4: usize = N4 * N4 * N4 * N4;
const FLUX: [i64; 4] = [2, 2, 1, -1];

fn idx4(x1: usize, y1: usize, x2: usize, y2: usize) -> usize {
    x1 + N4 * (y1 + N4 * (x2 + N4 * y2))
}

fn link_phase4(x1: usize, y1: usize, x2: usize, y2: usize, dir: usize) -> f64 {
    let nn = (N4 * N4) as f64;
    let two_pi = 2.0 * std::f64::consts::PI;
    let (p1, p2, pt, ps) = (
        two_pi * FLUX[0] as f64 / nn,
        two_pi * FLUX[1] as f64 / nn,
        two_pi * FLUX[2] as f64 / nn,
        two_pi * FLUX[3] as f64 / nn,
    );
    let nf = N4 as f64;
    match dir {
        0 => {
            if x1 == N4 - 1 {
                -(p1 * nf * y1 as f64 + pt * nf * y2 as f64)
            } else {
                0.0
            }
        }
        1 => {
            let mut th = p1 * x1 as f64;
            if y1 == N4 - 1 {
                th += -(ps * nf * x2 as f64);
            }
            th
        }
        2 => {
            let mut th = ps * y1 as f64;
            if x2 == N4 - 1 {
                th += -(p2 * nf * y2 as f64);
            }
            th
        }
        3 => p2 * x2 as f64 + pt * x1 as f64,
        _ => unreachable!(),
    }
}

fn hops_t4() -> Vec<(usize, usize, f64)> {
    let mut hops = Vec::new();
    for x1 in 0..N4 {
        for y1 in 0..N4 {
            for x2 in 0..N4 {
                for y2 in 0..N4 {
                    let i = idx4(x1, y1, x2, y2);
                    hops.push((i, idx4((x1 + 1) % N4, y1, x2, y2), link_phase4(x1, y1, x2, y2, 0)));
                    hops.push((i, idx4(x1, (y1 + 1) % N4, x2, y2), link_phase4(x1, y1, x2, y2, 1)));
                    hops.push((i, idx4(x1, y1, (x2 + 1) % N4, y2), link_phase4(x1, y1, x2, y2, 2)));
                    hops.push((i, idx4(x1, y1, x2, (y2 + 1) % N4), link_phase4(x1, y1, x2, y2, 3)));
                }
            }
        }
    }
    hops
}

/// hop リスト (i→j, θ) から複素 matvec を作る: H[i][j] = −e^{−iθ}, H[j][i] = −e^{+iθ}
fn matvec_hops(hops: &[(usize, usize, f64)], n: usize, v: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut o = vec![(0.0f64, 0.0f64); n];
    for &(i, j, th) in hops {
        let (c, s) = (th.cos(), th.sin());
        let (br, bi) = v[j];
        o[i].0 += -(c * br + s * bi);
        o[i].1 += -(c * bi - s * br);
        let (ar, ai) = v[i];
        o[j].0 += -(c * ar - s * ai);
        o[j].1 += -(c * ai + s * ar);
    }
    o
}

fn main() {
    self_test();
    println!("=== v13.1 Lanczos 疎ソルバの検証: 大格子への道具 ===\n");

    // ---- [1] 小行列: ヤコビ法と厳密一致 ----
    println!("[1] 稠密ランダムエルミート (n=40): 最低 6 固有値がヤコビ法と一致するか");
    let n = 40usize;
    let mut rng = Rng::new(13101);
    let mut hre = vec![0.0f64; n * n];
    let mut him = vec![0.0f64; n * n];
    for i in 0..n {
        for j in 0..=i {
            let a = rng.gauss();
            hre[i + j * n] = a;
            hre[j + i * n] = a;
            if i != j {
                let b = rng.gauss();
                him[i + j * n] = b;
                him[j + i * n] = -b;
            }
        }
    }
    // 実埋め込みでヤコビ
    let m2 = 2 * n;
    let mut emb = vec![0.0f64; m2 * m2];
    for i in 0..n {
        for j in 0..n {
            emb[i + j * m2] = hre[i + j * n];
            emb[i + (j + n) * m2] = -him[i + j * n];
            emb[(i + n) + j * m2] = him[i + j * n];
            emb[(i + n) + (j + n) * m2] = hre[i + j * n];
        }
    }
    let (w_j, _) = jacobi_eigh(&emb, m2);
    let dense_mv = |v: &[(f64, f64)]| -> Vec<(f64, f64)> {
        let mut o = vec![(0.0f64, 0.0f64); n];
        for i in 0..n {
            for j in 0..n {
                let (br, bi) = v[j];
                let (a, b) = (hre[i + j * n], him[i + j * n]);
                o[i].0 += a * br - b * bi;
                o[i].1 += a * bi + b * br;
            }
        }
        o
    };
    let (ev, _, res1) = lanczos_lowest_herm(&dense_mv, n, 6, 40, 4242);
    let mut max_d1: f64 = 0.0;
    for e in 0..6 {
        max_d1 = max_d1.max((ev[e] - w_j[2 * e]).abs());
    }
    let ok1 = max_d1 < 1e-9 && res1 < 1e-8;
    println!(
        "    max|Δλ| = {:.2e}, 残差 {:.2e}  {}",
        max_d1,
        res1,
        pass(ok1)
    );

    // ---- [2] 2D 磁束トーラス: 縮退 3 の最低バンド ----
    println!("\n[2] 2D 磁束トーラス (N=18, Q=3): 縮退クラスタの分解能");
    let hops2 = hops_2d();
    let mv2 = |v: &[(f64, f64)]| matvec_hops(&hops2, NS2, v);
    let t0 = std::time::Instant::now();
    let (ev2, _, res2) = lanczos_lowest_herm(&mv2, NS2, 4, 324, 777); // n=324 は全次元で厳密 (低磁束密度の 2D は低準位が密で m=120 では未収束 — 残差検査が検出した)
    let spread2 = ev2[2] - ev2[0];
    let gap2 = ev2[3] - ev2[2];
    // 一次ソース: results/v72_geomfn.txt [1] — 縮退幅 ~1e-13 ≪ ギャップ 0.115
    let ok2 = spread2.abs() < 1e-8 && (gap2 - 0.115).abs() < 0.005 && res2 < 1e-8;
    println!(
        "    幅 {:.2e} / ギャップ {:.4} (一次ソース 0.115) / 残差 {:.2e}  {}  ({} ms)",
        spread2,
        gap2,
        res2,
        pass(ok2),
        t0.elapsed().as_millis()
    );

    // ---- [3] N=6 傾き T⁴: v12.1 への回帰 + 速度 ----
    println!("\n[3] N=6 傾き T⁴ (2,2,1,−1): v12.1 (ヤコビ ~10 分) への回帰");
    let hops4 = hops_t4();
    let mv4 = |v: &[(f64, f64)]| matvec_hops(&hops4, NS4, v);
    let t1 = std::time::Instant::now();
    let (ev4, _, res4) = lanczos_lowest_herm(&mv4, NS4, 4, 200, 999);
    let ms = t1.elapsed().as_millis();
    let spread4 = ev4[2] - ev4[0];
    let gap4 = ev4[3] - ev4[2];
    let ok3 = spread4.abs() < 1e-8 && (gap4 - 0.3195).abs() < 5e-4 && res4 < 1e-8;
    println!(
        "    幅 {:.2e} / ギャップ {:.4} (v12.1: 0.3195) / 残差 {:.2e}  {}  ({} ms — ヤコビは ~594000 ms)",
        spread4,
        gap4,
        res4,
        pass(ok3),
        ms
    );

    // ---- [4] 決定論と初期値独立性 ----
    println!("\n[4] 決定論 (同シード) と初期値独立性 (異シード)");
    let (ev4b, _, _) = lanczos_lowest_herm(&mv4, NS4, 4, 200, 999);
    let same = (0..4).all(|e| ev4[e] == ev4b[e]);
    let (ev4c, _, _) = lanczos_lowest_herm(&mv4, NS4, 4, 200, 31337);
    let mut max_d4: f64 = 0.0;
    for e in 0..4 {
        max_d4 = max_d4.max((ev4[e] - ev4c[e]).abs());
    }
    let ok4 = same && max_d4 < 1e-9;
    println!(
        "    同シード bit 一致: {} / 異シード max|Δλ| = {:.2e}  {}",
        same,
        max_d4,
        pass(ok4)
    );

    let all_ok = ok1 && ok2 && ok3 && ok4;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-TOOL-001".into())),
        ("dense_dev".into(), Json::Num(max_d1)),
        ("torus2d_gap".into(), Json::Num(gap2)),
        ("t4_gap".into(), Json::Num(gap4)),
        ("t4_ms".into(), Json::Int(ms as i64)),
        ("jacobi_ms_reference".into(), Json::Int(594000)),
        ("max_residual".into(), Json::Num(res1.max(res2).max(res4))),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v131_lanczos.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {}", pass(all_ok));
    println!("\n結論: Lanczos 疎ソルバ (完全再直交化・残差検証つき) が 3 つの物理系で");
    println!("      ヤコビ法・一次ソースと一致した。N=6 T⁴ の最低バンドが ~10 分 → 数秒。");
    println!("      N=18 の T⁴ (10 万サイト) が射程に入った — v13.2 (深い磁束族) へ。");
    if !all_ok {
        std::process::exit(1);
    }
}
