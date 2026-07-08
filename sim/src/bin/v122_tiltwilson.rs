//! v12.2 傾き T⁴ の証拠検定 — アンザッツなしの 3 世代は機構のはしごを超えるか
//!
//! v12.1 は傾き磁束 (2,2,1,−1) の T⁴ が指数 = 3 の厳密 3 重縮退を持つことを実証した。
//! 本バイナリはその 9 量証拠を計算し、機構のはしご (T²×T² 系, v11.3):
//!   対角対 −23.61 / 一様オービフォールド −21.83 / 向き+並進 −20.86 / S₃ 対 −19.86
//! と比較する。
//!
//! **この構成の原理的な利点**: 3 つのゼロモードがそのまま 3 世代であり、
//! 「9 個から選ぶ射影」も「トーラス間の対」も存在しない。したがって
//! v9.2 のラベル綱渡りも v10.1 の σ も**原理的に発生しない** — 特異値と |CKM| は
//! モードの順序・位相の取り方に完全に不変で、規約が入り込む場所がない。
//!
//! Wilson 線: セクターごとに定数ホロノミー (w_{y1}, w_{y2}) = 2π(k₁,k₂)/6, k ∈ Z₆²
//! (36 配置/場 — はしごの模型と同じ事前体積 5·ln36)。傾きがあると並進では代用できない
//! (v12.1 の開発記録) ので、**36 配置を個別に対角化**する (std::thread で並列)。
//!
//! 検証: 全 36 配置で縮退幅 < 1e-8 (指数 3 が Wilson 線に依らない)、
//!       (0,0) 配置のギャップが v12.1 と一致 (回帰)。

use uft_sim::*;

const N: usize = 6;
const NS4: usize = N * N * N * N;
const FLUX: [i64; 4] = [2, 2, 1, -1];
const NK: usize = 6; // Wilson 格子 Z₆ × Z₆
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const REF_GAP_00: f64 = 0.3195; // v12.1 (results/v121_tiltflux.json tilt_gap, 表示 4 桁)

type M3 = [[(f64, f64); 3]; 3];

fn idx(x1: usize, y1: usize, x2: usize, y2: usize) -> usize {
    x1 + N * (y1 + N * (x2 + N * y2))
}

fn link_phase(x1: usize, y1: usize, x2: usize, y2: usize, dir: usize, w1: f64, w2: f64) -> f64 {
    let nn = (N * N) as f64;
    let two_pi = 2.0 * std::f64::consts::PI;
    let (p1, p2, pt, ps) = (
        two_pi * FLUX[0] as f64 / nn,
        two_pi * FLUX[1] as f64 / nn,
        two_pi * FLUX[2] as f64 / nn,
        two_pi * FLUX[3] as f64 / nn,
    );
    let nf = N as f64;
    match dir {
        0 => {
            if x1 == N - 1 {
                -(p1 * nf * y1 as f64 + pt * nf * y2 as f64)
            } else {
                0.0
            }
        }
        1 => {
            let mut th = p1 * x1 as f64 + w1;
            if y1 == N - 1 {
                th += -(ps * nf * x2 as f64);
            }
            th
        }
        2 => {
            let mut th = ps * y1 as f64;
            if x2 == N - 1 {
                th += -(p2 * nf * y2 as f64);
            }
            th
        }
        3 => p2 * x2 as f64 + pt * x1 as f64 + w2,
        _ => unreachable!(),
    }
}

/// Wilson (w1,w2) 付きで最低 3 複素モードを返す
fn t4_modes_w(w1: f64, w2: f64) -> (Vec<Vec<(f64, f64)>>, f64, f64) {
    let m = 2 * NS4;
    let mut a = vec![0.0f64; m * m];
    let mut addhop = |i: usize, j: usize, th: f64| {
        let (c, s) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + NS4) + (i + NS4) * m] += -c;
        a[(i + NS4) + (j + NS4) * m] += -c;
        a[j + (i + NS4) * m] += s;
        a[(j + NS4) + i * m] += -s;
        a[i + (j + NS4) * m] += -s;
        a[(i + NS4) + j * m] += s;
    };
    for x1 in 0..N {
        for y1 in 0..N {
            for x2 in 0..N {
                for y2 in 0..N {
                    let i = idx(x1, y1, x2, y2);
                    addhop(
                        i,
                        idx((x1 + 1) % N, y1, x2, y2),
                        link_phase(x1, y1, x2, y2, 0, w1, w2),
                    );
                    addhop(
                        i,
                        idx(x1, (y1 + 1) % N, x2, y2),
                        link_phase(x1, y1, x2, y2, 1, w1, w2),
                    );
                    addhop(
                        i,
                        idx(x1, y1, (x2 + 1) % N, y2),
                        link_phase(x1, y1, x2, y2, 2, w1, w2),
                    );
                    addhop(
                        i,
                        idx(x1, y1, x2, (y2 + 1) % N),
                        link_phase(x1, y1, x2, y2, 3, w1, w2),
                    );
                }
            }
        }
    }
    let (w, v) = jacobi_eigh(&a, m);
    let spread = w[5] - w[0];
    let gap = w[6] - w[5];
    let mut modes: Vec<Vec<(f64, f64)>> = Vec::new();
    for kk in 0..6 {
        let mut psi: Vec<(f64, f64)> = (0..NS4)
            .map(|i| (v[i + kk * m], v[(i + NS4) + kk * m]))
            .collect();
        for pm in &modes {
            let (mut pr, mut pi) = (0.0, 0.0);
            for i in 0..NS4 {
                pr += pm[i].0 * psi[i].0 + pm[i].1 * psi[i].1;
                pi += pm[i].0 * psi[i].1 - pm[i].1 * psi[i].0;
            }
            for i in 0..NS4 {
                let (ar, ai) = pm[i];
                psi[i].0 -= pr * ar - pi * ai;
                psi[i].1 -= pr * ai + pi * ar;
            }
        }
        let nrm: f64 = psi.iter().map(|&(r, i)| r * r + i * i).sum::<f64>().sqrt();
        if nrm > 1e-6 {
            for p in psi.iter_mut() {
                p.0 /= nrm;
                p.1 /= nrm;
            }
            modes.push(psi);
            if modes.len() == 3 {
                break;
            }
        }
    }
    (modes, spread, gap)
}

fn gram(y: &M3) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let mut hre = [[0.0f64; 3]; 3];
    let mut him = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let (a, b) = y[i][k];
                let (c, d) = y[j][k];
                hre[i][j] += a * c + b * d;
                him[i][j] += b * c - a * d;
            }
        }
    }
    (hre, him)
}

fn eig_herm3(hre: &[[f64; 3]; 3], him: &[[f64; 3]; 3]) -> ([f64; 3], M3) {
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
        let nrm: f64 = vecs[k]
            .iter()
            .map(|&(a, b)| a * a + b * b)
            .sum::<f64>()
            .sqrt();
        for i in 0..3 {
            vecs[k][i].0 /= nrm;
            vecs[k][i].1 /= nrm;
        }
    }
    (lam, vecs)
}

fn mass_and_vecs(y: &M3) -> ([f64; 2], M3) {
    let (hre, him) = gram(y);
    let (lam, vecs) = eig_herm3(&hre, &him);
    let sv = [
        lam[0].max(0.0).sqrt(),
        lam[1].max(0.0).sqrt(),
        lam[2].max(0.0).sqrt(),
    ];
    (
        [
            (sv[0].max(1e-300) / sv[2].max(1e-300)).ln(),
            (sv[1].max(1e-300) / sv[2].max(1e-300)).ln(),
        ],
        vecs,
    )
}

fn ckm3(vu: &M3, vd: &M3) -> [f64; 3] {
    let mut ckm = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = vu[i][k];
                let (c, d) = vd[j][k];
                re += a * c + b * d;
                im += a * d - b * c;
            }
            ckm[i][j] = (re * re + im * im).sqrt();
        }
    }
    [ckm[0][1], ckm[1][2], ckm[0][2]]
}

fn lse(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

fn main() {
    self_test();
    println!("=== v12.2 傾き T⁴ の証拠検定: アンザッツなしの 3 世代 vs 機構のはしご ===\n");
    let nw = 12usize.min(
        std::thread::available_parallelism()
            .map(|v| v.get())
            .unwrap_or(4),
    );
    println!(
        "[1] Wilson 36 配置の個別対角化 (2592², {} スレッド並列, ~30-60 分)",
        nw
    );
    let t0 = std::time::Instant::now();
    let two_pi = 2.0 * std::f64::consts::PI;
    // ワーカー: 配置 k = k1 + 6·k2 を処理
    let mut results: Vec<Option<(Vec<Vec<(f64, f64)>>, f64, f64)>> =
        (0..36).map(|_| None).collect();
    let mut next = 0usize;
    while next < 36 {
        let batch: Vec<usize> = (next..(next + nw).min(36)).collect();
        let handles: Vec<_> = batch
            .iter()
            .map(|&k| {
                std::thread::spawn(move || {
                    let (k1, k2) = (k % NK, k / NK);
                    let w1 = two_pi * k1 as f64 / NK as f64;
                    let w2 = two_pi * k2 as f64 / NK as f64;
                    (k, t4_modes_w(w1, w2))
                })
            })
            .collect();
        for h in handles {
            let (k, r) = h.join().expect("worker panic");
            results[k] = Some(r);
        }
        next += nw;
        println!(
            "    … {}/36 完了 ({} ms)",
            next.min(36),
            t0.elapsed().as_millis()
        );
    }
    let configs: Vec<(Vec<Vec<(f64, f64)>>, f64, f64)> =
        results.into_iter().map(|r| r.unwrap()).collect();
    let max_spread = configs.iter().map(|c| c.1).fold(0.0f64, f64::max);
    let min_gap = configs.iter().map(|c| c.2).fold(f64::INFINITY, f64::min);
    let ok_idx = max_spread < 1e-8 && min_gap > 0.02;
    println!(
        "    全 36 配置: 縮退幅 max {:.2e} / ギャップ min {:.4} — 指数 3 は Wilson に依らない  {}",
        max_spread,
        min_gap,
        pass(ok_idx)
    );
    let ok_reg = (configs[0].2 - REF_GAP_00).abs() < 5e-4;
    println!(
        "    (0,0) 配置のギャップ {:.4} = v12.1 の {:.4}  {}",
        configs[0].2,
        REF_GAP_00,
        pass(ok_reg)
    );

    // ---- [2] 9 量の証拠 ----
    println!("\n[2] 9 量の証拠 (事前: Wilson 36/場 × 5 場 + σ_H 4 点 — はしごと同体積)");
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let sig_grid = [0.8f64, 1.2, 1.6, 2.0];
    let nc = 36usize;
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let t2 = std::time::Instant::now();
    let mut terms = Vec::new();
    for &sh in &sig_grid {
        let mut phih = vec![0.0f64; NS4];
        for x1 in 0..N {
            for y1 in 0..N {
                for x2 in 0..N {
                    for y2 in 0..N {
                        let d = |a: usize| {
                            let v = a as f64;
                            v.min(N as f64 - v)
                        };
                        let r2 = d(x1) * d(x1) + d(y1) * d(y1) + d(x2) * d(x2) + d(y2) * d(y2);
                        phih[idx(x1, y1, x2, y2)] = (-r2 / (2.0 * sh * sh)).exp();
                    }
                }
            }
        }
        let pair: Vec<([f64; 2], M3)> = (0..nc * nc)
            .map(|ab| {
                let (a, b) = (ab % nc, ab / nc);
                let mut y = [[(0.0f64, 0.0f64); 3]; 3];
                for i in 0..3 {
                    for j in 0..3 {
                        let (mut sr, mut si) = (0.0, 0.0);
                        for s in 0..NS4 {
                            let (ar, ai) = configs[a].0[i][s];
                            let (br, bi) = configs[b].0[j][s];
                            sr += (ar * br + ai * bi) * phih[s];
                            si += (ar * bi - ai * br) * phih[s];
                        }
                        y[i][j] = (sr, si);
                    }
                }
                mass_and_vecs(&y)
            })
            .collect();
        let le: Vec<f64> = pair.iter().map(|(r, _)| ll2(r, tgt[4], tgt[5])).collect();
        let mut per_q = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut acc = (f64::NEG_INFINITY, 0.0f64);
            for ku in 0..nc {
                let (ru, vu) = &pair[kq + ku * nc];
                let llu = ll2(ru, tgt[0], tgt[1]);
                for kd in 0..nc {
                    let (rd, vd) = &pair[kq + kd * nc];
                    let lld = ll2(rd, tgt[2], tgt[3]);
                    let c = ckm3(vu, vd);
                    let mut ll = llu + lld;
                    for m in 0..3 {
                        let d = c[m].max(1e-300).ln() - tgt[6 + m];
                        ll += -d * d / (2.0 * sigma * sigma) + norm1;
                    }
                    if ll > acc.0 {
                        acc.1 = acc.1 * (acc.0 - ll).exp() + 1.0;
                        acc.0 = ll;
                    } else {
                        acc.1 += (ll - acc.0).exp();
                    }
                }
            }
            per_q.push(acc.0 + acc.1.ln());
        }
        terms.push(lse(&per_q) + lse(&le));
    }
    let lnz = lse(&terms) - (5.0 * (nc as f64).ln() + (sig_grid.len() as f64).ln());
    println!(
        "    lnZ₉(傾き T⁴) = {:.4}  ({} ms)",
        lnz,
        t2.elapsed().as_millis()
    );
    println!("\n    機構のはしご (T²×T² 系):");
    println!("      対角対      −23.61");
    println!("      一様オービ  −21.83");
    println!("      向き+並進   −20.86");
    println!("      S₃ 対       −19.86");
    println!(
        "      傾き T⁴     {:+.2}   ← 本バイナリ (アンザッツ・ラベル・対が原理的に不在)",
        lnz
    );
    let beats = lnz > -19.86;
    println!(
        "\n    => 傾き T⁴ は S₃ 対を{} (差 {:+.2})",
        if beats {
            "上回る — 「σ は因子化近似の影」が証拠でも支持された"
        } else {
            "上回らない — 指数 3 は実現するが、この磁束・格子では階層が足りない"
        },
        lnz - (-19.86)
    );

    let all_ok = ok_idx && ok_reg;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-013".into())),
        ("configs".into(), Json::Int(36)),
        ("max_spread".into(), Json::Num(max_spread)),
        ("min_gap".into(), Json::Num(min_gap)),
        ("lnZ_nine".into(), Json::Num(lnz)),
        ("beats_s3_ladder".into(), Json::Bool(beats)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v122_tiltwilson.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!(
        "\n総合判定: {} (PASS = 指数の Wilson 不変性と回帰 — はしご比較は [2] が本体)",
        pass(all_ok)
    );
    if !all_ok {
        std::process::exit(1);
    }
}
