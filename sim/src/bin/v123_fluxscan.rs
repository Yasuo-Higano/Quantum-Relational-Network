//! v12.3 指数 3 の磁束族の系統走査 — 深い傾き T⁴ を探す
//!
//! v12.2 の陰性 (磁束 (2,2,1,−1) は lnZ₉ = −94.5) の敗因は階層の浅さだった。
//! 同じ N=6 格子で指数 Pf(F) = Q₁Q₂ + ts = 3 を保つ磁束は他にもある:
//!   (3,1,0,0) 因子化参照 / (3,2,3,−1) / (3,2,1,−3) / (3,3,2,−3) / (3,3,3,−2)
//! Q=3 のタワーは深い階層を作れる (T² 系での実績) — 傾きと深さは両立するか。
//!
//! 手順:
//!  [1] 各磁束を 1 配置 (Wilson 0) 対角化 (並列): 指数 3 の確認 + 深さ代理
//!      (σ_H 走査での最小質量比 r1) + トーラス間 MI
//!  [2] 深さで勝者を選ぶ
//!  [3] 勝者だけ全 36 Wilson 配置の証拠計算 (v12.2 と同一の機械) — はしごと比較
//!
//! 検証: 各磁束の縮退幅 < 1e-8 (指数 3)、(2,2,1,−1) の参照値が v12.2 と一致。

use uft_sim::*;


const N: usize = 6;
const NS4: usize = N * N * N * N;
type Flux = [i64; 4];
const NK: usize = 6; // Wilson 格子 Z₆ × Z₆
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const REF_GAP_221: f64 = 0.3195; // v12.1/v12.2 の (2,2,1,−1) ギャップ (回帰用)

type M3 = [[(f64, f64); 3]; 3];

fn idx(x1: usize, y1: usize, x2: usize, y2: usize) -> usize {
    x1 + N * (y1 + N * (x2 + N * y2))
}

fn link_phase(f: &Flux, x1: usize, y1: usize, x2: usize, y2: usize, dir: usize, w1: f64, w2: f64) -> f64 {
    let nn = (N * N) as f64;
    let two_pi = 2.0 * std::f64::consts::PI;
    let (p1, p2, pt, ps) = (
        two_pi * f[0] as f64 / nn,
        two_pi * f[1] as f64 / nn,
        two_pi * f[2] as f64 / nn,
        two_pi * f[3] as f64 / nn,
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
fn t4_modes_w(f: &Flux, w1: f64, w2: f64) -> (Vec<Vec<(f64, f64)>>, f64, f64) {
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
                    addhop(i, idx((x1 + 1) % N, y1, x2, y2), link_phase(f, x1, y1, x2, y2, 0, w1, w2));
                    addhop(i, idx(x1, (y1 + 1) % N, x2, y2), link_phase(f, x1, y1, x2, y2, 1, w1, w2));
                    addhop(i, idx(x1, y1, (x2 + 1) % N, y2), link_phase(f, x1, y1, x2, y2, 2, w1, w2));
                    addhop(i, idx(x1, y1, x2, (y2 + 1) % N), link_phase(f, x1, y1, x2, y2, 3, w1, w2));
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

/// 深さ代理: 1 配置の湯川 (cfg0 × cfg0) の最小対数質量比を σ_H 走査で取る
fn depth_proxy(modes: &[Vec<(f64, f64)>]) -> f64 {
    let mut best = f64::INFINITY;
    for &sh in &[0.8f64, 1.2, 1.6, 2.0] {
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
        let mut y = [[(0.0f64, 0.0f64); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (mut sr, mut si) = (0.0, 0.0);
                for st in 0..NS4 {
                    let (ar, ai) = modes[i][st];
                    let (br, bi) = modes[j][st];
                    sr += (ar * br + ai * bi) * phih[st];
                    si += (ar * bi - ai * br) * phih[st];
                }
                y[i][j] = (sr, si);
            }
        }
        let (r, _) = mass_and_vecs(&y);
        best = best.min(r[0]);
    }
    best // ln(σ1/σ3) の最小 (深いほど負)
}

/// トーラス間 MI (バンド平均, v12.1 と同じ)
fn band_mi(modes: &[Vec<(f64, f64)>]) -> f64 {
    let mut p = vec![0.0f64; N * N];
    for m in modes {
        for x1 in 0..N {
            for y1 in 0..N {
                for x2 in 0..N {
                    for y2 in 0..N {
                        let (r, i) = m[idx(x1, y1, x2, y2)];
                        p[x1 + x2 * N] += r * r + i * i;
                    }
                }
            }
        }
    }
    let tot: f64 = p.iter().sum();
    for v in p.iter_mut() {
        *v /= tot;
    }
    let mut p1 = vec![0.0f64; N];
    let mut p2 = vec![0.0f64; N];
    for x1 in 0..N {
        for x2 in 0..N {
            p1[x1] += p[x1 + x2 * N];
            p2[x2] += p[x1 + x2 * N];
        }
    }
    let mut mi = 0.0;
    for x1 in 0..N {
        for x2 in 0..N {
            let q = p[x1 + x2 * N];
            if q > 1e-300 {
                mi += q * (q / (p1[x1] * p2[x2])).ln();
            }
        }
    }
    mi
}

fn main() {
    self_test();
    println!("=== v12.3 指数 3 の磁束族の系統走査: 深い傾き T⁴ を探す ===\n");
    let family: Vec<(Flux, &str)> = vec![
        ([3, 1, 0, 0], "因子化参照 (3·1)"),
        ([2, 2, 1, -1], "v12.2 の点 (4−1)"),
        ([3, 2, 3, -1], "6−3"),
        ([3, 2, 1, -3], "6−3'"),
        ([3, 3, 2, -3], "9−6"),
        ([3, 3, 3, -2], "9−6'"),
    ];
    println!("[1] 各磁束の 1 配置対角化 (並列, ~15 分) — 指数・深さ代理・MI");
    let t0 = std::time::Instant::now();
    let handles: Vec<_> = family
        .iter()
        .map(|&(f, _)| std::thread::spawn(move || t4_modes_w(&f, 0.0, 0.0)))
        .collect();
    let outs: Vec<(Vec<Vec<(f64, f64)>>, f64, f64)> =
        handles.into_iter().map(|h| h.join().expect("worker")).collect();
    println!("    完了 ({} ms)\n", t0.elapsed().as_millis());
    println!("    磁束              指数(幅/ギャップ)         深さ代理 ln r1   MI(x₁:x₂)");
    let mut ok_engine = true;
    let mut best_i = 0usize;
    let mut best_depth = f64::INFINITY;
    let mut rows = Vec::new();
    for (i, ((f, tag), (modes, spread, gap))) in family.iter().zip(outs.iter()).enumerate() {
        let ok = *spread < 1e-8 && *gap > 0.02;
        ok_engine &= ok;
        let d = depth_proxy(modes);
        let mi = band_mi(modes);
        println!(
            "    {:?} {:10} {:.1e}/{:.4} {}   {:+.2}        {:.4}",
            f, tag, spread, gap, pass(ok), d, mi
        );
        if d < best_depth {
            best_depth = d;
            best_i = i;
        }
        rows.push((*f, d, mi, *spread, *gap));
    }
    // 回帰: (2,2,1,−1) のギャップ
    let ok_reg = (outs[1].2 - REF_GAP_221).abs() < 5e-4;
    println!(
        "\n    回帰: (2,2,1,−1) ギャップ {:.4} = v12.1/v12.2 の {:.4}  {}",
        outs[1].2,
        REF_GAP_221,
        pass(ok_reg)
    );
    let (best_f, _) = family[best_i];
    println!(
        "\n[2] 深さの勝者: {:?} (深さ代理 ln r1 = {:+.2})",
        best_f, best_depth
    );

    // ---- [3] 勝者の全 36 Wilson 配置の証拠 ----
    println!("\n[3] 勝者の全 36 Wilson 配置の証拠 (v12.2 と同一の機械, ~45 分)");
    let nw = 12usize.min(std::thread::available_parallelism().map(|v| v.get()).unwrap_or(4));
    let two_pi = 2.0 * std::f64::consts::PI;
    let t1 = std::time::Instant::now();
    let mut results: Vec<Option<(Vec<Vec<(f64, f64)>>, f64, f64)>> = (0..36).map(|_| None).collect();
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
                    (k, t4_modes_w(&best_f, w1, w2))
                })
            })
            .collect();
        for h in handles {
            let (k, r) = h.join().expect("worker panic");
            results[k] = Some(r);
        }
        next += nw;
        println!("    … {}/36 完了 ({} ms)", next.min(36), t1.elapsed().as_millis());
    }
    let configs: Vec<(Vec<Vec<(f64, f64)>>, f64, f64)> =
        results.into_iter().map(|r| r.unwrap()).collect();
    let max_spread = configs.iter().map(|c| c.1).fold(0.0f64, f64::max);
    let ok_idx36 = max_spread < 1e-8;
    println!(
        "    全 36 配置の縮退幅 max {:.2e} (< 1e-8)  {}",
        max_spread,
        pass(ok_idx36)
    );

    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let sig_grid = [0.8f64, 1.2, 1.6, 2.0];
    let nc = 36usize;
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
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
                        for st in 0..NS4 {
                            let (ar, ai) = configs[a].0[i][st];
                            let (br, bi) = configs[b].0[j][st];
                            sr += (ar * br + ai * bi) * phih[st];
                            si += (ar * bi - ai * br) * phih[st];
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
    println!("\n    lnZ₉({:?}) = {:.4}", best_f, lnz);
    println!("\n    機構のはしご: 対角 −23.61 / 一様 −21.83 / 向き −20.86 / S₃ 対 −19.86 / (2,2,1,−1) −94.51");
    let beats = lnz > -19.86;
    println!(
        "    => 勝者磁束は S₃ 対を{} (差 {:+.2})",
        if beats {
            "上回る — 深い傾きが対の問いを解いた"
        } else {
            "上回らない — N=6 窓の傾き族では階層が届かない (正直な境界)"
        },
        lnz - (-19.86)
    );

    let all_ok = ok_engine && ok_reg && ok_idx36;
    let mut frows = Vec::new();
    for (f, d, mi, spread, gap) in &rows {
        frows.push(Json::Obj(vec![
            ("flux".into(), Json::Arr(f.iter().map(|&x| Json::Int(x)).collect())),
            ("depth_ln_r1".into(), Json::Num(*d)),
            ("band_mi".into(), Json::Num(*mi)),
            ("spread".into(), Json::Num(*spread)),
            ("gap".into(), Json::Num(*gap)),
        ]));
    }
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-014".into())),
        ("family".into(), Json::Arr(frows)),
        ("winner".into(), Json::Arr(best_f.iter().map(|&x| Json::Int(x)).collect())),
        ("winner_lnZ_nine".into(), Json::Num(lnz)),
        ("beats_s3_ladder".into(), Json::Bool(beats)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v123_fluxscan.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {} (PASS = 装置 — 走査の答えは [1]-[3] の表が本体)", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
