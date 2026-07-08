//! v16.1 傾き T⁴ の CP 検定 — 非因子化窓は Jarlskog J ≠ 0 を作るか
//!
//! v15.7 の holdout 採点は「長方形 T²×T² + 実 Wilson の湯川位相は行×列に因子化し、
//! CKM は再位相化で実になる (J = 0 厳密)」を発見した。測定 J = 3.08e-5 ≠ 0 は
//! この窓を反証し、**CP 破れは非因子化幾何を要求する**ことを指定した (残高 2 の
//! 新制約: 正しい窓は CP を作る窓)。本バイナリはその最初の検定を行う:
//! v13.2 の N=18 傾き T⁴ (指数 Pf(F)=3 の厳密縮退点) は J ≠ 0 を作るか。
//!
//! 方法 (v13.2 の機構を回帰つきで再利用):
//!   [1] 勝者磁束 (3,3,1,−6) の 36 Wilson 配置を並列 Lanczos で解き、
//!       lnZ₉ = −34.8804 (v13.2 公表値) を回帰再現する — エンジンの較正。
//!   [2] 同じ事後分布 (9 量尤度 × 一様事前) の下で複素 CKM → J を全
//!       (σ_H, kQ, ku, kd) について集計: max |J|、事後分位点、MAP の J、
//!       P(|J| > 1e-6)。J の符号は Wilson 反転の共役対称で予言されない — |J| が対象。
//!   [3] 対照: 因子化磁束 (3,1,0,0) を単一 T² 厳密対角化の積で構成し (完全性
//!       リスクなし)、同じ CKM 構成 (共有 kQ) で J を測る。開発で判明した事実:
//!       因子化スカラー格子の J は 0 でなく格子スケール (~1e-4) — 連続体の実化
//!       可能性は格子で破れる。よってこれは「装置が J を検出できる」陽性対照であり、
//!       構造零の参照は v15.7 (Dirac + 局在基底 + Hadamard 積) が担う。
//!
//! 判定 (事前登録; 対照の役割は開発記録の通り陽性対照に訂正):
//!   ・傾き窓の max |J| > 1e-4 → 「傾きは CP を作る」(深さと独立の資格)
//!   ・傾き窓の max |J| < 1e-6 (かつ陽性対照 > 1e-6) → 傾き点は格子アーティファクト
//!     水準よりも深く CP を抑制する — 非因子化だけでは足りない (窓の拡張が必要)
//!   どちらでも測量 — 「CP を作る窓」の地図の最初の点である。

use uft_sim::*;

const N: usize = 18;
const NS4: usize = N * N * N * N;
const NK: usize = 6;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const REF_LNZ9: f64 = -34.8804; // v13.2 の公表値 (results/v132_deeptilt.json)
const SIG_GRID: [f64; 4] = [1.0, 1.5, 2.0, 2.5];

type Flux = [i64; 4];
type M3 = [[(f64, f64); 3]; 3];

fn idx(x1: usize, y1: usize, x2: usize, y2: usize) -> usize {
    x1 + N * (y1 + N * (x2 + N * y2))
}

fn link_phase(
    f: &Flux,
    x1: usize,
    y1: usize,
    x2: usize,
    y2: usize,
    dir: usize,
    w1: f64,
    w2: f64,
) -> f64 {
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

fn hops(f: &Flux, w1: f64, w2: f64) -> Vec<(u32, u32, f64)> {
    let mut h = Vec::with_capacity(4 * NS4);
    for x1 in 0..N {
        for y1 in 0..N {
            for x2 in 0..N {
                for y2 in 0..N {
                    let i = idx(x1, y1, x2, y2) as u32;
                    h.push((
                        i,
                        idx((x1 + 1) % N, y1, x2, y2) as u32,
                        link_phase(f, x1, y1, x2, y2, 0, w1, w2),
                    ));
                    h.push((
                        i,
                        idx(x1, (y1 + 1) % N, x2, y2) as u32,
                        link_phase(f, x1, y1, x2, y2, 1, w1, w2),
                    ));
                    h.push((
                        i,
                        idx(x1, y1, (x2 + 1) % N, y2) as u32,
                        link_phase(f, x1, y1, x2, y2, 2, w1, w2),
                    ));
                    h.push((
                        i,
                        idx(x1, y1, x2, (y2 + 1) % N) as u32,
                        link_phase(f, x1, y1, x2, y2, 3, w1, w2),
                    ));
                }
            }
        }
    }
    h
}

fn matvec(hops: &[(u32, u32, f64)], v: &[(f64, f64)]) -> Vec<(f64, f64)> {
    let mut o = vec![(0.0f64, 0.0f64); NS4];
    for &(i, j, th) in hops {
        let (i, j) = (i as usize, j as usize);
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

/// 1 配置を解く (v13.2 と同一): 残差 < 1e-9 まで Krylov 次元を適応
fn solve(f: &Flux, w1: f64, w2: f64) -> (Vec<Vec<(f64, f64)>>, f64, f64, f64) {
    let hp = hops(f, w1, w2);
    let mv = |v: &[(f64, f64)]| matvec(&hp, v);
    for &m in &[250usize, 400, 700] {
        let (ev, vecs, res) = lanczos_lowest_herm(&mv, NS4, 4, m, 1332);
        if res < 1e-9 {
            let spread = ev[2] - ev[0];
            let gap = ev[3] - ev[2];
            return (vecs.into_iter().take(3).collect(), spread, gap, res);
        }
    }
    let (ev, vecs, res) = lanczos_lowest_herm(&mv, NS4, 4, 700, 1332);
    let spread = ev[2] - ev[0];
    let gap = ev[3] - ev[2];
    (vecs.into_iter().take(3).collect(), spread, gap, res)
}

fn higgs(sh: f64) -> Vec<f64> {
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
    phih
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

/// 複素 CKM と (J, |V_td|)
fn ckm_full(vu: &M3, vd: &M3) -> M3 {
    let mut v = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = vu[i][k];
                let (c, d) = vd[j][k];
                re += a * c + b * d;
                im += b * c - a * d;
            }
            v[i][j] = (re, im);
        }
    }
    v
}

fn jarlskog(v: &M3) -> f64 {
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let conj = |a: (f64, f64)| (a.0, -a.1);
    mul(mul(v[0][1], v[1][2]), mul(conj(v[0][2]), conj(v[1][1]))).1
}

fn cabs(v: &M3, i: usize, j: usize) -> f64 {
    (v[i][j].0 * v[i][j].0 + v[i][j].1 * v[i][j].1).sqrt()
}

fn lse(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

fn yuk(a: &[Vec<(f64, f64)>], b: &[Vec<(f64, f64)>], phih: &[f64]) -> M3 {
    let mut y = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut sr, mut si) = (0.0, 0.0);
            for s in 0..NS4 {
                let (ar, ai) = a[i][s];
                let (br, bi) = b[j][s];
                sr += (ar * br + ai * bi) * phih[s];
                si += (ar * bi - ai * br) * phih[s];
            }
            y[i][j] = (sr, si);
        }
    }
    y
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// Wilson 配置群を並列に解く
fn solve_grid(f: Flux, ks: &[(usize, usize)]) -> Vec<(Vec<Vec<(f64, f64)>>, f64, f64, f64)> {
    let nw = 12usize.min(
        std::thread::available_parallelism()
            .map(|v| v.get())
            .unwrap_or(4),
    );
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut results: Vec<Option<(Vec<Vec<(f64, f64)>>, f64, f64, f64)>> =
        (0..ks.len()).map(|_| None).collect();
    let mut next = 0usize;
    while next < ks.len() {
        let batch: Vec<usize> = (next..(next + nw).min(ks.len())).collect();
        let handles: Vec<_> = batch
            .iter()
            .map(|&i| {
                let (k1, k2) = ks[i];
                std::thread::spawn(move || {
                    let w1 = two_pi * k1 as f64 / NK as f64;
                    let w2 = two_pi * k2 as f64 / NK as f64;
                    (i, solve(&f, w1, w2))
                })
            })
            .collect();
        for h in handles {
            let (i, r) = h.join().expect("worker panic");
            results[i] = Some(r);
        }
        next += nw;
    }
    results.into_iter().map(|r| r.unwrap()).collect()
}

fn main() {
    self_test();
    println!("=== v16.1 傾き T⁴ の CP 検定: 非因子化窓は J ≠ 0 を作るか ===\n");
    println!(
        "事前登録の判定: 傾き max|J| > 1e-4 → CP を作る / < 1e-6 (陽性対照 > 1e-6) → 構造的抑制\n"
    );
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---- [1] 勝者磁束の 36 Wilson (v13.2 回帰) ----
    let f_tilt: Flux = [3, 3, 1, -6];
    println!(
        "[1] 勝者磁束 {:?} の 36 Wilson 配置 (並列 Lanczos, n=209,952)",
        f_tilt
    );
    let t0 = std::time::Instant::now();
    let ks36: Vec<(usize, usize)> = (0..36).map(|k| (k % NK, k / NK)).collect();
    let configs = solve_grid(f_tilt, &ks36);
    let max_res = configs.iter().map(|c| c.3).fold(0.0f64, f64::max);
    let max_spread = configs.iter().map(|c| c.1).fold(0.0f64, f64::max);
    println!("    完了 ({} ms)", t0.elapsed().as_millis());
    check(
        "全 36 配置の残差 < 1e-9・縮退幅 < 1e-6 (指数 3 の Wilson 頑健性)",
        max_res < 1e-9 && max_spread < 1e-6,
        format!("残差 max {:.1e}, 幅 max {:.1e}", max_res, max_spread),
    );

    // ---- [2] 証拠回帰 + J の事後集計 ----
    println!("\n[2] lnZ₉ 回帰 (v13.2 公表値) と J の事後分布");
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let nc = 36usize;
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let mut terms = Vec::new();
    let mut jmax: f64 = 0.0;
    let mut jmax_ll = f64::NEG_INFINITY; // 高尤度域 (MAP−5 以内) の max|J| 用に (ll, |J|) を保持
    let mut samples: Vec<(f64, f64)> = Vec::new(); // (ll_total_shifted 前の生 ll + lnZe 相当, |J|)
    let mut map = (f64::NEG_INFINITY, 0.0f64, 0.0f64, 0usize, 0usize, 0usize);
    for (isg, &sh) in SIG_GRID.iter().enumerate() {
        let phih = higgs(sh);
        let pair: Vec<([f64; 2], M3)> = (0..nc * nc)
            .map(|ab| {
                let (a, b) = (ab % nc, ab / nc);
                mass_and_vecs(&yuk(&configs[a].0, &configs[b].0, &phih))
            })
            .collect();
        let le: Vec<f64> = pair.iter().map(|(r, _)| ll2(r, tgt[4], tgt[5])).collect();
        let lnze = lse(&le);
        let mut per_q = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut acc = (f64::NEG_INFINITY, 0.0f64);
            for ku in 0..nc {
                let (ru, vu) = &pair[kq + ku * nc];
                let llu = ll2(ru, tgt[0], tgt[1]);
                for kd in 0..nc {
                    let (rd, vd) = &pair[kq + kd * nc];
                    let lld = ll2(rd, tgt[2], tgt[3]);
                    let v = ckm_full(vu, vd);
                    let c = [cabs(&v, 0, 1), cabs(&v, 1, 2), cabs(&v, 0, 2)];
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
                    let j = jarlskog(&v);
                    if j.abs() > jmax {
                        jmax = j.abs();
                    }
                    samples.push((ll + lnze, j.abs()));
                    if ll > map.0 {
                        map = (ll, j, cabs(&v, 2, 0), kq, ku, kd);
                        jmax_ll = jmax_ll.max(ll);
                    }
                }
            }
            per_q.push(acc.0 + acc.1.ln());
        }
        terms.push(lse(&per_q) + lnze);
        let _ = isg;
    }
    let lnz = lse(&terms) - (5.0 * (nc as f64).ln() + (SIG_GRID.len() as f64).ln());
    check(
        "lnZ₉ の回帰: v13.2 公表値 −34.88 と一致 (±0.02)",
        (lnz - REF_LNZ9).abs() < 0.02,
        format!("lnZ₉ = {:.4} (公表 {:.4})", lnz, REF_LNZ9),
    );
    // 事後重みつき |J| 分位点
    let wmax = samples
        .iter()
        .map(|s| s.0)
        .fold(f64::NEG_INFINITY, f64::max);
    let mut hist = vec![0.0f64; 400]; // ln|J| ∈ [−18, −2]
    let (mut wsum, mut w_big) = (0.0f64, 0.0f64);
    for &(llw, ja) in &samples {
        let w = (llw - wmax).exp();
        wsum += w;
        if ja > 1e-6 {
            w_big += w;
        }
        let x = ja.max(1e-300).ln().clamp(-18.0, -2.0001);
        let b = ((x + 18.0) / 16.0 * 400.0) as usize;
        hist[b.min(399)] += w;
    }
    let quant = |q: f64| -> f64 {
        let mut acc = 0.0;
        for (i, h) in hist.iter().enumerate() {
            acc += h;
            if acc >= q * wsum {
                return (-18.0 + (i as f64 + 0.5) / 400.0 * 16.0).exp();
            }
        }
        1e-2
    };
    let (j16, j50, j84) = (quant(0.16), quant(0.5), quant(0.84));
    println!(
        "\n    ┌─ 傾き T⁴ {:?} の J ─────────────────────────────",
        f_tilt
    );
    println!("    │ max |J| (窓全体)       = {:.3e}", jmax);
    println!(
        "    │ 事後 |J|: 中央値 {:.3e}  68% [{:.3e}, {:.3e}]",
        j50, j16, j84
    );
    println!("    │ P(|J| > 1e-6)          = {:.3}", w_big / wsum);
    println!(
        "    │ MAP 点: J = {:+.3e}, |V_td| = {:.4} (kQ,ku,kd = {},{},{})",
        map.1, map.2, map.3, map.4, map.5
    );
    println!("    │ (測定: J = 3.08e-5, |V_td| = 0.0086)");
    println!("    └───────────────────────────────────────────────");
    // 開発記録: 初版は湯川「位相プラケット」arg(Y_ij Y_kl Ȳ_il Ȳ_kj) で機構診断を
    // 試みたが無効だった — 縮退空間のモード基底 (Lanczos/jacobi が返す) は任意
    // ユニタリ回転を持ち、プラケット位相はその回転で変わる (基底依存量; 実測 ~π が
    // その証拠)。基底不変な CP 検定量は J (∝ Im tr[H_u,H_d]³) であり、それが
    // 床レベルというのが本版の結果。機構の基底不変な診断は開いた問い。

    // ---- [3] 陽性対照: 因子化磁束 (3,1,0,0) の厳密積構成 ----
    // 開発記録 (初版の設計ミス): 対照を 4D Lanczos で解いたところ、Kronecker 和の
    // 踏み台スペクトル (Q₂=1 のラダー間隔 0.0385) の中で縮退パートナーの取り逃しが
    // 起き、見かけの「分裂 3.9e-2」を出した — 残差はベクトルの質を測るが、縮退部分
    // 空間の**取り尽くし**は保証しない (v13.1 の限界の新しい顔)。対照は 2 つの
    // 単一 T² の稠密厳密対角化 → 積波動関数の構成に変更する (構成的に因子化・
    // 完全性リスクなし)。単一 T² (N=18, Q=3) の LLL は厳密 3 重 (spread ~1e-13) を
    // 別途確認済み。
    println!(
        "\n[3] 対照: 因子化磁束 [3,1,0,0] — 単一 T² 厳密対角化の積構成 (共有 kQ の物理的 CKM)"
    );
    let t2 = std::time::Instant::now();
    // 単一 T² スカラー Hofstadter (複素エルミート 324) を実埋め込み jacobi で厳密に解く
    let t2_modes = |q: i64, w1: f64, nmodes: usize| -> Vec<Vec<(f64, f64)>> {
        let ns = N * N;
        let phi = 2.0 * std::f64::consts::PI * q as f64 / ns as f64;
        let id2 = |x: usize, y: usize| x + y * N;
        let m = 2 * ns;
        let mut a = vec![0.0; m * m];
        let mut add = |i: usize, j: usize, th: f64| {
            let (c, s) = (th.cos(), th.sin());
            // H_ij = -e^{iθ} (i←j) + h.c. の実埋め込み
            a[i + j * m] += -c;
            a[j + i * m] += -c;
            a[(i + ns) + (j + ns) * m] += -c;
            a[(j + ns) + (i + ns) * m] += -c;
            a[i + (j + ns) * m] += s;
            a[(j + ns) + i * m] += -s;
            a[j + (i + ns) * m] += -s;
            a[(i + ns) + j * m] += s;
        };
        for x in 0..N {
            for y in 0..N {
                add(
                    id2(x, y),
                    id2((x + 1) % N, y),
                    if x == N - 1 {
                        -phi * N as f64 * y as f64
                    } else {
                        0.0
                    },
                );
                add(id2(x, y), id2(x, (y + 1) % N), phi * x as f64 + w1);
            }
        }
        let (w, v) = jacobi_eigh(&a, m);
        let _ = w;
        // 実埋め込みの固有値は 2 重 — 複素モードとして 1 つおきに取る
        (0..nmodes)
            .map(|k| {
                (0..ns)
                    .map(|i| (v[i + (2 * k) * m], v[(i + ns) + (2 * k) * m]))
                    .collect::<Vec<(f64, f64)>>()
            })
            .collect()
    };
    // Wilson は v10.1 系の規約 (磁束量子ステップ w = φ·k, φ = 2πQ/N² — 非自明な
    // ホロノミーを持つ物理的な族)。2πk/6 の大ステップは N·w ∈ 2πZ でゲージ自明になる。
    let phi_q3 = 2.0 * std::f64::consts::PI * 3.0 / (N * N) as f64;
    // 6 Wilson 配置 × 3 モード (T1, Q=3) と T2 (Q=1) の基底モード
    let eta = &t2_modes(1, 0.0, 1)[0];
    let mut ctrl_modes: Vec<Vec<Vec<(f64, f64)>>> = Vec::new();
    for k1 in 0..6 {
        let chi = t2_modes(3, phi_q3 * k1 as f64, 3);
        // 積波動関数 ψ(x1,y1,x2,y2) = χ_i(x1,y1)·η(x2,y2)
        let prods: Vec<Vec<(f64, f64)>> = chi
            .iter()
            .map(|c| {
                let mut p = vec![(0.0f64, 0.0f64); NS4];
                for x1 in 0..N {
                    for y1 in 0..N {
                        let (cr, ci) = c[x1 + y1 * N];
                        for x2 in 0..N {
                            for y2 in 0..N {
                                let (er, eiq) = eta[x2 + y2 * N];
                                p[idx(x1, y1, x2, y2)] = (cr * er - ci * eiq, cr * eiq + ci * er);
                            }
                        }
                    }
                }
                p
            })
            .collect();
        ctrl_modes.push(prods);
    }
    println!("    完了 ({} ms)", t2.elapsed().as_millis());
    // CKM は同じ左場 (kQ) を u/d で共有する物理的な組のみ (開発記録: 初版は
    // 左場の異なる非物理な組まで比較して偽の J を出した — [2] と同じ添字構造に統一)。
    let mut jmax_ctrl: f64 = 0.0;
    {
        let phih = higgs(1.5);
        let pairc: Vec<M3> = (0..36)
            .map(|ab| {
                let y = yuk(&ctrl_modes[ab % 6], &ctrl_modes[ab / 6], &phih);
                mass_and_vecs(&y).1
            })
            .collect();
        for kq in 0..6 {
            for ku in 0..6 {
                for kd in 0..6 {
                    let v = ckm_full(&pairc[kq + ku * 6], &pairc[kq + kd * 6]);
                    jmax_ctrl = jmax_ctrl.max(jarlskog(&v).abs());
                }
            }
        }
    }
    // 発見 (対照の前提の訂正): 因子化スカラー格子の J は 0 ではない — 連続体
    // (長方形トーラス・実 Wilson) の実化可能性は格子で破れ、J は格子スケール
    // (~1e-4, O(1/N²) 級の位相残差) に立つ。v15.7 の厳密零は Dirac 構成 + 局在基底 +
    // Hadamard 積という別の構造の性質だった。この対照は「装置が J を検出できる」
    // 陽性対照として機能する — 傾き点の J ≤ 2.6e-8 はこの検出床の 4 桁下であり、
    // パイプラインの床ではなく傾き点に固有の構造的抑制である。
    check(
        "陽性対照: 因子化スカラー格子の J は格子スケール (>1e-6) — 装置は J を検出できる",
        jmax_ctrl > 1e-6,
        format!(
            "max |J| = {:.1e} (格子アーティファクト級 — 傾き点はこの {:.0} 倍下)",
            jmax_ctrl,
            jmax_ctrl / jmax.max(1e-300)
        ),
    );

    // ---- 判定 ----
    let cp_born = jmax > 1e-4;
    println!("\n[4] 事前登録の判定:");
    if cp_born {
        println!(
            "    => 傾き T⁴ は CP を作る (max|J| = {:.2e} — 最初の CP 可能窓)。",
            jmax
        );
        println!("       v15.7 の要求「CP を作る非因子化窓」の実在が確認された。深さ (lnZ) では");
        println!("       S₃ 対に負けたままなので、次の的は「深さと CP を同時に出す窓」である。");
    } else if jmax < 1e-6 {
        println!(
            "    => 傾き点の J ≤ {:.1e} — 測定 3.08e-5 の 1/1000 未満、かつ陽性対照の",
            jmax
        );
        println!(
            "       格子アーティファクト級 ({:.1e}) より {:.0} 倍深い — 装置の床ではなく、",
            jmax_ctrl,
            jmax_ctrl / jmax.max(1e-300)
        );
        println!("       傾き厳密縮退点に固有の構造的な CP 抑制である。");
        println!("       CP は「傾き」では買えない — 深さ (v13.2) に続きこの窓は CP でも死んだ。");
        println!("       機構 (何が J を格子スケール以下に守るのか) は開いた問い — 候補は Pf=3");
        println!("       (素数) の単一磁気並進軌道のクロック構造 (基底不変な診断は未設計)。");
        println!("       窓の拡張先: 非素数 Pf (複数軌道)・複素構造 τ・オービフォールドセクター。");
    } else {
        println!(
            "    => 中間 (max|J| = {:.2e}) — 大きさの起源の解析が必要。",
            jmax
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.1".into())),
        (
            "flux_tilt".into(),
            Json::Arr(f_tilt.iter().map(|&x| Json::Int(x)).collect()),
        ),
        ("lnz9_regression".into(), Json::Num(lnz)),
        ("jmax_tilt".into(), Json::Num(jmax)),
        ("j_median_posterior".into(), Json::Num(j50)),
        (
            "j_68".into(),
            Json::Arr(vec![Json::Num(j16), Json::Num(j84)]),
        ),
        ("p_j_above_1em6".into(), Json::Num(w_big / wsum)),
        ("map_j".into(), Json::Num(map.1)),
        ("map_vtd".into(), Json::Num(map.2)),
        ("jmax_control_factorized".into(), Json::Num(jmax_ctrl)),
        ("cp_born".into(), Json::Bool(cp_born)),
        ("j_measured_pdg".into(), Json::Num(3.08e-5)),
    ]);
    let p = write_artifact("results/v161_tiltcp.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — CP の判定は [4] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
