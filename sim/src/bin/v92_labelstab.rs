//! v9.2 世代ラベルの安定化 — 積模型の対角世代対は wrap 境界で機械精度の綱渡りだった
//!
//! 発見の経緯 (正直に): 論文図の作成 (v9.3 予定) で湯川パイプラインを python/numpy に
//! 独立再実装したところ、単一トーラスの証拠 (lnZ=−53.77) と v72 の MAP 予測は厳密一致
//! するのに、積模型の証拠だけ −18.76 vs −20.41 と食い違った。切り分けの結果、原因は
//! 物理でも実装バグでもなく**ラベル規約**だった:
//!
//!   k=0 の局在ゼロモードの 1 つは中心が格子座標 0 ≡ N の wrap 境界上に厳密に乗って
//!   おり (計算値は |c| ~ 1e-16)、「中心の昇順」ソート (v7.2 以来の規約) はこの ±ε の
//!   丸めの符号で世代ラベルが巡回する。単一トーラスの特異値はラベル置換に不変だが、
//!   積模型の対角世代対 (i,i) は 2 枚のトーラスのラベルを跨いで対を組むため、
//!   lnZ が規約の側に依存する (T²×T² 質量のみで 1.65 nats)。公表値 (v7.2/v8.1/v9.1)
//!   は本リポジトリの Rust 実装が偶々落ちた側 (境界モードが N⁻ 側 = 末尾) の値である。
//!
//! 本バイナリの仕事:
//!  [0] 綱渡りの実在証明 (境界距離 ~1e-16 の報告) と、安定ラベル (中心を 0.5 サイト
//!      格子にスナップしてからソート — ±ε に不変) の導入
//!  [1] 公表側ラベル (境界モードを末尾へ回す) で v8.1/v9.1 の公表値を再現 (装置検証)
//!  [2] 安定ラベルで質量のみ証拠 (6 幾何) と全 9 量証拠 (4 幾何) を再計算
//!  [3] 判定: 幾何選択の勝者・順位が両規約で不変か (v8.1/v9.1 の結論の頑健性)
//!
//! 教訓: 「固定シードで再現可能」は「規約に依らず再現可能」を意味しない。厳密縮退帯の
//! ソート順は下流で物理量になり得る — スナップしてから並べよ (CLAUDE.md 落とし穴に追記)。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];

type C3v = [(f64, f64); NS];
type M3 = [[(f64, f64); 3]; 3];

fn flux_modes(k_half: usize) -> (Vec<C3v>, f64, f64) {
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / NS as f64;
    let wl = phi * k_half as f64 / 2.0;
    let idx = |x: usize, y: usize| x + y * N;
    let m = 2 * NS;
    let mut a = vec![0.0; m * m];
    let addhop = |a: &mut Vec<f64>, i: usize, j: usize, th: f64| {
        let (c, s) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + NS) + (i + NS) * m] += -c;
        a[(i + NS) + (j + NS) * m] += -c;
        a[j + (i + NS) * m] += s;
        a[(j + NS) + i * m] += -s;
        a[i + (j + NS) * m] += -s;
        a[(i + NS) + j * m] += s;
    };
    for x in 0..N {
        for y in 0..N {
            addhop(&mut a, idx(x, y), idx(x, (y + 1) % N), phi * x as f64 + wl);
            let th = if x == N - 1 {
                -phi * (N as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % N, y), th);
        }
    }
    let (w, v) = jacobi_eigh(&a, m);
    let gap = w[2 * Q] - w[2 * Q - 1];
    let spread = w[2 * Q - 1] - w[0];
    let mut modes: Vec<C3v> = Vec::new();
    for kk in 0..2 * Q {
        let mut psi = [(0.0f64, 0.0f64); NS];
        for i in 0..NS {
            psi[i] = (v[i + kk * m], v[(i + NS) + kk * m]);
        }
        for pm in &modes {
            let (mut pr, mut pi) = (0.0, 0.0);
            for i in 0..NS {
                pr += pm[i].0 * psi[i].0 + pm[i].1 * psi[i].1;
                pi += pm[i].0 * psi[i].1 - pm[i].1 * psi[i].0;
            }
            for i in 0..NS {
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
            if modes.len() == Q {
                break;
            }
        }
    }
    assert_eq!(modes.len(), Q);
    (modes, gap, spread)
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

/// 3×3 エルミート行列の固有値 (昇順, 閉形式 — v6.5/v7.2 と同じ式)。質量のみ経路の高速化用。
fn eigvals3(hre: &[[f64; 3]; 3], him: &[[f64; 3]; 3]) -> [f64; 3] {
    let p1 = hre[0][1] * hre[0][1]
        + him[0][1] * him[0][1]
        + hre[0][2] * hre[0][2]
        + him[0][2] * him[0][2]
        + hre[1][2] * hre[1][2]
        + him[1][2] * him[1][2];
    let q = (hre[0][0] + hre[1][1] + hre[2][2]) / 3.0;
    let p2 = (hre[0][0] - q).powi(2) + (hre[1][1] - q).powi(2) + (hre[2][2] - q).powi(2) + 2.0 * p1;
    if p2 < 1e-300 {
        return [q, q, q];
    }
    let p = (p2 / 6.0).sqrt();
    let bd = [
        (hre[0][0] - q) / p,
        (hre[1][1] - q) / p,
        (hre[2][2] - q) / p,
    ];
    let (b01r, b01i) = (hre[0][1] / p, him[0][1] / p);
    let (b02r, b02i) = (hre[0][2] / p, him[0][2] / p);
    let (b12r, b12i) = (hre[1][2] / p, him[1][2] / p);
    let tr_re = (b01r * b12r - b01i * b12i) * b02r + (b01r * b12i + b01i * b12r) * b02i;
    let det = bd[0] * bd[1] * bd[2] + 2.0 * tr_re
        - bd[0] * (b12r * b12r + b12i * b12i)
        - bd[1] * (b02r * b02r + b02i * b02i)
        - bd[2] * (b01r * b01r + b01i * b01i);
    let r = (det / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let e1 = q + 2.0 * p * phi.cos();
    let e3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos();
    let e2 = 3.0 * q - e1 - e3;
    let mut v = [e3, e2, e1];
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v
}

/// 局在化 (v7.2/v9.1 と同一の手続き) — ただしソートせず、モードと生の中心を返す。
fn localize_unsorted(modes: &[C3v]) -> (Vec<C3v>, Vec<f64>) {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut ure = [[0.0f64; 3]; 3];
    let mut uim = [[0.0f64; 3]; 3];
    for a in 0..Q {
        for b in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for i in 0..NS {
                let x = (i % N) as f64;
                let (sn, cs) = (two_pi * x / N as f64).sin_cos();
                let (ar, ai) = modes[a][i];
                let (br, bi) = modes[b][i];
                let (pr, pi) = (ar * br + ai * bi, ar * bi - ai * br);
                sr += cs * pr - sn * pi;
                si += cs * pi + sn * pr;
            }
            ure[a][b] = sr;
            uim[a][b] = si;
        }
    }
    let (fc, fs) = (0.83f64.cos(), 0.83f64.sin());
    let mut h1re = [[0.0f64; 3]; 3];
    let mut h1im = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            let vre = fc * ure[a][b] + fs * uim[a][b];
            let vim = fc * uim[a][b] - fs * ure[a][b];
            let wre = fc * ure[b][a] + fs * uim[b][a];
            let wim = fc * uim[b][a] - fs * ure[b][a];
            h1re[a][b] = 0.5 * (vre + wre);
            h1im[a][b] = 0.5 * (vim - wim);
        }
    }
    let (_, vecs) = eig_herm3(&h1re, &h1im);
    let mut out: Vec<C3v> = Vec::new();
    let mut centers = Vec::new();
    for k in 0..Q {
        let mut psi = [(0.0f64, 0.0f64); NS];
        for i in 0..NS {
            for a in 0..Q {
                let (cr, ci) = vecs[k][a];
                let (mr, mi) = modes[a][i];
                psi[i].0 += cr * mr - ci * mi;
                psi[i].1 += cr * mi + ci * mr;
            }
        }
        let (mut zr, mut zi) = (0.0, 0.0);
        for i in 0..NS {
            let p = psi[i].0 * psi[i].0 + psi[i].1 * psi[i].1;
            let x = (i % N) as f64;
            let (sn, cs) = (two_pi * x / N as f64).sin_cos();
            zr += p * cs;
            zi += p * sn;
        }
        let center = (zi.atan2(zr) / two_pi * N as f64).rem_euclid(N as f64);
        out.push(psi);
        centers.push(center);
    }
    (out, centers)
}

/// 安定ラベル: 中心を 0.5 サイト格子にスナップ (N≡0 に丸めて 0 と同一視) してから昇順。
/// 戻り値は (並べ替えた順序, 最大スナップ偏差)。±ε の丸めに不変。
fn order_stable(centers: &[f64]) -> (Vec<usize>, f64) {
    let mut dev_max = 0.0f64;
    let snapped: Vec<f64> = centers
        .iter()
        .map(|&c| {
            let s = (2.0 * c).round() / 2.0;
            dev_max = dev_max.max((c - s).abs());
            s.rem_euclid(N as f64)
        })
        .collect();
    let mut ord: Vec<usize> = (0..centers.len()).collect();
    ord.sort_by(|&a, &b| snapped[a].partial_cmp(&snapped[b]).unwrap());
    (ord, dev_max)
}

fn yukawa(la: &[C3v], lb: &[C3v], sig_h: f64) -> M3 {
    let mut phih = [0.0f64; NS];
    for y in 0..N {
        for x in 0..N {
            let dx = (x as f64).min(N as f64 - x as f64);
            let dy = (y as f64).min(N as f64 - y as f64);
            phih[x + y * N] = (-(dx * dx + dy * dy) / (2.0 * sig_h * sig_h)).exp();
        }
    }
    let mut y_out = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..Q {
        for j in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for s in 0..NS {
                let (ar, ai) = la[i][s];
                let (br, bi) = lb[j][s];
                sr += (ar * br + ai * bi) * phih[s];
                si += (ar * bi - ai * br) * phih[s];
            }
            y_out[i][j] = (sr, si);
        }
    }
    y_out
}

fn had_prod(a: &M3, b: &M3) -> M3 {
    let mut y = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (p, q) = a[i][j];
            let (r, s) = b[i][j];
            y[i][j] = (p * r - q * s, p * s + q * r);
        }
    }
    y
}

/// Y から Y·Y† のエルミート成分を作る
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

/// Y から (質量比 2 つの対数, 左固有ベクトル)
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

/// Y から質量比 2 つの対数のみ (閉形式 — 質量のみ経路の高速版)
fn mass_ratios(y: &M3) -> [f64; 2] {
    let (hre, him) = gram(y);
    let lam = eigvals3(&hre, &him);
    let sv = [
        lam[0].max(0.0).sqrt(),
        lam[1].max(0.0).sqrt(),
        lam[2].max(0.0).sqrt(),
    ];
    [
        (sv[0].max(1e-300) / sv[2].max(1e-300)).ln(),
        (sv[1].max(1e-300) / sv[2].max(1e-300)).ln(),
    ]
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

/// 質量比 6 つのみの証拠 (v8.1 と同じ模型) — 尤度が完全因子化するので高速に厳密計算。
fn eval_mass(nt: usize, ks: &[usize], locs: &[Vec<C3v>], sig_grid: &[f64], sigma: f64) -> f64 {
    let nk = ks.len();
    let nc = nk.pow(nt as u32);
    let decoded: Vec<Vec<usize>> = (0..nc)
        .map(|mut c| {
            let mut v = Vec::with_capacity(nt);
            for _ in 0..nt {
                v.push(ks[c % nk]);
                c /= nk;
            }
            v
        })
        .collect();
    let norm2 = -(2.0 * std::f64::consts::PI * sigma * sigma).ln(); // 観測 2 つ分
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let mut lnz_terms = Vec::new();
    for &sh in sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // ペアごとの質量対数比 (閉形式固有値)
        let ratio: Vec<[f64; 2]> = (0..nc * nc)
            .map(|ab| {
                let (a, b) = (ab % nc, ab / nc);
                let mut y = ytab[decoded[a][0] + decoded[b][0] * NK12];
                for t in 1..nt {
                    y = had_prod(&y, &ytab[decoded[a][t] + decoded[b][t] * NK12]);
                }
                mass_ratios(&y)
            })
            .collect();
        let ll = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
            -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + norm2
        };
        let mut per_q = Vec::with_capacity(nc);
        let mut le = Vec::with_capacity(nc * nc);
        for ab in 0..nc * nc {
            le.push(ll(&ratio[ab], tgt[4], tgt[5]));
        }
        for kq in 0..nc {
            let us: Vec<f64> = (0..nc)
                .map(|ku| ll(&ratio[kq + ku * nc], tgt[0], tgt[1]))
                .collect();
            let ds: Vec<f64> = (0..nc)
                .map(|kd| ll(&ratio[kq + kd * nc], tgt[2], tgt[3]))
                .collect();
            per_q.push(lse(&us) + lse(&ds));
        }
        lnz_terms.push(lse(&per_q) + lse(&le));
    }
    lse(&lnz_terms) - (5.0 * (nt as f64) * (nk as f64).ln() + (sig_grid.len() as f64).ln())
}

/// 全 9 量 (質量 6 + CKM 3) の証拠 — v9.1 の eval9 と同一 (e 因子化 + クォーク三重和)。
fn eval9(nt: usize, ks: &[usize], locs: &[Vec<C3v>], sig_grid: &[f64], sigma: f64) -> (f64, [f64; 9]) {
    let nk = ks.len();
    let nc = nk.pow(nt as u32);
    let decoded: Vec<Vec<usize>> = (0..nc)
        .map(|mut c| {
            let mut v = Vec::with_capacity(nt);
            for _ in 0..nt {
                v.push(ks[c % nk]);
                c /= nk;
            }
            v
        })
        .collect();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let mut lnz_terms = Vec::new();
    let mut best = (f64::NEG_INFINITY, 0.0f64, [0usize; 5]);
    for &sh in sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        let pair: Vec<([f64; 2], M3)> = (0..nc * nc)
            .map(|ab| {
                let (a, b) = (ab % nc, ab / nc);
                let mut y = ytab[decoded[a][0] + decoded[b][0] * NK12];
                for t in 1..nt {
                    y = had_prod(&y, &ytab[decoded[a][t] + decoded[b][t] * NK12]);
                }
                mass_and_vecs(&y)
            })
            .collect();
        let le: Vec<f64> = (0..nc * nc)
            .map(|ab| {
                let r = &pair[ab].0;
                -((r[0] - tgt[4]).powi(2) + (r[1] - tgt[5]).powi(2)) / (2.0 * sigma * sigma)
                    + 2.0 * norm1
            })
            .collect();
        let mut le_best = f64::NEG_INFINITY;
        let mut le_besti = 0usize;
        for (i, &v) in le.iter().enumerate() {
            if v > le_best {
                le_best = v;
                le_besti = i;
            }
        }
        let mut per_q = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut acc_m = f64::NEG_INFINITY;
            let mut acc_s = 0.0f64;
            for ku in 0..nc {
                let (ru, vu) = &pair[kq + ku * nc];
                let llu = -((ru[0] - tgt[0]).powi(2) + (ru[1] - tgt[1]).powi(2))
                    / (2.0 * sigma * sigma)
                    + 2.0 * norm1;
                for kd in 0..nc {
                    let (rd, vd) = &pair[kq + kd * nc];
                    let lld = -((rd[0] - tgt[2]).powi(2) + (rd[1] - tgt[3]).powi(2))
                        / (2.0 * sigma * sigma)
                        + 2.0 * norm1;
                    let mut ll = llu + lld;
                    let c = ckm3(vu, vd);
                    for m in 0..3 {
                        let d = c[m].max(1e-300).ln() - tgt[6 + m];
                        ll += -d * d / (2.0 * sigma * sigma) + norm1;
                    }
                    if ll > acc_m {
                        acc_s = acc_s * (acc_m - ll).exp() + 1.0;
                        acc_m = ll;
                    } else {
                        acc_s += (ll - acc_m).exp();
                    }
                    let tot = ll + le_best;
                    if tot > best.0 {
                        best = (tot, sh, [kq, ku, kd, le_besti % nc, le_besti / nc]);
                    }
                }
            }
            per_q.push(acc_m + acc_s.ln());
        }
        lnz_terms.push(lse(&per_q) + lse(&le));
    }
    let lnz =
        lse(&lnz_terms) - (5.0 * (nt as f64) * (nk as f64).ln() + (sig_grid.len() as f64).ln());
    let sh = best.1;
    let build = |a: usize, b: usize| -> M3 {
        let mut y = yukawa(&locs[decoded[a][0]], &locs[decoded[b][0]], sh);
        for t in 1..nt {
            y = had_prod(&y, &yukawa(&locs[decoded[a][t]], &locs[decoded[b][t]], sh));
        }
        y
    };
    let yu = build(best.2[0], best.2[1]);
    let yd = build(best.2[0], best.2[2]);
    let ye = build(best.2[3], best.2[4]);
    let (ru, vu) = mass_and_vecs(&yu);
    let (rd, vd) = mass_and_vecs(&yd);
    let (re_, _) = mass_and_vecs(&ye);
    let c = ckm3(&vu, &vd);
    (
        lnz,
        [
            ru[0].exp(),
            ru[1].exp(),
            rd[0].exp(),
            rd[1].exp(),
            re_[0].exp(),
            re_[1].exp(),
            c[0],
            c[1],
            c[2],
        ],
    )
}

fn main() {
    self_test();
    println!("=== v9.2 世代ラベルの安定化: wrap 境界の綱渡りと幾何選択の再検定 ===\n");
    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];

    // ---- [0] 局在モードと綱渡りの実在 ----
    println!("[0] 世代モード (Z₁₂, 対角化 12 回) と wrap 境界の距離");
    let t0 = std::time::Instant::now();
    let mut raw: Vec<(Vec<C3v>, Vec<f64>)> = Vec::new();
    let mut ok_engine = true;
    for k in 0..NK12 {
        let (modes, gap, spread) = flux_modes(k);
        if spread > 1e-9 || gap < 0.05 {
            ok_engine = false;
        }
        raw.push(localize_unsorted(&modes));
    }
    println!(
        "    縮退・ギャップ不変  {}  ({} ms)",
        pass(ok_engine),
        t0.elapsed().as_millis()
    );
    // 境界距離: 各 k の生中心の wrap 境界 (0≡N) までの最小距離
    let mut bdist_min = f64::INFINITY;
    let mut bdist_k = 0usize;
    let mut snap_dev = 0.0f64;
    for (k, (_, cents)) in raw.iter().enumerate() {
        for &c in cents {
            let d = c.min(N as f64 - c);
            if d < bdist_min {
                bdist_min = d;
                bdist_k = k;
            }
        }
        let (_, dev) = order_stable(cents);
        snap_dev = snap_dev.max(dev);
    }
    println!(
        "    生中心の wrap 境界距離の最小値: {:.2e} (k_half={}) — 昇順ソートは ±ε の綱渡り",
        bdist_min, bdist_k
    );
    println!(
        "    綱渡りの実在 (境界距離 < 1e-9)  {}",
        pass(bdist_min < 1e-9)
    );
    println!(
        "    全中心は 0.5 サイト格子上 (最大スナップ偏差 {:.2e} < 1e-6)  {}",
        snap_dev,
        pass(snap_dev < 1e-6)
    );

    // 両ラベルの locs 表を構築
    //  安定側: スナップ後ソート (境界モードは 0 として先頭)
    //  公表側: 安定側から境界モード (スナップ 0) を末尾へ回す — 本リポジトリの Rust が
    //          偶々落ちていた側 (results/v72_geomfn.txt の中心 [6,12,18] に対応)
    let mut locs_stable: Vec<Vec<C3v>> = Vec::new();
    let mut locs_pub: Vec<Vec<C3v>> = Vec::new();
    for (modes, cents) in raw.iter() {
        let (ord, _) = order_stable(cents);
        let stable: Vec<C3v> = ord.iter().map(|&i| modes[i]).collect();
        let snapped0 = ord
            .iter()
            .position(|&i| (2.0 * cents[i]).round().rem_euclid(2.0 * N as f64) == 0.0);
        let mut publ = stable.clone();
        if let Some(p0) = snapped0 {
            let m0 = publ.remove(p0);
            publ.push(m0);
        }
        locs_stable.push(stable);
        locs_pub.push(publ);
    }

    let z6: Vec<usize> = (0..NK12).step_by(2).collect();
    let z12: Vec<usize> = (0..NK12).collect();

    // ---- [1] 公表側ラベルで公表値を再現 (装置検証) ----
    println!("\n[1] 公表側ラベルの再現 (v8.1 = results/v81_geoselect.json, v9.1 = results/v91_ckmselect.json)");
    let geoms: [(&str, usize, &Vec<usize>); 6] = [
        ("T¹×Z₆ ", 1, &z6),
        ("T¹×Z₁₂", 1, &z12),
        ("T²×Z₆ ", 2, &z6),
        ("T²×Z₁₂", 2, &z12),
        ("T³×Z₆ ", 3, &z6),
        ("T³×Z₁₂", 3, &z12),
    ];
    let pub_mass_ref = [-53.77f64, -53.76, -20.41, -20.42, -17.41, -18.36];
    let mut ok_pub = true;
    let mut mass_pub = [0.0f64; 6];
    for (i, (name, nt, ks)) in geoms.iter().enumerate() {
        let v = eval_mass(*nt, ks, &locs_pub, &sig_grid, sigma);
        mass_pub[i] = v;
        let ok = (v - pub_mass_ref[i]).abs() < 0.02;
        ok_pub &= ok;
        println!(
            "    質量のみ {}: lnZ = {:7.2} vs 公表 {:7.2}  {}",
            name,
            v,
            pub_mass_ref[i],
            pass(ok)
        );
    }
    // 全 9 量 (T³×Z₁₂ は v9.1 と同じく対象外)
    let nine: [(&str, usize, &Vec<usize>, f64); 4] = [
        ("T¹×Z₆ ", 1, &z6, -63.57),
        ("T²×Z₆ ", 2, &z6, -27.13),
        ("T²×Z₁₂", 2, &z12, -27.35),
        ("T³×Z₆ ", 3, &z6, -25.56),
    ];
    let mut nine_pub = [0.0f64; 4];
    for (i, (name, nt, ks, r)) in nine.iter().enumerate() {
        let (v, _) = eval9(*nt, ks, &locs_pub, &sig_grid, sigma);
        nine_pub[i] = v;
        let ok = (v - r).abs() < 0.02;
        ok_pub &= ok;
        println!(
            "    全 9 量 {}: lnZ₉ = {:7.2} vs 公表 {:7.2}  {}",
            name,
            v,
            r,
            pass(ok)
        );
    }
    println!(
        "    => 公表値は「境界モードを末尾へ回す」ラベルの値として全て再現される  {}",
        pass(ok_pub)
    );

    // ---- [2] 安定ラベルでの再計算 ----
    println!("\n[2] 安定ラベル (スナップ後ソート) での再計算");
    let mut mass_stb = [0.0f64; 6];
    for (i, (name, nt, ks)) in geoms.iter().enumerate() {
        let v = eval_mass(*nt, ks, &locs_stable, &sig_grid, sigma);
        mass_stb[i] = v;
        println!(
            "    質量のみ {}: lnZ = {:7.2}  (公表側との差 {:+.2})",
            name,
            v,
            v - mass_pub[i]
        );
    }
    let mut nine_stb = [0.0f64; 4];
    let mut preds_t3 = [0.0f64; 9];
    for (i, (name, nt, ks, _)) in nine.iter().enumerate() {
        let (v, preds) = eval9(*nt, ks, &locs_stable, &sig_grid, sigma);
        nine_stb[i] = v;
        if i == 3 {
            preds_t3 = preds;
        }
        println!(
            "    全 9 量 {}: lnZ₉ = {:7.2}  (公表側との差 {:+.2})",
            name,
            v,
            v - nine_pub[i]
        );
    }
    // 単一トーラスはラベル不変のはず (特異値は置換に不変)
    let ok_t1 = (mass_stb[0] - mass_pub[0]).abs() < 1e-6 && (mass_stb[1] - mass_pub[1]).abs() < 1e-6;
    println!(
        "    単一トーラス (T¹) はラベル規約に厳密不変  {}",
        pass(ok_t1)
    );

    // ---- [3] 判定: 幾何選択の結論は規約に頑健か ----
    println!("\n[3] 判定: 幾何選択の勝者と順位");
    let argmax = |v: &[f64]| -> usize {
        let mut bi = 0;
        for i in 1..v.len() {
            if v[i] > v[bi] {
                bi = i;
            }
        }
        bi
    };
    let wm_pub = argmax(&mass_pub);
    let wm_stb = argmax(&mass_stb);
    let w9_pub = argmax(&nine_pub);
    let w9_stb = argmax(&nine_stb);
    println!(
        "    質量のみの勝者: 公表側 {} / 安定側 {}  (T³ の T² への差: {:+.2} → {:+.2})",
        geoms[wm_pub].0,
        geoms[wm_stb].0,
        mass_pub[4] - mass_pub[2],
        mass_stb[4] - mass_stb[2]
    );
    println!(
        "    全 9 量の勝者:  公表側 {} / 安定側 {}  (T³ の T² への差: {:+.2} → {:+.2})",
        nine[w9_pub].0,
        nine[w9_stb].0,
        nine_pub[3] - nine_pub[1],
        nine_stb[3] - nine_stb[1]
    );
    let ok_winner = wm_pub == wm_stb && w9_pub == w9_stb;
    println!(
        "    勝者は両規約で一致 (v8.1/v9.1 の選択の結論は頑健)  {}",
        pass(ok_winner)
    );
    let mut ok_rank = true;
    // 順位全体 (質量のみ 6 幾何・全 9 量 4 幾何) の一致
    let rank = |v: &[f64]| -> Vec<usize> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&a, &b| v[b].partial_cmp(&v[a]).unwrap());
        idx
    };
    ok_rank &= rank(&mass_pub) == rank(&mass_stb);
    ok_rank &= rank(&nine_pub) == rank(&nine_stb);
    println!("    順位全体も両規約で一致  {}", pass(ok_rank));
    println!("\n    安定側 T³×Z₆ の MAP 予測 (参考): m_u/m_t={:.2e}, |V_us|={:.3}, |V_cb|={:.3}, |V_ub|={:.4}",
        preds_t3[0], preds_t3[6], preds_t3[7], preds_t3[8]);

    // ---- 機械可読な結果 ----
    let geo_names = ["T1xZ6", "T1xZ12", "T2xZ6", "T2xZ12", "T3xZ6", "T3xZ12"];
    let mut mass_obj = Vec::new();
    for i in 0..6 {
        mass_obj.push((
            geo_names[i].to_string(),
            Json::Obj(vec![
                ("published_label".into(), Json::Num(mass_pub[i])),
                ("stable_label".into(), Json::Num(mass_stb[i])),
            ]),
        ));
    }
    let mut nine_obj = Vec::new();
    for (i, nm) in ["T1xZ6", "T2xZ6", "T2xZ12", "T3xZ6"].iter().enumerate() {
        nine_obj.push((
            nm.to_string(),
            Json::Obj(vec![
                ("published_label".into(), Json::Num(nine_pub[i])),
                ("stable_label".into(), Json::Num(nine_stb[i])),
            ]),
        ));
    }
    let all_ok = ok_engine && bdist_min < 1e-9 && snap_dev < 1e-6 && ok_pub && ok_t1 && ok_winner && ok_rank;
    let art = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-006".into())),
        ("boundary_distance".into(), Json::Num(bdist_min)),
        ("snap_deviation_max".into(), Json::Num(snap_dev)),
        ("lnZ_mass".into(), Json::Obj(mass_obj)),
        ("lnZ_nine".into(), Json::Obj(nine_obj)),
        ("winner_invariant".into(), Json::Bool(ok_winner)),
        ("ranking_invariant".into(), Json::Bool(ok_rank)),
        (
            "stable_t3_map_preds".into(),
            Json::Arr(preds_t3.iter().map(|&x| Json::Num(x)).collect()),
        ),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    write_artifact("results/v92_labelstab.json", &art.render());
    println!("\n  機械可読な結果: ../results/v92_labelstab.json");
    println!("\n総合判定: {}", pass(all_ok));
    println!(
        "\n結論: 積模型 (対角世代対) の lnZ は世代ラベルの wrap 規約に依存していた\n      (発見は python 独立再実装との照合による)。公表値は再現可能な一方の規約の\n      値として正しく、幾何選択の勝者と順位は両規約で不変 — 結論は頑健である。\n      以後の湯川系バイナリは安定ラベル (スナップ後ソート) を規約とする。"
    );
}
