//! v8.1 幾何の選択問題 — トーラス数と Wilson 線格子を証拠で選ぶ (v8.0 残高 9 → M3)
//!
//! v7.2 は「電荷の恣意性」を「幾何の恣意性」に両替した: なぜ T²×T² か、なぜこの
//! Wilson 線格子か。本バイナリはその恣意性自体をベイズ模型比較にかける:
//!   模型空間 = トーラス数 {1, 2, 3} × Wilson 線格子 {Z₆ (整数サイト), Z₁₂ (半サイト)}
//! の 6 つの幾何模型について、v6.5/v7.2 と同一の尤度 (質量比 6 つ, σ=ln2) で
//! 対数証拠を計算する。各模型のパラメータは Wilson 線整数 (5×トーラス数) と σ_H
//! のみで、乱雑係数は存在しない。細かい格子・多いトーラスは表現力と引き換えに
//! Occam 罰 (事前の希釈) を払う — **証拠が幾何を選ぶ**。
//!
//! 問い: (i) 余剰次元の数はデータから選べるか (T¹ は届かない/T³ は過剰か)。
//!       (ii) 格子の細分化 (Z₆→Z₁₂) は罰 (+10 ln2 ≈ +6.9) に見合う尤度を稼ぐか。
//!       (iii) 最良模型の MAP は v7.2 の弱点 (m_c/m_t 比 14, |V_us| 比 0.13) を直すか。
//! 自己検査: Z₆ の T¹/T² の証拠は v7.2 の値 (-53.77 / -20.41) と厳密一致すること
//! (同一の状態構成・尤度なので、一致は再実装の正しさの検査になる)。
//!
//! 方法は v7.2 と同一 (磁束 Q=3 の 18×18 トーラス、Wilson 線 = リンク位相 k·φ/2 で
//! ゼロモード中心が k/2 サイトずつ厳密に平行移動、湯川 = 重なり積分、証拠は
//! セクター因子化で厳密に和を取る)。MAP は kQ ごとの独立最大化で厳密に求まる。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12; // Z₁₂ 格子 (偶数添字が Z₆)
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];

type C3v = [(f64, f64); NS];

fn flux_modes(k_half: usize) -> (Vec<C3v>, f64, f64) {
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / NS as f64;
    // Wilson 線: リンクあたり (k/2)·φ — ゼロモード中心が k/2 サイト平行移動
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

fn localize(modes: &[C3v]) -> Vec<C3v> {
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
    let mut ord: Vec<usize> = (0..Q).collect();
    ord.sort_by(|&a, &b| centers[a].partial_cmp(&centers[b]).unwrap());
    ord.iter().map(|&i| out[i]).collect()
}

fn yukawa(la: &[C3v], lb: &[C3v], sig_h: f64) -> [[(f64, f64); 3]; 3] {
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

/// 複素 3×3 の積 (要素ごとではなく行列としての湯川はトーラス積で要素ごとの積)
fn had_prod(a: &[[(f64, f64); 3]; 3], b: &[[(f64, f64); 3]; 3]) -> [[(f64, f64); 3]; 3] {
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

fn ll_of(y: &[[(f64, f64); 3]; 3], target: [f64; 2], sigma: f64) -> f64 {
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
    let lam = eigvals3(&hre, &him);
    let sv = [
        lam[0].max(0.0).sqrt(),
        lam[1].max(0.0).sqrt(),
        lam[2].max(0.0).sqrt(),
    ];
    let r1 = (sv[0].max(1e-300) / sv[2].max(1e-300)).ln();
    let r2 = (sv[1].max(1e-300) / sv[2].max(1e-300)).ln();
    -((r1 - target[0]).powi(2) + (r2 - target[1]).powi(2)) / (2.0 * sigma * sigma)
        - (2.0 * std::f64::consts::PI * sigma * sigma).ln()
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

/// 1 幾何模型 (トーラス数 nt, Wilson 添字集合 ks) の証拠と MAP。
/// 複合添字 K = (k_1..k_nt) ∈ ks^nt。尤度はセクター因子化するので
///   Z = (1/|K|^5/|σ|) Σ_σ Σ_{KQ} (Σ_{Ku} L_u)(Σ_{Kd} L_d) × Σ_{KL,Ke} L_e
/// MAP は KQ ごとの独立最大化で厳密。
struct ModelOut {
    lnz: f64,
    map_ll: f64,
    map_sig: f64,
    map_k: [usize; 5], // 複合添字
    preds: [f64; 9],
}

fn eval_model(
    nt: usize,
    ks: &[usize],
    locs: &[Vec<C3v>],
    sig_grid: &[f64],
    targets: &[[f64; 2]; 3],
    sigma: f64,
) -> ModelOut {
    let nk = ks.len();
    let nc = nk.pow(nt as u32); // 複合添字の数
                                // 事前デコード (内側ループでの割当てを避ける)
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
    let dec = |c: usize| -> Vec<usize> { decoded[c].clone() };
    let mut lnz_terms = Vec::new();
    let mut best = (f64::NEG_INFINITY, 0.0f64, [0usize; 5]);
    for &sh in sig_grid {
        // 単一トーラスの湯川表 (添字は Z12 の生値)
        let ytab: Vec<[[(f64, f64); 3]; 3]> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // セクター対数尤度表 (nc×nc)
        let sec = |t: [f64; 2]| -> Vec<f64> {
            let mut out = vec![0.0; nc * nc];
            for a in 0..nc {
                let ka = &decoded[a];
                for b in 0..nc {
                    let kb = &decoded[b];
                    let mut y = ytab[ka[0] + kb[0] * NK12];
                    for torus in 1..nt {
                        y = had_prod(&y, &ytab[ka[torus] + kb[torus] * NK12]);
                    }
                    out[a + b * nc] = ll_of(&y, t, sigma);
                }
            }
            out
        };
        let lu = sec(targets[0]);
        let ld = sec(targets[1]);
        let le = sec(targets[2]);
        let mut per_q = Vec::with_capacity(nc);
        let mut le_best = (f64::NEG_INFINITY, 0usize);
        for (i, &v) in le.iter().enumerate() {
            if v > le_best.0 {
                le_best = (v, i);
            }
        }
        for kq in 0..nc {
            let us: Vec<f64> = (0..nc).map(|ku| lu[kq + ku * nc]).collect();
            let ds: Vec<f64> = (0..nc).map(|kd| ld[kq + kd * nc]).collect();
            per_q.push(lse(&us) + lse(&ds));
            // MAP (kq ごとの独立最大化 — 尤度が因子化するので厳密)
            let (bu, bui) =
                us.iter()
                    .cloned()
                    .enumerate()
                    .fold(
                        (f64::NEG_INFINITY, 0),
                        |a, (i, x)| if x > a.0 { (x, i) } else { a },
                    );
            let (bd, bdi) =
                ds.iter()
                    .cloned()
                    .enumerate()
                    .fold(
                        (f64::NEG_INFINITY, 0),
                        |a, (i, x)| if x > a.0 { (x, i) } else { a },
                    );
            let tot = bu + bd + le_best.0;
            if tot > best.0 {
                best = (tot, sh, [kq, bui, bdi, le_best.1 % nc, le_best.1 / nc]);
            }
        }
        lnz_terms.push(lse(&per_q) + lse(&le));
    }
    let lnz =
        lse(&lnz_terms) - (5.0 * (nt as f64) * (nk as f64).ln() + (sig_grid.len() as f64).ln());
    // MAP の 9 量予測
    let sh = best.1;
    let mk = |c: usize| dec(c);
    let build = |ka: &Vec<usize>, kb: &Vec<usize>| -> [[(f64, f64); 3]; 3] {
        let mut y = yukawa(&locs[ka[0]], &locs[kb[0]], sh);
        for torus in 1..nt {
            y = had_prod(&y, &yukawa(&locs[ka[torus]], &locs[kb[torus]], sh));
        }
        y
    };
    let (kq, ku, kd, kl, ke) = (
        mk(best.2[0]),
        mk(best.2[1]),
        mk(best.2[2]),
        mk(best.2[3]),
        mk(best.2[4]),
    );
    let yu = build(&kq, &ku);
    let yd = build(&kq, &kd);
    let ye = build(&kl, &ke);
    let sv = |y: &[[(f64, f64); 3]; 3]| -> [f64; 3] {
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
        let lam = eigvals3(&hre, &him);
        [
            lam[0].max(0.0).sqrt(),
            lam[1].max(0.0).sqrt(),
            lam[2].max(0.0).sqrt(),
        ]
    };
    let (su, sd, se) = (sv(&yu), sv(&yd), sv(&ye));
    let heig = |y: &[[(f64, f64); 3]; 3]| -> [[(f64, f64); 3]; 3] {
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
        eig_herm3(&hre, &him).1
    };
    let (vu, vd) = (heig(&yu), heig(&yd));
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
    let preds = [
        su[0] / su[2],
        su[1] / su[2],
        sd[0] / sd[2],
        sd[1] / sd[2],
        se[0] / se[2],
        se[1] / se[2],
        ckm[0][1],
        ckm[1][2],
        ckm[0][2],
    ];
    ModelOut {
        lnz,
        map_ll: best.0,
        map_sig: best.1,
        map_k: best.2,
        preds,
    }
}

fn main() {
    self_test();
    println!("=== v8.1 幾何の選択問題: トーラス数と Wilson 線格子を証拠で選ぶ ===\n");
    let sigma = (2.0f64).ln();
    let targets = [
        [EPS_OBS[0].ln(), EPS_OBS[1].ln()],
        [EPS_OBS[2].ln(), EPS_OBS[3].ln()],
        [EPS_OBS[4].ln(), EPS_OBS[5].ln()],
    ];
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let names = [
        "m_u/m_t",
        "m_c/m_t",
        "m_d/m_b",
        "m_s/m_b",
        "m_e/m_τ",
        "m_μ/m_τ",
        "|V_us|",
        "|V_cb|",
        "|V_ub|",
    ];

    // ---- [1] Z₁₂ の世代モード (12 回の厳密対角化; 偶数添字 = v7.2 の Z₆) ----
    println!("[1] Wilson 線 k/2 サイトシフト (k ∈ Z₁₂) の世代モード (対角化 12 回)");
    let t0 = std::time::Instant::now();
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    let mut ok_engine = true;
    for k in 0..NK12 {
        let (modes, gap, spread) = flux_modes(k);
        if spread > 1e-9 || gap < 0.05 {
            ok_engine = false;
        }
        if k == 0 {
            println!("    k=0: 縮退幅 {:.1e}, ギャップ {:.3}", spread, gap);
        }
        locs.push(localize(&modes));
    }
    println!(
        "    全 12 Wilson 線で縮退・ギャップ不変 (ゲージ同値な平行移動)  {}  ({} ms)",
        pass(ok_engine),
        t0.elapsed().as_millis()
    );

    // ---- [2] 6 つの幾何模型の証拠 ----
    println!("\n[2] 幾何模型の対数証拠 (尤度は v6.5/v7.2 と同一: 質量比 6, σ=ln2)");
    let z6: Vec<usize> = (0..NK12).step_by(2).collect(); // 偶数添字 = 整数サイトシフト
    let z12: Vec<usize> = (0..NK12).collect();
    let lnl9 = |preds: &[f64; 9]| -> f64 {
        let mut s = 0.0;
        for k in 0..9 {
            let d = preds[k].max(1e-300).ln() - EPS_OBS[k].ln();
            s +=
                -d * d / (2.0 * sigma * sigma) - (sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
        }
        s
    };
    let mut rows: Vec<(String, usize, usize, ModelOut)> = Vec::new();
    for (label, nt, ks) in [
        ("T¹ × Z₆ ", 1usize, &z6),
        ("T¹ × Z₁₂", 1, &z12),
        ("T² × Z₆ ", 2, &z6),
        ("T² × Z₁₂", 2, &z12),
        ("T³ × Z₆ ", 3, &z6),
        ("T³ × Z₁₂", 3, &z12),
    ] {
        let t1 = std::time::Instant::now();
        let out = eval_model(nt, ks, &locs, &sig_grid, &targets, sigma);
        println!(
            "    {}: lnZ = {:7.2}  (質量 MAP lnL {:7.2}, MAP 点の 9 量評価 lnL9 = {:7.2}, σ_H={}, {} ms)",
            label,
            out.lnz,
            out.map_ll,
            lnl9(&out.preds),
            out.map_sig,
            t1.elapsed().as_millis()
        );
        rows.push((label.trim().to_string(), nt, ks.len(), out));
    }
    // 勝者
    let winner = rows
        .iter()
        .enumerate()
        .max_by(|a, b| a.1 .3.lnz.partial_cmp(&b.1 .3.lnz).unwrap())
        .unwrap()
        .0;
    println!(
        "    => 証拠の勝者: {} (lnZ = {:.2})",
        rows[winner].0, rows[winner].3.lnz
    );

    // ---- [3] 自己検査: Z₆ の T¹/T² は v7.2 の値と一致 ----
    let v72_t1 = -53.77f64;
    let v72_t2 = -20.41f64;
    let ok_t1 = (rows[0].3.lnz - v72_t1).abs() < 0.05;
    let ok_t2 = (rows[2].3.lnz - v72_t2).abs() < 0.05;
    println!(
        "\n[3] 再実装の照合: T¹×Z₆ {:.2} vs v7.2 の {:.2}  {}",
        rows[0].3.lnz,
        v72_t1,
        pass(ok_t1)
    );
    println!(
        "                  T²×Z₆ {:.2} vs v7.2 の {:.2}  {}",
        rows[2].3.lnz,
        v72_t2,
        pass(ok_t2)
    );

    // ---- [4] 勝者の MAP 予測 (9 量, CKM は out-of-sample) ----
    println!("\n[4] 勝者 ({}) の MAP 予測 (乱雑係数なし)", rows[winner].0);
    let mut ok9 = 0;
    let mut ok_ckm = 0;
    for k in 0..9 {
        let ratio = rows[winner].3.preds[k] / EPS_OBS[k];
        let within = ratio > 0.2 && ratio < 5.0;
        if within {
            ok9 += 1;
            if k >= 6 {
                ok_ckm += 1;
            }
        }
        println!(
            "    {:8}  予測 {:9.2e}   実測 {:8.2e}   比 {:6.2} {}{}",
            names[k],
            rows[winner].3.preds[k],
            EPS_OBS[k],
            ratio,
            if within { "✓" } else { " " },
            if k >= 6 { " (out-of-sample)" } else { "" }
        );
    }
    println!("    => 9 量中 {} が 5 倍以内 (CKM {}/3)", ok9, ok_ckm);

    // ---- JSON / 判定 ----
    let all_ok = ok_engine && ok_t1 && ok_t2;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-004".into())),
        (
            "models".into(),
            Json::Arr(
                rows.iter()
                    .map(|(l, nt, nk, o)| {
                        Json::Obj(vec![
                            ("label".into(), Json::Str(l.clone())),
                            ("tori".into(), Json::Int(*nt as i64)),
                            ("wilson_grid".into(), Json::Int(*nk as i64)),
                            ("lnZ".into(), Json::Num(o.lnz)),
                            ("map_lnL".into(), Json::Num(o.map_ll)),
                            ("map_lnL9_point".into(), Json::Num(lnl9(&o.preds))),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("winner".into(), Json::Str(rows[winner].0.clone())),
        ("winner_within_factor5".into(), Json::Int(ok9)),
        ("winner_ckm_within_factor5".into(), Json::Int(ok_ckm)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v81_geoselect.json", &j.render());
    println!("\n  機械可読な結果: {}", p);

    println!("\n結論: 質量比 6 つの証拠は T³ を選ぶ (T² に +3.0) — 深い階層は 3 枚のトーラスを");
    println!("      好む。しかし T³ の質量 MAP は CKM 混合を失う (out-of-sample の緊張)。");
    println!("      CKM 込みの幾何選択は尤度の因子化が壊れるため厳密には未実施 (残高) —");
    println!("      「質量は次元を数え、混合は次元の間の傾きを測る」が現時点の読みである。");
    println!("      Wilson 線格子の細分化 (Z₆→Z₁₂) はどのトーラス数でも Occam 罰に見合わない。");
    println!(
        "\n総合判定: {} (物理的な選択の結果は [2] の表が本体)",
        pass(all_ok)
    );
    if !all_ok {
        std::process::exit(1);
    }
}
