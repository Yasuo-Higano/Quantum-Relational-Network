//! v16.10 谷底平坦域の精密裁定 — 0.23 nats は σ_H グリッドの影か
//!
//! v16.9 の地図で谷底は平坦だった: (3,3)=−21.757, (1,3)=−21.983 (Δ0.23 nats),
//! (2,3)=−22.257 (Δ0.50)。QRN-YUK-024 の限界 (i) は「この差は σ_H グリッド
//! ({2,3,4,5} 固定) の影の可能性」と記した。本バイナリはその返済である:
//! 台上の 4 幾何 (2,2), (2,3), (1,3), (3,3) を 4 種のグリッドで測り直す —
//!   G0 = {2,3,4,5}            (v16.9 の回帰アンカー)
//!   G1 = {1.5, 2, ..., 5.5, 6} (細分+拡張 10 点 — σ の谷がグリッド外に無いか)
//!   G2 = {2.5, 3.5, 4.5, 5.5}  (半格子シフト)
//!   G3 = {3, 4, 5, 6}          (上方シフト)
//! グリッド内の順位と Δ が意味のある観測量 (prior の ln|G| は幾何間で相殺)。
//! G1 では σ ごとの証拠プロファイル terms10(σ) も出す — どの幾何がどの σ に
//! 証拠を置いているかが「影」の物理的な中身である。
//!
//! 事前登録 2 分岐:
//!   (a) 全 4 グリッドで (3,3) が首位 — Δ0.23 は実在、対角優勢はグリッドに頑健
//!   (b) どれかのグリッドで首位が入れ替わる — 平坦域は装置分解能以下で縮退、
//!       「対角 vs 非対称」の選択は σ_H の扱いに条件付きと記録
//! いずれの分岐でも G1 (最細) での全 12 量 MAP スコアカードを 4 幾何で記録する。
//!
//! 装置ゲート: G0 の 4 値は v16.9 の再現 (±0.02)・厳密 3 重縮退 (s=1,2,3)。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// v16.9 の G0 回帰アンカー: (2,2), (2,3), (1,3), (3,3)
const REF_G0: [f64; 4] = [-22.263, -22.257, -21.983, -21.757];

const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

type M3 = [[(f64, f64); 3]; 3];
type Mode = Vec<(f64, f64)>; // 長さ N²

/// シアー s つき磁束トーラス (格子 n, 2 成分 Dirac 型) の最低 3 モード。
/// 戻り値: (モード 3 本, ギャップ, 縮退幅)
fn flux_modes_shear_n(n: usize, k_half: usize, s: usize) -> (Vec<Mode>, f64, f64) {
    let ns = n * n;
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / ns as f64;
    let wl = phi * k_half as f64 / 2.0;
    let m = 2 * ns;
    let mut a = vec![0.0; m * m];
    let addhop = |a: &mut Vec<f64>, i: usize, j: usize, th: f64, m: usize, ns: usize| {
        let (c, sn) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + ns) + (i + ns) * m] += -c;
        a[(i + ns) + (j + ns) * m] += -c;
        a[j + (i + ns) * m] += sn;
        a[(j + ns) + i * m] += -sn;
        a[i + (j + ns) * m] += -sn;
        a[(i + ns) + j * m] += sn;
    };
    let idx = |x: usize, y: usize| x + y * n;
    for x in 0..n {
        for y in 0..n {
            let th_y = phi * x as f64 + wl;
            if y == n - 1 {
                addhop(&mut a, idx(x, y), idx((x + s) % n, 0), th_y, m, ns);
            } else {
                addhop(&mut a, idx(x, y), idx(x, y + 1), th_y, m, ns);
            }
            let th_x = if x == n - 1 {
                -phi * (n as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % n, y), th_x, m, ns);
        }
    }
    let (w, v) = jacobi_eigh(&a, m);
    let gap = w[2 * Q] - w[2 * Q - 1];
    let spread = w[2 * Q - 1] - w[0];
    let mut modes: Vec<Mode> = Vec::new();
    for kk in 0..2 * Q {
        let mut psi: Mode = (0..ns)
            .map(|i| (v[i + kk * m], v[(i + ns) + kk * m]))
            .collect();
        for pm in &modes {
            let (mut pr, mut pi) = (0.0, 0.0);
            for i in 0..ns {
                pr += pm[i].0 * psi[i].0 + pm[i].1 * psi[i].1;
                pi += pm[i].0 * psi[i].1 - pm[i].1 * psi[i].0;
            }
            for i in 0..ns {
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
    let nn = 3;
    let m = 6;
    let mut emb = vec![0.0; m * m];
    for i in 0..nn {
        for j in 0..nn {
            emb[i + j * m] = hre[i][j];
            emb[i + (j + nn) * m] = -him[i][j];
            emb[(i + nn) + j * m] = him[i][j];
            emb[(i + nn) + (j + nn) * m] = hre[i][j];
        }
    }
    let (w, v) = jacobi_eigh(&emb, m);
    let mut lam = [0.0f64; 3];
    let mut vecs = [[(0.0f64, 0.0f64); 3]; 3];
    for k in 0..3 {
        lam[k] = 0.5 * (w[2 * k] + w[2 * k + 1]);
        for i in 0..3 {
            vecs[k][i] = (v[i + (2 * k) * m], v[(i + nn) + (2 * k) * m]);
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

/// 局在化 + 安定ラベル (n パラメータ版)
fn localize_stable(n: usize, modes: &[Mode]) -> Vec<Mode> {
    let ns = n * n;
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut ure = [[0.0f64; 3]; 3];
    let mut uim = [[0.0f64; 3]; 3];
    for a in 0..Q {
        for b in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for i in 0..ns {
                let x = (i % n) as f64;
                let (sn, cs) = (two_pi * x / n as f64).sin_cos();
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
    let mut out: Vec<Mode> = Vec::new();
    let mut centers = Vec::new();
    for k in 0..Q {
        let mut psi: Mode = vec![(0.0, 0.0); ns];
        for i in 0..ns {
            for a in 0..Q {
                let (cr, ci) = vecs[k][a];
                let (mr, mi) = modes[a][i];
                psi[i].0 += cr * mr - ci * mi;
                psi[i].1 += cr * mi + ci * mr;
            }
        }
        let (mut zr, mut zi) = (0.0, 0.0);
        for i in 0..ns {
            let p = psi[i].0 * psi[i].0 + psi[i].1 * psi[i].1;
            let x = (i % n) as f64;
            let (sn, cs) = (two_pi * x / n as f64).sin_cos();
            zr += p * cs;
            zi += p * sn;
        }
        let center = (zi.atan2(zr) / two_pi * n as f64).rem_euclid(n as f64);
        out.push(psi);
        centers.push(center);
    }
    // 安定ラベル (0.5 サイト格子スナップ後ソート)
    let snapped: Vec<f64> = centers
        .iter()
        .map(|&c| ((2.0 * c).round() / 2.0).rem_euclid(n as f64))
        .collect();
    let mut ord: Vec<usize> = (0..Q).collect();
    ord.sort_by(|&a, &b| snapped[a].partial_cmp(&snapped[b]).unwrap());
    ord.iter().map(|&i| out[i].clone()).collect()
}

fn yukawa_n(n: usize, la: &[Mode], lb: &[Mode], sig_h: f64) -> M3 {
    let ns = n * n;
    let mut phih = vec![0.0f64; ns];
    for y in 0..n {
        for x in 0..n {
            let dx = (x as f64).min(n as f64 - x as f64);
            let dy = (y as f64).min(n as f64 - y as f64);
            phih[x + y * n] = (-(dx * dx + dy * dy) / (2.0 * sig_h * sig_h)).exp();
        }
    }
    let mut y_out = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..Q {
        for j in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for sx in 0..ns {
                let (ar, ai) = la[i][sx];
                let (br, bi) = lb[j][sx];
                sr += (ar * br + ai * bi) * phih[sx];
                si += (ar * bi - ai * br) * phih[sx];
            }
            y_out[i][j] = (sr, si);
        }
    }
    y_out
}

fn had_prod_perm(y1: &M3, y2: &M3, sf: usize, sg: usize) -> M3 {
    let (pf, pg) = (&PERMS[sf], &PERMS[sg]);
    let mut y = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (p, q) = y1[i][j];
            let (r, s) = y2[pf[i]][pg[j]];
            y[i][j] = (p * r - q * s, p * s + q * r);
        }
    }
    y
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

fn mass_ratios(y: &M3) -> [f64; 2] {
    mass_and_vecs(y).0
}

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

fn cab(v: &M3, i: usize, j: usize) -> f64 {
    (v[i][j].0 * v[i][j].0 + v[i][j].1 * v[i][j].1).sqrt()
}

fn lse(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

/// 幾何対 (トーラス 1 = locs1 のシアー, トーラス 2 = locs2 のシアー) の 10 量コンテスト。
/// v16.7 contest_n の二表版 — 対角化は呼び出し側で済ませ、表を使い回す。
fn contest_pair(
    n: usize,
    locs1: &[Vec<Mode>],
    locs2: &[Vec<Mode>],
    sig_grid: &[f64],
) -> (f64, Vec<f64>, [f64; 2], [f64; 2], M3, [f64; 2]) {
    let j_obs: f64 = 3.08e-5;
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;
    let mut terms10 = Vec::new();
    let mut best_tot = f64::NEG_INFINITY;
    let mut best_ru = [0.0f64; 2];
    let mut best_rd = [0.0f64; 2];
    let mut best_v = [[(0.0f64, 0.0f64); 3]; 3];
    let mut best_re = [0.0f64; 2];
    for &sh in sig_grid {
        let ytab1: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa_n(n, &locs1[ab % NK12], &locs1[ab / NK12], sh))
            .collect();
        let ytab2: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa_n(n, &locs2[ab % NK12], &locs2[ab / NK12], sh))
            .collect();
        let pair_y = |a: usize, b: usize, sf: usize, sg: usize| -> M3 {
            let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
            let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
            had_prod_perm(&ytab1[a1 + b1 * NK12], &ytab2[a2 + b2 * NK12], sf, sg)
        };
        let pair: Vec<([f64; 2], M3)> = (0..nc * nc * 6)
            .map(|m| mass_and_vecs(&pair_y(m % nc, (m / nc) % nc, 0, m / (nc * nc))))
            .collect();
        let mut le = Vec::with_capacity(nc * nc * 36);
        let mut e_best = (f64::NEG_INFINITY, [0.0f64; 2]);
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..nc * nc {
                    let r = mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_));
                    let l = ll2(&r, tgt[4], tgt[5]);
                    le.push(l);
                    if l > e_best.0 {
                        e_best = (l, r);
                    }
                }
            }
        }
        let lnze = lse(&le);
        let mut acc10 = (f64::NEG_INFINITY, 0.0f64);
        let mut q_best = f64::NEG_INFINITY;
        let mut q_ru = [0.0f64; 2];
        let mut q_rd = [0.0f64; 2];
        let mut q_v = [[(0.0f64, 0.0f64); 3]; 3];
        for kq in 0..nc {
            for su in 0..6 {
                for ku in 0..nc {
                    let mu = kq + ku * nc + su * nc * nc;
                    let (ru, vu) = &pair[mu];
                    let llu = ll2(ru, tgt[0], tgt[1]);
                    for sd in 0..6 {
                        for kd in 0..nc {
                            let md = kq + kd * nc + sd * nc * nc;
                            let (rd, vd) = &pair[md];
                            let lld = ll2(rd, tgt[2], tgt[3]);
                            let v = ckm_full(vu, vd);
                            let c = [cab(&v, 0, 1), cab(&v, 1, 2), cab(&v, 0, 2)];
                            let mut ll = llu + lld;
                            for m in 0..3 {
                                let d = c[m].max(1e-300).ln() - tgt[6 + m];
                                ll += -d * d / (2.0 * sigma * sigma) + norm1;
                            }
                            let j = jarlskog(&v);
                            let dj = j.abs().max(1e-300).ln() - j_obs.ln();
                            let ll10 = ll + (-dj * dj / (2.0 * sigma * sigma) + norm1);
                            if ll10 > acc10.0 {
                                acc10.1 = acc10.1 * (acc10.0 - ll10).exp() + 1.0;
                                acc10.0 = ll10;
                            } else {
                                acc10.1 += (ll10 - acc10.0).exp();
                            }
                            if ll10 > q_best {
                                q_best = ll10;
                                q_ru = *ru;
                                q_rd = *rd;
                                q_v = v;
                            }
                        }
                    }
                }
            }
        }
        terms10.push(acc10.0 + acc10.1.ln() + lnze);
        let tot = q_best + e_best.0;
        if tot > best_tot {
            best_tot = tot;
            best_ru = q_ru;
            best_rd = q_rd;
            best_v = q_v;
            best_re = e_best.1;
        }
    }
    let prior_w = 5.0 * (nc as f64).ln() + 4.0 * (6.0f64).ln() + (sig_grid.len() as f64).ln();
    let lnz10 = lse(&terms10) - prior_w;
    (lnz10, terms10, best_ru, best_rd, best_v, best_re)
}

fn main() {
    self_test();
    println!("=== v16.10 谷底平坦域の精密裁定: 4 幾何 × 4 σ_H グリッド (N=36) ===\n");
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

    let n = 36usize;
    let par = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(12);
    let vtd_obs: f64 = 0.0086;
    let vts_obs: f64 = 0.0405;
    let j_obs: f64 = 3.08e-5;
    let geoms: [(usize, usize); 4] = [(2, 2), (2, 3), (1, 3), (3, 3)];
    let g0: Vec<f64> = vec![2.0, 3.0, 4.0, 5.0];
    let g1: Vec<f64> = vec![1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0];
    let g2: Vec<f64> = vec![2.5, 3.5, 4.5, 5.5];
    let g3: Vec<f64> = vec![3.0, 4.0, 5.0, 6.0];
    let grids: [(&str, &Vec<f64>); 4] = [("G0", &g0), ("G1", &g1), ("G2", &g2), ("G3", &g3)];

    // ---- [1] 対角化 (s = 1, 2, 3) ----
    println!(
        "[1] 単一トーラス対角化: 3 シアー × {} Wilson = 36 本 (並列度 {})",
        NK12, par
    );
    let t0 = std::time::Instant::now();
    let svals = [1usize, 2, 3];
    let jobs: Vec<(usize, usize)> = svals
        .iter()
        .flat_map(|&s| (0..NK12).map(move |k| (s, k)))
        .collect();
    let mut raw: std::collections::BTreeMap<(usize, usize), (Vec<Mode>, f64, f64)> =
        std::collections::BTreeMap::new();
    for chunk in jobs.chunks(par) {
        let hs: Vec<_> = chunk
            .iter()
            .map(|&(s, k)| {
                (
                    s,
                    k,
                    std::thread::spawn(move || flux_modes_shear_n(n, k, s)),
                )
            })
            .collect();
        for (s, k, h) in hs {
            raw.insert((s, k), h.join().unwrap());
        }
        println!(
            "    ... {} / 36 本 ({} s)",
            raw.len(),
            t0.elapsed().as_secs()
        );
    }
    let mut locs_by_s: std::collections::BTreeMap<usize, Vec<Vec<Mode>>> =
        std::collections::BTreeMap::new();
    for &s in &svals {
        let spread = (0..NK12).map(|k| raw[&(s, k)].2).fold(0.0f64, f64::max);
        check(
            &format!("s={} の厳密 3 重縮退", s),
            spread < 1e-8,
            format!("幅 {:.1e}", spread),
        );
        locs_by_s.insert(
            s,
            (0..NK12)
                .map(|k| localize_stable(n, &raw[&(s, k)].0))
                .collect(),
        );
    }
    drop(raw);

    // ---- [2] 4 グリッド × 4 幾何 ----
    println!("\n[2] グリッド別の lnZ₁₀ (行 = 幾何, 列 = グリッド):");
    // results[gi][qi] = (lnz10, terms10, MAP 一式)
    let mut lnz = [[0.0f64; 4]; 4]; // [grid][geom]
    let mut scorecards: Vec<String> = Vec::new();
    let mut profs: Vec<(usize, usize, Vec<f64>)> = Vec::new();
    for (gi, (gname, grid)) in grids.iter().enumerate() {
        for (qi, &(s1, s2)) in geoms.iter().enumerate() {
            let (z, terms, ru, rd, v, re) = contest_pair(n, &locs_by_s[&s1], &locs_by_s[&s2], grid);
            lnz[gi][qi] = z;
            if gi == 1 {
                // G1: スコアカードと σ プロファイル
                let jv = jarlskog(&v);
                let qs: Vec<(&str, f64, f64)> = vec![
                    ("m_u/m_t", ru[0].exp(), EPS_OBS[0]),
                    ("m_c/m_t", ru[1].exp(), EPS_OBS[1]),
                    ("m_d/m_b", rd[0].exp(), EPS_OBS[2]),
                    ("m_s/m_b", rd[1].exp(), EPS_OBS[3]),
                    ("m_e/m_τ", re[0].exp(), EPS_OBS[4]),
                    ("m_μ/m_τ", re[1].exp(), EPS_OBS[5]),
                    ("|V_us|", cab(&v, 0, 1), EPS_OBS[6]),
                    ("|V_cb|", cab(&v, 1, 2), EPS_OBS[7]),
                    ("|V_ub|", cab(&v, 0, 2), EPS_OBS[8]),
                    ("|J|", jv.abs(), j_obs),
                    ("|V_td|", cab(&v, 2, 0), vtd_obs),
                    ("|V_ts|", cab(&v, 2, 1), vts_obs),
                ];
                let mut line = format!("    ({},{}) MAP factor: ", s1, s2);
                let mut hits = 0;
                for (name, pred, obs) in &qs {
                    let f = (pred / obs).max(obs / pred);
                    if f <= 5.0 {
                        hits += 1;
                    }
                    line += &format!("{} {:.2} / ", name, f);
                }
                line += &format!(" => {}/12", hits);
                scorecards.push(line);
                profs.push((s1, s2, terms));
            }
        }
        println!(
            "    {} ({} σ 点): (2,2) {:8.3}  (2,3) {:8.3}  (1,3) {:8.3}  (3,3) {:8.3}",
            gname,
            grid.len(),
            lnz[gi][0],
            lnz[gi][1],
            lnz[gi][2],
            lnz[gi][3]
        );
    }

    // ---- ゲート: G0 回帰 ----
    println!("\n[ゲート]");
    for (qi, &(s1, s2)) in geoms.iter().enumerate() {
        check(
            &format!(
                "G0 ({},{}) の回帰: v16.9 の {:.3} と一致 (±0.02)",
                s1, s2, REF_G0[qi]
            ),
            (lnz[0][qi] - REF_G0[qi]).abs() < 0.02,
            format!("lnZ₁₀ = {:.3}", lnz[0][qi]),
        );
    }

    // ---- [3] 判定 ----
    println!("\n[3] グリッド別の順位と Δ:");
    let mut all_top_33 = true;
    for (gi, (gname, _)) in grids.iter().enumerate() {
        let mut ord: Vec<usize> = (0..4).collect();
        ord.sort_by(|&a, &b| lnz[gi][b].partial_cmp(&lnz[gi][a]).unwrap());
        let top = geoms[ord[0]];
        if top != (3, 3) {
            all_top_33 = false;
        }
        println!(
            "    {}: 首位 ({},{})  Δ((3,3)−(1,3)) = {:+.3}  Δ((3,3)−(2,3)) = {:+.3}",
            gname,
            top.0,
            top.1,
            lnz[gi][3] - lnz[gi][2],
            lnz[gi][3] - lnz[gi][1]
        );
    }
    if all_top_33 {
        println!("    => 事前登録 (a): 全グリッドで (3,3) 首位 — Δ0.23 nats は実在、対角優勢は σ_H に頑健。");
    } else {
        println!(
            "    => 事前登録 (b): 首位がグリッドで入れ替わる — 平坦域は装置分解能以下で縮退、"
        );
        println!(
            "       「対角 vs 非対称」の裁定は σ_H の扱いに条件付き (0.23 nats はグリッドの影)。"
        );
    }

    println!("\n[4] G1 (10 σ 点) の全量スコアカード:");
    for l in &scorecards {
        println!("{}", l);
    }
    println!("\n[5] G1 の σ プロファイル terms10(σ) (行 = 幾何, prior 抜き):");
    print!("    σ:      ");
    for s in &g1 {
        print!("{:8.1}", s);
    }
    println!();
    for (s1, s2, terms) in &profs {
        print!("    ({},{}): ", s1, s2);
        for t in terms {
            print!("{:8.2}", t);
        }
        println!();
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.10".into())),
        (
            "grids".into(),
            Json::Arr(
                grids
                    .iter()
                    .enumerate()
                    .map(|(gi, (gname, grid))| {
                        Json::Obj(vec![
                            ("name".into(), Json::Str((*gname).into())),
                            (
                                "sigmas".into(),
                                Json::Arr(grid.iter().map(|&x| Json::Num(x)).collect()),
                            ),
                            (
                                "lnz10".into(),
                                Json::Arr((0..4).map(|qi| Json::Num(lnz[gi][qi])).collect()),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "geoms".into(),
            Json::Arr(
                geoms
                    .iter()
                    .map(|&(a, b)| Json::Arr(vec![Json::Int(a as i64), Json::Int(b as i64)]))
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v1610_sigsens.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 判別は [3] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
