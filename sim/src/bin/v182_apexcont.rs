//! v18.2 収束点の連続極限検証 — τ = 1/12 + i/2 を固定して格子を倍にする
//!
//! 第十八期の収束点 (36×18; s=3) は τ = 1/12 + i/2 で 12/12・|V_us| 0.5%・γ 誤差内を
//! 達成した。本バイナリは第十九期の開幕として、その最初の連続極限テストを行う:
//! (n_x, n_y, s) = (72, 36, 6) — 複素構造 τ = 6/72 + i·36/72 = 1/12 + i/2 は厳密に
//! 同一、磁束は面積量子化 (φ = 2πQ/2592)、Higgs は物理幅固定 (σ = {2,3,4,5}·√2)。
//! 72×36 は 72² の半分の格子点 (5184² 行列) なので稠密 jacobi で届く —
//! v16.7 (N=18→36) と同じ判別を、今度は収束点そのものに対して行う。
//!
//! 事前登録 3 分岐:
//!   (a) lnZ₁₀ が維持/改善し MAP factor (|V_us|, |J|) が ±20% で不変
//!       → 収束点は格子の偶然でなく連続極限の構造 (v16.7 の再現をより高い山で)
//!   (b) factor が系統的に動く → apex の格子効果を記録 (どの量が動くかが情報)
//!   (c) 厳密縮退が割れる → 格子の縮退数論の新データ (72×36 は n_x=2n_y の有理点)
//! 比較アンカー (36×18; G0): lnZ₁₀ = −18.429, |J| factor 1.06, UT γ = ±66.8°。
//! 判定は lnZ の絶対値でなく「両格子の差」— 模型空間は同型なので直接比較可能。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const MODE_TAG: u64 = 1; // flux_modes_shear_n の構成タグ (v16.2 系ホッピング)

const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

type M3 = [[(f64, f64); 3]; 3];
type Mode = Vec<(f64, f64)>; // 長さ n_x·n_y

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

/// 矩形格子 n_x × n_y・シアー s の磁束トーラス最低 Q モード (2 成分 Dirac 型)。
/// 正方版 (flux_modes_shear_n) の一般化 — φ = 2πQ/(n_x n_y)。
fn flux_modes_rect(nx: usize, ny: usize, k_half: usize, s: usize) -> (Vec<Mode>, f64, f64) {
    let ns = nx * ny;
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
    let idx = |x: usize, y: usize| x + y * nx;
    for x in 0..nx {
        for y in 0..ny {
            let th_y = phi * x as f64 + wl;
            if y == ny - 1 {
                addhop(&mut a, idx(x, y), idx((x + s) % nx, 0), th_y, m, ns);
            } else {
                addhop(&mut a, idx(x, y), idx(x, y + 1), th_y, m, ns);
            }
            let th_x = if x == nx - 1 {
                -phi * (nx as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % nx, y), th_x, m, ns);
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

/// 局在化 + 安定ラベル (矩形版 — X̂ = e^{2πix/n_x} は x 方向のみ使う)
fn localize_stable_rect(nx: usize, ny: usize, modes: &[Mode]) -> Vec<Mode> {
    let ns = nx * ny;
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut ure = [[0.0f64; 3]; 3];
    let mut uim = [[0.0f64; 3]; 3];
    for a in 0..Q {
        for b in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for i in 0..ns {
                let x = (i % nx) as f64;
                let (sn, cs) = (two_pi * x / nx as f64).sin_cos();
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
            let x = (i % nx) as f64;
            let (sn, cs) = (two_pi * x / nx as f64).sin_cos();
            zr += p * cs;
            zi += p * sn;
        }
        let center = (zi.atan2(zr) / two_pi * nx as f64).rem_euclid(nx as f64);
        out.push(psi);
        centers.push(center);
    }
    let snapped: Vec<f64> = centers
        .iter()
        .map(|&c| ((2.0 * c).round() / 2.0).rem_euclid(nx as f64))
        .collect();
    let mut ord: Vec<usize> = (0..Q).collect();
    ord.sort_by(|&a, &b| snapped[a].partial_cmp(&snapped[b]).unwrap());
    ord.iter().map(|&i| out[i].clone()).collect()
}

/// 湯川重なり (矩形版 — 周期距離は方向別)
fn yukawa_rect(nx: usize, ny: usize, la: &[Mode], lb: &[Mode], sig_h: f64) -> M3 {
    let ns = nx * ny;
    let mut phih = vec![0.0f64; ns];
    for y in 0..ny {
        for x in 0..nx {
            let dx = (x as f64).min(nx as f64 - x as f64);
            let dy = (y as f64).min(ny as f64 - y as f64);
            phih[x + y * nx] = (-(dx * dx + dy * dy) / (2.0 * sig_h * sig_h)).exp();
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

/// UT 角 (β, γ) [rad] — β = arg(−V_cd V*_cb / (V_td V*_tb)), γ = arg(−V_ud V*_ub / (V_cd V*_cb))
fn ut_angles(v: &M3) -> (f64, f64) {
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let conj = |a: (f64, f64)| (a.0, -a.1);
    let div = |a: (f64, f64), b: (f64, f64)| {
        let d = b.0 * b.0 + b.1 * b.1;
        ((a.0 * b.0 + a.1 * b.1) / d, (a.1 * b.0 - a.0 * b.1) / d)
    };
    let neg = |a: (f64, f64)| (-a.0, -a.1);
    let beta = {
        let num = neg(mul(v[1][0], conj(v[1][2])));
        let den = mul(v[2][0], conj(v[2][2]));
        let r = div(num, den);
        r.1.atan2(r.0)
    };
    let gamma = {
        let num = neg(mul(v[0][0], conj(v[0][2])));
        let den = mul(v[1][0], conj(v[1][2]));
        let r = div(num, den);
        r.1.atan2(r.0)
    };
    (beta, gamma)
}

/// G0 コンテスト (12 量スコアカード込み): (lnZ₁₀, MAP の ru/rd/re/V)
fn contest_full(
    nx: usize,
    ny: usize,
    locs: &[Vec<Mode>],
) -> (f64, [f64; 2], [f64; 2], [f64; 2], M3) {
    let j_obs: f64 = 3.08e-5;
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;
    let scale = ((nx * ny) as f64 / (36.0 * 18.0)).sqrt() * ((18.0f64) / 36.0).sqrt();
    // 物理幅固定: σ = {2,3,4,5} · √(面積 / 36²) — (36×18) は ·√(1/2), (72×36) は ·√2·√(1/2)=·1
    // 一貫性のため基準を明示: σ_base · √(nx·ny) / 36
    let g0: Vec<f64> = [2.0f64, 3.0, 4.0, 5.0]
        .iter()
        .map(|s0| s0 * ((nx * ny) as f64).sqrt() / 36.0)
        .collect();
    let _ = scale;
    let mut terms = Vec::new();
    let mut best_tot = f64::NEG_INFINITY;
    let mut map_ru = [0.0f64; 2];
    let mut map_rd = [0.0f64; 2];
    let mut map_re = [0.0f64; 2];
    let mut map_v = [[(0.0f64, 0.0f64); 3]; 3];
    for &sh in &g0 {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa_rect(nx, ny, &locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        let pair_y = |a: usize, b: usize, sf: usize, sg: usize| -> M3 {
            let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
            let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
            had_prod_perm(&ytab[a1 + b1 * NK12], &ytab[a2 + b2 * NK12], sf, sg)
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
                            let jv = jarlskog(&v);
                            let dj = jv.abs().max(1e-300).ln() - j_obs.ln();
                            let ll10 = ll + (-dj * dj / (2.0 * sigma * sigma) + norm1);
                            if ll10 > acc10.0 {
                                acc10.1 = acc10.1 * (acc10.0 - ll10).exp() + 1.0;
                                acc10.0 = ll10;
                            } else {
                                acc10.1 += (ll10 - acc10.0).exp();
                            }
                            let tot = ll10 + e_best.0;
                            if tot > best_tot {
                                best_tot = tot;
                                map_ru = *ru;
                                map_rd = *rd;
                                map_re = e_best.1;
                                map_v = v;
                            }
                        }
                    }
                }
            }
        }
        terms.push(acc10.0 + acc10.1.ln() + lnze);
    }
    let prior_w = 5.0 * (nc as f64).ln() + 4.0 * (6.0f64).ln() + (g0.len() as f64).ln();
    (lse(&terms) - prior_w, map_ru, map_rd, map_re, map_v)
}

fn main() {
    self_test();
    println!("=== v18.2 収束点の連続極限検証: (72×36; s=6), τ = 1/12 + i/2 厳密固定 ===\n");
    println!("事前登録: (a) lnZ 維持/改善 + factor ±20% 不変 = 連続極限の構造 /");
    println!("          (b) factor 系統移動 = apex の格子効果 / (c) 縮退割れ = 数論\n");
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

    let par = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(12);
    let j_obs: f64 = 3.08e-5;
    let vtd_obs: f64 = 0.0086;
    let vts_obs: f64 = 0.0405;
    // 比較アンカー (36×18; 3,3) — v17.10/v17.13 の G0 値
    const REF_LNZ_APEX18: f64 = -18.429;

    let cases: [(usize, usize, usize); 2] = [(36, 18, 3), (72, 36, 6)];
    let mut rows = Vec::new();
    for &(nx, ny, s) in &cases {
        let t0 = std::time::Instant::now();
        let mut modes_k: Vec<(Vec<Mode>, f64, f64)> = Vec::new();
        let mut all_hit = true;
        for k in 0..NK12 {
            match cache_load_modes_rect(MODE_TAG, nx, ny, Q, s, k) {
                Some(v) => modes_k.push(v),
                None => {
                    all_hit = false;
                    break;
                }
            }
        }
        if !all_hit {
            modes_k.clear();
            let jobs: Vec<usize> = (0..NK12).collect();
            let mut got: std::collections::BTreeMap<usize, (Vec<Mode>, f64, f64)> =
                std::collections::BTreeMap::new();
            for chunk in jobs.chunks(par) {
                let hs: Vec<_> = chunk
                    .iter()
                    .map(|&k| (k, std::thread::spawn(move || flux_modes_rect(nx, ny, k, s))))
                    .collect();
                for (k, h) in hs {
                    let v = h.join().unwrap();
                    cache_save_modes_rect(MODE_TAG, nx, ny, Q, s, k, &v.0, v.1, v.2);
                    got.insert(k, v);
                }
            }
            modes_k = (0..NK12).map(|k| got.remove(&k).unwrap()).collect();
        }
        let spread = modes_k.iter().map(|r| r.2).fold(0.0f64, f64::max);
        let gap = modes_k.iter().map(|r| r.1).fold(f64::INFINITY, f64::min);
        check(
            &format!("({}×{}; s={}) の厳密 3 重縮退", nx, ny, s),
            spread < 1e-8,
            format!(
                "幅 {:.1e}, ギャップ {:.4} (キャッシュ {}) ({} s)",
                spread,
                gap,
                all_hit,
                t0.elapsed().as_secs()
            ),
        );
        let locs: Vec<Vec<Mode>> = modes_k
            .iter()
            .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
            .collect();
        let (lnz10, ru, rd, re, v) = contest_full(nx, ny, &locs);
        let jv = jarlskog(&v);
        let (bua, gua) = ut_angles(&v);
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
        let mut hits = 0;
        let mut fv = Vec::new();
        print!("    factor: ");
        for (name, pred, obs) in &qs {
            let f = (pred / obs).max(obs / pred);
            if f <= 5.0 {
                hits += 1;
            }
            fv.push(f);
            print!("{} {:.2} / ", name, f);
        }
        println!();
        println!(
            "    ({}×{})  lnZ₁₀ = {:.3}  {}/12  UT: β {:+.1}°, γ {:+.1}°  ({} s)",
            nx,
            ny,
            lnz10,
            hits,
            bua.to_degrees(),
            gua.to_degrees(),
            t0.elapsed().as_secs()
        );
        rows.push((nx, ny, lnz10, fv, hits, bua.to_degrees(), gua.to_degrees()));
    }

    // ゲート: 36×18 の回帰
    check(
        "(36×18) の回帰: v17.10 の −18.429 (±0.02)",
        (rows[0].2 - REF_LNZ_APEX18).abs() < 0.02,
        format!("lnZ₁₀ = {:.3}", rows[0].2),
    );

    // ---- [2] 判定 ----
    println!("\n[2] 連続極限判定 (τ = 1/12 + i/2 厳密固定):");
    let dz = rows[1].2 - rows[0].2;
    println!(
        "    lnZ₁₀: {:.3} (36×18) → {:.3} (72×36)  Δ = {:+.3}",
        rows[0].2, rows[1].2, dz
    );
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
        "|J|",
        "|V_td|",
        "|V_ts|",
    ];
    let mut max_shift = 0.0f64;
    let mut max_name = "";
    for i in 0..12 {
        let shift = (rows[1].3[i] / rows[0].3[i]).max(rows[0].3[i] / rows[1].3[i]);
        if shift > max_shift {
            max_shift = shift;
            max_name = names[i];
        }
    }
    println!(
        "    factor の最大移動: {} ({:.2}×) / 12 量 {}/12 → {}/12",
        max_name, max_shift, rows[0].4, rows[1].4
    );
    println!(
        "    UT 角: β {:+.1}° → {:+.1}°, γ {:+.1}° → {:+.1}°",
        rows[0].5, rows[1].5, rows[0].6, rows[1].6
    );
    if dz > -0.5 && max_shift < 1.2 && rows[1].4 == 12 {
        println!("    => 事前登録 (a): 収束点は連続極限の構造 — lnZ 維持/改善・factor ±20% 内・12/12 維持。");
    } else if rows[1].4 < 12 || max_shift >= 1.2 {
        println!("    => 事前登録 (b): factor の系統移動あり — apex の格子効果として記録 (最大: {} {:.2}×)。", max_name, max_shift);
    } else {
        println!(
            "    => 中間 — lnZ の低下 {:+.3} を記録 (factor は保存)。",
            dz
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v18.2".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|(nx, ny, z, fv, h, b, g)| {
                        Json::Obj(vec![
                            ("nx".into(), Json::Int(*nx as i64)),
                            ("ny".into(), Json::Int(*ny as i64)),
                            ("lnz10".into(), Json::Num(*z)),
                            (
                                "factors".into(),
                                Json::Arr(fv.iter().map(|&x| Json::Num(x)).collect()),
                            ),
                            ("hits".into(), Json::Int(*h as i64)),
                            ("ut_beta_deg".into(), Json::Num(*b)),
                            ("ut_gamma_deg".into(), Json::Num(*g)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v182_apexcont.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 判別は [2] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
