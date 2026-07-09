//! v17.13 orientation の検査 — 共役分枝で CP 角の符号は揃うか
//!
//! v17.9/v17.11 で UT 角の MAP が全窓で負符号 (大きさは測定に接近) と観測された。
//! 幾何の複素共役 (全ホッピング位相 θ → −θ、モード ψ → ψ*) は反ユニタリ対称で、
//! 数学的には V → V*・J → −J・全 CP 角 → 符号反転・全 |·| 厳密不変のはず。
//! 本バイナリは収束点 (36×18; 3,3) でこれを数値検証する:
//!   [1] 共役モードでの lnZ₁₀ (G0) が元と機械精度で一致 (|·| のみの尤度は盲目)
//!   [2] MAP の J・UT 角が符号反転し大きさ不変
//!   [3] 共役分枝の UT 角 (+13.7°, +44.6°) と測定 (+22.2°, +65.9°) の factor
//! 事前登録: (a) [1][2] が機械精度で成立 → orientation は Z₂ の離散データであり、
//! 測定された sign(J) > 0 が分枝を一意に固定する (10 量尤度は |J| なので見えない —
//! J の符号込み 11 量化は不要で、分枝選択として台帳に記す) / (b) 不成立 →
//! 実装か理解の誤り (要調査)。
//! 判定閾値: lnZ 一致 1e-6・|J| 一致 1e-12・角の和がゼロ 1e-9。

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

/// (36×18; 3,3) の G0 コンテスト — locs を与えて lnZ₁₀ と MAP (J, UT 角) を返す
fn contest_apex(nx: usize, ny: usize, locs: &[Vec<Mode>]) -> (f64, f64, f64, f64) {
    let j_obs: f64 = 3.08e-5;
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;
    let scale = ((ny as f64) / (nx as f64)).sqrt();
    let g0 = [2.0f64, 3.0, 4.0, 5.0];
    let mut terms = Vec::new();
    let mut best_tot = f64::NEG_INFINITY;
    let mut map_j = 0.0f64;
    let mut map_b = 0.0f64;
    let mut map_g = 0.0f64;
    for &s0 in &g0 {
        let sh = s0 * scale;
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
        let mut e_best = f64::NEG_INFINITY;
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..nc * nc {
                    let r = mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_));
                    let l = ll2(&r, tgt[4], tgt[5]);
                    le.push(l);
                    if l > e_best {
                        e_best = l;
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
                            let tot = ll10 + e_best;
                            if tot > best_tot {
                                best_tot = tot;
                                map_j = jv;
                                let (b, g) = ut_angles(&v);
                                map_b = b;
                                map_g = g;
                            }
                        }
                    }
                }
            }
        }
        terms.push(acc10.0 + acc10.1.ln() + lnze);
    }
    let prior_w = 5.0 * (nc as f64).ln() + 4.0 * (6.0f64).ln() + (g0.len() as f64).ln();
    (lse(&terms) - prior_w, map_j, map_b, map_g)
}

fn main() {
    self_test();
    println!("=== v17.13 orientation の検査: 共役分枝 (36×18; 3,3) ===\n");
    println!("事前登録: (a) lnZ 一致 (1e-6)・J/角の符号反転かつ大きさ不変 → orientation は");
    println!("          Z₂ 離散データ、測定 sign(J)>0 が分枝を固定 / (b) 不成立 → 要調査\n");
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

    let (nx, ny, s) = (36usize, 18usize, 3usize);
    let par = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(12);
    let (beta_obs, gamma_obs) = (22.2f64, 65.9f64);

    // モード表 (rect キャッシュ)
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
    let locs: Vec<Vec<Mode>> = modes_k
        .iter()
        .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
        .collect();
    // 共役分枝: ψ → ψ* (localize は共役集合上で同一手続き — ラベルの同型性も検証対象)
    let locs_conj: Vec<Vec<Mode>> = modes_k
        .iter()
        .map(|(m, _, _)| {
            let mc: Vec<Mode> = m
                .iter()
                .map(|psi| psi.iter().map(|&(re, im)| (re, -im)).collect())
                .collect();
            localize_stable_rect(nx, ny, &mc)
        })
        .collect();

    let t0 = std::time::Instant::now();
    let (z0, j0, b0, g0) = contest_apex(nx, ny, &locs);
    let (z1, j1, b1, g1) = contest_apex(nx, ny, &locs_conj);
    println!(
        "[1] 元分枝:   lnZ₁₀ = {:.6}, J_MAP = {:+.4e}, β = {:+.1}°, γ = {:+.1}° ({} s)",
        z0,
        j0,
        b0.to_degrees(),
        g0.to_degrees(),
        t0.elapsed().as_secs()
    );
    println!(
        "[2] 共役分枝: lnZ₁₀ = {:.6}, J_MAP = {:+.4e}, β = {:+.1}°, γ = {:+.1}°",
        z1,
        j1,
        b1.to_degrees(),
        g1.to_degrees()
    );

    println!("\n[ゲート]");
    check(
        "lnZ₁₀ の分枝不変 (1e-6)",
        (z0 - z1).abs() < 1e-6,
        format!("|Δ| = {:.1e}", (z0 - z1).abs()),
    );
    check(
        "J の符号反転・大きさ不変 (1e-12)",
        (j0 + j1).abs() < 1e-12 && (j0.abs() - j1.abs()).abs() < 1e-12,
        format!("J₀ = {:+.3e}, J₁ = {:+.3e}", j0, j1),
    );
    check(
        "UT 角の符号反転 (和 = 0, 1e-9)",
        (b0 + b1).abs() < 1e-9 && (g0 + g1).abs() < 1e-9,
        format!("β 和 {:.1e}, γ 和 {:.1e}", (b0 + b1).abs(), (g0 + g1).abs()),
    );

    // ---- [3] 測定分枝の holdout 採点 ----
    let (bp, gp) = if g1 > 0.0 {
        (b1.to_degrees(), g1.to_degrees())
    } else {
        (b0.to_degrees(), g0.to_degrees())
    };
    let fb = (bp / beta_obs).max(beta_obs / bp);
    let fg = (gp / gamma_obs).max(gamma_obs / gp);
    println!("\n[3] 測定 sign(J) > 0 が固定する分枝の UT 角 holdout:");
    println!(
        "    β = {:+.1}° (測定 +22.2°, factor {:.2}), γ = {:+.1}° (測定 +65.9°, factor {:.2})",
        bp, fb, gp, fg
    );
    println!(
        "    => 符号込みで {} — 12 量 + UT 2 角 = {} 級",
        if fb <= 5.0 && fg <= 5.0 && bp > 0.0 && gp > 0.0 {
            "両角とも factor 5 以内"
        } else {
            "factor 5 外あり"
        },
        if fb <= 5.0 && fg <= 5.0 && bp > 0.0 && gp > 0.0 {
            "14/14"
        } else {
            "—"
        }
    );
    if nfail == 0 {
        println!("\n    事前登録 (a) 成立: orientation は Z₂ の離散データ (共役分枝の対) であり、");
        println!("    10 量尤度 (|J|) はこれに盲目。測定された J > 0 が分枝を一意に固定する —");
        println!("    「J の符号込み 11 量化」は不要 (尤度でなく分枝選択として台帳に記す)。");
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v17.13".into())),
        ("lnz_orig".into(), Json::Num(z0)),
        ("lnz_conj".into(), Json::Num(z1)),
        ("j_orig".into(), Json::Num(j0)),
        ("j_conj".into(), Json::Num(j1)),
        ("beta_deg_measured_branch".into(), Json::Num(bp)),
        ("gamma_deg_measured_branch".into(), Json::Num(gp)),
    ]);
    let p = write_artifact("results/v1713_orientation.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 判別は [1]–[3] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
