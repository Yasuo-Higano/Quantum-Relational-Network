//! v17.5 τ_im の走査 — 最後の複素構造自由度と、矩形格子の縮退数論
//!
//! シアー s は τ の実部 (τ = (s + i·n_y)/n_x) を走査してきたが、虚部 (アスペクト比
//! n_y/n_x) は全版で 1 に固定されていた。本バイナリは最後の未走査自由度を開ける:
//! 両トーラス (3,3) シアー (τ_re 谷, v16.8)・n_x=36 固定で n_y ∈ {24, 30, 36, 42, 48}
//! → τ_im ∈ {2/3, 5/6, 1, 7/6, 4/3}。
//!
//! 事前設計 (結果より先に固定):
//!   ・磁束は面積量子化 φ = 2πQ/(n_x n_y) — 縮退 Q=3 は連続極限の指数だが、
//!     **矩形格子で厳密かは数論的に非自明** (v16.1 の教訓)。厳密縮退 (<1e-8) を
//!     満たすアスペクトだけがコンテスト資格を持つ (割れたら記録して除外 — それ自体が
//!     論文 §11 の開いた数論への追加データ)。
//!   ・Higgs 幅は面積比例 σ(n_y) = {2,3,4,5}·√(n_y/36) (等方 Gaussian の面積割合を
//!     固定 — 異方幅は限界に記す)。
//!   ・n_y=36 は正方格子の回帰点 (キャッシュから読み、(3,3) の lnZ₁₀ = −21.757)。
//! 事前登録 3 分岐:
//!   (a) 正方 τ_im=1 が最良 (±1 nat 以内で他が並ばない) — アスペクトは非レバー
//!   (b) ある τ_im≠1 が ≥1 nat 改善 — τ モジュライ空間が完全に開く (新レバー)
//!   (c) 中間 — 地図として記録
//! 併せて |V_us|・|J| の MAP factor を追跡 (張力の住所 [v17.4] は τ_im で動くか)。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// 正方格子 (τ_im=1) の回帰アンカー: v16.9 の (3,3)
const REF_SQ33: f64 = -21.756581;
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

fn main() {
    self_test();
    println!("=== v17.5 τ_im の走査: (3,3) シアー固定, n_x=36, n_y ∈ {{24,30,36,42,48}} ===\n");
    println!("事前登録: (a) 正方最良 (非レバー) / (b) τ_im≠1 が ≥1 nat 改善 (新レバー) / (c) 中間 (地図)。");
    println!(
        "資格条件: 厳密 3 重縮退 (<1e-8) — 割れたアスペクトは記録して除外 (格子の縮退数論)。\n"
    );
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

    let nx = 36usize;
    let s = 3usize;
    let nys = [24usize, 30, 36, 42, 48];
    let par = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(12);
    let j_obs: f64 = 3.08e-5;
    let vus_obs: f64 = 0.225;

    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;

    let t0 = std::time::Instant::now();
    let mut rows: Vec<(usize, f64, f64, f64, f64, f64, bool)> = Vec::new(); // (ny, τ_im, spread, lnz10, vus_f, j_f, eligible)
    for &ny in &nys {
        // ---- モード表 (正方はキャッシュ、矩形は計算) ----
        let mut modes_k: Vec<(Vec<Mode>, f64, f64)> = Vec::new();
        if ny == nx {
            let mut ok = true;
            for k in 0..NK12 {
                match cache_load_modes(MODE_TAG, nx, Q, s, k) {
                    Some(v) => modes_k.push(v),
                    None => {
                        ok = false;
                        break;
                    }
                }
            }
            if !ok {
                modes_k.clear();
            }
        }
        if modes_k.is_empty() {
            let mut got: std::collections::BTreeMap<usize, (Vec<Mode>, f64, f64)> =
                std::collections::BTreeMap::new();
            let jobs: Vec<usize> = (0..NK12).collect();
            for chunk in jobs.chunks(par) {
                let hs: Vec<_> = chunk
                    .iter()
                    .map(|&k| (k, std::thread::spawn(move || flux_modes_rect(nx, ny, k, s))))
                    .collect();
                for (k, h) in hs {
                    got.insert(k, h.join().unwrap());
                }
            }
            modes_k = (0..NK12).map(|k| got.remove(&k).unwrap()).collect();
        }
        let spread = modes_k.iter().map(|r| r.2).fold(0.0f64, f64::max);
        let gap = modes_k.iter().map(|r| r.1).fold(f64::INFINITY, f64::min);
        let eligible = spread < 1e-8;
        println!(
            "[n_y={}] τ_im = {:.3}: 縮退幅 {:.1e}, ギャップ {:.4} — {} ({} s)",
            ny,
            ny as f64 / nx as f64,
            spread,
            gap,
            if eligible {
                "資格あり"
            } else {
                "厳密縮退が割れた — コンテスト除外 (格子数論の記録)"
            },
            t0.elapsed().as_secs()
        );
        if !eligible {
            rows.push((
                ny,
                ny as f64 / nx as f64,
                spread,
                f64::NAN,
                f64::NAN,
                f64::NAN,
                false,
            ));
            continue;
        }
        let locs: Vec<Vec<Mode>> = modes_k
            .iter()
            .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
            .collect();
        // ---- 10 量コンテスト (σ は面積比例) ----
        let scale = ((ny as f64) / (nx as f64)).sqrt();
        let sig_grid: Vec<f64> = [2.0f64, 3.0, 4.0, 5.0]
            .iter()
            .map(|s0| s0 * scale)
            .collect();
        let mut terms10 = Vec::new();
        let mut best_ll = f64::NEG_INFINITY;
        let mut best_v = [[(0.0f64, 0.0f64); 3]; 3];
        for &sh in &sig_grid {
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
            for sl in 0..6 {
                for se_ in 0..6 {
                    for ab in 0..nc * nc {
                        let r = mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_));
                        le.push(ll2(&r, tgt[4], tgt[5]));
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
                                if ll10 > best_ll {
                                    best_ll = ll10;
                                    best_v = v;
                                }
                            }
                        }
                    }
                }
            }
            terms10.push(acc10.0 + acc10.1.ln() + lnze);
        }
        let prior_w = 5.0 * (nc as f64).ln() + 4.0 * (6.0f64).ln() + (sig_grid.len() as f64).ln();
        let lnz10 = lse(&terms10) - prior_w;
        let vus = cab(&best_v, 0, 1);
        let vus_f = (vus / vus_obs).max(vus_obs / vus);
        let jm = jarlskog(&best_v).abs();
        let j_f = (jm / j_obs).max(j_obs / jm);
        println!(
            "    lnZ₁₀ = {:.3}, |V_us| factor {:.2}, |J| factor {:.2} ({} s)",
            lnz10,
            vus_f,
            j_f,
            t0.elapsed().as_secs()
        );
        if ny == nx {
            check(
                "正方 (τ_im=1) の回帰: v16.9 の −21.757 (±0.02)",
                (lnz10 - REF_SQ33).abs() < 0.02,
                format!("lnZ₁₀ = {:.3}", lnz10),
            );
        }
        rows.push((ny, ny as f64 / nx as f64, spread, lnz10, vus_f, j_f, true));
    }

    // ---- [2] 地図と判定 ----
    println!("\n[2] τ_im の地図 (シアー (3,3), n_x=36, σ ∝ √面積):");
    println!("    n_y   τ_im    縮退幅     lnZ₁₀     |V_us| f   |J| f");
    for &(ny, ti, sp, z, vf, jf, el) in &rows {
        if el {
            println!(
                "    {:2}   {:.3}   {:.1e}  {:8.3}   {:6.2}   {:6.2}",
                ny, ti, sp, z, vf, jf
            );
        } else {
            println!("    {:2}   {:.3}   {:.1e}   (除外 — 縮退割れ)", ny, ti, sp);
        }
    }
    let elig: Vec<_> = rows.iter().filter(|r| r.6).collect();
    let n_split = rows.len() - elig.len();
    if n_split > 0 {
        println!("    格子数論の記録: {} アスペクトで厳密縮退が割れた (指数は連続極限の性質 — v16.1 と同型)", n_split);
    }
    if elig.len() >= 2 {
        let best = elig
            .iter()
            .max_by(|a, b| a.3.partial_cmp(&b.3).unwrap())
            .unwrap();
        let sq = elig.iter().find(|r| r.0 == nx).unwrap();
        println!(
            "\n    証拠最良: n_y={} (τ_im={:.3}, lnZ₁₀={:.3}) / 正方: {:.3}",
            best.0, best.1, best.3, sq.3
        );
        if best.0 == nx {
            println!("    => 事前登録 (a): 正方 τ_im=1 が最良 — アスペクト比は非レバー。τ 平面の走査が閉じた。");
        } else if best.3 - sq.3 >= 1.0 {
            println!("    => 事前登録 (b): τ_im={:.3} が +{:.2} nats — τ モジュライ空間が完全に開いた (新レバー)。", best.1, best.3 - sq.3);
        } else {
            println!(
                "    => 事前登録 (c): 中間 (最良差 {:+.2} nats < 1) — 地図として記録。",
                best.3 - sq.3
            );
        }
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v17.5".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(ny, ti, sp, z, vf, jf, el)| {
                        Json::Obj(vec![
                            ("ny".into(), Json::Int(ny as i64)),
                            ("tau_im".into(), Json::Num(ti)),
                            ("spread".into(), Json::Num(sp)),
                            ("lnz10".into(), Json::Num(if z.is_nan() { -1e9 } else { z })),
                            (
                                "vus_factor".into(),
                                Json::Num(if vf.is_nan() { -1.0 } else { vf }),
                            ),
                            (
                                "j_factor".into(),
                                Json::Num(if jf.is_nan() { -1.0 } else { jf }),
                            ),
                            ("eligible".into(), Json::Bool(el)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v175_tauim.json", &j.render());
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
