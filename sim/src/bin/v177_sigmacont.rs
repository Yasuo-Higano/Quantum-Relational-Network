//! v17.7 σ_H の連続周辺化 — グリッドの影を消す (PROMPT/3 v18.1)
//!
//! v16.10 は「谷底の順位は σ グリッドの置き方で反転する」ことを示した。測度判定
//! (v18.4) の前提として、σ_H を離散グリッドでなく**連続一様事前で周辺化**し、
//!   ・数値積分が収束するか (格子密度 2 倍で Δ < 0.05 nats)
//!   ・谷底 4 幾何の順位 (v16.10 で縮退) と τ_im 行 (v17.5 で +2.12) が
//!     連続化後も保たれるか
//! を測る。範囲は σ ∈ [1.2, 8.0] (v16.10 の σ プロファイルが両端で −20 nats 級に
//! 沈むことを確認済みの範囲 — 端の寄与は無視可能)。矩形格子は面積スケール
//! σ → σ·√(n_y/36) の同一物理範囲。
//!
//! 事前登録 2 分岐:
//!   (a) 谷底の縮退は連続化後も縮退 (順位差 <1 nat) かつ τ_im の順位は不変
//!       → v16.10/v17.5 の結論はグリッドの影の除去後も立つ
//!   (b) どれかの順位が ≥1 nat で入れ替わる → 当該結論はグリッド人工物と記録
//! 装置: 台形則、密度 2 段 (Δσ=0.4 / 0.2 — 粗い方は細かい方の部分集合なので
//! 一度の評価で両方出る)。回帰: G0 4 点の v16.9/v17.5 アンカー。
//! 矩形モード表は本版で cache_*_rect に格納する (以後の測度判定が分単位になる)。

use uft_sim::*;


const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// G0 回帰アンカー: v16.9 の谷底 4 幾何 + v17.5 の矩形 2 行
const REF_G0: [(usize, usize, usize, f64); 6] = [
    (36, 2, 2, -22.263460),
    (36, 2, 3, -22.256569),
    (36, 1, 3, -21.982785),
    (36, 3, 3, -21.756581),
    (24, 3, 3, -19.641),
    (30, 3, 3, -20.683),
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

fn main() {
    self_test();
    println!("=== v17.7 σ_H の連続周辺化: 谷底 4 幾何 + 矩形 2 行 (σ ∈ [1.2, 8.0] 一様) ===\n");
    println!("事前登録: (a) 谷底縮退は連続化後も縮退 (<1 nat) かつ τ_im 順位不変 / (b) ≥1 nat の順位交代 = グリッド人工物の記録\n");
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
    let par = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(12);
    let j_obs: f64 = 3.08e-5;

    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;

    // 幾何: (ny, s1, s2) — 正方 4 (谷底) + 矩形 2 (τ_im 行)
    let geoms: [(usize, usize, usize); 6] = [
        (36, 2, 2),
        (36, 2, 3),
        (36, 1, 3),
        (36, 3, 3),
        (24, 3, 3),
        (30, 3, 3),
    ];
    // 密グリッド: σ ∈ [1.2, 8.0], Δ=0.1 → 69 点 (D1 = 偶数番目 35 点 Δ=0.2)
    // 初走 (Δ=0.2/0.4) は収束ゲート 0.05 を 0.090 で外した — 基準は緩めず倍密で再走 (開発記録)
    let nsig = 69usize;
    let (s_lo, s_hi) = (1.2f64, 8.0f64);
    let sig_pts: Vec<f64> = (0..nsig)
        .map(|i| s_lo + (s_hi - s_lo) * i as f64 / (nsig - 1) as f64)
        .collect();

    let t0 = std::time::Instant::now();
    // モード表 (正方はキャッシュ、矩形は rect キャッシュ or 計算+保存)
    let mut locs_map: std::collections::BTreeMap<(usize, usize), Vec<Vec<Mode>>> =
        std::collections::BTreeMap::new(); // (ny, s) → 12 Wilson の局在モード
    let mut need: Vec<(usize, usize)> = Vec::new();
    for &(ny, s1, s2) in &geoms {
        for s in [s1, s2] {
            if !need.contains(&(ny, s)) {
                need.push((ny, s));
            }
        }
    }
    for &(ny, s) in &need {
        let mut modes_k: Vec<(Vec<Mode>, f64, f64)> = Vec::new();
        let mut all_hit = true;
        for k in 0..NK12 {
            let got = if ny == nx {
                cache_load_modes(MODE_TAG, nx, Q, s, k)
            } else {
                cache_load_modes_rect(MODE_TAG, nx, ny, Q, s, k)
            };
            match got {
                Some(v) => modes_k.push(v),
                None => {
                    all_hit = false;
                    break;
                }
            }
        }
        if !all_hit {
            modes_k.clear();
            let mut got: std::collections::BTreeMap<usize, (Vec<Mode>, f64, f64)> =
                std::collections::BTreeMap::new();
            let jobs: Vec<usize> = (0..NK12).collect();
            for chunk in jobs.chunks(par) {
                let hs: Vec<_> = chunk
                    .iter()
                    .map(|&k| (k, std::thread::spawn(move || flux_modes_rect(nx, ny, k, s))))
                    .collect();
                for (k, h) in hs {
                    let v = h.join().unwrap();
                    if ny == nx {
                        cache_save_modes(MODE_TAG, nx, Q, s, k, &v.0, v.1, v.2);
                    } else {
                        cache_save_modes_rect(MODE_TAG, nx, ny, Q, s, k, &v.0, v.1, v.2);
                    }
                    got.insert(k, v);
                }
            }
            modes_k = (0..NK12).map(|k| got.remove(&k).unwrap()).collect();
        }
        let spread = modes_k.iter().map(|r| r.2).fold(0.0f64, f64::max);
        check(
            &format!("(ny={}, s={}) の厳密 3 重縮退", ny, s),
            spread < 1e-8,
            format!("幅 {:.1e} ({} s)", spread, t0.elapsed().as_secs()),
        );
        locs_map.insert(
            (ny, s),
            modes_k
                .iter()
                .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
                .collect(),
        );
    }

    // ---- 幾何ごとの σ プロファイルと連続周辺化 ----
    println!("\n[1] σ プロファイル terms10(σ) と連続周辺化:");
    let mut rows: Vec<(usize, usize, usize, f64, f64, f64, f64, f64)> = Vec::new();
    // (ny,s1,s2, lnz_cont(D2), lnz_cont(D1), lnz_g0, σ_map, σ_med)
    for &(ny, s1, s2) in &geoms {
        let scale = ((ny as f64) / (nx as f64)).sqrt();
        let locs1 = &locs_map[&(ny, s1)];
        let locs2 = &locs_map[&(ny, s2)];
        let mut terms = vec![0.0f64; nsig];
        for (isg, &s0) in sig_pts.iter().enumerate() {
            let sh = s0 * scale;
            let ytab1: Vec<M3> = (0..NK12 * NK12)
                .map(|ab| yukawa_rect(nx, ny, &locs1[ab % NK12], &locs1[ab / NK12], sh))
                .collect();
            let ytab2: Vec<M3> = (0..NK12 * NK12)
                .map(|ab| yukawa_rect(nx, ny, &locs2[ab % NK12], &locs2[ab / NK12], sh))
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
                            }
                        }
                    }
                }
            }
            terms[isg] = acc10.0 + acc10.1.ln() + lnze;
        }
        // config 空間の規格化 (σ には連続一様密度)
        let prior_c = 5.0 * (nc as f64).ln() + 4.0 * (6.0f64).ln();
        // D2: 台形則 ln∫ = lse(ln terms + ln w) — w_i = Δ (端は Δ/2)
        let d2 = {
            let dx = (s_hi - s_lo) / (nsig - 1) as f64;
            let lw: Vec<f64> = (0..nsig)
                .map(|i| {
                    let w = if i == 0 || i == nsig - 1 {
                        dx / 2.0
                    } else {
                        dx
                    };
                    terms[i] + w.ln()
                })
                .collect();
            lse(&lw) - (s_hi - s_lo).ln() - prior_c
        };
        // D1: 偶数番目 (Δ=0.4)
        let d1 = {
            let pts: Vec<usize> = (0..nsig).step_by(2).collect();
            let dx = 2.0 * (s_hi - s_lo) / (nsig - 1) as f64;
            let lw: Vec<f64> = pts
                .iter()
                .enumerate()
                .map(|(ii, &i)| {
                    let w = if ii == 0 || ii == pts.len() - 1 {
                        dx / 2.0
                    } else {
                        dx
                    };
                    terms[i] + w.ln()
                })
                .collect();
            lse(&lw) - (s_hi - s_lo).ln() - prior_c
        };
        // G0 回帰 (4 点一様 — 従来定義)
        let g0 = {
            let g0pts = [2.0f64, 3.0, 4.0, 5.0];
            let lw: Vec<f64> = g0pts
                .iter()
                .map(|&s0| {
                    let i = ((s0 - s_lo) / ((s_hi - s_lo) / (nsig - 1) as f64)).round() as usize;
                    terms[i]
                })
                .collect();
            lse(&lw) - (4.0f64).ln() - prior_c
        };
        // σ 事後 (D2 重み)
        let (mut best_i, mut best_t) = (0usize, f64::NEG_INFINITY);
        for (i, &t) in terms.iter().enumerate() {
            if t > best_t {
                best_t = t;
                best_i = i;
            }
        }
        let sig_map = sig_pts[best_i];
        let sig_med = {
            let mx = best_t;
            let ws: Vec<f64> = terms.iter().map(|&t| (t - mx).exp()).collect();
            let tot: f64 = ws.iter().sum();
            let mut acc = 0.0;
            let mut med = sig_pts[nsig - 1];
            for i in 0..nsig {
                acc += ws[i];
                if acc >= 0.5 * tot {
                    med = sig_pts[i];
                    break;
                }
            }
            med
        };
        println!(
            "    ({:2},{},{})  cont(D2) {:8.3}  cont(D1) {:8.3}  G0 {:8.3}  σ* = {:.1} (中央値 {:.1}, 物理幅 {:.4}) ({} s)",
            ny,
            s1,
            s2,
            d2,
            d1,
            g0,
            sig_map,
            sig_med,
            sig_map * ((ny as f64) / (nx as f64)).sqrt() / (nx as f64 * ny as f64).sqrt(),
            t0.elapsed().as_secs()
        );
        rows.push((ny, s1, s2, d2, d1, g0, sig_map, sig_med));
    }

    // ---- ゲート ----
    println!("\n[ゲート]");
    for (ri, &(ny, s1, s2, _, _, g0, _, _)) in rows.iter().enumerate() {
        let (rny, rs1, rs2, refv) = REF_G0[ri];
        assert_eq!((rny, rs1, rs2), (ny, s1, s2));
        check(
            &format!("({},{},{}) の G0 回帰 (±0.05)", ny, s1, s2),
            (g0 - refv).abs() < 0.05,
            format!("{:.3} vs {:.3}", g0, refv),
        );
    }
    let max_conv = rows
        .iter()
        .map(|r| (r.3 - r.4).abs())
        .fold(0.0f64, f64::max);
    check(
        "台形則の収束 (D2 vs D1, 全幾何 <0.05 nats)",
        max_conv < 0.05,
        format!("最大 |Δ| = {:.3}", max_conv),
    );

    // ---- [2] 判定 ----
    println!("\n[2] 連続周辺化後の順位:");
    let mut sq: Vec<_> = rows.iter().filter(|r| r.0 == 36).collect();
    sq.sort_by(|a, b| b.3.partial_cmp(&a.3).unwrap());
    println!(
        "    正方谷底: {}",
        sq.iter()
            .map(|r| format!("({},{}) {:.2}", r.1, r.2, r.3))
            .collect::<Vec<_>>()
            .join(" > ")
    );
    let spread_sq = sq[0].3 - sq[sq.len() - 1].3;
    let r24 = rows.iter().find(|r| r.0 == 24).unwrap();
    let r30 = rows.iter().find(|r| r.0 == 30).unwrap();
    let r36 = rows
        .iter()
        .find(|r| r.0 == 36 && r.1 == 3 && r.2 == 3)
        .unwrap();
    println!(
        "    τ_im 行 (3,3): ny=24 {:.2} / ny=30 {:.2} / ny=36 {:.2}",
        r24.3, r30.3, r36.3
    );
    let tauim_ok = r24.3 > r30.3 && r30.3 > r36.3;
    println!("\n[3] 事前登録判定:");
    println!(
        "    谷底 4 幾何の全幅 = {:.2} nats ({})",
        spread_sq,
        if spread_sq < 1.0 {
            "縮退のまま"
        } else {
            "縮退が破れた"
        }
    );
    println!(
        "    τ_im 順位 (24 > 30 > 36): {}",
        if tauim_ok { "保存" } else { "交代" }
    );
    if spread_sq < 1.0 && tauim_ok {
        println!("    => (a) v16.10 の縮退と v17.5 の τ_im 単調性は、グリッドの影を除いても立つ。");
    } else {
        println!("    => (b) グリッド人工物の記録: 連続化で結論が動いた項目を上記のとおり記す。");
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v17.7".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(ny, s1, s2, d2, d1, g0, smap, smed)| {
                        Json::Obj(vec![
                            ("ny".into(), Json::Int(ny as i64)),
                            ("s1".into(), Json::Int(s1 as i64)),
                            ("s2".into(), Json::Int(s2 as i64)),
                            ("lnz_cont".into(), Json::Num(d2)),
                            ("lnz_cont_coarse".into(), Json::Num(d1)),
                            ("lnz_g0".into(), Json::Num(g0)),
                            ("sigma_map".into(), Json::Num(smap)),
                            ("sigma_med".into(), Json::Num(smed)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v177_sigmacont.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 判別は [2][3] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
