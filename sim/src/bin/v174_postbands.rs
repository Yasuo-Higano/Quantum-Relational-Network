//! v17.4 谷底の事後帯 — 縮退は分布のレベルでも縮退か
//!
//! v16.10 は谷底 4 幾何が「証拠 (lnZ) の分解能以下で縮退」することを示し、
//! QRN-YUK-024 の限界 (ii) に「事後帯の重なりは未評価」と記した。本バイナリは
//! その返済である: 谷底の代表 2 幾何 — 証拠最良の対角 (3,3) と、生存原理 (測度
//! 補正クラス, v16.11–13) が選ぶ非対称 (1,3) — について、10 量事後の下での
//! 5 つの CKM 系量 (|V_us|, |J|, |V_cb|, |V_td|, |V_ts|) の事後帯 (16/50/84%) を
//! v16.4 の WHist 方式 (シフト付き重み・対数ビン) で測る。
//!
//! 事前登録 2 分岐:
//!   (a) |V_us| の 68% 帯が 2 幾何で分離 → 幾何は事後分布で判別可能 —
//!       将来の |V_us| 側の精密化 (or 理論の精密化) がどちらかを落とす、という
//!       具体的な反証チャネルを登録する
//!   (b) 重なる → 縮退は分布のレベルでも縮退 (証拠だけでなく置き場所も同じ)
//! いずれでも holdout (|V_td|, |V_ts|) の 68% 帯が測定を含むかを両幾何で採点する。
//!
//! 装置: モード表キャッシュ (s=1,3 のみ読む)。回帰は v16.9 の lnZ₁₀ 2 本 (±0.02)。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// v16.9 の回帰アンカー
const REF_13: f64 = -21.982785;
const REF_33: f64 = -21.756581;
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

/// 対数ビンの重みつきヒストグラム (v16.4 WHist 方式)
struct WHist {
    lo: f64,
    hi: f64,
    bins: Vec<f64>,
    wsum: f64,
}
impl WHist {
    fn new(lo: f64, hi: f64, nb: usize) -> Self {
        WHist {
            lo,
            hi,
            bins: vec![0.0; nb],
            wsum: 0.0,
        }
    }
    fn add(&mut self, lnx: f64, w: f64) {
        let nb = self.bins.len();
        let x = lnx.clamp(self.lo, self.hi - 1e-9);
        let b = ((x - self.lo) / (self.hi - self.lo) * nb as f64) as usize;
        self.bins[b.min(nb - 1)] += w;
        self.wsum += w;
    }
    fn quant(&self, q: f64) -> f64 {
        let mut acc = 0.0;
        for (i, h) in self.bins.iter().enumerate() {
            acc += h;
            if acc >= q * self.wsum {
                let nb = self.bins.len() as f64;
                return (self.lo + (i as f64 + 0.5) / nb * (self.hi - self.lo)).exp();
            }
        }
        self.hi.exp()
    }
}

fn main() {
    self_test();
    println!("=== v17.4 谷底の事後帯: (3,3) vs (1,3) の 5 量事後 (10 量尤度, N=36) ===\n");
    println!("事前登録: (a) |V_us| の 68% 帯が分離 → 幾何は事後で判別可能 (反証チャネル登録) /");
    println!("          (b) 重なる → 縮退は分布レベルでも縮退。holdout 帯の被覆も両幾何で採点。\n");
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
    let sig_grid = [2.0f64, 3.0, 4.0, 5.0];
    let par = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(12);
    let j_obs: f64 = 3.08e-5;
    let vtd_obs: f64 = 0.0086;
    let vts_obs: f64 = 0.0405;
    let vus_obs: f64 = 0.225;
    let vcb_obs: f64 = 0.041;

    // ---- [0] モード表 (キャッシュ; s = 1, 3 のみ) ----
    let mut raw: std::collections::BTreeMap<(usize, usize), (Vec<Mode>, f64, f64)> =
        std::collections::BTreeMap::new();
    let mut misses: Vec<(usize, usize)> = Vec::new();
    for &s in &[1usize, 3] {
        for k in 0..NK12 {
            match cache_load_modes(MODE_TAG, n, Q, s, k) {
                Some(v) => {
                    raw.insert((s, k), v);
                }
                None => misses.push((s, k)),
            }
        }
    }
    println!(
        "[0] モード表: キャッシュ命中 {} / 24, 対角化 {} 本 (並列度 {})",
        24 - misses.len(),
        misses.len(),
        par
    );
    for chunk in misses.chunks(par) {
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
            let v = h.join().unwrap();
            cache_save_modes(MODE_TAG, n, Q, s, k, &v.0, v.1, v.2);
            raw.insert((s, k), v);
        }
    }
    let mut locs_by_s: std::collections::BTreeMap<usize, Vec<Vec<Mode>>> =
        std::collections::BTreeMap::new();
    for &s in &[1usize, 3] {
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

    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;

    // ---- [1] 幾何ごとの 5 量事後帯 ----
    let geoms: [(usize, usize, f64); 2] = [(3, 3, REF_33), (1, 3, REF_13)];
    let qnames = ["|V_us|", "|J|", "|V_cb|", "|V_td|", "|V_ts|"];
    let obs = [vus_obs, j_obs, vcb_obs, vtd_obs, vts_obs];
    // 事後帯の記録: [幾何][量] = (q16, q50, q84, MAP)
    let mut bands: Vec<[[f64; 4]; 5]> = Vec::new();
    for &(s1, s2, refz) in &geoms {
        let (locs1, locs2) = (&locs_by_s[&s1], &locs_by_s[&s2]);
        let mut terms10 = Vec::new();
        // 量ごとのヒストグラム (|J| だけ範囲が深い)
        let mut hists: Vec<WHist> = (0..5)
            .map(|qi| {
                if qi == 1 {
                    WHist::new(-20.0, -5.0, 600)
                } else {
                    WHist::new(-12.0, 0.0, 480)
                }
            })
            .collect();
        let shift = -12.0f64;
        let mut map_ll = f64::NEG_INFINITY;
        let mut map_vals = [0.0f64; 5];
        let t0 = std::time::Instant::now();
        for &sh in &sig_grid {
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
                                let w = ((ll10 + lnze) - shift).exp().min(1e30);
                                let vals = [c[0], jv.abs(), c[1], cab(&v, 2, 0), cab(&v, 2, 1)];
                                for (qi, &x) in vals.iter().enumerate() {
                                    hists[qi].add(x.max(1e-300).ln(), w);
                                }
                                if ll10 > map_ll {
                                    map_ll = ll10;
                                    map_vals = vals;
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
        check(
            &format!(
                "({},{}) の lnZ₁₀ 回帰: v16.9 の {:.3} (±0.02)",
                s1, s2, refz
            ),
            (lnz10 - refz).abs() < 0.02,
            format!("lnZ₁₀ = {:.3} ({} s)", lnz10, t0.elapsed().as_secs()),
        );
        let mut b = [[0.0f64; 4]; 5];
        for qi in 0..5 {
            b[qi] = [
                hists[qi].quant(0.16),
                hists[qi].quant(0.5),
                hists[qi].quant(0.84),
                map_vals[qi],
            ];
        }
        bands.push(b);
    }

    // ---- [2] 事後帯の表と判定 ----
    println!("\n[2] 10 量事後の 5 量帯 (16% / 中央値 / 84% / MAP — 測定との比):");
    for (gi, &(s1, s2, _)) in geoms.iter().enumerate() {
        println!("    ({},{})", s1, s2);
        for qi in 0..5 {
            let [q16, q50, q84, m] = bands[gi][qi];
            let cover = q16 <= obs[qi] && obs[qi] <= q84;
            println!(
                "      {:7}  [{:.4e}, {:.4e}, {:.4e}]  MAP {:.4e}  測定 {:.3e}  68% 帯{}",
                qnames[qi],
                q16,
                q50,
                q84,
                m,
                obs[qi],
                if cover {
                    "は測定を含む ✓"
                } else {
                    "は測定を外す ✗"
                }
            );
        }
    }
    // 事前登録判定: |V_us| 帯の分離
    let (a16, _a50, a84) = (bands[0][0][0], bands[0][0][1], bands[0][0][2]);
    let (b16, _b50, b84) = (bands[1][0][0], bands[1][0][1], bands[1][0][2]);
    let overlap = a16.max(b16) <= a84.min(b84);
    println!("\n[3] 事前登録判定 — |V_us| の 68% 帯:");
    println!(
        "    (3,3): [{:.4}, {:.4}] / (1,3): [{:.4}, {:.4}]",
        a16, a84, b16, b84
    );
    if overlap {
        println!(
            "    => (b) 帯は重なる — 縮退は分布のレベルでも縮退 (証拠だけでなく置き場所も同じ)。"
        );
    } else {
        println!("    => (a) 帯は分離 — 幾何は事後で判別可能。|V_us| 側の精密化がどちらかを落とす");
        println!("       (反証チャネルとして登録: 測定 0.225 がどちらの帯に居るかが裁く)。");
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v17.4".into())),
        (
            "geoms".into(),
            Json::Arr(
                geoms
                    .iter()
                    .enumerate()
                    .map(|(gi, &(s1, s2, _))| {
                        Json::Obj(vec![
                            ("s1".into(), Json::Int(s1 as i64)),
                            ("s2".into(), Json::Int(s2 as i64)),
                            (
                                "bands".into(),
                                Json::Arr(
                                    (0..5)
                                        .map(|qi| {
                                            Json::Obj(vec![
                                                ("q".into(), Json::Str(qnames[qi].into())),
                                                ("q16".into(), Json::Num(bands[gi][qi][0])),
                                                ("q50".into(), Json::Num(bands[gi][qi][1])),
                                                ("q84".into(), Json::Num(bands[gi][qi][2])),
                                                ("map".into(), Json::Num(bands[gi][qi][3])),
                                            ])
                                        })
                                        .collect(),
                                ),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v174_postbands.json", &j.render());
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
