//! v16.12 β の走査 — 最深原理の崩壊点を探す
//!
//! v16.11 で階層最大化 (Depth) 原理が初めて生存した (+2.31 nats, β 2 点で単調増加)。
//! 素朴な予想は「β を上げれば、いつか事前が観測より深い config に集中しすぎて
//! 崩壊する」(v15.5 の rect ではそうなった)。本バイナリは β ∈ {0.25..8} の 9 点と
//! 解析的な β→∞ 端点 (幾何・σ ごとの最深 config への hard selection) で
//! 生存曲線 Δ(β) を測る。
//!
//! 事前登録 2 分岐:
//!   (a) 崩壊あり — ある β* で Δ(β) < +1 (生存線割れ)。峰 β_peak と β* を記録し、
//!       Depth は「有限 β の窓の中の物語」として台帳に置く
//!   (b) 崩壊なし — 走査上限 β=8 と β→∞ 端点まで Δ ≥ +1 を維持。
//!       「シアー族の最深 config は観測に座る」という強い構造主張に昇格候補
//!   いずれでも谷底首位 ((1,3) — v16.11) の β 安定性を記録する。
//!
//! 装置ゲート: uniform 合計 = v16.9 marginal −23.121・β=0.25/0.5 の合計が v16.11 の
//! −21.605 / −20.370 を再現 (±0.02)・厳密縮退。モード表は disk キャッシュ (v16.11
//! 導入) — 初のフル命中運転になるはず (対角化 0 本・分単位)。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// v16.11 の回帰アンカー (二層合計)
const REF_UNIFORM: f64 = -23.1214;
const REF_B025: f64 = -21.6053;
const REF_B05: f64 = -20.3702;
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

/// 流し込み log-sum-exp 蓄積器 (v15.5 から移植)
#[derive(Clone, Copy)]
struct Acc {
    m: f64,
    s: f64,
}
impl Acc {
    fn new() -> Self {
        Acc {
            m: f64::NEG_INFINITY,
            s: 0.0,
        }
    }
    fn add(&mut self, x: f64) {
        if x > self.m {
            self.s = self.s * (self.m - x).exp() + 1.0;
            self.m = x;
        } else {
            self.s += (x - self.m).exp();
        }
    }
    fn val(&self) -> f64 {
        self.m + self.s.ln()
    }
}

fn main() {
    self_test();
    println!("=== v16.12 β の走査: 最深原理の崩壊点 (N=36, 21 幾何, 10 量) ===\n");
    println!("事前登録: (a) ある β* で Δ<+1 → 崩壊あり (峰と β* を記録) /");
    println!("          (b) β=8 と β→∞ 端点まで Δ≥+1 → 最深原理は hard selection まで生存\n");
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
    let betas = [0.25f64, 0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 6.0, 8.0];
    const NV: usize = 10; // 0: 一様 / 1..=9: β

    // ---- [0] モード表 (キャッシュ) ----
    let _t0 = std::time::Instant::now();
    let mut raw: std::collections::BTreeMap<(usize, usize), (Vec<Mode>, f64, f64)> =
        std::collections::BTreeMap::new();
    let mut misses: Vec<(usize, usize)> = Vec::new();
    for s in 0..6usize {
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
        "[0] モード表: キャッシュ命中 {} / 72, 対角化 {} 本 (並列度 {})",
        72 - misses.len(),
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
    let mut locs_by_s: Vec<Vec<Vec<Mode>>> = Vec::new();
    for s in 0..6usize {
        let spread = (0..NK12).map(|k| raw[&(s, k)].2).fold(0.0f64, f64::max);
        check(
            &format!("s={} の厳密 3 重縮退", s),
            spread < 1e-8,
            format!("幅 {:.1e}", spread),
        );
        locs_by_s.push(
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

    // ---- [1] 幾何ごとの Z_β(g) と hard-max 端点 ----
    println!("\n[1] 幾何ごとの評価 (β 9 点 + 一様 + β→∞ 端点):");
    let geoms: Vec<(usize, usize)> = (0..6usize)
        .flat_map(|a| (a..6usize).map(move |b| (a, b)))
        .collect();
    let mut zg: Vec<[f64; NV]> = Vec::new();
    let mut zg_inf: Vec<f64> = Vec::new();
    let t1 = std::time::Instant::now();
    for &(s1, s2) in &geoms {
        let (locs1, locs2) = (&locs_by_s[s1], &locs_by_s[s2]);
        let mut terms: Vec<[f64; NV]> = Vec::new();
        let mut terms_inf: Vec<f64> = Vec::new();
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
            let mut pair_r: Vec<[f64; 2]> = Vec::with_capacity(nc * nc * 6);
            let mut pair_v: Vec<M3> = Vec::with_capacity(nc * nc * 6);
            for m in 0..nc * nc * 6 {
                let y = pair_y(m % nc, (m / nc) % nc, 0, m / (nc * nc));
                let (r, v) = mass_and_vecs(&y);
                pair_r.push(r);
                pair_v.push(v);
            }
            // e セクター
            let mut ze_sh = [Acc::new(); NV];
            let mut ne_sh = [Acc::new(); NV];
            let mut e_inf = (f64::NEG_INFINITY, f64::NEG_INFINITY); // (深さ, その ll)
            for sl in 0..6 {
                for se_ in 0..6 {
                    for ab in 0..nc * nc {
                        let r = mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_));
                        let l = ll2(&r, tgt[4], tgt[5]);
                        let depth = -(r[0] + r[1]);
                        for v in 0..NV {
                            let s_p = if v == 0 { 0.0 } else { betas[v - 1] * depth };
                            ze_sh[v].add(s_p + l);
                            ne_sh[v].add(s_p);
                        }
                        if depth > e_inf.0 || (depth == e_inf.0 && l > e_inf.1) {
                            e_inf = (depth, l);
                        }
                    }
                }
            }
            // クォーク五重和 (10 量)
            let mut zq_sh = [Acc::new(); NV];
            let mut nq_sh = [Acc::new(); NV];
            let mut q_inf = (f64::NEG_INFINITY, f64::NEG_INFINITY);
            for kq in 0..nc {
                for su in 0..6 {
                    for ku in 0..nc {
                        let mu = kq + ku * nc + su * nc * nc;
                        let ru = &pair_r[mu];
                        let vu = &pair_v[mu];
                        let llu = ll2(ru, tgt[0], tgt[1]);
                        for sd in 0..6 {
                            for kd in 0..nc {
                                let md = kq + kd * nc + sd * nc * nc;
                                let rd = &pair_r[md];
                                let lld = ll2(rd, tgt[2], tgt[3]);
                                let v = ckm_full(vu, &pair_v[md]);
                                let c = [cab(&v, 0, 1), cab(&v, 1, 2), cab(&v, 0, 2)];
                                let mut ll = llu + lld;
                                for m in 0..3 {
                                    let d = c[m].max(1e-300).ln() - tgt[6 + m];
                                    ll += -d * d / (2.0 * sigma * sigma) + norm1;
                                }
                                let j = jarlskog(&v);
                                let dj = j.abs().max(1e-300).ln() - j_obs.ln();
                                let ll10 = ll + (-dj * dj / (2.0 * sigma * sigma) + norm1);
                                let depth = -(ru[0] + ru[1] + rd[0] + rd[1]);
                                for vv in 0..NV {
                                    let s_p = if vv == 0 { 0.0 } else { betas[vv - 1] * depth };
                                    zq_sh[vv].add(s_p + ll10);
                                    nq_sh[vv].add(s_p);
                                }
                                if depth > q_inf.0 || (depth == q_inf.0 && ll10 > q_inf.1) {
                                    q_inf = (depth, ll10);
                                }
                            }
                        }
                    }
                }
            }
            let mut row = [0.0f64; NV];
            for v in 0..NV {
                row[v] = (zq_sh[v].val() - nq_sh[v].val()) + (ze_sh[v].val() - ne_sh[v].val());
            }
            terms.push(row);
            terms_inf.push(q_inf.1 + e_inf.1);
        }
        let mut z = [0.0f64; NV];
        for v in 0..NV {
            let col: Vec<f64> = terms.iter().map(|r| r[v]).collect();
            z[v] = lse(&col) - (sig_grid.len() as f64).ln();
        }
        let zinf = lse(&terms_inf) - (sig_grid.len() as f64).ln();
        println!(
            "    ({},{})  一様 {:8.3}  β=1 {:8.3}  β=4 {:8.3}  β=8 {:9.3}  β→∞ {:9.3}   ({} s)",
            s1,
            s2,
            z[0],
            z[3],
            z[7],
            z[9],
            zinf,
            t1.elapsed().as_secs()
        );
        zg.push(z);
        zg_inf.push(zinf);
    }

    // ---- 合成 ----
    let total = |col: &[f64]| -> f64 { lse(col) - (col.len() as f64).ln() };
    let uni_tot = total(&zg.iter().map(|z| z[0]).collect::<Vec<_>>());
    let mut tot_b = [0.0f64; 9];
    for (bi, t) in tot_b.iter_mut().enumerate() {
        *t = total(&zg.iter().map(|z| z[bi + 1]).collect::<Vec<_>>());
    }
    let tot_inf = total(&zg_inf);

    // ---- ゲート ----
    println!("\n[ゲート]");
    check(
        "一様合計の v16.9 marginal 回帰 (±0.02)",
        (uni_tot - REF_UNIFORM).abs() < 0.02,
        format!("{:.4} vs {:.4}", uni_tot, REF_UNIFORM),
    );
    check(
        "β=0.25 の v16.11 回帰 (±0.02)",
        (tot_b[0] - REF_B025).abs() < 0.02,
        format!("{:.4} vs {:.4}", tot_b[0], REF_B025),
    );
    check(
        "β=0.5 の v16.11 回帰 (±0.02)",
        (tot_b[1] - REF_B05).abs() < 0.02,
        format!("{:.4} vs {:.4}", tot_b[1], REF_B05),
    );

    // ---- [2] 生存曲線と判定 ----
    println!("\n[2] 生存曲線 Δ(β) = lnZ_β − lnZ_一様 (生存線 +1):");
    let valley: [(usize, usize); 4] = [(2, 2), (2, 3), (1, 3), (3, 3)];
    let vidx: Vec<usize> = valley
        .iter()
        .map(|g| geoms.iter().position(|x| x == g).unwrap())
        .collect();
    println!("    β        lnZ_β      Δ       谷底首位");
    let mut curve: Vec<(f64, f64, f64, (usize, usize))> = Vec::new();
    for bi in 0..9 {
        let vals: Vec<f64> = vidx.iter().map(|&gi| zg[gi][bi + 1]).collect();
        let top = valley[vals
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0];
        println!(
            "    {:4.2}  {:9.3}  {:+7.3}    ({},{})",
            betas[bi],
            tot_b[bi],
            tot_b[bi] - uni_tot,
            top.0,
            top.1
        );
        curve.push((betas[bi], tot_b[bi], tot_b[bi] - uni_tot, top));
    }
    {
        let vals: Vec<f64> = vidx.iter().map(|&gi| zg_inf[gi]).collect();
        let top = valley[vals
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
            .unwrap()
            .0];
        println!(
            "    ∞     {:9.3}  {:+7.3}    ({},{})   (hard selection: 幾何・σ ごと最深 config)",
            tot_inf,
            tot_inf - uni_tot,
            top.0,
            top.1
        );
    }
    let peak = curve
        .iter()
        .cloned()
        .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
        .unwrap();
    let collapse = curve.iter().find(|c| c.0 > peak.0 && c.2 < 1.0);
    let d_inf = tot_inf - uni_tot;
    println!("\n    峰: β = {:.2} (Δ = {:+.3})", peak.0, peak.2);
    match collapse {
        Some(c) => {
            println!(
                "    => 事前登録 (a): 崩壊あり — β* = {:.2} で Δ = {:+.3} < +1。Depth は有限 β の窓の中の物語。",
                c.0, c.2
            );
        }
        None => {
            if d_inf >= 1.0 {
                println!(
                    "    => 事前登録 (b): 崩壊なし — β=8 と β→∞ 端点 (Δ = {:+.3}) まで生存線の上。",
                    d_inf
                );
                println!("       「シアー族の最深 config は観測に座る」— 強い構造主張に昇格候補。");
            } else {
                println!(
                    "    => 有限走査では崩壊なし (β≤8 で Δ≥+1) だが β→∞ 端点は Δ = {:+.3} < +1 —",
                    d_inf
                );
                println!("       崩壊点は β ∈ (8, ∞) にある。峰を過ぎた減衰の始まりとして記録。");
            }
        }
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.12".into())),
        ("uniform".into(), Json::Num(uni_tot)),
        (
            "betas".into(),
            Json::Arr(betas.iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "lnz_beta".into(),
            Json::Arr(tot_b.iter().map(|&x| Json::Num(x)).collect()),
        ),
        ("lnz_inf".into(), Json::Num(tot_inf)),
        (
            "zg_inf".into(),
            Json::Arr(zg_inf.iter().map(|&x| Json::Num(x)).collect()),
        ),
    ]);
    let p = write_artifact("results/v1612_depthscan.json", &j.render());
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
