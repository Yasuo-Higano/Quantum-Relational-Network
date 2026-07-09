//! v16.13 Jeffreys 判別 — β=1 の正体は「スケール不変測度」か
//!
//! v16.12 の峰 β=1 について、e^{β·depth} = Π rᵢ^{−β} なので β=1 の傾けは
//! 「質量比に対する対数一様 (スケール不変・Jeffreys 型) 測度への変換」と読める、
//! という解釈を残した。本バイナリはその直接検定である。データを見ない 3 つの
//! 事前を同じ 10 量証拠にかける:
//!   ・β=1 傾き (v16.12 の峰 — 回帰 −19.7834)
//!   ・1D 深さ平坦化: w = 1/n̂(d) — 配位密度の深さ方向ヒストグラムを平坦化
//!     (β 傾きと同じ 1 次元統計量を使う、最も公平な判別対手)
//!   ・2D 対数一様化: w = 1/n̂(ln r₁, ln r₂) — 両比を独立に平坦化 (完全な
//!     スケール不変化)。クォークは u 対 × d 対で因子化、e 対も同様
//! いずれも配位集団 (模型構造) のみから作る — 観測は使わない (データ盲目)。
//!
//! 事前登録 3 分岐 (Δ は一様基線 −23.1214 との差):
//!   (a) |Δ_flat − Δ_β1| ≤ 0.5 (どちらかの flatten) → β=1 の働きは測度変換と同一
//!       — 「生き残った原理はスケール不変性」の読みを採用
//!   (b) Δ_flat < Δ_β1 − 0.5 (両 flatten) → 傾きは測度変換を超える (深さ方向の
//!       選好が本体) — スケール不変説を棄却
//!   (c) Δ_flat > Δ_β1 + 0.5 → flatten の方が強い — 原理の座を flatten に譲る
//! 副検定: 配位密度 ln n(d) の勾配 (支持域中央) ≈ −1 なら β=1 = 平坦化の機構が
//! 直接確認される。装置変種: ビン幅 {0.25, 0.5, 1.0} で全結論の安定性を見る。
//!
//! 回帰: 一様 −23.1214 (v16.9 marginal)・β=1 −19.7834 (v16.12)。キャッシュ全命中前提。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// 回帰アンカー: 一様 (v16.9 marginal) と β=1 (v16.12)
const REF_UNIFORM: f64 = -23.1214;
const REF_B1: f64 = -19.7834;
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

/// 深さ 1D ヒストグラム (ビン幅 bw) から平坦化重み ln(1/n̂) を返す表を作る。
/// d の達成域のみ (空ビンの配位は存在しない)。
fn flatten_weights_1d(depths: &[f64], bw: f64) -> Vec<f64> {
    let dmin = depths.iter().cloned().fold(f64::INFINITY, f64::min);
    let bin = |d: f64| ((d - dmin) / bw).floor() as usize;
    let nb = depths.iter().map(|&d| bin(d)).max().unwrap() + 1;
    let mut cnt = vec![0usize; nb];
    for &d in depths {
        cnt[bin(d)] += 1;
    }
    depths
        .iter()
        .map(|&d| -((cnt[bin(d)] as f64).ln()))
        .collect()
}

/// 2D (ln r₁, ln r₂) ヒストグラム平坦化重み。
fn flatten_weights_2d(rs: &[[f64; 2]], bw: f64) -> Vec<f64> {
    let xmin = rs.iter().map(|r| r[0]).fold(f64::INFINITY, f64::min);
    let ymin = rs.iter().map(|r| r[1]).fold(f64::INFINITY, f64::min);
    let bx = |r: &[f64; 2]| ((r[0] - xmin) / bw).floor() as usize;
    let by = |r: &[f64; 2]| ((r[1] - ymin) / bw).floor() as usize;
    let nx = rs.iter().map(|r| bx(r)).max().unwrap() + 1;
    let ny = rs.iter().map(|r| by(r)).max().unwrap() + 1;
    let mut cnt = vec![0usize; nx * ny];
    for r in rs {
        cnt[bx(r) + by(r) * nx] += 1;
    }
    rs.iter()
        .map(|r| -((cnt[bx(r) + by(r) * nx] as f64).ln()))
        .collect()
}

fn main() {
    self_test();
    println!("=== v16.13 Jeffreys 判別: β=1 傾き vs 深さ平坦化 vs 対数一様化 (N=36, 21 幾何, 10 量) ===\n");
    println!("事前登録: (a) |Δ_flat−Δ_β1| ≤ 0.5 → 測度変換と同一 / (b) 両 flatten が −0.5 下 → 傾きは測度超え /");
    println!("          (c) flatten が +0.5 上 → flatten 昇格。副検定: ln n(d) 勾配 ≈ −1。ビン幅 {{0.25, 0.5, 1.0}}\n");
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
    let bws = [0.25f64, 0.5, 1.0];
    // 変種: 0 一様 / 1 β=1 / 2..5 flat1D(bw) / 5..8 flat2D(bw) — 計 8
    const NV: usize = 8;
    let vnames = [
        "一様",
        "β=1 傾き",
        "flat1D bw=.25",
        "flat1D bw=.5",
        "flat1D bw=1",
        "flat2D bw=.25",
        "flat2D bw=.5",
        "flat2D bw=1",
    ];

    // ---- [0] モード表 (キャッシュ) ----
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

    // ---- [1] 幾何ループ ----
    println!("\n[1] 幾何ごとの評価 (8 変種) + 配位密度勾配:");
    let geoms: Vec<(usize, usize)> = (0..6usize)
        .flat_map(|a| (a..6usize).map(move |b| (a, b)))
        .collect();
    let mut zg: Vec<[f64; NV]> = Vec::new();
    let mut slopes: Vec<f64> = Vec::new(); // (g,σ) ごとの ln n(d) 勾配 (bw=0.5)
    let t1 = std::time::Instant::now();
    for &(s1, s2) in &geoms {
        let (locs1, locs2) = (&locs_by_s[s1], &locs_by_s[s2]);
        let mut terms: Vec<[f64; NV]> = Vec::new();
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
            // 平坦化重み表 (クォーク対): 1D は深さ, 2D は (ln r₁, ln r₂)
            let depths: Vec<f64> = pair_r.iter().map(|r| -(r[0] + r[1])).collect();
            let w1: Vec<Vec<f64>> = bws
                .iter()
                .map(|&bw| flatten_weights_1d(&depths, bw))
                .collect();
            let w2: Vec<Vec<f64>> = bws
                .iter()
                .map(|&bw| flatten_weights_2d(&pair_r, bw))
                .collect();
            // 配位密度勾配 (bw=0.5, 支持域: 数 ≥ 10 のビン)
            {
                let bw = 0.5;
                let dmin = depths.iter().cloned().fold(f64::INFINITY, f64::min);
                let nb = depths
                    .iter()
                    .map(|&d| ((d - dmin) / bw) as usize)
                    .max()
                    .unwrap()
                    + 1;
                let mut cnt = vec![0usize; nb];
                for &d in &depths {
                    cnt[((d - dmin) / bw) as usize] += 1;
                }
                let pts: Vec<(f64, f64)> = cnt
                    .iter()
                    .enumerate()
                    .filter(|(_, &c)| c >= 10)
                    .map(|(i, &c)| (dmin + (i as f64 + 0.5) * bw, (c as f64).ln()))
                    .collect();
                if pts.len() >= 3 {
                    let xs: Vec<f64> = pts.iter().map(|p| p.0).collect();
                    let ys: Vec<f64> = pts.iter().map(|p| p.1).collect();
                    slopes.push(linfit(&xs, &ys).1);
                }
            }
            // e セクター: 比・平坦化重み
            let mut er: Vec<[f64; 2]> = Vec::with_capacity(nc * nc * 36);
            for sl in 0..6 {
                for se_ in 0..6 {
                    for ab in 0..nc * nc {
                        er.push(mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_)));
                    }
                }
            }
            let edepths: Vec<f64> = er.iter().map(|r| -(r[0] + r[1])).collect();
            let ew1: Vec<Vec<f64>> = bws
                .iter()
                .map(|&bw| flatten_weights_1d(&edepths, bw))
                .collect();
            let ew2: Vec<Vec<f64>> = bws.iter().map(|&bw| flatten_weights_2d(&er, bw)).collect();
            let mut ze_sh = [Acc::new(); NV];
            let mut ne_sh = [Acc::new(); NV];
            for (i, r) in er.iter().enumerate() {
                let l = ll2(r, tgt[4], tgt[5]);
                let depth = edepths[i];
                for v in 0..NV {
                    let s_p = match v {
                        0 => 0.0,
                        1 => depth,
                        2..=4 => ew1[v - 2][i],
                        _ => ew2[v - 5][i],
                    };
                    ze_sh[v].add(s_p + l);
                    ne_sh[v].add(s_p);
                }
            }
            // クォーク五重和 (10 量)
            let mut zq_sh = [Acc::new(); NV];
            let mut nq_sh = [Acc::new(); NV];
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
                                let jv = jarlskog(&v);
                                let dj = jv.abs().max(1e-300).ln() - j_obs.ln();
                                let ll10 = ll + (-dj * dj / (2.0 * sigma * sigma) + norm1);
                                for vv in 0..NV {
                                    let s_p = match vv {
                                        0 => 0.0,
                                        1 => depths[mu] + depths[md],
                                        2..=4 => w1[vv - 2][mu] + w1[vv - 2][md],
                                        _ => w2[vv - 5][mu] + w2[vv - 5][md],
                                    };
                                    zq_sh[vv].add(s_p + ll10);
                                    nq_sh[vv].add(s_p);
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
        }
        let mut z = [0.0f64; NV];
        for v in 0..NV {
            let col: Vec<f64> = terms.iter().map(|r| r[v]).collect();
            z[v] = lse(&col) - (sig_grid.len() as f64).ln();
        }
        println!(
            "    ({},{})  一様 {:8.3}  β=1 {:9.3}  f1D(.5) {:9.3}  f2D(.5) {:9.3}   ({} s)",
            s1,
            s2,
            z[0],
            z[1],
            z[3],
            z[6],
            t1.elapsed().as_secs()
        );
        zg.push(z);
    }

    // ---- 合成とゲート ----
    let total = |vi: usize| -> f64 {
        let col: Vec<f64> = zg.iter().map(|z| z[vi]).collect();
        lse(&col) - (col.len() as f64).ln()
    };
    let tots: Vec<f64> = (0..NV).map(total).collect();
    println!("\n[ゲート]");
    check(
        "一様合計の v16.9 marginal 回帰 (±0.02)",
        (tots[0] - REF_UNIFORM).abs() < 0.02,
        format!("{:.4} vs {:.4}", tots[0], REF_UNIFORM),
    );
    check(
        "β=1 の v16.12 回帰 (±0.02)",
        (tots[1] - REF_B1).abs() < 0.02,
        format!("{:.4} vs {:.4}", tots[1], REF_B1),
    );

    // ---- [2] 判定 ----
    println!("\n[2] 三つ巴 (Δ = lnZ − 一様):");
    for v in 0..NV {
        println!(
            "    {:14}  lnZ = {:9.3}   Δ = {:+7.3}",
            vnames[v],
            tots[v],
            tots[v] - tots[0]
        );
    }
    let d_b1 = tots[1] - tots[0];
    let d_f1: Vec<f64> = (2..5).map(|v| tots[v] - tots[0]).collect();
    let d_f2: Vec<f64> = (5..8).map(|v| tots[v] - tots[0]).collect();
    // 主判定はビン幅 0.5 (事前指定)、他は安定性検査
    let (df1, df2) = (d_f1[1], d_f2[1]);
    let near = |a: f64, b: f64| (a - b).abs() <= 0.5;
    println!(
        "\n    主判定 (bw=0.5): Δ_β1 = {:+.3}, Δ_flat1D = {:+.3}, Δ_flat2D = {:+.3}",
        d_b1, df1, df2
    );
    if near(df1, d_b1) || near(df2, d_b1) {
        println!("    => 事前登録 (a): flatten が β=1 と同じ働き (±0.5) — β=1 の正体は測度変換 (スケール不変性)。");
    } else if df1 > d_b1 + 0.5 || df2 > d_b1 + 0.5 {
        println!("    => 事前登録 (c): flatten が β=1 を上回る — 原理の座は対数一様化 (スケール不変測度) に移る。");
    } else {
        println!("    => 事前登録 (b): flatten は β=1 に届かない — 傾きは測度変換を超える (深さ方向の選好が本体)。");
    }
    println!(
        "    安定性: flat1D Δ = {:?} / flat2D Δ = {:?} (bw = {:?})",
        d_f1.iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Vec<_>>(),
        d_f2.iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Vec<_>>(),
        bws
    );

    // ---- [3] 副検定: 配位密度勾配 ----
    let mut ss = slopes.clone();
    ss.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let med = ss[ss.len() / 2];
    let (lo, hi) = (ss[ss.len() / 10], ss[ss.len() - 1 - ss.len() / 10]);
    println!("\n[3] 副検定: クォーク対の配位密度 ln n(d) の勾配 (84 (g,σ) 点, bw=0.5):");
    println!("    中央値 {:+.3}, 10–90% 帯 [{:+.3}, {:+.3}]", med, lo, hi);
    println!(
        "    => β=1 = 平坦化 なら勾配 ≈ −1 のはず: {}",
        if (-1.3..=-0.7).contains(&med) {
            "整合 (機構レベルで支持)"
        } else {
            "不整合 — β=1 の勝ちは密度勾配の相殺では説明できない"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.13".into())),
        (
            "totals".into(),
            Json::Arr(tots.iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "names".into(),
            Json::Arr(vnames.iter().map(|s| Json::Str((*s).into())).collect()),
        ),
        ("slope_median".into(), Json::Num(med)),
        ("slope_p10".into(), Json::Num(lo)),
        ("slope_p90".into(), Json::Num(hi)),
    ]);
    let p = write_artifact("results/v1613_jeffreys.json", &j.render());
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
