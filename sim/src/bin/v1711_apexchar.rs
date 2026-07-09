//! v17.11 収束点の全量特性 — τ = 1/12 + i/2 のスコアカード・事後帯・測度再検
//!
//! v17.10 で確定した谷の底 (36×18, 対角 (3,3), τ = 1/12 + i/2) を全量で特性化する:
//!   [A] 12 量スコアカード (10 量 MAP): 6 質量比 + 3 CKM + J + holdout |V_td|/|V_ts|
//!       + UT 角 β, γ (PRED-007/008)
//!   [B] 事後帯 (WHist, v17.4 方式): |V_us|, |J|, |V_cb|, |V_td|, |V_ts| の 68% 帯 —
//!       v17.4 (正方谷底) では |V_us| と |J| の帯が測定を外した。底では被覆するか
//!   [C] 測度再検 (v17.9 の正誤表チェック): 底でない窓で下した判定の構造発見
//!       「最良幾何ほど測度補正が不要」は、真の底でどうなるか
//!
//! 事前登録:
//!   [A] 12/12 維持なら 6 幾何連続。UT 角の符号を記録 (orientation の入力)
//!   [B] (a) |V_us|・|J| の 68% 帯が測定を被覆 = 張力は分布レベルでも解消 /
//!       (b) 外す = MAP は良いが事後質量は別 (張力は分布に残る)
//!   [C] (a) flatten Δ ≤ 0 = v17.9 の構造発見が底で確認・強化 /
//!       (b) flatten Δ ≥ +1 = 窓 3 の反転は τ_im=2/3 固有 — v17.9 の還元主張に正誤表
//! ゲート: G0 部分集合が v17.10 の −18.429 を再現 (±0.02)。σ は連続 69 点。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// v17.10 の G0 アンカー (36×18; 3,3)
const REF_G0_APEX: f64 = -18.429;
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

/// 流し込み log-sum-exp 蓄積器
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

/// 2D ヒストグラム平坦化重み (v16.13 から移植, bw 固定 0.5)
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

/// 標本モーメント (μ, Σ) と、二次傾き重み w = +½ (x−μ)ᵀ Σ⁻¹ (x−μ)
fn gauss_weights(rs: &[[f64; 2]]) -> (f64, f64, [f64; 3], Vec<f64>) {
    let n = rs.len() as f64;
    let mu0: f64 = rs.iter().map(|r| r[0]).sum::<f64>() / n;
    let mu1: f64 = rs.iter().map(|r| r[1]).sum::<f64>() / n;
    let (mut s00, mut s11, mut s01) = (0.0f64, 0.0, 0.0);
    for r in rs {
        let (a, b) = (r[0] - mu0, r[1] - mu1);
        s00 += a * a;
        s11 += b * b;
        s01 += a * b;
    }
    s00 /= n;
    s11 /= n;
    s01 /= n;
    let det = s00 * s11 - s01 * s01;
    let (i00, i11, i01) = (s11 / det, s00 / det, -s01 / det);
    let ws = rs
        .iter()
        .map(|r| {
            let (a, b) = (r[0] - mu0, r[1] - mu1);
            0.5 * (i00 * a * a + i11 * b * b + 2.0 * i01 * a * b)
        })
        .collect();
    (mu0 + mu1, s00 + s11 + 2.0 * s01, [s00, s11, s01], ws)
    // 戻り: (μ_d, σ_d² [深さ方向], Σ 成分, 重み)
}

/// 深さ周辺分布の QQ-RMS (99 分位, Gauss(μ_d, σ_d) 対比 — bin なし)
fn qq_rms_depth(depths: &mut [f64], mu_d: f64, var_d: f64) -> f64 {
    depths.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let n = depths.len();
    let sd = var_d.sqrt();
    // 標準正規分位の近似 (Acklam 型の簡易有理近似)
    let inv_norm = |p: f64| -> f64 {
        // Beasley-Springer-Moro
        let a = [
            2.50662823884f64,
            -18.61500062529,
            41.39119773534,
            -25.44106049637,
        ];
        let b = [
            -8.47351093090f64,
            23.08336743743,
            -21.06224101826,
            3.13082909833,
        ];
        let c = [
            0.3374754822726147f64,
            0.9761690190917186,
            0.1607979714918209,
            0.0276438810333863,
            0.0038405729373609,
            0.0003951896511919,
            0.0000321767881768,
            0.0000002888167364,
            0.0000003960315187,
        ];
        let y = p - 0.5;
        if y.abs() < 0.42 {
            let r = y * y;
            y * (((a[3] * r + a[2]) * r + a[1]) * r + a[0])
                / ((((b[3] * r + b[2]) * r + b[1]) * r + b[0]) * r + 1.0)
        } else {
            let mut r = if y > 0.0 { 1.0 - p } else { p };
            r = (-(r.ln())).ln();
            let mut x = c[0];
            let mut rp = 1.0;
            for ci in c.iter().skip(1) {
                rp *= r;
                x += ci * rp;
            }
            if y < 0.0 {
                -x
            } else {
                x
            }
        }
    };
    let mut acc = 0.0;
    for q in 1..100 {
        let p = q as f64 / 100.0;
        let emp = depths[((p * n as f64) as usize).min(n - 1)];
        let thr = mu_d + sd * inv_norm(p);
        acc += (emp - thr) * (emp - thr);
    }
    (acc / 99.0).sqrt()
}

/// 1D 深さヒストグラム平坦化重み (v16.13 から移植)
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

/// 対数ビンの重みつきヒストグラム (v16.4/v17.4 方式)
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

// 変種 (NV=6): 0 UNIFORM / 1 THERMO γ=1 / 2 DEPTH β=.5 / 3 FLAT1D .5 / 4 FLAT2D .5 / 5 GAUSSFLAT
const NV: usize = 6;

fn main() {
    self_test();
    println!("=== v17.11 収束点の全量特性: (36×18; 3,3), τ = 1/12 + i/2 ===\n");
    println!("事前登録: [A] 12/12 維持なら 6 幾何連続 / [B] |V_us|・|J| の 68% 帯被覆の有無 /");
    println!("          [C] flatten Δ≤0 = v17.9 の還元主張を底で確認 / Δ≥+1 = 正誤表\n");
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
    let j_obs: f64 = 3.08e-5;
    let vtd_obs: f64 = 0.0086;
    let vts_obs: f64 = 0.0405;
    let (_beta_obs, _gamma_obs) = (22.2f64.to_radians(), 65.9f64.to_radians());
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;
    let vnames = [
        "UNIFORM",
        "THERMO g=1",
        "DEPTH b=.5",
        "FLAT1D .5",
        "FLAT2D .5",
        "GAUSSFLAT",
    ];

    // ---- モード表 (rect キャッシュ — v17.10 が格納済み) ----
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
    check(
        "厳密 3 重縮退",
        spread < 1e-8,
        format!("幅 {:.1e} (キャッシュ命中 {})", spread, all_hit),
    );
    let locs: Vec<Vec<Mode>> = modes_k
        .iter()
        .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
        .collect();

    // σ 連続格子
    let nsig = 69usize;
    let (s_lo, s_hi) = (1.2f64, 8.0f64);
    let sig_pts: Vec<f64> = (0..nsig)
        .map(|i| s_lo + (s_hi - s_lo) * i as f64 / (nsig - 1) as f64)
        .collect();
    let scale = ((ny as f64) / (nx as f64)).sqrt();

    let t0 = std::time::Instant::now();
    let mut rowsig: Vec<[f64; NV]> = Vec::new();
    let mut qq_max = 0.0f64;
    // 事後帯 (一様事後): |V_us|, |J|, |V_cb|, |V_td|, |V_ts|
    let qnames = ["|V_us|", "|J|", "|V_cb|", "|V_td|", "|V_ts|"];
    let obs5 = [0.225f64, j_obs, 0.041, vtd_obs, vts_obs];
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
    // MAP (一様): クォーク部最良 + e 部最良 (σ ごと → 全体)
    let mut best_tot = f64::NEG_INFINITY;
    let mut map_ru = [0.0f64; 2];
    let mut map_rd = [0.0f64; 2];
    let mut map_re = [0.0f64; 2];
    let mut map_v = [[(0.0f64, 0.0f64); 3]; 3];
    for &s0 in &sig_pts {
        let sh = s0 * scale;
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa_rect(nx, ny, &locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        let pair_y = |a: usize, b: usize, sf: usize, sg: usize| -> M3 {
            let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
            let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
            had_prod_perm(&ytab[a1 + b1 * NK12], &ytab[a2 + b2 * NK12], sf, sg)
        };
        let mut pair_r: Vec<[f64; 2]> = Vec::with_capacity(nc * nc * 6);
        let mut pair_v: Vec<M3> = Vec::with_capacity(nc * nc * 6);
        let mut pair_f: Vec<f64> = Vec::with_capacity(nc * nc * 6);
        for m in 0..nc * nc * 6 {
            let y = pair_y(m % nc, (m / nc) % nc, 0, m / (nc * nc));
            let (r, v) = mass_and_vecs(&y);
            let fro: f64 = y
                .iter()
                .flatten()
                .map(|&(a, b)| a * a + b * b)
                .sum::<f64>()
                .sqrt();
            pair_r.push(r);
            pair_v.push(v);
            pair_f.push(fro.max(1e-300).ln());
        }
        let depths: Vec<f64> = pair_r.iter().map(|r| -(r[0] + r[1])).collect();
        let w1 = flatten_weights_1d(&depths, 0.5);
        let w2 = flatten_weights_2d(&pair_r, 0.5);
        let (mu_d, var_d, _c, wg) = gauss_weights(&pair_r);
        {
            let mut dcopy = depths.clone();
            let qq = qq_rms_depth(&mut dcopy, -mu_d, var_d);
            if qq > qq_max {
                qq_max = qq;
            }
        }
        let mut er: Vec<[f64; 2]> = Vec::with_capacity(nc * nc * 36);
        let mut ef: Vec<f64> = Vec::with_capacity(nc * nc * 36);
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..nc * nc {
                    let y = pair_y(ab % nc, ab / nc, sl, se_);
                    er.push(mass_ratios(&y));
                    let fro: f64 = y
                        .iter()
                        .flatten()
                        .map(|&(p, q)| p * p + q * q)
                        .sum::<f64>()
                        .sqrt();
                    ef.push(fro.max(1e-300).ln());
                }
            }
        }
        let edepths: Vec<f64> = er.iter().map(|r| -(r[0] + r[1])).collect();
        let ew1 = flatten_weights_1d(&edepths, 0.5);
        let ew2 = flatten_weights_2d(&er, 0.5);
        let (_, _, _, ewg) = gauss_weights(&er);
        let mut ze = [Acc::new(); NV];
        let mut ne = [Acc::new(); NV];
        let mut e_best = (f64::NEG_INFINITY, [0.0f64; 2]);
        for (i, r) in er.iter().enumerate() {
            let l = ll2(r, tgt[4], tgt[5]);
            if l > e_best.0 {
                e_best = (l, *r);
            }
            let sps = [0.0, ef[i], 0.5 * edepths[i], ew1[i], ew2[i], ewg[i]];
            for v in 0..NV {
                ze[v].add(sps[v] + l);
                ne[v].add(sps[v]);
            }
        }
        let lnze_uni = ze[0].val() - ((nc * nc * 36) as f64).ln() + ((nc * nc * 36) as f64).ln();
        let lnze = {
            // 一様の e 部 lse (ヒスト重みに使う)
            ze[0].val()
        };
        let _ = lnze_uni;
        let mut zq = [Acc::new(); NV];
        let mut nq = [Acc::new(); NV];
        let mut q_best = (f64::NEG_INFINITY, 0usize);
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
                            let sps = [
                                0.0,
                                pair_f[mu] + pair_f[md],
                                0.5 * (depths[mu] + depths[md]),
                                w1[mu] + w1[md],
                                w2[mu] + w2[md],
                                wg[mu] + wg[md],
                            ];
                            for vv in 0..NV {
                                zq[vv].add(sps[vv] + ll10);
                                nq[vv].add(sps[vv]);
                            }
                            if ll10 > q_best.0 {
                                q_best = (ll10, mu);
                                if ll10 + e_best.0 > best_tot {
                                    best_tot = ll10 + e_best.0;
                                    map_ru = *ru;
                                    map_rd = *rd;
                                    map_re = e_best.1;
                                    map_v = v;
                                }
                            }
                            // 事後帯 (一様重み)
                            let w = ((ll10 + lnze) - shift).exp().min(1e30);
                            let vals = [c[0], jv.abs(), c[1], cab(&v, 2, 0), cab(&v, 2, 1)];
                            for (qi, &x) in vals.iter().enumerate() {
                                hists[qi].add(x.max(1e-300).ln(), w);
                            }
                        }
                    }
                }
            }
        }
        let mut row = [0.0f64; NV];
        for v in 0..NV {
            row[v] = (zq[v].val() - nq[v].val()) + (ze[v].val() - ne[v].val());
        }
        rowsig.push(row);
    }
    println!(
        "[1] σ 連続評価 69 点完了, QQ最大 {:.2} ({} s)",
        qq_max,
        t0.elapsed().as_secs()
    );

    // 連続周辺化と G0 回帰
    let dx = (s_hi - s_lo) / (nsig - 1) as f64;
    let cont = |v: usize| -> f64 {
        let lw: Vec<f64> = (0..nsig)
            .map(|i| {
                let w = if i == 0 || i == nsig - 1 {
                    dx / 2.0
                } else {
                    dx
                };
                rowsig[i][v] + w.ln()
            })
            .collect();
        lse(&lw) - (s_hi - s_lo).ln()
    };
    let prior_c = 0.0; // rowsig は既に config 正規化済み (z−n)
    let _ = prior_c;
    let g0 = {
        let idx = [8usize, 18, 28, 38]; // σ0 = 2,3,4,5
        let lw: Vec<f64> = idx.iter().map(|&i| rowsig[i][0]).collect();
        lse(&lw) - (4.0f64).ln()
    };
    println!("\n[ゲート]");
    check(
        "G0 部分集合の v17.10 回帰 (±0.02)",
        (g0 - REF_G0_APEX).abs() < 0.02,
        format!("{:.3} vs {:.3}", g0, REF_G0_APEX),
    );

    // ---- [A] スコアカード ----
    println!("\n[A] 12 量スコアカード (10 量 MAP) + UT 角:");
    let jv = jarlskog(&map_v);
    let (bua, gua) = ut_angles(&map_v);
    let qs: Vec<(&str, f64, f64)> = vec![
        ("m_u/m_t", map_ru[0].exp(), EPS_OBS[0]),
        ("m_c/m_t", map_ru[1].exp(), EPS_OBS[1]),
        ("m_d/m_b", map_rd[0].exp(), EPS_OBS[2]),
        ("m_s/m_b", map_rd[1].exp(), EPS_OBS[3]),
        ("m_e/m_τ", map_re[0].exp(), EPS_OBS[4]),
        ("m_μ/m_τ", map_re[1].exp(), EPS_OBS[5]),
        ("|V_us|", cab(&map_v, 0, 1), EPS_OBS[6]),
        ("|V_cb|", cab(&map_v, 1, 2), EPS_OBS[7]),
        ("|V_ub|", cab(&map_v, 0, 2), EPS_OBS[8]),
        ("|J|", jv.abs(), j_obs),
        ("|V_td|", cab(&map_v, 2, 0), vtd_obs),
        ("|V_ts|", cab(&map_v, 2, 1), vts_obs),
    ];
    let mut hits = 0;
    let mut pct = 0;
    for (name, pred, obs) in &qs {
        let f = (pred / obs).max(obs / pred);
        if f <= 5.0 {
            hits += 1;
        }
        if f <= 1.1 {
            pct += 1;
        }
        println!(
            "    {:8} pred {:.4e}  obs {:.3e}  factor {:.2}",
            name, pred, obs, f
        );
    }
    println!("    => {}/12 (factor 5 以内), うち 10% 級 {} 個", hits, pct);
    println!(
        "    UT 角 (holdout): β = {:.1}° (測定 22.2°), γ = {:.1}° (測定 65.9°) — 符号 {}",
        bua.to_degrees(),
        gua.to_degrees(),
        if bua < 0.0 {
            "負 (orientation 問題継続)"
        } else {
            "正"
        }
    );

    // ---- [B] 事後帯 ----
    println!("\n[B] 事後帯 (68%) と被覆 — v17.4 の正方谷底では |V_us|・|J| が外れていた:");
    let mut cover_vus_j = [false; 2];
    for qi in 0..5 {
        let (q16, q50, q84) = (
            hists[qi].quant(0.16),
            hists[qi].quant(0.5),
            hists[qi].quant(0.84),
        );
        let cov = q16 <= obs5[qi] && obs5[qi] <= q84;
        if qi == 0 {
            cover_vus_j[0] = cov;
        }
        if qi == 1 {
            cover_vus_j[1] = cov;
        }
        println!(
            "    {:7}  [{:.4e}, {:.4e}, {:.4e}]  測定 {:.3e}  68% 帯{}",
            qnames[qi],
            q16,
            q50,
            q84,
            obs5[qi],
            if cov {
                "は測定を含む ✓"
            } else {
                "は測定を外す ✗"
            }
        );
    }
    println!(
        "    => 事前登録 [B]: {}",
        if cover_vus_j[0] && cover_vus_j[1] {
            "(a) |V_us| と |J| の両帯が被覆 — 張力は分布レベルでも解消"
        } else if cover_vus_j[0] || cover_vus_j[1] {
            "部分 — 片方のみ被覆 (どちらかの張力が分布に残る)"
        } else {
            "(b) 両方外す — MAP は良いが事後質量は別の場所"
        }
    );

    // ---- [C] 測度再検 ----
    println!("\n[C] 測度再検 (底の窓, 連続 σ) — Δ = variant − 一様:");
    let base = cont(0);
    for v in 1..NV {
        let d = cont(v) - base;
        let note = if v == 5 && qq_max >= 0.6 {
            " (資格外 — QQ ≥ 0.6)"
        } else {
            ""
        };
        println!("    {:12}  Δ = {:+.2}{}", vnames[v], d, note);
    }
    let d_f1 = cont(3) - base;
    let d_f2 = cont(4) - base;
    println!(
        "    => 事前登録 [C]: {}",
        if d_f1 <= 0.0 && d_f2 <= 0.0 {
            "(a) flatten は底で非正 — v17.9 の還元主張 (幾何が測度を不要にする) を底で確認"
        } else if d_f1 >= 1.0 || d_f2 >= 1.0 {
            "(b) flatten が底で ≥+1 — v17.9 の窓 3 反転は τ_im=2/3 固有 (正誤表対象)"
        } else {
            "中間 (0 < Δ < 1) — 補正の残余は僅少と記録"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v17.11".into())),
        ("lnz_cont_uniform".into(), Json::Num(base)),
        ("lnz_g0".into(), Json::Num(g0)),
        ("qq_max".into(), Json::Num(qq_max)),
        (
            "scorecard".into(),
            Json::Arr(
                qs.iter()
                    .map(|(n, p, o)| {
                        Json::Obj(vec![
                            ("q".into(), Json::Str((*n).into())),
                            ("pred".into(), Json::Num(*p)),
                            ("obs".into(), Json::Num(*o)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("ut_beta_deg".into(), Json::Num(bua.to_degrees())),
        ("ut_gamma_deg".into(), Json::Num(gua.to_degrees())),
        (
            "bands".into(),
            Json::Arr(
                (0..5)
                    .map(|qi| {
                        Json::Obj(vec![
                            ("q".into(), Json::Str(qnames[qi].into())),
                            ("q16".into(), Json::Num(hists[qi].quant(0.16))),
                            ("q50".into(), Json::Num(hists[qi].quant(0.5))),
                            ("q84".into(), Json::Num(hists[qi].quant(0.84))),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "measure_delta".into(),
            Json::Arr((1..NV).map(|v| Json::Num(cont(v) - base)).collect()),
        ),
    ]);
    let p = write_artifact("results/v1711_apexchar.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 判別は [A]–[C] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
