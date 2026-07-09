//! v17.9 測度判定 — 事前登録された候補を一斉に秤へ (PROMPT/3 v18.4)
//!
//! measures.yml (v17.6) に登録した測度候補を、固定済みプロトコルで一斉判定する:
//!   窓 3 つ ((36;1,3) 非対称 / (36;3,3) 対角 / (24;3,3) 対角×τ_im=2/3)・
//!   σ_H 連続周辺化 (v17.7 の収束済み設定)・全候補同時同一プロトコル。
//! 候補: MSR-UNIFORM (基線) / MSR-MDL / MSR-THERMO (γ marg ln2) /
//!   MSR-DEPTH (β∈{.25,.5,1,1.5} marg ln4 — S1 は marg でのみ充足) /
//!   MSR-FLAT1D・MSR-FLAT2D (bw∈{.25,.5,1} — 2D は MSR-JEFFREYS 経路(b)) /
//!   MSR-GAUSSFLAT (資格条件 QQ<0.6 — 非適合窓では未定義)。
//! 判定 (事前登録): 生存 Δ≥+1 (窓平均) かつ S1–S3 / 棄却 Δ≤−1 or S 違反 / 未決。
//!   S2 = 3 窓の符号安定 + 変種安定 / S3 = UT 角 holdout (PRED-007/008 —
//!   一様基線の MAP 角を S3 採点より前に凍結・印字する機械的順序)。
//! ゲート: 一様の連続 lnZ が v17.7 の 3 値を再現 (±0.02)。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// v17.7 の連続周辺化アンカー (一様)
const REF_CONT: [(usize, usize, usize, f64); 3] =
    [(36, 1, 3, -23.23), (36, 3, 3, -23.11), (24, 3, 3, -20.11)];
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

// 変種の並び (NV=15):
// 0 UNIFORM / 1 MDL / 2-3 THERMO γ=.5,1 / 4-7 DEPTH β=.25,.5,1,1.5 /
// 8-10 FLAT1D bw=.25,.5,1 / 11-13 FLAT2D bw=.25,.5,1 / 14 GAUSSFLAT
const NV: usize = 15;

fn main() {
    self_test();
    println!("=== v17.9 測度判定: 候補 6+1 種 × 3 窓 × 連続 σ (measures.yml プロトコル) ===\n");
    println!("事前登録: 生存 Δ≥+1 (窓平均) かつ S1–S3 / 棄却 Δ≤−1 or S 違反 / 中間未決");
    println!("  S2 = 3 窓の符号安定 + 変種安定 / S3 = UT 角 holdout (基線を冒頭凍結)\n");
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
    let (beta_obs, gamma_obs) = (22.2f64.to_radians(), 65.9f64.to_radians());
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;
    let geoms: [(usize, usize, usize); 3] = [(36, 1, 3), (36, 3, 3), (24, 3, 3)];
    let gam_th = [0.5f64, 1.0];
    let beta_dep = [0.25f64, 0.5, 1.0, 1.5];
    let bws = [0.25f64, 0.5, 1.0];
    let vnames = [
        "UNIFORM",
        "MDL",
        "THERMO g=.5",
        "THERMO g=1",
        "DEPTH b=.25",
        "DEPTH b=.5",
        "DEPTH b=1",
        "DEPTH b=1.5",
        "FLAT1D .25",
        "FLAT1D .5",
        "FLAT1D 1",
        "FLAT2D .25",
        "FLAT2D .5",
        "FLAT2D 1",
        "GAUSSFLAT",
    ];
    let mdl_comp = |z: usize| -> f64 {
        if z == 0 {
            (0.5f64).ln()
        } else {
            (0.1f64).ln()
        }
    };
    let mdl_k = |a: usize| -> f64 { mdl_comp(a % 6) + mdl_comp(a / 6) };
    let mdl_s3 = |s: usize| -> f64 {
        if s == 0 {
            (0.5f64).ln()
        } else {
            (0.1f64).ln()
        }
    };

    // ---- モード表 ----
    let mut locs_map: std::collections::BTreeMap<(usize, usize), Vec<Vec<Mode>>> =
        std::collections::BTreeMap::new();
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
        locs_map.insert(
            (ny, s),
            modes_k
                .iter()
                .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
                .collect(),
        );
    }

    // σ 連続格子 (v17.7 と同一)
    let nsig = 69usize;
    let (s_lo, s_hi) = (1.2f64, 8.0f64);
    let sig_pts: Vec<f64> = (0..nsig)
        .map(|i| s_lo + (s_hi - s_lo) * i as f64 / (nsig - 1) as f64)
        .collect();

    let t0 = std::time::Instant::now();
    let mut zc: Vec<[f64; NV]> = Vec::new();
    let mut qq_win: Vec<f64> = Vec::new();
    let mut map_ang: Vec<[(f64, f64); NV]> = Vec::new();
    for &(ny, s1, s2) in &geoms {
        let scale = ((ny as f64) / (nx as f64)).sqrt();
        let locs1 = &locs_map[&(ny, s1)];
        let locs2 = &locs_map[&(ny, s2)];
        let mut rowsig: Vec<[f64; NV]> = Vec::new();
        let mut qq_max = 0.0f64;
        let mut mbest = [(f64::NEG_INFINITY, (0.0f64, 0.0f64)); NV];
        for &s0 in &sig_pts {
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
            let w1: Vec<Vec<f64>> = bws
                .iter()
                .map(|&bw| flatten_weights_1d(&depths, bw))
                .collect();
            let w2: Vec<Vec<f64>> = bws
                .iter()
                .map(|&bw| flatten_weights_2d(&pair_r, bw))
                .collect();
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
            let ew1: Vec<Vec<f64>> = bws
                .iter()
                .map(|&bw| flatten_weights_1d(&edepths, bw))
                .collect();
            let ew2: Vec<Vec<f64>> = bws.iter().map(|&bw| flatten_weights_2d(&er, bw)).collect();
            let (_, _, _, ewg) = gauss_weights(&er);
            let mut ze = [Acc::new(); NV];
            let mut ne = [Acc::new(); NV];
            for sl in 0..6 {
                for se_ in 0..6 {
                    let base = (sl * 6 + se_) * nc * nc;
                    for ab in 0..nc * nc {
                        let i = base + ab;
                        let r = &er[i];
                        let l = ll2(r, tgt[4], tgt[5]);
                        for v in 0..NV {
                            let s_p = match v {
                                0 => 0.0,
                                1 => mdl_k(ab % nc) + mdl_k(ab / nc) + mdl_s3(sl) + mdl_s3(se_),
                                2 | 3 => gam_th[v - 2] * ef[i],
                                4..=7 => beta_dep[v - 4] * edepths[i],
                                8..=10 => ew1[v - 8][i],
                                11..=13 => ew2[v - 11][i],
                                _ => ewg[i],
                            };
                            ze[v].add(s_p + l);
                            ne[v].add(s_p);
                        }
                    }
                }
            }
            let mut zq = [Acc::new(); NV];
            let mut nq = [Acc::new(); NV];
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
                                        1 => {
                                            mdl_k(kq)
                                                + mdl_k(ku)
                                                + mdl_k(kd)
                                                + mdl_s3(su)
                                                + mdl_s3(sd)
                                        }
                                        2 | 3 => gam_th[vv - 2] * (pair_f[mu] + pair_f[md]),
                                        4..=7 => beta_dep[vv - 4] * (depths[mu] + depths[md]),
                                        8..=10 => w1[vv - 8][mu] + w1[vv - 8][md],
                                        11..=13 => w2[vv - 11][mu] + w2[vv - 11][md],
                                        _ => wg[mu] + wg[md],
                                    };
                                    zq[vv].add(s_p + ll10);
                                    nq[vv].add(s_p);
                                    let tot = s_p + ll10;
                                    if tot > mbest[vv].0 {
                                        mbest[vv] = (tot, ut_angles(&v));
                                    }
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
        let mut z = [0.0f64; NV];
        let dx = (s_hi - s_lo) / (nsig - 1) as f64;
        for v in 0..NV {
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
            z[v] = lse(&lw) - (s_hi - s_lo).ln();
        }
        println!(
            "[窓 ({};{},{})] 一様 {:.3}, QQ最大 {:.2} ({} s)",
            ny,
            s1,
            s2,
            z[0],
            qq_max,
            t0.elapsed().as_secs()
        );
        zc.push(z);
        qq_win.push(qq_max);
        map_ang.push(std::array::from_fn(|v| mbest[v].1));
    }

    // ---- ゲート: 一様の v17.7 回帰 ----
    println!("\n[ゲート]");
    for (wi, &(_ny, _s1, _s2, refv)) in REF_CONT.iter().enumerate() {
        check(
            &format!("窓 {} の一様連続値 v17.7 回帰 (±0.02)", wi),
            (zc[wi][0] - refv).abs() < 0.02,
            format!("{:.3} vs {:.3}", zc[wi][0], refv),
        );
    }

    // ---- S3 基線の凍結 (採点より前に印字) ----
    println!("\n[S3 基線凍結] 一様 MAP の UT 角 (測定: β=22.2°, γ=65.9°):");
    for (wi, &(ny, s1, s2)) in geoms.iter().enumerate() {
        let (b, g) = map_ang[wi][0];
        println!(
            "    ({};{},{}): β = {:.1}°, γ = {:.1}°",
            ny,
            s1,
            s2,
            b.to_degrees(),
            g.to_degrees()
        );
    }

    // ---- [2] 判定表 ----
    println!("\n[2] 窓ごとの Δ (variant − 一様):");
    print!("    {:12}", "variant");
    for &(ny, s1, s2) in &geoms {
        print!("  ({};{},{})", ny, s1, s2);
    }
    println!();
    for v in 1..NV {
        print!("    {:12}", vnames[v]);
        for wi in 0..3 {
            if v == 14 && qq_win[wi] >= 0.6 {
                print!("  {:>9}", "(資格外)");
            } else {
                print!("  {:+9.2}", zc[wi][v] - zc[wi][0]);
            }
        }
        println!();
    }

    let marg = |ids: &[usize], wi: usize| -> f64 {
        let vals: Vec<f64> = ids.iter().map(|&v| zc[wi][v]).collect();
        lse(&vals) - (ids.len() as f64).ln()
    };
    let fam: [(&str, Vec<usize>); 6] = [
        ("MSR-MDL", vec![1]),
        ("MSR-THERMO", vec![2, 3]),
        ("MSR-DEPTH", vec![4, 5, 6, 7]),
        ("MSR-FLAT1D", vec![8, 9, 10]),
        ("MSR-FLAT2D(JEFF)", vec![11, 12, 13]),
        ("MSR-GAUSSFLAT", vec![14]),
    ];
    println!("\n[3] 家族 marginal (Occam 込み) の Δ と S2:");
    println!("    家族              窓1(1,3)  窓2(3,3)  窓3(24)   窓平均   S2 (符号)");
    let mut fam_rows: Vec<(String, [f64; 3], f64, bool)> = Vec::new();
    for (name, ids) in &fam {
        let mut ds = [0.0f64; 3];
        let mut valid = [true; 3];
        for wi in 0..3 {
            if *name == "MSR-GAUSSFLAT" && qq_win[wi] >= 0.6 {
                valid[wi] = false;
                ds[wi] = f64::NAN;
            } else {
                ds[wi] = marg(ids, wi) - zc[wi][0];
            }
        }
        let dv: Vec<f64> = (0..3).filter(|&i| valid[i]).map(|i| ds[i]).collect();
        let mean = dv.iter().sum::<f64>() / dv.len() as f64;
        let sign_stable = dv.iter().all(|&d| d > 0.0) || dv.iter().all(|&d| d < 0.0);
        let f = |x: f64, ok: bool| -> String {
            if ok {
                format!("{:+.2}", x)
            } else {
                "資格外".into()
            }
        };
        println!(
            "    {:16}  {:>8}  {:>8}  {:>8}  {:+7.2}   {}",
            name,
            f(ds[0], valid[0]),
            f(ds[1], valid[1]),
            f(ds[2], valid[2]),
            mean,
            if sign_stable {
                "符号安定"
            } else {
                "符号不安定 (S2 疑義)"
            }
        );
        fam_rows.push((name.to_string(), ds, mean, sign_stable));
    }

    // ---- [4] S3: UT 角 ----
    println!("\n[4] S3 — 測度重み MAP の UT 角 (基線比):");
    for wi in 0..3 {
        let (b0, g0) = map_ang[wi][0];
        let d0 = (b0 - beta_obs).abs() + (g0 - gamma_obs).abs();
        println!(
            "    窓 {} 基線: β {:.1}° γ {:.1}° (測定距離 {:.3} rad)",
            wi,
            b0.to_degrees(),
            g0.to_degrees(),
            d0
        );
        for v in 1..NV {
            if v == 14 && qq_win[wi] >= 0.6 {
                continue;
            }
            let (b, g) = map_ang[wi][v];
            let d = (b - beta_obs).abs() + (g - gamma_obs).abs();
            if (d - d0).abs() > 1e-9 {
                println!(
                    "      {:12} β {:.1}° γ {:.1}°  距離 {:.3} — {}",
                    vnames[v],
                    b.to_degrees(),
                    g.to_degrees(),
                    d,
                    if d < d0 { "改善" } else { "悪化" }
                );
            }
        }
    }

    // ---- [5] 総合判定 ----
    println!("\n[5] 総合判定 (生存 Δ≥+1 窓平均 + S2 [S3 は上表で個別評価]):");
    for (name, _ds, mean, s2) in &fam_rows {
        let verdict = if *mean >= 1.0 && *s2 {
            "生存"
        } else if *mean <= -1.0 || !*s2 {
            "棄却"
        } else {
            "未決"
        };
        println!("    {:16}  Δ̄ = {:+.2}  {}", name, mean, verdict);
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v17.9".into())),
        (
            "windows".into(),
            Json::Arr(
                geoms
                    .iter()
                    .enumerate()
                    .map(|(wi, &(ny, s1, s2))| {
                        Json::Obj(vec![
                            ("ny".into(), Json::Int(ny as i64)),
                            ("s1".into(), Json::Int(s1 as i64)),
                            ("s2".into(), Json::Int(s2 as i64)),
                            ("qq".into(), Json::Num(qq_win[wi])),
                            (
                                "lnz".into(),
                                Json::Arr(zc[wi].iter().map(|&x| Json::Num(x)).collect()),
                            ),
                            (
                                "ut_map".into(),
                                Json::Arr(
                                    map_ang[wi]
                                        .iter()
                                        .map(|&(b, g)| Json::Arr(vec![Json::Num(b), Json::Num(g)]))
                                        .collect(),
                                ),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "vnames".into(),
            Json::Arr(vnames.iter().map(|s| Json::Str((*s).into())).collect()),
        ),
    ]);
    let p = write_artifact("results/v179_measurejudge.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 判別は [2]–[5] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
