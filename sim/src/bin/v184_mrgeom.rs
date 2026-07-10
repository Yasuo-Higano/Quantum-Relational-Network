//! v18.4 M_R の幾何化 — B−L 破れプロファイルは大混合を作るか
//!
//! v18.3 の情報量ある外れ (最小シーソー M∝I は θ23 大混合に敗北) の指差しを実装:
//! B−L 破れ場 φ_B (周期 Gaussian, 位置 (bx,by) 6×6 格子 × 幅 σ_B 4 種 = 144 種の
//! 新離散データ) との非共役 Majorana 重なり M_R[i][j] = M₀·∫ψ_νi ψ_νj φ_B を作り、
//! シーソー m_ν = Y_ν M_R⁻¹ Y_νᵀ の prior-predictive 帯を v18.3 と同一プロトコル
//! (e 部事後 × σ_H 事後 × [ν, B] 一様 — データ盲目) で測る。
//! M₀ は角度・比から脱落。特異 M_R (|det| < 1e-12‖M‖³) は skip し計数。
//! 事前登録: (a) sin²θ23 帯が 0.55 被覆 = M_R 幾何は生きた仮説 (証拠昇格の道) /
//!           (b) 依然外す = 不足はテクスチャ超え (別 ν 幾何/Weinberg 型を記録)。
//! ゲート: σ 重みの v18.3 回帰 (±0.02)・厳密縮退。σ_H 4 本を thread 並列。

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

/// 複素対称シーソー用の道具 -------------------------------------------------

/// B−L 破れ場の周期 Gaussian プロファイル (位置 (bx, by), 幅 σ_B)
fn phi_b(nx: usize, ny: usize, bx: usize, by: usize, sigb: f64) -> Vec<f64> {
    let mut v = vec![0.0f64; nx * ny];
    for y in 0..ny {
        for x in 0..nx {
            let dx = {
                let d = (x as f64 - bx as f64).abs();
                d.min(nx as f64 - d)
            };
            let dy = {
                let d = (y as f64 - by as f64).abs();
                d.min(ny as f64 - d)
            };
            v[x + y * nx] = (-(dx * dx + dy * dy) / (2.0 * sigb * sigb)).exp();
        }
    }
    v
}

/// Majorana 重なり (非共役積): M[i][j] = Σ_x ψ_i(x) ψ_j(x) φ_B(x)
fn majorana_rect(la: &[Mode], lb: &[Mode], phib: &[f64]) -> M3 {
    let mut m = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..Q {
        for j in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for (x, &p) in phib.iter().enumerate() {
                let (ar, ai) = la[i][x];
                let (br, bi) = lb[j][x];
                sr += (ar * br - ai * bi) * p;
                si += (ar * bi + ai * br) * p;
            }
            m[i][j] = (sr, si);
        }
    }
    m
}

/// 3×3 複素逆行列 (余因子/det)。相対条件 |det| < 1e-12·‖m‖³ で None (特異扱い)
fn inv3c(m: &M3) -> Option<M3> {
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let sub = |a: (f64, f64), b: (f64, f64)| (a.0 - b.0, a.1 - b.1);
    let c = |i: usize, j: usize| m[i][j];
    let cof = |i0: usize, i1: usize, j0: usize, j1: usize| {
        sub(mul(c(i0, j0), c(i1, j1)), mul(c(i0, j1), c(i1, j0)))
    };
    let a00 = cof(1, 2, 1, 2);
    let a01 = cof(1, 2, 0, 2);
    let a02 = cof(1, 2, 0, 1);
    let det = {
        let t0 = mul(c(0, 0), a00);
        let t1 = mul(c(0, 1), a01);
        let t2 = mul(c(0, 2), a02);
        (t0.0 - t1.0 + t2.0, t0.1 - t1.1 + t2.1)
    };
    let nrm: f64 = m
        .iter()
        .flatten()
        .map(|&(a, b)| a * a + b * b)
        .sum::<f64>()
        .sqrt();
    let dabs = (det.0 * det.0 + det.1 * det.1).sqrt();
    if dabs < 1e-12 * nrm * nrm * nrm {
        return None;
    }
    let dinv = (det.0 / (dabs * dabs), -det.1 / (dabs * dabs));
    // 余因子行列 (転置込み = adjugate)
    let mut adj = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (r0, r1) = match j {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            let (c0, c1) = match i {
                0 => (1, 2),
                1 => (0, 2),
                _ => (0, 1),
            };
            let v = sub(mul(c(r0, c0), c(r1, c1)), mul(c(r0, c1), c(r1, c0)));
            let sgn = if (i + j) % 2 == 0 { 1.0 } else { -1.0 };
            adj[i][j] = (sgn * v.0, sgn * v.1);
        }
    }
    let mut inv = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            inv[i][j] = mul(adj[i][j], dinv);
        }
    }
    Some(inv)
}

/// m_ν = Y M⁻¹ Yᵀ (複素, 転置 — 共役なし)
fn seesaw_full(y: &M3, minv: &M3) -> M3 {
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let mut t = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let p = mul(y[i][k], minv[k][j]);
                re += p.0;
                im += p.1;
            }
            t[i][j] = (re, im);
        }
    }
    let mut m = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let p = mul(t[i][k], y[j][k]); // Yᵀ: (k,j) 成分 = y[j][k]
                re += p.0;
                im += p.1;
            }
            m[i][j] = (re, im);
        }
    }
    m
}

/// m m† のエルミート固有分解 (質量² 昇順, 左 U)
fn takagi_like(m: &M3) -> ([f64; 3], M3) {
    let mut hre = [[0.0f64; 3]; 3];
    let mut him = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = m[i][k];
                let (c, d) = m[j][k];
                re += a * c + b * d;
                im += b * c - a * d;
            }
            hre[i][j] = re;
            him[i][j] = im;
        }
    }
    eig_herm3(&hre, &him)
}

/// PMNS = U_e† U_ν → (sin²θ12, sin²θ23, sin²θ13, J_lep)
fn pmns_angles(ue: &M3, unu: &M3) -> (f64, f64, f64, f64) {
    let mut u = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = ue[k][i];
                let (c, d) = unu[k][j];
                re += a * c + b * d;
                im += a * d - b * c;
            }
            u[i][j] = (re, im);
        }
    }
    let a2 = |z: (f64, f64)| z.0 * z.0 + z.1 * z.1;
    let s13sq = a2(u[0][2]);
    let d12 = a2(u[0][0]) + a2(u[0][1]);
    let s12sq = if d12 > 0.0 { a2(u[0][1]) / d12 } else { 0.0 };
    let d23 = a2(u[1][2]) + a2(u[2][2]);
    let s23sq = if d23 > 0.0 { a2(u[1][2]) / d23 } else { 0.0 };
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let conj = |a: (f64, f64)| (a.0, -a.1);
    let jl = mul(mul(u[0][1], u[1][2]), mul(conj(u[0][2]), conj(u[1][1]))).1;
    (s12sq, s23sq, s13sq, jl)
}

/// 重みつき線形ヒストグラム
#[derive(Clone)]
struct WH {
    lo: f64,
    hi: f64,
    bins: Vec<f64>,
    wsum: f64,
}
impl WH {
    fn new(lo: f64, hi: f64, nb: usize) -> Self {
        WH {
            lo,
            hi,
            bins: vec![0.0; nb],
            wsum: 0.0,
        }
    }
    fn add(&mut self, x: f64, w: f64) {
        let nb = self.bins.len();
        let xx = x.clamp(self.lo, self.hi - 1e-12);
        let b = ((xx - self.lo) / (self.hi - self.lo) * nb as f64) as usize;
        self.bins[b.min(nb - 1)] += w;
        self.wsum += w;
    }
    fn merge(&mut self, o: &WH) {
        for (a, b) in self.bins.iter_mut().zip(o.bins.iter()) {
            *a += b;
        }
        self.wsum += o.wsum;
    }
    fn quant(&self, q: f64) -> f64 {
        let mut acc = 0.0;
        for (i, h) in self.bins.iter().enumerate() {
            acc += h;
            if acc >= q * self.wsum {
                let nb = self.bins.len() as f64;
                return self.lo + (i as f64 + 0.5) / nb * (self.hi - self.lo);
            }
        }
        self.hi
    }
}

fn main() {
    self_test();
    println!("=== v18.4 M_R の幾何化: B−L 破れプロファイルは大混合を作るか ===\n");
    println!("設計 (データ盲目): M_R[i][j] = M₀·∫ψ_νi ψ_νj φ_B (非共役 Majorana 重なり)。");
    println!("  φ_B = 周期 Gaussian (位置 (bx,by) ∈ 6×6 格子, 幅 σ_B ∈ {{1.5,3,6,12}}) — 新離散データ、一様。");
    println!("  m_ν = Y_ν M_R⁻¹ Y_νᵀ (M₀ は角度・比から脱落)。特異 M_R は skip し計数。");
    println!(
        "事前登録: (a) sin²θ23 の 68% 帯が 0.55 を被覆 → M_R 幾何は生きた仮説 (証拠昇格の道) /"
    );
    println!("          (b) 依然外す → 不足はテクスチャ超え (別 ν 幾何/Weinberg 型を記録)\n");
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
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let nc = 36usize;

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
    let spread = modes_k.iter().map(|r| r.2).fold(0.0f64, f64::max);
    check("厳密 3 重縮退", spread < 1e-8, format!("幅 {:.1e}", spread));
    let locs: std::sync::Arc<Vec<Vec<Mode>>> = std::sync::Arc::new(
        modes_k
            .iter()
            .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
            .collect(),
    );

    // B プロファイル集合 (位置 6×6 × 幅 4 = 144)
    let bxs: Vec<usize> = (0..6).map(|i| i * 6).collect(); // 0,6,..,30 (nx=36)
    let bys: Vec<usize> = (0..6).map(|i| i * 3).collect(); // 0,3,..,15 (ny=18)
    let sigbs = [1.5f64, 3.0, 6.0, 12.0];
    let mut phib_list: Vec<Vec<f64>> = Vec::new();
    for &sb in &sigbs {
        for &bx in &bxs {
            for &by in &bys {
                phib_list.push(phi_b(nx, ny, bx, by, sb));
            }
        }
    }
    let phibs: std::sync::Arc<Vec<Vec<f64>>> = std::sync::Arc::new(phib_list);
    println!(
        "[0] B プロファイル {} 種 (位置 {}×{} × 幅 {:?})",
        phibs.len(),
        bxs.len(),
        bys.len(),
        sigbs
    );

    // σ_H ごとに並列評価
    let g0: [f64; 4] = [2.0, 3.0, 4.0, 5.0];
    let scale = ((ny as f64) / (nx as f64)).sqrt();
    let t0 = std::time::Instant::now();
    let handles: Vec<_> = g0
        .iter()
        .enumerate()
        .map(|(isg, &s0)| {
            let locs = locs.clone();
            let phibs = phibs.clone();
            let tgt = tgt.clone();
            std::thread::spawn(move || {
                let sh = s0 * scale;
                let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
                    -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma)
                        + 2.0 * norm1
                };
                let ytab: Vec<M3> = (0..NK12 * NK12)
                    .map(|ab| yukawa_rect(nx, ny, &locs[ab % NK12], &locs[ab / NK12], sh))
                    .collect();
                let pair_y = |a: usize, b: usize, sf: usize, sg: usize| -> M3 {
                    let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
                    let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
                    had_prod_perm(&ytab[a1 + b1 * NK12], &ytab[a2 + b2 * NK12], sf, sg)
                };
                // Majorana 表: mtab[phib][a + a'*12] — ただし M_R は a=a' (同一住所) しか使わないので対角のみ
                // mtabd[ip][a] = majorana(locs[a], locs[a], phib_ip)
                let mtabd: Vec<Vec<M3>> = phibs
                    .iter()
                    .map(|pb| {
                        (0..NK12)
                            .map(|a| majorana_rect(&locs[a], &locs[a], pb))
                            .collect()
                    })
                    .collect();
                // σ_H 事後重み (クォーク五重和 — v18.3 と同一)
                let pair: Vec<([f64; 2], M3)> = (0..nc * nc * 6)
                    .map(|m| mass_and_vecs(&pair_y(m % nc, (m / nc) % nc, 0, m / (nc * nc))))
                    .collect();
                let mut le_all = Vec::with_capacity(nc * nc * 36);
                for sl in 0..6 {
                    for se_ in 0..6 {
                        for ab in 0..nc * nc {
                            le_all.push(ll2(
                                &mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_)),
                                tgt[4],
                                tgt[5],
                            ));
                        }
                    }
                }
                let lnze = lse(&le_all);
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
                                    let v = ckm_full(vu, &pair[md].1);
                                    let _ = vd;
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
                let sig_w = acc10.0 + acc10.1.ln() + lnze;
                // e 事後の最大 l_e (カット用の 1 パス目)
                let mut le_max = f64::NEG_INFINITY;
                for &l in &le_all {
                    if l > le_max {
                        le_max = l;
                    }
                }
                // 主ループ: e 配位 (カット le > le_max − 12) × ν (kν, σν) × B (144)
                let mut h12 = WH::new(0.0, 1.0, 200);
                let mut h23 = WH::new(0.0, 1.0, 200);
                let mut h13 = WH::new(0.0, 1.0, 400);
                let mut hr = WH::new(-14.0, 0.0, 280);
                let mut jl_pos = 0.0f64;
                let mut jl_tot = 0.0f64;
                let mut n_sing = 0u64;
                let mut n_eval = 0u64;
                for sl in 0..6 {
                    for se_ in 0..6 {
                        for ab in 0..nc * nc {
                            let idx = (sl * 6 + se_) * nc * nc + ab;
                            let le = le_all[idx];
                            if le < le_max - 12.0 {
                                continue;
                            }
                            let ye = pair_y(ab % nc, ab / nc, sl, se_);
                            let (_re, ue) = mass_and_vecs(&ye);
                            let we = (le - le_max).exp();
                            let kl = ab % nc;
                            for snu in 0..6 {
                                for knu in 0..nc {
                                    let ynu = pair_y(kl, knu, sl, snu);
                                    // M_R: ν の住所 knu = (k1, k2), 対角 Majorana 表から
                                    let (a1, a2) = (2 * (knu % 6), 2 * (knu / 6));
                                    for (ip, _) in phibs.iter().enumerate() {
                                        let mr =
                                            had_prod_perm(&mtabd[ip][a1], &mtabd[ip][a2], snu, snu);
                                        let minv = match inv3c(&mr) {
                                            Some(m) => m,
                                            None => {
                                                n_sing += 1;
                                                continue;
                                            }
                                        };
                                        let mnu = seesaw_full(&ynu, &minv);
                                        let (m2, unu) = takagi_like(&mnu);
                                        let (s12, s23, s13, jl) = pmns_angles(&ue, &unu);
                                        let r = ((m2[1] - m2[0]) / (m2[2] - m2[0]).max(1e-300))
                                            .max(1e-300);
                                        h12.add(s12, we);
                                        h23.add(s23, we);
                                        h13.add(s13, we);
                                        hr.add(r.ln().clamp(-14.0, -1e-9), we);
                                        jl_tot += we;
                                        if jl > 0.0 {
                                            jl_pos += we;
                                        }
                                        n_eval += 1;
                                    }
                                }
                            }
                        }
                    }
                }
                (
                    isg, sig_w, h12, h23, h13, hr, jl_pos, jl_tot, n_sing, n_eval,
                )
            })
        })
        .collect();
    let mut rs: Vec<_> = handles.into_iter().map(|h| h.join().unwrap()).collect();
    rs.sort_by_key(|r| r.0);
    let sig_w: Vec<f64> = rs.iter().map(|r| r.1).collect();
    println!(
        "[1] σ_H 事後重み (G0, 10 量): {:?} ({} s)",
        sig_w
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Vec<_>>(),
        t0.elapsed().as_secs()
    );
    // v18.3 の回帰 (同一機構)
    const REF_SIGW: [f64; 4] = [5.76, 7.77, 6.02, 0.22];
    let max_dev = sig_w
        .iter()
        .zip(REF_SIGW.iter())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0f64, f64::max);
    check(
        "σ 重みの v18.3 回帰 (±0.02)",
        max_dev < 0.02,
        format!("最大偏差 {:.3}", max_dev),
    );

    // σ 事後重みで各スレッドのヒストを合成
    let smax = sig_w.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut h12 = WH::new(0.0, 1.0, 200);
    let mut h23 = WH::new(0.0, 1.0, 200);
    let mut h13 = WH::new(0.0, 1.0, 400);
    let mut hr = WH::new(-14.0, 0.0, 280);
    let mut jl_pos = 0.0f64;
    let mut jl_tot = 0.0f64;
    let mut n_sing = 0u64;
    let mut n_eval = 0u64;
    for (isg, sw, a, b, c, d, jp, jt, ns, ne) in rs.iter() {
        let _ = isg;
        let w = (sw - smax).exp();
        let mut aa = a.clone();
        for x in aa.bins.iter_mut() {
            *x *= w;
        }
        aa.wsum *= w;
        h12.merge(&aa);
        let mut bb = b.clone();
        for x in bb.bins.iter_mut() {
            *x *= w;
        }
        bb.wsum *= w;
        h23.merge(&bb);
        let mut cc = c.clone();
        for x in cc.bins.iter_mut() {
            *x *= w;
        }
        cc.wsum *= w;
        h13.merge(&cc);
        let mut dd = d.clone();
        for x in dd.bins.iter_mut() {
            *x *= w;
        }
        dd.wsum *= w;
        hr.merge(&dd);
        jl_pos += jp * w;
        jl_tot += jt * w;
        n_sing += ns;
        n_eval += ne;
    }
    println!(
        "[2] 評価 {} 配位 (特異 M_R skip {} = {:.2}%) ({} s)",
        n_eval,
        n_sing,
        100.0 * n_sing as f64 / (n_eval + n_sing).max(1) as f64,
        t0.elapsed().as_secs()
    );

    // ---- [3] 帯と採点 ----
    println!("\n[3] PMNS 帯 (幾何 M_R, 16/50/84%) と v18.3 (M∝I) との比較:");
    let rows = [
        ("sin²θ12", 0.307f64, &h12, [0.003, 0.013, 0.297]),
        ("sin²θ23", 0.55, &h23, [0.003, 0.003, 0.242]),
        ("sin²θ13", 0.022, &h13, [0.001, 0.036, 0.306]),
    ];
    let mut cov23 = false;
    for (name, ob, h, old) in &rows {
        let (a, m, b) = (h.quant(0.16), h.quant(0.5), h.quant(0.84));
        let cov = a <= *ob && *ob <= b;
        if *name == "sin²θ23" {
            cov23 = cov;
        }
        println!(
            "    {:9} [{:.3}, {:.3}, {:.3}]  測定 {:.3}  {}   (v18.3: [{:.3}, {:.3}, {:.3}])",
            name,
            a,
            m,
            b,
            ob,
            if cov { "✓" } else { "✗" },
            old[0],
            old[1],
            old[2]
        );
    }
    {
        let (a, m, b) = (
            hr.quant(0.16).exp(),
            hr.quant(0.5).exp(),
            hr.quant(0.84).exp(),
        );
        let ob = 0.030f64;
        let cov = a <= ob && ob <= b;
        println!(
            "    {:9} [{:.4}, {:.4}, {:.4}]  測定 {:.3}  {}   (v18.3: [~0, ~0, 0.0023])",
            "r=Δm²比",
            a,
            m,
            b,
            ob,
            if cov { "✓" } else { "✗" }
        );
    }
    println!(
        "    sign(J_lep) > 0 の確率: {:.2}",
        jl_pos / jl_tot.max(1e-300)
    );
    println!(
        "\n    => 事前登録: {}",
        if cov23 {
            "(a) θ23 の帯が 0.55 を被覆 — B−L 破れの幾何は大混合を作れる (M_R 幾何は生きた仮説)"
        } else {
            "(b) θ23 は依然外れる — 不足は M_R テクスチャを超える構造 (別 ν 幾何/Weinberg 型を記録)"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v18.4".into())),
        (
            "bands".into(),
            Json::Arr(
                vec![
                    ("s12", &h12, 0.307f64),
                    ("s23", &h23, 0.55),
                    ("s13", &h13, 0.022),
                ]
                .into_iter()
                .map(|(n, h, ob)| {
                    Json::Obj(vec![
                        ("q".into(), Json::Str(n.into())),
                        ("q16".into(), Json::Num(h.quant(0.16))),
                        ("q50".into(), Json::Num(h.quant(0.5))),
                        ("q84".into(), Json::Num(h.quant(0.84))),
                        ("obs".into(), Json::Num(ob)),
                    ])
                })
                .collect(),
            ),
        ),
        ("r_q16".into(), Json::Num(hr.quant(0.16).exp())),
        ("r_q50".into(), Json::Num(hr.quant(0.5).exp())),
        ("r_q84".into(), Json::Num(hr.quant(0.84).exp())),
        (
            "p_jlep_positive".into(),
            Json::Num(jl_pos / jl_tot.max(1e-300)),
        ),
        ("n_singular".into(), Json::Int(n_sing as i64)),
        ("n_eval".into(), Json::Int(n_eval as i64)),
    ]);
    let p = write_artifact("results/v184_mrgeom.json", &j.render());
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
