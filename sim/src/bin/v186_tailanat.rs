//! v18.6 θ23 尾の解剖 — 5.9% の正体と全量 hit 亜集団
//!
//! v18.5 で P(sin²θ23 ∈ [0.4, 0.7]) = 5.9% の尾が残った。本バイナリは尾の解剖を
//! 行う (σ_H は事後峰 σ₀=3 の 1 本に固定 — 解剖であって帯の再測定ではない):
//!   [A] 尾の特徴分布: σ_B / 巻き |w| / σν — どこに集中するか (全体分布との比)
//!   [B] M_R 構造: 非対角比 ‖off‖/‖diag‖ — 尾 vs 全体
//!   [C] 混合の帰属: 尾の s23 は U_ν 由来か U_e 由来か (各因子単独の s23 平均)
//!   [D] 全量 hit 亜集団: 尾のうち θ12/θ13/r も同時に測定圏に入る重み割合
//! 事前登録: (a) 尾が特定構造に集中 (どこかの特徴で全体比 ≥3 倍) → 射影設計の入力 /
//!           (b) 無構造 → 尾は逃げ道でない (局在幾何へ直行)。
//! ゲート: 厳密縮退・σ₀=3 の σ 重みが v18.5 の 7.77 と一致 (±0.02)。

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

/// B−L 破れ場: 周期 Gaussian × 巻き位相 e^{2πi(wx·x/nx + wy·y/ny)} (複素)
fn phi_b_c(
    nx: usize,
    ny: usize,
    bx: usize,
    by: usize,
    sigb: f64,
    wx: usize,
    wy: usize,
) -> Vec<(f64, f64)> {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut v = vec![(0.0f64, 0.0f64); nx * ny];
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
            let g = (-(dx * dx + dy * dy) / (2.0 * sigb * sigb)).exp();
            let th = two_pi * (wx as f64 * x as f64 / nx as f64 + wy as f64 * y as f64 / ny as f64);
            v[x + y * nx] = (g * th.cos(), g * th.sin());
        }
    }
    v
}

/// Majorana 重なり (非共役積, 複素プロファイル): M[i][j] = Σ_x ψ_i ψ_j φ_B
fn majorana_rect(la: &[Mode], lb: &[Mode], phib: &[(f64, f64)]) -> M3 {
    let mut m = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..Q {
        for j in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for (x, &(pr, pi)) in phib.iter().enumerate() {
                let (ar, ai) = la[i][x];
                let (br, bi) = lb[j][x];
                let (cr, ci) = (ar * br - ai * bi, ar * bi + ai * br);
                sr += cr * pr - ci * pi;
                si += cr * pi + ci * pr;
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

fn main() {
    self_test();
    println!("=== v18.6 θ23 尾の解剖: 5.9% の正体 (σ_H = 事後峰 1 本) ===\n");
    println!("事前登録: (a) 尾が特定構造に集中 (特徴比 ≥3 倍) = 射影設計の入力 / (b) 無構造 = 尾は逃げ道でない\n");
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
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;

    // モード表
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
    let locs: Vec<Vec<Mode>> = modes_k
        .iter()
        .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
        .collect();

    // B プロファイル (v18.5 と同一構成 — 特徴インデックスを保持)
    let bxs = [0usize, 12, 24];
    let bys = [0usize, 9];
    let sigbs = [2.0f64, 4.0, 8.0];
    let mut phib_list: Vec<Vec<(f64, f64)>> = Vec::new();
    let mut bmeta: Vec<(usize, usize, usize)> = Vec::new(); // (幅 idx, 巻き和 |w|, 位置 idx)
    for (isb, &sb) in sigbs.iter().enumerate() {
        for (ibx, &bx) in bxs.iter().enumerate() {
            for (iby, &by) in bys.iter().enumerate() {
                for wx in 0..4usize {
                    for wy in 0..4usize {
                        phib_list.push(phi_b_c(nx, ny, bx, by, sb, wx, wy));
                        bmeta.push((isb, wx + wy, ibx * 2 + iby));
                    }
                }
            }
        }
    }
    println!(
        "[0] B プロファイル {} 種 (v18.5 と同一構成)",
        phib_list.len()
    );

    // σ_H = 事後峰 (σ₀ = 3) の 1 本
    let scale = ((ny as f64) / (nx as f64)).sqrt();
    let sh = 3.0 * scale;
    let t0 = std::time::Instant::now();
    let ytab: Vec<M3> = (0..NK12 * NK12)
        .map(|ab| yukawa_rect(nx, ny, &locs[ab % NK12], &locs[ab / NK12], sh))
        .collect();
    let pair_y = |a: usize, b: usize, sf: usize, sg: usize| -> M3 {
        let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
        let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
        had_prod_perm(&ytab[a1 + b1 * NK12], &ytab[a2 + b2 * NK12], sf, sg)
    };
    let mtabd: Vec<Vec<M3>> = phib_list
        .iter()
        .map(|pb| {
            (0..NK12)
                .map(|a| majorana_rect(&locs[a], &locs[a], pb))
                .collect()
        })
        .collect();
    // σ 重み回帰 (クォーク五重和)
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
                        let (rd, _vd) = &pair[md];
                        let lld = ll2(rd, tgt[2], tgt[3]);
                        let v = ckm_full(vu, &pair[md].1);
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
    check(
        "σ₀=3 の σ 重みが v18.5 と一致 (±0.02)",
        (sig_w - 7.77).abs() < 0.02,
        format!("{:.3} vs 7.77 ({} s)", sig_w, t0.elapsed().as_secs()),
    );

    // ---- 解剖ループ ----
    let mut le_max = f64::NEG_INFINITY;
    for &l in &le_all {
        if l > le_max {
            le_max = l;
        }
    }
    // 集計器: 尾 vs 全体の特徴重み
    let nb = phib_list.len();
    let mut w_all = 0.0f64;
    let mut w_tail = 0.0f64;
    let mut f_sb = [[0.0f64; 3]; 2]; // [all/tail][幅]
    let mut f_w = [[0.0f64; 7]; 2]; // 巻き和 0..6
    let mut f_snu = [[0.0f64; 6]; 2];
    let mut od_all = (0.0f64, 0.0f64); // (Σ w·ratio, Σ w)
    let mut od_tail = (0.0f64, 0.0f64);
    let mut s23nu_all = (0.0f64, 0.0f64);
    let mut s23nu_tail = (0.0f64, 0.0f64);
    let mut s23e_all = (0.0f64, 0.0f64);
    let mut s23e_tail = (0.0f64, 0.0f64);
    let mut w_full = 0.0f64; // 全量 hit (θ12/θ13/r も測定圏)
    let mut best = (f64::NEG_INFINITY, 0usize, 0usize, 0usize, [0.0f64; 4]); // 展示用
    let a2 = |z: (f64, f64)| z.0 * z.0 + z.1 * z.1;
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
                // U_e 単独の s23 (PMNS = U_e† の (1,2)/(2,2) 構造)
                let s23e = {
                    let u12 = a2(ue[2][1]); // (U_e†)[1][2] = conj(ue[2][1])
                    let u22 = a2(ue[2][2]);
                    if u12 + u22 > 0.0 {
                        u12 / (u12 + u22)
                    } else {
                        0.0
                    }
                };
                let kl = ab % nc;
                for snu in 0..6 {
                    for knu in 0..nc {
                        let ynu = pair_y(kl, knu, sl, snu);
                        let (a1, a2i) = (2 * (knu % 6), 2 * (knu / 6));
                        for ip in 0..nb {
                            let mr = had_prod_perm(&mtabd[ip][a1], &mtabd[ip][a2i], snu, snu);
                            let minv = match inv3c(&mr) {
                                Some(m) => m,
                                None => continue,
                            };
                            // M_R 非対角比
                            let mut dg = 0.0f64;
                            let mut od = 0.0f64;
                            for i in 0..3 {
                                for j in 0..3 {
                                    let n = a2(mr[i][j]);
                                    if i == j {
                                        dg += n;
                                    } else {
                                        od += n;
                                    }
                                }
                            }
                            let odr = (od / dg.max(1e-300)).sqrt();
                            let mnu = seesaw_full(&ynu, &minv);
                            let (m2, unu) = takagi_like(&mnu);
                            let (s12, s23, s13, _jl) = pmns_angles(&ue, &unu);
                            let r = ((m2[1] - m2[0]) / (m2[2] - m2[0]).max(1e-300)).max(1e-300);
                            let s23nu = {
                                let u12 = a2(unu[1][2]);
                                let u22 = a2(unu[2][2]);
                                if u12 + u22 > 0.0 {
                                    u12 / (u12 + u22)
                                } else {
                                    0.0
                                }
                            };
                            let (isb, wsum, _ipos) = bmeta[ip];
                            w_all += we;
                            f_sb[0][isb] += we;
                            f_w[0][wsum] += we;
                            f_snu[0][snu] += we;
                            od_all.0 += we * odr;
                            od_all.1 += we;
                            s23nu_all.0 += we * s23nu;
                            s23nu_all.1 += we;
                            s23e_all.0 += we * s23e;
                            s23e_all.1 += we;
                            let in_tail = (0.4..=0.7).contains(&s23);
                            if in_tail {
                                w_tail += we;
                                f_sb[1][isb] += we;
                                f_w[1][wsum] += we;
                                f_snu[1][snu] += we;
                                od_tail.0 += we * odr;
                                od_tail.1 += we;
                                s23nu_tail.0 += we * s23nu;
                                s23nu_tail.1 += we;
                                s23e_tail.0 += we * s23e;
                                s23e_tail.1 += we;
                                let full = (0.25..=0.37).contains(&s12)
                                    && (0.01..=0.04).contains(&s13)
                                    && (0.02..=0.045).contains(&r);
                                if full {
                                    w_full += we;
                                }
                                // 展示: 尾内の全量スコア (粗い積 — フィットではない)
                                let sc = -((s12 - 0.307f64).powi(2) / 0.01
                                    + (s23 - 0.55f64).powi(2) / 0.01
                                    + (s13 - 0.022f64).powi(2) / 0.001
                                    + ((r.ln() - (0.03f64).ln()) / 1.0).powi(2));
                                if sc > best.0 {
                                    best = (sc, ip, knu, snu, [s12, s23, s13, r]);
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    println!(
        "[1] 解剖完了: 尾の重み割合 {:.4} ({} s)",
        w_tail / w_all.max(1e-300),
        t0.elapsed().as_secs()
    );

    // ---- [A] 特徴分布 (尾/全体の比) ----
    println!("\n[A] 尾の特徴集中 (尾内シェア ÷ 全体シェア — ≥3 で「集中」):");
    let mut max_conc = 0.0f64;
    let mut max_name = String::new();
    let show =
        |name: &str, tail: f64, all: f64, wt: f64, wa: f64, mc: &mut f64, mn: &mut String| {
            let ts = tail / wt.max(1e-300);
            let as_ = all / wa.max(1e-300);
            let ratio = ts / as_.max(1e-300);
            println!(
                "    {:12} 尾シェア {:.3} / 全体シェア {:.3} → 比 {:.2}",
                name, ts, as_, ratio
            );
            if ratio > *mc {
                *mc = ratio;
                *mn = name.to_string();
            }
        };
    for i in 0..3 {
        show(
            &format!("σ_B={}", sigbs[i]),
            f_sb[1][i],
            f_sb[0][i],
            w_tail,
            w_all,
            &mut max_conc,
            &mut max_name,
        );
    }
    for i in 0..7 {
        if f_w[0][i] > 0.0 {
            show(
                &format!("|w|={}", i),
                f_w[1][i],
                f_w[0][i],
                w_tail,
                w_all,
                &mut max_conc,
                &mut max_name,
            );
        }
    }
    for i in 0..6 {
        show(
            &format!("σν={}", i),
            f_snu[1][i],
            f_snu[0][i],
            w_tail,
            w_all,
            &mut max_conc,
            &mut max_name,
        );
    }

    // ---- [B][C] 構造と帰属 ----
    println!(
        "\n[B] M_R 非対角比 (重み平均): 全体 {:.4} / 尾 {:.4}",
        od_all.0 / od_all.1.max(1e-300),
        od_tail.0 / od_tail.1.max(1e-300)
    );
    println!(
        "[C] 混合の帰属 (s23 の因子別平均): U_ν 側 全体 {:.3} → 尾 {:.3} / U_e 側 全体 {:.3} → 尾 {:.3}",
        s23nu_all.0 / s23nu_all.1.max(1e-300),
        s23nu_tail.0 / s23nu_tail.1.max(1e-300),
        s23e_all.0 / s23e_all.1.max(1e-300),
        s23e_tail.0 / s23e_tail.1.max(1e-300)
    );

    // ---- [D] 全量 hit 亜集団 ----
    println!(
        "\n[D] 全量 hit 亜集団: P(θ12∧θ13∧r 測定圏 | 尾) = {:.4} (尾重みの内数)",
        w_full / w_tail.max(1e-300)
    );
    let (isb, wsum, ipos) = bmeta[best.1];
    println!(
        "    展示 (尾内の最良 config — フィットではない): σ_B={}, |w|={}, 位置#{}, kν={}, σν={} → s12={:.3}, s23={:.3}, s13={:.3}, r={:.4}",
        sigbs[isb], wsum, ipos, best.2, best.3, best.4[0], best.4[1], best.4[2], best.4[3]
    );

    // ---- 判定 ----
    println!(
        "\n[判定] 最大集中: {} (比 {:.2}) → {}",
        max_name,
        max_conc,
        if max_conc >= 3.0 {
            "事前登録 (a): 尾は特定構造に集中 — 射影設計の入力とする"
        } else {
            "事前登録 (b): 尾は無構造 (比 <3) — 尾は逃げ道でない (局在幾何へ直行)"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v18.6".into())),
        ("p_tail".into(), Json::Num(w_tail / w_all.max(1e-300))),
        (
            "p_full_given_tail".into(),
            Json::Num(w_full / w_tail.max(1e-300)),
        ),
        ("odr_all".into(), Json::Num(od_all.0 / od_all.1.max(1e-300))),
        (
            "odr_tail".into(),
            Json::Num(od_tail.0 / od_tail.1.max(1e-300)),
        ),
        (
            "s23nu_tail".into(),
            Json::Num(s23nu_tail.0 / s23nu_tail.1.max(1e-300)),
        ),
        (
            "s23e_tail".into(),
            Json::Num(s23e_tail.0 / s23e_tail.1.max(1e-300)),
        ),
        ("max_conc".into(), Json::Num(max_conc)),
        (
            "max_conc_feature".into(),
            Json::Str(max_name.clone().into()),
        ),
    ]);
    let p = write_artifact("results/v186_tailanat.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 解剖は [A]–[D] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
