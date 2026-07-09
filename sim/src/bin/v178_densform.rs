//! v17.8 config-density の解析形 — bin 規約を消す (PROMPT/3 v18.2)
//!
//! v16.13 で生存した「測度補正クラス」の弱点は、ヒストグラム平坦化の bin 規約
//! (S1 の残留自由度) だった。本バイナリは配位密度 n(x₁,x₂) (x = 質量対数比) に
//! **モーメント構成の解析形**を試す: 2D Gaussian (パラメータ = 標本平均と共分散 —
//! フィット自由度ゼロ、データ非参照)。これが立てば平坦化は
//!     w(c) = +½ (x−μ)ᵀ Σ⁻¹ (x−μ)   (二次傾き)
//! という bin なしの閉形式になる。
//!
//! 検査 2 段 (事前登録):
//!   [分布] 深さ周辺分布の QQ-RMS < 0.6 で「Gauss 近似は分布レベルでも良い」
//!     (裾は非 Gauss でも測度としては使える — 分布と測度の判定を分ける)
//!   [測度] 3 窓 ((36;1,3), (36;3,3), (24;3,3)) × 二次傾き測度の Δ が、同窓の
//!     binned flat2D (bw=0.5) の Δ と ±0.3 nats で一致 → bin 規約は消去可能
//! 分岐: (a) 分布も測度も合う = 解析形確定 / (b) 分布は外れるが測度は合う =
//! 「二次傾きで十分」(実用的解析形) / (c) 測度が合わない = 解析形未達 (binned 継続)。
//!
//! 装置: 全窓キャッシュ命中前提 (分単位)。一様の lnZ 回帰 (G0, v16.9/v17.5)。
//! 副検査: v16.13 の局所勾配 −0.19 が Gauss 形の予言 −(d−μ_d)/σ_d² と整合するか。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// G0 回帰アンカー (v16.9 / v17.5)
const REF_G0: [(usize, usize, usize, f64); 3] = [
    (36, 1, 3, -21.982785),
    (36, 3, 3, -21.756581),
    (24, 3, 3, -19.641),
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

fn main() {
    self_test();
    println!("=== v17.8 config-density の解析形: モーメント Gauss / 二次傾き測度 (3 窓) ===\n");
    println!("事前登録: [分布] 深さ QQ-RMS < 0.6 / [測度] 二次傾き Δ が binned flat2D Δ と ±0.3");
    println!("分岐: (a) 両方合う = 解析形確定 / (b) 測度だけ合う = 二次傾きで十分 / (c) 測度不一致 = 未達\n");
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
    let geoms: [(usize, usize, usize); 3] = [(36, 1, 3), (36, 3, 3), (24, 3, 3)];

    // モード表 (全キャッシュ命中前提)
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
    let mut hits = 0;
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
        } else {
            hits += 1;
        }
        locs_map.insert(
            (ny, s),
            modes_k
                .iter()
                .map(|(m, _, _)| localize_stable_rect(nx, ny, m))
                .collect(),
        );
    }
    println!(
        "[0] モード表: キャッシュ命中 {}/{} 系列\n",
        hits,
        need.len()
    );

    // ---- 窓ごとの密度解析と 3 測度評価 ----
    println!("[1] 窓ごとの密度モーメント・QQ・3 測度 (一様 / binned flat2D / 二次傾き):");
    let g0 = [2.0f64, 3.0, 4.0, 5.0];
    let mut rows = Vec::new();
    for &(ny, s1, s2) in &geoms {
        let scale = ((ny as f64) / (nx as f64)).sqrt();
        let locs1 = &locs_map[&(ny, s1)];
        let locs2 = &locs_map[&(ny, s2)];
        let mut trow: Vec<[f64; 3]> = Vec::new();
        let mut qq_all: Vec<f64> = Vec::new();
        let mut slope_pred = 0.0f64;
        for &s0 in &g0 {
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
            for m in 0..nc * nc * 6 {
                let y = pair_y(m % nc, (m / nc) % nc, 0, m / (nc * nc));
                let (r, v) = mass_and_vecs(&y);
                pair_r.push(r);
                pair_v.push(v);
            }
            // 密度解析 (クォーク対集団)
            let (mu_d, var_d, _cov, wq_gauss) = gauss_weights(&pair_r);
            let mut depths: Vec<f64> = pair_r.iter().map(|r| -(r[0] + r[1])).collect();
            let qq = qq_rms_depth(&mut depths, -mu_d, var_d);
            qq_all.push(qq);
            // v16.13 の勾配整合の予言: 局所勾配 at 中心+1σ ≈ −1σ_d/σ_d² = −1/σ_d
            slope_pred = -1.0 / var_d.sqrt();
            let wq_bin = flatten_weights_2d(&pair_r, 0.5);
            // e セクター
            let mut er: Vec<[f64; 2]> = Vec::with_capacity(nc * nc * 36);
            for sl in 0..6 {
                for se_ in 0..6 {
                    for ab in 0..nc * nc {
                        er.push(mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_)));
                    }
                }
            }
            let (_, _, _, we_gauss) = gauss_weights(&er);
            let we_bin = flatten_weights_2d(&er, 0.5);
            let mut ze = [Acc::new(); 3];
            let mut ne = [Acc::new(); 3];
            for (i, r) in er.iter().enumerate() {
                let l = ll2(r, tgt[4], tgt[5]);
                let sps = [0.0, we_bin[i], we_gauss[i]];
                for v in 0..3 {
                    ze[v].add(sps[v] + l);
                    ne[v].add(sps[v]);
                }
            }
            // クォーク五重和
            let mut zq = [Acc::new(); 3];
            let mut nq = [Acc::new(); 3];
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
                                let sps =
                                    [0.0, wq_bin[mu] + wq_bin[md], wq_gauss[mu] + wq_gauss[md]];
                                for vv in 0..3 {
                                    zq[vv].add(sps[vv] + ll10);
                                    nq[vv].add(sps[vv]);
                                }
                            }
                        }
                    }
                }
            }
            let mut row = [0.0f64; 3];
            for v in 0..3 {
                row[v] = (zq[v].val() - nq[v].val()) + (ze[v].val() - ne[v].val());
            }
            trow.push(row);
        }
        let mut z = [0.0f64; 3];
        for v in 0..3 {
            let col: Vec<f64> = trow.iter().map(|r| r[v]).collect();
            z[v] = lse(&col) - (g0.len() as f64).ln();
        }
        let qq_max = qq_all.iter().cloned().fold(0.0f64, f64::max);
        println!(
            "    ({:2};{},{})  一様 {:8.3}  flat2D {:8.3} (Δ{:+.2})  gauss {:8.3} (Δ{:+.2})  QQ最大 {:.2}  局所勾配予言 {:+.2}",
            ny, s1, s2, z[0], z[1], z[1] - z[0], z[2], z[2] - z[0], qq_max, slope_pred
        );
        rows.push((ny, s1, s2, z[0], z[1], z[2], qq_max));
    }

    // ---- ゲート ----
    println!("\n[ゲート]");
    for (ri, &(ny, s1, s2, u, _, _, _)) in rows.iter().enumerate() {
        let (rny, rs1, rs2, refv) = REF_G0[ri];
        assert_eq!((rny, rs1, rs2), (ny, s1, s2));
        check(
            &format!("({};{},{}) 一様 G0 回帰 (±0.05)", ny, s1, s2),
            (u - refv).abs() < 0.05,
            format!("{:.3} vs {:.3}", u, refv),
        );
    }

    // ---- 判定 ----
    println!("\n[2] 事前登録判定:");
    let mut meas_ok = true;
    let mut dist_ok = true;
    for &(ny, s1, s2, u, fb, fg, qq) in &rows {
        let (db, dg) = (fb - u, fg - u);
        let m_ok = (dg - db).abs() <= 0.3;
        let d_ok = qq < 0.6;
        println!(
            "    ({:2};{},{}): |Δ_gauss − Δ_binned| = {:.2} ({}), QQ = {:.2} ({})",
            ny,
            s1,
            s2,
            (dg - db).abs(),
            if m_ok {
                "測度一致"
            } else {
                "測度不一致"
            },
            qq,
            if d_ok { "分布良" } else { "分布外れ" }
        );
        meas_ok &= m_ok;
        dist_ok &= d_ok;
    }
    if meas_ok && dist_ok {
        println!(
            "    => (a) 解析形確定: モーメント Gauss — 平坦化測度は bin なしの二次傾きに置換可能。"
        );
    } else if meas_ok {
        println!("    => (b) 分布は非 Gauss だが測度としては二次傾きで十分 — MSR-FLATTEN の bin 規約は消去可能");
        println!("       (measures.yml に MSR-GAUSSFLAT を実用的解析形として追記する資格)。");
    } else {
        println!("    => (c) 解析形未達 — 二次傾きは binned 平坦化を再現しない。密度の高次構造が測度に効く。");
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v17.8".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(ny, s1, s2, u, fb, fg, qq)| {
                        Json::Obj(vec![
                            ("ny".into(), Json::Int(ny as i64)),
                            ("s1".into(), Json::Int(s1 as i64)),
                            ("s2".into(), Json::Int(s2 as i64)),
                            ("lnz_uniform".into(), Json::Num(u)),
                            ("lnz_flat2d_binned".into(), Json::Num(fb)),
                            ("lnz_gauss_quad".into(), Json::Num(fg)),
                            ("qq_rms_max".into(), Json::Num(qq)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v178_densform.json", &j.render());
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
