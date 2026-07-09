//! v18.3 レプトン混合の holdout — 収束点の PMNS prior-predictive (PROMPT/3 v19)
//!
//! 設計はバイナリ冒頭の println! に固定 (データ盲目性の宣言込み)。
//! 要点: ν_R は独自の Wilson 住所 (k_ν, σ_ν) を持ち 10 量尤度に拘束されない。
//! M_R = M₀·I (直交性による零パラメータ選択) → m_ν ∝ Y_ν Y_νᵀ (複素対称,
//! Takagi 相当は m m† のエルミート固有分解で実装)。PMNS = U_e† U_ν。
//! k_L は Y_e/Y_ν で共有 (v16.1 の共有 kQ の教訓)。帯 = e 部事後 × σ_H 事後 ×
//! ν 一様 — 幾何は収束点 (36×18; 3,3) 固定。
//! 採点対象: sin²θ12/23/13, r = Δm²21/Δm²31, |J_lep| 帯 + sign(J_lep) 確率
//! (orientation は quark sign(J) で固定済み — レプトン CP 符号は真の予言)。

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

/// 複素対称シーソー m ∝ Y Yᵀ (転置 — 共役なし) の 3×3 積
fn seesaw_mnu(y: &M3) -> M3 {
    let mut m = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = y[i][k];
                let (c, d) = y[j][k];
                re += a * c - b * d;
                im += a * d + b * c;
            }
            m[i][j] = (re, im);
        }
    }
    m
}

/// m m† (エルミート) の固有分解から (質量² 昇順, 左ユニタリ U)
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

/// PMNS = U_e† U_ν と (sin²θ12, sin²θ23, sin²θ13, J_lep)
fn pmns_angles(ue: &M3, unu: &M3) -> (f64, f64, f64, f64) {
    let mut u = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                // (U_e† U_ν)_ij = Σ_k conj(U_e[k? ]) ... eig_herm3 の vecs[k][i] = 固有ベクトル k の成分 i
                let (a, b) = ue[k][i];
                let (c, d) = unu[k][j];
                // conj(ue) · unu
                re += a * c + b * d;
                im += a * d - b * c;
            }
            u[i][j] = (re, im);
        }
    }
    let a2 = |z: (f64, f64)| z.0 * z.0 + z.1 * z.1;
    let s13sq = a2(u[0][2]);
    let denom12 = a2(u[0][0]) + a2(u[0][1]);
    let s12sq = if denom12 > 0.0 {
        a2(u[0][1]) / denom12
    } else {
        0.0
    };
    let denom23 = a2(u[1][2]) + a2(u[2][2]);
    let s23sq = if denom23 > 0.0 {
        a2(u[1][2]) / denom23
    } else {
        0.0
    };
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let conj = |a: (f64, f64)| (a.0, -a.1);
    let jl = mul(mul(u[0][1], u[1][2]), mul(conj(u[0][2]), conj(u[1][1]))).1;
    (s12sq, s23sq, s13sq, jl)
}

/// 重みつきヒストグラム (線形域 [0,1] or 対数域)
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
    println!(
        "=== v18.3 レプトン混合の holdout — 収束点の PMNS prior-predictive (PROMPT/3 v19) ===\n"
    );
    println!(
        "設計 (データ盲目): ν_R は独自の Wilson 住所 (k_ν, σ_ν) — 10 量尤度に一切拘束されない。"
    );
    println!("  M_R = M₀·I (モード直交性による零パラメータ選択) → m_ν ∝ Y_ν Y_νᵀ — 角度・質量比から M₀ 脱落。");
    println!("  帯 = e 部事後 exp(l_e) × σ_H 事後 × (k_ν, σ_ν) 一様の prior-predictive。");
    println!("  k_L は Y_e と Y_ν で共有 (クォークの共有 kQ と同じ要件 — v16.1 の教訓)。");
    println!("測定 (PDG/global fit 級): sin²θ12=0.307, sin²θ23=0.55, sin²θ13=0.022, r=Δm²21/|Δm²31|=0.030\n");
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

    // モード表 (rect キャッシュ — 収束点)
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

    // σ_H 事後の重み: G0 の terms10 (v17.10 の一次ソースから転記でなく再計算 — クォーク五重和)
    let g0: [f64; 4] = [2.0, 3.0, 4.0, 5.0];
    let scale = ((ny as f64) / (nx as f64)).sqrt();
    let t0 = std::time::Instant::now();
    let mut sig_w = [0.0f64; 4];
    let mut ytabs: Vec<Vec<M3>> = Vec::new();
    for (isg, &s0) in g0.iter().enumerate() {
        let sh = s0 * scale;
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
                    le.push(ll2(
                        &mass_ratios(&pair_y(ab % nc, ab / nc, sl, se_)),
                        tgt[4],
                        tgt[5],
                    ));
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
        sig_w[isg] = acc10.0 + acc10.1.ln() + lnze;
        ytabs.push(ytab);
    }
    println!(
        "[1] σ_H 事後重み (G0, 10 量): {:?} ({} s)",
        sig_w
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Vec<_>>(),
        t0.elapsed().as_secs()
    );

    // ---- [2] PMNS prior-predictive ----
    // 重み: exp(sig_w[isg] − max) × exp(l_e − max_e) × 一様 (k_ν, σ_ν)
    let smax = sig_w.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let mut h12 = WH::new(0.0, 1.0, 200);
    let mut h23 = WH::new(0.0, 1.0, 200);
    let mut h13 = WH::new(0.0, 1.0, 400);
    let mut hr = WH::new(-14.0, 0.0, 280); // ln r (質量² 比)
    let mut hjl = WH::new(-30.0, 0.0, 300); // ln |J_lep|
    let mut jl_pos = 0.0f64;
    let mut jl_tot = 0.0f64;
    let t1 = std::time::Instant::now();
    for (isg, ytab) in ytabs.iter().enumerate() {
        let wsig = (sig_w[isg] - smax).exp();
        let pair_y = |a: usize, b: usize, sf: usize, sg: usize| -> M3 {
            let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
            let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
            had_prod_perm(&ytab[a1 + b1 * NK12], &ytab[a2 + b2 * NK12], sf, sg)
        };
        // e 部: (kL⊗ke = ab, sl, se) — U_e と l_e。ν 部: (kL 共有, kν, sν)。
        // ab は kL (a) と ke (b) を含む: a = ab % nc, b = ab / nc
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..nc * nc {
                    let ye = pair_y(ab % nc, ab / nc, sl, se_);
                    let (re_, ue) = mass_and_vecs(&ye);
                    let le = ll2(&re_, tgt[4], tgt[5]);
                    let we = wsig * (le - 0.0).exp();
                    if we < 1e-8 {
                        continue; // 事後質量が無視可能な e 配位は飛ばす (高速化 — 帯には影響しない)
                    }
                    let kl = ab % nc;
                    for snu in 0..6 {
                        for knu in 0..nc {
                            let ynu = pair_y(kl, knu, sl, snu);
                            let mnu = seesaw_mnu(&ynu);
                            let (m2, unu) = takagi_like(&mnu);
                            let (s12, s23, s13, jl) = pmns_angles(&ue, &unu);
                            let r = if m2[2] > 0.0 {
                                ((m2[1] - m2[0]) / (m2[2] - m2[0]).max(1e-300)).max(1e-300)
                            } else {
                                1e-300
                            };
                            let w = we; // (kν, σν) 一様
                            h12.add(s12, w);
                            h23.add(s23, w);
                            h13.add(s13, w);
                            hr.add(r.ln().clamp(-14.0, -1e-9), w);
                            hjl.add(jl.abs().max(1e-300).ln().clamp(-30.0, -1e-9), w);
                            jl_tot += w;
                            if jl > 0.0 {
                                jl_pos += w;
                            }
                        }
                    }
                }
            }
        }
    }
    println!("[2] prior-predictive 完了 ({} s)", t1.elapsed().as_secs());

    // ---- [3] 帯と採点 ----
    let obs = [
        ("sin²θ12", 0.307f64, &h12),
        ("sin²θ23", 0.55, &h23),
        ("sin²θ13", 0.022, &h13),
    ];
    println!("\n[3] PMNS holdout 帯 (16/50/84%) と採点:");
    let mut score = Vec::new();
    for (name, ob, h) in &obs {
        let (a, m, b) = (h.quant(0.16), h.quant(0.5), h.quant(0.84));
        let cov = a <= *ob && *ob <= b;
        println!(
            "    {:9} [{:.3}, {:.3}, {:.3}]  測定 {:.3}  68% 帯{}",
            name,
            a,
            m,
            b,
            ob,
            if cov {
                "は測定を含む ✓"
            } else {
                "は測定を外す ✗"
            }
        );
        score.push((name.to_string(), a, m, b, *ob, cov));
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
            "    {:9} [{:.4}, {:.4}, {:.4}]  測定 {:.3}  68% 帯{}",
            "r=Δm²比",
            a,
            m,
            b,
            ob,
            if cov {
                "は測定を含む ✓"
            } else {
                "は測定を外す ✗"
            }
        );
        score.push(("r".into(), a, m, b, ob, cov));
    }
    {
        let (a, m, b) = (
            hjl.quant(0.16).exp(),
            hjl.quant(0.5).exp(),
            hjl.quant(0.84).exp(),
        );
        println!(
            "    {:9} [{:.1e}, {:.1e}, {:.1e}]  (測定は |J_lep|max≈0.033·sinδ — δ 未確定につき帯のみ登録)",
            "|J_lep|", a, m, b
        );
        println!(
            "    sign(J_lep) > 0 の事後確率: {:.2} (orientation は quark sign(J) で固定済み — レプトン δ_CP の符号予言)",
            jl_pos / jl_tot.max(1e-300)
        );
        score.push(("|J_lep|".into(), a, m, b, f64::NAN, true));
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v18.3".into())),
        (
            "bands".into(),
            Json::Arr(
                score
                    .iter()
                    .map(|(n, a, m, b, ob, cov)| {
                        Json::Obj(vec![
                            ("q".into(), Json::Str(n.clone().into())),
                            ("q16".into(), Json::Num(*a)),
                            ("q50".into(), Json::Num(*m)),
                            ("q84".into(), Json::Num(*b)),
                            (
                                "obs".into(),
                                Json::Num(if ob.is_nan() { -1.0 } else { *ob }),
                            ),
                            ("covered".into(), Json::Bool(*cov)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "p_jlep_positive".into(),
            Json::Num(jl_pos / jl_tot.max(1e-300)),
        ),
    ]);
    let p = write_artifact("results/v183_pmns.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 帯と採点は [3] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
