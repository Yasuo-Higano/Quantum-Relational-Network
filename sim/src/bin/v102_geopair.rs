//! v10.2 対込みの幾何選択 — T³ は勝ち続けるか (v10.1 の波及の解決)
//!
//! v10.1 は対角世代対がデータに否定されること、そして T²+自由対 (lnZ₉ −19.86) が
//! T³+対角対 (−21.84) に勝つことを示した。v8.1/v9.1 の「T³ が勝つ」は対仮定に
//! 条件付きだった。本バイナリは {T¹, T², T³}×Z₆ の証拠を**対 marginalize 込み**で
//! 再計算し、幾何選択の結論を対仮定から解放する。
//!
//! 対空間: 場 F の対は各追加トーラス t≥2 への置換 σᵗ_F ∈ S₃。大域relabel τᵗ で
//! σᵗ_F → τᵗ∘σᵗ_F (全場一斉) なので σ_Q = (e,…) にゲージ固定し、
//! (σ_u, σ_d, σ_L, σ_e) ∈ (S₃^(nt−1))⁴ を一様事前で marginalize
//! (Occam 罰 = 4(nt−1)ln6: T² で 7.17, T³ で 14.33)。
//!
//! T³ の全 9 量は (K_Q,K_u,K_d,σ_u,σ_d) の五重和 ~10¹³ で直接和は不可能。
//! **証明付き打ち切り**を使う: (ku,σu) と (kd,σd) のリストを質量尤度で降順に並べ、
//! llu+lld ≥ max−Δ (Δ=35) の対だけ CKM を計算。残りの寄与は
//!   Σ_残り e^{llu+lld+llck} ≤ e^{3·norm1} · (S_u·S_d − Σ_計算済み e^{llu+lld})
//! で厳密に上から抑える (CKM 項は各観測で ≤ norm1)。lnZ は区間 [下端, 上端] で
//! 報告し、幅 < 0.01 nats を PASS 条件とする。
//!
//! 回帰検査 (装置検証):
//!   T¹ (対なし):        質量 −53.7694 / 全9 −63.5703 (v9.2 安定側)
//!   T² 対角:            質量 −18.7619 (v9.2) — 対角経路
//!   T³ 対角:            質量 −16.4918 / 全9 −21.8436 (v9.2) — 対角経路
//!   T² 自由対:          質量 −15.0839 / 全9 −19.8633 (v10.1 の厳密五重和)
//!     — 特に全9は「打ち切り実装 vs v10.1 の全和実装」の相互検証になる。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];
const DELTA_CUT: f64 = 35.0; // 打ち切り深さ (nats)

// 回帰目標 (一次ソース: results/v92_labelstab.json, results/v101_pairing.json)
const REF_T1_MASS: f64 = -53.76939196753885;
const REF_T1_NINE: f64 = -63.57030452988765;
const REF_T2_MASS_DIAG: f64 = -18.761895471783593;
const REF_T3_MASS_DIAG: f64 = -16.4918220802609;
const REF_T3_NINE_DIAG: f64 = -21.843569420441305;
const REF_T2_MASS_PERM: f64 = -15.083852846457935;
const REF_T2_NINE_PERM: f64 = -19.86334559888438;

type C3v = [(f64, f64); NS];
type M3 = [[(f64, f64); 3]; 3];

fn flux_modes(k_half: usize) -> (Vec<C3v>, f64, f64) {
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / NS as f64;
    let wl = phi * k_half as f64 / 2.0;
    let idx = |x: usize, y: usize| x + y * N;
    let m = 2 * NS;
    let mut a = vec![0.0; m * m];
    let addhop = |a: &mut Vec<f64>, i: usize, j: usize, th: f64| {
        let (c, s) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + NS) + (i + NS) * m] += -c;
        a[(i + NS) + (j + NS) * m] += -c;
        a[j + (i + NS) * m] += s;
        a[(j + NS) + i * m] += -s;
        a[i + (j + NS) * m] += -s;
        a[(i + NS) + j * m] += s;
    };
    for x in 0..N {
        for y in 0..N {
            addhop(&mut a, idx(x, y), idx(x, (y + 1) % N), phi * x as f64 + wl);
            let th = if x == N - 1 {
                -phi * (N as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % N, y), th);
        }
    }
    let (w, v) = jacobi_eigh(&a, m);
    let gap = w[2 * Q] - w[2 * Q - 1];
    let spread = w[2 * Q - 1] - w[0];
    let mut modes: Vec<C3v> = Vec::new();
    for kk in 0..2 * Q {
        let mut psi = [(0.0f64, 0.0f64); NS];
        for i in 0..NS {
            psi[i] = (v[i + kk * m], v[(i + NS) + kk * m]);
        }
        for pm in &modes {
            let (mut pr, mut pi) = (0.0, 0.0);
            for i in 0..NS {
                pr += pm[i].0 * psi[i].0 + pm[i].1 * psi[i].1;
                pi += pm[i].0 * psi[i].1 - pm[i].1 * psi[i].0;
            }
            for i in 0..NS {
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
    let n = 3;
    let m = 6;
    let mut emb = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            emb[i + j * m] = hre[i][j];
            emb[i + (j + n) * m] = -him[i][j];
            emb[(i + n) + j * m] = him[i][j];
            emb[(i + n) + (j + n) * m] = hre[i][j];
        }
    }
    let (w, v) = jacobi_eigh(&emb, m);
    let mut lam = [0.0f64; 3];
    let mut vecs = [[(0.0f64, 0.0f64); 3]; 3];
    for k in 0..3 {
        lam[k] = 0.5 * (w[2 * k] + w[2 * k + 1]);
        for i in 0..3 {
            vecs[k][i] = (v[i + (2 * k) * m], v[(i + n) + (2 * k) * m]);
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

fn eigvals3(hre: &[[f64; 3]; 3], him: &[[f64; 3]; 3]) -> [f64; 3] {
    let p1 = hre[0][1] * hre[0][1]
        + him[0][1] * him[0][1]
        + hre[0][2] * hre[0][2]
        + him[0][2] * him[0][2]
        + hre[1][2] * hre[1][2]
        + him[1][2] * him[1][2];
    let q = (hre[0][0] + hre[1][1] + hre[2][2]) / 3.0;
    let p2 = (hre[0][0] - q).powi(2) + (hre[1][1] - q).powi(2) + (hre[2][2] - q).powi(2) + 2.0 * p1;
    if p2 < 1e-300 {
        return [q, q, q];
    }
    let p = (p2 / 6.0).sqrt();
    let bd = [
        (hre[0][0] - q) / p,
        (hre[1][1] - q) / p,
        (hre[2][2] - q) / p,
    ];
    let (b01r, b01i) = (hre[0][1] / p, him[0][1] / p);
    let (b02r, b02i) = (hre[0][2] / p, him[0][2] / p);
    let (b12r, b12i) = (hre[1][2] / p, him[1][2] / p);
    let tr_re = (b01r * b12r - b01i * b12i) * b02r + (b01r * b12i + b01i * b12r) * b02i;
    let det = bd[0] * bd[1] * bd[2] + 2.0 * tr_re
        - bd[0] * (b12r * b12r + b12i * b12i)
        - bd[1] * (b02r * b02r + b02i * b02i)
        - bd[2] * (b01r * b01r + b01i * b01i);
    let r = (det / 2.0).clamp(-1.0, 1.0);
    let phi = r.acos() / 3.0;
    let e1 = q + 2.0 * p * phi.cos();
    let e3 = q + 2.0 * p * (phi + 2.0 * std::f64::consts::PI / 3.0).cos();
    let e2 = 3.0 * q - e1 - e3;
    let mut v = [e3, e2, e1];
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v
}

fn localize_unsorted(modes: &[C3v]) -> (Vec<C3v>, Vec<f64>) {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut ure = [[0.0f64; 3]; 3];
    let mut uim = [[0.0f64; 3]; 3];
    for a in 0..Q {
        for b in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for i in 0..NS {
                let x = (i % N) as f64;
                let (sn, cs) = (two_pi * x / N as f64).sin_cos();
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
    let mut out: Vec<C3v> = Vec::new();
    let mut centers = Vec::new();
    for k in 0..Q {
        let mut psi = [(0.0f64, 0.0f64); NS];
        for i in 0..NS {
            for a in 0..Q {
                let (cr, ci) = vecs[k][a];
                let (mr, mi) = modes[a][i];
                psi[i].0 += cr * mr - ci * mi;
                psi[i].1 += cr * mi + ci * mr;
            }
        }
        let (mut zr, mut zi) = (0.0, 0.0);
        for i in 0..NS {
            let p = psi[i].0 * psi[i].0 + psi[i].1 * psi[i].1;
            let x = (i % N) as f64;
            let (sn, cs) = (two_pi * x / N as f64).sin_cos();
            zr += p * cs;
            zi += p * sn;
        }
        let center = (zi.atan2(zr) / two_pi * N as f64).rem_euclid(N as f64);
        out.push(psi);
        centers.push(center);
    }
    (out, centers)
}

fn order_stable(centers: &[f64]) -> Vec<usize> {
    let snapped: Vec<f64> = centers
        .iter()
        .map(|&c| ((2.0 * c).round() / 2.0).rem_euclid(N as f64))
        .collect();
    let mut ord: Vec<usize> = (0..centers.len()).collect();
    ord.sort_by(|&a, &b| snapped[a].partial_cmp(&snapped[b]).unwrap());
    ord
}

fn yukawa(la: &[C3v], lb: &[C3v], sig_h: f64) -> M3 {
    let mut phih = [0.0f64; NS];
    for y in 0..N {
        for x in 0..N {
            let dx = (x as f64).min(N as f64 - x as f64);
            let dy = (y as f64).min(N as f64 - y as f64);
            phih[x + y * N] = (-(dx * dx + dy * dy) / (2.0 * sig_h * sig_h)).exp();
        }
    }
    let mut y_out = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..Q {
        for j in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for s in 0..NS {
                let (ar, ai) = la[i][s];
                let (br, bi) = lb[j][s];
                sr += (ar * br + ai * bi) * phih[s];
                si += (ar * bi - ai * br) * phih[s];
            }
            y_out[i][j] = (sr, si);
        }
    }
    y_out
}

/// nt トーラスの対付き積湯川。a,b は Z₆^nt の複合 Wilson 添字、pf,pg は
/// (S₃)^(nt−1) の複合対添字 (トーラス 1 は常に恒等)。
fn pair_y(ytab: &[M3], nt: usize, a: usize, b: usize, pf: usize, pg: usize) -> M3 {
    let (mut aa, mut bb) = (a, b);
    let (a1, b1) = (2 * (aa % 6), 2 * (bb % 6));
    let mut y = ytab[a1 + b1 * NK12];
    let (mut pfc, mut pgc) = (pf, pg);
    for _ in 1..nt {
        aa /= 6;
        bb /= 6;
        let (at, bt) = (2 * (aa % 6), 2 * (bb % 6));
        let yt = &ytab[at + bt * NK12];
        let prow = &PERMS[pfc % 6];
        let pcol = &PERMS[pgc % 6];
        pfc /= 6;
        pgc /= 6;
        let mut z = [[(0.0f64, 0.0f64); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (p, q) = y[i][j];
                let (r, s) = yt[prow[i]][pcol[j]];
                z[i][j] = (p * r - q * s, p * s + q * r);
            }
        }
        y = z;
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

fn mass_ratios(y: &M3) -> [f64; 2] {
    let (hre, him) = gram(y);
    let lam = eigvals3(&hre, &him);
    let sv = [
        lam[0].max(0.0).sqrt(),
        lam[1].max(0.0).sqrt(),
        lam[2].max(0.0).sqrt(),
    ];
    [
        (sv[0].max(1e-300) / sv[2].max(1e-300)).ln(),
        (sv[1].max(1e-300) / sv[2].max(1e-300)).ln(),
    ]
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

fn ckm3(vu: &M3, vd: &M3) -> [f64; 3] {
    let mut ckm = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = vu[i][k];
                let (c, d) = vd[j][k];
                re += a * c + b * d;
                im += a * d - b * c;
            }
            ckm[i][j] = (re * re + im * im).sqrt();
        }
    }
    [ckm[0][1], ckm[1][2], ckm[0][2]]
}

fn lse(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    if m == f64::NEG_INFINITY {
        return m;
    }
    m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

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
    fn ln(&self) -> f64 {
        if self.m == f64::NEG_INFINITY {
            f64::NEG_INFINITY
        } else {
            self.m + self.s.ln()
        }
    }
}

/// 質量のみの証拠 (対 marginalize 有無を perm フラグで切替)。
fn eval_mass(nt: usize, locs: &[Vec<C3v>], sig_grid: &[f64], sigma: f64, perm: bool) -> f64 {
    let nc = 6usize.pow(nt as u32);
    let np = if perm { 6usize.pow(nt as u32 - 1) } else { 1 };
    let norm2 = -(2.0 * std::f64::consts::PI * sigma * sigma).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + norm2
    };
    let mut terms = Vec::new();
    for &sh in sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        let mut acc_u: Vec<Acc> = (0..nc).map(|_| Acc::new()).collect();
        let mut acc_d: Vec<Acc> = (0..nc).map(|_| Acc::new()).collect();
        let mut acc_e = Acc::new();
        for a in 0..nc {
            for b in 0..nc {
                for p in 0..np {
                    let r = mass_ratios(&pair_y(&ytab, nt, a, b, 0, p));
                    acc_u[a].add(ll2(&r, tgt[0], tgt[1]));
                    acc_d[a].add(ll2(&r, tgt[2], tgt[3]));
                }
                // e セクター: σL, σe が独立 (np² 通り)
                for pl in 0..np {
                    for pe in 0..np {
                        let r = mass_ratios(&pair_y(&ytab, nt, a, b, pl, pe));
                        acc_e.add(ll2(&r, tgt[4], tgt[5]));
                    }
                }
            }
        }
        let per_q: Vec<f64> = (0..nc).map(|k| acc_u[k].ln() + acc_d[k].ln()).collect();
        terms.push(lse(&per_q) + acc_e.ln());
    }
    let mut prior = 5.0 * (nt as f64) * (6.0f64).ln() + (sig_grid.len() as f64).ln();
    if perm {
        prior += 4.0 * (np as f64).ln();
    }
    lse(&terms) - prior
}

/// 全 9 量の証拠。perm=true では (ku,σu)×(kd,σd) を質量尤度で降順に並べ、
/// llu+lld ≥ max−Δ の対だけ CKM を計算。残差は e^{3norm1}·(SuSd−計算済み質量重み)
/// で厳密に抑え、lnZ の区間 (下端, 上端) を返す。
fn eval9(nt: usize, locs: &[Vec<C3v>], sig_grid: &[f64], sigma: f64, perm: bool) -> (f64, f64) {
    let nc = 6usize.pow(nt as u32);
    let np = if perm { 6usize.pow(nt as u32 - 1) } else { 1 };
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let llck_max = 3.0 * norm1;
    let mut terms_lo = Vec::new();
    let mut terms_hi = Vec::new();
    for &sh in sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // e セクター (streaming)
        let mut acc_e = Acc::new();
        for a in 0..nc {
            for b in 0..nc {
                for pl in 0..np {
                    for pe in 0..np {
                        let r = mass_ratios(&pair_y(&ytab, nt, a, b, pl, pe));
                        acc_e.add(ll2(&r, tgt[4], tgt[5]));
                    }
                }
            }
        }
        // クォーク部: kq ごとに (ku,p) リストを作り、打ち切り付き二重和
        let mut per_q_lo: Vec<f64> = Vec::with_capacity(nc);
        let mut per_q_hi: Vec<f64> = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut ulist: Vec<(f64, M3)> = Vec::with_capacity(nc * np);
            let mut dlist: Vec<(f64, M3)> = Vec::with_capacity(nc * np);
            for k2 in 0..nc {
                for p in 0..np {
                    let (r, v) = mass_and_vecs(&pair_y(&ytab, nt, kq, k2, 0, p));
                    ulist.push((ll2(&r, tgt[0], tgt[1]), v));
                    dlist.push((ll2(&r, tgt[2], tgt[3]), v));
                }
            }
            ulist.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
            dlist.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
            let su_all = lse(&ulist.iter().map(|x| x.0).collect::<Vec<_>>());
            let sd_all = lse(&dlist.iter().map(|x| x.0).collect::<Vec<_>>());
            let cut = ulist[0].0 + dlist[0].0 - DELTA_CUT;
            let mut acc = Acc::new(); // CKM 込みの計算済み和
            let mut accw = Acc::new(); // 計算済み対の質量重み (残差計算用)
            for (llu, vu) in &ulist {
                if llu + dlist[0].0 < cut {
                    break;
                }
                for (lld, vd) in &dlist {
                    if llu + lld < cut {
                        break;
                    }
                    let c = ckm3(vu, vd);
                    let mut ll = llu + lld;
                    for m in 0..3 {
                        let d = c[m].max(1e-300).ln() - tgt[6 + m];
                        ll += -d * d / (2.0 * sigma * sigma) + norm1;
                    }
                    acc.add(ll);
                    accw.add(llu + lld);
                }
            }
            // 残差上界: e^{llck_max} · (Su·Sd − 計算済み質量重み)
            let total_w = su_all + sd_all;
            let done_w = accw.ln();
            let rem = if done_w >= total_w {
                f64::NEG_INFINITY
            } else {
                llck_max + total_w + (1.0 - (done_w - total_w).exp()).max(0.0).ln()
            };
            let lo = acc.ln();
            let hi = lse(&[lo, rem]);
            per_q_lo.push(lo);
            per_q_hi.push(hi);
        }
        terms_lo.push(lse(&per_q_lo) + acc_e.ln());
        terms_hi.push(lse(&per_q_hi) + acc_e.ln());
    }
    let mut prior = 5.0 * (nt as f64) * (6.0f64).ln() + (sig_grid.len() as f64).ln();
    if perm {
        prior += 4.0 * (np as f64).ln();
    }
    (lse(&terms_lo) - prior, lse(&terms_hi) - prior)
}

fn main() {
    self_test();
    println!("=== v10.2 対込みの幾何選択: T³ は勝ち続けるか ===\n");
    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];

    println!("[0] 世代モード (Z₆ ⊂ Z₁₂, 対角化 12 回, 安定ラベル)");
    let t0 = std::time::Instant::now();
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    let mut ok_engine = true;
    for k in 0..NK12 {
        let (modes, gap, spread) = flux_modes(k);
        if spread > 1e-9 || gap < 0.05 {
            ok_engine = false;
        }
        let (raw, cents) = localize_unsorted(&modes);
        let ord = order_stable(&cents);
        locs.push(ord.iter().map(|&i| raw[i]).collect());
    }
    println!(
        "    縮退・ギャップ不変  {}  ({} ms)",
        pass(ok_engine),
        t0.elapsed().as_millis()
    );

    // ---- [1] 回帰検査 (装置検証) ----
    println!("\n[1] 回帰検査 — 対角経路 (v9.2) と T² 自由対 (v10.1) の再現");
    let mut all_reg = ok_engine;
    let mut check = |name: &str, got: f64, want: f64, tol: f64| {
        let ok = (got - want).abs() < tol;
        all_reg &= ok;
        println!(
            "    {}: {:+.4} vs {:+.4} (|Δ|<{})  {}",
            name,
            got,
            want,
            tol,
            pass(ok)
        );
    };
    let t1 = std::time::Instant::now();
    let m_t1 = eval_mass(1, &locs, &sig_grid, sigma, false);
    check("T¹ 質量        ", m_t1, REF_T1_MASS, 1e-6);
    let (n_t1_lo, n_t1_hi) = eval9(1, &locs, &sig_grid, sigma, false);
    check("T¹ 全9 (下端)   ", n_t1_lo, REF_T1_NINE, 0.01);
    let m_t2d = eval_mass(2, &locs, &sig_grid, sigma, false);
    check("T² 対角 質量    ", m_t2d, REF_T2_MASS_DIAG, 1e-6);
    let m_t3d = eval_mass(3, &locs, &sig_grid, sigma, false);
    check("T³ 対角 質量    ", m_t3d, REF_T3_MASS_DIAG, 1e-6);
    let (n_t3d_lo, n_t3d_hi) = eval9(3, &locs, &sig_grid, sigma, false);
    check("T³ 対角 全9 (下端)", n_t3d_lo, REF_T3_NINE_DIAG, 0.01);
    let m_t2p = eval_mass(2, &locs, &sig_grid, sigma, true);
    check("T² 自由対 質量  ", m_t2p, REF_T2_MASS_PERM, 1e-6);
    let (n_t2p_lo, n_t2p_hi) = eval9(2, &locs, &sig_grid, sigma, true);
    check("T² 自由対 全9   ", n_t2p_lo, REF_T2_NINE_PERM, 0.01);
    println!(
        "    打ち切り区間幅: T¹ {:.2e} / T³対角 {:.2e} / T²自由対 {:.2e}  ({} ms)",
        n_t1_hi - n_t1_lo,
        n_t3d_hi - n_t3d_lo,
        n_t2p_hi - n_t2p_lo,
        t1.elapsed().as_millis()
    );

    // ---- [2] 本番: T³ 自由対 ----
    println!(
        "\n[2] T³ 自由対 (対空間 (S₃×S₃)⁴ = 1.68e6, Occam 罰 4ln36 = {:.2})",
        4.0 * (36.0f64).ln()
    );
    let t2 = std::time::Instant::now();
    let m_t3p = eval_mass(3, &locs, &sig_grid, sigma, true);
    println!(
        "    質量のみ lnZ(T³ perm) = {:.4}  ({} ms)",
        m_t3p,
        t2.elapsed().as_millis()
    );
    let t3 = std::time::Instant::now();
    let (n_t3p_lo, n_t3p_hi) = eval9(3, &locs, &sig_grid, sigma, true);
    let width = n_t3p_hi - n_t3p_lo;
    let ok_width = width < 0.01;
    println!(
        "    全 9 量 lnZ₉(T³ perm) ∈ [{:.4}, {:.4}] (幅 {:.2e} < 0.01)  {}  ({} ms)",
        n_t3p_lo,
        n_t3p_hi,
        width,
        pass(ok_width),
        t3.elapsed().as_millis()
    );

    // ---- [3] 判定 ----
    println!("\n[3] 対込みの幾何選択 (全て安定ラベル)");
    println!("    幾何        質量のみ      全 9 量");
    println!("    T¹          {:+.2}        {:+.2}", m_t1, n_t1_lo);
    println!(
        "    T² 対角     {:+.2}        {:+.2}",
        m_t2d, REF_T2_NINE_DIAG_PRINT
    );
    println!("    T² 自由対   {:+.2}        {:+.2}", m_t2p, n_t2p_lo);
    println!("    T³ 対角     {:+.2}        {:+.2}", m_t3d, n_t3d_lo);
    println!(
        "    T³ 自由対   {:+.2}        [{:.2}, {:.2}]",
        m_t3p, n_t3p_lo, n_t3p_hi
    );
    let t3_wins_mass = m_t3p > m_t2p;
    let t3_wins_nine = n_t3p_lo > n_t2p_lo; // 区間下端で比較 (上端でも同じ側なら確定)
    let nine_decided = (n_t3p_lo > n_t2p_lo + 0.02) || (n_t3p_hi < n_t2p_lo - 0.02);
    println!(
        "\n    質量のみ: {} (T³−T² = {:+.2})",
        if t3_wins_mass {
            "T³ が勝つ"
        } else {
            "T² が勝つ"
        },
        m_t3p - m_t2p
    );
    println!(
        "    全 9 量:  {} (T³−T² = {:+.2}, 区間で{}確定)",
        if t3_wins_nine {
            "T³ が勝つ"
        } else {
            "T² が勝つ"
        },
        n_t3p_lo - n_t2p_lo,
        if nine_decided { "" } else { "未" }
    );

    let all_ok = all_reg && ok_width;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-008".into())),
        ("delta_cut_nats".into(), Json::Num(DELTA_CUT)),
        (
            "lnZ_mass".into(),
            Json::Obj(vec![
                ("T1".into(), Json::Num(m_t1)),
                ("T2_diag".into(), Json::Num(m_t2d)),
                ("T2_perm".into(), Json::Num(m_t2p)),
                ("T3_diag".into(), Json::Num(m_t3d)),
                ("T3_perm".into(), Json::Num(m_t3p)),
            ]),
        ),
        (
            "lnZ_nine".into(),
            Json::Obj(vec![
                ("T1".into(), Json::Num(n_t1_lo)),
                ("T2_perm".into(), Json::Num(n_t2p_lo)),
                ("T3_diag".into(), Json::Num(n_t3d_lo)),
                ("T3_perm_lo".into(), Json::Num(n_t3p_lo)),
                ("T3_perm_hi".into(), Json::Num(n_t3p_hi)),
            ]),
        ),
        ("t3_wins_mass".into(), Json::Bool(t3_wins_mass)),
        ("t3_wins_nine".into(), Json::Bool(t3_wins_nine)),
        ("nine_decided".into(), Json::Bool(nine_decided)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v102_geopair.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!(
        "\n総合判定: {} (PASS = 装置検証と区間幅 — 選択の答えは [3] の表が本体)",
        pass(all_ok)
    );
    if !all_ok {
        std::process::exit(1);
    }
}

// v9.2 の T² 対角 全9 (表示用; results/v92_labelstab.json)
const REF_T2_NINE_DIAG_PRINT: f64 = -23.611_871_688_966_34;
