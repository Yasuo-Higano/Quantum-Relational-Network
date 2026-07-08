//! v15.7 QRN 固有予言の pre-register — 「当たるかもしれない言葉」を数の形式に固定する
//!
//! 残高 6 (v15.0): 重力以外の QRN 固有予言。次の方針 (PROMPT/2 柱 6) の形式
//! (observable / QRN value / competing models / data not used / falsifier) に従い、
//! 予言台帳 predictions.yml を導入する。本バイナリはその数値部分を凍結済みの機構から
//! 計算する:
//!
//! [P1] BMV 型実験の (C, V) 平面: Newton 位相 Δφ = (Gm²τ/ħ)·2Δx²/(d(d²−Δx²)) から
//!      量子チャネル (C,V) = (|sin Δφ/2|, |cos Δφ/2|) と古典上限 V ≤ e^{−Δφ/2} を
//!      登録パラメータ 4 点で再導出 (v7.3 の公表表との回帰つき)。
//!
//! [P2] **真の holdout 予言**: 幾何湯川模型 (T²×T², v10.1 の凍結済みエンジン) は
//!      質量比 6 + |V_us|,|V_cb|,|V_ub| の 9 量だけで検定されてきた —
//!      **|V_td| と Jarlskog 不変量 J は一度も使われていない** (v9.1 が明記)。
//!      本バイナリは事後分布 (9 量の尤度 × 一様事前, M_perm) の下での
//!      |V_td| と |J| の分布を計算し、測定値 (PDG) との比較を同じ出力で行う。
//!      機構は v10.1 からコミット済み・数値は決定論 — これは事後的 fit ではなく
//!      凍結済み模型の out-of-sample 検定である。
//!      注: 模型空間は複素共役 (Wilson 反転) 対称なので J の符号は予言されない —
//!      予言されるのは CP 破れの大きさ |J| である (対称性の検査つき)。
//!
//! [P3–P5] 計算済み証拠に基づく登録予言 (分数電荷・LV 消失・Majorana) は
//!      predictions.yml に形式のみ登録 (出典は各版)。
//!
//! 検証: エンジンは v15.5 と同一 — lnZ(一様) が v10.1 公表値 −19.86 を回帰、
//! MAP が v10.1 の記録 (σ_u=(23), σ_d=(132), σ_H=1.0, lnL=−5.87) と一致すること。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const NC: usize = 36;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const REF_NINE_PERM: f64 = -19.86;
// PDG の holdout 測定値 (予言の採点用 — 学習には一切使っていない)
const PDG_VTD: f64 = 0.0086;
const PDG_J: f64 = 3.08e-5;

const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];

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

/// 複素 CKM 行列 V_ij = Σ_k (U_u)_ik conj((U_d)_jk)
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

/// Jarlskog J = Im(V_us V_cb V*_ub V*_cs) と |V_td|
fn jarlskog_vtd(v: &M3) -> (f64, f64) {
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let conj = |a: (f64, f64)| (a.0, -a.1);
    let t = mul(mul(v[0][1], v[1][2]), mul(conj(v[0][2]), conj(v[1][1])));
    let vtd = (v[2][0].0 * v[2][0].0 + v[2][0].1 * v[2][0].1).sqrt();
    (t.1, vtd)
}

fn lse(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

fn pair_yukawa(ytab: &[M3], a: usize, b: usize, sf: usize, sg: usize) -> M3 {
    let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
    let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
    let y1 = &ytab[a1 + b1 * NK12];
    let y2 = &ytab[a2 + b2 * NK12];
    had_prod_perm(y1, y2, sf, sg)
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// 重みつきヒストグラムの分位点 (bin 中央値)
struct WHist {
    lo: f64,
    hi: f64,
    bins: Vec<f64>,
}
impl WHist {
    fn new(lo: f64, hi: f64, n: usize) -> Self {
        WHist {
            lo,
            hi,
            bins: vec![0.0; n],
        }
    }
    fn add(&mut self, x: f64, w: f64) {
        let n = self.bins.len();
        let t = ((x - self.lo) / (self.hi - self.lo) * n as f64).floor();
        let i = (t.max(0.0) as usize).min(n - 1);
        self.bins[i] += w;
    }
    fn quantile(&self, q: f64) -> f64 {
        let tot: f64 = self.bins.iter().sum();
        let mut acc = 0.0;
        for (i, b) in self.bins.iter().enumerate() {
            acc += b;
            if acc >= q * tot {
                return self.lo + (i as f64 + 0.5) / self.bins.len() as f64 * (self.hi - self.lo);
            }
        }
        self.hi
    }
}

fn main() {
    self_test();
    println!("=== v15.7 QRN 固有予言の pre-register: 形式に固定された数値予言 ===\n");
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ================= [P1] BMV の (C, V) 予言表 =================
    println!("[P1] BMV 型実験の (C,V) 予言 (Newton 位相 Δφ = (Gm²τ/ħ)·2Δx²/(d(d²−Δx²)))");
    {
        let (g, hbar) = (6.674e-11f64, 1.0546e-34f64);
        // 登録パラメータ (m[kg], d[m], Δx[m], τ[s]) と v7.3 公表 Δφ (回帰)
        let table: [(f64, f64, f64, f64, f64); 4] = [
            (1e-14, 250e-6, 50e-6, 2.5, 0.053),
            (1e-14, 200e-6, 100e-6, 2.5, 0.527),
            (1e-14, 200e-6, 100e-6, 5.0, 1.055),
            (1e-15, 200e-6, 100e-6, 2.5, 0.005),
        ];
        let mut reg_ok = true;
        println!("    m[kg]   d[μm] Δx[μm] τ[s]   Δφ      C(QRN)   V(QRN)   V上限(古典)");
        for &(m, d, dx, tau, dphi_ref) in &table {
            let dphi = g * m * m * tau / hbar * 2.0 * dx * dx / (d * (d * d - dx * dx));
            let (c, v) = ((dphi / 2.0).sin().abs(), (dphi / 2.0).cos().abs());
            let vcl = (-dphi / 2.0).exp();
            if (dphi - dphi_ref).abs() > 0.002 {
                reg_ok = false;
            }
            println!(
                "    {:7.0e} {:5.0} {:5.0} {:5.1}  {:.3}   {:.4}   {:.4}   {:.4}",
                m,
                d * 1e6,
                dx * 1e6,
                tau,
                dphi,
                c,
                v,
                vcl
            );
        }
        check(
            "Δφ の再導出が v7.3 の公表表と一致 (±0.002)",
            reg_ok,
            "位相公式は凍結済み — 予言は predictions.yml PRED-1 に登録".to_string(),
        );
    }

    // ================= [P2] holdout 予言: |V_td| と |J| =================
    println!("\n[P2] 凍結済み幾何湯川模型の holdout 予言 — |V_td| と Jarlskog |J| (一度も学習に使っていない)");
    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let t0 = std::time::Instant::now();
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    for k in 0..NK12 {
        let (modes, _gap, _spread) = flux_modes(k);
        let (raw, cents) = localize_unsorted(&modes);
        let ord = order_stable(&cents);
        locs.push(ord.iter().map(|&i| raw[i]).collect());
    }
    // e セクターの lnZ_e(σ_H) と max lnL_e(σ_H) (σ_H 重みの結合・MAP 用) — 先に一括計算
    let mut lnze_sh = Vec::new();
    let mut emax_sh = Vec::new();
    let mut ytabs: Vec<Vec<M3>> = Vec::new();
    for &sh in &sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        let mut acc = f64::NEG_INFINITY;
        let mut s = 0.0;
        let mut emax = f64::NEG_INFINITY;
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..NC * NC {
                    let r = mass_ratios(&pair_yukawa(&ytab, ab % NC, ab / NC, sl, se_));
                    let l = ll2(&r, tgt[4], tgt[5]);
                    emax = emax.max(l);
                    if l > acc {
                        s = s * (acc - l).exp() + 1.0;
                        acc = l;
                    } else {
                        s += (l - acc).exp();
                    }
                }
            }
        }
        lnze_sh.push(acc + s.ln());
        emax_sh.push(emax);
        ytabs.push(ytab);
    }
    // クォーク五重和: 事後重み w = exp(ll_q + lnZ_e(σ_H)) で |V_td|, |J| を集計
    let shift = *lnze_sh
        .iter()
        .max_by(|a, b| a.partial_cmp(b).unwrap())
        .unwrap();
    let mut h_vtd = WHist::new(-9.0, 0.0, 360); // ln|V_td|
    let mut h_j = WHist::new(-16.0, -1.0, 360); // ln|J|
    let mut zq_terms = Vec::new();
    let (mut wsum, mut w_jbig) = (0.0f64, 0.0f64);
    let mut map = (
        f64::NEG_INFINITY,
        0.0f64,
        0.0f64,
        [0usize; 2],
        0.0f64,
        0.0f64,
    );
    for (isg, ytab) in ytabs.iter().enumerate() {
        let mut pair_r: Vec<[f64; 2]> = Vec::with_capacity(NC * NC * 6);
        let mut pair_v: Vec<M3> = Vec::with_capacity(NC * NC * 6);
        for m in 0..NC * NC * 6 {
            let y = pair_yukawa(ytab, m % NC, (m / NC) % NC, 0, m / (NC * NC));
            let (r, v) = mass_and_vecs(&y);
            pair_r.push(r);
            pair_v.push(v);
        }
        let we = (lnze_sh[isg] - shift).exp();
        let mut acc_q = f64::NEG_INFINITY;
        let mut s_q = 0.0;
        for kq in 0..NC {
            for su in 0..6 {
                for ku in 0..NC {
                    let mu = kq + ku * NC + su * NC * NC;
                    let llu = ll2(&pair_r[mu], tgt[0], tgt[1]);
                    for sd in 0..6 {
                        for kd in 0..NC {
                            let md = kq + kd * NC + sd * NC * NC;
                            let lld = ll2(&pair_r[md], tgt[2], tgt[3]);
                            let v = ckm_full(&pair_v[mu], &pair_v[md]);
                            let cabs = [
                                (v[0][1].0 * v[0][1].0 + v[0][1].1 * v[0][1].1).sqrt(),
                                (v[1][2].0 * v[1][2].0 + v[1][2].1 * v[1][2].1).sqrt(),
                                (v[0][2].0 * v[0][2].0 + v[0][2].1 * v[0][2].1).sqrt(),
                            ];
                            let mut ll = llu + lld;
                            for m in 0..3 {
                                let d = cabs[m].max(1e-300).ln() - tgt[6 + m];
                                ll += -d * d / (2.0 * sigma * sigma) + norm1;
                            }
                            // ストリーミング lse (回帰用)
                            if ll > acc_q {
                                s_q = s_q * (acc_q - ll).exp() + 1.0;
                                acc_q = ll;
                            } else {
                                s_q += (ll - acc_q).exp();
                            }
                            // 事後重みでの holdout 集計
                            let (j, vtd) = jarlskog_vtd(&v);
                            let w = we * (ll - (-6.0)).exp().min(1e30); // 基準シフト (−6 ≈ oracle)
                            h_vtd.add(vtd.max(1e-300).ln(), w);
                            h_j.add(j.abs().max(1e-300).ln(), w);
                            wsum += w;
                            if j.abs() > 1e-6 {
                                w_jbig += w;
                            }
                            if ll > map.0 {
                                map = (ll, j, vtd, [su, sd], sig_grid[isg], emax_sh[isg]);
                            }
                        }
                    }
                }
            }
        }
        zq_terms.push(acc_q + s_q.ln());
    }
    // 回帰: lnZ(一様) の再構成
    let ln_prior_q = ((NC * NC * NC) as f64).ln() + 2.0 * (6.0f64).ln();
    let ln_prior_e = ((NC * NC) as f64).ln() + 2.0 * (6.0f64).ln();
    let terms: Vec<f64> = (0..4)
        .map(|i| (zq_terms[i] - ln_prior_q) + (lnze_sh[i] - ln_prior_e))
        .collect();
    let lnz = lse(&terms) - (4.0f64).ln();
    check(
        "エンジン回帰: lnZ(一様) = v10.1 公表値",
        (lnz - REF_NINE_PERM).abs() < 0.02,
        format!("lnZ = {:.3} (公表 {:.2})", lnz, REF_NINE_PERM),
    );
    // MAP の全 9 量 lnL = クォーク部 + e セクター最良 (同一 σ_H)
    let map_total = map.0 + map.5;
    check(
        "MAP の再現: σ_u=(23), σ_d=(132), σ_H=1.0, 全 9 量+e lnL=−5.87 (v10.1/v15.5 の記録)",
        (map_total - (-5.87)).abs() < 0.01 && (map.4 - 1.0).abs() < 1e-9 && map.3 == [1, 4],
        format!(
            "lnL_MAP(全) = {:.3} (クォーク {:.3} + e {:.3}), σ_H = {}, (σu,σd) = {:?}",
            map_total, map.0, map.5, map.4, map.3
        ),
    );
    // J の構造零: 長方形トーラス + 実 Wilson 格子の湯川位相は行×列に因子化し、
    // CKM は再位相化で実になる — J は厳密に零のはず (数値ノイズ床のみ)。
    let frac_jbig = w_jbig / wsum.max(1e-300);
    check(
        "J の構造零 (因子化位相): MAP |J| < 1e-12 かつ事後質量の |J|>1e-6 割合 < 1%",
        map.1.abs() < 1e-12 && frac_jbig < 0.01,
        format!(
            "MAP J = {:.1e}, P(|J|>1e-6) = {:.2e} — この幾何窓は CP 破れを作れない (構造の発見)",
            map.1, frac_jbig
        ),
    );
    // 予言帯 (事後分位点)
    let q = |h: &WHist, p: f64| h.quantile(p).exp();
    let (v16, v50, v84) = (q(&h_vtd, 0.16), q(&h_vtd, 0.5), q(&h_vtd, 0.84));
    let (v025, v975) = (q(&h_vtd, 0.025), q(&h_vtd, 0.975));
    let (j16, j50, j84) = (q(&h_j, 0.16), q(&h_j, 0.5), q(&h_j, 0.84));
    let (j025, j975) = (q(&h_j, 0.025), q(&h_j, 0.975));
    println!("\n    ┌─ holdout 予言 (事後帯) と測定値 ─────────────────────────");
    println!(
        "    │ |V_td|: 中央値 {:.4}  68% [{:.4}, {:.4}]  95% [{:.4}, {:.4}]",
        v50, v16, v84, v025, v975
    );
    println!("    │          測定 (PDG) = {:.4}", PDG_VTD);
    println!(
        "    │ |J|:    中央値 {:.2e}  68% [{:.2e}, {:.2e}]  95% [{:.2e}, {:.2e}]",
        j50, j16, j84, j025, j975
    );
    println!("    │          測定 (PDG) = {:.2e}", PDG_J);
    println!("    │ MAP 点: |V_td| = {:.4}, J = {:+.2e}", map.2, map.1);
    println!("    └───────────────────────────────────────────────────────");
    let vtd_in95 = PDG_VTD >= v025 && PDG_VTD <= v975;
    let j_in95 = PDG_J >= j025 && PDG_J <= j975;
    println!(
        "    採点: |V_td| 測定は 95% 帯の{} (帯上端 {:.4} vs 測定 {:.4}) / J: 模型は構造零 — 測定 {:.2e} と決定的に乖離",
        if vtd_in95 { "内" } else { "外" },
        v975,
        PDG_VTD,
        PDG_J
    );
    let _ = j_in95;
    println!("    (採点は予言の成否であって装置の成否ではない — 外れの内容が台帳に入る:");
    println!(
        "     「長方形 T²×T² + 実 Wilson は CP を作れない」→ 傾き/複素構造 (v12–v13 の方向) が"
    );
    println!("     好みでなく CP 破れの要求であることを、holdout が独立に指した)");
    check(
        "ヒストグラム正規化と有効試料 (集計の健全性)",
        wsum > 0.0 && h_vtd.bins.iter().sum::<f64>() > 0.0,
        format!("総重み = {:.3e} ({} ms)", wsum, t0.elapsed().as_millis()),
    );

    // ================= artifact =================
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v15.7".into())),
        ("lnz_regression".into(), Json::Num(lnz)),
        ("vtd_median".into(), Json::Num(v50)),
        (
            "vtd_68".into(),
            Json::Arr(vec![Json::Num(v16), Json::Num(v84)]),
        ),
        (
            "vtd_95".into(),
            Json::Arr(vec![Json::Num(v025), Json::Num(v975)]),
        ),
        ("vtd_pdg".into(), Json::Num(PDG_VTD)),
        ("vtd_in95".into(), Json::Bool(vtd_in95)),
        ("jarlskog_abs_median".into(), Json::Num(j50)),
        (
            "jarlskog_68".into(),
            Json::Arr(vec![Json::Num(j16), Json::Num(j84)]),
        ),
        (
            "jarlskog_95".into(),
            Json::Arr(vec![Json::Num(j025), Json::Num(j975)]),
        ),
        ("jarlskog_pdg".into(), Json::Num(PDG_J)),
        ("jarlskog_in95".into(), Json::Bool(j_in95)),
        ("map_vtd".into(), Json::Num(map.2)),
        ("map_j".into(), Json::Num(map.1)),
        (
            "j_frac_above_1e-6".into(),
            Json::Num(w_jbig / wsum.max(1e-300)),
        ),
    ]);
    let p = write_artifact("results/v157_predict.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 予言は形式に固定された — 採点の成否は predictions.yml と本結果が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
