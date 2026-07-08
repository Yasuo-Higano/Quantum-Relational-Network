//! v16.5 勝者シアー (1,1) の全量スコアカード — 模型は何を当て、何を外すか
//!
//! v16.4 は 10 量証拠で最小対称シアー (1,1) を勝者とした (+306 nats の反転、
//! |V_td| holdout が hit)。本バイナリは v9.1 型の点検を行う: 勝者幾何の
//! 10 量 MAP 構成 (σ_H, kQ, ku, kd, σu, σd; e セクターは因子化して独立に MAP) で
//! **全観測量を実測と factor 比較**する:
//!
//!   尤度に入れた 10 量: 質量比 6 (m_u/m_t, m_c/m_t, m_d/m_b, m_s/m_b, m_e/m_τ, m_μ/m_τ)
//!                      + |V_us|, |V_cb|, |V_ub| + |J|
//!   holdout (尤度に入れない): |V_td| (v16.4 で帯の採点済 — ここでは MAP 点)
//!                            **|V_ts| (本プログラムで一度も使っていない新しい holdout)**
//!   導出量 (独立でない): CKM 位相 δ (J と混合角から従う — 参考表示)
//!
//! 検証: lnZ₁₀(1,1) = −24.29 (v16.4) の回帰。スコアカード自体は測定結果 (内容)。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const NC: usize = 36;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const REF_LNZ10_WINNER: f64 = -24.29; // v16.4 の勝者 (1,1) の公表値
const J_OBS: f64 = 3.08e-5; // 測定 Jarlskog (PDG)
const VTD_OBS: f64 = 0.0086; // |V_td| (holdout — 尤度に入れない)
const SIG_GRID: [f64; 4] = [1.0, 1.5, 2.0, 2.5];
const SHEARS: [usize; 1] = [1]; // 勝者幾何 (1,1) のみ

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

fn flux_modes_shear(k_half: usize, s: usize) -> (Vec<C3v>, f64, f64) {
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / NS as f64;
    let wl = phi * k_half as f64 / 2.0;
    let idx = |x: usize, y: usize| x + y * N;
    let m = 2 * NS;
    let mut a = vec![0.0; m * m];
    let addhop = |a: &mut Vec<f64>, i: usize, j: usize, th: f64| {
        let (c, sn) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + NS) + (i + NS) * m] += -c;
        a[(i + NS) + (j + NS) * m] += -c;
        a[j + (i + NS) * m] += sn;
        a[(j + NS) + i * m] += -sn;
        a[i + (j + NS) * m] += -sn;
        a[(i + NS) + j * m] += sn;
    };
    for x in 0..N {
        for y in 0..N {
            let th_y = phi * x as f64 + wl;
            if y == N - 1 {
                addhop(&mut a, idx(x, y), idx((x + s) % N, 0), th_y);
            } else {
                addhop(&mut a, idx(x, y), idx(x, y + 1), th_y);
            }
            let th_x = if x == N - 1 {
                -phi * (N as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % N, y), th_x);
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
            for sx in 0..NS {
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

/// 積湯川 (2 つの単一トーラス表から): Y[i][j] = Y¹[i][j] · Y²[σf(i)][σg(j)]
fn pair_yukawa2(yt1: &[M3], yt2: &[M3], a: usize, b: usize, sf: usize, sg: usize) -> M3 {
    let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
    let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
    had_prod_perm(&yt1[a1 + b1 * NK12], &yt2[a2 + b2 * NK12], sf, sg)
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

fn main() {
    self_test();
    println!("=== v16.5 勝者シアー (1,1) の全量スコアカード ===\n");
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // 実測値 (PDG): 尤度に入れた 10 量 + holdout 2 + UT 角 (導出量, 参考)
    let j_obs: f64 = 3.08e-5;
    let vtd_obs: f64 = 0.0086;
    let vts_obs: f64 = 0.0405;
    let (beta_obs, gamma_obs) = (22.2f64, 65.9f64); // 度

    // ---- [0] 勝者幾何 s=1 のモード (12 対角化) ----
    println!("[0] シアー s=1 のモード構築 (12 対角化, 安定ラベル)");
    let t0 = std::time::Instant::now();
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    let mut ok_deg = true;
    for k in 0..NK12 {
        let (modes, gap, spread) = flux_modes_shear(k, 1);
        if spread > 1e-9 || gap < 0.05 {
            ok_deg = false;
        }
        let (raw, cents) = localize_unsorted(&modes);
        let ord = order_stable(&cents);
        locs.push(ord.iter().map(|&i| raw[i]).collect());
    }
    println!("    完了 ({} ms)", t0.elapsed().as_millis());
    check(
        "厳密 3 重縮退・健全ギャップ",
        ok_deg,
        "幅 < 1e-9, ギャップ > 0.05".into(),
    );

    // ---- [1] 10 量証拠と MAP の完全な構成 ----
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let mut terms10 = Vec::new();
    // 全 MAP: (合計 ll10+lle, isg, [kq,ku,kd,su,sd], ru, rd, V, e 側 (lle, [kl,ke,sl,se], re))
    struct MapRec {
        tot: f64,
        sh: f64,
        q: [usize; 5],
        ru: [f64; 2],
        rd: [f64; 2],
        v: M3,
        e_cfg: [usize; 4],
        re: [f64; 2],
    }
    let mut best: Option<MapRec> = None;
    for &sh in SIG_GRID.iter() {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        let pair: Vec<([f64; 2], M3)> = (0..NC * NC * 6)
            .map(|m| pair_yukawa2(&ytab, &ytab, m % NC, (m / NC) % NC, 0, m / (NC * NC)))
            .map(|y| mass_and_vecs(&y))
            .collect();
        // e セクター: lnZe と MAP
        let mut le = Vec::with_capacity(NC * NC * 36);
        let mut e_best = (f64::NEG_INFINITY, [0usize; 4], [0.0f64; 2]);
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..NC * NC {
                    let r = mass_ratios(&pair_yukawa2(&ytab, &ytab, ab % NC, ab / NC, sl, se_));
                    let l = ll2(&r, tgt[4], tgt[5]);
                    le.push(l);
                    if l > e_best.0 {
                        e_best = (l, [ab % NC, ab / NC, sl, se_], r);
                    }
                }
            }
        }
        let lnze = lse(&le);
        // クォーク五重和 (10 量)
        let mut acc10 = (f64::NEG_INFINITY, 0.0f64);
        let mut q_best = (
            f64::NEG_INFINITY,
            [0usize; 5],
            [0.0f64; 2],
            [0.0f64; 2],
            [[(0.0f64, 0.0f64); 3]; 3],
        );
        for kq in 0..NC {
            for su in 0..6 {
                for ku in 0..NC {
                    let mu = kq + ku * NC + su * NC * NC;
                    let (ru, vu) = &pair[mu];
                    let llu = ll2(ru, tgt[0], tgt[1]);
                    for sd in 0..6 {
                        for kd in 0..NC {
                            let md = kq + kd * NC + sd * NC * NC;
                            let (rd, vd) = &pair[md];
                            let lld = ll2(rd, tgt[2], tgt[3]);
                            let v = ckm_full(vu, vd);
                            let c = [cab(&v, 0, 1), cab(&v, 1, 2), cab(&v, 0, 2)];
                            let mut ll = llu + lld;
                            for m in 0..3 {
                                let d = c[m].max(1e-300).ln() - tgt[6 + m];
                                ll += -d * d / (2.0 * sigma * sigma) + norm1;
                            }
                            let j = jarlskog(&v);
                            let dj = j.abs().max(1e-300).ln() - j_obs.ln();
                            let ll10 = ll + (-dj * dj / (2.0 * sigma * sigma) + norm1);
                            if ll10 > acc10.0 {
                                acc10.1 = acc10.1 * (acc10.0 - ll10).exp() + 1.0;
                                acc10.0 = ll10;
                            } else {
                                acc10.1 += (ll10 - acc10.0).exp();
                            }
                            if ll10 > q_best.0 {
                                q_best = (ll10, [kq, ku, kd, su, sd], *ru, *rd, v);
                            }
                        }
                    }
                }
            }
        }
        terms10.push(acc10.0 + acc10.1.ln() + lnze);
        let tot = q_best.0 + e_best.0;
        if best.as_ref().map_or(true, |b| tot > b.tot) {
            best = Some(MapRec {
                tot,
                sh,
                q: q_best.1,
                ru: q_best.2,
                rd: q_best.3,
                v: q_best.4,
                e_cfg: e_best.1,
                re: e_best.2,
            });
        }
    }
    let prior_w = 5.0 * (NC as f64).ln() + 4.0 * (6.0f64).ln() + (SIG_GRID.len() as f64).ln();
    let lnz10 = lse(&terms10) - prior_w;
    check(
        "lnZ₁₀(1,1) の回帰: v16.4 公表値 −24.29 と一致 (±0.02)",
        (lnz10 - REF_LNZ10_WINNER).abs() < 0.02,
        format!("lnZ₁₀ = {:.3}", lnz10),
    );
    let b = best.unwrap();

    // ---- [2] スコアカード ----
    let v = &b.v;
    let j_map = jarlskog(v);
    let quantities: Vec<(&str, f64, f64, bool)> = vec![
        // (名前, 予測, 実測, 尤度に入れたか)
        ("m_u/m_t", b.ru[0].exp(), EPS_OBS[0], true),
        ("m_c/m_t", b.ru[1].exp(), EPS_OBS[1], true),
        ("m_d/m_b", b.rd[0].exp(), EPS_OBS[2], true),
        ("m_s/m_b", b.rd[1].exp(), EPS_OBS[3], true),
        ("m_e/m_τ", b.re[0].exp(), EPS_OBS[4], true),
        ("m_μ/m_τ", b.re[1].exp(), EPS_OBS[5], true),
        ("|V_us|", cab(v, 0, 1), EPS_OBS[6], true),
        ("|V_cb|", cab(v, 1, 2), EPS_OBS[7], true),
        ("|V_ub|", cab(v, 0, 2), EPS_OBS[8], true),
        ("|J|", j_map.abs(), j_obs, true),
        ("|V_td|", cab(v, 2, 0), vtd_obs, false),
        ("|V_ts|", cab(v, 2, 1), vts_obs, false),
    ];
    println!(
        "\n[2] MAP 構成: σ_H = {}, (kQ,ku,kd) = ({},{},{}), (σu,σd) = ({},{}), e: (kL,ke,σL,σe) = {:?}",
        b.sh, b.q[0], b.q[1], b.q[2], b.q[3], b.q[4], b.e_cfg
    );
    println!("\n    量          予測         実測         factor   区分");
    let mut hits_in = 0;
    let mut n_in = 0;
    let mut hits_out = 0;
    for (name, pred, obs, in_lik) in &quantities {
        let f = (pred / obs).max(obs / pred);
        let hit = f <= 5.0;
        if *in_lik {
            n_in += 1;
            if hit {
                hits_in += 1;
            }
        } else if hit {
            hits_out += 1;
        }
        println!(
            "    {:8}  {:11.4e}  {:11.4e}  {:6.2}   {}{}",
            name,
            pred,
            obs,
            f,
            if *in_lik { "尤度" } else { "holdout" },
            if hit { "" } else { "  ← factor 5 超え" }
        );
    }
    println!(
        "\n    factor 5 以内: 尤度 10 量中 {} / holdout 2 量中 {}",
        hits_in, hits_out
    );

    // UT 角 (導出量 — 再位相化不変)
    let arg = |z: (f64, f64)| -> f64 { z.1.atan2(z.0) };
    let mul = |a: (f64, f64), bb: (f64, f64)| (a.0 * bb.0 - a.1 * bb.1, a.0 * bb.1 + a.1 * bb.0);
    let conj = |a: (f64, f64)| (a.0, -a.1);
    let neg = |a: (f64, f64)| (-a.0, -a.1);
    // γ = arg(−V_ud V*_ub / (V_cd V*_cb)), β = arg(−V_cd V*_cb / (V_td V*_tb))
    let num_g = neg(mul(v[0][0], conj(v[0][2])));
    let den_g = mul(v[1][0], conj(v[1][2]));
    let gamma = (arg(num_g) - arg(den_g)).rem_euclid(2.0 * std::f64::consts::PI);
    let num_b = neg(mul(v[1][0], conj(v[1][2])));
    let den_b = mul(v[2][0], conj(v[2][2]));
    let beta = (arg(num_b) - arg(den_b)).rem_euclid(2.0 * std::f64::consts::PI);
    let deg = 180.0 / std::f64::consts::PI;
    let norm_ang = |a: f64| -> f64 {
        let mut d = a * deg;
        if d > 180.0 {
            d -= 360.0;
        }
        d
    };
    println!(
        "    UT 角 (導出量, 参考): γ = {:+.1}° (実測 ±{:.1}°), β = {:+.1}° (実測 ±{:.1}°)",
        norm_ang(gamma),
        gamma_obs,
        norm_ang(beta),
        beta_obs
    );
    println!("    (J の符号/角の符号は共役対称で予言されない — 大きさが対象)");

    check(
        "|V_td| MAP の回帰: v16.4 の記録 0.0082 と一致 (±0.001)",
        (cab(v, 2, 0) - 0.0082).abs() < 0.001,
        format!("|V_td|_MAP = {:.4}", cab(v, 2, 0)),
    );

    // ---- artifact ----
    let mut qrows = Vec::new();
    for (name, pred, obs, in_lik) in &quantities {
        qrows.push(Json::Obj(vec![
            ("name".into(), Json::Str(name.to_string())),
            ("pred".into(), Json::Num(*pred)),
            ("obs".into(), Json::Num(*obs)),
            ("in_likelihood".into(), Json::Bool(*in_lik)),
        ]));
    }
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.5".into())),
        ("lnz10".into(), Json::Num(lnz10)),
        ("scorecard".into(), Json::Arr(qrows)),
        ("hits_in_likelihood".into(), Json::Int(hits_in as i64)),
        ("hits_holdout".into(), Json::Int(hits_out as i64)),
        ("gamma_deg".into(), Json::Num(norm_ang(gamma))),
        ("beta_deg".into(), Json::Num(norm_ang(beta))),
    ]);
    let p = write_artifact("results/v165_scorecard.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — スコアカードは [2] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
