//! v10.1 対角世代対の導出 — 仮定を選択結果に変える (v10.0 残高 2 の筆頭)
//!
//! T²×T² 積模型 (v7.2) は「世代 i はトーラス 1 のモード i とトーラス 2 のモード i の積」
//! という**対角世代対**を宣言された仮定として置いてきた。v9.2 はこの対の組み方が
//! lnZ を 1〜4 nats 動かす「模型の内容」であることを暴いた。本バイナリは対の組み方を
//! 仮定から外し、データに選ばせる:
//!
//!   場 F の世代対は置換 σ_F ∈ S₃ (トーラス 1 のモード i ↔ トーラス 2 のモード σ_F(i))。
//!   トーラス 2 のラベル付け替え τ (大域) で σ_F → τ∘σ_F (全場一斉) なので、
//!   物理的な対空間は S₃⁵/S₃ — ゲージ固定 σ_Q = e で (σ_u, σ_d, σ_L, σ_e) ∈ S₃⁴ が残る。
//!   湯川は Y[i][j] = Y¹[i][j] · Y²[σ_F(i)][σ_G(j)] (F=左場, G=右場)。
//!
//! 比較する模型 (世代ラベルは v9.2 の安定規約):
//!   M_diag: 全 σ = e (従来の対角対)。パラメータは Wilson 線 10 + σ_H のまま。
//!   M_perm: (σ_u,σ_d,σ_L,σ_e) を一様事前で marginalize (Occam 罰 4 ln 6 ≈ 7.17)。
//!
//! 問い:
//!   [A] lnZ(M_diag) ≥ lnZ(M_perm) か? — Yes なら「対角対はデータが Occam 選択する」
//!       となり、仮定は事後的に正当化される (導出の第一歩)。
//!   [B] MAP の対は e (= 安定ラベルの中心整列) か? — Yes なら「幾何 (どのモードが
//!       同じ住所に住むか) が対を決める」という規則がデータと整合する。
//!
//! 検証: σ 全て e の退化検査が v9.2 の安定側の値 (質量のみ -18.76 / 全 9 量 -23.61,
//! results/v92_labelstab.json) と厳密一致すること。対付き表の e 成分が無置換の
//! Hadamard 積と一致すること (実装の内部整合)。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];

// v9.2 の一次ソース (results/v92_labelstab.json, 安定ラベル) — 退化検査の目標値
const REF_MASS_DIAG: f64 = -18.76; // lnZ_mass.T2xZ6.stable_label (表示 2 桁; 判定は ±0.02)
const REF_NINE_DIAG: f64 = -23.61; // lnZ_nine.T2xZ6.stable_label

const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2], // e (恒等 = 対角対 = 安定ラベルの中心整列)
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

/// 安定ラベル (v9.2): 0.5 サイト格子にスナップしてから昇順
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

/// 対付き Hadamard 積: Y[i][j] = Y¹[i][j] · Y²[σf(i)][σg(j)]
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
    m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// T²×Z₆ の Wilson 複合添字 (a1,a2) と単一トーラス湯川表から、対 (σf,σg) 付きの
/// 積湯川を作る。ytab は 12×12 (半ステップ添字; Z₆ は偶数のみ使う)。
fn pair_yukawa(ytab: &[M3], a: usize, b: usize, sf: usize, sg: usize) -> M3 {
    let (a1, a2) = (2 * (a % 6), 2 * (a / 6)); // Z₆ → 半ステップ添字
    let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
    let y1 = &ytab[a1 + b1 * NK12];
    let y2 = &ytab[a2 + b2 * NK12];
    had_prod_perm(y1, y2, sf, sg)
}

fn main() {
    self_test();
    println!("=== v10.1 対角世代対の導出: 置換 marginalize と幾何的整列の検定 ===\n");
    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();

    // ---- [0] 世代モード (安定ラベル) ----
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

    // ---- 対付きセクター表の構築 (σ_H ごと) ----
    // lu[(kQ,ku,σu)], ld[(kQ,kd,σd)]: 36×36×6 / le[(kL,ke,σL,σe)]: 36×36×36
    let nc = 36usize;
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };

    // 内部整合検査: 対 e の pair_yukawa が無置換 Hadamard と bit 一致
    {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], 1.0))
            .collect();
        let y_pair = pair_yukawa(&ytab, 7, 13, 0, 0);
        let (a1, a2, b1, b2) = (2 * (7 % 6), 2 * (7 / 6), 2 * (13 % 6), 2 * (13 / 6));
        let y1 = &ytab[a1 + b1 * NK12];
        let y2 = &ytab[a2 + b2 * NK12];
        let mut ok_id = true;
        for i in 0..3 {
            for j in 0..3 {
                let (p, q) = y1[i][j];
                let (r, s) = y2[i][j];
                let want = (p * r - q * s, p * s + q * r);
                if (y_pair[i][j].0 - want.0).abs() > 0.0 || (y_pair[i][j].1 - want.1).abs() > 0.0 {
                    ok_id = false;
                }
            }
        }
        println!("    対 e = 無置換 Hadamard (bit 一致)  {}", pass(ok_id));
        if !ok_id {
            std::process::exit(1);
        }
    }

    // ---- [1] 質量のみ (6 比) の証拠: M_diag vs M_perm ----
    println!("\n[1] 質量のみ (6 比) の証拠");
    let t1 = std::time::Instant::now();
    let mut terms_diag = Vec::new();
    let mut terms_perm = Vec::new();
    let mut map_mass = (f64::NEG_INFINITY, [0usize; 4]);
    for &sh in &sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // セクター対数尤度表
        let mut lu = vec![f64::NEG_INFINITY; nc * nc * 6];
        let mut ld = vec![f64::NEG_INFINITY; nc * nc * 6];
        let mut le = vec![f64::NEG_INFINITY; nc * nc * 36];
        for a in 0..nc {
            for b in 0..nc {
                for s in 0..6 {
                    let r = mass_ratios(&pair_yukawa(&ytab, a, b, 0, s));
                    lu[a + b * nc + s * nc * nc] = ll2(&r, tgt[0], tgt[1]);
                    ld[a + b * nc + s * nc * nc] = ll2(&r, tgt[2], tgt[3]);
                }
                for sl in 0..6 {
                    for se in 0..6 {
                        let r = mass_ratios(&pair_yukawa(&ytab, a, b, sl, se));
                        le[a + b * nc + (sl + se * 6) * nc * nc] = ll2(&r, tgt[4], tgt[5]);
                    }
                }
            }
        }
        // M_diag (σ 全て e)
        let mut per_q = Vec::with_capacity(nc);
        for kq in 0..nc {
            let us: Vec<f64> = (0..nc).map(|ku| lu[kq + ku * nc]).collect();
            let ds: Vec<f64> = (0..nc).map(|kd| ld[kq + kd * nc]).collect();
            per_q.push(lse(&us) + lse(&ds));
        }
        let le_diag: Vec<f64> = (0..nc * nc).map(|ab| le[ab]).collect();
        terms_diag.push(lse(&per_q) + lse(&le_diag));
        // M_perm ((σu,σd,σL,σe) を marginalize)
        let mut per_q2 = Vec::with_capacity(nc);
        for kq in 0..nc {
            let us: Vec<f64> = (0..nc * 6)
                .map(|m| lu[kq + (m % nc) * nc + (m / nc) * nc * nc])
                .collect();
            let ds: Vec<f64> = (0..nc * 6)
                .map(|m| ld[kq + (m % nc) * nc + (m / nc) * nc * nc])
                .collect();
            per_q2.push(lse(&us) + lse(&ds));
        }
        terms_perm.push(lse(&per_q2) + lse(&le));
        // MAP σ (質量のみ; 参考) — kQ は u/d で共有するので kQ ごとに最大化してから結合
        for su in 0..6 {
            for sd in 0..6 {
                let mut bq = f64::NEG_INFINITY;
                for kq in 0..nc {
                    let mut bu = f64::NEG_INFINITY;
                    let mut bd = f64::NEG_INFINITY;
                    for k2 in 0..nc {
                        bu = bu.max(lu[kq + k2 * nc + su * nc * nc]);
                        bd = bd.max(ld[kq + k2 * nc + sd * nc * nc]);
                    }
                    bq = bq.max(bu + bd);
                }
                for sl in 0..6 {
                    for se_ in 0..6 {
                        let mut be = f64::NEG_INFINITY;
                        for ab in 0..nc * nc {
                            be = be.max(le[ab + (sl + se_ * 6) * nc * nc]);
                        }
                        let tot = bq + be;
                        if tot > map_mass.0 {
                            map_mass = (tot, [su, sd, sl, se_]);
                        }
                    }
                }
            }
        }
    }
    let prior_w = 10.0 * (6.0f64).ln() + (sig_grid.len() as f64).ln();
    let lnz_mass_diag = lse(&terms_diag) - prior_w;
    let lnz_mass_perm = lse(&terms_perm) - prior_w - 4.0 * (6.0f64).ln();
    println!(
        "    lnZ(M_diag) = {:.2}  vs v9.2 安定側 {:.2}  {}",
        lnz_mass_diag,
        REF_MASS_DIAG,
        pass((lnz_mass_diag - REF_MASS_DIAG).abs() < 0.02)
    );
    println!(
        "    lnZ(M_perm) = {:.2}  (Occam 罰 4ln6 = {:.2} 込み)",
        lnz_mass_perm,
        4.0 * (6.0f64).ln()
    );
    println!(
        "    => 質量のみ: {}  (差 {:+.2})   [MAP σ (質量のみ, 参考) = {:?}]",
        if lnz_mass_diag >= lnz_mass_perm {
            "M_diag (対角対) が勝つ"
        } else {
            "M_perm (自由対) が勝つ"
        },
        lnz_mass_diag - lnz_mass_perm,
        map_mass.1,
    );
    println!("    ({} ms)", t1.elapsed().as_millis());
    let ok_mass_reg = (lnz_mass_diag - REF_MASS_DIAG).abs() < 0.02;

    // ---- [2] 全 9 量 (質量 6 + CKM 3) の証拠: M_diag vs M_perm ----
    println!("\n[2] 全 9 量 (質量 6 + CKM 3) の証拠 — クォーク部は (K_Q,K_u,K_d,σ_u,σ_d) の五重和");
    let t2 = std::time::Instant::now();
    let mut terms9_diag = Vec::new();
    let mut terms9_perm = Vec::new();
    let mut map9 = (f64::NEG_INFINITY, 0.0f64, [0usize; 3], [0usize; 2], [0usize; 2]);
    for &sh in &sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // ペアキャッシュ: (a,b,σ) → (質量対数比, 左固有ベクトル)
        let pair_u: Vec<([f64; 2], M3)> = (0..nc * nc * 6)
            .map(|m| mass_and_vecs(&pair_yukawa(&ytab, m % nc, (m / nc) % nc, 0, m / (nc * nc))))
            .collect();
        // e セクター (σL, σe 込みで因子化)
        let mut le_diag = Vec::with_capacity(nc * nc);
        let mut le_all = Vec::with_capacity(nc * nc * 36);
        let mut le_best = (f64::NEG_INFINITY, [0usize; 2]);
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..nc * nc {
                    let r = mass_ratios(&pair_yukawa(&ytab, ab % nc, ab / nc, sl, se_));
                    let l = ll2(&r, tgt[4], tgt[5]);
                    le_all.push(l);
                    if sl == 0 && se_ == 0 {
                        le_diag.push(l);
                    }
                    if l > le_best.0 {
                        le_best = (l, [sl, se_]);
                    }
                }
            }
        }
        // クォーク部
        let mut per_q_diag = Vec::with_capacity(nc);
        let mut per_q_perm = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut acc_d = (f64::NEG_INFINITY, 0.0f64);
            let mut acc_p = (f64::NEG_INFINITY, 0.0f64);
            for su in 0..6 {
                for ku in 0..nc {
                    let (ru, vu) = &pair_u[kq + ku * nc + su * nc * nc];
                    let llu = ll2(ru, tgt[0], tgt[1]);
                    for sd in 0..6 {
                        for kd in 0..nc {
                            let (rd, vd) = &pair_u[kq + kd * nc + sd * nc * nc];
                            let lld = ll2(rd, tgt[2], tgt[3]);
                            let c = ckm3(vu, vd);
                            let mut ll = llu + lld;
                            for m in 0..3 {
                                let d = c[m].max(1e-300).ln() - tgt[6 + m];
                                ll += -d * d / (2.0 * sigma * sigma) + norm1;
                            }
                            // ストリーミング lse (perm は常に, diag は σ=e のみ)
                            if ll > acc_p.0 {
                                acc_p.1 = acc_p.1 * (acc_p.0 - ll).exp() + 1.0;
                                acc_p.0 = ll;
                            } else {
                                acc_p.1 += (ll - acc_p.0).exp();
                            }
                            if su == 0 && sd == 0 {
                                if ll > acc_d.0 {
                                    acc_d.1 = acc_d.1 * (acc_d.0 - ll).exp() + 1.0;
                                    acc_d.0 = ll;
                                } else {
                                    acc_d.1 += (ll - acc_d.0).exp();
                                }
                            }
                            let tot = ll + le_best.0;
                            if tot > map9.0 {
                                map9 = (tot, sh, [kq, ku, kd], [su, sd], le_best.1);
                            }
                        }
                    }
                }
            }
            per_q_diag.push(acc_d.0 + acc_d.1.ln());
            per_q_perm.push(acc_p.0 + acc_p.1.ln());
        }
        terms9_diag.push(lse(&per_q_diag) + lse(&le_diag));
        terms9_perm.push(lse(&per_q_perm) + lse(&le_all));
    }
    let lnz9_diag = lse(&terms9_diag) - prior_w;
    let lnz9_perm = lse(&terms9_perm) - prior_w - 4.0 * (6.0f64).ln();
    println!(
        "    lnZ₉(M_diag) = {:.2}  vs v9.2 安定側 {:.2}  {}",
        lnz9_diag,
        REF_NINE_DIAG,
        pass((lnz9_diag - REF_NINE_DIAG).abs() < 0.02)
    );
    println!(
        "    lnZ₉(M_perm) = {:.2}  (Occam 罰 4ln6 込み)",
        lnz9_perm
    );
    let sig_names = ["e", "(23)", "(12)", "(123)", "(132)", "(13)"];
    println!(
        "    => 全 9 量: {}  (差 {:+.2})",
        if lnz9_diag >= lnz9_perm {
            "M_diag (対角対) が勝つ"
        } else {
            "M_perm (自由対) が勝つ"
        },
        lnz9_diag - lnz9_perm
    );
    println!(
        "    MAP の対: σ_u={} σ_d={} σ_L={} σ_e={}  (σ_H={}, lnL={:.2})   ({} ms)",
        sig_names[map9.3[0]],
        sig_names[map9.3[1]],
        sig_names[map9.4[0]],
        sig_names[map9.4[1]],
        map9.1,
        map9.0,
        t2.elapsed().as_millis()
    );
    let ok_nine_reg = (lnz9_diag - REF_NINE_DIAG).abs() < 0.02;
    let map_is_diag = map9.3 == [0, 0] && map9.4 == [0, 0];
    // 物理の答え ([PASS]/[FAIL] は装置検査専用トークンなので使わない)
    println!(
        "    => [B] MAP の対は中心整列 (e) か: {}",
        if map_is_diag { "Yes" } else { "No" }
    );

    // ---- [3] 判定 ----
    println!("\n[3] 判定");
    let a_mass = lnz_mass_diag >= lnz_mass_perm;
    let a_nine = lnz9_diag >= lnz9_perm;
    println!(
        "    [A] 対角対は Occam 選択されるか: 質量のみ {} / 全 9 量 {}",
        if a_mass { "Yes" } else { "No" },
        if a_nine { "Yes" } else { "No" }
    );
    println!("    [B] MAP の対は中心整列 (e) か:  {}", if map_is_diag { "Yes" } else { "No" });

    let all_ok = ok_engine && ok_mass_reg && ok_nine_reg;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-007".into())),
        ("gauge".into(), Json::Str("sigma_Q = e (torus-2 relabeling fixed)".into())),
        ("lnZ_mass_diag".into(), Json::Num(lnz_mass_diag)),
        ("lnZ_mass_perm".into(), Json::Num(lnz_mass_perm)),
        ("lnZ_nine_diag".into(), Json::Num(lnz9_diag)),
        ("lnZ_nine_perm".into(), Json::Num(lnz9_perm)),
        (
            "map_sigma_nine".into(),
            Json::Arr(vec![
                Json::Int(map9.3[0] as i64),
                Json::Int(map9.3[1] as i64),
                Json::Int(map9.4[0] as i64),
                Json::Int(map9.4[1] as i64),
            ]),
        ),
        ("map_is_diagonal".into(), Json::Bool(map_is_diag)),
        ("diag_wins_mass".into(), Json::Bool(a_mass)),
        ("diag_wins_nine".into(), Json::Bool(a_nine)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v101_pairing.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {} (PASS 条件は装置の検証 — [A][B] の答えは上の表が本体)", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
