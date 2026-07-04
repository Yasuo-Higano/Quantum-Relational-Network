//! v11.3 場依存の磁化の向き — 混合パリティの物理的実現を検定する
//!
//! v11.2 の結論: データは「セクターごとに射影の向きが違う」ことを要求する。最有力の
//! 物理的実現は**場ごとに第 2 トーラスの磁束の向き (磁化の符号) が違う**構成である
//! (磁化 brane 模型で普通の状況)。磁束 −Q のゼロモードは +Q の複素共役なので、
//!
//!   M_orient: 場 F の第 2 トーラス因子 = (巡回オフセット c_F ∈ Z₃) × (向き ε_F ∈ Z₂)
//!             — 1 場 6 状態、事前 4 ln 6
//!
//! これは v10.1 の M_full (S₃ 対、1 場 6 状態、事前 4 ln 6) と**同じ事前での機構対決**。
//! v10.3 は「共役 ≠ 奇置換対 (指紋差 ~1e-3)」を示したので答えは自明でない。
//!
//! 構造の事実 (本バイナリで数値検証もする): svals は複素共役に不変なので
//! (ε_F, ε_G) には**相対の向きだけ**が効く — (C,C)≡(N,N), (C,N)≡(N,C)。
//! ゆえにゲージ固定 ε_Q = N, c_Q = 0 は厳密。必要な単一トーラス表は
//!   A[i][j] = Σ conj(ψ_i) φ ψ_j (通常),  B[i][j] = Σ ψ_i φ ψ_j (双線形; 相対 C 用)
//! の 2 つだけ: (N,C) = conj(B) だが svals には B で十分。
//!
//! 退化検査: (a) ε 全て N に制限 = v11.2 の M_uni (-21.8275) と一致,
//!           (b) 相対性の恒等式 svals(C,C)=svals(N,N) 等 (機械精度)。
//! 読み出し: MAP の (c, ε) パターン — 磁化符号の選択則の候補。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];

// v9.2 の一次ソース (results/v92_labelstab.json, 安定ラベル) — 退化検査の目標値
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

/// 双線形の単一トーラス表 B[i][j] = Σ_s ψ_i(s) φ_H(s) ψ_j(s) (共役なし)
fn yukawa_bilinear(la: &[C3v], lb: &[C3v], sig_h: f64) -> M3 {
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
                // conj なしの双線形
                sr += (ar * br - ai * bi) * phih[s];
                si += (ar * bi + ai * br) * phih[s];
            }
            y_out[i][j] = (sr, si);
        }
    }
    y_out
}

const EVEN: [usize; 3] = [0, 3, 4]; // e, (123), (132) — 巡回 (偶置換)

/// 向きつきペア湯川。行場 (cr, er) / 列場 (cg, eg)。表は εr⊕εg で選ぶ:
/// 相対 N → A (通常), 相対 C → B (双線形)。svals/|CKM| は共役に不変なので十分。
fn pair_y_orient(
    ta: &[M3],
    tb: &[M3],
    a: usize,
    b: usize,
    cr: usize,
    er: usize,
    cg: usize,
    eg: usize,
) -> M3 {
    let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
    let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
    let y1 = &ta[a1 + b1 * NK12];
    let t2 = if er == eg { &ta[a2 + b2 * NK12] } else { &tb[a2 + b2 * NK12] };
    let (pr, pg) = (&PERMS[EVEN[cr]], &PERMS[EVEN[cg]]);
    let mut y = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (p, q) = y1[i][j];
            let (r, s) = t2[pr[i]][pg[j]];
            y[i][j] = (p * r - q * s, p * s + q * r);
        }
    }
    y
}

fn main() {
    self_test();
    println!("=== v11.3 場依存の磁化の向き: 混合パリティの物理的実現の検定 ===\n");
    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();

    println!("[0] 世代モード (対角化 12 回, 安定ラベル)");
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
    println!("    縮退・ギャップ不変  {}  ({} ms)", pass(ok_engine), t0.elapsed().as_millis());

    // ---- [1] 相対性の恒等式 (svals は共役に不変) ----
    println!("\n[1] 恒等式: 向きは相対だけが効く (svals(C,C) = svals(N,N) を機械精度で)");
    // (C,C) 表 = conj(A): svals 同一は自明だが、実装経路 (共役モードで直接構成) と
    // 表経路 (A のまま) の一致を数点で確認する
    let mut max_dev: f64 = 0.0;
    {
        let sh = 1.0;
        let mk_conj = |v: &Vec<C3v>| -> Vec<C3v> {
            v.iter()
                .map(|m| {
                    let mut z = *m;
                    for p in z.iter_mut() {
                        p.1 = -p.1;
                    }
                    z
                })
                .collect()
        };
        for (ka, kb) in [(0usize, 2usize), (4, 8), (6, 10)] {
            let a_tab = yukawa(&locs[ka], &locs[kb], sh);
            let cc = yukawa(&mk_conj(&locs[ka]), &mk_conj(&locs[kb]), sh);
            let (ra, _) = mass_and_vecs(&a_tab);
            let (rc, _) = mass_and_vecs(&cc);
            max_dev = max_dev.max((ra[0] - rc[0]).abs()).max((ra[1] - rc[1]).abs());
            // (C,N) 直接 vs 双線形 B
            let cn = yukawa(&mk_conj(&locs[ka]), &locs[kb], sh);
            let bt = yukawa_bilinear(&locs[ka], &locs[kb], sh);
            let (r1, _) = mass_and_vecs(&cn);
            let (r2, _) = mass_and_vecs(&bt);
            max_dev = max_dev.max((r1[0] - r2[0]).abs()).max((r1[1] - r2[1]).abs());
        }
    }
    let ok_rel = max_dev < 1e-12;
    println!("    max|Δ質量対数比| = {:.2e} (< 1e-12)  {}", max_dev, pass(ok_rel));

    // ---- [2] 機構対決: M_orient vs M_full (同じ事前 4ln6) ----
    println!("\n[2] 機構対決 (T²×Z₆, ゲージ ε_Q=N, c_Q=0)");
    let nc = 36usize;
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let mut terms_or = Vec::new();
    let mut terms_n = Vec::new(); // ε 全 N 制限 (= v11.2 M_uni) — 退化検査
    let mut map_or = (f64::NEG_INFINITY, [9usize; 2], [9usize; 2]); // (cu,εu),(cd,εd)
    let mut map_le = (f64::NEG_INFINITY, [9usize; 2], [9usize; 2]); // (cL,εL),(ce,εe)
    for &sh in &sig_grid {
        let ta: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        let tb: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa_bilinear(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // クォーク列状態 m = c + 3ε (行 Q は N, c=0)
        let pair: Vec<([f64; 2], M3)> = (0..nc * nc * 6)
            .map(|q| {
                let (ab, m) = (q % (nc * nc), q / (nc * nc));
                mass_and_vecs(&pair_y_orient(&ta, &tb, ab % nc, ab / nc, 0, 0, m % 3, m / 3))
            })
            .collect();
        // e セクター: 行 L (mL) × 列 e (me) — 36 状態組
        let mut le_all = Vec::with_capacity(nc * nc * 36);
        let mut le_n = Vec::with_capacity(nc * nc * 9);
        let mut le_best = f64::NEG_INFINITY;
        for ml in 0..6usize {
            for me in 0..6usize {
                for ab in 0..nc * nc {
                    let r = mass_ratios(&pair_y_orient(
                        &ta,
                        &tb,
                        ab % nc,
                        ab / nc,
                        ml % 3,
                        ml / 3,
                        me % 3,
                        me / 3,
                    ));
                    let l = ll2(&r, tgt[4], tgt[5]);
                    le_all.push(l);
                    if ml / 3 == 0 && me / 3 == 0 {
                        le_n.push(l);
                    }
                    if l > le_best {
                        le_best = l;
                        map_le = (l, [ml % 3, ml / 3], [me % 3, me / 3]);
                    }
                }
            }
        }
        // クォーク部
        let mut per_q_or = Vec::with_capacity(nc);
        let mut per_q_n = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut ao = (f64::NEG_INFINITY, 0.0f64);
            let mut an = (f64::NEG_INFINITY, 0.0f64);
            for mu in 0..6usize {
                for ku in 0..nc {
                    let (ru, vu) = &pair[kq + ku * nc + mu * nc * nc];
                    let llu = ll2(ru, tgt[0], tgt[1]);
                    for md in 0..6usize {
                        for kd in 0..nc {
                            let (rd, vd) = &pair[kq + kd * nc + md * nc * nc];
                            let lld = ll2(rd, tgt[2], tgt[3]);
                            let c = ckm3(vu, vd);
                            let mut ll = llu + lld;
                            for m in 0..3 {
                                let d = c[m].max(1e-300).ln() - tgt[6 + m];
                                ll += -d * d / (2.0 * sigma * sigma) + norm1;
                            }
                            if ll > ao.0 {
                                ao.1 = ao.1 * (ao.0 - ll).exp() + 1.0;
                                ao.0 = ll;
                            } else {
                                ao.1 += (ll - ao.0).exp();
                            }
                            if mu / 3 == 0 && md / 3 == 0 {
                                if ll > an.0 {
                                    an.1 = an.1 * (an.0 - ll).exp() + 1.0;
                                    an.0 = ll;
                                } else {
                                    an.1 += (ll - an.0).exp();
                                }
                            }
                            let tot = ll + le_best;
                            if tot > map_or.0 {
                                map_or = (tot, [mu % 3, mu / 3], [md % 3, md / 3]);
                            }
                        }
                    }
                }
            }
            per_q_or.push(ao.0 + ao.1.ln());
            per_q_n.push(an.0 + an.1.ln());
        }
        terms_or.push(lse(&per_q_or) + lse(&le_all));
        terms_n.push(lse(&per_q_n) + lse(&le_n));
    }
    let prior_w = 10.0 * (6.0f64).ln() + (sig_grid.len() as f64).ln();
    let lnz_or = lse(&terms_or) - prior_w - 4.0 * (6.0f64).ln();
    let lnz_n = lse(&terms_n) - prior_w - 4.0 * (3.0f64).ln();
    let ref_uni = -21.82754; // v11.2 M_uni (results/v112_orbifold.json, 表示桁)
    let ref_full = -19.86334559888438; // v10.1 M_perm
    let ok_uni = (lnz_n - ref_uni).abs() < 0.01;
    println!("    退化検査: ε 全 N = v11.2 M_uni: {:.4} vs {:.4}  {}", lnz_n, ref_uni, pass(ok_uni));
    println!("    lnZ₉(M_orient) = {:.4}  (事前 4ln6 込み)", lnz_or);
    println!("    lnZ₉(M_full)   = {:.4}  (v10.1, S₃ 対)", ref_full);
    let verdict = if lnz_or > ref_full + 0.02 {
        "M_orient が勝つ — 混合パリティの正体は磁化の向きで実現される"
    } else if lnz_or > ref_full - 0.5 {
        "ほぼ互角 — 向き模型は S₃ 対の物理を実質的に実現する"
    } else {
        "M_full が勝つ — 共役 (向き) では S₃ 対の代わりにならない"
    };
    println!("\n    => {} (差 {:+.2})", verdict, lnz_or - ref_full);
    println!(
        "    MAP: u=(c{},{}) d=(c{},{}) / L=(c{},{}) e=(c{},{})",
        map_or.1[0],
        if map_or.1[1] == 0 { "N" } else { "C" },
        map_or.2[0],
        if map_or.2[1] == 0 { "N" } else { "C" },
        map_le.1[0],
        if map_le.1[1] == 0 { "N" } else { "C" },
        map_le.2[0],
        if map_le.2[1] == 0 { "N" } else { "C" }
    );

    let all_ok = ok_engine && ok_rel && ok_uni;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-011".into())),
        ("relativity_identity_dev".into(), Json::Num(max_dev)),
        ("lnZ_nine_orient".into(), Json::Num(lnz_or)),
        ("lnZ_nine_full_ref".into(), Json::Num(ref_full)),
        ("lnZ_nine_allN".into(), Json::Num(lnz_n)),
        (
            "map_orient".into(),
            Json::Arr(vec![
                Json::Int(map_or.1[0] as i64),
                Json::Int(map_or.1[1] as i64),
                Json::Int(map_or.2[0] as i64),
                Json::Int(map_or.2[1] as i64),
                Json::Int(map_le.1[0] as i64),
                Json::Int(map_le.1[1] as i64),
                Json::Int(map_le.2[0] as i64),
                Json::Int(map_le.2[1] as i64),
            ]),
        ),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v113_orient.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {} (PASS = 装置検証 — 機構対決は [2] が本体)", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
