//! v16.3 シアー入り T²×T² の三つ巴 — 深さ × Occam × CP
//!
//! v16.2 は「複素構造 (シアー τ_re) が CP を作り、世代数 3 と両立する」ことを
//! 単一トーラスで実証した。本バイナリはそれを物理的な土俵 — T²×T² 積模型の
//! 9 量証拠 (質量比 6 + |CKM| 3, v10.1 の M_perm) — に載せ、三つ巴を裁く:
//!
//!   深さ:  シアーは lnZ₉ (質量階層 + |CKM|) を助けるか壊すか
//!   Occam: シアー (s₁, s₂) ∈ {0..3}² は 2 つの離散パラメータ (罰 2 ln 4)
//!   CP:    シアー模型の事後 J 分布は測定 J = 3.08e-5 に届くか
//!          (J は依然として尤度に入れない — holdout の規律を維持)
//!
//! 方法: 各シアー s の単一トーラス表 (12 対角化, v16.2 の構成) から、
//! (s₁, s₂) の 16 幾何それぞれで v10.1 と同一の五重和 (kQ,σu,ku,σd,kd) ×
//! e 因子 × σ_H 4 点の M_perm 証拠を計算。(0,0) は v10.1 の公表値 −19.86 を
//! 回帰再現すること (エンジンの較正)。J は同じ事後重みで分布を集計する。
//!
//! 判定 (事前登録):
//!   (a) ある (s₁,s₂)≠(0,0) の lnZ₉ が rect + 1 nat を超える → 深さも CP も勝つ窓
//!   (b) rect ± 1 nat → 深さ中立: CP を出す自由は Occam 罰 (marginal 比較) が裁く
//!   (c) rect − 1 nat 未満 → 深さと CP の緊張を nats で記録 (どちらでも測量)

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const NC: usize = 36;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const REF_NINE_PERM: f64 = -19.86; // v10.1 の公表値 (rect 回帰)
const SIG_GRID: [f64; 4] = [1.0, 1.5, 2.0, 2.5];
const SHEARS: [usize; 4] = [0, 1, 2, 3];

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
    println!("=== v16.3 シアー入り T²×T² の三つ巴: 深さ × Occam × CP ===\n");
    println!(
        "事前登録の判定: (a) lnZ₉ > rect+1 で深さも勝つ窓 / (b) rect±1 で深さ中立 (Occam が裁く) /"
    );
    println!("               (c) rect−1 未満で深さと CP の緊張を nats で記録\n");
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---- [0] シアー別の単一トーラス局在モード (12 対角化 × 4 シアー) ----
    println!("[0] シアー別モード構築 (4 シアー × 12 対角化, 安定ラベル)");
    let t0 = std::time::Instant::now();
    let mut locs_s: Vec<Vec<Vec<C3v>>> = Vec::new();
    let mut ok_deg = true;
    for &s in &SHEARS {
        let mut locs: Vec<Vec<C3v>> = Vec::new();
        for k in 0..NK12 {
            let (modes, gap, spread) = flux_modes_shear(k, s);
            if spread > 1e-9 || gap < 0.05 {
                ok_deg = false;
            }
            let (raw, cents) = localize_unsorted(&modes);
            let ord = order_stable(&cents);
            locs.push(ord.iter().map(|&i| raw[i]).collect());
        }
        locs_s.push(locs);
    }
    println!("    完了 ({} ms)", t0.elapsed().as_millis());
    check(
        "全シアーで厳密 3 重縮退・健全ギャップ (v16.2 の再確認)",
        ok_deg,
        "縮退幅 < 1e-9, ギャップ > 0.05".to_string(),
    );

    // ---- [1] 16 幾何の証拠と J ----
    println!("\n[1] (s₁, s₂) ∈ {{0..3}}² の 16 幾何 — lnZ₉ (M_perm) と事後 J");
    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let t1 = std::time::Instant::now();
    // 単一トーラス湯川表: ytabs[s][isg][ab]
    let mut ytabs: Vec<Vec<Vec<M3>>> = Vec::new();
    for si in 0..SHEARS.len() {
        let mut per_sig = Vec::new();
        for &sh in &SIG_GRID {
            let tab: Vec<M3> = (0..NK12 * NK12)
                .map(|ab| yukawa(&locs_s[si][ab % NK12], &locs_s[si][ab / NK12], sh))
                .collect();
            per_sig.push(tab);
        }
        ytabs.push(per_sig);
    }
    println!("    湯川表 ({} ms)", t1.elapsed().as_millis());

    let t2 = std::time::Instant::now();
    let mut lnz_grid = [[0.0f64; 4]; 4];
    let mut jstats: Vec<(usize, usize, f64, f64, f64)> = Vec::new(); // (s1,s2,J_med,J16,J84)
    let mut jmap_grid = [[0.0f64; 4]; 4];
    for (i1, _s1) in SHEARS.iter().enumerate() {
        for (i2, _s2) in SHEARS.iter().enumerate() {
            let mut terms = Vec::new();
            // J 事後集計 (対数ヒストグラム)
            let mut hist = vec![0.0f64; 400];
            let (mut wsum, mut shift_used) = (0.0f64, f64::NEG_INFINITY);
            let mut map_ll = f64::NEG_INFINITY;
            let mut map_j = 0.0f64;
            for isg in 0..SIG_GRID.len() {
                let yt1 = &ytabs[i1][isg];
                let yt2 = &ytabs[i2][isg];
                // ペアキャッシュ (a,b,σ) → (比, 左固有ベクトル)
                let pair: Vec<([f64; 2], M3)> = (0..NC * NC * 6)
                    .map(|m| pair_yukawa2(yt1, yt2, m % NC, (m / NC) % NC, 0, m / (NC * NC)))
                    .map(|y| mass_and_vecs(&y))
                    .collect();
                // e セクター (σL, σe 込み)
                let mut le = Vec::with_capacity(NC * NC * 36);
                for sl in 0..6 {
                    for se_ in 0..6 {
                        for ab in 0..NC * NC {
                            let r = mass_ratios(&pair_yukawa2(yt1, yt2, ab % NC, ab / NC, sl, se_));
                            le.push(ll2(&r, tgt[4], tgt[5]));
                        }
                    }
                }
                let lnze = lse(&le);
                // クォーク五重和 + J 集計
                let mut acc_q = (f64::NEG_INFINITY, 0.0f64);
                if shift_used == f64::NEG_INFINITY {
                    shift_used = -8.0; // oracle 近傍の基準シフト
                }
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
                                    if ll > acc_q.0 {
                                        acc_q.1 = acc_q.1 * (acc_q.0 - ll).exp() + 1.0;
                                        acc_q.0 = ll;
                                    } else {
                                        acc_q.1 += (ll - acc_q.0).exp();
                                    }
                                    let j = jarlskog(&v);
                                    let w = ((ll + lnze) - shift_used).exp().min(1e30);
                                    wsum += w;
                                    let x = j.abs().max(1e-300).ln().clamp(-18.0, -0.0001);
                                    let b = ((x + 18.0) / 18.0 * 400.0) as usize;
                                    hist[b.min(399)] += w;
                                    if ll > map_ll {
                                        map_ll = ll;
                                        map_j = j;
                                    }
                                }
                            }
                        }
                    }
                }
                terms.push(acc_q.0 + acc_q.1.ln() + lnze);
            }
            let lnz = lse(&terms)
                - (5.0 * (NC as f64).ln() + 4.0 * (6.0f64).ln() + (SIG_GRID.len() as f64).ln());
            lnz_grid[i1][i2] = lnz;
            jmap_grid[i1][i2] = map_j;
            let quant = |q: f64| -> f64 {
                let mut acc = 0.0;
                for (i, h) in hist.iter().enumerate() {
                    acc += h;
                    if acc >= q * wsum {
                        return (-18.0 + (i as f64 + 0.5) / 400.0 * 18.0).exp();
                    }
                }
                1.0
            };
            jstats.push((i1, i2, quant(0.5), quant(0.16), quant(0.84)));
        }
    }
    println!("    16 幾何の証拠+J ({} ms)", t2.elapsed().as_millis());

    // ---- [2] 表示と判定 ----
    println!("\n    lnZ₉ (行 = s₁, 列 = s₂):");
    print!("         ");
    for s2 in SHEARS {
        print!("  s₂={}   ", s2);
    }
    println!();
    for (i1, s1) in SHEARS.iter().enumerate() {
        print!("    s₁={} ", s1);
        for i2 in 0..4 {
            print!(" {:8.2}", lnz_grid[i1][i2]);
        }
        println!();
    }
    println!("\n    事後 |J| 中央値 [68% 帯] と MAP の J (測定 J = 3.08e-5):");
    for &(i1, i2, jm, j16, j84) in &jstats {
        if (i1 == 0 && i2 == 0) || lnz_grid[i1][i2] > lnz_grid[0][0] - 3.0 {
            println!(
                "    (s₁,s₂)=({},{}):  |J| = {:.2e} [{:.2e}, {:.2e}]  J_MAP = {:+.2e}  lnZ₉ = {:.2}",
                SHEARS[i1], SHEARS[i2], jm, j16, j84, jmap_grid[i1][i2], lnz_grid[i1][i2]
            );
        }
    }

    let rect = lnz_grid[0][0];
    check(
        "rect (0,0) の回帰: v10.1 公表値 −19.86 と一致 (±0.02)",
        (rect - REF_NINE_PERM).abs() < 0.02,
        format!("lnZ₉ = {:.3}", rect),
    );
    // 最良シアー
    let mut best = (0usize, 0usize, f64::NEG_INFINITY);
    for i1 in 0..4 {
        for i2 in 0..4 {
            if (i1, i2) != (0, 0) && lnz_grid[i1][i2] > best.2 {
                best = (i1, i2, lnz_grid[i1][i2]);
            }
        }
    }
    // marginal 模型 (16 幾何一様 = Occam 罰 ln16) vs rect
    let all: Vec<f64> = (0..16).map(|k| lnz_grid[k / 4][k % 4]).collect();
    let lnz_marg = lse(&all) - (16.0f64).ln();
    println!(
        "\n    最良シアー: (s₁,s₂)=({},{}) lnZ₉ = {:.2} (rect との差 {:+.2})",
        SHEARS[best.0],
        SHEARS[best.1],
        best.2,
        best.2 - rect
    );
    println!(
        "    marginal 模型 (16 幾何一様, Occam 罰 ln16={:.2}): lnZ₉ = {:.2} (rect との差 {:+.2})",
        (16.0f64).ln(),
        lnz_marg,
        lnz_marg - rect
    );

    println!("\n[3] 事前登録の判定:");
    let d = best.2 - rect;
    if d > 1.0 {
        println!(
            "    => (a) 深さも CP も勝つ窓が在る — (s₁,s₂)=({},{}) が rect を {:+.2} nats 上回り、",
            SHEARS[best.0], SHEARS[best.1], d
        );
        println!(
            "       かつ J ≠ 0 を作る。幾何選択 (v9.1/v10.2) 以来の主結果候補 — 次版で頑健性検査。"
        );
    } else if d > -1.0 {
        println!(
            "    => (b) 深さ中立 (差 {:+.2} nats) — シアーは質量階層をほぼ損なわずに CP を作る。",
            d
        );
        println!(
            "       Occam の裁定: marginal − rect = {:+.2} nats (罰 ln16 込み)。",
            lnz_marg - rect
        );
        println!("       J という 10 個目の観測量を説明できるのはシアー側だけである点が次版の的");
        println!("       (J を尤度に入れた 10 量証拠での再戦 — ただし holdout を失う交換になる)。");
    } else {
        println!(
            "    => (c) 深さと CP の緊張: 最良シアーでも rect に {:+.2} nats 負ける。",
            d
        );
        println!(
            "       CP を出す幾何は質量階層を壊す — 緊張の構造 (どの観測量が壊れるか) が次の的。"
        );
    }

    // ---- artifact ----
    let mut rows = Vec::new();
    for i1 in 0..4 {
        for i2 in 0..4 {
            let js = jstats.iter().find(|r| r.0 == i1 && r.1 == i2).unwrap();
            rows.push(Json::Obj(vec![
                ("s1".into(), Json::Int(SHEARS[i1] as i64)),
                ("s2".into(), Json::Int(SHEARS[i2] as i64)),
                ("lnz9".into(), Json::Num(lnz_grid[i1][i2])),
                ("j_median".into(), Json::Num(js.2)),
                ("j_16".into(), Json::Num(js.3)),
                ("j_84".into(), Json::Num(js.4)),
                ("j_map".into(), Json::Num(jmap_grid[i1][i2])),
            ]));
        }
    }
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.3".into())),
        ("rect_lnz9".into(), Json::Num(rect)),
        ("best_shear_lnz9".into(), Json::Num(best.2)),
        (
            "best_shear".into(),
            Json::Arr(vec![
                Json::Int(SHEARS[best.0] as i64),
                Json::Int(SHEARS[best.1] as i64),
            ]),
        ),
        ("marginal_lnz9".into(), Json::Num(lnz_marg)),
        ("grid".into(), Json::Arr(rows)),
        ("j_measured_pdg".into(), Json::Num(3.08e-5)),
    ]);
    let p = write_artifact("results/v163_sheartriple.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 三つ巴の判定は [3] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
