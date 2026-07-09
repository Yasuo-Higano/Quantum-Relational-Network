//! v16.11 選択原理の再戦 — CP 制約つき 10 量・二層事前 (幾何 × config)
//!
//! v15.5 の有限候補戦 (rect N=18・9 量) は生存 0 だった: Robustness/Depth 棄却、
//! MDL 未決 (−0.33)、Thermo が唯一の正符号 (+0.32)。その後 v16.10 で幾何選択は
//! 証拠限界に到達した — 谷底 4 幾何は 10 量証拠で縮退し、「対角 vs 非対称」は
//! 決められない。本バイナリは両者を束ねる再戦である:
//!
//!   「データを見ない原理は、証拠が縮退させた幾何 (と config) を破れるか」
//!
//! 変更点 (v15.5 → 本版):
//!   ・窓: rect N=18 → シアー族 N=36 全 21 幾何 (s₁≤s₂ ∈ {0..5}²)
//!   ・証拠: 9 量 → 10 量 (J 込み — rect は尤度が構造零で殺す。これが「CP 制約」
//!     の入り方であり、原理への注入ではなく観測の側に置く)
//!   ・事前: 二層 π_P(g, c) = π_P(g)·π_P(c|g) — 幾何レベルは MDL のみ明示スコア
//!     (シアー整数の符号長)、他原理は config スコアが幾何ごとに違うことで
//!     幾何事後を誘導する
//!
//! 事前登録 (結果より先に固定):
//!   ・生存 ΔlnZ ≥ +1 / 棄却 ≤ −1 / 中間 = 未決 (v15.5 と同一)
//!   ・谷底の裁定: 一様で ≤0.5 nat だった谷底 4 幾何の対差が、原理の下で
//!     ≥1 nat 開けば「原理は縮退を破る」(向きも記録)
//!   ・グリッドは G0={2,3,4,5} — 結論が grid 条件付きであることは v16.10 で確立済み
//!
//! 装置ゲート: 一様の per-geometry 値が v16.9 の 21 値を再現 (±0.02)・一様の
//! 二層合計が v16.9 の marginal −23.121 を再現・厳密縮退 6 シアー。
//! モード表は lib.rs の disk キャッシュ (v16.11 で導入) を初めて使う — キャッシュは
//! 加速器であり一次ソースではない (削除すれば同一値が再計算される)。

use uft_sim::*;

const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
/// v16.9 の 21 幾何アンカー (s₁, s₂, lnZ₁₀) と marginal
const REF_V169: [(usize, usize, f64); 21] = [
    (0, 0, -269.103744),
    (0, 1, -24.360611),
    (0, 2, -22.608812),
    (0, 3, -23.523873),
    (0, 4, -25.032977),
    (0, 5, -25.784773),
    (1, 1, -24.520147),
    (1, 2, -22.004840),
    (1, 3, -21.982785),
    (1, 4, -23.928889),
    (1, 5, -24.391199),
    (2, 2, -22.263460),
    (2, 3, -22.256569),
    (2, 4, -22.570269),
    (2, 5, -25.284032),
    (3, 3, -21.756581),
    (3, 4, -24.530857),
    (3, 5, -25.270758),
    (4, 4, -25.893218),
    (4, 5, -26.160006),
    (5, 5, -26.475224),
];
const REF_MARGINAL: f64 = -23.121;
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
type Mode = Vec<(f64, f64)>; // 長さ N²

/// シアー s つき磁束トーラス (格子 n, 2 成分 Dirac 型) の最低 3 モード。
/// 戻り値: (モード 3 本, ギャップ, 縮退幅)
fn flux_modes_shear_n(n: usize, k_half: usize, s: usize) -> (Vec<Mode>, f64, f64) {
    let ns = n * n;
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
    let idx = |x: usize, y: usize| x + y * n;
    for x in 0..n {
        for y in 0..n {
            let th_y = phi * x as f64 + wl;
            if y == n - 1 {
                addhop(&mut a, idx(x, y), idx((x + s) % n, 0), th_y, m, ns);
            } else {
                addhop(&mut a, idx(x, y), idx(x, y + 1), th_y, m, ns);
            }
            let th_x = if x == n - 1 {
                -phi * (n as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % n, y), th_x, m, ns);
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

/// 局在化 + 安定ラベル (n パラメータ版)
fn localize_stable(n: usize, modes: &[Mode]) -> Vec<Mode> {
    let ns = n * n;
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut ure = [[0.0f64; 3]; 3];
    let mut uim = [[0.0f64; 3]; 3];
    for a in 0..Q {
        for b in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for i in 0..ns {
                let x = (i % n) as f64;
                let (sn, cs) = (two_pi * x / n as f64).sin_cos();
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
            let x = (i % n) as f64;
            let (sn, cs) = (two_pi * x / n as f64).sin_cos();
            zr += p * cs;
            zi += p * sn;
        }
        let center = (zi.atan2(zr) / two_pi * n as f64).rem_euclid(n as f64);
        out.push(psi);
        centers.push(center);
    }
    // 安定ラベル (0.5 サイト格子スナップ後ソート)
    let snapped: Vec<f64> = centers
        .iter()
        .map(|&c| ((2.0 * c).round() / 2.0).rem_euclid(n as f64))
        .collect();
    let mut ord: Vec<usize> = (0..Q).collect();
    ord.sort_by(|&a, &b| snapped[a].partial_cmp(&snapped[b]).unwrap());
    ord.iter().map(|&i| out[i].clone()).collect()
}

fn yukawa_n(n: usize, la: &[Mode], lb: &[Mode], sig_h: f64) -> M3 {
    let ns = n * n;
    let mut phih = vec![0.0f64; ns];
    for y in 0..n {
        for x in 0..n {
            let dx = (x as f64).min(n as f64 - x as f64);
            let dy = (y as f64).min(n as f64 - y as f64);
            phih[x + y * n] = (-(dx * dx + dy * dy) / (2.0 * sig_h * sig_h)).exp();
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

/// 流し込み log-sum-exp 蓄積器 (v15.5 から移植)
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

/// Wilson 成分空間 (Z6 × Z6) の ±1 近傍 4 点 (v15.5 から移植)
fn comp_neighbors(a: usize) -> [usize; 4] {
    let (a1, a2) = (a % 6, a / 6);
    [
        (a1 + 1) % 6 + a2 * 6,
        (a1 + 5) % 6 + a2 * 6,
        a1 + ((a2 + 1) % 6) * 6,
        a1 + ((a2 + 5) % 6) * 6,
    ]
}

fn main() {
    self_test();
    println!("=== v16.11 選択原理の再戦: CP 制約つき 10 量・二層事前 (N=36, 21 幾何) ===\n");
    println!("事前登録の判定規準 (結果より先に固定):");
    println!("  生存: ΔlnZ ≥ +1.0 nat / 棄却: ΔlnZ ≤ −1.0 nat / 中間: 未決 (v15.5 と同一)");
    println!(
        "  谷底の裁定: 谷底 4 幾何 (2,2)(2,3)(1,3)(3,3) の対差が原理下で ≥1 nat 開けば縮退を破る"
    );
    println!("  グリッド G0={{2,3,4,5}} 条件付き (v16.10)\n");
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

    let n = 36usize;
    let sig_grid = [2.0f64, 3.0, 4.0, 5.0];
    let par = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(12);
    let j_obs: f64 = 3.08e-5;

    // ---- [0] モード表 (disk キャッシュ → 無い分だけ対角化) ----
    let t0 = std::time::Instant::now();
    let mut raw: std::collections::BTreeMap<(usize, usize), (Vec<Mode>, f64, f64)> =
        std::collections::BTreeMap::new();
    let mut misses: Vec<(usize, usize)> = Vec::new();
    for s in 0..6usize {
        for k in 0..NK12 {
            match cache_load_modes(MODE_TAG, n, Q, s, k) {
                Some(v) => {
                    raw.insert((s, k), v);
                }
                None => misses.push((s, k)),
            }
        }
    }
    println!(
        "[0] モード表: キャッシュ命中 {} / 72, 対角化 {} 本 (並列度 {})",
        72 - misses.len(),
        misses.len(),
        par
    );
    for chunk in misses.chunks(par) {
        let hs: Vec<_> = chunk
            .iter()
            .map(|&(s, k)| {
                (
                    s,
                    k,
                    std::thread::spawn(move || flux_modes_shear_n(n, k, s)),
                )
            })
            .collect();
        for (s, k, h) in hs {
            let v = h.join().unwrap();
            cache_save_modes(MODE_TAG, n, Q, s, k, &v.0, v.1, v.2);
            raw.insert((s, k), v);
        }
        println!(
            "    ... {} / 72 本 ({} s)",
            raw.len(),
            t0.elapsed().as_secs()
        );
    }
    let mut locs_by_s: Vec<Vec<Vec<Mode>>> = Vec::new();
    for s in 0..6usize {
        let spread = (0..NK12).map(|k| raw[&(s, k)].2).fold(0.0f64, f64::max);
        check(
            &format!("s={} の厳密 3 重縮退", s),
            spread < 1e-8,
            format!("幅 {:.1e}", spread),
        );
        locs_by_s.push(
            (0..NK12)
                .map(|k| localize_stable(n, &raw[&(s, k)].0))
                .collect(),
        );
    }
    drop(raw);

    // ---- 原理スコアの定義 (v15.5 と同一の operationalization) ----
    let mdl_comp = |z: usize| -> f64 {
        if z == 0 {
            (0.5f64).ln()
        } else {
            (0.1f64).ln()
        }
    };
    let mdl_k = |a: usize| -> f64 { mdl_comp(a % 6) + mdl_comp(a / 6) };
    let mdl_s3 = |s: usize| -> f64 {
        if s == 0 {
            (0.5f64).ln()
        } else {
            (0.1f64).ln()
        }
    };
    // 幾何レベル (新層): シアー整数の符号長 — s=0 が「符号なし」
    let mdl_geo = |s1: usize, s2: usize| -> f64 { mdl_s3(s1) + mdl_s3(s2) };

    // 0: 一様 / 1: MDL / 2-3: Robust λ=0.5,1 / 4-5: Depth β=0.25,0.5 / 6-7: Thermo γ=0.5,1
    const NP: usize = 8;
    let lam_rob = [0.5f64, 1.0];
    let beta_dep = [0.25f64, 0.5];
    let gam_th = [0.5f64, 1.0];
    let pnames = [
        "一様 (基線)",
        "P1 MDL",
        "P2 Robust λ=.5",
        "P2 Robust λ=1",
        "P3 Depth β=.25",
        "P3 Depth β=.5",
        "P4 Thermo γ=.5",
        "P4 Thermo γ=1",
    ];

    let sigma = (2.0f64).ln();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let nc = 36usize;

    // ---- [1] 幾何ごとの原理つき証拠 zg[g][p] ----
    println!("\n[1] 幾何ごとの原理つき証拠 (config 層は幾何内で正規化):");
    let geoms: Vec<(usize, usize)> = (0..6usize)
        .flat_map(|a| (a..6usize).map(move |b| (a, b)))
        .collect();
    let mut zg: Vec<[f64; NP]> = Vec::new();
    let mut oracle = f64::NEG_INFINITY;
    let t1 = std::time::Instant::now();
    for &(s1, s2) in &geoms {
        let (locs1, locs2) = (&locs_by_s[s1], &locs_by_s[s2]);
        let mut terms: Vec<[f64; NP]> = Vec::new(); // σ_H ごとの正規化済み lnZq'+lnZe'
        for &sh in &sig_grid {
            let ytab1: Vec<M3> = (0..NK12 * NK12)
                .map(|ab| yukawa_n(n, &locs1[ab % NK12], &locs1[ab / NK12], sh))
                .collect();
            let ytab2: Vec<M3> = (0..NK12 * NK12)
                .map(|ab| yukawa_n(n, &locs2[ab % NK12], &locs2[ab / NK12], sh))
                .collect();
            let pair_y = |a: usize, b: usize, sf: usize, sg: usize| -> M3 {
                let (a1, a2) = (2 * (a % 6), 2 * (a / 6));
                let (b1, b2) = (2 * (b % 6), 2 * (b / 6));
                had_prod_perm(&ytab1[a1 + b1 * NK12], &ytab2[a2 + b2 * NK12], sf, sg)
            };
            // ペアキャッシュ: 質量対数比・左固有ベクトル・ln‖Y‖_F
            let mut pair_r: Vec<[f64; 2]> = Vec::with_capacity(nc * nc * 6);
            let mut pair_v: Vec<M3> = Vec::with_capacity(nc * nc * 6);
            let mut pair_f: Vec<f64> = Vec::with_capacity(nc * nc * 6);
            for m in 0..nc * nc * 6 {
                let y = pair_y(m % nc, (m / nc) % nc, 0, m / (nc * nc));
                let (r, v) = mass_and_vecs(&y);
                let fro: f64 = y
                    .iter()
                    .flatten()
                    .map(|&(a, b)| a * a + b * b)
                    .sum::<f64>()
                    .sqrt();
                pair_r.push(r);
                pair_v.push(v);
                pair_f.push(fro.max(1e-300).ln());
            }
            // 感度表 (Robustness)
            let mut sens: Vec<f64> = vec![0.0; nc * nc * 6];
            for s in 0..6 {
                for b in 0..nc {
                    for a in 0..nc {
                        let r0 = &pair_r[a + b * nc + s * nc * nc];
                        let mut acc = 0.0;
                        for &a2 in &comp_neighbors(a) {
                            let r = &pair_r[a2 + b * nc + s * nc * nc];
                            acc += (r[0] - r0[0]).abs() + (r[1] - r0[1]).abs();
                        }
                        for &b2 in &comp_neighbors(b) {
                            let r = &pair_r[a + b2 * nc + s * nc * nc];
                            acc += (r[0] - r0[0]).abs() + (r[1] - r0[1]).abs();
                        }
                        sens[a + b * nc + s * nc * nc] = acc / 8.0;
                    }
                }
            }
            // e セクター
            let mut ze_sh = [Acc::new(); NP];
            let mut ne_sh = [Acc::new(); NP];
            let mut ora_e = f64::NEG_INFINITY;
            let mut er: Vec<[f64; 2]> = Vec::with_capacity(nc * nc * 36);
            let mut ef: Vec<f64> = Vec::with_capacity(nc * nc * 36);
            for sl in 0..6 {
                for se_ in 0..6 {
                    for ab in 0..nc * nc {
                        let y = pair_y(ab % nc, ab / nc, sl, se_);
                        er.push(mass_ratios(&y));
                        let fro: f64 = y
                            .iter()
                            .flatten()
                            .map(|&(p, q)| p * p + q * q)
                            .sum::<f64>()
                            .sqrt();
                        ef.push(fro.max(1e-300).ln());
                    }
                }
            }
            for sl in 0..6 {
                for se_ in 0..6 {
                    let base = (sl * 6 + se_) * nc * nc;
                    for ab in 0..nc * nc {
                        let (a, b) = (ab % nc, ab / nc);
                        let r = &er[base + ab];
                        let l = ll2(r, tgt[4], tgt[5]);
                        if l > ora_e {
                            ora_e = l;
                        }
                        let mut acc = 0.0;
                        for &a2 in &comp_neighbors(a) {
                            let r2 = &er[base + a2 + b * nc];
                            acc += (r2[0] - r[0]).abs() + (r2[1] - r[1]).abs();
                        }
                        for &b2 in &comp_neighbors(b) {
                            let r2 = &er[base + a + b2 * nc];
                            acc += (r2[0] - r[0]).abs() + (r2[1] - r[1]).abs();
                        }
                        let sens_e = acc / 8.0;
                        let mdl = mdl_k(a) + mdl_k(b) + mdl_s3(sl) + mdl_s3(se_);
                        for p in 0..NP {
                            let s_p = match p {
                                0 => 0.0,
                                1 => mdl,
                                2 | 3 => -lam_rob[p - 2] * sens_e,
                                4 | 5 => -beta_dep[p - 4] * (r[0] + r[1]),
                                _ => gam_th[p - 6] * ef[base + ab],
                            };
                            ze_sh[p].add(s_p + l);
                            ne_sh[p].add(s_p);
                        }
                    }
                }
            }
            // クォーク五重和 (10 量: CKM 3 成分 + J)
            let mut zq_sh = [Acc::new(); NP];
            let mut nq_sh = [Acc::new(); NP];
            let mut ora_q = f64::NEG_INFINITY;
            for kq in 0..nc {
                let mdl_q = mdl_k(kq);
                for su in 0..6 {
                    for ku in 0..nc {
                        let mu = kq + ku * nc + su * nc * nc;
                        let ru = &pair_r[mu];
                        let vu = &pair_v[mu];
                        let llu = ll2(ru, tgt[0], tgt[1]);
                        let mdl_u = mdl_q + mdl_k(ku) + mdl_s3(su);
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
                                let j = jarlskog(&v);
                                let dj = j.abs().max(1e-300).ln() - j_obs.ln();
                                let ll10 = ll + (-dj * dj / (2.0 * sigma * sigma) + norm1);
                                if ll10 > ora_q {
                                    ora_q = ll10;
                                }
                                for p in 0..NP {
                                    let s_p = match p {
                                        0 => 0.0,
                                        1 => mdl_u + mdl_k(kd) + mdl_s3(sd),
                                        2 | 3 => -lam_rob[p - 2] * (sens[mu] + sens[md]),
                                        4 | 5 => -beta_dep[p - 4] * (ru[0] + ru[1] + rd[0] + rd[1]),
                                        _ => gam_th[p - 6] * (pair_f[mu] + pair_f[md]),
                                    };
                                    zq_sh[p].add(s_p + ll10);
                                    nq_sh[p].add(s_p);
                                }
                            }
                        }
                    }
                }
            }
            let mut row = [0.0f64; NP];
            for p in 0..NP {
                row[p] = (zq_sh[p].val() - nq_sh[p].val()) + (ze_sh[p].val() - ne_sh[p].val());
            }
            terms.push(row);
            if ora_q + ora_e > oracle {
                oracle = ora_q + ora_e;
            }
        }
        let mut z = [0.0f64; NP];
        for p in 0..NP {
            let col: Vec<f64> = terms.iter().map(|r| r[p]).collect();
            z[p] = lse(&col) - (sig_grid.len() as f64).ln();
        }
        println!(
            "    ({},{})  一様 {:8.3}  MDL {:8.3}  Rob(1) {:8.3}  Dep(.25) {:9.3}  Th(1) {:8.3}   ({} s)",
            s1,
            s2,
            z[0],
            z[1],
            z[3],
            z[4],
            z[7],
            t1.elapsed().as_secs()
        );
        zg.push(z);
    }

    // ---- ゲート: 一様の per-geometry が v16.9 を再現 ----
    println!("\n[ゲート]");
    let mut max_dev: f64 = 0.0;
    for (gi, &(s1, s2, r)) in REF_V169.iter().enumerate() {
        assert_eq!(geoms[gi], (s1, s2));
        max_dev = max_dev.max((zg[gi][0] - r).abs());
    }
    check(
        "一様 per-geometry の v16.9 回帰 (21 幾何, ±0.02)",
        max_dev < 0.02,
        format!("最大偏差 {:.1e}", max_dev),
    );
    let uni_marg = {
        let col: Vec<f64> = zg.iter().map(|z| z[0]).collect();
        lse(&col) - (geoms.len() as f64).ln()
    };
    check(
        "一様二層合計の v16.9 marginal 回帰 (±0.02)",
        (uni_marg - REF_MARGINAL).abs() < 0.02,
        format!("{:.3} vs {:.3}", uni_marg, REF_MARGINAL),
    );

    // ---- [2] 二層合成と生存判定 ----
    println!("\n[2] 二層合成 lnZ_P (幾何層: MDL のみ明示スコア、他は一様):");
    let mut lnz_tot = [0.0f64; NP];
    for p in 0..NP {
        let (mut zacc, mut wacc) = (Acc::new(), Acc::new());
        for (gi, &(s1, s2)) in geoms.iter().enumerate() {
            let w = if p == 1 { mdl_geo(s1, s2) } else { 0.0 };
            zacc.add(w + zg[gi][p]);
            wacc.add(w);
        }
        lnz_tot[p] = zacc.val() - wacc.val();
    }
    // marginalize (パラメータつき原理は 2 点グリッド, Occam 罰 ln2)
    let z_rob = lse(&[lnz_tot[2], lnz_tot[3]]) - (2.0f64).ln();
    let z_dep = lse(&[lnz_tot[4], lnz_tot[5]]) - (2.0f64).ln();
    let z_th = lse(&[lnz_tot[6], lnz_tot[7]]) - (2.0f64).ln();
    let verdict = |d: f64| -> &'static str {
        if d >= 1.0 {
            "生存"
        } else if d <= -1.0 {
            "棄却"
        } else {
            "未決"
        }
    };
    let base = lnz_tot[0];
    println!("    原理                lnZ        Δ vs 一様   判定");
    println!("    一様 (基線)      {:8.3}      —          —", base);
    println!(
        "    P1 MDL           {:8.3}   {:+8.3}     {}",
        lnz_tot[1],
        lnz_tot[1] - base,
        verdict(lnz_tot[1] - base)
    );
    println!(
        "    P2 Robustness    {:8.3}   {:+8.3}     {}   [λ=.5: {:+.2}, λ=1: {:+.2}]",
        z_rob,
        z_rob - base,
        verdict(z_rob - base),
        lnz_tot[2] - base,
        lnz_tot[3] - base
    );
    println!(
        "    P3 Depth         {:8.3}   {:+8.3}     {}   [β=.25: {:+.2}, β=.5: {:+.2}]",
        z_dep,
        z_dep - base,
        verdict(z_dep - base),
        lnz_tot[4] - base,
        lnz_tot[5] - base
    );
    println!(
        "    P4 Thermo        {:8.3}   {:+8.3}     {}   [γ=.5: {:+.2}, γ=1: {:+.2}]",
        z_th,
        z_th - base,
        verdict(z_th - base),
        lnz_tot[6] - base,
        lnz_tot[7] - base
    );
    println!(
        "    (oracle 上界)    {:8.3}   {:+8.3}     完全情報の理論上界",
        oracle,
        oracle - base
    );
    println!("    v15.5 (rect N=18, 9 量) との比較: Thermo は +0.32 だった — 符号の再現性が焦点");

    // ---- [3] 谷底の裁定 ----
    println!("\n[3] 谷底 4 幾何の原理別ランドスケープ (Δ = 各幾何 − その原理の谷底最良):");
    let valley: [(usize, usize); 4] = [(2, 2), (2, 3), (1, 3), (3, 3)];
    let vidx: Vec<usize> = valley
        .iter()
        .map(|g| geoms.iter().position(|x| x == g).unwrap())
        .collect();
    println!("    原理             (2,2)    (2,3)    (1,3)    (3,3)    最大対差");
    let mut broken: Vec<(String, f64, (usize, usize))> = Vec::new();
    for p in 0..NP {
        let vals: Vec<f64> = vidx
            .iter()
            .map(|&gi| {
                zg[gi][p]
                    + if p == 1 {
                        mdl_geo(geoms[gi].0, geoms[gi].1)
                    } else {
                        0.0
                    }
            })
            .collect();
        let best = vals.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let worst = vals.iter().cloned().fold(f64::INFINITY, f64::min);
        let top = valley[vals.iter().position(|&v| v == best).unwrap()];
        println!(
            "    {:14} {:+8.2} {:+8.2} {:+8.2} {:+8.2}   {:5.2}  (首位 ({},{}))",
            pnames[p],
            vals[0] - best,
            vals[1] - best,
            vals[2] - best,
            vals[3] - best,
            best - worst,
            top.0,
            top.1
        );
        if p > 0 && best - worst >= 1.0 {
            broken.push((pnames[p].to_string(), best - worst, top));
        }
    }
    if broken.is_empty() {
        println!("    => どの原理も谷底 4 幾何の縮退を ≥1 nat では破らない — 縮退は原理にも頑健。");
    } else {
        for (name, d, top) in &broken {
            println!(
                "    => {} は縮退を破る (最大対差 {:.2} nats, 首位 ({},{})) — 事前登録の裁定基準を満たす。",
                name, d, top.0, top.1
            );
        }
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.11".into())),
        (
            "lnz_total".into(),
            Json::Arr(lnz_tot.iter().map(|&x| Json::Num(x)).collect()),
        ),
        ("lnz_rob_marg".into(), Json::Num(z_rob)),
        ("lnz_dep_marg".into(), Json::Num(z_dep)),
        ("lnz_th_marg".into(), Json::Num(z_th)),
        ("oracle".into(), Json::Num(oracle)),
        ("uniform_marginal".into(), Json::Num(uni_marg)),
        (
            "zg".into(),
            Json::Arr(
                zg.iter()
                    .enumerate()
                    .map(|(gi, z)| {
                        Json::Obj(vec![
                            ("s1".into(), Json::Int(geoms[gi].0 as i64)),
                            ("s2".into(), Json::Int(geoms[gi].1 as i64)),
                            (
                                "lnz_p".into(),
                                Json::Arr(z.iter().map(|&x| Json::Num(x)).collect()),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v1611_selection2.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 判別は [2][3] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
