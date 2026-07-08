//! v15.5 選択原理の有限候補戦 — σ・Wilson 線を「選ぶ原理」を証拠にかける
//!
//! 残高 2/9 (v15.0): σ・磁束・Wilson 線・トーラス数を一括で決める上位原理は
//! 真空選択問題と同格の最深残高である。次の方針 (PROMPT/2 柱 4) は「もっと良い
//! fit を探す」ことではなく、**選択原理の候補を有限個にして全部落とすか残すかを
//! 決める**ことを求める。本バイナリは v10.1 の T²×T² 証拠エンジン (五重和・
//! σ marginalize・安定ラベル) の上で、データを見ない事前分布 π_P(config) を持つ
//! 原理候補 P を一様事前 (= 無知) と同じ土俵で比較する:
//!
//!     lnZ_P = ln Σ_c π_P(c) L(data|c)   (π_P は正規化 — 集中が当たれば得、外れれば損)
//!
//! 候補 (全て config のみから計算 — 実測値は使わない):
//!   P1 MDL:        接頭符号事前 — Wilson 成分 k=0 とσ=e に確率 1/2 (単純さの重み)。
//!                  自由パラメータ 0 (符号の設計は計数論から固定)。
//!   P2 Robustness: Wilson 隣接 (±1 目盛) シフトに対する質量対数比の感度が小さい
//!                  config を好む — 「自己安定な真空」。パラメータ λ (marginalize)。
//!   P3 Depth:      階層の深い config を好む (−β(ln r1 + ln r2) — 実測の深さは見ない)。
//!                  「真空は階層を最大化する」型。パラメータ β (marginalize)。
//!   P4 Thermo:     湯川の総結合 ln‖Y‖_F を凝縮エネルギーの代理として重み付け。
//!                  「真空選択の熱力学」の最小の operationalization。パラメータ γ。
//!   P5 計算可能性: この有限窓では全 config が無矛盾に符号化可能 — 判別力を持たない
//!                  (空な原理として記録し、C5 スローガンから削る)。
//!
//! 事前登録の判定規準 (結果を見る前に固定):
//!   生存: ΔlnZ ≥ +1.0 nat / 棄却: ΔlnZ ≤ −1.0 nat / それ以外: 未決。
//!   パラメータつき原理はグリッド marginalize (Occam 罰込み) で比較する。
//!
//! 検証: (i) 一様事前が v10.1 の公表値 (M_perm 全 9 量 −19.86, M_diag −23.61) を
//! ±0.02 で回帰再現 (エンジンの較正)。(ii) MDL の正規化子は構成から厳密に 0 —
//! 機械和が 1e-12 で 0 になること (積和機構の自己検査)。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const NC: usize = 36;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];
const REF_NINE_DIAG: f64 = -23.61; // v9.2/v10.1 の一次ソース
const REF_NINE_PERM: f64 = -19.86; // v10.1 の公表値

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

// ---- v10.1 と同一の物理エンジン (回帰検査つき再利用) ----

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

/// ストリーミング log-sum-exp
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

/// Z6×Z6 複合添字の隣接 (±1 目盛 × 2 成分 = 4 近傍)
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
    println!("=== v15.5 選択原理の有限候補戦: σ・Wilson 線を選ぶ原理を証拠にかける ===\n");
    println!("事前登録の判定規準 (結果より先に固定):");
    println!("  生存: ΔlnZ ≥ +1.0 nat / 棄却: ΔlnZ ≤ −1.0 nat / 中間: 未決");
    println!(
        "  パラメータつき原理はグリッド marginalize (Occam 罰込み)。基線 = 一様事前 (無知)。\n"
    );

    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };

    // ---- [0] 世代モード (安定ラベル) ----
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
        "[0] 世代モード構築 (12 対角化, 安定ラベル): 縮退・ギャップ  {}  ({} ms)",
        pass(ok_engine),
        t0.elapsed().as_millis()
    );

    // ---- 原理スコアの定義 ----
    // MDL: Z6 成分ごとに P(0)=1/2, P(k≠0)=1/10; σ ごとに P(e)=1/2, P(σ≠e)=1/10。
    let mdl_comp = |z: usize| -> f64 {
        if z == 0 {
            (0.5f64).ln()
        } else {
            (0.1f64).ln()
        }
    };
    let mdl_k = |a: usize| -> f64 { mdl_comp(a % 6) + mdl_comp(a / 6) };
    let mdl_s = |s: usize| -> f64 {
        if s == 0 {
            (0.5f64).ln()
        } else {
            (0.1f64).ln()
        }
    };

    // 原理の一覧 (名前, 変種ラベル)。lnZ は変種ごとに計算し、後で marginalize。
    // 0: 一様 (基線) / 1: MDL / 2-3: Robust λ=0.5,1 / 4-5: Depth β=0.25,0.5 / 6-7: Thermo γ=0.5,1
    const NP: usize = 8;
    let lam_rob = [0.5f64, 1.0];
    let beta_dep = [0.25f64, 0.5];
    let gam_th = [0.5f64, 1.0];

    // ---- σ_H ごとの証拠計算 ----
    // 事前は σ_H (共有ニュイサンス) の下で config に張る:
    //   Z_P = (1/|σ_H|) Σ_σH [Σ_q π_P L_q / Σ_q π_P] × [Σ_e π_P L_e / Σ_e π_P]
    // (クォーク部と e 部は添字を共有しないので σ_H ごとに厳密に因子化する)
    let mut terms: Vec<[f64; NP]> = Vec::new(); // σ_H ごとの lnZq'+lnZe' (正規化済み)
    let mut terms_diag: Vec<f64> = Vec::new();
    let mut nq_mdl_max: f64 = 0.0; // MDL 正規化子の解析値 0 からの最大偏差
    let mut oracle = f64::NEG_INFINITY; // max_c ln L (原理の理論上界)
    let t1 = std::time::Instant::now();
    for &sh in &sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // ペアキャッシュ: (a + b*NC + σ*NC*NC) → (質量対数比, 左固有ベクトル, ln‖Y‖_F)
        let mut pair_r: Vec<[f64; 2]> = Vec::with_capacity(NC * NC * 6);
        let mut pair_v: Vec<M3> = Vec::with_capacity(NC * NC * 6);
        let mut pair_f: Vec<f64> = Vec::with_capacity(NC * NC * 6);
        for m in 0..NC * NC * 6 {
            let y = pair_yukawa(&ytab, m % NC, (m / NC) % NC, 0, m / (NC * NC));
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
        // 感度表 (Robustness): (a,b,σ) → Wilson 近傍 8 点の |Δ(ln 比)|₁ の平均
        let mut sens: Vec<f64> = vec![0.0; NC * NC * 6];
        for s in 0..6 {
            for b in 0..NC {
                for a in 0..NC {
                    let r0 = &pair_r[a + b * NC + s * NC * NC];
                    let mut acc = 0.0;
                    for &a2 in &comp_neighbors(a) {
                        let r = &pair_r[a2 + b * NC + s * NC * NC];
                        acc += (r[0] - r0[0]).abs() + (r[1] - r0[1]).abs();
                    }
                    for &b2 in &comp_neighbors(b) {
                        let r = &pair_r[a + b2 * NC + s * NC * NC];
                        acc += (r[0] - r0[0]).abs() + (r[1] - r0[1]).abs();
                    }
                    sens[a + b * NC + s * NC * NC] = acc / 8.0;
                }
            }
        }

        // ---- e セクター (kL, ke, σL, σe — クォークと独立に因子化) ----
        let mut ze_sh = [Acc::new(); NP];
        let mut ne_sh = [Acc::new(); NP];
        let mut ze_diag_sh = Acc::new();
        let mut ora_e = f64::NEG_INFINITY;
        // 比のキャッシュ (σL, σe 込み): 46656 点
        let mut er: Vec<[f64; 2]> = Vec::with_capacity(NC * NC * 36);
        let mut ef: Vec<f64> = Vec::with_capacity(NC * NC * 36);
        for sl in 0..6 {
            for se_ in 0..6 {
                for ab in 0..NC * NC {
                    let y = pair_yukawa(&ytab, ab % NC, ab / NC, sl, se_);
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
                let base = (sl * 6 + se_) * NC * NC;
                for ab in 0..NC * NC {
                    let (a, b) = (ab % NC, ab / NC);
                    let r = &er[base + ab];
                    let l = ll2(r, tgt[4], tgt[5]);
                    if l > ora_e {
                        ora_e = l;
                    }
                    let mut acc = 0.0;
                    for &a2 in &comp_neighbors(a) {
                        let r2 = &er[base + a2 + b * NC];
                        acc += (r2[0] - r[0]).abs() + (r2[1] - r[1]).abs();
                    }
                    for &b2 in &comp_neighbors(b) {
                        let r2 = &er[base + a + b2 * NC];
                        acc += (r2[0] - r[0]).abs() + (r2[1] - r[1]).abs();
                    }
                    let sens_e = acc / 8.0;
                    let mdl = mdl_k(a) + mdl_k(b) + mdl_s(sl) + mdl_s(se_);
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
                    if sl == 0 && se_ == 0 {
                        ze_diag_sh.add(l);
                    }
                }
            }
        }

        // ---- クォーク部の五重和 (kQ, σu, ku, σd, kd) — CKM が因子化を壊す ----
        let mut zq_sh = [Acc::new(); NP];
        let mut nq_sh = [Acc::new(); NP];
        let mut zq_diag_sh = Acc::new();
        let mut ora_q = f64::NEG_INFINITY;
        for kq in 0..NC {
            let mdl_q = mdl_k(kq);
            for su in 0..6 {
                for ku in 0..NC {
                    let mu = kq + ku * NC + su * NC * NC;
                    let ru = &pair_r[mu];
                    let vu = &pair_v[mu];
                    let llu = ll2(ru, tgt[0], tgt[1]);
                    let mdl_u = mdl_q + mdl_k(ku) + mdl_s(su);
                    for sd in 0..6 {
                        for kd in 0..NC {
                            let md = kq + kd * NC + sd * NC * NC;
                            let rd = &pair_r[md];
                            let lld = ll2(rd, tgt[2], tgt[3]);
                            let c = ckm3(vu, &pair_v[md]);
                            let mut ll = llu + lld;
                            for m in 0..3 {
                                let d = c[m].max(1e-300).ln() - tgt[6 + m];
                                ll += -d * d / (2.0 * sigma * sigma) + norm1;
                            }
                            if ll > ora_q {
                                ora_q = ll;
                            }
                            for p in 0..NP {
                                let s_p = match p {
                                    0 => 0.0,
                                    1 => mdl_u + mdl_k(kd) + mdl_s(sd),
                                    2 | 3 => -lam_rob[p - 2] * (sens[mu] + sens[md]),
                                    4 | 5 => -beta_dep[p - 4] * (ru[0] + ru[1] + rd[0] + rd[1]),
                                    _ => gam_th[p - 6] * (pair_f[mu] + pair_f[md]),
                                };
                                zq_sh[p].add(s_p + ll);
                                nq_sh[p].add(s_p);
                            }
                            if su == 0 && sd == 0 {
                                zq_diag_sh.add(ll);
                            }
                        }
                    }
                }
            }
        }
        // σ_H ごとの項 (正規化済み積)
        let mut row = [0.0f64; NP];
        for p in 0..NP {
            row[p] = (zq_sh[p].val() - nq_sh[p].val()) + (ze_sh[p].val() - ne_sh[p].val());
        }
        terms.push(row);
        terms_diag.push(
            zq_diag_sh.val() - ((NC * NC * NC) as f64).ln() + ze_diag_sh.val()
                - ((NC * NC) as f64).ln(),
        );
        nq_mdl_max = nq_mdl_max
            .max(nq_sh[1].val().abs())
            .max(ne_sh[1].val().abs());
        oracle = oracle.max(ora_q + ora_e);
        let _ = sh;
    }
    println!(
        "    証拠エンジン (4 σ_H × 五重和 + e 因子): {} ms",
        t1.elapsed().as_millis()
    );

    // ---- [1] エンジンの回帰検査 (一様事前 = v10.1 の公表値) ----
    println!("\n[1] エンジンの回帰検査 (v10.1 の一次ソースとの照合)");
    let ln_sh = (sig_grid.len() as f64).ln();
    let lnz_p: Vec<f64> = (0..NP)
        .map(|p| lse(&terms.iter().map(|r| r[p]).collect::<Vec<_>>()) - ln_sh)
        .collect();
    let lnz_perm = lnz_p[0];
    let ok_perm = (lnz_perm - REF_NINE_PERM).abs() < 0.02;
    println!(
        "    lnZ(一様/M_perm) = {:.3}  (v10.1 公表 {:.2})  {}",
        lnz_perm,
        REF_NINE_PERM,
        pass(ok_perm)
    );
    let lnz_diag = lse(&terms_diag) - ln_sh;
    let ok_diag = (lnz_diag - REF_NINE_DIAG).abs() < 0.02;
    println!(
        "    lnZ(M_diag)      = {:.3}  (v9.2/v10.1 公表 {:.2})  {}",
        lnz_diag,
        REF_NINE_DIAG,
        pass(ok_diag)
    );
    // MDL の正規化子は構成上厳密に 0 (機械和の自己検査)
    let ok_mdl_norm = nq_mdl_max < 1e-9;
    println!(
        "    MDL 正規化子 (解析値 0): 最大偏差 {:.2e}  {}",
        nq_mdl_max,
        pass(ok_mdl_norm)
    );
    if !(ok_perm && ok_diag && ok_mdl_norm) {
        println!("\n総合判定: [FAIL] エンジンの較正が崩れている");
        std::process::exit(1);
    }

    // ---- [2] 候補戦の結果 ----
    println!(
        "\n[2] 候補戦の結果 (基線 = 一様事前 {:.2}, oracle 上界 = max ln L = {:.2})",
        lnz_perm, oracle
    );
    let verdict = |d: f64| -> &'static str {
        if d >= 1.0 {
            "生存"
        } else if d <= -1.0 {
            "棄却"
        } else {
            "未決"
        }
    };
    // marginalize (パラメータつき原理はグリッド 2 点の一様混合 = lse − ln2)
    let z_rob = lse(&[lnz_p[2], lnz_p[3]]) - (2.0f64).ln();
    let z_dep = lse(&[lnz_p[4], lnz_p[5]]) - (2.0f64).ln();
    let z_th = lse(&[lnz_p[6], lnz_p[7]]) - (2.0f64).ln();
    println!("    原理                          自由パラ  lnZ      ΔlnZ    判定");
    println!(
        "    P1 MDL (単純さ)                 0      {:7.2}  {:+6.2}   {}",
        lnz_p[1],
        lnz_p[1] - lnz_p[0],
        verdict(lnz_p[1] - lnz_p[0])
    );
    println!(
        "    P2 Robustness (自己安定)        1(λ)   {:7.2}  {:+6.2}   {}   [λ=0.5: {:+.2}, λ=1: {:+.2}]",
        z_rob,
        z_rob - lnz_p[0],
        verdict(z_rob - lnz_p[0]),
        lnz_p[2] - lnz_p[0],
        lnz_p[3] - lnz_p[0]
    );
    println!(
        "    P3 Depth (階層最大化)           1(β)   {:7.2}  {:+6.2}   {}   [β=0.25: {:+.2}, β=0.5: {:+.2}]",
        z_dep,
        z_dep - lnz_p[0],
        verdict(z_dep - lnz_p[0]),
        lnz_p[4] - lnz_p[0],
        lnz_p[5] - lnz_p[0]
    );
    println!(
        "    P4 Thermo (結合の凝縮)          1(γ)   {:7.2}  {:+6.2}   {}   [γ=0.5: {:+.2}, γ=1: {:+.2}]",
        z_th,
        z_th - lnz_p[0],
        verdict(z_th - lnz_p[0]),
        lnz_p[6] - lnz_p[0],
        lnz_p[7] - lnz_p[0]
    );
    println!("    P5 計算可能性                   —      —        —      空 (この窓では全 config が符号化可能 — 判別力なし, C5 から削除)");

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v15.5".into())),
        ("lnz_uniform".into(), Json::Num(lnz_p[0])),
        ("lnz_diag_regression".into(), Json::Num(lnz_diag)),
        ("oracle_maxll".into(), Json::Num(oracle)),
        ("lnz_mdl".into(), Json::Num(lnz_p[1])),
        ("lnz_robust_marg".into(), Json::Num(z_rob)),
        ("lnz_depth_marg".into(), Json::Num(z_dep)),
        ("lnz_thermo_marg".into(), Json::Num(z_th)),
        (
            "lnz_variants".into(),
            Json::Arr(lnz_p.iter().map(|&x| Json::Num(x)).collect()),
        ),
    ]);
    let p = write_artifact("results/v155_selection.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: [PASS] エンジンは較正済み — 判定は上の表が一次ソース (陰性でも棄却は前進)"
    );
}
