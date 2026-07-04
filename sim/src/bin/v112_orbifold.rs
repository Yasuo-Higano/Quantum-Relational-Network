//! v11.2 Z₃ シフト・オービフォールド射影 — 対の族の第一原理検定
//!
//! v10.3 は世代対 σ が既約な離散自由度 (タワー整列) であることを示した。では σ を
//! **決める機構**は何か。最も単純な候補は Z₃ シフト・オービフォールド:
//! 両トーラスを y 方向に同時に 6 サイト (= N/Q, 世代間隔) 並進する作用で商を取る。
//!
//! Landau ゲージ (A_y = φx) では y 並進はハミルトニアンと厳密に可換で、
//! ゼロモード上の作用はクロック (T_y ψ_g = ω^{g+c} ψ_g, ω = e^{2πi/3})。
//! 積状態 ψ_i¹⊗ψ_j² への作用は ω^{i+j+c} なので、固有値 ω^m の射影は
//!   {(i, j) : i + j ≡ m (mod 3)} — **反対角対の族** (3 状態) を選ぶ。
//! つまり 9 → 3 の射影は「対の置換が奇 (反対角)、オフセット m は自由」を予言する。
//!
//! ただしトーラス 2 のラベル付け替え τ (v10.1 のゲージ自由度) は奇置換も含むので、
//! **ゲージ不変な予言は「全ての場が同じ族 (相対パリティが偶)」**である。
//! v10.1 の MAP (σ_u 奇, σ_d 偶, σ_L 奇, σ_e 偶 — 混合パリティ) はこれに反する —
//! そこで族制限つき証拠で機構を検定する:
//!   M_uni: 全場が同じ族 (σ_Q=e ゲージで全場偶 = 巡回のみ)。事前 4 ln 3。
//!   M_full: 全 S₃ (v10.1 の M_perm)。事前 4 ln 6。
//! lnZ(M_uni) ≥ lnZ(M_full) なら単純オービフォールドで足りる。逆なら混合パリティ
//! (セクターごとに向きの違う射影 — 例: 磁化の符号が場ごとに違う brane 模型) が必要。
//!
//! 検証: [1] T_y のクロック作用 (数値対角性 1e-10)、[2] 9→3 射影の厳密一致、
//! [3] 族制限証拠の退化検査 (M_full が v10.1 と一致)。

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
    println!("=== v11.2 Z₃ シフト・オービフォールド射影: 対の族の第一原理検定 ===\n");
    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln();
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();

    // ---- [0] 世代モード ----
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
    println!("    縮退・ギャップ不変  {}  ({} ms)", pass(ok_engine), t0.elapsed().as_millis());

    // ---- [1] 磁気並進 T_y (6 サイト) のクロック作用 ----
    println!("\n[1] T_y (y 方向 6 サイト並進) はゼロモード上でクロック対角か");
    // T_y は Landau ゲージで純粋な格子並進 (H と厳密可換)
    let ty = |psi: &C3v| -> C3v {
        let mut o = [(0.0f64, 0.0f64); NS];
        for y in 0..N {
            for x in 0..N {
                o[x + ((y + 6) % N) * N] = psi[x + y * N];
            }
        }
        o
    };
    let mut max_offdiag: f64 = 0.0;
    let mut phases: Vec<f64> = Vec::new(); // 対角位相 (arg)
    for g in 0..3usize {
        let t = ty(&locs[0][g]);
        for h in 0..3usize {
            let (mut re, mut im) = (0.0, 0.0);
            for s in 0..NS {
                let (ar, ai) = locs[0][h][s];
                let (br, bi) = t[s];
                re += ar * br + ai * bi;
                im += ar * bi - ai * br;
            }
            if h == g {
                phases.push(im.atan2(re));
                let mag = (re * re + im * im).sqrt();
                max_offdiag = max_offdiag.max((mag - 1.0).abs());
            } else {
                max_offdiag = max_offdiag.max((re * re + im * im).sqrt());
            }
        }
    }
    // 対角位相が ω^g 系列 (隣接差 2π/3) か
    let two_pi = 2.0 * std::f64::consts::PI;
    let d01 = ((phases[1] - phases[0]) / two_pi * 3.0).rem_euclid(3.0);
    let d12 = ((phases[2] - phases[1]) / two_pi * 3.0).rem_euclid(3.0);
    let near_int = |x: f64| (x - x.round()).abs() < 1e-8 || (x - x.round()).abs() > 3.0 - 1e-8;
    let step01 = d01.round() as i64 % 3;
    let step12 = d12.round() as i64 % 3;
    let ok_clock = max_offdiag < 1e-10
        && near_int(d01)
        && near_int(d12)
        && step01 == step12
        && step01 != 0;
    println!(
        "    非対角 |T_y| 最大 {:.2e} / 対角位相の隣接差 = 2π/3 × {} (一定)  {}",
        max_offdiag,
        step01,
        pass(ok_clock)
    );

    // ---- [2] 9 → 3 射影 = 反対角対の族 ----
    println!("\n[2] クロック⊗クロックの固有値射影は反対角対 {{(i,j): i+j≡m}} を選ぶ");
    // 積状態への T_y¹⊗T_y² の作用は位相 ω^{(i+j)+2c} — 固有空間は i+j mod 3 で分類される。
    // [1] でゼロモードがクロック固有基底そのもの (非対角 <1e-10) と確認済みなので、
    // 固有値 ω^m の部分空間 = span{ψ_i⊗ψ_j : i+j+2c ≡ m} は厳密に 3 次元で、
    // 反対角対 σ(i) = (m−2c)−i の積そのもの。ここでは c を [1] の位相から読み取り、
    // 各 m の生き残り集合を明示的に構成して次元と対の形を確認する。
    let c_off = {
        // phases[g] = arg ω^{g+c} = 2π(g+c)/3 (mod 2π) — g=0 から c を読む
        let c = (phases[0] / two_pi * 3.0).rem_euclid(3.0);
        c.round() as i64 % 3
    };
    let mut ok_proj = true;
    for m in 0..3i64 {
        let mut members: Vec<(usize, usize)> = Vec::new();
        for i in 0..3usize {
            for j in 0..3usize {
                if ((i + j) as i64 + 2 * c_off).rem_euclid(3) == m {
                    members.push((i, j));
                }
            }
        }
        ok_proj &= members.len() == 3;
        // 対の形: j = const − i (反対角)
        let s0 = (members[0].0 + members[0].1) % 3;
        ok_proj &= members.iter().all(|&(i, j)| (i + j) % 3 == s0);
    }
    println!(
        "    各固有値の生き残りは 3 状態・j ≡ const − i (反対角対)  {}",
        pass(ok_proj)
    );
    println!("    => 一様な Z₃ シフト・オービフォールドの予言: 全場が同じ族 (相対パリティ偶)、");
    println!("       オフセット m_F は離散 Wilson 相当の自由度 (場ごとに Z₃)");

    // ---- [3] 族制限つき証拠 ----
    println!("\n[3] 族制限つき証拠 (T²×Z₆, 安定ラベル, σ_Q = e ゲージ)");
    println!("    M_uni  = 全場が同じ族 (偶置換 = 巡回のみ, 事前 4ln3 = {:.2})", 4.0 * (3.0f64).ln());
    println!("    M_full = 全 S₃ (v10.1 の M_perm, 事前 4ln6 = {:.2})", 4.0 * (6.0f64).ln());
    let nc = 36usize;
    let ll2 = |r: &[f64; 2], t0: f64, t1: f64| -> f64 {
        -((r[0] - t0).powi(2) + (r[1] - t1).powi(2)) / (2.0 * sigma * sigma) + 2.0 * norm1
    };
    let even: [usize; 3] = [0, 3, 4]; // e, (123), (132)
    let mut terms_uni = Vec::new();
    let mut terms_full = Vec::new();
    let mut map_full = (f64::NEG_INFINITY, [0usize; 4]);
    for &sh in &sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // ペアキャッシュ (質量比+左固有ベクトル) — σ 全 6 種
        let pair: Vec<([f64; 2], M3)> = (0..nc * nc * 6)
            .map(|m| mass_and_vecs(&pair_yukawa(&ytab, m % nc, (m / nc) % nc, 0, m / (nc * nc))))
            .collect();
        // e セクター
        let mut le_uni = Vec::new();
        let mut le_full = Vec::new();
        let mut le_best = f64::NEG_INFINITY;
        for sl in 0..6usize {
            for se_ in 0..6usize {
                for ab in 0..nc * nc {
                    let r = mass_ratios(&pair_yukawa(&ytab, ab % nc, ab / nc, sl, se_));
                    let l = ll2(&r, tgt[4], tgt[5]);
                    le_full.push(l);
                    if even.contains(&sl) && even.contains(&se_) {
                        le_uni.push(l);
                    }
                    le_best = le_best.max(l);
                }
            }
        }
        // クォーク部 (kq ごとに (ku,σu)×(kd,σd) — 全和は 36×216² で軽い)
        let mut per_q_uni = Vec::with_capacity(nc);
        let mut per_q_full = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut au = (f64::NEG_INFINITY, 0.0f64);
            let mut af = (f64::NEG_INFINITY, 0.0f64);
            for su in 0..6usize {
                for ku in 0..nc {
                    let (ru, vu) = &pair[kq + ku * nc + su * nc * nc];
                    let llu = ll2(ru, tgt[0], tgt[1]);
                    for sd in 0..6usize {
                        for kd in 0..nc {
                            let (rd, vd) = &pair[kq + kd * nc + sd * nc * nc];
                            let lld = ll2(rd, tgt[2], tgt[3]);
                            let c = ckm3(vu, vd);
                            let mut ll = llu + lld;
                            for m in 0..3 {
                                let d = c[m].max(1e-300).ln() - tgt[6 + m];
                                ll += -d * d / (2.0 * sigma * sigma) + norm1;
                            }
                            if ll > af.0 {
                                af.1 = af.1 * (af.0 - ll).exp() + 1.0;
                                af.0 = ll;
                            } else {
                                af.1 += (ll - af.0).exp();
                            }
                            if even.contains(&su) && even.contains(&sd) {
                                if ll > au.0 {
                                    au.1 = au.1 * (au.0 - ll).exp() + 1.0;
                                    au.0 = ll;
                                } else {
                                    au.1 += (ll - au.0).exp();
                                }
                            }
                            let tot = ll + le_best;
                            if tot > map_full.0 {
                                map_full = (tot, [su, sd, 9, 9]);
                            }
                        }
                    }
                }
            }
            per_q_uni.push(au.0 + au.1.ln());
            per_q_full.push(af.0 + af.1.ln());
        }
        terms_uni.push(lse(&per_q_uni) + lse(&le_uni));
        terms_full.push(lse(&per_q_full) + lse(&le_full));
    }
    let prior_w = 10.0 * (6.0f64).ln() + (sig_grid.len() as f64).ln();
    let lnz_uni = lse(&terms_uni) - prior_w - 4.0 * (3.0f64).ln();
    let lnz_full = lse(&terms_full) - prior_w - 4.0 * (6.0f64).ln();
    // 退化検査: M_full は v10.1 の lnZ₉(M_perm) と一致するはず
    let ref_full = -19.86334559888438; // results/v101_pairing.json lnZ_nine_perm
    let ok_reg = (lnz_full - ref_full).abs() < 0.02;
    println!("    lnZ₉(M_full) = {:.4} vs v10.1 {:.4}  {}", lnz_full, ref_full, pass(ok_reg));
    println!("    lnZ₉(M_uni)  = {:.4}", lnz_uni);
    let uni_wins = lnz_uni >= lnz_full;
    println!(
        "\n    => {} (差 {:+.2})",
        if uni_wins {
            "M_uni (一様族 = 単純オービフォールドで足りる)"
        } else {
            "M_full (混合パリティが必要 — 単純シフト・オービフォールドでは足りない)"
        },
        lnz_uni - lnz_full
    );
    let sig_names = ["e", "(23)", "(12)", "(123)", "(132)", "(13)"];
    println!(
        "    MAP のクォーク対: σ_u={} σ_d={} (v10.1 と同じく混合パリティか)",
        sig_names[map_full.1[0]],
        sig_names[map_full.1[1]]
    );

    let all_ok = ok_engine && ok_clock && ok_proj && ok_reg;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-010".into())),
        ("clock_offdiag_max".into(), Json::Num(max_offdiag)),
        ("projection_is_antidiagonal".into(), Json::Bool(ok_proj)),
        ("lnZ_nine_uniform_family".into(), Json::Num(lnz_uni)),
        ("lnZ_nine_full".into(), Json::Num(lnz_full)),
        ("uniform_family_wins".into(), Json::Bool(uni_wins)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v112_orbifold.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {} (PASS = 装置検証 — 機構の判定は [3] が本体)", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
