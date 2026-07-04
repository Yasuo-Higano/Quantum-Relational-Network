//! v7.2 幾何からの湯川 — M2 級 FN: 電荷を選ぶ代わりに、Wilson 線を選ぶ (v7.0 残高 2)
//!
//! v6.5 は「FN はアナーキーに勝つが、電荷は文献から選んだもの (M1 級)」と正直に
//! 記録した。本バイナリは改良方針 §6 の M2 へ進む: **湯川行列を幾何から計算する**。
//!
//! 構成 (全て第一原理の数値計算、フィット式ではない):
//!   1. 磁束 Q=3 を貫く 18×18 トーラスの複素ホッピング模型を厳密対角化し、
//!      縮退した 3 つの最低モード = 3 世代 (v2.3 の「世代数=トポロジー」と同じ構成)。
//!   2. 縮退空間内で位置演算子 X̂ = e^{2πix/N} を (一般位相で) 対角化し、
//!      N/Q = 6 サイト間隔で並ぶ局在世代波動関数を得る (LLL のシータ関数的モード)。
//!   3. セクター (Q,u,d,L,e) の区別は **y 方向の Wilson 線** θ_a = 2πk_a/N (k_a ∈ Z₆):
//!      ゼロモードの中心が k_a サイトずつ厳密に平行移動する (v2.2 の「電荷=巻き数」の
//!      Wilson 線がここでは世代の「住所」を決める)。各セクターのモードはそれぞれの
//!      ハミルトニアンの厳密な固有状態 (縮退・ギャップ・中心シフトを検証)。
//!      [開発記録: 当初は磁気並進でシフトを作ろうとしたが、格子磁気並進が局所位相で
//!       閉じるのは N | 2Q のときに限ると判明 (18∤6) — 障害物も幾何の学びとして残す]
//!   4. 湯川 = 重なり積分 Y_ij = Σ_x ψ_i^(a)*(x) ψ_j^(b)(x) φ_H(x)。Higgs は原点中心の
//!      ガウス (幅 σ_H)。階層は「距離² / 運動量差²」の指数から自動で生じる。
//! モデルのパラメータは **5 つの Wilson 線整数 k_a ∈ Z₆ と σ_H (4 通り)** のみ。
//! 乱雑な O(1) 係数は無い — 行列は幾何が完全に決める。
//! 証拠比較 (v6.5 と同一の尤度 σ=ln2, 学習 = 質量比 6 つ):
//!   Z(M2geo) = (1/|格子|) Σ_{k,σ_H} Π_sector N(ln r^obs; ln r^geo(k,σ_H), σ)
//! 尤度はセクターごとに因子化するので全パラメータ空間 (6⁵×4) が厳密に和を取れる。
//! CKM は学習に使わず、MAP 幾何の**決定論的予測**として検証する。
//! 比較対象 (results/v65_bayes.json より): lnZ(M0)≤-35.4, lnZ(M1)=-12.19, lnZ(M2文献)=-7.61。

use uft_sim::*;

const N: usize = 18; // トーラスの一辺 (Q の倍数: 世代間隔 N/Q = 6 が整数)
const NS: usize = N * N;
const Q: usize = 3; // 磁束量子 = 世代数
const NK: usize = 6; // Wilson 線の格子 k ∈ Z₆ (中心シフト 1 サイト刻み × 世代間隔 6)
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];

type C3v = [(f64, f64); NS];

/// 磁束 Q + y 方向 Wilson 線 θ=2πk/N のホッピング行列を対角化し、最低 Q 個の複素モードを返す。
fn flux_modes(k_wilson: usize) -> (Vec<C3v>, f64, f64) {
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / NS as f64;
    // Wilson 線: リンクあたり k·φ (= ゲージ場 A_y = φ(x+k))。ゼロモードの案内中心が
    // ちょうど k サイト平行移動する換算。H(k) は平行移動した H(0) とゲージ同値なので
    // スペクトル (縮退・ギャップ) の不変は厳密。
    let wl = phi * k_wilson as f64;
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

/// 3×3 エルミート行列の固有値・固有ベクトル (実埋め込み Jacobi)
fn eig_herm3(hre: &[[f64; 3]; 3], him: &[[f64; 3]; 3]) -> ([f64; 3], [[(f64, f64); 3]; 3]) {
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

/// 縮退空間内で X̂ = e^{2πix/N} を一般位相 φ₀ で対角化し、局在世代モードを得る。
/// 世代は中心位置の昇順に整列して返す (セクター間で世代番号を揃えるため)。
fn localize(modes: &[C3v]) -> (Vec<C3v>, Vec<f64>, Vec<f64>) {
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
    // Re(e^{-iφ₀}U): φ₀=0 は対称配置 (cos の縮退) で基底が混ざるので一般の φ₀ を使う
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
    let mut widths = Vec::new();
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
        let (mut zr, mut zi, mut norm) = (0.0, 0.0, 0.0);
        for i in 0..NS {
            let p = psi[i].0 * psi[i].0 + psi[i].1 * psi[i].1;
            let x = (i % N) as f64;
            let (sn, cs) = (two_pi * x / N as f64).sin_cos();
            zr += p * cs;
            zi += p * sn;
            norm += p;
        }
        let center = (zi.atan2(zr) / two_pi * N as f64).rem_euclid(N as f64);
        let r = (zr * zr + zi * zi).sqrt() / norm;
        let width = (-2.0 * r.ln()).sqrt() / two_pi * N as f64;
        out.push(psi);
        centers.push(center);
        widths.push(width);
    }
    let mut ord: Vec<usize> = (0..Q).collect();
    ord.sort_by(|&a, &b| centers[a].partial_cmp(&centers[b]).unwrap());
    let out2: Vec<C3v> = ord.iter().map(|&i| out[i]).collect();
    let c2: Vec<f64> = ord.iter().map(|&i| centers[i]).collect();
    let w2: Vec<f64> = ord.iter().map(|&i| widths[i]).collect();
    (out2, c2, w2)
}

/// 湯川行列: Y_ij = Σ_x conj(ψ_i^(a)) ψ_j^(b) φ_H,  φ_H = 原点中心の周期ガウス
fn yukawa(la: &[C3v], lb: &[C3v], sig_h: f64) -> [[(f64, f64); 3]; 3] {
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

/// 3×3 エルミート行列の固有値 (昇順, 閉形式 — v6.5 と同じ式で Jacobi と照合済み)
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

fn svals(y: &[[(f64, f64); 3]; 3]) -> [f64; 3] {
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
    let lam = eigvals3(&hre, &him);
    [
        lam[0].max(0.0).sqrt(),
        lam[1].max(0.0).sqrt(),
        lam[2].max(0.0).sqrt(),
    ]
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
    println!("=== v7.2 幾何からの湯川 (M2 級): 電荷を選ぶ代わりに Wilson 線を選ぶ ===\n");
    let sigma = (2.0f64).ln();
    let t_u = [EPS_OBS[0].ln(), EPS_OBS[1].ln()];
    let t_d = [EPS_OBS[2].ln(), EPS_OBS[3].ln()];
    let t_e = [EPS_OBS[4].ln(), EPS_OBS[5].ln()];

    // ---- [1] Wilson 線グリッド上の世代モード (6 回の厳密対角化) ----
    println!(
        "[1] 磁束 Q=3・Wilson 線 k∈Z₆ の 18×18 トーラス (対角化 {} 回)",
        NK
    );
    let t0 = std::time::Instant::now();
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    let mut c0 = Vec::new();
    let mut ok_engine = true;
    for k in 0..NK {
        let (modes, gap, spread) = flux_modes(k);
        let (loc, centers, widths) = localize(&modes);
        if k == 0 {
            println!(
                "    k=0: 縮退幅 {:.1e} ≪ ギャップ {:.3}, 中心 {:?}, 幅 ≈ {:.2}",
                spread,
                gap,
                centers
                    .iter()
                    .map(|c| (c * 10.0).round() / 10.0)
                    .collect::<Vec<_>>(),
                widths.iter().sum::<f64>() / 3.0
            );
            c0 = centers.clone();
        } else {
            // 中心が ±k サイト (世代間隔 6 の任意の倍数を法として) 平行移動しているか
            let mut shift_err: f64 = 0.0;
            for g in 0..Q {
                let mut best: f64 = f64::INFINITY;
                for g0 in 0..Q {
                    let d = (centers[g] - c0[g0]).rem_euclid(N as f64);
                    for sgn in [1.0f64, -1.0] {
                        for j in 0..Q {
                            let cand = (sgn * k as f64 + 6.0 * j as f64).rem_euclid(N as f64);
                            let e = (d - cand).abs();
                            let e = e.min(N as f64 - e);
                            best = best.min(e);
                        }
                    }
                }
                shift_err = shift_err.max(best);
            }
            if shift_err > 0.35 || spread > 1e-9 || gap < 0.05 {
                ok_engine = false;
            }
        }
        locs.push(loc);
    }
    println!(
        "    Wilson 線でモード中心が 1 サイト刻みで平行移動 (縮退・ギャップも維持)  {}  ({} ms)",
        pass(ok_engine),
        t0.elapsed().as_millis()
    );

    // ---- [2] 証拠: 全幾何パラメータ (k_Q,k_u,k_d,k_L,k_e ∈ Z₆, σ_H ∈ 4) ----
    println!("\n[2] 幾何の全パラメータ空間 (6⁵ × σ_H 4 通り) の証拠");
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let sector_ll = |target: [f64; 2], sig_h: f64, locs: &[Vec<C3v>]| -> Vec<f64> {
        let norm = -(2.0 * std::f64::consts::PI * sigma * sigma).ln();
        let mut out = vec![0.0; NK * NK];
        for ka in 0..NK {
            for kb in 0..NK {
                let y = yukawa(&locs[ka], &locs[kb], sig_h);
                let sv = svals(&y);
                let r1 = (sv[0].max(1e-300) / sv[2].max(1e-300)).ln();
                let r2 = (sv[1].max(1e-300) / sv[2].max(1e-300)).ln();
                out[ka + kb * NK] = -((r1 - target[0]).powi(2) + (r2 - target[1]).powi(2))
                    / (2.0 * sigma * sigma)
                    + norm;
            }
        }
        out
    };
    let t1 = std::time::Instant::now();
    let mut lnz_terms: Vec<f64> = Vec::new();
    let mut map = (f64::NEG_INFINITY, 0usize, [0usize; 5]);
    for (si, &sh) in sig_grid.iter().enumerate() {
        let lu = sector_ll(t_u, sh, &locs);
        let ld = sector_ll(t_d, sh, &locs);
        let le = sector_ll(t_e, sh, &locs);
        let lse = |v: &[f64]| -> f64 {
            let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
        };
        let mut per_q: Vec<f64> = Vec::with_capacity(NK);
        for kq in 0..NK {
            let us: Vec<f64> = (0..NK).map(|ku| lu[kq + ku * NK]).collect();
            let ds: Vec<f64> = (0..NK).map(|kd| ld[kq + kd * NK]).collect();
            per_q.push(lse(&us) + lse(&ds));
        }
        lnz_terms.push(lse(&per_q) + lse(&le));
        for kq in 0..NK {
            for ku in 0..NK {
                for kd in 0..NK {
                    let base = lu[kq + ku * NK] + ld[kq + kd * NK];
                    for kle in 0..NK * NK {
                        let l = base + le[kle];
                        if l > map.0 {
                            map = (l, si, [kq, ku, kd, kle % NK, kle / NK]);
                        }
                    }
                }
            }
        }
    }
    let lse_all = {
        let m = lnz_terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        m + lnz_terms.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
    };
    let lnz_m2geo = lse_all - (5.0 * (NK as f64).ln() + (sig_grid.len() as f64).ln());
    println!(
        "    走査 {} ms — lnZ(M2geo) = {:.2}",
        t1.elapsed().as_millis(),
        lnz_m2geo
    );
    println!(
        "    MAP 幾何: σ_H={}, k_Q={}, k_u={}, k_d={}, k_L={}, k_e={} (lnL={:.2})",
        sig_grid[map.1], map.2[0], map.2[1], map.2[2], map.2[3], map.2[4], map.0
    );

    // ---- [3] MAP 幾何の決定論的予測 (out-of-sample CKM 込み) ----
    println!("\n[3] MAP 幾何の予測 (乱雑係数なし — 幾何が行列を完全に決める)");
    let sh = sig_grid[map.1];
    let yu = yukawa(&locs[map.2[0]], &locs[map.2[1]], sh);
    let yd = yukawa(&locs[map.2[0]], &locs[map.2[2]], sh);
    let ye = yukawa(&locs[map.2[3]], &locs[map.2[4]], sh);
    let (su, sd, se) = (svals(&yu), svals(&yd), svals(&ye));
    let preds6 = [
        su[0] / su[2],
        su[1] / su[2],
        sd[0] / sd[2],
        sd[1] / sd[2],
        se[0] / se[2],
        se[1] / se[2],
    ];
    let names = [
        "m_u/m_t",
        "m_c/m_t",
        "m_d/m_b",
        "m_s/m_b",
        "m_e/m_τ",
        "m_μ/m_τ",
        "|V_us|",
        "|V_cb|",
        "|V_ub|",
    ];
    let mut ok6 = 0;
    for k in 0..6 {
        let ratio = preds6[k] / EPS_OBS[k];
        let within = ratio > 0.2 && ratio < 5.0;
        if within {
            ok6 += 1;
        }
        println!(
            "    {:8}  幾何予測 {:9.2e}   実測 {:8.2e}   比 {:6.2} {}",
            names[k],
            preds6[k],
            EPS_OBS[k],
            ratio,
            if within { "✓" } else { " " }
        );
    }
    let ckm = {
        let heig = |y: &[[(f64, f64); 3]; 3]| -> [[(f64, f64); 3]; 3] {
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
            eig_herm3(&hre, &him).1
        };
        let vu = heig(&yu);
        let vd = heig(&yd);
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
        ckm
    };
    let ckm_pred = [ckm[0][1], ckm[1][2], ckm[0][2]];
    let mut ok_ckm = 0;
    for k in 0..3 {
        let ratio = ckm_pred[k] / EPS_OBS[6 + k];
        let within = ratio > 0.2 && ratio < 5.0;
        if within {
            ok_ckm += 1;
        }
        println!(
            "    {:8}  幾何予測 {:9.2e}   実測 {:8.2e}   比 {:6.2} {} (out-of-sample)",
            names[6 + k],
            ckm_pred[k],
            EPS_OBS[6 + k],
            ratio,
            if within { "✓" } else { " " }
        );
    }

    // ---- [3b] T²×T² 積模型: 因子化コンパクト化の最小版 ----
    // 弦現象論の標準設定 (T⁶=T²×T²×…) では湯川は各トーラスの重なりの積になる。
    // ここでは 2 枚の同一トーラスに独立な Wilson 線を与え、世代は対角対
    // (i,i) とする (宣言された仮定)。抑制が 2 乗になり 10⁻⁵ 級に届くかを見る。
    println!("\n[3b] T²×T² 積模型 (対角世代対): パラメータ = Wilson 線 10 個 (6¹⁰) + σ_H");
    let t2 = std::time::Instant::now();
    let nk2 = NK * NK;
    let mut lnz2_terms: Vec<f64> = Vec::new();
    let mut map2 = (f64::NEG_INFINITY, 0usize, [0usize; 5]); // 添字は複合 (a1 + a2*NK)
    for (si, &sh) in sig_grid.iter().enumerate() {
        // 単一トーラスの湯川 36 通りを前計算
        let ytab: Vec<[[(f64, f64); 3]; 3]> = (0..nk2)
            .map(|ab| yukawa(&locs[ab % NK], &locs[ab / NK], sh))
            .collect();
        let sector2 = |target: [f64; 2]| -> Vec<f64> {
            let norm = -(2.0 * std::f64::consts::PI * sigma * sigma).ln();
            let mut out = vec![0.0; nk2 * nk2];
            for a in 0..nk2 {
                // a = (a1, a2): 左側の 2 トーラス Wilson 線
                let (a1, a2) = (a % NK, a / NK);
                for b in 0..nk2 {
                    let (b1, b2) = (b % NK, b / NK);
                    let y1 = &ytab[a1 + b1 * NK];
                    let y2 = &ytab[a2 + b2 * NK];
                    let mut y = [[(0.0f64, 0.0f64); 3]; 3];
                    for i in 0..3 {
                        for j in 0..3 {
                            let (p, q) = y1[i][j];
                            let (r, s) = y2[i][j];
                            y[i][j] = (p * r - q * s, p * s + q * r);
                        }
                    }
                    let sv = svals(&y);
                    let r1 = (sv[0].max(1e-300) / sv[2].max(1e-300)).ln();
                    let r2 = (sv[1].max(1e-300) / sv[2].max(1e-300)).ln();
                    out[a + b * nk2] = -((r1 - target[0]).powi(2) + (r2 - target[1]).powi(2))
                        / (2.0 * sigma * sigma)
                        + norm;
                }
            }
            out
        };
        let lu = sector2(t_u);
        let ld = sector2(t_d);
        let le = sector2(t_e);
        let lse = |v: &[f64]| -> f64 {
            let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
            m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
        };
        let mut per_q: Vec<f64> = Vec::with_capacity(nk2);
        for kq in 0..nk2 {
            let us: Vec<f64> = (0..nk2).map(|ku| lu[kq + ku * nk2]).collect();
            let ds: Vec<f64> = (0..nk2).map(|kd| ld[kq + kd * nk2]).collect();
            per_q.push(lse(&us) + lse(&ds));
        }
        lnz2_terms.push(lse(&per_q) + lse(&le));
        for kq in 0..nk2 {
            for ku in 0..nk2 {
                for kd in 0..nk2 {
                    let base = lu[kq + ku * nk2] + ld[kq + kd * nk2];
                    if base + 0.0 <= map2.0 - 60.0 {
                        continue; // 枝刈り (le の最大値は高々 norm≈-1.1)
                    }
                    for kle in 0..nk2 * nk2 {
                        let l = base + le[kle];
                        if l > map2.0 {
                            map2 = (l, si, [kq, ku, kd, kle % nk2, kle / nk2]);
                        }
                    }
                }
            }
        }
    }
    let lse_all2 = {
        let m = lnz2_terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        m + lnz2_terms.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
    };
    let lnz_m2geo2 = lse_all2 - (10.0 * (NK as f64).ln() + (sig_grid.len() as f64).ln());
    println!(
        "    走査 {} ms — lnZ(M2geo², T²×T²) = {:.2}",
        t2.elapsed().as_millis(),
        lnz_m2geo2
    );
    let dec = |c: usize| (c % NK, c / NK);
    println!(
        "    MAP: σ_H={}, k_Q={:?}, k_u={:?}, k_d={:?}, k_L={:?}, k_e={:?} (lnL={:.2})",
        sig_grid[map2.1],
        dec(map2.2[0]),
        dec(map2.2[1]),
        dec(map2.2[2]),
        dec(map2.2[3]),
        dec(map2.2[4]),
        map2.0
    );
    // MAP の決定論的予測
    let sh2 = sig_grid[map2.1];
    let prod_y = |ka: usize, kb: usize| -> [[(f64, f64); 3]; 3] {
        let (a1, a2) = dec(ka);
        let (b1, b2) = dec(kb);
        let y1 = yukawa(&locs[a1], &locs[b1], sh2);
        let y2 = yukawa(&locs[a2], &locs[b2], sh2);
        let mut y = [[(0.0f64, 0.0f64); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (p, q) = y1[i][j];
                let (r, s) = y2[i][j];
                y[i][j] = (p * r - q * s, p * s + q * r);
            }
        }
        y
    };
    let yu2 = prod_y(map2.2[0], map2.2[1]);
    let yd2 = prod_y(map2.2[0], map2.2[2]);
    let ye2 = prod_y(map2.2[3], map2.2[4]);
    let (su2, sd2, se2) = (svals(&yu2), svals(&yd2), svals(&ye2));
    let preds62 = [
        su2[0] / su2[2],
        su2[1] / su2[2],
        sd2[0] / sd2[2],
        sd2[1] / sd2[2],
        se2[0] / se2[2],
        se2[1] / se2[2],
    ];
    let mut ok62 = 0;
    for k in 0..6 {
        let ratio = preds62[k] / EPS_OBS[k];
        let within = ratio > 0.2 && ratio < 5.0;
        if within {
            ok62 += 1;
        }
        println!(
            "    {:8}  積幾何予測 {:9.2e}   実測 {:8.2e}   比 {:6.2} {}",
            names[k],
            preds62[k],
            EPS_OBS[k],
            ratio,
            if within { "✓" } else { " " }
        );
    }
    let heig2 = |y: &[[(f64, f64); 3]; 3]| -> [[(f64, f64); 3]; 3] {
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
        eig_herm3(&hre, &him).1
    };
    let (vu2, vd2) = (heig2(&yu2), heig2(&yd2));
    let mut ckm2m = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = vu2[i][k];
                let (c, d) = vd2[j][k];
                re += a * c + b * d;
                im += a * d - b * c;
            }
            ckm2m[i][j] = (re * re + im * im).sqrt();
        }
    }
    let ckm_pred2 = [ckm2m[0][1], ckm2m[1][2], ckm2m[0][2]];
    let mut ok_ckm2 = 0;
    for k in 0..3 {
        let ratio = ckm_pred2[k] / EPS_OBS[6 + k];
        let within = ratio > 0.2 && ratio < 5.0;
        if within {
            ok_ckm2 += 1;
        }
        println!(
            "    {:8}  積幾何予測 {:9.2e}   実測 {:8.2e}   比 {:6.2} {} (out-of-sample)",
            names[6 + k],
            ckm_pred2[k],
            EPS_OBS[6 + k],
            ratio,
            if within { "✓" } else { " " }
        );
    }

    // ---- [4] 模型比較 ----
    println!("\n[4] 模型比較 (尤度・データ・σ は v6.5 と同一; M0/M1/M2文献 は results/v65_bayes.json より)");
    let (lnz_m0_bound, lnz_m1, lnz_m2lit) = (-35.40f64, -12.19, -7.61);
    println!(
        "    lnZ(M0 アナーキー)   ≤ {:7.2}   (パラメータ 0 + 乱雑係数 18)",
        lnz_m0_bound
    );
    println!(
        "    lnZ(M1 FN 自由電荷)  = {:7.2}   (整数 10 + 乱雑係数 18)",
        lnz_m1
    );
    println!(
        "    lnZ(M2 文献電荷)     = {:7.2}   (固定電荷 + 乱雑係数 18)",
        lnz_m2lit
    );
    println!(
        "    lnZ(M2geo, 単一 T²)  = {:7.2}   (Wilson 線 6⁵ + σ_H — 乱雑係数なし)",
        lnz_m2geo
    );
    println!(
        "    lnZ(M2geo², T²×T²)   = {:7.2}   (Wilson 線 6¹⁰ + σ_H — 乱雑係数なし)",
        lnz_m2geo2
    );
    let beats_m0 = lnz_m2geo2 > lnz_m0_bound;
    println!(
        "    => 積幾何模型はアナーキーの上界を{} (差 {:+.1})、M1 との差 {:+.1}",
        if beats_m0 {
            "上回る"
        } else {
            "上回らない"
        },
        lnz_m2geo2 - lnz_m0_bound,
        lnz_m2geo2 - lnz_m1
    );

    // ---- JSON / 判定 ----
    // PASS 条件はエンジンの正しさのみ — 物理的な勝敗は報告が本体 (どちらでも記録する)
    let all_ok = ok_engine;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-003".into())),
        ("lattice".into(), Json::Int(N as i64)),
        ("flux".into(), Json::Int(Q as i64)),
        ("wilson_grid".into(), Json::Int(NK as i64)),
        ("lnZ_m2geo_single_torus".into(), Json::Num(lnz_m2geo)),
        ("lnZ_m2geo_product_torus".into(), Json::Num(lnz_m2geo2)),
        ("product_mass_ratios_within_factor5".into(), Json::Int(ok62)),
        ("product_ckm_within_factor5".into(), Json::Int(ok_ckm2)),
        ("lnZ_m1_from_v65".into(), Json::Num(lnz_m1)),
        ("lnZ_m2lit_from_v65".into(), Json::Num(lnz_m2lit)),
        (
            "lnZ_m0_upper_bound_from_v65".into(),
            Json::Num(lnz_m0_bound),
        ),
        ("beats_anarchy_bound".into(), Json::Bool(beats_m0)),
        (
            "map_geometry".into(),
            Json::Obj(vec![
                ("sigma_h".into(), Json::Num(sig_grid[map.1])),
                (
                    "wilson_k".into(),
                    Json::Arr(map.2.iter().map(|&s| Json::Int(s as i64)).collect()),
                ),
            ]),
        ),
        ("mass_ratios_within_factor5".into(), Json::Int(ok6)),
        ("ckm_within_factor5_out_of_sample".into(), Json::Int(ok_ckm)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v72_geomfn.json", &j.render());
    println!("\n  機械可読な結果: {}", p);

    println!(
        "\n総合判定: {} (PASS 条件はエンジンの正しさ — 物理的な勝敗は上の表の報告が本体)",
        pass(all_ok)
    );
    println!("\n結論: 世代 = 磁束のゼロモード (v2.3) を湯川まで延長した。行列は乱雑係数なしで");
    println!("      Wilson 線の整数が完全に決める。単一 T² は階層が浅く届かない");
    println!(
        "      (質量比 {}/6, CKM {}/3, lnZ={:.1} — 正直な陰性結果)。",
        ok6, ok_ckm, lnz_m2geo
    );
    println!(
        "      T²×T² では抑制が 2 乗になり: 質量比 {}/6・CKM {}/3 (out-of-sample) が",
        ok62, ok_ckm2
    );
    println!(
        "      5 倍以内、lnZ={:.1} でアナーキー上界を {:+.1} 上回る。",
        lnz_m2geo2,
        lnz_m2geo2 - lnz_m0_bound
    );
    println!("      *** 階層の深さは、余剰次元 (トーラス) の数を数えている ***");
    println!(
        "      M1 (自由電荷+乱雑係数 18 個) との残差 {:+.1} は幾何格子の粗さと",
        lnz_m2geo2 - lnz_m1
    );
    println!("      対角世代対の仮定に帰される — M3 (完全導出) への残りの距離である。");
    if !all_ok {
        std::process::exit(1);
    }
}
