//! v19.5 3+1D QNEC — 格子ヌル線に沿う S″ ≤ 2π⟨T_kk⟩ (Einstein 橋の第四段)
//!
//! 領域 A(λ) = {x ≤ n/2−1+λ} を時刻 t = λ と同時に伸ばす (格子ヌル変形 Δx=Δt=1)。
//! 純粋状態なので S(A) = S(Aᶜ) — λ=0 で両側を計算して一致をゲート化 (純粋性+機構検証)、
//! λ≥1 は補集合で評価 (縮む側 — 計算量が激減)。励起は波束 φ = G e^{ikx·x} の
//! 占有/非占有分割 (ũ, w̃) の完全反転 C̃′ = diag(occ) + ũ*ũᵀ − w̃*w̃ᵀ — rank-2 は
//! 時間発展でも rank-2 のまま (U(t) = V e^{−iεt}ũ, O(ns²))。モジュラー核は一切使わない。
//! T_kk = T_tt + 2T_tx + T_xx: T_tt/T_xx は面のボンドエネルギー、T_tx はエネルギー流
//! −dE_A/dt (中心差分) — 全成分をエネルギー保存 (厳密) に規約固定。
//! S″ の一次推定量は真空差引 ΔS(λ) = S(λ) − S_vac(λ) の曲率 (格子の薄板・staggered
//! 偶奇系統が差分で消える — run1 で確認)。S(λ) フィットは a + bλ + cλ² + d(−1)^λ。
//! 事前登録: N=12 の 2 構成 (静止 kx=π/2 / 移動 kx=π/4) で窓平均 Q_Δ = 2πT̄_kk − S″_Δ:
//! (a) 両構成 Q_Δ ≥ 0 (ε なし) = 第四段成立 / (a′) 静止のみ / (b) 静止で破れ。
//! 開発記録 (run1, results/v195_qnec3d_run1.txt): 電流の向きが未較正 (|ΔN|=|∫j| 5 桁一致・
//! 符号逆) + 「+x 伝播」ゲートは粒子・正孔対でエネルギー重心が静止する物理を誤解。
//! v2 較正: 電流の向き固定・判別量は対の相対分離 (移動 0.495 vs 静止 0.007 — 共動
//! ドリフトを除いた形)。物理数値 (S, T_kk) は run1 と不変。

use uft_sim::*;

/// 3+1D staggered ハミルトニアン (実対称, x 開放・y,z 反周期)
fn build_h3d(n: usize) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let mut add = |i: usize, j: usize, v: f64| {
        h[j + i * ns] += v;
        h[i + j * ns] += v;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                if x + 1 < n {
                    add(i, idx(x + 1, y, z), 0.5);
                }
                let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                if y + 1 < n {
                    add(i, idx(x, y + 1, z), ey);
                } else {
                    add(i, idx(x, 0, z), -ey);
                }
                let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                if z + 1 < n {
                    add(i, idx(x, y, z + 1), ez);
                } else {
                    add(i, idx(x, y, 0), -ez);
                }
            }
        }
    }
    h
}

/// フェルミオン・エントロピー (実対称 C の固有値から)
fn entropy_sym(c: &[f64], m: usize) -> f64 {
    let (w, _v) = jacobi_eigh(c, m);
    let mut s = 0.0;
    for &ck in &w {
        let x = ck.clamp(1e-14, 1.0 - 1e-14);
        s += -x * x.ln() - (1.0 - x) * (1.0 - x).ln();
    }
    s
}

/// 領域 A の C_A
fn restrict(c: &[f64], ns: usize, sel: &[usize]) -> Vec<f64> {
    let m = sel.len();
    let mut ca = vec![0.0f64; m * m];
    for (i, &si) in sel.iter().enumerate() {
        for (j, &sj) in sel.iter().enumerate() {
            ca[j + i * m] = c[sj + si * ns];
        }
    }
    ca
}

// 4×4 正規方程式の Gauss 消去 (部分ピボット)
fn gauss4(a: &mut [f64; 16], y: &mut [f64; 4]) -> [f64; 4] {
    for col in 0..4 {
        let mut piv = col;
        for r in col + 1..4 {
            if a[r * 4 + col].abs() > a[piv * 4 + col].abs() {
                piv = r;
            }
        }
        if piv != col {
            for c in 0..4 {
                a.swap(col * 4 + c, piv * 4 + c);
            }
            y.swap(col, piv);
        }
        let d = a[col * 4 + col];
        for r in col + 1..4 {
            let f = a[r * 4 + col] / d;
            for c in col..4 {
                a[r * 4 + c] -= f * a[col * 4 + c];
            }
            y[r] -= f * y[col];
        }
    }
    let mut b = [0.0f64; 4];
    for col in (0..4).rev() {
        let mut acc = y[col];
        for c in col + 1..4 {
            acc -= a[col * 4 + c] * b[c];
        }
        b[col] = acc / a[col * 4 + col];
    }
    b
}

// S(λ) = a + bλ + cλ² + d(−1)^λ の最小二乗 → (S″ = 2c, 残差 rms)
fn fit_s(lams: &[f64], ys: &[f64]) -> (f64, f64) {
    let basis =
        |l: f64| -> [f64; 4] { [1.0, l, l * l, if (l as i64) % 2 == 0 { 1.0 } else { -1.0 }] };
    let mut ata = [0.0f64; 16];
    let mut aty = [0.0f64; 4];
    for (i, &l) in lams.iter().enumerate() {
        let ph = basis(l);
        for r in 0..4 {
            for c in 0..4 {
                ata[r * 4 + c] += ph[r] * ph[c];
            }
            aty[r] += ph[r] * ys[i];
        }
    }
    let beta = gauss4(&mut ata, &mut aty);
    let mut ss = 0.0;
    for (i, &l) in lams.iter().enumerate() {
        let ph = basis(l);
        let pred: f64 = (0..4).map(|r| beta[r] * ph[r]).sum();
        ss += (ys[i] - pred).powi(2);
    }
    (2.0 * beta[2], (ss / lams.len() as f64).sqrt())
}

fn main() {
    self_test();
    println!("=== v19.5 3+1D QNEC: 格子ヌル線に沿う S″ ≤ 2π⟨T_kk⟩ (Einstein 橋の第四段) ===\n");
    println!("事前登録: 窓平均形 Q_Δ = 2π T̄_kk − S″_Δ を N=12 の 2 構成 (静止 kx=π/2・移動 kx=π/4) で測る。");
    println!(
        "  S″_Δ は真空差引 ΔS(λ) = S(λ) − S_vac(λ) の曲率 — 格子の薄板・偶奇系統は差分で消える"
    );
    println!("  (run1 で確認: ΔS 列は滑らか・S_vac 単独は S″_vac ~ +0.76 の格子系統を持つ)。");
    println!("  (a) 両構成で Q_Δ ≥ 0 (ε なし) = 3+1D QNEC (平均形) 成立 — 第四段成立 /");
    println!("  (a′) 静止のみ成立 = 移動構成の逸脱を記録 / (b) 静止で破れ = 破れの構造を記録。");
    println!("  T_tx はエネルギー流 −dE_A/dt (中心差分) — エネルギー保存に規約を固定。");
    println!(
        "  N=8 は較正記録。装置ゲート: 補集合エントロピー一致 (純粋性)・E_bond=E_eig (厳密)・"
    );
    println!("  エネルギー保存・粒子連続の式 (電流の向きは run1 の連続の式で較正済み)・対の分離\n");
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
    let two_pi = 2.0 * std::f64::consts::PI;

    // (n, config, s2d, tkk_mean, qd, s2raw, svac2)
    let mut summary: Vec<(usize, String, f64, f64, f64, f64, f64)> = Vec::new();
    for &n in &[8usize, 12] {
        let ns = n * n * n;
        let lmax = if n == 8 { 3usize } else { 4 };
        let t0 = std::time::Instant::now();
        let h = build_h3d(n);
        let (eps, v) = jacobi_eigh(&h, ns);
        let nocc = ns / 2;
        let gap = eps[nocc] - eps[nocc - 1];
        check(
            &format!("N={} の閉殻ギャップ", n),
            gap > 1e-6,
            format!("gap = {:.4} ({} s)", gap, t0.elapsed().as_secs()),
        );
        let mut c0 = vec![0.0f64; ns * ns];
        for k in 0..nocc {
            for i in 0..ns {
                let vi = v[i + k * ns];
                if vi == 0.0 {
                    continue;
                }
                for j in 0..ns {
                    c0[j + i * ns] += vi * v[j + k * ns];
                }
            }
        }
        let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
        // 全ボンド一覧: (i, j, t, xmid)。x ボンドは xmid = x+0.5、y/z ボンドは xmid = x
        let mut bonds: Vec<(usize, usize, f64, f64)> = Vec::new();
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let i = idx(x, y, z);
                    if x + 1 < n {
                        bonds.push((i, idx(x + 1, y, z), 0.5, x as f64 + 0.5));
                    }
                    let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                    if y + 1 < n {
                        bonds.push((i, idx(x, y + 1, z), ey, x as f64));
                    } else {
                        bonds.push((i, idx(x, 0, z), -ey, x as f64));
                    }
                    let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                    if z + 1 < n {
                        bonds.push((i, idx(x, y, z + 1), ez, x as f64));
                    } else {
                        bonds.push((i, idx(x, y, 0), -ez, x as f64));
                    }
                }
            }
        }
        // 領域 A(λ) = {x ≤ n/2−1+λ} の補集合
        let comp = |lam: usize| -> Vec<usize> {
            let mut sel = Vec::new();
            for z in 0..n {
                for y in 0..n {
                    for x in n / 2 + lam..n {
                        sel.push(idx(x, y, z));
                    }
                }
            }
            sel
        };
        let direct = |lam: usize| -> Vec<usize> {
            let mut sel = Vec::new();
            for z in 0..n {
                for y in 0..n {
                    for x in 0..n / 2 + lam {
                        sel.push(idx(x, y, z));
                    }
                }
            }
            sel
        };
        // ---- 真空基線 S_vac(λ) (静的・実対称) ----
        let mut svac = Vec::new();
        for lam in 0..=lmax {
            let sel = comp(lam);
            let ca = restrict(&c0, ns, &sel);
            svac.push(entropy_sym(&ca, sel.len()));
        }
        {
            let sel = direct(0);
            let ca = restrict(&c0, ns, &sel);
            let sd = entropy_sym(&ca, sel.len());
            check(
                &format!("N={} 真空の補集合一致 S(A)=S(Aᶜ)", n),
                (sd - svac[0]).abs() < 1e-6,
                format!("|ΔS| = {:.2e}", (sd - svac[0]).abs()),
            );
        }
        let lams: Vec<f64> = (0..=lmax).map(|l| l as f64).collect();
        let (svac2, svac_res) = fit_s(&lams, &svac);
        println!(
            "    N={} 真空: S_vac(λ) = {:?} → S″_vac = {:+.5e} (残差 rms {:.1e}) — 格子系統の目盛り",
            n,
            svac.iter().map(|s| (s * 1e4).round() / 1e4).collect::<Vec<_>>(),
            svac2,
            svac_res
        );
        // ---- 励起構成 ----
        let sig = n as f64 / 8.0;
        let x0 = n as f64 / 2.0 - 0.5;
        let mid = n as f64 / 2.0;
        for &(kx, tag) in &[
            (std::f64::consts::PI / 2.0, "静止 kx=π/2"),
            (std::f64::consts::PI / 4.0, "移動 kx=π/4"),
        ] {
            // 複素波束 φ = G(r−r₀) e^{i kx x} → 固有基底で占有/非占有に分割
            let mut phr = vec![0.0f64; ns];
            let mut phi = vec![0.0f64; ns];
            for x in 0..n {
                for y in 0..n {
                    for z in 0..n {
                        let d2 = (x as f64 - x0).powi(2)
                            + (y as f64 - mid).powi(2)
                            + (z as f64 - mid).powi(2);
                        let g = (-d2 / (2.0 * sig * sig)).exp();
                        phr[idx(x, y, z)] = g * (kx * x as f64).cos();
                        phi[idx(x, y, z)] = g * (kx * x as f64).sin();
                    }
                }
            }
            // φ̃ = Vᵀ φ
            let mut ptr = vec![0.0f64; ns];
            let mut pti = vec![0.0f64; ns];
            for k in 0..ns {
                let mut ar = 0.0;
                let mut ai = 0.0;
                for i in 0..ns {
                    let vv = v[i + k * ns];
                    ar += vv * phr[i];
                    ai += vv * phi[i];
                }
                ptr[k] = ar;
                pti[k] = ai;
            }
            let nrm_u: f64 = (nocc..ns).map(|k| ptr[k] * ptr[k] + pti[k] * pti[k]).sum();
            let nrm_w: f64 = (0..nocc).map(|k| ptr[k] * ptr[k] + pti[k] * pti[k]).sum();
            let tot = nrm_u + nrm_w;
            check(
                &format!("N={} {} の占有/非占有分割が非退化", n, tag),
                nrm_u / tot > 0.15 && nrm_u / tot < 0.85,
                format!("非占有比 = {:.3}", nrm_u / tot),
            );
            let (su, sw) = (nrm_u.sqrt(), nrm_w.sqrt());
            let mut ut = vec![(0.0f64, 0.0f64); ns]; // ũ (非占有部, 正規化)
            let mut wt = vec![(0.0f64, 0.0f64); ns]; // w̃ (占有部, 正規化)
            for k in 0..ns {
                if k >= nocc {
                    ut[k] = (ptr[k] / su, pti[k] / su);
                } else {
                    wt[k] = (ptr[k] / sw, pti[k] / sw);
                }
            }
            // ΔE (固有基底の厳密値): Σ ε (|ũ|² − |w̃|²)
            let de_eig: f64 = (0..ns)
                .map(|k| {
                    eps[k]
                        * (ut[k].0 * ut[k].0 + ut[k].1 * ut[k].1
                            - wt[k].0 * wt[k].0
                            - wt[k].1 * wt[k].1)
                })
                .sum();
            // U(t) = V e^{−iεt} ũ, W(t) = V e^{−iεt} w̃ (サイト基底, O(ns²))
            let evolve = |t: f64| -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
                let mut ur = vec![0.0f64; ns];
                let mut ui = vec![0.0f64; ns];
                let mut wr = vec![0.0f64; ns];
                let mut wi = vec![0.0f64; ns];
                for k in 0..ns {
                    let (cs, sn) = ((eps[k] * t).cos(), -(eps[k] * t).sin());
                    let (ar, ai) = (ut[k].0 * cs - ut[k].1 * sn, ut[k].0 * sn + ut[k].1 * cs);
                    let (br, bi) = (wt[k].0 * cs - wt[k].1 * sn, wt[k].0 * sn + wt[k].1 * cs);
                    if ar == 0.0 && ai == 0.0 && br == 0.0 && bi == 0.0 {
                        continue;
                    }
                    for i in 0..ns {
                        let vv = v[i + k * ns];
                        ur[i] += vv * ar;
                        ui[i] += vv * ai;
                        wr[i] += vv * br;
                        wi[i] += vv * bi;
                    }
                }
                (ur, ui, wr, wi)
            };
            // ΔC = U*Uᵀ − W*Wᵀ: dcre(i,j) = Re(U*_i U_j) − …, dcim(i,j) = Im(U*_i U_j) − …
            let dcre = |q: &(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), i: usize, j: usize| -> f64 {
                q.0[i] * q.0[j] + q.1[i] * q.1[j] - (q.2[i] * q.2[j] + q.3[i] * q.3[j])
            };
            let dcim = |q: &(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), i: usize, j: usize| -> f64 {
                q.0[i] * q.1[j] - q.1[i] * q.0[j] - (q.2[i] * q.3[j] - q.3[i] * q.2[j])
            };
            // ボンド読み出し: 全励起エネルギー・E_A(境界 xb, 跨ぎ半分)・跨ぎボンド量
            let e_tot = |q: &(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>)| -> f64 {
                bonds
                    .iter()
                    .map(|&(i, j, t, _)| 2.0 * t * dcre(q, i, j))
                    .sum()
            };
            let e_a = |q: &(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), xb: f64| -> f64 {
                bonds
                    .iter()
                    .map(|&(i, j, t, xm)| {
                        let w = if xm < xb - 0.01 {
                            1.0
                        } else if xm < xb + 0.01 {
                            0.5
                        } else {
                            0.0
                        };
                        w * 2.0 * t * dcre(q, i, j)
                    })
                    .sum()
            };
            let cross = |q: &(Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>), xb: f64| -> (f64, f64) {
                // (E_cross, E_trans±)
                let mut ec = 0.0;
                let mut et = 0.0;
                for &(i, j, t, xm) in &bonds {
                    if (xm - xb).abs() < 0.01 {
                        ec += 2.0 * t * dcre(q, i, j);
                    } else if (xm - (xb - 0.5)).abs() < 0.01 || (xm - (xb + 0.5)).abs() < 0.01 {
                        et += 2.0 * t * dcre(q, i, j);
                    }
                }
                (ec, et)
            };
            // ---- λ 掃引 ----
            let q0 = evolve(0.0);
            let e0 = e_tot(&q0);
            check(
                &format!("N={} {} の E_bond = E_eig (厳密)", n, tag),
                (e0 - de_eig).abs() < 1e-9 * (1.0 + de_eig.abs()),
                format!("E_bond = {:.9}, E_eig = {:.9}", e0, de_eig),
            );
            let mut s_exc = Vec::new();
            let mut ds_exc = Vec::new();
            let mut tkk = Vec::new();
            let mut elast = 0.0;
            for lam in 0..=lmax {
                let t = lam as f64;
                let q = evolve(t);
                elast = e_tot(&q);
                let xb = n as f64 / 2.0 - 0.5 + t;
                let (ec, et) = cross(&q, xb);
                // T_tx = −dE_A/dt (中心差分 h=0.25)
                let qp = evolve(t + 0.25);
                let qm = evolve(t - 0.25);
                let ttx = -(e_a(&qp, xb) - e_a(&qm, xb)) / 0.5;
                let t_tt = ec + 0.5 * et;
                let t_kk = t_tt + 2.0 * ttx + ec;
                tkk.push(t_kk);
                // エントロピー (補集合; λ=0 は直接側 too = 純粋性ゲート)
                let sel = comp(lam);
                let m = sel.len();
                let mut cre = vec![0.0f64; m * m];
                let mut cim = vec![0.0f64; m * m];
                for a in 0..m {
                    for b in 0..m {
                        cre[a + b * m] = c0[sel[a] + sel[b] * ns] + dcre(&q, sel[a], sel[b]);
                        cim[a + b * m] = dcim(&q, sel[a], sel[b]);
                    }
                }
                let s = entropy_corr_herm(&cre, &cim, m);
                if lam == 0 {
                    let seld = direct(0);
                    let md = seld.len();
                    let mut cred = vec![0.0f64; md * md];
                    let mut cimd = vec![0.0f64; md * md];
                    for a in 0..md {
                        for b in 0..md {
                            cred[a + b * md] =
                                c0[seld[a] + seld[b] * ns] + dcre(&q, seld[a], seld[b]);
                            cimd[a + b * md] = dcim(&q, seld[a], seld[b]);
                        }
                    }
                    let sd = entropy_corr_herm(&cred, &cimd, md);
                    check(
                        &format!("N={} {} の補集合一致 (純粋性+機構)", n, tag),
                        (sd - s).abs() < 1e-6,
                        format!("|ΔS| = {:.2e}", (sd - s).abs()),
                    );
                }
                println!(
                    "    N={} {} λ={}: S = {:.6} (ΔS_vac {:+.5}), 2πT_kk = {:+.5} ({} s)",
                    n,
                    tag,
                    lam,
                    s,
                    s - svac[lam],
                    two_pi * t_kk,
                    t0.elapsed().as_secs()
                );
                s_exc.push(s);
                ds_exc.push(s - svac[lam]);
            }
            check(
                &format!("N={} {} のエネルギー保存", n, tag),
                (elast - e0).abs() < 1e-9 * (1.0 + e0.abs()),
                format!("|ΔE| = {:.2e}", (elast - e0).abs()),
            );
            // 粒子連続の式 (固定領域 A(0), t ∈ [0,2], Simpson)。
            // 電流の向きは run1 の連続の式で較正: j(+x) = −Σ 2t·dcim (跨ぎボンド)
            {
                let a0 = direct(0);
                let xb = n as f64 / 2.0 - 0.5;
                let mut nvals = Vec::new();
                let mut jvals = Vec::new();
                for st in 0..5 {
                    let t = st as f64 * 0.5;
                    let q = evolve(t);
                    let na: f64 = a0.iter().map(|&i| dcre(&q, i, i)).sum();
                    let mut jw = 0.0;
                    for &(i, j, t_b, xm) in &bonds {
                        if (xm - xb).abs() < 0.01 {
                            jw -= 2.0 * t_b * dcim(&q, i, j);
                        }
                    }
                    nvals.push(na);
                    jvals.push(jw);
                }
                let dn = nvals[4] - nvals[0];
                let integ = 0.5 / 3.0
                    * (jvals[0] + 4.0 * jvals[1] + 2.0 * jvals[2] + 4.0 * jvals[3] + jvals[4]);
                check(
                    &format!("N={} {} の粒子連続の式 (Simpson)", n, tag),
                    (dn + integ).abs() < (0.02 * dn.abs()).max(2e-3),
                    format!("ΔN_A = {:+.5}, −∫j dt = {:+.5}", dn, -integ),
                );
            }
            // 対の分離: 粒子 U と正孔 W の密度重心が離れること (向きは記録 — 粒子・正孔対は
            // エネルギー重心が静止するのが正しい [run1 の教訓])
            {
                let q2 = evolve(2.0);
                let xden = |r: &[f64], im: &[f64]| -> f64 {
                    let mut sx = 0.0;
                    let mut sn = 0.0;
                    for x in 0..n {
                        for y in 0..n {
                            for z in 0..n {
                                let i = idx(x, y, z);
                                let d = r[i] * r[i] + im[i] * im[i];
                                sx += x as f64 * d;
                                sn += d;
                            }
                        }
                    }
                    sx / sn
                };
                let (xu0, xw0) = (xden(&q0.0, &q0.1), xden(&q0.2, &q0.3));
                let (xu2, xw2) = (xden(&q2.0, &q2.1), xden(&q2.2, &q2.3));
                // 判別量は相対分離 (共動ドリフトを除く): run1/v2 初走の較正 —
                // 移動 0.495 vs 静止 0.007 (70 倍) — 絶対分離は共動ドリフトで汚染される
                let sep = ((xu2 - xw2) - (xu0 - xw0)).abs();
                let need = kx < 1.0; // 移動構成のみゲート (静止は記録)
                let ok = !need || sep > 0.25;
                check(
                    &format!("N={} {} の対の相対分離 (移動構成のみゲート)", n, tag),
                    ok,
                    format!(
                        "x_U: {:.3}→{:.3}, x_W: {:.3}→{:.3} (相対分離 {:.3})",
                        xu0, xu2, xw0, xw2, sep
                    ),
                );
            }
            let (s2raw, res_raw) = fit_s(&lams, &s_exc);
            let (s2d, res_d) = fit_s(&lams, &ds_exc);
            let tbar = tkk.iter().sum::<f64>() / tkk.len() as f64;
            let qd = two_pi * tbar - s2d;
            println!(
                "    N={} {}: S″_Δ = {:+.5e} (残差 {:.1e}), S″_raw = {:+.5e} (残差 {:.1e}), 2πT̄_kk = {:+.5e} → Q_Δ = {:+.5e}",
                n, tag, s2d, res_d, s2raw, res_raw, two_pi * tbar, qd
            );
            summary.push((n, tag.to_string(), s2d, tbar, qd, s2raw, svac2));
        }
    }

    // ---- 判定 (記録) ----
    let n12: Vec<&(usize, String, f64, f64, f64, f64, f64)> =
        summary.iter().filter(|r| r.0 == 12).collect();
    let stand_ok = n12.iter().any(|r| r.1.contains("静止") && r.4 >= 0.0);
    let move_ok = n12.iter().any(|r| r.1.contains("移動") && r.4 >= 0.0);
    println!("\n[判定] N=12 の窓平均 QNEC (真空差引形, ε なし):");
    for r in &n12 {
        println!(
            "    {}: Q_Δ = {:+.5e} {}",
            r.1,
            r.4,
            if r.4 >= 0.0 { "≥ 0 ✓" } else { "< 0" }
        );
    }
    println!(
        "    => {}",
        if stand_ok && move_ok {
            "事前登録 (a): 3+1D QNEC (平均形) が両構成で成立 — 第四段成立"
        } else if stand_ok {
            "事前登録 (a′): 静止構成で成立 — 移動構成の逸脱を記録"
        } else {
            "事前登録 (b): 静止構成で破れ — 構造の記録"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v19.5".into())),
        (
            "rows".into(),
            Json::Arr(
                summary
                    .iter()
                    .map(|(n, tag, s2d, tbar, qd, s2raw, svac2)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(*n as i64)),
                            ("config".into(), Json::Str(tag.clone())),
                            ("s2_delta".into(), Json::Num(*s2d)),
                            ("tkk_mean".into(), Json::Num(*tbar)),
                            ("q_delta".into(), Json::Num(*qd)),
                            ("s2_raw".into(), Json::Num(*s2raw)),
                            ("svac2".into(), Json::Num(*svac2)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(stand_ok && move_ok)),
    ]);
    let p = write_artifact("results/v195_qnec3d.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(a′)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
