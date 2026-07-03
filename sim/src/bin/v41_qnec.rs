//! v4.1 量子ヌルエネルギー条件 (QNEC) — 単調性原理の局所版を実時間で検証
//!
//! 古典的なヌルエネルギー条件 T_kk ≥ 0 は量子論で破れる (Casimir 等)。その量子版:
//!   QNEC:  2π ⟨T_kk⟩ ≥ S''   (領域の端点を光的方向 k に変形したときのエントロピーの二階微分)
//! (2d CFT では強い形 S'' + (6/c)(S')² ≤ 2π⟨T_kk⟩ も成立)
//! これは v2.1 の相対エントロピー正値性の「微分版」であり、Wall による一般化第二法則の
//! 証明の要 — 単調性原理 M を時空の各点・各光的方向に局所化した命題である。
//!
//! 検証系: フェルミオン鎖の基底状態 + 右向きカイラル波束 (厳密な実時間発展)。
//! 領域 A = [a(λ), b] の端点を光的軌道 a = a₀ - v_F t (x + v_F t = 一定) で動かし、
//! 波束が端点を横切る瞬間の S(λ) と T₋₋ を比較する。

use uft_sim::*;

fn h2(z: f64) -> f64 {
    let z = z.clamp(1e-14, 1.0 - 1e-14);
    -z * z.ln() - (1.0 - z) * (1.0 - z).ln()
}
fn entropy_herm(cre: &[f64], cim: &[f64], n: usize) -> f64 {
    let m = 2 * n;
    let mut a = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            a[i + j * m] = cre[i + j * n];
            a[i + (j + n) * m] = -cim[i + j * n];
            a[(i + n) + j * m] = cim[i + j * n];
            a[(i + n) + (j + n) * m] = cre[i + j * n];
        }
    }
    let (w, _) = jacobi_eigh(&a, m);
    0.5 * w.iter().map(|&z| h2(z)).sum::<f64>()
}

fn main() {
    let n = 302usize;
    let nocc = 151usize; // 閉殻半充填
    let two_pi = 2.0 * std::f64::consts::PI;
    let vf = 2.0f64;
    println!("=== v4.1 QNEC: 光的変形に沿うエントロピーとエネルギー (N={}, 厳密) ===\n", n);

    // 基底状態相関 (実)
    let c0 = |d: isize| -> f64 {
        let d = d.unsigned_abs();
        if d == 0 {
            return nocc as f64 / n as f64;
        }
        (std::f64::consts::PI * d as f64 / 2.0).sin()
            / (n as f64 * (std::f64::consts::PI * d as f64 / n as f64).sin())
    };
    // 右向きカイラル波束 (占有側 j∈[55,75] 中心67 / 空側 j∈[76,97] 中心84, σ=5, α=0.5)
    let alpha = 0.5f64;
    let xc = 90.0f64;
    let (mut c0re, mut c0im) = (vec![0.0; n * n], vec![0.0; n * n]);
    for x in 0..n {
        for y in 0..n {
            c0re[x + y * n] = c0(x as isize - y as isize);
        }
    }
    {
        let (jh, jp, sig) = (67.0, 84.0, 5.0);
        let mut hre = vec![0.0; n];
        let mut him = vec![0.0; n];
        let mut pre = vec![0.0; n];
        let mut pim = vec![0.0; n];
        let (mut nh, mut np) = (0.0, 0.0);
        for j in 55..=75 {
            let wj = (-((j as f64 - jh) * (j as f64 - jh)) / (2.0 * sig * sig)).exp();
            nh += wj * wj;
            for x in 0..n {
                let ph = two_pi * j as f64 * (x as f64 - xc) / n as f64;
                hre[x] += wj * ph.cos();
                him[x] += wj * ph.sin();
            }
        }
        for j in 76..=97 {
            let wj = (-((j as f64 - jp) * (j as f64 - jp)) / (2.0 * sig * sig)).exp();
            np += wj * wj;
            for x in 0..n {
                let ph = two_pi * j as f64 * (x as f64 - xc) / n as f64;
                pre[x] += wj * ph.cos();
                pim[x] += wj * ph.sin();
            }
        }
        let (nh, np) = ((nh * n as f64).sqrt(), (np * n as f64).sqrt());
        for x in 0..n {
            hre[x] /= nh;
            him[x] /= nh;
            pre[x] /= np;
            pim[x] /= np;
        }
        let (s, c) = (alpha.sin(), alpha.cos());
        for x in 0..n {
            for y in 0..n {
                let hh_re = hre[x] * hre[y] + him[x] * him[y];
                let hh_im = him[x] * hre[y] - hre[x] * him[y];
                let pp_re = pre[x] * pre[y] + pim[x] * pim[y];
                let pp_im = pim[x] * pre[y] - pre[x] * pim[y];
                let hp_re = hre[x] * pre[y] + him[x] * pim[y];
                let hp_im = him[x] * pre[y] - hre[x] * pim[y];
                let ph_re = pre[x] * hre[y] + pim[x] * him[y];
                let ph_im = pim[x] * hre[y] - pre[x] * him[y];
                c0re[x + y * n] += -s * s * hh_re + s * s * pp_re + s * c * (hp_re + ph_re);
                c0im[x + y * n] += -s * s * hh_im + s * s * pp_im + s * c * (hp_im + ph_im);
            }
        }
    }

    // 実時間発展: U(t)_{xy} = u((x-y) mod N), u(d) = (1/N)Σ_j e^{i(k_j d - ε_j t)}
    let evolve = |t: f64| -> (Vec<f64>, Vec<f64>) {
        let mut ur = vec![0.0; n * n];
        let mut ui = vec![0.0; n * n];
        let mut ud = vec![(0.0f64, 0.0f64); n];
        for d in 0..n {
            let (mut a, mut b) = (0.0, 0.0);
            for j in 0..n {
                let k = two_pi * j as f64 / n as f64;
                let e = -2.0 * k.cos();
                let ph = k * d as f64 - e * t;
                a += ph.cos();
                b += ph.sin();
            }
            ud[d] = (a / n as f64, b / n as f64);
        }
        for x in 0..n {
            for y in 0..n {
                let d = (x + n - y) % n;
                ur[x + y * n] = ud[d].0;
                ui[x + y * n] = ud[d].1;
            }
        }
        // C(t) = U C0 U†
        let t1r = matmul(&ur, &c0re, n);
        let t1r2 = matmul(&ui, &c0im, n);
        let t1i = matmul(&ur, &c0im, n);
        let t1i2 = matmul(&ui, &c0re, n);
        let ar: Vec<f64> = t1r.iter().zip(&t1r2).map(|(a, b)| a - b).collect();
        let ai: Vec<f64> = t1i.iter().zip(&t1i2).map(|(a, b)| a + b).collect();
        // B = A U†: (U†)_{xy} = conj(U_{yx})
        let mut vr = vec![0.0; n * n];
        let mut vi = vec![0.0; n * n];
        for x in 0..n {
            for y in 0..n {
                vr[x + y * n] = ur[y + x * n];
                vi[x + y * n] = -ui[y + x * n];
            }
        }
        let b1 = matmul(&ar, &vr, n);
        let b2 = matmul(&ai, &vi, n);
        let b3 = matmul(&ar, &vi, n);
        let b4 = matmul(&ai, &vr, n);
        (
            b1.iter().zip(&b2).map(|(a, b)| a - b).collect(),
            b3.iter().zip(&b4).map(|(a, b)| a + b).collect(),
        )
    };

    // エネルギー密度と流れ (基底状態差し引き)
    let h_gs = -2.0 * c0(1);
    let t_minus = |cre: &Vec<f64>, cim: &Vec<f64>, x: usize| -> f64 {
        // 3 ボンド平均の h, j
        let (mut h, mut j) = (0.0, 0.0);
        for b in [x - 1, x, x + 1] {
            h += -2.0 * cre[b + (b + 1) * n] - h_gs;
            j += 2.0 * cim[(b - 1) + (b + 1) * n];
        }
        h /= 3.0;
        j /= 3.0;
        // 変形の接ベクトル k=(Δτ,Δx)=(1,-1) (c=1 単位, σ=サイト変位) に対する
        // T_kk = T_ττ - 2T_τx + T_xx = 2h/v_F + 2j/v_F²  (= 4·T₋₋)
        2.0 * (h / vf + j / (vf * vf))
    };

    // ---- 掃引: 端点 a_m = 190 - m, 時刻 t_m = m/2 (x + v_F t = 190 の光的線) ----
    let b_end = 250usize;
    let nsteps = 60usize;
    let mut svals = Vec::new();
    let mut tmm = Vec::new();
    let mut chir = (0.0f64, 0.0f64);
    for m in 0..=nsteps {
        let t = m as f64 / 2.0;
        let (cre, cim) = evolve(t);
        let a = 190 - m;
        let l = b_end - a + 1;
        let mut are = vec![0.0; l * l];
        let mut aim = vec![0.0; l * l];
        for i in 0..l {
            for jj in 0..l {
                are[i + jj * l] = cre[(a + i) + (a + jj) * n];
                aim[i + jj * l] = cim[(a + i) + (a + jj) * n];
            }
        }
        svals.push(entropy_herm(&are, &aim, l));
        tmm.push(t_minus(&cre, &cim, a));
        if m == 0 {
            // 較正: カイラル性 (エネルギー流 j ≈ v_F × エネルギー密度 h)
            let (mut sh, mut sj) = (0.0, 0.0);
            for x in 2..n - 2 {
                sh += -2.0 * cre[x + (x + 1) * n] - h_gs;
                sj += 2.0 * cim[(x - 1) + (x + 1) * n];
            }
            chir = (sh, sj);
        }
    }
    println!("[較正] 波束の全エネルギー Σh = {:.4}, 全エネルギー流 Σj = {:.4}", chir.0, chir.1);
    println!("       比 Σj/(v_F·Σh) = {:.3} (右向きカイラルなら +1)\n", chir.1 / (vf * chir.0));

    // ---- QNEC 検定 (同パリティ差分 Δσ=2) ----
    println!("[A] 光的変形 (x+v_Ft=const, 波束が端点を横切る) に沿う QNEC");
    println!("  σ    S(σ)      S''       6(S')²    2πT_kk    基本ギャップ  強ギャップ");
    let mut min_gap_basic = f64::INFINITY;
    let mut min_gap_strong = f64::INFINITY;
    for m in (2..nsteps - 1).step_by(2) {
        let sp = (svals[m + 2] - svals[m - 2]) / 4.0;
        let spp = (svals[m + 2] - 2.0 * svals[m] + svals[m - 2]) / 4.0;
        // S'' のステンシル [m-2, m+2] に合わせ T も同じ窓で加重平均 (1/4,1/2,1/4)
        let t_avg = 0.25 * tmm[m - 2] + 0.5 * tmm[m] + 0.25 * tmm[m + 2];
        let rhs = two_pi * t_avg;
        let gap_b = rhs - spp;
        let gap_s = rhs - spp - 6.0 * sp * sp;
        min_gap_basic = min_gap_basic.min(gap_b);
        min_gap_strong = min_gap_strong.min(gap_s);
        if m % 4 == 0 {
            println!(
                "  {:3}  {:.5}  {:+.5}  {:+.5}  {:+.5}  {:+.5}    {:+.5}",
                m, svals[m], spp, 6.0 * sp * sp, rhs, gap_b, gap_s
            );
        }
    }
    let tol = 2e-4;
    println!("\n  => 基本 QNEC (2πT₋₋ - S'' ≥ 0):   最小ギャップ {:+.5}  {}",
        min_gap_basic, pass(min_gap_basic > -tol));
    println!("     強い QNEC (- 6(S')² も含む):    最小ギャップ {:+.5}  {}",
        min_gap_strong, pass(min_gap_strong > -tol));

    // ---- 対照: 波束と並走する光的変形 (x - v_F t = const) — T₊₊ ≈ 0, S も不変のはず ----
    println!("\n[B] 並走する光的変形 (x-v_Ft=const): カイラル状態は静止して見える");
    let mut sco = Vec::new();
    for m in [0usize, 10, 20, 30] {
        let t = m as f64 / 2.0;
        let (cre, cim) = evolve(t);
        let a = 60 + m;
        let l = b_end - a + 1;
        let mut are = vec![0.0; l * l];
        let mut aim = vec![0.0; l * l];
        for i in 0..l {
            for jj in 0..l {
                are[i + jj * l] = cre[(a + i) + (a + jj) * n];
                aim[i + jj * l] = cim[(a + i) + (a + jj) * n];
            }
        }
        sco.push(entropy_herm(&are, &aim, l));
    }
    // 端点位置が変わると真空エントロピー自体は変わる (区間長が違う)…ではなく b は固定なので
    // 区間長も変わる。カイラル性の検定は S の「波束寄与」が一定であること:
    // 真空値を引いた ΔS を比較する
    let mut svac = Vec::new();
    for m in [0usize, 10, 20, 30] {
        let a = 60 + m;
        let l = b_end - a + 1;
        let mut are = vec![0.0; l * l];
        for i in 0..l {
            for jj in 0..l {
                are[i + jj * l] = c0((a + i) as isize - (a + jj) as isize);
            }
        }
        let (w, _) = jacobi_eigh(&are, l);
        svac.push(w.iter().map(|&z| h2(z)).sum::<f64>());
    }
    print!("  ΔS(波束寄与) = ");
    let mut dsmax: f64 = 0.0;
    let ds0 = sco[0] - svac[0];
    for i in 0..4 {
        let ds = sco[i] - svac[i];
        dsmax = dsmax.max((ds - ds0).abs());
        print!("{:.5}  ", ds);
    }
    println!("\n  => 並走系での ΔS の変動 {:.1e} (カイラル状態は x-v_Ft の関数 = 凍結)  {}",
        dsmax, pass(dsmax < 5e-3));
    println!("\n結論: 量子論では局所エネルギーは負にもなれるが、「情報の増え方の加速度 S''」が");
    println!("      常にエネルギー流でキャップされる (QNEC)。単調性原理 M の局所版が成立し、");
    println!("      これが一般化第二法則 (BH面積+外部エントロピーの単調性) の証明の要になる。");
    println!("      v2.1 (二次) → v4.1 (各点・各方向) — 非線形重力の情報側の部品が揃った。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
