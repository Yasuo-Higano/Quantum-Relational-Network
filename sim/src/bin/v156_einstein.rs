//! v15.6 非線形 Einstein への梯子 — 一次の等式・二次の正準エネルギー・全次数の BTZ
//!
//! 残高 3 (v15.0): 「QNEC → Einstein 全次数」。次の方針 (PROMPT/2 柱 5) は、これを
//! 引用で語るのではなく**次数ごとの residual を測る計算問題**に変換することを求める。
//! 1+1 次元では Einstein テンソルが恒等的に零 (2D 重力は自明) なので、非線形検定は
//! AdS₃/CFT₂ 辞書で行う: CFT₂ の区間エントロピー ⇔ AdS₃ の測地線 (RT 公式)。
//! 熱状態の幾何側は **BTZ ブラックホール = 非線形 Einstein 方程式の厳密解**であり、
//! その測地線公式 (Calabrese–Cardy の sinh 型) が格子データに合うかは
//! 「創発幾何が線形化を超えて Einstein 方程式の解になっているか」の検定になる。
//!
//! 梯子 (次数別 residual — 全て [PASS]/[FAIL] 内蔵):
//!   [一次]   δS = δ⟨K⟩ (エンタングルメント第一法則 = 線形化 Einstein; v0.7 の再検証を
//!            無限小ボンド摂動と Richardson 外挿で厳密化)
//!   [二次]   S_rel(ε) = Δ⟨K⟩ − ΔS = ½F ε² + c₃ ε³ + …
//!            F (Fisher 情報 = 正準エネルギー) > 0、冪 2.00、三次係数 c₃ を測定 —
//!            「三次に何が残るか」を言葉でなく数で残す
//!   [全次数] 熱状態 (温度 T = BTZ 質量に対応) の S(ℓ, β) が、自由パラメータ **0** で
//!            BTZ 測地線公式に従う (c=1, v_F=2 は模型から固定; UV 定数は真空との差で
//!            相殺)。線形化 (x² 打ち切り) との残差比を測り、非線形項の必要性を定量化。
//!            小 x 展開の x⁴ 係数 (−c/540) は「計量の二次応答」の係数 — 一致検査。
//!
//! 方法: 円環自由フェルミオン (N=402, 半充填)。熱相関 C_th = f(H₁粒子) は厳密。
//! モジュラー K は区間相関行列の一体スペクトルから厳密に構成 (κ = ln((1−ν)/ν);
//! ℓ ≲ 30 に留めて f64 の限界 [CLAUDE.md] を避ける)。
//! 既知の教訓: 第一法則の検証は無限小のハミルトニアン摂動で行う (有限 1 量子励起は
//! 不等式になる — v0.7 の記録)。

use std::f64::consts::PI;
use uft_sim::*;

const N: usize = 402;

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// 単一粒子ハミルトニアン (環, ホッピング -1, ボンド b の強さを 1+ε に) の基底相関行列
fn gs_correlation(eps: f64, bond: usize) -> Vec<f64> {
    let mut h = vec![0.0; N * N];
    for i in 0..N {
        let j = (i + 1) % N;
        let t = if i == bond { -(1.0 + eps) } else { -1.0 };
        h[i + j * N] += t;
        h[j + i * N] += t;
    }
    let (w, v) = jacobi_eigh(&h, N);
    // 半充填: 下から N/2 個を占有 (N=402 ≡ 2 mod 4 で閉殻 — 縮退の綱渡りなし)
    let nocc = N / 2;
    let _ = w;
    let mut c = vec![0.0; N * N];
    for k in 0..nocc {
        for i in 0..N {
            let vik = v[i + k * N];
            for j in 0..N {
                c[i + j * N] += vik * v[j + k * N];
            }
        }
    }
    c
}

/// 温度 1/β の熱相関行列 (一様環; 波数空間で厳密)
fn thermal_corr_g(beta: f64) -> Vec<f64> {
    // g(d) = (1/N) Σ_k f(ε_k) cos(k d),  ε_k = -2 cos k, f = 1/(1+e^{βε})
    let mut g = vec![0.0; N];
    for (n, gd) in g.iter_mut().enumerate() {
        let mut s = 0.0;
        for k in 0..N {
            let kk = 2.0 * PI * k as f64 / N as f64;
            let e = -2.0 * kk.cos();
            let f = 1.0 / (1.0 + (beta * e).exp());
            s += f * (kk * n as f64).cos();
        }
        *gd = s / N as f64;
    }
    g
}

/// 区間 [0, l) のエントロピー (実対称相関行列の部分行列)
fn interval_entropy_from_g(g: &[f64], l: usize) -> f64 {
    let mut ca = vec![0.0; l * l];
    for i in 0..l {
        for j in 0..l {
            let d = (i as isize - j as isize).unsigned_abs();
            let d = d.min(N - d);
            ca[i + j * l] = g[d];
        }
    }
    entropy_corr_real(&ca, l)
}

/// モジュラーハミルトニアン ⟨K₀⟩_ρ と S(ρ) (K₀ は基準状態 c0 の区間モジュラー演算子)。
/// 摂動状態の相関は複素エルミート (cp_re, cp_im) を許す (波束回転など)。
fn modular_energy_c(c0: &[f64], cp_re: &[f64], cp_im: &[f64], a: usize, l: usize) -> (f64, f64) {
    let mut c0a = vec![0.0; l * l];
    let mut cre = vec![0.0; l * l];
    let mut cim = vec![0.0; l * l];
    for i in 0..l {
        for j in 0..l {
            let (gi, gj) = ((a + i) % N, (a + j) % N);
            c0a[i + j * l] = c0[gi + gj * N];
            cre[i + j * l] = cp_re[gi + gj * N];
            cim[i + j * l] = cp_im[gi + gj * N];
        }
    }
    let (nu, vv) = jacobi_eigh(&c0a, l);
    let mut k_exp = 0.0;
    for k in 0..l {
        let v_k = nu[k].clamp(1e-14, 1.0 - 1e-14);
        let kappa = ((1.0 - v_k) / v_k).ln();
        // ñ_k = v_k† C' v_k (v 実 → 実部のみ寄与)
        let mut nk = 0.0;
        for i in 0..l {
            for j in 0..l {
                nk += vv[i + k * l] * cre[i + j * l] * vv[j + k * l];
            }
        }
        k_exp += kappa * nk - (1.0 - v_k).ln();
    }
    // 相関行列のエントロピーは Σ h2(ν) — 複素エルミートは実埋め込みで (lib)
    let s_p = entropy_corr_herm(&cre, &cim, l);
    (k_exp, s_p)
}

fn modular_energy(c0: &[f64], cp: &[f64], a: usize, l: usize) -> (f64, f64) {
    let zeros = vec![0.0; N * N];
    modular_energy_c(c0, cp, &zeros, a, l)
}

fn main() {
    self_test();
    println!(
        "=== v15.6 非線形 Einstein への梯子: 一次の等式・二次の正準エネルギー・全次数の BTZ ===\n"
    );
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ================= [一次] エンタングルメント第一法則 =================
    // 摂動は v0.7 と同じ粒子-正孔波束の無限小回転 α (状態の一次変分が健全な方向)。
    // 発見の記録: ボンド強度の摂動は半充填の粒子正孔対称性で一次応答がほぼ消える
    // 縮退方向だった (r が発散気味) — その性質は逆に二次 (Fisher) の精密測定に使える
    // ため、[二次] はボンド族で行う (方向の選択が次数を選ぶ)。
    println!("[一次] δS = δ⟨K⟩ (線形化 Einstein の情報側) — 波束回転 α + Richardson (厳密 K)");
    let l_int = 30usize;
    let a0 = 186usize; // 区間 [186, 216) — v0.7 [2a] と同じ数値安全域 (ℓ=30)
    let bond = a0 + 10; // [二次] 用: 区間内部のボンド
    let c0 = gs_correlation(0.0, 0); // 一様環 (ボンド摂動なし) の基底状態
    {
        let two_pi = 2.0 * PI;
        let lf = N as f64;
        let xc = a0 as f64 + 15.0; // 波束中心 = 区間中心
                                   // 正孔側 j∈[nocc/2−21, nocc/2−1] / 粒子側 j∈[nocc/2, nocc/2+21] 型の
                                   // k_F 近傍ガウス波束 (v0.7 の構成を N=402, nocc=201 で使用)
        let make_c = |alpha: f64| -> (Vec<f64>, Vec<f64>) {
            let (jh, jp, sig) = (92.0, 110.0, 5.0);
            let mut hre = vec![0.0; N];
            let mut him = vec![0.0; N];
            let mut pre = vec![0.0; N];
            let mut pim = vec![0.0; N];
            let (mut nh, mut np) = (0.0, 0.0);
            for j in 80..=100 {
                let wj = (-((j as f64 - jh) * (j as f64 - jh)) / (2.0 * sig * sig)).exp();
                nh += wj * wj;
                for x in 0..N {
                    let ph = two_pi * j as f64 * (x as f64 - xc) / lf;
                    hre[x] += wj * ph.cos();
                    him[x] += wj * ph.sin();
                }
            }
            for j in 101..=122 {
                let wj = (-((j as f64 - jp) * (j as f64 - jp)) / (2.0 * sig * sig)).exp();
                np += wj * wj;
                for x in 0..N {
                    let ph = two_pi * j as f64 * (x as f64 - xc) / lf;
                    pre[x] += wj * ph.cos();
                    pim[x] += wj * ph.sin();
                }
            }
            let (nh, np) = ((nh * lf).sqrt(), (np * lf).sqrt());
            for x in 0..N {
                hre[x] /= nh;
                him[x] /= nh;
                pre[x] /= np;
                pim[x] /= np;
            }
            let (s, c) = (alpha.sin(), alpha.cos());
            let mut dre = c0.clone();
            let mut dim = vec![0.0; N * N];
            for x in 0..N {
                for y in 0..N {
                    let hh_re = hre[x] * hre[y] + him[x] * him[y];
                    let hh_im = him[x] * hre[y] - hre[x] * him[y];
                    let pp_re = pre[x] * pre[y] + pim[x] * pim[y];
                    let pp_im = pim[x] * pre[y] - pre[x] * pim[y];
                    let hp_re = hre[x] * pre[y] + him[x] * pim[y];
                    let hp_im = him[x] * pre[y] - hre[x] * pim[y];
                    let ph_re = pre[x] * hre[y] + pim[x] * him[y];
                    let ph_im = pim[x] * hre[y] - pre[x] * him[y];
                    dre[x + y * N] += -s * s * hh_re + s * s * pp_re + s * c * (hp_re + ph_re);
                    dim[x + y * N] += -s * s * hh_im + s * s * pp_im + s * c * (hp_im + ph_im);
                }
            }
            (dre, dim)
        };
        // 射影子検査: 回転後の C' は純粋状態の射影子 (C'² = C') のまま
        {
            let (cre, cim) = make_c(0.05);
            let mut dmax: f64 = 0.0;
            for i in 0..N {
                for j in 0..N {
                    let (mut sr, mut si) = (0.0, 0.0);
                    for k in 0..N {
                        let (a1, b1) = (cre[i + k * N], cim[i + k * N]);
                        let (a2, b2) = (cre[k + j * N], cim[k + j * N]);
                        sr += a1 * a2 - b1 * b2;
                        si += a1 * b2 + b1 * a2;
                    }
                    dmax = dmax
                        .max((sr - cre[i + j * N]).abs())
                        .max((si - cim[i + j * N]).abs());
                }
            }
            check(
                "波束回転の健全性: C'² = C' (純粋 Slater 状態のまま)",
                dmax < 1e-10,
                format!("‖C'²−C'‖_max = {:.1e}", dmax),
            );
        }
        // 一次の推定は奇部で行う: [f(α) − f(−α)]/2 は偶数次 (α² の励起エネルギー項と
        // モジュラー核端のクランプ由来の偶バイアス) を厳密に消し、r_odd = 1 + O(α²)。
        let mut ratios = Vec::new();
        for &al in &[0.04f64, 0.02] {
            let (cre_p, cim_p) = make_c(al);
            let (cre_m, cim_m) = make_c(-al);
            let (kp, sp) = modular_energy_c(&c0, &cre_p, &cim_p, a0, l_int);
            let (km, sm) = modular_energy_c(&c0, &cre_m, &cim_m, a0, l_int);
            ratios.push((sp - sm) / (kp - km));
        }
        // Richardson: r(α) = 1 + bα² + O(α⁴) → [4r(α/2) − r(α)]/3 = 1 + O(α⁴)
        let extrap = (4.0 * ratios[1] - ratios[0]) / 3.0;
        check(
            "第一法則 δS/δ⟨K⟩ → 1 (奇部 + Richardson, 厳密 K)",
            (extrap - 1.0).abs() < 5e-3,
            format!(
                "r_odd(0.04) = {:.5}, r_odd(0.02) = {:.5}, 外挿 = {:.5}",
                ratios[0], ratios[1], extrap
            ),
        );
        // 規格の自己検査: ⟨K₀⟩_ρ0 = S₀
        let (k0, s0b) = modular_energy(&c0, &c0, a0, l_int);
        check(
            "モジュラー規格 ⟨K₀⟩_ρ₀ = S₀ (構成の自己検査)",
            (k0 - s0b).abs() < 1e-9,
            format!("|Δ| = {:.1e}", (k0 - s0b).abs()),
        );
    }

    // ================= [二次] 相対エントロピー = 正準エネルギー =================
    println!("\n[二次] S_rel(ε) = ½F ε² + c₃ ε³ — Fisher 情報 (正準エネルギー) と三次係数の測定");
    {
        let eps_grid = [-0.08f64, -0.04, -0.02, 0.02, 0.04, 0.08];
        let mut srels = Vec::new();
        let mut all_nonneg = true;
        for &eps in &eps_grid {
            let cp = gs_correlation(eps, bond);
            let (k1, s1) = modular_energy(&c0, &cp, a0, l_int);
            let srel = k1 - s1;
            if srel < -1e-12 {
                all_nonneg = false;
            }
            srels.push(srel);
        }
        check(
            "S_rel ≥ 0 (全 ε — 相対エントロピーの正値性)",
            all_nonneg,
            format!(
                "S_rel = {:?}",
                srels
                    .iter()
                    .map(|x| (x * 1e8).round() / 1e8)
                    .collect::<Vec<_>>()
            ),
        );
        // 冪の検査: 偶部から ε² 係数、奇部から ε³ 係数
        let f2 = (srels[3] + srels[2]) / (0.02f64.powi(2)); // ≈ F + O(ε²)
        let f2b = (srels[4] + srels[1]) / (0.04f64.powi(2));
        let p_est = ((srels[4] + srels[1]) / (srels[3] + srels[2])).ln() / (2.0f64).ln();
        check(
            "ε 冪 = 2.00 (二次性) と Fisher F > 0",
            (p_est - 2.0).abs() < 0.05 && f2 > 0.0,
            format!(
                "p̂ = {:.4}, F(0.02) = {:.6}, F(0.04) = {:.6}",
                p_est, f2, f2b
            ),
        );
        // 三次係数 (奇部): c₃ = [S_rel(ε) − S_rel(−ε)] / (2ε³)
        let c3_a = (srels[3] - srels[2]) / (2.0 * 0.02f64.powi(3));
        let c3_b = (srels[4] - srels[1]) / (2.0 * 0.04f64.powi(3));
        println!(
            "    三次係数 c₃ ≈ {:.4} (ε=0.02), {:.4} (ε=0.04) — 「三次に残る項」の大きさが数になった",
            c3_a, c3_b
        );
        check(
            "次数の階層 |c₃ ε³| ≪ ½F ε² (ε=0.04 で 1/10 以下)",
            (c3_b * 0.04).abs() < 0.1 * f2b / 2.0,
            format!("比 = {:.3}", (c3_b * 0.04).abs() / (f2b / 2.0)),
        );
    }

    // ================= [全次数] BTZ 測地線公式 (自由パラメータ 0) =================
    println!(
        "\n[全次数] 熱状態の S(ℓ, β) vs BTZ (非線形 Einstein の厳密解) — c=1, v_F=2 固定・定数なし"
    );
    {
        let c_th = 1.0; // Dirac CFT (v0.5/v15.4 で測定済み)
        let v_f = 2.0; // 半充填の Fermi 速度 (分散 -2cos k の傾き)
        let g_vac = thermal_corr_g(1e6); // β→∞ = 真空
        let betas = [8.0f64, 12.0, 16.0];
        let ls: Vec<usize> = vec![10, 20, 30, 40, 50, 60];
        let mut max_resid: f64 = 0.0;
        let mut max_resid_b12: f64 = 0.0;
        let mut max_rel_b8: f64 = 0.0;
        let mut max_x: f64 = 0.0;
        let mut lin_fail_ratio: f64 = f64::INFINITY;
        println!("    β     ℓ    x=πℓ/(vβ)   D_格子      D_BTZ       D_線形化    |BTZ 残差|");
        for &beta in &betas {
            let g_th = thermal_corr_g(beta);
            for &l in &ls {
                let s_th = interval_entropy_from_g(&g_th, l);
                let s_vac = interval_entropy_from_g(&g_vac, l);
                let d_lat = s_th - s_vac;
                let x = PI * l as f64 / (v_f * beta);
                let ring = PI * l as f64 / N as f64;
                // 理論: D = (c/3)[ln(sinh x / x) − ln(sin(ring)/ring)] … 真空側は環の弦長
                let d_btz = (c_th / 3.0) * ((x.sinh() / x).ln() - (ring.sin() / ring).ln());
                let d_lin = (c_th / 3.0) * (x * x / 6.0 - (ring.sin() / ring).ln());
                let resid = (d_lat - d_btz).abs();
                max_resid = max_resid.max(resid);
                if beta >= 12.0 {
                    max_resid_b12 = max_resid_b12.max(resid);
                } else {
                    max_rel_b8 = max_rel_b8.max(resid / d_lat.abs().max(1e-12));
                }
                if x > max_x {
                    max_x = x;
                }
                if x > 2.0 {
                    let r_lin = (d_lat - d_lin).abs();
                    lin_fail_ratio = lin_fail_ratio.min(r_lin / resid.max(1e-12));
                }
                if l % 20 == 10 || l == 60 {
                    println!(
                        "    {:4.0}  {:3}   {:6.3}     {:8.4}   {:8.4}   {:8.4}    {:.1e}",
                        beta, l, x, d_lat, d_btz, d_lin, resid
                    );
                }
            }
        }
        check(
            "BTZ 公式との一致 (β ≥ 12): |残差| < 0.01 nats — 自由パラメータ 0 で x ≈ 7.9 の深い非線形域まで",
            max_resid_b12 < 0.01,
            format!("max|残差| (β≥12) = {:.4}", max_resid_b12),
        );
        check(
            "β = 8 は相対 2% 以内 — 温度が格子分散を照らし始める境界 (創発幾何の UV カットオフの位置)",
            max_rel_b8 < 0.02,
            format!(
                "max|残差|/D (β=8) = {:.4} (絶対 {:.4} — 分散補正 O((πT/v)²) の域)",
                max_rel_b8, max_resid
            ),
        );
        check(
            "線形化 (x² 打ち切り) は x>2 で BTZ の 10 倍以上外れる — 非線形項は必要",
            lin_fail_ratio > 10.0,
            format!("最小比 (線形化残差 / BTZ 残差) = {:.0}", lin_fail_ratio),
        );
        // 小 x 域の係数: D_th − 環補正 = a₂x² + a₄x⁴ + a₆x⁶; a₂ = c/18, a₄ = −c/540。
        // x ≤ 1.6 に制限し x⁶ 項も入れて高次汚染を制御する (初版は x ≤ 3.9 の 2 項 fit で
        // a₄ が 2.5 倍ずれた — 展開の適用域を破った設計ミスとして記録)。
        let beta_c = 40.0;
        let g_c = thermal_corr_g(beta_c);
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for &l in &[10usize, 16, 22, 28, 34, 40] {
            let x = PI * l as f64 / (v_f * beta_c);
            let ring = PI * l as f64 / N as f64;
            let d = interval_entropy_from_g(&g_c, l) - interval_entropy_from_g(&g_vac, l)
                + (c_th / 3.0) * (ring.sin() / ring).ln();
            xs.push(x);
            ys.push(d);
        }
        // 最小二乗 3 項: y = a2 x² + a4 x⁴ + a6 x⁶ (正規方程式 3×3 を jacobi で解く)
        let mut ata = vec![0.0; 9];
        let mut aty = [0.0f64; 3];
        for (x, y) in xs.iter().zip(&ys) {
            let b = [x.powi(2), x.powi(4), x.powi(6)];
            for i in 0..3 {
                for j in 0..3 {
                    ata[i + j * 3] += b[i] * b[j];
                }
                aty[i] += b[i] * y;
            }
        }
        let (w3, v3) = jacobi_eigh(&ata, 3);
        let mut coef = [0.0f64; 3];
        for k in 0..3 {
            let mut proj = 0.0;
            for i in 0..3 {
                proj += v3[i + k * 3] * aty[i];
            }
            for i in 0..3 {
                coef[i] += v3[i + k * 3] * proj / w3[k];
            }
        }
        let (a2, a4, a6) = (coef[0], coef[1], coef[2]);
        let (a2_th, a4_th) = (c_th / 18.0, -c_th / 540.0);
        check(
            "小 x 展開 (x≤1.6, 3 項): a₂ = c/18 (±3%), a₄ = −c/540 (±15%) — 計量の一次・二次応答",
            (a2 / a2_th - 1.0).abs() < 0.03 && (a4 / a4_th - 1.0).abs() < 0.15,
            format!(
                "a₂ = {:.5} (理論 {:.5}), a₄ = {:.6} (理論 {:.6}), a₆ = {:.2e} (理論 {:.2e})",
                a2,
                a2_th,
                a4,
                a4_th,
                a6,
                c_th / 2835.0 / 3.0
            ),
        );
    }

    // ================= artifact =================
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v15.6".into())),
        ("n".into(), Json::Num(N as f64)),
        (
            "note".into(),
            Json::Str(
                "一次=第一法則 / 二次=Fisher / 全次数=BTZ (数値は stdout が一次ソース)".into(),
            ),
        ),
    ]);
    let p = write_artifact("results/v156_einstein.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 一次の等式・二次の正準エネルギー・全次数の BTZ が同じ格子データの上に立った"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
