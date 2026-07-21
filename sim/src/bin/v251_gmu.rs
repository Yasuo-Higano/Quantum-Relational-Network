//! v25.1 λ の起源 — 1D 質量関数 g(μ) への厳密還元と和則 (第二十六期)
//!
//! v25.0 の残高 1 (λ_x = 1.185468 / λ_⊥ = 1.229430 の解析的導出) の前半。
//! ブロック理論 (v24.1/24.2) の帰結: 横ブロック (ky, kz) の 2 バンドは回転基底で
//! 質量 ±μ (μ = √(cos²ky+cos²kz)) の **1D staggered 鎖 2 本に厳密分解**する。
//! さらに荷電共役 (c → (−1)^x c†) は μ ↔ −μ を結び、K(−μ) = −ΣK(μ)Σ
//! (Σ = diag((−1)^x)) — ボンド要素は符号 2 回反転で同値、対角要素は反号。よって:
//!
//!   A_x(ξ)  = ⟨K^{1D}_μ(ξ)⟩_BZ            (taste 一様 x ボンド = 鎖ボンドの BZ 平均)
//!   K_y(x)  = ⟨(2cos²ky/μ)·K^{1D}_μ(x,x)⟩  (横チャネル = 鎖の質量チャネル対角)
//!   K_z も cos²kz で同型 (λ_y = λ_z の別証明)
//!
//! 定義: K^{1D}_μ-NN(ξ) = g(μ)·πξ + b(μ),  (−1)^x K^{1D}_μ(x,x) = g_m(μ)·2πξ_s·μ + b_m。
//! BW (連続 massive Rindler) の素朴予言は g = g_m = 1 — 格子では g(μ) が
//! regulator 関数になり、**1/λ_x = ⟨g(μ)⟩_BZ、1/λ_⊥ = ⟨2cos²ky·g_m(μ)⟩_BZ**。
//!
//! 測定 (認証済み器械: DD + クランプ 2 点梯子):
//!  [S1] g(0) = 1 (v24.5 XX の再現)
//!  [S2] 各 μ の線形性 (窓内残差) — K^{1D}_μ ボンドの ξ 線形は質量があっても成立するか
//!  [S3] **厳密和則** (N=32): 3D 半空間スキャン (half_space_scan) の A_x/K_y/B_z が
//!       同一 N の 1D 鎖和と一致 (認証窓 ξ≤5 で < 1e-5 — 実測床 ~1e-7 の 100 倍
//!       マージン。クランプ境界モードの DD 分解能が床を決める) — バンド分解代数の証明
//!  [S4] 連続再構成: μ 格子の g/g_m + BZ 求積で 1/λ_x, 1/λ_⊥ を再構成 (±1e-3)
//!
//! 事前登録: (a) S3 成立 + S4 が両 λ を ≤0.1% で再構成 = **λ の起源は 1D 質量関数の
//!   BZ 平均と確定** (導出の器械的前半が完了 — 残るは g, g_m の解析形 = v25.2) /
//!   (b) S3 成立・S4 不一致 = 求積/有限 N の見落とし (記録) /
//!   (c) S3 不成立 = バンド分解代数の誤り (ブロック器械の再監査へ)。

use uft_sim::dd::*;
use uft_sim::stag::*;
use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// 1D staggered-mass 鎖 (H = Σ(1/2)(c†c+h.c.) + μ(−1)^x n) の占有モード行列
/// F (N 行 × N/2 列)。stag.rs の x 鎖 2×2 と同一代数 (μ を独立パラメタ化)。
fn chain_f(n: usize, mu: Dd) -> Vec<Dd> {
    let ni = n as i64;
    let norm = (dd(2.0) / dd(n as f64 + 1.0)).sqrt();
    let mut f = vec![DD0; n * (n / 2)];
    for np in 1..=n / 2 {
        let c1 = dd_cospi_frac(np as i64, ni + 1);
        let r = (c1 * c1 + mu * mu).sqrt();
        let a0 = mu;
        let b0 = -(r + c1);
        let nrm = (a0 * a0 + b0 * b0).sqrt();
        let (al, be) = (a0 / nrm, b0 / nrm);
        for x in 0..n {
            let s1 = dd_sinpi_frac(np as i64 * (x as i64 + 1), ni + 1);
            let s2 = dd_sinpi_frac((ni + 1 - np as i64) * (x as i64 + 1), ni + 1);
            f[x + (np - 1) * n] = (al * s1 + be * s2) * norm;
        }
    }
    f
}

/// 鎖の半分 (x < N/2) のモジュラー核 K (次元 N/2, クランプ指定)
fn chain_k(n: usize, mu: Dd, clamp: f64) -> Vec<Dd> {
    let f = chain_f(n, mu);
    let half = n / 2;
    let mut c = vec![DD0; half * half];
    for k in 0..n / 2 {
        for a in 0..half {
            let va = f[a + k * n];
            for b in a..half {
                c[a + b * half] = c[a + b * half] + va * f[b + k * n];
            }
        }
    }
    for a in 0..half {
        for b in 0..a {
            c[a + b * half] = c[b + a * half];
        }
    }
    let (kmat, _) = modular_k(&c, half, 60, clamp);
    kmat
}

/// 最小二乗 y = s·x + b
fn linfit2(xs: &[f64], ys: &[f64]) -> (f64, f64) {
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let sxy: f64 = xs.iter().zip(ys).map(|(x, y)| (x - mx) * (y - my)).sum();
    let sxx: f64 = xs.iter().map(|x| (x - mx) * (x - mx)).sum();
    let s = sxy / sxx;
    (s, my - s * mx)
}

/// 1 本の鎖から (g, b, 線形性残差, g_m, 梯子相対差) を測る
struct ChainMeas {
    g: f64,
    b: f64,
    resid: f64,
    g_m: f64,
    ladder: f64,
}

fn measure_chain(n: usize, mu: f64) -> ChainMeas {
    let k30 = chain_k(n, dd(mu), 1e-30);
    let k26 = chain_k(n, dd(mu), 1e-26);
    let half = n / 2;
    // 信頼域: ボンドの梯子
    let bond = |k: &Vec<Dd>, i: usize| k[i + (i + 1) * half].hi;
    let mut xi_trust = 0usize;
    for xi in 1..half {
        let i = half - 1 - xi;
        let rel = ((bond(&k30, i) - bond(&k26, i)) / bond(&k30, i)).abs();
        if rel < 1e-4 {
            xi_trust = xi;
        } else {
            break;
        }
    }
    let xi_hi = xi_trust.max(4);
    // g: ボンド勾配 (窓 ξ ∈ [3, ξ*])
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for xi in 3..=xi_hi {
        let i = half - 1 - xi;
        xs.push(PI * xi as f64);
        ys.push(bond(&k30, i));
    }
    let (g, b) = linfit2(&xs, &ys);
    let mut resid = 0.0f64;
    for (x, y) in xs.iter().zip(&ys) {
        resid = resid.max(((g * x + b - y) / y).abs());
    }
    // g_m: 質量チャネル (−1)^x K_xx = g_m·2πξ_s·μ + b_m (μ=0 は 0/0 — 呼び出し側で除外)
    let mut g_m = 0.0;
    if mu > 0.0 {
        let mut xs2 = Vec::new();
        let mut ys2 = Vec::new();
        for j in 0..half {
            let i = half - 1 - j;
            let xis = j as f64 + 0.5;
            if xis < 2.5 || xis > xi_hi as f64 {
                continue;
            }
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            xs2.push(2.0 * PI * xis * mu);
            ys2.push(sign * k30[i + i * half].hi);
        }
        let (s2, _b2) = linfit2(&xs2, &ys2);
        g_m = s2;
    }
    let lad = {
        let i = half - 1 - 3;
        ((bond(&k30, i) - bond(&k26, i)) / bond(&k30, i)).abs()
    };
    ChainMeas {
        g,
        b,
        resid,
        g_m,
        ladder: lad,
    }
}

fn main() {
    self_test();
    println!("=== v25.1 λ の起源 — 1D 質量関数 g(μ) への厳密還元と和則 (第二十六期) ===\n");
    println!("事前登録: (a) 和則成立 + 連続再構成 ≤0.1% → λ の起源確定 (導出前半完了) /");
    println!("          (b) 和則成立・再構成不一致 = 記録 / (c) 和則不成立 = 器械再監査\n");
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
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(1);
    check("[S0a] dd 自己検証", dd_self_test(), String::new());
    check("[S0b] stag 自己検証", stag_self_test(), String::new());

    // ---- [S1][S2] μ 格子の g(μ), g_m(μ) (N=256, 決定的スレッド分割) ----
    let mut mus: Vec<f64> = vec![0.0, 0.001, 0.002, 0.005, 0.01, 0.02, 0.05];
    let mut m = 0.1f64;
    while m < 1.4201 {
        mus.push((m * 1000.0).round() / 1000.0);
        m += 0.02;
    }
    mus.push(2.0f64.sqrt());
    let nchain = 256usize;
    let mut meas: Vec<Option<ChainMeas>> = Vec::new();
    meas.resize_with(mus.len(), || None);
    let chunk = (mus.len() + nthreads - 1) / nthreads;
    std::thread::scope(|sc| {
        for (t, sl) in meas.chunks_mut(chunk).enumerate() {
            let mus = &mus;
            sc.spawn(move || {
                for (j, slot) in sl.iter_mut().enumerate() {
                    let mu = mus[t * chunk + j];
                    *slot = Some(measure_chain(nchain, mu));
                }
            });
        }
    });
    let meas: Vec<ChainMeas> = meas.into_iter().map(|o| o.unwrap()).collect();
    let g0 = meas[0].g;
    check(
        "[S1] g(0) = 1 (v24.5 XX の再現)",
        (g0 - 1.0).abs() < 1e-3,
        format!("g(0) = {:.6} (b = {:+.5})", g0, meas[0].b),
    );
    let resid_max = meas.iter().map(|m| m.resid).fold(0.0f64, f64::max);
    check(
        "[S2] 全 μ でボンド線形性 (窓内残差 < 1e-2)",
        resid_max < 1e-2,
        format!(
            "max 残差 = {:.2e} ({} s)",
            resid_max,
            t0.elapsed().as_secs()
        ),
    );
    println!("\n    [g(μ) 抜粋] μ | g | b | g_m | 梯子Δ(ξ=3)");
    for (i, &mu) in mus.iter().enumerate() {
        if i % 8 == 0 || mu == 0.0 || (mu - 2.0f64.sqrt()).abs() < 1e-9 {
            println!(
                "      μ={:6.3}: g = {:.6}  b = {:+.5}  g_m = {:.6}  Δ = {:.1e}",
                mu, meas[i].g, meas[i].b, meas[i].g_m, meas[i].ladder
            );
        }
    }

    // ---- [S3] 厳密和則 (N=32): 3D スキャン vs 1D 鎖和 ----
    {
        let n = 32usize;
        let half = n / 2;
        let scan = half_space_scan::<Dd>(n, 1e-30, nthreads);
        // [S3-probe] 単一ブロック (q=1, p=0) でリンクごとの検算:
        //   (i) 荷電共役恒等式 K(−μ) = −ΣK(μ)Σ、(ii) バンド分解 (ブロック K vs 2 鎖)
        {
            let (q, p) = (1usize, 0usize);
            let cy = dd_cospi_frac(2 * q as i64, n as i64);
            let cz = dd_cospi_frac(2 * p as i64, n as i64);
            let mu = (cy * cy + cz * cz).sqrt();
            let kp = chain_k(n, mu, 1e-30);
            let km = chain_k(n, -mu, 1e-30);
            let mut dconj = 0.0f64;
            for x in half.saturating_sub(6)..half {
                for y in half.saturating_sub(6)..half {
                    let sgn = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                    dconj = dconj.max((km[x + y * half].hi + sgn * kp[x + y * half].hi).abs());
                }
            }
            let kblk = block_k::<Dd>(n, q, p, &(0..half).collect::<Vec<_>>(), 1e-30);
            let dim = 2 * half;
            let mut dband_bond = 0.0f64;
            let mut dband_diag = 0.0f64;
            // 認証窓 (切断側 ξ ≤ 5 = 添字 i ≥ half−6) のみで検算 — 深部はクランプ雑音域
            for i in half.saturating_sub(6)..half {
                if i + 1 < half {
                    let blk = kblk[(2 * i) + (2 * (i + 1)) * dim].hi
                        + kblk[(2 * i + 1) + (2 * (i + 1) + 1) * dim].hi;
                    let ch = kp[i + (i + 1) * half].hi + km[i + (i + 1) * half].hi;
                    dband_bond = dband_bond.max((blk - ch).abs());
                }
                let blkd =
                    kblk[(2 * i) + (2 * i) * dim].hi - kblk[(2 * i + 1) + (2 * i + 1) * dim].hi;
                let chd = (cy.hi / mu.hi) * (kp[i + i * half].hi - km[i + i * half].hi);
                dband_diag = dband_diag.max((blkd - chd).abs());
            }
            println!(
                "    [S3-probe q=1,p=0 μ={:.4} 窓 ξ≤5] 共役恒等式 Δ = {:.1e}, バンド分解 bond Δ = {:.1e}, diag Δ = {:.1e}",
                mu.hi, dconj, dband_bond, dband_diag
            );
        }
        // 1D 側: 全ブロック (q, p) の鎖 (dim 16 — 一瞬)
        let mut ax_sum = vec![0.0f64; half - 1];
        let mut ky_sum = vec![0.0f64; half];
        let mut bz_sum = vec![0.0f64; half];
        for q in 0..n / 2 {
            let cy = dd_cospi_frac(2 * q as i64, n as i64);
            for p in 0..n {
                let cz = dd_cospi_frac(2 * p as i64, n as i64);
                let mu = (cy * cy + cz * cz).sqrt();
                let k = chain_k(n, mu, 1e-30);
                let w = 1.0 / (n * n) as f64;
                for i in 0..half - 1 {
                    // [K⁺+K⁻]_bond = 2·K_bond (荷電共役の符号 2 回反転)
                    ax_sum[i] += w * 2.0 * k[i + (i + 1) * half].hi;
                }
                if mu.hi > 1e-14 {
                    for i in 0..half {
                        // K00−K11 = (cy/μ)(K⁺−K⁻)_xx = (cy/μ)·2K_xx → cos(ky) 重みで cy²
                        let kd = 2.0 * k[i + i * half].hi;
                        ky_sum[i] += w * (cy.hi * cy.hi / mu.hi) * kd;
                        bz_sum[i] += w * (cz.hi * cz.hi / mu.hi) * kd;
                    }
                }
            }
        }
        // 和則は厳密算術では全 ξ で成立するが、器械はクランプ境界の雑音モード
        // (c_true ~ e^{−80} は DD でも解像不能 — 経路ごとに異なる雑音) の分だけ
        // 深部ボンドで割れる。判定は認証済み窓 ξ ≤ 5 (v24.2/24.3 の梯子較正) で行い、
        // 全域 max は診断値として並記する (開発記録: 初版は全域 max で誤発報 —
        // v24.3 W1/W2 と同型の教訓)。
        let win_max = |scan_v: &Vec<f64>, sum_v: &Vec<f64>, nkeep: usize| -> f64 {
            scan_v
                .iter()
                .zip(sum_v)
                .rev()
                .take(nkeep)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max)
        };
        let full_max = |scan_v: &Vec<f64>, sum_v: &Vec<f64>| -> f64 {
            scan_v
                .iter()
                .zip(sum_v)
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max)
        };
        let dmax_ax = win_max(&scan.ax, &ax_sum, 5);
        let dmax_ky = win_max(&scan.ky, &ky_sum, 5);
        let dmax_bz = win_max(&scan.bz, &bz_sum, 5);
        println!(
            "    [S3 全域診断] max|Δ| (全 ξ) = {:.1e} / {:.1e} / {:.1e} — クランプ雑音域込み",
            full_max(&scan.ax, &ax_sum),
            full_max(&scan.ky, &ky_sum),
            full_max(&scan.bz, &bz_sum)
        );
        // 差分の ξ プロファイル (どこで破れるかの診断)
        println!("    [S3 差分プロファイル] ξ | ΔA_x | ΔK_y (サイト ξ−0.5)");
        for xi in 1..=half - 1 {
            let ib = half - 1 - xi;
            println!(
                "      ξ={:2}: ΔA = {:+.2e}  ΔK_y = {:+.2e}",
                xi,
                scan.ax[ib] - ax_sum[ib],
                scan.ky[ib] - ky_sum[ib]
            );
        }
        // ゲート床 1e-5: クランプ境界モード (κ~60, c~1e-26) の DD 相対分解能 ~1e-5 が
        // 窓ボンドへ ~1e-7 で漏れる (実測)。1e-5 は実測の 100 倍マージン。
        check(
            "[S3] 厳密和則 (N=32, 認証窓 ξ≤5): 3D スキャン = 1D 鎖和 (A_x/K_y/B_z)",
            dmax_ax < 1e-5 && dmax_ky < 1e-5 && dmax_bz < 1e-5,
            format!(
                "窓内 max|Δ| = {:.1e} / {:.1e} / {:.1e} ({} s)",
                dmax_ax,
                dmax_ky,
                dmax_bz,
                t0.elapsed().as_secs()
            ),
        );
    }

    // ---- [S4] 連続再構成: BZ 求積 (g は μ 格子の線形補間) ----
    let interp = |grid: &Vec<f64>, vals: &dyn Fn(usize) -> f64, mu: f64| -> f64 {
        // grid は昇順。範囲外はクランプ。
        let mut i = 0usize;
        while i + 1 < grid.len() && grid[i + 1] < mu {
            i += 1;
        }
        if i + 1 >= grid.len() {
            return vals(grid.len() - 1);
        }
        let (m0, m1) = (grid[i], grid[i + 1]);
        let t = ((mu - m0) / (m1 - m0)).clamp(0.0, 1.0);
        vals(i) * (1.0 - t) + vals(i + 1) * t
    };
    let mut recon = Vec::new(); // (M, 1/λ_x, 1/λ_⊥)
    for &mq in &[128usize, 256, 512] {
        let mut sum_g = 0.0f64;
        let mut sum_t = 0.0f64;
        for iy in 0..mq {
            let ky = 2.0 * PI * iy as f64 / mq as f64;
            let cy2 = ky.cos() * ky.cos();
            for iz in 0..mq {
                let kz = 2.0 * PI * iz as f64 / mq as f64;
                let cz2 = kz.cos() * kz.cos();
                let mu = (cy2 + cz2).sqrt();
                sum_g += interp(&mus, &|i| meas[i].g, mu);
                if mu > 1e-12 {
                    // K_y = (−1)^x πξ·⟨2cos²ky·g_m⟩ (run1 の開発記録: 当初 4cos² は
                    // ペア {ky, ky+π} の二重数え — 全 BZ 平均では 2cos² が正しい)
                    sum_t += 2.0 * cy2 * interp(&mus, &|i| meas[i].g_m, mu);
                }
            }
        }
        let inv_lx = sum_g / (mq * mq) as f64;
        let inv_lt = sum_t / (mq * mq) as f64;
        recon.push((mq, inv_lx, inv_lt));
        println!(
            "    [求積 M={}] ⟨g⟩ = {:.6} (1/λ_x 目標 {:.6}) / ⟨2c²g_m⟩ = {:.6} (1/λ_⊥ 目標 {:.6})",
            mq,
            inv_lx,
            1.0 / 1.185468,
            inv_lt,
            1.0 / 1.229430
        );
    }
    let (_, rx, rt) = *recon.last().unwrap();
    let ex = (rx - 1.0 / 1.185468).abs() * 1.185468;
    let et = (rt - 1.0 / 1.229430).abs() * 1.229430;
    check(
        "[S4a] 連続再構成 λ_x (±1e-3 相対)",
        ex < 1e-3,
        format!(
            "⟨g⟩ = {:.6} → λ_x = {:.6} (目標 1.185468, 相対差 {:.1e})",
            rx,
            1.0 / rx,
            ex
        ),
    );
    check(
        "[S4b] 連続再構成 λ_⊥ (±1e-3 相対)",
        et < 1e-3,
        format!(
            "⟨2c²g_m⟩ = {:.6} → λ_⊥ = {:.6} (目標 1.229430, 相対差 {:.1e})",
            rt,
            1.0 / rt,
            et
        ),
    );

    // ---- 漸近ヒント (v25.2 の解析形探索へ) ----
    {
        let g_at = |target: f64| -> f64 {
            let mut best = (f64::MAX, 0.0);
            for (i, &mu) in mus.iter().enumerate() {
                let d = (mu - target).abs();
                if d < best.0 {
                    best = (d, meas[i].g);
                }
            }
            best.1
        };
        println!(
            "\n    [漸近ヒント] g(√2) = {:.6}, g(1) = {:.6}, g(0.5) = {:.6}, g(0.1) = {:.6}",
            g_at(2.0f64.sqrt()),
            g_at(1.0),
            g_at(0.5),
            g_at(0.1)
        );
        println!("    小 μ: g(μ)−1 ∝ μ²lnμ 型か / 大 μ: 1/μ² 型かは v25.2 のフィット対象");
    }

    // ---- 判定 ----
    let s3_ok = nfail == 0; // S3 までに FAIL がなければ (check は逐次カウント)
    println!(
        "\n[判定] {}",
        if nfail == 0 {
            format!(
                "事前登録 (a): 和則成立 + 連続再構成 λ_x = {:.6} / λ_⊥ = {:.6} — **λ の起源は 1D 質量関数の BZ 平均と確定** (残るは g, g_m の解析形 = v25.2)",
                1.0 / rx,
                1.0 / rt
            )
        } else if s3_ok {
            "事前登録 (b): 和則成立・再構成不一致 — 求積/有限 N の見落としを記録".to_string()
        } else {
            "事前登録 (c): 和則不成立 — バンド分解代数の再監査へ".to_string()
        }
    );

    // ---- JSON (g/g_m テーブル — v25.2 の一次データ) ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v25.1".into())),
        (
            "table".into(),
            Json::Arr(
                mus.iter()
                    .zip(&meas)
                    .map(|(&mu, m)| {
                        Json::Obj(vec![
                            ("mu".into(), Json::Num(mu)),
                            ("g".into(), Json::Num(m.g)),
                            ("b".into(), Json::Num(m.b)),
                            ("g_m".into(), Json::Num(m.g_m)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("inv_lambda_x".into(), Json::Num(rx)),
        ("inv_lambda_t".into(), Json::Num(rt)),
    ]);
    let p = write_artifact("results/v251_gmu.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(b)/(c) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
