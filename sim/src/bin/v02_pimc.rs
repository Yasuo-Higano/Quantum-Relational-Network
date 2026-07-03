//! v0.2 Wick 回転の橋 — ユークリッド経路積分モンテカルロ (PIMC)
//! ħ = m = 1。虚時間 τ = it に回すと量子力学の経路積分は
//! 1次元統計力学 (弾性鎖の正準分布 e^{-S_E}) と厳密に同型になる:
//!   Z = Tr e^{-βH} = ∫Dx e^{-S_E},  S_E = Σ [ (x_{i+1}-x_i)²/2Δτ + Δτ V(x_i) ]
//!
//! A: 調和振動子の ⟨x²⟩(β) — 量子(β→∞)から古典(β→0)へのクロスオーバー
//! B: 基底エネルギー E0 (ビリアル)
//! C: 二重井戸のトンネル分裂 ΔE — 厳密対角化と比較(インスタントンの実演)

use uft_sim::*;

/// 単サイト Metropolis + 区間反転更新の 1 スイープ
fn sweep(
    x: &mut [f64],
    dtau: f64,
    v: &impl Fn(f64) -> f64,
    step: f64,
    rng: &mut Rng,
    seg_flip: bool,
) -> f64 {
    let n = x.len();
    let mut acc = 0usize;
    for _ in 0..n {
        let i = rng.range(n);
        let (l, r) = (x[(i + n - 1) % n], x[(i + 1) % n]);
        let xo = x[i];
        let xn = xo + step * (2.0 * rng.f64() - 1.0);
        let ds = ((xn - l).powi(2) + (xn - r).powi(2) - (xo - l).powi(2) - (xo - r).powi(2))
            / (2.0 * dtau)
            + dtau * (v(xn) - v(xo));
        if ds <= 0.0 || rng.f64() < (-ds).exp() {
            x[i] = xn;
            acc += 1;
        }
    }
    if seg_flip {
        // 区間 [i0, i0+len) を x → -x 反転 (偶ポテンシャル用): 境界 2 ボンドのみ ΔS
        let i0 = rng.range(n);
        let len = 1 + rng.range(n / 2);
        let jl = (i0 + n - 1) % n;
        let jr = (i0 + len) % n;
        let last = (i0 + len - 1) % n;
        let ds = 2.0 * (x[jl] * x[i0] + x[jr] * x[last]) / dtau;
        if ds <= 0.0 || rng.f64() < (-ds).exp() {
            for k in 0..len {
                let idx = (i0 + k) % n;
                x[idx] = -x[idx];
            }
        }
    }
    acc as f64 / n as f64
}

/// 1D ハミルトニアンの厳密対角化 (グリッド上, Jacobi 法)
fn exact_levels(v: &impl Fn(f64) -> f64, x0: f64, x1: f64, n: usize, nlevels: usize) -> Vec<f64> {
    let dx = (x1 - x0) / (n - 1) as f64;
    let mut h = vec![0.0; n * n];
    for i in 0..n {
        let x = x0 + i as f64 * dx;
        h[i + i * n] = 1.0 / (dx * dx) + v(x);
        if i + 1 < n {
            h[i + (i + 1) * n] = -0.5 / (dx * dx);
            h[(i + 1) + i * n] = -0.5 / (dx * dx);
        }
    }
    let (w, _) = jacobi_eigh(&h, n);
    w.into_iter().take(nlevels).collect()
}

fn main() {
    self_test();
    let mut rng = Rng::new(20260703);
    println!("=== v0.2 Wick 回転: 量子力学 = 虚時間の統計力学 ===\n");

    // ---- A: 調和振動子 V = x²/2 の ⟨x²⟩(β) ----
    println!("[A] 調和振動子 ⟨x²⟩(β): 量子↔古典クロスオーバー");
    println!("    連続理論: ⟨x²⟩ = (1/2)coth(β/2)   古典極限: kT = 1/β");
    println!("  β     MC結果        格子厳密解  連続理論  古典値");
    let dtau: f64 = 0.0625;
    let mut all_pass = true;
    for &beta in &[0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0] {
        let n = (beta / dtau).round() as usize;
        let mut x = vec![0.0f64; n];
        let vho = |y: f64| 0.5 * y * y;
        let step = 2.0 * (dtau.sqrt().min(1.0));
        for _ in 0..20_000 {
            sweep(&mut x, dtau, &vho, step, &mut rng, false);
        }
        let mut xs = Vec::with_capacity(200_000);
        for _ in 0..200_000 {
            sweep(&mut x, dtau, &vho, step, &mut rng, false);
            xs.push(x.iter().map(|a| a * a).sum::<f64>() / n as f64);
        }
        let (m, e) = mean_err(&xs);
        // 格子厳密解: ⟨x²⟩ = (1/N) Σ_q [ (2-2cos q)/Δτ + Δτ ]^{-1}
        let lat: f64 = (0..n)
            .map(|k| {
                let q = 2.0 * std::f64::consts::PI * k as f64 / n as f64;
                1.0 / ((2.0 - 2.0 * q.cos()) / dtau + dtau)
            })
            .sum::<f64>()
            / n as f64;
        let cont = 0.5 / (beta / 2.0).tanh();
        let ok = (m - lat).abs() < 4.0 * e.max(1e-9);
        all_pass &= ok;
        println!(
            "  {:4.1}  {:.4}±{:.4}  {:.4}     {:.4}    {:.4}  {}",
            beta,
            m,
            e,
            lat,
            cont,
            1.0 / beta,
            pass(ok)
        );
    }
    println!("  => 高温(β小)で古典値 kT に、低温で量子零点揺らぎ 1/2 に収束 {}\n", pass(all_pass));

    // ---- B: 基底エネルギー (ビリアル: E = ⟨V⟩ + ⟨xV'⟩/2 = ⟨x²⟩ for HO) ----
    println!("[B] 基底エネルギー E0 (β=32, ビリアル定理)");
    {
        let beta = 32.0;
        let n = (beta / dtau).round() as usize;
        let mut x = vec![0.0f64; n];
        let vho = |y: f64| 0.5 * y * y;
        for _ in 0..20_000 {
            sweep(&mut x, dtau, &vho, 0.5, &mut rng, false);
        }
        let mut es = Vec::new();
        for _ in 0..200_000 {
            sweep(&mut x, dtau, &vho, 0.5, &mut rng, false);
            es.push(x.iter().map(|a| a * a).sum::<f64>() / n as f64);
        }
        let (m, e) = mean_err(&es);
        println!(
            "  E0(MC) = {:.4} ± {:.4}   厳密値 ħω/2 = 0.5  {}",
            m,
            e,
            pass((m - 0.5).abs() < 0.01)
        );
        println!("  (零点エネルギー = 純粋に量子的な「揺らぎの熱力学」として得られる)\n");
    }

    // ---- C: 二重井戸のトンネル分裂 ----
    println!("[C] 二重井戸 V = (x²-η²)², η=1.2 — インスタントンとトンネル分裂");
    let eta = 1.2f64;
    let vdw = move |y: f64| (y * y - eta * eta).powi(2);
    let lv = exact_levels(&vdw, -4.5, 4.5, 361, 4);
    let gap_exact = lv[1] - lv[0];
    println!(
        "  厳密対角化: E0={:.5}, E1={:.5}, E2={:.5} → 分裂 ΔE = {:.5}",
        lv[0], lv[1], lv[2], gap_exact
    );
    println!("  (井戸の振動量子 ω = √(8η²) = {:.3} に比べ ΔE が桁違いに小さい = 障壁下の禁止過程)", (8.0 * eta * eta).sqrt());

    let beta: f64 = 24.0;
    let dtau_dw: f64 = 0.05;
    let n = (beta / dtau_dw).round() as usize;
    let mut x = vec![eta; n];
    let step = 0.55;
    for _ in 0..40_000 {
        sweep(&mut x, dtau_dw, &vdw, step, &mut rng, true);
    }
    let ncorr = n / 2;
    let mut corr = vec![0.0f64; ncorr + 1];
    let mut nm = 0usize;
    let mut kinks = 0usize;
    let nmeas = 300_000;
    for it in 0..nmeas {
        sweep(&mut x, dtau_dw, &vdw, step, &mut rng, true);
        if it % 10 == 0 {
            for d in 0..=ncorr {
                let mut c = 0.0;
                for i in 0..n {
                    c += x[i] * x[(i + d) % n];
                }
                corr[d] += c / n as f64;
            }
            nm += 1;
            kinks += (0..n).filter(|&i| x[i] * x[(i + 1) % n] < 0.0).count();
        }
    }
    for c in corr.iter_mut() {
        *c /= nm as f64;
    }
    // C(τ) ∝ cosh(Δ(β/2-τ)) から有効ギャップ: cosh比 (C(τ-Δτ)+C(τ+Δτ))/(2C(τ))
    let mut gaps = Vec::new();
    let d0 = (0.6 / dtau_dw) as usize;
    let d1 = (2.4 / dtau_dw) as usize;
    for d in d0..=d1 {
        let r = (corr[d - 1] + corr[d + 1]) / (2.0 * corr[d]);
        if r > 1.0 {
            gaps.push((r.acosh()) / dtau_dw);
        }
    }
    let (gap_mc, gap_err) = mean_err(&gaps);
    let kink_density = kinks as f64 / (nm as f64 * n as f64);
    println!("  PIMC: β={}, Δτ={}, 経路上のキンク(井戸間の壁)密度 = {:.4} /サイト", beta, dtau_dw, kink_density);
    println!(
        "  相関関数 C(τ)~cosh(ΔE(β/2-τ)) から: ΔE(MC) = {:.4} ± {:.4}  (厳密 {:.4})  {}",
        gap_mc,
        gap_err,
        gap_exact,
        pass((gap_mc - gap_exact).abs() < 5.0 * gap_err.max(0.01))
    );
    println!("  => 古典的には縮退した二つの真空が、虚時間の「キンク気体」(インスタントン) で混ざり");
    println!("     エネルギー分裂が生じる。トンネル = ユークリッド統計力学の位相欠陥。");
    println!("\n結論: 量子力学と統計力学は Wick 回転で厳密に同型 (時間 ↔ 逆温度)。");
    println!("      これは v0.5 以降の「時間・温度・エンタングルメントの三位一体」の最初の証拠。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
