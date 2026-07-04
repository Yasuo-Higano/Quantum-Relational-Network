//! v3.4 加速 = 温度 — Unruh 効果の検出器応答 (数値検証)
//!
//! v0.7 は「真空を区間に制限すると熱 (モジュラー温度) が現れる」ことを静的に示した。
//! その操作的・動的な顔が Unruh 効果 (1976): ミンコフスキー真空を一様加速 a で走る
//! 検出器は、温度 T = a/2π の熱浴を見る。これは Hawking 輻射 (BH 地平線) と
//! Gibbons–Hawking 温度 (dS 地平線, v2.4 の H/2π) の共通の数学的核心である。
//!
//! 方法: 検出器が結合する場の相関関数 (Wightman) を世界線に沿って評価し、
//! 応答関数 F(ω) = ∫ dΔτ e^{-iωΔτ} W(Δτ) を数値フーリエ変換。
//!   慣性系:  W ∝ -1/4π(Δτ-iε)²        → F(ω>0) = 0 (真空は検出器を励起しない)
//!   加速系:  W ∝ -a²/16π sinh²(a(Δτ-iε)/2) → F(ω) ∝ ω/(e^{2πω/a}-1) (プランク分布!)
//! 検証: 詳細釣り合い F(ω)/F(-ω) = e^{-2πω/a} と有効温度 T = a/2π。

fn main() {
    let a = 1.0f64;
    let eps = 2e-3;
    let big_l = 120.0;
    let dt = 4e-4;
    println!(
        "=== v3.4 Unruh 効果: 加速検出器の応答 (a={}, 数値FT: L={}, ε={}) ===\n",
        a, big_l, eps
    );

    // W(Δ) 複素値 (iε 処方)
    let w_acc = |d: f64| -> (f64, f64) {
        // s = sinh((a/2)(Δ - iε)) = sinh(x)cos(y')... x=(a/2)Δ, y=-(a/2)ε
        let x = 0.5 * a * d;
        let y = -0.5 * a * eps;
        let (sre, sim) = (x.sinh() * y.cos(), x.cosh() * y.sin());
        // W = -a²/(16π s²)
        let s2re = sre * sre - sim * sim;
        let s2im = 2.0 * sre * sim;
        let den = s2re * s2re + s2im * s2im;
        let pref = -a * a / (16.0 * std::f64::consts::PI);
        (pref * s2re / den, -pref * s2im / den)
    };
    let w_inert = |d: f64| -> (f64, f64) {
        // W = -1/(4π (Δ-iε)²)
        let s2re = d * d - eps * eps;
        let s2im = -2.0 * d * eps;
        let den = s2re * s2re + s2im * s2im;
        let pref = -1.0 / (4.0 * std::f64::consts::PI);
        (pref * s2re / den, -pref * s2im / den)
    };
    // F(ω) = ∫ e^{-iωΔ} W(Δ) dΔ (台形)
    let response = |wf: &dyn Fn(f64) -> (f64, f64), om: f64| -> f64 {
        let nn = (2.0 * big_l / dt) as usize;
        let mut sum_re = 0.0;
        for i in 0..=nn {
            let d = -big_l + i as f64 * dt;
            let (wr, wi) = wf(d);
            let (c, s) = ((om * d).cos(), (om * d).sin());
            // Re[e^{-iωΔ} W] = c·wr + s·wi
            let val = c * wr + s * wi;
            let weight = if i == 0 || i == nn { 0.5 } else { 1.0 };
            sum_re += weight * val;
        }
        sum_re * dt
    };

    println!("[A] 慣性検出器: 真空は検出器を励起しない");
    for &om in &[0.5f64, 1.0, 2.0] {
        let f_pos = response(&w_inert, om);
        let f_neg = response(&w_inert, -om);
        println!(
            "  ω={:.1}: F(+ω) = {:+.2e} (→0),  F(-ω) = {:+.4} (自発放出は有限)",
            om, f_pos, f_neg
        );
    }

    println!("\n[B] 加速検出器: プランク分布と詳細釣り合い");
    println!("  ω     F(ω)        F(ω)/F(-ω)   e^{{-2πω/a}}   T_eff = ω/ln[F(-ω)/F(ω)]");
    let mut max_terr = 0.0f64;
    for &om in &[0.25f64, 0.5, 0.75, 1.0, 1.25] {
        let f_pos = response(&w_acc, om);
        let f_neg = response(&w_acc, -om);
        let ratio = f_pos / f_neg;
        let boltz = (-2.0 * std::f64::consts::PI * om / a).exp();
        let teff = om / (f_neg / f_pos).ln();
        max_terr = max_terr.max(
            (teff - a / (2.0 * std::f64::consts::PI)).abs() / (a / (2.0 * std::f64::consts::PI)),
        );
        println!(
            "  {:.2}  {:+.3e}  {:.5e}  {:.5e}  {:.5}",
            om, f_pos, ratio, boltz, teff
        );
    }
    let t_unruh = a / (2.0 * std::f64::consts::PI);
    println!(
        "  => 有効温度の理論値 T = a/2π = {:.5}, 最大相対誤差 {:.1e}  {}",
        t_unruh,
        max_terr,
        pass(max_terr < 0.01)
    );

    // プランク形状そのものの検証: F(ω)·(e^{2πω/a}-1)/ω = 一定?
    println!("\n[C] スペクトル形状: F(ω)(e^{{2πω/a}}-1)/ω は一定か (プランク分布の検証)");
    let mut vals = Vec::new();
    for &om in &[0.25f64, 0.5, 0.75, 1.0] {
        let f_pos = response(&w_acc, om);
        let v = f_pos * ((2.0 * std::f64::consts::PI * om / a).exp() - 1.0) / om;
        vals.push(v);
        println!("  ω={:.2}: {:.5}", om, v);
    }
    let vmean: f64 = vals.iter().sum::<f64>() / vals.len() as f64;
    let vdev = vals
        .iter()
        .map(|v| (v / vmean - 1.0).abs())
        .fold(0.0f64, f64::max);
    println!("  => 一定性のずれ {:.1e}  {}", vdev, pass(vdev < 0.01));
    println!(
        "\n結論: 同じミンコフスキー真空が、慣性系では「無」、加速系では「温度 a/2π の熱浴」。"
    );
    println!("      温度・粒子・真空は観測者に依存する — 「何が在るか」は運動状態の関数である。");
    println!("      地平線 (加速・BH・dS) はどれもこの機構で熱を持つ: Hawking 温度 κ/2π (v0.8)、");
    println!("      Gibbons–Hawking 温度 H/2π (v2.4 の (H/2π)² の起源)、モジュラー温度 (v0.7)。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
