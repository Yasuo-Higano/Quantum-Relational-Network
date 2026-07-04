//! v2.4 宇宙構造の量子起源 — インフレーション摂動の厳密数値
//!
//! 銀河・CMB の温度ゆらぎは「引き伸ばされた真空の量子ゆらぎ」である (現代宇宙論の中心結果)。
//! Mukhanov–Sasaki 方程式  v_k'' + (k² - a''/a) v_k = 0  (共形時間 η)
//! を Bunch–Davies 真空から数値積分し:
//!   [A] de Sitter (a=-1/Hη): 厳密解と比較 (コード検証)
//!   [B] スペクトル P(k) = k³|v/a|²/(2π²) → (H/2π)² : スケール不変性
//!   [C] 準 dS (ν=3/2+δ): 傾き n_s-1 = 3-2ν — 観測される赤い傾き (Planck: n_s=0.965)
//!   [D] スクイージング: ホライズン通過後にゆらぎが凍結・増幅され古典化へ (v1.6 の機構)

use uft_sim::*;

/// v'' = -(k² - c/η²) v を η0→η1 まで RK4 積分 (複素 v)
fn evolve(k: f64, c: f64, eta0: f64, eta1: f64) -> (C64, C64) {
    // BD 初期条件 (平面波 + dS では厳密形)
    let x0 = k * eta0;
    let sq = 1.0 / (2.0 * k).sqrt();
    let mut v = if (c - 2.0).abs() < 1e-12 {
        // dS 厳密: v = e^{-ikη}(1 - i/(kη))/√(2k)
        (C64::expi(-x0) * C64::new(1.0, -1.0 / x0)).scale(sq)
    } else {
        C64::expi(-x0).scale(sq)
    };
    let mut u = if (c - 2.0).abs() < 1e-12 {
        (C64::expi(-x0)
            * (C64::new(0.0, -k) * C64::new(1.0, -1.0 / x0)
                + C64::new(0.0, 1.0 / (k * eta0 * eta0))))
        .scale(sq)
    } else {
        (C64::expi(-x0) * C64::new(0.0, -k)).scale(sq)
    };
    let mut eta = eta0;
    let acc = |v: C64, eta: f64| -> C64 { v.scale(-(k * k - c / (eta * eta))) };
    while eta < eta1 {
        let h = (0.015 / k).min((-eta) / 400.0).min(eta1 - eta);
        let k1 = (u, acc(v, eta));
        let k2 = (
            u + k1.1.scale(0.5 * h),
            acc(v + k1.0.scale(0.5 * h), eta + 0.5 * h),
        );
        let k3 = (
            u + k2.1.scale(0.5 * h),
            acc(v + k2.0.scale(0.5 * h), eta + 0.5 * h),
        );
        let k4 = (u + k3.1.scale(h), acc(v + k3.0.scale(h), eta + h));
        v = v + (k1.0 + k2.0.scale(2.0) + k3.0.scale(2.0) + k4.0).scale(h / 6.0);
        u = u + (k1.1 + k2.1.scale(2.0) + k3.1.scale(2.0) + k4.1).scale(h / 6.0);
        eta += h;
    }
    (v, u)
}

fn main() {
    println!("=== v2.4 宇宙構造の量子起源: インフレーションのゆらぎ (H=1) ===\n");

    // ---- [A] de Sitter 厳密解との比較 ----
    println!("[A] dS (a''/a = 2/η²): kη = -100 → -0.01 の数値積分 vs 厳密解");
    let mut max_rel = 0.0f64;
    for &k in &[1.0f64, 5.0, 20.0] {
        let (v, _) = evolve(k, 2.0, -100.0 / k, -0.01 / k);
        let x = -0.01;
        let exact = (1.0 / (2.0 * k)) * (1.0 + 1.0 / (x * x));
        let rel = ((v.norm2() - exact) / exact).abs();
        max_rel = max_rel.max(rel);
        println!(
            "  k={:5.1}: |v|²(数値)={:.6e}  厳密={:.6e}  相対誤差 {:.1e}",
            k,
            v.norm2(),
            exact,
            rel
        );
    }
    println!(
        "  => 最大相対誤差 {:.1e}  {}\n",
        max_rel,
        pass(max_rel < 1e-6)
    );

    // ---- [B] スケール不変スペクトル ----
    // 全モードを「インフレーション終了」の共通時刻 η_end で評価する (物理的に正しい比較。
    // 固定 kη での評価は自己相似で任意の ν が平坦に見えてしまう)
    let eta_end = -0.002f64;
    println!(
        "[B] スペクトル (共通の終了時刻 η_end={} で評価): P_δφ(k)/(H/2π)²",
        eta_end
    );
    println!("  k        P/(H/2π)²");
    let mut lnk = Vec::new();
    let mut lnp = Vec::new();
    for i in 0..10 {
        let k = 0.5 * (10.0f64).powf(i as f64 / 9.0); // 0.5 → 5 (1 桁)
        let (v, _) = evolve(k, 2.0, -100.0 / k, eta_end);
        let ratio = 2.0 * k.powi(3) * v.norm2() * eta_end * eta_end; // P/(H/2π)²
        println!("  {:7.3}  {:.6}", k, ratio);
        lnk.push(k.ln());
        lnp.push(ratio.ln());
    }
    let (_, ns1) = linfit(&lnk, &lnp);
    println!(
        "  => 傾き n_s - 1 = {:+.5} (dS の理論値 0: スケール不変)  {}\n",
        ns1,
        pass(ns1.abs() < 2e-3)
    );

    // ---- [C] 準 de Sitter: 赤い傾き ----
    println!("[C] 準 dS (a''/a = (ν²-¼)/η², ν = 1.55): 理論 n_s-1 = 3-2ν = -0.10");
    let nu: f64 = 1.55;
    let c = nu * nu - 0.25;
    let mut lnk = Vec::new();
    let mut lnp = Vec::new();
    for i in 0..10 {
        let k = 0.5 * (10.0f64).powf(i as f64 / 9.0);
        let (v, _) = evolve(k, c, -300.0 / k, eta_end);
        lnk.push(k.ln());
        lnp.push((2.0 * k.powi(3) * v.norm2() * eta_end * eta_end).ln());
    }
    let (_, ns1) = linfit(&lnk, &lnp);
    println!(
        "  数値: n_s - 1 = {:+.4}  (理論 {:+.4})  {}",
        ns1,
        3.0 - 2.0 * nu,
        pass((ns1 - (3.0 - 2.0 * nu)).abs() < 5e-3)
    );
    println!("  観測 (Planck 2018): n_s = 0.9649 ± 0.0042 — わずかに赤い = 膨張がわずかに減速");
    println!("  => 10^26 倍の引き伸ばしの「わずかな非一様さ」が、CMB の傾きとして今も見えている\n");

    // ---- [D] スクイージングと古典化 ----
    println!("[D] スクイージング: 真空ゆらぎとの振幅比 √(2k)|v| の成長 (k=1, dS)");
    for &etaf in &[-2.0f64, -0.5, -0.1, -0.02] {
        let (v, _) = evolve(1.0, 2.0, -100.0, etaf);
        println!(
            "  kη = {:5.2}: √(2k)|v| = {:8.2} (ホライズン外で ∝ 1/|kη| に凍結・増幅)",
            etaf,
            (2.0f64 * v.norm2()).sqrt()
        );
    }
    println!("  => 増幅されたゆらぎは強くスクイーズされた状態 = 実質的に古典確率場。");
    println!("     デコヒーレンス (v1.6) が引き継ぎ、密度ゆらぎ→重力崩壊→銀河へ (v0.1 の重力)。");
    println!("\n結論: 宇宙で最大の構造 (銀河・CMB の斑) は、最小のもの (真空のもつれ) の拡大写真である。");
    println!("      QRN の A2 (ゆらぎ=情報) と A6 (宇宙論) を繋ぐ、観測済みの量子重力現象。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
