//! v4.4 de Sitter 地平線と宇宙の情報収支
//!
//! [A] Gibbons–Hawking 温度: dS 空間の共動観測者は地平線から T = H/2π の熱を見る。
//!     共形場では応答関数の核が Unruh (v3.4) と同型 sinh² になる — 数値検証で三位一体を完結:
//!     加速地平線 a/2π (v3.4) = BH 地平線 κ/2π (v0.8) = 宇宙論的地平線 H/2π (本節)
//! [B] 観測可能宇宙のエントロピー台帳: 何がどれだけ情報を持っているか (プランク単位)
//! [C] ホログラフィック束縛: すべての内容物は dS 地平線の面積/4 に収まっているか

fn main() {
    println!("=== v4.4 de Sitter: 地平線の熱と宇宙の情報収支 ===\n");

    // ---- [A] dS 検出器 (共形真空, 核は sinh² — v3.4 と同じ数値 FT) ----
    let h = 1.0f64;
    let eps = 2e-3;
    let big_l = 120.0;
    let dt = 4e-4;
    let w_ds = |d: f64| -> (f64, f64) {
        let x = 0.5 * h * d;
        let y = -0.5 * h * eps;
        let (sre, sim) = (x.sinh() * y.cos(), x.cosh() * y.sin());
        let s2re = sre * sre - sim * sim;
        let s2im = 2.0 * sre * sim;
        let den = s2re * s2re + s2im * s2im;
        let pref = -h * h / (16.0 * std::f64::consts::PI);
        (pref * s2re / den, -pref * s2im / den)
    };
    let response = |om: f64| -> f64 {
        let nn = (2.0 * big_l / dt) as usize;
        let mut sum = 0.0;
        for i in 0..=nn {
            let d = -big_l + i as f64 * dt;
            let (wr, wi) = w_ds(d);
            let val = (om * d).cos() * wr + (om * d).sin() * wi;
            sum += if i == 0 || i == nn { 0.5 * val } else { val };
        }
        sum * dt
    };
    println!("[A] dS 共動検出器の応答 (H={}):", h);
    let mut max_err = 0.0f64;
    for &om in &[0.5f64, 1.0] {
        let fp = response(om);
        let fm = response(-om);
        let teff = om / (fm / fp).ln();
        max_err = max_err
            .max((teff - h / (2.0 * std::f64::consts::PI)).abs() * 2.0 * std::f64::consts::PI / h);
        println!(
            "  ω={:.1}: F(ω)/F(-ω) = {:.4e} (理論 e^{{-2πω/H}} = {:.4e}), T_eff = {:.5}",
            om,
            fp / fm,
            (-2.0 * std::f64::consts::PI * om / h).exp(),
            teff
        );
    }
    println!(
        "  => T = H/2π = {:.5} を相対誤差 {:.1e} で確認  {}",
        h / (2.0 * std::f64::consts::PI),
        max_err,
        pass(max_err < 0.01)
    );
    println!("     三位一体完結: 加速 a/2π (v3.4) = BH κ/2π (v0.8) = dS H/2π (本節)。");
    println!("     どの場合も「地平線 = 情報の遮蔽 = 熱」。\n");

    // ---- [B] エントロピー台帳 ----
    println!("[B] 観測可能宇宙のエントロピー台帳 (単位: k_B, 桁の見積り)");
    let h0 = 2.18e-18f64; // s^-1
    let t_p = 5.39e-44f64;
    let s_ds = std::f64::consts::PI / (h0 * t_p).powi(2);
    // CMB 光子
    let r_obs = 4.4e26f64; // m (共動半径)
    let vol = 4.0 / 3.0 * std::f64::consts::PI * r_obs.powi(3);
    let n_gamma = 4.11e8 * vol; // 411 /cm³
    let s_cmb = 3.6 * n_gamma;
    // ニュートリノ背景 ~ 同桁
    let s_nu = 5.2 * 3.0 / 4.0 * n_gamma;
    // 恒星質量 BH: ~10^8 個/銀河 × 10^11 銀河 × S(10 M_sun ~ 1e79)
    let s_stellar_bh = 1e19 * 1e79;
    // 超大質量 BH: ~10^11 個 × S(10^7 M_sun) (S ∝ M²: 1e79×(1e6)² = 1e91)
    let s_smbh = 1e11 * 1e93;
    println!("  CMB 光子           : {:.0e}", s_cmb);
    println!("  宇宙背景ニュートリノ : {:.0e}", s_nu);
    println!("  恒星質量ブラックホール: {:.0e}", s_stellar_bh);
    println!(
        "  超大質量ブラックホール: {:.0e}  ← 物質側の圧倒的王者",
        s_smbh
    );
    println!("  dS 地平線 (=宇宙の情報容量): {:.0e}", s_ds);
    println!(
        "  => 使用率 = S(BH)/S(dS) ~ {:.0e} — 宇宙は情報的にほぼ空 (まだ若い)",
        s_smbh / s_ds
    );
    println!("     エントロピーの主生産者はブラックホール = 「面積=情報」(v0.5) の天文学的実演\n");

    // ---- [C] ホログラフィック束縛 ----
    println!("[C] ホログラフィック束縛の検査: S(内容物) ≤ A/4 か");
    println!(
        "  S(全物質+BH) ~ {:.0e}  ≤  S_dS = {:.0e}  {}",
        s_smbh,
        s_ds,
        pass(s_smbh < s_ds)
    );
    println!("  ペンローズの初期条件問題: ビッグバンのエントロピー ~10^88 は最大値 ~10^122 の");
    println!(
        "  10^-34。初期宇宙は途方もなく「整った」状態だった — 時間の矢 (v1.1) の宇宙論的根源。"
    );
    println!("  これは QRN 最大の未解決問題として残る: なぜ初期状態のもつれは少なかったのか。");
    println!("\n結論: 宇宙論的地平線も熱を持ち (H/2π)、宇宙の情報収支はホログラフィック束縛に");
    println!("      従う。宇宙は情報容量の ~10^-18 しか使っておらず、進化の余地は膨大である。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
