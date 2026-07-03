//! v1.5 宇宙論への接続 — ΛCDM の再現と「常在 Λ」ゆらぎ模型
//!
//! [A] 確立部分: Friedmann 方程式で標準宇宙論の膨張史を再現 (較正)
//! [B] 観測的事実: Λ(観測) ≈ 1/√V₄(宇宙の4体積, プランク単位) — 因果集合が
//!     「事前に」オーダーを当てた唯一の予言 (Sorkin 1990s)
//! [C] 投機的トイ模型 (明確にラベル): 常在Λ (everpresent Λ, Ahmed-Dodelson-Greene-Sorkin
//!     2004 に触発): Λ が ±1/√V₄ で揺らぐと、Λ は全時代で物質密度と同オーダーを
//!     保ち続け「なぜ今問題」が消える。符号は揺らぐ → w(z) の揺らぎという検証可能な帰結。

use uft_sim::*;

fn main() {
    println!("=== v1.5 宇宙論: ΛCDM の再現と宇宙定数のゆらぎ仮説 ===\n");

    // ---- [A] ΛCDM 較正 (確立された物理) ----
    println!("[A] ΛCDM 膨張史 (Ω_m=0.315, Ω_r=9.15e-5, Ω_Λ=残り, h=0.674)");
    let om = 0.315f64;
    let or_ = 9.15e-5f64;
    let ol = 1.0 - om - or_;
    let e_of_a = |a: f64| (om / (a * a * a) + or_ / (a * a * a * a) + ol).sqrt();
    // t0 H0 = ∫0^1 da / (a E(a))  (シンプソン, ln a 変数)
    let nn = 200_000;
    let (la0, la1) = ((1e-10f64).ln(), 0.0f64);
    let dla = (la1 - la0) / nn as f64;
    let mut t0h0 = 0.0;
    for i in 0..=nn {
        let w = if i == 0 || i == nn {
            1.0
        } else if i % 2 == 1 {
            4.0
        } else {
            2.0
        };
        let a = (la0 + i as f64 * dla).exp();
        t0h0 += w * 1.0 / e_of_a(a);
    }
    t0h0 *= dla / 3.0;
    let hubble_time_gyr = 9.7779 / 0.674;
    println!("  宇宙年齢: t0·H0 = {:.4} → t0 = {:.2} Gyr   (Planck 観測: 13.80 Gyr)  {}",
        t0h0, t0h0 * hubble_time_gyr, pass((t0h0 * hubble_time_gyr - 13.80).abs() < 0.1));
    // 加速開始 (q=0): Ω_m a^-3 + 2Ω_r a^-4 = 2Ω_Λ
    let mut a_acc = 0.5f64;
    for _ in 0..100 {
        let f = om / a_acc.powi(3) + 2.0 * or_ / a_acc.powi(4) - 2.0 * ol;
        let df = -3.0 * om / a_acc.powi(4) - 8.0 * or_ / a_acc.powi(5);
        a_acc -= f / df;
    }
    println!("  減速→加速の転換: z = {:.3}   (SNe Ia 観測: z ≈ 0.6±0.1)  {}",
        1.0 / a_acc - 1.0, pass((1.0 / a_acc - 1.0 - 0.63).abs() < 0.1));
    println!("  物質-輻射等密度: z_eq = {:.0}   (Planck: ≈3400)  {}\n",
        om / or_ - 1.0, pass(((om / or_ - 1.0) - 3400.0).abs() < 100.0));

    // ---- [B] 数のゆらぎ仮説の算術 (観測的事実の指摘) ----
    println!("[B] 宇宙定数の値と 4 体積ゆらぎ (プランク単位の算術)");
    let h0_planck = 1.18e-61f64; // H0 in 1/t_P
    let lam_obs = 3.0 * ol * h0_planck * h0_planck; // Λ = 3 Ω_Λ H0²
    let t0_planck = t0h0 / h0_planck;
    let v4 = t0_planck.powi(4); // 過去光円錐 4 体積 ~ t0^4 (係数 O(1) は略)
    println!("  Λ(観測)              = {:.2e}  (プランク単位)", lam_obs);
    println!("  1/√V₄ (V₄ ~ t0⁴)     = {:.2e}", 1.0 / v4.sqrt());
    println!("  比                    = {:.1}", lam_obs * v4.sqrt());
    println!("  => 122 桁の算術で一致 (オーダー)。「Λ は 4 体積 N 個の要素の √N ゆらぎ」という");
    println!("     因果集合の予言 (Sorkin) は、加速膨張発見 *前* にこの値を指していた。\n");

    // ---- [C] 常在 Λ トイ模型 (投機的 — QRN の A2 と因果集合的離散性の帰結の候補) ----
    println!("[C] 常在Λトイ模型: Λ(t) = ξ·W(N)/V₄, W はランダムウォーク (⟨Λ²⟩^½ = ξ/√V₄)");
    println!("    物質優勢宇宙で t を 12 桁進め、比 r = |Λ/3| / (8πρ_m/3) を追う:");
    let mut rng = Rng::new(42);
    let xi = 1.0f64;
    let mut t = 1.0f64;
    let mut a = 1.0f64;
    let mut w = 0.0f64;
    let mut v4t = 1.0f64; // V4 ~ t^4 (過去光円錐)
    let rho0 = 1.0f64;
    let steps_per_efold = 200;
    let mut ratios = Vec::new();
    let mut neg_frac = 0usize;
    let mut cnt = 0usize;
    println!("  log10(t)   log10(ρ_m)   Λ/(8πρ_m)  (時刻サンプル)");
    let mut next_print = 1.0f64;
    while t < 1e12 {
        let dt = t / steps_per_efold as f64;
        let v4_new = (t + dt).powi(4);
        let dn = v4_new - v4t;
        w += rng.gauss() * dn.sqrt();
        v4t = v4_new;
        let lam = xi * w / v4t;
        let rho_m = rho0 / (a * a * a);
        let h2 = (8.0 * std::f64::consts::PI / 3.0) * rho_m + lam / 3.0;
        let h = h2.max(1e-30).sqrt();
        a *= 1.0 + h * dt;
        t += dt;
        let r = (lam / 3.0) / ((8.0 * std::f64::consts::PI / 3.0) * rho_m);
        ratios.push(r.abs());
        if lam < 0.0 {
            neg_frac += 1;
        }
        cnt += 1;
        if t >= next_print {
            println!("  {:7.2}    {:8.2}     {:+8.4}", t.log10(), rho_m.log10(), r);
            next_print *= 100.0;
        }
    }
    let mut sorted = ratios.clone();
    sorted.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let med = sorted[sorted.len() / 2];
    let in_band = sorted.iter().filter(|&&x| x > 0.001 && x < 10.0).count();
    println!("  => 密度が 36 桁変わる間、|Λ|/ρ_m の中央値 = {:.3}、[0.001,10] に留まる割合 {:.0}%",
        med, 100.0 * in_band as f64 / sorted.len() as f64);
    println!("     Λ が負の時間帯 {:.0}% (加速と減速が交代)", 100.0 * neg_frac as f64 / cnt as f64);
    println!("     対照: 定数Λなら比は 10^-36 → 1 を通過する一瞬しか同程度にならない (なぜ今問題)");
    println!("\n結論: [A] 標準宇宙論は再現できる (較正)。[B] Λ の観測値は 1/√V₄ と一致 — 偶然か手がかりか。");
    println!("      [C] Λ が離散時空の統計ゆらぎなら「なぜ今」は消え、w(z) の揺らぎという検証可能な");
    println!("      予言を持つ。ただし CMB 時代の制約など課題は多い — 投機的段階と明記する。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
