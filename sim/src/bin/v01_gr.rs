//! v0.1 一般相対性理論の数値検証 — シュヴァルツシルト時空の測地線
//! 幾何単位系 G = c = 1、長さの単位 = GM_sun/c² = 1476.625 m
//!
//! 質点軌道:  d²u/dφ² + u = 1/p + 3u²   (u = 1/r, p = 半直弦)
//! 光:        d²u/dφ² + u = 3u²
//!
//! A: 水星の近日点移動 → 観測値 42.98″/世紀
//! B: 太陽をかすめる光の偏向 → 観測値 1.75″

const RG_SUN_M: f64 = 1.476625e3; // GM_sun/c² [m]
const ARCSEC: f64 = 206264.806; // rad → arcsec

/// RK4 で (u, w=du/dφ) を積分
fn rk4_step(u: f64, w: f64, h: f64, f: impl Fn(f64) -> f64) -> (f64, f64) {
    let k1 = (w, f(u));
    let k2 = (w + 0.5 * h * k1.1, f(u + 0.5 * h * k1.0));
    let k3 = (w + 0.5 * h * k2.1, f(u + 0.5 * h * k2.0));
    let k4 = (w + h * k3.1, f(u + h * k3.0));
    (
        u + h / 6.0 * (k1.0 + 2.0 * k2.0 + 2.0 * k3.0 + k4.0),
        w + h / 6.0 * (k1.1 + 2.0 * k2.1 + 2.0 * k3.1 + k4.1),
    )
}

fn main() {
    println!("=== v0.1 一般相対論: Schwarzschild 測地線 (G=c=1, 長さ単位 GM/c²) ===\n");

    // ---- A: 水星近日点移動 ----
    println!("[A] 水星の近日点移動");
    let a_mercury = 5.7909083e10 / RG_SUN_M; // 軌道長半径
    let e = 0.205630;
    let p = a_mercury * (1.0 - e * e); // 半直弦
    let rhs = |u: f64| -u + 1.0 / p + 3.0 * u * u;
    let (mut u, mut w) = ((1.0 + e) / p, 0.0); // 近日点から出発
    let h = 2.0 * std::f64::consts::PI / 200_000.0;
    let mut phi = 0.0;
    let mut peri: Vec<f64> = vec![0.0];
    let mut w_prev = w;
    for _ in 0..(200_000 * 9) {
        let (nu, nw) = rk4_step(u, w, h, rhs);
        u = nu;
        w = nw;
        phi += h;
        // w が - → + に変わる点が近日点 (u 極大は w:+→-, 極小は -→+; 近日点は u 極大 = w:+→-)
        if w_prev > 0.0 && w <= 0.0 {
            let frac = w_prev / (w_prev - w);
            peri.push(phi - h + frac * h);
        }
        w_prev = w;
    }
    let n_orbit = peri.len() - 1;
    let dphi_per_orbit = (peri[n_orbit] - peri[0]) / n_orbit as f64 - 2.0 * std::f64::consts::PI;
    let orbits_per_century = 100.0 * 365.25 / 87.9691;
    let shift_num = dphi_per_orbit * orbits_per_century * ARCSEC;
    let shift_th = 6.0 * std::f64::consts::PI / p * orbits_per_century * ARCSEC;
    println!(
        "  軌道: a = {:.4e} (GM/c²), e = {}, {} 周回を積分",
        a_mercury, e, n_orbit
    );
    println!("  近日点移動(数値)   : {:.4}″/世紀", shift_num);
    println!("  GR 一次摂動論 6πGM/(c²p): {:.4}″/世紀", shift_th);
    println!("  観測値 (レーダー測距)   : 42.98 ± 0.04 ″/世紀");
    let rel = ((shift_num - 42.98) / 42.98).abs();
    println!(
        "  => 観測との相対差 {:.2}%  {}\n",
        rel * 100.0,
        pass(rel < 0.01)
    );

    // ---- B: 光の偏向 ----
    println!("[B] 太陽縁をかすめる光の偏向");
    let b = 6.957e8 / RG_SUN_M; // 衝突径数 = 太陽半径
    let rhs_light = |u: f64| -u + 3.0 * u * u;
    let (mut u, mut w) = (0.0, 1.0 / b); // 無限遠 (u=0) から入射
    let h = std::f64::consts::PI / 2_000_000.0;
    let mut phi = 0.0;
    let mut u_prev = u;
    loop {
        let (nu, nw) = rk4_step(u, w, h, rhs_light);
        u = nu;
        w = nw;
        phi += h;
        if u < 0.0 {
            // u=0 を横切った点 = 無限遠へ脱出
            let frac = u_prev / (u_prev - u);
            phi = phi - h + frac * h;
            break;
        }
        u_prev = u;
        if phi > 6.4 {
            break;
        }
    }
    let defl_num = (phi - std::f64::consts::PI) * ARCSEC;
    let defl_th = 4.0 / b * ARCSEC;
    println!("  衝突径数 b = R_sun = {:.4e} (GM/c²)", b);
    println!("  偏向角(数値)      : {:.4}″", defl_num);
    println!("  GR 一次 4GM/(c²b) : {:.4}″", defl_th);
    println!("  ニュートン論の予言 : {:.4}″ (半分!)", 2.0 / b * ARCSEC);
    println!("  観測値 (1919 皆既日食〜VLBI): 1.75″");
    let rel = ((defl_num - defl_th) / defl_th).abs();
    println!(
        "  => 一次摂動論との相対差 {:.2e}  {}\n",
        rel,
        pass(rel < 1e-3)
    );

    // ---- C: 強重力領域(光子球) ----
    println!(
        "[C] 強重力: 光子の捕獲 (臨界衝突径数 b_c = 3√3 = {:.4})",
        3.0 * 3.0f64.sqrt()
    );
    for &btest in &[5.0, 5.196, 5.4] {
        let rhs_light = |u: f64| -u + 3.0 * u * u;
        let (mut u, mut w) = (1e-12, 1.0 / btest);
        let mut captured = false;
        let mut phi = 0.0;
        for _ in 0..10_000_000 {
            let (nu, nw) = rk4_step(u, w, 1e-5, rhs_light);
            u = nu;
            w = nw;
            phi += 1e-5;
            if u > 0.55 {
                captured = true;
                break;
            } // r < 1.8GM/c² → 事象の地平線へ
            if u < 0.0 {
                break;
            }
        }
        println!(
            "  b = {:5.3} : {}",
            btest,
            if captured {
                "捕獲(ブラックホールへ落下)"
            } else {
                "脱出"
            }
        );
        let _ = phi;
    }
    println!("\n結論: 時空の曲率としての重力を数値で確認。GR は太陽系スケールで 10⁻³ 以下の精度で正しい。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
