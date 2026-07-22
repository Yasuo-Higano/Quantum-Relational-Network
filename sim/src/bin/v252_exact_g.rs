//! v25.2 g(μ) の厳密閉形式 — 既知厳密解 (Eisler) の取り込みと事前登録判定 (第二十六期)
//!
//! v25.1 は 3+1D staggered 半空間の λ を 1D staggered 鎖のモジュラー核 prefactor
//! g(μ), g_m(μ) の BZ 平均に厳密還元した (判定 a)。本版はその解析形を確定する。
//!
//! **先行研究 (C0)**: V. Eisler, "On the Bisognano-Wichmann entanglement Hamiltonian
//! of nonrelativistic fermions", J. Stat. Mech. (2025) 013101 [arXiv:2410.16433] は、
//! 半無限 staggered 鎖 Ĥ = −(1/2)Σ(c†c+h.c.) + μΣ(−1)ⁿ c†c の entanglement
//! Hamiltonian を厳密に
//!   H_EH = 4κ K(κ′) T,  T: t_m = −m/2, d_m = (−1)^m μ(m−1/2),
//!   κ = 1/√(1+μ²),  κ′ = |μ|/√(1+μ²)   (同論文 式 (64)(67)(72))
//! と与えた ([C,T]=0 の可換三重対角行列 + 留数和による Fermi 関数の再現)。
//! 本鎖 (hopping +1/2, v251 の chain_k) とは局所ユニタリ cₙ → (−1)ⁿcₙ で同値
//! (staggered 項は不変、ボンド符号のみ反転 — K の勾配絶対値は不変)。
//!
//! **QRN 規格化への翻訳** (v25.1 定義 K_NN = g·πξ + b, (−1)^x K_xx = g_m·2πξ_s·μ + b_m):
//!   πg = 2κK(κ′) かつ 2π·g_m·μ = 4κK(κ′)·μ  ⟹  **g(μ) = g_m(μ) = (2/π)κK(κ′)**。
//! 運動項と質量項の規格化が同一関数であること (v25.1 の観測 g_m ≈ g) は定理になる。
//!
//! 三表示 (数学的同値 — 導出は docs/uft-v25.2.md §1):
//!   [楕円]  g = (2/π)·K(κ′)/√(1+μ²)
//!   [AGM]   g = 1/AGM(1, √(1+μ²))          (K(k) = π/(2 AGM(1,k′)) と AGM 同次性)
//!   [積分]  g = (2/π)∫₀^{π/2} dθ/√(1+μ²cos²θ)  (Legendre 標準形 a=√(1+μ²), b=1)
//! 漸近形 (v25.1 §4 の開問への解答):
//!   小 μ: g = 1 − μ²/4 + 9μ⁴/64 − 25μ⁶/256 + 1225μ⁸/16384 − …  — μ² の解析級数
//!         (係数 (−1)ⁿ[(2n−1)!!/(2n)!!]²)。**μ²lnμ 型の非解析項は存在しない**。
//!   大 μ: g ≃ (2/(π|μ|))·ln(4|μ|) — 1/μ² 型でも純 1/μ 型でもなく log 補正つき 1/μ。
//!
//! 検査 (設計ノートの事前登録プロトコル: 判定は μ ≥ 0.5 部分集合 — 小 μ は
//! ξ_corr = 1/arsinh(μ) が N=256 の認証窓を超え crossover 系統が乗るため):
//!  [S1] 三表示の相互一致 (μ ∈ [1e-4, 1e4]) + 外部照合 g(1) = ガウス定数 1/AGM(1,√2)
//!  [S2] 小 μ 級数: 残差 g−S₄ = c₅μ¹⁰(1−(121/144)μ²+…) の両側検証 + 残差スロープ ≈ 10
//!  [S3] 旧仮説の棄却: (i) (g−1)/μ² → −1/4 (μ²lnμ なら log 発散) /
//!       (ii) πμg/(2ln(4√(1+μ²))) → 1 (1/μ² なら → 0, 純 c/μ なら μg が定数のはず)
//!  [S4] 事前登録判定 (μ≥0.5, v251_gmu.json): max|ĝ/g_exact − 1| < 5e-5
//!  [S5] g_m ≡ g: max|ĝ_m/g_exact − 1| < 5e-5 (μ≥0.5) — 同一 prefactor の定理の実測照合
//!  [S6] 性質: g(0)=1 (bit 一致)・偶関数・0<g≤1・狭義単調減少・g′<0 (積分表示と数値微分)
//!  [S7] 変異検出層 (破壊層): κ↔κ′ / 規格 4/π / 母数 κ′² / AGM(1,1+μ²) / 旧最良候補
//!       π/(2K(κ′)) — 全てゲートが検出すること (設計ノートの棄却 0.16%→3.8% の再現込み)
//!
//! 事前登録分岐: (a) S1–S7 全 PASS → **g の解析形確定** (楕円積分 — CTM 構造の見込みが的中、
//!   ただし閉形式自体は Eisler の既知結果 = 新規性は同定と接続のみ) /
//!   (b) S4/S5 のみ FAIL → 閉形式は自己整合だが実測と不一致 = 有限 N 系統の見積り誤り (記録) /
//!   (c) S1 FAIL → 実装または同定の誤り (器械再監査)。
//!
//! 新規性境界: 1D 閉形式 = Eisler (2024) の既知結果 (C0)。本版の寄与 = v25.1 実測との
//! 同定 (C1)・旧仮説 (μ²lnμ / 1/μ²) の棄却・BZ moment 解析化 (後続コミット) への接続。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;
/// ガウス定数 1/AGM(1,√2) (外部照合値 — OEIS A014549)
const GAUSS: f64 = 0.834_626_841_674_073_186_3;

/// 周期 2π の解析的周期関数の一様格子平均 (中点則、倍増収束 — 幾何収束)。
/// 戻り値 (平均, 直前レベルとの差, 使用 n)。fail-closed: 収束せねば est が残る。
fn circle_mean(f: impl Fn(f64) -> f64, tol: f64) -> (f64, f64, usize) {
    let mut n = 32usize;
    let mut prev = f64::NAN;
    loop {
        let mut s = 0.0;
        for i in 0..n {
            s += f(2.0 * PI * (i as f64 + 0.5) / n as f64);
        }
        let cur = s / n as f64;
        let est = (cur - prev).abs();
        if est <= tol * cur.abs().max(1e-300) || n >= (1 << 22) {
            return (cur, est, n);
        }
        prev = cur;
        n *= 2;
    }
}

/// 表示 1 [積分]: g = ⟨(1+μ²cos²θ)^{−1/2}⟩_circle (全周平均 = (2/π)∫₀^{π/2})
fn g_repr_int(mu: f64) -> f64 {
    let m2 = mu * mu;
    circle_mean(|t| 1.0 / (1.0 + m2 * t.cos() * t.cos()).sqrt(), 1e-15).0
}

/// 表示 2 [楕円]: g = (2/π)κK(κ′)、K は Legendre 形の直接求積 (AGM 不使用)
fn g_repr_ell(mu: f64) -> f64 {
    let kap = 1.0 / (1.0 + mu * mu).sqrt();
    let kp2 = mu * mu / (1.0 + mu * mu); // κ′²
                                         // K(κ′) = ∫₀^{π/2} dθ/√(1−κ′²sin²θ) = (π/2)·⟨…⟩_circle
    let kk = (PI / 2.0) * circle_mean(|t| 1.0 / (1.0 - kp2 * t.sin() * t.sin()).sqrt(), 1e-15).0;
    (2.0 / PI) * kap * kk
}

/// 旧最良候補 (設計ノートで棄却済み): π/(2K(κ′)) — T の固有値間隔関数であり
/// 演算子 prefactor ではない。S7 で棄却の再現に使う。
fn g_old_candidate(mu: f64) -> f64 {
    let kp = (mu * mu / (1.0 + mu * mu)).sqrt();
    PI / (2.0 * ellip_k(kp).unwrap())
}

/// 小 μ 級数 S₄ (μ⁸ まで)
fn g_series4(mu: f64) -> f64 {
    let x = mu * mu;
    1.0 - x / 4.0 + 9.0 * x * x / 64.0 - 25.0 * x * x * x / 256.0 + 1225.0 * x * x * x * x / 16384.0
}

/// v251_gmu.json の (mu, g, g_m) テーブルを読む (自作最小パーサ — 自リポジトリの
/// write_artifact 出力形式に限定)。読めなければ None (fail-closed)。
fn load_v251_table() -> Option<Vec<(f64, f64, f64)>> {
    let txt = std::fs::read_to_string("results/v251_gmu.json")
        .or_else(|_| std::fs::read_to_string("../results/v251_gmu.json"))
        .ok()?;
    let num_after = |chunk: &str, key: &str| -> Option<f64> {
        let p = chunk.find(key)? + key.len();
        let rest = &chunk[p..];
        let end = rest
            .find(|c: char| c == ',' || c == '}' || c == '\n')
            .unwrap_or(rest.len());
        rest[..end].trim().parse::<f64>().ok()
    };
    let mut out = Vec::new();
    for chunk in txt.split("\"mu\":").skip(1) {
        // chunk 先頭が mu の値。同一オブジェクト内に g, b, g_m が続く。
        let end = chunk
            .find(|c: char| c == ',' || c == '}')
            .unwrap_or(chunk.len());
        let mu = chunk[..end].trim().parse::<f64>().ok()?;
        let obj_end = chunk.find('}').unwrap_or(chunk.len());
        let obj = &chunk[..obj_end];
        let g = num_after(obj, "\"g\":")?;
        let gm = num_after(obj, "\"g_m\":")?;
        out.push((mu, g, gm));
    }
    if out.len() >= 10 {
        Some(out)
    } else {
        None
    }
}

fn main() {
    self_test();
    println!(
        "=== v25.2 g(μ) の厳密閉形式 — Eisler 厳密解の取り込みと事前登録判定 (第二十六期) ===\n"
    );
    println!(
        "事前登録: (a) S1–S7 全 PASS → 解析形確定 (閉形式は既知 [Eisler 2024]、新規性は同定) /"
    );
    println!("          (b) S4/S5 のみ FAIL → 有限 N 系統の見積り誤り (記録) / (c) S1 FAIL → 器械再監査\n");

    let mut stage_fail = std::collections::BTreeMap::<String, usize>::new();
    let mut checks: Vec<(String, bool)> = Vec::new();
    let mut check = |name: &str, ok: bool, detail: String| {
        println!(
            "  [{}] {}  {}",
            if ok { "PASS" } else { "FAIL" },
            name,
            detail
        );
        let stage = name
            .trim_start_matches('[')
            .split(|c| c == ']' || c == '-')
            .next()
            .unwrap_or("?")
            .to_string();
        if !ok {
            *stage_fail.entry(stage).or_insert(0) += 1;
        }
        checks.push((name.to_string(), ok));
    };

    // ---- [S1] 三表示の相互一致 + 外部照合 ----
    {
        let mut mus: Vec<f64> = Vec::new();
        let mut e = -4.0f64;
        while e <= 4.0001 {
            mus.push(10f64.powf(e));
            e += 0.25;
        }
        mus.extend_from_slice(&[0.5, 1.0, 2.0f64.sqrt()]);
        let (mut dmax_lo, mut dmax_hi) = (0.0f64, 0.0f64); // μ≤100 / μ>100
        for &mu in &mus {
            let a = g_exact(mu).unwrap();
            let b = g_repr_ell(mu);
            let c = g_repr_int(mu);
            let d = ((a - b).abs() / a)
                .max((a - c).abs() / a)
                .max((b - c).abs() / a);
            if mu <= 100.0 {
                dmax_lo = dmax_lo.max(d);
            } else {
                dmax_hi = dmax_hi.max(d);
            }
        }
        check(
            "[S1a] 三表示一致 (μ ≤ 100)",
            dmax_lo < 1e-12,
            format!(
                "max 相対差 = {:.1e} (AGM/楕円/積分, {} 点)",
                dmax_lo,
                mus.len()
            ),
        );
        check(
            "[S1b] 三表示一致 (μ > 100, Legendre 形 1−κ′²sin²θ の桁落ち許容域)",
            dmax_hi < 1e-8,
            format!("max 相対差 = {:.1e}", dmax_hi),
        );
        let g1 = g_exact(1.0).unwrap();
        check(
            "[S1c] 外部照合: g(1) = ガウス定数 1/AGM(1,√2)",
            (g1 - GAUSS).abs() < 5e-16,
            format!(
                "g(1) = {:.16} (照合値 {:.16}, Δ = {:.1e})",
                g1,
                GAUSS,
                g1 - GAUSS
            ),
        );
    }

    // ---- [S2] 小 μ 級数 (μ¹⁰ 残差の両側検証) ----
    {
        let c5 = -(63.0f64 / 256.0) * (63.0 / 256.0); // (−1)⁵[(9!!)/(10!!)]² = −(63/256)²
        let mus = [0.10, 0.15, 0.20, 0.25, 0.30];
        let mut ok_ratio = true;
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        let mut worst = 0.0f64;
        for &mu in &mus {
            let r = g_exact(mu).unwrap() - g_series4(mu);
            let ratio = r / (c5 * mu.powi(10));
            let pred = 1.0 - (121.0 / 144.0) * mu * mu; // 次係数比 c₆/c₅ = −(11/12)²
            let dev = (ratio - pred).abs();
            worst = worst.max(dev);
            if dev > 0.01 {
                ok_ratio = false;
            }
            xs.push(mu.ln());
            ys.push(r.abs().ln());
        }
        check(
            "[S2a] 残差 = c₅μ¹⁰(1−(121/144)μ²) の両側一致 (μ ∈ [0.1,0.3])",
            ok_ratio,
            format!("max|ratio−pred| = {:.1e} (< 0.01)", worst),
        );
        let (_, slope) = linfit_checked(&xs, &ys).unwrap();
        check(
            "[S2b] 残差スロープ d ln|g−S₄|/d ln μ ≈ 10",
            (9.7..10.05).contains(&slope),
            format!("slope = {:.3} (期待 10 − O(μ²))", slope),
        );
    }

    // ---- [S3] 旧仮説の棄却 ----
    {
        let q = |mu: f64| (g_exact(mu).unwrap() - 1.0) / (mu * mu);
        let (qa, qb) = (q(1e-3), q(1e-4));
        check(
            "[S3a] (g−1)/μ² → −1/4 (μ² の解析級数 — μ²lnμ 型は不在)",
            (qb + 0.25).abs() < 1e-6,
            format!(
                "(g−1)/μ²|₁ₑ₋₄ = {:.9} (−0.25 との差 {:.1e})",
                qb,
                (qb + 0.25).abs()
            ),
        );
        check(
            "[S3b] log 走行なし: |q(1e-3) − q(1e-4)| < 1e-5 → μ²lnμ 係数 |a| < 5e-6",
            (qa - qb).abs() < 1e-5,
            format!(
                "Δq = {:.2e} (μ²lnμ 仮説なら a·ln10 = 2.3a が残る)",
                (qa - qb).abs()
            ),
        );
        let t =
            |mu: f64| PI * mu * g_exact(mu).unwrap() / (2.0 * (4.0 * (1.0 + mu * mu).sqrt()).ln());
        let t3 = t(1e3);
        check(
            "[S3c] 大 μ: πμg/(2ln(4√(1+μ²))) → 1 (log 補正つき 1/μ)",
            (t3 - 1.0).abs() < 1e-5,
            format!("t(1e3) = 1{:+.1e}", t3 - 1.0),
        );
        let rho = (1e4 * g_exact(1e4).unwrap()) / (1e3 * g_exact(1e3).unwrap());
        let pred = (4.0 * (1.0f64 + 1e8).sqrt()).ln() / (4.0 * (1.0f64 + 1e6).sqrt()).ln();
        check(
            "[S3d] μg の log 成長 (純 c/μ なら 1.000、c/μ² なら 0.100)",
            (rho - pred).abs() < 1e-4,
            format!("μg 比 (1e4/1e3) = {:.5} (予言 ln 比 = {:.5})", rho, pred),
        );
    }

    // ---- [S4][S5] 事前登録判定: v25.1 一次データ (μ ≥ 0.5) との照合 ----
    let table = load_v251_table();
    match &table {
        None => {
            check(
                "[S4] 事前登録判定 (一次データ読込)",
                false,
                "results/v251_gmu.json が読めない (fail-closed)".into(),
            );
        }
        Some(tab) => {
            let sub: Vec<_> = tab.iter().filter(|(mu, _, _)| *mu >= 0.5 - 1e-12).collect();
            let mut max_g = 0.0f64;
            let mut max_gm = 0.0f64;
            let mut max_g_all = 0.0f64;
            for &&(mu, gh, gmh) in &sub {
                let ge = g_exact(mu).unwrap();
                max_g = max_g.max((gh / ge - 1.0).abs());
                max_gm = max_gm.max((gmh / ge - 1.0).abs());
            }
            for &(mu, gh, _) in tab.iter().filter(|(mu, _, _)| *mu > 0.0) {
                max_g_all = max_g_all.max((gh / g_exact(mu).unwrap() - 1.0).abs());
            }
            check(
                "[S4] 事前登録判定: max|ĝ/g_exact − 1| < 5e-5 (μ ≥ 0.5)",
                max_g < 5e-5,
                format!(
                    "max = {:.2e} ({} 点; 全域 max = {:.1e} — 小 μ crossover は登録済みの除外域)",
                    max_g,
                    sub.len(),
                    max_g_all
                ),
            );
            check(
                "[S5] g_m ≡ g (定理) の実測照合: max|ĝ_m/g_exact − 1| < 5e-5 (μ ≥ 0.5)",
                max_gm < 5e-5,
                format!("max = {:.2e} — 運動項と質量項の規格化は同一関数", max_gm),
            );
        }
    }

    // ---- [S6] 性質 ----
    {
        let g0 = g_exact(0.0).unwrap();
        check(
            "[S6a] g(0) = 1 (bit 一致)",
            g0 == 1.0,
            format!("g(0) = {}", g0),
        );
        let par = (0.5f64, 2.0f64.sqrt());
        let ok_par = g_exact(-par.0) == g_exact(par.0) && g_exact(-par.1) == g_exact(par.1);
        check("[S6b] 偶関数 g(−μ) = g(μ)", ok_par, String::new());
        let mut ok_mono = true;
        let mut ok_range = true;
        let mut prev = f64::INFINITY;
        let mut mu = 0.0;
        while mu <= 10.0001 {
            let g = g_exact(mu).unwrap();
            if g >= prev {
                ok_mono = false;
            }
            if !(0.0 < g && g <= 1.0) {
                ok_range = false;
            }
            prev = g;
            mu += 0.01;
        }
        check("[S6c] 0 < g ≤ 1 (μ ∈ [0,10] 格子)", ok_range, String::new());
        check("[S6d] 狭義単調減少 (同格子)", ok_mono, String::new());
        // g′(μ) = −μ·⟨cos²θ(1+μ²cos²θ)^{−3/2}⟩_circle < 0 と数値微分の照合
        let mut worst = 0.0f64;
        for &mu in &[0.5, 1.0, 2.0f64.sqrt()] {
            let m2 = mu * mu;
            let gp_int = -mu
                * circle_mean(
                    |t| {
                        let c2 = t.cos() * t.cos();
                        c2 / (1.0 + m2 * c2).powf(1.5)
                    },
                    1e-15,
                )
                .0;
            let h = 1e-5;
            let gp_num = (g_exact(mu + h).unwrap() - g_exact(mu - h).unwrap()) / (2.0 * h);
            worst = worst.max((gp_int / gp_num - 1.0).abs());
            if gp_int >= 0.0 {
                worst = f64::INFINITY;
            }
        }
        check(
            "[S6e] g′ < 0 (積分表示) と数値微分の一致",
            worst < 1e-7,
            format!("max 相対差 = {:.1e}", worst),
        );
    }

    // ---- [S7] 変異検出層 (破壊層): 誤った式は必ずゲートに掛かるか ----
    {
        let probe = [0.5, 1.0, 2.0f64.sqrt()];
        let dev = |f: &dyn Fn(f64) -> f64| -> f64 {
            probe
                .iter()
                .map(|&mu| (f(mu) / g_exact(mu).unwrap() - 1.0).abs())
                .fold(0.0, f64::max)
        };
        let m1 = dev(&|mu| {
            // κ↔κ′ 交換 (μ=1 の自己双対点では一致するため 3 点で判定)
            let kap = 1.0 / (1.0 + mu * mu).sqrt();
            (2.0 / PI) * (mu * kap) * ellip_k(kap).unwrap()
        });
        let m2 = dev(&|mu| 2.0 * g_exact(mu).unwrap()); // 規格 4/π
        let m3 = dev(&|mu| {
            // モジュラス/母数の混同 K(κ′²)
            let kap = 1.0 / (1.0 + mu * mu).sqrt();
            (2.0 / PI) * kap * ellip_k(mu * mu / (1.0 + mu * mu)).unwrap()
        });
        let m4 = dev(&|mu| 1.0 / agm(1.0, 1.0 + mu * mu).unwrap()); // √ 忘れ
        for (name, d) in [
            ("m1 κ↔κ′", m1),
            ("m2 規格 4/π", m2),
            ("m3 母数混同", m3),
            ("m4 √忘れ", m4),
        ] {
            check(
                &format!("[S7-{}] 変異検出", name),
                d > 1e-3,
                format!("逸脱 {:.1e} > 1e-3 (S4 ゲートが検出)", d),
            );
        }
        // 旧最良候補 π/(2K(κ′)) — 設計ノートの棄却 (0.16% → 3.8%) の再現
        let d05 = (g_old_candidate(0.5) / g_exact(0.5).unwrap() - 1.0).abs();
        let d14 = (g_old_candidate(2.0f64.sqrt()) / g_exact(2.0f64.sqrt()).unwrap() - 1.0).abs();
        check(
            "[S7-m5] 旧最良候補 π/(2K(κ′)) の棄却再現 (固有値間隔 ≠ prefactor)",
            (1e-3..3e-3).contains(&d05) && (0.025..0.055).contains(&d14),
            format!(
                "ドリフト {:.2}% @μ=0.5 → {:.2}% @√2 (設計ノート: 0.16% → 3.8%)",
                100.0 * d05,
                100.0 * d14
            ),
        );
    }

    // ---- 参照値の印字 (文書・後続コミットの一次ソース) ----
    println!("\n    [厳密値 (15 桁)] μ | κ′ | K(κ′) | g(μ) = g_m(μ)");
    for &mu in &[0.5, 1.0, 2.0f64.sqrt()] {
        let kap = 1.0 / (1.0 + mu * mu).sqrt();
        let kp = mu * kap;
        println!(
            "      μ = {:<8.6}  κ′ = {:.15}  K = {:.15}  g = {:.15}",
            mu,
            kp,
            ellip_k(kp).unwrap(),
            g_exact(mu).unwrap()
        );
    }
    println!("    [漸近] 小 μ: 1 − μ²/4 + 9μ⁴/64 − 25μ⁶/256 + 1225μ⁸/16384 − …");
    println!(
        "           大 μ: (2/(πμ))·ln(4μ)·(1 + O(1/μ²)) — v25.1 の 1/μ² 候補・μ²lnμ 候補は棄却\n"
    );

    // ---- 判定 ----
    let nfail: usize = stage_fail.values().sum();
    let s1_ok = stage_fail.get("S1").copied().unwrap_or(0) == 0;
    let s45_only = nfail > 0 && stage_fail.keys().all(|s| s == "S4" || s == "S5");
    println!(
        "[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **g(μ) = g_m(μ) = (2/π)κK(κ′) = 1/AGM(1,√(1+μ²)) と確定** — 閉形式は既知 (Eisler 2024)、本版の寄与は同定と旧仮説の棄却。λ の解析化は BZ moment (後続) へ".to_string()
        } else if s45_only {
            "事前登録 (b): 閉形式は自己整合だが v25.1 実測と不一致 — 有限 N 系統の再検討へ (記録)"
                .to_string()
        } else if !s1_ok {
            "事前登録 (c): 三表示が割れた — 実装/同定の誤り (器械再監査)".to_string()
        } else {
            format!("部分 FAIL {:?} — 個別に記録", stage_fail)
        }
    );

    let npass = checks.iter().filter(|(_, ok)| *ok).count();
    println!(
        "\n総合判定: {} ({} PASS / {} FAIL)",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" },
        npass,
        nfail
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
