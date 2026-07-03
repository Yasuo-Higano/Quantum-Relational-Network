//! v5.3 なぜ低エネルギーには「小さな理論」が見えるのか — 脱結合定理の定量検証
//!
//! v4.3 の残された恣意性:「なぜ ≤24 成分規模の理論を考えれば良いのか」。
//! 答えは 2 つの定理の合わせ技:
//!  (1) 脱結合 (Appelquist–Carazzone): 質量 M の場は低エネルギー観測量に 1/M² でしか効かない
//!  (2) カイラル保護: 質量項を書けない (=カイラルな) 場だけが低エネルギーに残る
//!      → 残る内容は v4.3/v5.2 の全数探索が「SM 一意」と決めた
//! ここでは (1) を数値で: 軽い鎖 a + 重い鎖 b (質量 M, 混合 g) の基底状態で、
//! 軽い場の相関関数への重い場の痕跡が M^{-2} で消えることをフィットで確認する。

use uft_sim::*;

fn light_corr(m_heavy: f64, g: f64, d: usize) -> f64 {
    // 2 バンド: h(k) = [[ε(k), g],[g, ε(k)+M]], ε = -2cos k。E<0 の状態を占有。
    // C^aa(d) = (1/N) Σ_k Σ_{band: E<0} |u_a|² cos(kd)
    let n = 4000usize;
    let mut c = 0.0;
    for i in 0..n {
        let k = 2.0 * std::f64::consts::PI * (i as f64 + 0.5) / n as f64;
        let e = -2.0 * k.cos();
        // 2×2 対角化
        let (h11, h22) = (e, e + m_heavy);
        let tr = 0.5 * (h11 + h22);
        let dd = 0.5 * (h11 - h22);
        let r = (dd * dd + g * g).sqrt();
        let (e1, e2) = (tr - r, tr + r);
        // 固有ベクトル: 下バンドの a 成分 = (dd-r)/denom, 上バンドの a 成分 = g/denom
        let denom = (dd - r).hypot(g);
        let (u1a, u2a) = if denom > 1e-15 {
            ((dd - r) / denom, g / denom)
        } else {
            (1.0, 0.0)
        };
        if e1 < 0.0 {
            c += u1a * u1a * (k * d as f64).cos();
        }
        if e2 < 0.0 {
            c += u2a * u2a * (k * d as f64).cos();
        }
    }
    c / n as f64
}

fn main() {
    println!("=== v5.3 脱結合定理: 重い場の痕跡は 1/M² で消える ===\n");
    let g = 0.5f64;
    let d = 5usize;
    let c_inf = light_corr(1e9, g, d); // M→∞ の基準 (事実上純粋な軽鎖)
    println!("軽い場の相関 C_aa(d=5)。混合 g={}。M→∞ 基準値 = {:.6}\n", g, c_inf);
    println!("  M      C_aa(M)     |ΔC| = |C(M)-C(∞)|");
    let mut lnm = Vec::new();
    let mut lnd = Vec::new();
    for &m in &[4.0f64, 8.0, 16.0, 32.0, 64.0, 128.0] {
        let c = light_corr(m, g, d);
        let dc = (c - c_inf).abs();
        println!("  {:5.0}  {:+.6}   {:.3e}", m, c, dc);
        lnm.push(m.ln());
        lnd.push(dc.ln());
    }
    let (_, slope) = linfit(&lnm, &lnd);
    println!("\n  => |ΔC| ∝ M^({:.3})   (脱結合定理の予言: -2)  {}", slope, pass((slope + 2.0).abs() < 0.15));
    println!("\n結論: 重い場は低エネルギーから M^-2 で退場する (数値実証)。よって:");
    println!("      - 低エネルギーに見える物質 = 質量項を持てないカイラルな場だけ");
    println!("      - その内容は無矛盾性で一意に決まる (v4.3/v5.2: SM 世代のみ)");
    println!("      - UV に何があろうと (弦・網・未知の巨大構造)、我々に見える理論は小さい");
    println!("      「なぜ小さな理論か」は謎ではなく、EFT + カイラル保護 + アノマリー一意性の定理連鎖。");
    println!("      裏返せば: プランクスケールの直接観測が絶望的に難しい理由でもある (正直に)。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
