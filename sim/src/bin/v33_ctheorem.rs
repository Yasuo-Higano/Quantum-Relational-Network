//! v3.3 RG の不可逆性 = 情報の単調性 — エントロピック c 定理 (厳密)
//!
//! くりこみ群 (v0.3) は「ミクロの詳細が洗い流される」一方向の流れに見える。
//! Zamolodchikov (1986): 2d QFT には RG 流に沿って単調減少する関数 c が存在する。
//! Casini–Huerta (2004) はこれをエンタングルメントで再導出した:
//!   c(ℓ) ≡ 3 ℓ dS/dℓ は強劣加法性 (SSA) とローレンツ対称性だけから単調減少。
//! *** RG の不可逆性は、量子情報の単調性 (データ処理不等式) の帰結である ***
//!
//! 検証: 交替質量 m のフェルミオン鎖 (UV: c=1 のディラック CFT → IR: ギャップ相 c=0)。
//! c(ℓ) が 1 → 0 へ単調に流れ、異なる m の曲線が x = ℓ·m でスケーリング収束するか。

use uft_sim::*;

fn h2(z: f64) -> f64 {
    let z = z.clamp(1e-14, 1.0 - 1e-14);
    -z * z.ln() - (1.0 - z) * (1.0 - z).ln()
}
fn entropy_real(c: &[f64], n: usize) -> f64 {
    let (w, _) = jacobi_eigh(c, n);
    w.iter().map(|&z| h2(z)).sum()
}

fn main() {
    let n = 400usize;
    println!(
        "=== v3.3 エントロピック c 定理: RG 流 c: 1 → 0 (N={} 鎖, 厳密) ===\n",
        n
    );
    println!("模型: H = -Σ(c†c+h.c.) + m Σ(-1)^x c†c  (m=0 で c=1 の CFT, m>0 でギャップ)");
    println!("c関数: c(ℓ) = 3·[S(ℓ+2)-S(ℓ-2)] / [ln chord(ℓ+2) - ln chord(ℓ-2)]\n");
    let chord =
        |l: f64| (n as f64 / std::f64::consts::PI) * (std::f64::consts::PI * l / n as f64).sin();
    let mut all_curves: Vec<(f64, Vec<(f64, f64)>)> = Vec::new();
    for &m in &[0.0f64, 0.05, 0.15, 0.5] {
        // 単一粒子ハミルトニアン
        let mut a = vec![0.0; n * n];
        for x in 0..n {
            let y = (x + 1) % n;
            a[x + y * n] = -1.0;
            a[y + x * n] = -1.0;
            a[x + x * n] = m * if x % 2 == 0 { 1.0 } else { -1.0 };
        }
        let (_, v) = jacobi_eigh(&a, n);
        let mut cf = vec![0.0; n * n];
        for k in 0..n / 2 {
            for i in 0..n {
                let vi = v[i + k * n];
                if vi == 0.0 {
                    continue;
                }
                for j in 0..n {
                    cf[i + j * n] += vi * v[j + k * n];
                }
            }
        }
        // S(ℓ) を偶数 ℓ で計算
        let ells: Vec<usize> = (4..=124).step_by(4).collect();
        let mut svals = Vec::new();
        for &l in &ells {
            let mut ca = vec![0.0; l * l];
            for i in 0..l {
                for j in 0..l {
                    ca[i + j * l] = cf[(100 + i) + (100 + j) * n];
                }
            }
            svals.push(entropy_real(&ca, l));
        }
        // c(ℓ) 中心差分
        let mut curve = Vec::new();
        for w2 in 1..ells.len() - 1 {
            let l = ells[w2] as f64;
            let ds = svals[w2 + 1] - svals[w2 - 1];
            let dln = (chord(ells[w2 + 1] as f64) / chord(ells[w2 - 1] as f64)).ln();
            curve.push((l, 3.0 * ds / dln));
        }
        all_curves.push((m, curve));
    }
    println!("  ℓ      c(m=0)   c(m=0.05)  c(m=0.15)  c(m=0.5)");
    let npts = all_curves[0].1.len();
    for i in 0..npts {
        let l = all_curves[0].1[i].0;
        print!("  {:4.0}  ", l);
        for (_, curve) in &all_curves {
            print!("  {:7.4} ", curve[i].1);
        }
        println!();
    }
    // 検定 1: m=0 は c ≈ 1 で一定
    let c0: Vec<f64> = all_curves[0].1.iter().map(|p| p.1).collect();
    let c0max = c0.iter().cloned().fold(0.0f64, f64::max);
    let c0min = c0.iter().cloned().fold(2.0f64, f64::min);
    println!(
        "\n  [1] UV 固定点 (m=0): c = {:.3}..{:.3} — 中心電荷 1 のプラトー  {}",
        c0min,
        c0max,
        pass((c0max - 1.0).abs() < 0.05 && (c0min - 1.0).abs() < 0.05)
    );
    // 検定 2: m>0 で単調非増加
    let mut mono = true;
    for (m, curve) in all_curves.iter().skip(1) {
        for w2 in 1..curve.len() {
            if curve[w2].1 > curve[w2 - 1].1 + 1e-3 {
                mono = false;
                println!("  違反: m={} ℓ={}", m, curve[w2].0);
            }
        }
    }
    println!("  [2] RG 単調性 (m>0 の全曲線が非増加): {}", pass(mono));
    // 検定 3: スケーリング収束 c(ℓ, m) = f(ℓ·m)
    println!("  [3] スケーリング: c を x=ℓ·m で比較 (異なる m が同一曲線に乗るか)");
    println!("      x=ℓm    c(m=0.05)   c(m=0.15)");
    let interp = |curve: &Vec<(f64, f64)>, x: f64, m: f64| -> Option<f64> {
        let l = x / m;
        for w2 in 1..curve.len() {
            if curve[w2 - 1].0 <= l && l <= curve[w2].0 {
                let t = (l - curve[w2 - 1].0) / (curve[w2].0 - curve[w2 - 1].0);
                return Some(curve[w2 - 1].1 * (1.0 - t) + curve[w2].1 * t);
            }
        }
        None
    };
    let mut maxdiff = 0.0f64;
    for &x in &[1.0f64, 2.0, 3.0, 4.0, 6.0] {
        let c1 = interp(&all_curves[1].1, x, 0.05);
        let c2 = interp(&all_curves[2].1, x, 0.15);
        if let (Some(c1), Some(c2)) = (c1, c2) {
            maxdiff = maxdiff.max((c1 - c2).abs());
            println!("      {:4.1}    {:.4}      {:.4}", x, c1, c2);
        }
    }
    println!(
        "      => 最大差 {:.3} (普遍スケーリング関数に収束)  {}",
        maxdiff,
        pass(maxdiff < 0.06)
    );
    println!("\n結論: RG 流は c(UV)=1 → c(IR)=0 へ単調に流れ、途中経過は普遍関数 f(ℓm) に乗る。");
    println!("      Casini–Huerta の証明はこの単調性が強劣加法性 (情報の一般定理) の帰結である");
    println!("      ことを示す。時間の矢 (v1.1)・データ処理 (v2.1)・RG (v3.3) — 物理の 3 つの");
    println!("      「不可逆性」は、すべて同じ情報単調性の異なる顔である。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
