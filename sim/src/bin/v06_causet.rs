//! v0.6a 因果集合 — 「時空 = 因果順序 + 数え上げ」の検証
//! d 次元ミンコフスキー時空の因果ダイヤモンドに一様ランダムに N 点を撒き (sprinkling)、
//! 座標を全て忘れて「どの点がどの点の因果的未来か」という半順序だけを残す。
//! その順序関係の統計 (順序比 r = 関係数/全ペア数) から Myrheim–Meyer 推定で
//! 時空の次元を復元する:  f(d) = Γ(d+1)Γ(d/2) / (4Γ(3d/2)) = r
//!
//! ポイント: 格子と違い、ポアソン sprinkling はローレンツブースト不変
//! (どの慣性系から見ても一様) — 「離散だがローレンツ対称」は可能である。

use uft_sim::*;

fn myrheim_meyer_f(d: f64) -> f64 {
    (ln_gamma(d + 1.0) + ln_gamma(d / 2.0) - ln_gamma(3.0 * d / 2.0)).exp() / 4.0
}

/// 順序比 r から次元を二分法で解く
fn dim_from_r(r: f64) -> f64 {
    let (mut lo, mut hi) = (1.0f64, 8.0f64);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if myrheim_meyer_f(mid) > r {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

/// d 次元因果ダイヤモンド (t∈[-1,1], |x|≤1-|t|) に N 点を撒く
fn sprinkle(n: usize, d: usize, rng: &mut Rng) -> Vec<Vec<f64>> {
    let mut pts = Vec::with_capacity(n);
    while pts.len() < n {
        let t = 2.0 * rng.f64() - 1.0;
        let x: Vec<f64> = (0..d - 1).map(|_| 2.0 * rng.f64() - 1.0).collect();
        let r2: f64 = x.iter().map(|a| a * a).sum();
        if r2.sqrt() <= 1.0 - t.abs() {
            let mut p = vec![t];
            p.extend(x);
            pts.push(p);
        }
    }
    pts
}

/// 因果関係の数 (x ≺ y: Δt ≥ |Δx|)
fn count_relations(pts: &[Vec<f64>]) -> u64 {
    let n = pts.len();
    let mut cnt = 0u64;
    for i in 0..n {
        for j in 0..n {
            if i == j {
                continue;
            }
            let dt = pts[j][0] - pts[i][0];
            if dt <= 0.0 {
                continue;
            }
            let dr2: f64 = (1..pts[i].len())
                .map(|k| (pts[j][k] - pts[i][k]).powi(2))
                .sum();
            if dt * dt >= dr2 {
                cnt += 1;
            }
        }
    }
    cnt
}

fn main() {
    let mut rng = Rng::new(20260703);
    println!("=== v0.6a 因果集合: 因果順序だけから時空次元を復元する ===\n");
    println!("[A] Myrheim–Meyer 次元推定 (N=2000, 5回平均)");
    println!("  真の次元   推定次元(±std)   順序比 r     理論値 f(d)");
    for &d in &[2usize, 3, 4] {
        let mut dims = Vec::new();
        for _ in 0..5 {
            let pts = sprinkle(2000, d, &mut rng);
            let rel = count_relations(&pts);
            let r = rel as f64 / (2000.0 * 1999.0); // 順序対 N(N-1) で規格化
            dims.push(dim_from_r(r));
        }
        let mean: f64 = dims.iter().sum::<f64>() / dims.len() as f64;
        let sd: f64 =
            (dims.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / (dims.len() - 1) as f64).sqrt();
        let pts = sprinkle(2000, d, &mut rng);
        let r = count_relations(&pts) as f64 / (2000.0 * 1999.0);
        println!(
            "     {}       {:.3} ± {:.3}     {:.4}      f({})={:.4}",
            d,
            mean,
            sd,
            r,
            d,
            myrheim_meyer_f(d as f64)
        );
    }
    println!("  => 座標も計量も使わず、「AはBの過去か?」の表だけで次元が測れる。");
    println!(
        "     時空の幾何情報は (因果順序, 要素数) に完全に含まれる (Malament の定理 + 体積=数)。\n"
    );

    println!("[B] ローレンツ不変性: sprinkling をブーストしても統計は不変か (d=2, v=0.6c)");
    {
        let n = 3000;
        let pts = sprinkle(n, 2, &mut rng);
        let r0 = count_relations(&pts) as f64 / (n as f64 * (n - 1) as f64);
        // ブースト後、同じダイヤモンドに入る点だけで再計測するため、
        // ブースト不変な部分ダイヤモンド (中心の小ダイヤモンド) 内の点で比較
        let v = 0.6;
        let g = 1.0 / (1.0 - v * v as f64).sqrt();
        let boosted: Vec<Vec<f64>> = pts
            .iter()
            .map(|p| vec![g * (p[0] - v * p[1]), g * (p[1] - v * p[0])])
            .collect();
        let sub = |ps: &[Vec<f64>]| -> Vec<Vec<f64>> {
            ps.iter()
                .filter(|p| p[0].abs() + p[1].abs() <= 0.35)
                .cloned()
                .collect()
        };
        let s0 = sub(&pts);
        let s1 = sub(&boosted);
        let rr0 = count_relations(&s0) as f64 / (s0.len() as f64 * (s0.len() - 1) as f64);
        let rr1 = count_relations(&s1) as f64 / (s1.len() as f64 * (s1.len() - 1) as f64);
        println!(
            "  静止系の小ダイヤモンド: N={}, r={:.4} → d={:.2}",
            s0.len(),
            rr0,
            dim_from_r(rr0)
        );
        println!(
            "  ブースト系の同領域   : N={}, r={:.4} → d={:.2}",
            s1.len(),
            rr1,
            dim_from_r(rr1)
        );
        println!(
            "  (因果関係そのものはブーストで厳密に不変: r_full={:.4} → {:.4})",
            r0,
            count_relations(&boosted) as f64 / (n as f64 * (n - 1) as f64)
        );
        println!("  => 正方格子なら特定の慣性系が最初から刻印されるが、ポアソン撒布には");
        println!("     どの方向も刻印されていない。「離散 かつ ローレンツ不変」は両立する。");
    }
    println!(
        "\n結論: 時空の候補となる最小データは (集合, 半順序)。距離・次元・体積は導出量である。"
    );
    println!(
        "      ただし因果集合だけでは量子重ね合わせがない — 量子データ(もつれ)が必要 (→v0.7)。"
    );
}
