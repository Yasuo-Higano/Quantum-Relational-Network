//! v24.4 演算子整合 — 傾き Z_T と方向・距離分解残差 (第二十五期の本丸 II)
//!
//! PROMPT/6 PR5。固定 ξ の比 λ(ξ) = πξ/A(ξ) は UV 境界層の加法的ずれ b に汚染される:
//! A(ξ) = Z·πξ + b なら λ(ξ) = 1/(Z + b/(πξ)) — ξ 有限では Z と b が縮退する。
//! 演算子整合はこれを分離する:
//!
//!   min_{Z, b} Σ_窓 |K_A のボンド成分(ξ) − Z·(2π/v)·ξ·h(ξ) − b·h(ξ)|²
//!
//! ここで v = 1 (v24.1 で厳密確定)、h は格子 H 自身のボンド演算子 (自由場では
//! T⁰⁰ の規格化は H で固定され、Z_T の独立決定に相当)。b は切断面付近に許される
//! 境界演算子 (ξ 非依存の NN 係数)。**λ_R = 1/Z が演算子整合後の温度比**。
//!
//! 測定 (認証済みブロック器械, DD + クランプ 2 点梯子の信頼域のみ):
//!   - x チャネル: A_x(ξ) を Z_x·πξ + b_x でフィット (窓 ξ ∈ [3, ξ*])
//!   - y チャネル: (−1)^x K_y(ξ_s) を Z_y·πξ_s + b_y で
//!   - z チャネル: (−1)^x B_z(ξ_s) を Z_z·πξ_s + b_z で
//!   - 頑健性: フィット族 {Zπξ+b, Zπξ+b+c/ξ}・holdout (窓上端 1 点予測)・N=32 vs 64
//!   - 残差プロファイル r(ξ) と、BW が禁止する成分の大きさ:
//!     NNN (A_x2/A_x)・on-site 一様 (on_uni の ξ 依存)
//!
//! 事前登録:
//!   (a) Z_x = Z_y = Z_z (±1%) かつ λ_R = 1/Z = 1 ± 2% = **BW は演算子整合で素朴回復**
//!       — 32/27 は「切断面直近の境界演算子 b を Z に混入させた固定 ξ 推定器の値」と確定。
//!   (a′) Z 等方だが 1/Z = 32/27 ± 1% = 有限繰り込みが整合後も生存 (重大 — v24.5 へ)。
//!   (b) Z が方向依存 (>1%) = Lorentz 等方性の破れ候補 / (c) その他 = 保留・記録。
//! 器械ゲート: [M1] 窓内の残差 |r(ξ)|/A_x(ξ) < 1% (線形性の成立)、[M2] 族間 Z 差 < 0.5%。

use uft_sim::dd::*;
use uft_sim::stag::*;
use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// 重み無し最小二乗 (2〜3 基底, 正規方程式)
fn lsq(xs: &[f64], ys: &[f64], basis: &dyn Fn(f64) -> Vec<f64>) -> Vec<f64> {
    let k = basis(xs[0]).len();
    let mut ata = vec![0.0f64; k * k];
    let mut atb = vec![0.0f64; k];
    for (&x, &y) in xs.iter().zip(ys) {
        let b = basis(x);
        for i in 0..k {
            atb[i] += b[i] * y;
            for j in 0..k {
                ata[i * k + j] += b[i] * b[j];
            }
        }
    }
    let mut m = ata;
    let mut v = atb;
    for col in 0..k {
        let mut piv = col;
        for r in col + 1..k {
            if m[r * k + col].abs() > m[piv * k + col].abs() {
                piv = r;
            }
        }
        for j in 0..k {
            m.swap(col * k + j, piv * k + j);
        }
        v.swap(col, piv);
        let d = m[col * k + col];
        for r in 0..k {
            if r == col {
                continue;
            }
            let f = m[r * k + col] / d;
            for j in 0..k {
                m[r * k + j] -= f * m[col * k + j];
            }
            v[r] -= f * v[col];
        }
    }
    (0..k).map(|i| v[i] / m[i * k + i]).collect()
}

struct Chan {
    name: &'static str,
    xs: Vec<f64>, // ξ (信頼域)
    ys: Vec<f64>, // 測定値 (符号整流済み)
}

fn main() {
    self_test();
    println!("=== v24.4 演算子整合 — 傾き Z_T と方向・距離分解残差 (第二十五期) ===\n");
    println!("事前登録: (a) Z 等方 & 1/Z = 1±2% → BW 素朴回復 (32/27 は境界演算子の混入) /");
    println!("          (a′) Z 等方 & 1/Z = 32/27±1% → 繰り込み生存 / (b) Z 方向依存 / (c) 保留\n");
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(1);
    check("[M0a] dd 自己検証", dd_self_test(), String::new());
    check("[M0b] stag 自己検証", stag_self_test(), String::new());

    let mut all_z: Vec<(usize, &'static str, f64, f64, f64)> = Vec::new(); // (N, ch, Z, b, Z_alt)
    let mut resid_max_rel = 0.0f64;
    let mut fam_dev_max = 0.0f64;

    for &n in &[32usize, 64] {
        let half = n / 2;
        let s30 = half_space_scan::<Dd>(n, 1e-30, nthreads);
        let s26 = half_space_scan::<Dd>(n, 1e-26, nthreads);
        // 信頼域 (x チャネルの梯子)
        let mut xi_trust = 0usize;
        for xi in 1..half {
            let i = half - 1 - xi;
            if ((s30.ax[i] - s26.ax[i]) / s30.ax[i]).abs() < 1e-4 {
                xi_trust = xi;
            } else {
                break;
            }
        }
        println!(
            "\n    N={}: 信頼域 ξ ≤ {} (κ_max = {:.1}, {} s)",
            n,
            xi_trust,
            s30.kappa_max,
            t0.elapsed().as_secs()
        );
        // チャネル構築 (符号整流): x はボンド中心 ξ、y/z はサイト ξ_s
        let mut chans: Vec<Chan> = Vec::new();
        {
            let mut xs = Vec::new();
            let mut ys = Vec::new();
            for xi in 3..=xi_trust {
                let i = half - 1 - xi;
                xs.push(xi as f64);
                ys.push(s30.ax[i]);
            }
            chans.push(Chan { name: "x", xs, ys });
        }
        for (name, arr) in [("y", &s30.ky), ("z", &s30.bz)] {
            let mut xs = Vec::new();
            let mut ys = Vec::new();
            for j in 0..half {
                let i = half - 1 - j;
                let xis = j as f64 + 0.5;
                if xis < 2.5 || xis > xi_trust as f64 {
                    continue;
                }
                let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
                xs.push(xis);
                ys.push(sign * arr[i]);
            }
            chans.push(Chan { name, xs, ys });
        }
        for ch in &chans {
            // 族 1: Z·πξ + b
            let c1 = lsq(&ch.xs, &ch.ys, &|x| vec![PI * x, 1.0]);
            // 族 2: Z·πξ + b + c/ξ
            let c2 = lsq(&ch.xs, &ch.ys, &|x| vec![PI * x, 1.0, 1.0 / x]);
            // 残差 (族 1)
            let mut rmax = 0.0f64;
            for (&x, &y) in ch.xs.iter().zip(&ch.ys) {
                let pred = c1[0] * PI * x + c1[1];
                rmax = rmax.max(((y - pred) / y).abs());
            }
            // holdout: 上端 1 点
            let nfit = ch.xs.len() - 1;
            let ch1 = lsq(&ch.xs[..nfit], &ch.ys[..nfit], &|x| vec![PI * x, 1.0]);
            let pred_h = ch1[0] * PI * ch.xs[nfit] + ch1[1];
            let herr = ((pred_h - ch.ys[nfit]) / ch.ys[nfit]).abs();
            println!(
                "      ch={} : Z = {:.5} (b = {:+.4}, 族2 Z = {:.5}, holdout {:.1e}), 1/Z = {:.5}, 残差max {:.2e}",
                ch.name,
                c1[0],
                c1[1],
                c2[0],
                herr,
                1.0 / c1[0],
                rmax
            );
            resid_max_rel = resid_max_rel.max(rmax);
            fam_dev_max = fam_dev_max.max(((c2[0] - c1[0]) / c1[0]).abs());
            all_z.push((n, ch.name, c1[0], c1[1], c2[0]));
        }
        // BW が禁止する成分の大きさ (窓内)
        let mut nnn_rel = 0.0f64;
        for xi in 3..=xi_trust.min(half.saturating_sub(3)) {
            let i = half - 1 - xi;
            if i + 2 < half && i < s30.ax2.len() {
                nnn_rel = nnn_rel.max((s30.ax2[i] / s30.ax[i]).abs());
            }
        }
        let on_drift = {
            // on-site 一様成分の ξ 依存 (窓の端点差)
            let i1 = half - 1 - 3;
            let i2 = half - 1 - xi_trust;
            (s30.on_uni[i2] - s30.on_uni[i1]).abs()
        };
        println!(
            "      禁止成分: NNN/NN = {:.2e} (窓内 max), on-site 一様のドリフト = {:.2e}",
            nnn_rel, on_drift
        );
    }

    // ---- ゲートと判定 ----
    check(
        "[M1] 窓内残差 < 1% (BW 線形形の成立)",
        resid_max_rel < 0.01,
        format!("max |r|/A = {:.2e}", resid_max_rel),
    );
    check(
        "[M2] フィット族間の Z 安定性 < 0.5%",
        fam_dev_max < 0.005,
        format!("max 族間差 = {:.2e}", fam_dev_max),
    );
    // N=64 の 3 チャネルで判定
    let z64: Vec<(&str, f64)> = all_z
        .iter()
        .filter(|r| r.0 == 64)
        .map(|r| (r.1, r.2))
        .collect();
    let zx = z64.iter().find(|r| r.0 == "x").unwrap().1;
    let zy = z64.iter().find(|r| r.0 == "y").unwrap().1;
    let zz = z64.iter().find(|r| r.0 == "z").unwrap().1;
    let iso_dev = ((zy / zx - 1.0).abs()).max((zz / zx - 1.0).abs());
    let lam_r = 1.0 / zx;
    let target = 32.0 / 27.0;
    println!(
        "\n    [N=64 整合結果] Z_x = {:.5}, Z_y = {:.5}, Z_z = {:.5} (等方偏差 {:.2e})",
        zx, zy, zz, iso_dev
    );
    println!(
        "    λ_R = 1/Z_x = {:.5}  (1 との差 {:+.4}, 32/27 との差 {:+.4})",
        lam_r,
        lam_r - 1.0,
        lam_r - target
    );
    let branch_a = iso_dev < 0.01 && (lam_r - 1.0).abs() < 0.02;
    let branch_ap = iso_dev < 0.01 && (lam_r - target).abs() < 0.01 * target;
    let branch_b = iso_dev >= 0.01;
    println!(
        "\n[判定] {}",
        if nfail > 0 {
            "装置ゲート故障 — 記録".to_string()
        } else if branch_a {
            format!(
                "事前登録 (a): Z 等方・λ_R = {:.4} → 1 — BW は演算子整合で素朴回復。32/27 は境界演算子 b の混入と確定",
                lam_r
            )
        } else if branch_ap {
            format!("事前登録 (a′): 整合後も λ_R = {:.5} ≈ 32/27 — 繰り込み生存 (重大)", lam_r)
        } else if branch_b {
            format!("事前登録 (b): Z の方向依存 {:.3} — 等方性の破れ候補", iso_dev)
        } else {
            format!("事前登録 (c): 保留 — λ_R = {:.5}, 全表を一次記録", lam_r)
        }
    );

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v24.4".into())),
        (
            "z_rows".into(),
            Json::Arr(
                all_z
                    .iter()
                    .map(|&(n, ch, z, b, z2)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("ch".into(), Json::Str(ch.into())),
                            ("z".into(), Json::Num(z)),
                            ("b".into(), Json::Num(b)),
                            ("z_fam2".into(), Json::Num(z2)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("lam_r".into(), Json::Num(lam_r)),
        ("iso_dev".into(), Json::Num(iso_dev)),
        ("branch_a".into(), Json::Bool(branch_a)),
        ("branch_ap".into(), Json::Bool(branch_ap)),
    ]);
    let p = write_artifact("results/v244_opmatch.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(a′)/(b)/(c) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
