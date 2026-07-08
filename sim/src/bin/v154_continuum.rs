//! v15.4 連続極限の監査 — residual(N) = A·N^(−p) + B を readout 群に課す
//!
//! ASM-LATTICE (格子正則化) は依存グラフ上で最大の影響範囲 (51 主張) を持つのに、
//! 系統的な連続極限監査を受けていなかった。次の方針 (PROMPT/2 柱 3) に従い、
//! 「ある N で PASS」を「N → ∞ で何が残るか」に置き換える:
//!
//!     observable(N) = B + A·N^(−p)   を fit し、
//!     B が連続極限の値 (既知なら比較、未知なら外挿値 ± 収束証拠) になる。
//!     消えるべき残差の B ≠ 0 が確立すれば FAL-CONTINUUM が発火する。
//!
//! 構成 (全て [PASS]/[FAIL] 内蔵):
//!   [0] fit 装置の自己検証 — 合成データ (B=0 / B≠0 / p の回復)。装置を信じない。
//!   [1] 厳密自由系列での較正 — Z2 core の h=0 厳密級数 e(L) (磁束 2 セクター min)
//!       の連続極限は解析値 −2w/π。fit がそれを復元するか (Casimir 補正 p≈2 込み)。
//!   [2] Z2GaugeRing (拘束つき相互作用 core, v15.3) の L 掃引 L=6..16:
//!       基底エネルギー密度 e(L)・ゲージ誘起凝縮 χ(L)・ギャップ Δ(L)・
//!       半分割エントロピー S(L)・弦張力 σ(L, 重物質)・幾何読み出し (全 L で安定)。
//!       各 L で h=0 厳密照合を先に走らせる (掃引自体の較正)。
//!       収束判定は fit と独立に Cauchy 差 |f(L)−f(L−2)| の単調減少でも行う。
//!   [3] 前線速度 v(L): 推定器の L 依存ドリフトの定量化 + 閉じ込め抑制の L 安定性。
//!   [4] RingChain (ガウス core, v6.7) の中心電荷 c(N), N=66..402:
//!       連続極限の既知値 c=1 に対する B=0 型の成功例 (肯定対照)。
//!
//! 出力: readout × (A, p, B, 収束判定) の「連続極限監査地図」+ JSON artifact。
//! 陰性結果も価値を持つ: 前線推定器のドリフトは推定器の限界として記録される。

use std::f64::consts::PI;
use uft_sim::*;

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// 3 パラメータ fit y = b + a·x^(−p): p をグリッド走査し (a, b) は閉形式最小二乗。
/// 戻り値 (a, p, b, rms)。
fn fit3(xs: &[f64], ys: &[f64]) -> (f64, f64, f64, f64) {
    let n = xs.len() as f64;
    let mut best = (0.0, 0.0, 0.0, f64::INFINITY);
    let mut p = 0.2;
    while p <= 4.0 + 1e-9 {
        let us: Vec<f64> = xs.iter().map(|&x| x.powf(-p)).collect();
        let su: f64 = us.iter().sum();
        let suu: f64 = us.iter().map(|u| u * u).sum();
        let sy: f64 = ys.iter().sum();
        let suy: f64 = us.iter().zip(ys).map(|(u, y)| u * y).sum();
        let det = n * suu - su * su;
        if det.abs() > 1e-14 {
            let a = (n * suy - su * sy) / det;
            let b = (sy - a * su) / n;
            let rms = (xs
                .iter()
                .zip(ys)
                .map(|(&x, &y)| (y - b - a * x.powf(-p)).powi(2))
                .sum::<f64>()
                / n)
                .sqrt();
            if rms < best.3 {
                best = (a, p, b, rms);
            }
        }
        p += 0.02;
    }
    best
}

/// h=0 の厳密自由エネルギー (磁束 2 セクターの最小) — v15.3 と同じ規約
fn free_gs_energy_min(l: usize, nf: usize, w: f64) -> f64 {
    let mut emin = f64::INFINITY;
    for t in [1.0f64, -1.0] {
        let twist = t * if (nf as i32 - 1) % 2 == 0 { 1.0 } else { -1.0 };
        let theta = if twist > 0.0 { 0.0 } else { PI };
        let mut ek: Vec<f64> = (0..l)
            .map(|n| -2.0 * w * ((2.0 * PI * n as f64 + theta) / l as f64).cos())
            .collect();
        ek.sort_by(|a, b| a.partial_cmp(b).unwrap());
        emin = emin.min(ek[..nf].iter().sum::<f64>());
    }
    emin
}

/// Cauchy 収束: 連続する差の絶対値が (緩みつきで) 減少しているか
fn cauchy_ok(ys: &[f64], slack: f64) -> bool {
    let d: Vec<f64> = ys.windows(2).map(|w| (w[1] - w[0]).abs()).collect();
    d.windows(2).all(|p| p[1] <= p[0] * (1.0 + slack) + 1e-12)
}

fn main() {
    self_test();
    println!("=== v15.4 連続極限の監査: residual(N) = A·N^(−p) + B ===\n");
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };
    let mut map_rows: Vec<(String, f64, f64, f64, String)> = Vec::new(); // (readout, A, p, B, 判定)

    // ================= [0] fit 装置の自己検証 =================
    println!("[0] fit 装置の自己検証 (合成データ — 装置を信じない)");
    {
        let xs: Vec<f64> = [6.0, 8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 26.0].to_vec();
        // (i) B=0, A=2, p=1.5
        let ys: Vec<f64> = xs.iter().map(|&x| 2.0 * x.powf(-1.5)).collect();
        let (a, p, b, _) = fit3(&xs, &ys);
        check(
            "合成 (B=0, A=2, p=1.5) の回復",
            b.abs() < 2e-3 && (a - 2.0).abs() < 0.05 && (p - 1.5).abs() < 0.05,
            format!("Â={:.3}, p̂={:.2}, B̂={:.1e}", a, p, b),
        );
        // (ii) B=0.5 (消えない残差) の検出
        let ys2: Vec<f64> = xs.iter().map(|&x| 0.5 + 1.0 * x.powf(-2.0)).collect();
        let (a2, p2, b2, _) = fit3(&xs, &ys2);
        check(
            "合成 (B=0.5, A=1, p=2) の検出 — 消えない残差を見逃さない",
            (b2 - 0.5).abs() < 2e-2 && (p2 - 2.0).abs() < 0.2 && (a2 - 1.0).abs() < 0.1,
            format!("Â={:.3}, p̂={:.2}, B̂={:.4}", a2, p2, b2),
        );
        // (iii) ノイズ耐性 (1e-4 の摂動で B の誤検出をしない)
        let mut rng = Rng::new(154);
        let ys3: Vec<f64> = xs
            .iter()
            .map(|&x| 3.0 * x.powf(-1.0) + 1e-4 * rng.gauss())
            .collect();
        let (_a3, _p3, b3, _) = fit3(&xs, &ys3);
        check(
            "合成 (B=0 + ノイズ 1e-4) で偽の B を出さない",
            b3.abs() < 5e-3,
            format!("B̂={:.1e}", b3),
        );
    }

    // ================= [1] 厳密自由系列での較正 =================
    println!("\n[1] 厳密自由系列 e(L) = E₀/L (h=0) — 連続極限の解析値 −2w/π への収束");
    {
        let ls: Vec<usize> = vec![10, 14, 18, 22, 26, 30, 38, 46];
        let xs: Vec<f64> = ls.iter().map(|&l| l as f64).collect();
        let ys: Vec<f64> = ls
            .iter()
            .map(|&l| free_gs_energy_min(l, l / 2, 1.0) / l as f64)
            .collect();
        let (a, p, b, rms) = fit3(&xs, &ys);
        let exact = -2.0 / PI;
        check(
            "e∞ の復元: B̂ = −2/π (解析値)、Casimir 型 p ≈ 2",
            (b - exact).abs() < 2e-4 && (p - 2.0).abs() < 0.25,
            format!(
                "B̂ = {:.6} (厳密 {:.6}, |Δ|={:.1e}), p̂ = {:.2}, rms = {:.1e}",
                b,
                exact,
                (b - exact).abs(),
                p,
                rms
            ),
        );
        map_rows.push(("e_free(L) [較正]".into(), a, p, b, "B=−2/π 一致".into()));
    }

    // ================= [2] Z2 core の L 掃引 =================
    let (w0, h0, m0) = (1.0, 0.6, 0.2);
    println!(
        "\n[2] Z2GaugeRing の L 掃引 (L=6..16, w={}, h={}, m={}) — 相互作用 readout の収束",
        w0, h0, m0
    );
    let ls: Vec<usize> = vec![6, 8, 10, 12, 14, 16];
    let mut e_l = Vec::new();
    let mut chi_l = Vec::new();
    let mut gap_l = Vec::new();
    let mut s_l = Vec::new();
    let mut adj_all = true;
    let mut free_check_max: f64 = 0.0;
    for &l in &ls {
        let nf = l / 2;
        // 掃引自体の較正: この L で h=0 が厳密自由と一致するか
        let gfree = Z2GaugeRing::try_new(l, nf, w0, 0.0, 0.0, vec![]).unwrap();
        let (ef, _sf, _rf) = gfree.ground_state(20260708);
        free_check_max = free_check_max.max((ef - free_gs_energy_min(l, nf, w0)).abs());
        // 相互作用点
        let g = Z2GaugeRing::try_new(l, nf, w0, h0, m0, vec![]).unwrap();
        let mv = |v: &[(f64, f64)]| g.matvec_c(v);
        let (ev, vecs, res) = lanczos_lowest_herm(&mv, g.dim(), 2, 160, 20260708);
        assert!(res < 1e-6, "Lanczos 未収束 L={}", l);
        e_l.push(ev[0] / l as f64);
        gap_l.push(ev[1] - ev[0]);
        let st = Z2CoreState {
            l,
            masks: (0u32..(1 << l))
                .filter(|x| x.count_ones() as usize == nf)
                .collect(),
            psi: vecs[0].clone(),
        };
        let chi: f64 = (0..l)
            .map(|x| st.density(x) * if x % 2 == 0 { 1.0 } else { -1.0 })
            .sum::<f64>()
            / l as f64;
        chi_l.push(chi);
        let region: Vec<usize> = (0..l / 2).collect();
        s_l.push(v2_entropy(&st, &region));
        if l >= 8 {
            let (mi, mim) = v2_mi_matrix(&st);
            let geo = readout_ring_geometry(&mi, l, mim);
            if geo.adjacency < 0.999 {
                adj_all = false;
            }
        }
    }
    check(
        "掃引の較正: 全 L で h=0 が厳密自由と一致",
        free_check_max < 1e-8,
        format!("max|Δ| = {:.1e}", free_check_max),
    );
    let xs: Vec<f64> = ls.iter().map(|&l| l as f64).collect();
    let quantities: Vec<(&str, &Vec<f64>, f64)> = vec![
        ("e(L) 基底エネルギー密度", &e_l, 0.05),
        ("χ(L) ゲージ誘起凝縮", &chi_l, 0.30),
        ("Δ(L) ギャップ", &gap_l, 0.30),
        ("S(L) 半分割エントロピー", &s_l, 0.40),
    ];
    println!("   L 掃引の生データ:");
    println!("     L    = {:?}", ls);
    println!(
        "     e    = {:?}",
        e_l.iter()
            .map(|x| (x * 1e4).round() / 1e4)
            .collect::<Vec<_>>()
    );
    println!(
        "     χ    = {:?}",
        chi_l
            .iter()
            .map(|x| (x * 1e4).round() / 1e4)
            .collect::<Vec<_>>()
    );
    println!(
        "     Δ    = {:?}",
        gap_l
            .iter()
            .map(|x| (x * 1e4).round() / 1e4)
            .collect::<Vec<_>>()
    );
    println!(
        "     S    = {:?}",
        s_l.iter()
            .map(|x| (x * 1e4).round() / 1e4)
            .collect::<Vec<_>>()
    );
    for (name, ys, slack) in &quantities {
        let (a, p, b, _rms) = fit3(&xs, ys);
        let ok = cauchy_ok(ys, *slack);
        let last_step = (ys[ys.len() - 1] - ys[ys.len() - 2]).abs();
        check(
            &format!("{} の収束 (Cauchy 差の減少)", name),
            ok,
            format!(
                "B̂(外挿) = {:.4}, p̂ = {:.2}, 最終差 |f(16)−f(14)| = {:.1e}",
                b, p, last_step
            ),
        );
        map_rows.push((
            name.to_string(),
            a,
            p,
            b,
            if ok { "収束" } else { "非収束" }.into(),
        ));
    }
    check(
        "幾何読み出しの安定性: 隣接復元 100% が全 L (8..16) で成立",
        adj_all,
        "スケールを変えても円環は円環".into(),
    );

    // 弦張力 σ(L) (重物質 m=3, h=0.4)
    {
        let hh = 0.4;
        let ls2: Vec<usize> = vec![10, 12, 14, 16];
        let mut sig_l = Vec::new();
        for &l in &ls2 {
            let nf = l / 2;
            let rmax = 4.min(l / 2 - 2);
            let mut es = Vec::new();
            for r in 0..=rmax {
                let ext = if r == 0 { vec![] } else { vec![3, 3 + r] };
                let g = Z2GaugeRing::try_new(l, nf, w0, hh, 3.0, ext).unwrap();
                let (er, _s, resr) = g.ground_state(5);
                assert!(resr < 1e-6);
                es.push(er);
            }
            let rs: Vec<f64> = (0..=rmax).map(|r| r as f64).collect();
            sig_l.push(linfit(&rs, &es).1);
        }
        let xs2: Vec<f64> = ls2.iter().map(|&l| l as f64).collect();
        let (a, p, b, _) = fit3(&xs2, &sig_l);
        println!(
            "   σ(L) = {:?} (2h = {})",
            sig_l
                .iter()
                .map(|x| (x * 1e4).round() / 1e4)
                .collect::<Vec<_>>(),
            2.0 * hh
        );
        let ok = cauchy_ok(&sig_l, 0.4) && (sig_l[sig_l.len() - 1] / (2.0 * hh) - 1.0).abs() < 0.1;
        check(
            "σ(L) 弦張力の収束と 2h 近傍への安定",
            ok,
            format!("B̂ = {:.4} (2h = {:.2}), p̂ = {:.2}", b, 2.0 * hh, p),
        );
        map_rows.push((
            "σ(L) 弦張力 (m=3)".into(),
            a,
            p,
            b,
            if ok { "収束" } else { "非収束" }.into(),
        ));
    }

    // ================= [3] 前線速度の L 依存 =================
    println!("\n[3] 前線速度推定器の L 依存 — ドリフトの定量化 (推定器の限界も記録)");
    {
        let mut v0_l = Vec::new();
        let mut v12_l = Vec::new();
        for &l in &[10usize, 12, 14, 16] {
            let nf = l / 2;
            for (hh, dst) in [(0.0, &mut v0_l), (1.2, &mut v12_l)] {
                let g = Z2GaugeRing::try_new(l, nf, w0, hh, m0, vec![]).unwrap();
                let (_e, gs, _r) = g.ground_state(20260708);
                let d0 = v2_density_profile(&gs);
                // 到着時刻法 (v15.3 と同一の推定器をインライン再実装せず、
                // ここでは同じ手続きを直接書く — 推定器の定義は v153 参照)
                let mut st = g.apply_bond_op(&gs, 3);
                let nstep = 60;
                let dt = 0.05;
                let mut dev = vec![vec![0.0f64; l]; nstep];
                for s in 0..nstep {
                    st = g.step(&st, dt);
                    let prof = v2_density_profile(&st);
                    for j in 0..l {
                        dev[s][j] = (prof[j] - d0[j]).abs();
                    }
                }
                let center = 3.5;
                let mut pts: Vec<(f64, f64)> = Vec::new();
                for j in 0..l {
                    let mut d = (j as f64 - center).abs();
                    d = d.min(l as f64 - d);
                    if !(1.5..=l as f64 / 2.0 - 1.0).contains(&d) {
                        continue;
                    }
                    let peak = (0..nstep).map(|s| dev[s][j]).fold(0.0f64, f64::max);
                    if peak < 1e-3 {
                        continue;
                    }
                    for s in 0..nstep {
                        if dev[s][j] > 0.5 * peak {
                            pts.push((dt * (s + 1) as f64, d));
                            break;
                        }
                    }
                }
                let v = if pts.len() >= 3 {
                    let ts: Vec<f64> = pts.iter().map(|p| p.0).collect();
                    let ds: Vec<f64> = pts.iter().map(|p| p.1).collect();
                    linfit(&ts, &ds).1.max(0.0)
                } else {
                    0.0
                };
                dst.push(v);
            }
        }
        println!(
            "   v(h=0)   = {:?}",
            v0_l.iter()
                .map(|x| (x * 1e3).round() / 1e3)
                .collect::<Vec<_>>()
        );
        println!(
            "   v(h=1.2) = {:?}",
            v12_l
                .iter()
                .map(|x| (x * 1e3).round() / 1e3)
                .collect::<Vec<_>>()
        );
        // 陰性結果の記録: v̂(L) は L=10..16 で 2w へ向かって単調にドリフトする —
        // 前線幅 ~t^{1/3} と系サイズの競合で、ピーク 50% 到着の推定器はこの窓では
        // 連続極限に達していない。これは物理の破れではなく推定器の未収束であり、
        // 検査は (i) 厳密な LR 上界と (ii) ドリフトの定量記録で行う。
        let vspread = v0_l.iter().cloned().fold(0.0f64, f64::max)
            - v0_l.iter().cloned().fold(f64::INFINITY, f64::min);
        check(
            "v̂(h=0): 全 L で LR 上界 2w 以下 (厳密不等式) — ドリフトは未収束として記録",
            v0_l.iter().all(|&v| v <= 2.0 * w0 * 1.05),
            format!(
                "v̂(L) = {:?}, 幅 {:.2} — 2w へ向かう単調ドリフト (推定器の限界, 陰性記録)",
                v0_l.iter()
                    .map(|x| (x * 1e3).round() / 1e3)
                    .collect::<Vec<_>>(),
                vspread
            ),
        );
        let supp_ok = v0_l.iter().zip(&v12_l).all(|(&v0, &v12)| v12 < 0.8 * v0);
        check(
            "閉じ込め抑制 v(1.2) < 0.8 v(0) が全 L で成立 (物理は L に頑健)",
            supp_ok,
            format!(
                "比 = {:?}",
                v0_l.iter()
                    .zip(&v12_l)
                    .map(|(&a, &b)| ((b / a) * 100.0).round() / 100.0)
                    .collect::<Vec<_>>()
            ),
        );
        let (a, p, b, _) = fit3(&[10.0, 12.0, 14.0, 16.0], &v0_l);
        map_rows.push((
            "v(L) 前線速度 (h=0)".into(),
            a,
            p,
            b,
            "未収束 — 推定器の限界 (外挿 B は無意味; LR 上界 2w が唯一の制約)".into(),
        ));
    }

    // ================= [4] RingChain の中心電荷 c(N) — B=0 の肯定対照 =================
    println!("\n[4] ガウス core (RingChain) の中心電荷 c(N) → 1 (連続極限の既知値)");
    {
        let ns = [66usize, 102, 202, 402];
        let mut cs = Vec::new();
        for &n in &ns {
            let model = RingChain { n };
            let st = model.init();
            // S(ℓ) = (c/3) ln[(N/π) sin(πℓ/N)] + const を 4 点で fit
            let fracs = [0.125, 0.25, 0.375, 0.5];
            let mut xs = Vec::new();
            let mut ys = Vec::new();
            for &f in &fracs {
                let len = ((n as f64 * f) as usize).max(2);
                let chord = (n as f64 / PI) * (PI * len as f64 / n as f64).sin();
                xs.push(chord.ln() / 3.0);
                ys.push(st.readout_entropy(0, len));
            }
            cs.push(linfit(&xs, &ys).1);
        }
        println!(
            "   c(N) = {:?}",
            cs.iter()
                .map(|x| (x * 1e4).round() / 1e4)
                .collect::<Vec<_>>()
        );
        let xs: Vec<f64> = ns.iter().map(|&n| n as f64).collect();
        let (a, p, b, _) = fit3(&xs, &cs);
        check(
            "c(N) の連続極限: B̂ = 1 (Dirac CFT, |B̂−1| < 0.005)",
            (b - 1.0).abs() < 0.005,
            format!("B̂ = {:.5}, p̂ = {:.2}, c(402) = {:.5}", b, p, cs[3]),
        );
        map_rows.push((
            "c(N) 中心電荷 (ガウス core)".into(),
            a,
            p,
            b,
            "B=1 一致".into(),
        ));
    }

    // ================= 監査地図と artifact =================
    println!("\n---- 連続極限監査地図 (readout × [A, p, B]) ----");
    for (name, a, p, b, verdict) in &map_rows {
        println!(
            "  {:34} A={:+.4}  p={:.2}  B={:+.6}  {}",
            name, a, p, b, verdict
        );
    }
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v15.4".into())),
        (
            "rows".into(),
            Json::Arr(
                map_rows
                    .iter()
                    .map(|(n, a, p, b, v)| {
                        Json::Obj(vec![
                            ("readout".into(), Json::Str(n.clone())),
                            ("a".into(), Json::Num(*a)),
                            ("p".into(), Json::Num(*p)),
                            ("b".into(), Json::Num(*b)),
                            ("verdict".into(), Json::Str(v.clone())),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v154_continuum.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 監査した readout 群に消えない残差 (B≠0 型の破れ) は見つからなかった — FAL-CONTINUUM は発火せず"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
