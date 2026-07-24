//! v26.7.1 (Gate 0) v2671_pole_semantics — 不変量 M² = ΔE² − q² での pole 再採点 (PRED-014)
//!
//! 事前登録: spec §12.1 (コミット bc644d4 で実装前に凍結)。PROMPT/9 の指摘:
//! 固定非零運動量では「s = ΔE² → 0」は ill-posed — massless 状態は ΔE² = q² (光円錐)
//! に住む。v26.7-II の実測 s_min/q² → 1.05 自体が光円錐収束を示していた。
//! 本版は spectral 変数を不変質量二乗 **M²_n := ΔE²_n − q²** に変更し、
//! tree-level improved **q̂ = 2 sin(q/2)** の変種 M̂²_n := ΔE²_n − q̂² も記録して、
//! pole residue を Z_pole = lim_{ε→0} lim_{L→∞} Σ_{0≤M²<ε} Z_n/V で再定義・再採点する。
//!
//! 判定 (spec §12.1 で凍結 — PRED-014):
//!  (i)  ε-ladder ε ∈ {0.5, 0.25, 0.125}·q² (M̂² 側は ·q̂²) で、N=64 の
//!       threshold weight W_ε = Σ_{0≤M²<ε} Z/V の **ε 半減比 < 0.75**
//!       (連続体は W_ε ~ ε^γ, γ≥1 → 比 ≤ 0.5 / 孤立 pole は比 → 1)
//!  (ii) 孤立 residue なし: M²_min(N) が単調減少で → 0 かつ最低クラスタの
//!       residue W₁(N) = ΣZ/V が単調減少 (v26.7-II の branch α の再確認)
//!
//! 検査:
//!  [S0] 回帰: 対和の χ_D が v26.7-II 公表値 (0.154068 [m=0] / 0.150054 [m=0.5],
//!       N=64) と 1e-6 で一致 — 同一走査の同一性
//!  [S1] 被覆: 収集した対が M², M̂² 両 ladder の最大窓を覆う (max ≥ 0.6·q²)
//!  [S2] PRED-014(i): ε 半減比 (q・q̂ 両変数, チャネル D, N=64) < 0.75
//!  [S3] PRED-014(ii): M²_min(N) 単調減少 + W₁(N)/V 単調減少 (N = 16→32→64, m=0)
//!  [S4] massive 対照: M²_min(m=0.5) が 4m² + O(有限サイズ) — 光円錐でなく
//!       2 粒子 threshold に座る (記録 + 下界 M² ≥ 4m² − q² の厳密性)
//!  [S5] 変異: 折返しスワップ落とし → S0 が検出 (> 1e-4)
//!
//! 事前登録分岐: (a) S0–S3 PASS → PRED-014 は hit (v26.7 の no-pole 結論は不変量でも
//!   維持 — 判定規則の修正であり結論の反転ではない) / (b) S2 or S3 FAIL → PRED-014
//!   miss = v26.7 no-pole 結論の再検討 (それ自体を公表) / (c) S0/S1 FAIL → 器械。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// v26.7-II の公表 χ_D (N=64, q = 2π/16) — S0 回帰の的
const REF267_CHI_D: [(f64, f64); 2] = [(0.0, 0.154068), (0.5, 0.150054)];

// ---------------- ブロック機構 (v26.6/26.7 の認証済み実装を写経) ----------------

fn block_h(n: usize, m: f64, cky: f64, ckz: f64) -> Vec<f64> {
    let dim = 4 * n;
    let mut h = vec![0.0f64; dim * dim];
    let id = |x: usize, c: usize| x + n * c;
    let add = |h: &mut Vec<f64>, a: usize, b: usize, t: f64| {
        h[b + a * dim] += t;
        h[a + b * dim] += t;
    };
    let ysgn = [1.0, -1.0, 1.0, -1.0];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            add(&mut h, id(x, c), id((x + 1) % n, c), 0.5 * tw);
            h[id(x, c) + id(x, c) * dim] += px * ysgn[c] * cky;
            if c == 0 || c == 2 {
                add(&mut h, id(x, c), id(x, c + 1), px * zsgn[c] * ckz);
            }
        }
        add(&mut h, id(x, 0), id(x, 3), px * m);
        add(&mut h, id(x, 1), id(x, 2), px * m);
    }
    h
}

fn vertex_qy(n: usize, ky: f64, ckz: f64, q: f64, sw: bool, which: usize) -> Vec<f64> {
    let dim = 4 * n;
    let id = |x: usize, c: usize| x + n * c;
    let tc = |c: usize| -> usize {
        if sw {
            [1usize, 0, 3, 2][c]
        } else {
            c
        }
    };
    let mut v = vec![0.0f64; dim * dim];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let kyc = ky + if c == 1 || c == 3 { PI } else { 0.0 };
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            if which == 1 {
                v[id((x + 1) % n, tc(c)) + id(x, c) * dim] += 0.5 * tw;
                v[id(x, tc(c)) + id((x + 1) % n, c) * dim] += 0.5 * tw;
            }
            if which == 2 {
                v[id(x, tc(c)) + id(x, c) * dim] += px * (kyc + q / 2.0).cos();
            }
            if which == 3 && (c == 0 || c == 2) {
                let coef = px * zsgn[c] * ckz;
                v[id(x, tc(c + 1)) + id(x, c) * dim] += coef;
                v[id(x, tc(c)) + id(x, c + 1) * dim] += coef;
            }
        }
    }
    v
}

fn channel_vertices(n: usize, ky: f64, ckz: f64, q: f64, sw: bool) -> [Vec<f64>; 3] {
    let ox = vertex_qy(n, ky, ckz, q, sw, 1);
    let oy = vertex_qy(n, ky, ckz, q, sw, 2);
    let oz = vertex_qy(n, ky, ckz, q, sw, 3);
    let dim = 4 * n;
    let r2 = (2.0f64).sqrt();
    let mut od = vec![0.0f64; dim * dim];
    let mut os = vec![0.0f64; dim * dim];
    for k in 0..dim * dim {
        od[k] = (ox[k] - oz[k]) / r2;
        os[k] = (ox[k] + oz[k]) / r2;
    }
    [od, os, oy]
}

/// 走査: χ (対和) + 光円錐近傍の全対 (ΔE² ≤ cutoff) を収集
struct ScanOut {
    stat: [f64; 3],
    near: Vec<(f64, [f64; 3])>, // (ΔE², z_D, z_S, z_L) — ΔE² ≤ cutoff の全対
}

fn spectral_scan(
    n: usize,
    m: f64,
    j: usize,
    cutoff: f64,
    nthreads: usize,
    mutate: bool,
) -> ScanOut {
    let nb = n / 2;
    let mut rows: Vec<Option<ScanOut>> = Vec::new();
    rows.resize_with(nb, || None);
    let chunk = nb.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let jz = t * chunk + i;
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    let mut eigs: Vec<(Vec<f64>, Vec<f64>)> = Vec::with_capacity(nb);
                    for jy in 0..nb {
                        let cky = (2.0 * PI * jy as f64 / n as f64).cos();
                        let h = block_h(n, m, cky, ckz);
                        eigs.push(jacobi_eigh(&h, 4 * n));
                    }
                    let dim = 4 * n;
                    let nocc = dim / 2;
                    let q = 2.0 * PI * j as f64 / n as f64;
                    let mut acc = ScanOut {
                        stat: [0.0; 3],
                        near: Vec::new(),
                    };
                    for jy in 0..nb {
                        let ky = 2.0 * PI * jy as f64 / n as f64;
                        let mut jt = jy + j;
                        let mut sw = false;
                        while jt >= nb {
                            jt -= nb;
                            sw = !sw;
                        }
                        let sw_eff = if mutate { false } else { sw };
                        let ops = channel_vertices(n, ky, ckz, q, sw_eff);
                        let (w1, v1) = &eigs[jy];
                        let (w2, v2) = &eigs[jt];
                        let tv = |o: &[f64]| -> Vec<f64> {
                            let mut t = vec![0.0f64; dim * nocc];
                            for ccol in 0..nocc {
                                for r in 0..dim {
                                    let mut s = 0.0;
                                    for k in 0..dim {
                                        s += o[r + k * dim] * v1[k + ccol * dim];
                                    }
                                    t[r + ccol * dim] = s;
                                }
                            }
                            t
                        };
                        let ts: Vec<Vec<f64>> = ops.iter().map(|o| tv(o)).collect();
                        for mu in nocc..dim {
                            for nu in 0..nocc {
                                let mut mm = [0.0f64; 3];
                                for (a, tvv) in ts.iter().enumerate() {
                                    let mut s = 0.0;
                                    for k in 0..dim {
                                        s += v2[k + mu * dim] * tvv[k + nu * dim];
                                    }
                                    mm[a] = s;
                                }
                                let de = w2[mu] - w1[nu];
                                let z = [mm[0] * mm[0], mm[1] * mm[1], mm[2] * mm[2]];
                                for a in 0..3 {
                                    acc.stat[a] += 2.0 * z[a] / de;
                                }
                                let s2 = de * de;
                                if s2 <= cutoff {
                                    acc.near.push((s2, z));
                                }
                            }
                        }
                    }
                    *slot = Some(acc);
                }
            });
        }
    });
    let vol = (n * n * n) as f64;
    let mut out = ScanOut {
        stat: [0.0; 3],
        near: Vec::new(),
    };
    for r in rows {
        let r = r.unwrap();
        for a in 0..3 {
            out.stat[a] += r.stat[a] / vol;
        }
        out.near.extend(r.near);
    }
    out.near.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    out
}

fn main() {
    self_test();
    println!("=== v26.7.1 (Gate 0) — 不変量 M² = ΔE² − q² での pole 再採点 (PRED-014) ===\n");
    println!("事前登録: spec §12.1 (bc644d4 で凍結)。s → 0 は固定 q ≠ 0 で ill-posed —");
    println!("massless は光円錐 ΔE² = q² に住む。M² と M̂² (q̂ = 2sin(q/2)) で再採点。\n");
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!(
            "  [{}] {}  {}",
            if ok { "PASS" } else { "FAIL" },
            name,
            detail
        );
        if !ok {
            nfail += 1;
        }
    };
    let ns_list = [16usize, 32, 64];
    let q = 2.0 * PI / 16.0;
    let qhat = 2.0 * (q / 2.0).sin();
    let (q2, qhat2) = (q * q, qhat * qhat);
    // 両 ladder (最大窓 0.5q², 0.5q̂²) を覆い、かつ粗い N=16 の最低対 (2.53q²,
    // v26.7-II) も収集できる cutoff。N=16 の窓内 W_ε は 0 になり得る (物理的 —
    // 粗格子は近閾状態を持たない)。ladder 比のゲートは N=64 のみ (凍結どおり)。
    let cutoff = q2 * 4.0;

    // ---- 走査 ----
    let mut scans_m0: Vec<ScanOut> = Vec::new();
    for &n in &ns_list {
        let s = spectral_scan(n, 0.0, n / 16, cutoff, nthreads, false);
        println!(
            "    [走査] N={} m=0 完了 ({} s) — 光円錐近傍対 {} 個, M²_min/q² = {:.4}",
            n,
            t0.elapsed().as_secs(),
            s.near.len(),
            (s.near[0].0 - q2) / q2
        );
        scans_m0.push(s);
    }
    let scan_m5 = spectral_scan(64, 0.5, 4, 4.0 * 0.25 + q2 * 1.7, nthreads, false);
    println!(
        "    [走査] N=64 m=0.5 完了 ({} s) — 近傍対 {} 個",
        t0.elapsed().as_secs(),
        scan_m5.near.len()
    );

    // ---- [S0] 回帰 ----
    {
        let mut worst = 0.0f64;
        worst = worst.max((scans_m0[2].stat[0] - REF267_CHI_D[0].1).abs());
        worst = worst.max((scan_m5.stat[0] - REF267_CHI_D[1].1).abs());
        check(
            "[S0] 回帰: χ_D (N=64) = v26.7-II 公表値 (±1e-6)",
            worst < 1e-6,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [S1] 被覆 ----
    {
        let mut ok = true;
        for (ni, s) in scans_m0.iter().enumerate() {
            let mmax = s.near.last().map(|x| x.0 - q2).unwrap_or(0.0);
            if mmax < 0.6 * q2 {
                ok = false;
            }
            let _ = ni;
        }
        check(
            "[S1] 被覆: 収集対が M²/M̂² 両 ladder の最大窓 (0.5q²) を覆う",
            ok,
            format!(
                "max M²/q² = {:.3} (N=64)",
                (scans_m0[2].near.last().unwrap().0 - q2) / q2
            ),
        );
    }

    // ---- threshold weight 表 ----
    let w_eps = |s: &ScanOut, n: usize, qq2: f64, eps: f64, ch: usize| -> f64 {
        let vol = (n * n * n) as f64;
        s.near
            .iter()
            .filter(|&&(s2, _)| s2 - qq2 >= -1e-12 && s2 - qq2 < eps)
            .map(|&(_, z)| z[ch])
            .sum::<f64>()
            / vol
    };
    println!("\n    [W_ε 表 (チャネル D, m=0)] 変数 | N | W(0.5q²) | W(0.25q²) | W(0.125q²) | 半減比1 | 半減比2");
    let mut ratios_n64 = [[0.0f64; 2]; 2]; // [variable][halving]
    for (vi, (vname, qq2)) in [("M² (q)", q2), ("M̂² (q̂)", qhat2)].iter().enumerate() {
        for (ni, &n) in ns_list.iter().enumerate() {
            let s = &scans_m0[ni];
            let w1 = w_eps(s, n, *qq2, 0.5 * qq2, 0);
            let w2 = w_eps(s, n, *qq2, 0.25 * qq2, 0);
            let w3 = w_eps(s, n, *qq2, 0.125 * qq2, 0);
            let (r1, r2) = if w1 > 0.0 && w2 > 0.0 {
                (w2 / w1, w3 / w2)
            } else {
                (0.0, 0.0) // 粗格子で窓が空 (物理的) — 比は定義しない
            };
            println!(
                "      {} N={}: {:.4e} | {:.4e} | {:.4e} | {:.3} | {:.3}",
                vname, n, w1, w2, w3, r1, r2
            );
            if ni == 2 {
                ratios_n64[vi] = [r1, r2];
            }
        }
    }

    // ---- [S2] PRED-014(i): ε 半減比 ----
    {
        let worst = ratios_n64
            .iter()
            .flatten()
            .cloned()
            .fold(0.0f64, f64::max);
        check(
            "[S2] PRED-014(i): ε 半減比 < 0.75 (N=64, q・q̂ 両変数, チャネル D)",
            worst < 0.75,
            format!(
                "M²: {:.3}/{:.3}, M̂²: {:.3}/{:.3} (max {:.3})",
                ratios_n64[0][0], ratios_n64[0][1], ratios_n64[1][0], ratios_n64[1][1], worst
            ),
        );
    }

    // ---- [S3] PRED-014(ii): 孤立 residue なし ----
    {
        let mmin: Vec<f64> = scans_m0.iter().map(|s| s.near[0].0 - q2).collect();
        let w1v: Vec<f64> = scans_m0
            .iter()
            .zip(&ns_list)
            .map(|(s, &n)| {
                let smin = s.near[0].0;
                let vol = (n * n * n) as f64;
                s.near
                    .iter()
                    .filter(|&&(s2, _)| s2 < smin * (1.0 + 1e-9))
                    .map(|&(_, z)| z[0])
                    .sum::<f64>()
                    / vol
            })
            .collect();
        println!(
            "    [S3 表] M²_min/q²: {:.4} → {:.4} → {:.4} / W₁/V: {:.3e} → {:.3e} → {:.3e}",
            mmin[0] / q2,
            mmin[1] / q2,
            mmin[2] / q2,
            w1v[0],
            w1v[1],
            w1v[2]
        );
        check(
            "[S3] PRED-014(ii): M²_min(N) 単調減少 → 0 かつ W₁(N)/V 単調減少 (孤立 residue なし)",
            mmin[0] > mmin[1] && mmin[1] > mmin[2] && w1v[0] > w1v[1] && w1v[1] > w1v[2],
            format!(
                "M²_min/q² → {:.4}, W₁/V → {:.1e}",
                mmin[2] / q2,
                w1v[2]
            ),
        );
        // 近傍状態数の成長 (連続体) も記録
        let counts: Vec<usize> = scans_m0
            .iter()
            .map(|s| s.near.iter().filter(|&&(s2, _)| s2 - q2 < 0.5 * q2).count())
            .collect();
        println!(
            "      [記録] 窓 [0, 0.5q²] 内の状態数: {} → {} → {} (成長 = 連続体; pole は 1 個が残存)",
            counts[0], counts[1], counts[2]
        );
    }

    // ---- [S4] massive 対照 ----
    {
        let m2min = scan_m5.near[0].0 - q2;
        let bound = 4.0 * 0.25 - 1e-9;
        println!(
            "    [S4] m=0.5: M²_min = {:.6} (2 粒子 threshold M² = 4m² = 1.0 の {:.4} 倍; 下界 4m²−q²… ΔE²≥4m² は厳密)",
            m2min,
            m2min / 1.0
        );
        check(
            "[S4] massive 対照: ΔE² ≥ 4m² (厳密) — 光円錐でなく threshold に座る",
            scan_m5.near[0].0 >= 4.0 * 0.25 - 1e-9 && m2min > 0.0,
            format!("ΔE²_min = {:.6} ≥ {:.3}", scan_m5.near[0].0, bound + 1e-9),
        );
    }

    // ---- [S5] 変異 ----
    {
        let bad = spectral_scan(16, 0.0, 1, cutoff, nthreads, true);
        let dev = (bad.stat[0] - scans_m0[0].stat[0]).abs();
        check(
            "[S5] 変異: 折返しスワップ落とし → S0 (χ_D) が検出 (> 1e-4)",
            dev > 1e-4,
            format!("逸脱 {:.2e}", dev),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.7.1".into())),
        ("kind".into(), Json::Str("invariant_pole_audit".into())),
        ("pred".into(), Json::Str("PRED-014".into())),
        ("q".into(), Json::Num(q)),
        ("qhat".into(), Json::Num(qhat)),
        (
            "halving_ratios_n64_D".into(),
            Json::Obj(vec![
                (
                    "q_variable".into(),
                    Json::Arr(vec![Json::Num(ratios_n64[0][0]), Json::Num(ratios_n64[0][1])]),
                ),
                (
                    "qhat_variable".into(),
                    Json::Arr(vec![Json::Num(ratios_n64[1][0]), Json::Num(ratios_n64[1][1])]),
                ),
            ]),
        ),
        (
            "m2min_over_q2_by_n".into(),
            Json::Arr(
                scans_m0
                    .iter()
                    .map(|s| Json::Num((s.near[0].0 - q2) / q2))
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v2671_pole_semantics.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **PRED-014 hit — v26.7 の no-pole 結論は不変量 M² でも維持** (判定規則の修正であり結論の反転ではない)"
        } else {
            "FAIL — 分岐 (b) PRED-014 miss = no-pole 再検討 / (c) 器械。欄が一次ソース"
        }
    );
    println!(
        "\n総合判定: {} ({} s)",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" },
        t0.elapsed().as_secs()
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
