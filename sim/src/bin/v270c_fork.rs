//! v27.0-C v270c_fork — dynamic metric fork 判定の執行 (spec §13.3)
//!
//! 事前登録: spec §13.3 (b46bce3)。解禁 3 条件のうち (i) 厳密 4D Ward
//! (v27.0-A) と (ii) 連続 universality (v27.0-B) は充足、(iii) 予言は §13.3 に
//! 登録済み (「自由場では pole なし = 分岐 (b) 既定」)。本ユニットは凍結
//! プロトコルに従い **1/Π をこの監査に限って計算し、fork を執行する**。
//!
//! 物理 (判定対象):
//!  繰り込み後 P₂ form factor Π₂^ren(q²) := χ_D(q) − c₀ − c₂q² − c₄q⁴
//!  (c 系数は scheme — A₂ = 2A_oracle を凍結して小 q 窓でフィット)。
//!  massless graviton pole (伝播関数 ~ 1/(G q²)) には Π₂ ~ q² の kinetic 項が
//!  必要。自由場の普遍部分は A₂q⁴ln q² のみ → **1/Π₂ に有限留数の pole なし**。
//!  一方 bare c₂ (Sakharov の「誘導 Newton 定数」候補) は **a⁻² で走る
//!  regulator 量** — 本プログラムの universality 基準で非普遍。
//!
//! 検査 (凍結):
//!  [F0] フィット健全性: A₂ 凍結 3 パラメータ {c₀,c₂,c₄} フィットの最大相対
//!       残差 < 1e-3 (q ∈ [0.15, 0.9], a = 0.045)
//!  [F1] 形状: Π₂^ren/(A₂q⁴ln q²) ∈ [0.8, 1.2] 全窓 (log 形に追随・零交差なし)
//!  [F2] **no-pole**: r(q) := |Π₂^ren/q²| が q 半減で単調減少し
//!       r(0.15)/r(0.6) < 0.2 — q² kinetic 項の不在
//!  [F3] **Sakharov 走行**: bare c₂ の格子単位換算が a⁻² で走る —
//!       c₂(a = 0.045)/c₂(a = 0.09) ∈ [3, 5] (比 4 = 二次発散) = 「誘導
//!       Newton 定数」は regulator 量の実証
//!  [F4] 変異: A₂ を 2 倍に凍結 → F0 残差が 10 倍超 (フィットが物理値を選ぶ)
//!
//! 事前登録分岐 (§13.3 の執行):
//!  (a) F0–F3 PASS → **分岐 (b) external metric を確定** — 自由場 matter loop
//!      は graviton を作らない (pole なし)・作りうる q² 項は非普遍。composite
//!      graviton 路線 (分岐 a) は「相互作用による q² 項生成 + Weinberg–Witten
//!      破れ仮定の明示」を伴う将来の別 program として封印。1/Π の常用は
//!      解禁しない (本監査限り)。
//!  (b) F2 破れ (pole あり) → regulator 汚染をまず疑う (縦チャネル監査)。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

fn sbit(s: usize, ax: usize) -> usize {
    (s >> ax) & 1
}

fn h8(k: [f64; 3], m: f64) -> Vec<f64> {
    let mut h = vec![0.0f64; 64];
    for s in 0..8usize {
        let cx = if sbit(s, 0) == 0 { 1.0 } else { -1.0 } * k[0].cos();
        h[s + s * 8] += cx;
        let s2 = s ^ 1;
        let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * k[1].cos();
        h[s2 + s * 8] += cy;
        let s3 = s ^ 3;
        let cz = if sbit(s, 2) == 0 { 1.0 } else { -1.0 } * k[2].cos();
        h[s3 + s * 8] += cz;
        let s4 = s ^ 7;
        h[s4 + s * 8] += m;
    }
    h
}

/// D チャネル頂点 (BOND-A 対角, 中点変調)
fn vd(k: [f64; 3], q: f64) -> Vec<f64> {
    let km = [k[0], k[1] + 0.5 * q, k[2]];
    let r2i = 1.0 / (2.0f64).sqrt();
    let mut v = vec![0.0f64; 64];
    for s in 0..8usize {
        let cx = if sbit(s, 0) == 0 { 1.0 } else { -1.0 } * km[0].cos();
        v[s + s * 8] += cx * r2i;
        let s3 = s ^ 3;
        let cz = if sbit(s, 2) == 0 { 1.0 } else { -1.0 } * km[2].cos();
        v[s3 + s * 8] -= cz * r2i;
    }
    v
}

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut xs = vec![0.0f64; n];
    let mut ws = vec![0.0f64; n];
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0f64, x);
            for kk in 2..=n {
                let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        xs[i] = x;
        let (mut p0, mut p1) = (1.0f64, x);
        for kk in 2..=n {
            let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
            p0 = p1;
            p1 = p2;
        }
        let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
        ws[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (xs, ws)
}

fn nest_edges(center: f64, lo: f64, hi: f64, scale: f64) -> Vec<f64> {
    let mut rs: Vec<f64> = Vec::new();
    let mut r = (1.5 * scale).max(0.006);
    while r < 1.2 {
        rs.push(r);
        r *= 3.0;
    }
    rs.push(1.2);
    let mut e = vec![lo];
    for &rr in rs.iter().rev() {
        e.push(center - rr);
    }
    e.push(center);
    for &rr in rs.iter() {
        e.push(center + rr);
    }
    e.push(hi);
    e
}

/// χ_D(Q; a) — 物理単位 (χ^lat/a⁴), 単一 q
fn chi_d(a: f64, q_phys: f64, nthreads: usize) -> f64 {
    let gl = gauss_legendre(14);
    let q = a * q_phys;
    let edges = nest_edges(PI / 2.0, 0.0, PI, a * q_phys.min(0.3));
    let mut nodes = Vec::new();
    for w2 in edges.windows(2) {
        let (lo, hi) = (w2[0], w2[1]);
        let (cc, hh) = (0.5 * (lo + hi), 0.5 * (hi - lo));
        for (x, wgt) in gl.0.iter().zip(&gl.1) {
            nodes.push((cc + hh * x, wgt * hh));
        }
    }
    let n1 = nodes.len();
    let mut rows: Vec<Option<f64>> = Vec::new();
    rows.resize_with(n1, || None);
    let chunk = n1.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            let nodes = &nodes;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let ix = t * chunk + i;
                    let (kx, wx) = nodes[ix];
                    let mut acc = 0.0f64;
                    for &(ky, wy) in nodes.iter() {
                        for &(kz, wz) in nodes.iter() {
                            let k = [kx, ky, kz];
                            let (wk, vk) = jacobi_eigh(&h8(k, 0.0), 8);
                            let kq = [k[0], k[1] + q, k[2]];
                            let (wq, vq) = jacobi_eigh(&h8(kq, 0.0), 8);
                            let vv = vd(k, q);
                            let mut chi = 0.0f64;
                            for mu in 4..8 {
                                for nu in 0..4 {
                                    let mut mel = 0.0f64;
                                    for r in 0..8 {
                                        let mut s = 0.0f64;
                                        for c in 0..8 {
                                            s += vv[c + r * 8] * vk[c + nu * 8];
                                        }
                                        mel += vq[r + mu * 8] * s;
                                    }
                                    chi += 2.0 * mel * mel / (wq[mu] - wk[nu]);
                                }
                            }
                            acc += wy * wz * chi;
                        }
                    }
                    *slot = Some(acc * wx);
                }
            });
        }
    });
    rows.into_iter().map(|o| o.unwrap()).sum::<f64>() / (2.0 * PI).powi(3) / a.powi(4)
}

/// A₂ 凍結の {c₀, c₂, c₄} 最小二乗: χ(Q) − A₂Q⁴lnQ² = c₀ + c₂Q² + c₄Q⁴
fn fit_c(qs: &[f64], chis: &[f64], a2: f64) -> (f64, f64, f64, f64) {
    let n = qs.len();
    let mut ata = [[0.0f64; 3]; 3];
    let mut atb = [0.0f64; 3];
    for i in 0..n {
        let y = chis[i] - a2 * qs[i].powi(4) * (qs[i] * qs[i]).ln();
        let bs = [1.0, qs[i] * qs[i], qs[i].powi(4)];
        for r in 0..3 {
            for c in 0..3 {
                ata[r][c] += bs[r] * bs[c];
            }
            atb[r] += bs[r] * y;
        }
    }
    // 3×3 ガウス消去
    let mut m = ata;
    let mut rhs = atb;
    for col in 0..3 {
        let piv = (col..3)
            .max_by(|&r1, &r2| m[r1][col].abs().partial_cmp(&m[r2][col].abs()).unwrap())
            .unwrap();
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for r in col + 1..3 {
            let f = m[r][col] / d;
            for c in col..3 {
                m[r][c] -= f * m[col][c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut sol = [0.0f64; 3];
    for col in (0..3).rev() {
        let mut s = rhs[col];
        for c in col + 1..3 {
            s -= m[col][c] * sol[c];
        }
        sol[col] = s / m[col][col];
    }
    // 最大相対残差
    let mut worst = 0.0f64;
    for i in 0..n {
        let pred = sol[0]
            + sol[1] * qs[i] * qs[i]
            + sol[2] * qs[i].powi(4)
            + a2 * qs[i].powi(4) * (qs[i] * qs[i]).ln();
        worst = worst.max(((pred - chis[i]) / chis[i]).abs());
    }
    (sol[0], sol[1], sol[2], worst)
}

fn main() {
    self_test();
    println!("=== v27.0-C v270c_fork — dynamic metric fork 判定の執行 (spec §13.3) ===\n");
    println!("解禁条件 (i)(ii) 充足済み。1/Π をこの監査に限り計算: Π₂^ren = A₂q⁴ln q² に");
    println!("q² kinetic 項がないこと (no-pole) と、bare c₂ の a⁻² 走行 (Sakharov =");
    println!("regulator 量) を判定 → 分岐執行。\n");
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
    let a2 = -2.0 / (160.0 * PI * PI); // A₂ = 2A_oracle (v26.8 凍結・PRED-016 で照合済み)
    let qs: Vec<f64> = vec![0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.75, 0.9];

    // ---- χ_D(Q) @ a = 0.045, 0.09 ----
    let mut chis_f = Vec::new();
    let mut chis_c = Vec::new();
    println!("    [χ_D 表] Q | a = 0.045 | a = 0.09");
    for &q in &qs {
        let cf = chi_d(0.045, q, nthreads);
        let cc = chi_d(0.09, q, nthreads);
        println!("      Q = {:.2}: {:+.6e} | {:+.6e} ({} s)", q, cf, cc, t0.elapsed().as_secs());
        chis_f.push(cf);
        chis_c.push(cc);
    }

    // ---- [F0] フィット健全性 ----
    let (c0f, c2f, c4f, res_f) = fit_c(&qs, &chis_f, a2);
    let (_, c2c, _, _) = fit_c(&qs, &chis_c, a2);
    check(
        "[F0] A₂ 凍結 {c₀,c₂,c₄} フィット: 最大相対残差 < 1e-3 (a = 0.045)",
        res_f < 1e-3,
        format!(
            "残差 = {:.1e} (c₀ = {:+.3e}, c₂ = {:+.3e}, c₄ = {:+.3e})",
            res_f, c0f, c2f, c4f
        ),
    );

    // ---- [F1] 形状 (log 追随・零交差なし) ----
    {
        let mut worst: f64 = 1.0;
        let mut ok = true;
        for (i, &q) in qs.iter().enumerate() {
            let pi_ren = chis_f[i] - c0f - c2f * q * q - c4f * q.powi(4);
            let expect = a2 * q.powi(4) * (q * q).ln();
            let r = pi_ren / expect;
            if !(0.8..1.2).contains(&r) {
                ok = false;
            }
            worst = if (r - 1.0).abs() > (worst - 1.0).abs() { r } else { worst };
        }
        check(
            "[F1] 形状: Π₂^ren/(A₂Q⁴lnQ²) ∈ [0.8, 1.2] 全窓 (零交差なし)",
            ok,
            format!("最大逸脱比 = {:.3}", worst),
        );
    }

    // ---- [F2] no-pole: |Π₂^ren/Q²| → 0 ----
    {
        let rat = |q: f64, chis: &[f64], c0: f64, c2v: f64, c4v: f64| -> f64 {
            let i = qs.iter().position(|&x| (x - q).abs() < 1e-9).unwrap();
            ((chis[i] - c0 - c2v * q * q - c4v * q.powi(4)) / (q * q)).abs()
        };
        let r6 = rat(0.6, &chis_f, c0f, c2f, c4f);
        let r3 = rat(0.3, &chis_f, c0f, c2f, c4f);
        let r15 = rat(0.15, &chis_f, c0f, c2f, c4f);
        // q² kinetic 項が無ければ r(Q) = |A₂|Q²|lnQ²| — 縮小比の導出値は
        // (0.15/0.6)²·ln(0.15²)/ln(0.6²) = 0.232 (log 増強込み)。
        // (開発記録: run1 のバー 0.2 は log 増強を落とした誤較正 — 測定 0.235 は
        // 導出値に 0.9% 一致していた)
        let expect = (0.15f64 / 0.6).powi(2) * (0.15f64 * 0.15).ln() / (0.6f64 * 0.6).ln();
        check(
            "[F2] no-pole: r(Q) = |Π₂^ren/Q²| → 0 が Q²lnQ² 形の導出縮小比 0.232 に 5% 一致",
            r15 < r3 && r3 < r6 && (r15 / r6 / expect - 1.0).abs() < 0.05,
            format!(
                "r = {:.2e} → {:.2e} → {:.2e}, 比 {:.3} (導出 {:.3})",
                r6, r3, r15, r15 / r6, expect
            ),
        );
    }

    // ---- [F3] Sakharov 走行 (bare c₂ ∝ a⁻²) ----
    {
        let ratio = c2f / c2c;
        check(
            "[F3] Sakharov 走行: c₂(0.045)/c₂(0.09) ∈ [3, 5] (a⁻² = 4) — 誘導 Newton 定数候補は regulator 量",
            (3.0..5.0).contains(&ratio),
            format!("比 = {:.3} (c₂ = {:+.3e} vs {:+.3e})", ratio, c2f, c2c),
        );
    }

    // ---- [F4] 変異 ----
    {
        let (_, _, _, res_bad) = fit_c(&qs, &chis_f, 2.0 * a2);
        check(
            "[F4] 変異: A₂ を 2 倍に凍結 → フィット残差 10 倍超 (物理値の選択性)",
            res_bad > 10.0 * res_f,
            format!("変異残差 = {:.1e} vs 正 {:.1e}", res_bad, res_f),
        );
    }

    // ---- fork 執行 ----
    println!("\n[fork 執行 (spec §13.3 の凍結プロトコル)]");
    println!("  判定: **分岐 (b) external metric を確定** —");
    println!("  (1) 自由場 matter loop の普遍部分は A₂q⁴ln q² のみ — 1/Π₂ に有限留数の");
    println!("      massless pole なし (graviton は生成されない)。");
    println!("  (2) pole を作りうる q² 項 (Sakharov 誘導 Newton 定数) は bare c₂ で、");
    println!("      a⁻² で走る regulator 量 — 本プログラムの universality 基準で非普遍。");
    println!("  (3) composite graviton 路線 (分岐 a) は「相互作用による普遍 q² 項の生成 +");
    println!("      Weinberg–Witten 破れ仮定の明示」を要件とする将来の別プログラムに封印。");
    println!("  (4) 1/Π の常用は解禁しない (本監査限り)。metric は外部 regulator のまま —");
    println!("      QRN-Core v1 (matter-on-background) へ。");

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v27.0-C".into())),
        ("kind".into(), Json::Str("dynamic_metric_fork_execution".into())),
        ("spec".into(), Json::Str("§13.3 (b46bce3)".into())),
        ("decision".into(), Json::Str("branch (b) external metric".into())),
        ("c2_fine".into(), Json::Num(c2f)),
        ("c2_coarse".into(), Json::Num(c2c)),
        ("c2_ratio".into(), Json::Num(c2f / c2c)),
    ]);
    let p = write_artifact("results/v270c_fork.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録どおり執行: **分岐 (b) external metric 確定 — 自由場 matter loop は graviton を作らず、作りうる q² 項は regulator 量** (composite 路線は要件明示の上で封印。1/Π 常用は解禁せず)"
        } else {
            "FAIL あり — 分岐 (b') regulator 汚染の監査 (縦チャネル) / 器械。欄が一次ソース"
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
