//! v26.8-S v268s_sumrule — スカラー spectral 和則の格子側 (PRED-017 採点)
//!
//! 事前登録: spec §12.7 (bc644d4) + PRED-017。taste-singlet scalar 演算子
//! θ = m ψ̄ψ の spectral 和則 (v26.8-A で自前規約導出・三重検証済み):
//!   ∫ ρ_θ(s)/s³ ds = 1/(80π²) per Dirac flavor
//! (UV 有限・local counterterm 非依存・質量非依存 — 文献規約 σ_f = ρ_θ/(3s³) では
//! 1/(240π²)/flavor, 2 taste で 1/(120π²) = PRED-017 の的)。**null 結合すら不要の
//! 単一 BZ 積分**で格子側を測る。凍結不変量は S = 2∫σ_θ(E)/E⁵ dE (v268a が
//! 3.5e-12 で三重検証した形。σ_θ(E) は E = ΔE の遷移密度) — q = 0 では対は
//! (±k)、ΔE = 2E(k) なので
//!   S^lat(a) = ∫_cell d³k/(2π)³ · 2·m²·T(k) / (2E(k))⁵
//! T(k) = Σ_pairs |⟨unocc|V_θ|occ⟩|² は射影子トレースの閉形式:
//!   staggered (8 成分, V_θ = Γm):  T = 4 Σcos²kᵢ / E²   (E² = Σcos²+ m²)
//!   Wilson    (4 成分, V_θ = β):   T = 2 Σsin²kᵢ / E²   (E² = Σsin²+ M(k)²)
//! (導出: tr[P₊ V P₋ V] = (1/4E²)[dim·E² − tr(HVHV)] — 本バイナリが明示行列
//! トレースと機械照合する)。trajectory: m^lat = a·m_phys。S は無次元かつ
//! m_phys 非依存 (連続極限) — a → 0 で
//!   **staggered → 2/(80π²) (2 tastes) / Wilson → 1/(80π²) (1 flavor)**
//! に収束するかが判定。Wilson の doubler (M ~ 2r) は 1/E⁵ 抑制で a² 級に脱結合。
//!
//! 外挿モデル (S 自身の漸近形から導出 — run2 の教訓を見よ):
//!   staggered: ノードの sin²p = p²(1 − p²/3 + …) 補正が ∫ m²p⁶/E⁹ d³p ~
//!     m²ln(1/(am)) を生む → δS = c·a²ln(1/a) + d·a² (ln(1/m_phys) は d に吸収)。
//!     数値でも Δ/(a²ln(1/a)) = 2.44–2.47 で定数 — モデル {1, a²ln(1/a), a²}。
//!   Wilson: r 項の交差項 2 r m p² (質量の p² 依存シフト) が O(m^lat) = O(a) を
//!     生む → δS = c·a + d·a² — モデル {1, a, a²}。
//! 各モデルを全域と尾部窓 (a ≤ 0.125) でフィットし相互 spread も偏差に数える。
//!
//! 検査 (凍結):
//!  [S0] T の閉形式 = 明示行列トレース (stag 8×8 / Wil 4×4, 決定的 k 点): 1e-12
//!  [S1] 求積自己整合: a_min で GL 14 → 18 < 1e-3
//!  [S2] **PRED-017 staggered**: 導出モデル {1, a²ln(1/a), a²} の全域/尾部窓
//!       2 フィット + spread で S₀·(80π²)/2 = 1 ± 0.02
//!  [S3] **Wilson**: 導出モデル {1, a, a²} の全域/尾部窓 2 フィットで
//!       S₀·(80π²) = 1 ± 0.02 (taste 1 — staggered の 2 との対比が taste 数の三たびの検証)
//!  [S4] m_phys 非依存: m_phys = 0.5 と 1.0 の外挿 S₀ が 1% で一致
//!  [S5] 変異: staggered の V_θ を Γm → 恒等 (taste-nonsinglet 密度) に置換 →
//!       S が > 10% 変化
//!
//! 事前登録分岐: (a) S0–S4 PASS → **PRED-017 scored-hit** (辞書換算で 1/(120π²)
//!   の的に一致 — operator benchmark として。gravitational response とは呼ばない) /
//!   (b) S2/S3 miss → scalar 規格化または taste count の誤り (PRED-017 の falsifier) /
//!   (c) S0/S1 FAIL → 器械。
//!
//! 開発記録 (物理は不変・器械の修正のみ):
//!   run1: ∫ρ(s)/s³ の s-measure を直書きし E-Jacobian 分 (×4E) を落とした →
//!     1/m 発散 (ladder が発散、S4 の m 依存 FAIL が正直に検出)。凍結不変量は
//!     σ-規約 S = 2∫σ/E⁵ dE の方 — 「和則は密度と測度の組で一つ」。
//!   run2: 外挿モデル (a², a+a²)/(a+a², a²+a⁴) は null 結合観測量 A (v268b/c) の
//!     漸近形の転記で 1/E⁵ 観測量 S には合わない (stag spread 0.07, Wil 0.13) →
//!     上の導出モデルに置換。教訓: **外挿モデルは観測量ごとに導出する — 別観測量
//!     からの転記は「凍結」ではなく汚染**。副産物: a²ln(1/a) 項の導出は PRED-016
//!     精密化のモデル正当化にそのまま使える。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

// ---------------- staggered 8 成分 (閉形式 + 照合用明示行列) ----------------

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

/// staggered T(k) 明示: Γm 頂点 (mutate: 恒等頂点 = taste-nonsinglet 密度) の
/// 占有→非占有 |M|² 和 (jacobi による)
fn t_stag_explicit(k: [f64; 3], m: f64, mutate: bool) -> f64 {
    let h = h8(k, m);
    let (w, v) = jacobi_eigh(&h, 8);
    // Γm = 全 flip (s ^ 7), 恒等 = δ
    let mut t = 0.0f64;
    for mu in 4..8 {
        for nu in 0..4 {
            let mut mel = 0.0f64;
            for s in 0..8usize {
                let s2 = if mutate { s } else { s ^ 7 };
                mel += v[s2 + mu * 8] * v[s + nu * 8];
            }
            t += mel * mel;
        }
    }
    let _ = w;
    t
}

/// staggered T(k) 閉形式: 4 Σcos²kᵢ / E²
fn t_stag_closed(k: [f64; 3], m: f64) -> f64 {
    let c2 = k[0].cos().powi(2) + k[1].cos().powi(2) + k[2].cos().powi(2);
    4.0 * c2 / (c2 + m * m)
}

fn e_stag(k: [f64; 3], m: f64) -> f64 {
    (k[0].cos().powi(2) + k[1].cos().powi(2) + k[2].cos().powi(2) + m * m).sqrt()
}

// ---------------- Wilson 4 成分 (閉形式 + 照合用明示行列) ----------------

type M4 = [(f64, f64); 16];

fn mzero() -> M4 {
    [(0.0, 0.0); 16]
}

fn mmul(a: &M4, b: &M4) -> M4 {
    let mut o = mzero();
    for r in 0..4 {
        for k in 0..4 {
            let av = a[k + r * 4];
            if av.0 == 0.0 && av.1 == 0.0 {
                continue;
            }
            for c in 0..4 {
                let bv = b[c + k * 4];
                o[c + r * 4].0 += av.0 * bv.0 - av.1 * bv.1;
                o[c + r * 4].1 += av.0 * bv.1 + av.1 * bv.0;
            }
        }
    }
    o
}

fn mtrace_re(a: &M4) -> f64 {
    (0..4).map(|i| a[i + i * 4].0).sum()
}

fn alphas() -> [M4; 4] {
    let mut ax = mzero();
    let mut ay = mzero();
    let mut az = mzero();
    let mut b = mzero();
    let put = |m: &mut M4, r: usize, c: usize, re: f64, im: f64| {
        m[c + r * 4] = (re, im);
    };
    for blk in 0..2 {
        let (ro, co) = if blk == 0 { (0, 2) } else { (2, 0) };
        put(&mut ax, ro, co + 1, 1.0, 0.0);
        put(&mut ax, ro + 1, co, 1.0, 0.0);
        put(&mut ay, ro, co + 1, 0.0, -1.0);
        put(&mut ay, ro + 1, co, 0.0, 1.0);
        put(&mut az, ro, co, 1.0, 0.0);
        put(&mut az, ro + 1, co + 1, -1.0, 0.0);
    }
    for i in 0..4 {
        put(&mut b, i, i, if i < 2 { 1.0 } else { -1.0 }, 0.0);
    }
    [ax, ay, az, b]
}

fn mwil(k: [f64; 3], m: f64, r: f64) -> f64 {
    m + r * (3.0 - k[0].cos() - k[1].cos() - k[2].cos())
}

fn e_wil(k: [f64; 3], m: f64, r: f64) -> f64 {
    let s2 = k[0].sin().powi(2) + k[1].sin().powi(2) + k[2].sin().powi(2);
    let mk = mwil(k, m, r);
    (s2 + mk * mk).sqrt()
}

fn t_wil_explicit(al: &[M4; 4], k: [f64; 3], m: f64, r: f64) -> f64 {
    // P± を閉形式で作り tr[P₊ β P₋ β]
    let mut h = mzero();
    for ax in 0..3 {
        let s = k[ax].sin();
        for i in 0..16 {
            h[i].0 += s * al[ax][i].0;
            h[i].1 += s * al[ax][i].1;
        }
    }
    let mk = mwil(k, m, r);
    for i in 0..16 {
        h[i].0 += mk * al[3][i].0;
        h[i].1 += mk * al[3][i].1;
    }
    let e = e_wil(k, m, r);
    let mut pp = mzero();
    let mut pm = mzero();
    for i in 0..16 {
        pp[i].0 = h[i].0 / (2.0 * e);
        pp[i].1 = h[i].1 / (2.0 * e);
        pm[i].0 = -h[i].0 / (2.0 * e);
        pm[i].1 = -h[i].1 / (2.0 * e);
    }
    for i in 0..4 {
        pp[i + i * 4].0 += 0.5;
        pm[i + i * 4].0 += 0.5;
    }
    mtrace_re(&mmul(&mmul(&pp, &al[3]), &mmul(&pm, &al[3])))
}

/// Wilson T(k) 閉形式: 2 Σsin²kᵢ / E²
fn t_wil_closed(k: [f64; 3], m: f64, r: f64) -> f64 {
    let s2 = k[0].sin().powi(2) + k[1].sin().powi(2) + k[2].sin().powi(2);
    2.0 * s2 / (e_wil(k, m, r).powi(2))
}

// ---------------- GL / BZ ----------------

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

fn make_nodes(edges: &[f64], gl: &(Vec<f64>, Vec<f64>)) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for w2 in edges.windows(2) {
        let (a, b) = (w2[0], w2[1]);
        let (cc, hh) = (0.5 * (a + b), 0.5 * (b - a));
        for (x, wgt) in gl.0.iter().zip(&gl.1) {
            out.push((cc + hh * x, wgt * hh));
        }
    }
    out
}

fn bz_sum<F: Fn([f64; 3]) -> f64 + Sync>(nodes: &[(f64, f64)], nthreads: usize, f: F) -> f64 {
    let n1 = nodes.len();
    let mut rows: Vec<Option<f64>> = Vec::new();
    rows.resize_with(n1, || None);
    let chunk = n1.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            let f = &f;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let ix = t * chunk + i;
                    let (kx, wx) = nodes[ix];
                    let mut acc = 0.0f64;
                    for &(ky, wy) in nodes.iter() {
                        for &(kz, wz) in nodes.iter() {
                            acc += wy * wz * f([kx, ky, kz]);
                        }
                    }
                    *slot = Some(acc * wx);
                }
            });
        }
    });
    rows.into_iter().map(|o| o.unwrap()).sum::<f64>() / (2.0 * PI).powi(3)
}

fn fit_a0(avals: &[(f64, f64)], basis: &dyn Fn(f64) -> Vec<f64>) -> f64 {
    let p = basis(avals[0].0).len();
    let mut ata = vec![0.0f64; p * p];
    let mut atb = vec![0.0f64; p];
    for &(a, y) in avals {
        let bs = basis(a);
        for i in 0..p {
            for j in 0..p {
                ata[j + i * p] += bs[i] * bs[j];
            }
            atb[i] += bs[i] * y;
        }
    }
    let mut m = ata;
    let mut rr = atb;
    for col in 0..p {
        let piv = (col..p)
            .max_by(|&r1, &r2| m[col + r1 * p].abs().partial_cmp(&m[col + r2 * p].abs()).unwrap())
            .unwrap();
        for c in 0..p {
            m.swap(c + col * p, c + piv * p);
        }
        rr.swap(col, piv);
        let d = m[col + col * p];
        for row in col + 1..p {
            let f = m[col + row * p] / d;
            for c in col..p {
                m[c + row * p] -= f * m[c + col * p];
            }
            rr[row] -= f * rr[col];
        }
    }
    let mut sol = vec![0.0f64; p];
    for col in (0..p).rev() {
        let mut s = rr[col];
        for c in col + 1..p {
            s -= m[c + col * p] * sol[c];
        }
        sol[col] = s / m[col + col * p];
    }
    sol[0]
}

fn main() {
    self_test();
    println!("=== v26.8-S v268s_sumrule — スカラー spectral 和則の格子側 (PRED-017) ===\n");
    println!("事前登録: spec §12.7 (bc644d4)。的 (自前規約): staggered → 2/(80π²) / Wilson →");
    println!("1/(80π²)。辞書 (×1/3s³) で文献規約の 1/(120π²) (2t) に対応。UV 有限・m 非依存。\n");
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
    let al = alphas();
    let target = 1.0 / (80.0 * PI * PI);

    // ---- [S0] 閉形式 = 明示行列トレース ----
    {
        let mut worst = 0.0f64;
        for &(k, m) in &[
            ([0.3f64, -1.1, 2.2], 0.2f64),
            ([1.4, 0.5, -0.7], 0.05),
            ([0.9, 2.6, 0.4], 0.5),
        ] {
            worst = worst.max((t_stag_explicit(k, m, false) - t_stag_closed(k, m)).abs());
            worst = worst.max((t_wil_explicit(&al, k, m, 1.0) - t_wil_closed(k, m, 1.0)).abs());
        }
        check(
            "[S0] T の閉形式 (stag: 4Σcos²/E² / Wil: 2Σsin²/E²) = 明示行列トレース",
            worst < 1e-12,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- ladder ----
    let ladder = [0.35f64, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.032, 0.022];
    let gl14 = gauss_legendre(14);
    // 入れ子半径 (×3 等比, 最内 1.5·a·m_phys): staggered セル [0,π)³ (ノード中心
    // π/2)、Wilson [−π,π]³ (ノード 0)
    let nest_edges = |center: f64, lo: f64, hi: f64, scale: f64| -> Vec<f64> {
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
    };
    let stag_edges = |a: f64, m_phys: f64| nest_edges(PI / 2.0, 0.0, PI, a * m_phys);
    let wil_edges = |a: f64, m_phys: f64| nest_edges(0.0, -PI, PI, a * m_phys);
    let run_stag = |m_phys: f64, ngl: usize| -> Vec<(f64, f64)> {
        let gl = gauss_legendre(ngl);
        ladder
            .iter()
            .map(|&a| {
                let m = a * m_phys;
                let nodes = make_nodes(&stag_edges(a, m_phys), &gl);
                let s = bz_sum(&nodes, nthreads, |k| {
                    let e = e_stag(k, m);
                    2.0 * m * m * t_stag_closed(k, m) / (2.0 * e).powi(5)
                });
                (a, s)
            })
            .collect()
    };
    let run_wil = |m_phys: f64, ngl: usize| -> Vec<(f64, f64)> {
        let gl = gauss_legendre(ngl);
        ladder
            .iter()
            .map(|&a| {
                let m = a * m_phys;
                let nodes = make_nodes(&wil_edges(a, m_phys), &gl);
                let s = bz_sum(&nodes, nthreads, |k| {
                    let e = e_wil(k, m, 1.0);
                    2.0 * m * m * t_wil_closed(k, m, 1.0) / (2.0 * e).powi(5)
                });
                (a, s)
            })
            .collect()
    };
    let s_stag = run_stag(1.0, 14);
    let s_wil = run_wil(1.0, 14);
    println!("    [S(a) 表 (m_phys = 1)] a | stag·80π²/2 | Wil·80π²");
    for i in 0..ladder.len() {
        println!(
            "      a = {:.3}: {:.4} | {:.4} ({} s)",
            ladder[i],
            s_stag[i].1 * 80.0 * PI * PI / 2.0,
            s_wil[i].1 * 80.0 * PI * PI,
            t0.elapsed().as_secs()
        );
    }

    // ---- [S1] 求積自己整合 ----
    {
        let s_fine = run_stag(1.0, 18);
        let w_fine = run_wil(1.0, 18);
        let rel = (s_fine.last().unwrap().1 / s_stag.last().unwrap().1 - 1.0)
            .abs()
            .max((w_fine.last().unwrap().1 / s_wil.last().unwrap().1 - 1.0).abs());
        check(
            "[S1] 求積自己整合: a_min で GL 14 → 18 の変化 < 1e-3 (両離散化)",
            rel < 1e-3,
            format!("max 相対変化 = {:.1e}", rel),
        );
    }

    // ---- [S2] staggered → 2/(80π²) (導出モデル {1, a²ln(1/a), a²}) ----
    let b_stag = |a: f64| vec![1.0, a * a * (1.0 / a).ln(), a * a];
    let st_full = fit_a0(&s_stag, &b_stag);
    let st_tail = fit_a0(&s_stag[3..], &b_stag);
    {
        let (rf, rt) = (st_full / (2.0 * target), st_tail / (2.0 * target));
        let dev = (rf - 1.0).abs().max((rt - 1.0).abs()).max((rf - rt).abs());
        println!(
            "    [S2 外挿 stag] S₀·80π²/2 = {:.4} (全域) / {:.4} (尾部 a ≤ 0.125), spread {:.4}",
            rf, rt, (rf - rt).abs()
        );
        check(
            "[S2] PRED-017 staggered: S₀·(80π²)/2 = 1 ± 0.02 (導出モデル, 全域/尾部窓)",
            dev < 0.02,
            format!("偏差 max = {:.4}", dev),
        );
    }

    // ---- [S3] Wilson → 1/(80π²) (導出モデル {1, a, a²}) ----
    let b_wil = |a: f64| vec![1.0, a, a * a];
    let wl_full = fit_a0(&s_wil, &b_wil);
    let wl_tail = fit_a0(&s_wil[3..], &b_wil);
    {
        let (rf, rt) = (wl_full / target, wl_tail / target);
        let dev = (rf - 1.0).abs().max((rt - 1.0).abs()).max((rf - rt).abs());
        println!(
            "    [S3 外挿 Wil] S₀·80π² = {:.4} (全域) / {:.4} (尾部 a ≤ 0.125), spread {:.4}",
            rf, rt, (rf - rt).abs()
        );
        check(
            "[S3] PRED-017 Wilson: S₀·80π² = 1 ± 0.02 (導出モデル {1,a,a²} — taste 1 の検証)",
            dev < 0.02,
            format!("偏差 max = {:.4}", dev),
        );
    }

    // ---- [S4] m_phys 非依存 ----
    {
        let s_m05 = run_stag(0.5, 14);
        let st05 = fit_a0(&s_m05, &b_stag);
        let rel = (st05 / st_full - 1.0).abs();
        check(
            "[S4] m_phys 非依存: staggered の外挿 S₀ が m_phys = 0.5 と 1.0 で 1% 一致",
            rel < 0.01,
            format!("相対差 = {:.3}", rel),
        );
    }

    // ---- [S5] 変異 ----
    {
        let a = 0.125f64;
        let m = a * 1.0;
        let nodes = make_nodes(&stag_edges(a, 1.0), &gl14);
        let bad = bz_sum(&nodes, nthreads, |k| {
            let e = e_stag(k, m);
            2.0 * m * m * t_stag_explicit(k, m, true) / (2.0 * e).powi(5)
        });
        let good = s_stag[3].1;
        check(
            "[S5] 変異: V_θ を Γm → 恒等 (taste-nonsinglet 密度) に置換 → S が > 10% 変化",
            (bad / good - 1.0).abs() > 0.10,
            format!("相対変化 = {:.3}", (bad / good - 1.0).abs()),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.8-S".into())),
        ("kind".into(), Json::Str("scalar_sum_rule_lattice".into())),
        ("pred".into(), Json::Str("PRED-017".into())),
        ("target_per_flavor".into(), Json::Num(target)),
        ("stag_s0_ratio_2t".into(), Json::Num(st_full / (2.0 * target))),
        ("wil_s0_ratio_1f".into(), Json::Num(wl_full / target)),
        (
            "stag_ladder_ratio".into(),
            Json::Arr(
                s_stag
                    .iter()
                    .map(|&(a, v)| {
                        Json::Obj(vec![
                            ("a".into(), Json::Num(a)),
                            ("ratio".into(), Json::Num(v * 80.0 * PI * PI / 2.0)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v268s_sumrule.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **PRED-017 scored-hit — 両離散化のスカラー和則が taste 数どおり連続値に収束** (operator benchmark — gravitational response とは呼ばない)"
        } else {
            "FAIL — 分岐 (b) scalar 規格化/taste count の誤り / (c) 器械。欄が一次ソース"
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
