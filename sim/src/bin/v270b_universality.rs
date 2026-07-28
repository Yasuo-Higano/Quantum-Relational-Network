//! v27.0-B v270b_universality — 4D kernel の連続 universality (spec §13.2-B)
//!
//! 事前登録: spec §13.2 (b46bce3)。P₂ 幾何 + Lean 証明済みチャネル辞書
//! (Projector.lean: L = yy = 純ゲージ / ProjectorND: D/X = P₂ 固有) から、
//! 静的 (q₀ = 0, q_L = (0, Qŷ) 空間的) の kernel チャネルの q⁴ln q² 係数は
//! 全て P₂ form factor A₂ 一つで決まる (massless — ρ₀ = 0):
//!   A_I = A₂ · ⟨I|P₂|I⟩,  θ₀₀ = 1・θ_yy = θ₀y = 0 (静的) より
//!   **A_00/A_D = P₂_{00,00} = 2/3 (厳密) / A_0y = 0 / A_yy = 0**
//! (A_D = A₂ は PRED-016 で oracle 照合済み)。温度セクター (h₀₀/h₀y source)
//! の連続極限が空間セクターと同一の form factor に流れるかの判定 —
//! spec §13.3 の 1/Π 解禁条件 (ii)。
//!
//! 器械: v268p (§12.9 凍結プロトコル) の null 結合 ladder を channels
//! {D, 00, 0y, yy} に適用。source は v26.9/v27.0-A の認証済み構成
//! (V₀₀ = h(k+q/2ŷ)・V₀y = −(1/2)sin(2k_y+q)𝟙・BOND-A yy・D = (xx−zz)/√2)。
//! ladder a ∈ {0.35 … 0.022}・導出モデル {1, a²ln(1/a), a²} 全域+尾部窓。
//!
//! 検査 (凍結):
//!  [B0] 回帰: A_D の a = 0.125 rung が v268p 公表 JSON (ratio 1.022197) と
//!       2e-3 相対一致
//!  [B1] **A_00/(⅔·A_D) = 1 ± 2%** (中心 = 尾部窓, 系統 = |全域 − 尾部| ≤ 1%)
//!  [B2] |A_0y/A_D| の外挿 → 0 (最終 |比| < 0.05)
//!  [B3] |A_yy/A_D| の外挿 → 0 (最終 |比| < 0.05 — Lean の「L = 純ゲージ」の
//!       one-loop 版)
//!  [B4] 変異: V₀₀ ×1.02 → B1 が 4% 逸脱 (< 2% ゲートの外)
//!
//! 事前登録分岐: (a) 全 PASS → **温度セクター universality 成立 — §13.3
//!   条件 (ii) 充足** (fork 判定 v27.0-C へ) / (b) B1 破れ → 温度 source の
//!   正規化または P₂ 幾何の誤り (公表) / (c) B0 FAIL → 器械。

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

/// チャネル頂点 (実行列): 0 = D, 1 = 00, 2 = 0y, 3 = yy。mutate: V₀₀ ×1.02
fn channel(ch: usize, k: [f64; 3], q: f64, mutate: bool) -> Vec<f64> {
    let km = [k[0], k[1] + 0.5 * q, k[2]];
    match ch {
        0 => {
            // D = (V_xx − V_zz)/√2 (BOND-A 対角: h8 片そのもの)
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
        1 => {
            let v = h8(km, 0.0);
            if mutate {
                v.iter().map(|x| 1.02 * x).collect()
            } else {
                v
            }
        }
        2 => {
            let val = -0.5 * (2.0 * km[1]).sin();
            let mut v = vec![0.0f64; 64];
            for s in 0..8usize {
                v[s + s * 8] = val;
            }
            v
        }
        _ => {
            let mut v = vec![0.0f64; 64];
            for s in 0..8usize {
                let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * km[1].cos();
                v[(s ^ 1) + s * 8] = cy;
            }
            v
        }
    }
}

/// χ_ch(q^lat) の被積分 (k 点, null 重み和) — v268p の構造
fn chi_integrand(ch: usize, k: [f64; 3], qs_lat: &[f64], w_null: &[f64], mutate: bool) -> f64 {
    let (wk, vk) = jacobi_eigh(&h8(k, 0.0), 8);
    let mut acc = 0.0f64;
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let (wq, vq) = jacobi_eigh(&h8(kq, 0.0), 8);
        let vv = channel(ch, k, q, mutate);
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
        acc += w_null[qi] * chi;
    }
    acc
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

fn null_weights(qs: &[f64; 4]) -> [f64; 4] {
    let mut a = [[0.0f64; 4]; 4];
    let mut rhs = [0.0f64; 4];
    for (i, &q) in qs.iter().enumerate() {
        a[0][i] = 1.0;
        a[1][i] = q * q;
        a[2][i] = q.powi(4);
        a[3][i] = q.powi(4) * (q * q).ln();
    }
    rhs[3] = 1.0;
    let mut m = a;
    for col in 0..4 {
        let piv = (col..4)
            .max_by(|&r1, &r2| m[r1][col].abs().partial_cmp(&m[r2][col].abs()).unwrap())
            .unwrap();
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for r in col + 1..4 {
            let f = m[r][col] / d;
            for c in col..4 {
                m[r][c] -= f * m[col][c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut w = [0.0f64; 4];
    for col in (0..4).rev() {
        let mut s = rhs[col];
        for c in col + 1..4 {
            s -= m[col][c] * w[c];
        }
        w[col] = s / m[col][col];
    }
    w
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
    println!("=== v27.0-B v270b_universality — 4D kernel の連続 universality (spec §13.2-B) ===\n");
    println!("予言 (P₂ 幾何, 静的 θ₀₀ = 1・θ_yy = θ₀y = 0): A_00/A_D = 2/3 厳密 /");
    println!("A_0y = A_yy = 0 (縦・時間 = 純ゲージ — Lean 定理の one-loop 版)。\n");
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
    let qset = [0.3f64, 0.6, 0.9, 1.2];
    let w_null = null_weights(&qset);
    let ladder = [0.35f64, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.032, 0.022];
    let names = ["D", "00", "0y", "yy"];

    // ---- 4 チャネルの A(a) ladder ----
    let run_ch = |ch: usize, a: f64, mutate: bool| -> f64 {
        let gl = gauss_legendre(14);
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let nodes = make_nodes(&nest_edges(PI / 2.0, 0.0, PI, a * qset[0]), &gl);
        bz_sum(&nodes, nthreads, |k| chi_integrand(ch, k, &qs_lat, &w_null, mutate))
            / a.powi(4)
    };
    let mut lads: Vec<Vec<(f64, f64)>> = vec![Vec::new(); 4];
    println!("    [A(a) 表] a | D | 00 | 0y | yy");
    for &a in &ladder {
        let mut row = String::new();
        for ch in 0..4 {
            let v = run_ch(ch, a, false);
            lads[ch].push((a, v));
            row = format!("{} {:+.4e} |", row, v);
        }
        println!("      a = {:.4}: {} ({} s)", a, row, t0.elapsed().as_secs());
    }

    // ---- [B0] 回帰: A_D rung a = 0.125 vs v268p 公表 ----
    {
        let aor2 = -2.0 / (160.0 * PI * PI);
        let got = lads[0][3].1 / aor2;
        let refv = 1.022196610668463;
        check(
            "[B0] 回帰: A_D (a = 0.125) が v268p 公表 JSON と 2e-3 相対一致",
            (got / refv - 1.0).abs() < 2e-3,
            format!("比 = {:.6} vs {:.6}", got, refv),
        );
    }

    // ---- 外挿 (§12.9 導出モデル) ----
    let b_stag = |a: f64| vec![1.0, a * a * (1.0 / a).ln(), a * a];
    let a0_full: Vec<f64> = (0..4).map(|ch| fit_a0(&lads[ch], &b_stag)).collect();
    let a0_tail: Vec<f64> = (0..4).map(|ch| fit_a0(&lads[ch][3..], &b_stag)).collect();
    println!("    [外挿] ch | 全域 | 尾部 (a ≤ 0.125)");
    for ch in 0..4 {
        println!(
            "      {}: {:+.4e} | {:+.4e}",
            names[ch], a0_full[ch], a0_tail[ch]
        );
    }

    // ---- [B1] A_00/(⅔ A_D) ----
    {
        let central = a0_tail[1] / (2.0 / 3.0 * a0_tail[0]);
        let full = a0_full[1] / (2.0 / 3.0 * a0_full[0]);
        let spread = (central - full).abs();
        check(
            "[B1] 温度チャネル universality: A_00/(⅔·A_D) = 1 ± 2% (系統 ≤ 1%)",
            (central - 1.0).abs() < 0.02 && spread < 0.01,
            format!("中心 = {:.4} (全域 {:.4}, spread {:.4})", central, full, spread),
        );
    }

    // ---- [B2][B3] 純ゲージチャネル ----
    {
        let r0y = (a0_tail[2] / a0_tail[0]).abs();
        let ryy = (a0_tail[3] / a0_tail[0]).abs();
        check(
            "[B2] |A_0y/A_D| < 0.05 (時間縦チャネルは純ゲージ)",
            r0y < 0.05,
            format!("|比| = {:.4}", r0y),
        );
        check(
            "[B3] |A_yy/A_D| < 0.05 (L = yy = 純ゲージ — Lean 定理の one-loop 版)",
            ryy < 0.05,
            format!("|比| = {:.4}", ryy),
        );
    }

    // ---- [B4] 変異 ----
    {
        // V₀₀ ×1.02 → A_00 ×1.0404 → B1 中心が +4% 逸脱
        let mut lad_mut = Vec::new();
        for &a in &ladder[3..] {
            lad_mut.push((a, run_ch(1, a, true)));
        }
        let a0m = fit_a0(&lad_mut, &b_stag);
        let central = a0m / (2.0 / 3.0 * a0_tail[0]);
        check(
            "[B4] 変異: V₀₀ ×1.02 → A_00/(⅔A_D) が +4% 逸脱 (> 2% ゲート外)",
            (central - 1.0).abs() > 0.03,
            format!("変異中心 = {:.4}", central),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v27.0-B".into())),
        ("kind".into(), Json::Str("kernel_4d_universality".into())),
        ("spec".into(), Json::Str("§13.2-B (b46bce3)".into())),
        ("a00_over_twothirds_ad".into(), Json::Num(a0_tail[1] / (2.0 / 3.0 * a0_tail[0]))),
        ("a0y_over_ad".into(), Json::Num(a0_tail[2] / a0_tail[0])),
        ("ayy_over_ad".into(), Json::Num(a0_tail[3] / a0_tail[0])),
    ]);
    let p = write_artifact("results/v270b_universality.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **温度セクター universality 成立 — A_00 = (2/3)A_D・縦/時間チャネルは純ゲージ (spec §13.3 条件 (ii) 充足)** — fork 判定 v27.0-C へ"
        } else {
            "FAIL あり — 分岐 (b) 温度 source 正規化 / P₂ 幾何の誤り (公表) / (c) 器械。欄が一次ソース"
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
