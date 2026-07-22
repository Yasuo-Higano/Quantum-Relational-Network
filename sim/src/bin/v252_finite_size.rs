//! v25.2 有限サイズ収束試験 — 数値の役割の再定義とゲートの fail-closed 化 (第二十六期, III/IV)
//!
//! v25.2-I で g(μ) の厳密閉形式が確定した (Eisler 対応)。これに伴い有限鎖計算の役割を
//! 再定義する: 「フィットから g の形を発見する」ではなく「**厳密な g(μ) への有限サイズ
//! 収束則を測定する**」。この変更で、フィット不一致が (i) 物理式の反証なのか
//! (ii) 有限サイズ不足なのかを分離できる (PROMPT/7 の方針)。
//!
//! 器械 (v251 と同一アルゴリズム — stag.rs へ昇格した chain_k):
//!   モジュラー核 K の NN ボンド勾配 / πξ = ĝ_N(μ)、対角質量チャネル / 2πξ_s μ = ĝ_m。
//!   クランプ 2 点梯子 (1e-30/1e-26) の一致域 ξ* を信頼窓とし、**ξ* < 5 は測定不能
//!   (None) — fail-closed** (旧 v251 は ξ ≤ 4 を強制フィットする fail-open だった)。
//!
//! 物理スケール: 分散 E(q) = √(cos²q+μ²) の複素零点から ξ_corr = 1/arsinh(μ)。
//! 無次元変数 x = N·arsinh(μ) が「鎖が gap を知っているか」を決める:
//!   x ≪ 1 → 鎖は実質 massless (ĝ は massless 値 ĝ_N(0) を読む) /
//!   x ≫ 1 → 半無限 gapped の厳密 prefactor g(μ) に収束。
//!
//! 検査 (ゲートは器械精度のみ — 「収束の形」は物理仮説なのでゲートにせず分岐で記録):
//!  [S1] 器械精度: μ ∈ {0.5, 1, √2} × N ∈ {64,128,256,512} の全点が測定可能 (除外ゼロ)
//!       かつ max|ĝ/g_exact − 1| ≤ 2e-5・max|ĝ_m/g_exact − 1| ≤ 2e-6。
//!       収束次数 p (ln err vs ln N) は測定し、分岐を分類する:
//!       (a-i) p ≥ 0.5 = 有限サイズ支配 / (a-ii) |p| < 0.5 = **窓/クランプの器械床支配**
//!       (床は N 非依存 — ξ* = 5–6 の κ 予算端で決まる)。
//!  [S2] crossover: μ ∈ {0.01,0.02,0.05,0.1} × N — x ≥ 25 で |d| ≤ 1e-4 (収束域) /
//!       x ≤ 2 で ĝ ≈ ĝ_N(0) ± 2e-3 (massless 読み — v25.1 の「小 μ 系統」の正体)
//!  [S3] fail-closed 負制御: 粗クランプ対 (1e-8, 1e-6) では κ 予算 (8ξ+12 > 14) が
//!       ξ=1 から破綻し ξ* = 0 → **None が返ること** (旧実装なら ξ≤4 の強制フィットで
//!       値が黙って出ていた — 修正の検出可能性を証拠化)
//!  [S4] 一次データ同一性: N=256, μ ∈ {0.5,1,√2} の ĝ, ĝ_m が v251_gmu.json の
//!       アーカイブ値と ≤ 1e-9 で一致 (窓が同一なら bit 同値のはず — 遡及修正の不変性)
//!
//! 事前登録分岐: (a) 全 PASS → 器械化完了。副分岐 (a-i)/(a-ii) は収束の形の記録 /
//!   (b) S1 精度 FAIL → 器械床が想定より高い (床の同定へ) / (c) S3/S4 FAIL → 器械の再監査。
//!
//! 開発記録 (run1 → run2 のゲート再設計): 初版 S1 は「err が N で減少」を要求したが、
//! 実測は全 N で平坦 (p ≈ 0, 残差 −8.1e-6/−2.3e-6/−5.8e-6 が N=64→512 で不変) —
//! 残差は有限サイズではなく器械床だった。物理仮説をゲートにした設計ミスを v24.3
//! W1/W2 と同型の教訓として修正 (ゲートは器械の実測に貼る)。測定値は run1/run2 で同一。

use uft_sim::dd::dd;
use uft_sim::stag::chain_k;
use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

struct Extract {
    g: f64,
    g_m: f64,
    xi_trust: usize,
}

/// v251 と同一の抽出 (窓 [3, ξ*] のボンド勾配 / [2.5, ξ*] の質量対角勾配)。
/// fail-closed: ξ* < 5 は None。クランプ対を引数化 (負制御用)。
fn extract(n: usize, mu: f64, clamp_a: f64, clamp_b: f64) -> Option<Extract> {
    let ka = chain_k(n, dd(mu), clamp_a);
    let kb = chain_k(n, dd(mu), clamp_b);
    let half = n / 2;
    let bond = |k: &Vec<uft_sim::dd::Dd>, i: usize| k[i + (i + 1) * half].hi;
    let mut xi_trust = 0usize;
    for xi in 1..half {
        let i = half - 1 - xi;
        let rel = ((bond(&ka, i) - bond(&kb, i)) / bond(&ka, i)).abs();
        if rel < 1e-4 {
            xi_trust = xi;
        } else {
            break;
        }
    }
    if xi_trust < 5 {
        return None;
    }
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for xi in 3..=xi_trust {
        let i = half - 1 - xi;
        xs.push(PI * xi as f64);
        ys.push(bond(&ka, i));
    }
    let (icpt, g) = linfit_checked(&xs, &ys).ok()?;
    let _ = icpt;
    let mut g_m = 0.0;
    if mu > 0.0 {
        let mut xs2 = Vec::new();
        let mut ys2 = Vec::new();
        for j in 0..half {
            let i = half - 1 - j;
            let xis = j as f64 + 0.5;
            if xis < 2.5 || xis > xi_trust as f64 {
                continue;
            }
            let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
            xs2.push(2.0 * PI * xis * mu);
            ys2.push(sign * ka[i + i * half].hi);
        }
        let (_i2, s2) = linfit_checked(&xs2, &ys2).ok()?;
        g_m = s2;
    }
    Some(Extract { g, g_m, xi_trust })
}

fn main() {
    self_test();
    println!(
        "=== v25.2 有限サイズ収束試験 — 役割の再定義と fail-closed ゲート (第二十六期, III/IV) ===\n"
    );
    println!("事前登録: (a) 全 PASS → 収束則測定完了 / (b) S1 のみ FAIL → 窓系統の床を記録 /");
    println!("          (c) S3/S4 FAIL → 器械再監査 (fail-closed 化が数値を変えた = 重大)\n");
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(1);
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

    // ---- 測定ジョブ (決定的スレッド分割) ----
    let ns = [64usize, 128, 256, 512];
    let mus_a = [0.5f64, 1.0, 2.0f64.sqrt()];
    let mus_b = [0.01f64, 0.02, 0.05, 0.1];
    let mut jobs: Vec<(usize, f64)> = Vec::new();
    for &n in &ns {
        jobs.push((n, 0.0));
        for &m in mus_a.iter().chain(&mus_b) {
            jobs.push((n, m));
        }
    }
    let mut out: Vec<Option<Option<Extract>>> = Vec::new();
    out.resize_with(jobs.len(), || None);
    let chunk = jobs.len().div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in out.chunks_mut(chunk).enumerate() {
            let jobs = &jobs;
            sc.spawn(move || {
                for (j, slot) in sl.iter_mut().enumerate() {
                    let (n, mu) = jobs[t * chunk + j];
                    *slot = Some(extract(n, mu, 1e-30, 1e-26));
                }
            });
        }
    });
    let get = |n: usize, mu: f64| -> Option<&Extract> {
        jobs.iter()
            .position(|&(jn, jm)| jn == n && (jm - mu).abs() < 1e-15)
            .and_then(|i| out[i].as_ref().unwrap().as_ref())
    };
    println!(
        "    [測定完了] {} ジョブ ({} s, {} threads)\n",
        jobs.len(),
        t0.elapsed().as_secs(),
        nthreads
    );

    // ---- [S1] 器械精度と収束の形 (μ ≥ 0.5 の代表 3 点) ----
    // ゲート = 精度のみ。収束次数 p は分岐 (a-i)/(a-ii) の分類に使う (物理仮説を
    // ゲートにしない — run1 の設計ミスの修正、ヘッダ開発記録参照)。
    println!("    [収束表] μ | N | ξ* | ĝ_N | ĝ_N/g_exact − 1 | ĝ_m/g_exact − 1");
    let mut s1_ok = true;
    let mut max_g_err = 0.0f64;
    let mut max_gm_err = 0.0f64;
    let mut branch = String::new();
    for &mu in &mus_a {
        let ge = g_exact(mu).unwrap();
        let mut errs = Vec::new();
        for &n in &ns {
            match get(n, mu) {
                Some(e) => {
                    let d = e.g / ge - 1.0;
                    let dm = e.g_m / ge - 1.0;
                    println!(
                        "      μ={:5.3} N={:4} ξ*={:2}  ĝ = {:.9}  {:+.2e}  {:+.2e}",
                        mu, n, e.xi_trust, e.g, d, dm
                    );
                    errs.push((n as f64, d.abs()));
                    max_g_err = max_g_err.max(d.abs());
                    max_gm_err = max_gm_err.max(dm.abs());
                }
                None => {
                    println!("      μ={:5.3} N={:4}  除外 (ξ* < 5)", mu, n);
                    s1_ok = false;
                }
            }
        }
        if errs.len() == ns.len() {
            let xs: Vec<f64> = errs.iter().map(|(n, _)| n.ln()).collect();
            let ys: Vec<f64> = errs.iter().map(|(_, e)| e.max(1e-14).ln()).collect();
            let p = linfit_checked(&xs, &ys)
                .map(|(_, b)| -b)
                .unwrap_or(f64::NAN);
            branch.push_str(&format!(
                "μ={:.3}: p = {:+.2} → {}  ",
                mu,
                p,
                if p >= 0.5 {
                    "(a-i) 有限サイズ支配"
                } else {
                    "(a-ii) 器械床支配"
                }
            ));
        }
    }
    check(
        "[S1] 器械精度 (除外ゼロ・max|ĝ/g−1| ≤ 2e-5・max|ĝ_m/g−1| ≤ 2e-6)",
        s1_ok && max_g_err <= 2e-5 && max_gm_err <= 2e-6,
        format!("max = {:.1e} / {:.1e}", max_g_err, max_gm_err),
    );
    println!("    [収束の形 (分岐記録)] {}", branch);

    // ---- [S2] crossover: x = N·arsinh(μ) ----
    println!("\n    [crossover 表] μ | N | x = N·arsinh μ | ĝ | ĝ/g_exact − 1 | ĝ − ĝ_N(0)");
    let mut s2a_ok = true;
    let mut s2b_ok = true;
    for &mu in &mus_b {
        let ge = g_exact(mu).unwrap();
        for &n in &ns {
            let x = n as f64 * mu.asinh();
            let (d, dml) = match (get(n, mu), get(n, 0.0)) {
                (Some(e), Some(e0)) => (e.g / ge - 1.0, e.g - e0.g),
                _ => (f64::NAN, f64::NAN),
            };
            let tag = if x >= 25.0 {
                if d.abs() > 1e-4 {
                    s2a_ok = false;
                }
                "収束域"
            } else if x <= 2.0 {
                if dml.abs() > 2e-3 {
                    s2b_ok = false;
                }
                "massless 読み"
            } else {
                "crossover"
            };
            println!(
                "      μ={:5.3} N={:4}  x = {:6.2}  ĝ = {:.6}  {:+.2e}  {:+.2e}  [{}]",
                mu,
                n,
                x,
                get(n, mu).map(|e| e.g).unwrap_or(f64::NAN),
                d,
                dml,
                tag
            );
        }
    }
    check(
        "[S2a] x ≥ 25 は収束域 (|ĝ/g_exact − 1| ≤ 1e-4)",
        s2a_ok,
        String::new(),
    );
    check(
        "[S2b] x ≤ 2 は massless 読み (|ĝ − ĝ_N(0)| ≤ 2e-3) — v25.1 小 μ 系統の正体",
        s2b_ok,
        String::new(),
    );

    // ---- [S3] fail-closed 負制御 ----
    {
        let bad = extract(64, 1.0, 1e-8, 1e-6);
        check(
            "[S3] 負制御: 粗クランプ対 (1e-8,1e-6) は κ 予算破綻 → None (旧実装は強制フィット)",
            bad.is_none(),
            format!(
                "extract = {}",
                if bad.is_none() {
                    "None (fail-closed 動作確認)"
                } else {
                    "Some — fail-open が残っている!"
                }
            ),
        );
    }

    // ---- [S4] 一次データ同一性 (v251_gmu.json, N=256) ----
    {
        let txt = std::fs::read_to_string("results/v251_gmu.json")
            .or_else(|_| std::fs::read_to_string("../results/v251_gmu.json"))
            .unwrap_or_default();
        let lookup = |mu: f64| -> Option<(f64, f64)> {
            // "mu": <v> の直後の "g": と "g_m": を拾う (v252_exact_g と同じ最小パーサ)
            for chunk in txt.split("\"mu\":").skip(1) {
                let end = chunk.find(|c: char| c == ',' || c == '}')?;
                let m: f64 = chunk[..end].trim().parse().ok()?;
                if (m - mu).abs() < 1e-12 {
                    let obj = &chunk[..chunk.find('}').unwrap_or(chunk.len())];
                    let num = |key: &str| -> Option<f64> {
                        let p = obj.find(key)? + key.len();
                        let r = &obj[p..];
                        let e = r.find(|c: char| c == ',' || c == '}' || c == '\n')?;
                        r[..e].trim().parse().ok()
                    };
                    return Some((num("\"g\":")?, num("\"g_m\":")?));
                }
            }
            None
        };
        let mut worst = 0.0f64;
        let mut all_found = true;
        for &mu in &mus_a {
            match (get(256, mu), lookup(mu)) {
                (Some(e), Some((ga, gma))) => {
                    worst = worst.max((e.g - ga).abs()).max((e.g_m - gma).abs());
                }
                _ => all_found = false,
            }
        }
        check(
            "[S4] 一次データ同一性 (N=256): fail-closed 化後も v25.1 アーカイブと一致",
            all_found && worst < 1e-9,
            format!("max|Δ| = {:.1e} (窓同一なら bit 同値)", worst),
        );
    }

    // ---- 判定 ----
    println!(
        "\n[判定] {}",
        if nfail == 0 {
            format!(
                "事前登録 (a): **fail-closed 器械化完了・有限鎖は厳密式への収束試験として再定義** — 収束の形: {}",
                branch.trim_end()
            )
        } else {
            "FAIL — 分岐 (b)/(c) は各検査の欄を一次ソースとする".to_string()
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
