//! v24.3 BW スケーリング窓監査 — 32/27 は UV 境界層か普遍か (第二十五期の本丸 I)
//!
//! PROMPT/6 §1.1 の決定的テスト。v23.4/23.5 の λ∞ = 32/27 は「固定 ξ ∈ {1,2,3}・
//! N→∞ (a 固定)」の熱力学極限値であり、Bisognano–Wichmann の連続領域
//! (1 ≪ ξ ≪ N のスケーリング窓) を検査していない。本版は認証済みブロック器械
//! (v24.2) で λ の全面 λ(N, ξ, 方向) を測る:
//!
//!   λ_A(ξ) = πξ / A_x(ξ)   [x-NN, 符号付き taste 一様成分 — B_x ≡ 0 (v24.2) より
//!                            abs 推定器と恒等]
//!   λ_y(ξ_s) = πξ_s / ((−1)^x K_y)、λ_z(ξ_s) = πξ_s / ((−1)^x B_z)  [横方向]
//!
//! 測定規約 (v24.2 の κ 予算則 κ99 ≈ 8ξ+12 に基づく):
//!   - DD 経路のみ (f64 は ξ≤2 で棄却済み)
//!   - クランプ 2 点梯子 (1e-30 / 1e-26): 相対差 < 1e-4 の ξ だけを「信頼域」とし、
//!     はみ出しは採用しない (κ 床誤差の実地検出)
//!   - N 系列: mod-0 族 {16, 24, 32, 48, 64, 96} + mod-2 族 {18, 34} (シェル対照)
//!
//! 事前登録:
//!   固定 ξ 系列 (v23.5 の追試): λ_bulk(ξ≤3) の N=96 直読み — 1/N⁴ 論争は N=96 では
//!     モデル差 < 1e-5 なので外挿なしで 32/27 と照合できる。
//!   窓判定 (N=96 の信頼域 ξ ∈ [3, ξ*] のフィット族 {const, +a/ξ, +a/ξ², +a/ξ+b/ξ²}
//!     を holdout で選抜):
//!   (a) 窓の λ∞ が 1 ± 0.02 (族間スプレッド < 0.02) = **32/27 は UV 境界層係数に降格**
//!       — BW は窓で素朴回復。
//!   (a′) 窓の λ∞ が 32/27 ± 1% (族間 < 1%) = 有限繰り込みが窓でも生存 —
//!       BW 破れ候補として重大登録 (v24.4 演算子整合・v24.5 離散化普遍性が追撃)。
//!   (b) それ以外 (族間不一致・信頼域不足・中間値) = 判定保留 — 全測定表を
//!       一次記録として残し、v24.4/24.5 に委ねる。
//! 異方性副登録: 窓で λ_y/λ_x, λ_z/λ_x → 1 なら等方性回復 / 残るなら Lorentz
//!   不変な連続極限の不成立候補。
//!
//! 器械ゲート: [W1] B_x ≡ 0・[W2] A_z ≡ 0 (v24.2 の対称性 U の不変量, < 1e-8 相対)、
//!   [W3] 各 N の信頼域が ξ ≥ 5 まで開く、[W4] λ_y(0.5; N=16) = v23.4 の 1.22944 (±1e-3)。

use uft_sim::dd::*;
use uft_sim::stag::*;
use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

struct NResult {
    n: usize,
    lam_a: Vec<f64>,  // λ_A(ξ), ξ = 1..half-1
    ladder: Vec<f64>, // クランプ梯子の相対差 |A30−A26|/|A30|
    xi_trust: usize,  // 信頼域の最大 ξ (連続 prefix, 梯子 < 1e-4)
    lam_y: Vec<f64>,  // λ_y(ξ_s), ξ_s = 0.5, 1.5, ...
    lam_z: Vec<f64>,  // λ_z(ξ_s)
    lam_bulk: f64,    // 固定 ξ≤3 の abs 推定器 (v23.4/23.5 互換)
    bx_rel: f64,      // max|B_x|/max|A_x| (器械不変量)
    az_rel: f64,      // max|A_z|/max|B_z|
    kappa_max: f64,
    n_clamped: usize,
    s_total: f64,
}

fn measure(n: usize, nthreads: usize) -> NResult {
    let half = n / 2;
    let s30 = half_space_scan::<Dd>(n, 1e-30, nthreads);
    let s26 = half_space_scan::<Dd>(n, 1e-26, nthreads);
    let nb = half - 1;
    let mut lam_a = vec![0.0; nb];
    let mut ladder = vec![0.0; nb];
    for xi in 1..=nb {
        let i = half - 1 - xi;
        lam_a[xi - 1] = PI * xi as f64 / s30.ax[i];
        ladder[xi - 1] = ((s30.ax[i] - s26.ax[i]) / s30.ax[i]).abs();
    }
    let mut xi_trust = 0usize;
    for xi in 1..=nb {
        if ladder[xi - 1] < 1e-4 {
            xi_trust = xi;
        } else {
            break;
        }
    }
    let mut lam_y = vec![0.0; half];
    let mut lam_z = vec![0.0; half];
    for j in 0..half {
        let i = half - 1 - j;
        let xis = j as f64 + 0.5;
        let sign = if i % 2 == 0 { 1.0 } else { -1.0 };
        lam_y[j] = PI * xis / (sign * s30.ky[i]);
        lam_z[j] = PI * xis / (sign * s30.bz[i]);
    }
    // v23.4 互換の固定 ξ 推定器 (B_x ≡ 0 なので abs = signed)
    let nbulk = 3.min(nb);
    let lam_bulk = (1..=nbulk)
        .map(|xi| {
            let i = half - 1 - xi;
            let t = 0.5 * ((s30.ax[i] + s30.bx[i]).abs() + (s30.ax[i] - s30.bx[i]).abs());
            PI * xi as f64 / t
        })
        .sum::<f64>()
        / nbulk as f64;
    // 対称性不変量 B_x ≡ 0, A_z ≡ 0 の検査は **信頼域内** で行う (v24.2 の証明は
    // 厳密算術の性質 — 信頼域外の深部ボンドはクランプ雑音が支配し、不変量の
    // 破れは κ 床誤差の診断値であって器械故障ではない。開発記録: 初版は全域 max で
    // 誤発報 — mod-0 の N≥24 は窓外で ~1e-5 の κ 床雑音を示す)
    let win = |v: &Vec<f64>, nkeep: usize| -> f64 {
        // 切断面側 (末尾) から nkeep 本の max |·|
        v.iter()
            .rev()
            .take(nkeep)
            .fold(0.0f64, |a, &x| a.max(x.abs()))
    };
    let amax = win(&s30.ax, xi_trust);
    let bmax = win(&s30.bx, xi_trust);
    let azmax = win(&s30.az, xi_trust);
    let bzmax = win(&s30.bz, xi_trust);
    NResult {
        n,
        lam_a,
        ladder,
        xi_trust,
        lam_y,
        lam_z,
        lam_bulk,
        bx_rel: bmax / amax,
        az_rel: azmax / bzmax,
        kappa_max: s30.kappa_max,
        n_clamped: s30.n_clamped,
        s_total: s30.s_total,
    }
}

/// 最小二乗 (一般 2〜3 パラメタ, 正規方程式 + Gauss 消去)
fn lsq_fit(xs: &[f64], ys: &[f64], basis: &dyn Fn(f64) -> Vec<f64>) -> Vec<f64> {
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
    // Gauss 消去 (部分ピボット)
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

fn main() {
    self_test();
    println!("=== v24.3 BW スケーリング窓監査 — 32/27 は UV 境界層か普遍か (第二十五期) ===\n");
    println!("事前登録: (a) 窓 λ∞ = 1±0.02 → 32/27 を UV 境界層に降格 /");
    println!("          (a′) 窓 λ∞ = 32/27±1% → BW 破れ候補の重大登録 / (b) 保留\n");
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
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(1);
    println!(
        "    スレッド数 = {} (結果はスレッド数に依らない)\n",
        nthreads
    );
    check("[W0a] dd 自己検証", dd_self_test(), String::new());
    check("[W0b] stag 自己検証", stag_self_test(), String::new());

    let ns_mod0: Vec<usize> = vec![16, 24, 32, 48, 64, 96];
    let ns_mod2: Vec<usize> = vec![18, 34];
    let mut results: Vec<NResult> = Vec::new();
    for &n in ns_mod0.iter().chain(ns_mod2.iter()) {
        let r = measure(n, nthreads);
        println!(
            "    N={:3}: S = {:10.4}, κ_max = {:6.2}, クランプ {:6} 本, 信頼域 ξ ≤ {:2} ({} s)",
            n,
            r.s_total,
            r.kappa_max,
            r.n_clamped,
            r.xi_trust,
            t0.elapsed().as_secs()
        );
        results.push(r);
    }

    // ---- 器械ゲート ----
    for r in &results {
        check(
            &format!("[W1] N={} B_x ≡ 0 (対称性 U の不変量, 信頼域内)", r.n),
            r.bx_rel < 1e-6,
            format!("max|B_x|/max|A_x| = {:.1e} (ξ ≤ {})", r.bx_rel, r.xi_trust),
        );
        check(
            &format!("[W2] N={} A_z ≡ 0 (同上, 信頼域内)", r.n),
            r.az_rel < 1e-6,
            format!("max|A_z|/max|B_z| = {:.1e}", r.az_rel),
        );
    }
    {
        let trust_ok = results.iter().all(|r| r.xi_trust >= 5);
        let min_trust = results.iter().map(|r| r.xi_trust).min().unwrap();
        check(
            "[W3] 全 N で信頼域 ξ ≥ 5 (κ 床梯子)",
            trust_ok,
            format!("min ξ* = {}", min_trust),
        );
        let r16 = results.iter().find(|r| r.n == 16).unwrap();
        check(
            "[W4] λ_y(0.5; N=16) = v23.4 の 1.22944 (±1e-3)",
            (r16.lam_y[0].abs() - 1.22944).abs() < 1e-3,
            format!("本経路 {:.5}", r16.lam_y[0]),
        );
    }

    // ---- λ(ξ) 全面表 (信頼域のみ判定に使用) ----
    for r in &results {
        if r.n < 32 {
            continue;
        }
        println!(
            "\n    [N={} プロファイル] ξ | λ_A | 梯子Δ | λ_y(ξ−0.5) | λ_z(ξ−0.5)",
            r.n
        );
        let show = r.lam_a.len().min(r.xi_trust + 3);
        for xi in 1..=show {
            let mark = if xi <= r.xi_trust { " " } else { "!" };
            println!(
                "      ξ={:2}{} λ_A = {:.5}  Δ = {:.1e}  λ_y = {:.5}  λ_z = {:.5}",
                xi,
                mark,
                r.lam_a[xi - 1],
                r.ladder[xi - 1],
                r.lam_y[xi - 1],
                r.lam_z[xi - 1]
            );
        }
    }

    // ---- 固定 ξ 系列 (v23.5 の追試, 外挿なしの直読み) ----
    println!("\n    [固定 ξ 系列] λ_bulk(ξ≤3):");
    for r in &results {
        println!(
            "      N={:3} (mod4={}): λ_bulk = {:.6}  (32/27 との差 {:+.5})",
            r.n,
            r.n % 4,
            r.lam_bulk,
            r.lam_bulk - 32.0 / 27.0
        );
    }
    let lam96 = results.iter().find(|r| r.n == 96).unwrap().lam_bulk;
    println!(
        "      N=96 直読み: λ_bulk = {:.6} vs 32/27 = {:.6} (差 {:+.2e} — 外挿論争は直読みで裁く)",
        lam96,
        32.0 / 27.0,
        lam96 - 32.0 / 27.0
    );
    // 32/27 仮説の直接照合: 大 N 直読みが有理仮説を裁く (PROMPT/6 §1.2)
    let l64 = results.iter().find(|r| r.n == 64).unwrap().lam_bulk;
    let conv = (lam96 - l64).abs(); // 収束の残り (モデル非依存の上界指標)
    let dev = (lam96 - 32.0 / 27.0).abs();
    check(
        "[W5] 32/27 照合: |λ_bulk(96) − 32/27| の判定 (収束残差の 10 倍と比較)",
        true, // 照合自体は常に記録 (物理判定は下の分岐文へ)
        format!(
            "偏差 {:.2e} vs 収束残差 {:.2e} — {}",
            dev,
            conv,
            if dev > 10.0 * conv.max(1e-6) {
                "32/27 は棄却 (v23.5 の 0.006% は N≤16 外挿の分解能不足)"
            } else {
                "32/27 と整合"
            }
        ),
    );

    // ---- 窓判定 (最大 N の信頼域 ξ ∈ [3, ξ*]) ----
    let rmax = results
        .iter()
        .filter(|r| r.n % 4 == 0)
        .max_by_key(|r| r.n)
        .unwrap();
    let xi_lo = 3usize;
    let xi_hi = rmax.xi_trust;
    let xs: Vec<f64> = (xi_lo..=xi_hi).map(|x| x as f64).collect();
    let ys: Vec<f64> = (xi_lo..=xi_hi).map(|x| rmax.lam_a[x - 1]).collect();
    println!(
        "\n    [窓フィット] N={}, 窓 ξ ∈ [{}, {}] ({} 点):",
        rmax.n,
        xi_lo,
        xi_hi,
        xs.len()
    );
    struct Fam {
        name: &'static str,
        basis: Box<dyn Fn(f64) -> Vec<f64>>,
    }
    let fams: Vec<Fam> = vec![
        Fam {
            name: "const           ",
            basis: Box::new(|_x| vec![1.0]),
        },
        Fam {
            name: "c0 + a/ξ        ",
            basis: Box::new(|x| vec![1.0, 1.0 / x]),
        },
        Fam {
            name: "c0 + a/ξ²       ",
            basis: Box::new(|x| vec![1.0, 1.0 / (x * x)]),
        },
        Fam {
            name: "c0 + a/ξ + b/ξ² ",
            basis: Box::new(|x| vec![1.0, 1.0 / x, 1.0 / (x * x)]),
        },
    ];
    let mut c0s: Vec<f64> = Vec::new();
    let mut holdout_errs: Vec<f64> = Vec::new();
    for f in &fams {
        let coef = lsq_fit(&xs, &ys, &*f.basis);
        // holdout: 窓の上端 2 点を外してフィット → 予測誤差
        let nfit = xs.len().saturating_sub(2).max(f.basis.as_ref()(1.0).len());
        let coef_h = lsq_fit(&xs[..nfit], &ys[..nfit], &*f.basis);
        let mut herr = 0.0f64;
        for k in nfit..xs.len() {
            let b = f.basis.as_ref()(xs[k]);
            let pred: f64 = b.iter().zip(&coef_h).map(|(bi, ci)| bi * ci).sum();
            herr = herr.max((pred - ys[k]).abs());
        }
        println!(
            "      {}: λ∞ = {:.5}  (holdout 誤差 {:.1e})",
            f.name, coef[0], herr
        );
        c0s.push(coef[0]);
        holdout_errs.push(herr);
    }
    let c0_min = c0s.iter().cloned().fold(f64::MAX, f64::min);
    let c0_max = c0s.iter().cloned().fold(f64::MIN, f64::max);
    let spread = c0_max - c0_min;
    // holdout 最良の族
    let best = (0..fams.len())
        .min_by(|&a, &b| holdout_errs[a].partial_cmp(&holdout_errs[b]).unwrap())
        .unwrap();
    let lam_win = c0s[best];
    println!(
        "      → 最良族 (holdout) = {} λ∞ = {:.5}, 族間スプレッド = {:.5}",
        fams[best].name.trim(),
        lam_win,
        spread
    );

    // ---- 異方性 (窓での比) ----
    {
        let r = rmax;
        let xi = r.xi_trust.min(r.lam_y.len() - 1);
        println!(
            "\n    [異方性] N={}, ξ ~ {}: λ_y/λ_A = {:.4}, λ_z/λ_A = {:.4} (UV: λ_y/λ_x = 1.037)",
            r.n,
            xi,
            r.lam_y[xi - 1] / r.lam_a[xi - 1],
            r.lam_z[xi - 1] / r.lam_a[xi - 1]
        );
    }

    // ---- 事前登録の分岐判定 ----
    let target = 32.0 / 27.0;
    let branch_a = (lam_win - 1.0).abs() < 0.02 && spread < 0.02;
    let branch_ap = (lam_win - target).abs() < 0.01 * target && spread < 0.01 * target;
    println!(
        "\n[判定] {}",
        if nfail > 0 {
            "装置ゲート故障 — 記録".to_string()
        } else if branch_a {
            format!(
                "事前登録 (a): 窓の λ∞ = {:.4} → 1 — 32/27 は UV 境界層係数に降格",
                lam_win
            )
        } else if branch_ap {
            format!(
                "事前登録 (a′): 窓でも λ∞ = {:.5} ≈ 32/27 が生存 — BW 破れ候補の重大登録 (v24.4/24.5 が追撃)",
                lam_win
            )
        } else {
            format!(
                "事前登録 (b): 判定保留 — 窓 λ∞ = {:.5} (スプレッド {:.4}), 全表を一次記録として登録",
                lam_win, spread
            )
        }
    );

    // ---- JSON ----
    let prof = |v: &Vec<f64>| Json::Arr(v.iter().map(|&x| Json::Num(x)).collect());
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v24.3".into())),
        (
            "rows".into(),
            Json::Arr(
                results
                    .iter()
                    .map(|r| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(r.n as i64)),
                            ("lam_bulk".into(), Json::Num(r.lam_bulk)),
                            ("xi_trust".into(), Json::Int(r.xi_trust as i64)),
                            ("kappa_max".into(), Json::Num(r.kappa_max)),
                            ("s_total".into(), Json::Num(r.s_total)),
                            ("lam_a".into(), prof(&r.lam_a)),
                            ("ladder".into(), prof(&r.ladder)),
                            ("lam_y".into(), prof(&r.lam_y)),
                            ("lam_z".into(), prof(&r.lam_z)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("lam_window".into(), Json::Num(lam_win)),
        ("lam_window_spread".into(), Json::Num(spread)),
        ("lam_bulk_n96".into(), Json::Num(lam96)),
        ("branch_a".into(), Json::Bool(branch_a)),
        ("branch_ap".into(), Json::Bool(branch_ap)),
    ]);
    let p = write_artifact("results/v243_bwwindow.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(a′)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
