//! v24.6 エントロピー普遍性 — 点ノード面積則の正しい有限サイズ形と G^bare の再分類
//!
//! PROMPT/6 §1.3–1.4。v23.6 の「対数 vs 定数」二択は Fermi 面 (余次元 1) 前提の
//! Gioev–Klich 形 (S ~ N²lnN) を対照にしていたが、本模型の零集合は点ノード
//! (余次元 3, v24.1) — 正しい候補族は S = αN² + 抑制補正。本版は分散から予告される
//! 族を**事前登録**し、holdout N で裁く:
//!
//!   F1: s = α                    (定数)
//!   F2: s = α + c/N²             (点ノード予想の主形)
//!   F3: s = α + (β·lnN + γ)/N²   (曲面対数の平面残滓の検査 — β の要不要)
//!   F4: s = α + c/N              (境界項対照 — 周期横方向では不要のはず)
//!   F5: s = α + β·lnN            (Gioev–Klich 対照 — 点ノードでは誤形のはず)
//!
//!   フィット: mod-0 族 N ∈ {16,24,32,48,64} / holdout: N ∈ {96,128,192}
//!   (シェル系統は mod-2 族 {18,34,50,66} を並記)
//!
//! 追加測定:
//!   - UV 有限量: 相互情報 I(A:B) (x スラブ対, 距離 s) — 将来の離散化間比較の
//!     基準値として登録 (leading area 発散が相殺する連続比較用観測量)
//!   - G^bare の再分類: G = a²/(4·λ·s∞) は λ の方向依存 (v24.4) により
//!     一意でない — G_entropy^bare(a; regulator, 方向スキーム) と正しく命名し、
//!     x/y 両スキームの値域を記録 (v23.6 の 1.06a² の再解釈)
//!
//! 器械ゲート: [E1] ブロック経路の S(N=16)/N² が v23.6 公表値と一致 (歴史照合)。
//! 事前登録: (a) F5 が holdout 最悪 & F3 の β 不要 (F2 で十分) = 定数面積則 +
//!   1/N² 補正の確定 — v23.6 の二択は前提から誤りだったが結論 (定数則) は生存 /
//!   (b) β·lnN/N² が holdout で必要 = 平面対数残滓の発見 / (c) その他 = 記録。

use uft_sim::dd::*;
use uft_sim::stag::*;
use uft_sim::*;

fn json_num(s: &str, key: &str, occurrence: usize) -> f64 {
    let pat = format!("\"{}\":", key);
    let mut idx = 0usize;
    let mut from = 0usize;
    loop {
        let i = s[from..]
            .find(&pat)
            .unwrap_or_else(|| panic!("json key {}", key))
            + from;
        if idx == occurrence {
            let rest = &s[i + pat.len()..];
            let end = rest
                .find(|c| c == ',' || c == '}' || c == '\n')
                .unwrap_or(rest.len());
            return rest[..end].trim().parse().expect("parse");
        }
        idx += 1;
        from = i + pat.len();
    }
}

/// 最小二乗 (正規方程式 + Gauss 消去)
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

fn main() {
    self_test();
    println!("=== v24.6 エントロピー普遍性 — 点ノード面積則と G^bare の再分類 (第二十五期) ===\n");
    println!("事前登録: (a) F5 (Gioev–Klich) が holdout 最悪 & F2 で十分 → 定数面積則 + 1/N² /");
    println!("          (b) β·lnN/N² が必要 → 平面対数残滓 / (c) その他 = 記録\n");
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
    check("[E0a] dd 自己検証", dd_self_test(), String::new());
    check("[E0b] stag 自己検証", stag_self_test(), String::new());

    // ---- s_area(N) 系列 (f64 — S はクランプ鈍感) ----
    let ns_fit: Vec<usize> = vec![16, 24, 32, 48, 64];
    let ns_hold: Vec<usize> = vec![96, 128, 192];
    let ns_mod2: Vec<usize> = vec![18, 34, 50, 66];
    let mut s_area = std::collections::BTreeMap::<usize, f64>::new();
    for &n in ns_fit.iter().chain(ns_hold.iter()).chain(ns_mod2.iter()) {
        let xsel: Vec<usize> = (0..n / 2).collect();
        let s = region_entropy::<f64>(n, &xsel, nthreads);
        let sa = s / ((n * n) as f64);
        s_area.insert(n, sa);
        println!(
            "    N={:3} (mod4={}): S = {:12.4}, s_area = {:.6} ({} s)",
            n,
            n % 4,
            s,
            sa,
            t0.elapsed().as_secs()
        );
    }

    // ---- [E1] 歴史照合 (v23.6 の s_area(N=16)) ----
    {
        let js = std::fs::read_to_string("results/v236_gread.json").expect("v236 json");
        // rows は N 昇順 (8,10,12,14,16) — s_area の 5 番目 (N=16)
        let s16_v236 = json_num(&js, "s_area", 4);
        let d = (s_area[&16] - s16_v236).abs();
        check(
            "[E1] s_area(16) = v23.6 公表値 (歴史照合)",
            d < 1e-6,
            format!(
                "ブロック {:.6} vs v23.6 {:.6} (Δ {:.1e})",
                s_area[&16], s16_v236, d
            ),
        );
    }

    // ---- [E2] 事前登録フィット族 + holdout ----
    struct Fam {
        name: &'static str,
        basis: Box<dyn Fn(f64) -> Vec<f64>>,
    }
    let fams: Vec<Fam> = vec![
        Fam {
            name: "F1: α                ",
            basis: Box::new(|_| vec![1.0]),
        },
        Fam {
            name: "F2: α + c/N²         ",
            basis: Box::new(|n| vec![1.0, 1.0 / (n * n)]),
        },
        Fam {
            name: "F3: α + (βlnN + γ)/N²",
            basis: Box::new(|n| vec![1.0, n.ln() / (n * n), 1.0 / (n * n)]),
        },
        Fam {
            name: "F4: α + c/N          ",
            basis: Box::new(|n| vec![1.0, 1.0 / n]),
        },
        Fam {
            name: "F5: α + β·lnN        ",
            basis: Box::new(|n| vec![1.0, n.ln()]),
        },
    ];
    let xs: Vec<f64> = ns_fit.iter().map(|&n| n as f64).collect();
    let ys: Vec<f64> = ns_fit.iter().map(|&n| s_area[&n]).collect();
    println!("\n    [フィット (mod-0, N=16..64) → holdout (96/128/192) 予測誤差]");
    let mut hold_errs: Vec<f64> = Vec::new();
    let mut alphas: Vec<f64> = Vec::new();
    for f in &fams {
        let coef = lsq(&xs, &ys, &*f.basis);
        let mut herr = 0.0f64;
        for &nh in &ns_hold {
            let b = f.basis.as_ref()(nh as f64);
            let pred: f64 = b.iter().zip(&coef).map(|(bi, ci)| bi * ci).sum();
            herr = herr.max((pred - s_area[&nh]).abs());
        }
        println!(
            "      {}: α = {:.6}  holdout max|Δ| = {:.2e}",
            f.name, coef[0], herr
        );
        hold_errs.push(herr);
        alphas.push(coef[0]);
    }
    let best = (0..fams.len())
        .min_by(|&a, &b| hold_errs[a].partial_cmp(&hold_errs[b]).unwrap())
        .unwrap();
    let worst = (0..fams.len())
        .max_by(|&a, &b| hold_errs[a].partial_cmp(&hold_errs[b]).unwrap())
        .unwrap();
    println!(
        "      → 最良 = {} / 最悪 = {}",
        fams[best].name.trim(),
        fams[worst].name.trim()
    );
    // シェル系統: mod-2 と mod-0 の大 N 差
    let shell_dev = (s_area[&66] - s_area[&64]).abs();
    println!(
        "      シェル系統 (mod2 N=66 vs mod0 N=64): Δs = {:.2e}",
        shell_dev
    );
    // s∞ = 最良族の α (mod-0)。全 N (holdout 込み) での再フィットで確定値
    let xs_all: Vec<f64> = ns_fit
        .iter()
        .chain(ns_hold.iter())
        .map(|&n| n as f64)
        .collect();
    let ys_all: Vec<f64> = ns_fit
        .iter()
        .chain(ns_hold.iter())
        .map(|&n| s_area[&n])
        .collect();
    let coef_best = lsq(&xs_all, &ys_all, &*fams[best].basis);
    let s_inf = coef_best[0];
    println!("      s∞ (最良族, 全 mod-0) = {:.6}", s_inf);

    // ---- [E4] UV 有限量: 相互情報 I(A:B) (N=64, スラブ対) ----
    {
        let n = 64usize;
        let a: Vec<usize> = (0..16).collect();
        let s_a = region_entropy::<f64>(n, &a, nthreads);
        println!("\n    [相互情報 (UV 有限量), N=64, A = [0,16)]");
        for &sep in &[2usize, 4, 8] {
            let b: Vec<usize> = (16 + sep..32 + sep).collect();
            let ab: Vec<usize> = a.iter().cloned().chain(b.iter().cloned()).collect();
            let s_b = region_entropy::<f64>(n, &b, nthreads);
            let s_ab = region_entropy::<f64>(n, &ab, nthreads);
            let mi = s_a + s_b - s_ab;
            println!(
                "      s={}: I = {:.6}  (I/N² = {:.3e}, 種数規格化 I/2 = {:.6})",
                sep,
                mi,
                mi / (n * n) as f64,
                mi / 2.0
            );
        }
    }

    // ---- [E5] G^bare の再分類 ----
    let lam_x = 1.185468; // v24.3 (N=96 直読み)
    let lam_y = 1.229430; // v24.3
    let g_x = 1.0 / (4.0 * lam_x * s_inf);
    let g_y = 1.0 / (4.0 * lam_y * s_inf);
    println!("\n    [G^bare 再分類] G_entropy^bare(a; regulator, 方向スキーム)/a²:");
    println!(
        "      x スキーム (λ_x = {:.6}): G/a² = {:.5} / y,z スキーム (λ_y = {:.6}): G/a² = {:.5}",
        lam_x, g_x, lam_y, g_y
    );
    println!(
        "      → 方向スキーム分だけで {:.1}% の定義依存 — 物理的 Newton 定数ではない",
        ((g_x / g_y) - 1.0).abs() * 100.0
    );
    println!("      (v23.6 の G∞ ≈ 1.06a² は λ = 32/27 [棄却済み]・x スキーム相当の一点)");

    // ---- 判定 ----
    let beta_needed = {
        // F3 の holdout が F2 より 3 倍以上良ければ β 必要と判定
        let f2 = hold_errs[1];
        let f3 = hold_errs[2];
        f3 < f2 / 3.0
    };
    let branch_a = worst == 4 && !beta_needed;
    let branch_b = beta_needed;
    println!(
        "\n[判定] {}",
        if nfail > 0 {
            "装置ゲート故障 — 記録".to_string()
        } else if branch_a {
            format!(
                "事前登録 (a): F5 (Gioev–Klich) が holdout 最悪・F2 で十分 — **定数面積則 s∞ = {:.5} + O(1/N²)** の確定。v23.6 の二択は前提から誤り (点ノードに Fermi 面形を対照した) だが、結論 (定数則) は生存",
                s_inf
            )
        } else if branch_b {
            "事前登録 (b): β·lnN/N² が必要 — 平面対数残滓の発見 (記録)".to_string()
        } else {
            format!(
                "事前登録 (c): 記録 — 最良 {} / 最悪 {}",
                fams[best].name.trim(),
                fams[worst].name.trim()
            )
        }
    );

    // ---- JSON ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v24.6".into())),
        (
            "s_area".into(),
            Json::Arr(
                s_area
                    .iter()
                    .map(|(&n, &s)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("s_area".into(), Json::Num(s)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("s_inf".into(), Json::Num(s_inf)),
        ("g_x".into(), Json::Num(g_x)),
        ("g_y".into(), Json::Num(g_y)),
        (
            "holdout_errs".into(),
            Json::Arr(hold_errs.iter().map(|&x| Json::Num(x)).collect()),
        ),
    ]);
    let p = write_artifact("results/v246_entropyscale.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(b)/(c) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
