//! v25.2 BZ moment の区間証明書 — λ_x, λ_⊥ の包含区間と異方性の証明 (第二十六期, II/IV)
//!
//! v25.1 (厳密還元) + v25.2-I (厳密閉形式 g = 1/AGM(1,√(1+μ²))) の帰結として、
//!   1/λ_x = r_x = ⟨g(μ)⟩_BZ,   1/λ_⊥ = r_⊥ = ⟨2cos²ky·g(μ)⟩_BZ,
//!   μ² = cos²ky + cos²kz,  (ky,kz) は BZ 一様
//! は**引数不要の確定した数**になった。本版はこれを「倍精度で合った」から一段上げ、
//! (i) 独立 2 実装の相互照合 (1e-13 級) と (ii) 区間演算による厳密包含 (証明) を与える。
//!
//! 実装 (共通の補間表・離散化を共有しない 3 経路 — 共通故障モードの排除):
//!   A: 3D トーラス中点則 — g の積分表示を展開した単一被積分関数
//!      F(θ,ky,kz) = [1+(cos²ky+cos²kz)cos²θ]^{−1/2} の (θ,ky,kz) 平均。AGM 不使用。
//!   B: 2D BZ 中点則 + AGM (lib::agm) — θ 積分を AGM で厳密に畳んだもの。
//!   R: **区間 Riemann** (iv.rs, 外向き丸め・cos は剰余項つき Taylor・π は区間定数):
//!      [0,π/2]² を n=32768 (2^15) 分割、各セルで cos の単調性から Y=cos²ky の厳密範囲
//!      → g の区間 → 階層和。**全ての真値がセル区間に含まれる構成**なので
//!      Σ/n² は r の厳密な包含区間 (n² = 2^30 での除算は正確)。幅は O(h) —
//!      モノトーン包の Riemann 和の宿命で、証明幅 ~1e-5 と浮動小数 15 桁の分業。
//!
//! 異方性定理 (解析 — docs/uft-v25.2.md §2 で証明): Y = cos²ky, Z = cos²kz は独立、
//! E[Y] = 1/2、h(s) = g(√s) は狭義単調減少 (g′<0, v252_exact_g S6e) なので
//!   r_⊥ − r_x = 2Cov(Y, h(Y+Z)) < 0  (Chebyshev の相関不等式を Z で条件付けて適用)
//! ⟹ **λ_⊥ > λ_x は幾何の必然** (Lorentz 破れではなくチャネル重みと質量依存
//! prefactor の負の共分散)。本バイナリは同じ不等式を区間演算で数値的に**証明**する
//! (r_x.lo − r_⊥.hi > 0)。
//!
//! 検査:
//!  [S0] iv.rs 自己検証 (π/cos(π/3)=1/2/AGM(1,√2)=ガウス定数逆数 等の外部照合)
//!  [S1] A の倍増収束 (中点則は解析的周期関数で幾何収束)
//!  [S2] B の倍増収束
//!  [S3] A = B (≤ 5e-13 相対) — 独立実装の相互照合
//!  [S4] R ∋ A かつ R ∋ B (厳密区間が浮動小数推定を包含)
//!  [S5] 区間幅 ≤ 1e-4 (n=32768 の O(h) 予算内)
//!  [S6] **異方性の区間証明**: (r_x − r_⊥).lo > 0 と λ_⊥ − λ_x の証明付き下界
//!  [S7] v24.3 直接測定 (λ_x = 1.185468 / λ_⊥ = 1.229430) との整合 ≤ 5e-5
//!  [S8] 変異検出層: 重み sin²ky / 重み係数 1 / μ² の kz 落とし — 全て検出
//!
//! 事前登録分岐: (a) 全 PASS → λ の解析的完結 (証明付き数値 + 異方性定理) /
//!   (b) S4 FAIL → 区間か浮動小数の実装欠陥 (監査) / (c) S6 FAIL → n 不足 (増強して再走) /
//!   (d) S7 FAIL → 3D 測定と解析値の不一致 = v25.1 還元の見落とし (重大 — 再監査)。
//!
//! 証明書: results/v252_bz_certificate.json (式・領域・規格化・格子・区間・SHA-256)。
//! 決定性のため日付・commit hash は含めない (それらは凍結マニフェスト = IV が持つ)。

use uft_sim::iv::*;
use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// 実装 A: 3D トーラス中点則 (n³ 格子)。戻り値 (r_x, r_⊥)。
/// 行 (θ 添字) ごとの部分和 → 添字順の畳み込み — スレッド数に依らず決定的。
fn torus3d(n: usize, nthreads: usize) -> (f64, f64) {
    let ct2: Vec<f64> = (0..n)
        .map(|i| {
            let t = 2.0 * PI * (i as f64 + 0.5) / n as f64;
            t.cos() * t.cos()
        })
        .collect();
    let mut rows_x = vec![0.0f64; n];
    let mut rows_p = vec![0.0f64; n];
    let chunk = n.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for ((tid, rx), rp) in rows_x
            .chunks_mut(chunk)
            .enumerate()
            .zip(rows_p.chunks_mut(chunk))
        {
            let ct2 = &ct2;
            sc.spawn(move || {
                for (j, (sx, sp)) in rx.iter_mut().zip(rp.iter_mut()).enumerate() {
                    let i = tid * chunk + j; // θ 添字
                    let c2t = ct2[i];
                    let (mut ax, mut ap) = (0.0f64, 0.0f64);
                    for &cy2 in ct2 {
                        let mut row = 0.0f64;
                        for &cz2 in ct2 {
                            row += 1.0 / (1.0 + (cy2 + cz2) * c2t).sqrt();
                        }
                        ax += row;
                        ap += 2.0 * cy2 * row;
                    }
                    *sx = ax;
                    *sp = ap;
                }
            });
        }
    });
    let n3 = (n * n * n) as f64;
    (
        rows_x.iter().sum::<f64>() / n3,
        rows_p.iter().sum::<f64>() / n3,
    )
}

/// 実装 B: 2D BZ 中点則 + AGM (θ 方向を厳密に畳んだ g を直接평均)
fn bz2d(n: usize) -> (f64, f64) {
    let c2: Vec<f64> = (0..n)
        .map(|i| {
            let t = 2.0 * PI * (i as f64 + 0.5) / n as f64;
            t.cos() * t.cos()
        })
        .collect();
    let g_of_mu2 = |m2: f64| 1.0 / agm(1.0, (1.0 + m2).sqrt()).unwrap();
    let (mut sx, mut sp) = (0.0f64, 0.0f64);
    for &cy2 in &c2 {
        let (mut rx, mut rp) = (0.0f64, 0.0f64);
        for &cz2 in &c2 {
            let g = g_of_mu2(cy2 + cz2);
            rx += g;
            rp += 2.0 * cy2 * g;
        }
        sx += rx;
        sp += rp;
    }
    let n2 = (n * n) as f64;
    (sx / n2, sp / n2)
}

/// 実装 R: 区間 Riemann。[0,π/2]² の n×n セル (n = 2^13)、厳密包含を返す。
fn interval_riemann(n: usize, nthreads: usize) -> (Iv, Iv) {
    assert!(n.is_power_of_two());
    // 格子点 t_i = (π/2)·i/n の cos の区間 (i = 0..=n)。cos は [0,π/2] で単調減少。
    let half_pi = iv_pi().half();
    let cgrid: Vec<Iv> = (0..=n)
        .map(|i| cos_iv(half_pi.mul(iv(i as f64)).div(iv(n as f64))))
        .collect();
    // セル i の Y = cos²k の厳密範囲 (真の範囲 [cos²t_{i+1}, cos²t_i] を包含)
    let ycell: Vec<Iv> = (0..n)
        .map(|i| {
            let lo = cgrid[i + 1].lo.max(0.0);
            Iv {
                lo: fnext_down(lo * lo),
                hi: fnext_up(cgrid[i].hi * cgrid[i].hi),
            }
        })
        .collect();
    let mut rows_x = vec![IV_ZERO; n];
    let mut rows_p = vec![IV_ZERO; n];
    let chunk = n.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for ((tid, rx), rp) in rows_x
            .chunks_mut(chunk)
            .enumerate()
            .zip(rows_p.chunks_mut(chunk))
        {
            let ycell = &ycell;
            sc.spawn(move || {
                for (j, (sx, sp)) in rx.iter_mut().zip(rp.iter_mut()).enumerate() {
                    let iy = tid * chunk + j;
                    let y = ycell[iy];
                    let mut row = IV_ZERO;
                    for z in ycell {
                        row = row.add(g_iv_mu2(y.add(*z)));
                    }
                    *sx = row;
                    // 重み 2Y は行内で共通 — 因子化した積の方が区間は厳密かつタイト
                    *sp = y.dbl().mul(row);
                }
            });
        }
    });
    let (mut tx, mut tp) = (IV_ZERO, IV_ZERO);
    for i in 0..n {
        tx = tx.add(rows_x[i]);
        tp = tp.add(rows_p[i]);
    }
    let k = 2 * n.trailing_zeros(); // ÷ n² は 2 冪除算で正確
    (tx.div_pow2(k), tp.div_pow2(k))
}

fn main() {
    self_test();
    println!(
        "=== v25.2 BZ moment 区間証明書 — λ_x, λ_⊥ の包含と異方性の証明 (第二十六期, II/IV) ===\n"
    );
    println!("事前登録: (a) 全 PASS → λ の解析的完結 / (b) S4 FAIL → 実装欠陥 /");
    println!("          (c) S6 FAIL → n 不足 / (d) S7 FAIL → v25.1 還元の見落とし (重大)\n");
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

    check("[S0] 区間演算の自己検証", iv_self_test(), String::new());

    // ---- [S1] 実装 A (3D トーラス, AGM 不使用) ----
    let mut a_seq = Vec::new();
    for &n in &[64usize, 128, 256, 512] {
        a_seq.push((n, torus3d(n, nthreads)));
    }
    let (arx, arp) = a_seq.last().unwrap().1;
    let (prx, prp) = a_seq[a_seq.len() - 2].1;
    let da = ((arx - prx).abs() / arx).max((arp - prp).abs() / arp);
    for (n, (x, p)) in &a_seq {
        println!("    [A n={:4}] r_x = {:.15}  r_⊥ = {:.15}", n, x, p);
    }
    check(
        "[S1] A の倍増収束 (256→512)",
        da < 1e-13,
        format!("|Δ| = {:.1e} ({} s)", da, t0.elapsed().as_secs()),
    );

    // ---- [S2] 実装 B (2D BZ + AGM) ----
    let mut b_seq = Vec::new();
    for &n in &[512usize, 1024, 2048, 4096] {
        b_seq.push((n, bz2d(n)));
    }
    let (brx, brp) = b_seq.last().unwrap().1;
    let (qrx, qrp) = b_seq[b_seq.len() - 2].1;
    let db = ((brx - qrx).abs() / brx).max((brp - qrp).abs() / brp);
    for (n, (x, p)) in &b_seq {
        println!("    [B n={:4}] r_x = {:.15}  r_⊥ = {:.15}", n, x, p);
    }
    check(
        "[S2] B の倍増収束 (2048→4096)",
        db < 1e-13,
        format!("|Δ| = {:.1e} ({} s)", db, t0.elapsed().as_secs()),
    );

    // ---- [S3] A = B ----
    let dab = ((arx - brx).abs() / brx).max((arp - brp).abs() / brp);
    check(
        "[S3] 独立実装の相互照合 A = B",
        dab < 5e-13,
        format!("max 相対差 = {:.1e}", dab),
    );

    // ---- [R] 区間 Riemann ----
    let n_r = 32768usize;
    let (rx_iv, rp_iv) = interval_riemann(n_r, nthreads);
    println!(
        "\n    [R n={}] r_x ∈ [{:.15}, {:.15}] (幅 {:.1e})",
        n_r,
        rx_iv.lo,
        rx_iv.hi,
        rx_iv.width()
    );
    println!(
        "               r_⊥ ∈ [{:.15}, {:.15}] (幅 {:.1e})  ({} s)",
        rp_iv.lo,
        rp_iv.hi,
        rp_iv.width(),
        t0.elapsed().as_secs()
    );
    check(
        "[S4] 厳密包含 R ∋ A かつ R ∋ B",
        rx_iv.contains(arx) && rx_iv.contains(brx) && rp_iv.contains(arp) && rp_iv.contains(brp),
        String::new(),
    );
    check(
        "[S5] 区間幅 ≤ 1e-4 (O(h) 予算, h = (π/2)/32768)",
        rx_iv.width() < 1e-4 && rp_iv.width() < 1e-4,
        format!("幅 = {:.1e} / {:.1e}", rx_iv.width(), rp_iv.width()),
    );

    // ---- [S6] 異方性の区間証明 ----
    let diff = rx_iv.sub(rp_iv);
    let lam_x = IV_ONE.div(rx_iv);
    let lam_p = IV_ONE.div(rp_iv);
    let dl = IV_ONE.div(rp_iv).sub(IV_ONE.div(rx_iv));
    println!(
        "\n    [証明付き数値] λ_x ∈ [{:.15}, {:.15}]",
        lam_x.lo, lam_x.hi
    );
    println!(
        "                   λ_⊥ ∈ [{:.15}, {:.15}]",
        lam_p.lo, lam_p.hi
    );
    println!(
        "                   r_x − r_⊥ ≥ {:.6e},  λ_⊥ − λ_x ≥ {:.6e}",
        diff.lo, dl.lo
    );
    check(
        "[S6] 異方性の区間証明: (r_x − r_⊥).lo > 0 ⟹ λ_⊥ > λ_x (厳密)",
        diff.lo > 0.0,
        format!("証明付き下界 λ_⊥ − λ_x ≥ {:.5}", dl.lo),
    );

    // ---- [S7] v24.3 直接測定との整合 ----
    let (lx_b, lp_b) = (1.0 / brx, 1.0 / brp);
    let (e_x, e_p) = ((lx_b / 1.185468 - 1.0).abs(), (lp_b / 1.229430 - 1.0).abs());
    check(
        "[S7] v24.3 直接測定との整合 (λ_x 1.185468 / λ_⊥ 1.229430)",
        e_x < 5e-5 && e_p < 5e-5,
        format!(
            "λ_x = {:.15} (相対差 {:.1e}) / λ_⊥ = {:.15} ({:.1e})",
            lx_b, e_x, lp_b, e_p
        ),
    );

    // ---- [S8] 変異検出層 (B 経路, n=1024) ----
    {
        let n = 1024usize;
        let c2: Vec<f64> = (0..n)
            .map(|i| {
                let t = 2.0 * PI * (i as f64 + 0.5) / n as f64;
                t.cos() * t.cos()
            })
            .collect();
        let g_of_mu2 = |m2: f64| 1.0 / agm(1.0, (1.0 + m2).sqrt()).unwrap();
        let (mut m1, mut m2v, mut m3) = (0.0f64, 0.0f64, 0.0f64);
        for &cy2 in &c2 {
            for &cz2 in &c2 {
                let g = g_of_mu2(cy2 + cz2);
                m1 += 2.0 * (1.0 - cy2) * g; // 重みを sin²ky に取り違え
                m2v += cy2 * g; // 係数 2 を忘れ
                m3 += 2.0 * cy2 * g_of_mu2(cy2); // μ² の kz 落とし
            }
        }
        let n2 = (n * n) as f64;
        let devs = [
            ("重み sin²", (m1 / n2 - brp).abs()),
            ("係数 2 忘れ", (m2v / n2 - brp).abs()),
            ("μ² の kz 落とし", (m3 / n2 - brp).abs()),
        ];
        for (nm, d) in devs {
            check(
                &format!("[S8-{}] 変異検出", nm),
                d > 0.01,
                format!("逸脱 {:.3} > 0.01 (S3/S4/S7 ゲートが検出)", d),
            );
        }
    }

    // ---- 証明書 JSON (決定的 — 日付・commit は凍結マニフェスト IV が持つ) ----
    let src = std::fs::read("sim/src/bin/v252_bz_certificate.rs")
        .or_else(|_| std::fs::read("src/bin/v252_bz_certificate.rs"))
        .unwrap_or_default();
    let libiv = std::fs::read("sim/src/iv.rs")
        .or_else(|_| std::fs::read("src/iv.rs"))
        .unwrap_or_default();
    let spec = "g(mu)=1/AGM(1,sqrt(1+mu^2)) [=(2/pi)*kappa*K(kappa'), Eisler arXiv:2410.16433 (72) -> QRN K_NN=g*pi*xi]; mu^2=cos^2(ky)+cos^2(kz); r_x=<g>_BZ; r_perp=<2cos^2(ky)*g>_BZ; lambda=1/r; BZ uniform, reduced to [0,pi/2]^2 by cos^2 symmetry";
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v25.2".into())),
        ("kind".into(), Json::Str("bz_moment_interval_certificate".into())),
        ("spec".into(), Json::Str(spec.into())),
        ("spec_sha256".into(), Json::Str(sha256_hex(spec.as_bytes()))),
        (
            "method".into(),
            Json::Obj(vec![
                ("A".into(), Json::Str("3D torus midpoint n=512 (integral form, no AGM)".into())),
                ("B".into(), Json::Str("2D BZ midpoint n=4096 + f64 AGM".into())),
                (
                    "R".into(),
                    Json::Str(format!(
                        "interval Riemann n={} on [0,pi/2]^2, outward rounding, cos via remainder-padded Taylor, interval AGM, hierarchical sums, exact /2^k",
                        n_r
                    )),
                ),
            ]),
        ),
        (
            "results".into(),
            Json::Obj(vec![
                ("r_x_float".into(), Json::Num(brx)),
                ("r_perp_float".into(), Json::Num(brp)),
                ("lambda_x_float".into(), Json::Num(lx_b)),
                ("lambda_perp_float".into(), Json::Num(lp_b)),
                ("r_x_lo".into(), Json::Num(rx_iv.lo)),
                ("r_x_hi".into(), Json::Num(rx_iv.hi)),
                ("r_perp_lo".into(), Json::Num(rp_iv.lo)),
                ("r_perp_hi".into(), Json::Num(rp_iv.hi)),
                ("lambda_x_lo".into(), Json::Num(lam_x.lo)),
                ("lambda_x_hi".into(), Json::Num(lam_x.hi)),
                ("lambda_perp_lo".into(), Json::Num(lam_p.lo)),
                ("lambda_perp_hi".into(), Json::Num(lam_p.hi)),
                ("anisotropy_lower_bound".into(), Json::Num(dl.lo)),
            ]),
        ),
        (
            "anisotropy_theorem".into(),
            Json::Str(
                "r_perp - r_x = 2 Cov(Y, h(Y+Z)) < 0 with Y=cos^2 ky, Z=cos^2 kz independent, E[Y]=1/2, h(s)=g(sqrt(s)) strictly decreasing (g'<0) — Chebyshev correlation inequality conditioned on Z; interval proof: r_x.lo - r_perp.hi > 0".into(),
            ),
        ),
        (
            "sha256".into(),
            Json::Obj(vec![
                ("src".into(), Json::Str(sha256_hex(&src))),
                ("iv_rs".into(), Json::Str(sha256_hex(&libiv))),
            ]),
        ),
    ]);
    let p = write_artifact("results/v252_bz_certificate.json", &j.render());
    println!("\n[artifact] {}", p);

    // ---- 判定 ----
    println!(
        "\n[判定] {}",
        if nfail == 0 {
            format!(
                "事前登録 (a): **λ の解析的完結** — λ_x ∈ [{:.13}, {:.13}], λ_⊥ ∈ [{:.13}, {:.13}] (証明付き), λ_⊥ − λ_x ≥ {:.5} (異方性は負共分散の定理)",
                lam_x.lo, lam_x.hi, lam_p.lo, lam_p.hi, dl.lo
            )
        } else {
            "FAIL — 分岐 (b)/(c)/(d) は各検査の欄を一次ソースとする".to_string()
        }
    );
    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
