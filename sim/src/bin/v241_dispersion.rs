//! v24.1 分散の解析 — 点ノード・taste 数・開放境界の厳密スペクトル (第二十五期)
//!
//! PROMPT/6 の最短経路の第一歩「分散の解析」。v23.x 系格子 (3+1D staggered,
//! x 開放・y/z 周期, ホップ ±1/2) の解析構造を確定し、dense jacobi と照合する:
//!
//!   (1) 反交換: {D_i, D_j} = 0 (i≠j) — staggered 位相の帰結 (厳密, 半整数演算)
//!   (2) H² = D_x² + D_y² + D_z² — 交差項が反交換で消える (厳密)
//!   (3) バルク分散 E(k) = ±√(cos²kx + cos²ky + cos²kz):
//!       零集合は cos kx = cos ky = cos kz = 0 の **8 点ノード (余次元 3)**。
//!       Fermi 面 (余次元 1) ではない → Gioev–Klich 対数面積則 (S ~ N²·lnN) の
//!       前提が不成立 — v23.6 の「対数 vs 定数」二択は正しい対ではない (v24.6 で再裁定)。
//!   (4) taste 計数: 8 ノード × 1 成分 = 8 = **2 taste × 4 成分 Dirac**
//!       (3+1D staggered の標準結果)。離散化間比較 (v24.5) の種数規格化は 2。
//!   (5) Fermi 速度 v = 1 (等方, ノード近傍で E ≈ |q|)。
//!   (6) x 開放でも厳密: H = D_x + (−1)^x M と (−1)^x φ_n = φ_{N+1−n} により
//!       E = ±√(cos²(πn/(N+1)) + μ²), μ² = cos²ky + cos²kz — **開放境界の
//!       スペクトルが閉形式** (stag.rs のブロック理論 = v24.2 以降の器械の土台)。
//!   (7) シェル gap 公式: gap(N) = 2·min E — N mod 4 のシェル構造 (v23.x の観察) の起源。
//!
//! 事前登録: (a) 全恒等式が dense jacobi と一致 (< 1e-10) = ブロック理論は厳密 —
//!   v24.2 (N=128 級ブロック経路) を解禁 / (b) いずれか破れ = 理論の誤りを記録し、
//!   ブロック経路は使用禁止。

use uft_sim::dd::*;
use uft_sim::stag::*;
use uft_sim::*;

fn main() {
    self_test();
    println!("=== v24.1 分散の解析 — 点ノード・taste・開放境界の厳密スペクトル (第二十五期) ===\n");
    println!("事前登録: (a) 全恒等式が dense と一致 → ブロック理論解禁 / (b) 破れ = 記録\n");
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

    // [A0] 共有モジュールの自己検証
    check(
        "[A0a] dd 自己検証 (EFT/sqrt/ln/exp/有理角)",
        dd_self_test(),
        String::new(),
    );
    check(
        "[A0b] stag 自己検証 (ペア恒等式 DD < 1e-30, ブロック直交規格)",
        stag_self_test(),
        String::new(),
    );

    // [A1][A2] 反交換と H² = ΣD² — 半整数演算なので厳密ゼロを要求
    for &(open_x, label) in &[(false, "全周期"), (true, "x開放")] {
        let n = 6usize;
        let ns = n * n * n;
        let ds = build_d_matrices(n, open_x);
        let mut max_anti = 0.0f64;
        for i in 0..3 {
            for j in i + 1..3 {
                let ab = matmul(&ds[i], &ds[j], ns);
                let ba = matmul(&ds[j], &ds[i], ns);
                for k in 0..ns * ns {
                    max_anti = max_anti.max((ab[k] + ba[k]).abs());
                }
            }
        }
        check(
            &format!("[A1] {{D_i,D_j}} = 0 ({}, N=6)", label),
            max_anti == 0.0,
            format!("max |{{D_i,D_j}}| = {:.1e}", max_anti),
        );
        let mut h = vec![0.0f64; ns * ns];
        for d in &ds {
            for k in 0..ns * ns {
                h[k] += d[k];
            }
        }
        let h2 = matmul(&h, &h, ns);
        let mut max_dev = 0.0f64;
        let d2: Vec<Vec<f64>> = ds.iter().map(|d| matmul(d, d, ns)).collect();
        for k in 0..ns * ns {
            let s = d2[0][k] + d2[1][k] + d2[2][k];
            max_dev = max_dev.max((h2[k] - s).abs());
        }
        check(
            &format!("[A2] H² = ΣD_i² ({}, N=6)", label),
            max_dev == 0.0,
            format!("max |H² − ΣD²| = {:.1e}", max_dev),
        );
    }

    // [A3] 全周期 N=8: H² のスペクトル = {Σcos²(2πq/N)} 多重集合, H は ± 対称
    {
        let n = 8usize;
        let ns = n * n * n;
        let ds = build_d_matrices(n, false);
        let mut h = vec![0.0f64; ns * ns];
        for d in &ds {
            for k in 0..ns * ns {
                h[k] += d[k];
            }
        }
        let (ev, _) = jacobi_eigh(&h, ns);
        let mut dev_sym = 0.0f64;
        for i in 0..ns {
            dev_sym = dev_sym.max((ev[i] + ev[ns - 1 - i]).abs());
        }
        check(
            "[A3a] H スペクトルの ± 対称 (全周期 N=8)",
            dev_sym < 1e-10,
            format!(
                "max |E_i + E_(rev)| = {:.1e} ({} s)",
                dev_sym,
                t0.elapsed().as_secs()
            ),
        );
        let analytic = h2_spectrum_periodic(n);
        let mut e2: Vec<f64> = ev.iter().map(|e| e * e).collect();
        e2.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let dev = e2
            .iter()
            .zip(&analytic)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        check(
            "[A3b] H² スペクトル = {Σcos²k} 多重集合 (全周期 N=8)",
            dev < 1e-10,
            format!("多重集合 max Δ = {:.1e}", dev),
        );
        // [A4] ノード計数: 零モード数 = 8 (k = ±π/2 が格子上にある 4|N)
        let nzero = ev.iter().filter(|e| e.abs() < 1e-10).count();
        check(
            "[A4] 零モード数 = 8 (8 点ノード, 全周期 N=8)",
            nzero == 8,
            format!("零モード {} 本", nzero),
        );
    }

    // [A5] x 開放 N ∈ {8,10,12}: dense jacobi vs 解析多重集合 ±√(cos²k_n + μ²)
    for &n in &[8usize, 10, 12] {
        let ns = n * n * n;
        let h = build_h3d_open_x(n);
        let (ev, _) = jacobi_eigh(&h, ns);
        let analytic = spectrum_analytic_open(n);
        let dev = ev
            .iter()
            .zip(&analytic)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        check(
            &format!("[A5] x開放 N={} スペクトル = 解析閉形式", n),
            dev < 1e-10,
            format!(
                "多重集合 max Δ = {:.1e} ({} s)",
                dev,
                t0.elapsed().as_secs()
            ),
        );
        // [A8] gap の照合
        let nocc = ns / 2;
        let gap_dense = ev[nocc] - ev[nocc - 1];
        let gap_th = gap_analytic(n);
        check(
            &format!("[A8] N={} gap = 解析式", n),
            (gap_dense - gap_th).abs() < 1e-10,
            format!("dense {:.6} vs 解析 {:.6}", gap_dense, gap_th),
        );
    }

    // [A6] Fermi 速度 v = 1 等方 (解析分散のノード近傍展開)
    {
        let e3 = |kx: f64, ky: f64, kz: f64| -> f64 {
            (kx.cos().powi(2) + ky.cos().powi(2) + kz.cos().powi(2)).sqrt()
        };
        let hp = std::f64::consts::FRAC_PI_2;
        let q = 1e-4;
        let vx = e3(hp + q, hp, hp) / q;
        let vy = e3(hp, hp + q, hp) / q;
        let vz = e3(hp, hp, hp + q) / q;
        let s3 = (3.0f64).sqrt();
        let vd = e3(hp + q / s3, hp + q / s3, hp + q / s3) / q;
        let dev = [vx, vy, vz, vd]
            .iter()
            .map(|v| (v - 1.0).abs())
            .fold(0.0f64, f64::max);
        check(
            "[A6] Fermi 速度 v = 1 (x/y/z/対角, 等方)",
            dev < 1e-6,
            format!("v = {:.8}/{:.8}/{:.8}/{:.8}", vx, vy, vz, vd),
        );
    }

    // [A7] 点ノード分類: 80³ BZ 格子で E=0 は厳密に 8 点、ノード球外で E は下に有界
    {
        let m = 80usize;
        let nodes: Vec<[f64; 3]> = {
            let hp = std::f64::consts::FRAC_PI_2;
            let mut v = Vec::new();
            for &sx in &[-1.0f64, 1.0] {
                for &sy in &[-1.0f64, 1.0] {
                    for &sz in &[-1.0f64, 1.0] {
                        v.push([sx * hp, sy * hp, sz * hp]);
                    }
                }
            }
            v
        };
        let mut nzero = 0usize;
        let mut emin_far = f64::MAX;
        for ix in 0..m {
            let kx = -std::f64::consts::PI + 2.0 * std::f64::consts::PI * ix as f64 / m as f64;
            for iy in 0..m {
                let ky = -std::f64::consts::PI + 2.0 * std::f64::consts::PI * iy as f64 / m as f64;
                for iz in 0..m {
                    let kz =
                        -std::f64::consts::PI + 2.0 * std::f64::consts::PI * iz as f64 / m as f64;
                    let e = (kx.cos().powi(2) + ky.cos().powi(2) + kz.cos().powi(2)).sqrt();
                    if e < 1e-12 {
                        nzero += 1;
                    }
                    // 最近接ノードへの距離 (2π 周期)
                    let mut dmin = f64::MAX;
                    for nd in &nodes {
                        let per = |a: f64| {
                            let mut d = (a).abs() % (2.0 * std::f64::consts::PI);
                            if d > std::f64::consts::PI {
                                d = 2.0 * std::f64::consts::PI - d;
                            }
                            d
                        };
                        let d2 = per(kx - nd[0]).powi(2)
                            + per(ky - nd[1]).powi(2)
                            + per(kz - nd[2]).powi(2);
                        dmin = dmin.min(d2.sqrt());
                    }
                    if dmin > 0.3 {
                        emin_far = emin_far.min(e);
                    }
                }
            }
        }
        check(
            "[A7] 零集合 = 孤立 8 点 (余次元 3) — Fermi 面ではない",
            nzero == 8 && emin_far > 0.15,
            format!(
                "BZ 格子の零点 {} 個, ノード球外の min E = {:.4}",
                nzero, emin_far
            ),
        );
    }

    // ---- シェル gap 表 (解析式) と帰結 ----
    println!("\n    シェル gap 表 (解析式 2·min√(cos²k_n + cos²ky + cos²kz)):");
    let mut gap_rows: Vec<(usize, f64)> = Vec::new();
    for &n in &[8usize, 10, 12, 14, 16, 24, 32, 48, 64, 96, 128] {
        let g = gap_analytic(n);
        gap_rows.push((n, g));
        println!(
            "      N={:3} (mod4={}): gap = {:.6}  (≈ {:.3}/N{})",
            n,
            n % 4,
            g,
            g * n as f64,
            if n % 4 == 0 {
                " — 横 μ=0 系列"
            } else {
                " — 横ギャップ系列"
            }
        );
    }

    println!("\n[帰結]");
    println!("  1. 零集合は 8 点ノード (余次元 3) = 2 taste × 4 成分 Dirac, v = 1 等方。");
    println!("     Gioev–Klich 対数面積則 (S ~ N²lnN) は余次元 1 の Fermi 面の結果 —");
    println!("     v23.6 の「対数 vs 定数」二択は前提不成立。正しい有限サイズ族で再裁定 = v24.6。");
    println!("  2. x 開放スペクトルは閉形式 E = ±√(cos²(πn/(N+1)) + cos²ky + cos²kz) —");
    println!("     半空間 C_A は (ky ペア, kz) の N²/2 個の N×N 実対称ブロックに厳密分解。");
    println!("     dense N³ 対角化 (v23.4 の限界 N=16) → ブロック経路で N=128 級へ (v24.2)。");
    println!(
        "  3. シェル構造 (N mod 4) の起源 = 横運動量 π/2 の有無 (μ=0 系列) — gap 公式で閉じた。"
    );
    println!("  4. 離散化間比較 (v24.5) の種数規格化: staggered = 2 Dirac / Wilson = 1 Dirac。");

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v24.1".into())),
        (
            "gap_table".into(),
            Json::Arr(
                gap_rows
                    .iter()
                    .map(|&(n, g)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("gap".into(), Json::Num(g)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("nodes".into(), Json::Int(8)),
        ("tastes".into(), Json::Int(2)),
        ("fermi_velocity".into(), Json::Num(1.0)),
    ]);
    let p = write_artifact("results/v241_dispersion.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 事前登録 (a) — ブロック理論は厳密。v24.2 (ブロック経路) を解禁"
        } else {
            "[FAIL] 事前登録 (b) — 理論の誤りを記録、ブロック経路は使用禁止"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
