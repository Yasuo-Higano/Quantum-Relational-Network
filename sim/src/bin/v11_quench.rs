//! v1.1 実時間・因果・時間の矢の創発 — クエンチ動力学 (厳密)
//!
//! v1.0 の未解決問題 #2 (実時間の定式化) への一歩。格子上のシュレーディンガー動力学は
//! 非相対論的 (原理的には瞬時の影響が可能) なのに、実際には:
//!   (a) 有効な光円錐が創発する (Lieb-Robinson 束縛): 相関の前線は速度 2v_max で広がる
//!   (b) エンタングルメントは光円錐に従い線形成長 → 体積則へ飽和
//!   (c) 全系は純粋なまま、局所的には完全な熱平衡 (ここでは T=∞ Gibbs) に到達
//!       = 時間の矢はもつれの拡散である
//! 初期状態: 電荷密度波 |101010...⟩ (もつれゼロ)。H = -Σ(c†c+h.c.) で発展。
//! 厳密検証: 秩序変数 m(t) = |J₀(4t)| (無限鎖の解析解) との比較。

use uft_sim::*;

/// 長方行列積 C[m×n] = A[m×k] B[k×n] (列優先)
fn matmul_rect(a: &[f64], b: &[f64], m: usize, k: usize, n: usize) -> Vec<f64> {
    let mut c = vec![0.0; m * n];
    for j in 0..n {
        for l in 0..k {
            let blj = b[l + j * k];
            if blj == 0.0 {
                continue;
            }
            for i in 0..m {
                c[i + j * m] += a[i + l * m] * blj;
            }
        }
    }
    c
}

/// ベッセル J0 (積分表示, 台形則)
fn bessel_j0(x: f64) -> f64 {
    let n = 400;
    let mut s = 0.0;
    for i in 0..n {
        let th = std::f64::consts::PI * (i as f64 + 0.5) / n as f64;
        s += (x * th.sin()).cos();
    }
    s / n as f64
}

fn h2(z: f64) -> f64 {
    let z = z.clamp(1e-14, 1.0 - 1e-14);
    -z * z.ln() - (1.0 - z) * (1.0 - z).ln()
}

fn entropy_herm(cre: &[f64], cim: &[f64], n: usize) -> f64 {
    let m = 2 * n;
    let mut a = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            a[i + j * m] = cre[i + j * n];
            a[i + (j + n) * m] = -cim[i + j * n];
            a[(i + n) + j * m] = cim[i + j * n];
            a[(i + n) + (j + n) * m] = cre[i + j * n];
        }
    }
    let (w, _) = jacobi_eigh(&a, m);
    0.5 * w.iter().map(|&z| h2(z)).sum::<f64>()
}

fn main() {
    let n = 400usize;
    println!(
        "=== v1.1 実時間の創発: CDW クエンチの厳密動力学 (N={}) ===\n",
        n
    );
    // OBC 単一粒子モード: φ_k(x) = √(2/(N+1)) sin(πk x/(N+1)), ε_k = -2cos(πk/(N+1))
    let nf = (n + 1) as f64;
    let mut v = vec![0.0; n * n]; // v[x + k*n]
    let mut eps = vec![0.0; n];
    for k in 0..n {
        let kk = (k + 1) as f64;
        eps[k] = -2.0 * (std::f64::consts::PI * kk / nf).cos();
        for x in 0..n {
            v[x + k * n] =
                (2.0 / nf).sqrt() * (std::f64::consts::PI * kk * (x as f64 + 1.0) / nf).sin();
        }
    }
    // 占有サイト (CDW): x = 0,2,4,... (格子座標 1,3,5,.. と思ってよい)
    let occ: Vec<usize> = (0..n).step_by(2).collect();
    let nocc = occ.len();

    // C(t) を計算するクロージャ: U(t) = V e^{-iεt} Vᵀ, C = U_occ U_occ†
    let compute_c = |t: f64| -> (Vec<f64>, Vec<f64>) {
        // P_re[x,k] = V[x,k] cos(ε_k t), P_im = -V sin(ε_k t)
        let mut p_re = vec![0.0; n * n];
        let mut p_im = vec![0.0; n * n];
        for k in 0..n {
            let (c, s) = ((eps[k] * t).cos(), -(eps[k] * t).sin());
            for x in 0..n {
                p_re[x + k * n] = v[x + k * n] * c;
                p_im[x + k * n] = v[x + k * n] * s;
            }
        }
        // U_occ[x, z∈occ] = Σ_k P[x,k] V[z,k]  → まず Vocc[k, j] = V[occ_j, k]
        let mut vocc = vec![0.0; n * nocc];
        for (j, &z) in occ.iter().enumerate() {
            for k in 0..n {
                vocc[k + j * n] = v[z + k * n];
            }
        }
        let ur = matmul_rect(&p_re, &vocc, n, n, nocc);
        let ui = matmul_rect(&p_im, &vocc, n, n, nocc);
        // C = U U†: C_re = Ur Urᵀ + Ui Uiᵀ, C_im = Ui Urᵀ - Ur Uiᵀ
        let mut urt = vec![0.0; nocc * n];
        let mut uit = vec![0.0; nocc * n];
        for x in 0..n {
            for j in 0..nocc {
                urt[j + x * nocc] = ur[x + j * n];
                uit[j + x * nocc] = ui[x + j * n];
            }
        }
        let c1 = matmul_rect(&ur, &urt, n, nocc, n);
        let c2 = matmul_rect(&ui, &uit, n, nocc, n);
        let c3 = matmul_rect(&ui, &urt, n, nocc, n);
        let c4 = matmul_rect(&ur, &uit, n, nocc, n);
        let cre: Vec<f64> = c1.iter().zip(&c2).map(|(a, b)| a + b).collect();
        let cim: Vec<f64> = c3.iter().zip(&c4).map(|(a, b)| a - b).collect();
        (cre, cim)
    };

    // ---- (a)+(b): 秩序変数の厳密解比較 & 光円錐 ----
    println!("[A] 秩序変数 m(t) = |(2/N)Σ(-1)^x(n_x-½)| vs 解析解 |J₀(4t)| (無限鎖)");
    println!("  t     m(数値)   |J0(4t)|");
    let mut max_dev = 0.0f64;
    for &t in &[0.25f64, 0.5, 1.0, 1.5, 2.0, 3.0] {
        let (cre, _) = compute_c(t);
        let mut m = 0.0;
        for x in 100..300 {
            // 端の効果を避け中央で測る
            let sign = if x % 2 == 0 { 1.0 } else { -1.0 };
            m += sign * (cre[x + x * n] - 0.5);
        }
        m = (2.0 * m / 200.0).abs();
        let j0 = bessel_j0(4.0 * t).abs();
        max_dev = max_dev.max((m - j0).abs());
        println!("  {:4.2}  {:.5}   {:.5}", t, m, j0);
    }
    println!("  => 最大偏差 {:.1e}  {}", max_dev, pass(max_dev < 1e-3));
    println!("     (秩序の融解は可逆なユニタリー発展なのに実質的に不可逆 — 位相の脱干渉)\n");

    println!("[B] 創発する光円錐: 相関 |C(x₀,x₀+d)| の前線 (閾値 0.01) — Lieb-Robinson 束縛");
    println!("  t    前線 d*(t)   2v_max·t (v_max=2)");
    let x0 = 100usize;
    let mut ts = Vec::new();
    let mut ds = Vec::new();
    for &t in &[2.0f64, 4.0, 6.0, 8.0, 10.0, 12.0] {
        let (cre, cim) = compute_c(t);
        let mut front = 0usize;
        for d in (1..200).rev() {
            let cabs = (cre[x0 + (x0 + d) * n].powi(2) + cim[x0 + (x0 + d) * n].powi(2)).sqrt();
            if cabs > 0.01 {
                front = d;
                break;
            }
        }
        println!("  {:4.1}   {:4}         {:4.0}", t, front, 4.0 * t);
        ts.push(t);
        ds.push(front as f64);
    }
    let (_, slope) = linfit(&ts, &ds);
    println!(
        "  => 前線速度 = {:.2} (理論値 2v_max = 4: 対から生まれる準粒子が左右へ v_max=2 で走る)  {}",
        slope,
        pass((slope - 4.0).abs() < 0.4)
    );
    println!("     格子の動力学に光速は入っていないのに、有効な因果構造(光円錐)が創発する。\n");

    // ---- (c): もつれの成長と局所熱化 ----
    println!("[C] もつれの成長と局所熱化 (区間 A = 中央 100 サイト)");
    println!("  t     S_A(t)    S_A/(ℓ ln2)   ⟨n⟩(中央)  |⟨c†c⟩|(中央ボンド)");
    let (ia, la) = (150usize, 100usize);
    for &t in &[0.0f64, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0, 40.0, 60.0] {
        let (cre, cim) = compute_c(t);
        let mut are = vec![0.0; la * la];
        let mut aim = vec![0.0; la * la];
        for i in 0..la {
            for j in 0..la {
                are[i + j * la] = cre[(ia + i) + (ia + j) * n];
                aim[i + j * la] = cim[(ia + i) + (ia + j) * n];
            }
        }
        let s = entropy_herm(&are, &aim, la);
        let nb = cre[200 + 200 * n];
        let bond = (cre[200 + 201 * n].powi(2) + cim[200 + 201 * n].powi(2)).sqrt();
        println!(
            "  {:4.0}  {:8.3}   {:.3}        {:.4}     {:.4}",
            t,
            s,
            s / (la as f64 * (2.0f64).ln()),
            nb,
            bond
        );
    }
    println!("  T=∞ の熱平衡値: ⟨n⟩=0.5, ⟨c†c⟩=0, S_A = ℓ·ln2 (最大)");
    println!("\n結論: (a) 光円錐、(b) もつれの線形成長→体積則飽和、(c) 局所観測量の熱平衡化。");
    println!("      全系は S=0 の純粋状態のまま (情報は失われていない) なのに、局所的には");
    println!("      完全に熱的 — 「時間の矢」とは、もつれが広がって戻らないことである。");
    println!("      因果構造 (光円錐) すら動力学から創発する: 相対論的時空構造への実時間の橋。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
