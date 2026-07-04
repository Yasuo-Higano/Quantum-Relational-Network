//! v1.3 時空 = 量子誤り訂正符号 — 3 クォートリット完全符号 (Almheiri–Dong–Harlow)
//!
//! ホログラフィーの難問:「バルク(時空内部)の 1 点の情報は、境界のどこにあるのか?」
//! 答え: どこか特定の場所ではなく、冗長に符号化されている。
//! AdS/CFT の辞書は量子誤り訂正符号の構造を持つ (ADH 2015)。その最小模型:
//!
//!   論理クォートリット 1 個を 3 つの物理クォートリットに符号化:
//!   |0̄⟩=(|000⟩+|111⟩+|222⟩)/√3, |1̄⟩=(|012⟩+|120⟩+|201⟩)/√3, |2̄⟩=(|021⟩+|102⟩+|210⟩)/√3
//!
//! 参照系 (バルクの情報の「原本」) R と最大もつれさせ、
//!   I(R : 任意の1領域) = 0        … どの単一領域もバルクを全く知らない
//!   I(R : 任意の2領域) = 2ln3 (最大) … どの2領域からもバルクを完全復元できる
//! を厳密に検証する。「バルクの点」は境界のどこにもなく、相関の中に遍在する。

use uft_sim::*;

fn h_ent(w: &[f64]) -> f64 {
    w.iter()
        .map(|&l| if l > 1e-14 { -l * l.ln() } else { 0.0 })
        .sum()
}

/// 4 クォートリット純粋状態 (dim 81, 実振幅) の部分系エントロピー
/// keep: 残す因子のリスト (0=ref, 1..3=code qutrits)
fn subsystem_entropy(psi: &[f64], keep: &[usize]) -> f64 {
    let dims = [3usize, 3, 3, 3];
    let dk: usize = keep.iter().map(|&i| dims[i]).product();
    let env: Vec<usize> = (0..4).filter(|i| !keep.contains(i)).collect();
    let de: usize = env.iter().map(|&i| dims[i]).product();
    // インデックス変換
    let full_index = |ks: usize, es: usize| -> usize {
        let mut digits = [0usize; 4];
        let (mut ks, mut es) = (ks, es);
        for &i in keep.iter().rev() {
            digits[i] = ks % 3;
            ks /= 3;
        }
        for &i in env.iter().rev() {
            digits[i] = es % 3;
            es /= 3;
        }
        ((digits[0] * 3 + digits[1]) * 3 + digits[2]) * 3 + digits[3]
    };
    let mut rho = vec![0.0; dk * dk];
    for a in 0..dk {
        for b in 0..=a {
            let mut s = 0.0;
            for e in 0..de {
                s += psi[full_index(a, e)] * psi[full_index(b, e)];
            }
            rho[a + b * dk] = s;
            rho[b + a * dk] = s;
        }
    }
    let (w, _) = jacobi_eigh(&rho, dk);
    h_ent(&w)
}

fn main() {
    println!("=== v1.3 時空=量子誤り訂正: 3 クォートリット完全符号 ===\n");
    // 符号語 (27 次元)
    let idx3 = |a: usize, b: usize, c: usize| (a * 3 + b) * 3 + c;
    let mut code = vec![[0.0f64; 27]; 3];
    let s3 = 1.0 / (3.0f64).sqrt();
    for a in 0..3 {
        // |ā⟩ = (1/√3) Σ_c |c, c+a, c+2a mod 3⟩
        for c in 0..3 {
            code[a][idx3(c, (c + a) % 3, (c + 2 * a) % 3)] = s3;
        }
    }
    // 符号語の直交性確認
    for a in 0..3 {
        for b in 0..3 {
            let ip: f64 = (0..27).map(|i| code[a][i] * code[b][i]).sum();
            assert!((ip - if a == b { 1.0 } else { 0.0 }).abs() < 1e-14);
        }
    }
    println!("符号語の正規直交性: OK");

    // |Φ⟩ = (1/√3) Σ_a |a⟩_R |ā⟩ (参照系との最大もつれ = バルク情報の原本)
    let mut psi = vec![0.0f64; 81];
    for a in 0..3 {
        for i in 0..27 {
            psi[a * 27 + i] = s3 * code[a][i];
        }
    }
    let ln3 = 3.0f64.ln();
    println!("\n[A] バルク情報 (参照系 R) は境界のどこにあるか");
    let s_r = subsystem_entropy(&psi, &[0]);
    println!(
        "  S(R) = {:.6} (= ln3 = {:.6}: R は符号と最大もつれ)",
        s_r, ln3
    );
    println!("\n  領域        S(領域)   I(R:領域)     判定");
    let mut all_ok = true;
    for &i in &[1usize, 2, 3] {
        let s_i = subsystem_entropy(&psi, &[i]);
        let s_ri = subsystem_entropy(&psi, &[0, i]);
        let mi = s_r + s_i - s_ri;
        let ok = mi.abs() < 1e-10;
        all_ok &= ok;
        println!(
            "  クォートリット{}   {:.4}    {:+.2e}     単一領域は何も知らない {}",
            i,
            s_i,
            mi,
            pass(ok)
        );
    }
    for &(i, j) in &[(1usize, 2usize), (1, 3), (2, 3)] {
        let s_ij = subsystem_entropy(&psi, &[i, j]);
        let s_rij = subsystem_entropy(&psi, &[0, i, j]);
        let mi = s_r + s_ij - s_rij;
        let ok = (mi - 2.0 * ln3).abs() < 1e-10;
        all_ok &= ok;
        println!(
            "  ペア({},{})     {:.4}    {:.4}       任意の2領域が完全な情報 (=2ln3) {}",
            i,
            j,
            s_ij,
            mi,
            pass(ok)
        );
    }
    println!("  => {}", pass(all_ok));

    println!("\n[B] 消失誤り耐性: 異なる論理状態でも、単一クォートリットの密度行列は同一か");
    // 2つの論理状態を符号化し、site1 の縮約状態のトレース距離を測る
    let enc = |amp: [f64; 3]| -> Vec<f64> {
        let mut v = vec![0.0f64; 27];
        for a in 0..3 {
            for i in 0..27 {
                v[i] += amp[a] * code[a][i];
            }
        }
        v
    };
    let psi1 = enc([1.0, 0.0, 0.0]);
    let psi2 = enc([0.5f64.sqrt(), 0.5, -0.5]);
    let rho_site = |v: &[f64], site: usize| -> Vec<f64> {
        let mut rho = vec![0.0f64; 9];
        for q1 in 0..3 {
            for q2 in 0..3 {
                for q3 in 0..3 {
                    let q = [q1, q2, q3];
                    for r in 0..3 {
                        let mut q2v = q;
                        q2v[site - 1] = r;
                        rho[q[site - 1] + r * 3] +=
                            v[idx3(q[0], q[1], q[2])] * v[idx3(q2v[0], q2v[1], q2v[2])];
                    }
                }
            }
        }
        rho
    };
    for site in 1..=3 {
        let r1 = rho_site(&psi1, site);
        let r2 = rho_site(&psi2, site);
        let mut diff = vec![0.0; 9];
        for i in 0..9 {
            diff[i] = r1[i] - r2[i];
        }
        let (w, _) = jacobi_eigh(&diff, 3);
        let td: f64 = 0.5 * w.iter().map(|x| x.abs()).sum::<f64>();
        println!(
            "  site{}: トレース距離 = {:.2e} (0 なら区別不能 = 消失しても情報は漏れも失われもしない) {}",
            site,
            td,
            pass(td < 1e-12)
        );
    }

    println!("\n結論: バルクの 1 量子情報は、(i) どの単一境界領域にも存在せず (I=0)、");
    println!("      (ii) どの 2 領域の組にも完全に存在する (I=2ln3)。");
    println!("      「時空の内部の点」とは、境界自由度への冗長な符号化パターンである。");
    println!("      v0.7 の部分領域双対性・v0.8 の島公式はこの誤り訂正構造の帰結。");
    println!("      QRN 公理 A0 の精密化: テンソル分解は誤り訂正符号の構造を持つ。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
