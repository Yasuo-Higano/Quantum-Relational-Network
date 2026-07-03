//! v3.2 湯川階層の機構 — フロガット=ニールセン (FN)
//!
//! 観測されるフェルミオン質量は 6 桁にわたる (m_u/m_t ~ 10^-5)。基本定数が
//! 6 桁も「微調整」されているのか? FN 機構 (1979): 世代ごとの U(1) 電荷 q_i と
//! 小さいパラメータ ε (~0.22) があれば、湯川行列は Y_ij = c_ij ε^{q_i+q̄_j}
//! (c_ij は全て O(1)!) となり、階層は「電荷の算術」から自動的に生じる。
//! O(1) 係数をモンテカルロで振り、質量比と CKM 行列が実測と (対数スケールで)
//! 一致するかを検証する。

use uft_sim::*;

/// 3×3 複素行列 (行優先 [ [re,im]; 9 ])
type M3 = [[(f64, f64); 3]; 3];

/// H = Y Y† (エルミート) の固有値と複素固有ベクトル (実埋め込み Jacobi)
fn eig_yyd(y: &M3) -> ([f64; 3], [[(f64, f64); 3]; 3]) {
    // H_ij = Σ_k y_ik conj(y_jk)
    let mut hre = [[0.0f64; 3]; 3];
    let mut him = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let (a, b) = y[i][k];
                let (c, d) = y[j][k];
                hre[i][j] += a * c + b * d;
                him[i][j] += b * c - a * d;
            }
        }
    }
    let n = 3;
    let m = 6;
    let mut emb = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            emb[i + j * m] = hre[i][j];
            emb[i + (j + n) * m] = -him[i][j];
            emb[(i + n) + j * m] = him[i][j];
            emb[(i + n) + (j + n) * m] = hre[i][j];
        }
    }
    let (w, v) = jacobi_eigh(&emb, m);
    // 固有値は 2 重: 列 0,2,4 から複素ベクトル ψ = 上半分 + i 下半分
    let mut lam = [0.0f64; 3];
    let mut vecs = [[(0.0f64, 0.0f64); 3]; 3];
    for k in 0..3 {
        lam[k] = 0.5 * (w[2 * k] + w[2 * k + 1]);
        for i in 0..3 {
            vecs[k][i] = (v[i + (2 * k) * m], v[(i + n) + (2 * k) * m]);
        }
        // 正規化
        let nrm: f64 = vecs[k].iter().map(|&(a, b)| a * a + b * b).sum::<f64>().sqrt();
        for i in 0..3 {
            vecs[k][i].0 /= nrm;
            vecs[k][i].1 /= nrm;
        }
    }
    (lam, vecs)
}

fn main() {
    let mut rng = Rng::new(31415);
    let eps = 0.22f64; // ~ カビボ角
    // FN 電荷 (文献の標準的な割当て)
    let q_q = [3.0f64, 2.0, 0.0];
    let q_u = [4.0f64, 2.0, 0.0];
    let q_d = [1.0f64, 0.0, 0.0];
    let q_l = [1.0f64, 0.0, 0.0];
    let q_e = [4.0f64, 2.0, 0.0];
    println!("=== v3.2 湯川階層: フロガット=ニールセン機構 (ε = {}) ===", eps);
    println!("    FN電荷: q_Q={:?} q_u={:?} q_d={:?} q_L={:?} q_e={:?}\n", q_q, q_u, q_d, q_l, q_e);

    let ntrial = 2000;
    // 記録: [u/t, c/t, d/b, s/b, e/τ, μ/τ, Vus, Vcb, Vub]
    let mut samples: Vec<Vec<f64>> = vec![Vec::with_capacity(ntrial); 9];
    let o1 = |rng: &mut Rng| -> (f64, f64) {
        // |c| ∈ [1/3, 3] 対数一様、位相一様 — 「O(1) で構造なし」
        let r = (3.0f64).powf(2.0 * rng.f64() - 1.0);
        let th = 2.0 * std::f64::consts::PI * rng.f64();
        (r * th.cos(), r * th.sin())
    };
    for _ in 0..ntrial {
        let mut yu: M3 = [[(0.0, 0.0); 3]; 3];
        let mut yd: M3 = [[(0.0, 0.0); 3]; 3];
        let mut ye: M3 = [[(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (a, b) = o1(&mut rng);
                let f = eps.powf(q_q[i] + q_u[j]);
                yu[i][j] = (a * f, b * f);
                let (a, b) = o1(&mut rng);
                let f = eps.powf(q_q[i] + q_d[j]);
                yd[i][j] = (a * f, b * f);
                let (a, b) = o1(&mut rng);
                let f = eps.powf(q_l[i] + q_e[j]);
                ye[i][j] = (a * f, b * f);
            }
        }
        let (lu, vu) = eig_yyd(&yu);
        let (ld, vd) = eig_yyd(&yd);
        let (le, _) = eig_yyd(&ye);
        let mu: Vec<f64> = lu.iter().map(|x| x.max(0.0).sqrt()).collect();
        let md: Vec<f64> = ld.iter().map(|x| x.max(0.0).sqrt()).collect();
        let me: Vec<f64> = le.iter().map(|x| x.max(0.0).sqrt()).collect();
        samples[0].push(mu[0] / mu[2]);
        samples[1].push(mu[1] / mu[2]);
        samples[2].push(md[0] / md[2]);
        samples[3].push(md[1] / md[2]);
        samples[4].push(me[0] / me[2]);
        samples[5].push(me[1] / me[2]);
        // CKM = V_u† V_d (行 = u 固有ベクトル, 列 = d 固有ベクトル; 固有値昇順 → 3行3列目が第3世代)
        let mut ckm = [[0.0f64; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                let (mut re, mut im) = (0.0, 0.0);
                for k in 0..3 {
                    let (a, b) = vu[i][k]; // conj(vu) · vd
                    let (c, d) = vd[j][k];
                    re += a * c + b * d;
                    im += a * d - b * c;
                }
                ckm[i][j] = (re * re + im * im).sqrt();
            }
        }
        // 固有値昇順: 世代 1,2,3 = index 0,1,2 ✓
        samples[6].push(ckm[0][1]); // V_us
        samples[7].push(ckm[1][2]); // V_cb
        samples[8].push(ckm[0][2]); // V_ub
    }
    let obs = [1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037];
    let names = ["m_u/m_t", "m_c/m_t", "m_d/m_b", "m_s/m_b", "m_e/m_τ", "m_μ/m_τ", "|V_us|", "|V_cb|", "|V_ub|"];
    let naive = [
        eps.powi(7), eps.powi(4), eps.powi(4), eps.powi(2), eps.powi(5), eps.powi(2),
        eps, eps * eps, eps.powi(3),
    ];
    println!("  量        FN予言(中央値 [16%,84%])       ε冪の見積り   実測        中央値/実測");
    let mut ok_count = 0;
    for k in 0..9 {
        samples[k].sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = samples[k][ntrial / 2];
        let lo = samples[k][ntrial * 16 / 100];
        let hi = samples[k][ntrial * 84 / 100];
        let ratio = med / obs[k];
        let within = ratio > 0.2 && ratio < 5.0;
        if within {
            ok_count += 1;
        }
        println!(
            "  {:8}  {:9.2e} [{:8.2e},{:8.2e}]   {:8.2e}     {:8.2e}   {:5.2} {}",
            names[k], med, lo, hi, naive[k], obs[k], ratio,
            if within { "✓" } else { " " }
        );
    }
    println!("  => 9 量中 {} 量が実測の 5 倍以内 (対数で 6 桁にわたる階層を、O(1) 係数 + 電荷の算術で再現)", ok_count);

    // 対照実験: FN 構造なし (すべて O(1))
    let mut flat: Vec<f64> = Vec::new();
    for _ in 0..ntrial {
        let mut y: M3 = [[(0.0, 0.0); 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                y[i][j] = o1(&mut rng);
            }
        }
        let (l, _) = eig_yyd(&y);
        let m: Vec<f64> = l.iter().map(|x| x.max(0.0).sqrt()).collect();
        flat.push(m[0] / m[2]);
    }
    flat.sort_by(|a, b| a.partial_cmp(b).unwrap());
    println!("\n  対照 (FN構造なし, 全て O(1)): m_1/m_3 の中央値 = {:.3} — 階層は決して出ない", flat[ntrial / 2]);
    println!("     (実測 m_u/m_t = 1.3e-5 は 4 桁以上外れる)");
    println!("\n結論: 質量の 6 桁の階層は「6 桁の微調整」を要求しない。小さい ε の冪 (=電荷の算術)");
    println!("      と O(1) の乱雑な係数で、質量比と混合角の全パターンが桁で再現される。");
    println!("      ε の起源 (対称性の破れのスケール比) と電荷の割当てが残る問い (正直に)。");
    println!("      QRN: 階層は隠れ構造の「距離/電荷の整数」に由来しうる (v2.2-2.3 の幾何と同型)。");
}
