//! v31.4 不変 operator response atlas と factorization no-go (PROMPT/12 第三十一期)
//!
//! **Part A — whitened operator atlas**: 各ノードの parity-even 局所演算子基底
//! O_i^a = c†M^a c (2 軌道ノード: n₁, n₂, Re hop, Im hop) に対し
//!   静的:  X_ij^{ab} = ⟨O_i^a O_j^b⟩_c = Tr(M^a (I−C) M^b C)   (Wick)
//!   動的:  R_ij^{ab}(t) = i⟨[O_i^a(t), O_j^b]⟩ = i Tr(C [M^a(t), M^b])
//!   局所 Gram: G_i = X_ii → whitening X̂_ij = G_i^{−1/2} X_ij G_j^{−1/2}
//! の特異スペクトル・作用素/Frobenius/核ノルム・時間積分ノルムを不変量にする —
//! **局所演算子基底の任意の可逆再結合 L_i に対して不変** (unitary に限らない)。
//! B3 は atlas の (n,n) 成分・B4 は動的成分の到着時刻圧縮として位置づける。
//! rank 欠損 Gram には ObservableSupportCertificate を返し無条件擬似逆を禁止。
//!
//! **Part B — factorization no-go**: 同一の大域 Gaussian 状態に対し
//!   (i) site 因子分解 → ring 幾何 (12 辺)
//!   (ii) eigenmode 因子分解 → 幾何なし (モード間静的核は厳密 0 — C が対角化される)
//!   (iii) pair 回転因子分解 → 別の自己整合幾何
//! の 3 通りが全て「同じ状態の正しい読み」— **state-only の静的共分散では spatial
//! factorization を一意選択できない** (v29.5 [C5] の定理化)。さらに自然な state-only
//! 選択基準 (K の疎性) は mode 基底 (nnz 最小 = 対角) を選ぶ — 幾何は自明化する。
//! 選択には OperationalAlgebra (準備・介入・測定・両立性) が必要: site-local probe の
//! 応答は ring を返し、mode-local probe の応答は厳密 0 — **因子分解は状態でなく
//! 操作代数が運ぶ**。state-only で factorization が「読めた」場合は hidden basis
//! convention の流入を疑い成功扱いにしない (負制御を常設)。
//!
//! 検査: [A0] Wick 二重実装照合 [A1] 可逆再結合不変性 [A2] 動的 atlas 二重経路
//! [A3] B3/B4 = atlas の成分/圧縮 [A4] rank 欠損の支持証明書 [A5] mode 対角化 no-go +
//! 疎性基準の負制御 [A6] 3 因子分解の相異 [A7] OperationalAlgebra が幾何を選ぶ
//!
//! 実行: cargo run --release --bin v314_operator_atlas

use uft_sim::readout_contract::*;
use uft_sim::{jacobi_eigh, matfun_sym, matmul, Rng};

fn gibbs_c(h: &[f64], n: usize, beta: f64) -> Vec<f64> {
    matfun_sym(h, n, |x| 1.0 / (1.0 + (beta * x).exp()))
}

/// 静的 connected covariance X_AB = Tr(M_A (I−C) M_B C) (実対称 M, C)
fn cov_static(ma: &[f64], mb: &[f64], c: &[f64], n: usize) -> f64 {
    // 4 重ループ版 (独立照合用の素朴実装)
    let mut s = 0.0;
    for p in 0..n {
        for q in 0..n {
            if ma[p * n + q] == 0.0 {
                continue;
            }
            for r in 0..n {
                for t in 0..n {
                    if mb[r * n + t] == 0.0 {
                        continue;
                    }
                    let d_qr = if q == r { 1.0 } else { 0.0 };
                    s += ma[p * n + q] * mb[r * n + t] * c[p * n + t] * (d_qr - c[r * n + q]);
                }
            }
        }
    }
    s
}

/// 同じ量の行列積版 X_AB = Tr(M_A (I−C) M_B C)
fn cov_static_mat(ma: &[f64], mb: &[f64], c: &[f64], n: usize) -> f64 {
    let mut imc = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            imc[i * n + j] = (if i == j { 1.0 } else { 0.0 }) - c[i * n + j];
        }
    }
    let t1 = matmul(ma, &imc, n);
    let t2 = matmul(mb, c, n);
    let t3 = matmul(&t1, &t2, n);
    (0..n).map(|i| t3[i * n + i]).sum()
}

/// 動的応答 R_AB(t) = i⟨[A(t), B]⟩ = i Tr(C [M_A(t), M_B])。
/// M_A(t) = e^{iht} M_A e^{−iht} を h の固有系で構成 (実対称 h)。
/// 返り値は実数 (エルミート交換子の i 倍の期待値)。
fn response_r(
    vals: &[f64],
    vecs: &[f64],
    ma: &[f64],
    mb: &[f64],
    c: &[f64],
    n: usize,
    t: f64,
) -> f64 {
    // M̃_A = Vᵀ M_A V, 位相 e^{i(λa−λb)t} を掛けて戻す — 実部/虚部を分けて実装
    let mut mt = vec![0.0; n * n]; // Ṽ: Vᵀ M V
    for a in 0..n {
        for b in 0..n {
            let mut s = 0.0;
            for i in 0..n {
                for j in 0..n {
                    s += vecs[a * n + i] * ma[i * n + j] * vecs[b * n + j];
                }
            }
            mt[a * n + b] = s;
        }
    }
    // M_A(t)_{ij} = Σ_ab V_a(i) e^{i(λa−λb)t} M̃_ab V_b(j) → 実部 re, 虚部 im
    let mut mre = vec![0.0; n * n];
    let mut mim = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut sre = 0.0;
            let mut sim = 0.0;
            for a in 0..n {
                for b in 0..n {
                    let ph = (vals[a] - vals[b]) * t;
                    let w = vecs[a * n + i] * mt[a * n + b] * vecs[b * n + j];
                    sre += w * ph.cos();
                    sim += w * ph.sin();
                }
            }
            mre[i * n + j] = sre;
            mim[i * n + j] = sim;
        }
    }
    // i Tr(C [M_A(t), M_B]) — M_A(t) = mre + i·mim, M_B 実, C 実対称:
    // [M_A(t), M_B] = ([mre,M_B]) + i([mim,M_B]) → i Tr(C·..) の実部 = −Tr(C [mim, M_B])
    let mut s = 0.0;
    for i in 0..n {
        for j in 0..n {
            for k in 0..n {
                s += c[i * n + j] * (mim[j * n + k] * mb[k * n + i] - mb[j * n + k] * mim[k * n + i]);
            }
        }
    }
    -s
}

/// 2×2 対称行列の逆平方根 (G^{-1/2}) — whitening 用 (一般 m×m は jacobi で)
fn inv_sqrt_sym(g: &[f64], m: usize) -> Vec<f64> {
    matfun_sym(g, m, |x| {
        assert!(x > 1e-12, "Gram が正定でない — 支持証明書経路を使うこと");
        1.0 / x.sqrt()
    })
}

/// m_i × m_j の X block の特異値 (降順)
fn singvals(x: &[f64], mi: usize, mj: usize) -> Vec<f64> {
    // XᵀX (mj×mj) の固有値の平方根
    let mut g = vec![0.0; mj * mj];
    for a in 0..mj {
        for b in 0..mj {
            let mut s = 0.0;
            for r in 0..mi {
                s += x[r * mj + a] * x[r * mj + b];
            }
            g[a * mj + b] = s;
        }
    }
    let (ev, _) = jacobi_eigh(&g, mj);
    let mut sv: Vec<f64> = ev.iter().map(|&e| e.max(0.0).sqrt()).collect();
    sv.sort_by(|a, b| b.partial_cmp(a).unwrap());
    sv
}

/// ノード i の parity-even 演算子基底 (2 軌道): n₁, n₂, Re hop, Im hop → M 行列 4 本
fn node_ops(sites: &[usize], n: usize) -> Vec<Vec<f64>> {
    let (s1, s2) = (sites[0], sites[1]);
    let mut ops = Vec::new();
    let mut m1 = vec![0.0; n * n];
    m1[s1 * n + s1] = 1.0;
    ops.push(m1);
    let mut m2 = vec![0.0; n * n];
    m2[s2 * n + s2] = 1.0;
    ops.push(m2);
    let mut m3 = vec![0.0; n * n];
    m3[s1 * n + s2] = 1.0;
    m3[s2 * n + s1] = 1.0;
    ops.push(m3);
    // 注: 2 サイトノードの実対称 quadratic 空間は 3 次元 (n₁, n₂, Re hop) で完備。
    // Im hop はエルミートだが M が複素 — 実 Wick 枠の外 (複素側は v31.2 が担保)。
    ops
}

/// atlas block X_ij (m×m) を構成
fn atlas_block(
    ops_i: &[Vec<f64>],
    ops_j: &[Vec<f64>],
    c: &[f64],
    n: usize,
) -> Vec<f64> {
    let (mi, mj) = (ops_i.len(), ops_j.len());
    let mut x = vec![0.0; mi * mj];
    for a in 0..mi {
        for b in 0..mj {
            x[a * mj + b] = cov_static_mat(&ops_i[a], &ops_j[b], c, n);
        }
    }
    x
}

fn frob(x: &[f64]) -> f64 {
    x.iter().map(|v| v * v).sum::<f64>().sqrt()
}

fn main() {
    uft_sim::self_test();
    println!("=== v31.4 不変 operator response atlas / factorization no-go (PROMPT/12) ===\n");
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

    // 試験系: 8 サイト (4 ノード × 2 軌道) 重みつき環 + ring12 (1 サイトノード)
    let n8 = 8usize;
    let h8 = {
        let mut h = vec![0.0; n8 * n8];
        let mut rng = Rng::new(314);
        for b in 0..4 {
            let (i, j) = (2 * b, 2 * b + 1);
            let t = 0.8 + 0.4 * rng.f64();
            h[i * n8 + j] = -t;
            h[j * n8 + i] = -t;
        }
        for b in 0..4 {
            let (i, j) = (2 * b + 1, (2 * b + 2) % n8);
            let t = 0.9 + 0.3 * rng.f64();
            h[i * n8 + j] = -t;
            h[j * n8 + i] = -t;
        }
        h
    };
    let beta = 1.0;
    let c8 = gibbs_c(&h8, n8, beta);
    let nodes8: Vec<Vec<usize>> = (0..4).map(|b| vec![2 * b, 2 * b + 1]).collect();
    let ops8: Vec<Vec<Vec<f64>>> = nodes8.iter().map(|s| node_ops(s, n8)).collect();

    // ---- [A0] Wick 二重実装照合 + n-n 成分の解析値 ----
    {
        let mut worst: f64 = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                for a in 0..3 {
                    for b in 0..3 {
                        let v1 = cov_static(&ops8[i][a], &ops8[j][b], &c8, n8);
                        let v2 = cov_static_mat(&ops8[i][a], &ops8[j][b], &c8, n8);
                        worst = worst.max((v1 - v2).abs());
                    }
                }
            }
        }
        // n-n 成分の解析値: X(n_p, n_q) = −C_pq² (p ≠ q, 実 C)
        let mut worst_nn: f64 = 0.0;
        for p in 0..n8 {
            for q in 0..n8 {
                if p == q {
                    continue;
                }
                let mut mp = vec![0.0; n8 * n8];
                mp[p * n8 + p] = 1.0;
                let mut mq = vec![0.0; n8 * n8];
                mq[q * n8 + q] = 1.0;
                let x = cov_static_mat(&mp, &mq, &c8, n8);
                worst_nn = worst_nn.max((x + c8[p * n8 + q] * c8[p * n8 + q]).abs());
            }
        }
        check(
            "[A0] Wick 静的共分散: 4 重ループ vs 行列積の二重実装 ≤ 1e-13 + n-n 成分 = −C_pq² 厳密",
            worst <= 1e-13 && worst_nn <= 1e-14,
            format!("二重実装乖離 {:.2e} / n-n 解析値乖離 {:.2e}", worst, worst_nn),
        );
    }

    // ---- [A1] whitening の可逆再結合不変性 ----
    {
        let mut rng = Rng::new(3141);
        let gs: Vec<Vec<f64>> = (0..4)
            .map(|i| atlas_block(&ops8[i], &ops8[i], &c8, n8))
            .collect();
        let m = 3usize;
        let mut worst: f64 = 0.0;
        // 各ノードにランダム可逆 L_i (unitary ではない) を掛けた基底で再計算
        let mut draws = 0usize;
        while draws < 8 {
            let mut ls: Vec<Vec<f64>> = Vec::new();
            let mut ok = true;
            for _ in 0..4 {
                let l: Vec<f64> = (0..m * m).map(|_| rng.f64() * 2.0 - 1.0).collect();
                if det3m(&l).abs() < 0.1 {
                    ok = false;
                }
                ls.push(l);
            }
            if !ok {
                continue;
            }
            draws += 1;
            for i in 0..4 {
                for j in 0..4 {
                    if i == j {
                        continue;
                    }
                    let x = atlas_block(&ops8[i], &ops8[j], &c8, n8);
                    // 元の whitened 特異値
                    let wi = inv_sqrt_sym(&gs[i], m);
                    let wj = inv_sqrt_sym(&gs[j], m);
                    let xw = matmul4(&wi, &matmul4(&x, &wj, m), m);
                    let sv0 = singvals(&xw, m, m);
                    // 再結合基底: X' = L_i X L_jᵀ, G' = L G Lᵀ
                    let xp = lxl(&ls[i], &x, &ls[j], m);
                    let gi = lxl(&ls[i], &gs[i], &ls[i], m);
                    let gj = lxl(&ls[j], &gs[j], &ls[j], m);
                    let wip = inv_sqrt_sym(&gi, m);
                    let wjp = inv_sqrt_sym(&gj, m);
                    let xwp = matmul4(&wip, &matmul4(&xp, &wjp, m), m);
                    let sv1 = singvals(&xwp, m, m);
                    for (a, b) in sv0.iter().zip(sv1.iter()) {
                        worst = worst.max((a - b).abs());
                    }
                    // ノルム 3 種
                    worst = worst
                        .max((frob(&xw) - frob(&xwp)).abs())
                        .max((sv0[0] - sv1[0]).abs())
                        .max((sv0.iter().sum::<f64>() - sv1.iter().sum::<f64>()).abs());
                }
            }
        }
        check(
            "[A1] whitening 不変性: 特異スペクトル/作用素/Frobenius/核ノルムが任意の可逆再結合 L_i で不変 ≤ 1e-9",
            worst <= 1e-9,
            format!("max 乖離 {:.2e} (乱択 8 draws — unitary に限らない)", worst),
        );
    }

    // ---- [A2] 動的 atlas: 二重経路照合 + 時間積分ノルムの不変性 ----
    {
        let (vals, vecs) = jacobi_eigh(&h8, n8);
        // 二重経路: response_r (固有系位相) vs 微小 t の直接展開 i Tr(C[M+it[h,M], M']) + O(t²)
        let mut worst_dual: f64 = 0.0;
        let tsmall = 1e-5;
        for i in 0..2 {
            for j in 2..4 {
                let ma = &ops8[i][0];
                let mb = &ops8[j][0];
                let r = response_r(&vals, &vecs, ma, mb, &c8, n8, tsmall);
                // 展開: M_A(t) ≈ M + it[h,M] → i Tr(C[M_A(t),M_B]) の実部
                //   ≈ −t Tr(C [[h,M_A], M_B])
                let hm = matmul(&h8, ma, n8);
                let mh = matmul(ma, &h8, n8);
                let comm: Vec<f64> = (0..n8 * n8).map(|k| hm[k] - mh[k]).collect();
                let t1 = matmul(&comm, mb, n8);
                let t2 = matmul(mb, &comm, n8);
                let mut s = 0.0;
                for p in 0..n8 {
                    for q in 0..n8 {
                        s += c8[p * n8 + q] * (t1[q * n8 + p] - t2[q * n8 + p]);
                    }
                }
                let expect = -tsmall * s;
                worst_dual = worst_dual.max((r - expect).abs() / tsmall);
            }
        }
        // 時間積分ノルム ∫₀^T ‖R̂‖_F dt の再結合不変性 (1 draw)
        let mut rng = Rng::new(3142);
        let m = 3usize;
        let mut ls: Vec<Vec<f64>> = Vec::new();
        loop {
            ls.clear();
            let mut ok = true;
            for _ in 0..4 {
                let l: Vec<f64> = (0..m * m).map(|_| rng.f64() * 2.0 - 1.0).collect();
                if det3m(&l).abs() < 0.1 {
                    ok = false;
                }
                ls.push(l);
            }
            if ok {
                break;
            }
        }
        let gs: Vec<Vec<f64>> = (0..4)
            .map(|i| atlas_block(&ops8[i], &ops8[i], &c8, n8))
            .collect();
        let (i, j) = (0usize, 2usize);
        let integ = |use_l: bool| -> f64 {
            let mut acc = 0.0;
            for k in 1..=20 {
                let t = 0.1 * k as f64;
                let mut r = vec![0.0; m * m];
                for a in 0..m {
                    for b in 0..m {
                        r[a * m + b] =
                            response_r(&vals, &vecs, &ops8[i][a], &ops8[j][b], &c8, n8, t);
                    }
                }
                let (rw, giw, gjw) = if use_l {
                    (
                        lxl(&ls[i], &r, &ls[j], m),
                        lxl(&ls[i], &gs[i], &ls[i], m),
                        lxl(&ls[j], &gs[j], &ls[j], m),
                    )
                } else {
                    (r.clone(), gs[i].clone(), gs[j].clone())
                };
                let wi = inv_sqrt_sym(&giw, m);
                let wj = inv_sqrt_sym(&gjw, m);
                let rhat = matmul4(&wi, &matmul4(&rw, &wj, m), m);
                acc += frob(&rhat) * 0.1;
            }
            acc
        };
        let (i0, i1) = (integ(false), integ(true));
        check(
            "[A2] 動的 atlas: 短時間展開との二重経路 ≤ 1e-6 + 時間積分ノルム ∫‖R̂‖dt の再結合不変 ≤ 1e-9",
            worst_dual <= 1e-6 && (i0 - i1).abs() <= 1e-9,
            format!(
                "二重経路乖離 {:.2e} / ∫‖R̂‖dt = {:.6} (再結合後乖離 {:.2e})",
                worst_dual,
                i0,
                (i0 - i1).abs()
            ),
        );
    }

    // ---- [A3] B3/B4 = atlas の成分/圧縮 (ring12, 1 サイトノード) ----
    {
        let n = 12usize;
        let mut h = vec![0.0; n * n];
        for k in 0..n {
            let t = 1.0 + 0.3 * (std::f64::consts::TAU * k as f64 / n as f64 + 0.5).sin();
            h[k * n + (k + 1) % n] = -t;
            h[((k + 1) % n) * n + k] = -t;
        }
        let c = gibbs_c(&h, n, beta);
        // B3 核 = |X^{nn}|: X(n_i, n_j) = −C_ij² → |X| = C_ij² (B3-COV そのもの)
        let mut worst_b3: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let mut mi = vec![0.0; n * n];
                mi[i * n + i] = 1.0;
                let mut mj = vec![0.0; n * n];
                mj[j * n + j] = 1.0;
                let x = cov_static_mat(&mi, &mj, &c, n);
                worst_b3 = worst_b3.max((x.abs() - c[i * n + j] * c[i * n + j]).abs());
            }
        }
        // B4 = 動的成分の到着時刻圧縮: |R^{nn}_{0j}(t)| の初交差 τ_j と重みつき最短路の単調性
        let (vals, vecs) = jacobi_eigh(&h, n);
        let mut m0 = vec![0.0; n * n];
        m0[0] = 1.0;
        let mut tau = vec![f64::NAN; n];
        for k in 1..=800 {
            let t = 0.01 * k as f64;
            for j in 1..n {
                if !tau[j].is_nan() {
                    continue;
                }
                let mut mj = vec![0.0; n * n];
                mj[j * n + j] = 1.0;
                let r = response_r(&vals, &vecs, &m0, &mj, &c, n, t);
                if r.abs() > 1e-3 {
                    tau[j] = t;
                }
            }
            if (1..n).all(|j| !tau[j].is_nan()) {
                break;
            }
        }
        let dist = {
            // Dijkstra (v31.3 と同型)
            let mut dist = vec![f64::INFINITY; n];
            dist[0] = 0.0;
            let mut done = vec![false; n];
            for _ in 0..n {
                let mut u = usize::MAX;
                let mut best = f64::INFINITY;
                for v in 0..n {
                    if !done[v] && dist[v] < best {
                        best = dist[v];
                        u = v;
                    }
                }
                if u == usize::MAX {
                    break;
                }
                done[u] = true;
                for v in 0..n {
                    let w = h[u * n + v].abs();
                    if v != u && w > 1e-12 {
                        let nd = dist[u] + 1.0 / w;
                        if nd < dist[v] {
                            dist[v] = nd;
                        }
                    }
                }
            }
            dist
        };
        let idx: Vec<usize> = (1..n).collect();
        let tv: Vec<f64> = idx.iter().map(|&j| tau[j]).collect();
        let dv: Vec<f64> = idx.iter().map(|&j| dist[j]).collect();
        let rho = spearman(&tv, &dv);
        check(
            "[A3] B3 = atlas の (n,n) 成分 (|X| = C² 厳密) / B4 = 動的 (n,n) 成分の到着時刻圧縮 (Spearman ≥ 0.8)",
            worst_b3 <= 1e-14 && rho >= 0.8,
            format!("B3 成分乖離 {:.2e} / Spearman(τ_R, dist) = {:.3}", worst_b3, rho),
        );
    }

    // ---- [A4] rank 欠損 Gram → 支持証明書 (無条件擬似逆の禁止) ----
    {
        // 基底に厳密な線形従属を仕込む: O⁴ = O¹ + O² → G は rank 3 / nullspace 1
        let mut ops5 = ops8[0].clone();
        let dup: Vec<f64> = ops5[0]
            .iter()
            .zip(ops5[1].iter())
            .map(|(a, b)| a + b)
            .collect();
        ops5.push(dup);
        let m5 = 4usize;
        let g5 = {
            let mut g = vec![0.0; m5 * m5];
            for a in 0..m5 {
                for b in 0..m5 {
                    g[a * m5 + b] = cov_static_mat(&ops5[a], &ops5[b], &c8, n8);
                }
            }
            g
        };
        let (ev, evec) = jacobi_eigh(&g5, m5);
        let thr = 1e-10 * ev.iter().fold(0.0f64, |m, &e| m.max(e.abs()));
        let rank = ev.iter().filter(|&&e| e > thr).count();
        let cert = ObservableSupportCertificate {
            rank,
            threshold: thr,
            nullspace_dim: m5 - rank,
        };
        // 支持制限 whitening: 正固有空間に射影してから逆平方根 (擬似逆の無条件使用を禁止 —
        // 証明書の rank でのみ動く)
        let x05 = {
            let mut x = vec![0.0; m5 * 3];
            for a in 0..m5 {
                for b in 0..3 {
                    x[a * 3 + b] = cov_static_mat(&ops5[a], &ops8[2][b], &c8, n8);
                }
            }
            x
        };
        // 支持基底 (rank 本の固有ベクトル) に写像: X_s = Sᵀ X, G_s = diag(ev_pos)
        let mut sv_support = {
            let mut xs = vec![0.0; rank * 3];
            let mut row = 0usize;
            let mut gs_half = vec![0.0; rank];
            for e in 0..m5 {
                if ev[e] <= thr {
                    continue;
                }
                gs_half[row] = 1.0 / ev[e].sqrt();
                for b in 0..3 {
                    let mut s = 0.0;
                    for a in 0..m5 {
                        s += evec[e * m5 + a] * x05[a * 3 + b];
                    }
                    xs[row * 3 + b] = s;
                }
                row += 1;
            }
            // whitening (行側): diag(1/√ev) Xs、列側は ops8[2] の Gram
            let g2 = atlas_block(&ops8[2], &ops8[2], &c8, n8);
            let w2 = inv_sqrt_sym(&g2, 3);
            let mut xw = vec![0.0; rank * 3];
            for r in 0..rank {
                for b in 0..3 {
                    let mut s = 0.0;
                    for q in 0..3 {
                        s += xs[r * 3 + q] * w2[q * 3 + b];
                    }
                    xw[r * 3 + b] = gs_half[r] * s;
                }
            }
            singvals(&xw, rank, 3)
        };
        // 清浄基底 (4 本) の whitened 特異値と一致するはず
        let sv_clean = {
            let g0 = atlas_block(&ops8[0], &ops8[0], &c8, n8);
            let x = atlas_block(&ops8[0], &ops8[2], &c8, n8);
            let g2 = atlas_block(&ops8[2], &ops8[2], &c8, n8);
            let w0 = inv_sqrt_sym(&g0, 3);
            let w2 = inv_sqrt_sym(&g2, 3);
            let xw = matmul4(&w0, &matmul4(&x, &w2, 3), 3);
            singvals(&xw, 3, 3)
        };
        sv_support.resize(3, 0.0);
        let mut worst: f64 = 0.0;
        for (a, b) in sv_support.iter().zip(sv_clean.iter()) {
            worst = worst.max((a - b).abs());
        }
        check(
            "[A4] rank 欠損 Gram (従属演算子): ObservableSupportCertificate {rank 3, null 1} + 支持制限 whitening = 清浄基底の不変量 ≤ 1e-9",
            cert.rank == 3 && cert.nullspace_dim == 1 && worst <= 1e-9,
            format!(
                "rank {} / nullspace {} / 支持 vs 清浄の特異値乖離 {:.2e} (無条件擬似逆は不使用)",
                cert.rank, cert.nullspace_dim, worst
            ),
        );
    }

    // ---- Part B: factorization no-go (ring12) ----
    let n = 12usize;
    let h12 = {
        let mut h = vec![0.0; n * n];
        for k in 0..n {
            let t = 1.0 + 0.3 * (std::f64::consts::TAU * k as f64 / n as f64 + 0.5).sin();
            h[k * n + (k + 1) % n] = -t;
            h[((k + 1) % n) * n + k] = -t;
        }
        h
    };
    let c12 = gibbs_c(&h12, n, beta);
    let k12 = matfun_sym(&c12, n, |x| ((1.0 - x) / x).ln());

    // ---- [A5] mode 対角化 no-go + 疎性基準の負制御 ----
    {
        let (_, vecs) = jacobi_eigh(&c12, n);
        // mode 基底での C: 厳密対角 → モード間 B3 核は 0
        let mut worst_off: f64 = 0.0;
        for a in 0..n {
            for b in 0..n {
                if a == b {
                    continue;
                }
                let mut s = 0.0;
                for i in 0..n {
                    for j in 0..n {
                        s += vecs[a * n + i] * c12[i * n + j] * vecs[b * n + j];
                    }
                }
                worst_off = worst_off.max(s.abs());
            }
        }
        // 疎性基準の負制御: nnz(K_site) vs nnz(K_mode = 対角)
        let nnz = |m: &[f64], tol: f64| m.iter().filter(|v| v.abs() > tol).count();
        let nnz_site = nnz(&k12, 1e-8);
        let nnz_mode = n; // 対角のみ (worst_off ≤ 1e-12 が上で機械確認される)
        check(
            "[A5] mode factorization no-go: モード間静的相関は厳密 0 (幾何なし) + 疎性基準は mode 基底 (nnz 最小) を選ぶ — state-only 選択の負制御",
            worst_off <= 1e-12 && nnz_mode < nnz_site,
            format!(
                "モード間 |C̃| max {:.2e} / nnz: site {} vs mode {} — 「最も疎な基底」は幾何を自明化する",
                worst_off, nnz_site, nnz_mode
            ),
        );
    }

    // ---- [A6] 同一状態の 3 因子分解 → 3 つの異なる幾何 ----
    {
        // (i) site: B3 支持 = ring 12 辺
        let sup_site = {
            let mut w = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    if i != j {
                        w[i * n + j] = c12[i * n + j] * c12[i * n + j];
                    }
                }
            }
            support_from_weights(&w, n)
        };
        // (ii) mode: 支持なし ([A5])
        // (iii) pair 回転: W = ⊕ R(π/4) on (0,1),(2,3),… → C' = Wᵀ C W
        let cp = {
            let th = std::f64::consts::FRAC_PI_4;
            let (co, si) = (th.cos(), th.sin());
            let mut w = vec![0.0; n * n];
            for b in 0..n / 2 {
                let (i, j) = (2 * b, 2 * b + 1);
                w[i * n + i] = co;
                w[i * n + j] = -si;
                w[j * n + i] = si;
                w[j * n + j] = co;
            }
            let t1 = matmul(&transpose(&w, n), &c12, n);
            matmul(&t1, &w, n)
        };
        let sup_rot = {
            let mut wgt = vec![0.0; n * n];
            for i in 0..n {
                for j in 0..n {
                    if i != j {
                        wgt[i * n + j] = cp[i * n + j] * cp[i * n + j];
                    }
                }
            }
            support_from_weights(&wgt, n)
        };
        let ring: Vec<(usize, usize)> = (0..n)
            .map(|k| {
                let (a, b) = (k, (k + 1) % n);
                (a.min(b), a.max(b))
            })
            .collect();
        let site_is_ring =
            sup_site.len() == 12 && ring.iter().all(|e| sup_site.contains(e));
        let rot_differs = sup_rot != sup_site;
        check(
            "[A6] 同一状態の 3 因子分解 → 3 幾何 (site = ring 12 辺 / mode = 幾何なし / pair 回転 = 別の支持) — 一意選定は state-only では不可能",
            site_is_ring && rot_differs,
            format!(
                "site 支持 {} 辺 (= ring) / 回転基底 支持 {} 辺 (≠ site) / mode 0 辺 — 全て同じ大域状態の「正しい読み」",
                sup_site.len(),
                sup_rot.len()
            ),
        );
    }

    // ---- [A7] OperationalAlgebra が幾何を選ぶ ----
    {
        // probe 応答 (v31.2 曲率則) の代数計算: (n̈⁺−n̈⁻) = −2ε Tr(P_j [h,[h,P_i]])。
        // 二重交換子トレースを完全な行列演算で計算し (簡約を使わない)、
        // site 代数では 4ε·h_ji²・mode 代数では 0 になることを検査する。
        let eps = 0.3;
        let curvature_diff = |proj: &dyn Fn(usize) -> Vec<f64>, i: usize, j: usize| -> f64 {
            let pi = proj(i);
            let pj = proj(j);
            let hp = matmul(&h12, &pi, n);
            let ph = matmul(&pi, &h12, n);
            let c1: Vec<f64> = (0..n * n).map(|k| hp[k] - ph[k]) .collect();
            let hc = matmul(&h12, &c1, n);
            let ch = matmul(&c1, &h12, n);
            let dd: Vec<f64> = (0..n * n).map(|k| hc[k] - ch[k]).collect();
            let t = matmul(&pj, &dd, n);
            -2.0 * eps * (0..n).map(|a| t[a * n + a]).sum::<f64>()
        };
        // site 代数: P_i = e_i e_iᵀ → 応答 = 4ε h_ji²
        let site_proj = |i: usize| -> Vec<f64> {
            let mut p = vec![0.0; n * n];
            p[i * n + i] = 1.0;
            p
        };
        let mut worst_site: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                if i == j {
                    continue;
                }
                let lhs = curvature_diff(&site_proj, i, j);
                let rhs = 4.0 * eps * h12[j * n + i] * h12[j * n + i];
                worst_site = worst_site.max((lhs - rhs).abs());
            }
        }
        // mode 代数: P_m = v_m v_mᵀ → 応答 = 0 (h は自分の固有基底で対角)
        let (_, hvec) = jacobi_eigh(&h12, n);
        let mode_proj = |m: usize| -> Vec<f64> {
            let mut p = vec![0.0; n * n];
            for a in 0..n {
                for b in 0..n {
                    p[a * n + b] = hvec[m * n + a] * hvec[m * n + b];
                }
            }
            p
        };
        let mut worst_mode: f64 = 0.0;
        for a in 0..n {
            for b in 0..n {
                if a != b {
                    worst_mode = worst_mode.max(curvature_diff(&mode_proj, a, b).abs());
                }
            }
        }
        // OperationalAlgebra の値を構成 (契約型の実演)
        let site_algebra = OperationalAlgebra {
            preparations: vec!["C± = I/2 ± εP_i (site 射影)"],
            interventions: vec!["局所バイアス ±ε (v31.2)"],
            measurements: vec!["ノード密度時系列 n_j(t)"],
            compatibility: vec!["[P_i, P_j] = 0 (直交射影)"],
        };
        let mode_algebra = OperationalAlgebra {
            preparations: vec!["C± = I/2 ± εP_m (mode 射影)"],
            interventions: vec!["モードバイアス ±ε"],
            measurements: vec!["モード占有時系列 n_m(t)"],
            compatibility: vec!["[P_m, P_m'] = 0"],
        };
        println!(
            "     site 代数 {:?} → ring / mode 代数 {:?} → 幾何なし",
            site_algebra.preparations, mode_algebra.preparations
        );
        check(
            "[A7] OperationalAlgebra が幾何を選ぶ: site-local probe → ring (厳密) / mode-local probe → 応答 0 (厳密) — 因子分解は操作代数が運ぶ (状態ではない)",
            worst_site <= 1e-14 && worst_mode <= 1e-12,
            format!(
                "site 応答 = 隣接重み² (乖離 {:.2e}) / mode 間応答 max {:.2e}",
                worst_site, worst_mode
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "atlas と no-go が確定 — whitened 不変量は演算子基底の可逆再結合に不変で B3/B4 はその成分/圧縮。同一状態は複数の因子分解で異なる幾何を返し (mode 基底では幾何が消える)、state-only の選択基準 (疎性) は自明幾何を選ぶ。FactorizationGivenObservables には OperationalAlgebra (準備・介入・測定・両立性) の入力が必要 — v29.5 [C5] の空隙が定理化された"
        } else {
            "**atlas/no-go の破れ**"
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

// ---------------------------------------------------------------- 補助 (4×4 と支持)

fn det3m(l: &[f64]) -> f64 {
    // 3×3 行列式 (可逆性の粗い確認用)
    l[0] * (l[4] * l[8] - l[5] * l[7]) - l[1] * (l[3] * l[8] - l[5] * l[6])
        + l[2] * (l[3] * l[7] - l[4] * l[6])
}

fn matmul4(a: &[f64], b: &[f64], m: usize) -> Vec<f64> {
    let mut c = vec![0.0; m * m];
    for i in 0..m {
        for k in 0..m {
            let aik = a[i * m + k];
            for j in 0..m {
                c[i * m + j] += aik * b[k * m + j];
            }
        }
    }
    c
}

/// L X Rᵀ (m×m)
fn lxl(l: &[f64], x: &[f64], r: &[f64], m: usize) -> Vec<f64> {
    let mut t = vec![0.0; m * m];
    for i in 0..m {
        for j in 0..m {
            let mut s = 0.0;
            for p in 0..m {
                for q in 0..m {
                    s += l[i * m + p] * x[p * m + q] * r[j * m + q];
                }
            }
            t[i * m + j] = s;
        }
    }
    t
}

fn transpose(a: &[f64], n: usize) -> Vec<f64> {
    let mut t = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            t[j * n + i] = a[i * n + j];
        }
    }
    t
}

fn spearman(x: &[f64], y: &[f64]) -> f64 {
    let m = x.len();
    let rank = |v: &[f64]| -> Vec<f64> {
        let mut idx: Vec<usize> = (0..v.len()).collect();
        idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap());
        let mut r = vec![0.0; v.len()];
        for (k, &i) in idx.iter().enumerate() {
            r[i] = k as f64;
        }
        r
    };
    let rx = rank(x);
    let ry = rank(y);
    let mean = (m as f64 - 1.0) / 2.0;
    let mut num = 0.0;
    let mut dx = 0.0;
    let mut dy = 0.0;
    for k in 0..m {
        num += (rx[k] - mean) * (ry[k] - mean);
        dx += (rx[k] - mean).powi(2);
        dy += (ry[k] - mean).powi(2);
    }
    num / (dx * dy).sqrt()
}

fn support_from_weights(w: &[f64], n: usize) -> Vec<(usize, usize)> {
    let mut vals: Vec<f64> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            vals.push(w[i * n + j].abs().max(1e-300));
        }
    }
    let mut sorted = vals.clone();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
    // 最終規則 (v31.6 で凍結): スケールガード窓 (max·1e-3) の**内側**の最大対数段差が
    // 有意 (≥ ln 3) ならそこで切る。有意な窓内段差がなければ窓内は単一クラスタ =
    // 全て辺として窓境界で切る。旧規則 2 種の故障を両方閉じる:
    //   (i) ガードなし: f64 尾部の発散段差を拾う (v31.3 で訂正済み)
    //   (ii) 跨ぎ段差を本命にする: 多段階系で物理段差を尾部跨ぎが上書き
    let guard = sorted[0] * 1e-3;
    let mut cut: Option<usize> = None;
    let mut best_gap = 0.0;
    for k in 0..sorted.len() - 1 {
        if sorted[k + 1] < guard {
            break; // 窓内段差のみ (両端 ≥ guard)
        }
        let gap = (sorted[k] / sorted[k + 1]).ln();
        if gap > best_gap {
            best_gap = gap;
            cut = Some(k);
        }
    }
    let thr = match cut {
        Some(k) if best_gap >= 3.0f64.ln() => (sorted[k] * sorted[k + 1]).sqrt(),
        _ => guard,
    };
    let mut e = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if w[i * n + j].abs() > thr {
                e.push((i, j));
            }
        }
    }
    e
}
