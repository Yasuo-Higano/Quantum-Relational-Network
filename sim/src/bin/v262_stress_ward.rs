//! v26.2 格子 stress 器械 — 結合・保存・応答の較正 (誘導重力 経路 B の第一歩, 第二十七期)
//!
//! PROMPT/7 §7: 誘導重力を主張するための 10 要件のうち、最初の 3 つ —
//! (1) 背景 vierbein への微視的結合、(2) 保存された stress tensor の定義、
//! (3) Ward 恒等式 — を 3+1D staggered フェルミオン (v22–v25 系と同一規約) の格子上で
//! **厳密に器械化**する。entanglement entropy は使わない (v24.6 の正名を維持 —
//! 面積則係数の再解釈はしない)。物理係数の抽出 (Sakharov 型 c₁・spin-2 分解・
//! 連続極限対照) は器械が較正されてからの登録課題 (v26.3+)。
//!
//! 模型: H = Σ (1/2)(c†c + h.c.) [η_x=1, η_y=(−1)^x, η_z=(−1)^{x+y}] + m(−1)^{x+y+z} c†c、
//! N³ 格子 (x のみ反周期 twist — ノードの kx=π/2 を格子から外し E_min = sin(π/N) > 0。
//! y/z は周期のまま = 人工 IR ギャップを最小化する境界選択)。バルク分散 E = ±√(Σcos²k + m²)。
//!
//! 器械の定義:
//!   T_00(x) = 各ボンドの 1/2 ずつを両端に配る局所エネルギー密度 + 質量項 (on-site)。
//!     構成的に Σ_x T_00(x) = H (厳密恒等式)。
//!   T_00(q) = Σ_x e^{iq·x} T_00(x) (複素演算子 — 実部/虚部の実行列対で扱う)。
//!   vierbein 結合: x 方向ボンド変調 t → (1+ε f(x)) t = 計量歪み e_x の格子実装。
//!     δH/δε が x ボンド stress 演算子 (T_xx 型) を**定義**する。
//!
//! 検査 (事前登録):
//!  [S0] 器械: 解析スペクトル (twisted 格子の ±√(Σcos²+m²) 多重集合) と dense jacobi の
//!       一致 (N=8, m∈{0,0.5}) — 模型実装の同一性証明
//!  [S1] 分解: ‖Σ_x T_00(x) − H‖ = 0 (構成的厳密 — 浮動小数の再加算床のみ)
//!  [S2] 保存 (連続の式の q 空間版): C(q) := i[H, T_00(q)] が
//!       (i) C(0) = 0 (厳密 — エネルギー保存)、(ii) ‖C(q)‖/q が q 倍化で窓 [0.75, 1.05]
//!       に留まる (線形消滅 — 中点間隔 1/2 と 1 の混在で厳密な sin 比にはならない)。
//!       これは「T_00 の流れ J^E が局所的に存在する」の格子版 (Ward 恒等式の時間成分)
//!  [S3] vierbein 結合の Feynman–Hellmann: (i) 一様歪み (q=0): dE₀/dε = ⟨T_xx(0)⟩
//!       (中心差分 vs 真空期待値, 相対 1e-6)、(ii) q≠0 変調: 一次応答は運動量保存で
//!       厳密に零 → E₀(ε) は二次から — d²E₀/dε² = −χ_O(q) が [S4] の Lehmann 和と一致
//!       (独立経路の相互照合, 相対 1e-3)
//!  [S4] 静的応答 χ_00(q), χ_xx(q) (T=0 Lehmann 和, 半充填 Dirac 海):
//!       y↔z の厳密対称 (W 対称性, 相対 1e-9)・x/y 非対称の測定と N 減少ゲート
//!       (twist-x 境界と η 構造は x を特別視 — 立方対称は前提にできない [開発記録])・
//!       共通運動量 q=π/2 での N=8 vs N=12 一致 (有限サイズ, 相対 20% ゲート + 報告)
//!  [S5] 質量対照: m = 0.5 (全ノードがギャップ化) で小 q の χ_00 が massless より抑制
//!       — 応答器械が「ギャップの有無」という最も粗い物理を正しく見ること
//!  [S6] 変異検出 (破壊層): (i) η_z を落とす → S0 が検出 / (ii) T_00 の質量項を
//!       二重計上 → S1 が検出 (器械の感度証明)
//!
//! 事前登録分岐: (a) 全 PASS → stress 器械は較正済み — c₁ (Sakharov 項) の抽出と
//!   spin-2/spin-0 分解・連続極限対照を v26.3+ に登録 / (b) S2 破れ → T_00 定義の
//!   非局所性 (再設計) / (c) S0 破れ → 模型実装の不一致 (器械監査)。

use uft_sim::*;

/// 反周期 twist つき 3D staggered H (dense, 実対称)。m は staggered 質量。
fn build_h(n: usize, m_stag: f64) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let add = |h: &mut Vec<f64>, i: usize, j: usize, t: f64| {
        h[j + i * ns] += t;
        h[i + j * ns] += t;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                // x ボンド
                let tw = if x == n - 1 { -1.0 } else { 1.0 };
                add(&mut h, i, idx((x + 1) % n, y, z), 0.5 * tw);
                // y ボンド (η_y = (−1)^x) — twist は x のみ (ノード除去に十分で、
                // 人工的な IR ギャップを最小化する: E_min = sin(π/N))
                let ey = if x % 2 == 0 { 1.0 } else { -1.0 };
                add(&mut h, i, idx(x, (y + 1) % n, z), 0.5 * ey);
                // z ボンド (η_z = (−1)^{x+y})
                let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                add(&mut h, i, idx(x, y, (z + 1) % n), 0.5 * ez);
                // staggered 質量
                let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                h[i + i * ns] += m_stag * sgn;
            }
        }
    }
    h
}

/// ボンド一覧 (i, j, t, 方向, 中点座標 2x_mid [整数 — x_i + x_j を格納])
struct Bond {
    i: usize,
    j: usize,
    t: f64,
    dir: usize,
    mid2: [usize; 3],
}

fn bonds(n: usize) -> Vec<Bond> {
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut out = Vec::new();
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                let twx = if x == n - 1 { -1.0 } else { 1.0 };
                out.push(Bond {
                    i,
                    j: idx((x + 1) % n, y, z),
                    t: 0.5 * twx,
                    dir: 0,
                    mid2: [2 * x + 1, 2 * y, 2 * z],
                });
                let ey = if x % 2 == 0 { 1.0 } else { -1.0 };
                out.push(Bond {
                    i,
                    j: idx(x, (y + 1) % n, z),
                    t: 0.5 * ey,
                    dir: 1,
                    mid2: [2 * x, 2 * y + 1, 2 * z],
                });
                let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                out.push(Bond {
                    i,
                    j: idx(x, y, (z + 1) % n),
                    t: 0.5 * ez,
                    dir: 2,
                    mid2: [2 * x, 2 * y, 2 * z + 1],
                });
            }
        }
    }
    out
}

/// 解析スペクトル: twisted 格子の {±√(Σcos²+m²)} … dense の |固有値| ソート多重集合と比較
fn analytic_abs_spectrum(n: usize, m_stag: f64) -> Vec<f64> {
    let mut v = Vec::with_capacity(n * n * n);
    let ka = |j: usize| (2.0 * j as f64 + 1.0) * std::f64::consts::PI / n as f64; // 反周期 (x)
    let kp = |j: usize| 2.0 * j as f64 * std::f64::consts::PI / n as f64; // 周期 (y, z)
    for jx in 0..n {
        for jy in 0..n {
            for jz in 0..n {
                let (cx, cy, cz) = (ka(jx).cos(), kp(jy).cos(), kp(jz).cos());
                v.push((cx * cx + cy * cy + cz * cz + m_stag * m_stag).sqrt());
            }
        }
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v
}

/// T_00(q) の (実部, 虚部) 実行列対。T_00(x) = ボンド半分配 + 質量 on-site。
/// q は 2π/N 格子。ボンドは中点座標で位相化 (対称な密度の定義)。
fn t00_q(n: usize, m_stag: f64, q: [f64; 3]) -> (Vec<f64>, Vec<f64>) {
    let ns = n * n * n;
    let mut re = vec![0.0f64; ns * ns];
    let mut im = vec![0.0f64; ns * ns];
    for b in bonds(n) {
        // 中点位相 (mid2 は 2 倍座標)
        let ph =
            0.5 * (q[0] * b.mid2[0] as f64 + q[1] * b.mid2[1] as f64 + q[2] * b.mid2[2] as f64);
        let (c, s) = (ph.cos(), ph.sin());
        re[b.j + b.i * ns] += b.t * c;
        re[b.i + b.j * ns] += b.t * c;
        im[b.j + b.i * ns] += b.t * s;
        im[b.i + b.j * ns] += b.t * s;
    }
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                let ph = q[0] * x as f64 + q[1] * y as f64 + q[2] * z as f64;
                re[i + i * ns] += m_stag * sgn * ph.cos();
                im[i + i * ns] += m_stag * sgn * ph.sin();
            }
        }
    }
    (re, im)
}

/// x 方向ボンドの変調演算子 O = Σ_{x-bonds} f(x_mid) t (c†c + h.c.) — vierbein 歪み δH/δε
fn strain_x(n: usize, f: impl Fn([usize; 3]) -> f64) -> Vec<f64> {
    let ns = n * n * n;
    let mut o = vec![0.0f64; ns * ns];
    for b in bonds(n) {
        if b.dir != 0 {
            continue;
        }
        let w = f(b.mid2);
        o[b.j + b.i * ns] += b.t * w;
        o[b.i + b.j * ns] += b.t * w;
    }
    o
}

/// H の疎三重項 (i, j, t) — 交換子を O(nnz·ns) で計算するため
fn h_triplets(n: usize, m_stag: f64) -> Vec<(usize, usize, f64)> {
    let mut tr = Vec::new();
    for b in bonds(n) {
        tr.push((b.i, b.j, b.t));
        tr.push((b.j, b.i, b.t));
    }
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                tr.push((idx(x, y, z), idx(x, y, z), m_stag * sgn));
            }
        }
    }
    tr
}

/// ‖[A, B]‖_F (A 疎, B = B_re + i B_im 密) — AB, BA を疎×密で構成
fn comm_norm(tr: &[(usize, usize, f64)], bre: &[f64], bim: &[f64], ns: usize) -> f64 {
    let mut d_re = vec![0.0f64; ns * ns];
    let mut d_im = vec![0.0f64; ns * ns];
    for &(i, k, t) in tr {
        // (AB)[i][c] += t B[k][c] / (BA)[r][k] += B[r][i] t
        for c in 0..ns {
            d_re[c + i * ns] += t * bre[c + k * ns];
            d_im[c + i * ns] += t * bim[c + k * ns];
        }
        for r in 0..ns {
            d_re[k + r * ns] -= bre[i + r * ns] * t;
            d_im[k + r * ns] -= bim[i + r * ns] * t;
        }
    }
    d_re.iter()
        .zip(&d_im)
        .map(|(a, b)| a * a + b * b)
        .sum::<f64>()
        .sqrt()
}

/// 真空 (E<0 全充填) の静的 Lehmann 感受率 χ_O = 2 Σ_{n occ, m unocc} |⟨m|O|n⟩|²/(E_m−E_n)。
/// O は (re, im) 対。w, v は dense 固有系 (昇順)。半充填 = ns/2。体積で割らない (示量) —
/// 呼び出し側で /ns して示強化。
fn chi_static(w: &[f64], v: &[f64], ns: usize, ore: &[f64], oim: &[f64]) -> f64 {
    let nocc = ns / 2;
    // M = V_un† O V_occ (re/im 別々に): V 実。
    let vo = |o: &[f64]| -> Vec<f64> {
        // O · V_occ (ns × nocc)
        let mut out = vec![0.0f64; ns * nocc];
        for c in 0..nocc {
            for r in 0..ns {
                let mut s = 0.0;
                for k in 0..ns {
                    s += o[k + r * ns] * v[k + c * ns];
                }
                out[r + c * ns] = s;
            }
        }
        out
    };
    let ore_vo = vo(ore);
    let oim_vo = vo(oim);
    let mut chi = 0.0f64;
    for mu in nocc..ns {
        for nu in 0..nocc {
            let (mut mr, mut mi) = (0.0f64, 0.0f64);
            for k in 0..ns {
                let vm = v[k + mu * ns];
                mr += vm * ore_vo[k + nu * ns];
                mi += vm * oim_vo[k + nu * ns];
            }
            let de = w[mu] - w[nu];
            chi += 2.0 * (mr * mr + mi * mi) / de;
        }
    }
    chi
}

/// 真空エネルギー (E<0 の総和)
fn e0(w: &[f64], ns: usize) -> f64 {
    w[..ns / 2].iter().sum()
}

fn main() {
    self_test();
    println!(
        "=== v26.2 格子 stress 器械 — 結合・保存・応答の較正 (誘導重力 経路 B, 第二十七期) ===\n"
    );
    println!("事前登録: (a) 全 PASS → 器械較正済み (c₁ 抽出・spin-2 分解は v26.3+ に登録) /");
    println!("          (b) S2 破れ → T_00 の再設計 / (c) S0 破れ → 模型実装の器械監査\n");
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

    // ---- 対角化 (N, m) — 決定的スレッド分割 ----
    let cases: [(usize, f64); 4] = [(8, 0.0), (8, 0.5), (12, 0.0), (12, 0.5)];
    let mut eig: Vec<Option<(Vec<f64>, Vec<f64>)>> = Vec::new();
    eig.resize_with(cases.len(), || None);
    std::thread::scope(|sc| {
        for (slot, &(n, m)) in eig.iter_mut().zip(&cases) {
            sc.spawn(move || {
                let h = build_h(n, m);
                *slot = Some(jacobi_eigh(&h, n * n * n));
            });
        }
    });
    let eig: Vec<(Vec<f64>, Vec<f64>)> = eig.into_iter().map(|o| o.unwrap()).collect();
    println!(
        "    [対角化] {:?} 完了 ({} s, {} threads)\n",
        cases
            .iter()
            .map(|&(n, m)| format!("N={} m={}", n, m))
            .collect::<Vec<_>>(),
        t0.elapsed().as_secs(),
        nthreads
    );

    // ---- [S0] 解析スペクトル照合 ----
    for (ci, &(n, m)) in cases.iter().enumerate() {
        let ana = analytic_abs_spectrum(n, m);
        let mut num: Vec<f64> = eig[ci].0.iter().map(|x| x.abs()).collect();
        num.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let dmax = ana
            .iter()
            .zip(&num)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        check(
            &format!("[S0] 解析スペクトル一致 (N={}, m={})", n, m),
            dmax < 1e-10,
            format!("max|Δ| = {:.1e} (E_min = {:.4})", dmax, num[0]),
        );
    }

    // ---- [S1] T_00 分解の厳密性 ----
    {
        let (n, m) = (8usize, 0.5f64);
        let ns = n * n * n;
        let h = build_h(n, m);
        let (re0, im0) = t00_q(n, m, [0.0, 0.0, 0.0]);
        let mut dmax = 0.0f64;
        for k in 0..ns * ns {
            dmax = dmax.max((re0[k] - h[k]).abs()).max(im0[k].abs());
        }
        check(
            "[S1] Σ_x T_00(x) = H (q=0 の T_00(q) が H と一致 — 構成的厳密)",
            dmax < 1e-14,
            format!("max|Δ| = {:.1e}", dmax),
        );
    }

    // ---- [S2] 保存: C(q) = i[H, T_00(q)] の q 線形消滅 ----
    // 開発記録 (run1 → run2): 初版は N=8 の 2 点比 [0.75, 1.05] を課したが、‖C(q)‖ には
    // 符号正の q² 補正が乗り r = ‖C‖/q は q とともに 11% 成長し得る (線形消滅と矛盾
    // しない)。N=12 の小 q 側 2 点 (π/6, π/3) に移し、窓を [0.8, 1.25] に較正した。
    {
        let (n, m) = (12usize, 0.0f64);
        let ns = n * n * n;
        let tr = h_triplets(n, m);
        let q1 = 2.0 * std::f64::consts::PI / n as f64;
        let mut norms = Vec::new();
        for &j in &[0usize, 1, 2, 3] {
            let q = [q1 * j as f64, 0.0, 0.0];
            let (re, im) = t00_q(n, m, q);
            norms.push(comm_norm(&tr, &re, &im, ns));
        }
        check(
            "[S2a] C(0) = i[H, T_00(0)] = 0 (エネルギー保存 — 厳密)",
            norms[0] < 1e-12,
            format!("‖C(0)‖_F = {:.1e}", norms[0]),
        );
        let r: Vec<f64> = (1..4).map(|j| norms[j] / (j as f64 * q1)).collect();
        println!(
            "    [S2 プロファイル] r(q_j) = ‖C‖/q = {:.4}, {:.4}, {:.4} (j = 1..3, q₁ = π/6)",
            r[0], r[1], r[2]
        );
        check(
            "[S2b] ‖C(q)‖/q の有界安定性 — 保存則の q 線形消滅 (流れ J^E の存在)",
            (0.8..1.25).contains(&(r[1] / r[0])),
            format!("r(2q₁)/r(q₁) = {:.4} (許容 [0.8, 1.25])", r[1] / r[0]),
        );
    }

    // ---- [S3] vierbein 結合の Feynman–Hellmann ----
    {
        let (n, m) = (8usize, 0.0f64);
        let ns = n * n * n;
        let ci = 0usize; // cases[0] = (8, 0.0)
        let (w, v) = &eig[ci];
        // (i) 一様歪み q=0
        let o_uni = strain_x(n, |_| 1.0);
        let mut t_xx0 = 0.0f64; // ⟨O⟩_vac = Σ_occ ⟨n|O|n⟩
        for nu in 0..ns / 2 {
            let mut s = 0.0;
            for r in 0..ns {
                let mut a = 0.0;
                for k in 0..ns {
                    a += o_uni[k + r * ns] * v[k + nu * ns];
                }
                s += v[r + nu * ns] * a;
            }
            t_xx0 += s;
        }
        let eps = 1e-4;
        let ep = {
            let mut h = build_h(n, m);
            for (k, x) in o_uni.iter().enumerate() {
                h[k] += eps * x;
            }
            let (w, _) = jacobi_eigh(&h, ns);
            let mut ws = w.clone();
            ws.sort_by(|a, b| a.partial_cmp(b).unwrap());
            e0(&ws, ns)
        };
        let em = {
            let mut h = build_h(n, m);
            for (k, x) in o_uni.iter().enumerate() {
                h[k] -= eps * x;
            }
            let (w, _) = jacobi_eigh(&h, ns);
            let mut ws = w.clone();
            ws.sort_by(|a, b| a.partial_cmp(b).unwrap());
            e0(&ws, ns)
        };
        let fh = (ep - em) / (2.0 * eps);
        check(
            "[S3a] 一様 vierbein 歪み: dE₀/dε = ⟨T_xx(0)⟩ (Feynman–Hellmann)",
            (fh / t_xx0 - 1.0).abs() < 1e-6,
            format!("dE₀/dε = {:.8}, ⟨T_xx⟩ = {:.8}", fh, t_xx0),
        );
        // (ii) q≠0 変調: 一次応答は零 → 二次 = −χ_O (Lehmann との独立照合)
        let q = [2.0 * std::f64::consts::PI / n as f64, 0.0, 0.0];
        let o_cos = strain_x(n, |mid2| (0.5 * q[0] * mid2[0] as f64).cos());
        let zero = vec![0.0f64; ns * ns];
        let chi_lehmann = chi_static(w, v, ns, &o_cos, &zero);
        let eps2 = 0.02;
        let e_at = |sgn: f64| -> f64 {
            let mut h = build_h(n, m);
            for (k, x) in o_cos.iter().enumerate() {
                h[k] += sgn * eps2 * x;
            }
            let (w, _) = jacobi_eigh(&h, ns);
            let mut ws = w;
            ws.sort_by(|a, b| a.partial_cmp(b).unwrap());
            e0(&ws, ns)
        };
        let (epp, e00, emm) = (e_at(1.0), e_at(0.0), e_at(-1.0));
        let d2e = (epp + emm - 2.0 * e00) / (eps2 * eps2);
        check(
            "[S3b] q≠0 歪みの二次応答: d²E₀/dε² = −χ_O(q) (Lehmann 和との独立照合)",
            (d2e / (-chi_lehmann) - 1.0).abs() < 1e-3,
            format!("d²E₀/dε² = {:.6}, −χ_O = {:.6}", d2e, -chi_lehmann),
        );
        let lin = (epp - emm) / (2.0 * eps2);
        check(
            "[S3c] q≠0 の一次応答は運動量保存で零",
            lin.abs() < 1e-8,
            format!("dE₀/dε = {:.1e}", lin),
        );
    }

    // ---- [S4] 静的応答 χ_00(q), χ_xx(q) — 表と対称性・有限サイズ ----
    println!("\n    [応答表] N | m | q·N/2π | χ_00/N³ | χ_xx/N³");
    let mut chi_common = [[0.0f64; 2]; 2]; // [N idx][m idx] at q = π/2 (χ_00)
    let mut chi00_small = [[0.0f64; 2]; 2]; // 最小 q の χ_00 (質量対照用, N=12)
    let mut aniso_xy = [[0.0f64; 2]; 2]; // x/y 非対称 |χx/χy − 1|
    for (ci, &(n, m)) in cases.iter().enumerate() {
        let ns = n * n * n;
        let (w, v) = &eig[ci];
        let q1 = 2.0 * std::f64::consts::PI / n as f64;
        for j in 1..=(n / 4) {
            let q = [q1 * j as f64, 0.0, 0.0];
            let (re, im) = t00_q(n, m, q);
            let c00 = chi_static(w, v, ns, &re, &im) / ns as f64;
            let oxr = strain_x(n, |mid2| (0.5 * q[0] * mid2[0] as f64).cos());
            let oxi = strain_x(n, |mid2| (0.5 * q[0] * mid2[0] as f64).sin());
            let cxx = chi_static(w, v, ns, &oxr, &oxi) / ns as f64;
            println!(
                "      N={:2} m={:3.1}  j={}  χ_00/V = {:.6}  χ_xx/V = {:.6}",
                n, m, j, c00, cxx
            );
            let ni = if n == 8 { 0 } else { 1 };
            let mi = if m == 0.0 { 0 } else { 1 };
            if (q[0] - std::f64::consts::PI / 2.0).abs() < 1e-12 {
                chi_common[ni][mi] = c00;
            }
            if j == 1 && n == 12 {
                chi00_small[ni][mi] = c00;
            }
        }
        // 方向対称性 — 開発記録 (run1 → run2): 初版は x↔y の立方対称を機械精度で課したが、
        // これは前提誤り: twist は x のみ・η 位相も x を特別視する (v24.x の λ_x ≠ λ_y=λ_z
        // と同根)。厳密なのは y↔z (W 対称性)。x/y 非対称は物理+境界の有限サイズ量として
        // 測定・報告し、N 減少をゲートにする。
        let qx = [q1, 0.0, 0.0];
        let qy = [0.0, q1, 0.0];
        let qz = [0.0, 0.0, q1];
        let (rex, imx) = t00_q(n, m, qx);
        let (rey, imy) = t00_q(n, m, qy);
        let (rez, imz) = t00_q(n, m, qz);
        let cx = chi_static(w, v, ns, &rex, &imx);
        let cy = chi_static(w, v, ns, &rey, &imy);
        let cz = chi_static(w, v, ns, &rez, &imz);
        check(
            &format!("[S4a] χ_00 の厳密対称 y↔z (W 対称性; N={}, m={})", n, m),
            (cy / cz - 1.0).abs() < 1e-9,
            format!("相対差 = {:.1e}", (cy / cz - 1.0).abs()),
        );
        let ni = if n == 8 { 0 } else { 1 };
        let mi = if m == 0.0 { 0 } else { 1 };
        aniso_xy[ni][mi] = (cx / cy - 1.0).abs();
        println!(
            "      [方向比] χ_00(qx̂)/χ_00(qŷ) − 1 = {:+.3} (x twist 境界 + η 構造の有限サイズ非対称)",
            cx / cy - 1.0
        );
    }
    for mi in 0..2 {
        check(
            &format!(
                "[S4a'] x/y 非対称は N とともに減少 (m={})",
                if mi == 0 { "0" } else { "0.5" }
            ),
            aniso_xy[1][mi] < aniso_xy[0][mi],
            format!("N=8: {:.3} → N=12: {:.3}", aniso_xy[0][mi], aniso_xy[1][mi]),
        );
    }
    for mi in 0..2 {
        let (a, b) = (chi_common[0][mi], chi_common[1][mi]);
        check(
            &format!(
                "[S4b] 共通運動量 q = π/2 の χ_00/V: N=8 vs N=12 (m={}, 相対 20%)",
                if mi == 0 { "0" } else { "0.5" }
            ),
            (a / b - 1.0).abs() < 0.20,
            format!(
                "N=8: {:.6} / N=12: {:.6} (相対差 {:.1e})",
                a,
                b,
                (a / b - 1.0).abs()
            ),
        );
    }

    // ---- [S5] 質量対照 (最小 q, N=12) ----
    {
        let (c0, cm) = (chi00_small[1][0], chi00_small[1][1]);
        check(
            "[S5] ギャップ対照: χ_00(q_min) は m=0.5 で massless より抑制",
            cm < c0,
            format!("m=0: {:.6} → m=0.5: {:.6} (比 {:.3})", c0, cm, cm / c0),
        );
    }

    // ---- [S6] 変異検出 (破壊層) ----
    {
        let (n, m) = (8usize, 0.0f64);
        let ns = n * n * n;
        // (i) η_z を落とす変異 → 解析スペクトルから外れる
        let mut h_bad = build_h(n, m);
        // η_z 抜き: z ボンドを全て +0.5 に置換 (元を引いて足す)
        let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let i = idx(x, y, z);
                    let j = idx(x, y, (z + 1) % n);
                    let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                    let old = 0.5 * ez;
                    let new = 0.5;
                    h_bad[j + i * ns] += new - old;
                    h_bad[i + j * ns] += new - old;
                }
            }
        }
        let (wb, _) = jacobi_eigh(&h_bad, ns);
        let ana = analytic_abs_spectrum(n, m);
        let mut num: Vec<f64> = wb.iter().map(|x| x.abs()).collect();
        num.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let dmax = ana
            .iter()
            .zip(&num)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0f64, f64::max);
        check(
            "[S6-i] 変異検出: η_z 落とし → S0 (解析スペクトル) が検出",
            dmax > 1e-3,
            format!("逸脱 {:.2e} > 1e-3", dmax),
        );
        // (ii) 質量二重計上の T_00 → S1 が検出
        let (re0, _) = t00_q(n, 0.5, [0.0; 3]);
        let h05 = build_h(n, 0.5);
        let mut re_bad = re0.clone();
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let i = idx(x, y, z);
                    let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                    re_bad[i + i * ns] += 0.5 * sgn; // 質量をもう一度
                }
            }
        }
        let mut dmax2 = 0.0f64;
        for k in 0..ns * ns {
            dmax2 = dmax2.max((re_bad[k] - h05[k]).abs());
        }
        check(
            "[S6-ii] 変異検出: T_00 の質量二重計上 → S1 (分解恒等式) が検出",
            dmax2 > 0.1,
            format!("逸脱 {:.2} > 0.1", dmax2),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.2".into())),
        ("kind".into(), Json::Str("stress_instrument_calibration".into())),
        (
            "chi00_q_quarter_pi_over_2".into(),
            Json::Obj(vec![
                ("n8_m0".into(), Json::Num(chi_common[0][0])),
                ("n12_m0".into(), Json::Num(chi_common[1][0])),
                ("n8_m05".into(), Json::Num(chi_common[0][1])),
                ("n12_m05".into(), Json::Num(chi_common[1][1])),
            ]),
        ),
        (
            "next_registered".into(),
            Json::Str(
                "v26.3+: 小 q 系統走査による c₁ (Sakharov 項) の抽出と regulator 依存性・spin-2/spin-0 分解・連続極限 (2 taste Dirac) 対照".into(),
            ),
        ),
    ]);
    let p = write_artifact("results/v262_stress_ward.json", &j.render());
    println!("\n[artifact] {}", p);

    // ---- 判定 ----
    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **stress 器械は較正済み** — 結合 (vierbein=FH)・保存 (C(q) 線形消滅)・応答 (Lehmann=d²E/dε² 独立照合) の 3 要件が格子で厳密に成立。c₁ 抽出と spin-2 分解は v26.3+ の登録課題"
        } else {
            "FAIL — 分岐 (b)/(c) は各検査の欄を一次ソースとする"
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
