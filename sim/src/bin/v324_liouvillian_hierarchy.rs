//! v32.4 Liouvillian 応答階層 — 一階/二階恒等式・中心を除く決定・H ↔ −H no-go (PROMPT/13 §3)
//!
//! 相関ホッピングを「曲率則が壊れるか」で採点する前に、**何がどの次数の応答に
//! 読めるか**を一般の有限次元系で恒等式として固定する。基準状態対
//! ρ± = ρ₀ ± εA_α と測定 B_β に対し (ρ(t) = e^{−iHt} ρ e^{iHt}):
//!
//!   R⁽¹⁾_βα = (ḃ_β⁺(0) − ḃ_β⁻(0))/(2ε) = −i·Tr(B_β [H, A_α])
//!   R⁽²⁾_βα = (b̈_β⁺(0) − b̈_β⁻(0))/(2ε) = −Tr(B_β [H,[H,A_α]]) = Tr([H,B_β][H,A_α])
//!   (Schrödinger 規約 ρ(t) = e^{−iHt} ρ e^{iHt} — v31.2 の C(t) と同一。
//!    PROMPT/13 の「i Tr」は逆向き発展の規約差で、階層の結論は不変)
//!
//! 能力階層 (機械化):
//!   一階は H に**線形** — 情報完全基底上で ad_H が決まれば H は**中心 (スカラー
//!   恒等項) を除いて**決まる。二階は ad_H² — **H ↔ −H を区別できない** (厳密)。
//!   密度だけに制限すると同値類はさらに大きい (**ring の磁束は密度二階に不可視・
//!   coherent 一階は分離**)。v31.2 の曲率則は本階層の**特殊化** (密度観測 ×
//!   積状態接ベクトル の R⁽²⁾ 対角) であることを厳密照合する。
//!
//!   [L0] 恒等式の測定照合: dim 8 の一般エルミート H — 時間発展 (複素 Taylor
//!        expm) + 4 次 stencil の測定値が R⁽¹⁾/R⁽²⁾ 公式と一致 (rel ≤ 1e-6)・
//!        ε ∈ {0.02, 0.10} で差分商が厳密一致 (恒等式 — 線形応答近似ではない)
//!   [L1] **一階応答は H を中心を除いて決める**: dim 4・情報完全 Pauli 基底 15 ×
//!        15 の R⁽¹⁾ から線形逆問題で Ĥ を復元 — traceless 部で rel ≤ 1e-10・
//!        H + cI は R⁽¹⁾ を厳密に変えない (中心の不可視性)
//!   [L2] **no-go の対**: (a) R⁽²⁾[H] = R⁽²⁾[−H] 厳密 (ad_H² の偶性) — 一階は
//!        符号を分離 (R⁽¹⁾[−H] = −R⁽¹⁾[H])。(b) **磁束 ring (3 モード JW)**:
//!        密度二階核は θ に不可視 (|h_ij|² = t²)・coherent 一階 (電流観測) は分離
//!   [L3] **保存則と PSD**: 数保存 H で B = N̂ の応答は全次数で厳密 0 (和則)・
//!        応答 Gram 核 K_αα' = Tr([H,A_α]†[H,A_α']) は PSD (固有値 ≥ −1e-12)
//!   [L4] **v31.2 曲率則の統一**: 3 モード自由鎖 (Fock dim 8) の many-body
//!        R⁽²⁾_ji = Tr([H, n_j][H, A_i]) (A_i = 積状態接ベクトル (n_i − 1/2)/4) が
//!        one-body の ‖P_j h P_i‖² = |h_ij|² と厳密一致 — 曲率則 = 本階層の対角
//!   [L5] 文書アンカー — uft-v32.4.md
//!
//! 実行: cargo run --release --bin v324_liouvillian_hierarchy

use std::fs;
use std::path::Path;
use uft_sim::operational_net::{cdag, cmul, commutator, hs_norm};
use uft_sim::{jacobi_eigh, C64};

// ---------------------------------------------------------------- 素子

fn pauli(which: char) -> Vec<C64> {
    let (o, l) = (C64::new(0.0, 0.0), C64::new(1.0, 0.0));
    match which {
        'I' => vec![l, o, o, l],
        'X' => vec![o, l, l, o],
        'Y' => vec![o, C64::new(0.0, -1.0), C64::new(0.0, 1.0), o],
        'Z' => vec![l, o, o, C64::new(-1.0, 0.0)],
        _ => panic!("未知の Pauli"),
    }
}

fn kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
    let n = na * nb;
    let mut out = vec![C64::new(0.0, 0.0); n * n];
    for i1 in 0..na {
        for j1 in 0..na {
            for i2 in 0..nb {
                for j2 in 0..nb {
                    out[(i1 * nb + i2) * n + (j1 * nb + j2)] = a[i1 * na + j1] * b[i2 * nb + j2];
                }
            }
        }
    }
    out
}

/// Pauli 文字列 (長さ任意) の行列
fn op_str(s: &str) -> Vec<C64> {
    let cs: Vec<char> = s.chars().collect();
    let mut m = pauli(cs[0]);
    let mut dim = 2;
    for &c in &cs[1..] {
        m = kron(&m, dim, &pauli(c), 2);
        dim *= 2;
    }
    m
}

fn ident(n: usize) -> Vec<C64> {
    let mut m = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

fn add_scaled(a: &mut [C64], b: &[C64], s: C64) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x = *x + s * *y;
    }
}

fn norm1(m: &[C64], n: usize) -> f64 {
    (0..n)
        .map(|i| (0..n).map(|j| m[i * n + j].abs()).sum::<f64>())
        .fold(0.0f64, f64::max)
}

/// exp(M) — スケーリング + Taylor (‖M‖ 小の領域用, 決定的)
fn expm(m: &[C64], n: usize) -> Vec<C64> {
    let mut k = 0u32;
    let mut nrm = norm1(m, n);
    while nrm > 0.1 {
        nrm *= 0.5;
        k += 1;
    }
    let scale = 1.0 / (1u64 << k) as f64;
    let ms: Vec<C64> = m.iter().map(|c| c.scale(scale)).collect();
    let mut out = ident(n);
    let mut term = ident(n);
    for j in 1..=14 {
        term = cmul(&term, &ms, n);
        let inv = 1.0 / (j as f64);
        for t in term.iter_mut() {
            *t = t.scale(inv);
        }
        for (o, t) in out.iter_mut().zip(term.iter()) {
            *o = *o + *t;
        }
    }
    for _ in 0..k {
        out = cmul(&out, &out, n);
    }
    out
}

/// b(t) = Re Tr(B ρ(t)), ρ(t) = e^{−iHt} ρ e^{iHt}
fn observe(h: &[C64], rho: &[C64], b: &[C64], t: f64, n: usize) -> f64 {
    let a: Vec<C64> = h.iter().map(|c| C64::new(c.im * t, -c.re * t)).collect(); // −iHt
    let u = expm(&a, n);
    let rt = cmul(&u, &cmul(rho, &cdag(&u, n), n), n);
    let mut s = C64::new(0.0, 0.0);
    for i in 0..n {
        for k in 0..n {
            s = s + b[i * n + k] * rt[k * n + i];
        }
    }
    s.re
}

/// 解析 R⁽¹⁾ = −i Tr(B[H,A]) = Im Tr(B[H,A]) (実数) —
/// Schrödinger 発展 ρ(t) = e^{−iHt} ρ e^{iHt} (v31.2 の C(t) と同じ規約) の一階
fn r1_exact(h: &[C64], b: &[C64], a: &[C64], n: usize) -> f64 {
    let c = commutator(h, a, n);
    let mut s = C64::new(0.0, 0.0);
    for i in 0..n {
        for k in 0..n {
            s = s + b[i * n + k] * c[k * n + i];
        }
    }
    s.im
}

/// 解析 R⁽²⁾ = Tr([H,B][H,A]) (実数)
fn r2_exact(h: &[C64], b: &[C64], a: &[C64], n: usize) -> f64 {
    let hb = commutator(h, b, n);
    let ha = commutator(h, a, n);
    let mut s = C64::new(0.0, 0.0);
    for i in 0..n {
        for k in 0..n {
            s = s + hb[i * n + k] * ha[k * n + i];
        }
    }
    s.re
}

/// 測定 lane: 4 次 stencil で (ḃ⁺−ḃ⁻)/(2ε), (b̈⁺−b̈⁻)/(2ε)
fn measure_r12(
    h: &[C64],
    rho0: &[C64],
    a: &[C64],
    b: &[C64],
    eps: f64,
    n: usize,
) -> (f64, f64) {
    let dt = 0.02 / norm1(h, n).max(1.0);
    let hh = dt / 2.0;
    let mut rp = rho0.to_vec();
    add_scaled(&mut rp, a, C64::new(eps, 0.0));
    let mut rm = rho0.to_vec();
    add_scaled(&mut rm, a, C64::new(-eps, 0.0));
    let f = |rho: &[C64], t: f64| observe(h, rho, b, t, n);
    let stencil = |rho: &[C64]| -> (f64, f64) {
        let (fm2, fm1, f0, f1, f2) = (
            f(rho, -2.0 * hh),
            f(rho, -hh),
            f(rho, 0.0),
            f(rho, hh),
            f(rho, 2.0 * hh),
        );
        let d1 = (fm2 - 8.0 * fm1 + 8.0 * f1 - f2) / (12.0 * hh);
        let d2 = (-fm2 + 16.0 * fm1 - 30.0 * f0 + 16.0 * f1 - f2) / (12.0 * hh * hh);
        (d1, d2)
    };
    let (d1p, d2p) = stencil(&rp);
    let (d1m, d2m) = stencil(&rm);
    ((d1p - d1m) / (2.0 * eps), (d2p - d2m) / (2.0 * eps))
}

/// エルミート結合の構築: Σ c_k · P_k (Pauli 文字列)
fn build_h(terms: &[(&str, f64)]) -> (Vec<C64>, usize) {
    let n = 1usize << terms[0].0.len();
    let mut h = vec![C64::new(0.0, 0.0); n * n];
    for (s, c) in terms {
        add_scaled(&mut h, &op_str(s), C64::new(*c, 0.0));
    }
    (h, n)
}

/// 15×15 実対称連立の Gauss 解 (部分 pivot)
fn gauss_solve(a: &mut [f64], b: &mut [f64], n: usize) -> Vec<f64> {
    for col in 0..n {
        let mut piv = col;
        for r in col + 1..n {
            if a[r * n + col].abs() > a[piv * n + col].abs() {
                piv = r;
            }
        }
        if piv != col {
            for k in 0..n {
                a.swap(col * n + k, piv * n + k);
            }
            b.swap(col, piv);
        }
        let d = a[col * n + col];
        for r in col + 1..n {
            let f = a[r * n + col] / d;
            for k in col..n {
                a[r * n + k] -= f * a[col * n + k];
            }
            b[r] -= f * b[col];
        }
    }
    let mut x = vec![0.0; n];
    for r in (0..n).rev() {
        let mut s = b[r];
        for k in r + 1..n {
            s -= a[r * n + k] * x[k];
        }
        x[r] = s / a[r * n + r];
    }
    x
}

// ---------------------------------------------------------------- JW フェルミオン (3 モード)

/// 3 モードの消滅演算子 (JW): a₁ = σ⁻II, a₂ = Zσ⁻I, a₃ = ZZσ⁻
fn jw_annihilators() -> Vec<Vec<C64>> {
    let sm = {
        // σ⁻ = (X + iY)/2 = [[0,1],[0,0]] は c|1⟩ = |0⟩ の規約 (占有 = bit 1)
        let mut m = vec![C64::new(0.0, 0.0); 4];
        m[1] = C64::new(1.0, 0.0);
        m
    };
    let z = pauli('Z');
    let i2 = ident(2);
    let a1 = kron(&kron(&sm, 2, &i2, 2), 4, &i2, 2);
    let a2 = kron(&kron(&z, 2, &sm, 2), 4, &i2, 2);
    let a3 = kron(&kron(&z, 2, &z, 2), 4, &sm, 2);
    vec![a1, a2, a3]
}

/// hopping H(θ) = Σ_bonds t(e^{iθ} a†_i a_j + h.c.)
fn ring_h(theta: f64, bonds: &[(usize, usize)], t: f64) -> Vec<C64> {
    let ann = jw_annihilators();
    let n = 8;
    let mut h = vec![C64::new(0.0, 0.0); n * n];
    for &(i, j) in bonds {
        let hop = cmul(&cdag(&ann[i], n), &ann[j], n);
        let phase = C64::expi(theta).scale(t);
        // e^{iθ} a†_i a_j + e^{−iθ} a†_j a_i
        for (kx, hv) in hop.iter().enumerate() {
            h[kx] = h[kx] + phase * *hv;
        }
        let hop2 = cmul(&cdag(&ann[j], n), &ann[i], n);
        let phase2 = C64::expi(-theta).scale(t);
        for (kx, hv) in hop2.iter().enumerate() {
            h[kx] = h[kx] + phase2 * *hv;
        }
    }
    h
}

fn number_op(k: usize) -> Vec<C64> {
    let ann = jw_annihilators();
    cmul(&cdag(&ann[k], 8), &ann[k], 8)
}

/// 積状態接ベクトル A_i = (⊗_{k≠i} I/2) ⊗ (n_i − 1/2) = (n_i − I/2)/4 (8 次元)
fn product_tangent(i: usize) -> Vec<C64> {
    let n = 8;
    let mut a = number_op(i);
    for d in 0..n {
        a[d * n + d] = a[d * n + d] - C64::new(0.5, 0.0);
    }
    for x in a.iter_mut() {
        *x = x.scale(0.25);
    }
    a
}

fn main() {
    uft_sim::self_test();
    println!(
        "=== v32.4 Liouvillian 応答階層 — 恒等式・中心を除く決定・H ↔ −H no-go (PROMPT/13 §3) ===\n"
    );
    let root = if Path::new("core.schema.yml").exists() {
        "."
    } else {
        ".."
    };
    let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p));
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---- [L0] 恒等式の測定照合 (dim 8) ----
    {
        let (h8, n) = build_h(&[
            ("XXI", 0.9),
            ("IYZ", -0.7),
            ("ZIZ", 0.5),
            ("YII", 0.3),
            ("IZI", -0.45),
            ("XYX", 0.62),
        ]);
        let rho0: Vec<C64> = ident(n).iter().map(|c| c.scale(1.0 / n as f64)).collect();
        // 非自明対の選別: A 4 種 × 全 traceless 三重 63 種から |解析値| > 0.1 の対を
        // 解析式で拾い (Pauli 直交の偶然ゼロを避ける)、先頭 8 対を測定と照合する
        let a_set = ["ZII", "IXI", "ZXI", "IIZ"];
        let cs = ['I', 'X', 'Y', 'Z'];
        let mut b_labels = Vec::new();
        for &c1 in &cs {
            for &c2 in &cs {
                for &c3 in &cs {
                    let s: String = [c1, c2, c3].iter().collect();
                    if s != "III" {
                        b_labels.push(s);
                    }
                }
            }
        }
        // R⁽¹⁾ の非自明対 (奇積) と R⁽²⁾ の非自明対 (偶積) は一般に別 — 各 4 対
        let mut p1 = Vec::new();
        let mut p2 = Vec::new();
        for sa in a_set {
            for sb in &b_labels {
                let a = op_str(sa);
                let b = op_str(sb);
                let r1 = r1_exact(&h8, &b, &a, n);
                let r2 = r2_exact(&h8, &b, &a, n);
                if r1.abs() > 0.1 && p1.len() < 4 {
                    p1.push((a.clone(), b.clone(), r1));
                }
                if r2.abs() > 0.1 && p2.len() < 4 {
                    p2.push((a, b, r2));
                }
            }
        }
        let mut max_dev1 = 0.0f64;
        let mut max_dev2 = 0.0f64;
        let mut max_eps = 0.0f64;
        let mut sample = (0.0, 0.0, 0.0, 0.0);
        for (i, (a, b, r1)) in p1.iter().enumerate() {
            let (m1a, _) = measure_r12(&h8, &rho0, a, b, 0.02, n);
            max_dev1 = max_dev1.max((m1a - r1).abs() / r1.abs());
            if i == 0 {
                let (m1b, _) = measure_r12(&h8, &rho0, a, b, 0.10, n);
                max_eps = max_eps.max((m1a - m1b).abs());
                sample.0 = m1a;
                sample.1 = *r1;
            }
        }
        for (i, (a, b, r2)) in p2.iter().enumerate() {
            let (_, m2a) = measure_r12(&h8, &rho0, a, b, 0.02, n);
            max_dev2 = max_dev2.max((m2a - r2).abs() / r2.abs());
            if i == 0 {
                let (_, m2b) = measure_r12(&h8, &rho0, a, b, 0.10, n);
                max_eps = max_eps.max((m2a - m2b).abs());
                sample.2 = m2a;
                sample.3 = *r2;
            }
        }
        let ok = p1.len() >= 3 && p2.len() >= 3 && max_dev1 <= 1e-6 && max_dev2 <= 1e-6 && max_eps <= 1e-7;
        check(
            "[L0] 恒等式 — 測定 (expm + 4 次 stencil) が R⁽¹⁾ = −i Tr(B[H,A])・R⁽²⁾ = Tr([H,B][H,A]) と一致・ε 非依存",
            ok,
            format!(
                "非自明対 R⁽¹⁾ {} 本 / R⁽²⁾ {} 本: max rel 偏差 {:.1e} / {:.1e}・ε 0.02 vs 0.10 差 {:.1e}・例 R⁽¹⁾ {:.6} = {:.6}, R⁽²⁾ {:.6} = {:.6}",
                p1.len(),
                p2.len(),
                max_dev1,
                max_dev2,
                max_eps,
                sample.0,
                sample.1,
                sample.2,
                sample.3
            ),
        );
    }

    // ---- [L1] 一階応答は H を中心を除いて決める (dim 4, 情報完全基底) ----
    {
        let n = 4;
        let (h4, _) = build_h(&[("XX", 0.8), ("ZI", -0.5), ("YZ", 0.35), ("IY", 0.2)]);
        // 情報完全 traceless 基底: II 以外の Pauli 対 15 種 (正規化 /2)
        let labels: Vec<String> = {
            let cs = ['I', 'X', 'Y', 'Z'];
            let mut v = Vec::new();
            for &c1 in &cs {
                for &c2 in &cs {
                    let s: String = [c1, c2].iter().collect();
                    if s != "II" {
                        v.push(s);
                    }
                }
            }
            v
        };
        let basis: Vec<Vec<C64>> = labels
            .iter()
            .map(|s| {
                op_str(s)
                    .iter()
                    .map(|c| c.scale(0.5))
                    .collect::<Vec<C64>>()
            })
            .collect();
        let k = basis.len();
        // R⁽¹⁾ 行列と設計テンソル T_{(βα),γ} = i Tr(B_β [E_γ, A_α])
        let mut gram = vec![0.0; k * k];
        let mut rhs = vec![0.0; k];
        for bi in 0..k {
            for ai in 0..k {
                let m = r1_exact(&h4, &basis[bi], &basis[ai], n);
                let mut t_row = vec![0.0; k];
                for gi in 0..k {
                    t_row[gi] = r1_exact(&basis[gi], &basis[bi], &basis[ai], n);
                }
                for g1 in 0..k {
                    rhs[g1] += t_row[g1] * m;
                    for g2 in 0..k {
                        gram[g1 * k + g2] += t_row[g1] * t_row[g2];
                    }
                }
            }
        }
        let x = gauss_solve(&mut gram, &mut rhs, k);
        // 復元 Ĥ = Σ x_γ E_γ と真の traceless 部の比較
        let mut hhat = vec![C64::new(0.0, 0.0); n * n];
        for (gi, e) in basis.iter().enumerate() {
            add_scaled(&mut hhat, e, C64::new(x[gi], 0.0));
        }
        let diff: f64 = hhat
            .iter()
            .zip(h4.iter())
            .map(|(a, b)| (*a - *b).norm2())
            .sum::<f64>()
            .sqrt();
        let scale = hs_norm(&h4);
        // 中心の不可視性: H + cI は R⁽¹⁾ を厳密に変えない
        let mut hc = h4.clone();
        for i in 0..n {
            hc[i * n + i] = hc[i * n + i] + C64::new(0.37, 0.0);
        }
        let mut center_dev = 0.0f64;
        for bi in 0..k {
            for ai in 0..k {
                center_dev = center_dev.max(
                    (r1_exact(&h4, &basis[bi], &basis[ai], n)
                        - r1_exact(&hc, &basis[bi], &basis[ai], n))
                    .abs(),
                );
            }
        }
        let ok = diff / scale <= 1e-10 && center_dev <= 1e-12;
        check(
            "[L1] 一階の完全性 — 情報完全 15 基底の R⁽¹⁾ から線形逆問題で Ĥ 復元 (中心を除く)・H + cI は不可視",
            ok,
            format!(
                "‖Ĥ − H_traceless‖/‖H‖ = {:.2e}・H + 0.37I の R⁽¹⁾ 差 = {:.2e} (中心はスカラー恒等項)",
                diff / scale,
                center_dev
            ),
        );
    }

    // ---- [L2] no-go の対: H ↔ −H (二階) と磁束 (密度 lane) ----
    {
        let mut bad = Vec::new();
        // (a) R⁽²⁾[H] = R⁽²⁾[−H] 厳密・R⁽¹⁾ は符号分離
        let (h8, n) = build_h(&[("XXI", 0.9), ("IYZ", -0.7), ("ZIZ", 0.5), ("YII", 0.3)]);
        let hneg: Vec<C64> = h8.iter().map(|c| c.scale(-1.0)).collect();
        // A は 4 種・B は全 traceless Pauli 三重 (63) — 偶然ゼロを避けて網羅
        let a_probes = [op_str("ZII"), op_str("IXI"), op_str("ZXI"), op_str("IIZ")];
        let cs = ['I', 'X', 'Y', 'Z'];
        let mut b_probes = Vec::new();
        for &c1 in &cs {
            for &c2 in &cs {
                for &c3 in &cs {
                    let s: String = [c1, c2, c3].iter().collect();
                    if s != "III" {
                        b_probes.push(op_str(&s));
                    }
                }
            }
        }
        let mut max_r2_diff = 0.0f64;
        let mut max_r1_sum = 0.0f64;
        let mut max_r1_abs = 0.0f64;
        for b in &b_probes {
            for a in &a_probes {
                let r2p = r2_exact(&h8, b, a, n);
                let r2m = r2_exact(&hneg, b, a, n);
                max_r2_diff = max_r2_diff.max((r2p - r2m).abs());
                let r1p = r1_exact(&h8, b, a, n);
                let r1m = r1_exact(&hneg, b, a, n);
                max_r1_sum = max_r1_sum.max((r1p + r1m).abs());
                max_r1_abs = max_r1_abs.max(r1p.abs());
            }
        }
        if max_r2_diff > 1e-12 {
            bad.push(format!("R⁽²⁾[H] ≠ R⁽²⁾[−H] ({:.1e})", max_r2_diff));
        }
        if !(max_r1_sum <= 1e-12 && max_r1_abs > 0.1) {
            bad.push("R⁽¹⁾ が符号を分離しない".into());
        }
        // (b) 磁束 ring: 密度二階核は θ 不可視・coherent 一階は分離
        let bonds = [(0usize, 1usize), (1, 2), (2, 0)];
        let h0 = ring_h(0.0, &bonds, 1.0);
        let hf = ring_h(0.4, &bonds, 1.0);
        let mut dens_diff = 0.0f64;
        for j in 0..3 {
            for i in 0..3 {
                if i == j {
                    continue;
                }
                let r0 = r2_exact(&h0, &number_op(j), &product_tangent(i), 8);
                let rf = r2_exact(&hf, &number_op(j), &product_tangent(i), 8);
                dens_diff = dens_diff.max((r0 - rf).abs());
            }
        }
        if dens_diff > 1e-12 {
            bad.push(format!("密度二階核が磁束を見た ({:.1e})", dens_diff));
        }
        // coherent 一階: 電流観測 J₀₁ = i(a†₀a₁ − a†₁a₀)
        let ann = jw_annihilators();
        let hop = cmul(&cdag(&ann[0], 8), &ann[1], 8);
        let hop_d = cdag(&hop, 8);
        let mut j01 = vec![C64::new(0.0, 0.0); 64];
        for kx in 0..64 {
            let d = hop[kx] - hop_d[kx];
            j01[kx] = C64::new(-d.im, d.re); // i·(hop − hop†)
        }
        let c0 = r1_exact(&h0, &j01, &product_tangent(0), 8);
        let cf = r1_exact(&hf, &j01, &product_tangent(0), 8);
        if (c0 - cf).abs() < 0.05 {
            bad.push(format!("coherent 一階が磁束を分離しない (|Δ| = {:.3})", (c0 - cf).abs()));
        }
        check(
            "[L2] no-go の対 — R⁽²⁾ は H ↔ −H 不可視 (厳密)・磁束は密度二階に不可視で coherent 一階が分離",
            bad.is_empty(),
            format!(
                "max|R⁽²⁾[H]−R⁽²⁾[−H]| = {:.1e}・R⁽¹⁾ 符号分離 (max|和| = {:.1e})・密度核の θ 差 = {:.1e}・電流 R⁽¹⁾: θ=0 → {:.4} / θ=0.4 → {:.4}",
                max_r2_diff, max_r1_sum, dens_diff, c0, cf
            ),
        );
    }

    // ---- [L3] 保存則と PSD ----
    {
        let mut bad = Vec::new();
        let bonds = [(0usize, 1usize), (1, 2)];
        let h = ring_h(0.0, &bonds, 1.0); // 数保存 (開鎖 hopping)
        let ntot = {
            let mut m = number_op(0);
            for k in 1..3 {
                let nk = number_op(k);
                for (x, y) in m.iter_mut().zip(nk.iter()) {
                    *x = *x + *y;
                }
            }
            m
        };
        let probes: Vec<Vec<C64>> = (0..3).map(product_tangent).collect();
        let mut max_r = 0.0f64;
        for a in &probes {
            max_r = max_r.max(r1_exact(&h, &ntot, a, 8).abs());
            max_r = max_r.max(r2_exact(&h, &ntot, a, 8).abs());
        }
        if max_r > 1e-12 {
            bad.push(format!("N̂ の応答が 0 でない ({:.1e})", max_r));
        }
        // PSD Gram: K_αα' = Tr([H,A_α]†[H,A_α']) — dim-4 の情報完全基底で
        let n4 = 4;
        let (h4, _) = build_h(&[("XX", 0.8), ("ZI", -0.5), ("YZ", 0.35), ("IY", 0.2)]);
        let cs = ['I', 'X', 'Y', 'Z'];
        let mut probes4 = Vec::new();
        for &c1 in &cs {
            for &c2 in &cs {
                let s: String = [c1, c2].iter().collect();
                if s != "II" {
                    probes4.push(op_str(&s));
                }
            }
        }
        let k = probes4.len();
        let mut gram = vec![0.0; k * k];
        for (i, ai) in probes4.iter().enumerate() {
            let hai = commutator(&h4, ai, n4);
            for (j, aj) in probes4.iter().enumerate() {
                let haj = commutator(&h4, aj, n4);
                let mut s = C64::new(0.0, 0.0);
                for (x, y) in cdag(&hai, n4).iter().zip(haj.iter()) {
                    // tr(hai† haj) = Σ (hai†)_{ik} (haj)_{ki}: 直接は cmul のトレースだが
                    // HS 内積 tr(A†B) = Σ conj(A_k) B_k で足す
                    let _ = (x, y);
                }
                for kx in 0..n4 * n4 {
                    let a = hai[kx];
                    let b = haj[kx];
                    s = s + C64::new(a.re * b.re + a.im * b.im, a.re * b.im - a.im * b.re);
                }
                gram[i * k + j] = s.re;
            }
        }
        let (evals, _) = jacobi_eigh(&gram, k);
        let emin = evals.iter().cloned().fold(f64::INFINITY, f64::min);
        let emax = evals.iter().cloned().fold(0.0f64, f64::max);
        if emin < -1e-12 * emax.max(1.0) {
            bad.push(format!("応答 Gram 核が PSD でない (λ_min = {:.1e})", emin));
        }
        check(
            "[L3] 保存則と PSD — 数保存 H で N̂ の応答は全次数 0 (和則)・Gram 核 K = Tr([H,A]†[H,A']) は PSD",
            bad.is_empty(),
            format!("max|R(N̂)| = {:.1e}・K の λ_min = {:.2e} (λ_max = {:.2})", max_r, emin, emax),
        );
    }

    // ---- [L4] v31.2 曲率則 = 本階層の対角特殊化 ----
    {
        let bonds = [(0usize, 1usize), (1, 2)];
        let t = 0.8;
        let h = ring_h(0.0, &bonds, t);
        // one-body 行列 h_ij (3×3): 隣接 t
        let mut max_dev = 0.0f64;
        let mut vals = Vec::new();
        for j in 0..3 {
            for i in 0..3 {
                if i == j {
                    continue;
                }
                let r2 = r2_exact(&h, &number_op(j), &product_tangent(i), 8);
                let h_ij = if bonds.contains(&(i, j)) || bonds.contains(&(j, i)) {
                    t
                } else {
                    0.0
                };
                max_dev = max_dev.max((r2 - h_ij * h_ij).abs());
                if h_ij > 0.0 {
                    vals.push(r2);
                }
            }
        }
        let ok = max_dev <= 1e-12;
        check(
            "[L4] 統一 — many-body R⁽²⁾_ji (密度 × 積状態接) = one-body |h_ij|² (v31.2 曲率則の特殊化)",
            ok,
            format!(
                "max|R⁽²⁾_ji − |h_ij|²| = {:.1e} (辺の読み {:?} vs t² = {:.2})",
                max_dev,
                vals.iter().map(|v| (v * 100.0).round() / 100.0).collect::<Vec<f64>>(),
                t * t
            ),
        );
    }

    // ---- [L5] 文書アンカー ----
    {
        let mut bad = Vec::new();
        let doc = rd("docs/uft-v32.4.md").unwrap_or_default();
        for needle in [
            "R⁽¹⁾",
            "R⁽²⁾",
            "H ↔ −H",
            "中心",
            "磁束",
            "保存則",
            "曲率則の特殊化",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v32.4.md: 「{}」が無い", needle));
            }
        }
        check(
            "[L5] 文書アンカー — uft-v32.4.md の恒等式・no-go・統一",
            bad.is_empty(),
            if bad.is_empty() {
                "応答階層が恒等式・no-go・統一の三点で凍結された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "一階は中心を除いて H を読み・二階は符号を読めない — 「破れ」ではなく観測契約ごとの読める同値類が確定した"
        } else {
            "**応答階層の破れ** — 恒等式と文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
