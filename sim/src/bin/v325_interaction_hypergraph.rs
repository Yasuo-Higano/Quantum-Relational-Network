//! v32.5 Interaction hypergraph — H_S 直交分解・相関ホッピング・条件付き probe (PROMPT/13 §3)
//!
//! 因子分解が確定した後 (v32.3)、Hamiltonian を**局所 Hilbert–Schmidt 条件期待値**で
//! H = Σ_{S⊆V} H_S と直交分解し、w_S = ‖H_S‖²_F を S-体相互作用の重みと定義する。
//! これは局所演算子基底の選び方に依存せず、block-local unitary で不変 —
//! |S| = 1: on-site / 2: graph edge / 3: correlated-hopping 級 hyperedge / ≥4: 高体。
//!
//! 相関ホッピング V·n₃(c†₁c₂ + h.c.) は中心化 n₃ = 1/2 + (2n₃−1)/2 により
//! **二体成分 (V/2)h₁₂ と真の三体成分 −(V/2)Z₃h₁₂ に正準分離**される。v31.5 の
//! 「密度対角 V は厳密転移・相関 hopping は未走査」の未走査側を、v32.4 の応答階層で
//! **「破れ」でなく観測契約ごとの読みとして採点**する:
//!
//!   [G0] 正準分解の二重定義一致: Möbius 条件期待値 H_S = Σ_{T⊆S} (−1)^{|S\T|} E_T(H)
//!        = Pauli 支持射影 (≤ 1e-12)・直交性 ⟨H_S, H_S'⟩ = 0・完全性 Σ H_S = H・
//!        w_S ≥ 0
//!   [G1] **局所 unitary 不変性**: w_S(UHU†) = w_S(H) (U = u₁⊗u₂⊗u₃, ≤ 1e-12)。
//!        負制御: 非局所 unitary (DFT₈) は w_S を変える (> 0.1)
//!   [G2] **中心化分離の実例**: V·n₃h₁₂ → w_{12} = w_{123} = V² (等重み) —
//!        「三体項」の半分は二体に住む (半充填平均の dressing)
//!   [G3] **条件付き密度 probe (Möbius/Boolean-Fourier lane)**: 補助ノード 3 を
//!        P(n₃=v) に条件付けた曲率 K(v) = |t + vV|² (厳密)・K(1)−K(0) は
//!        {1,2,3} hyperedge 検出器 (V = 0 で厳密 0 = 負制御)・非条件付き曲率 =
//!        (K(0)+K(1))/2 (混合恒等式)
//!   [G4] **coherent parity-even probe が符号を回復**: t ↔ −t は密度曲率に不可視
//!        (厳密 0)・coherent 一階 (電流) が分離 — 密度 lane 単独の正答は
//!        符号同値類 (強制回答しない)
//!   [G5] **hypergraph の組み立てと遷移率和則**: 正準 w から InteractionHypergraph
//!        (v32.2 先凍結型) を構成 — 支持 {1}/{1,2}/{2,3}/{1,2,3}・**非条件付き
//!        密度曲率 = (w_{12} + w_{123})/4 (対を含む hyperedge 重みの和を読む)** —
//!        「条件付き遷移率の Gram 核」の等式化
//!   [G6] 文書アンカー — uft-v32.5.md
//!
//! 実行: cargo run --release --bin v325_interaction_hypergraph

use std::fs;
use std::path::Path;
use uft_sim::operational_net::{cdag, cmul, commutator, hs_inner, InteractionHypergraph};
use uft_sim::C64;

// ---------------------------------------------------------------- 素子 (v322/v324 と同一)

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

fn op3(s: &str) -> Vec<C64> {
    let cs: Vec<char> = s.chars().collect();
    let a = kron(&pauli(cs[0]), 2, &pauli(cs[1]), 2);
    kron(&a, 4, &pauli(cs[2]), 2)
}

fn ident(n: usize) -> Vec<C64> {
    let mut m = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

fn add_scaled(a: &mut [C64], b: &[C64], s: f64) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x = *x + y.scale(s);
    }
}

fn dft8() -> Vec<C64> {
    let n = 8;
    let inv = 1.0 / (n as f64).sqrt();
    let mut f = vec![C64::new(0.0, 0.0); n * n];
    for j in 0..n {
        for k in 0..n {
            f[j * n + k] =
                C64::expi(2.0 * std::f64::consts::PI * (j * k) as f64 / n as f64).scale(inv);
        }
    }
    f
}

fn conj_by(v: &[C64], a: &[C64], n: usize) -> Vec<C64> {
    cmul(&cmul(v, a, n), &cdag(v, n), n)
}

/// JW: 3 モードの消滅演算子
fn jw_annihilators() -> Vec<Vec<C64>> {
    let sm = {
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

fn number_op(k: usize) -> Vec<C64> {
    let ann = jw_annihilators();
    cmul(&cdag(&ann[k], 8), &ann[k], 8)
}

/// hopping h_{ij} = c†_i c_j + c†_j c_i
fn hop(i: usize, j: usize) -> Vec<C64> {
    let ann = jw_annihilators();
    let a = cmul(&cdag(&ann[i], 8), &ann[j], 8);
    let b = cmul(&cdag(&ann[j], 8), &ann[i], 8);
    a.iter().zip(b.iter()).map(|(x, y)| *x + *y).collect()
}

/// R⁽²⁾ = Tr([H,B][H,A]) (v32.4 と同一)
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

/// R⁽¹⁾ = −i Tr(B[H,A]) (v32.4 と同一)
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

// ---------------------------------------------------------------- 正準分解 (3 qubit / 3 モード)

/// Pauli 展開係数 c_P = ⟨P̂, H⟩ (P̂ = P/√8) — 64 文字列
fn pauli_coeffs(h: &[C64]) -> Vec<(String, f64, f64)> {
    let cs = ['I', 'X', 'Y', 'Z'];
    let inv = 1.0 / (8.0f64).sqrt();
    let mut out = Vec::new();
    for &c1 in &cs {
        for &c2 in &cs {
            for &c3 in &cs {
                let s: String = [c1, c2, c3].iter().collect();
                let p: Vec<C64> = op3(&s).iter().map(|c| c.scale(inv)).collect();
                let x = hs_inner(&p, h);
                out.push((s, x.re, x.im));
            }
        }
    }
    out
}

fn supp_of(s: &str) -> Vec<u32> {
    s.chars()
        .enumerate()
        .filter(|(_, c)| *c != 'I')
        .map(|(i, _)| i as u32)
        .collect()
}

/// Pauli 支持射影による H_S (S を昇順ラベルで指定)
fn h_s_pauli(h: &[C64], s: &[u32]) -> Vec<C64> {
    let inv = 1.0 / (8.0f64).sqrt();
    let mut out = vec![C64::new(0.0, 0.0); 64];
    for (label, re, im) in pauli_coeffs(h) {
        if supp_of(&label) == s {
            let p: Vec<C64> = op3(&label).iter().map(|c| c.scale(inv)).collect();
            for (o, pv) in out.iter_mut().zip(p.iter()) {
                *o = *o + C64::new(re, im) * *pv;
            }
        }
    }
    out
}

/// 条件期待値 E_T(H) = (Tr_{T^c} H) ⊗ I_{T^c} / 2^{|T^c|} (T = qubit 集合)
fn cond_exp(h: &[C64], t_set: &[u32]) -> Vec<C64> {
    let n = 8;
    let in_t = |q: u32| t_set.contains(&q);
    // インデックス bit: qubit 0 = 最上位 (op3 の kron 順)
    let bit = |idx: usize, q: u32| (idx >> (2 - q)) & 1;
    let mut out = vec![C64::new(0.0, 0.0); n * n];
    let comp_size = 3 - t_set.len() as u32;
    let scale = 1.0 / (1u32 << comp_size) as f64;
    for r in 0..n {
        for c in 0..n {
            // T 上の成分が (r,c)・T^c 上は対角で足す
            let mut acc = C64::new(0.0, 0.0);
            let mut compatible = true;
            for q in 0..3u32 {
                if !in_t(q) && bit(r, q) != bit(c, q) {
                    compatible = false;
                    break;
                }
            }
            if !compatible {
                continue;
            }
            // T^c の対角自由度を走る: r, c の T^c bit を共通に置換
            for d in 0..n {
                let mut r2 = 0usize;
                let mut c2 = 0usize;
                for q in 0..3u32 {
                    let (br, bc) = if in_t(q) {
                        (bit(r, q), bit(c, q))
                    } else {
                        (bit(d, q), bit(d, q))
                    };
                    r2 |= br << (2 - q);
                    c2 |= bc << (2 - q);
                }
                acc = acc + h[r2 * n + c2];
            }
            // d ループは T の bit も走ってしまう (2^3) — T^c のみ数えるため補正
            let over = (1u32 << t_set.len()) as f64;
            out[r * n + c] = acc.scale(scale / over);
        }
    }
    out
}

/// Möbius: H_S = Σ_{T ⊆ S} (−1)^{|S|−|T|} E_T(H)
fn h_s_moebius(h: &[C64], s: &[u32]) -> Vec<C64> {
    let k = s.len();
    let mut out = vec![C64::new(0.0, 0.0); 64];
    for mask in 0..(1u32 << k) {
        let t: Vec<u32> = (0..k).filter(|&i| mask & (1 << i) != 0).map(|i| s[i]).collect();
        let sign = if (k - t.len()) % 2 == 0 { 1.0 } else { -1.0 };
        let e = cond_exp(h, &t);
        add_scaled(&mut out, &e, sign);
    }
    out
}

/// 全 S の正準重み w_S = ‖H_S‖²̂ (‖·‖̂ = HS/√8 正規化ノルムの二乗 = Σ|c_P|², supp = S)
fn canonical_weights(h: &[C64]) -> std::collections::BTreeMap<Vec<u32>, f64> {
    let mut w: std::collections::BTreeMap<Vec<u32>, f64> = Default::default();
    for (label, re, im) in pauli_coeffs(h) {
        let s = supp_of(&label);
        if s.is_empty() {
            continue;
        }
        *w.entry(s).or_insert(0.0) += re * re + im * im;
    }
    w.retain(|_, v| *v > 1e-20);
    w
}

fn main() {
    uft_sim::self_test();
    println!(
        "=== v32.5 Interaction hypergraph — H_S 直交分解・相関ホッピング・条件付き probe (PROMPT/13 §3) ===\n"
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
    let n = 8usize;

    // 相関ホッピング模型: H = t(h₁₂ + h₂₃) + V·n₃h₁₂ + μ·n₁
    let (t, vv, mu) = (0.8, 0.6, 0.35);
    let build_model = |t: f64, vv: f64, mu: f64| -> Vec<C64> {
        let mut h = vec![C64::new(0.0, 0.0); n * n];
        add_scaled(&mut h, &hop(0, 1), t);
        add_scaled(&mut h, &hop(1, 2), t);
        let n3h12 = cmul(&number_op(2), &hop(0, 1), n); // n₃ は h₁₂ と可換 (積はエルミート)
        add_scaled(&mut h, &n3h12, vv);
        add_scaled(&mut h, &number_op(0), mu);
        h
    };
    let h_model = build_model(t, vv, mu);

    // ---- [G0] 正準分解の二重定義一致 ----
    {
        let mut bad = Vec::new();
        // 全 2^3 − 1 の S 候補
        let all_s: Vec<Vec<u32>> = (1u32..8)
            .map(|m| (0..3u32).filter(|&q| m & (1 << q) != 0).collect())
            .collect();
        let mut max_dev = 0.0f64;
        let mut sum = vec![C64::new(0.0, 0.0); 64];
        let e_empty = cond_exp(&h_model, &[]);
        add_scaled(&mut sum, &e_empty, 1.0); // S = ∅ (トレース部)
        let mut parts: Vec<Vec<C64>> = vec![e_empty];
        for s in &all_s {
            let hp = h_s_pauli(&h_model, s);
            let hm = h_s_moebius(&h_model, s);
            let dev: f64 = hp
                .iter()
                .zip(hm.iter())
                .map(|(a, b)| (*a - *b).norm2())
                .sum::<f64>()
                .sqrt();
            max_dev = max_dev.max(dev);
            add_scaled(&mut sum, &hp, 1.0);
            parts.push(hp);
        }
        if max_dev > 1e-12 {
            bad.push(format!("Möbius ≠ Pauli 支持射影 ({:.1e})", max_dev));
        }
        // 完全性
        let comp: f64 = sum
            .iter()
            .zip(h_model.iter())
            .map(|(a, b)| (*a - *b).norm2())
            .sum::<f64>()
            .sqrt();
        if comp > 1e-12 {
            bad.push(format!("Σ H_S ≠ H ({:.1e})", comp));
        }
        // 直交性
        let mut max_cross = 0.0f64;
        for i in 0..parts.len() {
            for j in i + 1..parts.len() {
                max_cross = max_cross.max(hs_inner(&parts[i], &parts[j]).abs());
            }
        }
        if max_cross > 1e-12 {
            bad.push(format!("⟨H_S, H_S'⟩ ≠ 0 ({:.1e})", max_cross));
        }
        check(
            "[G0] 正準分解 — Möbius 条件期待値 = Pauli 支持射影・直交・完全 (w_S ≥ 0 は構成上)",
            bad.is_empty(),
            format!(
                "max|H_S^Möbius − H_S^Pauli| = {:.1e}・Σ H_S − H = {:.1e}・max 直交残差 = {:.1e}",
                max_dev, comp, max_cross
            ),
        );
    }

    // ---- [G1] 局所 unitary 不変性 + 非局所の負制御 ----
    {
        let rot = |theta: f64, nx: f64, ny: f64, nz: f64| -> Vec<C64> {
            let (c, s) = (theta.cos(), theta.sin());
            vec![
                C64::new(c, s * nz),
                C64::new(s * ny, s * nx),
                C64::new(-s * ny, s * nx),
                C64::new(c, -s * nz),
            ]
        };
        let u12 = kron(&rot(0.4, 1.0, 0.0, 0.0), 2, &rot(0.9, 0.0, 0.6, 0.8), 2);
        let u = kron(&u12, 4, &rot(1.3, 0.0, 1.0, 0.0), 2);
        let w0 = canonical_weights(&h_model);
        let w_loc = canonical_weights(&conj_by(&u, &h_model, n));
        let mut max_dev = 0.0f64;
        for (s, w) in &w0 {
            max_dev = max_dev.max((w - w_loc.get(s).copied().unwrap_or(0.0)).abs());
        }
        let w_nl = canonical_weights(&conj_by(&dft8(), &h_model, n));
        let mut max_change = 0.0f64;
        for (s, w) in &w0 {
            max_change = max_change.max((w - w_nl.get(s).copied().unwrap_or(0.0)).abs());
        }
        let ok = max_dev <= 1e-12 && max_change > 0.1;
        check(
            "[G1] w_S は block-local unitary 不変・非局所 unitary (DFT₈) は変える (負制御)",
            ok,
            format!(
                "局所 U での max|Δw_S| = {:.1e}・DFT₈ での max|Δw_S| = {:.3} (> 0.1)",
                max_dev, max_change
            ),
        );
    }

    // ---- [G2] 中心化分離の実例 ----
    {
        // 純粋な相関項 V·n₃h₁₂ のみ
        let mut h_corr = vec![C64::new(0.0, 0.0); n * n];
        add_scaled(&mut h_corr, &cmul(&number_op(2), &hop(0, 1), n), vv);
        let w = canonical_weights(&h_corr);
        let w12 = w.get(&vec![0u32, 1]).copied().unwrap_or(0.0);
        let w123 = w.get(&vec![0u32, 1, 2]).copied().unwrap_or(0.0);
        let expect = (vv / 2.0) * (vv / 2.0) * 4.0; // (V/2)²·‖h₁₂‖_F² = V² (‖h₁₂‖_F² = 4)
        let ok = (w12 - expect).abs() <= 1e-12
            && (w123 - expect).abs() <= 1e-12
            && w.len() == 2;
        check(
            "[G2] 中心化分離 — V·n₃h₁₂ = 二体 (V/2)h₁₂ ⊕ 三体 −(V/2)Z₃h₁₂ (等重み (V/2)²‖h₁₂‖²_F = V²)",
            ok,
            format!(
                "w_{{12}} = {:.10} / w_{{123}} = {:.10} (= V² = {:.10})・他の S なし = {}",
                w12,
                w123,
                expect,
                w.len() == 2
            ),
        );
    }

    // ---- [G3] 条件付き密度 probe (Möbius lane) ----
    {
        let mut bad = Vec::new();
        // 条件付き接ベクトル: A_v = (n₁ − 1/2) ⊗ (I/2) ⊗ P_v(qubit3)
        let tangent_cond = |v: usize| -> Vec<C64> {
            let mut a = vec![C64::new(0.0, 0.0); n * n];
            for idx in 0..n {
                let b1 = (idx >> 2) & 1;
                let b3 = idx & 1;
                if b3 != v {
                    continue;
                }
                let val = (b1 as f64 - 0.5) * 0.5; // (n₁−1/2)·(I/2 on site2)
                a[idx * n + idx] = C64::new(val, 0.0);
            }
            a
        };
        let tangent_uncond = {
            let mut a = vec![C64::new(0.0, 0.0); n * n];
            for idx in 0..n {
                let b1 = (idx >> 2) & 1;
                let val = (b1 as f64 - 0.5) * 0.25;
                a[idx * n + idx] = C64::new(val, 0.0);
            }
            a
        };
        let n2op = number_op(1);
        let k0 = r2_exact(&h_model, &n2op, &tangent_cond(0), n);
        let k1 = r2_exact(&h_model, &n2op, &tangent_cond(1), n);
        let ku = r2_exact(&h_model, &n2op, &tangent_uncond, n);
        // 期待: K(v) = |t + vV|² (占有 v の条件付き遷移率)
        if (k0 - t * t).abs() > 1e-12 {
            bad.push(format!("K(0) = {:.10} ≠ t² = {:.10}", k0, t * t));
        }
        if (k1 - (t + vv) * (t + vv)).abs() > 1e-12 {
            bad.push(format!("K(1) = {:.10} ≠ (t+V)² = {:.10}", k1, (t + vv) * (t + vv)));
        }
        if (ku - 0.5 * (k0 + k1)).abs() > 1e-12 {
            bad.push("混合恒等式 K_uncond = (K(0)+K(1))/2 が破れた".into());
        }
        // 負制御: V = 0 で条件差は厳密 0
        let h_free = build_model(t, 0.0, mu);
        let d_free = r2_exact(&h_free, &n2op, &tangent_cond(1), n)
            - r2_exact(&h_free, &n2op, &tangent_cond(0), n);
        if d_free.abs() > 1e-13 {
            bad.push(format!("V=0 の条件差が 0 でない ({:.1e})", d_free));
        }
        check(
            "[G3] 条件付き密度 probe — K(v) = |t + vV|² 厳密・K(1)−K(0) は hyperedge 検出器 (V=0 負制御)・混合恒等式",
            bad.is_empty(),
            format!(
                "K(0) = {:.6} = t²・K(1) = {:.6} = (t+V)²・K_uncond = {:.6} = 平均・V=0 条件差 = {:.1e}",
                k0, k1, ku, d_free
            ),
        );
    }

    // ---- [G4] coherent parity-even probe が符号を回復 ----
    {
        // 純 hopping の符号 t ↔ −t を比較する (V ≠ 0 では K(1) = |±t+V|² が
        // 異なるため、密度不可視性の対象は V = 0 の符号自由度)
        let hp0 = build_model(t, 0.0, mu);
        let hm0 = build_model(-t, 0.0, mu);
        let mut dens_diff = 0.0f64;
        for j in 0..3 {
            for i in 0..3 {
                if i == j {
                    continue;
                }
                let mut a = vec![C64::new(0.0, 0.0); n * n];
                for idx in 0..n {
                    let bi = (idx >> (2 - i)) & 1;
                    let val = (bi as f64 - 0.5) * 0.25;
                    a[idx * n + idx] = C64::new(val, 0.0);
                }
                let d = r2_exact(&hp0, &number_op(j), &a, n)
                    - r2_exact(&hm0, &number_op(j), &a, n);
                dens_diff = dens_diff.max(d.abs());
            }
        }
        // coherent 一階: B = i(c†₀c₁ − c†₁c₀), A = hop 実部演算子 (parity-even probe)
        let ann = jw_annihilators();
        let hopc = cmul(&cdag(&ann[0], n), &ann[1], n);
        let hopd = cdag(&hopc, n);
        let jcur: Vec<C64> = hopc
            .iter()
            .zip(hopd.iter())
            .map(|(x, y)| {
                let d = *x - *y;
                C64::new(-d.im, d.re)
            })
            .collect();
        // 密度不均衡の準備 (parity-even) → 電流の一階応答は t に奇
        let a_dens: Vec<C64> = {
            let mut a = vec![C64::new(0.0, 0.0); n * n];
            for idx in 0..n {
                let b0 = (idx >> 2) & 1;
                a[idx * n + idx] = C64::new((b0 as f64 - 0.5) * 0.25, 0.0);
            }
            a
        };
        let c_plus = r1_exact(&hp0, &jcur, &a_dens, n);
        let c_minus = r1_exact(&hm0, &jcur, &a_dens, n);
        let ok = dens_diff <= 1e-12 && (c_plus - c_minus).abs() > 0.1 && (c_plus + c_minus).abs() <= 1e-12;
        check(
            "[G4] 符号の観測契約 — t ↔ −t は密度曲率に厳密不可視・coherent 一階 (電流) が分離 (密度単独は符号同値類)",
            ok,
            format!(
                "密度曲率の t↔−t 差 = {:.1e}・coherent R⁽¹⁾: +t → {:.4} / −t → {:.4} (奇 = 符号読み出し)",
                dens_diff, c_plus, c_minus
            ),
        );
    }

    // ---- [G5] hypergraph の組み立てと遷移率和則 ----
    {
        let mut bad = Vec::new();
        let w = canonical_weights(&h_model);
        let hg = InteractionHypergraph {
            n_nodes: 3,
            weights: w.clone(),
        };
        let thr = 1e-12;
        let s1 = hg.support_of_order(1, thr);
        let s2 = hg.support_of_order(2, thr);
        let s3 = hg.support_of_order(3, thr);
        // 期待: on-site {0} (μ n₁ の Z₁ 成分) と {2} は? n₃h₁₂ の I 部? no —
        // n₃h₁₂ 中心化は S = {0,1} と {0,1,2} のみ。on-site は {0} (μ)
        if s1 != vec![vec![0u32]] {
            bad.push(format!("order-1 支持 {:?} ≠ [[0]]", s1));
        }
        if s2 != vec![vec![0u32, 1], vec![1, 2]] {
            bad.push(format!("order-2 支持 {:?} ≠ [[0,1],[1,2]]", s2));
        }
        if s3 != vec![vec![0u32, 1, 2]] {
            bad.push(format!("order-3 支持 {:?} ≠ [[0,1,2]]", s3));
        }
        // 数値 (‖h_bond‖²_F = 4, ‖Z₁‖²_F = 8): w_{01} = 4(t+V/2)², w_{12} = 4t²,
        // w_{012} = 4(V/2)² = V², w_{0} = 2μ² (μn₁ = μ(I−Z₁)/2 の Z₁ 成分)
        let w01 = w.get(&vec![0u32, 1]).copied().unwrap_or(0.0);
        let w12b = w.get(&vec![1u32, 2]).copied().unwrap_or(0.0);
        let w012 = w.get(&vec![0u32, 1, 2]).copied().unwrap_or(0.0);
        let w0 = w.get(&vec![0u32]).copied().unwrap_or(0.0);
        for (name, got, want) in [
            ("w_{01}", w01, 4.0 * (t + vv / 2.0) * (t + vv / 2.0)),
            ("w_{12}", w12b, 4.0 * t * t),
            ("w_{012}", w012, vv * vv),
            ("w_{0}", w0, 2.0 * mu * mu),
        ] {
            if (got - want).abs() > 1e-12 {
                bad.push(format!("{} = {:.10} ≠ {:.10}", name, got, want));
            }
        }
        // 遷移率和則: 非条件付き密度曲率 K_uncond(2←1) = 2·(w_{01} + w_{012})
        let mut a = vec![C64::new(0.0, 0.0); n * n];
        for idx in 0..n {
            let b1 = (idx >> 2) & 1;
            a[idx * n + idx] = C64::new((b1 as f64 - 0.5) * 0.25, 0.0);
        }
        let ku = r2_exact(&h_model, &number_op(1), &a, n);
        if (ku - (w01 + w012) / 4.0).abs() > 1e-12 {
            bad.push(format!(
                "和則 K_uncond = (w_{{01}} + w_{{012}})/4 が破れた ({:.10} vs {:.10})",
                ku,
                (w01 + w012) / 4.0
            ));
        }
        check(
            "[G5] hypergraph 組み立て + 遷移率和則 — 支持 {0}/{0,1}/{1,2}/{0,1,2}・K_uncond = (w₂ + w₃)/4",
            bad.is_empty(),
            format!(
                "w = {{0}}: {:.4} / {{0,1}}: {:.4} / {{1,2}}: {:.4} / {{0,1,2}}: {:.4}・K_uncond(2←1) = {:.6} = (w_{{01}}+w_{{012}})/4 — 密度曲率は対を含む hyperedge 重みの和を読む",
                w0, w01, w12b, w012, ku
            ),
        );
    }

    // ---- [G6] 文書アンカー ----
    {
        let mut bad = Vec::new();
        let doc = rd("docs/uft-v32.5.md").unwrap_or_default();
        for needle in [
            "Hilbert–Schmidt 条件期待値",
            "局所 unitary 不変",
            "中心化",
            "hyperedge 検出器",
            "遷移率和則",
            "符号同値類",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v32.5.md: 「{}」が無い", needle));
            }
        }
        check(
            "[G6] 文書アンカー — uft-v32.5.md の正準分解・lane・和則",
            bad.is_empty(),
            if bad.is_empty() {
                "hypergraph の定義・観測契約・和則が文書に凍結された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "相関ホッピングは「破れ」ではない — 密度曲率は hyperedge 重みの和 (条件付き遷移率の Gram 核) を読み、条件付き probe が次数を分離し、coherent probe が符号を読む"
        } else {
            "**hypergraph 分解の破れ** — 正準分解と文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
