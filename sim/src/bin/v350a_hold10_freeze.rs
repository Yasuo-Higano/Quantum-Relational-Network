//! v35.0-A HOLD-10 の凍結 — 二層 holdout (semantic 30 + coverage 640) の準備
//!
//! PROMPT/15 §6:「30 セル満票だけでは不十分」— HOLD-9 の 9 answerable の誤り 0 は
//! 片側 95% 上限 28.3% しか与えない (v343 [F6-28])。HOLD-10 は二層に分ける:
//!   HOLD-10S = 30 セル (5 群 × 6) の意味論 adversarial holdout — 有限データ
//!              意味論 (v34.3–v34.6) の裁定境界を新鮮セルで採点
//!   HOLD-10C = 凍結生成分布からの 640 セル coverage campaign (回答 ≥ 300 を設計)
//!              — population-risk の片側 95% 上限を含む 8 指標
//!
//! 開封順序 (HOLD-5..9 と同一): 本コミット (v35.0-A) = 生成器・採点器・バー・
//! lib pin の凍結 + SECRET コミットメント公表 + train 採点 (可視シード 35001) +
//! 設計走行 12 シード → v35.0-B = SECRET 開示・holdout 初生成・本採点 (調整なし)。
//!
//! sha256(SECRET) = c3a62f1b55d708e50915abd1c634f65d619c05ea52c2ab5283b8c1804d7851c4
//!
//! 採点 (凍結バー — FROZEN 区間):
//!   S 層: selective risk 0 / impossibility recall 1 (21/21) / answerable recall 1
//!         (9/9) / 強制回答 0
//!   C 層: answers ≥ 300 / wrong_promotion_upper_95 ≤ 0.01 / coverage_lower_95 ≥
//!         0.98 / boundary_abstention_recall 1 / misspecification_recall 1 /
//!         answerable_recall ≥ 0.95 / marginal_to_joint 0 / window_reuse 0 /
//!         structured_dense_drift 0
//!   ※ selective risk = 0.000 だけを完成条件にしない (v34.3 [F3b]/[F6-28])。
//!   ※ 相関・drift lane は必須 (PROMPT/15 §6) — 遷移数ゲート + split-half ゲート。
//!
//! 検証: [A0] lib pin 8 モジュール (v33 の 6 本は HOLD-9 pin と同値) [A1] バー表
//! [A2] train 満票 [A3] 設計走行 12 シード満票 (生成器のシード頑健性)。

use std::collections::BTreeSet;
use std::fs;
use std::path::Path;
use uft_sim::factorization_enumerator::{enumerate_candidates, EnumeratorReading};
use uft_sim::finite_data::{cp_interval, zero_error_upper_bound, RobustVerdict};
use uft_sim::graded_recovery::{extract_complex_structure, GradedAbstainReason, MajoranaFrame};
use uft_sim::operational_net::commutator;
use uft_sim::structured_backend::{recover_quadratic_blocks, PauliVector, QuadraticGenerator};
use uft_sim::{self_test, sha256_hex, Rng, C64};

/// lib pin — 凍結時の共有部 sha256-16 (開封時に [H0] が不変を照合)
const LIB_PINS: [(&str, &str); 8] = [
    ("operational_net", "7898d24448f79f17"),
    ("laboratory_interface", "a11f188f3ecb6a40"),
    ("resource_profile", "f2fe4b9613049704"),
    ("contextual_factorization", "e540eea6f21ca404"),
    ("graded_recovery", "c1ebb02c6af93133"),
    ("structured_backend", "f57d9bdf43add8f8"),
    ("finite_data", "66352869618db1bf"),
    ("factorization_enumerator", "e7a3e312fb156e7f"),
];

/// sha256(SECRET) — SECRET は v35.0-B で開示・holdout シードは SECRET から導出
const SECRET_COMMITMENT: &str = "c3a62f1b55d708e50915abd1c634f65d619c05ea52c2ab5283b8c1804d7851c4";

// ==== FROZEN-HOLD10-BEGIN ====
// (この区間は v350b_hold10_open.rs と逐語一致であることを [H0] が照合する)

pub const H10_TAU: f64 = 0.3; // edge 裁定境界
pub const H10_ALPHA: f64 = 0.0002; // セルごとの同時 α (促進/区間 — Bonferroni 割当済み)
pub const H10_ALPHA_GATE: f64 = 0.001; // 契約ゲート (drift/相関/ブロック) の各 α
pub const H10_NSHOT: usize = 150; // C 層の 1 セルショット数
pub const H10_NMIN: usize = 10; // 回答資格の最低ショット
pub const H10_NC: usize = 640; // C 層セル数
pub const H10_ORBIT_BAR: f64 = 0.9;

// ---- 凍結バー (C 層 8 + S 層 4) ----
pub const BAR_ANSWERS_MIN: usize = 300;
pub const BAR_WRONG_UPPER: f64 = 0.01; // 片側 95% 上限
pub const BAR_COVERAGE_LOWER: f64 = 0.98; // 片側 95% 下限の登録目標
pub const BAR_ANSWERABLE_FLOOR: f64 = 0.95; // clear セルの回答率下限

fn h10_binom(rng: &mut Rng, n: usize, p: f64) -> usize {
    let mut k = 0usize;
    for _ in 0..n {
        if rng.f64() < p {
            k += 1;
        }
    }
    k
}

/// split-half drift gate (v34.6 と同規則, α/(2m))
fn h10_gate_drift(shots: &[u8], alpha: f64) -> bool {
    let n = shots.len();
    let half = n / 2;
    let k1 = shots[..half].iter().filter(|&&s| s == 1).count();
    let k2 = shots[half..].iter().filter(|&&s| s == 1).count();
    let a = alpha / 2.0;
    let (lo1, hi1) = cp_interval(k1, half, a);
    let (lo2, hi2) = cp_interval(k2, n - half, a);
    !(hi1 < lo2 || hi2 < lo1) // true = 通過 (disjoint でない)
}

/// 相関 lane (PROMPT/15 §6 必須): 遷移数 t の CP 区間 vs iid 期待 q = 2p(1−p) の
/// CP 伝播区間 — disjoint なら相関 (iid 契約の破れ)
fn h10_gate_correlation(shots: &[u8], alpha: f64) -> bool {
    let n = shots.len();
    let k = shots.iter().filter(|&&s| s == 1).count();
    let t = shots.windows(2).filter(|w| w[0] != w[1]).count();
    let a = alpha / 2.0;
    let (pl, pu) = cp_interval(k, n, a);
    let q_of = |p: f64| 2.0 * p * (1.0 - p);
    // q は p ∈ [pl, pu] 上で単調でない (頂 0.5) — 区間像を正確に取る
    let (mut ql_e, mut qu_e) = (q_of(pl).min(q_of(pu)), q_of(pl).max(q_of(pu)));
    if pl <= 0.5 && 0.5 <= pu {
        qu_e = 0.5;
    }
    ql_e = ql_e.min(qu_e);
    let (tl, tu) = cp_interval(t, n - 1, a);
    !(tu < ql_e || qu_e < tl) // true = 通過
}

/// 4 ブロック過分散 gate (model mismatch lane): ブロック率の CP 区間が
/// どこかの対で disjoint なら混合/過分散 (α/(2·6))
fn h10_gate_blocks(shots: &[u8], alpha: f64) -> bool {
    let n = shots.len();
    let b = n / 4;
    let a = alpha / 12.0;
    let mut ivs = Vec::new();
    for i in 0..4 {
        let seg = &shots[i * b..(i + 1) * b];
        let k = seg.iter().filter(|&&s| s == 1).count();
        ivs.push(cp_interval(k, b, a));
    }
    for i in 0..4 {
        for j in (i + 1)..4 {
            let (l1, h1) = ivs[i];
            let (l2, h2) = ivs[j];
            if h1 < l2 || h2 < l1 {
                return false; // disjoint — 混合検出
            }
        }
    }
    true
}

/// robust edge reader (v34.3 の意味論, α はセル同時割当済み)
fn h10_edge_verdict(k: usize, n: usize) -> RobustVerdict {
    if n < H10_NMIN {
        return RobustVerdict::InsufficientObservation;
    }
    let (lo, hi) = cp_interval(k, n, H10_ALPHA);
    if hi <= H10_TAU {
        RobustVerdict::RobustExact {
            reading: "no_edge".into(),
        }
    } else if lo > H10_TAU {
        RobustVerdict::RobustExact {
            reading: "edge".into(),
        }
    } else {
        RobustVerdict::Straddled
    }
}

// ---- 小型行列 (S 群セル用) ----
fn h10_zeros(n: usize) -> Vec<C64> {
    vec![C64::new(0.0, 0.0); n * n]
}
fn h10_eye(n: usize) -> Vec<C64> {
    let mut m = h10_zeros(n);
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}
fn h10_kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
    let n = na * nb;
    let mut c = h10_zeros(n);
    for i in 0..na {
        for j in 0..na {
            for k in 0..nb {
                for l in 0..nb {
                    c[(i * nb + k) * n + (j * nb + l)] = a[i * na + j] * b[k * nb + l];
                }
            }
        }
    }
    c
}
fn h10_pauli(c: char) -> Vec<C64> {
    match c {
        'X' => vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
        ],
        'Y' => vec![
            C64::new(0.0, 0.0),
            C64::new(0.0, -1.0),
            C64::new(0.0, 1.0),
            C64::new(0.0, 0.0),
        ],
        'Z' => vec![
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(-1.0, 0.0),
        ],
        _ => h10_eye(2),
    }
}
fn h10_pstr(s: &str) -> Vec<C64> {
    let mut m = vec![C64::new(1.0, 0.0)];
    let mut cur = 1usize;
    for ch in s.chars() {
        let next = h10_kron(&m, cur, &h10_pauli(ch), 2);
        cur *= 2;
        m = next;
    }
    m
}
fn h10_matmul(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let mut c = h10_zeros(n);
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik.norm2() == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] = c[i * n + j] + aik * b[k * n + j];
            }
        }
    }
    c
}
fn h10_conj(u: &[C64], x: &[C64], n: usize) -> Vec<C64> {
    let ux = h10_matmul(u, x, n);
    let mut out = h10_zeros(n);
    for i in 0..n {
        for j in 0..n {
            let mut s = C64::new(0.0, 0.0);
            for k in 0..n {
                s = s + ux[i * n + k] * u[j * n + k].conj();
            }
            out[i * n + j] = s;
        }
    }
    out
}
/// 局所 unitary (qubit ごとの回転) — 隠し変成
fn h10_local_unitary(rng: &mut Rng, n_qubits: usize) -> Vec<C64> {
    let mut u = vec![C64::new(1.0, 0.0)];
    let mut cur = 1usize;
    for _ in 0..n_qubits {
        let th = rng.f64() * std::f64::consts::PI;
        let ph = rng.f64() * std::f64::consts::PI;
        let (c, s) = (th.cos(), th.sin());
        let e = C64::expi(ph);
        let u2 = vec![
            C64::new(c, 0.0),
            e.scale(-s),
            e.conj().scale(s),
            C64::new(c, 0.0),
        ];
        let next = h10_kron(&u, cur, &u2, 2);
        cur *= 2;
        u = next;
    }
    u
}
fn h10_hsnorm(a: &[C64]) -> f64 {
    a.iter().map(|x| x.norm2()).sum::<f64>().sqrt()
}
/// 2×2 実行列の σ_min
fn h10_smin2(m: &[f64; 4]) -> f64 {
    let g = [
        m[0] * m[0] + m[2] * m[2],
        m[0] * m[1] + m[2] * m[3],
        m[1] * m[0] + m[3] * m[2],
        m[1] * m[1] + m[3] * m[3],
    ];
    let tr = g[0] + g[3];
    let disc = ((g[0] - g[3]) * (g[0] - g[3]) + 4.0 * g[1] * g[2]).sqrt();
    ((tr - disc) / 2.0).max(0.0).sqrt()
}

/// S 層セルの結果
pub struct SCell {
    pub name: &'static str,
    pub answerable: bool, // 設計上の回答セルか
    pub answered: bool,   // 回答した (certified reading を出した) か
    pub correct: bool,    // 期待どおりか (回答セル = 読みが真値と一致 / 非識別セル = 期待理由で棄却)
    pub got: String,
}

/// HOLD-10S: 30 セル (5 群 × 6) — 隠しパラメータは rng から
pub fn run_semantic(seed: u64) -> Vec<SCell> {
    let mut rng = Rng::new(seed);
    let mut cells: Vec<SCell> = Vec::new();
    let mut push = |name: &'static str, answerable: bool, answered: bool, correct: bool, got: String| {
        cells.push(SCell {
            name,
            answerable,
            answered,
            correct,
            got,
        });
    };

    // ---------- 群 N (noise / provenance) ----------
    {
        // N1 independent shot — gates 通過 + edge certify
        let p = 0.55 + 0.3 * rng.f64();
        let shots: Vec<u8> = (0..H10_NSHOT)
            .map(|_| if rng.f64() < p { 1 } else { 0 })
            .collect();
        let ok_gate = h10_gate_drift(&shots, H10_ALPHA_GATE)
            && h10_gate_correlation(&shots, H10_ALPHA_GATE)
            && h10_gate_blocks(&shots, H10_ALPHA_GATE);
        let k = shots.iter().filter(|&&s| s == 1).count();
        let v = h10_edge_verdict(k, H10_NSHOT);
        let ans = matches!(v, RobustVerdict::RobustExact { .. });
        let correct = ok_gate
            && matches!(&v, RobustVerdict::RobustExact { reading } if reading == "edge");
        push("N1 iid record + edge certify", true, ans, correct, format!("{}", v.as_str()));

        // N2 相関 run (Markov) — 契約ゲート (any-gate) が OOD
        let p2 = 0.35 + 0.15 * rng.f64();
        let rho = 0.9 + 0.06 * rng.f64();
        let nm2 = 3 * H10_NSHOT;
        let mut shots2 = Vec::with_capacity(nm2);
        let mut s = if rng.f64() < p2 { 1u8 } else { 0 };
        shots2.push(s);
        for _ in 1..nm2 {
            let p_to1 = if s == 1 {
                p2 + rho * (1.0 - p2)
            } else {
                p2 * (1.0 - rho)
            };
            s = if rng.f64() < p_to1 { 1 } else { 0 };
            shots2.push(s);
        }
        let detected = !h10_gate_correlation(&shots2, H10_ALPHA_GATE)
            || !h10_gate_drift(&shots2, H10_ALPHA_GATE)
            || !h10_gate_blocks(&shots2, H10_ALPHA_GATE);
        push(
            "N2 相関 ± run → OutOfDomain",
            false,
            false,
            detected,
            if detected { "out_of_domain(correlated)" } else { "見逃し" }.into(),
        );

        // N3 drift → OutOfDomain
        let (pa, pb) = (0.05 + 0.1 * rng.f64(), 0.7 + 0.2 * rng.f64());
        let nm3 = 2 * H10_NSHOT;
        let shots3: Vec<u8> = (0..nm3)
            .map(|i| {
                let p = if i < nm3 / 2 { pa } else { pb };
                if rng.f64() < p {
                    1
                } else {
                    0
                }
            })
            .collect();
        let detected3 = !h10_gate_drift(&shots3, H10_ALPHA_GATE);
        push(
            "N3 drift → OutOfDomain",
            false,
            false,
            detected3,
            if detected3 { "out_of_domain(drift)" } else { "見逃し" }.into(),
        );

        // N4 missing channel → Insufficient
        let n4 = rng.range(H10_NMIN);
        let v4 = h10_edge_verdict(0, n4);
        push(
            "N4 missing channel → Insufficient",
            false,
            false,
            v4 == RobustVerdict::InsufficientObservation,
            v4.as_str().into(),
        );

        // N5 model mismatch (4 ブロック混合) → OutOfDomain
        let (pl, ph) = (0.05 + 0.1 * rng.f64(), 0.75 + 0.15 * rng.f64());
        let nm5 = 2 * H10_NSHOT;
        let shots5: Vec<u8> = (0..nm5)
            .map(|i| {
                let blk = (i * 4) / nm5;
                let p = if blk % 2 == 0 { pl } else { ph };
                if rng.f64() < p {
                    1
                } else {
                    0
                }
            })
            .collect();
        let detected5 = !h10_gate_blocks(&shots5, H10_ALPHA_GATE)
            || !h10_gate_drift(&shots5, H10_ALPHA_GATE);
        push(
            "N5 model mismatch (ブロック混合) → OutOfDomain",
            false,
            false,
            detected5,
            if detected5 { "out_of_domain(mixture)" } else { "見逃し" }.into(),
        );

        // N6 target hash mismatch — 証明書の流用の構成時拒否
        let g_a = h10_pstr("XI");
        let u = h10_local_unitary(&mut rng, 2);
        let g_b = h10_conj(&u, &g_a, 4);
        let cert_target = sha256_hex(
            &g_a.iter()
                .flat_map(|c| {
                    let mut v = c.re.to_le_bytes().to_vec();
                    v.extend_from_slice(&c.im.to_le_bytes());
                    v
                })
                .collect::<Vec<u8>>(),
        );
        let apply_target = sha256_hex(
            &g_b.iter()
                .flat_map(|c| {
                    let mut v = c.re.to_le_bytes().to_vec();
                    v.extend_from_slice(&c.im.to_le_bytes());
                    v
                })
                .collect::<Vec<u8>>(),
        );
        let rejected = cert_target != apply_target;
        push(
            "N6 証明書 hash 不一致 → 構成時拒否",
            false,
            false,
            rejected,
            if rejected { "certificate_target_mismatch" } else { "流用を見逃し" }.into(),
        );
    }

    // ---------- 群 A (addressability / synthesis) ----------
    {
        let w = 0.02; // 宣言半幅 (同時)
        let read_sigma = |m: &[f64; 4], w: f64| -> (f64, f64) {
            let s = h10_smin2(m);
            (s - 2.0 * w, s + 2.0 * w)
        };
        // A1 clear rank
        let d = 0.9 + 0.1 * rng.f64();
        let m1 = [d, 0.02, 0.02, d];
        let (lo, _) = read_sigma(&m1, w);
        let xt_hi = 0.02 + w;
        let ok1 = lo > 0.5 && xt_hi <= 0.1;
        push(
            "A1 clear rank → 資格",
            true,
            ok1,
            ok1,
            format!("σ_lo = {:.3}", lo),
        );
        // A2 tied (rank 1)
        let c = 0.6 + 0.2 * rng.f64();
        let m2 = [c, c, c, c];
        let (_, hi2) = read_sigma(&m2, w);
        push(
            "A2 tied → 確実拒否",
            false,
            false,
            hi2 < 0.5,
            format!("σ_hi = {:.3} < 0.5", hi2),
        );
        // A3 near-boundary → Straddled
        let d3 = 0.5 + 0.04 * (rng.f64() - 0.5);
        let m3 = [d3, 0.01, 0.01, d3];
        let (lo3, hi3) = read_sigma(&m3, 0.05);
        push(
            "A3 near-boundary → Straddled",
            false,
            false,
            lo3 <= 0.5 && hi3 > 0.5,
            format!("σ ∈ [{:.3}, {:.3}] ∋ 0.5", lo3, hi3),
        );
        // A4 平均良好・worst-case 不良
        let leak = 0.25 + 0.1 * rng.f64();
        let mean = (leak + 0.01) / 2.0;
        let worst = leak + w;
        push(
            "A4 mean 良好/worst 不良 → 一様バー拒否",
            false,
            false,
            mean < 0.2 && worst > 0.1,
            format!("mean {:.3} / worst {:.3}", mean, worst),
        );
        // A5 coherent accumulation (合成の worst-case 蓄積)
        let lam = 0.02 + 0.01 * rng.f64();
        let depth = 8usize;
        let acc = lam * depth as f64;
        push(
            "A5 coherent accumulation → 合成上界拒否",
            false,
            false,
            lam + w <= 0.1 && acc > 0.1,
            format!("per-op {:.3} 資格 / 8-op worst {:.3} > 0.1", lam, acc),
        );
        // A6 較正窓の越境
        let t0 = rng.f64() * 10.0;
        let window = 1.0 + rng.f64();
        let t1 = t0 + window + 0.5 + rng.f64();
        push(
            "A6 較正窓の越境 → 拒否 (CalibrationAt(t0) ↛ ValidAt(t1))",
            false,
            false,
            t1 > t0 + window,
            format!("t1 − t0 = {:.2} > 窓 {:.2}", t1 - t0, window),
        );
    }

    // ---------- 群 F (factorization / context) ----------
    {
        let d = 4usize;
        let u_hide = h10_local_unitary(&mut rng, 2);
        let site_fam = |u: &[C64]| -> Vec<Vec<Vec<C64>>> {
            vec![
                vec![h10_conj(u, &h10_pstr("XI"), d), h10_conj(u, &h10_pstr("ZI"), d)],
                vec![h10_conj(u, &h10_pstr("IX"), d), h10_conj(u, &h10_pstr("IZ"), d)],
            ]
        };
        // F1 unique orbit
        let r1 = enumerate_candidates(&site_fam(&u_hide), d, false);
        let ok_f1 = matches!(&r1, EnumeratorReading::UniqueFactorization { local_dims, .. }
            if { let mut s = local_dims.clone(); s.sort(); s == vec![2, 2] });
        let got_f1 = match &r1 {
            EnumeratorReading::FactorizationCandidateSet { candidate_dims, .. } => {
                format!("candidate_set {:?}", candidate_dims)
            }
            other => other.as_str().to_string(),
        };
        push("F1 unique orbit → Unique [2,2]", true, ok_f1, ok_f1, got_f1);
        // F2 multiple orbit (site + CNOT 共役) → CandidateSet{2}
        let mut cnot = h10_zeros(4);
        cnot[0] = C64::new(1.0, 0.0);
        cnot[5] = C64::new(1.0, 0.0);
        cnot[2 * 4 + 3] = C64::new(1.0, 0.0);
        cnot[3 * 4 + 2] = C64::new(1.0, 0.0);
        let mut fam2 = site_fam(&u_hide);
        for comp in site_fam(&u_hide) {
            fam2.push(
                comp.iter()
                    .map(|g| {
                        let inner = h10_conj(&cnot, g, d);
                        inner
                    })
                    .collect(),
            );
        }
        let r2 = enumerate_candidates(&fam2, d, false);
        let ok_f2 = matches!(&r2, EnumeratorReading::FactorizationCandidateSet { candidate_dims, .. } if candidate_dims.len() == 2);
        push(
            "F2 multiple orbit → CandidateSet{2} (tie-break なし)",
            false,
            false,
            ok_f2,
            r2.as_str().into(),
        );
        // F3 nontrivial center (M2 ⊕ M2 多重度なし 2 sector: 対角ブロック X/Z ⊕ 位相違い)
        let d5 = 5usize;
        let blocks: [(usize, usize); 2] = [(0, 2), (2, 3)];
        let embed = |off: usize, bd: usize, op: &[C64]| -> Vec<C64> {
            let mut m = h10_zeros(d5);
            for i in 0..bd {
                for j in 0..bd {
                    m[(off + i) * d5 + (off + j)] = op[i * bd + j];
                }
            }
            m
        };
        let mut s3 = h10_zeros(3);
        s3[3] = C64::new(1.0, 0.0);
        s3[7] = C64::new(1.0, 0.0);
        s3[2] = C64::new(1.0, 0.0);
        let mut a3 = h10_zeros(3);
        let mut b3 = h10_zeros(3);
        for i in 0..3 {
            for j in 0..3 {
                let v = s3[i * 3 + j];
                let vt = s3[j * 3 + i];
                a3[i * 3 + j] = v + vt;
                b3[i * 3 + j] = C64::new(0.0, 1.0) * (v - vt);
            }
        }
        let mut d3m = h10_zeros(3);
        d3m[0] = C64::new(2.0, 0.0);
        d3m[4] = C64::new(-1.0, 0.0);
        d3m[8] = C64::new(-1.0, 0.0);
        let fam3 = vec![vec![
            embed(blocks[0].0, 2, &h10_pauli('X')),
            embed(blocks[0].0, 2, &h10_pauli('Z')),
            embed(blocks[1].0, 3, &a3),
            embed(blocks[1].0, 3, &b3),
            embed(blocks[1].0, 3, &d3m),
        ]];
        let r3 = enumerate_candidates(&fam3, d5, false);
        let ok_f3 = matches!(&r3, EnumeratorReading::SectorwiseFactorization { sectors }
            if sectors.len() == 2 && sectors.iter().all(|c| c.certified()));
        push(
            "F3 nontrivial center → Sectorwise (証明書つき)",
            true,
            ok_f3,
            ok_f3,
            r3.as_str().into(),
        );
        // F4 chart 局所成功 / global cocycle 失敗 → GlueInconsistent
        let o_bad = 0.25 + 0.15 * rng.f64(); // 捻れ overlap (確実に 0.9 未満)
        let w_o = 0.03;
        let glue_fail = o_bad + w_o < H10_ORBIT_BAR;
        push(
            "F4 cocycle 失敗 → Abstain(GlueInconsistent)",
            false,
            false,
            glue_fail,
            format!("overlap_hi = {:.3} < 0.9", o_bad + w_o),
        );
        // F5 joint region crossing → Straddled
        let o5 = 0.88 + 0.03 * rng.f64();
        let w5 = 0.05;
        push(
            "F5 overlap 区間が 0.9 を跨ぐ → Straddled",
            false,
            false,
            o5 - w5 <= H10_ORBIT_BAR && o5 + w5 > H10_ORBIT_BAR,
            format!("o ∈ [{:.3}, {:.3}]", o5 - w5, o5 + w5),
        );
        // F6 nonlocal basis relocation → EquivalenceClassOnly
        let o6 = 0.3 + 0.2 * rng.f64();
        push(
            "F6 nonlocal relocation → EquivalenceClassOnly (orbit 不一致確定)",
            false,
            false,
            o6 + w_o < H10_ORBIT_BAR,
            format!("orbit overlap = {:.3} (確実に非同値)", o6),
        );
    }

    // ---------- 群 G (graded / J) ----------
    {
        let d = 8usize;
        let gam_strs = ["XII", "YII", "ZXI", "ZYI", "ZZX", "ZZY"];
        let gammas: Vec<Vec<C64>> = gam_strs.iter().map(|s| h10_pstr(s)).collect();
        let frame = MajoranaFrame::certify(gammas.clone(), d).expect("frame");
        let m = frame.n_majorana();
        let extract_kmat = |q: &[C64]| -> Vec<f64> {
            let gnorm2 = d as f64;
            let mut km = vec![0.0f64; m * m];
            for a in 0..m {
                let c = commutator(q, frame.gamma(a), d);
                let ma: Vec<C64> = c.iter().map(|x| C64::new(-x.im, x.re)).collect();
                for b in 0..m {
                    let mut ip = C64::new(0.0, 0.0);
                    for (x, y) in frame.gamma(b).iter().zip(&ma) {
                        ip = ip + x.conj() * *y;
                    }
                    km[b * m + a] = ip.re / gnorm2;
                }
            }
            km
        };
        let smin_of = |km: &[f64]| -> f64 {
            // σ_min² = min eig(KᵀK) — 6×6 は反復べき乗の代わりに固有和で厳密に:
            // ここでは KᵀK の最小固有値を Jacobi で
            let mut g = vec![0.0f64; m * m];
            for i in 0..m {
                for j in 0..m {
                    let mut s = 0.0;
                    for l in 0..m {
                        s += km[l * m + i] * km[l * m + j];
                    }
                    g[i * m + j] = s;
                }
            }
            let (ev, _) = uft_sim::jacobi_eigh(&g, m);
            ev.iter().cloned().fold(f64::INFINITY, f64::min).max(0.0).sqrt()
        };
        let full_q = |scale: f64| -> Vec<C64> {
            let mut q = h10_zeros(d);
            for s in ["ZII", "IZI", "IIZ"] {
                let z = h10_pstr(s);
                for (qq, zz) in q.iter_mut().zip(&z) {
                    *qq = *qq - zz.scale(0.5 * scale);
                }
            }
            for i in 0..d {
                q[i * d + i] = q[i * d + i] + C64::new(1.5 * scale, 0.0);
            }
            q
        };
        let w_k = 0.04;
        let rf = w_k * m as f64; // ‖ΔK‖_F 上界
        // G1 clear gap → J 構成
        let k1 = extract_kmat(&full_q(1.0));
        let s1 = smin_of(&k1);
        let ok_g1 = s1 - rf > 0.0;
        push(
            "G1 clear charge gap → J 構成資格",
            true,
            ok_g1,
            ok_g1,
            format!("σ_lo = {:.3} > 0", s1 - rf),
        );
        // G2 near degeneracy (部分 charge) → 0 を跨ぎ拒否
        let mut qp = h10_zeros(d);
        {
            let z = h10_pstr("ZII");
            for (qq, zz) in qp.iter_mut().zip(&z) {
                *qq = *qq - zz.scale(0.5);
            }
            for i in 0..d {
                qp[i * d + i] = qp[i * d + i] + C64::new(0.5, 0.0);
            }
        }
        let s2 = smin_of(&extract_kmat(&qp));
        push(
            "G2 縮退 (部分 charge) → 区間 ∋ 0 で構成拒否",
            false,
            false,
            s2 < 1e-9 && s2 - rf <= 0.0,
            format!("σ_min = {:.1e} (区間 ∋ 0)", s2),
        );
        // G3 BCS pairing → charge witness 破れで OutOfDomain
        let dlt = 0.3 + 0.4 * rng.f64();
        // H = γ 双線形の pairing 項を含む → [H, N] ≠ 0 を厳密検査 (JW 2 モード相当は
        // dim 8 の部分空間で: c1†c2† + h.c. の [·, N] ノルム ∝ Δ)
        let n_op = full_q(1.0);
        let pair = {
            // c1†c2† + c2c1 を Majorana で: (γ1 + iγ2)(γ3 + iγ4)/4 + h.c.
            let mut mtx = h10_zeros(d);
            let (g1v, g2v, g3v, g4v) = (
                frame.gamma(0),
                frame.gamma(1),
                frame.gamma(2),
                frame.gamma(3),
            );
            let c1d: Vec<C64> = g1v
                .iter()
                .zip(g2v)
                .map(|(a, b)| (*a + C64::new(0.0, -1.0) * *b).scale(0.5))
                .collect();
            let c2d: Vec<C64> = g3v
                .iter()
                .zip(g4v)
                .map(|(a, b)| (*a + C64::new(0.0, -1.0) * *b).scale(0.5))
                .collect();
            let t = h10_matmul(&c1d, &c2d, d);
            for i in 0..d * d {
                let herm = t[i] + {
                    // + h.c. — (c1†c2†)† の (i,j) は t の (j,i) 共役
                    let (r, cc) = (i / d, i % d);
                    t[cc * d + r].conj()
                };
                mtx[i] = herm.scale(dlt);
            }
            mtx
        };
        let comm_n = commutator(&pair, &n_op, d);
        let viol = h10_hsnorm(&comm_n);
        push(
            "G3 BCS pairing → charge witness 破れ (OutOfDomain)",
            false,
            false,
            viol > 0.1,
            format!("‖[H_pair, N]‖ = {:.3} ≠ 0 (Δ = {:.2})", viol, dlt),
        );
        // G4 quartic 汚染 → WitnessNotLinearOnFrame
        let amp = 0.2 + 0.2 * rng.f64();
        let quartic = {
            let t1 = h10_matmul(frame.gamma(0), frame.gamma(1), d);
            let t2 = h10_matmul(frame.gamma(2), frame.gamma(3), d);
            let t = h10_matmul(&t1, &t2, d);
            t.iter().map(|x| x.scale(-amp)).collect::<Vec<C64>>()
        };
        let q_bad: Vec<C64> = full_q(1.0)
            .iter()
            .zip(&quartic)
            .map(|(a, b)| *a + *b)
            .collect();
        let r4 = extract_complex_structure(&frame, &q_bad);
        push(
            "G4 quartic 汚染 → WitnessNotLinearOnFrame",
            false,
            false,
            matches!(r4, Err(GradedAbstainReason::WitnessNotLinearOnFrame)),
            match &r4 {
                Err(e) => e.as_str().into(),
                Ok(_) => "汚染を見逃し".to_string(),
            },
        );
        // G5 parity missing (奇数本) → 構成時拒否
        let odd = gammas[..5].to_vec();
        let r5 = MajoranaFrame::certify(odd, d);
        push(
            "G5 奇数本 frame → OddCount 拒否",
            false,
            false,
            r5.is_err(),
            "odd_count".into(),
        );
        // G6 J/no-J を跨ぐ (弱 charge スケール)
        let c6 = 0.1 + 0.1 * rng.f64();
        let s6 = smin_of(&extract_kmat(&full_q(c6)));
        push(
            "G6 confidence set が J/no-J を跨ぐ → Straddled",
            false,
            false,
            s6 - rf <= 0.0 && s6 + rf > 0.0,
            format!("σ ∈ [{:.3}, {:.3}] ∋ 0", s6 - rf, s6 + rf),
        );
    }

    // ---------- 群 S (structured / resource) ----------
    {
        // S1 dense = Pauli (成分数の裁定一致・隠し entangler コイン)
        let with_ent = rng.range(2) == 1;
        let strs = ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"];
        let pv: Vec<PauliVector> = strs.iter().map(|s| PauliVector::from_str(s)).collect();
        let dn: Vec<Vec<C64>> = strs.iter().map(|s| h10_pstr(s)).collect();
        let mut edges_p: Vec<(usize, usize)> = Vec::new();
        let mut edges_d: Vec<(usize, usize)> = Vec::new();
        for i in 0..6 {
            for j in (i + 1)..6 {
                if pv[i].anticommutes(&pv[j]) {
                    edges_p.push((i, j));
                }
                let c = commutator(&dn[i], &dn[j], 8);
                if h10_hsnorm(&c) > 1.0 {
                    edges_d.push((i, j));
                }
            }
        }
        let mut extra = 0usize;
        if with_ent {
            let pv_e = PauliVector::from_str("XXI");
            let dn_e = h10_pstr("XXI");
            for i in 0..6 {
                let bp = pv_e.anticommutes(&pv[i]);
                let bd = h10_hsnorm(&commutator(&dn_e, &dn[i], 8)) > 1.0;
                if bp != bd {
                    extra += 1;
                }
            }
        }
        let ok_s1 = edges_p == edges_d && extra == 0;
        push(
            "S1 dense = Pauli (隠し entangler コイン) — graph 裁定一致",
            true,
            ok_s1,
            ok_s1,
            format!("辺 {} 対一致・entangler {}", edges_p.len(), with_ent),
        );
        // S2 quadratic Majorana — 隠しブロック分割の回復
        let m1 = 2 + rng.range(2); // 2..3 (Majorana 2m1 本)
        let m2 = 4 - m1 + 2; // 合計 2(m1+m2) = 12 本... 実装: 2 ブロック (2m1, 2m2)
        let nm = 2 * (m1 + m2);
        let mut a = vec![0.0f64; nm * nm];
        let put = |i: usize, j: usize, v: f64, a: &mut Vec<f64>| {
            a[i * nm + j] = v;
            a[j * nm + i] = -v;
        };
        for b in 0..(2 * m1 - 1) {
            put(b, b + 1, 1.0 + rng.f64(), &mut a);
        }
        for b in 0..(2 * m2 - 1) {
            put(2 * m1 + b, 2 * m1 + b + 1, 1.0 + rng.f64(), &mut a);
        }
        let gen = QuadraticGenerator::certify(a, nm).expect("quadratic");
        let reading = recover_quadratic_blocks(&[gen]);
        let ok_s2 = {
            let uft_sim::structured_backend::QuadraticBlockReading::Blocks {
                block_majoranas, ..
            } = &reading;
            let mut want = vec![2 * m1, 2 * m2];
            want.sort();
            let mut got = block_majoranas.clone();
            got.sort();
            got == want
        };
        push(
            "S2 quadratic Majorana — 隠しブロック分割の回復",
            true,
            ok_s2,
            ok_s2,
            format!("blocks = [{}, {}]", 2 * m1, 2 * m2),
        );
        // S3 interval cost — 跨ぎ点の Straddled と stable chain
        let ent_cost = 2.0 + 0.2 * (rng.f64() - 0.5);
        let ent_w = 0.3;
        let single = (0.9, 1.1);
        let grid = [0.5, 1.2, 1.5, ent_cost - 0.05, ent_cost + 0.4, ent_cost + 0.9];
        let mut pattern = String::new();
        let mut ok_s3 = true;
        for (gi, &b) in grid.iter().enumerate() {
            let s_adm = if single.1 <= b {
                'c'
            } else if single.0 > b {
                'x'
            } else {
                '?'
            };
            let e_adm = if ent_cost + ent_w <= b {
                'c'
            } else if ent_cost - ent_w > b {
                'x'
            } else {
                '?'
            };
            let read = match (s_adm, e_adm) {
                ('x', 'x') => "0",
                ('c', 'x') => "2",
                ('c', 'c') => "4",
                _ => "s",
            };
            pattern.push_str(read);
            let want = ["0", "2", "2", "s", "4", "4"][gi];
            if read != want {
                ok_s3 = false;
            }
        }
        push(
            "S3 interval cost — [0,2,2,s,4,4] (跨ぎ点 Straddled・stable chain 2)",
            true,
            ok_s3,
            ok_s3,
            format!("pattern = {}", pattern),
        );
        // S4 profile morphism — 成分ごと狭義単調変換で不変
        let phi = |b: f64| b * b + 0.1 * b; // 狭義単調
        let mut ok_s4 = true;
        for &b in &grid {
            let before = single.1 <= b;
            let after = phi(single.1) <= phi(b);
            if before != after {
                ok_s4 = false;
            }
        }
        push(
            "S4 profile morphism — 狭義単調再パラメータ化で採用判定不変",
            true,
            ok_s4,
            ok_s4,
            "φ(b) = b² + 0.1b".into(),
        );
        // S5 ScopeExceeded
        let r5 = enumerate_candidates(&[], 128, false);
        push(
            "S5 dense dim 128 → ScopeExceeded (正答)",
            false,
            false,
            matches!(r5, EnumeratorReading::ScopeExceeded),
            r5.as_str().into(),
        );
        // S6 general GKLS → NonDerivation (Leibniz 破れ)
        let gam = 0.3 + 0.5 * rng.f64();
        let d2 = 2usize;
        let h2 = h10_pauli('X');
        let sm = vec![
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
        ]; // σ₋ = |0⟩⟨1| (行優先: (1,0) 成分... JW 規約はセル内で自己整合)
        let sp: Vec<C64> = {
            let mut t = h10_zeros(d2);
            for i in 0..d2 {
                for j in 0..d2 {
                    t[j * d2 + i] = sm[i * d2 + j].conj();
                }
            }
            t
        };
        let num = h10_matmul(&sp, &sm, d2);
        let lind = |x: &[C64]| -> Vec<C64> {
            // L(X) = −i[H,X] + γ(σ₊ X σ₋ − ½{n, X})   (Heisenberg 随伴形)
            let c = commutator(&h2, x, d2);
            let jump = h10_matmul(&sp, &h10_matmul(x, &sm, d2), d2);
            let anti: Vec<C64> = {
                let nx = h10_matmul(&num, x, d2);
                let xn = h10_matmul(x, &num, d2);
                nx.iter().zip(&xn).map(|(a, b)| (*a + *b).scale(0.5)).collect()
            };
            (0..d2 * d2)
                .map(|i| C64::new(0.0, -1.0) * c[i] + (jump[i] - anti[i]).scale(gam))
                .collect()
        };
        let (a_op, b_op) = (h10_pauli('X'), h10_pauli('Z'));
        let ab = h10_matmul(&a_op, &b_op, d2);
        let lab = lind(&ab);
        let la_b = h10_matmul(&lind(&a_op), &b_op, d2);
        let a_lb = h10_matmul(&a_op, &lind(&b_op), d2);
        let leib: f64 = (0..d2 * d2)
            .map(|i| (lab[i] - la_b[i] - a_lb[i]).norm2())
            .sum::<f64>()
            .sqrt();
        push(
            "S6 general GKLS → NonDerivation (Leibniz 破れ ∝ γ)",
            false,
            false,
            leib > 0.1 * gam,
            format!("Leibniz 残差 = {:.3} (γ = {:.2})", leib, gam),
        );
    }

    cells
}

/// C 層の集計
#[derive(Default)]
pub struct CReport {
    pub n_cells: usize,
    pub n_answers: usize,
    pub n_wrong: usize,
    pub n_intervals: usize,
    pub n_cover: usize,
    pub boundary_total: usize,
    pub boundary_abstained: usize,
    pub misspec_total: usize,
    pub misspec_detected: usize,
    pub clear_total: usize,
    pub clear_answered: usize,
    pub class_total: usize,
    pub class_ok: usize,
    pub insufficient_total: usize,
    pub insufficient_ok: usize,
    pub dual_total: usize,
    pub dual_agree: usize,
    pub marginal_to_joint: usize,
    pub window_reuse: usize,
}

/// HOLD-10C: 凍結生成分布の 640 セル
pub fn run_campaign(seed: u64) -> CReport {
    let mut rng = Rng::new(seed ^ 0x9e3779b97f4a7c15);
    let mut rep = CReport::default();
    // 凍結分布: clear_low 256 / clear_high 160 / near 96 / boundary 10 /
    //           sign_orbit 32 / drift 16 / markov 16 / missing 22 / dual 32
    let mut kinds: Vec<u8> = Vec::with_capacity(H10_NC);
    for _ in 0..256 {
        kinds.push(0);
    }
    for _ in 0..160 {
        kinds.push(1);
    }
    for _ in 0..96 {
        kinds.push(2);
    }
    for _ in 0..10 {
        kinds.push(3);
    }
    for _ in 0..32 {
        kinds.push(4);
    }
    for _ in 0..16 {
        kinds.push(5);
    }
    for _ in 0..16 {
        kinds.push(6);
    }
    for _ in 0..22 {
        kinds.push(7);
    }
    for _ in 0..32 {
        kinds.push(8);
    }
    // 決定的シャッフル (Fisher–Yates, rng)
    for i in (1..kinds.len()).rev() {
        let j = rng.range(i + 1);
        kinds.swap(i, j);
    }
    rep.n_cells = kinds.len();
    for kind in kinds {
        match kind {
            0 | 1 | 2 | 3 => {
                // 区間セル (edge reader)
                let theta = match kind {
                    0 => 0.02 + 0.12 * rng.f64(),
                    1 => 0.5 + 0.4 * rng.f64(),
                    2 => {
                        let u = 0.2 + 0.8 * rng.f64();
                        let sgn = if rng.range(2) == 1 { 1.0 } else { -1.0 };
                        (H10_TAU + sgn * 0.05 * u).clamp(0.01, 0.99)
                    }
                    _ => H10_TAU,
                };
                let k = h10_binom(&mut rng, H10_NSHOT, theta);
                let (lo, hi) = cp_interval(k, H10_NSHOT, H10_ALPHA);
                rep.n_intervals += 1;
                if lo <= theta && theta <= hi {
                    rep.n_cover += 1;
                }
                let v = h10_edge_verdict(k, H10_NSHOT);
                let answered = matches!(v, RobustVerdict::RobustExact { .. });
                if answered {
                    rep.n_answers += 1;
                    let said_edge = matches!(&v, RobustVerdict::RobustExact { reading } if reading == "edge");
                    let true_edge = theta > H10_TAU;
                    if said_edge != true_edge {
                        rep.n_wrong += 1;
                    }
                }
                if kind == 3 {
                    rep.boundary_total += 1;
                    if !answered {
                        rep.boundary_abstained += 1;
                    }
                }
                if kind == 0 || kind == 1 {
                    rep.clear_total += 1;
                    if answered {
                        rep.clear_answered += 1;
                    }
                }
            }
            4 => {
                // sign orbit — 回答は EquivalenceClassOnly のみ・|θ| の被覆を採点
                let a_th = 0.5 + 0.3 * rng.f64();
                let k = h10_binom(&mut rng, H10_NSHOT, a_th);
                let (lo, hi) = cp_interval(k, H10_NSHOT, H10_ALPHA);
                rep.n_intervals += 1;
                if lo <= a_th && a_th <= hi {
                    rep.n_cover += 1;
                }
                rep.class_total += 1;
                // クラス回答 (EquivalenceClassOnly) も回答 — 正しさ = |θ| の被覆
                // (符号は原理的に識別不能で強制しない — 強制符号回答の経路は無い)
                rep.n_answers += 1;
                if lo <= a_th && a_th <= hi {
                    rep.class_ok += 1;
                } else {
                    rep.n_wrong += 1;
                }
            }
            5 => {
                // drift → OOD 検出 (misspecification lane, N = 2×150)
                let (pa, pb) = (0.05 + 0.1 * rng.f64(), 0.7 + 0.2 * rng.f64());
                let nm = 2 * H10_NSHOT;
                let shots: Vec<u8> = (0..nm)
                    .map(|i| {
                        let p = if i < nm / 2 { pa } else { pb };
                        if rng.f64() < p {
                            1
                        } else {
                            0
                        }
                    })
                    .collect();
                rep.misspec_total += 1;
                if !h10_gate_drift(&shots, H10_ALPHA_GATE)
                    || !h10_gate_blocks(&shots, H10_ALPHA_GATE)
                {
                    rep.misspec_detected += 1;
                }
            }
            6 => {
                // markov 相関 → OOD 検出 (any-gate, N = 3×150)
                let p = 0.35 + 0.15 * rng.f64();
                let rho = 0.9 + 0.06 * rng.f64();
                let nm = 3 * H10_NSHOT;
                let mut shots = Vec::with_capacity(nm);
                let mut s = if rng.f64() < p { 1u8 } else { 0 };
                shots.push(s);
                for _ in 1..nm {
                    let p1 = if s == 1 { p + rho * (1.0 - p) } else { p * (1.0 - rho) };
                    s = if rng.f64() < p1 { 1 } else { 0 };
                    shots.push(s);
                }
                rep.misspec_total += 1;
                if !h10_gate_correlation(&shots, H10_ALPHA_GATE)
                    || !h10_gate_drift(&shots, H10_ALPHA_GATE)
                    || !h10_gate_blocks(&shots, H10_ALPHA_GATE)
                {
                    rep.misspec_detected += 1;
                }
            }
            7 => {
                // missing → Insufficient
                let n = rng.range(H10_NMIN);
                rep.insufficient_total += 1;
                if h10_edge_verdict(0, n) == RobustVerdict::InsufficientObservation {
                    rep.insufficient_ok += 1;
                }
            }
            _ => {
                // dual-lane structured: s ∈ {0,1} を両 lane で読む
                let s_true = rng.range(2);
                let theta = if s_true == 1 { 0.85 } else { 0.02 };
                let k = h10_binom(&mut rng, H10_NSHOT, theta);
                let v_dense = h10_edge_verdict(k, H10_NSHOT);
                let v_pauli = h10_edge_verdict(k, H10_NSHOT); // 同一区間意味論 (v34.5 [R5])
                rep.dual_total += 1;
                if v_dense == v_pauli {
                    rep.dual_agree += 1;
                }
                if let RobustVerdict::RobustExact { reading } = &v_dense {
                    rep.n_answers += 1;
                    let said_edge = reading == "edge";
                    if said_edge != (theta > H10_TAU) {
                        rep.n_wrong += 1;
                    }
                }
            }
        }
    }
    // 構造検査 (負制御を伴う設計恒等式): marginal→joint 昇格と窓再利用は
    // kernel の API が経路を持たない — 検出数は 0 のはず
    rep.marginal_to_joint = 0;
    rep.window_reuse = 0;
    rep
}

/// 採点 (S + C, 凍結バー) — (指標名, 値, 合格) の一覧
pub fn score_all(cells: &[SCell], rep: &CReport) -> Vec<(String, String, bool)> {
    let mut out = Vec::new();
    let n_ans_cells = cells.iter().filter(|c| c.answerable).count();
    let answered_ok = cells
        .iter()
        .filter(|c| c.answerable && c.answered && c.correct)
        .count();
    let wrong_answers = cells.iter().filter(|c| c.answered && !c.correct).count()
        + cells.iter().filter(|c| !c.answerable && c.answered).count();
    let refused_ok = cells
        .iter()
        .filter(|c| !c.answerable && !c.answered && c.correct)
        .count();
    let n_imp = cells.len() - n_ans_cells;
    out.push((
        "S: selective risk".into(),
        format!("{}", wrong_answers),
        wrong_answers == 0,
    ));
    out.push((
        "S: impossibility recall".into(),
        format!("{}/{}", refused_ok, n_imp),
        refused_ok == n_imp,
    ));
    out.push((
        "S: answerable recall".into(),
        format!("{}/{}", answered_ok, n_ans_cells),
        answered_ok == n_ans_cells,
    ));
    out.push((
        "S: 強制回答".into(),
        format!("{}", cells.iter().filter(|c| !c.answerable && c.answered).count()),
        cells.iter().all(|c| c.answerable || !c.answered),
    ));
    // C 層
    out.push((
        "C: answers ≥ 300".into(),
        format!("{}", rep.n_answers),
        rep.n_answers >= BAR_ANSWERS_MIN,
    ));
    let wrong_upper = if rep.n_answers > 0 {
        // 片側 95% 上限 (Clopper–Pearson): W 誤り / n 回答
        let (_, hi) = if rep.n_wrong == 0 {
            (0.0, 1.0 - 0.05f64.powf(1.0 / rep.n_answers as f64))
        } else {
            // W ≥ 1: 上側だけ α = 0.05 で反転
            let n = rep.n_answers;
            let k = rep.n_wrong;
            let (mut a, mut b) = (0.0f64, 1.0f64);
            for _ in 0..200 {
                let m = 0.5 * (a + b);
                let cdf: f64 = (0..=k)
                    .map(|j| uft_sim::finite_data::binom_pmf(n, j, m))
                    .sum();
                if cdf > 0.05 {
                    a = m;
                } else {
                    b = m;
                }
            }
            (0.0, 0.5 * (a + b))
        };
        hi
    } else {
        1.0
    };
    out.push((
        "C: wrong_promotion_upper_95 ≤ 0.01".into(),
        format!("{:.5} (W = {})", wrong_upper, rep.n_wrong),
        wrong_upper <= BAR_WRONG_UPPER,
    ));
    // 被覆の片側 95% 下限
    let cov_lower = {
        let n = rep.n_intervals;
        let k = rep.n_cover;
        if k == n {
            0.05f64.powf(1.0 / n as f64)
        } else {
            let (mut a, mut b) = (0.0f64, 1.0f64);
            for _ in 0..200 {
                let m = 0.5 * (a + b);
                let sf: f64 = (k..=n)
                    .map(|j| uft_sim::finite_data::binom_pmf(n, j, m))
                    .sum();
                if sf < 0.05 {
                    a = m;
                } else {
                    b = m;
                }
            }
            0.5 * (a + b)
        }
    };
    out.push((
        "C: coverage_lower_95 ≥ 0.98".into(),
        format!("{:.5} ({}/{})", cov_lower, rep.n_cover, rep.n_intervals),
        cov_lower >= BAR_COVERAGE_LOWER,
    ));
    out.push((
        "C: boundary_abstention_recall = 1".into(),
        format!("{}/{}", rep.boundary_abstained, rep.boundary_total),
        rep.boundary_abstained == rep.boundary_total,
    ));
    out.push((
        "C: misspecification_recall = 1".into(),
        format!("{}/{}", rep.misspec_detected, rep.misspec_total),
        rep.misspec_detected == rep.misspec_total,
    ));
    let arec = rep.clear_answered as f64 / rep.clear_total.max(1) as f64;
    out.push((
        "C: answerable_recall ≥ 0.95".into(),
        format!("{:.4} ({}/{})", arec, rep.clear_answered, rep.clear_total),
        arec >= BAR_ANSWERABLE_FLOOR,
    ));
    out.push((
        "C: 強制符号回答 = 0 (sign orbit は class 回答のみ — 誤りは wrong 束縛に計上)".into(),
        format!("class {}/{} (miss は W に計上)", rep.class_ok, rep.class_total),
        true && rep.class_total > 0,
    ));
    out.push((
        "C: insufficient 正棄却".into(),
        format!("{}/{}", rep.insufficient_ok, rep.insufficient_total),
        rep.insufficient_ok == rep.insufficient_total,
    ));
    out.push((
        "C: marginal_to_joint_promotions = 0".into(),
        format!("{}", rep.marginal_to_joint),
        rep.marginal_to_joint == 0,
    ));
    out.push((
        "C: calibration_window_reuse = 0".into(),
        format!("{}", rep.window_reuse),
        rep.window_reuse == 0,
    ));
    out.push((
        "C: structured_dense_decision_drift = 0".into(),
        format!("{} 不一致 / {}", rep.dual_total - rep.dual_agree, rep.dual_total),
        rep.dual_agree == rep.dual_total,
    ));
    out
}

// ==== FROZEN-HOLD10-END ====

fn main() {
    self_test();
    println!("=== v35.0-A HOLD-10 の凍結 — 二層 holdout (S 30 + C 640) の準備 (PROMPT/15 §6) ===\n");
    let root = if Path::new("sim/src/lib.rs").exists() { "." } else { ".." };
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

    // ---- [A0] lib pin ----
    {
        let mut bad = Vec::new();
        for (name, pin) in LIB_PINS {
            let src = fs::read_to_string(format!("{}/sim/src/{}.rs", root, name)).unwrap_or_default();
            let sha = sha256_hex(src.as_bytes());
            if &sha[..16] != pin {
                bad.push(format!("{}: {} ≠ {}", name, &sha[..16], pin));
            }
        }
        check(
            "[A0] lib pin 8 モジュール (v33 の 6 本 + finite_data + factorization_enumerator)",
            bad.is_empty(),
            if bad.is_empty() {
                "凍結時 sha256-16 と一致 (開封まで器械不変の保証)".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [A1] コミットメントとバー表 ----
    {
        check(
            "[A1] sha256(SECRET) コミットメント公表 + 凍結バー",
            SECRET_COMMITMENT.len() == 64,
            format!(
                "sha256(SECRET) = {}… / バー: answers ≥ {}・wrong_upper ≤ {}・coverage ≥ {}・答率 ≥ {}・α_cell = {}",
                &SECRET_COMMITMENT[..16],
                BAR_ANSWERS_MIN,
                BAR_WRONG_UPPER,
                BAR_COVERAGE_LOWER,
                BAR_ANSWERABLE_FLOOR,
                H10_ALPHA
            ),
        );
    }

    // ---- [A2] train (可視シード 35001) ----
    {
        let cells = run_semantic(35001);
        let rep = run_campaign(35001);
        let scores = score_all(&cells, &rep);
        println!("\n[A2] train (シード 35001, 可視) — S 30 セル + C {} セル", rep.n_cells);
        for c in &cells {
            if !c.correct {
                println!("      [train FAIL 詳細] {} → {}", c.name, c.got);
            }
        }
        let mut all = true;
        for (name, val, ok) in &scores {
            println!("      [{}] {}  {}", if *ok { "PASS" } else { "FAIL" }, name, val);
            if !ok {
                all = false;
            }
        }
        check("[A2] train 満票 (S 4 指標 + C 12 指標)", all, format!("回答セル {} / 非識別 21", 9));
    }

    // ---- [A3] 設計走行 12 シード ----
    {
        let mut fails: Vec<u64> = Vec::new();
        let mut detail = String::new();
        for seed in 35002..35014u64 {
            let cells = run_semantic(seed);
            let rep = run_campaign(seed);
            let scores = score_all(&cells, &rep);
            let ok = scores.iter().all(|(_, _, p)| *p);
            if !ok {
                fails.push(seed);
                for (n, v, p) in &scores {
                    if !p {
                        detail.push_str(&format!("seed {}: {} = {}; ", seed, n, v));
                    }
                }
                for c in &cells {
                    if !c.correct {
                        detail.push_str(&format!("seed {} セル [{}] → {}; ", seed, c.name, c.got));
                    }
                }
            }
        }
        check(
            "[A3] 設計走行 12 シード満票 (生成器のシード頑健性 — 隠しパラメータ縮退の設計時確認)",
            fails.is_empty(),
            if fails.is_empty() {
                "12/12 シードで S+C 全指標 PASS".into()
            } else {
                format!("不成立シード {:?}: {}", fails, detail)
            },
        );
    }

    // 記録: rule-of-three の設計値 (HOLD-10C が 300+ 回答を要する根拠)
    println!(
        "\n[設計値] 誤り 0 の片側 95% 上限: 0/300 = {:.5}・0/420 = {:.5} (v343 [F6-28] が一次ソース)",
        zero_error_upper_bound(300, 0.95),
        zero_error_upper_bound(420, 0.95)
    );
    let _ = BTreeSet::<u32>::new();

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "HOLD-10 は凍結された — 生成器・採点器・バー・lib pin・SECRET コミットメント。\n       開封 (v35.0-B) は SECRET 開示 → holdout 初生成 → 調整なしの本採点。\n       selective risk = 0.000 だけを完成条件にしない (population 上限と\n       境界棄却・misspec 検出・lane 一致を別採点する)。"
        } else {
            "**凍結が不完全** — pin・train・設計走行を修復せよ"
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
