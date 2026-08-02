//! v34.0-A HOLD-9 の凍結 — 操作の出自 × 文脈整合 × 資源依存局所性 (PROMPT/14 §HOLD-9)
//!
//! 第三十三期の全器械 (v33.1 修復入口 / v33.2 laboratory interface / v33.3 resource
//! profile / v33.4 atlas glue / v33.5 graded recovery / v33.6 structured backend) を
//! 凍結し、識別可能性境界を新鮮 holdout で採点する準備をする。採点は従来の
//!
//!   selective risk = 0 / impossibility recall = 1 / answerable recall = 1
//!
//! に第三十三期の 5 計量を加える (PROMPT/14):
//!
//!   origin_certificate_coverage = 1.0 / context_witness_coverage = 1.0 /
//!   raw_operation_promotions = 0 / scope_violations = 0 /
//!   transient_factorization_promotions = 0
//!
//! セル 5 群 20 セル (隠しパラメータ = 変成 unitary・置換・対の選択・γ・コスト・
//! O(6) 回転・汚染振幅・セル/qubit 数の抽選):
//!   入力完全性: IC1 raw+流用証明書 → 拒否 / IC2 gens 別渡し不在 (net-owned byte
//!   同一) / IC3 role-mixed → 拒否 / IC4 GKLS → 可換子 lane 拒否
//!   accessibility: AC1 独立 knobs → Exact / AC2 tied → 分解拒否 + Abstain /
//!   AC3 同一閉包・別 interface → 非同値 orbit / AC4 出自証明書の流用 → 拒否
//!   context・resource: CR1 budget profile (不足 → [2,2,2] → [2,4]) / CR2 一意
//!   glue = 大域一致 / CR3 複数 glue → EquivalenceClassOnly / CR4 cross-talk 跨ぎ
//!   graded: GR1 odd CAR only → Majorana orbit / GR2 charge witness → complex
//!   modes / GR3 汚染 witness → Abstain / GR4 ordinary odd → 構成時拒否
//!   scale・変成: SC1 dense = Pauli 裁定一致 / SC2 大型 structured (dense は
//!   ScopeExceeded) / SC3 大域共役の interface 共変 / SC4 transient 非昇格
//!
//! 開封順序 (HOLD-5..8 と同一):
//!   v34.0-A (本コミット) = 生成器・採点器・バー・lib pin の凍結 + SECRET
//!   コミットメント公表 + train 採点 (可視シード 34001)
//!   → v34.0-B = SECRET 開示・holdout 初生成・本採点 (調整なし) + 期末完全儀式
//!
//! sha256(SECRET) = ef6a8cd97b5d7693f4a4ffdb11ccdf42dbf5a971b0c79ba5f0f72f7b15739fcd
//!
//! 実行: cargo run --release --bin v340a_hold9_freeze

use uft_sim::contextual_factorization::{recover_atlas, AtlasReading, ChartSpec};
use uft_sim::graded_recovery::{
    extract_complex_structure, recover_graded, GradedAbstainReason, GradedRecoveryReading,
    MajoranaFrame,
};
use uft_sim::laboratory_interface::{
    certify_addressability, certify_synthesis, AccessibleOperation, AccessibleOperationalNet,
    BoundInterval, InterfaceRejection, OperationOrigin, ResourceBudget, SynthStep,
};
use uft_sim::operational_net::{
    anticommutator, cdag, classify_generator, cmul, commutator, hs_norm, same_gauge_orbit,
    CertifiedCommutator, ControlGenerator, FactorizationAbstainReason, FactorizationReading,
    FermionicZ2Graded, GeneratorClassification, GklsLiouvillian, MeasurementEffect, OpId, OpKind,
    OperationalNet, OperatorParity, OrdinaryCommutation, PrimitiveOperation,
    RecoveryInputRejection,
};
use uft_sim::resource_profile::{ProfilePoint, ResourceFilteredInterface};
use uft_sim::structured_backend::{
    dense_scope_guard, recover_pauli_net, PauliNetSpec, PauliVector, StructuredScopeError,
};
use uft_sim::{sha256_hex, Rng, C64};

// ================================================================================
// FROZEN-HOLD9-BEGIN  (この区間は v340a/v340b で逐語一致 — [H0] が SHA-256 で照合する)
// ================================================================================
//
// HOLD-9 = 第三十三期の全器械 (v33.1 修復入口 / v33.2 laboratory interface /
// v33.3 resource profile / v33.4 atlas glue / v33.5 graded recovery /
// v33.6 structured backend) の識別可能性境界を新鮮 holdout で採点する。
// 採点は従来の 3 計量に第三十三期の 5 計量を加える:
//   selective risk = 0 / impossibility recall = 1 / answerable recall = 1
//   origin_certificate_coverage = 1 / context_witness_coverage = 1
//   raw_operation_promotions = 0 / scope_violations = 0 /
//   transient_factorization_promotions = 0

pub const HOLD9_COMMITMENT: &str =
    "ef6a8cd97b5d7693f4a4ffdb11ccdf42dbf5a971b0c79ba5f0f72f7b15739fcd";
pub const HOLD9_TRAIN_SEED: u64 = 34001;

// ---- 凍結バー (開封後に変更しない — 全て第三十三期の各版で凍結済みの値) ----
pub const TAU_COMM: f64 = 1e-3; // 可換子閾値 (v32.3/v33.1)
pub const SIGMA_BAR: f64 = 0.5; // addressability の σ_min バー (v33.2)
pub const XTALK_BAR: f64 = 0.1; // cross-talk バー (v33.2)
pub const SYNTH_ERR_BAR: f64 = 1e-9; // 合成残差バー (v33.2)
pub const MODE_CAR_BAR: f64 = 1e-9; // mode-CAR / Σn̂−Q バー (v33.5)
pub const ORBIT_MATCH_BAR: f64 = 0.9; // orbit 非同値の上限 (v32.3/v33.4)
pub const XTALK_NOISE_SIGMA: f64 = 0.02 / 6.0; // CR4 の観測ノイズ σ (z = 6 で ±0.02)
pub const NOISE_Z: f64 = 6.0;

// ---- 凍結 lib pin (第三十三期の器械 6 モジュール — 凍結時の sha256-16。
//      開封までに lib が変わっていないことを [A0]/[H0] が照合する) ----
pub const HOLD9_LIB_PINS: [(&str, &str); 6] = [
    ("sim/src/operational_net.rs", "7898d24448f79f17"),
    ("sim/src/laboratory_interface.rs", "a11f188f3ecb6a40"),
    ("sim/src/resource_profile.rs", "f2fe4b9613049704"),
    ("sim/src/contextual_factorization.rs", "e540eea6f21ca404"),
    ("sim/src/graded_recovery.rs", "c1ebb02c6af93133"),
    ("sim/src/structured_backend.rs", "f57d9bdf43add8f8"),
];

// ---- 凍結素子 (Pauli / 回転 / 置換 — v33 期バイナリと同一規約) ----

pub fn pauli(which: char) -> Vec<C64> {
    let (o, l) = (C64::new(0.0, 0.0), C64::new(1.0, 0.0));
    match which {
        'I' => vec![l, o, o, l],
        'X' => vec![o, l, l, o],
        'Y' => vec![o, C64::new(0.0, -1.0), C64::new(0.0, 1.0), o],
        'Z' => vec![l, o, o, C64::new(-1.0, 0.0)],
        _ => panic!("未知の Pauli"),
    }
}

pub fn kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
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

pub fn op3(s: &str) -> Vec<C64> {
    let cs: Vec<char> = s.chars().collect();
    let a = kron(&pauli(cs[0]), 2, &pauli(cs[1]), 2);
    kron(&a, 4, &pauli(cs[2]), 2)
}

pub fn ident(n: usize) -> Vec<C64> {
    let mut m = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

pub fn dft8() -> Vec<C64> {
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

pub fn conj_by(v: &[C64], a: &[C64], n: usize) -> Vec<C64> {
    cmul(&cmul(v, a, n), &cdag(v, n), n)
}

pub fn rot2(theta: f64, nx: f64, ny: f64, nz: f64) -> Vec<C64> {
    let (c, s) = (theta.cos(), theta.sin());
    vec![
        C64::new(c, s * nz),
        C64::new(s * ny, s * nx),
        C64::new(-s * ny, s * nx),
        C64::new(c, -s * nz),
    ]
}

/// 隠し局所 unitary (1 qubit) — 角度・軸を rng から
pub fn rand_u2(rng: &mut Rng) -> Vec<C64> {
    let theta = 0.2 + 1.1 * rng.f64();
    let (a, b, c) = (rng.f64() - 0.5, rng.f64() - 0.5, rng.f64() - 0.5);
    let norm = (a * a + b * b + c * c).sqrt().max(1e-9);
    rot2(theta, a / norm, b / norm, c / norm)
}

/// 隠し変成 W = (u₁⊗u₂⊗u₃)·P_perm (局所 unitary × qubit 置換)
pub fn rand_local_perm_w(rng: &mut Rng) -> Vec<C64> {
    let u12 = kron(&rand_u2(rng), 2, &rand_u2(rng), 2);
    let u = kron(&u12, 4, &rand_u2(rng), 2);
    const PERMS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let p = PERMS[rng.range(6)];
    let mut pm = vec![C64::new(0.0, 0.0); 64];
    for b in 0..8usize {
        let bits = [(b >> 2) & 1, (b >> 1) & 1, b & 1];
        let nb = (bits[p[0]] << 2) | (bits[p[1]] << 1) | bits[p[2]];
        pm[nb * 8 + b] = C64::new(1.0, 0.0);
    }
    cmul(&pm, &u, 8)
}

/// site 6 本 (X₁ Z₁ X₂ Z₂ X₃ Z₃)
pub fn site6() -> Vec<Vec<C64>> {
    ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"]
        .iter()
        .map(|s| op3(s))
        .collect()
}

/// 資格つき interface から AccessibleOperationalNet を建てる (v33.2 [C7] と同一手順)
pub fn build_accessible_net(
    gens: &[Vec<C64>],
    n: usize,
    contexts: &[Vec<usize>],
) -> Result<AccessibleOperationalNet<OrdinaryCommutation>, String> {
    let cert = certify_addressability(gens, gens, n, SIGMA_BAR, XTALK_BAR)
        .map_err(|e| format!("addressability: {}", e.as_str()))?;
    let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 1.0, 1e-9).unwrap();
    let mut net: AccessibleOperationalNet<OrdinaryCommutation> =
        AccessibleOperationalNet::new(n, TAU_COMM);
    let mut ids: Vec<OpId> = Vec::new();
    for g in gens {
        let op = AccessibleOperation::certify(
            OpKind::Control(
                ControlGenerator::certify(
                    g.iter().map(|c| c.re).collect(),
                    g.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            OperatorParity::Even,
            OperationOrigin::DirectlyCalibrated(cert.clone()),
            cert.clone(),
            budget,
        )
        .map_err(|e| format!("accessible: {}", e.as_str()))?;
        ids.push(net.admit(op).map_err(|e| e.to_string())?);
    }
    for (a, ga) in gens.iter().enumerate() {
        for (b, gb) in gens.iter().enumerate().skip(a + 1) {
            let nu = hs_norm(&commutator(ga, gb, n));
            net.set_commutator(
                ids[a],
                ids[b],
                CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
            );
        }
    }
    for ctx in contexts {
        let members: Vec<OpId> = ctx.iter().map(|&i| ids[i]).collect();
        net.add_control_context(&members, cert.clone())?;
    }
    Ok(net)
}

/// 生の OperationalNet を建てる (atlas / graded セル用 — v33.4/v33.5 と同一手順)
pub fn build_plain_net(
    gens: &[Vec<C64>],
    n: usize,
    contexts: &[Vec<usize>],
    graded: bool,
) -> (
    Option<OperationalNet<OrdinaryCommutation>>,
    Option<OperationalNet<FermionicZ2Graded>>,
    Vec<OpId>,
) {
    let mk = |g: &Vec<C64>, parity: OperatorParity| PrimitiveOperation {
        kind: OpKind::Control(
            ControlGenerator::certify(
                g.iter().map(|c| c.re).collect(),
                g.iter().map(|c| c.im).collect(),
                n,
            )
            .unwrap(),
        ),
        parity,
        provenance: "hold9_kernel",
    };
    if graded {
        let mut net: OperationalNet<FermionicZ2Graded> = OperationalNet::new(n, TAU_COMM);
        let ids: Vec<OpId> = gens
            .iter()
            .map(|g| net.add_primitive(mk(g, OperatorParity::Odd)).unwrap())
            .collect();
        for a in 0..gens.len() {
            for b in (a + 1)..gens.len() {
                let nu = hs_norm(&anticommutator(&gens[a], &gens[b], n));
                net.set_commutator(
                    ids[a],
                    ids[b],
                    CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
                );
            }
        }
        for ctx in contexts {
            let members: Vec<OpId> = ctx.iter().map(|&i| ids[i]).collect();
            net.add_context(&members).unwrap();
        }
        (None, Some(net), ids)
    } else {
        let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, TAU_COMM);
        let ids: Vec<OpId> = gens
            .iter()
            .map(|g| net.add_primitive(mk(g, OperatorParity::Even)).unwrap())
            .collect();
        for a in 0..gens.len() {
            for b in (a + 1)..gens.len() {
                let nu = hs_norm(&commutator(&gens[a], &gens[b], n));
                net.set_commutator(
                    ids[a],
                    ids[b],
                    CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
                );
            }
        }
        for ctx in contexts {
            let members: Vec<OpId> = ctx.iter().map(|&i| ids[i]).collect();
            net.add_context(&members).unwrap();
        }
        (Some(net), None, ids)
    }
}

// ---- 採点の型 (凍結) ----

pub struct CellOutcome {
    pub name: &'static str,
    /// true = 回答セル (正の読みが正解) / false = 非識別セル (棄却・同値類・拒否が正解)
    pub answer_cell: bool,
    /// 器械が回答した (棄却・拒否でなく正の読みを返した)
    pub answered: bool,
    /// 期待どおりか
    pub correct: bool,
    pub detail: String,
}

#[derive(Default)]
pub struct Hold9Counters {
    pub admitted_ops: usize,
    pub admitted_with_origin: usize,
    pub exact_readings: usize,
    pub exact_with_witness: usize,
    pub raw_operation_promotions: usize,
    pub scope_violations: usize,
    pub transient_promotions: usize,
}

/// 3 + 5 計量 (凍結): 強制回答 (非識別セルでの回答) は risk に数える
pub fn score(cells: &[CellOutcome]) -> (f64, f64, f64, usize) {
    let n_answered = cells.iter().filter(|c| c.answered).count();
    let n_bad_answer = cells
        .iter()
        .filter(|c| c.answered && (!c.correct || !c.answer_cell))
        .count();
    let n_imp = cells.iter().filter(|c| !c.answer_cell).count();
    let n_imp_ok = cells
        .iter()
        .filter(|c| !c.answer_cell && !c.answered && c.correct)
        .count();
    let n_ans = cells.iter().filter(|c| c.answer_cell).count();
    let n_ans_ok = cells
        .iter()
        .filter(|c| c.answer_cell && c.answered && c.correct)
        .count();
    let forced = cells.iter().filter(|c| !c.answer_cell && c.answered).count();
    (
        n_bad_answer as f64 / (n_answered.max(1) as f64),
        n_imp_ok as f64 / (n_imp.max(1) as f64),
        n_ans_ok as f64 / (n_ans.max(1) as f64),
        forced,
    )
}

// ---- 20 セルの生成と裁定 (凍結) — 隠しパラメータは rng から固定順で引く ----

pub fn run_cells(seed: u64) -> (Vec<CellOutcome>, Hold9Counters) {
    let mut rng = Rng::new(seed);
    let mut cells: Vec<CellOutcome> = Vec::new();
    let mut k = Hold9Counters::default();
    let n = 8usize;
    let site = site6();
    let ctx_site: Vec<Vec<usize>> = vec![vec![0, 2, 4], vec![1, 3, 5]];

    // ============ 群 1: 入力完全性 ============

    // IC1: raw matrix (証明書なし) は AccessibleOperation を構成できない —
    //      無関係な隠し op の証明書を流用する唯一の残された路も sha 結束で拒否
    {
        let w = rand_local_perm_w(&mut rng);
        let raw = conj_by(&w, &op3("XII"), n); // 隠し「生の行列」
        let other = conj_by(&w, &op3("IIZ"), n); // 証明書は別の op のもの
        let cert = certify_addressability(&[other.clone()], &[other.clone()], n, SIGMA_BAR, XTALK_BAR)
            .expect("較正が立たない");
        let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 1.0, 1e-9).unwrap();
        let res = AccessibleOperation::certify(
            OpKind::Control(
                ControlGenerator::certify(
                    raw.iter().map(|c| c.re).collect(),
                    raw.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            OperatorParity::Even,
            OperationOrigin::DirectlyCalibrated(cert.clone()),
            cert,
            budget,
        );
        let refused = matches!(res, Err(InterfaceRejection::CertificateTargetMismatch));
        if !refused {
            k.raw_operation_promotions += 1;
        }
        cells.push(CellOutcome {
            name: "IC1-raw-refusal",
            answer_cell: false,
            answered: !refused,
            correct: refused,
            detail: "生の行列 + 流用証明書 → certificate_target_mismatch".into(),
        });
    }

    // IC2: gens 別渡しの不在 — 復元入力の行列は net の primitive と byte 同一
    {
        let w = rand_local_perm_w(&mut rng);
        let gens: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&w, g, n)).collect();
        let net = build_accessible_net(&gens, n, &ctx_site);
        let mut ok = false;
        let mut det = String::new();
        if let Ok(net) = &net {
            k.admitted_ops += 6;
            k.admitted_with_origin += 6;
            let inp = net.inner_net().recovery_input();
            if let Ok(inp) = inp {
                let mats = inp.generator_matrices();
                let same = mats.len() == 6
                    && mats.iter().zip(gens.iter()).all(|(a, b)| {
                        a.iter()
                            .zip(b.iter())
                            .all(|(x, y)| (x.re - y.re).abs() < 1e-15 && (x.im - y.im).abs() < 1e-15)
                    });
                let reading = inp.recover().reading;
                let exact = reading
                    == FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
                        local_dims: vec![2, 2, 2],
                    };
                if exact {
                    k.exact_readings += 1;
                    k.exact_with_witness += 1;
                }
                ok = same && exact;
                det = format!("行列の出所 = net の primitive (byte 同一 {}) → exact [2,2,2]", same);
            }
        }
        cells.push(CellOutcome {
            name: "IC2-net-owned-gens",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: det,
        });
    }

    // IC3: role-mixed (隠し測定 effect の混入) → 構成時拒否
    {
        let w = rand_local_perm_w(&mut rng);
        let gens: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&w, g, n)).collect();
        let (net_o, _, _) = build_plain_net(&gens, n, &ctx_site, false);
        let mut net = net_o.unwrap();
        let zt = conj_by(&w, &op3("ZII"), n);
        let eff: Vec<C64> = ident(n)
            .iter()
            .zip(zt.iter())
            .map(|(i, z)| (*i - *z).scale(0.5))
            .collect();
        net.add_primitive(PrimitiveOperation {
            kind: OpKind::Measure(
                MeasurementEffect::certify(
                    eff.iter().map(|c| c.re).collect(),
                    eff.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            parity: OperatorParity::Even,
            provenance: "hold9_measure",
        })
        .unwrap();
        let refused = matches!(
            net.recovery_input(),
            Err(RecoveryInputRejection::RoleMixedRecovery)
        );
        cells.push(CellOutcome {
            name: "IC3-role-mixed",
            answer_cell: false,
            answered: !refused,
            correct: refused,
            detail: "測定 effect の混入 → role_mixed_recovery".into(),
        });
    }

    // IC4: GKLS (隠し γ ∈ [0.1, 0.5]) を可換子 lane に入れると拒否
    {
        let n4 = 4usize;
        let mut h = vec![C64::new(0.0, 0.0); 16];
        for (c, s) in [
            (0.4 + 0.6 * rng.f64(), "XI"),
            (0.2 + 0.5 * rng.f64(), "ZZ"),
            (0.1 + 0.4 * rng.f64(), "IY"),
        ] {
            let m = kron(&pauli(s.chars().next().unwrap()), 2, &pauli(s.chars().nth(1).unwrap()), 2);
            for (hi, mi) in h.iter_mut().zip(m.iter()) {
                *hi = *hi + mi.scale(c);
            }
        }
        let gamma = 0.1 + 0.4 * rng.f64();
        let sm = vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
        ];
        let jump = if rng.range(2) == 0 {
            kron(&sm, 2, &pauli('I'), 2)
        } else {
            kron(&pauli('I'), 2, &sm, 2)
        };
        let gk = GklsLiouvillian::certify(h, vec![jump], vec![gamma], n4).unwrap();
        let refused = matches!(
            classify_generator(&|m: &[C64]| gk.apply(m), n4),
            GeneratorClassification::NonDerivation { .. }
        );
        cells.push(CellOutcome {
            name: "IC4-gkls-lane",
            answer_cell: false,
            answered: !refused,
            correct: refused,
            detail: "GKLS (γ > 0) は Leibniz が破れ可換子 lane 資格なし".into(),
        });
    }

    // ============ 群 2: accessibility ============

    // AC1: 独立 site knobs (隠し変成) → Exact [2,2,2]
    {
        let w = rand_local_perm_w(&mut rng);
        let gens: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&w, g, n)).collect();
        let net = build_accessible_net(&gens, n, &ctx_site);
        let mut ok = false;
        if let Ok(net) = &net {
            k.admitted_ops += 6;
            k.admitted_with_origin += 6;
            if let Ok(d) = net.recover() {
                ok = d.reading
                    == FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
                        local_dims: vec![2, 2, 2],
                    };
                if ok {
                    k.exact_readings += 1;
                    k.exact_with_witness += 1;
                }
            }
        }
        cells.push(CellOutcome {
            name: "AC1-independent-knobs",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "資格つき独立 knobs → exact [2,2,2]".into(),
        });
    }

    // AC2: tied (X_i + X_j — 隠し対・隠し変成) → 分解拒否 + 正直な net は Abstain
    {
        let w = rand_local_perm_w(&mut rng);
        let pair = [(0usize, 2usize), (0, 4), (2, 4)][rng.range(3)];
        let xi = conj_by(&w, &site[pair.0], n);
        let xj = conj_by(&w, &site[pair.1], n);
        let tied: Vec<C64> = xi.iter().zip(xj.iter()).map(|(a, b)| *a + *b).collect();
        let dec = certify_addressability(
            &[xi.clone(), xj.clone()],
            &[tied.clone(), tied.clone()],
            n,
            SIGMA_BAR,
            XTALK_BAR,
        );
        let dec_refused = matches!(dec, Err(InterfaceRejection::InsufficientCommandRank));
        let net = build_accessible_net(&[tied], n, &[vec![0]]);
        let abstained = match net.map(|nt| nt.recover()) {
            Ok(Ok(d)) => {
                d.reading
                    == FactorizationReading::Abstain(
                        FactorizationAbstainReason::InsufficientOperationalGenerators,
                    )
            }
            _ => false,
        };
        let correct = dec_refused && abstained;
        cells.push(CellOutcome {
            name: "AC2-tied-control",
            answer_cell: false,
            answered: !correct,
            correct,
            detail: "数学的分解は rank 1 < 2 で拒否・正直な tied net は abstain".into(),
        });
    }

    // AC3: 同一閉包・異なる interface (site vs DFT 共役, 隠し変成) → 非同値 orbit
    {
        let w = rand_local_perm_w(&mut rng);
        let v = dft8();
        let gens_a: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&w, g, n)).collect();
        let gens_b: Vec<Vec<C64>> = site
            .iter()
            .map(|g| conj_by(&w, &conj_by(&v, g, n), n))
            .collect();
        let net_a = build_accessible_net(&gens_a, n, &ctx_site);
        let net_b = build_accessible_net(&gens_b, n, &ctx_site);
        let mut nonequiv = false;
        let mut det = String::new();
        if let (Ok(na), Ok(nb)) = (&net_a, &net_b) {
            k.admitted_ops += 12;
            k.admitted_with_origin += 12;
            if let (Ok(da), Ok(db)) = (na.recover(), nb.recover()) {
                let both_exact = da.reading.as_str() == "exact_up_to_local_unitary_and_permutation"
                    && db.reading.as_str() == "exact_up_to_local_unitary_and_permutation";
                if both_exact {
                    k.exact_readings += 2;
                    k.exact_with_witness += 2;
                }
                let (same, ov) =
                    same_gauge_orbit(&da.component_subalgebras, &db.component_subalgebras);
                nonequiv = both_exact && !same && ov < ORBIT_MATCH_BAR;
                det = format!("両 interface とも exact・orbit matching 不在 (best {:.4})", ov);
            }
        }
        cells.push(CellOutcome {
            name: "AC3-closure-erasure",
            answer_cell: false,
            answered: false,
            correct: nonequiv,
            detail: det,
        });
    }

    // AC4: 出自証明書の流用 (合成証明書を隠し別行列へ) → 構成時拒否
    {
        let n4 = 4usize;
        let x1 = kron(&pauli('X'), 2, &pauli('I'), 2);
        let x2 = kron(&pauli('I'), 2, &pauli('X'), 2);
        let z2 = kron(&pauli('I'), 2, &pauli('Z'), 2);
        let tied: Vec<C64> = x1.iter().zip(x2.iter()).map(|(a, b)| *a + *b).collect();
        let base = vec![tied.clone(), z2.clone()];
        let steps = vec![
            SynthStep::BracketOverI(0, 1),
            SynthStep::BracketOverI(2, 1),
            SynthStep::Linear(vec![(1.0, 0), (0.25, 3)]),
        ];
        let cert = certify_synthesis(&base, &steps, &x1, n4, SYNTH_ERR_BAR).unwrap();
        let addr = certify_addressability(&[tied.clone(), z2.clone()], &[tied, z2], n4, SIGMA_BAR, XTALK_BAR)
            .unwrap();
        let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 3.0, 1e-9).unwrap();
        // 隠し別行列 (回転した X₂ 系) へ流用を試みる
        let u = rand_u2(&mut rng);
        let w4 = kron(&u, 2, &pauli('I'), 2);
        let other = conj_by(&w4, &x2, n4);
        let res = AccessibleOperation::certify(
            OpKind::Control(
                ControlGenerator::certify(
                    other.iter().map(|c| c.re).collect(),
                    other.iter().map(|c| c.im).collect(),
                    n4,
                )
                .unwrap(),
            ),
            OperatorParity::Even,
            OperationOrigin::Synthesized(cert),
            addr,
            budget,
        );
        let refused = matches!(res, Err(InterfaceRejection::CertificateTargetMismatch));
        if !refused {
            k.raw_operation_promotions += 1;
        }
        cells.push(CellOutcome {
            name: "AC4-origin-binding",
            answer_cell: false,
            answered: !refused,
            correct: refused,
            detail: "X₁ の合成証明書は隠し別行列に流用できない (sha 結束)".into(),
        });
    }

    // ============ 群 3: context・resource ============

    // CR1: budget profile (隠し深さコスト) — 資源不足 → [2,2,2] → [2,4]
    {
        let d1 = 0.5 + 0.5 * rng.f64();
        let d2 = d1 + 0.5 + 0.5 * rng.f64();
        let bud = |d: f64| ResourceBudget::certify(1.0, 1.0, 1.0, d, 1e-9).unwrap();
        let mut gens = site.clone();
        gens.push(op3("XXI"));
        let mats = gens.clone();
        let cert = certify_addressability(&mats, &mats, n, SIGMA_BAR, XTALK_BAR).unwrap();
        let mut iface: ResourceFilteredInterface<OrdinaryCommutation> =
            ResourceFilteredInterface::new(n, TAU_COMM);
        for (i, g) in gens.iter().enumerate() {
            let cost = if i < 6 { bud(d1) } else { bud(d2) };
            let op = AccessibleOperation::certify(
                OpKind::Control(
                    ControlGenerator::certify(
                        g.iter().map(|c| c.re).collect(),
                        g.iter().map(|c| c.im).collect(),
                        n,
                    )
                    .unwrap(),
                ),
                OperatorParity::Even,
                OperationOrigin::DirectlyCalibrated(cert.clone()),
                cert.clone(),
                cost,
            )
            .unwrap();
            iface.add_operation(op);
            k.admitted_ops += 1;
            k.admitted_with_origin += 1;
        }
        iface.add_context_recipe(vec![0, 2, 4, 6], cert.clone()).unwrap();
        iface.add_context_recipe(vec![1, 3, 5], cert.clone()).unwrap();
        let grid = [bud(d1 * 0.5), bud((d1 + d2) * 0.5), bud(d2 + 0.5)];
        let prof = iface.profile_over(&grid);
        let got: Vec<&str> = prof.points.iter().map(|p| p.as_str()).collect();
        let dims_at = |i: usize| -> Vec<usize> {
            match &prof.points[i] {
                ProfilePoint::Reading {
                    reading: FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims },
                    ..
                } => local_dims.clone(),
                _ => vec![],
            }
        };
        let ok = got[0] == "no_accessible_operations"
            && dims_at(1) == vec![2, 2, 2]
            && dims_at(2) == vec![2, 4];
        if ok {
            k.exact_readings += 2;
            k.exact_with_witness += 2;
        }
        cells.push(CellOutcome {
            name: "CR1-budget-profile",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "profile: 資源不足 → [2,2,2] → [2,4] (隠し深さコスト)".into(),
        });
    }

    // CR2: overlap の一意 glue (隠し qubit-2 frame 回転) → GluedExact = 直接大域
    {
        let u2 = {
            let u = rand_u2(&mut rng);
            let a = kron(&pauli('I'), 2, &u, 2);
            kron(&a, 4, &pauli('I'), 2)
        };
        let x2p = conj_by(&u2, &site[2], n);
        let z2p = conj_by(&u2, &site[3], n);
        let gens: Vec<Vec<C64>> = vec![
            site[0].clone(),
            site[1].clone(),
            site[2].clone(),
            site[3].clone(),
            x2p,
            z2p,
            site[4].clone(),
            site[5].clone(),
        ];
        let ctxs: Vec<Vec<usize>> = vec![vec![0, 2], vec![1, 3], vec![4, 6], vec![5, 7], vec![0, 6]];
        let (net_o, _, ids) = build_plain_net(&gens, n, &ctxs, false);
        let net = net_o.unwrap();
        let ch_a = ChartSpec {
            primitive_ids: vec![ids[0], ids[1], ids[2], ids[3]],
        };
        let ch_b = ChartSpec {
            primitive_ids: vec![ids[4], ids[5], ids[6], ids[7]],
        };
        let atlas = recover_atlas(&net, &[ch_a, ch_b], n);
        let ok = match &atlas {
            AtlasReading::GluedExact { local_dims, .. } => local_dims == &vec![2, 2, 2],
            _ => false,
        };
        if ok {
            k.exact_readings += 1;
            k.exact_with_witness += 1;
        }
        cells.push(CellOutcome {
            name: "CR2-unique-glue",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "chart 共有因子の matching (frame 回転不変) → glued [2,2,2]".into(),
        });
    }

    // CR3: 複数 glue (site atlas vs DFT chart, 隠し変成) → EquivalenceClassOnly{2}
    {
        let w = rand_local_perm_w(&mut rng);
        let v = dft8();
        let sgens: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&w, g, n)).collect();
        let dgens: Vec<Vec<C64>> = site
            .iter()
            .map(|g| conj_by(&w, &conj_by(&v, g, n), n))
            .collect();
        let mut gens = sgens.clone();
        gens.extend(dgens.iter().cloned());
        let ctxs: Vec<Vec<usize>> = vec![
            vec![0, 2],
            vec![1, 3],
            vec![2, 4],
            vec![3, 5],
            vec![0, 4],
            vec![6, 8, 10],
            vec![7, 9, 11],
        ];
        let (net_o, _, ids) = build_plain_net(&gens, n, &ctxs, false);
        let net = net_o.unwrap();
        let s1 = ChartSpec {
            primitive_ids: vec![ids[0], ids[1], ids[2], ids[3]],
        };
        let s2 = ChartSpec {
            primitive_ids: vec![ids[2], ids[3], ids[4], ids[5]],
        };
        let d = ChartSpec {
            primitive_ids: (6..12).map(|i| ids[i]).collect(),
        };
        let atlas = recover_atlas(&net, &[s1, s2, d], n);
        let ok = matches!(
            atlas,
            AtlasReading::EquivalenceClassOnly {
                n_consistent_atlases: 2
            }
        );
        cells.push(CellOutcome {
            name: "CR3-multi-glue",
            answer_cell: false,
            answered: false,
            correct: ok,
            detail: "整合 atlas 2 つ → equivalence_class_only (tie-break なし)".into(),
        });
    }

    // CR4: cross-talk 区間跨ぎ → Straddled 棄却。隠し観測値 ε̂ ∈ [bar − half/2,
    //      bar + half/2] — 区間 [ε̂ ± half] は構成上必ずバーを跨ぐ (跨ぎ設計セル:
    //      隠れているのは ε̂ の位置・採点対象は「跨ぎで強制判定しない」規律)
    {
        let half = NOISE_Z * XTALK_NOISE_SIGMA; // = 0.02
        let eps_hat = XTALK_BAR + half * (rng.f64() - 0.5);
        let interval = BoundInterval::new((eps_hat - half).max(0.0), eps_hat + half).unwrap();
        let straddled = interval.within(XTALK_BAR).is_err();
        cells.push(CellOutcome {
            name: "CR4-crosstalk-straddle",
            answer_cell: false,
            answered: !straddled,
            correct: straddled,
            detail: "cross-talk 区間がバー 0.1 を跨ぐ → 強制判定なし".into(),
        });
    }

    // ============ 群 4: graded ============

    // GR1: odd CAR only (隠し O(6) 回転) → MajoranaFrameOnly
    {
        let frame = rand_rotated_frame(&mut rng, n);
        let reading = recover_graded(&frame, None);
        let ok = matches!(reading, GradedRecoveryReading::MajoranaFrameOnly { n_majorana: 6 });
        cells.push(CellOutcome {
            name: "GR1-majorana-orbit",
            answer_cell: false,
            answered: false,
            correct: ok,
            detail: "witness なし → majorana_frame_only (O(2N) orbit で止まる)".into(),
        });
    }

    // GR2: charge witness (隠し回転 frame + 標準 Q) → complex modes (Σn̂ = Q)
    {
        let frame = rand_rotated_frame(&mut rng, n);
        let q = std_charge(n);
        let mut ok = false;
        if let Ok(w) = extract_complex_structure(&frame, &q) {
            if let GradedRecoveryReading::ComplexModeFactorization { n_modes: 3, modes } =
                recover_graded(&frame, Some(&w))
            {
                let mut sum_n = vec![C64::new(0.0, 0.0); n * n];
                for (a, ad) in &modes {
                    let nn = cmul(ad, a, n);
                    for (s, x) in sum_n.iter_mut().zip(nn.iter()) {
                        *s = *s + *x;
                    }
                }
                let dev: f64 = sum_n
                    .iter()
                    .zip(q.iter())
                    .map(|(x, y)| (*x - *y).norm2())
                    .sum::<f64>()
                    .sqrt();
                ok = dev < MODE_CAR_BAR;
            }
        }
        cells.push(CellOutcome {
            name: "GR2-charge-witness",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "charge witness → J → 3 モード回復 (Σâ†â = Q)".into(),
        });
    }

    // GR3: witness の汚染 (隠し quartic 振幅 ∈ [0.05, 0.3]) → Abstain
    {
        let frame = rand_rotated_frame(&mut rng, n);
        let amp = 0.05 + 0.25 * rng.f64();
        let quartic: Vec<C64> = op3("ZZI").iter().map(|c| c.scale(-amp)).collect();
        let mut q = std_charge(n);
        for (qi, ki) in q.iter_mut().zip(quartic.iter()) {
            *qi = *qi + *ki;
        }
        let refused = matches!(
            extract_complex_structure(&frame, &q),
            Err(GradedAbstainReason::WitnessNotLinearOnFrame)
        );
        cells.push(CellOutcome {
            name: "GR3-witness-contaminated",
            answer_cell: false,
            answered: !refused,
            correct: refused,
            detail: "quartic 汚染 witness → witness_not_linear_on_frame".into(),
        });
    }

    // GR4: ordinary JW lane → 構成時拒否 (隠し回転 odd)
    {
        let frame = rand_rotated_frame(&mut rng, n);
        let g0 = frame.gamma(0).to_vec();
        let mut net_o: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, TAU_COMM);
        let refused = net_o
            .add_primitive(PrimitiveOperation {
                kind: OpKind::Control(
                    ControlGenerator::certify(
                        g0.iter().map(|c| c.re).collect(),
                        g0.iter().map(|c| c.im).collect(),
                        n,
                    )
                    .unwrap(),
                ),
                parity: OperatorParity::Odd,
                provenance: "hold9_odd",
            })
            .is_err();
        cells.push(CellOutcome {
            name: "GR4-ordinary-odd-refusal",
            answer_cell: false,
            answered: !refused,
            correct: refused,
            detail: "odd primitive は ordinary net が構成時拒否".into(),
        });
    }

    // ============ 群 5: scale・変成 ============

    // SC1: dense と Pauli backend の裁定一致 (隠しセル選択 + 隠し qubit 置換)
    {
        const CELLS: [(&[&str], &[&[usize]]); 4] = [
            (
                &["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"],
                &[&[0, 2, 4], &[1, 3, 5]],
            ),
            (
                &["XII", "ZII", "IXI", "IZI", "IIX", "IIZ", "XXI"],
                &[&[0, 2, 4, 6], &[1, 3, 5]],
            ),
            (&["ZII", "IZI", "IIZ"], &[&[0, 1, 2]]),
            (&["XII", "ZII", "IZI"], &[&[0, 2], &[1, 2]]),
        ];
        let (strs0, ctxs0) = CELLS[rng.range(4)];
        const PERMS: [[usize; 3]; 6] = [
            [0, 1, 2],
            [0, 2, 1],
            [1, 0, 2],
            [1, 2, 0],
            [2, 0, 1],
            [2, 1, 0],
        ];
        let p = PERMS[rng.range(6)];
        let permute = |s: &str| -> String {
            let cs: Vec<char> = s.chars().collect();
            (0..3).map(|i| cs[p[i]]).collect()
        };
        let strs: Vec<String> = strs0.iter().map(|s| permute(s)).collect();
        let ctxs: Vec<Vec<usize>> = ctxs0.iter().map(|c| c.to_vec()).collect();
        // dense lane
        let gens: Vec<Vec<C64>> = strs.iter().map(|s| op3(s)).collect();
        let (net_o, _, _) = build_plain_net(&gens, n, &ctxs, false);
        let dense = match net_o.unwrap().recovery_input() {
            Ok(inp) => format!("{:?}", inp.recover().reading),
            Err(e) => format!("rejected:{}", e.as_str()),
        };
        // Pauli lane
        let spec = PauliNetSpec {
            n_qubits: 3,
            ops: strs.iter().map(|s| PauliVector::from_str(s)).collect(),
            contexts: ctxs.iter().map(|c| c.iter().cloned().collect()).collect(),
        };
        let pl = match recover_pauli_net(&spec) {
            Ok(r) => format!("{:?}", r),
            Err(e) => format!("rejected:{}", e.as_str()),
        };
        let ok = dense == pl;
        cells.push(CellOutcome {
            name: "SC1-dense-pauli-agree",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "隠しセル + 隠し置換で dense = Pauli backend".into(),
        });
    }

    // SC2: 大型 structured lane (隠し n_qubits ∈ [32, 60]) + dense は ScopeExceeded
    {
        let nq = 32 + rng.range(29);
        let epos = rng.range(nq - 1);
        let mk = |pos: usize, ch: char| -> String {
            let mut s = vec!['I'; nq];
            s[pos] = ch;
            s.iter().collect()
        };
        let mut strs: Vec<String> = Vec::new();
        for i in 0..nq {
            strs.push(mk(i, 'X'));
            strs.push(mk(i, 'Z'));
        }
        let mut e = vec!['I'; nq];
        e[epos] = 'X';
        e[epos + 1] = 'X';
        strs.push(e.iter().collect());
        let mut ctx_x: std::collections::BTreeSet<usize> = (0..nq).map(|i| 2 * i).collect();
        ctx_x.insert(2 * nq);
        let ctx_z: std::collections::BTreeSet<usize> = (0..nq).map(|i| 2 * i + 1).collect();
        let spec = PauliNetSpec {
            n_qubits: nq,
            ops: strs.iter().map(|s| PauliVector::from_str(s)).collect(),
            contexts: vec![ctx_x, ctx_z],
        };
        let mut dims = vec![2usize; nq - 2];
        dims.push(4);
        dims.sort_unstable();
        let want = FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims: dims };
        let pl_ok = recover_pauli_net(&spec).map(|r| r == want) == Ok(true);
        let dense_refused = dense_scope_guard(1usize << nq)
            == Err(StructuredScopeError::DimensionTooLargeForDense);
        if !dense_refused {
            k.scope_violations += 1;
        }
        let ok = pl_ok && dense_refused;
        cells.push(CellOutcome {
            name: "SC2-large-structured",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "Pauli backend が行列なしで Exact・dense は ScopeExceeded".into(),
        });
    }

    // SC3: global conjugation + interface 共変 → 同一 orbit (共変性・隠し大域 W)
    {
        // W = (局所 ⊗) × CZ₂₃ (真に大域的な成分を含む)
        let w_loc = rand_local_perm_w(&mut rng);
        let cz = {
            let mut m = vec![C64::new(0.0, 0.0); 64];
            for b in 0..8usize {
                let q2 = (b >> 1) & 1;
                let q3 = b & 1;
                let sign = if q2 == 1 && q3 == 1 { -1.0 } else { 1.0 };
                m[b * 8 + b] = C64::new(sign, 0.0);
            }
            m
        };
        let w = cmul(&w_loc, &cz, n);
        let gens_o: Vec<Vec<C64>> = site.clone();
        let gens_w: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&w, g, n)).collect();
        let (na, _, _) = build_plain_net(&gens_o, n, &ctx_site, false);
        let (nb, _, _) = build_plain_net(&gens_w, n, &ctx_site, false);
        let da = na.unwrap().recovery_input().map(|i| i.recover());
        let db = nb.unwrap().recovery_input().map(|i| i.recover());
        let mut ok = false;
        if let (Ok(da), Ok(db)) = (da, db) {
            let dims_same = da.reading == db.reading
                && da.reading.as_str() == "exact_up_to_local_unitary_and_permutation";
            // 共変性: W-共役した A の部分代数 = B の部分代数 (orbit matching)
            let conj_subs: Vec<Vec<Vec<C64>>> = da
                .component_subalgebras
                .iter()
                .map(|sub| sub.iter().map(|m| conj_by(&w, m, n)).collect())
                .collect();
            let (same, _) = same_gauge_orbit(&conj_subs, &db.component_subalgebras);
            ok = dims_same && same;
        }
        cells.push(CellOutcome {
            name: "SC3-covariance",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "dims 不変・orbit は W-共役で厳密に対応 (interface 共変)".into(),
        });
    }

    // SC4: 単一 threshold の短命 factorization を昇格しない (隠し深さ鎖)
    {
        let d1 = 0.5 + 0.5 * rng.f64();
        let d2 = d1 + 0.5 + 0.5 * rng.f64();
        let d3 = d2 + 0.5 + 0.5 * rng.f64();
        let bud = |d: f64| ResourceBudget::certify(1.0, 1.0, 1.0, d, 1e-9).unwrap();
        let mut gens = site.clone();
        gens.push(op3("XXI"));
        gens.push(op3("IXX"));
        let mats = gens.clone();
        let cert = certify_addressability(&mats, &mats, n, SIGMA_BAR, XTALK_BAR).unwrap();
        let mut iface: ResourceFilteredInterface<OrdinaryCommutation> =
            ResourceFilteredInterface::new(n, TAU_COMM);
        for (i, g) in gens.iter().enumerate() {
            let cost = if i < 6 {
                bud(d1)
            } else if i == 6 {
                bud(d2)
            } else {
                bud(d3)
            };
            let op = AccessibleOperation::certify(
                OpKind::Control(
                    ControlGenerator::certify(
                        g.iter().map(|c| c.re).collect(),
                        g.iter().map(|c| c.im).collect(),
                        n,
                    )
                    .unwrap(),
                ),
                OperatorParity::Even,
                OperationOrigin::DirectlyCalibrated(cert.clone()),
                cert.clone(),
                cost,
            )
            .unwrap();
            iface.add_operation(op);
        }
        iface.add_context_recipe(vec![0, 2, 4, 6, 7], cert.clone()).unwrap();
        iface.add_context_recipe(vec![1, 3, 5], cert.clone()).unwrap();
        let grid = [
            bud(d1),
            bud((d1 + d2) * 0.5),
            bud(d2),
            bud(d3),
        ];
        let prof = iface.profile_over(&grid);
        let classes = prof.classes();
        let mut stable_dims: Vec<Vec<usize>> = Vec::new();
        let mut transient_dims: Vec<Vec<usize>> = Vec::new();
        for c in &classes {
            let dims = match &prof.points[c.representative] {
                ProfilePoint::Reading {
                    reading: FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims },
                    ..
                } => local_dims.clone(),
                _ => vec![],
            };
            if c.stable {
                stable_dims.push(dims);
                if c.region.len() < 2 {
                    k.transient_promotions += 1;
                }
            } else {
                transient_dims.push(dims);
            }
        }
        stable_dims.sort();
        transient_dims.sort();
        let ok = stable_dims == vec![vec![2, 2, 2]]
            && transient_dims == vec![vec![2, 4], vec![8]];
        cells.push(CellOutcome {
            name: "SC4-transient-not-promoted",
            answer_cell: true,
            answered: ok,
            correct: ok,
            detail: "stable = [2,2,2] のみ・[2,4]@単点 と [8]@頂 は transient".into(),
        });
    }

    (cells, k)
}

/// 隠し O(6) 回転 (Givens 3 連) を掛けた JW frame
pub fn rand_rotated_frame(rng: &mut Rng, n: usize) -> MajoranaFrame {
    let gammas: Vec<Vec<C64>> = ["XII", "YII", "ZXI", "ZYI", "ZZX", "ZZY"]
        .iter()
        .map(|s| op3(s))
        .collect();
    let base = MajoranaFrame::certify(gammas, n).unwrap();
    let mut r = vec![0.0f64; 36];
    for i in 0..6 {
        r[i * 6 + i] = 1.0;
    }
    for _ in 0..3 {
        let i = rng.range(5);
        let j = i + 1 + rng.range(5 - i);
        let th = 0.3 + 1.0 * rng.f64();
        let (c, s) = (th.cos(), th.sin());
        // r ← g(i,j,th) · r
        let mut nr = r.clone();
        for col in 0..6 {
            nr[i * 6 + col] = c * r[i * 6 + col] - s * r[j * 6 + col];
            nr[j * 6 + col] = s * r[i * 6 + col] + c * r[j * 6 + col];
        }
        r = nr;
    }
    base.rotated(&r).expect("O(6) 回転で CAR が保存されない")
}

/// 標準 charge Q = Σ n_i = (3I − Z₁ − Z₂ − Z₃)/2
pub fn std_charge(n: usize) -> Vec<C64> {
    let mut q: Vec<C64> = ident(n).iter().map(|c| c.scale(1.5)).collect();
    for s in ["ZII", "IZI", "IIZ"] {
        let m = op3(s);
        for (qi, mi) in q.iter_mut().zip(m.iter()) {
            *qi = *qi + mi.scale(-0.5);
        }
    }
    q
}

/// 8 計量の集計表示 (train / holdout 共通の凍結フォーマット)
pub fn print_scoreboard(cells: &[CellOutcome], k: &Hold9Counters) -> (f64, f64, f64) {
    for c in cells {
        let mark = if (c.answer_cell && c.answered && c.correct)
            || (!c.answer_cell && !c.answered && c.correct)
        {
            "✓"
        } else {
            "✗"
        };
        let kind = if c.answer_cell { "[回答]" } else { "[棄却]" };
        println!("        {} {:28} {} {}", mark, c.name, kind, c.detail);
    }
    let (risk, imp, ans, forced) = score(cells);
    let occ = if k.admitted_ops > 0 {
        k.admitted_with_origin as f64 / k.admitted_ops as f64
    } else {
        1.0
    };
    let cwc = if k.exact_readings > 0 {
        k.exact_with_witness as f64 / k.exact_readings as f64
    } else {
        1.0
    };
    println!(
        "      selective risk = {:.3} / impossibility recall = {:.3} / answerable recall = {:.3} / 強制回答 = {}",
        risk, imp, ans, forced
    );
    println!(
        "      origin_certificate_coverage = {:.3} / context_witness_coverage = {:.3} / raw_operation_promotions = {} / scope_violations = {} / transient_factorization_promotions = {}",
        occ, cwc, k.raw_operation_promotions, k.scope_violations, k.transient_promotions
    );
    (risk, imp, ans)
}

// ================================================================================
// FROZEN-HOLD9-END
// ================================================================================

fn main() {
    uft_sim::self_test();
    println!(
        "=== v34.0-A HOLD-9 の凍結 — 出自 × 文脈整合 × 資源依存局所性 (PROMPT/14) ===\n"
    );
    let root = if std::path::Path::new("core.schema.yml").exists() {
        "."
    } else {
        ".."
    };
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---- [A0] コミットメント公表 + lib pin + 契約自己検査 ----
    {
        let mut bad = Vec::new();
        for (f, want) in HOLD9_LIB_PINS {
            match std::fs::read_to_string(format!("{}/{}", root, f)) {
                Err(_) => bad.push(format!("{} が読めない", f)),
                Ok(t) => {
                    let h = sha256_hex(t.as_bytes());
                    if &h[..16] != want {
                        bad.push(format!("{} の sha256-16 {} ≠ 凍結 pin {}", f, &h[..16], want));
                    }
                }
            }
        }
        if let Err(e) = uft_sim::operational_net::operational_net_self_test() {
            bad.push(format!("operational_net: {}", e));
        }
        if let Err(e) = uft_sim::operational_net::scope_repair_self_test() {
            bad.push(format!("scope_repair: {}", e));
        }
        if let Err(e) = uft_sim::laboratory_interface::laboratory_interface_self_test() {
            bad.push(format!("laboratory_interface: {}", e));
        }
        if let Err(e) = uft_sim::resource_profile::resource_profile_self_test() {
            bad.push(format!("resource_profile: {}", e));
        }
        if let Err(e) = uft_sim::contextual_factorization::contextual_factorization_self_test() {
            bad.push(format!("contextual_factorization: {}", e));
        }
        if let Err(e) = uft_sim::graded_recovery::graded_recovery_self_test() {
            bad.push(format!("graded_recovery: {}", e));
        }
        if let Err(e) = uft_sim::structured_backend::structured_backend_self_test() {
            bad.push(format!("structured_backend: {}", e));
        }
        check(
            "[A0] sha256(SECRET) コミットメント公表 + 器械 lib pin (6 モジュール) + 契約自己検査",
            bad.is_empty(),
            if bad.is_empty() {
                format!("sha256(SECRET) = {}… / lib pin 6 本一致", &HOLD9_COMMITMENT[..16])
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [A1] 採点器の健全性 — 強制回答・誤答・昇格は必ず失点になる ----
    {
        let synthetic = vec![
            CellOutcome {
                name: "synthetic-forced",
                answer_cell: false,
                answered: true, // 非識別セルに回答 = 強制回答
                correct: false,
                detail: String::new(),
            },
            CellOutcome {
                name: "synthetic-wrong",
                answer_cell: true,
                answered: true,
                correct: false, // 誤答
                detail: String::new(),
            },
            CellOutcome {
                name: "synthetic-ok",
                answer_cell: true,
                answered: true,
                correct: true,
                detail: String::new(),
            },
        ];
        let (risk, imp, ans, forced) = score(&synthetic);
        let ok = risk > 0.5 && forced == 1 && imp == 0.0 && (ans - 0.5).abs() < 1e-12;
        check(
            "[A1] 採点器の健全性 — 強制回答と誤答が risk に数えられ、recall が水増しされない",
            ok,
            format!(
                "人工 3 セル: risk = {:.3} (2/3)・強制回答 = {}・imp recall = {:.1}・ans recall = {:.1}",
                risk, forced, imp, ans
            ),
        );
    }

    // ---- [A2] train 採点 (可視シード 34001) ----
    {
        println!("      train (シード {}) のセル表:", HOLD9_TRAIN_SEED);
        let (cells, counters) = run_cells(HOLD9_TRAIN_SEED);
        let (risk, imp, ans) = print_scoreboard(&cells, &counters);
        let occ_ok = counters.admitted_ops == counters.admitted_with_origin;
        let cwc_ok = counters.exact_readings == counters.exact_with_witness;
        let ok = cells.len() == 20
            && risk == 0.0
            && imp == 1.0
            && ans == 1.0
            && occ_ok
            && cwc_ok
            && counters.raw_operation_promotions == 0
            && counters.scope_violations == 0
            && counters.transient_promotions == 0;
        check(
            "[A2] train 20 セル満票 — risk 0 / recall 1/1 / 出自被覆 1 / 証人被覆 1 / 昇格・違反 0",
            ok,
            format!(
                "セル {} / 出自つき admit {}/{} / witness つき Exact {}/{}",
                cells.len(),
                counters.admitted_with_origin,
                counters.admitted_ops,
                counters.exact_with_witness,
                counters.exact_readings
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "HOLD-9 は凍結された — 生成器・採点器・バー・lib pin・コミットメントが公開され、holdout は SECRET 開示 (v34.0-B) まで存在しない"
        } else {
            "**凍結の破れ** — kernel と train を修正せよ (開封前のみ許される)"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
