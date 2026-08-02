//! v33.1 境界監査と型スコープ修復 — contexts は定理入力・accessible primitive は未証明 (PROMPT/14)
//!
//! 第三十三期の開幕。期テーゼ:
//!
//! > **可アクセス性は作用素単体の属性ではなく、系・制御器・測定器・資源・誤差
//! > 証明書の関係である。局所性は、証明付き laboratory interface が生成する
//! > role-typed context atlas 上で整合する因子分解の、資源スケールにわたる安定な
//! > 同値類としてのみ識別される。**
//!
//! その第一歩は新機能ではなく**境界監査**である。v32.3 の復元器 (v323 — 凍結原本)
//! は (net, gens) を並行に受け、`contexts` を一切参照しなかった。これは定理の誤り
//! ではない — 定理の仮定「各ノードの局所生成子だけが選別された primitive family」
//! が型に存在しなかった、という型スコープの空隙である。本版はこの空隙を反例・型・
//! 禁止変換で固定する (黙って直さない — 旧署名は v323 に凍結保存し、修復は
//! operational_net.rs の新入口として並置):
//!
//!   [B0] 契約自己検査 — scope_repair_self_test (新契約) + operational_net_self_test
//!        (v32.2 契約の不変) + qrn_core/readout_contract の封鎖不変
//!   [B1] **contexts 盲目性の機械実証**: 同一 primitive・同一証明書で contexts だけ
//!        が異なる (∅ vs atlas) 2 つの net に対し、v32.3 参照手順は成分・読みとも
//!        完全一致 (= contexts は入力でなかった)。修復入口 recovery_input は ∅ を
//!        構成時拒否 (NoDeclaredContexts) — contexts が load-bearing になった
//!   [B2] **entangler 負制御 (primitive 選別の循環)**: site 6 primitive → Exact
//!        [2,2,2]。independently accessible な entangler X₁X₂ を 1 本加えるだけで
//!        [2,4] に併合 — 大域閉包は同一 (M₈) のまま、**選別が答えを入力している**。
//!        この循環は型修復では解けない (v33.2 Certified Laboratory Interface の動機)
//!   [B3] **修復入口の資格ゲート**: role-mixed (測定混入) / 文脈 0 / 被覆不完全は
//!        構成時拒否 (Abstain でなく型エラー)。資格を満たす net では旧解を再現
//!        (site → Exact [2,2,2]・qutrit×qubit → Exact [2,3])
//!   [B4] **代数的可換 ↛ 操作的両立 (禁止変換 12)**: singleton 文脈のみの net は
//!        全対に Commuting 証明書があっても成分間の共同 addressability の証人が
//!        なく Abstain(OperationalCompatibilityUnwitnessed) — 参照手順 (旧) は同じ
//!        素材を Exact [2,2,2] と読む (修復が変えた点の対照)。JointContextWitness
//!        の唯一の構成は宣言済み文脈の共有
//!   [B5] **Liouvillian lane の型分離 (禁止変換 13)**: 導分 (Leibniz) 証明書 —
//!        L = −i[H,·] は Leibniz/†-共変/unital 全成立 + Ĥ 復元 (rel ≤ 1e-10)・
//!        GKLS (γ > 0) は Leibniz が γ 比例で破れ NonDerivation・v32.4 の R⁽¹⁾
//!        公式は GKLS 測定と不一致 (γ = 0.3 で 0.6・γ = 0 で ≤ 1e-6 に回復) —
//!        HamiltonianCommutatorLiouvillian の証明書は GklsLiouvillian に昇格しない
//!   [B6] 封鎖の source/schema/文書検査 — 禁止 impl From の不在・禁止変換 12/13 の
//!        登録・新概念の登録・v323 旧署名の凍結保存・uft-v33.1.md アンカー
//!
//! 実行: cargo run --release --bin v331_scope_repair

use std::fs;
use std::path::Path;
use uft_sim::operational_net::*;
use uft_sim::C64;

// ---------------------------------------------------------------- Pauli / kron 素子 (v322/v323 と同一)

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

fn op2(s: &str) -> Vec<C64> {
    let cs: Vec<char> = s.chars().collect();
    kron(&pauli(cs[0]), 2, &pauli(cs[1]), 2)
}

fn ident(n: usize) -> Vec<C64> {
    let mut m = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

// ---------------------------------------------------------------- net 構築 (exact ノルム証明書)

fn build_net(
    gens: &[Vec<C64>],
    n: usize,
    tau: f64,
) -> (OperationalNet<OrdinaryCommutation>, Vec<OpId>) {
    let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
    let mut ids = Vec::new();
    for g in gens {
        let p = PrimitiveOperation {
            kind: OpKind::Control(
                ControlGenerator::certify(
                    g.iter().map(|c| c.re).collect(),
                    g.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            parity: OperatorParity::Even,
            provenance: "v331_control",
        };
        ids.push(net.add_primitive(p).unwrap());
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
    (net, ids)
}

// ---------------------------------------------------------------- v32.3 参照手順 (旧署名の監査対象)

/// v32.3 の凍結決定手順の参照再現 — **旧署名 (net と別渡しの gens・contexts 非参照)
/// をそのまま持つ**。歴史的原本は v323_factorization_recovery.rs (凍結)。ここでの
/// 役割は「contexts が入力でなかったこと」を機械実証する監査対象であり、新規復元は
/// operational_net::MarkedRecoveryInput::recover (v33.1 の唯一の型付き入口) を使う。
fn reference_recover_v323(
    net: &OperationalNet<OrdinaryCommutation>,
    gens: &[Vec<C64>],
    n: usize,
) -> FactorizationReading {
    let abstain = FactorizationReading::Abstain;
    let comps = match net.noncommutation_components() {
        Ok(c) => c,
        Err(r) => return abstain(r),
    };
    let joint = algebra_closure(gens, n);
    let mut commutative = true;
    'outer: for a in gens {
        for b in gens {
            if hs_norm(&commutator(a, b, n)) > 1e-9 {
                commutative = false;
                break 'outer;
            }
        }
    }
    if commutative {
        return abstain(FactorizationAbstainReason::InsufficientOperationalGenerators);
    }
    let center = closure_center_basis(&joint, gens, n);
    if center.is_empty() {
        return abstain(FactorizationAbstainReason::ComponentNotFactor);
    }
    if center.len() == 1 {
        if joint.len() != n * n {
            return abstain(FactorizationAbstainReason::InsufficientOperationalGenerators);
        }
        let mut dims = Vec::new();
        for comp in &comps {
            let sub: Vec<Vec<C64>> = comp.iter().map(|&i| gens[i as usize].clone()).collect();
            let cl = algebra_closure(&sub, n);
            let d2 = cl.len();
            let d = (d2 as f64).sqrt().round() as usize;
            if d * d != d2 || d < 2 {
                return abstain(FactorizationAbstainReason::ComponentNotFactor);
            }
            if closure_center_basis(&cl, &sub, n).len() != 1 {
                return abstain(FactorizationAbstainReason::ComponentNotFactor);
            }
            dims.push(d);
        }
        if dims.iter().product::<usize>() != n {
            return abstain(FactorizationAbstainReason::ComponentNotFactor);
        }
        dims.sort_unstable();
        return FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims: dims };
    }
    let projs = match closure_central_projectors(&center, n) {
        Some(p) => p,
        None => return abstain(FactorizationAbstainReason::ComponentNotFactor),
    };
    let mut sectors = Vec::new();
    for p in &projs {
        let tr: f64 = (0..n).map(|i| p[i * n + i].re).sum();
        let b_dim = tr.round() as usize;
        if b_dim == 0 || (tr - b_dim as f64).abs() > 1e-7 {
            return abstain(FactorizationAbstainReason::ComponentNotFactor);
        }
        let mut restricted: Vec<Vec<C64>> = Vec::new();
        for b in &joint {
            let pbp = cmul(p, &cmul(b, p, n), n);
            if hs_norm(&pbp) < 1e-9 {
                continue;
            }
            push_ortho(&mut restricted, &pbp, 1e-8);
        }
        let m2 = restricted.len();
        let m = (m2 as f64).sqrt().round() as usize;
        if m * m != m2 || b_dim % m != 0 {
            return abstain(FactorizationAbstainReason::ComponentNotFactor);
        }
        sectors.push((m, b_dim / m));
    }
    sectors.sort_unstable();
    FactorizationReading::SuperselectionSectors { sectors }
}

// ---------------------------------------------------------------- GKLS 応答測定 (RK4 + 4 次 stencil)

/// dρ/dt = L(ρ) の RK4 固定刻み積分 (決定的)。t < 0 も ODE として同じに扱う。
fn rk4_evolve(l: &GklsLiouvillian, rho0: &[C64], t: f64, steps: usize) -> Vec<C64> {
    let n2 = rho0.len();
    let dt = t / steps as f64;
    let mut rho = rho0.to_vec();
    for _ in 0..steps {
        let k1 = l.apply(&rho);
        let mut tmp: Vec<C64> = (0..n2).map(|i| rho[i] + k1[i].scale(dt / 2.0)).collect();
        let k2 = l.apply(&tmp);
        tmp = (0..n2).map(|i| rho[i] + k2[i].scale(dt / 2.0)).collect();
        let k3 = l.apply(&tmp);
        tmp = (0..n2).map(|i| rho[i] + k3[i].scale(dt)).collect();
        let k4 = l.apply(&tmp);
        for i in 0..n2 {
            rho[i] = rho[i]
                + (k1[i] + k2[i].scale(2.0) + k3[i].scale(2.0) + k4[i]).scale(dt / 6.0);
        }
    }
    rho
}

fn tr_prod_re(b: &[C64], m: &[C64], n: usize) -> f64 {
    let mut s = C64::new(0.0, 0.0);
    for i in 0..n {
        for k in 0..n {
            s = s + b[i * n + k] * m[k * n + i];
        }
    }
    s.re
}

/// 測定 lane: GKLS 発展下の (ḃ⁺ − ḃ⁻)/(2ε) — 4 次 stencil
fn measure_r1_gkls(
    l: &GklsLiouvillian,
    rho0: &[C64],
    a: &[C64],
    b: &[C64],
    eps: f64,
    n: usize,
) -> f64 {
    let h = 0.005;
    let steps = 400usize;
    let mut d1 = [0.0f64; 2];
    for (s, sign) in [(0usize, 1.0f64), (1usize, -1.0f64)] {
        let rho: Vec<C64> = rho0
            .iter()
            .zip(a.iter())
            .map(|(r, x)| *r + x.scale(sign * eps))
            .collect();
        let f = |t: f64| -> f64 {
            if t == 0.0 {
                tr_prod_re(b, &rho, n)
            } else {
                tr_prod_re(b, &rk4_evolve(l, &rho, t, steps), n)
            }
        };
        let (fm2, fm1, f1, f2) = (f(-2.0 * h), f(-h), f(h), f(2.0 * h));
        d1[s] = (fm2 - 8.0 * fm1 + 8.0 * f1 - f2) / (12.0 * h);
    }
    (d1[0] - d1[1]) / (2.0 * eps)
}

fn main() {
    uft_sim::self_test();
    println!("=== v33.1 境界監査と型スコープ修復 — 第三十三期 開幕 (PROMPT/14) ===\n");
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
    let tau = 1e-3; // 可換子閾値 (v32.3 と同一の凍結値)

    let site_gens: Vec<Vec<C64>> = (0..3)
        .flat_map(|i| {
            let mut s = ['I', 'I', 'I'];
            s[i] = 'X';
            let x = op3(&s.iter().collect::<String>());
            s[i] = 'Z';
            let z = op3(&s.iter().collect::<String>());
            vec![x, z]
        })
        .collect();

    // site net の文脈 atlas (可換分解の 2 区画): {X₁,X₂,X₃} と {Z₁,Z₂,Z₃}
    let add_site_atlas = |net: &mut OperationalNet<OrdinaryCommutation>, ids: &[OpId]| {
        net.add_context(&[ids[0], ids[2], ids[4]]).unwrap();
        net.add_context(&[ids[1], ids[3], ids[5]]).unwrap();
    };

    // ---- [B0] 契約自己検査 ----
    {
        let mut bad = Vec::new();
        if let Err(e) = operational_net_self_test() {
            bad.push(format!("operational_net_self_test (v32.2 契約): {}", e));
        }
        if let Err(e) = scope_repair_self_test() {
            bad.push(format!("scope_repair_self_test (v33.1 契約): {}", e));
        }
        if let Err(e) = uft_sim::qrn_core::qrn_core_self_test() {
            bad.push(format!("qrn_core_self_test: {}", e));
        }
        if let Err(e) = uft_sim::readout_contract::readout_contract_self_test() {
            bad.push(format!("readout_contract_self_test: {}", e));
        }
        check(
            "[B0] 契約自己検査 — v33.1 新契約 + v32.2 契約の不変 + qrn_core/readout_contract 封鎖不変",
            bad.is_empty(),
            if bad.is_empty() {
                "拒否は型エラー・棄却は裁定 — 二つを混ぜない契約が通った".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [B1] contexts 盲目性の機械実証 ----
    {
        let mut bad = Vec::new();
        let (net_a, _) = build_net(&site_gens, n, tau); // 文脈 ∅
        let (mut net_b, ids_b) = build_net(&site_gens, n, tau); // 文脈 atlas
        add_site_atlas(&mut net_b, &ids_b);
        // 旧手順 (v32.3 参照再現): contexts が違っても成分・読みは完全一致
        let comps_a = net_a.noncommutation_components();
        let comps_b = net_b.noncommutation_components();
        if comps_a != comps_b {
            bad.push("成分が contexts に依存した (旧手順の前提が破れ)".into());
        }
        let ref_a = reference_recover_v323(&net_a, &site_gens, n);
        let ref_b = reference_recover_v323(&net_b, &site_gens, n);
        let want = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 2, 2],
        };
        if !(ref_a == want && ref_b == want) {
            bad.push(format!(
                "参照手順の読みが不一致 (∅: {} / atlas: {})",
                ref_a.as_str(),
                ref_b.as_str()
            ));
        }
        // 修復入口: ∅ は構成時拒否・atlas は同じ正解
        let new_a = net_a.recovery_input();
        let refused = matches!(new_a, Err(RecoveryInputRejection::NoDeclaredContexts));
        if !refused {
            bad.push("文脈 ∅ の復元入力が拒否されない".into());
        }
        let new_b = net_b
            .recovery_input()
            .map(|inp| inp.recover().reading)
            .unwrap_or(FactorizationReading::Abstain(
                FactorizationAbstainReason::ComponentNotFactor,
            ));
        if new_b != want {
            bad.push(format!("修復入口の読みが {} (期待 Exact [2,2,2])", new_b.as_str()));
        }
        check(
            "[B1] contexts 盲目性 — 旧手順は ∅/atlas で読み完全一致 (入力でなかった)・修復入口は ∅ を構成時拒否",
            bad.is_empty(),
            format!(
                "参照手順: ∅ → {} = atlas → {} (一致)・修復入口: ∅ → no_declared_contexts (拒否) / atlas → {}",
                ref_a.as_str(),
                ref_b.as_str(),
                new_b.as_str()
            ),
        );
    }

    // ---- [B2] entangler 負制御 — primitive 選別の循環 ----
    {
        let mut bad = Vec::new();
        // 6 primitive の site net
        let (mut net6, ids6) = build_net(&site_gens, n, tau);
        add_site_atlas(&mut net6, &ids6);
        let det6 = net6.recovery_input().map(|i| i.recover().reading);
        // + independently accessible な entangler X₁X₂
        let mut gens7 = site_gens.clone();
        gens7.push(op3("XXI"));
        let (mut net7, ids7) = build_net(&gens7, n, tau);
        net7.add_context(&[ids7[0], ids7[2], ids7[4], ids7[6]]).unwrap(); // {X₁,X₂,X₃,X₁X₂}
        net7.add_context(&[ids7[1], ids7[3], ids7[5]]).unwrap(); // {Z₁,Z₂,Z₃}
        let det7 = net7.recovery_input().map(|i| i.recover().reading);
        // 旧手順でも同じ併合が起こる (選別の循環は lane に依らない)
        let ref7 = reference_recover_v323(&net7, &gens7, n);
        // 大域閉包は同一 (erasure no-go の対照)
        let cl6 = closure_of(&site_gens, n);
        let cl7 = closure_of(&gens7, n);
        let want222 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 2, 2],
        };
        let want24 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 4],
        };
        if det6 != Ok(want222.clone()) {
            bad.push(format!("site 6 本の読みが {:?}", det6.map(|r| r.as_str().to_string())));
        }
        if det7 != Ok(want24.clone()) {
            bad.push(format!("+X₁X₂ の読みが {:?}", det7.map(|r| r.as_str().to_string())));
        }
        if ref7 != want24 {
            bad.push(format!("旧手順の +X₁X₂ 読みが {}", ref7.as_str()));
        }
        if !(cl6.dim_algebra == 64 && cl7.dim_algebra == 64) {
            bad.push(format!("閉包次元 {} / {}", cl6.dim_algebra, cl7.dim_algebra));
        }
        check(
            "[B2] entangler 負制御 — X₁X₂ を 1 本加えるだけで Exact [2,2,2] → [2,4] (閉包は同一 M₈)",
            bad.is_empty(),
            format!(
                "6 本 → [2,2,2] / 7 本 → [2,4] (新旧 lane 一致)・閉包 dim {} = {} — 選別が答えを入力する循環は型修復では解けない (v33.2 Certified Laboratory Interface の動機)",
                cl6.dim_algebra, cl7.dim_algebra
            ),
        );
    }

    // ---- [B3] 修復入口の資格ゲート ----
    {
        let mut bad = Vec::new();
        // (a) role-mixed: 測定 effect の混入 → 構成時拒否
        let (mut net_m, ids_m) = build_net(&site_gens, n, tau);
        add_site_atlas(&mut net_m, &ids_m);
        let n1: Vec<C64> = ident(n)
            .iter()
            .zip(op3("ZII").iter())
            .map(|(i, z)| (*i - *z).scale(0.5))
            .collect();
        net_m
            .add_primitive(PrimitiveOperation {
                kind: OpKind::Measure(
                    MeasurementEffect::certify(
                        n1.iter().map(|c| c.re).collect(),
                        n1.iter().map(|c| c.im).collect(),
                        n,
                    )
                    .unwrap(),
                ),
                parity: OperatorParity::Even,
                provenance: "v331_measurement",
            })
            .unwrap();
        let rej_role = net_m.recovery_input().err();
        if rej_role != Some(RecoveryInputRejection::RoleMixedRecovery) {
            bad.push(format!("role-mixed の拒否が {:?}", rej_role.map(|r| r.as_str())));
        }
        // (b) 被覆不完全: Z₃ がどの文脈にも属さない → 構成時拒否
        let (mut net_c, ids_c) = build_net(&site_gens, n, tau);
        net_c.add_context(&[ids_c[0], ids_c[2], ids_c[4]]).unwrap();
        net_c.add_context(&[ids_c[1], ids_c[3]]).unwrap(); // {Z₁,Z₂} — Z₃ 欠落
        let rej_cov = net_c.recovery_input().err();
        if rej_cov != Some(RecoveryInputRejection::ContextCoverageIncomplete) {
            bad.push(format!("被覆不完全の拒否が {:?}", rej_cov.map(|r| r.as_str())));
        }
        // (c) qutrit × qubit (C⁶) — 資格を満たす非自明次元セルで旧解 [2,3] を再現
        let n6 = 6;
        let mut a3 = vec![C64::new(0.0, 0.0); 9];
        a3[1] = C64::new(1.0, 0.0);
        a3[3] = C64::new(1.0, 0.0);
        a3[5] = C64::new(1.0, 0.0);
        a3[7] = C64::new(1.0, 0.0);
        let mut b3 = vec![C64::new(0.0, 0.0); 9];
        b3[0] = C64::new(1.0, 0.0);
        b3[8] = C64::new(-1.0, 0.0);
        let i3 = ident(3);
        let i2 = ident(2);
        let gens6 = vec![
            kron(&a3, 3, &i2, 2),
            kron(&b3, 3, &i2, 2),
            kron(&i3, 3, &pauli('X'), 2),
            kron(&i3, 3, &pauli('Z'), 2),
        ];
        let (mut net_q, ids_q) = build_net(&gens6, n6, tau);
        net_q.add_context(&[ids_q[0], ids_q[2]]).unwrap(); // {path₃, X}
        net_q.add_context(&[ids_q[1], ids_q[3]]).unwrap(); // {diag₃, Z}
        let det_q = net_q.recovery_input().map(|i| i.recover().reading);
        let want23 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 3],
        };
        if det_q != Ok(want23) {
            bad.push(format!("C⁶ の読みが {:?}", det_q.map(|r| r.as_str().to_string())));
        }
        check(
            "[B3] 資格ゲート — role-mixed/文脈 0/被覆不完全は構成時拒否・資格を満たせば旧解を再現 (C⁶ → [2,3])",
            bad.is_empty(),
            format!(
                "拒否 3 種 = 型エラー (Abstain と混ぜない)・qutrit×qubit → {} — 生成子行列の出所は net の primitive のみ (gens 別渡しの経路なし)",
                if bad.is_empty() { "exact [2,3]" } else { "?" }
            ),
        );
    }

    // ---- [B4] 代数的可換 ↛ 操作的両立 (禁止変換 12) ----
    {
        let mut bad = Vec::new();
        // singleton 文脈のみ: 被覆は満たすが成分間の共同 addressability の証人がない
        let (mut net_s, ids_s) = build_net(&site_gens, n, tau);
        for id in &ids_s {
            net_s.add_context(&[*id]).unwrap();
        }
        let det_s = net_s.recovery_input().map(|i| i.recover().reading);
        let want_abst = FactorizationReading::Abstain(
            FactorizationAbstainReason::OperationalCompatibilityUnwitnessed,
        );
        if det_s != Ok(want_abst) {
            bad.push(format!(
                "singleton 文脈の読みが {:?} (期待 unwitnessed 棄却)",
                det_s.map(|r| r.as_str().to_string())
            ));
        }
        // 同じ素材を旧手順は Exact と読む (修復が変えた点の対照)
        let ref_s = reference_recover_v323(&net_s, &site_gens, n);
        if ref_s.as_str() != "exact_up_to_local_unitary_and_permutation" {
            bad.push(format!("旧手順の読みが {} (対照が立たない)", ref_s.as_str()));
        }
        // JointContextWitness の唯一の構成 = 宣言済み文脈の共有
        if net_s.joint_context_witness(ids_s[0], ids_s[2]).is_some() {
            bad.push("共有文脈なしに witness が構成できた".into());
        }
        let (mut net_w, ids_w) = build_net(&site_gens, n, tau);
        add_site_atlas(&mut net_w, &ids_w);
        let wit = net_w.joint_context_witness(ids_w[0], ids_w[2]);
        if wit.is_none() {
            bad.push("共有文脈があるのに witness が構成できない".into());
        }
        // 可換子証明書は全対に存在し definite — それでも昇格しない
        let all_certified = (0..ids_s.len()).all(|a| {
            ((a + 1)..ids_s.len())
                .all(|b| net_s.commutator_verdict(ids_s[a], ids_s[b]).is_some())
        });
        if !all_certified {
            bad.push("可換子証明書が欠けている (負制御の前提が破れ)".into());
        }
        check(
            "[B4] 代数的可換 ↛ 操作的両立 — singleton 文脈は全対証明書ありでも Abstain(unwitnessed)・旧手順は Exact (対照)",
            bad.is_empty(),
            format!(
                "修復入口 → operational_compatibility_unwitnessed / 旧手順 → {} — JointContextWitness の唯一の構成は共有文脈 (禁止変換 12)",
                ref_s.as_str()
            ),
        );
    }

    // ---- [B5] Liouvillian lane の型分離 (禁止変換 13) ----
    {
        let mut bad = Vec::new();
        let n4 = 4usize;
        // 決定的な非自明 H (traceless・非可換な項の和)
        let mut h = vec![C64::new(0.0, 0.0); 16];
        for (c, s) in [(0.7, "XI"), (0.4, "ZZ"), (0.3, "IY"), (0.2, "XX")] {
            for (hi, oi) in h.iter_mut().zip(op2(s).iter()) {
                *hi = *hi + oi.scale(c);
            }
        }
        let lane = HamiltonianCommutatorLiouvillian::certify(h.clone(), n4).unwrap();
        // (a) 可換子 lane は導分証明書を通り Ĥ を復元する
        let mut resid_h = f64::INFINITY;
        let mut leib_h = f64::INFINITY;
        match classify_generator(&|m: &[C64]| lane.apply(m), n4) {
            GeneratorClassification::HamiltonianCommutator {
                h_hat,
                leibniz_defect,
                reconstruction_residual,
            } => {
                leib_h = leibniz_defect;
                resid_h = reconstruction_residual;
                // Ĥ = H (traceless gauge) の直接照合
                let diff: f64 = h_hat
                    .iter()
                    .zip(h.iter())
                    .map(|(a, b)| (*a - *b).norm2())
                    .sum::<f64>()
                    .sqrt();
                let hn: f64 = h.iter().map(|c| c.norm2()).sum::<f64>().sqrt();
                if diff / hn > 1e-10 {
                    bad.push(format!("Ĥ 復元の相対誤差 {:.1e}", diff / hn));
                }
            }
            GeneratorClassification::NonDerivation { leibniz_defect } => {
                bad.push(format!("可換子 lane が NonDerivation ({:.1e})", leibniz_defect));
            }
        }
        if resid_h > 1e-10 || leib_h > 1e-12 {
            bad.push(format!("導分証明書 leib = {:.1e} / resid = {:.1e}", leib_h, resid_h));
        }
        // (b) GKLS (γ > 0) は Leibniz が破れる — γ 比例
        let sm = vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
        ];
        let jump = kron(&sm, 2, &pauli('I'), 2);
        let mk_gkls = |g: f64| {
            GklsLiouvillian::certify(h.clone(), vec![jump.clone()], vec![g], n4).unwrap()
        };
        let mut leib_gamma = Vec::new();
        for g in [0.15, 0.3] {
            let gk = mk_gkls(g);
            match classify_generator(&|m: &[C64]| gk.apply(m), n4) {
                GeneratorClassification::NonDerivation { leibniz_defect } => {
                    leib_gamma.push(leibniz_defect)
                }
                GeneratorClassification::HamiltonianCommutator { .. } => {
                    bad.push(format!("GKLS (γ = {}) が可換子 lane の資格を通った", g));
                    leib_gamma.push(0.0);
                }
            }
        }
        let ratio = if leib_gamma.len() == 2 && leib_gamma[0] > 0.0 {
            leib_gamma[1] / leib_gamma[0]
        } else {
            0.0
        };
        if leib_gamma.iter().any(|&l| l < 1e-3) || (ratio - 2.0).abs() > 1e-9 {
            bad.push(format!("Leibniz 破れが γ 比例でない ({:?}, 比 {:.6})", leib_gamma, ratio));
        }
        // (c) R⁽¹⁾ 公式は GKLS 測定と一致しない (γ → 0 で回復)
        let rho0: Vec<C64> = ident(n4).iter().map(|c| c.scale(0.25)).collect();
        let probes = [op2("XI"), op2("YI"), op2("ZI")];
        let eps = 0.05;
        let gk03 = mk_gkls(0.3);
        let gk00 = mk_gkls(0.0);
        let mut dev03 = 0.0f64;
        let mut dev00 = 0.0f64;
        let mut meas_vs_gkls = 0.0f64;
        for a in &probes {
            for b in &probes {
                let m03 = measure_r1_gkls(&gk03, &rho0, a, b, eps, n4);
                let m00 = measure_r1_gkls(&gk00, &rho0, a, b, eps, n4);
                let formula = lane.r1(b, a);
                dev03 = dev03.max((m03 - formula).abs());
                dev00 = dev00.max((m00 - formula).abs());
                // 測定は GKLS 自身の線形応答 Tr(B L(A)) とは一致する (乖離は数値でない)
                let gkls_pred = tr_prod_re(b, &gk03.apply(a), n4);
                meas_vs_gkls = meas_vs_gkls.max((m03 - gkls_pred).abs());
            }
        }
        if !(dev03 > 0.05 && dev00 <= 1e-6 && meas_vs_gkls <= 1e-6) {
            bad.push(format!(
                "応答負制御 dev(γ=0.3) = {:.3} / dev(γ=0) = {:.1e} / 測定 vs GKLS 予測 = {:.1e}",
                dev03, dev00, meas_vs_gkls
            ));
        }
        // (d) γ = 0 の縮退点は可換子 lane の資格を通る (dissipator_strength = 0)
        let degen_ok = matches!(
            classify_generator(&|m: &[C64]| gk00.apply(m), n4),
            GeneratorClassification::HamiltonianCommutator { .. }
        ) && gk00.dissipator_strength() == 0.0;
        if !degen_ok {
            bad.push("γ = 0 の縮退点が可換子 lane の資格を通らない".into());
        }
        check(
            "[B5] Liouvillian lane 分離 — 導分証明書 (Ĥ 復元 rel ≤ 1e-10)・GKLS は Leibniz γ 比例破れ・R⁽¹⁾ 公式は GKLS で不成立",
            bad.is_empty(),
            format!(
                "可換子 lane: leib {:.1e}/resid {:.1e}・GKLS leib 破れ {:.4} (γ 比 {:.6} = 2)・R⁽¹⁾ 乖離: γ=0.3 → {:.3} / γ=0 → {:.1e} (測定 vs GKLS 予測 {:.1e})",
                leib_h,
                resid_h,
                leib_gamma.last().cloned().unwrap_or(0.0),
                ratio,
                dev03,
                dev00,
                meas_vs_gkls
            ),
        );
    }

    // ---- [B6] 封鎖の source/schema/文書検査 ----
    {
        let mut bad = Vec::new();
        // (a) 禁止 impl From の不在 (sim/src 全走査 — 本監査自身はパターン定義を含むため対象外)
        const FORBIDDEN_IMPLS: [&str; 5] = [
            "impl From<CertifiedCommutator",
            "impl From<HamiltonianCommutatorLiouvillian",
            "impl From<GklsLiouvillian",
            "impl From<MarkedRecoveryInput",
            "impl From<JointContextWitness",
        ];
        const EXEMPT: [&str; 1] = ["v331_scope_repair.rs"];
        let mut scan = |dir: &str| {
            if let Ok(rdir) = fs::read_dir(format!("{}/{}", root, dir)) {
                for e in rdir.filter_map(|e| e.ok()) {
                    let p = e.path();
                    let name = p
                        .file_name()
                        .map(|f| f.to_string_lossy().to_string())
                        .unwrap_or_default();
                    if EXEMPT.contains(&name.as_str())
                        || !p.extension().map(|x| x == "rs").unwrap_or(false)
                    {
                        continue;
                    }
                    if let Ok(text) = fs::read_to_string(&p) {
                        for f in FORBIDDEN_IMPLS {
                            if text.contains(f) {
                                bad.push(format!("{}: {}", p.display(), f));
                            }
                        }
                    }
                }
            }
        };
        scan("sim/src");
        scan("sim/src/bin");
        // (b) operational_net.rs の v33.1 ゲート実在 + v32.2 ゲート不変
        let src = rd("sim/src/operational_net.rs").unwrap_or_default();
        for needle in [
            "禁止変換 12",
            "禁止変換 13",
            "RoleMixedRecovery",
            "NoDeclaredContexts",
            "ContextCoverageIncomplete",
            "OperationalCompatibilityUnwitnessed",
            "pub fn recovery_input",
            "禁止変換 11",
            "ACCEPTS_ODD",
        ] {
            if !src.contains(needle) {
                bad.push(format!("operational_net.rs: 「{}」が無い", needle));
            }
        }
        // (c) v323 旧署名の凍結保存 (黙って直さない — 歴史的原本の確認)
        let v323 = rd("sim/src/bin/v323_factorization_recovery.rs").unwrap_or_default();
        if !(v323.contains("gens: &[Vec<C64>]") && v323.contains("fn recover_factorization")) {
            bad.push("v323 の旧署名 (gens 別渡し) が凍結保存されていない".into());
        }
        // (d) schema: 新概念 + 禁止変換 12/13
        let schema = rd("core.schema.yml").unwrap_or_default();
        for needle in [
            "- name: MarkedRecoveryInput",
            "- name: JointContextWitness",
            "- name: HamiltonianCommutatorLiouvillian",
            "- name: GklsLiouvillian",
        ] {
            if !schema.contains(needle) {
                bad.push(format!("core.schema.yml: 「{}」が無い", needle));
            }
        }
        if !schema.contains("- from: CertifiedCommutator\n  to: JointContextWitness\n  reason:") {
            bad.push("禁止変換 12 (CertifiedCommutator → JointContextWitness) が未登録".into());
        }
        if !schema
            .contains("- from: HamiltonianCommutatorLiouvillian\n  to: GklsLiouvillian\n  reason:")
        {
            bad.push(
                "禁止変換 13 (HamiltonianCommutatorLiouvillian → GklsLiouvillian) が未登録".into(),
            );
        }
        // (e) 文書アンカー
        let doc = rd("docs/uft-v33.1.md").unwrap_or_default();
        for needle in [
            "境界監査",
            "primitive 選別の循環",
            "OperationalCompatibilityUnwitnessed",
            "禁止変換 12",
            "禁止変換 13",
            "Leibniz",
            "role-mixed recovery",
            "gens 別渡しの廃止",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v33.1.md: 「{}」が無い", needle));
            }
        }
        check(
            "[B6] 封鎖の source/schema/文書 — 禁止 impl From 不在・禁止変換 12/13 登録・v323 凍結保存・文書アンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "修復は型・schema・文書の三点で凍結され、旧署名は歴史として保存された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "contexts は定理入力になり、代数的可換と操作的両立・可換子 lane と GKLS は型で分離された — 残る循環 (accessible primitive の選別) が v33.2 の主題"
        } else {
            "**型スコープの破れ** — operational_net と schema/文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
