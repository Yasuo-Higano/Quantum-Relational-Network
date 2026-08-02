//! v33.2 Certified Laboratory Interface — 操作の出自証明と controller-free no-go (PROMPT/14)
//!
//! v33.1 が固定した循環 (primitive 選別が答えを入力する) の正側修復。目標は
//! 「物理的に可能な操作を状態や Hamiltonian だけから自動導出する」ことでは**ない**
//! (それは原理的に不可能 — 本版 [C5] が no-go として機械化)。目標は
//!
//! > 「accessible operations を入力として与えた」から、「各 operation が、どの
//! > command・校正・合成列・誤差・資源によって accessible なのかを証明した」へ
//!
//! 進めることである。共有契約 = `sim/src/laboratory_interface.rs`。
//!
//!   [C0] 契約自己検査 — laboratory_interface_self_test (新) +
//!        operational_net/scope_repair/qrn_core/readout_contract の不変
//!   [C1] **DeclaredOperation は資格なし (禁止変換 14)**: 宣言 (行列 + 意図) から
//!        AccessibleOperation への直接変換は存在しない — 門は較正・合成・トモグラフィ
//!        の 3 証明書のみ。証明書は対象行列の sha256 に**結束** (X₁ の合成証明書を
//!        X₂ に流用 → 構成時拒否)。資源予算は成分半順序 (比較不能対の実在 — 単一
//!        スカラー化 Ord なし)
//!   [C2] **IndependentAddressabilityCertificate**: 独立 site knobs (6 command 1:1)
//!        → rank 6・σ_min ≈ 1・cross-talk 0 で資格。cross-coupling ε = 0.05 は
//!        バー内で記録つき資格・ε = 0.3 は CrosstalkExcess で拒否
//!   [C3] **tied control no-go (禁止変換 15)**: 装置が u(t)(X₁+X₂) しか持たないとき
//!        {X₁, X₂} への数学的分解は rank 1 < 2 で構成時拒否。正直な interface
//!        (tied 1 本) の net は Abstain(InsufficientOperationalGenerators) — HOLD-9
//!        tied セルの器械
//!   [C4] **可アクセス性は interface との関係** (SynthesisCertificate): 同じ X₁ が
//!        interface A = {X₁+X₂} では合成路なし (Lie 閉包残差 0.71)・interface B =
//!        {X₁+X₂, Z₂} では bracket 3 手 + 線形 1 手で Synthesized 資格 (残差 0・
//!        機械実行で検証)
//!   [C5] **controller-free decomposition no-go (E3-A)**: 同一の (H, drift H = 0,
//!        ρ = I/8) に対し、4 つの certified interface が非同値な読みを生成 —
//!        site → Exact [2,2,2] (orbit α)・DFT → Exact [2,2,2] (orbit β ≠ α,
//!        matching 不在)・+entangler → Exact [2,4]・tied → Abstain。ゆえに状態・
//!        Hamiltonian・大域代数だけから accessible operations を選ぶ写像は存在しない
//!        (第三の no-go: 状態単独 ✗ [v31.4]・閉包 ✗ [v32.2]・controller-free ✗ [本版])
//!   [C6] **role-typed 文脈 4 型 (禁止変換 16)**: MeasurementContext = joint
//!        measurability (**可換性より広い** — 非可換 unsharp 対 η = 0.6 は明示的
//!        joint POVM で資格・η = 0.8 は正値性破れで拒否 [不偏 qubit 対の canonical
//!        構成 — Busch の iff の器械化])・トモグラフィ出自 (情報完全 6 状態から
//!        effect 再構成 — 偏りデータは残差バーで拒否)・PreparationFamily = 凸可達
//!        (重み機械検証・到達不能標的は拒否)・DriftRegime = 安定性 (変動バー)
//!   [C7] **end-to-end**: 全 primitive が出自証明書つきの AccessibleOperationalNet
//!        (DeclaredOperation を受ける口が無い) から v33.1 修復入口で Exact [2,2,2]
//!   [C8] 封鎖の source/schema/文書検査 — 禁止変換 14/15/16 登録・impl From 不在・
//!        ResourceBudget の Ord/PartialOrd 不在・uft-v33.2.md アンカー
//!
//! 実行: cargo run --release --bin v332_certified_interface

use std::fs;
use std::path::Path;
use uft_sim::laboratory_interface::*;
use uft_sim::operational_net::{
    commutator, hs_inner, hs_norm, CertifiedCommutator, ControlGenerator, DriftGenerator,
    FactorizationAbstainReason, FactorizationReading, OpId, OpKind, OperatorParity,
    OrdinaryCommutation, Preparation,
};
use uft_sim::C64;

// ---------------------------------------------------------------- Pauli / kron 素子 (v322/v331 と同一)

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
    use uft_sim::operational_net::{cdag, cmul};
    cmul(&cmul(v, a, n), &cdag(v, n), n)
}

/// gauge orbit 照合 (v32.3 [F3] と同一の判定器 — 成分 traceless 部分代数の置換 matching)
fn same_gauge_orbit(a: &[Vec<Vec<C64>>], b: &[Vec<Vec<C64>>]) -> (bool, f64) {
    if a.len() != b.len() {
        return (false, 0.0);
    }
    let k = a.len();
    let overlap = |u: &Vec<Vec<C64>>, w: &Vec<Vec<C64>>| -> f64 {
        if u.len() != w.len() {
            return 0.0;
        }
        let mut acc = 0.0;
        for x in w {
            for y in u {
                acc += hs_inner(y, x).norm2();
            }
        }
        acc / (u.len() as f64)
    };
    let mut perm: Vec<usize> = (0..k).collect();
    let mut best = 0.0f64;
    let mut found = false;
    loop {
        let mut minov = f64::INFINITY;
        for i in 0..k {
            minov = minov.min(overlap(&a[i], &b[perm[i]]));
        }
        best = best.max(minov);
        if minov >= 1.0 - 1e-9 {
            found = true;
            break;
        }
        let mut i = k as isize - 2;
        while i >= 0 && perm[i as usize] >= perm[(i + 1) as usize] {
            i -= 1;
        }
        if i < 0 {
            break;
        }
        let mut j = k - 1;
        while perm[j] <= perm[i as usize] {
            j -= 1;
        }
        perm.swap(i as usize, j);
        perm[(i as usize + 1)..].reverse();
    }
    (found, best)
}

// ---------------------------------------------------------------- certified interface の組立て

fn mk_ctrl(g: &[C64], n: usize) -> OpKind {
    OpKind::Control(
        ControlGenerator::certify(
            g.iter().map(|c| c.re).collect(),
            g.iter().map(|c| c.im).collect(),
            n,
        )
        .unwrap(),
    )
}

/// commands = targets (1:1 直接較正) の interface から資格つき net を建てる
fn build_certified_net(
    gens: &[Vec<C64>],
    n: usize,
    tau: f64,
    contexts: &[Vec<usize>],
) -> Result<AccessibleOperationalNet<OrdinaryCommutation>, String> {
    let cert = certify_addressability(gens, gens, n, 0.5, 0.1)
        .map_err(|e| format!("addressability: {}", e.as_str()))?;
    let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 1.0, 1e-9).unwrap();
    let mut net: AccessibleOperationalNet<OrdinaryCommutation> =
        AccessibleOperationalNet::new(n, tau);
    let mut ids: Vec<OpId> = Vec::new();
    for g in gens {
        let op = AccessibleOperation::certify(
            mk_ctrl(g, n),
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

fn main() {
    uft_sim::self_test();
    println!("=== v33.2 Certified Laboratory Interface — 操作の出自証明と controller-free no-go (PROMPT/14) ===\n");
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
    let tau = 1e-3;

    // ---- [C0] 契約自己検査 ----
    {
        let mut bad = Vec::new();
        if let Err(e) = laboratory_interface_self_test() {
            bad.push(format!("laboratory_interface_self_test: {}", e));
        }
        if let Err(e) = uft_sim::operational_net::operational_net_self_test() {
            bad.push(format!("operational_net_self_test: {}", e));
        }
        if let Err(e) = uft_sim::operational_net::scope_repair_self_test() {
            bad.push(format!("scope_repair_self_test: {}", e));
        }
        if let Err(e) = uft_sim::qrn_core::qrn_core_self_test() {
            bad.push(format!("qrn_core_self_test: {}", e));
        }
        if let Err(e) = uft_sim::readout_contract::readout_contract_self_test() {
            bad.push(format!("readout_contract_self_test: {}", e));
        }
        check(
            "[C0] 契約自己検査 — laboratory_interface (新) + operational_net/scope_repair/qrn_core/readout_contract 不変",
            bad.is_empty(),
            if bad.is_empty() {
                "出自は型 + sha256 結束・拒否は構成時 — 既存契約は不変".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    let n4 = 4usize;
    let x1 = op2("XI");
    let x2 = op2("IX");
    let z2 = op2("IZ");
    let x3_8 = op3("IIX");
    let tied: Vec<C64> = x1.iter().zip(x2.iter()).map(|(a, b)| *a + *b).collect();

    // ---- [C1] DeclaredOperation は資格なし (禁止変換 14) ----
    {
        let mut bad: Vec<String> = Vec::new();
        // 宣言 (行列 + 意図) — これだけでは AccessibleOperation を構成する路がない
        let declared = DeclaredOperation {
            re: x1.iter().map(|c| c.re).collect(),
            im: x1.iter().map(|c| c.im).collect(),
            dim: n4,
            intent: RoleIntent::Control,
        };
        let m = declared.matrix_c64();
        // 門 1 (較正) を通ると資格が立つ
        let cert = certify_addressability(&[m.clone()], &[m.clone()], n4, 0.5, 0.1);
        if cert.is_err() {
            bad.push("較正の門が開かない".into());
        }
        let cert = cert.unwrap();
        let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 1.0, 1e-9).unwrap();
        let ok_op = AccessibleOperation::certify(
            mk_ctrl(&m, n4),
            OperatorParity::Even,
            OperationOrigin::DirectlyCalibrated(cert.clone()),
            cert.clone(),
            budget,
        );
        if ok_op.is_err() {
            bad.push("較正済み操作の受理が拒否された".into());
        }
        // 証明書の流用 (X₁ の較正証明書を X₂ に付ける) は構成時拒否
        match AccessibleOperation::certify(
            mk_ctrl(&x2, n4),
            OperatorParity::Even,
            OperationOrigin::DirectlyCalibrated(cert.clone()),
            cert.clone(),
            budget,
        ) {
            Err(InterfaceRejection::CertificateTargetMismatch) => {}
            _ => bad.push("X₁ の証明書が X₂ に流用できた".into()),
        }
        // 資源予算は成分半順序 — 比較不能対
        let b1 = ResourceBudget::certify(1.0, 2.0, 1.0, 1.0, 1e-3).unwrap();
        let b2 = ResourceBudget::certify(2.0, 1.0, 1.0, 1.0, 1e-3).unwrap();
        if b1.comparable(&b2) {
            bad.push("資源予算が全順序化されている".into());
        }
        check(
            "[C1] 宣言 ≠ 資格 (禁止変換 14) — 門は 3 証明書のみ・sha256 結束 (流用拒否)・予算は成分半順序",
            bad.is_empty(),
            format!(
                "declared(X₁, intent=control) → 較正の門で資格 / X₂ への流用 → certificate_target_mismatch / 予算 (1,2,·) vs (2,1,·) 比較不能 = {}",
                !b1.comparable(&b2)
            ),
        );
    }

    // ---- [C2] IndependentAddressabilityCertificate ----
    {
        let mut bad = Vec::new();
        let n8 = 8usize;
        let site: Vec<Vec<C64>> = ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"]
            .iter()
            .map(|s| op3(s))
            .collect();
        let cert = certify_addressability(&site, &site, n8, 0.5, 0.1);
        let (rank, sigma, xtalk) = match &cert {
            Ok(c) => (
                c.command_jacobian_rank.rank,
                c.smallest_singular_value.hi(),
                c.crosstalk_bound.hi(),
            ),
            Err(_) => (0, 0.0, 9.9),
        };
        if !(cert.is_ok() && rank == 6 && (sigma - 1.0).abs() < 1e-9 && xtalk < 1e-9) {
            bad.push(format!("独立 knobs: rank {} σ {} xtalk {}", rank, sigma, xtalk));
        }
        // cross-coupling ε: G_k = X_k + ε X_{k+1} (標的は X_k)
        let mk_eps = |eps: f64| -> Vec<Vec<C64>> {
            vec![
                site[0].iter().zip(site[2].iter()).map(|(a, b)| *a + b.scale(eps)).collect(),
                site[2].iter().zip(site[4].iter()).map(|(a, b)| *a + b.scale(eps)).collect(),
                site[4].clone(),
            ]
        };
        let targets3 = vec![site[0].clone(), site[2].clone(), site[4].clone()];
        let c_ok = certify_addressability(&targets3, &mk_eps(0.05), n8, 0.5, 0.1);
        let xt_rec = c_ok.as_ref().map(|c| c.crosstalk_bound.hi()).unwrap_or(9.9);
        if !(c_ok.is_ok() && xt_rec > 0.04 && xt_rec < 0.06) {
            bad.push(format!("ε = 0.05 の記録つき資格が立たない (xtalk {})", xt_rec));
        }
        let c_bad = certify_addressability(&targets3, &mk_eps(0.3), n8, 0.5, 0.1);
        match c_bad {
            Err(InterfaceRejection::CrosstalkExcess) => {}
            _ => bad.push("ε = 0.3 が拒否されない".into()),
        }
        check(
            "[C2] addressability 証明書 — 独立 6 knobs: rank 6・σ_min 1・xtalk 0 / ε = 0.05 記録つき資格 / ε = 0.3 拒否",
            bad.is_empty(),
            format!(
                "rank = {}・σ_min = {:.9}・xtalk = {:.1e} / ε = 0.05 → 記録 {:.6} (バー 0.1 内) / ε = 0.3 → crosstalk_excess",
                rank, sigma, xtalk, xt_rec
            ),
        );
    }

    // ---- [C3] tied control no-go (禁止変換 15) ----
    {
        let mut bad = Vec::new();
        // 数学的分解の申告: 標的 {X₁, X₂}・command は同一の X₁+X₂
        let dec = certify_addressability(
            &[x1.clone(), x2.clone()],
            &[tied.clone(), tied.clone()],
            n4,
            0.5,
            0.1,
        );
        match dec {
            Err(InterfaceRejection::InsufficientCommandRank) => {}
            _ => bad.push("数学的分解が rank 不足で拒否されない".into()),
        }
        // 正直な interface: tied 1 本を 1 標的として較正 → net は Abstain
        let net_tied = build_certified_net(&[tied.clone()], n4, tau, &[vec![0]]);
        let reading = net_tied
            .and_then(|net| net.recover().map_err(|e| e.as_str().to_string()))
            .map(|d| d.reading);
        let want = FactorizationReading::Abstain(
            FactorizationAbstainReason::InsufficientOperationalGenerators,
        );
        if reading.as_ref() != Ok(&want) {
            bad.push(format!("tied net の読みが {:?}", reading.map(|r| r.as_str().to_string())));
        }
        check(
            "[C3] tied control no-go (禁止変換 15) — {X₁,X₂} への数学的分解は rank 1 < 2 で拒否・正直な tied net は Abstain",
            bad.is_empty(),
            "u(t)(X₁+X₂) しか無い装置に「二つの独立 primitive」は立たない — 読みは insufficient_operational_generators (HOLD-9 tied セルの器械)".into(),
        );
    }

    // ---- [C4] 可アクセス性は interface との関係 (SynthesisCertificate) ----
    {
        let mut bad = Vec::new();
        // interface A = {X₁+X₂}: X₁ への合成路なし
        let resid_a = synthesis_path_residual(&[tied.clone()], &x1, n4);
        if resid_a < 0.5 {
            bad.push(format!("interface A の Lie 閉包残差が小さすぎる ({:.3})", resid_a));
        }
        // interface B = {X₁+X₂, Z₂}: bracket 3 手 + 線形 1 手で X₁
        let base = vec![tied.clone(), z2.clone()];
        let steps = vec![
            SynthStep::BracketOverI(0, 1),
            SynthStep::BracketOverI(2, 1),
            SynthStep::Linear(vec![(1.0, 0), (0.25, 3)]),
        ];
        let synth = certify_synthesis(&base, &steps, &x1, n4, 1e-9);
        let (depth, resid_b) = match &synth {
            Ok(c) => (c.depth, c.residual.hi()),
            Err(_) => (0, 9.9),
        };
        if synth.is_err() || depth != 3 || resid_b > 1e-9 {
            bad.push(format!("interface B の合成資格が立たない (depth {} resid {:.1e})", depth, resid_b));
        }
        // 資格つき Synthesized 操作として受理される (sha256 結束)
        if let Ok(cert) = synth {
            let addr = certify_addressability(&[tied.clone(), z2.clone()], &[tied.clone(), z2.clone()], n4, 0.5, 0.1).unwrap();
            let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 3.0, 1e-9).unwrap();
            if AccessibleOperation::certify(
                mk_ctrl(&x1, n4),
                OperatorParity::Even,
                OperationOrigin::Synthesized(cert),
                addr,
                budget,
            )
            .is_err()
            {
                bad.push("Synthesized X₁ の受理が拒否された".into());
            }
        }
        check(
            "[C4] 可アクセス性は関係 — 同じ X₁: interface A = {X₁+X₂} は路なし・B = {X₁+X₂, Z₂} は Synthesized 資格",
            bad.is_empty(),
            format!(
                "A: Lie 閉包への相対残差 {:.6} (> 0.5 = 路なし) / B: bracket 2 手 + 線形 1 手 (depth 3)・残差 ≤ {:.1e} — 作用素単体の属性ではない",
                resid_a, resid_b
            ),
        );
    }

    // ---- [C5] controller-free decomposition no-go (E3-A) ----
    {
        let mut bad = Vec::new();
        let n8 = 8usize;
        // 同一の (H = (C²)⊗³, drift H = 0, ρ = I/8) — interface だけを変える
        let rho = Preparation::certify(
            (0..64).map(|k| if k % 9 == 0 { 0.125 } else { 0.0 }).collect(),
            vec![0.0; 64],
            n8,
        );
        if rho.is_err() {
            bad.push("共通 ρ = I/8 が準備資格を通らない".into());
        }
        let site: Vec<Vec<C64>> = ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"]
            .iter()
            .map(|s| op3(s))
            .collect();
        let v = dft8();
        let mode: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&v, g, n8)).collect();
        let mut ent = site.clone();
        ent.push(op3("XXI"));
        let ctx_site = vec![vec![0usize, 2, 4], vec![1usize, 3, 5]];
        let ctx_ent = vec![vec![0usize, 2, 4, 6], vec![1usize, 3, 5]];
        let net_site = build_certified_net(&site, n8, tau, &ctx_site);
        let net_mode = build_certified_net(&mode, n8, tau, &ctx_site);
        let net_ent = build_certified_net(&ent, n8, tau, &ctx_ent);
        let tied8: Vec<C64> = op3("XII")
            .iter()
            .zip(op3("IXI").iter())
            .map(|(a, b)| *a + *b)
            .collect();
        let net_tied = build_certified_net(&[tied8], n8, tau, &[vec![0]]);
        let recover =
            |net: Result<AccessibleOperationalNet<OrdinaryCommutation>, String>| match net {
                Ok(net) => net
                    .recover()
                    .map(|d| (d.reading, d.component_subalgebras))
                    .map_err(|e| e.as_str().to_string()),
                Err(e) => Err(e),
            };
        let r_site = recover(net_site);
        let r_mode = recover(net_mode);
        let r_ent = recover(net_ent);
        let r_tied = recover(net_tied);
        let want222 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 2, 2],
        };
        let want24 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 4],
        };
        let want_abst = FactorizationReading::Abstain(
            FactorizationAbstainReason::InsufficientOperationalGenerators,
        );
        let mut orbit_overlap = 1.0f64;
        let mut orbit_same = true;
        match (&r_site, &r_mode, &r_ent, &r_tied) {
            (Ok((a, sa)), Ok((b, sb)), Ok((c, _)), Ok((d, _))) => {
                if a != &want222 || b != &want222 {
                    bad.push(format!("site/DFT の読みが {} / {}", a.as_str(), b.as_str()));
                }
                let (same, ov) = same_gauge_orbit(sa, sb);
                orbit_same = same;
                orbit_overlap = ov;
                if same || ov > 0.9 {
                    bad.push(format!("site×DFT の orbit が一致してしまう (overlap {:.4})", ov));
                }
                if c != &want24 {
                    bad.push(format!("+entangler の読みが {}", c.as_str()));
                }
                if d != &want_abst {
                    bad.push(format!("tied の読みが {}", d.as_str()));
                }
            }
            _ => bad.push(format!(
                "interface の構成に失敗: {:?} {:?} {:?} {:?}",
                r_site.as_ref().err(),
                r_mode.as_ref().err(),
                r_ent.as_ref().err(),
                r_tied.as_ref().err()
            )),
        }
        check(
            "[C5] controller-free no-go (E3-A) — 同一 (H, drift 0, ρ = I/8) で 4 interface が非同値な読み",
            bad.is_empty(),
            format!(
                "site → [2,2,2] (orbit α) / DFT → [2,2,2] (orbit β — matching {} overlap {:.4}) / +entangler → [2,4] / tied → abstain — 状態・Hamiltonian・大域代数だけから accessible operations を選ぶ写像は存在しない",
                if orbit_same { "あり" } else { "不在" },
                orbit_overlap
            ),
        );
    }

    // ---- [C6] role-typed 文脈 4 型 (禁止変換 16) ----
    {
        let mut bad = Vec::new();
        let n2 = 2usize;
        let i2 = pauli('I');
        let px = pauli('X');
        let pz = pauli('Z');
        // (a) MeasurementContext: 非可換 unsharp 対の joint measurability
        let mk_pair = |eta: f64| -> (Vec<Vec<C64>>, Vec<Vec<C64>>) {
            let e_x: Vec<C64> = (0..4)
                .map(|k| i2[k].scale(0.5) + px[k].scale(0.5 * eta))
                .collect();
            let e_z: Vec<C64> = (0..4)
                .map(|k| i2[k].scale(0.5) + pz[k].scale(0.5 * eta))
                .collect();
            let mut joint = Vec::new();
            for s in [1.0f64, -1.0] {
                for t in [1.0f64, -1.0] {
                    joint.push(
                        (0..4)
                            .map(|k| {
                                i2[k].scale(0.25)
                                    + px[k].scale(0.25 * s * eta)
                                    + pz[k].scale(0.25 * t * eta)
                            })
                            .collect::<Vec<C64>>(),
                    );
                }
            }
            (vec![e_x, e_z], joint)
        };
        let marg = vec![vec![0usize, 1], vec![0usize, 2]];
        let (eff6, joint6) = mk_pair(0.6);
        let comm_norm = hs_norm(&commutator(&eff6[0], &eff6[1], n2));
        let j_ok = certify_joint_measurement(&eff6, &joint6, &marg, n2);
        if !(j_ok.is_ok() && comm_norm > 0.2) {
            bad.push(format!(
                "非可換 unsharp 対 (‖[E,F]‖ = {:.4}) の joint 資格が立たない: {:?}",
                comm_norm,
                j_ok.err().map(|e| e.as_str())
            ));
        }
        let (eff8, joint8) = mk_pair(0.8);
        match certify_joint_measurement(&eff8, &joint8, &marg, n2) {
            Err(InterfaceRejection::JointCandidateNotPositive) => {}
            r => bad.push(format!("η = 0.8 が拒否されない: {:?}", r.err().map(|e| e.as_str()))),
        }
        // (b) トモグラフィ出自: 情報完全 6 状態から effect 再構成
        let bloch = |x: f64, y: f64, z: f64| -> Vec<C64> {
            (0..4)
                .map(|k| {
                    i2[k].scale(0.5)
                        + px[k].scale(0.5 * x)
                        + pauli('Y')[k].scale(0.5 * y)
                        + pz[k].scale(0.5 * z)
                })
                .collect()
        };
        let states = vec![
            bloch(0.0, 0.0, 1.0),
            bloch(0.0, 0.0, -1.0),
            bloch(1.0, 0.0, 0.0),
            bloch(-1.0, 0.0, 0.0),
            bloch(0.0, 1.0, 0.0),
            bloch(0.0, -1.0, 0.0),
        ];
        let e_true: Vec<C64> = (0..4).map(|k| i2[k].scale(0.5) + px[k].scale(0.25)).collect();
        let probs: Vec<f64> = states
            .iter()
            .map(|s| {
                let mut acc = 0.0;
                for i in 0..n2 {
                    for k in 0..n2 {
                        let a = e_true[i * n2 + k];
                        let b = s[k * n2 + i];
                        acc += a.re * b.re - a.im * b.im;
                    }
                }
                acc
            })
            .collect();
        let tomo = certify_effect_tomography(&states, &probs, n2, 1e-3);
        let resid_ok = tomo.as_ref().map(|(_, c)| c.residual.hi()).unwrap_or(9.9);
        if tomo.is_err() || resid_ok > 1e-9 {
            bad.push(format!("正確データのトモグラフィが立たない (resid {:.1e})", resid_ok));
        }
        let mut probs_biased = probs.clone();
        probs_biased[2] += 0.02;
        match certify_effect_tomography(&states, &probs_biased, n2, 1e-3) {
            Err(InterfaceRejection::TomographyResidualExcess) => {}
            r => bad.push(format!(
                "偏りデータが拒否されない: {:?}",
                r.err().map(|e| e.as_str())
            )),
        }
        // (c) PreparationFamily: 凸可達 (w = 0.625) / 到達不能 (z = 0.9)
        let prep = |z: f64| {
            Preparation::certify(
                vec![0.5 * (1.0 + z), 0.0, 0.0, 0.5 * (1.0 - z)],
                vec![0.0; 4],
                n2,
            )
            .unwrap()
        };
        let fam = [prep(-0.3), prep(0.5)];
        let reach = certify_convex_reachability(&fam, &[0.375, 0.625], &prep(0.2));
        if reach.is_err() {
            bad.push("z = 0.2 の凸可達が立たない".into());
        }
        match certify_convex_reachability(&fam, &[0.0, 1.0], &prep(0.9)) {
            Err(InterfaceRejection::MixtureMismatch) => {}
            r => bad.push(format!(
                "z = 0.9 (到達不能) が拒否されない: {:?}",
                r.err().map(|e| e.as_str())
            )),
        }
        // (d) DriftRegime: 安定 (0.02 ≤ 0.05) / 不安定 (0.5)
        let drift = |amp: f64| {
            DriftGenerator::certify(
                pz.iter().map(|c| c.re * amp).collect(),
                vec![0.0; 4],
                n2,
            )
            .unwrap()
        };
        let stab = certify_drift_stability(&[drift(1.0), drift(1.0 + 0.02 / hs_norm(&pz))], 0.05);
        if stab.is_err() {
            bad.push("安定 regime の資格が立たない".into());
        }
        match certify_drift_stability(&[drift(1.0), drift(1.0 + 0.5 / hs_norm(&pz))], 0.05) {
            Err(InterfaceRejection::DriftRegimeUnstable) => {}
            r => bad.push(format!(
                "不安定 drift が拒否されない: {:?}",
                r.err().map(|e| e.as_str())
            )),
        }
        check(
            "[C6] role-typed 文脈 — joint measurability は可換性より広い (η=0.6 資格/0.8 拒否)・トモグラフィ/凸可達/安定性",
            bad.is_empty(),
            format!(
                "‖[E^X,E^Z]‖ = {:.4} ≠ 0 でも joint POVM 資格 (禁止変換 16 の裏側)・η=0.8 は正値性破れ・トモグラフィ resid {:.1e} (偏り +0.02 は拒否)・凸重み (0.375, 0.625)・drift 変動バー 0.05",
                comm_norm, resid_ok
            ),
        );
    }

    // ---- [C7] end-to-end: 全 primitive 出自つきの復元 ----
    {
        let mut bad = Vec::new();
        let n8 = 8usize;
        let site: Vec<Vec<C64>> = ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"]
            .iter()
            .map(|s| op3(s))
            .collect();
        let net = build_certified_net(&site, n8, tau, &[vec![0, 2, 4], vec![1, 3, 5]]);
        let mut origin_kinds = Vec::new();
        let reading = match &net {
            Ok(net) => {
                for i in 0..6u32 {
                    origin_kinds.push(net.origin(OpId(i)).as_str());
                }
                net.recover().map(|d| d.reading).map_err(|e| e.as_str().to_string())
            }
            Err(e) => Err(e.clone()),
        };
        let want = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 2, 2],
        };
        if reading.as_ref() != Ok(&want) {
            bad.push(format!("end-to-end の読みが {:?}", reading));
        }
        if !origin_kinds.iter().all(|k| *k == "directly_calibrated") {
            bad.push(format!("出自の欠落: {:?}", origin_kinds));
        }
        // DeclaredOperation を net に入れる口は存在しない — 型 (source 検査は [C8])。
        // 較正されていない行列 (X₃ 標的外) を無理に admit する路が拒否されること:
        let cert1 = certify_addressability(&[x1.clone()], &[x1.clone()], n4, 0.5, 0.1).unwrap();
        let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 1.0, 1e-9).unwrap();
        match AccessibleOperation::certify(
            mk_ctrl(&x3_8[..].to_vec(), 8),
            OperatorParity::Even,
            OperationOrigin::DirectlyCalibrated(cert1.clone()),
            cert1,
            budget,
        ) {
            Err(InterfaceRejection::CertificateTargetMismatch) => {}
            _ => bad.push("未較正の行列が受理された".into()),
        }
        check(
            "[C7] end-to-end — 全 primitive が出自証明書つきの AccessibleOperationalNet から v33.1 入口で Exact [2,2,2]",
            bad.is_empty(),
            format!(
                "出自 = {:?} (文字列 provenance なし — sha256 結束)・読み = exact [2,2,2]・未較正行列は構成時拒否",
                origin_kinds.first().unwrap_or(&"?")
            ),
        );
    }

    // ---- [C8] 封鎖の source/schema/文書検査 ----
    {
        let mut bad = Vec::new();
        // パターンは分割リテラルで組む (本版以降の規約): 走査パターン自身が他の
        // 走査器 (v331 [B6] 等) の誤検出源にならないため。除外リストは歴史的
        // リテラルを含む v331 のみで足りる。
        let forbidden_impls: [String; 6] = [
            format!("impl From{}", "<DeclaredOperation"),
            format!("impl From{}", "<AccessibleOperation"),
            format!("impl From{}", "<CertifiedCommutator"),
            format!("impl Ord {}", "for ResourceBudget"),
            format!("impl PartialOrd {}", "for ResourceBudget"),
            format!("impl From{}", "<OperationOrigin"),
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
                        for f in &forbidden_impls {
                            if text.contains(f.as_str()) {
                                bad.push(format!("{}: {}", p.display(), f));
                            }
                        }
                    }
                }
            }
        };
        scan("sim/src");
        scan("sim/src/bin");
        let src = rd("sim/src/laboratory_interface.rs").unwrap_or_default();
        for needle in [
            "禁止変換 14",
            "禁止変換 15",
            "禁止変換 16",
            "DeclaredOperation",
            "CertificateTargetMismatch",
            "InsufficientCommandRank",
            "componentwise_le",
            "certified_accessible_operation",
        ] {
            if !src.contains(needle) {
                bad.push(format!("laboratory_interface.rs: 「{}」が無い", needle));
            }
        }
        let schema = rd("core.schema.yml").unwrap_or_default();
        for needle in [
            "- name: DeclaredOperation",
            "- name: AccessibleOperation",
            "- name: OperationOrigin",
            "- name: IndependentAddressabilityCertificate",
            "- name: SynthesisCertificate",
            "- name: TomographyCertificate",
            "- name: ResourceBudget",
            "- name: JointMeasurementCertificate",
            "- name: AccessibleOperationalNet",
            "- name: ControlContext",
            "- name: MeasurementContext",
            "- name: PreparationFamily",
            "- name: DriftRegime",
        ] {
            if !schema.contains(needle) {
                bad.push(format!("core.schema.yml: 「{}」が無い", needle));
            }
        }
        if !schema.contains("- from: DeclaredOperation\n  to: AccessibleOperation\n  reason:") {
            bad.push("禁止変換 14 が未登録".into());
        }
        if !schema
            .contains("- from: OperatorDecomposition\n  to: IndependentlyAddressablePrimitives\n  reason:")
        {
            bad.push("禁止変換 15 が未登録".into());
        }
        if !schema
            .contains("- from: CertifiedCommutator\n  to: JointMeasurementCertificate\n  reason:")
        {
            bad.push("禁止変換 16 が未登録".into());
        }
        let doc = rd("docs/uft-v33.2.md").unwrap_or_default();
        for needle in [
            "Certified Laboratory Interface",
            "DeclaredOperation",
            "禁止変換 14",
            "禁止変換 15",
            "禁止変換 16",
            "joint measurability",
            "文字列 provenance の廃止",
            "controller-free",
            "成分半順序",
            "sha256 結束",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v33.2.md: 「{}」が無い", needle));
            }
        }
        check(
            "[C8] 封鎖の source/schema/文書 — 禁止変換 14/15/16 登録・impl From/Ord 不在・アンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "宣言から資格への直接路はコンパイル不能・予算の全順序化も型で不在".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "操作の出自が証明書になった — 可アクセス性は interface との関係であり、controller-free な因子分解写像は存在しない (第三の no-go)。資源の semantics (budget profile) が v33.3 の主題"
        } else {
            "**interface 契約の破れ** — laboratory_interface と schema/文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
