//! v33.5 Graded recovery の正しい境界 — Majorana locality ≠ Dirac locality (PROMPT/14)
//!
//! 現行の graded lane の課題表現は「odd 演算子からモード構造を復元」だったが、
//! **odd CAR だけから通常の複素 fermion mode を一意に得ることはできない**: 2N 本の
//! Majorana の CAR は O(2N) に不変で、二本ずつの組 (Dirac pairing) には直交複素構造
//! J (J² = −I) が追加で要る。共有契約 = `sim/src/graded_recovery.rs`。
//!
//! > **Majorana locality と Dirac-mode locality は同じ識別問題ではない。後者は
//! > 追加の U(1) charge / complex structure witness を必要とする。**
//!
//!   [M0] 契約自己検査 — graded_recovery_self_test (新) + 既存 5 契約の不変
//!   [M1] **odd CAR は O(2N) 不変・graded graph は空**: JW 6 本の CAR 資格・
//!        Givens 回転 (γ₂,γ₃ 混合) 後も CAR 完全一致・graded bracket ノルムは全対
//!        厳密 0 (ordinary では 2√8 の K₆ に見える — v32.2 [N5] の罠の再確認) —
//!        「どの二本が一組か」の情報は CAR に無い
//!   [M2] **Dirac pairing no-go (禁止変換 20)**: witness なしの読みは
//!        MajoranaFrameOnly (O(2N) orbit)。標準 pairing (12)(34)(56) と回転 pairing
//!        の両方が完全な mode-CAR を満たす — CAR データからは区別不能で、
//!        MajoranaFrame → ComplexModeFactorization の witness なし昇格は存在しない
//!   [M3] **正側: charge witness → 複素構造 → U(N) を除く回復**: Q = Σ n_i から
//!        J (実・反対称・J² = −I 厳密) を抽出 → 3 モード回復 — 全 mode-CAR 厳密・
//!        Σ â†â = Q 厳密再現。回転 frame 座標でも同じ physical content (Σ n̂ = Q)
//!   [M4] **縮退/非線形 witness → Abstain**: 部分 charge n₁ は J² ≠ −I →
//!        ComplexStructureUnresolved・quartic 汚染 (Q + 0.3 γ₁γ₂γ₃γ₄) は adjoint
//!        作用が frame 上で閉じず WitnessNotLinearOnFrame・微小汚染 (1e-12) は資格
//!        (バー規律)
//!   [M5] **既存復元器は捏造しない**: ordinary net は odd を構成時拒否 (v32.2 ゲート
//!        継承)・graded net (反可換子証明書 — 全対 0) の marked recovery は成分が
//!        全て単本で閉包 dim 2 → Abstain(ComponentNotFactor) — モード構造は graded
//!        graph から出ない (witness 経路が唯一)
//!   [M6] 封鎖の schema/文書検査 — 概念登録 + 禁止変換 20・impl From 不在・アンカー
//!
//! 実行: cargo run --release --bin v335_graded_recovery

use std::fs;
use std::path::Path;
use uft_sim::graded_recovery::*;
use uft_sim::operational_net::{
    anticommutator, commutator, hs_norm, CertifiedCommutator, ControlGenerator,
    FactorizationAbstainReason, FactorizationReading, FermionicZ2Graded, OpKind, OperationalNet,
    OperatorParity, OrdinaryCommutation, PrimitiveOperation,
};
use uft_sim::C64;

// ---------------------------------------------------------------- Pauli / kron 素子

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

fn add_m(a: &[C64], b: &[C64], ca: f64, cb: f64) -> Vec<C64> {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| x.scale(ca) + y.scale(cb))
        .collect()
}

/// mode-CAR の最大偏差: {a_i, a_j†} = δ_ij I・{a_i, a_j} = 0
fn mode_car_defect(modes: &[(Vec<C64>, Vec<C64>)], n: usize) -> f64 {
    let mut dev = 0.0f64;
    for (i, (ai, _)) in modes.iter().enumerate() {
        for (j, (aj, ajd)) in modes.iter().enumerate() {
            let ac1 = anticommutator(ai, ajd, n);
            for r in 0..n {
                for c in 0..n {
                    let want = if i == j && r == c { 1.0 } else { 0.0 };
                    let x = ac1[r * n + c];
                    dev = dev.max((x.re - want).hypot(x.im));
                }
            }
            let ac2 = anticommutator(ai, aj, n);
            dev = dev.max(hs_norm(&ac2));
        }
    }
    dev
}

fn main() {
    uft_sim::self_test();
    println!("=== v33.5 Graded recovery — Majorana locality ≠ Dirac locality (PROMPT/14) ===\n");
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

    // JW Majorana 6 本 (3 モード)
    let gammas: Vec<Vec<C64>> = ["XII", "YII", "ZXI", "ZYI", "ZZX", "ZZY"]
        .iter()
        .map(|s| op3(s))
        .collect();
    // 標準 charge Q = Σ n_i = (3I − Z₁ − Z₂ − Z₃)/2
    let q_std = {
        let i8 = op3("III");
        let mut q = i8.iter().map(|c| c.scale(1.5)).collect::<Vec<C64>>();
        for s in ["ZII", "IZI", "IIZ"] {
            q = add_m(&q, &op3(s), 1.0, -0.5);
        }
        q
    };
    // Givens 回転 (γ₂, γ₃ を θ = 0.7 で混合) の O(6) 行列
    let theta = 0.7f64;
    let r6 = {
        let mut r = vec![0.0f64; 36];
        for i in 0..6 {
            r[i * 6 + i] = 1.0;
        }
        let (c, s) = (theta.cos(), theta.sin());
        r[1 * 6 + 1] = c;
        r[1 * 6 + 2] = -s;
        r[2 * 6 + 1] = s;
        r[2 * 6 + 2] = c;
        r
    };

    // ---- [M0] 契約自己検査 ----
    {
        let mut bad = Vec::new();
        if let Err(e) = uft_sim::graded_recovery::graded_recovery_self_test() {
            bad.push(format!("graded_recovery_self_test: {}", e));
        }
        if let Err(e) = uft_sim::contextual_factorization::contextual_factorization_self_test() {
            bad.push(format!("contextual_factorization_self_test: {}", e));
        }
        if let Err(e) = uft_sim::resource_profile::resource_profile_self_test() {
            bad.push(format!("resource_profile_self_test: {}", e));
        }
        if let Err(e) = uft_sim::laboratory_interface::laboratory_interface_self_test() {
            bad.push(format!("laboratory_interface_self_test: {}", e));
        }
        if let Err(e) = uft_sim::operational_net::operational_net_self_test() {
            bad.push(format!("operational_net_self_test: {}", e));
        }
        check(
            "[M0] 契約自己検査 — graded_recovery (新) + 既存 4 契約の不変",
            bad.is_empty(),
            if bad.is_empty() {
                "Dirac pairing は witness の関数 — CAR だけの昇格は型に無い".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [M1] odd CAR は O(2N) 不変・graded graph は空 ----
    {
        let mut bad = Vec::new();
        let frame = MajoranaFrame::certify(gammas.clone(), n);
        if frame.is_err() {
            bad.push("JW 6 本が CAR 資格を通らない".into());
        }
        let frame = frame.unwrap();
        let rot = frame.rotated(&r6);
        if rot.is_err() {
            bad.push("O(6) 回転 frame が CAR 資格を通らない (CAR が不変でない?)".into());
        }
        // graded graph: 全対の反可換子ノルム 0・ordinary は 2√8
        let mut max_graded = 0.0f64;
        let mut min_ord = f64::INFINITY;
        for a in 0..6 {
            for b in (a + 1)..6 {
                max_graded = max_graded.max(hs_norm(&anticommutator(&gammas[a], &gammas[b], n)));
                min_ord = min_ord.min(hs_norm(&commutator(&gammas[a], &gammas[b], n)));
            }
        }
        if !(max_graded < 1e-12 && (min_ord - 2.0 * (8.0f64).sqrt()).abs() < 1e-9) {
            bad.push(format!("graded max {:.1e} / ordinary min {:.6}", max_graded, min_ord));
        }
        check(
            "[M1] odd CAR は O(2N) 不変・graded graph は空 (ordinary では K₆ に見える罠)",
            bad.is_empty(),
            format!(
                "CAR 資格: 原 frame ✓・Givens 回転 frame ✓ (完全一致)・graded 全対 ‖{{γ,γ'}}‖ ≤ {:.1e}・ordinary ‖[γ,γ']‖ = {:.6} — pairing の情報は CAR に無い",
                max_graded, min_ord
            ),
        );
    }

    // ---- [M2] Dirac pairing no-go (禁止変換 20) ----
    {
        let mut bad = Vec::new();
        let frame = MajoranaFrame::certify(gammas.clone(), n).unwrap();
        // witness なし → MajoranaFrameOnly
        match recover_graded(&frame, None) {
            GradedRecoveryReading::MajoranaFrameOnly { n_majorana: 6 } => {}
            r => bad.push(format!("witness なしの読みが {}", r.as_str())),
        }
        // 標準 pairing と回転 pairing はどちらも完全な mode-CAR を満たす
        let mk_modes = |fr: &MajoranaFrame| -> Vec<(Vec<C64>, Vec<C64>)> {
            (0..3)
                .map(|i| {
                    let gv = fr.gamma(2 * i).to_vec();
                    let gw = fr.gamma(2 * i + 1).to_vec();
                    let a: Vec<C64> = gv
                        .iter()
                        .zip(gw.iter())
                        .map(|(x, y)| C64::new((x.re - y.im) * 0.5, (x.im + y.re) * 0.5))
                        .collect();
                    let ad: Vec<C64> = gv
                        .iter()
                        .zip(gw.iter())
                        .map(|(x, y)| C64::new((x.re + y.im) * 0.5, (x.im - y.re) * 0.5))
                        .collect();
                    (a, ad)
                })
                .collect()
        };
        let frame_rot = frame.rotated(&r6).unwrap();
        let car_std = mode_car_defect(&mk_modes(&frame), n);
        let car_rot = mode_car_defect(&mk_modes(&frame_rot), n);
        if !(car_std < 1e-9 && car_rot < 1e-9) {
            bad.push(format!("pairing CAR: 標準 {:.1e} / 回転 {:.1e}", car_std, car_rot));
        }
        // 2 つの pairing のモード数演算子は実際に異なる (n̂₁ の差のノルム > 0.1)
        let n_of = |modes: &Vec<(Vec<C64>, Vec<C64>)>| -> Vec<C64> {
            uft_sim::operational_net::cmul(&modes[0].1, &modes[0].0, n)
        };
        let d: f64 = n_of(&mk_modes(&frame))
            .iter()
            .zip(n_of(&mk_modes(&frame_rot)).iter())
            .map(|(x, y)| (*x - *y).norm2())
            .sum::<f64>()
            .sqrt();
        if d < 0.1 {
            bad.push(format!("2 つの pairing が同じモードを与えた ({:.3})", d));
        }
        check(
            "[M2] Dirac pairing no-go — witness なしは MajoranaFrameOnly・非同値な 2 pairing がともに完全な mode-CAR",
            bad.is_empty(),
            format!(
                "標準 (12)(34)(56) と Givens 回転 pairing: mode-CAR 偏差 {:.1e} / {:.1e}・n̂₁ の差 {:.4} — CAR データは選べない (禁止変換 20: witness なし昇格の門は無い)",
                car_std, car_rot, d
            ),
        );
    }

    // ---- [M3] 正側: charge witness → 複素構造 → U(N) を除く回復 ----
    {
        let mut bad = Vec::new();
        let frame = MajoranaFrame::certify(gammas.clone(), n).unwrap();
        let w_res = extract_complex_structure(&frame, &q_std);
        let mut struct_resid = f64::INFINITY;
        let mut sum_n: Vec<C64> = vec![C64::new(0.0, 0.0); n * n];
        let mut car_dev = f64::INFINITY;
        let mut n_modes = 0usize;
        match &w_res {
            Err(e) => bad.push(format!("標準 charge の抽出が {}", e.as_str())),
            Ok(w) => {
                struct_resid = w.structure_residual;
                match recover_graded(&frame, Some(w)) {
                    GradedRecoveryReading::ComplexModeFactorization { n_modes: k, modes } => {
                        n_modes = k;
                        car_dev = mode_car_defect(&modes, n);
                        for (a, ad) in &modes {
                            let nn = uft_sim::operational_net::cmul(ad, a, n);
                            for (s, x) in sum_n.iter_mut().zip(nn.iter()) {
                                *s = *s + *x;
                            }
                        }
                    }
                    r => bad.push(format!("witness つきの読みが {}", r.as_str())),
                }
            }
        }
        let q_dev: f64 = sum_n
            .iter()
            .zip(q_std.iter())
            .map(|(x, y)| (*x - *y).norm2())
            .sum::<f64>()
            .sqrt();
        if !(n_modes == 3 && car_dev < 1e-9 && q_dev < 1e-9) {
            bad.push(format!(
                "モード {} / CAR {:.1e} / Σn̂ − Q {:.1e}",
                n_modes, car_dev, q_dev
            ));
        }
        // 回転 frame 座標でも同じ physical content (Σ n̂ = Q)
        let frame_rot = frame.rotated(&r6).unwrap();
        let w_rot = extract_complex_structure(&frame_rot, &q_std);
        let mut q_dev_rot = f64::INFINITY;
        if let Ok(wr) = &w_rot {
            if let GradedRecoveryReading::ComplexModeFactorization { modes, .. } =
                recover_graded(&frame_rot, Some(wr))
            {
                let mut s2: Vec<C64> = vec![C64::new(0.0, 0.0); n * n];
                for (a, ad) in &modes {
                    let nn = uft_sim::operational_net::cmul(ad, a, n);
                    for (s, x) in s2.iter_mut().zip(nn.iter()) {
                        *s = *s + *x;
                    }
                }
                q_dev_rot = s2
                    .iter()
                    .zip(q_std.iter())
                    .map(|(x, y)| (*x - *y).norm2())
                    .sum::<f64>()
                    .sqrt();
            }
        }
        if q_dev_rot > 1e-9 {
            bad.push(format!("回転 frame 座標の Σn̂ − Q = {:.1e}", q_dev_rot));
        }
        check(
            "[M3] 正側 — Q = Σn_i から J (実・反対称・J² = −I) を抽出し 3 モード回復 (U(N) gauge を除く)",
            bad.is_empty(),
            format!(
                "J² + I 残差 {:.1e}・mode-CAR {:.1e}・Σâ†â = Q 残差 {:.1e} (回転 frame 座標でも {:.1e}) — witness が O(6) orbit から pairing を選ぶ",
                struct_resid, car_dev, q_dev, q_dev_rot
            ),
        );
    }

    // ---- [M4] 縮退/非線形 witness → Abstain ----
    {
        let mut bad = Vec::new();
        let frame = MajoranaFrame::certify(gammas.clone(), n).unwrap();
        // 部分 charge n₁ = (I − Z₁)/2 → J² ≠ −I
        let q_part = add_m(&op3("III"), &op3("ZII"), 0.5, -0.5);
        match extract_complex_structure(&frame, &q_part) {
            Err(GradedAbstainReason::ComplexStructureUnresolved) => {}
            r => bad.push(format!(
                "部分 charge が棄却されない: {:?}",
                r.err().map(|e| e.as_str())
            )),
        }
        // quartic 汚染 (γ₁γ₂γ₃γ₄ = (iγ₁γ₂)(iγ₃γ₄)·(−1) = −Z₁Z₂) → 線形性破れ
        let quartic = op3("ZZI").iter().map(|c| c.scale(-1.0)).collect::<Vec<C64>>();
        let q_dirty = add_m(&q_std, &quartic, 1.0, 0.3);
        match extract_complex_structure(&frame, &q_dirty) {
            Err(GradedAbstainReason::WitnessNotLinearOnFrame) => {}
            r => bad.push(format!(
                "quartic 汚染が棄却されない: {:?}",
                r.err().map(|e| e.as_str())
            )),
        }
        // 微小汚染 (1e-12) はバー内で資格
        let q_tiny = add_m(&q_std, &quartic, 1.0, 1e-12);
        if extract_complex_structure(&frame, &q_tiny).is_err() {
            bad.push("微小汚染 (1e-12) が資格を通らない (バー規律の破れ)".into());
        }
        check(
            "[M4] 縮退/非線形 witness — 部分 charge は ComplexStructureUnresolved・quartic 汚染は WitnessNotLinearOnFrame・微小汚染は資格",
            bad.is_empty(),
            "n₁ だけでは残り 4 本の pairing が決まらない (J² ≠ −I)・0.3 γ₁γ₂γ₃γ₄ は adjoint 作用が frame 外へ漏れる・1e-12 はバー内".into(),
        );
    }

    // ---- [M5] 既存復元器は捏造しない ----
    {
        let mut bad = Vec::new();
        // ordinary net は odd を構成時拒否 (v32.2 ゲートの継承)
        let mk_odd = |g: &Vec<C64>| PrimitiveOperation {
            kind: OpKind::Control(
                ControlGenerator::certify(
                    g.iter().map(|c| c.re).collect(),
                    g.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            parity: OperatorParity::Odd,
            provenance: "v335_majorana",
        };
        let mut net_o: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, 1e-3);
        if net_o.add_primitive(mk_odd(&gammas[0])).is_ok() {
            bad.push("ordinary net が odd を受理した".into());
        }
        // graded net: 反可換子証明書 (全対 0 = graded-commuting) → marked recovery は
        // 成分単本・閉包 dim 2 → Abstain(ComponentNotFactor)
        let mut net_g: OperationalNet<FermionicZ2Graded> = OperationalNet::new(n, 1e-3);
        let ids: Vec<_> = gammas
            .iter()
            .map(|g| net_g.add_primitive(mk_odd(g)).unwrap())
            .collect();
        for a in 0..6 {
            for b in (a + 1)..6 {
                let nu = hs_norm(&anticommutator(&gammas[a], &gammas[b], n));
                net_g.set_commutator(
                    ids[a],
                    ids[b],
                    CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
                );
            }
        }
        net_g.add_context(&ids).unwrap();
        let reading = net_g.recovery_input().map(|i| i.recover().reading);
        let want = FactorizationReading::Abstain(FactorizationAbstainReason::ComponentNotFactor);
        if reading.as_ref().map(|r| r == &want) != Ok(true) {
            bad.push(format!(
                "graded net の読みが {:?} (期待 Abstain(ComponentNotFactor))",
                reading.map(|r| r.as_str().to_string())
            ));
        }
        check(
            "[M5] 既存復元器は捏造しない — ordinary は odd 構成時拒否・graded net の marked recovery は Abstain(ComponentNotFactor)",
            bad.is_empty(),
            "graded graph (全対 0) の成分は単本 ×6・閉包 dim 2 は factor でない — モード構造は graded graph から出ない (witness 経路が唯一)".into(),
        );
    }

    // ---- [M6] 封鎖の schema/文書検査 ----
    {
        let mut bad = Vec::new();
        let forbidden_impls: [String; 2] = [
            format!("impl From{}", "<MajoranaFrame"),
            format!("impl From{}", "<ComplexStructureWitness"),
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
        let src = rd("sim/src/graded_recovery.rs").unwrap_or_default();
        for needle in [
            "禁止変換 20",
            "MajoranaFrameOnly",
            "ComplexModeFactorization",
            "ComplexStructureUnresolved",
            "WitnessNotLinearOnFrame",
            "O(2N)",
        ] {
            if !src.contains(needle) {
                bad.push(format!("graded_recovery.rs: 「{}」が無い", needle));
            }
        }
        let schema = rd("core.schema.yml").unwrap_or_default();
        for needle in [
            "- name: MajoranaFrame",
            "- name: ComplexStructureWitness",
            "- name: GradedRecoveryReading",
        ] {
            if !schema.contains(needle) {
                bad.push(format!("core.schema.yml: 「{}」が無い", needle));
            }
        }
        if !schema.contains("- from: MajoranaFrame\n  to: ComplexModeFactorization\n  reason:") {
            bad.push("禁止変換 20 が未登録".into());
        }
        let doc = rd("docs/uft-v33.5.md").unwrap_or_default();
        for needle in [
            "Majorana locality",
            "Dirac",
            "O(2N)",
            "複素構造",
            "禁止変換 20",
            "MajoranaFrameOnly",
            "ComplexStructureUnresolved",
            "witness",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v33.5.md: 「{}」が無い", needle));
            }
        }
        check(
            "[M6] 封鎖の schema/文書 — 概念登録 + 禁止変換 20・impl From 不在・アンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "CAR から pairing への直接路はコンパイル不能 — witness が唯一の門".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "Majorana locality と Dirac locality は型で分離された — pairing は charge witness の関数であり、CAR 単独の読みは O(2N) orbit で止まる。構造化 backend のスケーリングが v33.6 の主題"
        } else {
            "**graded 契約の破れ** — graded_recovery と schema/文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
