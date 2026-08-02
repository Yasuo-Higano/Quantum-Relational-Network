//! v33.4 Contextual factorization — chart 局所復元と overlap glue (PROMPT/14)
//!
//! 実際の実験室は持ち場 (chart) の集まりであり、全系を一望する制御器はどこにも
//! ない。v33.4 は「context ごとの局所復元 → overlap 上の algebra matching →
//! glue 裁定」を機械化する。共有契約 = `sim/src/contextual_factorization.rs`。
//!
//!   [A0] 契約自己検査 — contextual_factorization_self_test (新) + 既存 5 契約不変
//!   [A1] **整合 atlas → GluedExact + glue 定理**: chart A = {X₁,Z₁,X₂,Z₂}・
//!        B = {X₂,Z₂,X₃,Z₃} (qubit 2 を共有)。局所復元は各 [2,2]・共有因子の
//!        matching (overlap 1) で束ね、被覆 Π d = 8・大域証人 (bridge {X₁,X₃} 込み)
//!        → GluedExact [2,2,2]。**直接大域復元 (v33.1 入口) と読み・gauge orbit が
//!        一致** — atlas glue = 大域復元の機械定理
//!   [A2] **変成不変**: chart B の qubit-2 操作を局所 unitary u₂ で回した net でも
//!        因子部分代数は集合として不変 (u M₂ u† = M₂) — glue は同じ [2,2,2]・
//!        同じ orbit (frame 変換は glue を破らない)
//!   [A3] **cocycle 不整合 → Abstain**: chart B を entangler W = CZ₂₃ で捻ると
//!        B の因子は W M₂ W† (overlap 1/3) — matching が破れ、全被覆の整合 chart 群
//!        も無い (A も B も単独では Π d = 4 ≠ 8) → Abstain(GlueInconsistent)。
//!        chart 単独の局所 Exact は大域主張にならない (**禁止変換 19**)
//!   [A4] **複数 glue → EquivalenceClassOnly**: site charts (S1, S2) + DFT chart D
//!        の atlas — site 群と D 群がそれぞれ全被覆で内部整合・相互 matching は
//!        破れ (overlap ≈ 0.56) → 整合 atlas 2 つ = EquivalenceClassOnly{2}
//!        (無制約 tie-break の禁止 — v32.3 [F3] の atlas 版)
//!   [A5] **witness 境界の両 lane 一致**: bridge 文脈 {X₁,X₃} を外すと glue は
//!        Abstain(CompatibilityUnwitnessed)・直接大域復元も
//!        Abstain(OperationalCompatibilityUnwitnessed) — v33.1 の証人規律を atlas が
//!        継承する (glue で規律が緩まない)
//!   [A6] 封鎖の schema/文書検査 — 概念登録 + 禁止変換 19・impl From 不在・アンカー
//!
//! 実行: cargo run --release --bin v334_contextual_factorization

use std::fs;
use std::path::Path;
use uft_sim::contextual_factorization::*;
use uft_sim::operational_net::{
    cdag, cmul, commutator, hs_norm, same_gauge_orbit, CertifiedCommutator, ControlGenerator,
    FactorizationAbstainReason, FactorizationReading, OpId, OpKind, OperationalNet,
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

/// net の組立て: 全対 exact 証明書 + 文脈 (id 添字で指定)
fn build_net(
    gens: &[Vec<C64>],
    n: usize,
    tau: f64,
    contexts: &[Vec<usize>],
) -> (OperationalNet<OrdinaryCommutation>, Vec<OpId>) {
    let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
    let mut ids = Vec::new();
    for g in gens {
        ids.push(
            net.add_primitive(PrimitiveOperation {
                kind: OpKind::Control(
                    ControlGenerator::certify(
                        g.iter().map(|c| c.re).collect(),
                        g.iter().map(|c| c.im).collect(),
                        n,
                    )
                    .unwrap(),
                ),
                parity: OperatorParity::Even,
                provenance: "v334_control",
            })
            .unwrap(),
        );
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
        net.add_context(&members).unwrap();
    }
    (net, ids)
}

fn chart(ids: &[OpId], sel: &[usize]) -> ChartSpec {
    ChartSpec {
        primitive_ids: sel.iter().map(|&i| ids[i]).collect(),
    }
}

fn main() {
    uft_sim::self_test();
    println!("=== v33.4 Contextual factorization — chart 局所復元と overlap glue (PROMPT/14) ===\n");
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
    let tau = 1e-3;

    // 素材: site 6 本 (並び: X₁ Z₁ X₂ Z₂ X₃ Z₃)
    let site: Vec<Vec<C64>> = ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"]
        .iter()
        .map(|s| op3(s))
        .collect();
    // [A1] の文脈: A 内 {X₁,X₂}/{Z₁,Z₂}・B 内 {X₂,X₃}/{Z₂,Z₃}・bridge {X₁,X₃}
    let ctx_a1: Vec<Vec<usize>> = vec![
        vec![0, 2],
        vec![1, 3],
        vec![2, 4],
        vec![3, 5],
        vec![0, 4],
    ];

    // ---- [A0] 契約自己検査 ----
    {
        let mut bad = Vec::new();
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
        if let Err(e) = uft_sim::operational_net::scope_repair_self_test() {
            bad.push(format!("scope_repair_self_test: {}", e));
        }
        check(
            "[A0] 契約自己検査 — contextual_factorization (新) + 既存 4 契約の不変",
            bad.is_empty(),
            if bad.is_empty() {
                "chart は持ち場しか語らない・大域は glue だけが与える".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [A1] 整合 atlas → GluedExact + glue 定理 ----
    {
        let mut bad = Vec::new();
        let (net, ids) = build_net(&site, n, tau, &ctx_a1);
        let ch_a = chart(&ids, &[0, 1, 2, 3]);
        let ch_b = chart(&ids, &[2, 3, 4, 5]);
        // chart 局所復元 (それぞれ因子 2 つ)
        let la = recover_chart(&net, &ch_a, n);
        let lb = recover_chart(&net, &ch_b, n);
        let dims_of = |l: &Result<ChartLocalFactorization, ChartFailure>| -> Vec<usize> {
            l.as_ref()
                .map(|f| {
                    let mut d: Vec<usize> = f.factors.iter().map(|(d, _, _)| *d).collect();
                    d.sort_unstable();
                    d
                })
                .unwrap_or_default()
        };
        if dims_of(&la) != vec![2, 2] || dims_of(&lb) != vec![2, 2] {
            bad.push(format!("chart 局所因子 {:?} / {:?}", dims_of(&la), dims_of(&lb)));
        }
        // atlas glue
        let atlas = recover_atlas(&net, &[ch_a, ch_b], n);
        let (glue_dims, glue_subs) = match &atlas {
            AtlasReading::GluedExact {
                local_dims,
                factor_subalgebras,
            } => (local_dims.clone(), factor_subalgebras.clone()),
            r => {
                bad.push(format!("glue が {}", r.as_str()));
                (Vec::new(), Vec::new())
            }
        };
        if glue_dims != vec![2, 2, 2] {
            bad.push(format!("glue dims {:?}", glue_dims));
        }
        // glue 定理: 直接大域復元と読み・orbit 一致
        let direct = net.recovery_input().map(|i| i.recover());
        let mut orbit_ok = false;
        if let Ok(det) = &direct {
            if det.reading
                == (FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
                    local_dims: vec![2, 2, 2],
                })
            {
                let (same, _) = same_gauge_orbit(&glue_subs, &det.component_subalgebras);
                orbit_ok = same;
            }
        }
        if !orbit_ok {
            bad.push("直接大域復元との orbit 一致が立たない".into());
        }
        check(
            "[A1] 整合 atlas → GluedExact [2,2,2]・glue 定理 (直接大域復元と読み・orbit 一致)",
            bad.is_empty(),
            "chart A/B とも局所 [2,2]・共有 qubit 2 の因子 matching (overlap 1) で束ね・被覆 Π = 8・証人完備 — atlas glue = 大域復元".into(),
        );
    }

    // ---- [A2] 変成不変 (局所 unitary は glue を破らない) ----
    {
        let mut bad = Vec::new();
        // u₂ = exp(iθ n·σ) on qubit 2
        let rot = |theta: f64, nx: f64, ny: f64, nz: f64| -> Vec<C64> {
            let (c, s) = (theta.cos(), theta.sin());
            vec![
                C64::new(c, s * nz),
                C64::new(s * ny, s * nx),
                C64::new(-s * ny, s * nx),
                C64::new(c, -s * nz),
            ]
        };
        let u2 = {
            let u = rot(0.6, 0.6, 0.0, 0.8);
            let a = kron(&pauli('I'), 2, &u, 2);
            kron(&a, 4, &pauli('I'), 2)
        };
        // B の qubit-2 操作だけ u₂ 共役 (chart B' は別 frame で qubit 2 を制御)
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
        // 文脈: A 内 {X₁,X₂}/{Z₁,Z₂}・B' 内 {X₂',X₃}/{Z₂',Z₃}・bridge {X₁,X₃}
        let ctxs: Vec<Vec<usize>> = vec![
            vec![0, 2],
            vec![1, 3],
            vec![4, 6],
            vec![5, 7],
            vec![0, 6],
        ];
        let (net, ids) = build_net(&gens, n, tau, &ctxs);
        let ch_a = chart(&ids, &[0, 1, 2, 3]);
        let ch_bp = chart(&ids, &[4, 5, 6, 7]);
        let atlas = recover_atlas(&net, &[ch_a, ch_bp], n);
        let ok = match &atlas {
            AtlasReading::GluedExact { local_dims, .. } => local_dims == &vec![2, 2, 2],
            _ => false,
        };
        if !ok {
            bad.push(format!("変成 atlas の glue が {}", atlas.as_str()));
        }
        check(
            "[A2] 変成不変 — chart B の qubit-2 frame を u₂ で回しても glue は同じ [2,2,2] (u M₂ u† = M₂)",
            bad.is_empty(),
            "因子部分代数は集合として frame 非依存 — matching は局所 unitary に不変 (変成対照)".into(),
        );
    }

    // ---- [A3] cocycle 不整合 → Abstain (禁止変換 19 の執行) ----
    {
        let mut bad = Vec::new();
        // W = CZ₂₃: X₂ → X₂Z₃, Z₂ → Z₂, X₃ → Z₂X₃, Z₃ → Z₃
        let cz23 = {
            let mut m = vec![C64::new(0.0, 0.0); 64];
            for b in 0..8usize {
                let q2 = (b >> 1) & 1;
                let q3 = b & 1;
                let sign = if q2 == 1 && q3 == 1 { -1.0 } else { 1.0 };
                m[b * 8 + b] = C64::new(sign, 0.0);
            }
            m
        };
        let bpp: Vec<Vec<C64>> = [&site[2], &site[3], &site[4], &site[5]]
            .iter()
            .map(|g| conj_by(&cz23, g, n))
            .collect();
        let gens: Vec<Vec<C64>> = vec![
            site[0].clone(),
            site[1].clone(),
            site[2].clone(),
            site[3].clone(),
            bpp[0].clone(), // X₂Z₃
            bpp[1].clone(), // Z₂
            bpp[2].clone(), // Z₂X₃
            bpp[3].clone(), // Z₃
        ];
        // 文脈: A 内 ×2・B'' 内 {X₂Z₃, Z₂X₃} (可換) と {Z₂, Z₃}
        let ctxs: Vec<Vec<usize>> = vec![vec![0, 2], vec![1, 3], vec![4, 6], vec![5, 7]];
        let (net, ids) = build_net(&gens, n, tau, &ctxs);
        let ch_a = chart(&ids, &[0, 1, 2, 3]);
        let ch_bpp = chart(&ids, &[4, 5, 6, 7]);
        // 各 chart は局所 Exact (因子 2 つずつ) — だが大域は立たない
        let la = recover_chart(&net, &ch_a, n).is_ok();
        let lb = recover_chart(&net, &ch_bpp, n).is_ok();
        let atlas = recover_atlas(&net, &[ch_a, ch_bpp], n);
        let ok = la
            && lb
            && matches!(
                atlas,
                AtlasReading::Abstain(AtlasAbstainReason::GlueInconsistent)
            );
        if !ok {
            bad.push(format!(
                "局所 A = {} / B'' = {} / atlas = {}",
                la,
                lb,
                atlas.as_str()
            ));
        }
        check(
            "[A3] cocycle 不整合 — entangler W = CZ₂₃ で捻れた chart は matching が破れ Abstain(GlueInconsistent)",
            bad.is_empty(),
            "A の M₂(2) と B'' の W M₂(2) W† は overlap 1/3 — 両 chart とも局所 Exact なのに大域は無い (chart 局所 Exact ↛ 大域 = 禁止変換 19)".into(),
        );
    }

    // ---- [A4] 複数 glue → EquivalenceClassOnly ----
    {
        let mut bad = Vec::new();
        let v = dft8();
        let mode: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&v, g, n)).collect();
        let mut gens = site.clone();
        gens.extend(mode.iter().cloned());
        // 文脈: site 側 (A 内・B 内・bridge) + DFT 側 ({VX_iV†}₃, {VZ_iV†}₃)
        let ctxs: Vec<Vec<usize>> = vec![
            vec![0, 2],
            vec![1, 3],
            vec![2, 4],
            vec![3, 5],
            vec![0, 4],
            vec![6, 8, 10],
            vec![7, 9, 11],
        ];
        let (net, ids) = build_net(&gens, n, tau, &ctxs);
        let s1 = chart(&ids, &[0, 1, 2, 3]);
        let s2 = chart(&ids, &[2, 3, 4, 5]);
        let d = chart(&ids, &[6, 7, 8, 9, 10, 11]);
        let atlas = recover_atlas(&net, &[s1, s2, d], n);
        let ok = matches!(
            atlas,
            AtlasReading::EquivalenceClassOnly {
                n_consistent_atlases: 2
            }
        );
        if !ok {
            bad.push(format!("atlas 裁定が {}", atlas.as_str()));
        }
        check(
            "[A4] 複数 glue — site atlas (S1+S2) と DFT chart が各々全被覆で内部整合・相互 matching 破れ → EquivalenceClassOnly{2}",
            bad.is_empty(),
            "整合 atlas が 2 つ — 無制約 tie-break で 1 つを選ばない (v32.3 [F3] site×DFT の atlas 版)".into(),
        );
    }

    // ---- [A5] witness 境界の両 lane 一致 ----
    {
        let mut bad = Vec::new();
        // [A1] から bridge {X₁,X₃} を外す
        let ctxs: Vec<Vec<usize>> = vec![vec![0, 2], vec![1, 3], vec![2, 4], vec![3, 5]];
        let (net, ids) = build_net(&site, n, tau, &ctxs);
        let ch_a = chart(&ids, &[0, 1, 2, 3]);
        let ch_b = chart(&ids, &[2, 3, 4, 5]);
        let atlas = recover_atlas(&net, &[ch_a, ch_b], n);
        let glue_abst = matches!(
            atlas,
            AtlasReading::Abstain(AtlasAbstainReason::CompatibilityUnwitnessed)
        );
        let direct = net.recovery_input().map(|i| i.recover().reading);
        let direct_abst = direct
            == Ok(FactorizationReading::Abstain(
                FactorizationAbstainReason::OperationalCompatibilityUnwitnessed,
            ));
        if !(glue_abst && direct_abst) {
            bad.push(format!(
                "glue = {} / 直接 = {:?}",
                atlas.as_str(),
                direct.map(|r| r.as_str().to_string())
            ));
        }
        check(
            "[A5] witness 境界 — bridge 文脈を外すと glue も直接大域復元も unwitnessed で棄却 (両 lane 一致)",
            bad.is_empty(),
            "v33.1 の証人規律は atlas でも緩まない — 因子対 (1,3) の共同 addressability の証人がない".into(),
        );
    }

    // ---- [A6] 封鎖の schema/文書検査 ----
    {
        let mut bad = Vec::new();
        // 禁止 impl From (分割リテラル規約 — v332 [C8] と同じ)
        let forbidden_impls: [String; 2] = [
            format!("impl From{}", "<ChartLocalFactorization"),
            format!("impl From{}", "<AtlasReading"),
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
        let src = rd("sim/src/contextual_factorization.rs").unwrap_or_default();
        for needle in [
            "禁止変換 19",
            "ChartLocalFactorization",
            "GlueInconsistent",
            "CompatibilityUnwitnessed",
            "黙って捨てない",
        ] {
            if !src.contains(needle) {
                bad.push(format!("contextual_factorization.rs: 「{}」が無い", needle));
            }
        }
        let schema = rd("core.schema.yml").unwrap_or_default();
        for needle in [
            "- name: ChartSpec",
            "- name: ChartLocalFactorization",
            "- name: AtlasReading",
        ] {
            if !schema.contains(needle) {
                bad.push(format!("core.schema.yml: 「{}」が無い", needle));
            }
        }
        if !schema.contains("- from: ChartLocalFactorization\n  to: GlobalFactorization\n  reason:")
        {
            bad.push("禁止変換 19 が未登録".into());
        }
        let doc = rd("docs/uft-v33.4.md").unwrap_or_default();
        for needle in [
            "Contextual factorization",
            "chart",
            "overlap",
            "glue 定理",
            "禁止変換 19",
            "GluedExact",
            "EquivalenceClassOnly",
            "cocycle",
            "witness 境界",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v33.4.md: 「{}」が無い", needle));
            }
        }
        check(
            "[A6] 封鎖の schema/文書 — 概念登録 + 禁止変換 19・impl From 不在・アンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "chart の局所 Exact から大域への直接路はコンパイル不能".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "局所性は chart atlas の glue として識別されるようになった — 整合すれば大域復元と一致し、捻れれば棄却し、複数の整合 atlas は選ばない。graded lane (Majorana/Dirac) が v33.5 の主題"
        } else {
            "**atlas 契約の破れ** — contextual_factorization と schema/文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
