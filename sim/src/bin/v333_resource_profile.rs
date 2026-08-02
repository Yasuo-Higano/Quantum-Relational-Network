//! v33.3 Resource-Filtered OperationalNet — poset 上の factorization profile (PROMPT/14)
//!
//! 可アクセス性は resource budget に依存する。ゆえに v33.3 の中心対象は「単一
//! 因子分解」ではなく **profile**: ResourceBudget ↦ FactorizationReading。
//! 共有契約 = `sim/src/resource_profile.rs`。
//!
//!   [R0] 契約自己検査 — resource_profile_self_test (新) + 既存 4 契約の不変
//!   [R1] **budget chain**: 資源不足 → NoAccessibleOperations・局所制御 →
//!        Exact [2,2,2]・entangler 可能 → [2,4]・完全 → [8] (深さ軸の鎖で機械化)。
//!        資格のない操作は budget をいくら積んでも現れない (禁止変換 14 の維持)
//!   [R2] **poset は barcode ではない**: 比較不能な budget 対 (amp 2, depth 1) と
//!        (amp 1, depth 2) が**同じ dims [2,4] で別の gauge orbit** ({12|3} vs
//!        {1|23} — matching 不在) を持つ。分裂・併合が一次元の出生死滅で書けない —
//!        有限 poset 上の constructible profile として扱う (zigzag 等への昇格は
//!        写像と安定性定理の後)
//!   [R3] **昇格規則 (禁止変換 17)**: stable = 領域に比較可能な対 (chain ≥ 2) が
//!        存在。深さ鎖 {1, 1.5, 2, 3} で [2,2,2] は {1, 1.5} の chain で stable・
//!        [2,4]@2 と [8]@3 は単点で transient (grid 相対の正直な記録) —
//!        transient_factorization_promotions = 0
//!   [R4] **command 再パラメータ化不変性 + スカラー潰しの負制御 (禁止変換 18)**:
//!        成分ごとの狭義単調再パラメータ化 (amp×3, depth²) で profile は点ごとに
//!        不変。恣意的な重み付き和 (amp + depth) への全順序化は accessibility 集合を
//!        変え、(2,1) の読みを [2,4] から [8] へ反転させる — 潰しは不変でない
//!   [R5] **頂は経路を消す (erasure 対照)**: コスト構造だけが異なる 2 interface
//!        (site-first vs entangler-first) は top budget で同じ [8] に合流するが、
//!        低予算の読み ([2,2,2] vs Abstain) が異なる — 最終 budget だけを見ると
//!        v32.2 の erasure no-go に戻る。profile が経路の情報を運ぶ
//!   [R6] 封鎖の schema/文書検査 — 概念登録 + 禁止変換 17/18・uft-v33.3.md アンカー
//!
//! 実行: cargo run --release --bin v333_resource_profile

use std::fs;
use std::path::Path;
use uft_sim::laboratory_interface::{
    certify_addressability, AccessibleOperation, IndependentAddressabilityCertificate,
    OperationOrigin, ResourceBudget,
};
use uft_sim::operational_net::{
    same_gauge_orbit, ControlGenerator, FactorizationReading, OpKind, OperatorParity,
    OrdinaryCommutation,
};
use uft_sim::resource_profile::*;
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

/// budget/コスト: (amplitude, depth) だけを動かし、他 3 成分は固定
fn bud(amp: f64, depth: f64) -> ResourceBudget {
    ResourceBudget::certify(1.0, amp, 1.0, depth, 1e-9).unwrap()
}

/// interface の組立て: ops = (行列, コスト) の列 (1:1 直接較正)・recipes = 文脈レシピ
fn build_interface(
    ops: &[(Vec<C64>, ResourceBudget)],
    n: usize,
    tau: f64,
    recipes: &[Vec<usize>],
) -> (
    ResourceFilteredInterface<OrdinaryCommutation>,
    IndependentAddressabilityCertificate,
) {
    let mats: Vec<Vec<C64>> = ops.iter().map(|(m, _)| m.clone()).collect();
    let cert = certify_addressability(&mats, &mats, n, 0.5, 0.1).expect("直接較正が立たない");
    let mut iface: ResourceFilteredInterface<OrdinaryCommutation> =
        ResourceFilteredInterface::new(n, tau);
    for (m, cost) in ops {
        let op = AccessibleOperation::certify(
            OpKind::Control(
                ControlGenerator::certify(
                    m.iter().map(|c| c.re).collect(),
                    m.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            OperatorParity::Even,
            OperationOrigin::DirectlyCalibrated(cert.clone()),
            cert.clone(),
            *cost,
        )
        .expect("資格つき操作の構成が失敗");
        iface.add_operation(op);
    }
    for r in recipes {
        iface
            .add_context_recipe(r.clone(), cert.clone())
            .expect("文脈レシピの登録が失敗");
    }
    (iface, cert)
}

fn reading_of(p: &ProfilePoint) -> String {
    match p {
        ProfilePoint::Reading { reading, .. } => match reading {
            FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims } => {
                format!("exact {:?}", local_dims)
            }
            r => r.as_str().to_string(),
        },
        p => p.as_str().to_string(),
    }
}

fn main() {
    uft_sim::self_test();
    println!("=== v33.3 Resource-Filtered OperationalNet — poset 上の factorization profile (PROMPT/14) ===\n");
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

    // 共通素材: site 6 本 + entangler 2 本 ((C²)⊗³)
    let site: Vec<Vec<C64>> = ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"]
        .iter()
        .map(|s| op3(s))
        .collect();
    let e12 = op3("XXI");
    let e23 = op3("IXX");
    // 文脈レシピ (ops の並び: 0..5 = site, 6 = E12, 7 = E23):
    // X 文脈 {X₁, X₂, X₃, E12, E23} (全対可換)・Z 文脈 {Z₁, Z₂, Z₃}
    let recipes: Vec<Vec<usize>> = vec![vec![0, 2, 4, 6, 7], vec![1, 3, 5]];
    let mk_ops = |c_site: ResourceBudget,
                  c_e12: ResourceBudget,
                  c_e23: ResourceBudget|
     -> Vec<(Vec<C64>, ResourceBudget)> {
        let mut v: Vec<(Vec<C64>, ResourceBudget)> =
            site.iter().map(|m| (m.clone(), c_site)).collect();
        v.push((e12.clone(), c_e12));
        v.push((e23.clone(), c_e23));
        v
    };

    // ---- [R0] 契約自己検査 ----
    {
        let mut bad = Vec::new();
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
        if let Err(e) = uft_sim::qrn_core::qrn_core_self_test() {
            bad.push(format!("qrn_core_self_test: {}", e));
        }
        check(
            "[R0] 契約自己検査 — resource_profile (新) + laboratory_interface/operational_net/scope_repair/qrn_core 不変",
            bad.is_empty(),
            if bad.is_empty() {
                "profile は poset 上の constructible データ・昇格規則は chain ≥ 2".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [R1] budget chain: 資源不足 → 局所 → entangler → global ----
    {
        let mut bad = Vec::new();
        let (iface, _) = build_interface(
            &mk_ops(bud(1.0, 1.0), bud(1.0, 2.0), bud(1.0, 3.0)),
            n,
            tau,
            &recipes,
        );
        let grid = [bud(1.0, 0.5), bud(1.0, 1.0), bud(1.0, 2.0), bud(1.0, 3.0)];
        let prof = iface.profile_over(&grid);
        let got: Vec<String> = prof.points.iter().map(reading_of).collect();
        let want = [
            "no_accessible_operations",
            "exact [2, 2, 2]",
            "exact [2, 4]",
            "exact [8]",
        ];
        if got != want {
            bad.push(format!("chain の読みが {:?} (期待 {:?})", got, want));
        }
        check(
            "[R1] budget chain — 資源不足 → 局所 [2,2,2] → entangler [2,4] → global [8] (深さ鎖)",
            bad.is_empty(),
            format!(
                "depth 0.5 → {} / 1 → {} / 2 → {} / 3 → {} — 読みは budget の関数であり単一の「正解因子分解」は存在しない",
                got[0], got[1], got[2], got[3]
            ),
        );
    }

    // ---- [R2] poset は barcode ではない ----
    {
        let mut bad = Vec::new();
        let (iface, _) = build_interface(
            &mk_ops(bud(1.0, 1.0), bud(2.0, 1.0), bud(1.0, 2.0)),
            n,
            tau,
            &recipes,
        );
        let b_c = bud(2.0, 1.0);
        let b_d = bud(1.0, 2.0);
        if b_c.comparable(&b_d) {
            bad.push("対が比較可能になっている (poset の設計ミス)".into());
        }
        let p_c = iface.reading_at(&b_c);
        let p_d = iface.reading_at(&b_d);
        let extract = |p: &ProfilePoint| -> Option<(Vec<usize>, Vec<Vec<Vec<C64>>>)> {
            match p {
                ProfilePoint::Reading {
                    reading:
                        FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims },
                    component_subalgebras,
                } => Some((local_dims.clone(), component_subalgebras.clone())),
                _ => None,
            }
        };
        let mut ov = 0.0f64;
        match (extract(&p_c), extract(&p_d)) {
            (Some((dims_c, sub_c)), Some((dims_d, sub_d))) => {
                let (same, best) = same_gauge_orbit(&sub_c, &sub_d);
                ov = best;
                if !(dims_c == vec![2, 4]
                    && dims_d == vec![2, 4]
                    && !same
                    && !p_c.same_class(&p_d))
                {
                    bad.push(format!(
                        "dims {:?}/{:?}・orbit matching {} (期待: 同 dims・別 orbit)",
                        dims_c, dims_d, same
                    ));
                }
            }
            _ => bad.push(format!(
                "(2,1)/(1,2) の読みが {}/{} (期待: 両方 Exact)",
                reading_of(&p_c),
                reading_of(&p_d)
            )),
        }
        let p_top = iface.reading_at(&bud(2.0, 2.0));
        if reading_of(&p_top) != "exact [8]" {
            bad.push(format!("(2,2) の読みが {}", reading_of(&p_top)));
        }
        check(
            "[R2] poset は barcode ではない — 比較不能対 (2,1)/(1,2) が同 dims [2,4] で別 gauge orbit",
            bad.is_empty(),
            format!(
                "orbit {{12|3}} vs {{1|23}}: matching 不在 (最良 min-overlap {:.4})・join (2,2) → [8] — 分裂/併合は 1 次元の出生死滅で書けず、profile は poset 上の constructible データ",
                ov
            ),
        );
    }

    // ---- [R3] 昇格規則 (禁止変換 17) ----
    {
        let mut bad = Vec::new();
        let (iface, _) = build_interface(
            &mk_ops(bud(1.0, 1.0), bud(1.0, 2.0), bud(1.0, 3.0)),
            n,
            tau,
            &recipes,
        );
        let grid = [bud(1.0, 1.0), bud(1.0, 1.5), bud(1.0, 2.0), bud(1.0, 3.0)];
        let prof = iface.profile_over(&grid);
        let classes = prof.classes();
        // [2,2,2] = {depth 1, 1.5} chain → stable / [2,4]@2, [8]@3 → transient
        let mut stable_readings = Vec::new();
        let mut transient_readings = Vec::new();
        for c in &classes {
            let r = reading_of(&prof.points[c.representative]);
            if c.stable {
                stable_readings.push((r, c.region.len()));
            } else {
                transient_readings.push((r, c.region.len()));
            }
        }
        stable_readings.sort();
        transient_readings.sort();
        let ok = stable_readings == vec![("exact [2, 2, 2]".to_string(), 2)]
            && transient_readings
                == vec![
                    ("exact [2, 4]".to_string(), 1),
                    ("exact [8]".to_string(), 1),
                ];
        if !ok {
            bad.push(format!(
                "stable {:?} / transient {:?}",
                stable_readings, transient_readings
            ));
        }
        // transient の昇格 0 (HOLD-9 の採点語彙)
        let promoted: Vec<&ProfileClass> = classes.iter().filter(|c| c.stable).collect();
        let transient_promotions = promoted
            .iter()
            .filter(|c| c.region.len() < 2)
            .count();
        if transient_promotions != 0 {
            bad.push(format!("transient の昇格 {}", transient_promotions));
        }
        check(
            "[R3] 昇格規則 (禁止変換 17) — stable は chain ≥ 2 の領域のみ・単点の読みは transient (昇格 0)",
            bad.is_empty(),
            format!(
                "stable = {:?} / transient = {:?} — 単一閾値で一瞬だけ出現した因子分解は局所性に昇格しない (grid の頂の [8] も調査 grid 相対で transient と正直に記録)",
                stable_readings, transient_readings
            ),
        );
    }

    // ---- [R4] 再パラメータ化不変性 + スカラー潰しの負制御 (禁止変換 18) ----
    {
        let mut bad = Vec::new();
        // 元の [R2] interface と grid
        let costs = [
            (bud(1.0, 1.0), bud(2.0, 1.0), bud(1.0, 2.0)),
        ][0];
        let (iface, _) = build_interface(&mk_ops(costs.0, costs.1, costs.2), n, tau, &recipes);
        let grid = [bud(1.0, 1.0), bud(2.0, 1.0), bud(1.0, 2.0), bud(2.0, 2.0)];
        let prof = iface.profile_over(&grid);
        // 成分ごとの狭義単調再パラメータ化 φ: amp → 3·amp, depth → depth²
        let phi = |b: &ResourceBudget| bud(3.0 * b.max_amplitude, b.max_depth * b.max_depth);
        let (iface2, _) = build_interface(
            &mk_ops(phi(&costs.0), phi(&costs.1), phi(&costs.2)),
            n,
            tau,
            &recipes,
        );
        let grid2: Vec<ResourceBudget> = grid.iter().map(|b| phi(b)).collect();
        let prof2 = iface2.profile_over(&grid2);
        let same_pointwise = prof
            .points
            .iter()
            .zip(prof2.points.iter())
            .all(|(a, b)| a.same_class(b));
        if !same_pointwise {
            bad.push("狭義単調再パラメータ化で profile が変わった".into());
        }
        // スカラー潰し (amp + depth): (2,1) の読みが [2,4] → [8] に反転
        let scalar_cost = |b: &ResourceBudget| b.max_amplitude + b.max_depth;
        let costs_all = [costs.0, costs.0, costs.0, costs.0, costs.0, costs.0, costs.1, costs.2];
        let b_c = bud(2.0, 1.0);
        let budget_scalar = scalar_cost(&b_c); // = 3
        let accessible_scalar: Vec<usize> = (0..8)
            .filter(|&i| scalar_cost(&costs_all[i]) <= budget_scalar)
            .collect();
        // 成分半順序では E23 (1,2) は (2,1) に不可 — スカラーでは 1+2 = 3 ≤ 3 で可
        let flip = accessible_scalar.contains(&7);
        let componentwise = iface.reading_at(&b_c);
        if !(flip && reading_of(&componentwise) == "exact [2, 4]") {
            bad.push(format!(
                "スカラー潰しの反転が立たない (flip = {}, 成分 = {})",
                flip,
                reading_of(&componentwise)
            ));
        }
        // スカラー化した読み: 全 8 op accessible → [8]
        let (iface_s, _) = build_interface(
            &mk_ops(bud(1.0, 1.0), bud(1.0, 1.0), bud(1.0, 1.0)),
            n,
            tau,
            &recipes,
        );
        let p_s = iface_s.reading_at(&bud(1.0, 1.0));
        if reading_of(&p_s) != "exact [8]" {
            bad.push(format!("スカラー同値の読みが {}", reading_of(&p_s)));
        }
        check(
            "[R4] 再パラメータ化不変 (成分ごと狭義単調) + スカラー潰しは裁定を反転 (禁止変換 18)",
            bad.is_empty(),
            format!(
                "φ = (amp×3, depth²) で profile 点ごと不変 = {} / スカラー amp+depth は (2,1) で E23 を accessible にし読みが [2,4] → [8] に反転 — 恣意的な重み付き和は新しい選択バイアス",
                same_pointwise
            ),
        );
    }

    // ---- [R5] 頂は経路を消す (erasure 対照) ----
    {
        let mut bad = Vec::new();
        // P: site-first (site 安い)・Q: entangler-first (entangler 安い)
        let (iface_p, _) = build_interface(
            &mk_ops(bud(1.0, 1.0), bud(2.0, 2.0), bud(2.0, 2.0)),
            n,
            tau,
            &recipes,
        );
        let (iface_q, _) = build_interface(
            &mk_ops(bud(2.0, 2.0), bud(1.0, 1.0), bud(1.0, 1.0)),
            n,
            tau,
            &recipes,
        );
        let lo = bud(1.0, 1.0);
        let hi = bud(2.0, 2.0);
        let p_lo = iface_p.reading_at(&lo);
        let q_lo = iface_q.reading_at(&lo);
        let p_hi = iface_p.reading_at(&hi);
        let q_hi = iface_q.reading_at(&hi);
        let top_same = p_hi.same_class(&q_hi) && reading_of(&p_hi) == "exact [8]";
        let path_differs = !p_lo.same_class(&q_lo)
            && reading_of(&p_lo) == "exact [2, 2, 2]"
            && reading_of(&q_lo) == "abstain";
        if !top_same {
            bad.push(format!("頂が一致しない ({} vs {})", reading_of(&p_hi), reading_of(&q_hi)));
        }
        if !path_differs {
            bad.push(format!(
                "経路が区別されない (P@lo = {} / Q@lo = {})",
                reading_of(&p_lo),
                reading_of(&q_lo)
            ));
        }
        check(
            "[R5] 頂は経路を消す — top budget は両 interface とも [8] に合流・低予算の読みが profile を区別",
            bad.is_empty(),
            format!(
                "P (site-first): lo → {} / Q (entangler-first): lo → {} / 両者 hi → exact [8] — 最終 budget だけを見ると v32.2 erasure no-go に戻る。profile が経路 (どの操作から局所性が組み上がったか) を運ぶ",
                reading_of(&p_lo),
                reading_of(&q_lo)
            ),
        );
    }

    // ---- [R6] 封鎖の schema/文書検査 ----
    {
        let mut bad = Vec::new();
        let src = rd("sim/src/resource_profile.rs").unwrap_or_default();
        for needle in [
            "禁止変換 17",
            "禁止変換 18",
            "barcode",
            "chain ≥ 2",
            "NoAccessibleOperations",
            "constructible",
        ] {
            if !src.contains(needle) {
                bad.push(format!("resource_profile.rs: 「{}」が無い", needle));
            }
        }
        let schema = rd("core.schema.yml").unwrap_or_default();
        for needle in [
            "- name: ResourceFilteredInterface",
            "- name: OperationalFactorizationProfile",
            "- name: ProfilePoint",
            "- name: StableFactorizationRegion",
        ] {
            if !schema.contains(needle) {
                bad.push(format!("core.schema.yml: 「{}」が無い", needle));
            }
        }
        if !schema
            .contains("- from: TransientFactorizationPoint\n  to: StableFactorizationRegion\n  reason:")
        {
            bad.push("禁止変換 17 が未登録".into());
        }
        if !schema.contains("- from: ResourceBudget\n  to: ScalarResourceCost\n  reason:") {
            bad.push("禁止変換 18 が未登録".into());
        }
        let doc = rd("docs/uft-v33.3.md").unwrap_or_default();
        for needle in [
            "Resource-Filtered",
            "profile",
            "barcode",
            "禁止変換 17",
            "禁止変換 18",
            "chain ≥ 2",
            "transient",
            "頂は経路を消す",
            "再パラメータ化",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v33.3.md: 「{}」が無い", needle));
            }
        }
        check(
            "[R6] 封鎖の schema/文書 — 概念登録 + 禁止変換 17/18・アンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "profile・昇格規則・スカラー潰しの禁止が schema/文書の三点で凍結された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "局所性の読みは budget の関数になった — 昇格できるのは perturbation に生き残る安定同値類だけで、単点の読みと頂の読みは経路を語らない。文脈の overlap 整合 (glue) が v33.4 の主題"
        } else {
            "**profile 契約の破れ** — resource_profile と schema/文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
