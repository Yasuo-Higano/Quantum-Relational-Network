//! v27.2 型付き Core Contract の機械検査 (PROMPT/10 §4)
//!
//! qrn_core (sim/src/qrn_core.rs) が定義する型境界・bridge law の門・状態表示が、
//! core.schema.yml / docs の凍結内容と整合し、かつ「水増し」がないことを検査する:
//!   [T0] qrn_core_self_test — bridge law 登録簿が空・証明書発行不能・
//!        dynamics/bridge に Defined の混入なし・昇格先が居住不能型 (サイズ 0)
//!   [T1] QRN_CORE_V1 (Rust const) ↔ core.schema.yml の status_* の一致
//!   [T2] 旧名の根絶 — sim/src に QrnState / QrnModel トークンが残存しない
//!        (GaussianFermionState / GaussianToyModel / ConstrainedToy* へ改名済み)
//!   [T3] 外部時間の型分離 — GaussianToyModel::evolve と ConstrainedToyDynamicsV2::step
//!        の引数が EvolutionParameter 型 (f64 の裸渡しは型エラー)
//!   [T4] 証明書の門 — 未登録 claim id への register は全て None
//!   [T5] Lean 形式化 (proofs/QrnPromotion.lean) の存在と定理数
//!
//! 「Lean 証明済み」が意味するのは昇格グラフ・状態表示の形式的性質であって、
//! 物理的正しさではない (PROMPT/10 §4 の明示要件)。

use std::fs;
use std::path::Path;
use uft_sim::qrn_core::*;

fn main() {
    uft_sim::self_test();
    println!("=== v27.2 型付き Core Contract の機械検査 (PROMPT/10 §4) ===\n");
    let root = if Path::new("core.schema.yml").exists() {
        "."
    } else {
        ".."
    };
    let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p));
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

    // ---- [T0] qrn_core の不変条件 ----
    {
        let r = qrn_core_self_test();
        check(
            "[T0] qrn_core_self_test — 登録簿空・証明書発行不能・状態表示の水増しなし・居住不能型",
            r.is_ok(),
            match &r {
                Ok(()) => format!(
                    "REGISTERED_BRIDGE_LAWS = {:?} / QrnEvidence 等 4 型は構成不能 (空 enum)",
                    REGISTERED_BRIDGE_LAWS
                ),
                Err(e) => e.clone(),
            },
        );
    }

    // ---- [T1] QRN_CORE_V1 ↔ core.schema.yml ----
    {
        let schema = rd("core.schema.yml").unwrap_or_default();
        let scalar = |key: &str| -> String {
            schema
                .lines()
                .find(|l| l.starts_with(&format!("{}:", key)))
                .and_then(|l| l.split(':').nth(1))
                .map(|v| v.trim().trim_matches('"').to_string())
                .unwrap_or_default()
        };
        let all_layer = |layer: &str, want: ContractStatus| -> bool {
            QRN_CORE_V1
                .iter()
                .filter(|c| c.layer == layer)
                .all(|c| c.status == want)
        };
        let comp = |name: &str| QRN_CORE_V1.iter().find(|c| c.name == name).unwrap();
        let mut bad = Vec::new();
        if !(scalar("status_kinematics") == "defined" && all_layer("core", ContractStatus::Defined))
        {
            bad.push("kinematics");
        }
        if !(scalar("status_dynamics") == "model_family_only"
            && all_layer("dynamics", ContractStatus::ModelFamilyOnly))
        {
            bad.push("dynamics");
        }
        if !(scalar("status_geometry_bridge") == "conjectural"
            && comp("GeometryBridge").status == ContractStatus::Conjectural)
        {
            bad.push("geometry_bridge");
        }
        if !(scalar("status_gravity_bridge") == "unsupported"
            && comp("GravityBridge").status == ContractStatus::Unsupported)
        {
            bad.push("gravity_bridge");
        }
        if scalar("status_empirical_prediction") != "none" {
            bad.push("empirical_prediction");
        }
        check(
            "[T1] QRN_CORE_V1 (Rust const 10 成分) ↔ core.schema.yml の状態表示",
            bad.is_empty(),
            if bad.is_empty() {
                format!(
                    "kinematics: defined / dynamics: model_family_only / geometry: {} / gravity: {} / clock: {}",
                    comp("GeometryBridge").status.as_str(),
                    comp("GravityBridge").status.as_str(),
                    comp("ClockBridge").status.as_str()
                )
            } else {
                format!("不一致: {:?}", bad)
            },
        );
    }

    // ---- [T2] 旧名の根絶 ----
    {
        // 本バイナリ自身は検査パターンを文字列として含むため対象外 (v271 [S12] と同型)
        const EXEMPT: [&str; 1] = ["v272_core_contract.rs"];
        const OLD_NAMES: [&str; 2] = ["QrnState", "QrnModel"];
        let mut hits = Vec::new();
        for dir in ["sim/src", "sim/src/bin"] {
            if let Ok(entries) = fs::read_dir(format!("{}/{}", root, dir)) {
                for e in entries.filter_map(|e| e.ok()) {
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
                        // コメント行 (// /// //!) は改名台帳への歴史的言及を許す —
                        // 検査対象は識別子としての使用 (コード行) のみ
                        for (lno, line) in text.lines().enumerate() {
                            if line.trim_start().starts_with("//") {
                                continue;
                            }
                            for o in OLD_NAMES {
                                if line.contains(o) {
                                    hits.push(format!("{}:{}: {}", name, lno + 1, o));
                                }
                            }
                        }
                    }
                }
            }
        }
        check(
            "[T2] 旧名 QrnState/QrnModel の根絶 (GaussianFermionState/GaussianToyModel へ)",
            hits.is_empty(),
            if hits.is_empty() {
                "sim/src 全走査 — 残存 0 (toy の名から Qrn を外した)".into()
            } else {
                format!("{:?}", hits)
            },
        );
    }

    // ---- [T3] 外部時間の型分離 (署名の source 検査) ----
    {
        let lib = rd("sim/src/lib.rs").unwrap_or_default();
        let ok1 = lib
            .contains("fn evolve(&self, s: &GaussianFermionState, t: EvolutionParameter) -> GaussianFermionState");
        let ok2 = lib.contains("fn step(&self, s: &S, dt: EvolutionParameter) -> S");
        // 型そのものの実演: EvolutionParameter は ProperTime と別型。v30.0 以降、
        // ProperTime へ至る能力 (ClockCalibration) は BridgeCapability 未実装で
        // 証明書型が構成不能 — 門となる関数も削除済み (qrn_core の source 検査)
        let t = EvolutionParameter(1.5);
        let core_src = rd("sim/src/qrn_core.rs").unwrap_or_default();
        let ok3 = !core_src.contains("fn promote_evolution_to_proper_time(")
            && core_src.contains("impl sealed_cap::Sealed for ClockCalibration {}")
            && !core_src.contains("impl_capability!(ClockCalibration");
        check(
            "[T3] 外部時間の型分離 — evolve/step の t は EvolutionParameter (ProperTime への門は関数ごと不在)",
            ok1 && ok2 && ok3,
            format!(
                "GaussianToyModel::evolve {} / ConstrainedToyDynamicsV2::step {} / 門は閉 (t = {:?})",
                ok1, ok2, t
            ),
        );
    }

    // ---- [T4] 証明書の門 (v30.0: 能力別) ----
    {
        let ids = [
            "QRN-GRAV-001",
            "QRN-META-029",
            "QRN-CORE-001",
            "QRN-BRIDGE-004",
            "MutualInformationGeometry",
        ];
        fn locked<C: BridgeCapability>(ids: &[&'static str]) -> bool {
            C::REGISTERED.is_empty()
                && ids
                    .iter()
                    .all(|id| BridgeLawCertificate::<C>::register(id).is_none())
        }
        let all_none = locked::<FactorizationGivenObservables>(&ids)
            && locked::<SpatialTopologyGivenFactorization>(&ids)
            && locked::<SpatialMetricUpToGlobalScale>(&ids)
            && locked::<CausalOrderGivenExternalClock>(&ids)
            && locked::<ConformalLorentzianStructure>(&ids)
            && locked::<VolumeMeasure>(&ids);
        check(
            "[T4] BridgeLawCertificate<能力>::register — 全 6 実装能力 × 未登録 id で None (bridge law 0 件・ClockCalibration/FullLorentzianMetric は型レベル構成不能)",
            all_none,
            format!("{} 個の id × 6 能力で発行拒否を確認", ids.len()),
        );
    }

    // ---- [T5] Lean 形式化の存在 ----
    {
        let lean = rd("proofs/QrnPromotion.lean").unwrap_or_default();
        let n_thm = lean
            .lines()
            .filter(|l| l.trim_start().starts_with("theorem "))
            .count();
        check(
            "[T5] proofs/QrnPromotion.lean — 昇格グラフ・状態表示の形式化 (定理 12 本)",
            n_thm == 12,
            format!(
                "theorem 宣言 {} 本 (証明対象は形式的性質であり物理的正しさではない)",
                n_thm
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "型境界は成立 — regulator/存在論・外部/創発計量・発展/固有時間の暗黙変換はコンパイル不能、bridge law の門は閉"
        } else {
            "**契約の破れ** — qrn_core と schema/docs の整合を修正せよ"
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
