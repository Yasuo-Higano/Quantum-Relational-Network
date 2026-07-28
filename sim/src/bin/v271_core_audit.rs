//! v27.1 核意味論監査 (Core Semantic Audit) — 第二十八期の最初の道具 (PROMPT/10 §2–§3)
//!
//! 「何が QRN 本体で、何が既知 QFT の測定器か」をコード上も文書上も混同不能にする:
//!   [S0] core.schema.yml の書式・層語彙・6 軸語彙・QRN-Core v1 状態表示の解析
//!   [S1] 全概念の層分類 — 未分類 (unclassified) 0 件・必須概念の登録 (ゲート:
//!        未分類が 1 件でもあれば v27.2 に進まない)
//!   [S2] 型レベル禁止変換 7 種の登録 (RegulatorSite→RelationalNode 等)
//!   [S3] claims.yml 全主張の 6 軸完備・語彙適合 (v61_ledger [7] と独立に再検査)
//!   [S4–S9] 昇格禁止規則 R1–R7 の機械検査:
//!        R1/R2 natural_observation・external_replication = 0 件
//!        R3 independent_author は C0 (外部の確立結果) のみ
//!        R4 core/dynamics/bridge 層 (C0 以外) は physical_scope ∈ {toy, effective_model}
//!        R5 data_relation = future_observation の主張は存在しない
//!        R6 continuum_universal は adapter/instrument 層 (または C0) の専有
//!        R7 internal_holdout ⇒ preregistered_holdout
//!   [S10] predictions.yml との照合 — 全 PRED が分類済み・future_observation の PRED は
//!        全て未採点 (registered) ⇒「自然の観測量の的中 0」の機械化・PRED-019 不在
//!   [S11] 文書アンカー — 状態表示 (spec §7)・監査注記 (uft-v0.7/v1.0/v27.0)・
//!        QRN-C0-001 の限定・QRN-GRAV-001 の降格・README の残高行
//!   [S12] 禁止 impl From の不在 (sim/src 全走査 — v27.2 の型実装後も常設)
//!
//! 本監査はスイートの常時実行層 (ALWAYS_RUN)。分類の意味論は ASM-LAYER-SEMANTICS
//! (規約) — 監査が保証するのは完備性・語彙適合・禁止昇格の不在であり、
//! 分類自体の物理的正しさではない。

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;

// ---------------------------------------------------------------- 共通

fn unquote(s: &str) -> String {
    let t = s.trim();
    if t.len() >= 2 && t.starts_with('"') && t.ends_with('"') {
        t[1..t.len() - 1].to_string()
    } else {
        t.to_string()
    }
}

fn split_kv(s: &str) -> Option<(String, String)> {
    let idx = s.find(':')?;
    Some((s[..idx].trim().to_string(), unquote(&s[idx + 1..])))
}

fn parse_inline_list(v: &str) -> Option<Vec<String>> {
    let t = v.trim();
    if !t.starts_with('[') || !t.ends_with(']') {
        return None;
    }
    let inner = &t[1..t.len() - 1];
    if inner.trim().is_empty() {
        return Some(Vec::new());
    }
    Some(inner.split(',').map(|x| x.trim().to_string()).collect())
}

// ---------------------------------------------------------------- schema

#[derive(Default)]
struct Schema {
    scalars: BTreeMap<String, String>,
    lists: BTreeMap<String, Vec<String>>,
    concepts: Vec<BTreeMap<String, String>>,
    conversions: Vec<BTreeMap<String, String>>,
    predictions: Vec<BTreeMap<String, String>>,
}

fn parse_schema(text: &str) -> Result<Schema, String> {
    let mut s = Schema::default();
    #[derive(PartialEq, Clone, Copy)]
    enum Sec {
        Top,
        Concepts,
        Conversions,
        Predictions,
    }
    let mut sec = Sec::Top;
    for (lno, raw) in text.lines().enumerate() {
        let lno = lno + 1;
        let line = raw.trim_end();
        if line.trim_start().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        if !line.starts_with(' ') && !line.starts_with('-') {
            // トップレベル
            let (k, v) = split_kv(line).ok_or(format!("{}行目: トップレベルの書式", lno))?;
            match k.as_str() {
                "concepts" => sec = Sec::Concepts,
                "forbidden_conversions" => sec = Sec::Conversions,
                "predictions" => sec = Sec::Predictions,
                _ => {
                    sec = Sec::Top;
                    if let Some(lst) = parse_inline_list(&v) {
                        s.lists.insert(k, lst);
                    } else {
                        s.scalars.insert(k, v);
                    }
                }
            }
        } else if let Some(rest) = line.strip_prefix("- ") {
            let (k, v) = split_kv(rest).ok_or(format!("{}行目: 項目先頭の書式", lno))?;
            let mut m = BTreeMap::new();
            m.insert(k, v);
            match sec {
                Sec::Concepts => s.concepts.push(m),
                Sec::Conversions => s.conversions.push(m),
                Sec::Predictions => s.predictions.push(m),
                Sec::Top => return Err(format!("{}行目: セクション外の項目", lno)),
            }
        } else if let Some(rest) = line.strip_prefix("  ") {
            let (k, v) = split_kv(rest).ok_or(format!("{}行目: フィールドの書式", lno))?;
            let cur = match sec {
                Sec::Concepts => s.concepts.last_mut(),
                Sec::Conversions => s.conversions.last_mut(),
                Sec::Predictions => s.predictions.last_mut(),
                Sec::Top => None,
            }
            .ok_or(format!("{}行目: 項目外のフィールド", lno))?;
            cur.insert(k, v);
        } else {
            return Err(format!("{}行目: 解釈できない行", lno));
        }
    }
    Ok(s)
}

// ---------------------------------------------------------------- claims

#[derive(Default)]
struct Claim {
    id: String,
    level: String,
    claim: String,
    status: String,
    axes: BTreeMap<String, String>,
}

fn parse_claims(text: &str) -> Result<Vec<Claim>, String> {
    let mut out: Vec<Claim> = Vec::new();
    let mut in_block = false; // evidence/inputs/limitations ブロック内
    for (lno, raw) in text.lines().enumerate() {
        let lno = lno + 1;
        let line = raw.trim_end();
        if line.trim_start().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        if let Some(rest) = line.strip_prefix("- id:") {
            out.push(Claim {
                id: unquote(rest),
                ..Default::default()
            });
            in_block = false;
        } else if line.starts_with("    ") {
            if !in_block {
                return Err(format!("{}行目: ブロック外の 4 スペース行", lno));
            }
        } else if let Some(rest) = line.strip_prefix("  ") {
            let cur = out.last_mut().ok_or(format!("{}行目: エントリ外", lno))?;
            match rest {
                "evidence:" | "inputs:" | "limitations:" => in_block = true,
                _ => {
                    in_block = false;
                    let (k, v) = split_kv(rest).ok_or(format!("{}行目: フィールド書式", lno))?;
                    match k.as_str() {
                        "level" => cur.level = v,
                        "claim" => cur.claim = v,
                        "status" => cur.status = v,
                        "version" => {}
                        "layer" | "evidence_kind" | "independence" | "universality"
                        | "data_relation" | "physical_scope" => {
                            cur.axes.insert(k, v);
                        }
                        _ => return Err(format!("{}行目: 未知フィールド '{}'", lno, k)),
                    }
                }
            }
        } else {
            return Err(format!("{}行目: 解釈できない行", lno));
        }
    }
    Ok(out)
}

/// predictions.yml — id と status だけ拾う
fn parse_predictions(text: &str) -> Vec<(String, String)> {
    let mut out: Vec<(String, String)> = Vec::new();
    for raw in text.lines() {
        let line = raw.trim_end();
        if let Some(rest) = line.strip_prefix("- id:") {
            out.push((unquote(rest), String::new()));
        } else if let Some(rest) = line.strip_prefix("  status:") {
            if let Some(last) = out.last_mut() {
                last.1 = unquote(rest);
            }
        }
    }
    out
}

// ---------------------------------------------------------------- main

fn main() {
    println!("=== v27.1 核意味論監査 — 層分類・多軸台帳・昇格禁止 (第二十八期, PROMPT/10) ===\n");
    let root = if Path::new("core.schema.yml").exists() {
        "."
    } else if Path::new("../core.schema.yml").exists() {
        ".."
    } else {
        println!("core.schema.yml が見つからない  [FAIL]");
        std::process::exit(1);
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

    const LAYERS: [&str; 7] = [
        "core",
        "dynamics",
        "adapter",
        "instrument",
        "bridge",
        "phenomenology",
        "meta",
    ];

    // ---- [S0] schema の解析と語彙 ----
    let schema_text = rd("core.schema.yml").expect("core.schema.yml が読めない");
    let schema = match parse_schema(&schema_text) {
        Ok(s) => s,
        Err(e) => {
            println!("  [FAIL] [S0] core.schema.yml: {}", e);
            std::process::exit(1);
        }
    };
    {
        let layers_ok = schema.lists.get("layers").map(|l| l.as_slice())
            == Some(
                LAYERS
                    .iter()
                    .map(|s| s.to_string())
                    .collect::<Vec<_>>()
                    .as_slice(),
            );
        let axes_ok = [
            "axis_evidence_kind",
            "axis_independence",
            "axis_universality",
            "axis_data_relation",
            "axis_physical_scope",
        ]
        .iter()
        .all(|k| schema.lists.contains_key(*k));
        let status_ok = [
            ("status_kinematics", "defined"),
            ("status_dynamics", "model_family_only"),
            ("status_geometry_bridge", "conjectural"),
            ("status_gravity_bridge", "unsupported"),
            ("status_empirical_prediction", "none"),
        ]
        .iter()
        .all(|(k, v)| schema.scalars.get(*k).map(|s| s.as_str()) == Some(*v));
        check(
            "[S0] schema 解析・層 7 値・軸語彙 5 本・QRN-Core v1 状態表示",
            layers_ok && axes_ok && status_ok,
            format!(
                "concepts {} / conversions {} / predictions {}",
                schema.concepts.len(),
                schema.conversions.len(),
                schema.predictions.len()
            ),
        );
    }

    // ---- [S1] 概念の層分類 (ゲート: 未分類 0 件) ----
    {
        const REQUIRED: [&str; 35] = [
            "StateSpace",
            "ObservableAlgebra",
            "RelationalDecomposition",
            "EvolutionLaw",
            "ConstraintAlgebra",
            "InitialConditionRule",
            "RelationalNodeId",
            "RegulatorSiteId",
            "ContinuumPoint",
            "ExternalMetricSource",
            "EmergentMetricCandidate",
            "EvolutionParameter",
            "ModularParameter",
            "OperationalClockReading",
            "ProperTime",
            "GaussianFermionState",
            "GaussianToyModel",
            "BondACoupling",
            "MidpointModulation",
            "BelinfanteImprovement",
            "ContactTermCompletion",
            "StaggeredDiscretization",
            "WilsonDiscretization",
            "NullCombinationLadder",
            "ShellIntegration",
            "MatsubaraWard",
            "DerivedAsymptoticExtrapolator",
            "AnalyticOracle",
            "MutualInformationGeometry",
            "ModularFlowGeometry",
            "FisherResponseBridge",
            "CommutatorCausalGeometry",
            "GeometryReadout",
            "CausalReadout",
            "GravityBridge",
        ];
        let mut bad = Vec::new();
        let mut hist: BTreeMap<String, usize> = BTreeMap::new();
        let mut names = Vec::new();
        for c in &schema.concepts {
            let name = c.get("name").cloned().unwrap_or_default();
            let layer = c.get("layer").cloned().unwrap_or_default();
            if name.is_empty() {
                bad.push("name 欠落".to_string());
            }
            if layer == "unclassified" || layer.is_empty() {
                bad.push(format!("{}: 未分類", name));
            } else if !LAYERS.contains(&layer.as_str()) {
                bad.push(format!("{}: 層語彙外 '{}'", name, layer));
            }
            if c.get("status").is_none() || c.get("note").is_none() {
                bad.push(format!("{}: status/note 欠落", name));
            }
            *hist.entry(layer).or_default() += 1;
            names.push(name);
        }
        for r in REQUIRED {
            if !names.iter().any(|n| n == r) {
                bad.push(format!("必須概念 {} が未登録", r));
            }
        }
        let hists: Vec<String> = hist.iter().map(|(k, v)| format!("{}:{}", k, v)).collect();
        check(
            "[S1] 全概念の層分類 — 未分類 0 件・必須 35 概念の登録 (v27.2 への進行ゲート)",
            bad.is_empty(),
            if bad.is_empty() {
                format!("{} 概念 [{}]", schema.concepts.len(), hists.join(" "))
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S2] 型レベル禁止変換の登録 ----
    {
        const CONV: [(&str, &str); 7] = [
            ("RegulatorSiteId", "RelationalNodeId"),
            ("EvolutionParameter", "ProperTime"),
            ("ExternalMetricSource", "EmergentMetricCandidate"),
            ("CalibrationEvidence", "QrnEvidence"),
            ("InternalPrediction", "NaturalObservation"),
            ("SameAuthorReplication", "IndependentReplication"),
            ("RegulatorQuantity", "UniversalQuantity"),
        ];
        let mut bad = Vec::new();
        for (f, t) in CONV {
            let found = schema.conversions.iter().any(|c| {
                c.get("from").map(|s| s.as_str()) == Some(f)
                    && c.get("to").map(|s| s.as_str()) == Some(t)
            });
            if !found {
                bad.push(format!("{} -> {}", f, t));
            }
        }
        for c in &schema.conversions {
            if c.get("reason").is_none() {
                bad.push(format!("{:?}: reason 欠落", c.get("from")));
            }
        }
        check(
            "[S2] 型レベル禁止変換 7 種の登録 (reason つき)",
            bad.is_empty(),
            if bad.is_empty() {
                "regulator/時間/計量/証拠の 4 部門".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- claims の解析 ----
    let claims_text = rd("claims.yml").expect("claims.yml が読めない");
    let claims = match parse_claims(&claims_text) {
        Ok(c) => c,
        Err(e) => {
            println!("  [FAIL] [S3] claims.yml: {}", e);
            std::process::exit(1);
        }
    };

    // ---- [S3] 6 軸の完備・語彙適合 ----
    {
        let vocab = |axis: &str| -> Vec<String> {
            schema
                .lists
                .get(&format!("axis_{}", axis))
                .cloned()
                .unwrap_or_default()
        };
        let mut bad = Vec::new();
        for c in &claims {
            for axis in [
                "evidence_kind",
                "independence",
                "universality",
                "data_relation",
                "physical_scope",
            ] {
                match c.axes.get(axis) {
                    None => bad.push(format!("{}: {} 欠落", c.id, axis)),
                    Some(v) if !vocab(axis).contains(v) => {
                        bad.push(format!("{}: {}='{}' 語彙外", c.id, axis, v))
                    }
                    _ => {}
                }
            }
            match c.axes.get("layer") {
                None => bad.push(format!("{}: layer 欠落", c.id)),
                Some(v) if !LAYERS.contains(&v.as_str()) => {
                    bad.push(format!("{}: layer='{}' 語彙外", c.id, v))
                }
                _ => {}
            }
        }
        bad.truncate(10);
        check(
            "[S3] claims.yml 全主張の 6 軸完備・語彙適合",
            bad.is_empty(),
            if bad.is_empty() {
                format!("{} 主張", claims.len())
            } else {
                format!("{:?}", bad)
            },
        );
    }

    let ax = |c: &Claim, k: &str| c.axes.get(k).cloned().unwrap_or_default();

    // ---- [S4] R1/R2 ----
    {
        let n_nat = claims
            .iter()
            .filter(|c| ax(c, "evidence_kind") == "natural_observation")
            .count();
        let n_ext = claims
            .iter()
            .filter(|c| ax(c, "evidence_kind") == "external_replication")
            .count();
        check(
            "[S4] R1/R2 自然観測・外部再現の証拠 0 件 (正直な残高の機械化)",
            n_nat == 0 && n_ext == 0,
            format!(
                "natural_observation {} / external_replication {}",
                n_nat, n_ext
            ),
        );
    }

    // ---- [S5] R3 ----
    {
        let bad: Vec<&str> = claims
            .iter()
            .filter(|c| ax(c, "independence") == "independent_author" && c.level != "C0")
            .map(|c| c.id.as_str())
            .collect();
        check(
            "[S5] R3 independent_author は C0 (外部の確立結果) のみ — 独立外部再現 0",
            bad.is_empty(),
            if bad.is_empty() {
                format!(
                    "C0 の {} 件のみ",
                    claims
                        .iter()
                        .filter(|c| ax(c, "independence") == "independent_author")
                        .count()
                )
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S6] R4 ----
    {
        let bad: Vec<&str> = claims
            .iter()
            .filter(|c| {
                ["core", "dynamics", "bridge"].contains(&ax(c, "layer").as_str())
                    && c.level != "C0"
                    && !["toy", "effective_model"].contains(&ax(c, "physical_scope").as_str())
            })
            .map(|c| c.id.as_str())
            .collect();
        check(
            "[S6] R4 core/dynamics/bridge 層 (C0 以外) は toy/effective_model のみ",
            bad.is_empty(),
            if bad.is_empty() {
                "toy mechanism → theory of nature の昇格なし".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S7] R5 ----
    {
        let bad: Vec<&str> = claims
            .iter()
            .filter(|c| ax(c, "data_relation") == "future_observation")
            .map(|c| c.id.as_str())
            .collect();
        check(
            "[S7] R5 future_observation を根拠にする主張は存在しない",
            bad.is_empty(),
            if bad.is_empty() {
                "未来のデータは主張を支えない".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S8] R6 ----
    {
        let bad: Vec<&str> = claims
            .iter()
            .filter(|c| {
                ax(c, "universality") == "continuum_universal"
                    && c.level != "C0"
                    && !["adapter", "instrument"].contains(&ax(c, "layer").as_str())
            })
            .map(|c| c.id.as_str())
            .collect();
        check(
            "[S8] R6 continuum_universal は adapter/instrument (または C0) の専有",
            bad.is_empty(),
            if bad.is_empty() {
                "連続普遍性の bridge/core への漏出なし".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S9] R7 ----
    {
        let bad: Vec<&str> = claims
            .iter()
            .filter(|c| {
                ax(c, "evidence_kind") == "internal_holdout"
                    && ax(c, "data_relation") != "preregistered_holdout"
            })
            .map(|c| c.id.as_str())
            .collect();
        check(
            "[S9] R7 internal_holdout ⇒ preregistered_holdout",
            bad.is_empty(),
            if bad.is_empty() {
                format!(
                    "{} 件の holdout 主張が整合",
                    claims
                        .iter()
                        .filter(|c| ax(c, "evidence_kind") == "internal_holdout")
                        .count()
                )
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S10] predictions.yml との照合 ----
    {
        let pred_text = rd("predictions.yml").expect("predictions.yml が読めない");
        let preds = parse_predictions(&pred_text);
        let mut bad = Vec::new();
        if preds.iter().any(|(id, _)| id == "PRED-019") {
            bad.push("PRED-019 が登録されている (解析的導出まで登録禁止 — PROMPT/10 §8)".into());
        }
        for (id, status) in &preds {
            match schema
                .predictions
                .iter()
                .find(|p| p.get("id").map(|s| s.as_str()) == Some(id.as_str()))
            {
                None => bad.push(format!("{} が schema 未分類", id)),
                Some(p) => {
                    let dr = p.get("data_relation").cloned().unwrap_or_default();
                    if dr == "future_observation" && status != "registered" {
                        bad.push(format!(
                            "{}: future_observation なのに status='{}' (自然の的中 0 の破れ)",
                            id, status
                        ));
                    }
                }
            }
        }
        for p in &schema.predictions {
            let id = p.get("id").cloned().unwrap_or_default();
            if !preds.iter().any(|(i, _)| *i == id) {
                bad.push(format!("schema の {} が predictions.yml に無い", id));
            }
        }
        let n_hit_nat = preds
            .iter()
            .filter(|(id, status)| {
                status == "scored-hit"
                    && schema.predictions.iter().any(|p| {
                        p.get("id").map(|s| s.as_str()) == Some(id.as_str())
                            && p.get("physical_scope").map(|s| s.as_str()) == Some("natural")
                            && p.get("data_relation").map(|s| s.as_str())
                                == Some("future_observation")
                    })
            })
            .count();
        check(
            "[S10] predictions 照合 — 自然の観測量の的中 (natural×future×hit) = 0・PRED-019 不在",
            bad.is_empty() && n_hit_nat == 0,
            if bad.is_empty() {
                format!(
                    "{} PRED 分類済み / 自然の的中 {} (公知測定値 holdout の hit は別枠)",
                    preds.len(),
                    n_hit_nat
                )
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S11] 文書アンカー ----
    {
        let mut bad = Vec::new();
        let anchors: [(&str, &[&str]); 6] = [
            (
                "docs/qrn-core-v1-spec.md",
                &[
                    "kinematics: defined",
                    "dynamics: model_family_only",
                    "geometry_bridge: conjectural",
                    "gravity_bridge: unsupported",
                    "empirical_prediction: none",
                ],
            ),
            (
                "docs/qrn-terminology.md",
                &[
                    "RegulatorSiteId -> RelationalNodeId",
                    "EvolutionParameter -> ProperTime",
                    "GaussianFermionState",
                ],
            ),
            ("docs/uft-v0.7.md", &["監査注記 (v27.1"]),
            ("docs/uft-v1.0.md", &["監査注記 (v27.1"]),
            ("docs/uft-v27.0.md", &["監査注記 (v27.1"]),
            ("README.md", &["自然の観測量の的中", "独立外部再現 0"]),
        ];
        for (path, needles) in anchors {
            match rd(path) {
                Ok(text) => {
                    for n in needles {
                        if !text.contains(n) {
                            bad.push(format!("{}: 「{}」が無い", path, n));
                        }
                    }
                }
                Err(_) => bad.push(format!("{} が読めない", path)),
            }
        }
        // claims の限定・降格アンカー
        if let Some(c) = claims.iter().find(|c| c.id == "QRN-C0-001") {
            if !c.claim.contains("無条件の『同値性』ではない") {
                bad.push("QRN-C0-001 に限定条項が無い".into());
            }
        } else {
            bad.push("QRN-C0-001 が無い".into());
        }
        if let Some(c) = claims.iter().find(|c| c.id == "QRN-GRAV-001") {
            if !c.status.contains("降格") {
                bad.push("QRN-GRAV-001 が降格されていない".into());
            }
        } else {
            bad.push("QRN-GRAV-001 が無い".into());
        }
        check(
            "[S11] 文書アンカー — 状態表示・監査注記 (v0.7/v1.0/v27.0)・C0-001 限定・GRAV-001 降格・README 残高",
            bad.is_empty(),
            if bad.is_empty() {
                "6 文書 + 2 主張".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S12] 禁止 impl From の不在 (sim/src 全走査) ----
    {
        const FORBIDDEN_IMPLS: [&str; 3] = [
            "impl From<RegulatorSiteId> for RelationalNodeId",
            "impl From<EvolutionParameter> for ProperTime",
            "impl From<ExternalMetricSource> for EmergentMetricCandidate",
        ];
        // 監査基盤自身はパターン定義を含むため対象外 (v61_ledger [6] の exempt と同型)
        const EXEMPT: [&str; 1] = ["v271_core_audit.rs"];
        let mut hits = Vec::new();
        let mut scan = |dir: &str| {
            if let Ok(rd) = fs::read_dir(format!("{}/{}", root, dir)) {
                for e in rd.filter_map(|e| e.ok()) {
                    let p = e.path();
                    let name = p
                        .file_name()
                        .map(|f| f.to_string_lossy().to_string())
                        .unwrap_or_default();
                    if EXEMPT.contains(&name.as_str()) {
                        continue;
                    }
                    if p.extension().map(|x| x == "rs").unwrap_or(false) {
                        if let Ok(text) = fs::read_to_string(&p) {
                            for f in FORBIDDEN_IMPLS {
                                if text.contains(f) {
                                    hits.push(format!("{}: {}", p.display(), f));
                                }
                            }
                        }
                    }
                }
            }
        };
        scan("sim/src");
        scan("sim/src/bin");
        check(
            "[S12] 禁止 impl From の不在 (regulator→ontology / 発展→固有時間 / 外部→創発計量)",
            hits.is_empty(),
            if hits.is_empty() {
                "sim/src 全走査 — 検出 0".into()
            } else {
                format!("{:?}", hits)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "層分類は完備・昇格禁止は保たれている — 未分類 0 件 (v27.2 への進行ゲート開通)"
        } else {
            "**意味論の破れ** — 分類・軸・文書アンカーを修正せよ (v27.2 に進まない)"
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
