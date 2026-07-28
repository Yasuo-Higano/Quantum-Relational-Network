//! v27.3 器械台帳 (instruments.yml) の機械検証 (PROMPT/10 §5)
//!
//! uft-v27.0.md §D.3 の器械群を QRN-Matter-on-Background Adapter v1 (family:
//! adapter) と QRN-Metrology Suite v1 (family: metrology) に分離・凍結した台帳を
//! 検査する:
//!   [I0] instruments.yml の書式解析 (器械 22 + 常設回帰 5)
//!   [I1] 必須 16 フィールドの完備と family 語彙 (adapter | metrology)
//!   [I2] concept ↔ core.schema.yml の層一致 (adapter → adapter / metrology → instrument)
//!   [I3] 較正記録・認証ファイルの実在 + certificate SHA-256 (先頭 16 桁) の一致
//!        — 器械の無断変更 (認証なしの頂点・規約の書き換え) を FAIL として検出
//!   [I4] allowed_claims の id が claims.yml に実在
//!   [I5] forbidden_interpretations が全器械で非空 (凍結解釈の明示)
//!   [I6] 常設回帰 (4 比 universality・和則・Ward 64/kernel・fork) の記録照合 —
//!        result が must_contain (凍結判定文) を含む
//!   [I7] §D.3 要件の被覆 — PROMPT/10 §5 指定の 9 器械が全て登録済み
//!
//! 本監査はスイートの常時実行層 (ALWAYS_RUN)。恒久解釈: 全器械の成果は
//! 「測定器が正しい」ことの証明であり QRN・創発重力の証拠ではない。

use std::collections::BTreeMap;
use std::fs;
use std::path::Path;
use uft_sim::sha256_hex;

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

/// instruments.yml — セクション (instruments / regressions) 付き制限 YAML
fn parse_registry(
    text: &str,
) -> Result<(Vec<BTreeMap<String, String>>, Vec<BTreeMap<String, String>>), String> {
    let mut instruments = Vec::new();
    let mut regressions = Vec::new();
    let mut sec = 0u8; // 0 = top, 1 = instruments, 2 = regressions
    for (lno, raw) in text.lines().enumerate() {
        let lno = lno + 1;
        let line = raw.trim_end();
        if line.trim_start().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        if !line.starts_with(' ') && !line.starts_with('-') {
            let (k, _) = split_kv(line).ok_or(format!("{}行目: トップレベル書式", lno))?;
            sec = match k.as_str() {
                "instruments" => 1,
                "regressions" => 2,
                _ => 0,
            };
        } else if let Some(rest) = line.strip_prefix("- ") {
            let (k, v) = split_kv(rest).ok_or(format!("{}行目: 項目書式", lno))?;
            let mut m = BTreeMap::new();
            m.insert(k, v);
            match sec {
                1 => instruments.push(m),
                2 => regressions.push(m),
                _ => return Err(format!("{}行目: セクション外の項目", lno)),
            }
        } else if let Some(rest) = line.strip_prefix("  ") {
            let (k, v) = split_kv(rest).ok_or(format!("{}行目: フィールド書式", lno))?;
            let cur = match sec {
                1 => instruments.last_mut(),
                2 => regressions.last_mut(),
                _ => None,
            }
            .ok_or(format!("{}行目: 項目外のフィールド", lno))?;
            cur.insert(k, v);
        } else {
            return Err(format!("{}行目: 解釈できない行", lno));
        }
    }
    Ok((instruments, regressions))
}

/// core.schema.yml の concepts から (name → layer) を拾う
fn schema_concept_layers(text: &str) -> BTreeMap<String, String> {
    let mut out = BTreeMap::new();
    let mut in_concepts = false;
    let mut cur_name = String::new();
    for raw in text.lines() {
        let line = raw.trim_end();
        if line.trim_start().starts_with('#') {
            continue;
        }
        if !line.starts_with(' ') && !line.starts_with('-') && line.contains(':') {
            in_concepts = line.starts_with("concepts:");
            continue;
        }
        if !in_concepts {
            continue;
        }
        if let Some(rest) = line.strip_prefix("- name:") {
            cur_name = rest.trim().to_string();
        } else if let Some(rest) = line.strip_prefix("  layer:") {
            if !cur_name.is_empty() {
                out.insert(cur_name.clone(), rest.trim().to_string());
            }
        }
    }
    out
}

fn main() {
    uft_sim::self_test();
    println!("=== v27.3 器械台帳の機械検証 — Adapter v1 / Metrology Suite v1 の凍結 (PROMPT/10 §5) ===\n");
    let root = if Path::new("instruments.yml").exists() {
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

    // ---- [I0] 解析 ----
    let text = rd("instruments.yml").expect("instruments.yml が読めない");
    let (instruments, regressions) = match parse_registry(&text) {
        Ok(x) => x,
        Err(e) => {
            println!("  [FAIL] [I0] instruments.yml: {}", e);
            std::process::exit(1);
        }
    };
    let n_adapter = instruments
        .iter()
        .filter(|i| i.get("family").map(|s| s.as_str()) == Some("adapter"))
        .count();
    let n_metrology = instruments
        .iter()
        .filter(|i| i.get("family").map(|s| s.as_str()) == Some("metrology"))
        .count();
    check(
        "[I0] instruments.yml の解析",
        !instruments.is_empty() && !regressions.is_empty(),
        format!(
            "器械 {} (adapter {} + metrology {}) / 常設回帰 {}",
            instruments.len(),
            n_adapter,
            n_metrology,
            regressions.len()
        ),
    );

    // ---- [I1] 必須フィールドと family ----
    {
        const REQUIRED: [&str; 16] = [
            "id",
            "family",
            "concept",
            "name",
            "input_type",
            "output_type",
            "normalization",
            "regulator",
            "continuum_contract",
            "calibration_source",
            "negative_controls",
            "known_failure_modes",
            "allowed_claims",
            "forbidden_interpretations",
            "certificate_code",
            "certificate_sha256_16",
        ];
        let mut bad = Vec::new();
        for inst in &instruments {
            let id = inst.get("id").cloned().unwrap_or_default();
            for f in REQUIRED {
                if inst.get(f).map(|v| v.is_empty()).unwrap_or(true) {
                    bad.push(format!("{}: {} 欠落", id, f));
                }
            }
            let fam = inst.get("family").cloned().unwrap_or_default();
            if fam != "adapter" && fam != "metrology" {
                bad.push(format!("{}: family '{}' 語彙外", id, fam));
            }
        }
        bad.truncate(8);
        check(
            "[I1] 必須 16 フィールドの完備・family 語彙",
            bad.is_empty(),
            if bad.is_empty() {
                format!("{} 器械 × 16 フィールド", instruments.len())
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [I2] concept ↔ schema 層一致 ----
    {
        let schema = rd("core.schema.yml").expect("core.schema.yml が読めない");
        let layers = schema_concept_layers(&schema);
        let mut bad = Vec::new();
        for inst in &instruments {
            let id = inst.get("id").cloned().unwrap_or_default();
            let concept = inst.get("concept").cloned().unwrap_or_default();
            let fam = inst.get("family").cloned().unwrap_or_default();
            match layers.get(&concept) {
                None => bad.push(format!("{}: concept {} が schema に無い", id, concept)),
                Some(layer) => {
                    let want = if fam == "adapter" {
                        "adapter"
                    } else {
                        "instrument"
                    };
                    if layer != want {
                        bad.push(format!(
                            "{}: {} は layer {} (family {} と不整合)",
                            id, concept, layer, fam
                        ));
                    }
                }
            }
        }
        check(
            "[I2] concept ↔ core.schema.yml の層一致 (adapter → adapter / metrology → instrument)",
            bad.is_empty(),
            if bad.is_empty() {
                "台帳間の分類が単一の層辞書に従う".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [I3] ファイル実在と certificate ハッシュ ----
    {
        let mut bad = Vec::new();
        for inst in &instruments {
            let id = inst.get("id").cloned().unwrap_or_default();
            for src in inst
                .get("calibration_source")
                .map(|s| {
                    s.split(',')
                        .map(|x| x.trim().to_string())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
            {
                if !Path::new(&format!("{}/{}", root, src)).exists() {
                    bad.push(format!("{}: {} が実在しない", id, src));
                }
            }
            let code = inst.get("certificate_code").cloned().unwrap_or_default();
            let want = inst
                .get("certificate_sha256_16")
                .cloned()
                .unwrap_or_default();
            match fs::read(format!("{}/{}", root, code)) {
                Err(_) => bad.push(format!("{}: 認証ファイル {} が読めない", id, code)),
                Ok(bytes) => {
                    let got = &sha256_hex(&bytes)[..16];
                    if got != want {
                        bad.push(format!(
                            "{}: {} の sha256 先頭 16 桁が不一致 (記録 {} / 現物 {}) — 器械の無断変更",
                            id, code, want, got
                        ));
                    }
                }
            }
        }
        check(
            "[I3] 較正記録・認証ファイルの実在 + certificate SHA-256 の一致",
            bad.is_empty(),
            if bad.is_empty() {
                format!("{} 器械の認証ハッシュが凍結どおり", instruments.len())
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [I4] allowed_claims の実在 ----
    {
        let claims = rd("claims.yml").expect("claims.yml が読めない");
        let mut bad = Vec::new();
        for inst in &instruments {
            let id = inst.get("id").cloned().unwrap_or_default();
            for cid in inst
                .get("allowed_claims")
                .map(|s| {
                    s.split(',')
                        .map(|x| x.trim().to_string())
                        .collect::<Vec<_>>()
                })
                .unwrap_or_default()
            {
                if !claims.contains(&format!("- id: {}", cid)) {
                    bad.push(format!("{}: {} が claims.yml に無い", id, cid));
                }
            }
        }
        check(
            "[I4] allowed_claims の id が claims.yml に実在",
            bad.is_empty(),
            if bad.is_empty() {
                "器械 → 主張の対応が台帳で閉じる".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [I5] 禁止解釈の明示 ----
    {
        let n = instruments
            .iter()
            .filter(|i| {
                !i.get("forbidden_interpretations")
                    .map(|v| v.is_empty())
                    .unwrap_or(true)
            })
            .count();
        check(
            "[I5] forbidden_interpretations が全器械で非空 (凍結解釈の明示)",
            n == instruments.len(),
            format!("{}/{} 器械", n, instruments.len()),
        );
    }

    // ---- [I6] 常設回帰の記録照合 ----
    {
        let mut bad = Vec::new();
        for reg in &regressions {
            let id = reg.get("id").cloned().unwrap_or_default();
            let result = reg.get("result").cloned().unwrap_or_default();
            let needle = reg.get("must_contain").cloned().unwrap_or_default();
            match rd(&result) {
                Err(_) => bad.push(format!("{}: {} が読めない", id, result)),
                Ok(text) => {
                    if !text.contains(&needle) {
                        bad.push(format!("{}: 凍結判定文が {} に無い", id, result));
                    }
                    if text.contains("[FAIL]") {
                        bad.push(format!("{}: {} に FAIL がある", id, result));
                    }
                }
            }
        }
        check(
            "[I6] 常設回帰 (4 比 universality・和則・Ward 64/kernel・fork) の記録照合",
            bad.is_empty(),
            if bad.is_empty() {
                format!("{} 回帰の凍結判定文と FAIL 不在を確認", regressions.len())
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [I7] §D.3 / PROMPT/10 §5 指定 9 器械の被覆 ----
    {
        const REQUIRED_CONCEPTS: [&str; 9] = [
            "BondACoupling",
            "MidpointModulation",
            "BelinfanteImprovement",
            "NullCombinationLadder",
            "ShellIntegration",
            "MatsubaraWard",
            "ContactTermCompletion",
            "DerivedAsymptoticExtrapolator",
            "StaggeredDiscretization",
        ];
        let mut missing = Vec::new();
        for c in REQUIRED_CONCEPTS {
            if !instruments
                .iter()
                .any(|i| i.get("concept").map(|s| s.as_str()) == Some(c))
            {
                missing.push(c);
            }
        }
        // Wilson は staggered との対 (対照離散化) として明示要求
        if !instruments
            .iter()
            .any(|i| i.get("concept").map(|s| s.as_str()) == Some("WilsonDiscretization"))
        {
            missing.push("WilsonDiscretization");
        }
        check(
            "[I7] PROMPT/10 §5 指定の器械 (BOND-A〜staggered/Wilson matching) の全登録",
            missing.is_empty(),
            if missing.is_empty() {
                "9 + Wilson 対照の 10 概念を被覆".into()
            } else {
                format!("未登録: {:?}", missing)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "Adapter v1 / Metrology Suite v1 は凍結された — 器械の変更は認証ハッシュの更新 (意識的な再認証) を要する"
        } else {
            "**台帳の破れ** — 器械・認証・回帰記録を修正せよ"
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
