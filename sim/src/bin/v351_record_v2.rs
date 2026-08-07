//! v35.1 実験前仕様修復の常設監査 (PROMPT/16 §3–4)
//!
//! 対象:
//!   (1) RECORD v2 三段 schema (preregistration / acquisition / disclosure) —
//!       v1 は不変のまま、**データ取得前に**別 hash で凍結 (取得後の schema 変更は
//!       「結果を見た後の変更」に見えるため禁止 — MANIFEST drift として検出)。
//!   (2) 相関粒度 3 型 (`sim/src/record_v2.rs`): OrderedShots / TimestampedBatches /
//!       AggregateCounts — iid 証明書は OrderedShots からのみ (禁止変換 30/31)。
//!   (3) D2-R outreach kit (`campaigns/d2r-v1-outreach/`) — d2r-v1 凍結物は不変・
//!       停止規則 failed_at_current_burden の事前登録・数は実記録のみ。
//!   (4) FollowUp 状態の主リポジトリ側整合 (replications.yml REP-001 —
//!       partially_replicated / shared human operator / external 0 維持)。
//!
//! 検査: [P0] モジュール自己検証 [P1] v1 凍結ピン [P2] v2 MANIFEST [P3] schema 構造
//! (additionalProperties:false 全 object・oneOf 3 粒度) [P4] fixture 正 4 負 3 +
//! 三段結束 (sha256 連鎖・予測/topology コミットメント照合・held-out 整合・± 対の
//! program hash 規則) [P5] 相関裁定の較正 (正例/Markov/drift/batch drift/過分散/
//! aggregate + level/power) [P6] outreach kit (凍結ピン・コードなし・停止規則)
//! [P7] replications.yml の内部整合。
//!
//! これは instrument maintenance であり物理的前進に数えない (PROMPT/16 §2)。

use std::fs;
use std::path::Path;
use uft_sim::record_v2::{
    assess, chi2_quantile, iid_certificate, record_v2_self_test, BatchRecord, CorrelationVerdict,
    RecordData, MIN_ORDERED_SHOTS,
};
use uft_sim::{self_test, sha256_hex, Rng};

// ---------------------------------------------------------------- 凍結ピン (v35.1 で確定)

/// v1 schema — 不変 (v34.6 凍結)
const PIN_RECORD_V1: &str = "780cd728f2ffed3fbc0aa48d6cee930f81231cf007a81c4558b0edc862991674";
/// d2r-v1 campaign MANIFEST — 不変 (protocol/schema/判定規則の凍結の外周ピン)
const PIN_D2R_MANIFEST: &str = "af967b3eb9a34511bd93785c05b4b6ca1899cd9459f0f9e2680219767e22de28";

// ---------------------------------------------------------------- 最小 JSON パーサ (v346 と同形)

#[derive(Debug, Clone)]
enum Json {
    Obj(Vec<(String, Json)>),
    Arr(Vec<Json>),
    Str(String),
    Num(f64),
    Bool(bool),
    Null,
}

impl Json {
    fn get(&self, k: &str) -> Option<&Json> {
        if let Json::Obj(fields) = self {
            fields.iter().find(|(key, _)| key == k).map(|(_, v)| v)
        } else {
            None
        }
    }
    fn as_str(&self) -> Option<&str> {
        if let Json::Str(s) = self {
            Some(s)
        } else {
            None
        }
    }
    fn keys(&self) -> Vec<&str> {
        if let Json::Obj(fields) = self {
            fields.iter().map(|(k, _)| k.as_str()).collect()
        } else {
            Vec::new()
        }
    }
}

struct Jp<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Jp<'a> {
    fn ws(&mut self) {
        while self.i < self.b.len() && (self.b[self.i] as char).is_whitespace() {
            self.i += 1;
        }
    }
    fn parse(&mut self) -> Result<Json, String> {
        self.ws();
        if self.i >= self.b.len() {
            return Err("途中終端".into());
        }
        match self.b[self.i] {
            b'{' => {
                self.i += 1;
                let mut fields = Vec::new();
                self.ws();
                if self.i < self.b.len() && self.b[self.i] == b'}' {
                    self.i += 1;
                    return Ok(Json::Obj(fields));
                }
                loop {
                    self.ws();
                    let key = match self.parse()? {
                        Json::Str(s) => s,
                        _ => return Err("オブジェクトのキーが文字列でない".into()),
                    };
                    self.ws();
                    if self.i >= self.b.len() || self.b[self.i] != b':' {
                        return Err("':' 欠落".into());
                    }
                    self.i += 1;
                    let v = self.parse()?;
                    fields.push((key, v));
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b'}') => {
                            self.i += 1;
                            return Ok(Json::Obj(fields));
                        }
                        _ => return Err("オブジェクトの区切り".into()),
                    }
                }
            }
            b'[' => {
                self.i += 1;
                let mut items = Vec::new();
                self.ws();
                if self.i < self.b.len() && self.b[self.i] == b']' {
                    self.i += 1;
                    return Ok(Json::Arr(items));
                }
                loop {
                    items.push(self.parse()?);
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => self.i += 1,
                        Some(b']') => {
                            self.i += 1;
                            return Ok(Json::Arr(items));
                        }
                        _ => return Err("配列の区切り".into()),
                    }
                }
            }
            b'"' => {
                self.i += 1;
                let mut s = String::new();
                while self.i < self.b.len() {
                    match self.b[self.i] {
                        b'"' => {
                            self.i += 1;
                            return Ok(Json::Str(s));
                        }
                        b'\\' => {
                            self.i += 1;
                            match self.b.get(self.i) {
                                Some(b'n') => s.push('\n'),
                                Some(b't') => s.push('\t'),
                                Some(&c) => s.push(c as char),
                                None => return Err("エスケープ途中終端".into()),
                            }
                            self.i += 1;
                        }
                        c => {
                            s.push(c as char);
                            self.i += 1;
                        }
                    }
                }
                Err("文字列途中終端".into())
            }
            _ => {
                let start = self.i;
                while self.i < self.b.len()
                    && !matches!(self.b[self.i], b',' | b'}' | b']')
                    && !(self.b[self.i] as char).is_whitespace()
                {
                    self.i += 1;
                }
                let tok = std::str::from_utf8(&self.b[start..self.i]).map_err(|e| e.to_string())?;
                match tok {
                    "true" => Ok(Json::Bool(true)),
                    "false" => Ok(Json::Bool(false)),
                    "null" => Ok(Json::Null),
                    t => t
                        .parse::<f64>()
                        .map(Json::Num)
                        .map_err(|_| format!("不明トークン '{}'", t)),
                }
            }
        }
    }
}

fn parse_json(text: &str) -> Result<Json, String> {
    let mut p = Jp {
        b: text.as_bytes(),
        i: 0,
    };
    let v = p.parse()?;
    p.ws();
    if p.i != p.b.len() {
        return Err("末尾に余分な内容".into());
    }
    Ok(v)
}

// ---------------------------------------------------------------- v2 validator (schema の器械化)

fn is_hex(s: &str, n: usize) -> bool {
    s.len() == n && s.chars().all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
}

/// 許可キー集合との厳密一致 (additionalProperties:false の器械化)
fn only_keys(j: &Json, allowed: &[&str], required: &[&str], ctx: &str) -> Result<(), String> {
    for k in j.keys() {
        if !allowed.contains(&k) {
            return Err(format!("{}: 追加プロパティ '{}' (additionalProperties:false)", ctx, k));
        }
    }
    for k in required {
        if j.get(k).is_none() {
            return Err(format!("{}: 必須 '{}' 欠落", ctx, k));
        }
    }
    Ok(())
}

fn str_arr(j: &Json) -> Option<Vec<String>> {
    if let Json::Arr(a) = j {
        a.iter()
            .map(|v| v.as_str().map(|s| s.to_string()))
            .collect::<Option<Vec<_>>>()
    } else {
        None
    }
}

/// A. preregistration-v2 の意味論
fn validate_prereg_v2(j: &Json) -> Result<(), String> {
    let req = [
        "schema_version",
        "observation_contract_hash",
        "analysis_commit",
        "calibration_channels",
        "held_out_channels",
        "prediction_commitment_sha256",
        "score_rule",
        "model_misspecification_gates",
        "topology_and_mapping_commitment",
        "committed_at",
    ];
    only_keys(j, &req, &req, "prereg")?;
    if j.get("schema_version").and_then(|v| v.as_str()) != Some("record-v2-preregistration") {
        return Err("schema_version 不一致".into());
    }
    for k in [
        "observation_contract_hash",
        "prediction_commitment_sha256",
        "topology_and_mapping_commitment",
    ] {
        let s = j.get(k).and_then(|v| v.as_str()).ok_or(format!("{} 型", k))?;
        if !is_hex(s, 64) {
            return Err(format!("{} が sha256 hex でない", k));
        }
    }
    let ac = j
        .get("analysis_commit")
        .and_then(|v| v.as_str())
        .ok_or("analysis_commit 型")?;
    if !is_hex(ac, 40) {
        return Err("analysis_commit が git sha (hex40) でない".into());
    }
    let cal = j
        .get("calibration_channels")
        .and_then(str_arr)
        .ok_or("calibration_channels 型")?;
    let held = j
        .get("held_out_channels")
        .and_then(str_arr)
        .ok_or("held_out_channels 型")?;
    if cal.is_empty() || held.is_empty() {
        return Err("channels が空".into());
    }
    if held.iter().any(|h| cal.contains(h)) {
        return Err("held_out と calibration が交差 (held-out ではない)".into());
    }
    let sr = j.get("score_rule").ok_or("score_rule 欠落")?;
    only_keys(sr, &["type", "alpha", "notes"], &["type", "alpha"], "score_rule")?;
    if sr.get("type").and_then(|v| v.as_str()) != Some("held_out_interval_hit") {
        return Err("score_rule.type が登録語彙にない".into());
    }
    match sr.get("alpha") {
        Some(Json::Num(a)) if *a > 0.0 && *a <= 0.1 => {}
        _ => return Err("score_rule.alpha ∉ (0, 0.1]".into()),
    }
    let g = j.get("model_misspecification_gates").ok_or("gates 欠落")?;
    only_keys(
        g,
        &[
            "drift_gate_alpha",
            "correlation_gate_alpha",
            "signed_linearity_epsilon_pair",
            "negative_control_channels",
        ],
        &[
            "drift_gate_alpha",
            "correlation_gate_alpha",
            "signed_linearity_epsilon_pair",
        ],
        "gates",
    )?;
    match g.get("signed_linearity_epsilon_pair") {
        Some(Json::Arr(a)) if a.len() == 2 => {
            for v in a {
                match v {
                    Json::Num(x) if *x > 0.0 => {}
                    _ => return Err("epsilon_pair に非正値".into()),
                }
            }
        }
        _ => return Err("signed_linearity_epsilon_pair は長さ 2".into()),
    }
    Ok(())
}

/// data oneOf — 3 粒度のどれか厳密に 1 つに適合
fn validate_data_oneof(d: &Json) -> Result<&'static str, String> {
    let mut matched = Vec::new();
    // ordered_shots
    let m1 = (|| -> Result<(), String> {
        only_keys(
            d,
            &["kind", "acquired_at", "shots"],
            &["kind", "acquired_at", "shots"],
            "ordered_shots",
        )?;
        if d.get("kind").and_then(|v| v.as_str()) != Some("ordered_shots") {
            return Err("kind".into());
        }
        match d.get("shots") {
            Some(Json::Arr(a)) if !a.is_empty() => {
                for s in a {
                    match s {
                        Json::Num(x) if *x == 0.0 || *x == 1.0 => {}
                        _ => return Err("shots に {0,1} 以外".into()),
                    }
                }
                Ok(())
            }
            _ => Err("shots 欠落/空".into()),
        }
    })();
    if m1.is_ok() {
        matched.push("ordered_shots");
    }
    // timestamped_batches
    let m2 = (|| -> Result<(), String> {
        only_keys(d, &["kind", "batches"], &["kind", "batches"], "batches")?;
        if d.get("kind").and_then(|v| v.as_str()) != Some("timestamped_batches") {
            return Err("kind".into());
        }
        match d.get("batches") {
            Some(Json::Arr(a)) if a.len() >= 2 => {
                for b in a {
                    only_keys(
                        b,
                        &["batch_id", "started_at", "ended_at", "n_shots", "n_ones"],
                        &["batch_id", "started_at", "ended_at", "n_shots", "n_ones"],
                        "batch",
                    )?;
                    let ns = match b.get("n_shots") {
                        Some(Json::Num(x)) if *x >= 1.0 => *x,
                        _ => return Err("n_shots".into()),
                    };
                    match b.get("n_ones") {
                        Some(Json::Num(x)) if *x >= 0.0 && *x <= ns => {}
                        _ => return Err("n_ones".into()),
                    }
                }
                Ok(())
            }
            _ => Err("batches < 2".into()),
        }
    })();
    if m2.is_ok() {
        matched.push("timestamped_batches");
    }
    // aggregate_counts
    let m3 = (|| -> Result<(), String> {
        only_keys(
            d,
            &["kind", "acquired_at", "n_shots", "n_ones"],
            &["kind", "acquired_at", "n_shots", "n_ones"],
            "aggregate",
        )?;
        if d.get("kind").and_then(|v| v.as_str()) != Some("aggregate_counts") {
            return Err("kind".into());
        }
        Ok(())
    })();
    if m3.is_ok() {
        matched.push("aggregate_counts");
    }
    match matched.len() {
        1 => Ok(matched[0]),
        0 => Err("data がどの粒度にも適合しない (oneOf 失敗)".into()),
        _ => Err("data が複数粒度に適合 (oneOf 失敗)".into()),
    }
}

/// B. acquisition-v2 の意味論 (± probe 対の program hash 規則を含む)
fn validate_acquisition_v2(j: &Json) -> Result<(), String> {
    let req = [
        "schema_version",
        "provenance",
        "evolution_provenance",
        "preregistration_sha256",
        "device",
        "operator",
        "acquisition_randomization_seed",
        "channels",
        "results_including_failures_public",
    ];
    only_keys(j, &req, &req, "acquisition")?;
    if j.get("schema_version").and_then(|v| v.as_str()) != Some("record-v2-acquisition") {
        return Err("schema_version 不一致".into());
    }
    if j.get("provenance").and_then(|v| v.as_str()) != Some("recorded_experimental") {
        return Err("provenance が recorded_experimental でない (lane 分離)".into());
    }
    let ep = j
        .get("evolution_provenance")
        .and_then(|v| v.as_str())
        .ok_or("evolution_provenance 欠落")?;
    if !["native_analog_evolution", "compiled_digital_channel_family", "emulator"].contains(&ep) {
        return Err(format!("evolution_provenance '{}' は語彙にない", ep));
    }
    let ps = j
        .get("preregistration_sha256")
        .and_then(|v| v.as_str())
        .ok_or("preregistration_sha256 欠落")?;
    if !is_hex(ps, 64) {
        return Err("preregistration_sha256 が hex64 でない".into());
    }
    let dev = j.get("device").ok_or("device 欠落")?;
    only_keys(
        dev,
        &["vendor", "device_id", "calibration_snapshot_hash"],
        &["vendor", "device_id"],
        "device",
    )?;
    let op = j.get("operator").ok_or("operator 欠落")?;
    only_keys(
        op,
        &["name", "organization", "organizationally_external"],
        &["name", "organization", "organizationally_external"],
        "operator",
    )?;
    match op.get("organizationally_external") {
        Some(Json::Bool(_)) => {}
        _ => return Err("organizationally_external が boolean でない".into()),
    }
    match j.get("results_including_failures_public") {
        Some(Json::Bool(true)) => {}
        _ => return Err("results_including_failures_public = true が必要".into()),
    }
    let chs = match j.get("channels") {
        Some(Json::Arr(a)) if !a.is_empty() => a,
        _ => return Err("channels が空/欠落".into()),
    };
    let ch_req = [
        "channel_id",
        "source_id",
        "target_id",
        "probe_sign",
        "epsilon",
        "evolution_time",
        "evolution_family_id",
        "job_id",
        "program_hash",
        "compiler_or_pulse_version",
        "calibration_snapshot_hash",
        "used_in_calibration",
        "raw_artifact_sha256",
        "data",
    ];
    // (source, target, ε, t, family) → (sign, program_hash) — ± 対の同一 evolution program 検査
    let mut pairs: Vec<((String, String, String, String, String), (String, String))> = Vec::new();
    for ch in chs {
        only_keys(ch, &ch_req, &ch_req, "channel")?;
        let sign = ch
            .get("probe_sign")
            .and_then(|v| v.as_str())
            .ok_or("probe_sign 欠落")?;
        if !["plus", "minus"].contains(&sign) {
            return Err(format!("probe_sign '{}' は語彙にない", sign));
        }
        match ch.get("epsilon") {
            Some(Json::Num(e)) if *e > 0.0 => {}
            _ => return Err("epsilon ≤ 0".into()),
        }
        match ch.get("evolution_time") {
            Some(Json::Num(t)) if *t >= 0.0 => {}
            _ => return Err("evolution_time < 0".into()),
        }
        for k in ["program_hash", "calibration_snapshot_hash", "raw_artifact_sha256"] {
            let s = ch.get(k).and_then(|v| v.as_str()).ok_or(format!("{} 欠落", k))?;
            if !is_hex(s, 64) {
                return Err(format!("{} が hex64 でない", k));
            }
        }
        match ch.get("used_in_calibration") {
            Some(Json::Bool(_)) => {}
            _ => return Err("used_in_calibration (必須 boolean) 欠落".into()),
        }
        validate_data_oneof(ch.get("data").ok_or("data 欠落")?)?;
        let key = (
            ch.get("source_id").and_then(|v| v.as_str()).unwrap_or("").to_string(),
            ch.get("target_id").and_then(|v| v.as_str()).unwrap_or("").to_string(),
            format!("{:?}", ch.get("epsilon").map(|v| if let Json::Num(x) = v { *x } else { 0.0 })),
            format!("{:?}", ch.get("evolution_time").map(|v| if let Json::Num(x) = v { *x } else { 0.0 })),
            ch.get("evolution_family_id").and_then(|v| v.as_str()).unwrap_or("").to_string(),
        );
        let ph = ch.get("program_hash").and_then(|v| v.as_str()).unwrap_or("").to_string();
        pairs.push((key, (sign.to_string(), ph)));
    }
    // ± 対: 同一 (source,target,ε,t,family) で符号が異なる 2 チャネルは同一 evolution program
    for i in 0..pairs.len() {
        for k in (i + 1)..pairs.len() {
            if pairs[i].0 == pairs[k].0
                && pairs[i].1 .0 != pairs[k].1 .0
                && pairs[i].1 .1 != pairs[k].1 .1
            {
                return Err("± probe 対で program_hash 不一致 (evolution program が同一でない)".into());
            }
        }
    }
    Ok(())
}

/// C. disclosure-v2 の意味論
fn validate_disclosure_v2(j: &Json) -> Result<(), String> {
    let req = [
        "schema_version",
        "preregistration_sha256",
        "acquisition_sha256",
        "disclosed_prediction",
        "disclosed_topology",
        "disclosed_label_mapping",
        "commitment_verification",
        "score",
        "final_verdict",
        "failures_public",
    ];
    only_keys(j, &req, &req, "disclosure")?;
    if j.get("schema_version").and_then(|v| v.as_str()) != Some("record-v2-disclosure") {
        return Err("schema_version 不一致".into());
    }
    for k in ["preregistration_sha256", "acquisition_sha256"] {
        let s = j.get(k).and_then(|v| v.as_str()).ok_or(format!("{} 欠落", k))?;
        if !is_hex(s, 64) {
            return Err(format!("{} が hex64 でない", k));
        }
    }
    let cv = j.get("commitment_verification").ok_or("commitment_verification 欠落")?;
    only_keys(
        cv,
        &["prediction_sha256_matches", "topology_mapping_sha256_matches"],
        &["prediction_sha256_matches", "topology_mapping_sha256_matches"],
        "commitment_verification",
    )?;
    let sc = j.get("score").ok_or("score 欠落")?;
    only_keys(
        sc,
        &["held_out_total", "hits", "misses", "straddled", "notes"],
        &["held_out_total", "hits", "misses", "straddled"],
        "score",
    )?;
    let num = |k: &str| -> Result<f64, String> {
        match sc.get(k) {
            Some(Json::Num(x)) if *x >= 0.0 && x.fract() == 0.0 => Ok(*x),
            _ => Err(format!("score.{} が非負整数でない", k)),
        }
    };
    let (tot, h, m, s) = (num("held_out_total")?, num("hits")?, num("misses")?, num("straddled")?);
    if h + m + s != tot {
        return Err("score の会計不一致 (hits+misses+straddled ≠ total)".into());
    }
    let fv = j.get("final_verdict").and_then(|v| v.as_str()).ok_or("final_verdict 欠落")?;
    if !["hit", "miss", "straddled", "out_of_domain", "aborted"].contains(&fv) {
        return Err(format!("final_verdict '{}' は語彙にない", fv));
    }
    match j.get("failures_public") {
        Some(Json::Bool(true)) => {}
        _ => return Err("failures_public = true が必要".into()),
    }
    Ok(())
}

// ---------------------------------------------------------------- main

fn main() {
    self_test();
    record_v2_self_test().expect("record_v2 self test");
    println!("=== v35.1 実験前仕様修復 — RECORD v2・相関粒度型・outreach kit (PROMPT/16 §3–4) ===\n");
    let root = if Path::new("reproducer/real_data/RECORDED_LANE.md").exists() {
        "."
    } else {
        ".."
    };
    let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p)).unwrap_or_default();
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---------------- [P0] モジュール自己検証 ----------------
    {
        let q = chi2_quantile(11.0, 0.995);
        check(
            "[P0] record_v2 自己検証 (χ² 分位・粒度到達集合) + lib self_test",
            (q - 26.7568).abs() < 1e-2,
            format!("χ²₁₁(0.995) = {:.4} (期待 26.7568)", q),
        );
    }

    // ---------------- [P1] v1 凍結ピン ----------------
    {
        let v1 = rd("reproducer/real_data/RECORD.schema.json");
        let h = sha256_hex(v1.as_bytes());
        let lane = rd("reproducer/real_data/RECORDED_LANE.md");
        let anchors_ok = ["RECORD v2 — 三段 schema", "禁止変換 30", "禁止変換 31", "CorrelationUnassessed"]
            .iter()
            .all(|a| lane.contains(a));
        check(
            "[P1] v1 凍結: RECORD.schema.json 不変 (v34.6 ピン) + lane 文書に v2 節",
            h == PIN_RECORD_V1 && anchors_ok,
            format!("sha256 = {}…", &h[..16]),
        );
    }

    // ---------------- [P2] v2 MANIFEST (凍結の drift 検出) ----------------
    {
        let man = rd("reproducer/real_data/RECORD_V2/MANIFEST.sha256");
        let mut bad = Vec::new();
        let mut n = 0usize;
        for line in man.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut it = line.split_whitespace();
            let (Some(want), Some(fname)) = (it.next(), it.next()) else {
                bad.push(format!("MANIFEST 行不明: {}", line));
                continue;
            };
            n += 1;
            let t = rd(&format!("reproducer/real_data/RECORD_V2/{}", fname));
            if t.is_empty() {
                bad.push(format!("{} が無い", fname));
            } else if sha256_hex(t.as_bytes()) != want {
                bad.push(format!("{} の sha256 が MANIFEST と不一致 (取得前凍結の drift)", fname));
            }
        }
        if n != 3 {
            bad.push(format!("schema は 3 段 (実測 {})", n));
        }
        check(
            "[P2] v2 MANIFEST: 三段 schema の凍結 hash 一致 (取得後の変更は drift として検出)",
            bad.is_empty(),
            if bad.is_empty() {
                "prereg/acquisition/disclosure 3 本認証".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [P3] schema 構造の自己検査 ----------------
    {
        let mut bad = Vec::new();
        // 全 object ノードに additionalProperties:false (テキスト検査 + parse 検査)
        for f in [
            "preregistration-v2.schema.json",
            "acquisition-v2.schema.json",
            "disclosure-v2.schema.json",
        ] {
            let t = rd(&format!("reproducer/real_data/RECORD_V2/{}", f));
            match parse_json(&t) {
                Err(e) => bad.push(format!("{}: parse 失敗 {}", f, e)),
                Ok(_) => {
                    let n_obj = t.matches("\"type\": \"object\"").count();
                    let n_ap = t.matches("\"additionalProperties\": false").count();
                    if n_ap < n_obj {
                        bad.push(format!(
                            "{}: object {} 個に additionalProperties:false が {} 個",
                            f, n_obj, n_ap
                        ));
                    }
                }
            }
        }
        // acquisition の oneOf 3 粒度
        let acq = rd("reproducer/real_data/RECORD_V2/acquisition-v2.schema.json");
        for k in ["ordered_shots", "timestamped_batches", "aggregate_counts", "oneOf"] {
            if !acq.contains(k) {
                bad.push(format!("acquisition schema に {} が無い", k));
            }
        }
        // 必須フィールドの表明 (PROMPT/16 §3.1 の指定)
        for k in [
            "used_in_calibration",
            "probe_sign",
            "epsilon",
            "evolution_time",
            "evolution_family_id",
            "program_hash",
            "acquisition_randomization_seed",
            "raw_artifact_sha256",
            "evolution_provenance",
        ] {
            if !acq.contains(k) {
                bad.push(format!("acquisition schema に {} が無い", k));
            }
        }
        let pre = rd("reproducer/real_data/RECORD_V2/preregistration-v2.schema.json");
        for k in [
            "held_out_channels",
            "prediction_commitment_sha256",
            "score_rule",
            "model_misspecification_gates",
            "topology_and_mapping_commitment",
        ] {
            if !pre.contains(k) {
                bad.push(format!("prereg schema に {} が無い", k));
            }
        }
        check(
            "[P3] schema 構造: 全 object に additionalProperties:false・oneOf 3 粒度・必須フィールド",
            bad.is_empty(),
            if bad.is_empty() {
                "3 schema とも構造規約を満たす".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [P4] fixture 較正 (正 4 負 3) + 三段結束 ----------------
    {
        let fx = |name: &str| rd(&format!("reproducer/real_data/RECORD_V2/FIXTURES_V2/{}", name));
        let mut bad = Vec::new();
        // 正例
        let prereg_t = fx("prereg_valid.json");
        let prereg = parse_json(&prereg_t).expect("prereg parse");
        if let Err(e) = validate_prereg_v2(&prereg) {
            bad.push(format!("prereg_valid が不適合: {}", e));
        }
        let acq_ord_t = fx("acquisition_ordered_valid.json");
        let acq_ord = parse_json(&acq_ord_t).expect("acq parse");
        if let Err(e) = validate_acquisition_v2(&acq_ord) {
            bad.push(format!("acquisition_ordered_valid が不適合: {}", e));
        }
        let acq_bat_t = fx("acquisition_batches_valid.json");
        if let Err(e) = validate_acquisition_v2(&parse_json(&acq_bat_t).expect("parse")) {
            bad.push(format!("acquisition_batches_valid が不適合: {}", e));
        }
        let disc_t = fx("disclosure_valid.json");
        let disc = parse_json(&disc_t).expect("disc parse");
        if let Err(e) = validate_disclosure_v2(&disc) {
            bad.push(format!("disclosure_valid が不適合: {}", e));
        }
        // 負例
        for (name, why) in [
            ("acquisition_invalid_extra_field.json", "追加プロパティ"),
            ("acquisition_invalid_shots_and_counts.json", "oneOf"),
        ] {
            match parse_json(&fx(name)) {
                Ok(j) => {
                    if validate_acquisition_v2(&j).is_ok() {
                        bad.push(format!("{} が拒否されない ({})", name, why));
                    }
                }
                Err(e) => bad.push(format!("{} parse 失敗 {}", name, e)),
            }
        }
        match parse_json(&fx("prereg_invalid_no_heldout.json")) {
            Ok(j) => {
                if validate_prereg_v2(&j).is_ok() {
                    bad.push("prereg_invalid_no_heldout が拒否されない".into());
                }
            }
            Err(e) => bad.push(format!("prereg_invalid parse 失敗 {}", e)),
        }
        // 三段結束: sha256 連鎖
        let h_prereg = sha256_hex(prereg_t.as_bytes());
        for (t, name) in [(&acq_ord_t, "acq_ordered"), (&acq_bat_t, "acq_batches")] {
            let j = parse_json(t).unwrap();
            if j.get("preregistration_sha256").and_then(|v| v.as_str()) != Some(h_prereg.as_str()) {
                bad.push(format!("{} の preregistration_sha256 が prereg 実 hash と不一致", name));
            }
        }
        let h_acq = sha256_hex(acq_ord_t.as_bytes());
        if disc.get("acquisition_sha256").and_then(|v| v.as_str()) != Some(h_acq.as_str()) {
            bad.push("disclosure の acquisition_sha256 が実 hash と不一致".into());
        }
        // コミットメント照合 (予測 / topology+mapping)
        let pred = disc.get("disclosed_prediction").and_then(|v| v.as_str()).unwrap_or("");
        let want_pred = prereg
            .get("prediction_commitment_sha256")
            .and_then(|v| v.as_str())
            .unwrap_or("");
        if sha256_hex(pred.as_bytes()) != want_pred {
            bad.push("sha256(disclosed_prediction) ≠ prediction_commitment".into());
        }
        let topo = disc.get("disclosed_topology").and_then(|v| v.as_str()).unwrap_or("");
        let map = disc.get("disclosed_label_mapping").and_then(|v| v.as_str()).unwrap_or("");
        let combined = format!("{}\n{}", topo, map);
        let want_topo = prereg
            .get("topology_and_mapping_commitment")
            .and_then(|v| v.as_str())
            .unwrap_or("");
        if sha256_hex(combined.as_bytes()) != want_topo {
            bad.push("sha256(topology+mapping) ≠ commitment".into());
        }
        // held-out 整合: held_out_channels のチャネルは used_in_calibration = false
        let held = prereg.get("held_out_channels").and_then(str_arr).unwrap_or_default();
        if let Some(Json::Arr(chs)) = acq_ord.get("channels") {
            for ch in chs {
                let cid = ch.get("channel_id").and_then(|v| v.as_str()).unwrap_or("");
                let uic = matches!(ch.get("used_in_calibration"), Some(Json::Bool(true)));
                if held.contains(&cid.to_string()) && uic {
                    bad.push(format!("held-out {} が used_in_calibration = true", cid));
                }
            }
        }
        // ± 対規則の較正 (inline — 同一 (source,target,ε,t,family) で符号違い・hash 違い → 拒否)
        let mk = |ph1: &str, ph2: &str| -> String {
            format!(
                r#"{{"schema_version": "record-v2-acquisition", "provenance": "recorded_experimental",
"evolution_provenance": "native_analog_evolution",
"preregistration_sha256": "{h}", "device": {{"vendor": "v", "device_id": "d"}},
"operator": {{"name": "n", "organization": "o", "organizationally_external": false}},
"acquisition_randomization_seed": "s",
"channels": [
 {{"channel_id": "c+", "source_id": "n1", "target_id": "n2", "probe_sign": "plus",
   "epsilon": 0.1, "evolution_time": 1.0, "evolution_family_id": "F", "job_id": "j1",
   "program_hash": "{p1}", "compiler_or_pulse_version": "v",
   "calibration_snapshot_hash": "{h}", "used_in_calibration": true,
   "raw_artifact_sha256": "{h}",
   "data": {{"kind": "aggregate_counts", "acquired_at": "t", "n_shots": 10, "n_ones": 5}}}},
 {{"channel_id": "c-", "source_id": "n1", "target_id": "n2", "probe_sign": "minus",
   "epsilon": 0.1, "evolution_time": 1.0, "evolution_family_id": "F", "job_id": "j2",
   "program_hash": "{p2}", "compiler_or_pulse_version": "v",
   "calibration_snapshot_hash": "{h}", "used_in_calibration": true,
   "raw_artifact_sha256": "{h}",
   "data": {{"kind": "aggregate_counts", "acquired_at": "t", "n_shots": 10, "n_ones": 5}}}}
], "results_including_failures_public": true}}"#,
                h = h_prereg,
                p1 = ph1,
                p2 = ph2
            )
        };
        let same = mk(&h_prereg, &h_prereg);
        let diff = mk(&h_prereg, &sha256_hex(b"other-program"));
        if validate_acquisition_v2(&parse_json(&same).unwrap()).is_err() {
            bad.push("± 対 (同一 program_hash) が拒否された".into());
        }
        if validate_acquisition_v2(&parse_json(&diff).unwrap()).is_ok() {
            bad.push("± 対 (program_hash 不一致) が拒否されない".into());
        }
        check(
            "[P4] fixture 較正 正 4/負 3 + 三段 sha256 結束 + コミットメント照合 + held-out 整合 + ± 対規則",
            bad.is_empty(),
            if bad.is_empty() {
                format!("prereg {}… → acq {}… → disclosure (連鎖照合 OK)", &h_prereg[..8], &h_acq[..8])
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [P5] 相関粒度の裁定較正 ----------------
    {
        let alpha = 0.01;
        let mut bad = Vec::new();
        // (a) fixture の ordered (0110 反復, n=40) → IidConsistent
        let fx_t = rd("reproducer/real_data/RECORD_V2/FIXTURES_V2/acquisition_ordered_valid.json");
        let fx_j = parse_json(&fx_t).expect("parse");
        let mut fixture_shots: Vec<u8> = Vec::new();
        if let Some(Json::Arr(chs)) = fx_j.get("channels") {
            if let Some(Json::Arr(a)) = chs[0].get("data").and_then(|d| d.get("shots")) {
                fixture_shots = a
                    .iter()
                    .map(|v| if let Json::Num(x) = v { *x as u8 } else { 0 })
                    .collect();
            }
        }
        let v = assess(&RecordData::OrderedShots(fixture_shots.clone()), alpha);
        if iid_certificate(&v).is_none() {
            bad.push(format!("fixture ordered が {} (期待 iid_consistent)", v.as_str()));
        }
        // (b) 持続的 Markov 鎖 (同一周辺分布 0.5, stay=0.9, n=400) → serial correlation
        let mut rng = Rng::new(35101);
        let mut s = 0u8;
        let markov: Vec<u8> = (0..400)
            .map(|_| {
                if rng.f64() > 0.9 {
                    s ^= 1;
                }
                s
            })
            .collect();
        let v = assess(&RecordData::OrderedShots(markov), alpha);
        if !matches!(v, CorrelationVerdict::SerialCorrelationDetected { .. }) {
            bad.push(format!("Markov が {} (期待 serial_correlation)", v.as_str()));
        }
        // (c) drift (前半 p=0.15 / 後半 p=0.85, n=400) → drift
        let mut rng = Rng::new(35102);
        let drift: Vec<u8> = (0..400)
            .map(|i| {
                let p = if i < 200 { 0.15 } else { 0.85 };
                (rng.f64() < p) as u8
            })
            .collect();
        let v = assess(&RecordData::OrderedShots(drift), alpha);
        if !matches!(v, CorrelationVerdict::DriftDetected { .. }) {
            bad.push(format!("drift 列が {} (期待 drift)", v.as_str()));
        }
        // (d) 短列は資格を出さない
        let short: Vec<u8> = (0..20).map(|i| (i % 2) as u8).collect();
        let v = assess(&RecordData::OrderedShots(short), alpha);
        if !matches!(v, CorrelationVerdict::CorrelationUnresolved) {
            bad.push(format!("短列 (n=20 < {}) が {} (期待 unresolved)", MIN_ORDERED_SHOTS, v.as_str()));
        }
        // (e) 良性バッチ → CorrelationUnresolved (iid 証明書は出ない — 禁止変換 31)
        let benign = RecordData::TimestampedBatches(vec![
            BatchRecord { n_shots: 50, n_ones: 12 },
            BatchRecord { n_shots: 50, n_ones: 13 },
            BatchRecord { n_shots: 50, n_ones: 11 },
            BatchRecord { n_shots: 50, n_ones: 14 },
        ]);
        let v = assess(&benign, alpha);
        if !matches!(v, CorrelationVerdict::CorrelationUnresolved) || iid_certificate(&v).is_some() {
            bad.push(format!("良性バッチが {} (期待 unresolved・証明書なし)", v.as_str()));
        }
        // (f) バッチ間 drift (5/50 vs 45/50) → batch_drift
        let bd = RecordData::TimestampedBatches(vec![
            BatchRecord { n_shots: 50, n_ones: 5 },
            BatchRecord { n_shots: 50, n_ones: 45 },
        ]);
        let v = assess(&bd, alpha);
        if !matches!(v, CorrelationVerdict::BatchDriftDetected { .. }) {
            bad.push(format!("batch drift が {} (期待 batch_drift)", v.as_str()));
        }
        // (g) 過分散 (m=12, n_b=100, k = 38/62 交互 — 対では disjoint しない) → overdispersion
        let od = RecordData::TimestampedBatches(
            (0..12)
                .map(|i| BatchRecord {
                    n_shots: 100,
                    n_ones: if i % 2 == 0 { 38 } else { 62 },
                })
                .collect(),
        );
        let v = assess(&od, alpha);
        if !matches!(v, CorrelationVerdict::OverdispersionDetected { .. }) {
            bad.push(format!("過分散バッチが {} (期待 overdispersion)", v.as_str()));
        }
        // (h) AggregateCounts は常に Unassessed (禁止変換 30 — ショット数によらず)
        let agg = RecordData::AggregateCounts {
            n_shots: 1_000_000,
            n_ones: 500_000,
        };
        let v = assess(&agg, alpha);
        if !matches!(v, CorrelationVerdict::CorrelationUnassessed) || iid_certificate(&v).is_some() {
            bad.push(format!("aggregate が {} (期待 unassessed)", v.as_str()));
        }
        check(
            "[P5a] 裁定較正: iid/Markov/drift/短列/良性バッチ/batch drift/過分散/aggregate の 8 セル",
            bad.is_empty(),
            if bad.is_empty() {
                "全セル期待裁定 (IidConsistent は OrderedShots のみ)".into()
            } else {
                format!("{:?}", bad)
            },
        );

        // (i) level: iid (p=0.3, n=200) R=2000 — 偽検出率 ≤ 1.5α
        let mut rng = Rng::new(35103);
        let mut false_det = 0usize;
        let reps = 2000;
        for _ in 0..reps {
            let shots: Vec<u8> = (0..200).map(|_| (rng.f64() < 0.3) as u8).collect();
            match assess(&RecordData::OrderedShots(shots), alpha) {
                CorrelationVerdict::IidConsistent(_) => {}
                _ => false_det += 1,
            }
        }
        let rate = false_det as f64 / reps as f64;
        check(
            "[P5b] ordered level 較正: iid 偽検出率 ≤ 1.5α (CP ゲートの保守性)",
            rate <= 1.5 * alpha,
            format!("偽検出 {}/{} = {:.4} (α = {})", false_det, reps, rate, alpha),
        );
        // (j) power: Markov (stay 0.9, n=800) R=200 — 検出率 ≥ 0.95。
        // 検出 = 非 IidConsistent (持続鎖は遷移数ゲートだけでなく split-half でも
        // 正しく棄却され得る — どちらのモードでも iid 証明書は発行されない)。
        let mut rng = Rng::new(35104);
        let (mut det_serial, mut det_drift, mut miss) = (0usize, 0usize, 0usize);
        for _ in 0..200 {
            let mut s = (rng.f64() < 0.5) as u8;
            let shots: Vec<u8> = (0..800)
                .map(|_| {
                    if rng.f64() > 0.9 {
                        s ^= 1;
                    }
                    s
                })
                .collect();
            match assess(&RecordData::OrderedShots(shots), alpha) {
                CorrelationVerdict::SerialCorrelationDetected { .. } => det_serial += 1,
                CorrelationVerdict::DriftDetected { .. } => det_drift += 1,
                CorrelationVerdict::IidConsistent(_) => miss += 1,
                _ => miss += 1,
            }
        }
        check(
            "[P5c] ordered power 較正: Markov (stay 0.9, n=800) 非 IidConsistent 率 ≥ 0.95",
            det_serial + det_drift >= 190,
            format!(
                "検出 {}/200 (遷移数 {} + split-half {}), 誤資格 {}",
                det_serial + det_drift,
                det_serial,
                det_drift,
                miss
            ),
        );
        // (k) batches level: iid (m=8, n_b=100, p=0.3) R=2000 — 偽検出率 ≤ 2α (χ² 近似の登録 slack)
        let mut rng = Rng::new(35105);
        let mut false_det = 0usize;
        for _ in 0..2000 {
            let batches: Vec<BatchRecord> = (0..8)
                .map(|_| {
                    let k = (0..100).filter(|_| rng.f64() < 0.3).count();
                    BatchRecord {
                        n_shots: 100,
                        n_ones: k,
                    }
                })
                .collect();
            match assess(&RecordData::TimestampedBatches(batches), alpha) {
                CorrelationVerdict::CorrelationUnresolved => {}
                _ => false_det += 1,
            }
        }
        let rate = false_det as f64 / 2000.0;
        check(
            "[P5d] batches level 較正: iid 偽検出率 ≤ 2α (χ²_{m−1} は近似と登録)",
            rate <= 2.0 * alpha,
            format!("偽検出 {}/2000 = {:.4}", false_det, rate),
        );
        // (l) batches power: 交互 p = 0.38/0.62 (m=12, n_b=100) R=200 — 検出率 ≥ 0.9
        let mut rng = Rng::new(35106);
        let mut det = 0usize;
        for _ in 0..200 {
            let batches: Vec<BatchRecord> = (0..12)
                .map(|i| {
                    let p = if i % 2 == 0 { 0.38 } else { 0.62 };
                    let k = (0..100).filter(|_| rng.f64() < p).count();
                    BatchRecord {
                        n_shots: 100,
                        n_ones: k,
                    }
                })
                .collect();
            match assess(&RecordData::TimestampedBatches(batches), alpha) {
                CorrelationVerdict::OverdispersionDetected { .. }
                | CorrelationVerdict::BatchDriftDetected { .. } => det += 1,
                _ => {}
            }
        }
        check(
            "[P5e] batches power 較正: 交互 0.38/0.62 の検出率 ≥ 0.9",
            det >= 180,
            format!("検出 {}/200", det),
        );
    }

    // ---------------- [P6] outreach kit (d2r-v1 は不変) ----------------
    {
        let mut bad = Vec::new();
        // d2r-v1 の凍結物は不変 (MANIFEST 自体のピン)
        let d2r_man = rd("reproducer/campaigns/d2r-v1/MANIFEST.sha256");
        let h = sha256_hex(d2r_man.as_bytes());
        if h != PIN_D2R_MANIFEST {
            // 台帳行の正当な更新は許す — その場合は本ピンを PROMPT/16 §4 の手続きで更新
            bad.push(format!("d2r-v1 MANIFEST が変化 ({}… ≠ ピン)", &h[..16]));
        }
        // kit の MANIFEST 一致
        let kit_man = rd("reproducer/campaigns/d2r-v1-outreach/MANIFEST.sha256");
        let mut n = 0usize;
        for line in kit_man.lines() {
            let line = line.trim();
            if line.is_empty() {
                continue;
            }
            let mut it = line.split_whitespace();
            let (Some(want), Some(fname)) = (it.next(), it.next()) else {
                bad.push("kit MANIFEST 行不明".into());
                continue;
            };
            n += 1;
            let t = rd(&format!("reproducer/campaigns/d2r-v1-outreach/{}", fname));
            if t.is_empty() {
                bad.push(format!("kit {} が無い", fname));
            } else if sha256_hex(t.as_bytes()) != want {
                bad.push(format!("kit {} の sha256 不一致", fname));
            }
        }
        if n != 4 {
            bad.push(format!("kit は 4 成果物 (実測 {})", n));
        }
        // kit の凍結文言 (事前登録された設計と停止規則)
        let kit = rd("reproducer/campaigns/d2r-v1-outreach/OUTREACH_KIT.md");
        for needle in [
            "3 cohort × 各 10 件 = 合計 30 件",
            "failed_at_current_burden",
            "v1 の条件は変更しない",
            "数値 kernel・実装骨格・翻訳コードは渡さない",
            "予定・シミュレーション・水増しの記載は禁止",
        ] {
            if !kit.contains(needle) {
                bad.push(format!("OUTREACH_KIT: 「{}」が無い", needle));
            }
        }
        // kit 内にコードなし (campaign と同じ規律)
        if let Ok(dir) = fs::read_dir(format!("{}/reproducer/campaigns/d2r-v1-outreach", root)) {
            for e in dir.filter_map(|e| e.ok()) {
                if let Some(ext) = e.path().extension().and_then(|x| x.to_str()) {
                    if ["rs", "py", "c", "cpp", "js", "jl", "f90"].contains(&ext) {
                        bad.push(format!("kit にコード {:?} が混入", e.path()));
                    }
                }
            }
        }
        // funnel は実記録のみ (単調性は v310 が検査 — ここでは段の実在)
        let funnel = rd("reproducer/campaigns/d2r-v1/REPLICATION_FUNNEL.yml");
        if !funnel.contains("preregistered:") || !funnel.contains("not_instrumented") {
            bad.push("funnel の段/非計測宣言が無い".into());
        }
        check(
            "[P6] outreach kit: d2r-v1 不変ピン・kit MANIFEST 4 本・停止規則の事前登録・コードなし",
            bad.is_empty(),
            if bad.is_empty() {
                "配布物は凍結・d2r-v1 は逐語不変・30 件 0 なら failed_at_current_burden".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [P7] FollowUp 状態の主リポジトリ側整合 ----------------
    {
        let rep = rd("replications.yml");
        let mut bad = Vec::new();
        for needle in [
            "external_replications: 0",
            "id: REP-001",
            "status: \"partially_replicated\"",
            "human_operator: shared",
            "different_author: false",
            "replication_kind: cross_model_clean_room",
        ] {
            if !rep.contains(needle) {
                bad.push(format!("replications.yml: 「{}」が無い", needle));
            }
        }
        check(
            "[P7] FollowUp 整合: REP-001 = partially_replicated / shared operator / external 0 維持",
            bad.is_empty(),
            if bad.is_empty() {
                "cross-model clean-room は organizationally external に数えない".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "全検査 PASS — RECORD v2 は取得前に凍結・相関粒度は型で分離・d2r-v1 は不変".to_string()
        } else {
            format!("FAIL {} 件", nfail)
        }
    );
}
