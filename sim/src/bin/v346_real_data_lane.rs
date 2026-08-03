//! v34.6 real-data lane — synthetic と実測の構造的分離と受け皿 (PROMPT/15 §7)
//!
//! 背景: v34.5 までの robust 器械は synthetic lane (登録契約 = iid Bernoulli の
//! 決定的代表カウント) で走った。**synthetic shot noise を「実測ノイズ」と呼ばない**
//! (PROMPT/15 非交渉 6) を型と受け皿で守るのが本版:
//!   - `DataProvenance` (finite_data.rs): SyntheticCoverage ↛ RecordedExperimental
//!   - `reproducer/real_data/`: 実装置の repeated finite-shot record の受け皿
//!     (RECORD.schema.json・事前登録コミットメント・vendor topology commitment・
//!     drift gate・台帳 REAL_DATA_LEDGER.yml [数は実記録のみ])
//!   - `reproducer/D2R_PACKET.md`: D2-R 実施者向けの配布パケット (campaign layer は
//!     凍結のまま — PROMPT/15 §8「schema の改良ではなく配布と実施者の獲得」)
//!
//! 検証 (全て [PASS]/[FAIL] 内蔵):
//!   [L0] lane 型の存在と分離 (DataProvenance — 変換不在の文書アンカー)
//!   [L1] スキーマとレーン文書の規範アンカー (二 lane 表・OutOfDomain 正答・
//!        v35.0 完成条件 3 択・捏造禁止)
//!   [L2] fixture の構造検証 — 適合 2 (stationary/drifting)・不適合 2
//!        (synthetic の提出 / 事前登録欠落) が正しく分別される
//!   [L3] 事前登録コミットメントの機械照合 — sha256(開示予測) = commitment
//!   [L4] drift gate — split-half Clopper–Pearson 区間の disjoint 検査:
//!        stationary (6/20 vs 6/20) は通過・drifting (2/20 vs 14/20) は
//!        OutOfDomain (iid 契約の破れの正検出 — 禁止変換 25/29 の運用形)
//!   [L5] 台帳の正直さ — recorded_runs = 0 = entries 件数・fixture は数えない
//!   [L6] 配布パケットの規範アンカー (一件で足りる・6 条件・OCS-1.0 参照)
//!
//! 実記録の受理・externally operated run は本バイナリでは作れない (それは外部の
//! 人間の仕事) — 受け皿と採点規則を凍結し、台帳 0 を正直に監査する。

use std::fs;
use std::path::Path;
use uft_sim::finite_data::{cp_interval, DataProvenance};
use uft_sim::{self_test, sha256_hex};

// ---------------------------------------------------------------- 最小 JSON パーサ

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
            fields.iter().find(|(kk, _)| kk == k).map(|(_, v)| v)
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
        match self.b.get(self.i) {
            Some(b'{') => {
                self.i += 1;
                let mut fields = Vec::new();
                self.ws();
                if self.b.get(self.i) == Some(&b'}') {
                    self.i += 1;
                    return Ok(Json::Obj(fields));
                }
                loop {
                    self.ws();
                    let key = match self.parse()? {
                        Json::Str(s) => s,
                        _ => return Err("キーが文字列でない".into()),
                    };
                    self.ws();
                    if self.b.get(self.i) != Some(&b':') {
                        return Err(format!("位置 {}: ':' が必要", self.i));
                    }
                    self.i += 1;
                    let v = self.parse()?;
                    fields.push((key, v));
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => {
                            self.i += 1;
                        }
                        Some(b'}') => {
                            self.i += 1;
                            return Ok(Json::Obj(fields));
                        }
                        _ => return Err(format!("位置 {}: ',' か '}}' が必要", self.i)),
                    }
                }
            }
            Some(b'[') => {
                self.i += 1;
                let mut items = Vec::new();
                self.ws();
                if self.b.get(self.i) == Some(&b']') {
                    self.i += 1;
                    return Ok(Json::Arr(items));
                }
                loop {
                    items.push(self.parse()?);
                    self.ws();
                    match self.b.get(self.i) {
                        Some(b',') => {
                            self.i += 1;
                        }
                        Some(b']') => {
                            self.i += 1;
                            return Ok(Json::Arr(items));
                        }
                        _ => return Err(format!("位置 {}: ',' か ']' が必要", self.i)),
                    }
                }
            }
            Some(b'"') => {
                self.i += 1;
                let mut s = String::new();
                while let Some(&c) = self.b.get(self.i) {
                    self.i += 1;
                    match c {
                        b'"' => return Ok(Json::Str(s)),
                        b'\\' => {
                            if let Some(&e) = self.b.get(self.i) {
                                self.i += 1;
                                s.push(match e {
                                    b'n' => '\n',
                                    b't' => '\t',
                                    other => other as char,
                                });
                            }
                        }
                        _ => {
                            // UTF-8 継続バイトも素通し (キー照合は ASCII)
                            s.push(c as char);
                        }
                    }
                }
                Err("文字列が閉じない".into())
            }
            Some(_) => {
                let start = self.i;
                while let Some(&c) = self.b.get(self.i) {
                    if (c as char).is_whitespace() || c == b',' || c == b'}' || c == b']' {
                        break;
                    }
                    self.i += 1;
                }
                let tok = std::str::from_utf8(&self.b[start..self.i]).unwrap_or("");
                match tok {
                    "true" => Ok(Json::Bool(true)),
                    "false" => Ok(Json::Bool(false)),
                    "null" => Ok(Json::Null),
                    _ => tok
                        .parse::<f64>()
                        .map(Json::Num)
                        .map_err(|_| format!("不正なトークン '{}'", tok)),
                }
            }
            None => Err("空の JSON".into()),
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

// ---------------------------------------------------------------- 構造検証と drift gate

fn is_sha256_hex(s: &str) -> bool {
    s.len() == 64 && s.chars().all(|c| c.is_ascii_hexdigit() && !c.is_ascii_uppercase())
}

/// RECORD.schema.json の必須意味論 (スキーマの器械化) — Ok(理由なし) / Err(理由)
fn validate_record(j: &Json) -> Result<(), String> {
    let prov = j
        .get("provenance")
        .and_then(|v| v.as_str())
        .ok_or("provenance 欠落")?;
    if prov != DataProvenance::RecordedExperimental.as_str() {
        return Err(format!(
            "provenance '{}' は本 lane に提出できない (lane 分離)",
            prov
        ));
    }
    for k in ["device", "operator", "channels", "preregistration"] {
        if j.get(k).is_none() {
            return Err(format!("{} 欠落", k));
        }
    }
    let prereg = j.get("preregistration").unwrap();
    let commit = prereg
        .get("commitment_sha256")
        .and_then(|v| v.as_str())
        .ok_or("preregistration.commitment_sha256 欠落")?;
    if !is_sha256_hex(commit) {
        return Err("commitment_sha256 が sha256 hex でない".into());
    }
    let topo = j
        .get("vendor_topology_commitment")
        .and_then(|v| v.as_str())
        .ok_or("vendor_topology_commitment 欠落")?;
    if !is_sha256_hex(topo) {
        return Err("vendor_topology_commitment が sha256 hex でない".into());
    }
    match j.get("results_including_failures_public") {
        Some(Json::Bool(true)) => {}
        _ => return Err("results_including_failures_public = true が必要".into()),
    }
    // channels: shots ∈ {0,1}
    if let Some(Json::Arr(chs)) = j.get("channels") {
        if chs.is_empty() {
            return Err("channels が空".into());
        }
        for ch in chs {
            let shots = match ch.get("shots") {
                Some(Json::Arr(a)) if !a.is_empty() => a,
                _ => return Err("channel.shots 欠落/空".into()),
            };
            for s in shots {
                match s {
                    Json::Num(x) if *x == 0.0 || *x == 1.0 => {}
                    _ => return Err("shots に {0,1} 以外".into()),
                }
            }
        }
    } else {
        return Err("channels が配列でない".into());
    }
    Ok(())
}

/// drift gate: 各チャネルの split-half CP 区間 (Bonferroni α/(2m)) が disjoint なら
/// OutOfDomain (iid 契約の破れ)。Ok(true) = 通過 / Ok(false) = OutOfDomain。
fn drift_gate(j: &Json, alpha: f64) -> Result<bool, String> {
    let chs = match j.get("channels") {
        Some(Json::Arr(a)) => a,
        _ => return Err("channels 不在".into()),
    };
    let m = chs.len();
    for ch in chs {
        let shots: Vec<u8> = match ch.get("shots") {
            Some(Json::Arr(a)) => a
                .iter()
                .map(|v| if let Json::Num(x) = v { *x as u8 } else { 0 })
                .collect(),
            _ => return Err("shots 不在".into()),
        };
        let n = shots.len();
        let half = n / 2;
        let k1 = shots[..half].iter().filter(|&&s| s == 1).count();
        let k2 = shots[half..].iter().filter(|&&s| s == 1).count();
        let a_eff = alpha / (2.0 * m as f64);
        let (lo1, hi1) = cp_interval(k1, half, a_eff);
        let (lo2, hi2) = cp_interval(k2, n - half, a_eff);
        if hi1 < lo2 || hi2 < lo1 {
            return Ok(false); // disjoint — drift 検出
        }
    }
    Ok(true)
}

fn main() {
    self_test();
    println!("=== v34.6 real-data lane — synthetic と実測の構造的分離 (PROMPT/15 §7) ===\n");
    let root = if Path::new("reproducer/real_data/RECORDED_LANE.md").exists() {
        "."
    } else {
        ".."
    };
    let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p)).unwrap_or_default();
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

    // ---------------- [L0] lane 型 ----------------
    {
        let s = DataProvenance::SyntheticCoverage;
        let r = DataProvenance::RecordedExperimental;
        let src = rd("sim/src/finite_data.rs");
        check(
            "[L0] DataProvenance 型 (synthetic_coverage / recorded_experimental) — 変換不在の宣言",
            s != r
                && s.as_str() == "synthetic_coverage"
                && r.as_str() == "recorded_experimental"
                && src.contains("二つの間に変換は存在しない"),
            format!("{} ≠ {}", s.as_str(), r.as_str()),
        );
    }

    // ---------------- [L1] 文書アンカー ----------------
    {
        let lane = rd("reproducer/real_data/RECORDED_LANE.md");
        let mut bad = Vec::new();
        for a in [
            "Synthetic coverage lane",
            "Recorded experimental lane",
            "synthetic ↛ experimental",
            "OutOfDomain",
            "externally operated D2-R report",
            "未使用応答チャネル予測",
            "新しい厳密 no-go",
            "instrumental closure",
            "実記録のみ",
        ] {
            if !lane.contains(a) {
                bad.push(a);
            }
        }
        let schema = rd("reproducer/real_data/RECORD.schema.json");
        for a in [
            "recorded_experimental",
            "commitment_sha256",
            "vendor_topology_commitment",
            "results_including_failures_public",
        ] {
            if !schema.contains(a) {
                bad.push(a);
            }
        }
        check(
            "[L1] レーン文書とスキーマの規範アンカー (二 lane 分離・v35.0 完成条件 3 択・捏造禁止)",
            bad.is_empty(),
            if bad.is_empty() {
                "13 アンカー".into()
            } else {
                format!("欠落 {:?}", bad)
            },
        );
    }

    // ---------------- [L2] fixture の分別 ----------------
    let fixtures_dir = "reproducer/real_data/FIXTURES";
    let stationary = parse_json(&rd(&format!("{}/example_stationary.json", fixtures_dir)));
    let drifting = parse_json(&rd(&format!("{}/example_drifting.json", fixtures_dir)));
    {
        let bad_prov = parse_json(&rd(&format!(
            "{}/invalid_provenance_synthetic.json",
            fixtures_dir
        )));
        let bad_prereg = parse_json(&rd(&format!(
            "{}/invalid_missing_prereg.json",
            fixtures_dir
        )));
        let ok_pos = matches!(&stationary, Ok(j) if validate_record(j).is_ok())
            && matches!(&drifting, Ok(j) if validate_record(j).is_ok());
        let neg1 = matches!(&bad_prov, Ok(j) if validate_record(j).is_err());
        let neg2 = matches!(&bad_prereg, Ok(j) if validate_record(j).is_err());
        let neg1_reason = if let Ok(j) = &bad_prov {
            validate_record(j).err().unwrap_or_default()
        } else {
            String::new()
        };
        check(
            "[L2] fixture 分別: 適合 2 (stationary/drifting)・不適合 2 (synthetic 提出/事前登録欠落)",
            ok_pos && neg1 && neg2,
            format!("synthetic 提出の拒否理由: {}", neg1_reason),
        );
    }

    // ---------------- [L3] 事前登録コミットメントの照合 ----------------
    {
        let ok = if let Ok(j) = &stationary {
            let prereg = j.get("preregistration").unwrap();
            let commit = prereg
                .get("commitment_sha256")
                .and_then(|v| v.as_str())
                .unwrap_or("");
            let disclosed = prereg
                .get("disclosed_prediction")
                .and_then(|v| v.as_str())
                .unwrap_or("");
            sha256_hex(disclosed.as_bytes()) == commit
        } else {
            false
        };
        check(
            "[L3] 事前登録: sha256(開示予測) = commitment (HOLD の SECRET 機構と同型)",
            ok,
            "コミット → データ → 開示 の順序が hash で検証可能".into(),
        );
    }

    // ---------------- [L4] drift gate ----------------
    {
        let g_stat = stationary.as_ref().ok().and_then(|j| drift_gate(j, 0.05).ok());
        let g_drift = drifting.as_ref().ok().and_then(|j| drift_gate(j, 0.05).ok());
        check(
            "[L4] drift gate: stationary (6/20 vs 6/20) 通過・drifting (2/20 vs 14/20) は OutOfDomain",
            g_stat == Some(true) && g_drift == Some(false),
            "split-half CP 区間の disjoint = iid 契約の破れの正検出 (禁止変換 25/29 の運用形)".into(),
        );
    }

    // ---------------- [L5] 台帳の正直さ ----------------
    {
        let ledger = rd("reproducer/real_data/REAL_DATA_LEDGER.yml");
        let n_entries = ledger.matches("\n- id: RDR-").count();
        let zero = ["recorded_runs: 0", "externally_operated_runs: 0", "preregistered_prediction_hits: 0"]
            .iter()
            .all(|a| ledger.contains(a));
        let fixture_note = rd(&format!("{}/README.md", fixtures_dir));
        check(
            "[L5] 台帳: recorded_runs = 0 = entries 件数・fixture は実記録に数えない (捏造禁止)",
            zero && n_entries == 0
                && ledger.contains("実記録のみ")
                && fixture_note.contains("実記録ではない"),
            format!("entries = {} / 全計数 0", n_entries),
        );
    }

    // ---------------- [L6] 配布パケット ----------------
    {
        let packet = rd("reproducer/D2R_PACKET.md");
        let mut bad = Vec::new();
        for a in [
            "一件で足りる",
            "different_author",
            "no_shared_numerical_kernel",
            "protocol_frozen_before_run",
            "operational-core-spec.md",
            "失敗・不一致の報告も",
            "結果に合わせず版を上げる",
        ] {
            if !packet.contains(a) {
                bad.push(a);
            }
        }
        check(
            "[L6] D2R_PACKET (配布用): 約束の凍結・6 条件・OCS-1.0 参照 — campaign layer は不変",
            bad.is_empty(),
            if bad.is_empty() {
                "配布パケット完備 (schema 改良ではなく配布と実施者の獲得へ)".into()
            } else {
                format!("欠落 {:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "実測と synthetic が型と受け皿で分離された — 実記録の受理・externally operated\n       run は外部の人間の仕事であり、本版は受け皿と採点規則の凍結と、台帳 0 の\n       正直な監査を提供する。v35.0 の科学的完成条件 (外部報告 / 実データ予測 /\n       新 no-go のいずれか) はこの分離の上で判定される。"
        } else {
            "**real-data lane の破れ** — 分離・fixture・台帳を修復せよ"
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
