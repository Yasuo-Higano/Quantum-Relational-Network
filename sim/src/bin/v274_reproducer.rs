//! v27.4 外部再現単位 (reproducer/) と再現台帳 (replications.yml) の機械検証
//! (PROMPT/10 §6)
//!
//! 「独立外部再現 0」を正直に維持しながら、その 0 を破るための単位を外部へ固定する:
//!   [P0] reproducer/ 一式 (SPEC / INPUTS×3 / EXPECTED_SCHEMA / TOLERANCES /
//!        CLAIMS / NO_SHARED_CODE) + replications.yml の実在
//!   [P1] JSON 4 ファイルの構文検証 (自前の最小 JSON パーサ — 外部クレート不使用)
//!   [P2] TOLERANCES ↔ 凍結一次ソースの一致 — 単位 A の SHA-256 =
//!        certificates/v62_sha256.txt (R1_v31)・単位 B の区間 =
//!        results/v252_bz_certificate.json (v25.2 凍結値)・単位 C の oracle = −1/10
//!   [P3] 単位 A の領域・期待解 ↔ certificates (v62_domains / v62_solutions) の一致
//!   [P4] replications.yml — 6 条件の定義・entries の整合・
//!        external_replications = 0 = claims.yml の external_replication 件数
//!   [P5] CLAIMS.md の昇格対象 id が claims.yml に実在
//!   [P6] クリーンルーム条項 (同一 AI は独立でない等) の文書アンカー
//!
//! 本監査はスイートの常時実行層 (ALWAYS_RUN)。独立外部再現が成立するときは、
//! 台帳と本監査の期待値を同一コミットで意識的に更新する (CLAIMS.md の手続き)。

use std::fs;
use std::path::Path;

// ---------------------------------------------------------------- 最小 JSON パーサ (構文検証)

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
    fn peek(&mut self) -> Option<u8> {
        self.ws();
        self.b.get(self.i).copied()
    }
    fn eat(&mut self, c: u8) -> Result<(), String> {
        self.ws();
        if self.b.get(self.i) == Some(&c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("位置 {}: '{}' が必要", self.i, c as char))
        }
    }
    fn string(&mut self) -> Result<(), String> {
        self.eat(b'"')?;
        while let Some(&c) = self.b.get(self.i) {
            self.i += 1;
            match c {
                b'"' => return Ok(()),
                b'\\' => {
                    self.i += 1;
                }
                _ => {}
            }
        }
        Err("文字列が閉じない".into())
    }
    fn value(&mut self) -> Result<(), String> {
        match self.peek() {
            Some(b'{') => {
                self.eat(b'{')?;
                if self.peek() == Some(b'}') {
                    return self.eat(b'}');
                }
                loop {
                    self.string()?;
                    self.eat(b':')?;
                    self.value()?;
                    match self.peek() {
                        Some(b',') => {
                            self.eat(b',')?;
                        }
                        _ => return self.eat(b'}'),
                    }
                }
            }
            Some(b'[') => {
                self.eat(b'[')?;
                if self.peek() == Some(b']') {
                    return self.eat(b']');
                }
                loop {
                    self.value()?;
                    match self.peek() {
                        Some(b',') => {
                            self.eat(b',')?;
                        }
                        _ => return self.eat(b']'),
                    }
                }
            }
            Some(b'"') => self.string(),
            Some(_) => {
                // number / true / false / null
                let start = self.i;
                while let Some(&c) = self.b.get(self.i) {
                    if (c as char).is_whitespace() || c == b',' || c == b'}' || c == b']' {
                        break;
                    }
                    self.i += 1;
                }
                let tok = std::str::from_utf8(&self.b[start..self.i]).unwrap_or("");
                if tok == "true" || tok == "false" || tok == "null" || tok.parse::<f64>().is_ok() {
                    Ok(())
                } else {
                    Err(format!("位置 {}: 不正なトークン '{}'", start, tok))
                }
            }
            None => Err("空の JSON".into()),
        }
    }
}

fn json_ok(text: &str) -> Result<(), String> {
    let mut p = Jp {
        b: text.as_bytes(),
        i: 0,
    };
    p.value()?;
    p.ws();
    if p.i != p.b.len() {
        return Err(format!("位置 {}: 末尾に余分な内容", p.i));
    }
    Ok(())
}

/// JSON 中の `"key": <num>` を拾う (v252_manifest と同型)
fn num_after(txt: &str, key: &str) -> Option<f64> {
    let pat = format!("\"{}\":", key);
    let p = txt.find(&pat)? + pat.len();
    let rest = &txt[p..];
    let end = rest.find(|c: char| c == ',' || c == '}' || c == '\n')?;
    rest[..end].trim().parse().ok()
}

/// TOLERANCES.yml (平坦 key: value) から値を拾う
fn tol_value(txt: &str, key: &str) -> Option<String> {
    for line in txt.lines() {
        if let Some(rest) = line.strip_prefix(&format!("{}:", key)) {
            let v = rest.trim().trim_matches('"');
            return Some(v.to_string());
        }
    }
    None
}

fn main() {
    uft_sim::self_test();
    println!(
        "=== v27.4 外部再現単位の機械検証 — reproducer/ と replications.yml (PROMPT/10 §6) ===\n"
    );
    let root = if Path::new("replications.yml").exists() {
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

    // ---- [P0] 実在 ----
    const FILES: [&str; 9] = [
        "reproducer/SPEC.md",
        "reproducer/INPUTS/unit_a_domain.json",
        "reproducer/INPUTS/unit_b_params.json",
        "reproducer/INPUTS/unit_c_params.json",
        "reproducer/EXPECTED_SCHEMA.json",
        "reproducer/TOLERANCES.yml",
        "reproducer/CLAIMS.md",
        "reproducer/NO_SHARED_CODE.md",
        "replications.yml",
    ];
    {
        let missing: Vec<&str> = FILES
            .iter()
            .filter(|f| !Path::new(&format!("{}/{}", root, f)).exists())
            .copied()
            .collect();
        check(
            "[P0] reproducer/ 一式 + replications.yml の実在",
            missing.is_empty(),
            if missing.is_empty() {
                format!("{} ファイル", FILES.len())
            } else {
                format!("欠落: {:?}", missing)
            },
        );
    }

    // ---- [P1] JSON 構文 ----
    {
        let mut bad = Vec::new();
        for f in [
            "reproducer/INPUTS/unit_a_domain.json",
            "reproducer/INPUTS/unit_b_params.json",
            "reproducer/INPUTS/unit_c_params.json",
            "reproducer/EXPECTED_SCHEMA.json",
        ] {
            match rd(f) {
                Err(_) => bad.push(format!("{} が読めない", f)),
                Ok(t) => {
                    if let Err(e) = json_ok(&t) {
                        bad.push(format!("{}: {}", f, e));
                    }
                }
            }
        }
        check(
            "[P1] JSON 4 ファイルの構文 (最小パーサ)",
            bad.is_empty(),
            if bad.is_empty() {
                "全て well-formed".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [P2] TOLERANCES ↔ 凍結一次ソース ----
    {
        let tol = rd("reproducer/TOLERANCES.yml").unwrap_or_default();
        let mut bad = Vec::new();
        // 単位 A: certificates/v62_sha256.txt の R1_v31 行
        let cert_sha = rd("certificates/v62_sha256.txt").unwrap_or_default();
        let r1_hash = cert_sha
            .lines()
            .find(|l| l.contains("R1_v31"))
            .and_then(|l| l.split_whitespace().next())
            .unwrap_or("")
            .to_string();
        if tol_value(&tol, "unit_a_sha256").as_deref() != Some(r1_hash.as_str()) {
            bad.push("unit_a_sha256 ≠ certificates/v62_sha256.txt (R1_v31)".to_string());
        }
        // 単位 B: v25.2 凍結証明書の区間
        let cert = rd("results/v252_bz_certificate.json").unwrap_or_default();
        for (k, jk) in [
            ("unit_b_lambda_x_lo", "lambda_x_lo"),
            ("unit_b_lambda_x_hi", "lambda_x_hi"),
            ("unit_b_lambda_perp_lo", "lambda_perp_lo"),
            ("unit_b_lambda_perp_hi", "lambda_perp_hi"),
        ] {
            let t: Option<f64> = tol_value(&tol, k).and_then(|v| v.parse().ok());
            let c = num_after(&cert, jk);
            if t.is_none() || c.is_none() || t != c {
                bad.push(format!("{} ≠ v252 証明書 {} ({:?} vs {:?})", k, jk, t, c));
            }
        }
        // 単位 C: oracle = −1/10・PRED-016 バー
        let ok_c = tol_value(&tol, "unit_c_oracle_16pi2_A").as_deref() == Some("-0.1")
            && tol_value(&tol, "unit_c_ratio_central_tol").as_deref() == Some("0.01")
            && tol_value(&tol, "unit_c_ratio_systematic_tol").as_deref() == Some("0.005");
        if !ok_c {
            bad.push("unit_c の oracle/バーが凍結値 (−0.1, 1%, 0.5%) と不一致".to_string());
        }
        check(
            "[P2] TOLERANCES ↔ 凍結一次ソース (v62 証明書・v25.2 凍結値・oracle −1/10・PRED-016 バー)",
            bad.is_empty(),
            if bad.is_empty() {
                "外部への約束と内部の凍結が同一の数".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [P3] 単位 A の領域と期待解 ↔ certificates ----
    {
        let dom_in = rd("reproducer/INPUTS/unit_a_domain.json").unwrap_or_default();
        let dom_src = rd("certificates/v62_domains.json").unwrap_or_default();
        let sol_src = rd("certificates/v62_solutions.json").unwrap_or_default();
        let tol = rd("reproducer/TOLERANCES.yml").unwrap_or_default();
        let mut bad = Vec::new();
        for anchor in [
            "\"hypercharge_6y_max\": 9",
            "\"max_multiplets\": 5",
            "\"max_components\": 15",
        ] {
            if !dom_in.contains(anchor) {
                bad.push(format!("INPUTS に {} が無い", anchor));
            }
            if !dom_src.contains(anchor) {
                bad.push(format!("certificates に {} が無い", anchor));
            }
        }
        // 期待多重項 = certificates の R1_v31 解と一致
        let expected = tol_value(&tol, "unit_a_multiplets").unwrap_or_default();
        for m in expected.split_whitespace() {
            if !sol_src.contains(&format!("\"{}\"", m)) {
                bad.push(format!("期待多重項 {} が v62_solutions.json に無い", m));
            }
        }
        if !sol_src.contains("\"run\": \"R1_v31\",\n    \"n_solutions\": 1,") {
            bad.push("R1_v31 の n_solutions = 1 が確認できない".to_string());
        }
        check(
            "[P3] 単位 A の領域 (D₁) と期待解 (SM 15 成分) ↔ certificates の一致",
            bad.is_empty(),
            if bad.is_empty() {
                "領域も答えも証明書から転写されている (転記ミスの機械検出)".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [P4] replications.yml と external_replication = 0 の整合 ----
    {
        let rep = rd("replications.yml").unwrap_or_default();
        let claims = rd("claims.yml").unwrap_or_default();
        let mut bad = Vec::new();
        for c in [
            "different_author",
            "independent_repository",
            "no_shared_numerical_kernel",
            "protocol_frozen_before_run",
            "commit_hash_recorded",
            "result_including_failures_public",
        ] {
            if !rep.contains(c) {
                bad.push(format!("条件 {} が未定義", c));
            }
        }
        let declared: usize = rep
            .lines()
            .find_map(|l| l.strip_prefix("external_replications:"))
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(999);
        // entries の passed × 6 条件充足の数
        let n_entries = rep.matches("\n- id: REP-").count();
        let n_claims_ext = claims
            .matches("evidence_kind: external_replication")
            .count();
        if declared != n_claims_ext {
            bad.push(format!(
                "external_replications ({}) ≠ claims の external_replication 件数 ({})",
                declared, n_claims_ext
            ));
        }
        if n_entries == 0 && declared != 0 {
            bad.push("entries が空なのに計数が 0 でない".to_string());
        }
        check(
            "[P4] replications.yml — 6 条件定義・計数 = claims 台帳 (現在どちらも 0)",
            bad.is_empty(),
            if bad.is_empty() {
                format!(
                    "external_replications = {} / entries = {} / claims 側 = {}",
                    declared, n_entries, n_claims_ext
                )
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [P5] CLAIMS.md の昇格対象 id の実在 ----
    {
        let cl_md = rd("reproducer/CLAIMS.md").unwrap_or_default();
        let claims = rd("claims.yml").unwrap_or_default();
        let mut ids = Vec::new();
        let bytes = cl_md.as_bytes();
        let mut i = 0;
        while let Some(p) = cl_md[i..].find("QRN-") {
            let s = i + p;
            let mut e = s;
            while e < bytes.len() && (bytes[e].is_ascii_alphanumeric() || bytes[e] == b'-') {
                e += 1;
            }
            ids.push(cl_md[s..e].to_string());
            i = e;
        }
        ids.sort();
        ids.dedup();
        let missing: Vec<String> = ids
            .iter()
            .filter(|id| !claims.contains(&format!("- id: {}", id)))
            .cloned()
            .collect();
        check(
            "[P5] CLAIMS.md の昇格対象 id が claims.yml に実在",
            missing.is_empty() && !ids.is_empty(),
            if missing.is_empty() {
                format!("{} id を確認", ids.len())
            } else {
                format!("不在: {:?}", missing)
            },
        );
    }

    // ---- [P6] クリーンルーム条項の文書アンカー ----
    {
        let nsc = rd("reproducer/NO_SHARED_CODE.md").unwrap_or_default();
        let spec = rd("reproducer/SPEC.md").unwrap_or_default();
        let mut bad = Vec::new();
        for (f, t, needle) in [
            ("NO_SHARED_CODE.md", &nsc, "同一 AI"),
            ("NO_SHARED_CODE.md", &nsc, "algorithmic diversity"),
            ("SPEC.md", &spec, "失敗・不一致も同じ形式で提出する"),
            ("SPEC.md", &spec, "algorithmic diversity"),
        ] {
            if !t.contains(needle) {
                bad.push(format!("{}: 「{}」が無い", f, needle));
            }
        }
        check(
            "[P6] クリーンルーム条項 (同一 AI 非独立・失敗も公開) の文書アンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "外部への約束が文書に凍結されている".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "外部再現単位は固定された — 独立外部再現 0 は「未達の目標」として機械監査下にある"
        } else {
            "**再現単位の破れ** — reproducer/ と台帳を修正せよ"
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
