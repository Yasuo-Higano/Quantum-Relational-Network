//! v15.1 監査版: 主張依存グラフ (claims.graph.yml) の機械検証
//!
//! 第十六期の最初の道具。claims.yml (v6.1) が「主張の等級と証拠」を機械可読に
//! したのに続き、本監査は「主張の依存構造」を機械可読にする:
//!   [1] 4 台帳 (claims / graph / assumptions / falsifiers) の書式解析
//!   [2] id の一意性と命名規則 (ASM-*, FAL-*)
//!   [3] graph が claims.yml と 1:1 対応 (過不足なし)
//!   [4] 全参照の実在 (deps → claims / asm → assumptions / fal → falsifiers)
//!   [5] deps の非循環性 (DAG) と深さの計測
//!   [6] 等級単調性: rank(dep) ≤ rank(claim)。順位は C0 < C1 < C2 < {C3,C4} < C5
//!       (機構 C3 と現象論 C4 は横並び)。C2 の計算定理が現象論・解釈に依存しない
//!       ことを機械的に保証する。C0 は依存を持たない (グラフの根)。
//!   [7] 被覆: C1 以上は反証条件 ≥1 / C2 は仮定 ≥1 / 孤児の仮定・反証条件なし /
//!       棄却済み仮定は棄却した主張を指す
//!   [8] Lean 定理数の照合: proofs/*.lean の theorem 宣言を機械計数し、
//!       期待値表・最新統合文書・README の記述と照合 (v15.0 の「9 本」誤記の再発防止)
//!   [9] README の網羅性: 全文書・全バイナリへの言及、最新版の明示
//!   [10] 生成物 (claims.graph.json / evidence_matrix.md) が committed 内容と一致
//!
//! 目的は「正しさの主張」ではなく「どの仮定を抜くとどの結論が落ちるか」を
//! 機械的に見えるようにすること。各仮定・反証条件の影響範囲 (直接 + 依存閉包)
//! を導出し、evidence_matrix.md に出力する。--write で生成物を書き出す。

use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::path::Path;

// ---------------------------------------------------------------- データ構造

#[derive(Default, Clone)]
struct Claim {
    id: String,
    version: String,
    level: String,
    evidence: Vec<(String, String)>,
}

#[derive(Default, Clone)]
struct GraphEntry {
    id: String,
    deps: Vec<String>,
    asm: Vec<String>,
    fal: Vec<String>,
}

#[derive(Default, Clone)]
struct Assumption {
    id: String,
    typ: String,
    scope: String,
    introduced: String,
    status: String,
    statement: String,
    falsified_by: String,
}

#[derive(Default, Clone)]
struct Falsifier {
    id: String,
    source: String,
    status: String,
    condition: String,
    effect: String,
}

// ---------------------------------------------------------------- 共通パーサ

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

/// インライン列挙 "[A, B, C]" → Vec (空 "[]" は空 Vec)
fn parse_list(v: &str) -> Result<Vec<String>, String> {
    let t = v.trim();
    if !t.starts_with('[') || !t.ends_with(']') {
        return Err(format!("インライン列挙 [..] でない: {}", t));
    }
    let inner = &t[1..t.len() - 1];
    if inner.trim().is_empty() {
        return Ok(Vec::new());
    }
    Ok(inner.split(',').map(|x| x.trim().to_string()).collect())
}

/// claims.yml (v6.1 書式) — id/version/level/evidence のみ使う
fn parse_claims(text: &str) -> Result<Vec<Claim>, String> {
    #[derive(PartialEq)]
    enum Mode {
        Top,
        Evidence,
        Inputs,
        Limitations,
    }
    let mut out: Vec<Claim> = Vec::new();
    let mut mode = Mode::Top;
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
            mode = Mode::Top;
        } else if let Some(rest) = line.strip_prefix("    ") {
            let cur = out.last_mut().ok_or(format!("{}行目: エントリ外", lno))?;
            match mode {
                Mode::Evidence => {
                    let (k, v) = split_kv(rest).ok_or(format!("{}行目: evidence 書式", lno))?;
                    cur.evidence.push((k, v));
                }
                Mode::Inputs | Mode::Limitations => {}
                Mode::Top => return Err(format!("{}行目: ブロック外の 4 スペース行", lno)),
            }
        } else if let Some(rest) = line.strip_prefix("  ") {
            let cur = out.last_mut().ok_or(format!("{}行目: エントリ外", lno))?;
            match rest {
                "evidence:" => mode = Mode::Evidence,
                "inputs:" => mode = Mode::Inputs,
                "limitations:" => mode = Mode::Limitations,
                _ => {
                    mode = Mode::Top;
                    let (k, v) = split_kv(rest).ok_or(format!("{}行目: フィールド書式", lno))?;
                    match k.as_str() {
                        "version" => cur.version = v,
                        "level" => cur.level = v,
                        "claim" | "status" => {}
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

/// claims.graph.yml
fn parse_graph(text: &str) -> Result<Vec<GraphEntry>, String> {
    let mut out: Vec<GraphEntry> = Vec::new();
    for (lno, raw) in text.lines().enumerate() {
        let lno = lno + 1;
        let line = raw.trim_end();
        if line.trim_start().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        if let Some(rest) = line.strip_prefix("- id:") {
            out.push(GraphEntry {
                id: unquote(rest),
                ..Default::default()
            });
        } else if let Some(rest) = line.strip_prefix("  ") {
            let cur = out.last_mut().ok_or(format!("{}行目: エントリ外", lno))?;
            let (k, v) = split_kv(rest).ok_or(format!("{}行目: フィールド書式", lno))?;
            let lst = parse_list(&v).map_err(|e| format!("{}行目: {}", lno, e))?;
            match k.as_str() {
                "deps" => cur.deps = lst,
                "asm" => cur.asm = lst,
                "fal" => cur.fal = lst,
                _ => return Err(format!("{}行目: 未知フィールド '{}'", lno, k)),
            }
        } else {
            return Err(format!("{}行目: 解釈できない行", lno));
        }
    }
    Ok(out)
}

/// assumptions.yml / falsifiers.yml 共通の平坦エントリパーサ
fn parse_flat(text: &str, known: &[&str]) -> Result<Vec<BTreeMap<String, String>>, String> {
    let mut out: Vec<BTreeMap<String, String>> = Vec::new();
    for (lno, raw) in text.lines().enumerate() {
        let lno = lno + 1;
        let line = raw.trim_end();
        if line.trim_start().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        if let Some(rest) = line.strip_prefix("- id:") {
            let mut m = BTreeMap::new();
            m.insert("id".to_string(), unquote(rest));
            out.push(m);
        } else if let Some(rest) = line.strip_prefix("  ") {
            let cur = out.last_mut().ok_or(format!("{}行目: エントリ外", lno))?;
            let (k, v) = split_kv(rest).ok_or(format!("{}行目: フィールド書式", lno))?;
            if !known.contains(&k.as_str()) {
                return Err(format!("{}行目: 未知フィールド '{}'", lno, k));
            }
            cur.insert(k, v);
        } else {
            return Err(format!("{}行目: 解釈できない行", lno));
        }
    }
    Ok(out)
}

fn tag_id_ok(id: &str, prefix: &str) -> bool {
    id.starts_with(prefix)
        && id.len() > prefix.len()
        && id
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit() || c == '-')
        && !id.ends_with('-')
        && !id.contains("--")
}

fn rank(level: &str) -> i32 {
    match level {
        "C0" => 0,
        "C1" => 1,
        "C2" => 2,
        "C3" | "C4" => 3, // 機構と現象論は横並び (ASM-EDGE-SEMANTICS)
        "C5" => 4,
        _ => 99,
    }
}

// ---------------------------------------------------------------- JSON 生成

fn jesc(s: &str) -> String {
    let mut o = String::new();
    for c in s.chars() {
        match c {
            '"' => o.push_str("\\\""),
            '\\' => o.push_str("\\\\"),
            '\n' => o.push_str("\\n"),
            c if (c as u32) < 0x20 => o.push_str(&format!("\\u{:04x}", c as u32)),
            c => o.push(c),
        }
    }
    o
}

fn jlist(v: &[String]) -> String {
    let items: Vec<String> = v.iter().map(|s| format!("\"{}\"", jesc(s))).collect();
    format!("[{}]", items.join(", "))
}

// ---------------------------------------------------------------- main

fn main() {
    println!("=== v15.1 監査版: 主張依存グラフの機械検証 ===\n");
    let write_mode = std::env::args().any(|a| a == "--write");
    let root = if Path::new("claims.yml").exists() {
        "."
    } else if Path::new("../claims.yml").exists() {
        ".."
    } else {
        println!("claims.yml が見つからない  [FAIL]");
        std::process::exit(1);
    };
    let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p));
    let mut nfail = 0;

    // [1] 4 台帳の解析
    let claims_text = rd("claims.yml").expect("claims.yml 読めない");
    let graph_text = rd("claims.graph.yml").expect("claims.graph.yml 読めない");
    let asm_text = rd("assumptions.yml").expect("assumptions.yml 読めない");
    let fal_text = rd("falsifiers.yml").expect("falsifiers.yml 読めない");
    let claims = match parse_claims(&claims_text) {
        Ok(c) => c,
        Err(e) => {
            println!("[1] claims.yml: {}  [FAIL]", e);
            std::process::exit(1);
        }
    };
    let graph = match parse_graph(&graph_text) {
        Ok(g) => g,
        Err(e) => {
            println!("[1] claims.graph.yml: {}  [FAIL]", e);
            std::process::exit(1);
        }
    };
    let asms: Vec<Assumption> = match parse_flat(
        &asm_text,
        &[
            "id",
            "type",
            "scope",
            "introduced",
            "status",
            "statement",
            "note",
            "falsified_by",
        ],
    ) {
        Ok(v) => v
            .into_iter()
            .map(|m| Assumption {
                id: m.get("id").cloned().unwrap_or_default(),
                typ: m.get("type").cloned().unwrap_or_default(),
                scope: m.get("scope").cloned().unwrap_or_default(),
                introduced: m.get("introduced").cloned().unwrap_or_default(),
                status: m.get("status").cloned().unwrap_or_default(),
                statement: m.get("statement").cloned().unwrap_or_default(),
                falsified_by: m.get("falsified_by").cloned().unwrap_or_default(),
            })
            .collect(),
        Err(e) => {
            println!("[1] assumptions.yml: {}  [FAIL]", e);
            std::process::exit(1);
        }
    };
    let fals: Vec<Falsifier> = match parse_flat(
        &fal_text,
        &["id", "source", "status", "condition", "effect"],
    ) {
        Ok(v) => v
            .into_iter()
            .map(|m| Falsifier {
                id: m.get("id").cloned().unwrap_or_default(),
                source: m.get("source").cloned().unwrap_or_default(),
                status: m.get("status").cloned().unwrap_or_default(),
                condition: m.get("condition").cloned().unwrap_or_default(),
                effect: m.get("effect").cloned().unwrap_or_default(),
            })
            .collect(),
        Err(e) => {
            println!("[1] falsifiers.yml: {}  [FAIL]", e);
            std::process::exit(1);
        }
    };
    println!(
        "[1] 書式解析: 主張 {} / グラフ {} / 仮定 {} / 反証条件 {}  [PASS]",
        claims.len(),
        graph.len(),
        asms.len(),
        fals.len()
    );

    // [2] id の一意性と命名規則
    {
        let mut bad = Vec::new();
        let mut seen = BTreeSet::new();
        for a in &asms {
            if !tag_id_ok(&a.id, "ASM-") {
                bad.push(format!("{} (命名規則)", a.id));
            }
            if !seen.insert(a.id.clone()) {
                bad.push(format!("{} (重複)", a.id));
            }
            if ![
                "framework",
                "window",
                "model",
                "design",
                "data",
                "observational",
                "convention",
                "definition",
                "trust",
                "ontology",
            ]
            .contains(&a.typ.as_str())
            {
                bad.push(format!("{} (type '{}')", a.id, a.typ));
            }
            if !["global", "local"].contains(&a.scope.as_str()) {
                bad.push(format!("{} (scope '{}')", a.id, a.scope));
            }
            if !["active", "falsified", "superseded"].contains(&a.status.as_str()) {
                bad.push(format!("{} (status '{}')", a.id, a.status));
            }
            if a.statement.is_empty() || a.introduced.is_empty() {
                bad.push(format!("{} (必須フィールド欠落)", a.id));
            }
        }
        let mut seen_f = BTreeSet::new();
        for f in &fals {
            if !tag_id_ok(&f.id, "FAL-") {
                bad.push(format!("{} (命名規則)", f.id));
            }
            if !seen_f.insert(f.id.clone()) {
                bad.push(format!("{} (重複)", f.id));
            }
            if !["open", "triggered", "retired"].contains(&f.status.as_str()) {
                bad.push(format!("{} (status '{}')", f.id, f.status));
            }
            if f.condition.is_empty() || f.effect.is_empty() || f.source.is_empty() {
                bad.push(format!("{} (必須フィールド欠落)", f.id));
            }
        }
        if bad.is_empty() {
            println!("[2] id 一意性・命名・語彙: OK  [PASS]");
        } else {
            println!("[2] id 検査: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // [3] graph ↔ claims の 1:1 対応
    let claim_ids: BTreeSet<&str> = claims.iter().map(|c| c.id.as_str()).collect();
    let graph_ids: BTreeSet<&str> = graph.iter().map(|g| g.id.as_str()).collect();
    {
        let missing: Vec<&&str> = claim_ids.difference(&graph_ids).collect();
        let extra: Vec<&&str> = graph_ids.difference(&claim_ids).collect();
        let dup = graph.len() != graph_ids.len();
        if missing.is_empty() && extra.is_empty() && !dup {
            println!(
                "[3] graph ↔ claims 1:1 対応: 全 {} 主張  [PASS]",
                claims.len()
            );
        } else {
            println!(
                "[3] 対応検査: graph 欠落 {:?} / 余剰 {:?} / 重複 {}  [FAIL]",
                missing, extra, dup
            );
            nfail += 1;
        }
    }

    // [4] 全参照の実在
    let asm_ids: BTreeSet<&str> = asms.iter().map(|a| a.id.as_str()).collect();
    let fal_ids: BTreeSet<&str> = fals.iter().map(|f| f.id.as_str()).collect();
    {
        let mut bad = Vec::new();
        for g in &graph {
            for d in &g.deps {
                if !claim_ids.contains(d.as_str()) {
                    bad.push(format!("{}: deps {}", g.id, d));
                }
                if d == &g.id {
                    bad.push(format!("{}: 自己依存", g.id));
                }
            }
            for a in &g.asm {
                if !asm_ids.contains(a.as_str()) {
                    bad.push(format!("{}: asm {}", g.id, a));
                }
            }
            for f in &g.fal {
                if !fal_ids.contains(f.as_str()) {
                    bad.push(format!("{}: fal {}", g.id, f));
                }
            }
        }
        for a in &asms {
            if a.status == "falsified" && !claim_ids.contains(a.falsified_by.as_str()) {
                bad.push(format!("{}: falsified_by '{}' 不在", a.id, a.falsified_by));
            }
        }
        if bad.is_empty() {
            let ne: usize = graph.iter().map(|g| g.deps.len()).sum();
            println!("[4] 参照の実在: 依存辺 {} 本ほか全て OK  [PASS]", ne);
        } else {
            println!("[4] 参照: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // 以降のためのインデックス
    let gidx: BTreeMap<&str, &GraphEntry> = graph.iter().map(|g| (g.id.as_str(), g)).collect();
    let cidx: BTreeMap<&str, &Claim> = claims.iter().map(|c| (c.id.as_str(), c)).collect();

    // [5] DAG 非循環性と深さ
    let mut depth: BTreeMap<String, i32> = BTreeMap::new();
    {
        // 反復的トポロジカル整列 (依存が全て解決した節点から深さ確定)
        let mut remaining: BTreeSet<&str> = graph_ids.clone();
        let mut changed = true;
        while changed && !remaining.is_empty() {
            changed = false;
            let snapshot: Vec<&str> = remaining.iter().cloned().collect();
            for id in snapshot {
                let g = gidx[id];
                if g.deps.iter().all(|d| depth.contains_key(d.as_str())) {
                    let dmax = g.deps.iter().map(|d| depth[d.as_str()]).max().unwrap_or(-1);
                    depth.insert(id.to_string(), dmax + 1);
                    remaining.remove(id);
                    changed = true;
                }
            }
        }
        if remaining.is_empty() {
            let maxd = depth.values().cloned().max().unwrap_or(0);
            let deepest: Vec<&String> = depth
                .iter()
                .filter(|(_, v)| **v == maxd)
                .map(|(k, _)| k)
                .collect();
            println!(
                "[5] DAG 非循環性: OK (最大深さ {} — 例 {})  [PASS]",
                maxd,
                deepest.first().map(|s| s.as_str()).unwrap_or("-")
            );
        } else {
            println!("[5] 循環に関与する主張: {:?}  [FAIL]", remaining);
            nfail += 1;
        }
    }

    // [6] 等級単調性と C0 の根性
    {
        let mut bad = Vec::new();
        for g in &graph {
            let c = match cidx.get(g.id.as_str()) {
                Some(c) => c,
                None => continue, // [3] で検出済み
            };
            if c.level == "C0" && (!g.deps.is_empty() || !g.asm.is_empty() || !g.fal.is_empty()) {
                bad.push(format!("{}: C0 が依存/仮定/反証条件を持つ", g.id));
            }
            for d in &g.deps {
                if let Some(dc) = cidx.get(d.as_str()) {
                    if rank(&dc.level) > rank(&c.level) {
                        bad.push(format!(
                            "{} ({}) → {} ({}): 等級単調性違反",
                            g.id, c.level, d, dc.level
                        ));
                    }
                }
            }
        }
        if bad.is_empty() {
            println!("[6] 等級単調性 (C0<C1<C2<{{C3,C4}}<C5): 全依存辺 OK  [PASS]");
        } else {
            println!("[6] 等級単調性: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // [7] 被覆: 反証条件・仮定・孤児
    {
        let mut bad = Vec::new();
        for g in &graph {
            let c = match cidx.get(g.id.as_str()) {
                Some(c) => c,
                None => continue,
            };
            if c.level != "C0" && g.fal.is_empty() {
                bad.push(format!("{}: 反証条件なし", g.id));
            }
            if c.level == "C2" && g.asm.is_empty() {
                bad.push(format!("{}: C2 に仮定なし", g.id));
            }
        }
        let used_asm: BTreeSet<&str> = graph
            .iter()
            .flat_map(|g| g.asm.iter().map(|s| s.as_str()))
            .collect();
        let used_fal: BTreeSet<&str> = graph
            .iter()
            .flat_map(|g| g.fal.iter().map(|s| s.as_str()))
            .collect();
        for a in &asms {
            if a.scope == "local" && !used_asm.contains(a.id.as_str()) {
                bad.push(format!("{}: どの主張からも参照されない", a.id));
            }
        }
        for f in &fals {
            if !used_fal.contains(f.id.as_str()) {
                bad.push(format!("{}: どの主張からも参照されない", f.id));
            }
        }
        if bad.is_empty() {
            println!("[7] 被覆 (C1+ に fal / C2 に asm / 孤児なし): OK  [PASS]");
        } else {
            println!("[7] 被覆: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // 影響範囲の導出 (逆依存の閉包)
    let mut rev: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
    for g in &graph {
        for d in &g.deps {
            rev.entry(d.as_str()).or_default().push(g.id.as_str());
        }
    }
    let closure = |start: &BTreeSet<&str>| -> BTreeSet<String> {
        let mut seen: BTreeSet<String> = start.iter().map(|s| s.to_string()).collect();
        let mut stack: Vec<&str> = start.iter().cloned().collect();
        while let Some(x) = stack.pop() {
            if let Some(ups) = rev.get(x) {
                for u in ups {
                    if seen.insert(u.to_string()) {
                        stack.push(u);
                    }
                }
            }
        }
        seen
    };
    // 主張ごとの被支持閉包 (自分を除く)
    let mut dependents: BTreeMap<&str, usize> = BTreeMap::new();
    for id in &graph_ids {
        let mut s = BTreeSet::new();
        s.insert(*id);
        dependents.insert(id, closure(&s).len() - 1);
    }
    // 仮定・反証条件ごとの影響範囲
    let mut asm_direct: BTreeMap<&str, usize> = BTreeMap::new();
    let mut asm_blast: BTreeMap<&str, BTreeSet<String>> = BTreeMap::new();
    for a in &asms {
        let direct: BTreeSet<&str> = graph
            .iter()
            .filter(|g| g.asm.iter().any(|x| x == &a.id))
            .map(|g| g.id.as_str())
            .collect();
        asm_direct.insert(a.id.as_str(), direct.len());
        asm_blast.insert(a.id.as_str(), closure(&direct));
    }
    let mut fal_direct: BTreeMap<&str, usize> = BTreeMap::new();
    let mut fal_blast: BTreeMap<&str, BTreeSet<String>> = BTreeMap::new();
    for f in &fals {
        let direct: BTreeSet<&str> = graph
            .iter()
            .filter(|g| g.fal.iter().any(|x| x == &f.id))
            .map(|g| g.id.as_str())
            .collect();
        fal_direct.insert(f.id.as_str(), direct.len());
        fal_blast.insert(f.id.as_str(), closure(&direct));
    }

    // [8] Lean 定理数の照合
    let lean_expected: [(&str, usize); 5] = [
        ("Anomaly.lean", 3),
        ("AnomalyArray.lean", 2),
        ("AnomalyArrayBig.lean", 2),
        ("AnomalyBig.lean", 3),
        ("AnomalyUpstream.lean", 13),
    ];
    let lean_total_expected: usize = lean_expected.iter().map(|(_, n)| n).sum();
    {
        let mut bad = Vec::new();
        let mut actual: BTreeMap<String, usize> = BTreeMap::new();
        let dir = fs::read_dir(format!("{}/proofs", root)).expect("proofs/ が読めない");
        for e in dir.filter_map(|e| e.ok()) {
            let name = e.file_name().to_string_lossy().to_string();
            if name.ends_with(".lean") {
                let text = fs::read_to_string(e.path()).unwrap_or_default();
                let n = text
                    .lines()
                    .filter(|l| l.trim_start().starts_with("theorem "))
                    .count();
                actual.insert(name, n);
            }
        }
        for (f, n) in &lean_expected {
            match actual.get(*f) {
                Some(m) if m == n => {}
                Some(m) => bad.push(format!("{}: 期待 {} 実測 {}", f, n, m)),
                None => bad.push(format!("{}: ファイル不在", f)),
            }
        }
        for f in actual.keys() {
            if !lean_expected.iter().any(|(e, _)| e == f) {
                bad.push(format!(
                    "{}: 台帳外の Lean ファイル (期待値表に追記せよ)",
                    f
                ));
            }
        }
        let total: usize = actual.values().sum();
        // 最新統合文書と README の記述照合
        let readme = rd("README.md").unwrap_or_default();
        let v150 = rd("docs/uft-v15.0.md").unwrap_or_default();
        let want = format!("Lean 定理 {} 本", lean_total_expected);
        let wrong = "Lean 定理 9 本";
        // 現在値の照合は README のみ (v21.6 で改訂: uft-v15.0.md は当時の統合アンカーで
        // 歴史文書 — 定理数が成長する現在は照合先にしない。誤記検査は両方に残す)
        if !readme.contains(&want) {
            bad.push(format!("README.md: 「{}」の記述がない", want));
        }
        for (name, text) in [("README.md", &readme), ("docs/uft-v15.0.md", &v150)] {
            if text.contains(wrong) {
                bad.push(format!("{}: 誤記「{}」が残存", name, wrong));
            }
        }
        if bad.is_empty() {
            println!(
                "[8] Lean 定理数: {} ファイル計 {} 本 = 期待値、文書記述とも一致  [PASS]",
                lean_expected.len(),
                total
            );
        } else {
            println!("[8] Lean 定理数: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // [9] README の網羅性
    {
        let mut bad = Vec::new();
        let readme = rd("README.md").unwrap_or_default();
        // 全文書がリンクされている
        let mut latest = (0u32, 0u32);
        let dir = fs::read_dir(format!("{}/docs", root)).expect("docs/ が読めない");
        for e in dir.filter_map(|e| e.ok()) {
            let name = e.file_name().to_string_lossy().to_string();
            if !name.ends_with(".md") {
                continue;
            }
            if !readme.contains(&name) {
                bad.push(format!("README に {} への言及がない", name));
            }
            if let Some(v) = name
                .strip_prefix("uft-v")
                .and_then(|s| s.strip_suffix(".md"))
            {
                let head = v.split('-').next().unwrap_or(v); // "5.1-5.3" → "5.1"
                let mut it = head.split('.');
                if let (Some(a), Some(b)) = (it.next(), it.next()) {
                    if let (Ok(a), Ok(b)) = (a.parse::<u32>(), b.parse::<u32>()) {
                        if (a, b) > latest {
                            latest = (a, b);
                        }
                    }
                }
            }
        }
        // 全バイナリが言及されている
        let dir = fs::read_dir(format!("{}/sim/src/bin", root)).expect("bin が読めない");
        for e in dir.filter_map(|e| e.ok()) {
            let name = e.file_name().to_string_lossy().to_string();
            if let Some(stem) = name.strip_suffix(".rs") {
                if !readme.contains(stem) {
                    bad.push(format!("README にバイナリ {} への言及がない", stem));
                }
            }
        }
        // 最新版が「現在の到達点」として明示されている
        let want = format!("現在の到達点: v{}.{}", latest.0, latest.1);
        if !readme.contains(&want) {
            bad.push(format!("README に「{}」の明示がない", want));
        }
        if bad.is_empty() {
            println!(
                "[9] README 網羅性: 全文書・全バイナリ言及、到達点 v{}.{} 明示  [PASS]",
                latest.0, latest.1
            );
        } else {
            println!("[9] README 網羅性: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // ---------------- 生成物: claims.graph.json / evidence_matrix.md ----------------
    let mut json = String::new();
    {
        json.push_str("{\n");
        json.push_str(&format!(
            "  \"generated_by\": \"v151_audit\",\n  \"n_claims\": {},\n  \"n_edges\": {},\n  \"n_assumptions\": {},\n  \"n_falsifiers\": {},\n",
            claims.len(),
            graph.iter().map(|g| g.deps.len()).sum::<usize>(),
            asms.len(),
            fals.len()
        ));
        json.push_str("  \"claims\": [\n");
        for (i, c) in claims.iter().enumerate() {
            let g = gidx[c.id.as_str()];
            let supports: Vec<String> = rev
                .get(c.id.as_str())
                .map(|v| v.iter().map(|s| s.to_string()).collect())
                .unwrap_or_default();
            let ev: Vec<String> = c
                .evidence
                .iter()
                .map(|(k, v)| format!("{{\"kind\": \"{}\", \"path\": \"{}\"}}", jesc(k), jesc(v)))
                .collect();
            json.push_str(&format!(
                "    {{\"id\": \"{}\", \"version\": \"{}\", \"level\": \"{}\", \"depth\": {}, \"deps\": {}, \"asm\": {}, \"fal\": {}, \"supported_by_closure\": {}, \"supports_direct\": {}, \"evidence\": [{}]}}{}\n",
                jesc(&c.id),
                jesc(&c.version),
                jesc(&c.level),
                depth.get(&c.id).cloned().unwrap_or(-1),
                jlist(&g.deps),
                jlist(&g.asm),
                jlist(&g.fal),
                dependents.get(c.id.as_str()).cloned().unwrap_or(0),
                jlist(&supports),
                ev.join(", "),
                if i + 1 < claims.len() { "," } else { "" }
            ));
        }
        json.push_str("  ],\n  \"assumptions\": [\n");
        for (i, a) in asms.iter().enumerate() {
            json.push_str(&format!(
                "    {{\"id\": \"{}\", \"type\": \"{}\", \"scope\": \"{}\", \"introduced\": \"{}\", \"status\": \"{}\", \"direct_users\": {}, \"blast_radius\": {}, \"statement\": \"{}\"}}{}\n",
                jesc(&a.id),
                jesc(&a.typ),
                jesc(&a.scope),
                jesc(&a.introduced),
                jesc(&a.status),
                asm_direct[a.id.as_str()],
                asm_blast[a.id.as_str()].len(),
                jesc(&a.statement),
                if i + 1 < asms.len() { "," } else { "" }
            ));
        }
        json.push_str("  ],\n  \"falsifiers\": [\n");
        for (i, f) in fals.iter().enumerate() {
            json.push_str(&format!(
                "    {{\"id\": \"{}\", \"source\": \"{}\", \"status\": \"{}\", \"direct_targets\": {}, \"blast_radius\": {}, \"condition\": \"{}\", \"effect\": \"{}\"}}{}\n",
                jesc(&f.id),
                jesc(&f.source),
                jesc(&f.status),
                fal_direct[f.id.as_str()],
                fal_blast[f.id.as_str()].len(),
                jesc(&f.condition),
                jesc(&f.effect),
                if i + 1 < fals.len() { "," } else { "" }
            ));
        }
        json.push_str("  ]\n}\n");
    }

    let mut md = String::new();
    {
        md.push_str("# QRN 証拠マトリクス (evidence matrix)\n\n");
        md.push_str("**このファイルは `v151_audit --write` が生成する。手で編集しない。**\n");
        md.push_str("機械可読版は [claims.graph.json](claims.graph.json)、辺の定義は [claims.graph.yml](claims.graph.yml)。\n\n");
        md.push_str(&format!(
            "主張 {} 件 / 依存辺 {} 本 / 仮定 {} 件 / 反証条件 {} 件。\n",
            claims.len(),
            graph.iter().map(|g| g.deps.len()).sum::<usize>(),
            asms.len(),
            fals.len()
        ));
        md.push_str(
            "等級順位 C0 < C1 < C2 < {C3,C4} < C5 の単調性・非循環性は CI で機械検証される。\n\n",
        );
        md.push_str("## 主張 × 証拠・依存\n\n");
        md.push_str("凡例: 証拠 C=code R=result D=doc L=Lean。「被支持」= この主張が落ちると (依存閉包で) 落ちる主張の数。\n\n");
        md.push_str("| ID | 等級 | 版 | 証拠 | deps | asm | fal | 被支持 |\n");
        md.push_str("|---|---|---|---|---|---|---|---|\n");
        for c in &claims {
            let g = gidx[c.id.as_str()];
            let mut kinds = BTreeSet::new();
            for (k, v) in &c.evidence {
                if v.starts_with("proofs/") {
                    kinds.insert("L");
                } else {
                    kinds.insert(match k.as_str() {
                        "code" => "C",
                        "result" => "R",
                        _ => "D",
                    });
                }
            }
            let ks: Vec<&str> = kinds.into_iter().collect();
            md.push_str(&format!(
                "| {} | {} | {} | {} | {} | {} | {} | {} |\n",
                c.id,
                c.level,
                c.version,
                ks.join(""),
                g.deps.len(),
                g.asm.len(),
                g.fal.len(),
                dependents[c.id.as_str()]
            ));
        }
        md.push_str("\n## 仮定の影響範囲 — これを抜くと何が落ちるか\n\n");
        md.push_str("「直接」= この仮定を明示的に使う主張の数。「閉包」= 依存を遡って落ちる主張の総数。\n\n");
        md.push_str("| 仮定 | type | status | 直接 | 閉包 | 閉包に含まれる主張 (抜粋) |\n");
        md.push_str("|---|---|---|---|---|---|\n");
        let mut sorted_asms: Vec<&Assumption> = asms.iter().collect();
        sorted_asms.sort_by_key(|a| std::cmp::Reverse(asm_blast[a.id.as_str()].len()));
        for a in sorted_asms {
            let blast = &asm_blast[a.id.as_str()];
            let sample: Vec<&str> = blast.iter().take(4).map(|s| s.as_str()).collect();
            md.push_str(&format!(
                "| {} | {} | {} | {} | {} | {}{} |\n",
                a.id,
                a.typ,
                a.status,
                asm_direct[a.id.as_str()],
                blast.len(),
                sample.join(", "),
                if blast.len() > 4 { ", …" } else { "" }
            ));
        }
        md.push_str("\n## 反証条件の射程 — これが発火すると何が落ちるか\n\n");
        md.push_str("| 反証条件 | status | 直接 | 閉包 | 条件 (要約) |\n");
        md.push_str("|---|---|---|---|---|\n");
        let mut sorted_fals: Vec<&Falsifier> = fals.iter().collect();
        sorted_fals.sort_by_key(|f| std::cmp::Reverse(fal_blast[f.id.as_str()].len()));
        for f in sorted_fals {
            md.push_str(&format!(
                "| {} | {} | {} | {} | {} |\n",
                f.id,
                f.status,
                fal_direct[f.id.as_str()],
                fal_blast[f.id.as_str()].len(),
                f.condition
            ));
        }
        md.push('\n');
    }

    // [10] 生成物の同一性 (--write で更新)
    {
        let jpath = format!("{}/claims.graph.json", root);
        let mpath = format!("{}/evidence_matrix.md", root);
        if write_mode {
            fs::write(&jpath, &json).expect("claims.graph.json 書き込み失敗");
            fs::write(&mpath, &md).expect("evidence_matrix.md 書き込み失敗");
            println!("[10] 生成物を書き出した (--write): claims.graph.json / evidence_matrix.md");
        }
        let jdisk = fs::read_to_string(&jpath).unwrap_or_default();
        let mdisk = fs::read_to_string(&mpath).unwrap_or_default();
        if jdisk == json && mdisk == md {
            println!("[10] 生成物の同一性: committed 内容 = 再生成内容  [PASS]");
        } else {
            println!("[10] 生成物が古い — `v151_audit --write` で更新せよ  [FAIL]");
            nfail += 1;
        }
    }

    // ---------------- 集計表示 ----------------
    println!("\n---- 依存グラフの概観 ----");
    let nroots = graph.iter().filter(|g| g.deps.is_empty()).count();
    println!(
        "  節点 {} / 辺 {} / 根 (deps なし) {} / 最大深さ {}",
        graph.len(),
        graph.iter().map(|g| g.deps.len()).sum::<usize>(),
        nroots,
        depth.values().cloned().max().unwrap_or(0)
    );
    let mut top_claims: Vec<(&str, usize)> = dependents.iter().map(|(k, v)| (*k, *v)).collect();
    top_claims.sort_by_key(|(id, n)| (std::cmp::Reverse(*n), *id));
    println!("  被支持数の上位 (この主張が落ちると連鎖する数):");
    for (id, n) in top_claims.iter().take(5) {
        println!("    {:22} {:2} 主張", id, n);
    }
    let mut top_asm: Vec<(&str, usize)> = asms
        .iter()
        .map(|a| (a.id.as_str(), asm_blast[a.id.as_str()].len()))
        .collect();
    top_asm.sort_by_key(|(id, n)| (std::cmp::Reverse(*n), *id));
    println!("  仮定の影響範囲の上位 (抜くと落ちる主張の閉包):");
    for (id, n) in top_asm.iter().take(5) {
        println!("    {:22} {:2} 主張", id, n);
    }

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 依存グラフは非循環・等級単調・完全被覆 (辺の意味論は ASM-EDGE-SEMANTICS)"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
