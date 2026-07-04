//! v6.1 主張台帳 (claims.yml) の機械検証
//!
//! 第七期 (監査期) の最初の道具。リポジトリの全主張を C0-C5 の等級に分類した
//! claims.yml が、形式的に健全であることを検査する:
//!   [1] 書式 (制限 YAML) が仕様どおり解析できる
//!   [2] id が一意で命名規則 QRN-<領域>-<3桁> に従う
//!   [3] 等級が C0..C5 のいずれかで、必須フィールドが揃っている
//!   [4] 証拠 (code/result/doc) のファイルが全て実在する
//!   [5] C0 以外の全主張が限界 (limitations) を最低 1 件持つ
//!   [6] 全シミュレーションバイナリが台帳から参照されている (取りこぼしなし)
//! これにより「どこまでが導出でどこからが解釈か分からない」という批判に
//! 機械検証可能な形で答える。台帳の内容そのものの正しさは各文書と results/ が担う。

use std::collections::{BTreeMap, HashSet};
use std::fs;
use std::path::Path;

#[derive(Default)]
struct Claim {
    id: String,
    version: String,
    level: String,
    claim: String,
    status: String,
    evidence: Vec<(String, String)>,
    limitations: Vec<String>,
    inputs: Vec<(String, String)>,
    line: usize,
}

fn unquote(s: &str) -> String {
    let t = s.trim();
    if t.len() >= 2 && t.starts_with('"') && t.ends_with('"') {
        t[1..t.len() - 1].to_string()
    } else {
        t.to_string()
    }
}

/// (key, value) に分割 (最初の ':' のみで切る — 値の中の ':' は保持)
fn split_kv(s: &str) -> Option<(String, String)> {
    let idx = s.find(':')?;
    Some((s[..idx].trim().to_string(), unquote(&s[idx + 1..])))
}

fn parse(text: &str) -> Result<Vec<Claim>, String> {
    #[derive(PartialEq)]
    enum Mode {
        Top,
        Evidence,
        Inputs,
        Limitations,
    }
    let mut claims: Vec<Claim> = Vec::new();
    let mut mode = Mode::Top;
    for (lno, raw) in text.lines().enumerate() {
        let lno = lno + 1;
        let line = raw.trim_end();
        if line.trim_start().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        if let Some(rest) = line.strip_prefix("- id:") {
            let mut c = Claim::default();
            c.id = unquote(rest);
            c.line = lno;
            claims.push(c);
            mode = Mode::Top;
        } else if let Some(rest) = line.strip_prefix("    ") {
            // 4 スペース: ブロック内要素
            let cur = claims
                .last_mut()
                .ok_or(format!("{}行目: エントリ外の要素", lno))?;
            match mode {
                Mode::Limitations => {
                    let item = rest.strip_prefix("- ").ok_or(format!(
                        "{}行目: limitations の要素は '- ' で始まる必要",
                        lno
                    ))?;
                    cur.limitations.push(unquote(item));
                }
                Mode::Evidence => {
                    let (k, v) =
                        split_kv(rest).ok_or(format!("{}行目: evidence 内の書式エラー", lno))?;
                    cur.evidence.push((k, v));
                }
                Mode::Inputs => {
                    let (k, v) =
                        split_kv(rest).ok_or(format!("{}行目: inputs 内の書式エラー", lno))?;
                    cur.inputs.push((k, v));
                }
                Mode::Top => return Err(format!("{}行目: ブロック外の 4 スペース行", lno)),
            }
        } else if let Some(rest) = line.strip_prefix("  ") {
            let cur = claims
                .last_mut()
                .ok_or(format!("{}行目: エントリ外のフィールド", lno))?;
            match rest {
                "evidence:" => mode = Mode::Evidence,
                "inputs:" => mode = Mode::Inputs,
                "limitations:" => mode = Mode::Limitations,
                _ => {
                    mode = Mode::Top;
                    let (k, v) =
                        split_kv(rest).ok_or(format!("{}行目: フィールドの書式エラー", lno))?;
                    match k.as_str() {
                        "version" => cur.version = v,
                        "level" => cur.level = v,
                        "claim" => cur.claim = v,
                        "status" => cur.status = v,
                        _ => return Err(format!("{}行目: 未知のフィールド '{}'", lno, k)),
                    }
                }
            }
        } else {
            return Err(format!("{}行目: 解釈できない行: {}", lno, line));
        }
    }
    Ok(claims)
}

fn id_ok(id: &str) -> bool {
    // QRN-<英数大文字>-<3桁>
    let parts: Vec<&str> = id.split('-').collect();
    parts.len() == 3
        && parts[0] == "QRN"
        && !parts[1].is_empty()
        && parts[1]
            .chars()
            .all(|c| c.is_ascii_uppercase() || c.is_ascii_digit())
        && parts[2].len() == 3
        && parts[2].chars().all(|c| c.is_ascii_digit())
}

fn main() {
    println!("=== v6.1 主張台帳 (claims.yml) の機械検証 ===\n");
    let root = if Path::new("claims.yml").exists() {
        "."
    } else if Path::new("../claims.yml").exists() {
        ".."
    } else {
        println!("claims.yml が見つからない  [FAIL]");
        std::process::exit(1);
    };
    let text = fs::read_to_string(format!("{}/claims.yml", root)).expect("claims.yml 読み込み失敗");
    let mut nfail = 0;

    // [1] 解析
    let claims = match parse(&text) {
        Ok(c) => {
            println!("[1] 書式解析: {} 件の主張を解析  [PASS]", c.len());
            c
        }
        Err(e) => {
            println!("[1] 書式解析: {}  [FAIL]", e);
            std::process::exit(1);
        }
    };

    // [2] id の一意性と命名規則
    {
        let mut seen = HashSet::new();
        let mut bad = Vec::new();
        for c in &claims {
            if !id_ok(&c.id) {
                bad.push(format!("{} (命名規則違反)", c.id));
            }
            if !seen.insert(c.id.clone()) {
                bad.push(format!("{} (重複)", c.id));
            }
        }
        if bad.is_empty() {
            println!("[2] id 一意性・命名規則: 全 {} 件 OK  [PASS]", claims.len());
        } else {
            println!("[2] id 検査: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // [3] 必須フィールドと等級
    {
        const LEVELS: [&str; 6] = ["C0", "C1", "C2", "C3", "C4", "C5"];
        let mut bad = Vec::new();
        for c in &claims {
            if !LEVELS.contains(&c.level.as_str()) {
                bad.push(format!("{}: 等級 '{}'", c.id, c.level));
            }
            if c.version.is_empty() || c.claim.is_empty() || c.status.is_empty() {
                bad.push(format!("{}: 必須フィールド欠落", c.id));
            }
        }
        if bad.is_empty() {
            println!("[3] 等級・必須フィールド: OK  [PASS]");
        } else {
            println!("[3] 等級・必須フィールド: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // [4] 証拠ファイルの実在
    {
        let mut bad = Vec::new();
        for c in &claims {
            if c.evidence.is_empty() {
                bad.push(format!("{}: 証拠なし", c.id));
            }
            for (k, p) in &c.evidence {
                if !["code", "result", "doc"].contains(&k.as_str()) {
                    bad.push(format!("{}: 未知の証拠種別 '{}'", c.id, k));
                }
                if !Path::new(&format!("{}/{}", root, p)).exists() {
                    bad.push(format!("{}: {} が実在しない", c.id, p));
                }
            }
        }
        if bad.is_empty() {
            let ne: usize = claims.iter().map(|c| c.evidence.len()).sum();
            println!("[4] 証拠ファイルの実在: 全 {} 参照 OK  [PASS]", ne);
        } else {
            println!("[4] 証拠ファイル: {:?}  [FAIL]", bad);
            nfail += 1;
        }
    }

    // [5] 限界の明示 (C0 以外)
    {
        let bad: Vec<String> = claims
            .iter()
            .filter(|c| c.level != "C0" && c.limitations.is_empty())
            .map(|c| c.id.clone())
            .collect();
        if bad.is_empty() {
            println!("[5] 限界の明示 (C0 以外は必須): OK  [PASS]");
        } else {
            println!("[5] 限界の明示: {:?} に limitations がない  [FAIL]", bad);
            nfail += 1;
        }
    }

    // [6] 全バイナリの被覆 (台帳の取りこぼし検査)
    {
        let exempt = ["v61_ledger.rs"]; // 監査基盤そのものは対象外
        let mut bins: Vec<String> = fs::read_dir(format!("{}/sim/src/bin", root))
            .expect("sim/src/bin が読めない")
            .filter_map(|e| e.ok())
            .map(|e| e.file_name().to_string_lossy().to_string())
            .filter(|n| n.ends_with(".rs") && !exempt.contains(&n.as_str()))
            .collect();
        bins.sort();
        let covered: HashSet<String> = claims
            .iter()
            .flat_map(|c| c.evidence.iter())
            .filter(|(k, _)| k == "code")
            .filter_map(|(_, p)| {
                Path::new(p)
                    .file_name()
                    .map(|f| f.to_string_lossy().to_string())
            })
            .collect();
        let missing: Vec<&String> = bins.iter().filter(|b| !covered.contains(*b)).collect();
        if missing.is_empty() {
            println!(
                "[6] バイナリ被覆: 全 {} 本が台帳から参照されている  [PASS]",
                bins.len()
            );
        } else {
            println!("[6] バイナリ被覆: 未参照 {:?}  [FAIL]", missing);
            nfail += 1;
        }
    }

    // 集計
    println!("\n---- 等級別の集計 ----");
    let mut by_level: BTreeMap<String, Vec<&str>> = BTreeMap::new();
    for c in &claims {
        by_level.entry(c.level.clone()).or_default().push(&c.id);
    }
    const DESC: [(&str, &str); 6] = [
        ("C0", "教科書的・既存定理 (前提)"),
        ("C1", "既存結果の再現実装"),
        ("C2", "制限付き計算定理"),
        ("C3", "toy model による機構提示"),
        ("C4", "現象論的 fit / オーダー推定"),
        ("C5", "解釈・哲学・スローガン"),
    ];
    for (lv, desc) in DESC {
        let n = by_level.get(lv).map_or(0, |v| v.len());
        println!("  {} {:28} {:2} 件", lv, desc, n);
    }
    println!("  計 {} 件", claims.len());

    println!("\n---- 期別の主張数 ----");
    let mut by_major: BTreeMap<String, usize> = BTreeMap::new();
    for c in &claims {
        let major = c.version.split('.').next().unwrap_or("?").to_string();
        *by_major.entry(major).or_default() += 1;
    }
    for (m, n) in &by_major {
        println!("  {}.x : {} 件", m, n);
    }

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 台帳は形式的に健全 (内容の正しさは各文書と results/ が担う)"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
