//! v25.2 凍結マニフェスト — 凍結成果物の SHA-256 台帳と不変性の常設監査 (第二十六期, IV/IV)
//!
//! v25 系列は v25.2 で停止する (PROMPT/7)。本バイナリは凍結を「宣言」ではなく
//! **機械検査可能な状態**にする:
//!   - 初走 (台帳なし): 凍結対象の SHA-256 と凍結値を results/v252_manifest.json に書く。
//!   - 以後の走行 (スイート含む): 台帳と現物を照合し、**凍結の破れを FAIL として検出**。
//!
//! 凍結対象 = 決定的なファイルのみ (文書・論文原稿・ソース・JSON 成果物・決定的 txt)。
//! タイミング行を含む結果 txt は対象外 (スイート再走の byte 決定性がない) — それらの
//! 数値内容は決定的な JSON と docs/uft-v25.2.md が固定する。
//!
//! 凍結値 (証明書 JSON と毎走行照合):
//!   λ_x = 1.185467287349258 / λ_⊥ = 1.229428764341310 (独立 2 実装 3.9e-16 一致)
//!   λ_x ∈ [1.185462071698428, 1.185472502995920] (区間証明)
//!   λ_⊥ ∈ [1.229385046338885, 1.229472484821826]
//!   異方性下界 λ_⊥ − λ_x ≥ 0.0439 (証明付き)
//!
//! 検査: [S1] 凍結対象の実在と SHA-256 / [S2] 凍結値 = 証明書 JSON /
//!       [S3] 台帳照合 (2 走目以降) または台帳作成 (初走)。

use uft_sim::*;

const FROZEN: &[&str] = &[
    "docs/uft-v25.2.md",
    "paper/modular-bw-full.md",
    "sim/src/iv.rs",
    "sim/src/bin/v251_gmu.rs",
    "sim/src/bin/v252_exact_g.rs",
    "sim/src/bin/v252_bz_certificate.rs",
    "sim/src/bin/v252_finite_size.rs",
    "sim/src/bin/v252_manifest.rs",
    "results/v251_gmu.json",
    "results/v252_exact_g.txt",
    "results/v252_bz_certificate.json",
];

// 証明書 JSON の完全精度値 (shortest-repr — 文書の {:.15} 表示より下位桁まで持つ)
const LAM_X: f64 = 1.185467287349258;
const LAM_P: f64 = 1.229_428_764_341_309_5;
const LAM_X_LO: f64 = 1.185462071698428;
const LAM_X_HI: f64 = 1.18547250299592;
const LAM_P_LO: f64 = 1.229_385_046_338_885_1;
const LAM_P_HI: f64 = 1.229_472_484_821_826_4;
const ANISO_LB: f64 = 0.0439;

fn read_repo(path: &str) -> Option<Vec<u8>> {
    std::fs::read(path)
        .or_else(|_| std::fs::read(format!("../{}", path)))
        .ok()
}

fn read_repo_str(path: &str) -> Option<String> {
    read_repo(path).and_then(|b| String::from_utf8(b).ok())
}

/// JSON 中の `"key": <num>` を拾う (自リポジトリの write_artifact 形式限定)
fn num_after(txt: &str, key: &str) -> Option<f64> {
    let pat = format!("\"{}\":", key);
    let p = txt.find(&pat)? + pat.len();
    let rest = &txt[p..];
    let end = rest.find(|c: char| c == ',' || c == '}' || c == '\n')?;
    rest[..end].trim().parse().ok()
}

/// 台帳から path の記録済みハッシュを拾う
fn recorded_hash(manifest: &str, path: &str) -> Option<String> {
    let p = manifest.find(&format!("\"{}\"", path))?;
    let rest = &manifest[p..];
    let q = rest.find("\"sha256\":")? + 9;
    let rest = &rest[q..];
    let a = rest.find('"')? + 1;
    let b = rest[a..].find('"')? + a;
    Some(rest[a..b].to_string())
}

/// 凍結時点の親コミット (best-effort — 台帳作成時のみ記録)
fn git_head() -> String {
    let head = read_repo_str(".git/HEAD").unwrap_or_default();
    let head = head.trim();
    if let Some(r) = head.strip_prefix("ref: ") {
        read_repo_str(&format!(".git/{}", r))
            .map(|s| s.trim().to_string())
            .unwrap_or_else(|| "unknown".into())
    } else if !head.is_empty() {
        head.to_string()
    } else {
        "unknown".into()
    }
}

fn main() {
    self_test();
    println!("=== v25.2 凍結マニフェスト — 不変性の常設監査 (第二十六期, IV/IV) ===\n");
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

    // ---- [S1] 凍結対象の実在とハッシュ ----
    let mut hashes: Vec<(String, String)> = Vec::new();
    let mut missing = Vec::new();
    for &p in FROZEN {
        match read_repo(p) {
            Some(b) => hashes.push((p.to_string(), sha256_hex(&b))),
            None => missing.push(p),
        }
    }
    check(
        "[S1] 凍結対象 11 ファイルの実在",
        missing.is_empty(),
        if missing.is_empty() {
            String::new()
        } else {
            format!("欠落: {:?}", missing)
        },
    );

    // ---- [S2] 凍結値 = 証明書 JSON ----
    {
        let cert = read_repo_str("results/v252_bz_certificate.json").unwrap_or_default();
        let eq = |k: &str, v: f64| num_after(&cert, k).map(|x| x == v).unwrap_or(false);
        let ok = eq("lambda_x_float", LAM_X)
            && eq("lambda_perp_float", LAM_P)
            && eq("lambda_x_lo", LAM_X_LO)
            && eq("lambda_x_hi", LAM_X_HI)
            && eq("lambda_perp_lo", LAM_P_LO)
            && eq("lambda_perp_hi", LAM_P_HI)
            && num_after(&cert, "anisotropy_lower_bound")
                .map(|x| x >= ANISO_LB)
                .unwrap_or(false);
        check(
            "[S2] 凍結値 = 証明書 (λ 15 桁・区間・異方性下界 ≥ 0.0439)",
            ok,
            format!("λ_x = {} / λ_⊥ = {}", LAM_X, LAM_P),
        );
    }

    // ---- [S3] 台帳の照合 (2 走目以降) または作成 (初走) ----
    match read_repo_str("results/v252_manifest.json") {
        Some(man) => {
            let mut drift = Vec::new();
            for (p, h) in &hashes {
                match recorded_hash(&man, p) {
                    Some(r) if &r == h => {}
                    Some(_) => drift.push(format!("{} (内容変更)", p)),
                    None => drift.push(format!("{} (台帳に無い)", p)),
                }
            }
            check(
                "[S3] 凍結不変性: 全 11 ファイルが台帳と一致",
                drift.is_empty(),
                if drift.is_empty() {
                    "凍結は保たれている".into()
                } else {
                    format!("**凍結の破れ**: {:?}", drift)
                },
            );
        }
        None => {
            let j = Json::Obj(vec![
                ("version".into(), Json::Str("v25.2".into())),
                ("kind".into(), Json::Str("freeze_manifest".into())),
                (
                    "policy".into(),
                    Json::Str(
                        "v25 系列は v25.2 で凍結 — v25.3 は開始しない。再開条件は FAL 台帳の発火または本台帳の照合破れのみ (PROMPT/7)".into(),
                    ),
                ),
                ("frozen_at_parent_commit".into(), Json::Str(git_head())),
                (
                    "novelty_boundary".into(),
                    Json::Str(
                        "1D closed form g(mu) = Eisler (J.Stat.Mech.(2025)013101) = prior art; this repo contributes: exact 3D->1D reduction + BZ moment formulas + interval certificates + anisotropy theorem + machine-validated claim boundary".into(),
                    ),
                ),
                (
                    "frozen_values".into(),
                    Json::Obj(vec![
                        ("lambda_x".into(), Json::Num(LAM_X)),
                        ("lambda_perp".into(), Json::Num(LAM_P)),
                        ("lambda_x_lo".into(), Json::Num(LAM_X_LO)),
                        ("lambda_x_hi".into(), Json::Num(LAM_X_HI)),
                        ("lambda_perp_lo".into(), Json::Num(LAM_P_LO)),
                        ("lambda_perp_hi".into(), Json::Num(LAM_P_HI)),
                        ("anisotropy_lower_bound_min".into(), Json::Num(ANISO_LB)),
                    ]),
                ),
                (
                    "files".into(),
                    Json::Arr(
                        hashes
                            .iter()
                            .map(|(p, h)| {
                                Json::Obj(vec![
                                    ("path".into(), Json::Str(p.clone())),
                                    ("sha256".into(), Json::Str(h.clone())),
                                ])
                            })
                            .collect(),
                    ),
                ),
            ]);
            let p = write_artifact("results/v252_manifest.json", &j.render());
            check(
                "[S3] 凍結台帳を作成 (初走 — 以後の走行は照合モード)",
                true,
                format!("{} ({} ファイル)", p, hashes.len()),
            );
        }
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "凍結は機械検査可能な状態にある — 破れは本バイナリの FAIL として現れる"
        } else {
            "**凍結の破れ、または凍結成果物の欠落** — 台帳と現物を照合せよ"
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
