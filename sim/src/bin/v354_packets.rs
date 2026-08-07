//! v35.4 FollowUp adversary packet A/B と OCS-2.0 の gating の常設監査 (PROMPT/16 §9)
//!
//! FollowUp は外部再現ではなく **paper-closed theorem adversary** として使う
//! (shared human operator は final report 自身が明示 — 組織的外部性はない)。
//! 論文更新は実装 commit ごとではなく atomic packet 単位:
//!   Packet A = Open Response (GQF 系 + 観測商 gauge の open 問題 OQ-0..OQ-8)
//!   Packet B = Resource Profile (RPF 系 + nesting 構成の open 問題 RP-1..RP-6)
//!   Packet C = OCS-2.0 — **A/B の反証と修正を受理した後にだけ凍結** (gating)。
//!
//! 検査:
//!   [K1] packet A/B の実在と sha256 凍結 pin
//!   [K2] 各 claim の必須フィールド (status/claim/falsifier/forbidden interpretation)
//!        と claim ID の全数
//!   [K3] paper-closed 規律: リポジトリ内部参照 (sim/src, results/, proofs/, .rs) と
//!        記録済み出力値の持ち込み禁止・expected_outputs_in_packet: false
//!   [K4] adversary 規律: RefutedAsStated / Inconclusive を一級の成果として受理する
//!        文言・独立導出の要求 (OQ-0 較正課題 — 規約の導出を渡さない)
//!   [K5] OCS-2.0 gating: paper/ に ocs-2.0 ファイルが**存在しない**こと (A/B の
//!        反証受理前の凍結は規律違反) + 文書側の gating 宣言
//!   [K6] scope discipline: no-go の契約 scope (curvature-only / 単一時点値) が
//!        packet 内で明示され、full-time 問題が open として攻撃面に出ていること

use std::fs;
use std::path::Path;
use uft_sim::{self_test, sha256_hex};

const PIN_PACKET_A: &str = "e23e82db5ca2128d07f0a115baf5cd6a4cabbcd1bc4f7d978cd0fb9eefd5abd4";
const PIN_PACKET_B: &str = "76d7267cb48b8e68b5c9a8b76b72f5e33ebbba982845a11a6eda2584d040dbea";

fn main() {
    self_test();
    println!("=== v35.4 adversary packet A/B と OCS-2.0 gating (PROMPT/16 §9) ===\n");
    let root = if Path::new("paper/packet-a-open-response.md").exists() {
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

    let pa = rd("paper/packet-a-open-response.md");
    let pb = rd("paper/packet-b-resource-profile.md");

    // ---------------- [K1] 実在と凍結 pin ----------------
    {
        let ha = sha256_hex(pa.as_bytes());
        let hb = sha256_hex(pb.as_bytes());
        check(
            "[K1] packet A/B の実在と sha256 凍結 pin",
            !pa.is_empty() && !pb.is_empty() && ha == PIN_PACKET_A && hb == PIN_PACKET_B,
            format!("A = {}…, B = {}…", &ha[..16], &hb[..16]),
        );
    }

    // ---------------- [K2] claim の必須フィールドと ID 全数 ----------------
    {
        let mut bad = Vec::new();
        let ids_a = [
            "OQ-0", "OQ-1", "OQ-2", "OQ-3", "OQ-4", "OQ-5", "OQ-6", "OQ-7", "OQ-8",
        ];
        for id in ids_a {
            if !pa.contains(id) {
                bad.push(format!("A に {} が無い", id));
            }
        }
        let ids_b = ["RP-1", "RP-2", "RP-3", "RP-4", "RP-5", "RP-6"];
        for id in ids_b {
            if !pb.contains(id) {
                bad.push(format!("B に {} が無い", id));
            }
        }
        // 必須フィールドの本数 (claim 節ごとに status/falsifier/forbidden)
        let count = |t: &str, needle: &str| t.matches(needle).count();
        if count(&pa, "- status:") < 8 {
            bad.push(format!("A の status: が {} 本 (< 8)", count(&pa, "- status:")));
        }
        if count(&pa, "- falsifier:") < 6 {
            bad.push(format!("A の falsifier: が {} 本 (< 6)", count(&pa, "- falsifier:")));
        }
        if count(&pa, "- forbidden interpretation:") < 7 {
            bad.push(format!(
                "A の forbidden interpretation: が {} 本 (< 7)",
                count(&pa, "- forbidden interpretation:")
            ));
        }
        if count(&pb, "- status:") < 6 || count(&pb, "- falsifier:") < 5 {
            bad.push("B の status/falsifier が不足".into());
        }
        check(
            "[K2] claim ID 全数 (A: OQ-0..8 / B: RP-1..6) と必須フィールド",
            bad.is_empty(),
            if bad.is_empty() {
                "A 9 claims / B 6 claims — status/falsifier/forbidden interpretation 完備".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [K3] paper-closed 規律 ----------------
    {
        let mut bad = Vec::new();
        for (name, t) in [("A", &pa), ("B", &pb)] {
            for forbidden in ["sim/src", "results/", "proofs/", ".rs", "cargo"] {
                if t.contains(forbidden) {
                    bad.push(format!("{} にリポジトリ内部参照 '{}'", name, forbidden));
                }
            }
            if !t.contains("expected_outputs_in_packet: false") {
                bad.push(format!("{} に expected_outputs_in_packet: false が無い", name));
            }
            if !t.contains("paper-closed") {
                bad.push(format!("{} に paper-closed 宣言が無い", name));
            }
        }
        check(
            "[K3] paper-closed: 内部参照なし・記録済み出力値なしの宣言",
            bad.is_empty(),
            if bad.is_empty() {
                "両 packet とも自己完結 (コード・出力値の持ち込みなし)".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [K4] adversary 規律 ----------------
    {
        let mut bad = Vec::new();
        for (name, t) in [("A", &pa), ("B", &pb)] {
            for needle in ["RefutedAsStated", "Inconclusive"] {
                if !t.contains(needle) {
                    bad.push(format!("{} に {} 受理宣言が無い", name, needle));
                }
            }
        }
        // OQ-0: 規約 (X, Y の導出・転置) を渡さず導出課題として出す
        if !pa.contains("deliberately left to the adversary") {
            bad.push("A の OQ-0 較正課題 (規約を渡さない) が無い".into());
        }
        check(
            "[K4] adversary 規律: RefutedAsStated/Inconclusive の一級受理・規約導出は渡さない",
            bad.is_empty(),
            if bad.is_empty() {
                "反証・不確定を成果として受理する宣言つき".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [K5] OCS-2.0 gating ----------------
    {
        let mut bad = Vec::new();
        // paper/ に ocs-2.0 系ファイルが存在しないこと
        if let Ok(dir) = fs::read_dir(format!("{}/paper", root)) {
            for e in dir.filter_map(|e| e.ok()) {
                let name = e.file_name().to_string_lossy().to_lowercase();
                if name.contains("ocs-2") || name.contains("ocs2") {
                    bad.push(format!("paper/{} が存在 (gating 違反)", name));
                }
            }
        }
        let doc = rd("docs/uft-v35.4.md");
        if !doc.contains("OCS-2.0") || !doc.contains("反証と修正を受理した後") {
            bad.push("uft-v35.4.md に OCS-2.0 gating 宣言が無い".into());
        }
        // OCS-1.0 は不変 (v342 が pin — ここでは存在のみ)
        if rd("paper/operational-core-spec.md").is_empty() {
            bad.push("OCS-1.0 が見つからない".into());
        }
        check(
            "[K5] OCS-2.0 gating: A/B の反証受理前に spec を凍結しない (ファイル不在 + 宣言)",
            bad.is_empty(),
            if bad.is_empty() {
                "OCS-2.0 は未凍結 (正しい状態) — OCS-1.0 は不変".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [K6] scope discipline ----------------
    {
        let mut bad = Vec::new();
        for needle in [
            "CurvatureOnlyOpenResponse",
            "full time series",
            "single-point value",
            "OQ-8",
        ] {
            if !pa.contains(needle) {
                bad.push(format!("A に scope 文言 '{}' が無い", needle));
            }
        }
        check(
            "[K6] scope discipline: no-go の契約 scope 明示 + full-time 問題が open 攻撃面",
            bad.is_empty(),
            if bad.is_empty() {
                "no-go は契約 scope つき — 誇張なし・攻撃面は open と宣言".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "全検査 PASS — packet A/B は凍結・OCS-2.0 は反証受理まで gating".to_string()
        } else {
            format!("FAIL {} 件", nfail)
        }
    );
}
