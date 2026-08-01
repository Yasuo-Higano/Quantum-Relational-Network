//! v31.0 幾何読み出しの識別可能性契約 — 意味論・型・protocol 凍結の機械検査 (PROMPT/12)
//!
//! 第三十一期の最初の道具。「何が状態に符号化され (E0)、完全な大域状態から何が
//! 逆転でき (E1)、制限された観測から何が読め (E2)、因子分解自体は選べるか (E3)、
//! 自然への橋か (E4)」を別能力として型で分離した契約 (sim/src/readout_contract.rs)
//! と、その文書・schema・台帳への波及を検査する:
//!   [R0] readout_contract_self_test — タグ 13 種の一意性・exact lane の床・
//!        生成子 2 型の門の棄却挙動・証明書コンストラクタの封鎖
//!   [R1] 型レベル封鎖の source 検査 — 禁止 impl From 3 種 (親→物理生成子 /
//!        正則化→exact / oracle→operational patch) の不在 + verdict 私有
//!   [R2] core.schema.yml — 新概念 15 種の登録・RelationalDecomposition の意味論差
//!        是正 (「入力でなく読み出し」→ 設計目標と入力の分離)・禁止変換 8–10
//!   [R3] replications.yml — claim/capability scoped 拡張 (Unit D schema)・
//!        external_replications = 0 維持・v27.4 の 6 条件不変
//!   [R4] 文書アンカー — docs/uft-v31.0.md の E0–E4 分離・絶対禁止・採点原則
//!   [R5] qrn_core 追加のみ検査 — v30.0 の封鎖 (ProperTime 門の不在・ClockCalibration
//!        構成不能・登録簿空) が第三十一期の変更後も破れていない
//!   [R6] 常設監査への自己登録 — tools/suite.sh の ALWAYS_RUN に本監査が載っている
//!   [R7] 熱的 Gaussian round-trip 実演 — P4 鎖 C = (I+e^{β(h−μI)})⁻¹ から
//!        K(C) = log[(I−C)C⁻¹] = β(h−μI) を経て門で h を復元 (v31.1 の先行最小例)
//!   [R9] 外部再現プロトコルの版分離 (v32.1, PROMPT/13 で追加) — 凍結プロトコル
//!        一式の sha256-16 認証・v27.4 版付き複製 = 原本の byte 一致・
//!        supersession 台帳 (ERR-D2-V1) の実在。凍結ファイルの変更は本監査の
//!        認証値の意識的更新 (= 版分離) を要する
//!
//! 本監査はスイートの常時実行層 (ALWAYS_RUN)。契約の分類自体の物理的正しさは
//! 保証しない (ASM-LAYER-SEMANTICS と同種の規約) — 保証するのは封鎖経路の不在と
//! 文書・schema・型の一致である。

use std::fs;
use std::path::Path;
use uft_sim::readout_contract::*;
use uft_sim::{jacobi_eigh, matfun_sym, sha256_hex};

fn main() {
    uft_sim::self_test();
    println!("=== v31.0 幾何読み出しの識別可能性契約 — 意味論・型・protocol 凍結 (PROMPT/12) ===\n");
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

    // ---- [R0] 契約の自己検査 ----
    {
        let r = readout_contract_self_test();
        check(
            "[R0] readout_contract_self_test — タグ 13 種・exact 床・門の棄却・証明書封鎖",
            r.is_ok(),
            match &r {
                Ok(()) => format!(
                    "状態領域 4 / 観測契約 6 / 因子分解 3 タグ・棄却理由 {} 種・δ 床 {:.0e}",
                    AbstainReason::ALL.len(),
                    DELTA_EXACT_FLOOR
                ),
                Err(e) => e.clone(),
            },
        );
    }

    // ---- [R1] 型レベル封鎖の source 検査 ----
    {
        const FORBIDDEN_IMPLS: [&str; 3] = [
            "impl From<ParentModularGenerator> for PhysicalGenerator",
            "impl From<RegularizedCorrelation> for ExactFullRankCorrelation",
            "impl From<OraclePatch> for OperationalPatch",
        ];
        const EXEMPT: [&str; 1] = ["v310_readout_semantics.rs"];
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
                        for f in FORBIDDEN_IMPLS {
                            if text.contains(f) {
                                hits.push(format!("{}: {}", name, f));
                            }
                        }
                    }
                }
            }
        }
        // verdict 私有 + ExactWitness に Regularized variant がないことの source 確認
        let contract = rd("sim/src/readout_contract.rs").unwrap_or_default();
        if !contract.contains("    verdict: IdentifiabilityVerdict,") {
            hits.push("readout_contract.rs: verdict 私有フィールドが見つからない".into());
        }
        if contract.contains("pub verdict") {
            hits.push("readout_contract.rs: verdict が pub になっている".into());
        }
        if contract.contains("FullRankCorrelation(&'a RegularizedCorrelation)") {
            hits.push("readout_contract.rs: ExactWitness に正則化相関が混入".into());
        }
        check(
            "[R1] 型レベル封鎖 — 禁止 impl From 3 種の不在 (sim/src 全走査) + verdict 私有",
            hits.is_empty(),
            if hits.is_empty() {
                "親→物理生成子 / 正則化→exact / oracle→operational patch の経路なし".into()
            } else {
                format!("{:?}", hits)
            },
        );
    }

    // ---- [R2] core.schema.yml の登録 ----
    {
        let schema = rd("core.schema.yml").unwrap_or_default();
        let mut bad = Vec::new();
        const NEW_CONCEPTS: [&str; 15] = [
            "ReadoutCertificate",
            "StateDomainTag",
            "ObservationContractTag",
            "FactorizationStatusTag",
            "IdentifiabilityVerdict",
            "ParentModularGenerator",
            "PhysicalGenerator",
            "GaussianityEvidence",
            "GibbsProvenance",
            "GivenNodeFactorization",
            "RelationalDecompositionGoal",
            "OraclePatch",
            "OperationalPatch",
            "ExactFullRankCorrelation",
            "RegularizedCorrelation",
        ];
        for c in NEW_CONCEPTS {
            if !schema.contains(&format!("- name: {}\n", c)) {
                bad.push(format!("概念 {} が未登録", c));
            }
        }
        // RelationalDecomposition の意味論差是正: 旧文言の根絶と新分離の明示
        if schema.contains("テンソル分解は入力でなく読み出し") {
            bad.push("RelationalDecomposition の旧文言「入力でなく読み出し」が残存".into());
        }
        for needle in ["RelationalDecompositionGoal, 未構成", "GivenNodeFactorization を入力に取る"] {
            if !schema.contains(needle) {
                bad.push(format!("意味論差是正の文言「{}」が無い", needle));
            }
        }
        // 禁止変換 8–10
        for (f, t) in [
            ("ParentModularGenerator", "PhysicalGenerator"),
            ("RegularizedCorrelation", "ExactFullRankCorrelation"),
            ("OraclePatch", "OperationalPatch"),
        ] {
            if !schema.contains(&format!("- from: {}\n  to: {}\n  reason:", f, t)) {
                bad.push(format!("禁止変換 {} → {} が未登録", f, t));
            }
        }
        check(
            "[R2] core.schema.yml — 新概念 15 種・意味論差是正・禁止変換 8–10 (reason つき)",
            bad.is_empty(),
            if bad.is_empty() {
                "「定義済み読み出し」と「設計目標」と「現行 bridge の入力」が分離された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [R3] replications.yml の capability scoped 拡張 ----
    {
        let rep = rd("replications.yml").unwrap_or_default();
        let mut bad = Vec::new();
        for field in [
            "capability_scope_version:",
            "claim_ids:",
            "capabilities:",
            "replication_level:",
            "independence_scope:",
            "protocol_commit:",
            "generator_hash:",
            "input_hash:",
        ] {
            if !rep.contains(field) {
                bad.push(format!("フィールド {} が無い", field));
            }
        }
        for needle in ["D1", "D2", "D3", "D4", "unit A/B/C は geometry 能力の blocker を解除しない"] {
            if !rep.contains(needle) {
                bad.push(format!("Unit D schema の「{}」が無い", needle));
            }
        }
        // v27.4 の既存契約の不変
        for c in [
            "different_author",
            "independent_repository",
            "no_shared_numerical_kernel",
            "protocol_frozen_before_run",
            "commit_hash_recorded",
            "result_including_failures_public",
        ] {
            if !rep.contains(c) {
                bad.push(format!("v27.4 の条件 {} が壊れた", c));
            }
        }
        let declared: usize = rep
            .lines()
            .find_map(|l| l.strip_prefix("external_replications:"))
            .and_then(|v| v.trim().parse().ok())
            .unwrap_or(999);
        if declared != 0 {
            bad.push(format!("external_replications = {} (0 のはず)", declared));
        }
        check(
            "[R3] replications.yml — claim/capability scoped 拡張・Unit D schema・external 0 維持",
            bad.is_empty(),
            if bad.is_empty() {
                "gauge 単位の成功は geometry 能力を解除しない (D2 のみ解除可)".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [R4] 文書アンカー ----
    {
        let doc = rd("docs/uft-v31.0.md").unwrap_or_default();
        let mut bad = Vec::new();
        for needle in [
            "E0",
            "E4",
            "ParentModularGenerator",
            "絶対禁止",
            "非識別セル",
            "識別可能性",
            "GaussianGibbsInverseOracle",
            "oracle ceiling",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v31.0.md: 「{}」が無い", needle));
            }
        }
        let core_src = rd("sim/src/qrn_core.rs").unwrap_or_default();
        if !core_src.contains("GivenNodeFactorization を入力に取る") {
            bad.push("qrn_core.rs の RelationalDecomposition note が未是正".into());
        }
        check(
            "[R4] 文書アンカー — uft-v31.0.md の E0–E4 分離・絶対禁止・採点原則・qrn_core note",
            bad.is_empty(),
            if bad.is_empty() {
                "意味論凍結が文書に固定されている".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [R5] qrn_core の封鎖の不変 (v30.0 → v31.0 で破れていない) ----
    {
        let core_src = rd("sim/src/qrn_core.rs").unwrap_or_default();
        let mut bad = Vec::new();
        if core_src.contains("fn promote_evolution_to_proper_time(") {
            bad.push("ProperTime への門が復活している".to_string());
        }
        if core_src.contains("impl_capability!(ClockCalibration") {
            bad.push("ClockCalibration に能力が実装されている".to_string());
        }
        if !core_src.contains("impl sealed_cap::Sealed for ClockCalibration {}") {
            bad.push("ClockCalibration の sealed が消えた".to_string());
        }
        if let Err(e) = uft_sim::qrn_core::qrn_core_self_test() {
            bad.push(format!("qrn_core_self_test: {}", e));
        }
        check(
            "[R5] qrn_core の封鎖不変 — ProperTime 門の不在・ClockCalibration 構成不能・登録簿空",
            bad.is_empty(),
            if bad.is_empty() {
                "v30.0 の封鎖は v31.0 の追加後も保たれている".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [R6] 常設監査への自己登録 ----
    {
        let suite = rd("tools/suite.sh").unwrap_or_default();
        let registered = suite
            .lines()
            .find(|l| l.starts_with("ALWAYS_RUN="))
            .map(|l| l.contains("v310_readout_semantics"))
            .unwrap_or(false);
        check(
            "[R6] tools/suite.sh の ALWAYS_RUN に本監査が登録されている",
            registered,
            if registered {
                "リポジトリ状態を読む監査層として常時実行".into()
            } else {
                "ALWAYS_RUN に v310_readout_semantics が無い".into()
            },
        );
    }

    // ---- [R7] 熱的 Gaussian round-trip 実演 (v31.1 の先行最小例) ----
    {
        // P4 鎖 (4 ノード) の一体生成子 h = −t A + 対角, β = 1.7, μ = 0.3
        let n = 4;
        let mut h = vec![0.0; n * n];
        for i in 0..n - 1 {
            h[i * n + i + 1] = -1.0;
            h[(i + 1) * n + i] = -1.0;
        }
        h[0] = 0.2; // 対角の不均一 (縮退回避)
        h[n * n - 1] = -0.1;
        let beta = 1.7;
        let mu = 0.3;
        // C = (I + e^{β(h−μI)})⁻¹ — matfun_sym で厳密構成
        let mut hm = h.clone();
        for i in 0..n {
            hm[i * n + i] -= mu;
        }
        let c = matfun_sym(&hm, n, |x| 1.0 / (1.0 + (beta * x).exp()));
        // 資格審査 (clamp なし) → K(C) → 門 → h 復元 → projector 対照
        let (ok, detail) = match ExactFullRankCorrelation::certify_real(&c, n) {
            Err(e) => (false, format!("資格審査に失敗: {:?}", e)),
            Ok(cert) => {
                // K(C) = log[(I−C)C⁻¹] = β(h−μI)
                let k = matfun_sym(cert.c_re(), n, |x| ((1.0 - x) / x).ln());
                let parent = ParentModularGenerator {
                    re: k,
                    im: vec![0.0; n * n],
                    n,
                };
                // 門: β, μ 既知 → h を厳密復元
                let (mut ok, mut detail) = match identify_physical_generator(
                    &parent,
                    GaussianityEvidence::ByConstruction,
                    GibbsProvenance::KnownBetaMu { beta, mu },
                ) {
                    Ok(PhysicalGeneratorReading::Exact(hr)) => {
                        let mut err_max: f64 = 0.0;
                        for (a, b) in hr.re.iter().zip(h.iter()) {
                            err_max = err_max.max((a - b).abs());
                        }
                        // 縮退のない固有系での logit round-trip: f64 で ~1e-12 級
                        (
                            err_max < 1e-10,
                            format!(
                                "δ = {:.3e} / 復元誤差 max|ĥ−h| = {:.2e} (バー 1e-10)",
                                cert.spectral_margin(),
                                err_max
                            ),
                        )
                    }
                    other => (false, format!("門が開かない: {:?}", other.err())),
                };
                // 対照: projector (β→∞ 極限) は資格を通らない
                let (evals, evecs) = jacobi_eigh(&hm, n);
                let mut proj = vec![0.0; n * n];
                for m in 0..n {
                    if evals[m] < 0.0 {
                        for i in 0..n {
                            for j in 0..n {
                                proj[i * n + j] += evecs[m * n + i] * evecs[m * n + j];
                            }
                        }
                    }
                }
                match ExactFullRankCorrelation::certify_real(&proj, n) {
                    Err(AbstainReason::RankDeficient) => {}
                    _ => {
                        ok = false;
                        detail.push_str(" / projector が RankDeficient にならない");
                    }
                }
                (ok, detail)
            }
        };
        check(
            "[R7] 熱的 Gaussian round-trip — C = (I+e^{β(h−μI)})⁻¹ → K(C) → 門 → h 復元・projector は棄却",
            ok,
            detail,
        );
    }

    // ---- [R8] 外部再現 Unit D の実在と整合 (Track X — v31 期で追加) ----
    {
        let mut bad = Vec::new();
        let unit_d = rd("reproducer/UNIT_D.md").unwrap_or_default();
        if unit_d.is_empty() {
            bad.push("reproducer/UNIT_D.md が無い".to_string());
        }
        for needle in [
            "D1",
            "D2",
            "D3",
            "D4",
            "geometry 能力の blocker を解除できるのはこの水準のみ",
            "同一 AI による再実装は independence を満たさない",
            "β₃ には ∂₄ が必須",
            "L ≥ 4",
        ] {
            if !unit_d.contains(needle) {
                bad.push(format!("UNIT_D.md: 「{}」が無い", needle));
            }
        }
        let d1 = rd("reproducer/INPUTS/unit_d1_frozen_c.json").unwrap_or_default();
        for needle in ["u693", "C_reference_1e-9", "\"beta\": 1.0"] {
            if !d1.contains(needle) {
                bad.push(format!("unit_d1_frozen_c.json: 「{}」が無い", needle));
            }
        }
        let claims_md = rd("reproducer/CLAIMS.md").unwrap_or_default();
        if !claims_md.contains("単位 D") || !claims_md.contains("QRN-BRIDGE-013") {
            bad.push("CLAIMS.md に単位 D 表が無い".to_string());
        }
        check(
            "[R8] 外部再現 Unit D (Track X): UNIT_D.md (D1–D4)・凍結入力 D1・CLAIMS.md の単位 D 表 — geometry 解除は D2 のみ",
            bad.is_empty(),
            if bad.is_empty() {
                "geometry 用 reproducer が公開された — external_replications = 0 は維持".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [R9] 外部再現プロトコルの版分離 (v32.1, PROMPT/13) — 凍結一式の常設監査 ----
    {
        let mut bad = Vec::new();
        // 凍結プロトコル一式の sha256-16 (変更は版分離 + 認証値の意識的更新のみ)
        let frozen: [(&str, &str); 6] = [
            (
                "reproducer/protocols/v32.1/unit-d-report.schema.json",
                "f578816c54db3d23",
            ),
            (
                "reproducer/protocols/v32.1/unit-d-tolerances.yml",
                "0ebc7098c6961355",
            ),
            ("reproducer/protocols/v32.1/d2-static.md", "f858fa4bdeaa3554"),
            ("reproducer/protocols/v32.1/d2-response.md", "ce27f04b56303110"),
            (
                "reproducer/protocols/v32.1/protocol-index.yml",
                "62d2ef94632ec5ac",
            ),
            (
                "reproducer/protocols/v31.7/d2-v1-superseded.md",
                "a4419dc8794558a6",
            ),
        ];
        for (f, want) in frozen {
            match rd(f) {
                Err(_) => bad.push(format!("{} が無い", f)),
                Ok(t) => {
                    let h = sha256_hex(t.as_bytes());
                    if &h[..16] != want {
                        bad.push(format!("{} の sha256-16 {} ≠ 凍結値 {}", f, &h[..16], want));
                    }
                }
            }
        }
        // v27.4 の版付き複製 = 凍結原本と byte 一致 (drift の禁止)
        for (copy, orig) in [
            (
                "reproducer/protocols/v27.4/abc-report.schema.json",
                "reproducer/EXPECTED_SCHEMA.json",
            ),
            (
                "reproducer/protocols/v27.4/abc-tolerances.yml",
                "reproducer/TOLERANCES.yml",
            ),
        ] {
            match (rd(copy), rd(orig)) {
                (Ok(c), Ok(o)) => {
                    if sha256_hex(c.as_bytes()) != sha256_hex(o.as_bytes()) {
                        bad.push(format!("{} が原本 {} と byte 不一致", copy, orig));
                    }
                }
                _ => bad.push(format!("{} または {} が読めない", copy, orig)),
            }
        }
        // supersession 台帳と UNIT_D.md の版分離ポインタ
        let rep = rd("replications.yml").unwrap_or_default();
        for needle in ["ERR-D2-V1", "superseded_before_external_run"] {
            if !rep.contains(needle) {
                bad.push(format!("replications.yml: 「{}」が無い", needle));
            }
        }
        let unit_d = rd("reproducer/UNIT_D.md").unwrap_or_default();
        for needle in ["D2-S", "D2-R", "protocols/v32.1/"] {
            if !unit_d.contains(needle) {
                bad.push(format!("UNIT_D.md: 「{}」が無い", needle));
            }
        }
        check(
            "[R9] 外部再現プロトコルの版分離 (v32.1): 凍結一式 sha256 認証・v27.4 複製 byte 一致・supersession 台帳",
            bad.is_empty(),
            if bad.is_empty() {
                "D2-v1 の supersession と D2-S/D2-R の凍結が常設監査下にある".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "意味論・型・protocol は凍結された — encoding / inversion / operational readout / factorization / bridge は別能力として型で分離されている"
        } else {
            "**契約の破れ** — readout_contract と schema/台帳/文書の整合を修正せよ"
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
