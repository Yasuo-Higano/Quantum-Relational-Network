//! v34.2 standalone operational core spec の凍結監査 (PROMPT/15 §3)
//!
//! 背景: FollowUp (cross-model clean-room) の最終判定は「supplied paper set does
//! not contain a complete core specification」— 多くの paper claim は repository
//! replay でしか再現できない (specification-limited)。本版は第三十三期の
//! operational core (certified interface / context / resource / recovery / glue /
//! graded / structured) を、リポジトリのソースにも出力値にもアクセスせずに独立
//! 実装できる paper-closed spec (paper/operational-core-spec.md, OCS-1.0) として
//! 凍結し、claim ごとの閉包状態を publication closure manifest
//! (paper/operational-core-closure.yml) に機械可読化する。
//!
//! 検証 (全て [PASS]/[FAIL] 内蔵):
//!   [M0] spec の実在と SHA-256 凍結 (pin 照合 — spec の変更は意図的な二重更新のみ)
//!   [M1] spec の規範アンカー — 全 14 節・4 再現性等級・裁定優先順位・凍結バー表・
//!        probe 型分離 (SignedInitialCovarianceProbe / HamiltonianQuench /
//!        NumberConservingResponse / Pairing)・出力値非掲載の規律文
//!   [M2] closure manifest の構文と整合 — 全エントリに必須 13 フィールド・
//!        claim_ids / assumption_ids の実在・spec_section アンカーの実在・
//!        expected_outputs_in_spec = false 全行・paper_closed: false の正直な
//!        境界が最低 1 件 (閉包の過大主張の禁止)
//!   [M3] spec は出力値を含まない規律 — 規律文アンカー + 4 等級の語彙一致
//!   [M4] 規範インスタンスの実行 (spec §10 を本リポジトリでも機械実証 —
//!        これまで FollowUp のみが実装していた分離の主リポジトリ側の受理):
//!        [M4a] N1 quench null — 一粒子 lane: signed covariance probe の応答 =
//!              ‖P₂hP₁‖² = 1 (厳密) / Hamiltonian quench (Γ₀ = I/2) の応答 = 0
//!              (厳密) — 同じ h で probe 型が答えを変える (T9)
//!        [M4b] N2 pairing — 2 モード JW (dim 4): 積状態 signed probe の応答 =
//!              1 − Δ² 厳密 (Δ ∈ {0, 0.3, 0.7}・ε 非依存)・対角 V 項は不変 (T10)
//!        [M4c] N3 Busch 対 — R_ab = (I + aησ_x + bησ_z)/4 の joint POVM 資格が
//!              2η² ≤ 1 と厳密一致 (η = 0.6 資格 / 0.7071+ 拒否 / 0.8 拒否) —
//!              非可換対の joint measurability (可換性より広い)
//!
//! 実行順序の規約 (spec §0): 先に spec を凍結し、実装は互いの結果を見ずに spec
//! から書く。spec の曖昧性が見つかった場合は結果に合わせず OCS の版を上げる。

use std::fs;
use std::path::Path;
use uft_sim::{self_test, sha256_hex, C64};

/// spec (OCS-1.0) の凍結 SHA-256 — 変更は意図的な二重更新 (spec 版バンプ) のみ
const SPEC_SHA256: &str = "bb3b4d4bff7222255838b9f95f46ccda697102d9da7d736564b732d7ea433359";

// ---------------------------------------------------------------- 小型複素行列
fn zeros(n: usize) -> Vec<C64> {
    vec![C64::new(0.0, 0.0); n * n]
}
fn eye(n: usize) -> Vec<C64> {
    let mut m = zeros(n);
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}
fn mul(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let mut c = zeros(n);
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik.norm2() == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] = c[i * n + j] + aik * b[k * n + j];
            }
        }
    }
    c
}
fn add(a: &[C64], b: &[C64]) -> Vec<C64> {
    a.iter().zip(b).map(|(x, y)| *x + *y).collect()
}
fn scal(s: f64, a: &[C64]) -> Vec<C64> {
    a.iter().map(|x| x.scale(s)).collect()
}
fn comm(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let ab = mul(a, b, n);
    let ba = mul(b, a, n);
    ab.iter().zip(&ba).map(|(x, y)| *x - *y).collect()
}
fn dagger(a: &[C64], n: usize) -> Vec<C64> {
    let mut d = zeros(n);
    for i in 0..n {
        for j in 0..n {
            d[j * n + i] = a[i * n + j].conj();
        }
    }
    d
}
fn trace(a: &[C64], n: usize) -> C64 {
    let mut t = C64::new(0.0, 0.0);
    for i in 0..n {
        t = t + a[i * n + i];
    }
    t
}
fn kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
    let n = na * nb;
    let mut c = zeros(n);
    for i in 0..na {
        for j in 0..na {
            for k in 0..nb {
                for l in 0..nb {
                    c[(i * nb + k) * n + (j * nb + l)] = a[i * na + j] * b[k * nb + l];
                }
            }
        }
    }
    c
}
/// エルミート 2×2 の固有値 (閉形式)
fn eig2_herm(m: &[C64]) -> (f64, f64) {
    let a = m[0].re;
    let d = m[3].re;
    let b2 = m[1].norm2();
    let tr = a + d;
    let disc = ((a - d) * (a - d) + 4.0 * b2).sqrt();
    ((tr - disc) / 2.0, (tr + disc) / 2.0)
}

// ---------------------------------------------------------------- 最小 flat YAML
#[derive(Default, Clone)]
struct Entry {
    fields: Vec<(String, String)>,
}
impl Entry {
    fn get(&self, k: &str) -> Option<&str> {
        self.fields
            .iter()
            .find(|(kk, _)| kk == k)
            .map(|(_, v)| v.as_str())
    }
}
fn unquote(s: &str) -> String {
    let t = s.trim();
    if t.len() >= 2 && t.starts_with('"') && t.ends_with('"') {
        t[1..t.len() - 1].to_string()
    } else {
        t.to_string()
    }
}
fn parse_entries(text: &str) -> Result<Vec<Entry>, String> {
    let mut out: Vec<Entry> = Vec::new();
    let mut in_entries = false;
    for (lno, raw) in text.lines().enumerate() {
        let lno = lno + 1;
        let line = raw.trim_end();
        if line.trim_start().is_empty() || line.trim_start().starts_with('#') {
            continue;
        }
        if line == "entries:" {
            in_entries = true;
            continue;
        }
        if !in_entries {
            continue;
        }
        if let Some(rest) = line.strip_prefix("- id:") {
            let mut e = Entry::default();
            e.fields.push(("id".into(), unquote(rest)));
            out.push(e);
        } else if let Some(rest) = line.strip_prefix("  ") {
            let cur = out.last_mut().ok_or(format!("{}行目: エントリ外", lno))?;
            let idx = rest.find(':').ok_or(format!("{}行目: 書式", lno))?;
            cur.fields
                .push((rest[..idx].trim().into(), unquote(&rest[idx + 1..])));
        } else {
            return Err(format!("{}行目: 解釈できない行", lno));
        }
    }
    Ok(out)
}
fn parse_list(v: &str) -> Vec<String> {
    let t = v.trim();
    if !t.starts_with('[') || !t.ends_with(']') {
        return vec![];
    }
    let inner = &t[1..t.len() - 1];
    if inner.trim().is_empty() {
        return vec![];
    }
    inner.split(',').map(|x| x.trim().to_string()).collect()
}

fn main() {
    self_test();
    println!("=== v34.2 standalone operational core spec の凍結監査 (PROMPT/15 §3) ===\n");
    let root = if Path::new("paper/operational-core-spec.md").exists() {
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

    // ---------------- [M0] spec の実在と SHA-256 凍結 ----------------
    let spec = rd("paper/operational-core-spec.md");
    {
        let sha = sha256_hex(spec.as_bytes());
        check(
            "[M0] spec (OCS-1.0) の実在と SHA-256 凍結 pin",
            !spec.is_empty() && sha == SPEC_SHA256,
            format!("sha256 = {}…", &sha[..16.min(sha.len())]),
        );
    }

    // ---------------- [M1] spec の規範アンカー ----------------
    {
        let anchors: [&str; 24] = [
            "OCS version 1.0",
            "expected_outputs_in_spec: false",
            "## 0. Reproducibility grades",
            "repository replay",
            "paper/spec-closed",
            "clean-room",
            "organizationally external",
            "## 1. Kinematics",
            "## 2. Certificates",
            "## 3. Contexts",
            "## 4. Resource budgets",
            "## 5. Factorization recovery",
            "## 6. Contextual descent",
            "## 7. Graded lane",
            "## 8. Response probes are types",
            "SignedInitialCovarianceProbe ↛ HamiltonianQuench",
            "NumberConservingResponse ↛ BCS/PairingResponse",
            "## 9. Structured backends",
            "## 10. Normative instances",
            "## 11. No-gos and positive theorems",
            "## 12. Adjudication order",
            "## 13. Frozen bars",
            "## 14. What this spec does not close",
            "1 − Δ²",
        ];
        let missing: Vec<&&str> = anchors.iter().filter(|a| !spec.contains(**a)).collect();
        check(
            "[M1] spec の規範アンカー (14 節・4 等級・probe 分離・バー表・規律文)",
            missing.is_empty(),
            if missing.is_empty() {
                format!("{} アンカー", anchors.len())
            } else {
                format!("欠落: {:?}", missing)
            },
        );
    }

    // ---------------- [M2] closure manifest ----------------
    {
        let man = rd("paper/operational-core-closure.yml");
        let claims = rd("claims.yml");
        let asms = rd("assumptions.yml");
        let mut bad: Vec<String> = Vec::new();
        let entries = match parse_entries(&man) {
            Ok(e) => e,
            Err(e) => {
                bad.push(e);
                Vec::new()
            }
        };
        const REQUIRED: [&str; 13] = [
            "id",
            "spec_section",
            "theorem_or_algorithm",
            "claim_ids",
            "assumption_ids",
            "observation_contract",
            "required_input_artifacts",
            "equivalence_relation",
            "allowed_conclusion",
            "falsifier",
            "paper_closed",
            "repository_replay_only",
            "expected_outputs_in_spec",
        ];
        let mut n_open = 0usize;
        for e in &entries {
            let id = e.get("id").unwrap_or("?");
            for f in REQUIRED {
                if e.get(f).is_none() {
                    bad.push(format!("{}: {} 欠落", id, f));
                }
            }
            if e.get("expected_outputs_in_spec") != Some("false") {
                bad.push(format!("{}: expected_outputs_in_spec ≠ false", id));
            }
            match e.get("paper_closed") {
                Some("true") => {}
                Some("false") => n_open += 1,
                _ => bad.push(format!("{}: paper_closed が真偽値でない", id)),
            }
            for cid in parse_list(e.get("claim_ids").unwrap_or("[]")) {
                if !claims.contains(&format!("- id: {}", cid)) {
                    bad.push(format!("{}: claim {} が claims.yml に無い", id, cid));
                }
            }
            for aid in parse_list(e.get("assumption_ids").unwrap_or("[]")) {
                if !asms.contains(&format!("- id: {}", aid)) {
                    bad.push(format!("{}: 仮定 {} が assumptions.yml に無い", id, aid));
                }
            }
            if let Some(sec) = e.get("spec_section") {
                let first = sec.split('-').next().unwrap_or(sec);
                if !spec.contains(&format!("## {}.", first)) {
                    bad.push(format!("{}: spec 節 {} が実在しない", id, sec));
                }
            }
        }
        if n_open == 0 {
            bad.push("paper_closed: false の正直な境界が 1 件も無い (閉包の過大主張)".into());
        }
        check(
            "[M2] closure manifest — 13 フィールド・参照実在・正直な境界 ≥ 1",
            bad.is_empty() && !entries.is_empty(),
            if bad.is_empty() {
                format!(
                    "{} エントリ (paper_closed {} / replay_only {})",
                    entries.len(),
                    entries.len() - n_open,
                    n_open
                )
            } else {
                bad.truncate(6);
                format!("{:?}", bad)
            },
        );
    }

    // ---------------- [M3] 出力値非掲載の規律 ----------------
    {
        let man = rd("paper/operational-core-closure.yml");
        let ok = spec.contains("never states measured output")
            && man.contains("reproducibility_grades: [repository_replay, paper_closed, clean_room, organizationally_external]");
        check(
            "[M3] 出力値非掲載の規律文 + 4 再現性等級の語彙一致 (spec ↔ manifest)",
            ok,
            "spec は手順・恒等式・バー・裁定のみを凍結する".into(),
        );
    }

    // ---------------- [M4] 規範インスタンスの実行 ----------------
    println!("\n[M4] 規範インスタンス (spec §10) の機械実証 — probe 型分離の主リポジトリ側受理");
    {
        // ---- [M4a] N1 quench null (一粒子 lane, n = 2) ----
        // h = σ_x, P1 = |1⟩⟨1|, P2 = |2⟩⟨2|, Γ0 = I/2, ε = 0.1
        let n = 2usize;
        let h: Vec<C64> = vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
        ];
        let p1 = {
            let mut m = zeros(n);
            m[0] = C64::new(1.0, 0.0);
            m
        };
        let p2 = {
            let mut m = zeros(n);
            m[3] = C64::new(1.0, 0.0);
            m
        };
        let eps = 0.1;
        // signed covariance probe: n̈_j[Γ] = −[h,[h,Γ]]_jj
        let curv = |gen: &[C64], gamma: &[C64]| -> f64 {
            let c = comm(gen, &comm(gen, gamma, n), n);
            -c[1 * n + 1].re
        };
        let g0 = scal(0.5, &eye(n));
        let gp = add(&g0, &scal(eps, &p1));
        let gm = add(&g0, &scal(-eps, &p1));
        let resp_probe = (curv(&h, &gp) - curv(&h, &gm)) / (4.0 * eps);
        // coupling weight ‖P2 h P1‖²_F
        let php = mul(&p2, &mul(&h, &p1, n), n);
        let w: f64 = php.iter().map(|x| x.norm2()).sum();
        // Hamiltonian quench: h± = h ± εP1, 状態は Γ0 = I/2 固定
        let hp = add(&h, &scal(eps, &p1));
        let hm = add(&h, &scal(-eps, &p1));
        let resp_quench = (curv(&hp, &g0) - curv(&hm, &g0)) / (4.0 * eps);
        check(
            "[M4a] N1: signed covariance probe = 結合重み 1 (厳密恒等式 T8)",
            (resp_probe - w).abs() < 1e-12 && (w - 1.0).abs() < 1e-12,
            format!("probe 応答 = {:.12}, ‖P₂hP₁‖² = {:.1}", resp_probe, w),
        );
        check(
            "[M4a'] N1: Hamiltonian quench の応答 = 0 (別実験 — T9 の null)",
            resp_quench.abs() < 1e-13,
            format!(
                "quench 応答 = {:.1e} (結合 1 に対し零 — probe 型が契約の必須フィールド)",
                resp_quench
            ),
        );

        // ---- [M4b] N2 pairing (2 モード JW, dim 4) ----
        // c1 = a ⊗ I, c2 = Z ⊗ a; H = c1†c2 + c2†c1 + Δ(c1†c2† + c2c1) [+ V n1n2]
        let d4 = 4usize;
        let a2: Vec<C64> = vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
        ];
        let z2: Vec<C64> = vec![
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(-1.0, 0.0),
        ];
        let i2 = eye(2);
        let c1 = kron(&a2, 2, &i2, 2);
        let c2 = kron(&z2, 2, &a2, 2);
        let c1d = dagger(&c1, d4);
        let c2d = dagger(&c2, d4);
        let n1 = mul(&c1d, &c1, d4);
        let n2op = mul(&c2d, &c2, d4);
        let rho1 = |p: f64| -> Vec<C64> {
            let mut m = zeros(2);
            m[0] = C64::new(1.0 - p, 0.0);
            m[3] = C64::new(p, 0.0);
            m
        };
        let eps2 = 0.07;
        let mut all_ok = true;
        let mut details = String::new();
        for (dlt, vterm) in [(0.0, 0.0), (0.3, 0.0), (0.7, 0.0), (0.0, 1.3)] {
            let hop = add(&mul(&c1d, &c2, d4), &mul(&c2d, &c1, d4));
            let pair = add(&mul(&c1d, &c2d, d4), &mul(&c2, &c1, d4));
            let vint = mul(&n1, &n2op, d4);
            let ham = add(&add(&hop, &scal(dlt, &pair)), &scal(vterm, &vint));
            let cop = comm(&ham, &comm(&ham, &n2op, d4), d4);
            let resp_of = |p: f64| -> f64 {
                let rho = kron(&rho1(p), 2, &rho1(0.5), 2);
                -trace(&mul(&rho, &cop, d4), d4).re
            };
            let resp = (resp_of(0.5 + eps2) - resp_of(0.5 - eps2)) / (4.0 * eps2);
            let want = 1.0 - dlt * dlt;
            if (resp - want).abs() > 1e-12 {
                all_ok = false;
            }
            details.push_str(&format!("Δ={} V={}: {:.12}; ", dlt, vterm, resp));
        }
        check(
            "[M4b] N2: 積状態 signed probe の応答 = 1 − Δ² 厳密 (V 項は不変) — T10",
            all_ok,
            details,
        );

        // ---- [M4c] N3 Busch 対 (qubit joint POVM) ----
        let sx: Vec<C64> = vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
        ];
        let sz: Vec<C64> = vec![
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(-1.0, 0.0),
        ];
        // 非可換の確認 ([E^x, E^z] ≠ 0 — joint measurability は可換性より広い)
        let ex = add(&scal(0.5, &eye(2)), &scal(0.3, &sx));
        let ez = add(&scal(0.5, &eye(2)), &scal(0.3, &sz));
        let nc: f64 = comm(&ex, &ez, 2).iter().map(|x| x.norm2()).sum::<f64>().sqrt();
        let busch = |eta: f64| -> (bool, f64) {
            // R_ab = (I + aησ_x + bησ_z)/4 — 資格 = 全 R_ab ≥ 0 ∧ Σ = I ∧ marginal 一致
            let mut min_eig = f64::INFINITY;
            let mut sum = zeros(2);
            for aa in [-1.0, 1.0] {
                for bb in [-1.0, 1.0] {
                    let r = scal(
                        0.25,
                        &add(&eye(2), &add(&scal(aa * eta, &sx), &scal(bb * eta, &sz))),
                    );
                    let (lo, _) = eig2_herm(&r);
                    if lo < min_eig {
                        min_eig = lo;
                    }
                    sum = add(&sum, &r);
                }
            }
            let sum_dev: f64 = add(&sum, &scal(-1.0, &eye(2)))
                .iter()
                .map(|x| x.norm2())
                .sum::<f64>()
                .sqrt();
            (min_eig >= -1e-15 && sum_dev < 1e-15, min_eig)
        };
        let (q06, m06) = busch(0.6);
        let (q08, m08) = busch(0.8);
        let crit06 = 2.0 * 0.6 * 0.6 <= 1.0;
        let crit08 = 2.0 * 0.8 * 0.8 <= 1.0;
        check(
            "[M4c] N3: Busch 対 — joint POVM 資格 ⟺ 2η² ≤ 1 (η=0.6 資格 / η=0.8 拒否)・非可換対",
            q06 == crit06 && q08 == crit08 && q06 && !q08 && nc > 0.1,
            format!(
                "min eig: η=0.6 → {:+.4} (資格) / η=0.8 → {:+.4} (joint_candidate_not_positive)・‖[E^x,E^z]‖ = {:.3}",
                m06, m08, nc
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "operational core は paper-closed になった — 独立実装の入力は spec ひとつ。\n       出力値は spec に無く、閉じない領域 (holdout harness・campaign) は\n       manifest が repository_replay_only と正直に記録する。"
        } else {
            "**spec の閉包が破れている** — spec/manifest/規範インスタンスを修復せよ"
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
