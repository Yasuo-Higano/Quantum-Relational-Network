//! v23.5 λ∞ の外挿判定 — 「→1」と「→λ*」を正しく判別する (第二十四期)
//!
//! v23.4 の登録基準 (|λ−1| の族内単調減少) は λ → 1 と λ → λ* ≠ 1 を区別できない
//! 設計不足だった (開発記録)。本版が判別基準を再登録する: λ(N) = λ∞ + c/N⁴ の
//! 族別 Richardson 外挿 (mod-4 shell 族: {8,12,16} と {10,14}) で λ∞ を決め、
//! 1 との比較で分岐する。一次データ = results/v234_lambda.json (v23.4 の測定)。
//!
//! 装置ゲート: [J1] 1/N⁴ 則の自己整合 — mod-0 族の 2 通りの外挿が ±0.002 で一致,
//!   [J2] 族間整合 — |λ∞(mod0) − λ∞(mod2)| < 0.003 (shell 系統の許容)。
//! 事前登録: (a) λ∞ ∈ 1 ± 0.02 = 連続極限で BW 素朴回復 (温度 2π) /
//!   (a′) λ∞ が 1 ± 0.02 の外 = **有限繰り込みの確定** — 格子エンタングルメント
//!   温度は 2π/λ∞。記録: 32/27 との照合 / (b) [J1/J2] 破れ = スケーリング不成立。

use uft_sim::*;

fn json_num(s: &str, key: &str, occurrence: usize) -> f64 {
    let pat = format!("\"{}\":", key);
    let mut idx = 0usize;
    let mut from = 0usize;
    loop {
        let i = s[from..]
            .find(&pat)
            .unwrap_or_else(|| panic!("json key {}", key))
            + from;
        if idx == occurrence {
            let rest = &s[i + pat.len()..];
            let end = rest
                .find(|c| c == ',' || c == '}' || c == '\n')
                .unwrap_or(rest.len());
            return rest[..end].trim().parse().expect("parse");
        }
        idx += 1;
        from = i + pat.len();
    }
}

fn main() {
    self_test();
    println!("=== v23.5 λ∞ の外挿判定 — 有限繰り込みか素朴回復か (第二十四期) ===\n");
    println!("事前登録: (a) λ∞ = 1 ± 0.02 = BW 素朴回復 / (a′) 外 = 有限繰り込み 2π/λ∞ /");
    println!("          (b) 1/N⁴ 自己整合・族間整合の破れ\n");
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
    let js = std::fs::read_to_string("results/v234_lambda.json").expect("results/v234_lambda.json");
    // rows は N 昇順 (8,10,12,14,16) — lambda の出現順で取得
    let ns = [8usize, 10, 12, 14, 16];
    let lam: Vec<f64> = (0..5).map(|i| json_num(&js, "lambda", i)).collect();
    for i in 0..5 {
        println!("    N={:2}: λ = {:.5}", ns[i], lam[i]);
    }
    // 1/N⁴ 外挿: λ∞ = (λ2·N2⁴ − λ1·N1⁴)/(N2⁴ − N1⁴)
    let extrap = |n1: usize, l1: f64, n2: usize, l2: f64| -> f64 {
        let p1 = (n1 as f64).powi(4);
        let p2 = (n2 as f64).powi(4);
        (l2 * p2 - l1 * p1) / (p2 - p1)
    };
    let e_8_12 = extrap(8, lam[0], 12, lam[2]);
    let e_12_16 = extrap(12, lam[2], 16, lam[4]);
    let e_10_14 = extrap(10, lam[1], 14, lam[3]);
    println!(
        "    外挿: mod-0 族 (8,12) → {:.5} / (12,16) → {:.5} / mod-2 族 (10,14) → {:.5}",
        e_8_12, e_12_16, e_10_14
    );
    check(
        "[J1] 1/N⁴ 自己整合 (mod-0 の 2 外挿が ± 0.002)",
        (e_8_12 - e_12_16).abs() < 0.002,
        format!("Δ = {:.4}", (e_8_12 - e_12_16).abs()),
    );
    let lam_inf0 = e_12_16; // 大 N 側を採用
    check(
        "[J2] 族間整合 |λ∞(mod0) − λ∞(mod2)| < 0.003",
        (lam_inf0 - e_10_14).abs() < 0.003,
        format!(
            "mod0 {:.5} vs mod2 {:.5} (Δ {:.4})",
            lam_inf0,
            e_10_14,
            (lam_inf0 - e_10_14).abs()
        ),
    );
    let lam_inf = 0.5 * (lam_inf0 + e_10_14);
    let naive = (lam_inf - 1.0).abs() < 0.02;
    println!(
        "\n    λ∞ = {:.5} (族平均) — 32/27 = {:.5} との比 = {:.5}",
        lam_inf,
        32.0 / 27.0,
        lam_inf / (32.0 / 27.0)
    );
    println!(
        "\n[判定] {}",
        if nfail > 0 {
            "事前登録 (b): スケーリング不成立 — 記録"
        } else if naive {
            "事前登録 (a): λ∞ = 1 — 連続極限で BW 素朴回復"
        } else {
            "事前登録 (a′): 有限繰り込みの確定 — 格子エンタングルメント温度は 2π/λ∞。G 読み出し (v23.6) はこの温度で行う"
        }
    );

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v23.5".into())),
        ("lam_inf".into(), Json::Num(lam_inf)),
        ("lam_inf_mod0".into(), Json::Num(lam_inf0)),
        ("lam_inf_mod2".into(), Json::Num(e_10_14)),
        ("ratio_32_27".into(), Json::Num(lam_inf / (32.0 / 27.0))),
        ("branch_a_naive".into(), Json::Bool(nfail == 0 && naive)),
        ("branch_ap_renorm".into(), Json::Bool(nfail == 0 && !naive)),
    ]);
    let p = write_artifact("results/v235_lambdainf.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(a′)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
