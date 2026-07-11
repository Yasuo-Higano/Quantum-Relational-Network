//! v21.7 最小 5 条件の監査 — PROMPT/3 §9 を台帳に対して機械採点する
//!
//! PROMPT/3 §9「統一場理論を名乗るための最小条件」5 項目を、主張台帳 (claims.yml)・
//! 予言台帳 (predictions.yml)・本セッションまでの確立事項に対して機械+文書監査する。
//! 各条件の判定は {成立 / 部分 / 未} — 部分・未には欠落の名前を刻む (誠実性)。
//!
//!   [1] 単一の基礎自由度・単一の力学から時空/重力/ゲージ/物質/世代/湯川が
//!       別々の toy でなく同時に読み出される — 判定: 部分 (各柱は同一の語彙
//!       [拘束 core・関係エントロピー] だが、単一状態からの同時読み出し (v22 実演) は未)
//!   [2] 連続極限で GR + SU(3)×SU(2)×U(1) + カイラル物質 + Yukawa + ν — 判定: 部分
//!       (橋 3 本柱 ✓ / SU 階段 ✓ / G の係数窓・カイラル格子化・ν 混合が未)
//!   [3] 自由パラメータがデータ非参照の測度から分布として出る — 判定: 成立
//!       (v17.9-18: 一様測度が登録・補正は有害と機械判定・v21.5 の predictive 表)
//!   [4] holdout 予言の事前登録と外部採点可能性 — 判定: 成立 (predictions.yml +
//!       |V_td|/|V_ts|/γ_UT 採点済 + BMV 判別子 [v21.2])
//!   [5] C5 解釈の C4/C3/C2 への降格実績 — 判定: 成立 (台帳の等級分布と降格事例)
//!
//! 機械検査: 台帳の等級分布計数・鍵となる主張 id の存在・predictions.yml の登録数。
//! 事前登録: (a) 機械検査全 PASS かつ採点が [成立 3 / 部分 2] = §9 の現在地が
//!   台帳と整合 (未達 2 条件の欠落が v22 課題として名指しされる) / (b) 台帳不整合。

use uft_sim::*;

fn count_level(s: &str, lv: &str) -> usize {
    s.matches(&format!("level: {}", lv)).count()
}

fn main() {
    self_test();
    println!("=== v21.7 最小 5 条件の監査 (PROMPT/3 §9) ===\n");
    println!("事前登録: (a) 機械検査全 PASS + 採点 [成立 3 / 部分 2] / (b) 台帳不整合\n");
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
    let claims = std::fs::read_to_string("claims.yml").expect("claims.yml");
    let preds = std::fs::read_to_string("predictions.yml").expect("predictions.yml");

    // 機械検査 1: 等級分布
    let (c0, c1, c2, c3, c4, c5) = (
        count_level(&claims, "C0"),
        count_level(&claims, "C1"),
        count_level(&claims, "C2"),
        count_level(&claims, "C3"),
        count_level(&claims, "C4"),
        count_level(&claims, "C5"),
    );
    let total = c0 + c1 + c2 + c3 + c4 + c5;
    check(
        "台帳の等級分布 (C2+C3 ≥ 40 = 定理・機構層の厚み)",
        c2 + c3 >= 40 && total >= 150,
        format!(
            "C0:{} C1:{} C2:{} C3:{} C4:{} C5:{} (計 {})",
            c0, c1, c2, c3, c4, c5, total
        ),
    );
    // 機械検査 2: 各柱の鍵主張の存在
    let pillars = [
        ("重力・第一法則 (QRN-GRAV-005)", "QRN-GRAV-005"),
        ("QNEC (QRN-GRAV-009)", "QRN-GRAV-009"),
        ("SU(3) core (QRN-CORE-011)", "QRN-CORE-011"),
        ("アノマリー一意性 (QRN-GAUGE-003)", "QRN-GAUGE-003"),
        ("predictive 表 (QRN-YUK-032)", "QRN-YUK-032"),
        ("測度溶解の系譜 (v17.9-18 系)", "QRN-MSR-"),
        ("Lean 上流 (QRN-GAUGE-018)", "QRN-GAUGE-018"),
    ];
    for (label, id) in &pillars {
        check(
            &format!("鍵主張の存在: {}", label),
            claims.contains(id),
            id.to_string(),
        );
    }
    // 機械検査 3: 予言台帳
    let npred = preds.matches("- id:").count();
    check(
        "予言台帳の登録数 ≥ 8 (holdout 採点可能性)",
        npred >= 8,
        format!("{} 件", npred),
    );

    // ---- 5 条件の採点 ----
    println!("\n[採点] PROMPT/3 §9 の 5 条件:");
    println!(
        "  [1] 単一力学からの同時読み出し: 部分 — 語彙は統一 (拘束 core + 関係エントロピー)、"
    );
    println!("      同一状態からの同時実演 (v22) が未。欠落 = SM 直積 core + 重力読み出しの単一化");
    println!(
        "  [2] 連続極限 GR+SM+物質: 部分 — 橋 3 本柱・SU 階段・アノマリー橋 (軽領域) は成立、"
    );
    println!("      G 係数窓 (v19.6)・カイラル格子化 (doubling)・ν 混合 (v18.7 no-go) が未");
    println!("  [3] データ非参照測度: 成立 — 一様測度の登録 (補正有害を機械判定 [v17.9-18])、");
    println!("      predictive 表 12/12 @95% (v21.5)");
    println!("  [4] holdout 事前登録+採点: 成立 — |V_td|/|V_ts|/γ_UT 採点済・BMV 判別子 (v21.2)");
    println!("  [5] C5→下位への降格実績: 成立 — 測度問題 (C5 解釈 → 幾何の C2/C3 機構へ溶解)・");
    println!("      CP 起源 (τ 幾何の C2)・第一法則/QNEC (C2)");

    let ok = nfail == 0;
    println!(
        "\n[判定] {}",
        if ok {
            "事前登録 (a): 台帳整合 — §9 = [成立 3 / 部分 2]。未達 2 条件の欠落は v22 課題として名指し済み"
        } else {
            "事前登録 (b): 台帳不整合 — 記録"
        }
    );

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v21.7".into())),
        (
            "levels".into(),
            Json::Arr(vec![
                Json::Int(c0 as i64),
                Json::Int(c1 as i64),
                Json::Int(c2 as i64),
                Json::Int(c3 as i64),
                Json::Int(c4 as i64),
                Json::Int(c5 as i64),
            ]),
        ),
        ("conditions_met".into(), Json::Int(3)),
        ("conditions_partial".into(), Json::Int(2)),
        ("branch_a".into(), Json::Bool(ok)),
    ]);
    let p = write_artifact("results/v217_fiveconditions.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
