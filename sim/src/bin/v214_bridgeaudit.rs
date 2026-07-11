//! v21.4 Einstein 橋の残項目監査 — 測れるものを測り、測れないものの依存連鎖を刻む (PROMPT/3 §2R)
//!
//! PROMPT/3 §2 の検査対象 7 項目の到達可能性監査:
//!   [測定 1] ローレンツ不変性の主要次数保護: staggered 分散 E²(k) = Σ sin²k_μ の
//!     方向異方性 A(k) = [E_axis(|k|) − E_diag(|k|)] / E は O(k²) で消えるか —
//!     事前登録: log-log 勾配 p = 2 ± 0.1 (異方性は O(k²a²) = 主要次数で保護)。
//!   [測定 2] ∇T = 0 の格子形: エネルギー保存 (v19.5 で ~1e-15 実測済) と粒子連続の式
//!     (v19.5 で Simpson 5 桁一致) の再掲照合 — 本監査はこの 2 つを保存則台帳として引用。
//!   [限界 3-7] Newton 1/r・重力波 2 偏極・S_eff (G,Λ,α,β) 読み出し・R² 上界・
//!     QNEC 連続極限: 全て G の読み出し (第二段) に依存 — v19.2/19.4/19.6 の三角測量で
//!     「測定可能窓 d≲2 の外」と定量化済み。依存連鎖: 1/r と 2 偏極は S_eff の
//!     Einstein-Hilbert 項の係数 (=G) の同定が前提 / R² 上界は S_eff の次の次数 /
//!     QNEC 連続極限は N≥24 級の格子系列。脱出路 (v19.6 §2): N≥24 (疎化 jacobi)・
//!     横質量ゼロ構成・多倍長核 — いずれも本監査時点で未着手 = 限界として明示。
//! 事前登録: (a) 測定 1 の p = 2 ± 0.1 = ローレンツ保護の主要次数成立 (§2 の測定可能分が
//!   閉じる) / (b) 外れ = 格子の等方性破れが低次に残る (中心テーゼへの反証材料)。

use uft_sim::*;

fn disp(kx: f64, ky: f64, kz: f64) -> f64 {
    (kx.sin().powi(2) + ky.sin().powi(2) + kz.sin().powi(2)).sqrt()
}

fn main() {
    self_test();
    println!("=== v21.4 Einstein 橋の残項目監査 — 測れるものと依存連鎖 (PROMPT/3 §2R) ===\n");
    println!("事前登録: (a) 分散異方性の log-log 勾配 p = 2 ± 0.1 = ローレンツ主要次数保護 /");
    println!("          (b) 外れ。G 依存項目 (1/r・2 偏極・S_eff・R²) は限界の依存連鎖を明示\n");
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

    // ---- [測定 1] 分散の方向異方性 ----
    // |k| を固定し、軸方向 (k,0,0) と対角方向 (k/√3, k/√3, k/√3) の E を比較
    let ks = [0.05f64, 0.1, 0.2, 0.4];
    let mut lnk = Vec::new();
    let mut lna = Vec::new();
    for &k in &ks {
        let e_axis = disp(k, 0.0, 0.0);
        let kd = k / 3.0f64.sqrt();
        let e_diag = disp(kd, kd, kd);
        let aniso = ((e_axis - e_diag) / (0.5 * (e_axis + e_diag))).abs();
        println!(
            "    |k| = {:.2}: E_axis = {:.6}, E_diag = {:.6}, 異方性 = {:.3e}",
            k, e_axis, e_diag, aniso
        );
        lnk.push(k.ln());
        lna.push(aniso.ln());
    }
    let (_ic, p) = linfit(&lnk, &lna);
    check(
        "[測定 1] 異方性の勾配 p = 2 ± 0.1 (O(k²a²) 保護)",
        (p - 2.0).abs() < 0.1,
        format!("p = {:.4}", p),
    );
    // 理論確認: E² = k² − Σk_μ⁴/3 + O(k⁶): 軸 Σk⁴ = k⁴, 対角 = k⁴/3 —
    // 異方性 ≈ (k²/6)(1 − 1/3) = k²/9 の係数チェック (記録)
    let a01 = ((disp(0.1, 0.0, 0.0)
        - disp(
            0.1 / 3.0f64.sqrt(),
            0.1 / 3.0f64.sqrt(),
            0.1 / 3.0f64.sqrt(),
        ))
        / disp(0.1, 0.0, 0.0))
    .abs();
    println!(
        "    [記録] |k|=0.1 の異方性 {:.3e} vs 解析 k²/9 = {:.3e} (小 k 展開)",
        a01,
        0.01 / 9.0
    );

    // ---- [測定 2] 保存則台帳 (v19.5 の実測を引用照合) ----
    // v19.5 (results/v195_qnec3d.txt): エネルギー保存 |ΔE| ≤ 2e-14, 連続の式 5 桁一致
    let txt = std::fs::read_to_string("results/v195_qnec3d.txt").unwrap_or_default();
    let has_e = txt.contains("のエネルギー保存");
    let has_c = txt.contains("の粒子連続の式 (Simpson)");
    check(
        "[測定 2] ∇T=0 の格子形 (v19.5 の保存則ゲートの存在照合)",
        has_e && has_c,
        format!("エネルギー保存 {}, 連続の式 {}", has_e, has_c),
    );

    // ---- 限界の依存連鎖 (機械照合: 該当文書の存在) ----
    let doc196 = std::fs::read_to_string("docs/uft-v19.6.md").unwrap_or_default();
    check(
        "[限界 3-7] G 窓の定量化文書 (v19.6 §2 の三角測量) の存在",
        doc196.contains("測定可能窓") && doc196.contains("脱出路"),
        "docs/uft-v19.6.md".into(),
    );
    println!("\n    依存連鎖 (限界の明示):");
    println!("      Newton 1/r ← S_eff の EH 項係数 G ← 第二段 (窓 d≲2 の外, v19.2/19.4/19.6)");
    println!("      重力波 2 偏極 ← 線形化 S_eff ← 同上 / R² 上界 ← S_eff 次次数 ← 同上");
    println!("      QNEC 連続極限 ← N≥24 格子系列 (jacobi 経済則: ~27h/対角化 — 要疎化)");
    println!("      脱出路 (登録済): N≥24 疎化・横質量ゼロ構成・多倍長核");

    // ---- 判定 ----
    let ok = nfail == 0;
    println!(
        "\n[判定] {}",
        if ok {
            "事前登録 (a): ローレンツ主要次数保護が成立 — §2 の測定可能分は閉じ、残 4 項目は G 窓の依存連鎖として限界化"
        } else {
            "事前登録 (b): 外れ — 記録"
        }
    );

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v21.4".into())),
        ("aniso_slope".into(), Json::Num(p)),
        ("branch_a".into(), Json::Bool(ok)),
    ]);
    let pth = write_artifact("results/v214_bridgeaudit.json", &j.render());
    println!("\n[artifact] {}", pth);

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
