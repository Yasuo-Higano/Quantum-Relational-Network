//! v9.4 QNEC 掃引の core 化 — v6.3 を共通 core (PacketRing + 読み出し) から出す
//!
//! 改良方針 §7 (単一の QRN core) の続き (v6.7 で 1 模型、v7.4 で 3 模型)。v6.3 の
//! QNEC 誤差予算バイナリは状態構築・エントロピー・T₋₋・カイラル度を全て局所実装して
//! いた。本バイナリは同じ物理を
//!   PacketRing (円環 GS + カイラル/定在波束; lib.rs の QrnModel) ×
//!   QrnState の読み出し (readout_entropy / readout_null_energy / readout_chirality)
//! だけで再構成し、v6.3 の一次ソース (results/v63_qnec_budget.json) の数値への
//! **回帰一致**を [PASS] 条件とする (N=302 のカイラル掃引・定在波対照・共動凍結比)。
//!
//! core 化の付加価値 (v6.3 に無かった検査):
//!   - purity_defect: 波束回転後も C²=C (ガウス純粋状態) が機械精度で保たれること。
//!     nocc バグ (v6.3 の正誤表) を一発で捕まえたのはこの検査だった (v6.7)。
//!   - readout_fermi_momentum: 波束は粒子数を保存するので k_F = π/2 のまま。
//!   - ユニタリー性: 掃引の最終時刻でも purity_defect < 1e-9。
//!
//! これで core 模型は 4 つ (RingChain / TfdPair / GrowingChain / PacketRing) になり、
//! 残高 7 の core 移行は「動的テンソル分解」を残すのみ。

use uft_sim::*;

const N: usize = 302;
const A0: usize = 190;
const B_END: usize = 250;
const NSTEPS: usize = 60;
const VF: f64 = 2.0;

// v6.3 の一次ソース (results/v63_qnec_budget.json, N=302) — 回帰一致の目標値
const REF_BASIC: f64 = 4.46187813099372e-05;
const REF_STRONG: f64 = -1.8673890336609765e-05;
const REF_CTRL_CHIR: f64 = 3.3690858680592092e-12;
const REF_CTRL_GAP: f64 = 6.593513202610355e-05;
const REF_COMOVING_RATIO: f64 = 9835.786533238072;
const REF_CHIR: f64 = 0.976; // results/v63_qnec_budget.txt (表示 3 桁)

struct Gaps {
    basic: f64,
    strong: f64,
}

/// 掃引 S(σ), T₋₋(σ) から最小ギャップ (v6.3 の gaps と同一の差分)
fn gaps(svals: &[f64], tmm: &[f64]) -> Gaps {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut basic = f64::INFINITY;
    let mut strong = f64::INFINITY;
    for m in (4..NSTEPS - 3).step_by(2) {
        let sp2 = (svals[m + 2] - svals[m - 2]) / 4.0;
        let spp2 = (svals[m + 2] - 2.0 * svals[m] + svals[m - 2]) / 4.0;
        let t_avg = 0.25 * tmm[m - 2] + 0.5 * tmm[m] + 0.25 * tmm[m + 2];
        let rhs = two_pi * t_avg;
        basic = basic.min(rhs - spp2);
        strong = strong.min(rhs - spp2 - 6.0 * sp2 * sp2);
    }
    Gaps { basic, strong }
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// 光的掃引: 各時刻 t=m/2 で区間 [A0−m, B_END] の S と左端の T₋₋ を core 読み出しで取る
fn sweep(model: &PacketRing) -> (Vec<f64>, Vec<f64>, f64, f64) {
    let s0 = model.init();
    let h_gs = model.ring.gs_bond_energy();
    let mut svals = Vec::new();
    let mut tmm = Vec::new();
    let mut chir = 0.0;
    let mut purity_last = 0.0;
    for m in 0..=NSTEPS {
        let st = model.evolve(&s0, m as f64 / 2.0);
        let a = A0 - m;
        let l = B_END - a + 1;
        svals.push(st.readout_entropy(a, l));
        tmm.push(st.readout_null_energy(a, VF, h_gs));
        if m == 0 {
            chir = st.readout_chirality(VF, h_gs);
        }
        if m == NSTEPS {
            purity_last = st.purity_defect();
        }
    }
    (svals, tmm, chir, purity_last)
}

/// 共動変形 (x − v_F t = 一定) での ΔS 変動幅 — カイラルなら凍結、定在波なら壊れる
fn comoving_var(model: &PacketRing) -> f64 {
    let s0 = model.init();
    let vac = model.ring.init();
    let mut ds = Vec::new();
    for &m in &[0usize, 10, 20, 30] {
        let st = model.evolve(&s0, m as f64 / 2.0);
        let a = 60 + m;
        let l = B_END - a + 1;
        ds.push(st.readout_entropy(a, l) - vac.readout_entropy(a, l));
    }
    let d0 = ds[0];
    ds.iter().map(|d| (d - d0).abs()).fold(0.0, f64::max)
}

fn main() {
    self_test();
    println!("=== v9.4 QNEC 掃引の core 化: PacketRing + 読み出しで v6.3 を再現 ===\n");
    let mk = |standing: bool| PacketRing {
        ring: RingChain { n: N },
        xc: 90.0,
        alpha: 0.5,
        sig: 5.0,
        standing,
    };
    let chiral = mk(false);
    let standing = mk(true);

    // ---- [1] core の健全性検査 (v6.3 には無かった付加価値) ----
    println!("[1] 状態の健全性 (core の語彙)");
    let s0 = chiral.init();
    let pd = s0.purity_defect();
    let kf = s0.readout_fermi_momentum();
    let kf_err = (kf - std::f64::consts::PI / 2.0).abs();
    println!("    波束回転後の purity_defect = {:.2e} (< 1e-9)  {}", pd, pass(pd < 1e-9));
    println!(
        "    k_F = π/2 (波束は粒子数を保存): 偏差 {:.2e} (< 1e-12)  {}",
        kf_err,
        pass(kf_err < 1e-12)
    );
    let ok_health = pd < 1e-9 && kf_err < 1e-12;
    println!(
        "    模型の前提 {} 項目 / 主張 {:?}",
        chiral.assumptions().len(),
        chiral.claims()
    );

    // ---- [2] カイラル掃引の回帰一致 (v6.3 N=302) ----
    println!("\n[2] カイラル掃引 (N=302, 60 步) — v6.3 の一次ソースへの回帰一致");
    let t0 = std::time::Instant::now();
    let (svals, tmm, chir, pur_end) = sweep(&chiral);
    let g = gaps(&svals, &tmm);
    let ok_basic = (g.basic - REF_BASIC).abs() < 1e-9;
    let ok_strong = (g.strong - REF_STRONG).abs() < 1e-9;
    let ok_chir = (chir - REF_CHIR).abs() < 5e-4;
    println!(
        "    基本 min gap {:+.6e} vs v6.3 {:+.6e} (|Δ|<1e-9)  {}",
        g.basic,
        REF_BASIC,
        pass(ok_basic)
    );
    println!(
        "    強形 min gap {:+.6e} vs v6.3 {:+.6e} (|Δ|<1e-9)  {}",
        g.strong,
        REF_STRONG,
        pass(ok_strong)
    );
    println!(
        "    カイラル度 {:+.4} vs v6.3 {:+.3} (±5e-4)  {}",
        chir,
        REF_CHIR,
        pass(ok_chir)
    );
    let ok_unitary = pur_end < 1e-9;
    println!(
        "    掃引最終時刻 (t=30) の purity_defect = {:.2e} (< 1e-9, ユニタリー性)  {}  ({} ms)",
        pur_end,
        pass(ok_unitary),
        t0.elapsed().as_millis()
    );

    // ---- [3] 定在波対照と共動凍結の回帰一致 ----
    println!("\n[3] 定在波対照 (非カイラル) と共動凍結");
    let (sv_st, tm_st, chir_st, _) = sweep(&standing);
    let g_st = gaps(&sv_st, &tm_st);
    let cv_ch = comoving_var(&chiral);
    let cv_st = comoving_var(&standing);
    let ratio = cv_st / cv_ch.max(1e-300);
    let ok_ctrl_chir = (chir_st - REF_CTRL_CHIR).abs() < 1e-9;
    let ok_ctrl_gap = (g_st.basic - REF_CTRL_GAP).abs() < 1e-9;
    let ok_ratio = (ratio / REF_COMOVING_RATIO - 1.0).abs() < 1e-4;
    println!(
        "    定在波のカイラル度 {:+.2e} vs v6.3 {:+.2e} (|Δ|<1e-9)  {}",
        chir_st,
        REF_CTRL_CHIR,
        pass(ok_ctrl_chir)
    );
    println!(
        "    定在波の基本 min gap {:+.6e} vs v6.3 {:+.6e} (|Δ|<1e-9)  {}",
        g_st.basic,
        REF_CTRL_GAP,
        pass(ok_ctrl_gap)
    );
    println!(
        "    共動 ΔS 変動比 (定在波/カイラル) {:.3} vs v6.3 {:.3} (±0.01%)  {}",
        ratio,
        REF_COMOVING_RATIO,
        pass(ok_ratio)
    );

    // ---- JSON artifact ----
    let all_ok = ok_health
        && ok_basic
        && ok_strong
        && ok_chir
        && ok_unitary
        && ok_ctrl_chir
        && ok_ctrl_gap
        && ok_ratio;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-CORE-003".into())),
        ("model".into(), Json::Str("PacketRing (RingChain + chiral/standing packet)".into())),
        ("regression_source".into(), Json::Str("results/v63_qnec_budget.json (N=302)".into())),
        ("min_gap_basic".into(), Json::Num(g.basic)),
        ("min_gap_strong".into(), Json::Num(g.strong)),
        ("chirality".into(), Json::Num(chir)),
        ("standing_chirality".into(), Json::Num(chir_st)),
        ("standing_min_gap".into(), Json::Num(g_st.basic)),
        ("comoving_ratio".into(), Json::Num(ratio)),
        ("purity_defect_init".into(), Json::Num(pd)),
        ("purity_defect_end".into(), Json::Num(pur_end)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v94_qneccore.json", &j.render());
    println!("\n  機械可読な結果: {}", p);

    println!("\n総合判定: {}", pass(all_ok));
    println!("\n結論: v6.3 の QNEC 掃引 (最小ギャップ・カイラル度・対照・共動凍結比) は");
    println!("      共通 core の PacketRing + 4 つの読み出しから桁落ちなく再現された。");
    println!("      core 模型は 4 つ目 — QNEC の物理も RingChain/TfdPair/GrowingChain と");
    println!("      同じ状態空間 (ガウスフェルミオン網) と同じ語彙で語れる。");
    if !all_ok {
        std::process::exit(1);
    }
}
