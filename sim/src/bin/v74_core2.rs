//! v7.4 core 移行 (2): TFD と成長する鎖を同じ QRN core から読み出す (v7.0 残高 7)
//!
//! v6.7 は RingChain 1 模型で「同じ状態から幾何・エントロピー・物質・因果」を実演した。
//! 本バイナリは v1.2 (ER=EPR) と v5.1 (成長する宇宙) を lib.rs の core (QrnState /
//! QrnModel / 共有読み出し) の実装として再現する — toy の集まりが一つの語彙に揃っていく。
//!
//! [A] TfdPair (v1.2 の再現):
//!     - 純粋性 (C²=C) — TFD は「もつれで熱を装う」純粋状態
//!     - 片側の局所観測量が厳密に熱的 (フェルミ分布)
//!     - 鏡像ブロック間 MI が β とともに単調に落ちる (もつれを絞ると接続が千切れる)
//!       v1.2 の値: β=0.02 で 10 サイトブロック MI ≈ 13.86 (≈ 20 ln 2)、β=20 で ≈ 0
//!     - H_L+H_R の発展で TFD が定常 (ΔC ≈ 0) — ER=EPR の「橋」は動力学的にも安定
//! [B] GrowingChain (v5.1 の再現):
//!     - シナリオ A (真空で到着): 窓のエントロピーが成長し、全系は純粋のまま
//!     - シナリオ B (熱的に到着): 全系の純粋性が壊れる (対照)
//!     - v6.7 の教訓を核に: purity_defect (C²=C) を全ステップで検査する

use uft_sim::*;

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

fn main() {
    self_test();
    println!("=== v7.4 core 移行 (2): TFD と成長する鎖 — 同じ core、別の宇宙 ===\n");
    let mut checks: Vec<(String, bool)> = Vec::new();
    let record = |name: &str, ok: bool, checks: &mut Vec<(String, bool)>| {
        println!("  => {}  {}", name, pass(ok));
        checks.push((name.to_string(), ok));
    };

    // ---- [A] TfdPair ----
    println!("[A] TfdPair (v1.2 = ER=EPR を core の語彙で)");
    let n = 64usize;
    println!("    仮定: {:?}", TfdPair { n, beta: 1.0 }.assumptions());
    let betas = [0.02f64, 0.5, 2.0, 20.0];
    let mut mis = Vec::new();
    let mut max_defect: f64 = 0.0;
    let mut max_thermal_err: f64 = 0.0;
    for &beta in &betas {
        let model = TfdPair { n, beta };
        let st = model.init();
        max_defect = max_defect.max(st.purity_defect());
        // 片側 (L) の局所密度: 厳密に 1/2 (半充填の熱状態)
        let dens: f64 = (0..n).map(|i| st.cre[i + i * st.n]).sum::<f64>() / n as f64;
        max_thermal_err = max_thermal_err.max((dens - 0.5).abs());
        // 鏡像 10 サイトブロック間の MI (v1.2 と同じ観測量)
        let ra: Vec<usize> = (10..20).collect();
        let rb: Vec<usize> = (n + 10..n + 20).collect();
        let mi = st.readout_mi_regions(&ra, &rb);
        mis.push(mi);
        println!(
            "    β={:5.2}: 鏡像ブロック MI = {:7.4}  (純粋性欠陥 {:.1e})",
            beta,
            mi,
            st.purity_defect()
        );
    }
    let max_mi_expected = 20.0 * (2.0f64).ln(); // 10 サイト × 2 ln 2
    record(
        &format!("TFD は全 β で厳密に純粋 (欠陥 {:.1e} < 1e-9)", max_defect),
        max_defect < 1e-9,
        &mut checks,
    );
    record(
        &format!(
            "高温極限で MI → 2 ln2/サイト (観測 {:.2} / 理論 {:.2})、β で単調減少、低温で ≈ 0",
            mis[0], max_mi_expected
        ),
        (mis[0] - max_mi_expected).abs() < 0.15
            && mis.windows(2).all(|w| w[1] < w[0])
            && mis[3] < 0.3,
        &mut checks,
    );
    // 動力学: TFD は H_L+H_R で定常
    {
        let model = TfdPair { n, beta: 2.0 };
        let st = model.init();
        let st2 = model.evolve(&st, 7.0);
        let mut dmax: f64 = 0.0;
        for i in 0..st.cre.len() {
            dmax = dmax.max((st.cre[i] - st2.cre[i]).abs());
            dmax = dmax.max((st.cim[i] - st2.cim[i]).abs());
        }
        record(
            &format!(
                "TFD は H_L+H_R の発展で定常 (t=7 で ΔC = {:.1e}) — 橋は動力学的に安定",
                dmax
            ),
            dmax < 1e-9,
            &mut checks,
        );
    }
    println!("    (v1.2 の結果と同型: 高温 β=0.02 で MI={:.2} ≈ v1.2 の 13.86 — もつれを絞ると接続が千切れる)", mis[0]);

    // ---- [B] GrowingChain ----
    println!("\n[B] GrowingChain (v5.1 = 成長する宇宙を core の語彙で)");
    let gc = GrowingChain { n_max: 96 };
    let window: Vec<usize> = (0..16).collect();
    // シナリオ A: 2 サイトずつ・局所基底状態 (純粋な「真空」) で到着 — v5.1 と同一プロトコル
    let mut st = gc.init(16);
    let mut s_hist = Vec::new();
    let mut defect_a: f64 = st.purity_defect();
    let s_of = |st: &QrnState| -> f64 {
        let k = window.len();
        let mut cre = vec![0.0; k * k];
        let mut cim = vec![0.0; k * k];
        for (a, &sa) in window.iter().enumerate() {
            for (b, &sb) in window.iter().enumerate() {
                cre[a + b * k] = st.cre[sa + sb * st.n];
                cim[a + b * k] = st.cim[sa + sb * st.n];
            }
        }
        entropy_corr_herm(&cre, &cim, k)
    };
    s_hist.push(s_of(&st));
    let mut active = 16usize;
    while active < 96 {
        st = gc.arrive_pair_vacuum(&st, active);
        active += 2;
        st = gc.evolve_active(&st, active, 3.0);
        defect_a = defect_a.max(st.purity_defect());
        s_hist.push(s_of(&st));
    }
    let ups = s_hist.windows(2).filter(|w| w[1] > w[0]).count();
    println!(
        "    シナリオ A (真空到着): 窓 S: {:.3} → {:.3} (上昇 {}/{} ステップ), 全系純粋性欠陥 {:.1e}",
        s_hist[0],
        s_hist.last().unwrap(),
        ups,
        s_hist.len() - 1,
        defect_a
    );
    record(
        "シナリオ A: 全系は純粋なまま (欠陥 < 1e-9) 窓のエントロピーが成長 (S > 4, 上昇 ≥ 70%)",
        defect_a < 1e-9
            && *s_hist.last().unwrap() > 4.0
            && ups as f64 >= 0.7 * (s_hist.len() - 1) as f64,
        &mut checks,
    );
    // シナリオ B: 熱的到着 (対照)
    let mut stb = gc.init(16);
    stb = gc.arrive_pair_thermal(&stb, 16);
    stb = gc.evolve_active(&stb, 18, 3.0);
    let defect_b = stb.purity_defect();
    println!(
        "    シナリオ B (熱的到着, 対照): 1 ステップ目で純粋性欠陥 {:.2} (混合状態)",
        defect_b
    );
    record(
        "シナリオ B (対照): 熱的到着は全系の純粋性を壊す (欠陥 > 0.1)",
        defect_b > 0.1,
        &mut checks,
    );
    println!(
        "    (v5.1 と同じ機構: 低エントロピーの過去 = 自由度が少ない過去 + 真空で生まれる新モード)"
    );

    // ---- JSON / 判定 ----
    let all_ok = checks.iter().all(|(_, ok)| *ok);
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-CORE-002".into())),
        (
            "tfd_mirror_mi_by_beta".into(),
            Json::Arr(mis.iter().map(|&m| Json::Num(m)).collect()),
        ),
        (
            "tfd_betas".into(),
            Json::Arr(betas.iter().map(|&b| Json::Num(b)).collect()),
        ),
        (
            "growing_window_entropy".into(),
            Json::Arr(s_hist.iter().map(|&s| Json::Num(s)).collect()),
        ),
        ("growing_purity_defect".into(), Json::Num(defect_a)),
        ("thermal_control_defect".into(), Json::Num(defect_b)),
        (
            "checks".into(),
            Json::Arr(
                checks
                    .iter()
                    .map(|(nm, ok)| {
                        Json::Obj(vec![
                            ("name".into(), Json::Str(nm.clone())),
                            ("pass".into(), Json::Bool(*ok)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v74_core2.json", &j.render());
    println!("\n  機械可読な結果: {}", p);

    println!("\n---- 検査一覧 ----");
    for (nm, ok) in &checks {
        println!("  {} {}", pass(*ok), nm);
    }
    println!("\n総合判定: {}", pass(all_ok));
    println!("\n結論: ER=EPR (v1.2) と成長する宇宙 (v5.1) が、RingChain (v6.7) と同じ");
    println!("      QrnState・同じ読み出し・同じ健全性検査 (C²=C) の上で動いた。");
    println!("      core 移行地図: RingChain ✓ / TfdPair ✓ / GrowingChain ✓ —");
    println!("      「toy の集まり」から「一つの網の異なる模型」への移行が進んでいる。");
    if !all_ok {
        std::process::exit(1);
    }
}
