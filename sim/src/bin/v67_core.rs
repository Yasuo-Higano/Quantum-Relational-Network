//! v6.7 単一 QRN core からの複数読み出し — 「同じ状態」から幾何・エントロピー・物質・因果
//!
//! これまでの 32 本は個別の toy model の集まりだった。統一理論としての説得力は
//! 新しいシミュレーションを増やすことではなく、**既存のシミュレーションを同じ core
//! から出すこと**で上がる (改良方針 §7)。lib.rs に共通状態 (QrnState = ガウス
//! フェルミオン網) と模型トレイト (QrnModel) と読み出し群を定義し、本バイナリは
//! 一つの具体模型 RingChain (円環自由フェルミオン・半充填基底状態) の **同一の状態**
//! から 4 つの読み出しを行う:
//!   [G] 幾何:      MI ブロック → MDS → 円環 (v0.7 の再現; v6.4 の判定規準)
//!   [S] エントロピー: S(ℓ) の対数則 → 中心電荷 c (理論値 1; v0.5 と同型)
//!   [M] 物質:      密度 → フェルミ運動量 k_F = π⟨n⟩ (ラッティンジャー) → v_F = 2 sin k_F
//!   [C] 因果:      局所クエンチ (1 ボンド切断の基底状態から発展) → 前線速度 (v1.1 と同型)
//! さらに「物質と因果の握手」: 前線速度 ≈ v_F = 2 sin(k_F) という読み出し間の
//! 整合性を検査する。幾何・物質・因果のどれも公理には入っていない — 全て
//! 一つの相関行列とユニタリー発展からの読み出しである。

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
    println!(
        "=== v6.7 単一 QRN core: 一つの状態から幾何・エントロピー・物質・因果を読み出す ===\n"
    );
    let model = RingChain { n: 202 };
    println!(
        "[模型] RingChain (N={}) の仮定 (QrnModel::assumptions):",
        202
    );
    for a in model.assumptions() {
        println!("   - {}", a);
    }
    println!("   支える主張: {:?}\n", model.claims());
    let st = model.init();

    // ---- [0] core の健全性: 基底状態は純粋なガウス状態 (C² = C) ----
    let n = st.n;
    {
        let c2 = matmul(&st.cre, &st.cre, n);
        let mut dmax: f64 = 0.0;
        for i in 0..n * n {
            dmax = dmax.max((c2[i] - st.cre[i]).abs());
        }
        println!(
            "[0] 状態の健全性: ‖C²−C‖_max = {:.1e} (純粋状態の射影子)  {}",
            dmax,
            pass(dmax < 1e-9)
        );
    }

    // ---- [G] 幾何の読み出し ----
    let (mi, nb, mi_max) = st.readout_mi_blocks();
    let g = readout_ring_geometry(&mi, nb, mi_max);
    println!(
        "\n[G] 幾何: 隣接復元率 {:.1}% / 半径ばらつき {:.1}% / λ2/λ1 = {:.2} → {}  {}",
        100.0 * g.adjacency,
        100.0 * g.rsd,
        g.lam21,
        if g.ring { "円環" } else { "失敗" },
        pass(g.ring && g.adjacency > 0.99)
    );

    // ---- [S] エントロピーの読み出し: 対数則 → 中心電荷 ----
    let (mut xs, mut ys) = (Vec::new(), Vec::new());
    for l in (6..=50).step_by(4) {
        let chord =
            (n as f64 / std::f64::consts::PI) * (std::f64::consts::PI * l as f64 / n as f64).sin();
        xs.push(chord.ln() / 3.0);
        ys.push(st.readout_entropy(10, l));
    }
    let (_, c_fit) = linfit(&xs, &ys);
    let ok_c = (c_fit - 1.0).abs() < 0.05;
    println!(
        "[S] エントロピー: S(ℓ) = (c/3)ln[(N/π)sin(πℓ/N)] + 定数 のフィット → c = {:.3} (理論値 1)  {}",
        c_fit,
        pass(ok_c)
    );

    // ---- [M] 物質の読み出し: 密度 → k_F → v_F ----
    let kf = st.readout_fermi_momentum();
    let vf = 2.0 * kf.sin();
    let ok_kf = (kf - std::f64::consts::PI / 2.0).abs() < 1e-6;
    println!(
        "[M] 物質: k_F = π⟨n⟩ = {:.6} (理論値 π/2 = {:.6}) → v_F = 2 sin k_F = {:.4}  {}",
        kf,
        std::f64::consts::PI / 2.0,
        vf,
        pass(ok_kf)
    );

    // ---- [C] 因果の読み出し: 局所クエンチの前線速度 ----
    // 1 ボンド (n-1, 0) を切った開鎖の基底状態から、円環 H でユニタリー発展させる。
    // 密度の乱れが切断点から広がる速度 = 光円錐 (Lieb–Robinson)。
    let st_open = {
        let mut h = vec![0.0; n * n];
        for x in 0..n - 1 {
            h[x + (x + 1) * n] = -1.0;
            h[(x + 1) + x * n] = -1.0;
        }
        let (_, v) = jacobi_eigh(&h, n);
        let nocc = n / 2;
        let mut cre = vec![0.0; n * n];
        for m in 0..nocc {
            for i in 0..n {
                let vi = v[i + m * n];
                if vi == 0.0 {
                    continue;
                }
                for j in 0..n {
                    cre[i + j * n] += vi * v[j + m * n];
                }
            }
        }
        QrnState {
            n,
            cre,
            cim: vec![0.0; n * n],
        }
    };
    // 観測量の選択 (物理ノート): 二部格子+半充填のカイラル対称性により、局所密度は
    // このクエンチで厳密に不変 (Δρ ≡ 0 を数値でも確認)。伝播を運ぶのはボンド運動
    // エネルギー −2C(x,x+1) なので、その変化で前線を測る。
    let bond0: Vec<f64> = (0..n - 1).map(|x| st_open.cre[x + (x + 1) * n]).collect();
    let (mut ts, mut fronts) = (Vec::new(), Vec::new());
    println!("[C] 因果: 局所クエンチ (切断点 x=0/N-1) のボンドエネルギー前線");
    println!("    t     前線距離");
    for &t in &[6.0f64, 10.0, 14.0, 18.0, 22.0] {
        let stt = model.evolve(&st_open, t);
        let mut front = 0usize;
        for x in 0..n - 1 {
            let d = (x + 1).min(n - 1 - x);
            if (stt.cre[x + (x + 1) * n] - bond0[x]).abs() > 1e-3 {
                front = front.max(d);
            }
        }
        println!("    {:4.0}  {:5}", t, front);
        ts.push(t);
        fronts.push(front as f64);
    }
    let (_, v_front) = linfit(&ts, &fronts);
    let ok_front = (v_front - 2.0).abs() < 0.2;
    println!(
        "    => 前線速度 = {:.3} (理論値 v_max = 2)  {}",
        v_front,
        pass(ok_front)
    );

    // ---- 読み出し間の握手: 因果の速度 = 物質の v_F ----
    let ok_hand = (v_front - vf).abs() < 0.2;
    println!(
        "\n[握手] 因果の前線速度 {:.3} ≈ 物質読み出しの v_F = 2 sin(k_F) = {:.3}  {}",
        v_front,
        vf,
        pass(ok_hand)
    );
    println!(
        "       (光円錐の速さは物質の分散の読み出しと一致する — 別々の読み出しが同じ網を見ている)"
    );

    // ---- JSON artifact ----
    let all_ok = g.ring && ok_c && ok_kf && ok_front && ok_hand;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-CORE-001".into())),
        ("model".into(), Json::Str("RingChain".into())),
        ("n".into(), Json::Int(n as i64)),
        (
            "geometry".into(),
            Json::Obj(vec![
                ("adjacency".into(), Json::Num(g.adjacency)),
                ("radius_scatter".into(), Json::Num(g.rsd)),
                ("ring".into(), Json::Bool(g.ring)),
            ]),
        ),
        ("central_charge".into(), Json::Num(c_fit)),
        ("fermi_momentum".into(), Json::Num(kf)),
        ("front_velocity".into(), Json::Num(v_front)),
        ("vf_from_matter".into(), Json::Num(vf)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v67_core.json", &j.render());
    println!("\n  機械可読な結果: {}", p);

    println!("\n総合判定: {}", pass(all_ok));
    println!("\n結論: 幾何 (円環)・エントロピー (c=1 の対数則)・物質 (k_F)・因果 (光円錐) が、");
    println!("      一つの相関行列と一つのユニタリー発展 — 同じ量子情報網 — の 4 つの読み出し");
    println!("      として得られ、読み出し同士が整合する (前線速度 = 2 sin k_F)。");
    println!("      既存の各シミュレーションを QrnModel 実装として core に移す道が開いた");
    println!("      (移行地図は docs/uft-v6.7.md)。");
    if !all_ok {
        std::process::exit(1);
    }
}
