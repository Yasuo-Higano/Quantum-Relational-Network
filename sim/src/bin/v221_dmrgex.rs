//! v22.1 監査層: 重電荷ボソンの DMRG 励起状態を全空間 ED アンカーで一次ソース化
//!
//! 第二十三期 課題 3 (v21.1 の重電荷窓を開く)。v21.1 の環 ED では巻き付き ε の
//! フラックス状態が低エネルギー帯を占拠し、電荷 q ≥ 2 のボソン質量が読めなかった。
//! 開鎖には ε が無い (厳密) — E₁ − E₀ が中性セクターのボソン質量そのもの。
//! 計算層 (v221_dmrg.rs — Rust DMRG, 直交罰則, χ=128) の N=64 ギャップを、
//! N=10 全空間 ED (dim 1024, セクター仮定なし — 罰則法の大域 E₁ と同土俵) で照合。
//!
//! 模型 (探索層と同一規約): H = −x Σ (c†_l c_{l+1} + h.c.) + Σ_links L_n² + λ_Q (Σ_k q_k)²,
//! L_n = Σ_{k≤n} q_k, q_k = qf·(n_k − [k odd])。隣接ホップは JW 符号なし (開鎖)。
//! 教科書値 (ボソン化): M(q)/g = q/√π — 比 M(2)/M(1) = 2 が連続極限の判定線。
//! 開鎖の箱運動量 p/g ≈ π√x/N による既知の下方偏差 (x=16, N=64 で比 ≈ 1.92) は
//! 窓 ±7% の内側に事前見積もりとして織り込む。
//!
//! 開発記録 (run1, 保存): λ_Q = 0 の run1 は判定 (b) — 全空間 E₁ が ⟨n⟩ = 4.000 の
//! 荷電状態だった。開鎖の端の電荷はどのリンクにも寄与せず自由 (最終サイトの q は
//! L_n に入らない) — 荷電端状態がボソンより下に潜り E₁ を偽装する。環の磁束状態
//! (v21.1) の開鎖版の教訓。修正 = 全電荷罰則 λ_Q (Σq)² (ペア重みに +λ_Q)。
//! 中性セクターのスペクトルは厳密に不変 — これ自体を [I3] ゲートにする。
//!
//! 装置ゲート (q ∈ {1,2}): [I1] E₀(ED) = anchor ± 1e-4, [I2] E₁(ED) = anchor ± 1e-3
//! (罰則法 λ_o=50 の系統許容), [I3] 罰則 H の E₀/E₁ = 中性セクター制限 ED (λ_Q=0) の
//! E₀/E₁ ± 1e-9 (罰則の中性不変性 = 厳密恒等式) ∧ E₁ の ⟨n⟩ = N/2。
//!
//! 開発記録 — 二層 (Rust/python) × 三推定器の完結 (計算層 v221_dmrg run1-10 全記録保存):
//!   q=2 の gap は両エンジン一致 (Δ0.1%)。q=1 の E₁ は低スペクトルが密 (多ボソン梯子)
//!   で Rust 罰則法 8 走が変分上界どまり — python (31h) が収束させた。E₀ は全 4 点で
//!   交差エンジン厳密一致 (Δ ≤ 4e-9)。第三推定器 = 接続電荷相関の有効質量 ([G8] で
//!   gap 法と ±20% 整合、バイアス +10-14% は q に一様) が独立に同じ比を与える。
//! 一次データ: results/v221_dmrg_gaps.json (Rust: アンカー・全 E₀・q=2 gap・m_eff) +
//!   explore/dmrg_heavy_out_pyq1.json (python: q=1 gap — 31h 走の収束値)。
//! 装置ゲート追加: [I4] q=1 E₀ の交差エンジン一致 ± 1e-6 (実測 Δ ≤ 4e-9)。
//! 事前登録 (判定量は当初から不変 = x=16 の質量比 2 ± 7%, gap 法):
//!   (a) 装置ゲート全 PASS ∧ x=16 の gap 比 (q=2)/(q=1) = 2 ± 7% = 重電荷窓の開通 /
//!   (b) 外れ = 記録。記録: x=9 比・相関質量比・q=1 の Rust E₁ 上界。

use uft_sim::*;

const N_ED: usize = 10;
const LAMQ: f64 = 10.0; // 全電荷罰則 λ_Q (荷電端状態の排除 — run1 の教訓)

// 全空間 (dim 2^N) の密 H を構築 — セクター仮定なし, 罰則 λ_Q (Σq)² 込み
fn build_h(x: f64, qf: f64, lamq: f64) -> Vec<f64> {
    let dim = 1usize << N_ED;
    let mut h = vec![0.0f64; dim * dim];
    for s in 0..dim {
        // 対角: 電場 Σ_{links} L² + λ_Q (Σ_k q_k)²
        let mut e = 0.0;
        let mut acc = 0.0;
        for k in 0..N_ED {
            let occ = ((s >> k) & 1) as f64;
            let bg = if k % 2 == 1 { 1.0 } else { 0.0 };
            acc += qf * (occ - bg);
            if k < N_ED - 1 {
                e += acc * acc;
            }
        }
        e += lamq * acc * acc;
        h[s * dim + s] += e;
        // ホップ: 隣接 JW 符号なし
        for k in 0..N_ED - 1 {
            let b0 = (s >> k) & 1;
            let b1 = (s >> (k + 1)) & 1;
            if b0 != b1 {
                let t = s ^ (1 << k) ^ (1 << (k + 1));
                h[s * dim + t] += -x;
            }
        }
    }
    h
}

fn json_num(s: &str, key: &str) -> f64 {
    let pat = format!("\"{}\":", key);
    let i = s.find(&pat).unwrap_or_else(|| panic!("json key {}", key));
    let rest = &s[i + pat.len()..];
    let end = rest
        .find(|c| c == ',' || c == '}' || c == '\n')
        .unwrap_or(rest.len());
    rest[..end].trim().parse().expect("parse")
}

fn main() {
    self_test();
    println!("=== v22.1 重電荷ボソンの DMRG 励起状態 — ED アンカー照合 (第二十三期 課題 3) ===\n");
    println!("事前登録: (a) 装置ゲート PASS ∧ x=16 の M/g 比 (q=2)/(q=1) = 2 ± 7% /");
    println!("          (b) 外れ = 記録。教科書線: M(q)/g = q/√π (ボソン化)\n");
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
    let t0 = std::time::Instant::now();
    // 一次データ = Rust 計算層 (v221_dmrg, PROMPT/4)。python 走行の生き残り
    // (run2 の q=2 対) はクロスエンジン照合として文書に記録。
    let js = std::fs::read_to_string("results/v221_dmrg_gaps.json")
        .expect("results/v221_dmrg_gaps.json");

    // ---- 装置ゲート: N=10, x=2 — 罰則 H の全空間 ED と DMRG アンカーの照合 + 中性不変性 ----
    let dim = 1usize << N_ED;
    for qf in [1.0f64, 2.0] {
        let h = build_h(2.0, qf, LAMQ);
        let (ev, vv) = jacobi_eigh(&h, dim);
        // E₀/E₁ の全粒子数 (罰則後は両方 N/2 のはず)
        let mut n0 = 0.0;
        let mut n1 = 0.0;
        for i in 0..dim {
            let c0 = vv[i];
            let c1 = vv[i + dim];
            let pc = i.count_ones() as f64;
            n0 += c0 * c0 * pc;
            n1 += c1 * c1 * pc;
        }
        // 中性セクター制限 ED (λ_Q = 0) — 罰則の中性不変性の対照
        let h0 = build_h(2.0, qf, 0.0);
        let neut: Vec<usize> = (0..dim)
            .filter(|s| s.count_ones() as usize == N_ED / 2)
            .collect();
        let m = neut.len();
        let mut hn = vec![0.0f64; m * m];
        for a in 0..m {
            for b in 0..m {
                hn[a * m + b] = h0[neut[a] * dim + neut[b]];
            }
        }
        let (evn, _) = jacobi_eigh(&hn, m);
        let a0 = json_num(&js, &format!("anchor_n10_q{}_e0", qf as i64));
        let a1 = json_num(&js, &format!("anchor_n10_q{}_e1", qf as i64));
        check(
            &format!("[I1] q={} E₀: ED = DMRG ± 1e-4", qf as i64),
            (ev[0] - a0).abs() < 1e-4,
            format!(
                "ED {:.8} vs DMRG {:.8} (Δ {:.1e}, {} s)",
                ev[0],
                a0,
                (ev[0] - a0).abs(),
                t0.elapsed().as_secs()
            ),
        );
        check(
            &format!("[I2] q={} E₁: ED = DMRG ± 1e-3 (罰則法許容)", qf as i64),
            (ev[1] - a1).abs() < 1e-3,
            format!(
                "ED {:.8} vs DMRG {:.8} (Δ {:.1e}, gap ED {:.6})",
                ev[1],
                a1,
                (ev[1] - a1).abs(),
                ev[1] - ev[0]
            ),
        );
        check(
            &format!(
                "[I3] q={} 罰則の中性不変性 (厳密恒等式) ∧ E₁ 中性",
                qf as i64
            ),
            (ev[0] - evn[0]).abs() < 1e-9
                && (ev[1] - evn[1]).abs() < 1e-9
                && (n1 - (N_ED / 2) as f64).abs() < 1e-6,
            format!(
                "ΔE₀ {:.1e}, ΔE₁ {:.1e} (中性 dim {}), ⟨n⟩(E₀/E₁) = {:.3}/{:.3}",
                (ev[0] - evn[0]).abs(),
                (ev[1] - evn[1]).abs(),
                m,
                n0,
                n1
            ),
        );
    }

    // ---- [I4] q=1 E₀ の交差エンジン一致 + 本測定 ----
    let jpy = std::fs::read_to_string("explore/dmrg_heavy_out_pyq1.json")
        .expect("explore/dmrg_heavy_out_pyq1.json");
    println!();
    for x in [9i64, 16] {
        let e0_rs = json_num(&js, &format!("n64_x{}_q1_e0", x));
        let e0_py = json_num(&jpy, &format!("n64_x{}_q1_e0", x));
        check(
            &format!("[I4] x={} q=1 E₀ 交差エンジン ± 1e-6", x),
            (e0_rs - e0_py).abs() < 1e-6,
            format!(
                "rust {:.9} vs python {:.9} (Δ {:.1e})",
                e0_rs,
                e0_py,
                (e0_rs - e0_py).abs()
            ),
        );
    }
    println!();
    let sqrt_pi_inv = 1.0 / std::f64::consts::PI.sqrt();
    let mut ratios = Vec::new();
    for x in [9.0f64, 16.0] {
        let g1 = json_num(&jpy, &format!("n64_x{}_q1_gap", x as i64));
        let g2 = json_num(&js, &format!("n64_x{}_q2_gap", x as i64));
        let m1 = json_num(&js, &format!("n64_x{}_q1_meff", x as i64));
        let m2 = json_num(&js, &format!("n64_x{}_q2_meff", x as i64));
        println!(
            "    x={:2}: M/g(1) = {:.4} [gap py], M/g(2) = {:.4} [gap rs], 比 = {:.4} (相関比 {:.4})",
            x as i64,
            g1 / (2.0 * x.sqrt()),
            g2 / (2.0 * x.sqrt()),
            g2 / g1,
            m2 / m1
        );
        ratios.push(g2 / g1);
    }
    println!(
        "    [記録] M/g(q=1, x=16) = {:.4} vs 1/√π = {:.4} (偏差 {:+.1}%; 相関法 {:.4})",
        json_num(&jpy, "n64_x16_q1_gap") / 8.0,
        sqrt_pi_inv,
        (json_num(&jpy, "n64_x16_q1_gap") / 8.0 / sqrt_pi_inv - 1.0) * 100.0,
        json_num(&js, "n64_x16_q1_meff") * 4.0
    );
    check(
        "x=16, N=64 の gap 質量比 (q=2)/(q=1) = 2 ± 7% (当初登録の判定量)",
        (ratios[1] - 2.0).abs() < 0.14,
        format!("比 = {:.4} (x=9: {:.4})", ratios[1], ratios[0]),
    );

    let ok = nfail == 0;
    println!(
        "\n[判定] {}",
        if ok {
            "事前登録 (a): 重電荷窓の開通 — 開鎖 (ε なし) + 直交罰則 DMRG で M(2)/M(1) = 2 を確認"
        } else {
            "事前登録 (b): 外れ — 記録"
        }
    );
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v22.1".into())),
        ("ratio_x9".into(), Json::Num(ratios[0])),
        ("ratio_x16".into(), Json::Num(ratios[1])),
        (
            "gap_q1_x16_py".into(),
            Json::Num(json_num(&jpy, "n64_x16_q1_gap")),
        ),
        ("branch_a".into(), Json::Bool(ok)),
    ]);
    let p = write_artifact("results/v221_dmrgex.json", &j.render());
    println!("\n[artifact] {}", p);
    println!("\n総合判定: {}", if ok { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
