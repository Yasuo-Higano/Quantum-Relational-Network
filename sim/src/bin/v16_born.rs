//! v1.6 測定問題 — デコヒーレンス (アインセレクション) とボルン則の導出
//!
//! QRN は A5 で「発展は常にユニタリー」と置いた。では「観測」で何が起きるのか。
//! [A] デコヒーレンス: 中心スピン + N 環境量子ビットの厳密動力学。
//!     重ね合わせの干渉項は r(t) = Π_k cos(2g_k t) で崩壊し、N 大で事実上不可逆。
//!     環境と可換な基底 (ポインター基底) だけが安定に生き残る = 古典性の起源。
//! [B] ボルン則 (Zurek のエンバリアンス論法): 等振幅の枝は「環境側で打ち消せる対称性」
//!     により等確率でなければならない。振幅 √(m/M) の枝を m 本の等振幅微細枝に
//!     分割すれば P = m/M = |振幅|² が **数え上げ** から出る。対称性を数値で厳密検証する。

use uft_sim::*;

fn main() {
    println!("=== v1.6 測定問題: デコヒーレンスとボルン則 ===\n");
    // ---- [A] 中心スピンモデル (厳密) ----
    println!("[A] 中心スピン S + 環境 N qubits:  H = σ_z ⊗ Σ g_k σ_x^(k)");
    println!("    重ね合わせ (α|0⟩+β|1⟩)|E⟩ の干渉項 ∝ r(t) = Π cos(2 g_k t)");
    let mut rng = Rng::new(7);
    println!("\n  |r(t)| の時間発展 (g_k は [0.5,1.5] の一様乱数):");
    print!("  t:      ");
    let times: Vec<f64> = (0..=10).map(|i| i as f64 * 0.4).collect();
    for &t in &times {
        print!("{:7.1}", t);
    }
    println!();
    for &nenv in &[4usize, 32] {
        let gs: Vec<f64> = (0..nenv).map(|_| 0.5 + rng.f64()).collect();
        print!("  N={:3}:  ", nenv);
        for &t in &times {
            let r: f64 = gs.iter().map(|&g| (2.0 * g * t).cos()).product();
            print!("{:7.3}", r.abs());
        }
        println!();
    }
    // 再帰の有無: 長時間の最大値
    for &nenv in &[4usize, 32] {
        let gs: Vec<f64> = (0..nenv).map(|_| 0.5 + rng.f64()).collect();
        let mut max_rev = 0.0f64;
        let mut t = 2.0;
        while t < 2000.0 {
            let r: f64 = gs.iter().map(|&g| (2.0 * g * t).cos()).product();
            max_rev = max_rev.max(r.abs());
            t += 0.01;
        }
        println!(
            "  N={:3}: t∈[2,2000] での再帰の最大 |r| = {:.2e} {}",
            nenv,
            max_rev,
            if nenv == 4 { "(小さい環境は情報を返す)" } else { "(大環境では事実上不可逆 = 古典性)" }
        );
    }
    // デコヒーレンス時間のスケーリング: 短時間 |r| ≈ exp(-2t² Σg²) → τ_D ∝ 1/√N
    println!("\n  デコヒーレンス時間 τ_D (|r|=1/e) のスケーリング:");
    print!("  ");
    for &nenv in &[8usize, 16, 32, 64, 128] {
        let gs: Vec<f64> = (0..nenv).map(|_| 0.5 + rng.f64()).collect();
        let mut t = 0.0;
        loop {
            t += 1e-4;
            let r: f64 = gs.iter().map(|&g| (2.0 * g * t).cos()).product();
            if r.abs() < (-1.0f64).exp() {
                break;
            }
        }
        print!("N={}: τ_D·√N={:.3}  ", nenv, t * (nenv as f64).sqrt());
    }
    println!("\n  => τ_D ∝ 1/√N: 環境が大きいほど速く「古典化」する。");
    println!("     ポインター基底 = 相互作用と可換な σ_z 基底のみが安定 (アインセレクション)。");
    println!("     v1.1 の熱化・v0.8 の Page と同じ機構: 情報は消えず環境へ拡散する。\n");

    // ---- [B] ボルン則: エンバリアンス + 数え上げ ----
    println!("[B] ボルン則の導出 (Zurek): |ψ⟩ = √(2/5)|0⟩|E₀⟩ + √(3/5)|1⟩|E₁⟩ → P(0) = ?");
    // 微細化: |0⟩ 枝を 2 本、|1⟩ 枝を 3 本の等振幅枝へ (環境⊗カウンター次元 5)
    // 状態ベクトル: sys(2) ⊗ fine(5): psi[s*5 + f]
    let mut psi = [0.0f64; 10];
    let a0 = (2.0f64 / 5.0).sqrt();
    let a1 = (3.0f64 / 5.0).sqrt();
    // |E0⟩→(|F1⟩+|F2⟩)/√2, |E1⟩→(|F3⟩+|F4⟩+|F5⟩)/√3
    psi[0 * 5 + 0] = a0 / 2.0f64.sqrt();
    psi[0 * 5 + 1] = a0 / 2.0f64.sqrt();
    psi[1 * 5 + 2] = a1 / 3.0f64.sqrt();
    psi[1 * 5 + 3] = a1 / 3.0f64.sqrt();
    psi[1 * 5 + 4] = a1 / 3.0f64.sqrt();
    let amps: Vec<f64> = (0..5)
        .map(|f| (psi[0 * 5 + f].powi(2) + psi[1 * 5 + f].powi(2)).sqrt())
        .collect();
    println!("  微細化後の 5 本の枝の振幅: {:?}", amps.iter().map(|x| format!("{:.4}", x)).collect::<Vec<_>>());
    let equal = amps.iter().all(|&x| (x - 1.0 / 5.0f64.sqrt()).abs() < 1e-12);
    println!("  => 全枝が等振幅 1/√5  {}", pass(equal));
    // エンバリアンス: 等振幅枝のスワップ (系側) は環境側のスワップで打ち消せる
    // 例: 枝1↔枝2 (どちらも s=0): 環境スワップ F1↔F2 で状態は厳密に不変
    let mut psi_swapped = psi;
    psi_swapped.swap(0 * 5 + 0, 0 * 5 + 1); // 環境側 F1↔F2
    let diff: f64 = psi.iter().zip(&psi_swapped).map(|(a, b)| (a - b).powi(2)).sum::<f64>().sqrt();
    println!("  エンバリアンス検証 (F1↔F2 スワップで不変): ‖Δψ‖ = {:.1e}  {}", diff, pass(diff < 1e-14));
    // 系の確率は環境ユニタリで変わらない (局所性): ρ_S の不変性を乱数ユニタリで検証
    let mut max_change = 0.0f64;
    for trial in 0..20 {
        // 5×5 乱数複素ユニタリ (グラム・シュミット)
        let mut rng2 = Rng::new(100 + trial);
        let mut u = vec![[CZERO; 5]; 5];
        for i in 0..5 {
            for j in 0..5 {
                u[i][j] = C64::new(rng2.gauss(), rng2.gauss());
            }
        }
        for i in 0..5 {
            for k in 0..i {
                let mut ip = CZERO;
                for j in 0..5 {
                    ip = ip + u[k][j].conj() * u[i][j];
                }
                for j in 0..5 {
                    u[i][j] = u[i][j] - ip * u[k][j];
                }
            }
            let nrm: f64 = u[i].iter().map(|z| z.norm2()).sum::<f64>().sqrt();
            for j in 0..5 {
                u[i][j] = u[i][j].scale(1.0 / nrm);
            }
        }
        // ψ' = (1 ⊗ U) ψ, ρ_S = Tr_E |ψ'⟩⟨ψ'|
        let mut psi2 = [[CZERO; 5]; 2];
        for s in 0..2 {
            for f in 0..5 {
                for f2 in 0..5 {
                    psi2[s][f] = psi2[s][f] + u[f][f2].scale(psi[s * 5 + f2]);
                }
            }
        }
        let mut rho = [[CZERO; 2]; 2];
        for s in 0..2 {
            for s2 in 0..2 {
                for f in 0..5 {
                    rho[s][s2] = rho[s][s2] + psi2[s][f] * psi2[s2][f].conj();
                }
            }
        }
        let p0_ref = 2.0 / 5.0;
        max_change = max_change.max((rho[0][0].re - p0_ref).abs());
    }
    println!("  局所性検証 (環境の任意ユニタリで ρ_S 不変): 20 個の乱数 U で |ΔP(0)| ≤ {:.1e}  {}",
        max_change, pass(max_change < 1e-12));
    println!("\n  論証: (i) 等振幅の枝は環境側で打ち消せるスワップ対称性を持つ (検証済)。");
    println!("        (ii) 系の確率は環境操作に依存しない (検証済)。");
    println!("        (i)+(ii) ⇒ 等振幅の枝は等確率 ⇒ 5 本中 2 本が「0」⇒ P(0) = 2/5 = |√(2/5)|²");
    println!("  *** ボルン則 P=|ψ|² は追加公理ではなく、もつれの対称性 + 数え上げの帰結 ***");
    println!("\n結論: 測定の見かけの非ユニタリー性はデコヒーレンス (もつれの一方向拡散) で、");
    println!("      確率則はエンバリアンスで説明できる。公理は「ユニタリー + テンソル構造」だけに減る。");
    println!("      (残る哲学的問い: 枝の「実在」の解釈。導出の前提 (局所性・微細化可能性) は明示した)");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
