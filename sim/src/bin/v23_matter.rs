//! v2.3 物質の内容はどこまで「必然」か — カイラリティ・世代数・超電荷
//!
//! 標準模型の最大の謎「なぜこの物質内容か」に、3 つの無矛盾性・トポロジーの機構で迫る:
//! [A] フェルミオン倍化 (Nielsen–Ninomiya): 格子(離散網)上のカイラルフェルミオンは
//!     単純には作れない — カイラリティは QRN への強い制約 (逃げ道: 高次元の境界モード)
//! [B] 世代数 = トポロジー: 磁束 Q を貫いたトーラス上の最低ランダウ準位の縮退度は
//!     厳密に Q (指数定理の最低次)。「同じ粒子が 3 コピー」は隠れ空間の位相不変量でありうる
//! [C] アノマリー消去: 標準模型 1 世代の 4 つのアノマリーが奇跡的に消えること、
//!     さらに逆に「消去を要求すると超電荷が一意に決まる」ことを数値で示す

use uft_sim::*;

fn main() {
    println!("=== v2.3 物質内容の必然性: カイラリティ・世代・電荷 ===\n");

    // ---- [A] フェルミオン倍化 ----
    println!("[A] 格子上のディラック分散 ε(k) = sin k のゼロ点 (質量ゼロのモード)");
    let nk = 100000usize;
    let kf = |i: usize| -> f64 {
        -std::f64::consts::PI + 2.0 * std::f64::consts::PI * (i as f64 + 0.5) / nk as f64
    };
    let mut zeros = Vec::new();
    for i in 0..nk {
        let (k1, k2) = (kf(i), kf((i + 1) % nk));
        // ブリルアンゾーンは周期的: 端 (±π) をまたぐ交差も数える
        if (k1.sin() < 0.0) != (k2.sin() < 0.0) {
            zeros.push((k1, k1.cos().signum()));
        }
    }
    println!("  ゼロ点: {} 個", zeros.len());
    for (k, s) in &zeros {
        println!(
            "   k = {:+.4}, 群速度の符号 = {:+.0} ({})",
            k,
            s,
            if *s > 0.0 { "右巻き" } else { "左巻き" }
        );
    }
    println!("  => 必ず左右ペアで現れる (Nielsen–Ninomiya の定理): 離散格子はカイラリティを消す。");
    println!("     Wilson 項 r(1-cos k) を足すと k=π の複製は質量 2r を得て消えるが、");
    println!("     カイラル対称性が壊れる。逃げ道 = 高次元のドメイン壁/境界モード (→[B] と同根)。");
    println!("     ***観測されるニュートリノのカイラリティは、下部構造への強い制約である***\n");

    // ---- [B] 世代数 = 磁束のトポロジー ----
    println!(
        "[B] 磁束 Q を貫いた {}×{} トーラス: 最低準位の縮退度 = Q (指数定理)",
        12, 12
    );
    let n = 12usize;
    let ns = n * n;
    for q in 1..=4usize {
        let phi = 2.0 * std::f64::consts::PI * q as f64 / ns as f64;
        // ランダウゲージ: y リンク位相 e^{iφx}, x 境界リンク (N-1→0) に e^{-iφNy}
        let idx = |x: usize, y: usize| x + y * n;
        let m = 2 * ns;
        let mut a = vec![0.0; m * m];
        let mut addhop = |i: usize, j: usize, th: f64| {
            // H[j][i] += -e^{iθ} (i→j), H[i][j] += -e^{-iθ}
            let (c, s) = (th.cos(), th.sin());
            a[j + i * m] += -c;
            a[i + j * m] += -c;
            a[(j + ns) + (i + ns) * m] += -c;
            a[(i + ns) + (j + ns) * m] += -c;
            // Him[j,i] = -s, Him[i,j] = +s → embed [[Re,-Im],[Im,Re]]
            a[j + (i + ns) * m] += s;
            a[(j + ns) + i * m] += -s;
            a[i + (j + ns) * m] += -s;
            a[(i + ns) + j * m] += s;
        };
        for x in 0..n {
            for y in 0..n {
                // y ホップ (x,y)→(x,y+1): e^{iφx}
                addhop(idx(x, y), idx(x, (y + 1) % n), phi * x as f64);
                // x ホップ: 境界のみ位相
                let th = if x == n - 1 {
                    -phi * n as f64 * y as f64
                } else {
                    0.0
                };
                addhop(idx(x, y), idx((x + 1) % n, y), th);
            }
        }
        let (w, _) = jacobi_eigh(&a, m);
        let evs: Vec<f64> = (0..ns).map(|i| 0.5 * (w[2 * i] + w[2 * i + 1])).collect();
        // 最低クラスタ: 最初の 10 準位のうち最大ギャップで切る
        let mut best = (0.0, 0usize);
        for i in 0..9 {
            let g = evs[i + 1] - evs[i];
            if g > best.0 {
                best = (g, i + 1);
            }
        }
        println!(
            "  Q={}: 最低準位群 {:?}... → 縮退度 {} (ギャップ {:.3})  {}",
            q,
            &evs[0..5.min(ns)]
                .iter()
                .map(|x| format!("{:.4}", x))
                .collect::<Vec<_>>(),
            best.1,
            best.0,
            pass(best.1 == q)
        );
    }
    println!(
        "  => 縮退度 = 磁束量子数 Q。ディラック場では「カイラルゼロモードの数 = Q」(指数定理)。"
    );
    println!("     隠れ空間に位相的な捻れ Q=3 があれば、同じ場が 3 世代に見える —");
    println!(
        "     ***「なぜ 3 世代?」は「隠れ空間のトポロジカル不変量はいくつか?」に翻訳される***\n"
    );

    // ---- [C] アノマリー消去と超電荷の一意性 ----
    println!("[C] 標準模型 1 世代のアノマリー (Y: 超電荷, 左巻きワイル基底)");
    // (名前, 色多重度, 弱多重度, Y)
    let gen: [(&str, f64, f64, f64); 6] = [
        ("Q (クォーク二重項)", 3.0, 2.0, 1.0 / 6.0),
        ("u^c", 3.0, 1.0, -2.0 / 3.0),
        ("d^c", 3.0, 1.0, 1.0 / 3.0),
        ("L (レプトン二重項)", 1.0, 2.0, -1.0 / 2.0),
        ("e^c", 1.0, 1.0, 1.0),
        ("ν^c (右巻きν)", 1.0, 1.0, 0.0),
    ];
    let a3: f64 = gen.iter().filter(|f| f.1 > 1.0).map(|f| f.2 * f.3).sum(); // [SU(3)]²U(1)
    let a2: f64 = gen.iter().filter(|f| f.2 > 1.0).map(|f| f.1 * f.3).sum(); // [SU(2)]²U(1)
    let ag: f64 = gen.iter().map(|f| f.1 * f.2 * f.3).sum(); // 重力²U(1)
    let ay: f64 = gen.iter().map(|f| f.1 * f.2 * f.3 * f.3 * f.3).sum(); // U(1)³
    let ndoublets: f64 = gen.iter().filter(|f| f.2 > 1.0).map(|f| f.1).sum();
    println!(
        "  [SU(3)]²U(1): {:+.2e}  [SU(2)]²U(1): {:+.2e}  [grav]²U(1): {:+.2e}  [U(1)]³: {:+.2e}",
        a3, a2, ag, ay
    );
    println!(
        "  Witten SU(2) 大域アノマリー: 二重項の数 = {} (偶数なら無矛盾)",
        ndoublets
    );
    let ok = a3.abs() < 1e-12
        && a2.abs() < 1e-12
        && ag.abs() < 1e-12
        && ay.abs() < 1e-12
        && (ndoublets as i64) % 2 == 0;
    println!(
        "  => 4 つの三角アノマリー + 大域アノマリーが全て消える  {}",
        ok_str(ok)
    );
    // 部分集合では消えない
    let ay_q: f64 = gen[..3].iter().map(|f| f.1 * f.2 * f.3.powi(3)).sum();
    let ay_l: f64 = gen[3..].iter().map(|f| f.1 * f.2 * f.3.powi(3)).sum();
    println!(
        "  クォークだけ: [U(1)]³ = {:+.4} ≠ 0 / レプトンだけ: {:+.4} ≠ 0",
        ay_q, ay_l
    );
    println!("  => クォークとレプトンは互いを必要とする (量子論的な整合が両者を縫い合わせる)\n");

    println!("  逆問題: 多重項構造だけ与えて、消去条件から超電荷を解く (規格化 Y_e^c = 1)");
    // 線形 3 条件: 2a+b+c=0, 3a+d=0, 6a+3b+3c+2d+e=0 → a=e/6, d=-e/2, b+c=-e/3
    // 立方条件: 6a³+3b³+3c³+2d³+e³=0 → bc が決まり b,c は二次方程式の根
    let e = 1.0f64;
    let a_ = e / 6.0;
    let d = -e / 2.0;
    let s = -e / 3.0; // b + c
                      // 3(b³+c³) = -(6a³+2d³+e³);  b³+c³ = s³-3bc·s → bc = (s³ - (b³+c³))/(3s)
    let b3c3 = -(6.0 * a_.powi(3) + 2.0 * d.powi(3) + e.powi(3)) / 3.0;
    let bc = (s.powi(3) - b3c3) / (3.0 * s);
    let disc = (s * s - 4.0 * bc).sqrt();
    let (b, c) = ((s + disc) / 2.0, (s - disc) / 2.0);
    println!(
        "  解: Y_Q = {:.4} (SM: 1/6), Y_u^c = {:.4} (SM: -2/3), Y_d^c = {:.4} (SM: 1/3),",
        a_, c, b
    );
    println!("      Y_L = {:.4} (SM: -1/2), Y_e^c = 1 (規格化)", d);
    let sm_ok = (a_ - 1.0 / 6.0).abs() < 1e-12
        && (c + 2.0 / 3.0).abs() < 1e-12
        && (b - 1.0 / 3.0).abs() < 1e-12
        && (d + 0.5).abs() < 1e-12;
    println!(
        "  => アノマリー消去だけで超電荷は標準模型の値に一意に決まる  {}",
        ok_str(sm_ok)
    );
    // 帰結: 水素原子の中性
    let q_u = 0.5 + a_;
    let q_d = -0.5 + a_;
    let q_e = -0.5 + d;
    println!(
        "  帰結: Q_p + Q_e = (2Q_u+Q_d) + Q_e = {:+.1e} — 陽子と電子の電荷が正確に相殺",
        2.0 * q_u + q_d + q_e
    );
    println!(
        "  ***水素原子が中性である (10^-21 の精度で検証済) のは、アノマリー消去の帰結でありうる***"
    );
    println!("\n結論: (a) カイラリティは離散構造への制約、(b) 世代数はトポロジーの候補、");
    println!("      (c) 電荷は無矛盾性が決める。「物質の内容は恣意的でなく、整合性が強く絞る」——");
    println!("      残る大穴は群 SU(3)×SU(2)×U(1) そのものの選択理由 (v3.0 の未解決問題へ)。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
fn ok_str(ok: bool) -> &'static str {
    pass(ok)
}
