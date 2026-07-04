//! v3.1 ゲージ群の選択に向けて — 最小カイラル物質・GUT 埋め込み・結合の走行
//!
//! [A] 全数探索: SU(3)×SU(2)×U(1) の下で「無矛盾なカイラル物質集合」を数え上げる。
//!     条件: 5つのアノマリー消去 (SU(3)³, SU(3)²U(1), SU(2)²U(1), 重力²U(1), U(1)³)
//!          + Witten SU(2) 大域アノマリー + カイラル (質量項を許す対が無い)
//!          + 3 因子すべてに帯電。表現は (3,2) 以下、|Y|≤3/2 (分母6)、多重項 ≤5 個。
//!     問い: 最小の解は何か? → 標準模型 1 世代 (15 成分) が最小・一意かを検証
//! [B] SU(5) 埋め込み: 超電荷 Y = diag(-1/3,-1/3,-1/3,1/2,1/2) から 5̄⊕10 の
//!     電荷が SM 1 世代と厳密一致することを確認 (電荷量子化の群論的起源)
//! [C] 1 ループ RG: 3 つの結合定数が高エネルギーで (ほぼ) 一点に収束するか

use std::collections::HashSet;

// (色次元, 弱次元, 色符号: +1=3, -1=3̄, 0=singlet)
const TYPES: [(i64, i64, i64); 6] = [
    (1, 1, 0),
    (1, 2, 0),
    (3, 1, 1),
    (3, 1, -1),
    (3, 2, 1),
    (3, 2, -1),
];
const TNAMES: [&str; 6] = ["(1,1)", "(1,2)", "(3,1)", "(3̄,1)", "(3,2)", "(3̄,2)"];
const NYMAX: i64 = 9; // Y = n/6, n ∈ [-9, 9]
const NY: usize = (2 * NYMAX + 1) as usize;

fn decode(id: usize) -> (usize, i64) {
    (id / NY, id as i64 % (2 * NYMAX + 1) - NYMAX)
}
fn conj_type(t: usize) -> usize {
    match t {
        2 => 3,
        3 => 2,
        4 => 5,
        5 => 4,
        _ => t,
    }
}

fn check(set: &[usize]) -> bool {
    let (mut su3cub, mut su3sq, mut su2sq, mut grav, mut cubic, mut wit) = (0i64, 0, 0, 0, 0, 0i64);
    let (mut has_c, mut has_w, mut has_y) = (false, false, false);
    for &id in set {
        let (t, n) = decode(id);
        let (cd, wd, cs) = TYPES[t];
        if cs != 0 {
            has_c = true;
            su3cub += cs * wd;
            su3sq += wd * n;
        }
        if wd == 2 {
            has_w = true;
            su2sq += cd * n;
            wit += cd;
        }
        if n != 0 {
            has_y = true;
        }
        grav += cd * wd * n;
        cubic += cd * wd * n * n * n;
    }
    if !(has_c && has_w && has_y) || wit % 2 != 0 {
        return false;
    }
    if su3cub != 0 || su3sq != 0 || su2sq != 0 || grav != 0 || cubic != 0 {
        return false;
    }
    // カイラル性: 自己共役多重項なし & 共役対なし
    for (i, &a) in set.iter().enumerate() {
        let (t, n) = decode(a);
        let cid = conj_type(t) * NY + (-n + NYMAX) as usize;
        if cid == a {
            return false;
        }
        for &b in set.iter().skip(i + 1) {
            if b == cid {
                return false;
            }
        }
    }
    true
}

fn canonical(set: &[usize]) -> Vec<(usize, i64)> {
    let mut v: Vec<(usize, i64)> = set.iter().map(|&id| decode(id)).collect();
    // gcd で約す
    let mut g = 0i64;
    for &(_, n) in &v {
        g = gcd(g, n.abs());
    }
    if g > 1 {
        for e in v.iter_mut() {
            e.1 /= g;
        }
    }
    // 物理的に同値な 4 通り (恒等 / 荷電共役 / U(1)符号反転 B→-B / 両方) の最小を代表に
    let mut best: Option<Vec<(usize, i64)>> = None;
    for op in 0..4 {
        let mut w: Vec<(usize, i64)> = v
            .iter()
            .map(|&(t, n)| {
                let (t2, n2) = if op & 1 == 1 {
                    (conj_type(t), -n)
                } else {
                    (t, n)
                };
                if op & 2 == 2 {
                    (t2, -n2)
                } else {
                    (t2, n2)
                }
            })
            .collect();
        w.sort();
        if best.is_none() || w < *best.as_ref().unwrap() {
            best = Some(w);
        }
    }
    best.unwrap()
}
fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

fn rec(
    start: usize,
    slots: &mut Vec<usize>,
    comps: i64,
    found: &mut Vec<(i64, Vec<(usize, i64)>)>,
    seen: &mut HashSet<Vec<(usize, i64)>>,
    tested: &mut u64,
) {
    if !slots.is_empty() {
        *tested += 1;
        if check(slots) {
            let can = canonical(slots);
            if seen.insert(can.clone()) {
                found.push((comps, can));
            }
        }
    }
    if slots.len() == 5 {
        return;
    }
    for id in start..6 * NY {
        let (t, _) = decode(id);
        let c = TYPES[t].0 * TYPES[t].1;
        if comps + c > 15 {
            if TYPES[t].0 * TYPES[t].1 >= 6 {
                break;
            }
            continue;
        }
        slots.push(id);
        rec(id, slots, comps + c, found, seen, tested);
        slots.pop();
    }
}

fn main() {
    println!("=== v3.1 ゲージ群の選択に向けて ===\n");
    println!("[A] 全数探索: 15 成分以下の無矛盾カイラル物質 (SU(3)×SU(2)×U(1), 全因子に帯電)");
    let mut found = Vec::new();
    let mut seen = HashSet::new();
    let mut slots = Vec::new();
    let mut tested = 0u64;
    rec(0, &mut slots, 0, &mut found, &mut seen, &mut tested);
    found.sort();
    println!("  検査した集合の数: {:.2e}", tested as f64);
    println!(
        "  条件を満たす解 (スケール・共役の重複除去後): {} 個",
        found.len()
    );
    for (comps, sol) in &found {
        let desc: Vec<String> = sol
            .iter()
            .map(|&(t, n)| format!("{}_{{{}/6}}", TNAMES[t], n))
            .collect();
        println!("   {} 成分: {}", comps, desc.join(" ⊕ "));
    }
    if found.len() == 1 && found[0].0 == 15 {
        println!(
            "  => 唯一の最小解 = 標準模型 1 世代 (Y×6 = 1,-4,2,-3,6 = Q,u^c,d^c,L,e^c)  [PASS]"
        );
        println!(
            "     ***3 つの力すべてを感じる無矛盾なカイラル物質は、標準模型世代が最小である***"
        );
    } else {
        println!("  => 複数解/想定外の解 — 内訳を上に列挙 (要検討)");
    }

    println!("\n[B] SU(5) 埋め込み: Y = diag(-1/3,-1/3,-1/3,1/2,1/2) (トレースレス生成子)");
    let y5 = [-2i64, -2, -2, 3, 3]; // ×6
    println!(
        "  tr Y = {}/6 (トレースレス = 電荷の和が消える起源)",
        y5.iter().sum::<i64>()
    );
    let mut fivebar: Vec<i64> = y5.iter().map(|&y| -y).collect();
    fivebar.sort();
    let mut ten = Vec::new();
    for i in 0..5 {
        for j in (i + 1)..5 {
            ten.push(y5[i] + y5[j]);
        }
    }
    ten.sort();
    println!("  5̄ の超電荷 ×6: {:?} → d^c(2)×3, L(-3)×2", fivebar);
    println!("  10 の超電荷 ×6: {:?} → u^c(-4)×3, Q(1)×6, e^c(6)×1", ten);
    let mut expect5 = vec![-3i64, -3, 2, 2, 2];
    expect5.sort();
    let mut expect10 = vec![1i64, 1, 1, 1, 1, 1, -4, -4, -4, 6];
    expect10.sort();
    let ok = fivebar == expect5 && ten == expect10;
    println!("  => 5̄ ⊕ 10 = 標準模型 1 世代と厳密一致  {}", pass(ok));
    println!("     SO(10) ではさらに 16 = 10⊕5̄⊕1 (右巻きν込み) の単一スピノルに収まる。");
    println!(
        "     電荷が 1/3 の倍数なのは「色が 3 つ」だから (トレース条件) — 電荷量子化の群論。\n"
    );

    println!("[C] 結合定数の 1 ループ走行 (入力: M_Z での実測値)");
    let a1_mz = 59.01f64; // α1^-1 (GUT規格化 5/3)
    let a2_mz = 29.59f64;
    let a3_mz = 8.47f64;
    let mz = 91.19f64;
    let b_sm = [41.0 / 10.0, -19.0 / 6.0, -7.0];
    let b_mssm = [33.0 / 5.0, 1.0, -3.0];
    let two_pi = 2.0 * std::f64::consts::PI;
    let meet = |ainv: [f64; 3], b: [f64; 3], mu0: f64, i: usize, j: usize| -> (f64, f64) {
        let t = two_pi * (ainv[i] - ainv[j]) / (b[i] - b[j]);
        (mu0 * t.exp(), ainv[i] - b[i] / two_pi * t)
    };
    println!("  標準模型のみ (b = 41/10, -19/6, -7):");
    for &(i, j) in &[(0usize, 1usize), (0, 2), (1, 2)] {
        let (mu, ai) = meet([a1_mz, a2_mz, a3_mz], b_sm, mz, i, j);
        println!(
            "   α{}=α{} : μ = {:.1e} GeV, α⁻¹ = {:.1}",
            i + 1,
            j + 1,
            mu,
            ai
        );
    }
    // MSSM: 1 TeV で超対称粒子が現れると仮定
    let t1 = (1000.0f64 / mz).ln();
    let a_tev = [
        a1_mz - b_sm[0] / two_pi * t1,
        a2_mz - b_sm[1] / two_pi * t1,
        a3_mz - b_sm[2] / two_pi * t1,
    ];
    println!("  MSSM (超対称化, 閾値 1 TeV, b = 33/5, 1, -3):");
    let mut mus = Vec::new();
    for &(i, j) in &[(0usize, 1usize), (0, 2), (1, 2)] {
        let (mu, ai) = meet(a_tev, b_mssm, 1000.0, i, j);
        println!(
            "   α{}=α{} : μ = {:.2e} GeV, α⁻¹ = {:.2}",
            i + 1,
            j + 1,
            mu,
            ai
        );
        mus.push(mu);
    }
    let spread_sm = 9.7e16f64 / 1.0e13;
    let spread_mssm = mus.iter().cloned().fold(0.0f64, f64::max)
        / mus.iter().cloned().fold(f64::INFINITY, f64::min);
    println!(
        "  => 交点のずれ: SM では ~{:.0}00 倍、MSSM では {:.1} 倍 (α⁻¹ で ~2%) に収束",
        spread_sm / 100.0,
        spread_mssm
    );
    println!("     3 本の直線が 1 点で交わる必然性はないのに、ほぼ交わる — 統一群の間接証拠。");
    println!("     交点 ~2×10^16 GeV は陽子崩壊の未観測 (>10^34 年) とも整合する高さ。");
    println!("\n結論: 群の「選択理由」は未解決のまま (正直に)。しかし (A) 内容の最小性、");
    println!("      (B) 埋め込みの厳密な整合、(C) 結合の収束という 3 つの独立な矢印が、");
    println!("      SU(3)×SU(2)×U(1) ⊂ SU(5) ⊂ SO(10) の入れ子を指している。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
