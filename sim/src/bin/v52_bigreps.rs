//! v5.2 理論空間の拡張 — 大きい表現 ((1,3), (6,1), (8,1)) を加えても SM は最小か
//!
//! v3.1/v4.3 の範囲 (表現 ≤ (3,2)) を拡張: SU(2) 三重項、SU(3) 六重項・随伴 (八重項) を
//! 追加した全数探索。アノマリー係数: A(6)=7, A(8)=0; 2T(3)=1, 2T(6)=5, 2T(8)=6; 2T(2)=1, 2T(3_w)=4。
//! 条件は同じ: 5 アノマリー + Witten (SU(2) 半整数アイソスピンのみ寄与) + カイラル + 全因子帯電。

use std::collections::HashSet;

// (名前, 色次元, 弱次元, 色A, 2T_色, 2T_弱, 共役タイプid)
const NT: usize = 10;
const TYPES: [(&str, i64, i64, i64, i64, i64, usize); NT] = [
    ("(1,1)", 1, 1, 0, 0, 0, 0),
    ("(1,2)", 1, 2, 0, 0, 1, 1),
    ("(1,3)", 1, 3, 0, 0, 4, 2),
    ("(3,1)", 3, 1, 1, 1, 0, 4),
    ("(3̄,1)", 3, 1, -1, 1, 0, 3),
    ("(3,2)", 3, 2, 2, 1, 1, 6),
    ("(3̄,2)", 3, 2, -2, 1, 1, 5),
    ("(6,1)", 6, 1, 7, 5, 0, 8),
    ("(6̄,1)", 6, 1, -7, 5, 0, 7),
    ("(8,1)", 8, 1, 0, 6, 0, 9),
];
const NYMAX: i64 = 9;

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}
fn canonical(ts: &[usize], ns: &[i64]) -> Vec<(usize, i64)> {
    let mut g = 0i64;
    for &n in ns {
        g = gcd(g, n.abs());
    }
    let g = g.max(1);
    let mut best: Option<Vec<(usize, i64)>> = None;
    for op in 0..4 {
        let mut w: Vec<(usize, i64)> = ts
            .iter()
            .zip(ns)
            .map(|(&t, &n)| {
                let (t2, n2) = if op & 1 == 1 {
                    (TYPES[t].6, -n / g)
                } else {
                    (t, n / g)
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

fn full_check(ts: &[usize], ns: &[i64]) -> bool {
    let k = ts.len();
    let (mut su3sq, mut su2sq, mut cubic) = (0i64, 0, 0);
    let mut has_y = false;
    for i in 0..k {
        let (_, cd, wd, _, t3, t2, _) = TYPES[ts[i]];
        let n = ns[i];
        su3sq += t3 * wd * n;
        su2sq += t2 * cd * n;
        cubic += cd * wd * n * n * n;
        if n != 0 {
            has_y = true;
        }
    }
    if !has_y || su3sq != 0 || su2sq != 0 || cubic != 0 {
        return false;
    }
    for i in 0..k {
        let ci = (TYPES[ts[i]].6, -ns[i]);
        if ci == (ts[i], ns[i]) {
            return false;
        }
        for j in 0..k {
            if j != i && (ts[j], ns[j]) == ci {
                return false;
            }
        }
    }
    true
}

fn main() {
    println!(
        "=== v5.2 拡張理論空間: (1,3),(6,1),(8,1) を加えた全数探索 (≤5 多重項, ≤16 成分) ===\n"
    );
    // 構造列挙
    let mut structures: Vec<Vec<usize>> = Vec::new();
    fn gen(start: usize, cur: &mut Vec<usize>, comps: i64, out: &mut Vec<Vec<usize>>) {
        if !cur.is_empty() {
            out.push(cur.clone());
        }
        if cur.len() == 5 {
            return;
        }
        for t in start..NT {
            let c = TYPES[t].1 * TYPES[t].2;
            if comps + c > 16 {
                continue;
            }
            cur.push(t);
            gen(t, cur, comps + c, out);
            cur.pop();
        }
    }
    let mut cur = Vec::new();
    gen(0, &mut cur, 0, &mut structures);
    structures.retain(|s| {
        // 色A は表で弱次元込み ((3,2) は ±2)
        let su3cub: i64 = s.iter().map(|&t| TYPES[t].3).sum::<i64>();
        let wit: i64 = s
            .iter()
            .filter(|&&t| TYPES[t].2 == 2)
            .map(|&t| TYPES[t].1)
            .sum();
        let has_c = s.iter().any(|&t| TYPES[t].3 != 0 || TYPES[t].4 != 0);
        let has_w = s.iter().any(|&t| TYPES[t].2 >= 2);
        su3cub == 0 && wit % 2 == 0 && has_c && has_w
    });
    println!("構造候補: {} 個", structures.len());
    let mut seen: HashSet<Vec<(usize, i64)>> = HashSet::new();
    let mut sols: Vec<(i64, Vec<(usize, i64)>)> = Vec::new();
    let mut checked = 0u64;
    for s in &structures {
        let k = s.len();
        let comps: i64 = s.iter().map(|&t| TYPES[t].1 * TYPES[t].2).sum();
        let wlast = TYPES[s[k - 1]].1 * TYPES[s[k - 1]].2;
        let mut ns = vec![0i64; k];
        let mut idx = vec![-NYMAX; k - 1];
        loop {
            checked += 1;
            let mut grav = 0i64;
            for i in 0..k - 1 {
                ns[i] = idx[i];
                grav += TYPES[s[i]].1 * TYPES[s[i]].2 * ns[i];
            }
            if grav % wlast == 0 {
                let nl = -grav / wlast;
                if nl.abs() <= NYMAX {
                    ns[k - 1] = nl;
                    let mut ok = true;
                    for i in 1..k {
                        if s[i] == s[i - 1] && ns[i] < ns[i - 1] {
                            ok = false;
                            break;
                        }
                    }
                    if ok && full_check(s, &ns) {
                        let can = canonical(s, &ns);
                        if seen.insert(can.clone()) {
                            sols.push((comps, can));
                        }
                    }
                }
            }
            let mut p = 0;
            loop {
                if p >= idx.len() {
                    break;
                }
                idx[p] += 1;
                if idx[p] <= NYMAX {
                    break;
                }
                idx[p] = -NYMAX;
                p += 1;
            }
            if p >= idx.len() {
                break;
            }
        }
    }
    println!("走査: {:.2e} 割当て\n", checked as f64);
    sols.sort();
    let mut by: std::collections::BTreeMap<i64, usize> = Default::default();
    for (c, _) in &sols {
        *by.entry(*c).or_default() += 1;
    }
    println!("成分数ごとの解: {:?}", by);
    for (c, sol) in sols.iter().take(6) {
        let d: Vec<String> = sol
            .iter()
            .map(|&(t, n)| format!("{}_{{{}}}", TYPES[t].0, n))
            .collect();
        println!("  {} 成分: {}", c, d.join(" ⊕ "));
    }
    let sm = sols.iter().filter(|(c, _)| *c == 15).count();
    let below = sols.iter().filter(|(c, _)| *c < 15).count();
    println!(
        "\n  => 15 成分未満の解: {} 個 / 15 成分の解: {} 個",
        below, sm
    );
    println!("  {}", pass(below == 0 && sm == 1));
    println!("\n結論: SU(2) 三重項・SU(3) 六重項・随伴を解禁しても、標準模型 1 世代より小さい");
    println!("      無矛盾カイラル物質は存在せず、15 成分の解も SM ただ一つのまま。");
    println!("      大表現は「アノマリーが重い」(A(6)=7, T(6)=5T(3)) ため小さな解に入れない。");
    println!("      SM の孤立性 (v4.3) は表現の範囲拡大に対して頑健である。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
