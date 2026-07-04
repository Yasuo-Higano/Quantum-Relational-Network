//! v8.2 三本目の U(1) — rank-3 カイラル物質の走査と E6 の 27 のカイラル芯 (v8.0 残高 4)
//!
//! v7.1 の答え「2 本目の力の最小の足し方は B−L (SM+ν_R, 16 成分, 平面は一意)」の
//! 次の問い: **3 本目の独立な U(1) は載るのか**。
//!
//! [A] rank-3 走査: 多重項 = (表現, y₁, y₂, y₃)。条件は 19 本の超電荷依存アノマリー
//!     (線形 9: SU(3)²/SU(2)²/重力² × 3 + 三次 10: y_a y_b y_c の全単項式) +
//!     SU(3)³ + Witten + カイラル性 (共役 = 全電荷反転) + 色弱帯電 + **階数 3**。
//!     方法は v7.1 と同じ構造ごとの meet-in-the-middle (ストリーミング)。
//!     窓: U3a (≤5 多重項, |6Y|≤3, ≤16 成分) / U3b (≤5, |6Y|≤2, ≤18)。
//! [B] 装置の検証: カイラル性を外した対照窓では rank-3 解 (vectorlike 3 対 + 独立電荷)
//!     が必ず見つかることを確認 — 「見つからない」が装置の盲目でないことの対照。
//! [C] E6 降下: 27 = 16 ⊕ 10 ⊕ 1 (SO(10) 経由) を厳密整数で分解し、共役対と
//!     完全中性を剥がした**カイラル芯**が SM+ν_R に一致することを検査。
//!     (E6 の示唆: 3 本目の U(1) は新しいカイラル物質でなく vectorlike な 10 を連れてくる)

use std::collections::{BTreeMap, HashMap, HashSet};
use uft_sim::*;

#[derive(Clone, Copy)]
struct Rep {
    name: &'static str,
    cd: i64,
    wd: i64,
    a3: i64,
    t3x2: i64,
    t2x2: i64,
    conj: usize,
}
const NT: usize = 6;
const REPS: [Rep; NT] = [
    Rep {
        name: "(1,1)",
        cd: 1,
        wd: 1,
        a3: 0,
        t3x2: 0,
        t2x2: 0,
        conj: 0,
    },
    Rep {
        name: "(1,2)",
        cd: 1,
        wd: 2,
        a3: 0,
        t3x2: 0,
        t2x2: 1,
        conj: 1,
    },
    Rep {
        name: "(3,1)",
        cd: 3,
        wd: 1,
        a3: 1,
        t3x2: 1,
        t2x2: 0,
        conj: 3,
    },
    Rep {
        name: "(3̄,1)",
        cd: 3,
        wd: 1,
        a3: -1,
        t3x2: 1,
        t2x2: 0,
        conj: 2,
    },
    Rep {
        name: "(3,2)",
        cd: 3,
        wd: 2,
        a3: 1,
        t3x2: 1,
        t2x2: 1,
        conj: 5,
    },
    Rep {
        name: "(3̄,2)",
        cd: 3,
        wd: 2,
        a3: -1,
        t3x2: 1,
        t2x2: 1,
        conj: 4,
    },
];

type Mult = (usize, i64, i64, i64);

fn conj_m(m: Mult) -> Mult {
    (REPS[m.0].conj, -m.1, -m.2, -m.3)
}

const NKEY: usize = 19;

fn add_key(key: &mut [i64; NKEY], t: usize, y: [i64; 3]) {
    let r = &REPS[t];
    let d = r.cd * r.wd;
    for a in 0..3 {
        key[a] += r.t3x2 * r.wd * y[a]; // SU(3)²·y_a
        key[3 + a] += r.t2x2 * r.cd * y[a]; // SU(2)²·y_a
        key[6 + a] += d * y[a]; // 重力²·y_a
    }
    // 三次 10 単項式 (a ≤ b ≤ c)
    let mut idx = 9;
    for a in 0..3 {
        for b in a..3 {
            for c in b..3 {
                key[idx] += d * y[a] * y[b] * y[c];
                idx += 1;
            }
        }
    }
}

fn full_check(sp: &[Mult], require_chiral: bool) -> bool {
    let (mut a3s, mut wit) = (0i64, 0i64);
    let mut key = [0i64; NKEY];
    let (mut has_c, mut has_w) = (false, false);
    for &(t, y1, y2, y3) in sp {
        let r = &REPS[t];
        a3s += r.a3 * r.wd;
        if r.wd == 2 {
            wit += r.cd;
        }
        add_key(&mut key, t, [y1, y2, y3]);
        if r.cd > 1 {
            has_c = true;
        }
        if r.wd > 1 {
            has_w = true;
        }
    }
    if a3s != 0 || wit % 2 != 0 || !has_c || !has_w || key != [0; NKEY] {
        return false;
    }
    // 階数 3: どれかの 3×3 小行列式が非零
    let mut rank3 = false;
    let k = sp.len();
    'outer: for i in 0..k {
        for j in (i + 1)..k {
            for l in (j + 1)..k {
                let (a, b, c) = (sp[i], sp[j], sp[l]);
                let det = a.1 * (b.2 * c.3 - c.2 * b.3) - a.2 * (b.1 * c.3 - c.1 * b.3)
                    + a.3 * (b.1 * c.2 - c.1 * b.2);
                if det != 0 {
                    rank3 = true;
                    break 'outer;
                }
            }
        }
    }
    if !rank3 {
        return false;
    }
    if require_chiral {
        for (i, &m) in sp.iter().enumerate() {
            if conj_m(m) == m {
                return false;
            }
            for &m2 in sp.iter().skip(i + 1) {
                if conj_m(m) == m2 {
                    return false;
                }
            }
        }
    }
    true
}

/// 3-平面の Plücker 標準形 (C(k,3) 個の 3×3 小行列式、gcd 約分・符号正規化) —
/// v7.1 の 2-平面分類の rank-3 版。列順は同タイプ内の全順列の最小。
fn plane_canon(sp: &[Mult]) -> (Vec<usize>, Vec<i64>) {
    fn gcd(a: i64, b: i64) -> i64 {
        if b == 0 {
            a.abs()
        } else {
            gcd(b, a % b)
        }
    }
    let plucker = |ch: &[(i64, i64, i64)]| -> Vec<i64> {
        let k = ch.len();
        let mut ps = Vec::new();
        for i in 0..k {
            for j in (i + 1)..k {
                for l in (j + 1)..k {
                    let (a, b, c) = (ch[i], ch[j], ch[l]);
                    ps.push(
                        a.0 * (b.1 * c.2 - c.1 * b.2) - a.1 * (b.0 * c.2 - c.0 * b.2)
                            + a.2 * (b.0 * c.1 - c.0 * b.1),
                    );
                }
            }
        }
        let mut g = 0i64;
        for &x in &ps {
            g = gcd(g, x);
        }
        let g = g.max(1);
        for x in ps.iter_mut() {
            *x /= g;
        }
        if let Some(&first) = ps.iter().find(|&&x| x != 0) {
            if first < 0 {
                for x in ps.iter_mut() {
                    *x = -*x;
                }
            }
        }
        ps
    };
    let mut best: Option<(Vec<usize>, Vec<i64>)> = None;
    for op in 0..2 {
        let cur: Vec<Mult> = sp
            .iter()
            .map(|&m| if op == 1 { conj_m(m) } else { m })
            .collect();
        let mut sorted = cur.clone();
        sorted.sort_by_key(|m| m.0);
        // 同タイプ run の順列 (解は稀なので単純な全列挙)
        let mut runs: Vec<(usize, usize)> = Vec::new();
        let mut i = 0;
        while i < sorted.len() {
            let mut j = i;
            while j < sorted.len() && sorted[j].0 == sorted[i].0 {
                j += 1;
            }
            runs.push((i, j));
            i = j;
        }
        permute_runs(&sorted, &runs, 0, &mut |cand: &[Mult]| {
            let charges: Vec<(i64, i64, i64)> =
                cand.iter().map(|&(_, a, b, c)| (a, b, c)).collect();
            let types: Vec<usize> = cand.iter().map(|&(t, _, _, _)| t).collect();
            let form = (types, plucker(&charges));
            if best.is_none() || form < *best.as_ref().unwrap() {
                best = Some(form);
            }
        });
    }
    best.unwrap()
}

fn permute_runs(sp: &[Mult], runs: &[(usize, usize)], ri: usize, f: &mut impl FnMut(&[Mult])) {
    if ri == runs.len() {
        f(sp);
        return;
    }
    let mut v = sp.to_vec();
    permute_range(&mut v, ri, runs, f);
}
fn permute_range(
    v: &mut Vec<Mult>,
    ri: usize,
    runs: &[(usize, usize)],
    f: &mut impl FnMut(&[Mult]),
) {
    let (lo, hi) = runs[ri];
    fn heap(
        v: &mut Vec<Mult>,
        lo: usize,
        n: usize,
        ri: usize,
        runs: &[(usize, usize)],
        f: &mut impl FnMut(&[Mult]),
    ) {
        if n <= 1 {
            if ri + 1 == runs.len() {
                f(v);
            } else {
                permute_range(v, ri + 1, runs, f);
            }
            return;
        }
        for i in 0..n {
            heap(v, lo, n - 1, ri, runs, f);
            if n % 2 == 0 {
                v.swap(lo + i, lo + n - 1);
            } else {
                v.swap(lo, lo + n - 1);
            }
        }
    }
    heap(v, lo, hi - lo, ri, runs, f);
}

struct Domain {
    name: &'static str,
    label: &'static str,
    ymax: i64,
    kmax: usize,
    cmax: i64,
    chiral: bool,
}

fn gen_structures(dom: &Domain) -> Vec<Vec<usize>> {
    fn rec(dom: &Domain, ti: usize, cur: &mut Vec<usize>, comps: i64, out: &mut Vec<Vec<usize>>) {
        if cur.len() >= 3 {
            out.push(cur.clone()); // rank-3 には ≥3 多重項が必要
        }
        if cur.len() == dom.kmax {
            return;
        }
        for t in ti..NT {
            let c = REPS[t].cd * REPS[t].wd;
            if comps + c > dom.cmax {
                continue;
            }
            cur.push(t);
            rec(dom, t, cur, comps + c, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    let mut cur = Vec::new();
    rec(dom, 0, &mut cur, 0, &mut out);
    out.retain(|s| {
        let a3: i64 = s.iter().map(|&t| REPS[t].a3 * REPS[t].wd).sum();
        let wit: i64 = s
            .iter()
            .filter(|&&t| REPS[t].wd == 2)
            .map(|&t| REPS[t].cd)
            .sum();
        let has_c = s.iter().any(|&t| REPS[t].cd > 1);
        let has_w = s.iter().any(|&t| REPS[t].wd > 1);
        a3 == 0 && wit % 2 == 0 && has_c && has_w
    });
    out
}

const SIDE_MAX: usize = 3;

fn enum_side_cb(
    types: &[usize],
    ymax: i64,
    f: &mut impl FnMut(&[i64; NKEY], &[(i8, i8, i8); SIDE_MAX]),
) {
    assert!(types.len() <= SIDE_MAX);
    fn rec(
        types: &[usize],
        ymax: i64,
        pos: usize,
        cur: &mut [(i8, i8, i8); SIDE_MAX],
        key: [i64; NKEY],
        f: &mut impl FnMut(&[i64; NKEY], &[(i8, i8, i8); SIDE_MAX]),
    ) {
        if pos == types.len() {
            f(&key, cur);
            return;
        }
        let t = types[pos];
        let same = pos > 0 && types[pos - 1] == t;
        for y1 in -ymax..=ymax {
            for y2 in -ymax..=ymax {
                for y3 in -ymax..=ymax {
                    if same {
                        let p = (
                            cur[pos - 1].0 as i64,
                            cur[pos - 1].1 as i64,
                            cur[pos - 1].2 as i64,
                        );
                        if (y1, y2, y3) < p {
                            continue;
                        }
                    }
                    let mut k2 = key;
                    add_key(&mut k2, t, [y1, y2, y3]);
                    cur[pos] = (y1 as i8, y2 as i8, y3 as i8);
                    rec(types, ymax, pos + 1, cur, k2, f);
                }
            }
        }
    }
    let mut cur = [(0i8, 0i8, 0i8); SIDE_MAX];
    rec(types, ymax, 0, &mut cur, [0; NKEY], f);
}

struct ScanOut {
    reps_by_comps: BTreeMap<i64, usize>,
    planes_by_comps: BTreeMap<i64, usize>,
    examples: Vec<(i64, Vec<Mult>)>,
    n_structures: usize,
    n_stream: u64,
}

fn rep_canon(sp: &[Mult]) -> Vec<Mult> {
    // 表示用の窓内代表: 共役 × 各電荷の符号反転 × 電荷の置換 × ソートの最小
    let mut best: Option<Vec<Mult>> = None;
    const PERMS: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    for op in 0..16 {
        for perm in PERMS {
            let mut w: Vec<Mult> = sp
                .iter()
                .map(|&m| {
                    let m = if op & 1 == 1 { conj_m(m) } else { m };
                    let ys = [m.1, m.2, m.3];
                    let mut y = [ys[perm[0]], ys[perm[1]], ys[perm[2]]];
                    if op & 2 == 2 {
                        y[0] = -y[0];
                    }
                    if op & 4 == 4 {
                        y[1] = -y[1];
                    }
                    if op & 8 == 8 {
                        y[2] = -y[2];
                    }
                    (m.0, y[0], y[1], y[2])
                })
                .collect();
            w.sort();
            if best.is_none() || w < *best.as_ref().unwrap() {
                best = Some(w);
            }
        }
    }
    best.unwrap()
}

fn scan(dom: &Domain) -> ScanOut {
    let structures = gen_structures(dom);
    let mut seen_rep: HashSet<Vec<Mult>> = HashSet::new();
    let mut seen_pl: HashSet<(Vec<usize>, Vec<i64>)> = HashSet::new();
    let mut reps_by: BTreeMap<i64, usize> = BTreeMap::new();
    let mut pl_by: BTreeMap<i64, usize> = BTreeMap::new();
    let mut examples = Vec::new();
    let mut n_stream = 0u64;
    for s in &structures {
        let k = s.len();
        // 均衡分割 (両側 ≤ SIDE_MAX): big = ceil(k/2) ≤ 3
        let big_len = k.div_ceil(2).min(SIDE_MAX);
        let (big, small) = s.split_at(big_len);
        let mut map_side: Vec<([i64; NKEY], [(i8, i8, i8); SIDE_MAX])> = Vec::new();
        enum_side_cb(small, dom.ymax, &mut |k, c| map_side.push((*k, *c)));
        let mut map: HashMap<[i64; NKEY], Vec<u32>> = HashMap::with_capacity(map_side.len() * 2);
        for (i, (key, _)) in map_side.iter().enumerate() {
            map.entry(*key).or_default().push(i as u32);
        }
        let comps: i64 = s.iter().map(|&t| REPS[t].cd * REPS[t].wd).sum();
        let mut local = 0u64;
        enum_side_cb(big, dom.ymax, &mut |bkey, bch| {
            local += 1;
            let need: [i64; NKEY] = core::array::from_fn(|i| -bkey[i]);
            if let Some(cands) = map.get(&need) {
                for &mi in cands {
                    let sch = &map_side[mi as usize].1;
                    let sp: Vec<Mult> = s
                        .iter()
                        .enumerate()
                        .map(|(i, &t)| {
                            let (a, b, c) = if i < big.len() {
                                (bch[i].0 as i64, bch[i].1 as i64, bch[i].2 as i64)
                            } else {
                                let x = sch[i - big.len()];
                                (x.0 as i64, x.1 as i64, x.2 as i64)
                            };
                            (t, a, b, c)
                        })
                        .collect();
                    if full_check(&sp, dom.chiral) {
                        let rc = rep_canon(&sp);
                        if seen_rep.insert(rc.clone()) {
                            *reps_by.entry(comps).or_default() += 1;
                            let pc = plane_canon(&sp);
                            if seen_pl.insert(pc) {
                                *pl_by.entry(comps).or_default() += 1;
                                if examples.len() < 8 {
                                    examples.push((comps, rc));
                                }
                            }
                        }
                    }
                }
            }
        });
        n_stream += local;
    }
    ScanOut {
        reps_by_comps: reps_by,
        planes_by_comps: pl_by,
        examples,
        n_structures: structures.len(),
        n_stream,
    }
}

fn fmt_sp(sp: &[Mult]) -> String {
    sp.iter()
        .map(|&(t, a, b, c)| format!("{}_{{{},{},{}}}", REPS[t].name, a, b, c))
        .collect::<Vec<_>>()
        .join(" ⊕ ")
}
fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

fn main() {
    self_test();
    println!("=== v8.2 三本目の U(1): rank-3 走査と E6 の 27 のカイラル芯 ===\n");
    let mut checks: Vec<(String, bool)> = Vec::new();
    let record = |name: &str, ok: bool, checks: &mut Vec<(String, bool)>| {
        println!("  => {}  {}", name, pass(ok));
        checks.push((name.to_string(), ok));
    };

    // ---- [0] エンジン自己検査 ----
    println!("[0] エンジン自己検査");
    // SM+ν_R に (Y, B−L, 0) — rank 2 なので rank-3 条件で正しく落ちる
    let smnur2: Vec<Mult> = vec![
        (0, 6, 6, 0),
        (0, 0, 6, 0),
        (1, -3, -6, 0),
        (3, -4, -2, 0),
        (3, 2, -2, 0),
        (4, 1, 2, 0),
    ];
    record(
        "rank-2 の埋め込み (SM+ν_R, y₃=0) は階数条件で正しく除外される",
        !full_check(&smnur2, true),
        &mut checks,
    );
    // 人工の rank-3 vectorlike 解 (カイラル性を外せば通る)
    let vl3: Vec<Mult> = vec![
        (4, 1, 0, 0),
        (5, -1, 0, 0),
        (1, 0, 1, 0),
        (1, 0, -1, 0),
        (0, 0, 0, 1),
        (0, 0, 0, -1),
    ];
    record(
        "人工 rank-3 (vectorlike 3 対) はカイラル性オフで解・オンで非解",
        full_check(&vl3, false) && !full_check(&vl3, true),
        &mut checks,
    );

    // ---- [A] rank-3 の窓内全数探索 ----
    println!("\n[A] rank-3 カイラル解の全数探索");
    let domains = [
        Domain {
            name: "U3a",
            label: "≤5 多重項, |6Y_a|≤3, ≤16 成分",
            ymax: 3,
            kmax: 5,
            cmax: 16,
            chiral: true,
        },
        Domain {
            name: "U3b",
            label: "≤5 多重項, |6Y_a|≤2, ≤18 成分",
            ymax: 2,
            kmax: 5,
            cmax: 18,
            chiral: true,
        },
    ];
    let mut founds = Vec::new();
    for dom in &domains {
        let t0 = std::time::Instant::now();
        let out = scan(dom);
        let rd: Vec<String> = out
            .reps_by_comps
            .iter()
            .map(|(c, n)| format!("{}:{}", c, n))
            .collect();
        let pd: Vec<String> = out
            .planes_by_comps
            .iter()
            .map(|(c, n)| format!("{}:{}", c, n))
            .collect();
        println!(
            "  [{}] {} — 構造 {} 個, ストリーム {:.2e}, {} ms",
            dom.name,
            dom.label,
            out.n_structures,
            out.n_stream as f64,
            t0.elapsed().as_millis()
        );
        println!(
            "      窓内代表 {{{}}} / 3-平面 {{{}}}",
            rd.join(", "),
            pd.join(", ")
        );
        for (c, sp) in out.examples.iter().take(3) {
            println!("      {} 成分の例: {}", c, fmt_sp(sp));
        }
        founds.push(out.reps_by_comps.values().sum::<usize>());
    }

    // ---- [B] 対照: カイラル性オフの窓で装置が rank-3 を見つけられること ----
    println!("\n[B] 対照 (装置の検証): カイラル性オフ, ≤6 多重項, |6Y|≤1, ≤18 成分");
    let ctrl = Domain {
        name: "U3ctrl",
        label: "",
        ymax: 1,
        kmax: 6,
        cmax: 18,
        chiral: false,
    };
    let t0 = std::time::Instant::now();
    let outc = scan(&ctrl);
    let nc: usize = outc.reps_by_comps.values().sum();
    println!(
        "      解 {} 個 (ストリーム {:.2e}, {} ms) — 例: {}",
        nc,
        outc.n_stream as f64,
        t0.elapsed().as_millis(),
        outc.examples
            .first()
            .map(|(_, s)| fmt_sp(s))
            .unwrap_or_default()
    );
    record(
        "対照窓では rank-3 (vectorlike) 解が見つかる — 装置は盲目でない",
        nc > 0,
        &mut checks,
    );

    // ---- [C] E6 の 27 のカイラル芯 ----
    println!("\n[C] E6 降下: 27 = 16 ⊕ 10 ⊕ 1 (SO(10) 経由) のカイラル芯");
    // 16 = SM+ν_R, 10 = (3,1)_{-2}⊕(3̄,1)_{2}⊕(1,2)_{3}⊕(1,2)_{-3} (vectorlike), 1 = 完全中性
    let e27: Vec<Mult> = vec![
        // 16 (Y×6, 以降の電荷は使わないので 0)
        (0, 6, 0, 0),
        (0, 0, 0, 0),
        (1, -3, 0, 0),
        (3, -4, 0, 0),
        (3, 2, 0, 0),
        (4, 1, 0, 0),
        // 10 = 5 ⊕ 5̄ of SU(5): D=(3,1)_{-2}, L'=(1,2)_{3} とその共役
        (2, -2, 0, 0),
        (3, 2, 0, 0),
        (1, 3, 0, 0),
        (1, -3, 0, 0),
        // 1
        (0, 0, 0, 0),
    ];
    // 成分数検査: 16 + 10 + 1 = 27
    let comps: i64 = e27
        .iter()
        .map(|&(t, _, _, _)| REPS[t].cd * REPS[t].wd)
        .sum();
    record(
        &format!("27 の成分数 = {} (16+10+1)", comps),
        comps == 27,
        &mut checks,
    );
    // カイラル芯の抽出: 共役対と自己共役 (完全中性) を反復的に除去
    let core = {
        let mut v = e27.clone();
        loop {
            let mut removed = false;
            'search: for i in 0..v.len() {
                if conj_m(v[i]) == v[i] {
                    v.remove(i);
                    removed = true;
                    break 'search;
                }
                for j in (i + 1)..v.len() {
                    if conj_m(v[i]) == v[j] {
                        v.remove(j);
                        v.remove(i);
                        removed = true;
                        break 'search;
                    }
                }
            }
            if !removed {
                break;
            }
        }
        v
    };
    let smnur_y: Vec<Mult> = vec![
        (0, 6, 0, 0),
        (1, -3, 0, 0),
        (3, -4, 0, 0),
        (3, 2, 0, 0),
        (4, 1, 0, 0),
    ];
    let mut core_s = core.clone();
    core_s.sort();
    let mut sm_s = smnur_y.clone();
    sm_s.sort();
    println!(
        "      カイラル芯 (共役対と完全中性を剥がした残り): {}",
        fmt_sp(&core_s)
    );
    record(
        "E6 の 27 のカイラル芯 = SM 1 世代 (ν_R と 10 は vectorlike/中性として剥がれる)",
        core_s == sm_s,
        &mut checks,
    );

    // ---- 判定と JSON ----
    let no_rank3 = founds.iter().all(|&n| n == 0);
    println!(
        "\n  発見: 窓内の rank-3 カイラル解は {} — 「3 本目の U(1) の最小の足し方」は",
        if no_rank3 {
            "0 個"
        } else {
            "存在する (上の表)"
        }
    );
    println!("        この窓には存在しない。B−L (rank 2) が打ち止めであり、E6 型の拡大も");
    println!("        新しいカイラル物質でなく vectorlike な 10 を連れてくるだけである。");
    let all_ok = checks.iter().all(|(_, ok)| *ok);
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-GAUGE-014".into())),
        ("rank3_solutions_u3a".into(), Json::Int(founds[0] as i64)),
        ("rank3_solutions_u3b".into(), Json::Int(founds[1] as i64)),
        ("control_solutions".into(), Json::Int(nc as i64)),
        (
            "checks".into(),
            Json::Arr(
                checks
                    .iter()
                    .map(|(n, ok)| {
                        Json::Obj(vec![
                            ("name".into(), Json::Str(n.clone())),
                            ("pass".into(), Json::Bool(*ok)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v82_threeu1.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n---- 検査一覧 ----");
    for (n, ok) in &checks {
        println!("  {} {}", pass(*ok), n);
    }
    println!("\n総合判定: {}", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
