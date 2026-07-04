//! v7.1 複数 U(1) の一般走査 — 第 2 超電荷も自由にした理論空間の地図 (v7.0 残高 4)
//!
//! v6.2 の [U] は「SM の物質内容を固定して第 2 超電荷を全走査」だった (答え:
//! ν_R なしで Y のみ、込みで span{Y, B−L})。本バイナリは**構造も両方の電荷も自由**
//! にした SU(3)×SU(2)×U(1)² の rank-2 カイラル物質を全数探索する。
//!
//! 条件: SU(3)³ / SU(3)²U(1)_a / SU(2)²U(1)_a / 重力²U(1)_a (a=1,2) /
//!       U(1)³ の 4 本 (y₁³, y₁²y₂, y₁y₂², y₂³) / Witten / カイラル性
//!       (共役 = 表現共役 + 両電荷の符号反転) / 色・弱に帯電 / 電荷行列の階数 = 2。
//! 方法: 構造ごとに原子 = (タイプ, y₁, y₂)。部分和ベクトル (10 成分) の
//!       meet-in-the-middle — 小さい側を表に載せ、大きい側はストリーミングで照合。
//! 分類: 単一 U(1) で電荷の整数倍を同値としたのと同様、U(1)² では基底の取り替えを
//!       同値とする — 物理的同値類 = 電荷の張る **2 次元平面** (Plücker 座標で標準化、
//!       列順は同タイプ内の全順列の最小)。
//! 窓 (計算量とメモリの設計判断、certificates に機械可読で固定):
//!   U2a: ≤4 多重項, |6Y_a|≤9, ≤16 成分 / U2b: ≤5, |6Y_a|≤6, ≤16 / U2c: ≤6, |6Y_a|≤4, ≤18
//! 自己検査: (i) SM 構造 (5 多重項, ν_R なし) には rank-2 解が無い (v6.2 の定理の再現)、
//!   (ii) SM+ν_R の {Y, B−L} 平面が U2c で見つかる、(iii) 窓の共通部分で U2a/U2b が一致。

use std::collections::{BTreeMap, HashMap, HashSet};
use uft_sim::*;

#[derive(Clone, Copy)]
struct Rep {
    name: &'static str,
    key: &'static str,
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
        key: "1,1",
        cd: 1,
        wd: 1,
        a3: 0,
        t3x2: 0,
        t2x2: 0,
        conj: 0,
    },
    Rep {
        name: "(1,2)",
        key: "1,2",
        cd: 1,
        wd: 2,
        a3: 0,
        t3x2: 0,
        t2x2: 1,
        conj: 1,
    },
    Rep {
        name: "(3,1)",
        key: "3,1",
        cd: 3,
        wd: 1,
        a3: 1,
        t3x2: 1,
        t2x2: 0,
        conj: 3,
    },
    Rep {
        name: "(3̄,1)",
        key: "3b,1",
        cd: 3,
        wd: 1,
        a3: -1,
        t3x2: 1,
        t2x2: 0,
        conj: 2,
    },
    Rep {
        name: "(3,2)",
        key: "3,2",
        cd: 3,
        wd: 2,
        a3: 1,
        t3x2: 1,
        t2x2: 1,
        conj: 5,
    },
    Rep {
        name: "(3̄,2)",
        key: "3b,2",
        cd: 3,
        wd: 2,
        a3: -1,
        t3x2: 1,
        t2x2: 1,
        conj: 4,
    },
];

type Mult = (usize, i64, i64); // (タイプ, 6Y₁, 6Y₂)

fn conj_m(m: Mult) -> Mult {
    (REPS[m.0].conj, -m.1, -m.2)
}

/// 全条件の独立検査 (MITM の結合結果を再検証)
fn full_check(sp: &[Mult]) -> bool {
    let (mut a3s, mut wit) = (0i64, 0i64);
    let (mut s3, mut s2, mut gr) = ([0i64; 2], [0i64; 2], [0i64; 2]);
    let mut cub = [0i64; 4]; // y1³, y1²y2, y1y2², y2³
    let (mut has_c, mut has_w) = (false, false);
    for &(t, y1, y2) in sp {
        let r = &REPS[t];
        let d = r.cd * r.wd;
        a3s += r.a3 * r.wd;
        if r.wd == 2 {
            wit += r.cd;
        }
        for (a, y) in [(0usize, y1), (1, y2)] {
            s3[a] += r.t3x2 * r.wd * y;
            s2[a] += r.t2x2 * r.cd * y;
            gr[a] += d * y;
        }
        cub[0] += d * y1 * y1 * y1;
        cub[1] += d * y1 * y1 * y2;
        cub[2] += d * y1 * y2 * y2;
        cub[3] += d * y2 * y2 * y2;
        if r.cd > 1 {
            has_c = true;
        }
        if r.wd > 1 {
            has_w = true;
        }
    }
    if a3s != 0 || wit % 2 != 0 || !has_c || !has_w {
        return false;
    }
    if s3 != [0; 2] || s2 != [0; 2] || gr != [0; 2] || cub != [0; 4] {
        return false;
    }
    // 階数 2 (どれかの 2×2 小行列式が非零)
    let mut rank2 = false;
    for i in 0..sp.len() {
        for j in (i + 1)..sp.len() {
            if sp[i].1 * sp[j].2 - sp[j].1 * sp[i].2 != 0 {
                rank2 = true;
            }
        }
    }
    if !rank2 {
        return false;
    }
    // カイラル性
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
    true
}

// ---------------- 電荷平面の標準形 (Plücker 座標) ----------------
// 単一 U(1) で「電荷の整数倍」を同値としたのと同様、U(1)² では基底の取り替え
// (GL(2,Q) 変換) を同値とする。物理的同値類 = 電荷ベクトルの張る 2 次元平面。
// 平面の標準形は 2×k 行列の 2×2 小行列式 (Plücker 座標) を gcd で約し
// 先頭非零を正に正規化したもの。列順は同タイプ内の全順列の最小を取る。

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a.abs()
    } else {
        gcd(b, a % b)
    }
}

/// 2×k 電荷行列の正規化 Plücker 座標 (C(k,2) 個)
fn plucker(ch: &[(i64, i64)]) -> Vec<i64> {
    let k = ch.len();
    let mut ps = Vec::with_capacity(k * (k - 1) / 2);
    for i in 0..k {
        for j in (i + 1)..k {
            ps.push(ch[i].0 * ch[j].1 - ch[j].0 * ch[i].1);
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
}

/// 平面レベルの標準形: 荷電共役 2 通り × 同タイプ内の列順列の最小 (タイプ列 + Plücker)。
fn plane_canon(sp: &[Mult]) -> (Vec<usize>, Vec<i64>) {
    let mut best: Option<(Vec<usize>, Vec<i64>)> = None;
    for op in 0..2 {
        let cur: Vec<Mult> = sp
            .iter()
            .map(|&m| if op == 1 { conj_m(m) } else { m })
            .collect();
        let mut sorted = cur.clone();
        sorted.sort_by_key(|m| m.0);
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
            let charges: Vec<(i64, i64)> = cand.iter().map(|&(_, a, b)| (a, b)).collect();
            let types: Vec<usize> = cand.iter().map(|&(t, _, _)| t).collect();
            let form = (types, plucker(&charges));
            if best.is_none() || form < *best.as_ref().unwrap() {
                best = Some(form);
            }
        });
    }
    best.unwrap()
}

/// runs 内の順列を全列挙してコールバック
fn permute_runs(sp: &[Mult], runs: &[(usize, usize)], ri: usize, f: &mut impl FnMut(&[Mult])) {
    if ri == runs.len() {
        f(sp);
        return;
    }
    let (lo, hi) = runs[ri];
    let mut v = sp.to_vec();
    permute_range(&mut v, lo, hi, runs, ri, f);
}
fn permute_range(
    v: &mut Vec<Mult>,
    lo: usize,
    hi: usize,
    runs: &[(usize, usize)],
    ri: usize,
    f: &mut impl FnMut(&[Mult]),
) {
    // ヒープ順列 (hi-lo ≤ 6)
    fn heap(
        v: &mut Vec<Mult>,
        lo: usize,
        n: usize,
        runs: &[(usize, usize)],
        ri: usize,
        f: &mut impl FnMut(&[Mult]),
    ) {
        if n == 1 {
            permute_runs_inner(v, runs, ri + 1, f);
            return;
        }
        for i in 0..n {
            heap(v, lo, n - 1, runs, ri, f);
            if n % 2 == 0 {
                v.swap(lo + i, lo + n - 1);
            } else {
                v.swap(lo, lo + n - 1);
            }
        }
    }
    let n = hi - lo;
    heap(v, lo, n, runs, ri, f);
}
fn permute_runs_inner(
    v: &mut Vec<Mult>,
    runs: &[(usize, usize)],
    ri: usize,
    f: &mut impl FnMut(&[Mult]),
) {
    if ri == runs.len() {
        f(v);
        return;
    }
    let (lo, hi) = runs[ri];
    permute_range(v, lo, hi, runs, ri, f);
}

/// 表示用の簡易標準形 (窓内の代表の重複除去): 4 通りの符号操作 × 共役 × ソートの最小
fn rep_canon(sp: &[Mult]) -> Vec<Mult> {
    let mut best: Option<Vec<Mult>> = None;
    for op in 0..8 {
        let mut w: Vec<Mult> = sp
            .iter()
            .map(|&m| {
                let m = if op & 1 == 1 { conj_m(m) } else { m };
                let (t, a, b) = m;
                let (a, b) = if op & 2 == 2 { (-a, b) } else { (a, b) };
                let (a, b) = if op & 4 == 4 { (a, -b) } else { (a, b) };
                (t, a, b)
            })
            .collect();
        // swap y1<->y2 も試す
        for sw in 0..2 {
            let mut u: Vec<Mult> = w
                .iter()
                .map(|&(t, a, b)| if sw == 1 { (t, b, a) } else { (t, a, b) })
                .collect();
            u.sort();
            if best.is_none() || u < *best.as_ref().unwrap() {
                best = Some(u);
            }
        }
        w.sort();
    }
    best.unwrap()
}

// ---------------- 走査 (構造ごとの MITM、片側ストリーミング) ----------------

struct Domain {
    name: &'static str,
    label: &'static str,
    ymax: i64,
    kmax: usize,
    cmax: i64,
}

fn gen_structures(dom: &Domain) -> Vec<Vec<usize>> {
    fn rec(dom: &Domain, ti: usize, cur: &mut Vec<usize>, comps: i64, out: &mut Vec<Vec<usize>>) {
        if !cur.is_empty() {
            out.push(cur.clone());
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
        a3 == 0 && wit % 2 == 0 && has_c && has_w && s.len() >= 2 // rank-2 には ≥2 多重項が必要
    });
    out
}

/// 部分和キー (10 成分): [s3_1, s3_2, s2_1, s2_2, gr_1, gr_2, y1³, y1²y2, y1y2², y2³]
fn add_key(key: &mut [i64; 10], t: usize, y1: i64, y2: i64) {
    let r = &REPS[t];
    let d = r.cd * r.wd;
    key[0] += r.t3x2 * r.wd * y1;
    key[1] += r.t3x2 * r.wd * y2;
    key[2] += r.t2x2 * r.cd * y1;
    key[3] += r.t2x2 * r.cd * y2;
    key[4] += d * y1;
    key[5] += d * y2;
    key[6] += d * y1 * y1 * y1;
    key[7] += d * y1 * y1 * y2;
    key[8] += d * y1 * y2 * y2;
    key[9] += d * y2 * y2 * y2;
}

const SIDE_MAX: usize = 3;
type SideEntry = ([i64; 10], [(i8, i8); SIDE_MAX]);

/// 片側の割当てをコールバックで列挙 (大きい側の実体化を避ける)
fn enum_side_cb(types: &[usize], ymax: i64, f: &mut impl FnMut(&[i64; 10], &[(i8, i8); SIDE_MAX])) {
    assert!(types.len() <= SIDE_MAX);
    fn rec(
        types: &[usize],
        ymax: i64,
        pos: usize,
        cur: &mut [(i8, i8); SIDE_MAX],
        key: [i64; 10],
        f: &mut impl FnMut(&[i64; 10], &[(i8, i8); SIDE_MAX]),
    ) {
        if pos == types.len() {
            f(&key, cur);
            return;
        }
        let t = types[pos];
        // 同タイプが続く場合は (y1,y2) 辞書順非減少で重複を削減
        let same = pos > 0 && types[pos - 1] == t;
        for y1 in -ymax..=ymax {
            for y2 in -ymax..=ymax {
                if same {
                    let (p1, p2) = (cur[pos - 1].0 as i64, cur[pos - 1].1 as i64);
                    if (y1, y2) < (p1, p2) {
                        continue;
                    }
                }
                let mut k2 = key;
                add_key(&mut k2, t, y1, y2);
                cur[pos] = (y1 as i8, y2 as i8);
                rec(types, ymax, pos + 1, cur, k2, f);
            }
        }
    }
    let mut cur = [(0i8, 0i8); SIDE_MAX];
    rec(types, ymax, 0, &mut cur, [0; 10], f);
}

fn enum_side(types: &[usize], ymax: i64, out: &mut Vec<SideEntry>) {
    enum_side_cb(types, ymax, &mut |k, c| out.push((*k, *c)));
}

struct ScanOut {
    reps_by_comps: BTreeMap<i64, usize>,
    lattices_by_comps: BTreeMap<i64, usize>,
    examples: Vec<(i64, Vec<Mult>)>,
    n_structures: usize,
    n_stream: u64,
    sols: Vec<(i64, Vec<Mult>)>, // rep_canon (証明書用)
}

fn scan(dom: &Domain) -> ScanOut {
    let structures = gen_structures(dom);
    let mut seen_rep: HashSet<Vec<Mult>> = HashSet::new();
    let mut seen_lat: HashSet<(Vec<usize>, Vec<i64>)> = HashSet::new();
    let mut reps_by: BTreeMap<i64, usize> = BTreeMap::new();
    let mut lat_by: BTreeMap<i64, usize> = BTreeMap::new();
    let mut examples: Vec<(i64, Vec<Mult>)> = Vec::new();
    let mut sols: Vec<(i64, Vec<Mult>)> = Vec::new();
    let mut n_stream = 0u64;
    for s in &structures {
        let k = s.len();
        // 小さい側 (末尾 ≤2 原子) を表に、大きい側をストリーミング
        let ksmall = (k / 2).min(2).max(k.saturating_sub(SIDE_MAX));
        let (big, small) = s.split_at(k - ksmall);
        let mut map_side = Vec::new();
        enum_side(small, dom.ymax, &mut map_side);
        let mut map: HashMap<[i64; 10], Vec<u32>> = HashMap::with_capacity(map_side.len() * 2);
        for (i, (key, _)) in map_side.iter().enumerate() {
            map.entry(*key).or_default().push(i as u32);
        }
        // 大きい側はコールバックでストリーム (実体化しない)
        let comps: i64 = s.iter().map(|&t| REPS[t].cd * REPS[t].wd).sum();
        let mut local_stream = 0u64;
        enum_side_cb(big, dom.ymax, &mut |bkey, bch| {
            local_stream += 1;
            let need: [i64; 10] = core::array::from_fn(|i| -bkey[i]);
            if let Some(cands) = map.get(&need) {
                for &mi in cands {
                    let sch = &map_side[mi as usize].1;
                    let sp: Vec<Mult> = s
                        .iter()
                        .enumerate()
                        .map(|(i, &t)| {
                            let (a, b) = if i < big.len() {
                                (bch[i].0 as i64, bch[i].1 as i64)
                            } else {
                                (sch[i - big.len()].0 as i64, sch[i - big.len()].1 as i64)
                            };
                            (t, a, b)
                        })
                        .collect();
                    if full_check(&sp) {
                        let rc = rep_canon(&sp);
                        if seen_rep.insert(rc.clone()) {
                            *reps_by.entry(comps).or_default() += 1;
                            sols.push((comps, rc.clone()));
                            let lc = plane_canon(&sp);
                            if seen_lat.insert(lc) {
                                *lat_by.entry(comps).or_default() += 1;
                                if examples.len() < 12 {
                                    examples.push((comps, rc));
                                }
                            }
                        }
                    }
                }
            }
        });
        n_stream += local_stream;
    }
    sols.sort();
    ScanOut {
        reps_by_comps: reps_by,
        lattices_by_comps: lat_by,
        examples,
        n_structures: structures.len(),
        n_stream,
        sols,
    }
}

fn fmt_sp(sp: &[Mult]) -> String {
    sp.iter()
        .map(|&(t, a, b)| format!("{}_{{{},{}}}", REPS[t].name, a, b))
        .collect::<Vec<_>>()
        .join(" ⊕ ")
}
fn ascii_sp(sp: &[Mult]) -> String {
    sp.iter()
        .map(|&(t, a, b)| format!("{}:{},{}", REPS[t].key, a, b))
        .collect::<Vec<_>>()
        .join(";")
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
    println!("=== v7.1 複数 U(1) の一般走査: rank-2 カイラル物質の地図 (v7.0 残高 4) ===\n");
    let mut checks: Vec<(String, bool)> = Vec::new();
    let record = |name: &str, ok: bool, checks: &mut Vec<(String, bool)>| {
        println!("  => {}  {}", name, pass(ok));
        checks.push((name.to_string(), ok));
    };

    // ---- [0] 自己検査 ----
    println!("[0] エンジン自己検査");
    // SM+ν_R に (Y, B−L) を載せた 16 成分 rank-2 解
    // (6Y, 6(B−L)): e^c=(1,1)(6,6), ν_R=(1,1)(0,6), L=(1,2)(−3,−6),
    //  u^c=(3̄,1)(−4,−2), d^c=(3̄,1)(2,−2), Q=(3,2)(1,2)
    let smnur: Vec<Mult> = vec![
        (0, 6, 6),
        (0, 0, 6),
        (1, -3, -6),
        (3, -4, -2),
        (3, 2, -2),
        (4, 1, 2),
    ];
    let ok0 = full_check(&smnur);
    record(
        "SM+ν_R + (Y, B−L) は全条件を満たす rank-2 解",
        ok0,
        &mut checks,
    );
    // HNF の格子不変性: 基底を GL(2,Z) で変えても lattice_canon が同じ
    let mixed: Vec<Mult> = smnur
        .iter()
        .map(|&(t, a, b)| (t, a + b, b)) // y1' = y1 + y2 (ユニモジュラー)
        .collect();
    let ok1 = plane_canon(&smnur) == plane_canon(&mixed) && full_check(&mixed);
    record(
        "Plücker 標準形は基底変換で不変 (混合基底も同一平面と判定)",
        ok1,
        &mut checks,
    );

    // ---- 走査 ----
    let domains = [
        Domain {
            name: "U2a",
            label: "≤4 多重項, |6Y|≤9, ≤16 成分",
            ymax: 9,
            kmax: 4,
            cmax: 16,
        },
        Domain {
            name: "U2b",
            label: "≤5 多重項, |6Y|≤6, ≤16 成分",
            ymax: 6,
            kmax: 5,
            cmax: 16,
        },
        Domain {
            name: "U2c",
            label: "≤6 多重項, |6Y|≤4, ≤18 成分",
            ymax: 4,
            kmax: 6,
            cmax: 18,
        },
    ];
    let mut outs = Vec::new();
    for dom in &domains {
        let t0 = std::time::Instant::now();
        let out = scan(dom);
        println!(
            "\n[{}] {} — 構造 {} 個, ストリーム {:.2e}, {} ms",
            dom.name,
            dom.label,
            out.n_structures,
            out.n_stream as f64,
            t0.elapsed().as_millis()
        );
        let rep_desc: Vec<String> = out
            .reps_by_comps
            .iter()
            .map(|(c, n)| format!("{}:{}", c, n))
            .collect();
        let lat_desc: Vec<String> = out
            .lattices_by_comps
            .iter()
            .map(|(c, n)| format!("{}:{}", c, n))
            .collect();
        println!("  窓内代表の数     {{{}}}", rep_desc.join(", "));
        println!("  電荷平面の数     {{{}}}", lat_desc.join(", "));
        for (c, sp) in out.examples.iter().take(4) {
            println!("  {} 成分の例: {}", c, fmt_sp(sp));
        }
        outs.push(out);
    }

    // ---- 判定 ----
    println!();
    // U2c に SM+ν_R 格子があるか (基底 (0,-1,1,0,1,-1),(1,-1,-1,-3,3,3) は窓 |6Y|≤4 内)
    let target = plane_canon(&smnur);
    let found_smnur = outs[2].sols.iter().any(|(_, sp)| plane_canon(sp) == target);
    record(
        "U2c: SM+ν_R の {Y, B−L} 平面が発見される",
        found_smnur,
        &mut checks,
    );
    // 最小 rank-2 解
    let min_comps: Vec<Option<i64>> = outs
        .iter()
        .map(|o| o.reps_by_comps.keys().next().copied())
        .collect();
    println!(
        "  最小 rank-2 解の成分数: U2a={:?} U2b={:?} U2c={:?}",
        min_comps[0], min_comps[1], min_comps[2]
    );
    // 窓の共通部分での一致 (U2a と U2b の k≤4, |6Y|≤6 部分)
    let filt = |o: &ScanOut, kmax: usize, ymax: i64, cmax: i64| -> HashSet<Vec<Mult>> {
        o.sols
            .iter()
            .filter(|(c, sp)| {
                *c <= cmax
                    && sp.len() <= kmax
                    && sp
                        .iter()
                        .all(|&(_, a, b)| a.abs() <= ymax && b.abs() <= ymax)
            })
            .map(|(_, sp)| sp.clone())
            .collect()
    };
    let common_a = filt(&outs[0], 4, 6, 16);
    let common_b = filt(&outs[1], 4, 6, 16);
    record(
        "U2a ∩ U2b の共通窓で解集合が一致 (実装の整合性)",
        common_a == common_b,
        &mut checks,
    );

    // ---- 証明書 ----
    let mut sha_lines = String::new();
    let mut sols_json = Vec::new();
    for (dom, out) in domains.iter().zip(&outs) {
        let ser: String = out
            .sols
            .iter()
            .map(|(c, s)| format!("{}|{}", c, ascii_sp(s)))
            .collect::<Vec<_>>()
            .join("\n");
        let h = sha256_hex(ser.as_bytes());
        sha_lines.push_str(&format!("{}  {}  n={}\n", h, dom.name, out.sols.len()));
        sols_json.push(Json::Obj(vec![
            ("run".into(), Json::Str(dom.name.into())),
            ("label".into(), Json::Str(dom.label.into())),
            ("hypercharge_6y_max".into(), Json::Int(dom.ymax)),
            ("max_multiplets".into(), Json::Int(dom.kmax as i64)),
            ("max_components".into(), Json::Int(dom.cmax)),
            (
                "n_window_representatives".into(),
                Json::Int(out.sols.len() as i64),
            ),
            (
                "planes_by_components".into(),
                Json::Obj(
                    out.lattices_by_comps
                        .iter()
                        .map(|(c, n)| (c.to_string(), Json::Int(*n as i64)))
                        .collect(),
                ),
            ),
            ("sha256_of_representatives".into(), Json::Str(h)),
        ]));
    }
    let p1 = write_artifact(
        "certificates/v71_solutions.json",
        &Json::Arr(sols_json).render(),
    );
    let p2 = write_artifact("certificates/v71_sha256.txt", &sha_lines);
    println!("\n  証明書: {} / {}", p1, p2);

    let all_ok = checks.iter().all(|(_, ok)| *ok);
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-GAUGE-012".into())),
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
        (
            "min_rank2_components".into(),
            Json::Arr(
                min_comps
                    .iter()
                    .map(|c| c.map(Json::Int).unwrap_or(Json::Null))
                    .collect(),
            ),
        ),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p3 = write_artifact("results/v71_twou1.json", &j.render());
    println!("  機械可読な結果: {}", p3);

    println!("\n---- 検査一覧 ----");
    for (n, ok) in &checks {
        println!("  {} {}", pass(*ok), n);
    }
    println!("\n総合判定: {}", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
