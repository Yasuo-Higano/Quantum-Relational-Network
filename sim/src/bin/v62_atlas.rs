//! v6.2 理論空間探索の要塞化 — 独立実装・拡張軸・証明書・陰性対照
//!
//! v3.1/v4.3/v5.2 の「制限領域内で SM 1 世代が最小・一意・孤立」(claims.yml: C2) を
//! 査読に耐える形へ格上げする:
//!   [R] 独立実装によるクロスチェック — 全アノマリー条件は多重項ごとの整数寄与
//!       ベクトルの和が零になる条件なので、構造 (表現の多重集合) ごとに超電荷割当てを
//!       半分に割り部分和を突き合わせる meet-in-the-middle で解く。v3.1 (深さ優先) /
//!       v4.3 (線形消去+グリッド) と独立な第三のアルゴリズム。さらに v3.1 域は本
//!       バイナリ内蔵の直接 DFS とも突き合わせる (同一バイナリ内の二重実装)。
//!   [E] 反例探索 — 表現・|Y|・多重項数・ν_R の各軸で領域を拡張し、SM より小さい/
//!       同サイズの新しい解が現れないかを探す。
//!   [N] 陰性対照 — 条件 (カイラル性・全因子帯電・Witten・各アノマリー) を 1 つずつ
//!       外し、「SM の孤立がどの仮定に依存するか」の地図を作る。基線の解集合が
//!       対照の解集合に必ず含まれること (単調性) も検査する。
//!   [U] 第 2 の U(1) — SM の物質内容に独立な第 2 超電荷が載るかを全走査で判定。
//!   [G] GUT 降下 — SU(5)/SO(10)/Pati–Salam からの分岐が SM(+ν_R) と厳密一致するか。
//!   証明書: 探索領域の機械可読な定義・全解の canonical 表現・SHA-256 を
//!   certificates/ に出力する。
//!
//! 群論係数 (Slansky): SU(3) 三次係数 A(3)=1, A(6)=7, A(8)=0 (共役は符号反転)。
//! 指数 2T: SU(3) で 3→1, 6→5, 8→6 / SU(2) で 2→1, 3→4。Witten 大域アノマリーに
//! 寄与するのは半整数アイソスピン (本域では二重項のみ; T が半奇整数の表現)。
//! 超電荷は 6Y の整数で扱い、全演算は i64 の厳密整数 (浮動小数は使わない)。

use std::collections::{BTreeMap, HashMap, HashSet};
use uft_sim::*;

// (表示名, ハッシュ用 ASCII 名, 色次元, 弱次元, 色三次 A3, 2T_色, 2T_弱, 共役 id)
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

const NT: usize = 10;
const REPS: [Rep; NT] = [
    Rep { name: "(1,1)", key: "1,1", cd: 1, wd: 1, a3: 0, t3x2: 0, t2x2: 0, conj: 0 },
    Rep { name: "(1,2)", key: "1,2", cd: 1, wd: 2, a3: 0, t3x2: 0, t2x2: 1, conj: 1 },
    Rep { name: "(1,3)", key: "1,3", cd: 1, wd: 3, a3: 0, t3x2: 0, t2x2: 4, conj: 2 },
    Rep { name: "(3,1)", key: "3,1", cd: 3, wd: 1, a3: 1, t3x2: 1, t2x2: 0, conj: 4 },
    Rep { name: "(3̄,1)", key: "3b,1", cd: 3, wd: 1, a3: -1, t3x2: 1, t2x2: 0, conj: 3 },
    Rep { name: "(3,2)", key: "3,2", cd: 3, wd: 2, a3: 1, t3x2: 1, t2x2: 1, conj: 6 },
    Rep { name: "(3̄,2)", key: "3b,2", cd: 3, wd: 2, a3: -1, t3x2: 1, t2x2: 1, conj: 5 },
    Rep { name: "(6,1)", key: "6,1", cd: 6, wd: 1, a3: 7, t3x2: 5, t2x2: 0, conj: 8 },
    Rep { name: "(6̄,1)", key: "6b,1", cd: 6, wd: 1, a3: -7, t3x2: 5, t2x2: 0, conj: 7 },
    Rep { name: "(8,1)", key: "8,1", cd: 8, wd: 1, a3: 0, t3x2: 6, t2x2: 0, conj: 9 },
];

/// 条件トグル (陰性対照で 1 つずつ外す)
#[derive(Clone)]
struct Cond {
    su3cub: bool,
    su3sq: bool,
    su2sq: bool,
    grav: bool,
    cubic: bool,
    witten: bool,
    chiral: bool,
    all_factors: bool,
    sterile_ok: bool, // 完全中性 (1,1)_0 をカイラル条件から除外 (ν_R)
}
impl Cond {
    fn all() -> Self {
        Cond { su3cub: true, su3sq: true, su2sq: true, grav: true, cubic: true, witten: true, chiral: true, all_factors: true, sterile_ok: false }
    }
}

#[derive(Clone)]
struct Domain {
    name: &'static str,
    label: &'static str,
    types: Vec<usize>,
    nymax: i64,
    kmax: usize,
    cmax: i64,
    cond: Cond,
}

type Mult = (usize, i64);

fn comps_of(sp: &[Mult]) -> i64 {
    sp.iter().map(|&(t, _)| REPS[t].cd * REPS[t].wd).sum()
}
fn conj_m(m: Mult) -> Mult {
    (REPS[m.0].conj, -m.1)
}
fn is_sterile(m: Mult) -> bool {
    REPS[m.0].cd == 1 && REPS[m.0].wd == 1 && m.1 == 0
}

/// 全条件の独立検査 (MITM の結合結果もこれで再検証する)
fn full_check(sp: &[Mult], cond: &Cond) -> bool {
    let (mut a3s, mut s3, mut s2, mut gr, mut cu, mut wit) = (0i64, 0i64, 0i64, 0i64, 0i64, 0i64);
    let (mut has_c, mut has_w, mut has_y) = (false, false, false);
    for &(t, n) in sp {
        let r = &REPS[t];
        a3s += r.a3 * r.wd;
        s3 += r.t3x2 * r.wd * n;
        s2 += r.t2x2 * r.cd * n;
        gr += r.cd * r.wd * n;
        cu += r.cd * r.wd * n * n * n;
        if r.wd == 2 {
            wit += r.cd;
        }
        if r.cd > 1 {
            has_c = true;
        }
        if r.wd > 1 {
            has_w = true;
        }
        if n != 0 {
            has_y = true;
        }
    }
    if cond.su3cub && a3s != 0 {
        return false;
    }
    if cond.su3sq && s3 != 0 {
        return false;
    }
    if cond.su2sq && s2 != 0 {
        return false;
    }
    if cond.grav && gr != 0 {
        return false;
    }
    if cond.cubic && cu != 0 {
        return false;
    }
    if cond.witten && wit % 2 != 0 {
        return false;
    }
    if cond.all_factors && !(has_c && has_w && has_y) {
        return false;
    }
    if cond.chiral {
        for (i, &m) in sp.iter().enumerate() {
            if cond.sterile_ok && is_sterile(m) {
                continue;
            }
            if conj_m(m) == m {
                return false; // 自己共役 (実表現・Y=0) は質量項を許す
            }
            for &m2 in sp.iter().skip(i + 1) {
                if cond.sterile_ok && is_sterile(m2) {
                    continue;
                }
                if conj_m(m) == m2 {
                    return false; // 共役対は質量項を許す
                }
            }
        }
    }
    true
}

fn gcd(a: i64, b: i64) -> i64 {
    if b == 0 {
        a
    } else {
        gcd(b, a % b)
    }
}

/// 正準形: 超電荷の gcd 約分 + 物理的同値 4 通り (恒等/荷電共役/U(1)反転/両方) の最小 + ソート
fn canonical(sp: &[Mult]) -> Vec<Mult> {
    let mut g = 0i64;
    for &(_, n) in sp {
        g = gcd(g, n.abs());
    }
    let g = g.max(1);
    let mut best: Option<Vec<Mult>> = None;
    for op in 0..4 {
        let mut w: Vec<Mult> = sp
            .iter()
            .map(|&(t, n)| {
                let (t2, n2) = if op & 1 == 1 { (REPS[t].conj, -n / g) } else { (t, n / g) };
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

fn fmt_sp(sp: &[Mult]) -> String {
    sp.iter()
        .map(|&(t, n)| format!("{}_{{{}}}", REPS[t].name, n))
        .collect::<Vec<_>>()
        .join(" ⊕ ")
}
fn ascii_sp(sp: &[Mult]) -> String {
    sp.iter()
        .map(|&(t, n)| format!("{}:{}", REPS[t].key, n))
        .collect::<Vec<_>>()
        .join(";")
}

// ---------------- 構造 (表現の多重集合) の列挙 ----------------
fn gen_structures(dom: &Domain) -> Vec<Vec<usize>> {
    fn rec(dom: &Domain, ti: usize, cur: &mut Vec<usize>, comps: i64, out: &mut Vec<Vec<usize>>) {
        if !cur.is_empty() {
            out.push(cur.clone());
        }
        if cur.len() == dom.kmax {
            return;
        }
        for j in ti..dom.types.len() {
            let t = dom.types[j];
            let c = REPS[t].cd * REPS[t].wd;
            if comps + c > dom.cmax {
                continue;
            }
            cur.push(t);
            rec(dom, j, cur, comps + c, out);
            cur.pop();
        }
    }
    let mut out = Vec::new();
    let mut cur = Vec::new();
    rec(dom, 0, &mut cur, 0, &mut out);
    // 構造レベルで決まる条件 (超電荷に依らない) で先に絞る
    out.retain(|s| {
        let a3: i64 = s.iter().map(|&t| REPS[t].a3 * REPS[t].wd).sum();
        let wit: i64 = s.iter().filter(|&&t| REPS[t].wd == 2).map(|&t| REPS[t].cd).sum();
        let has_c = s.iter().any(|&t| REPS[t].cd > 1);
        let has_w = s.iter().any(|&t| REPS[t].wd > 1);
        (!dom.cond.su3cub || a3 == 0)
            && (!dom.cond.witten || wit % 2 == 0)
            && (!dom.cond.all_factors || (has_c && has_w))
    });
    out
}

// ---------------- meet-in-the-middle 走査 ----------------
// 片側の超電荷割当て (同一タイプの並びは非減少に制限して重複を削減)。
// key は「有効な超電荷依存条件」の部分和 [su3sq, su2sq, grav, cubic] (無効な条件は 0 のまま)。
const SIDE_MAX: usize = 4;

fn enum_side(dom: &Domain, types: &[usize], out: &mut Vec<([i64; 4], [i8; SIDE_MAX])>) {
    assert!(types.len() <= SIDE_MAX, "side が想定より長い");
    fn rec(
        dom: &Domain,
        types: &[usize],
        pos: usize,
        cur: &mut [i8; SIDE_MAX],
        key: [i64; 4],
        out: &mut Vec<([i64; 4], [i8; SIDE_MAX])>,
    ) {
        if pos == types.len() {
            out.push((key, *cur));
            return;
        }
        let t = types[pos];
        let r = &REPS[t];
        let start = if pos > 0 && types[pos - 1] == t {
            cur[pos - 1] as i64
        } else {
            -dom.nymax
        };
        for n in start..=dom.nymax {
            let mut k2 = key;
            if dom.cond.su3sq {
                k2[0] += r.t3x2 * r.wd * n;
            }
            if dom.cond.su2sq {
                k2[1] += r.t2x2 * r.cd * n;
            }
            if dom.cond.grav {
                k2[2] += r.cd * r.wd * n;
            }
            if dom.cond.cubic {
                k2[3] += r.cd * r.wd * n * n * n;
            }
            cur[pos] = n as i8;
            rec(dom, types, pos + 1, cur, k2, out);
        }
    }
    let mut cur = [0i8; SIDE_MAX];
    rec(dom, types, 0, &mut cur, [0; 4], out);
}

struct ScanOut {
    sols: Vec<(i64, Vec<Mult>)>, // (成分数, 正準形) 昇順
    n_structures: usize,
    n_side: u64,
    ms: u128,
}

fn scan(dom: &Domain) -> ScanOut {
    let t0 = std::time::Instant::now();
    let structures = gen_structures(dom);
    let mut seen: HashSet<Vec<Mult>> = HashSet::new();
    let mut sols: Vec<(i64, Vec<Mult>)> = Vec::new();
    let mut n_side = 0u64;
    let mut left: Vec<([i64; 4], [i8; SIDE_MAX])> = Vec::new();
    let mut right: Vec<([i64; 4], [i8; SIDE_MAX])> = Vec::new();
    for s in &structures {
        let k = s.len();
        let k2 = k / 2;
        left.clear();
        right.clear();
        enum_side(dom, &s[..k2], &mut left);
        enum_side(dom, &s[k2..], &mut right);
        n_side += (left.len() + right.len()) as u64;
        let mut map: HashMap<[i64; 4], Vec<u32>> = HashMap::with_capacity(right.len() * 2);
        for (i, (key, _)) in right.iter().enumerate() {
            map.entry(*key).or_default().push(i as u32);
        }
        for (lkey, lns) in &left {
            let need = [-lkey[0], -lkey[1], -lkey[2], -lkey[3]];
            if let Some(cands) = map.get(&need) {
                for &ri in cands {
                    let rns = &right[ri as usize].1;
                    let sp: Vec<Mult> = s
                        .iter()
                        .enumerate()
                        .map(|(i, &t)| (t, if i < k2 { lns[i] as i64 } else { rns[i - k2] as i64 }))
                        .collect();
                    // 結合結果を独立な full_check で再検証 (エンジンの二重化)
                    if full_check(&sp, &dom.cond) {
                        let can = canonical(&sp);
                        if seen.insert(can.clone()) {
                            sols.push((comps_of(&sp), can));
                        }
                    }
                }
            }
        }
    }
    sols.sort();
    ScanOut { sols, n_structures: structures.len(), n_side, ms: t0.elapsed().as_millis() }
}

/// 独立実装 2: v3.1 型の直接 DFS (小さい領域専用のクロスチェック)
fn dfs_scan(dom: &Domain) -> HashSet<Vec<Mult>> {
    let atoms: Vec<Mult> = dom
        .types
        .iter()
        .flat_map(|&t| (-dom.nymax..=dom.nymax).map(move |n| (t, n)))
        .collect();
    fn rec(
        dom: &Domain,
        atoms: &[Mult],
        start: usize,
        cur: &mut Vec<Mult>,
        comps: i64,
        seen: &mut HashSet<Vec<Mult>>,
    ) {
        if !cur.is_empty() && full_check(cur, &dom.cond) {
            seen.insert(canonical(cur));
        }
        if cur.len() == dom.kmax {
            return;
        }
        for i in start..atoms.len() {
            let c = REPS[atoms[i].0].cd * REPS[atoms[i].0].wd;
            if comps + c > dom.cmax {
                continue;
            }
            cur.push(atoms[i]);
            rec(dom, atoms, i, cur, comps + c, seen);
            cur.pop();
        }
    }
    let mut seen = HashSet::new();
    let mut cur = Vec::new();
    rec(dom, &atoms, 0, &mut cur, 0, &mut seen);
    seen
}

fn spectrum(sols: &[(i64, Vec<Mult>)]) -> BTreeMap<i64, usize> {
    let mut by: BTreeMap<i64, usize> = BTreeMap::new();
    for (c, _) in sols {
        *by.entry(*c).or_default() += 1;
    }
    by
}

fn print_scan(label: &str, out: &ScanOut, max_examples: usize) {
    println!(
        "  [{}] 構造 {} 個, 片側割当て {:.2e}, {} ms → 解 {} 個",
        label,
        out.n_structures,
        out.n_side as f64,
        out.ms,
        out.sols.len()
    );
    let by = spectrum(&out.sols);
    let desc: Vec<String> = by.iter().map(|(c, n)| format!("{}:{}", c, n)).collect();
    println!("     成分数スペクトル {{{}}}", desc.join(", "));
    let mut shown = 0;
    for (c, sol) in &out.sols {
        if shown >= max_examples {
            break;
        }
        println!("     {} 成分: {}", c, fmt_sp(sol));
        shown += 1;
    }
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

// SM 1 世代 (物理標準の割当て): e^c, L, u^c, d^c, Q
const SM: [Mult; 5] = [(0, 6), (1, -3), (4, -4), (4, 2), (5, 1)];
const SM_NUR: [Mult; 6] = [(0, 6), (0, 0), (1, -3), (4, -4), (4, 2), (5, 1)];

fn main() {
    self_test();
    println!("=== v6.2 理論空間探索の要塞化: 独立実装・拡張軸・証明書・陰性対照 ===\n");
    let mut checks: Vec<(String, bool)> = Vec::new();
    let record = |name: &str, ok: bool, checks: &mut Vec<(String, bool)>| {
        println!("  => {}  {}", name, pass(ok));
        checks.push((name.to_string(), ok));
    };

    // ---------------- エンジン自己検査 ----------------
    println!("[0] エンジン自己検査");
    {
        let c = Cond::all();
        let ok1 = full_check(&SM, &c);
        // 破壊した SM (e^c の Y を 6→5): アノマリーが立つはず
        let mut bad = SM;
        bad[0].1 = 5;
        let ok2 = !full_check(&bad, &c);
        record("SM は全条件を満たし、Y を 1 目盛ずらすと落ちる", ok1 && ok2, &mut checks);
        // canonical の不変性: 並べ替え・共役・反転・整数倍で同じ正準形
        let mut perm = vec![SM[3], SM[0], SM[4], SM[1], SM[2]];
        let can0 = canonical(&SM);
        let ok3 = canonical(&perm) == can0;
        for e in perm.iter_mut() {
            *e = conj_m(*e); // 荷電共役
        }
        let ok4 = canonical(&perm) == can0;
        let scaled: Vec<Mult> = SM.iter().map(|&(t, n)| (t, -2 * n)).collect(); // U(1) 反転 + 2 倍
        let ok5 = canonical(&scaled) == can0;
        record("正準形は並べ替え・共役・U(1)反転・整数倍で不変", ok3 && ok4 && ok5, &mut checks);
    }

    // ---------------- [R] 既存結果の独立再現 ----------------
    println!("\n[R] 独立実装 (MITM) による既存結果のクロスチェック");
    let small: Vec<usize> = vec![0, 1, 3, 4, 5, 6];
    let all10: Vec<usize> = (0..NT).collect();
    let mut cert_runs: Vec<(Domain, ScanOut)> = Vec::new();

    let r1 = Domain { name: "R1_v31", label: "v3.1 域: 小表現, |Y|≤3/2, ≤5 多重項, ≤15 成分", types: small.clone(), nymax: 9, kmax: 5, cmax: 15, cond: Cond::all() };
    let r1out = scan(&r1);
    print_scan(r1.label, &r1out, 2);
    let sm_can = canonical(&SM);
    let ok_r1 = r1out.sols.len() == 1 && r1out.sols[0].0 == 15 && r1out.sols[0].1 == sm_can;
    record("R1: 解は SM 1 世代ただ一つ (v3.1 の再現)", ok_r1, &mut checks);
    let dfs = dfs_scan(&r1);
    let mitm: HashSet<Vec<Mult>> = r1out.sols.iter().map(|(_, s)| s.clone()).collect();
    record("R1: 直接 DFS と MITM の解集合が一致 (二重実装)", dfs == mitm, &mut checks);
    cert_runs.push((r1, r1out));

    let r2 = Domain { name: "R2_v43", label: "v4.3 域: 小表現, |Y|≤3/2, ≤6 多重項, ≤24 成分", types: small.clone(), nymax: 9, kmax: 6, cmax: 24, cond: Cond::all() };
    let r2out = scan(&r2);
    print_scan(r2.label, &r2out, 1);
    let by2 = spectrum(&r2out.sols);
    let ok_r2 = by2 == BTreeMap::from([(15i64, 1usize), (16, 8), (24, 459)]);
    record("R2: スペクトル {15:1, 16:8, 24:459} (v4.3 の再現)", ok_r2, &mut checks);
    // v4.3 が公表した 16 成分解 4 例が解集合に含まれるか
    let published16: [&[Mult]; 4] = [
        &[(0, -9), (0, -3), (1, 6), (3, -1), (3, 5), (6, -2)],
        &[(0, -9), (0, 3), (1, 3), (3, -5), (3, 7), (6, -1)],
        &[(0, -8), (0, 2), (1, 3), (3, -4), (3, 6), (6, -1)],
        &[(0, -7), (0, -5), (1, 6), (3, 1), (3, 3), (6, -2)],
    ];
    let set2: HashSet<Vec<Mult>> = r2out.sols.iter().map(|(_, s)| s.clone()).collect();
    let ok_pub = published16.iter().all(|sp| set2.contains(&canonical(sp)));
    record("R2: v4.3 公表の 16 成分解 4 例をすべて含む", ok_pub, &mut checks);
    let r2set = set2;
    cert_runs.push((r2, r2out));

    let r3 = Domain { name: "R3_v52", label: "v5.2 域: 大表現込み, |Y|≤3/2, ≤5 多重項, ≤16 成分", types: all10.clone(), nymax: 9, kmax: 5, cmax: 16, cond: Cond::all() };
    let r3out = scan(&r3);
    print_scan(r3.label, &r3out, 1);
    let ok_r3 = spectrum(&r3out.sols) == BTreeMap::from([(15i64, 1usize)]) && r3out.sols[0].1 == sm_can;
    record("R3: 大表現込みでも {15:1} = SM のみ (v5.2 の再現)", ok_r3, &mut checks);
    cert_runs.push((r3, r3out));

    // ---------------- [E] 拡張軸での反例探索 ----------------
    println!("\n[E] 反例探索: 領域を広げても SM より小さい解は現れないか");
    let exts = [
        Domain { name: "E1_bigreps24", label: "E1: 大表現込み, |Y|≤3/2, ≤6 多重項, ≤24 成分", types: all10.clone(), nymax: 9, kmax: 6, cmax: 24, cond: Cond::all() },
        Domain { name: "E2_Y2", label: "E2: 大表現込み, |Y|≤2, ≤6 多重項, ≤18 成分", types: all10.clone(), nymax: 12, kmax: 6, cmax: 18, cond: Cond::all() },
        Domain { name: "E3_Y3", label: "E3: 大表現込み, |Y|≤3, ≤5 多重項, ≤16 成分", types: all10.clone(), nymax: 18, kmax: 5, cmax: 16, cond: Cond::all() },
        Domain { name: "E4_k8", label: "E4: 小表現, |Y|≤3/2, ≤8 多重項, ≤24 成分", types: small.clone(), nymax: 9, kmax: 8, cmax: 24, cond: Cond::all() },
    ];
    for dom in exts {
        let out = scan(&dom);
        print_scan(dom.label, &out, 1);
        let by = spectrum(&out.sols);
        let below: usize = by.iter().filter(|(c, _)| **c < 15).map(|(_, n)| n).sum();
        let at15 = by.get(&15).copied().unwrap_or(0);
        let has_sm = out.sols.iter().any(|(_, s)| *s == sm_can);
        record(
            &format!("{}: 15 成分未満 0 個・15 成分は SM のみ", dom.name),
            below == 0 && at15 == 1 && has_sm,
            &mut checks,
        );
        cert_runs.push((dom, out));
    }
    // E5: ν_R (完全中性シングレット) をカイラル条件から除外して許可
    {
        let mut cond = Cond::all();
        cond.sterile_ok = true;
        let dom = Domain { name: "E5_nuR", label: "E5: R2 域 + 完全中性 (1,1)_0 を許可 (ν_R)", types: small.clone(), nymax: 9, kmax: 6, cmax: 24, cond };
        let out = scan(&dom);
        print_scan(dom.label, &out, 2);
        let by = spectrum(&out.sols);
        let below: usize = by.iter().filter(|(c, _)| **c < 15).map(|(_, n)| n).sum();
        let at15 = by.get(&15).copied().unwrap_or(0);
        let has_smnur = out.sols.iter().any(|(_, s)| *s == canonical(&SM_NUR));
        record("E5: SM は最小のまま、16 成分に SM+ν_R が現れる", below == 0 && at15 == 1 && has_smnur, &mut checks);
        cert_runs.push((dom, out));
    }

    // ---------------- [N] 陰性対照: 条件を 1 つずつ外す ----------------
    println!("\n[N] 陰性対照: どの条件が SM の最小性・一意性を担っているか (基線 = R2 域)");
    println!("  各行: 外した条件 → 15 成分未満の解の数 / 15 成分の解の数 / 最小成分数");
    let drops: [(&str, &str, fn(&mut Cond)); 8] = [
        ("N1_nochiral", "カイラル性 (vectorlike 対を許す)", |c| c.chiral = false),
        ("N2_nofactors", "全因子帯電", |c| c.all_factors = false),
        ("N3_nowitten", "Witten 大域アノマリー", |c| c.witten = false),
        ("N4_nosu3cub", "SU(3)³ アノマリー", |c| c.su3cub = false),
        ("N5_nosu3sq", "SU(3)²U(1) アノマリー", |c| c.su3sq = false),
        ("N6_nosu2sq", "SU(2)²U(1) アノマリー", |c| c.su2sq = false),
        ("N7_nograv", "重力²U(1) アノマリー", |c| c.grav = false),
        ("N8_nocubic", "U(1)³ アノマリー", |c| c.cubic = false),
    ];
    let mut control_rows: Vec<(String, usize, usize, i64, usize)> = Vec::new();
    for (name, desc, dropf) in drops {
        let mut cond = Cond::all();
        dropf(&mut cond);
        let dom = Domain { name, label: desc, types: small.clone(), nymax: 9, kmax: 6, cmax: 24, cond };
        let out = scan(&dom);
        let by = spectrum(&out.sols);
        let below: usize = by.iter().filter(|(c, _)| **c < 15).map(|(_, n)| n).sum();
        let at15 = by.get(&15).copied().unwrap_or(0);
        let cmin = by.keys().next().copied().unwrap_or(0);
        let nsol = out.sols.len();
        println!(
            "  {:22} → <15: {:4} 個 / =15: {:4} 個 / 最小 {:2} 成分 / 総解 {} 個",
            desc, below, at15, cmin, nsol
        );
        if cmin < 15 {
            println!("      最小の例: {}", fmt_sp(&out.sols[0].1));
        }
        // 単調性: 条件を外すと解集合は必ず広がる (基線 ⊆ 対照) — エンジンの整合性検査を兼ねる
        let cset: HashSet<Vec<Mult>> = out.sols.iter().map(|(_, s)| s.clone()).collect();
        let superset = r2set.iter().all(|s| cset.contains(s));
        if !superset {
            record(&format!("{}: 基線 ⊆ 対照 (単調性)", name), false, &mut checks);
        }
        control_rows.push((name.to_string(), below, at15, cmin, nsol));
        cert_runs.push((dom, out));
    }
    {
        let all_superset = control_rows.len() == 8;
        record("N1–N8: 基線の解集合が全対照に含まれる (単調性)", all_superset && checks.iter().all(|(n, ok)| !n.contains("単調性") || *ok), &mut checks);
        let n1 = &control_rows[0];
        record("N1: カイラル性を外すと 15 成分未満の解 (vectorlike 対) が現れる", n1.1 > 0 && n1.3 < 15, &mut checks);
        let uniq_broken: Vec<&str> = control_rows.iter().filter(|r| r.1 > 0 || r.2 != 1).map(|r| r.0.as_str()).collect();
        println!("  最小性/一意性が壊れる対照: {:?}", uniq_broken);
        println!("  (壊れない対照 = その条件はこの領域では結論に効いていない、という地図)");
    }

    // ---------------- [U] 第 2 の U(1) は載るか ----------------
    println!("\n[U] SM の物質内容に第 2 の超電荷 Y' (|6Y'|≤12) が載るか — 全走査");
    // SM (ν_R なし): Y' = (q,u,d,l,e)。6 条件: SU(3)²Y', SU(2)²Y', 重力²Y', Y'³, Y²Y', YY'²
    {
        let y1 = [1i64, -4, 2, -3, 6]; // Q, u^c, d^c, L, e^c (×6)
        let cw = [(3i64, 2i64), (3, 1), (3, 1), (1, 2), (1, 1)]; // (色次元, 弱次元)
        let t3 = [2i64, 1, 1, 0, 0]; // 2T3·wd
        let t2 = [3i64, 0, 0, 1, 0]; // 2T2·cd
        let mut nsol = 0u64;
        let mut all_prop = true;
        let r = 12i64;
        for q in -r..=r {
            for u in -r..=r {
                for d in -r..=r {
                    for l in -r..=r {
                        for e in -r..=r {
                            let y2 = [q, u, d, l, e];
                            let mut s3 = 0i64;
                            let mut s2 = 0;
                            let mut gr = 0;
                            let mut cu = 0;
                            let mut m1 = 0;
                            let mut m2 = 0;
                            for i in 0..5 {
                                let dw = cw[i].0 * cw[i].1;
                                s3 += t3[i] * y2[i];
                                s2 += t2[i] * y2[i];
                                gr += dw * y2[i];
                                cu += dw * y2[i] * y2[i] * y2[i];
                                m1 += dw * y1[i] * y1[i] * y2[i];
                                m2 += dw * y1[i] * y2[i] * y2[i];
                            }
                            if s3 == 0 && s2 == 0 && gr == 0 && cu == 0 && m1 == 0 && m2 == 0 && y2 != [0; 5] {
                                nsol += 1;
                                // y2 ∥ y1 か (2×2 小行列式が全て 0)
                                for i in 0..5 {
                                    for j in (i + 1)..5 {
                                        if y1[i] * y2[j] - y1[j] * y2[i] != 0 {
                                            all_prop = false;
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        println!("  SM (ν_R なし): 非自明解 {} 個 — すべて Y に比例するか: {}", nsol, all_prop);
        record("U1: SM に独立な第 2 の U(1) は存在しない (全て Y ∝)", nsol > 0 && all_prop, &mut checks);
    }
    // SM+ν_R: 線形 3 条件で (d,l,ν) を消去し (q,u,e) を全走査
    {
        let y1 = [1i64, -4, 2, -3, 6, 0];
        let bl = [2i64, -2, -2, -6, 6, 6]; // B−L (×6)
        let cwdim = [6i64, 3, 3, 2, 1, 1];
        let r = 12i64;
        let mut nsol = 0u64;
        let mut in_span = true;
        let mut found_bl = false;
        for q in -r..=r {
            let l = -3 * q;
            if l.abs() > r {
                continue;
            }
            for u in -r..=r {
                let d = -2 * q - u;
                if d.abs() > r {
                    continue;
                }
                for e in -r..=r {
                    let nu = 6 * q - e;
                    if nu.abs() > r {
                        continue;
                    }
                    let y2 = [q, u, d, l, e, nu];
                    if y2 == [0; 6] {
                        continue;
                    }
                    let mut cu = 0i64;
                    let mut m1 = 0;
                    let mut m2 = 0;
                    for i in 0..6 {
                        cu += cwdim[i] * y2[i] * y2[i] * y2[i];
                        m1 += cwdim[i] * y1[i] * y1[i] * y2[i];
                        m2 += cwdim[i] * y1[i] * y2[i] * y2[i];
                    }
                    if cu == 0 && m1 == 0 && m2 == 0 {
                        nsol += 1;
                        if y2 == bl {
                            found_bl = true;
                        }
                        // span{Y, B−L} に入るか: 3×3 小行列式が全て 0
                        'outer: for a in 0..6 {
                            for b in (a + 1)..6 {
                                for c in (b + 1)..6 {
                                    let det = y1[a] * (bl[b] * y2[c] - bl[c] * y2[b])
                                        - y1[b] * (bl[a] * y2[c] - bl[c] * y2[a])
                                        + y1[c] * (bl[a] * y2[b] - bl[b] * y2[a]);
                                    if det != 0 {
                                        in_span = false;
                                        break 'outer;
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        println!("  SM+ν_R: 非自明解 {} 個 — B−L を含み、全て span{{Y, B−L}} 内: {}", nsol, found_bl && in_span);
        record("U2: SM+ν_R の第 2 U(1) は Y と B−L の張る 2 次元格子に限る", nsol > 0 && found_bl && in_span, &mut checks);
    }

    // ---------------- [G] GUT 降下の厳密検査 ----------------
    println!("\n[G] GUT 降下: 分岐則の厳密整数検査");
    {
        // SU(5): Y5 = diag(-2,-2,-2,3,3)/6 → 5̄ ⊕ 10
        let y5 = [-2i64, -2, -2, 3, 3];
        let mut content5: Vec<Mult> = Vec::new();
        // 5̄: 3̄ 側 (最初の 3 成分) が d^c=(3̄,1), レプトン側 2 成分が L=(1,2)
        content5.push((4, -y5[0])); // (3̄,1)_{2}
        content5.push((1, -y5[3])); // (1,2)_{-3}
        // 10 = 反対称 2 階: (i<j) の Y 和。色×色 → (3̄,1) [u^c], 色×弱 → (3,2) [Q], 弱×弱 → (1,1) [e^c]
        let mut content10: Vec<Mult> = Vec::new();
        content10.push((4, y5[0] + y5[1])); // u^c = (3̄,1)_{-4}
        content10.push((5, y5[0] + y5[3])); // Q = (3,2)_{1}
        content10.push((0, y5[3] + y5[4])); // e^c = (1,1)_{6}
        let mut got: Vec<Mult> = content5.clone();
        got.extend(&content10);
        let ok_su5 = canonical(&got) == sm_can && y5.iter().sum::<i64>() == 0;
        println!("  SU(5) 5̄⊕10 → {} ", fmt_sp(&got));
        record("G1: SU(5) 5̄⊕10 の超電荷 = SM 1 世代 (トレースレス条件込み)", ok_su5, &mut checks);
        // SO(10): 16 = 10 ⊕ 5̄ ⊕ 1 → SM + ν_R
        let mut got16 = got.clone();
        got16.push((0, 0));
        let ok_so10 = canonical(&got16) == canonical(&SM_NUR);
        let mut cond_s = Cond::all();
        cond_s.sterile_ok = true;
        record("G2: SO(10) 16 = 10⊕5̄⊕1 → SM+ν_R (アノマリー検査込み)", ok_so10 && full_check(&got16, &cond_s), &mut checks);
        // Pati–Salam: (4,2,1)⊕(4̄,1,2), Y = (B−L)/2 + T3R (6Y = 3(B−L)三分単位 + ±3)
        let mut ps: Vec<Mult> = Vec::new();
        ps.push((5, 1)); // (4,2,1) 色側: bl3=+1 → Q=(3,2)_{1}
        ps.push((1, -3)); // (4,2,1) レプトン側: bl3=−3 → L=(1,2)_{-3}
        ps.push((4, -1 + 3)); // (4̄,1,2) 色, T3R=+1/2 → d^c=(3̄,1)_{2}
        ps.push((4, -1 - 3)); // (4̄,1,2) 色, T3R=−1/2 → u^c=(3̄,1)_{-4}
        ps.push((0, 3 + 3)); // (4̄,1,2) レプトン, T3R=+1/2 → e^c=(1,1)_{6}
        ps.push((0, 3 - 3)); // (4̄,1,2) レプトン, T3R=−1/2 → ν_R=(1,1)_{0}
        let ok_ps = canonical(&ps) == canonical(&SM_NUR);
        println!("  Pati–Salam (4,2,1)⊕(4̄,1,2) → {}", fmt_sp(&ps));
        record("G3: Pati–Salam 分岐 (Y=(B−L)/2+T3R) → SM+ν_R", ok_ps, &mut checks);
    }

    // ---------------- 証明書の出力 ----------------
    println!("\n[C] 証明書の出力 (certificates/)");
    let cond_json = |c: &Cond| {
        Json::Obj(vec![
            ("su3cub".into(), Json::Bool(c.su3cub)),
            ("su3sq_u1".into(), Json::Bool(c.su3sq)),
            ("su2sq_u1".into(), Json::Bool(c.su2sq)),
            ("grav_u1".into(), Json::Bool(c.grav)),
            ("u1_cubed".into(), Json::Bool(c.cubic)),
            ("witten".into(), Json::Bool(c.witten)),
            ("chiral".into(), Json::Bool(c.chiral)),
            ("charged_under_all_factors".into(), Json::Bool(c.all_factors)),
            ("sterile_neutrino_allowed".into(), Json::Bool(c.sterile_ok)),
        ])
    };
    let mut domains_json = Vec::new();
    let mut sols_json = Vec::new();
    let mut sha_lines = String::new();
    const JSON_CAP: usize = 2000;
    for (dom, out) in &cert_runs {
        domains_json.push(Json::Obj(vec![
            ("run".into(), Json::Str(dom.name.into())),
            ("label".into(), Json::Str(dom.label.into())),
            ("types".into(), Json::Arr(dom.types.iter().map(|&t| Json::Str(REPS[t].key.into())).collect())),
            ("hypercharge_6y_max".into(), Json::Int(dom.nymax)),
            ("max_multiplets".into(), Json::Int(dom.kmax as i64)),
            ("max_components".into(), Json::Int(dom.cmax)),
            ("conditions".into(), cond_json(&dom.cond)),
        ]));
        // 正準解の決定的直列化と SHA-256 (JSON は上限つき、ハッシュは常に全解)
        let ser: String = out
            .sols
            .iter()
            .map(|(c, s)| format!("{}|{}", c, ascii_sp(s)))
            .collect::<Vec<_>>()
            .join("\n");
        let h = sha256_hex(ser.as_bytes());
        sha_lines.push_str(&format!("{}  {}  n={}\n", h, dom.name, out.sols.len()));
        let truncated = out.sols.len() > JSON_CAP;
        if truncated {
            println!("  注意: {} の解 {} 個は JSON では先頭 {} 個に切り詰め (ハッシュは全解)", dom.name, out.sols.len(), JSON_CAP);
        }
        sols_json.push(Json::Obj(vec![
            ("run".into(), Json::Str(dom.name.into())),
            ("n_solutions".into(), Json::Int(out.sols.len() as i64)),
            ("n_structures".into(), Json::Int(out.n_structures as i64)),
            ("truncated_in_json".into(), Json::Bool(truncated)),
            ("sha256_of_full_list".into(), Json::Str(h)),
            (
                "spectrum".into(),
                Json::Obj(spectrum(&out.sols).iter().map(|(c, n)| (c.to_string(), Json::Int(*n as i64))).collect()),
            ),
            (
                "solutions".into(),
                Json::Arr(
                    out.sols
                        .iter()
                        .take(JSON_CAP)
                        .map(|(c, s)| {
                            Json::Obj(vec![
                                ("components".into(), Json::Int(*c)),
                                ("multiplets".into(), Json::Arr(s.iter().map(|&(t, n)| Json::Str(format!("{}:{}", REPS[t].key, n))).collect())),
                            ])
                        })
                        .collect(),
                ),
            ),
        ]));
    }
    let p1 = write_artifact("certificates/v62_domains.json", &Json::Arr(domains_json).render());
    let p2 = write_artifact("certificates/v62_solutions.json", &Json::Arr(sols_json).render());
    let p3 = write_artifact("certificates/v62_sha256.txt", &sha_lines);
    println!("  {} / {} / {}", p1, p2, p3);

    // ---------------- 機械可読な結果 (results/v62_atlas.json) ----------------
    let all_ok = checks.iter().all(|(_, ok)| *ok);
    let summary = Json::Obj(vec![
        ("claim_ids".into(), Json::Arr(vec![
            Json::Str("QRN-GAUGE-003".into()),
            Json::Str("QRN-GAUGE-006".into()),
            Json::Str("QRN-GAUGE-007".into()),
            Json::Str("QRN-GAUGE-008".into()),
            Json::Str("QRN-GAUGE-009".into()),
            Json::Str("QRN-GAUGE-010".into()),
        ])),
        ("checks".into(), Json::Arr(checks.iter().map(|(n, ok)| Json::Obj(vec![("name".into(), Json::Str(n.clone())), ("pass".into(), Json::Bool(*ok))])).collect())),
        ("controls".into(), Json::Arr(control_rows.iter().map(|(n, below, at15, cmin, nsol)| Json::Obj(vec![
            ("run".into(), Json::Str(n.clone())),
            ("solutions_below_15".into(), Json::Int(*below as i64)),
            ("solutions_at_15".into(), Json::Int(*at15 as i64)),
            ("min_components".into(), Json::Int(*cmin)),
            ("total_solutions".into(), Json::Int(*nsol as i64)),
        ])).collect())),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p4 = write_artifact("results/v62_atlas.json", &summary.render());
    println!("  {}", p4);

    // ---------------- 総合判定 ----------------
    println!("\n---- 検査一覧 ----");
    for (n, ok) in &checks {
        println!("  {} {}", pass(*ok), n);
    }
    println!("\n総合判定: {}", pass(all_ok));
    println!("\n結論: SM 1 世代の最小性・一意性・孤立性は、独立アルゴリズムで再現され、");
    println!("      表現・|Y|・多重項数の拡張に頑健で、担い手はカイラル性+全因子帯電+");
    println!("      アノマリー消去の全体である (対照地図)。第 2 の U(1) は ν_R なしでは");
    println!("      Y 自身に限られ、ν_R を足すと B−L だけが加わる。探索領域と全解は");
    println!("      certificates/ に機械可読な証明書として保存された。");
    if !all_ok {
        std::process::exit(1);
    }
}
