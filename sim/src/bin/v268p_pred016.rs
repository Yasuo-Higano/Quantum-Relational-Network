//! v26.8-P v268p_pred016 — PRED-016 の 4 比を導出モデルで再採点 (修正条項 1)
//!
//! 事前登録: spec §12.9 (08a1321 — 本バイナリの実装・走行より先にコミット済み)。
//! v26.8-X の interim 判定 (4/4 比 = 1 ± 1.3%、登録バー 1%/系統 0.5% 未達) を、
//! **観測量自身の漸近形から導出した外挿モデル**で再採点する。的・観測量・バーは
//! 不変。修正は 1 回限り (届かなければ interim を維持し原因分析は独立ユニット)。
//!
//! 導出モデル (§12.9):
//!   staggered: null 重みは Σw = ΣwQ² = ΣwQ⁴ = 0 — χ^lat の q⁴ln(aQ)² が持つ
//!     ln a 片は ΣwQ⁴ = 0 が消すが、a² 次補正 a²Q⁶ln(aQ)² は ΣwQ⁶ ≠ 0 で生存
//!     → δA = c·a²ln(1/a) + d·a² — **モデル {1, a²ln(1/a), a²}** (D・X 共通)
//!   Wilson: r 項 (カイラル破れ) の O(a) が先頭 → δA = b·a + d·a²
//!     — **モデル {1, a, a²}** (D・X 共通)
//! 中心値 = 尾部窓 (a ≤ 0.125) フィット、系統 = |全域 − 尾部| spread。
//!
//! 器械: 4 チャネルの積分核は認証済みバイナリの写経 —
//!   stag D: v268b (V_D = (V_xx − V_zz)/√2, 一般頂点公式, 8 成分 jacobi)
//!   stag X: v268x (v268z 認証の 4 隅 point-split, 厳密 taste-singlet, Z = 2)
//!   Wil D: v268c (bond_i = αᵢsin kᵢ − rβcos kᵢ, 射影子閉形式, r = 1)
//!   Wil X: v268x (V = ½(αx sin kz + αz sin kx), Z = 1)
//! 求積は v268s の入れ子半径方式 (×3 等比, 最内 ~1.5·a·Q_min)。
//! ladder: a ∈ {0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.032, 0.022}。
//! 的: stag D → 2A_or, stag X → 4A_or, Wil D → A_or, Wil X → A_or/2
//! (A_or = −1/(160π²), v268a 三重検証)。
//!
//! 検査 (凍結, §12.9):
//!  [P0] 転記なし回帰: 各チャネルの a = 0.125 rung が公表済み ladder (v268b/c/x
//!       の results JSON) と 2e-3 相対で一致 (求積方式差込み)
//!  [P1] 求積自己整合: a_min = 0.022 で GL 14 → 18 < 1e-3 (4 チャネル)
//!  [P2] モデル妥当性: 全域フィットの最大相対残差 < 0.3% (4 チャネル)
//!  [P3] **PRED-016 中心**: 4 比すべて |尾部窓 S₀/的 − 1| ≤ 0.01
//!  [P4] **PRED-016 系統**: 4 比すべて |全域 − 尾部| ≤ 0.005
//!
//! 事前登録分岐: (a) P0–P4 PASS → **PRED-016 scored-hit (登録バー到達)** /
//!   (b) P3/P4 のみ破れ → interim 維持 (再修正禁止 — §12.9) /
//!   (c) P0–P2 FAIL → 器械 (採点せず修理のみ公表)。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

fn a_oracle_1d() -> f64 {
    -1.0 / (160.0 * PI * PI)
}

/// 公表済み a = 0.125 rung (的に対する比) — P0 回帰の的 (一次ソース: results/
/// v268b_continuum.json, v268c_wilson.json, v268x_completion.json)
const REF_A0125: [(&str, f64); 4] = [
    ("stagD", 1.022196610668463),
    ("wilD", 1.0396001713192307),
    ("stagX", 0.9852942211046332),
    ("wilX", 1.0234937128826342),
];

// ================= staggered 8 成分エンジン (v268b/v268x の写経) =================

fn sbit(s: usize, ax: usize) -> usize {
    (s >> ax) & 1
}

fn h8(k: [f64; 3], m: f64) -> Vec<f64> {
    let mut h = vec![0.0f64; 64];
    for s in 0..8usize {
        let cx = if sbit(s, 0) == 0 { 1.0 } else { -1.0 } * k[0].cos();
        h[s + s * 8] += cx;
        let s2 = s ^ 1;
        let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * k[1].cos();
        h[s2 + s * 8] += cy;
        let s3 = s ^ 3;
        let cz = if sbit(s, 2) == 0 { 1.0 } else { -1.0 } * k[2].cos();
        h[s3 + s * 8] += cz;
        let s4 = s ^ 7;
        h[s4 + s * 8] += m;
    }
    h
}

struct Term {
    eps: usize,
    d: [i32; 3],
    w: f64,
}

fn vertex8(terms: &[Term], k: [f64; 3], q: f64) -> Vec<f64> {
    let mut v = vec![0.0f64; 64];
    for t in terms {
        for s in 0..8usize {
            let s2 = s ^ t.eps;
            let mut ph1 = 0.0f64;
            let mut ph2 = 0.0f64;
            for ax in 0..3 {
                let ka = k[ax] + PI * sbit(s, ax) as f64;
                let ka2 = k[ax] + if ax == 1 { q } else { 0.0 } + PI * sbit(s2, ax) as f64;
                ph1 += ka * t.d[ax] as f64;
                ph2 += ka2 * t.d[ax] as f64;
            }
            let mid = 0.5 * q * t.d[1] as f64;
            v[s + s2 * 8] += t.w * ((ph1 + mid).cos() + (-ph2 + mid).cos());
        }
    }
    v
}

/// D 頂点素材: V_xx / V_zz (v268b と同一)
fn t_xx() -> Vec<Term> {
    vec![Term { eps: 0, d: [1, 0, 0], w: 0.5 }]
}
fn t_zz() -> Vec<Term> {
    vec![Term { eps: 3, d: [0, 0, 1], w: 0.5 }]
}

/// X 頂点: v268z 認証済み 4 隅 point-split (両片 w = −σρ/16, v268x と同一)
fn t_x_split() -> Vec<Term> {
    let mut v = Vec::new();
    for sg in [1i32, -1] {
        for rh in [1i32, -1] {
            let c = (sg * rh) as f64 / 16.0;
            v.push(Term { eps: 0, d: [sg, 0, 2 * rh], w: -c });
            v.push(Term { eps: 3, d: [2 * rh, 0, sg], w: -c });
        }
    }
    v
}

/// staggered null 結合被積分 (k 点ごと)。xch: false = D ((vx−vz)/√2), true = X
fn stag_integrand(k: [f64; 3], qs_lat: &[f64], w_null: &[f64], xch: bool) -> f64 {
    let hk = h8(k, 0.0);
    let (wk, vk) = jacobi_eigh(&hk, 8);
    let r2i = 1.0 / (2.0f64).sqrt();
    let mut acc = 0.0f64;
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let hq = h8(kq, 0.0);
        let (wq, vq) = jacobi_eigh(&hq, 8);
        let vv: Vec<f64> = if xch {
            vertex8(&t_x_split(), k, q)
        } else {
            let vx = vertex8(&t_xx(), k, q);
            let vz = vertex8(&t_zz(), k, q);
            (0..64).map(|i| (vx[i] - vz[i]) * r2i).collect()
        };
        let mut chi = 0.0f64;
        for mu in 4..8 {
            for nu in 0..4 {
                let mut mel = 0.0f64;
                for r in 0..8 {
                    let mut s = 0.0f64;
                    for c in 0..8 {
                        s += vv[c + r * 8] * vk[c + nu * 8];
                    }
                    mel += vq[r + mu * 8] * s;
                }
                chi += 2.0 * mel * mel / (wq[mu] - wk[nu]);
            }
        }
        acc += w_null[qi] * chi;
    }
    acc
}

// ================= Wilson 4 成分エンジン (v268c/v268x の写経) =================

type M4 = [(f64, f64); 16];

fn mzero() -> M4 {
    [(0.0, 0.0); 16]
}

fn mmul(a: &M4, b: &M4) -> M4 {
    let mut o = mzero();
    for r in 0..4 {
        for k in 0..4 {
            let av = a[k + r * 4];
            if av.0 == 0.0 && av.1 == 0.0 {
                continue;
            }
            for c in 0..4 {
                let bv = b[c + k * 4];
                o[c + r * 4].0 += av.0 * bv.0 - av.1 * bv.1;
                o[c + r * 4].1 += av.0 * bv.1 + av.1 * bv.0;
            }
        }
    }
    o
}

fn mtrace_re(a: &M4) -> f64 {
    (0..4).map(|i| a[i + i * 4].0).sum()
}

fn alphas() -> [M4; 4] {
    let mut ax = mzero();
    let mut ay = mzero();
    let mut az = mzero();
    let mut b = mzero();
    let put = |m: &mut M4, r: usize, c: usize, re: f64, im: f64| {
        m[c + r * 4] = (re, im);
    };
    for blk in 0..2 {
        let (ro, co) = if blk == 0 { (0, 2) } else { (2, 0) };
        put(&mut ax, ro, co + 1, 1.0, 0.0);
        put(&mut ax, ro + 1, co, 1.0, 0.0);
        put(&mut ay, ro, co + 1, 0.0, -1.0);
        put(&mut ay, ro + 1, co, 0.0, 1.0);
        put(&mut az, ro, co, 1.0, 0.0);
        put(&mut az, ro + 1, co + 1, -1.0, 0.0);
    }
    for i in 0..4 {
        put(&mut b, i, i, if i < 2 { 1.0 } else { -1.0 }, 0.0);
    }
    [ax, ay, az, b]
}

fn mwil(k: [f64; 3], m: f64, r: f64) -> f64 {
    m + r * (3.0 - k[0].cos() - k[1].cos() - k[2].cos())
}

fn ewil(k: [f64; 3], m: f64, r: f64) -> f64 {
    let s2 = k[0].sin().powi(2) + k[1].sin().powi(2) + k[2].sin().powi(2);
    let mk = mwil(k, m, r);
    (s2 + mk * mk).sqrt()
}

fn projw(al: &[M4; 4], k: [f64; 3], m: f64, r: f64, sign: f64) -> M4 {
    let mut h = mzero();
    for ax in 0..3 {
        let s = k[ax].sin();
        for i in 0..16 {
            h[i].0 += s * al[ax][i].0;
            h[i].1 += s * al[ax][i].1;
        }
    }
    let mk = mwil(k, m, r);
    for i in 0..16 {
        h[i].0 += mk * al[3][i].0;
        h[i].1 += mk * al[3][i].1;
    }
    let e = ewil(k, m, r);
    let mut p = mzero();
    for i in 0..16 {
        p[i].0 = sign * h[i].0 / (2.0 * e);
        p[i].1 = sign * h[i].1 / (2.0 * e);
    }
    for i in 0..4 {
        p[i + i * 4].0 += 0.5;
    }
    p
}

/// bond_i(kᵢ) = αᵢ sin kᵢ − rβ cos kᵢ (v268c)
fn bondop(al: &[M4; 4], dir: usize, karg: f64, r: f64) -> M4 {
    let (s, c) = (karg.sin(), karg.cos());
    let mut v = mzero();
    for i in 0..16 {
        v[i].0 = s * al[dir][i].0 - r * c * al[3][i].0;
        v[i].1 = s * al[dir][i].1 - r * c * al[3][i].1;
    }
    v
}

/// Wilson 頂点: D = [bond_x − bond_z]/√2 / X = ½(αx sin kz + αz sin kx)
fn vertex_wil(al: &[M4; 4], k: [f64; 3], r: f64, xch: bool) -> M4 {
    if xch {
        let mut v = mzero();
        let (sz, sx) = (k[2].sin(), k[0].sin());
        for i in 0..16 {
            v[i].0 = 0.5 * (sz * al[0][i].0 + sx * al[2][i].0);
            v[i].1 = 0.5 * (sz * al[0][i].1 + sx * al[2][i].1);
        }
        v
    } else {
        let vx = bondop(al, 0, k[0], r);
        let vz = bondop(al, 2, k[2], r);
        let r2i = 1.0 / (2.0f64).sqrt();
        let mut v = mzero();
        for i in 0..16 {
            v[i].0 = (vx[i].0 - vz[i].0) * r2i;
            v[i].1 = (vx[i].1 - vz[i].1) * r2i;
        }
        v
    }
}

fn wil_integrand(al: &[M4; 4], k: [f64; 3], qs_lat: &[f64], w_null: &[f64], xch: bool) -> f64 {
    let r = 1.0f64;
    let pm = projw(al, k, 0.0, r, -1.0);
    let e2 = ewil(k, 0.0, r);
    let vd = vertex_wil(al, k, r, xch);
    let pmv = mmul(&pm, &vd);
    let mut acc = 0.0f64;
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let pp = projw(al, kq, 0.0, r, 1.0);
        let e1 = ewil(kq, 0.0, r);
        let tr = mtrace_re(&mmul(&mmul(&pp, &vd), &pmv));
        acc += w_null[qi] * tr * 2.0 / (e1 + e2);
    }
    acc
}

// ================= 共通: GL / 入れ子 edges / BZ / 重み / fit =================

fn gauss_legendre(n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut xs = vec![0.0f64; n];
    let mut ws = vec![0.0f64; n];
    for i in 0..n {
        let mut x = (PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        for _ in 0..100 {
            let (mut p0, mut p1) = (1.0f64, x);
            for kk in 2..=n {
                let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
                p0 = p1;
                p1 = p2;
            }
            let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
            let dx = p1 / dp;
            x -= dx;
            if dx.abs() < 1e-15 {
                break;
            }
        }
        xs[i] = x;
        let (mut p0, mut p1) = (1.0f64, x);
        for kk in 2..=n {
            let p2 = ((2 * kk - 1) as f64 * x * p1 - (kk - 1) as f64 * p0) / kk as f64;
            p0 = p1;
            p1 = p2;
        }
        let dp = n as f64 * (x * p1 - p0) / (x * x - 1.0);
        ws[i] = 2.0 / ((1.0 - x * x) * dp * dp);
    }
    (xs, ws)
}

/// 入れ子半径 edges (v268s 方式): ×3 等比, 最内 ~1.5·scale
fn nest_edges(center: f64, lo: f64, hi: f64, scale: f64) -> Vec<f64> {
    let mut rs: Vec<f64> = Vec::new();
    let mut r = (1.5 * scale).max(0.006);
    while r < 1.2 {
        rs.push(r);
        r *= 3.0;
    }
    rs.push(1.2);
    let mut e = vec![lo];
    for &rr in rs.iter().rev() {
        e.push(center - rr);
    }
    e.push(center);
    for &rr in rs.iter() {
        e.push(center + rr);
    }
    e.push(hi);
    e
}

fn make_nodes(edges: &[f64], gl: &(Vec<f64>, Vec<f64>)) -> Vec<(f64, f64)> {
    let mut out = Vec::new();
    for w2 in edges.windows(2) {
        let (a, b) = (w2[0], w2[1]);
        let (cc, hh) = (0.5 * (a + b), 0.5 * (b - a));
        for (x, wgt) in gl.0.iter().zip(&gl.1) {
            out.push((cc + hh * x, wgt * hh));
        }
    }
    out
}

fn bz_sum<F: Fn([f64; 3]) -> f64 + Sync>(nodes: &[(f64, f64)], nthreads: usize, f: F) -> f64 {
    let n1 = nodes.len();
    let mut rows: Vec<Option<f64>> = Vec::new();
    rows.resize_with(n1, || None);
    let chunk = n1.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            let f = &f;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let ix = t * chunk + i;
                    let (kx, wx) = nodes[ix];
                    let mut acc = 0.0f64;
                    for &(ky, wy) in nodes.iter() {
                        for &(kz, wz) in nodes.iter() {
                            acc += wy * wz * f([kx, ky, kz]);
                        }
                    }
                    *slot = Some(acc * wx);
                }
            });
        }
    });
    rows.into_iter().map(|o| o.unwrap()).sum::<f64>() / (2.0 * PI).powi(3)
}

fn null_weights(qs: &[f64; 4]) -> [f64; 4] {
    let mut a = [[0.0f64; 4]; 4];
    let mut rhs = [0.0f64; 4];
    for (i, &q) in qs.iter().enumerate() {
        a[0][i] = 1.0;
        a[1][i] = q * q;
        a[2][i] = q.powi(4);
        a[3][i] = q.powi(4) * (q * q).ln();
    }
    rhs[3] = 1.0;
    let mut m = a;
    for col in 0..4 {
        let piv = (col..4)
            .max_by(|&r1, &r2| m[r1][col].abs().partial_cmp(&m[r2][col].abs()).unwrap())
            .unwrap();
        m.swap(col, piv);
        rhs.swap(col, piv);
        let d = m[col][col];
        for r in col + 1..4 {
            let f = m[r][col] / d;
            for c in col..4 {
                m[r][c] -= f * m[col][c];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    let mut w = [0.0f64; 4];
    for col in (0..4).rev() {
        let mut s = rhs[col];
        for c in col + 1..4 {
            s -= m[col][c] * w[c];
        }
        w[col] = s / m[col][col];
    }
    w
}

/// 最小二乗 fit — (S₀, 最大相対残差) を返す
fn fit_a0_resid(avals: &[(f64, f64)], basis: &dyn Fn(f64) -> Vec<f64>) -> (f64, f64) {
    let p = basis(avals[0].0).len();
    let mut ata = vec![0.0f64; p * p];
    let mut atb = vec![0.0f64; p];
    for &(a, y) in avals {
        let bs = basis(a);
        for i in 0..p {
            for j in 0..p {
                ata[j + i * p] += bs[i] * bs[j];
            }
            atb[i] += bs[i] * y;
        }
    }
    let mut m = ata;
    let mut rr = atb;
    for col in 0..p {
        let piv = (col..p)
            .max_by(|&r1, &r2| m[col + r1 * p].abs().partial_cmp(&m[col + r2 * p].abs()).unwrap())
            .unwrap();
        for c in 0..p {
            m.swap(c + col * p, c + piv * p);
        }
        rr.swap(col, piv);
        let d = m[col + col * p];
        for row in col + 1..p {
            let f = m[col + row * p] / d;
            for c in col..p {
                m[c + row * p] -= f * m[c + col * p];
            }
            rr[row] -= f * rr[col];
        }
    }
    let mut sol = vec![0.0f64; p];
    for col in (0..p).rev() {
        let mut s = rr[col];
        for c in col + 1..p {
            s -= m[c + col * p] * sol[c];
        }
        sol[col] = s / m[col + col * p];
    }
    let mut worst = 0.0f64;
    for &(a, y) in avals {
        let bs = basis(a);
        let pred: f64 = (0..p).map(|i| sol[i] * bs[i]).sum();
        worst = worst.max(((pred - y) / y).abs());
    }
    (sol[0], worst)
}

fn main() {
    self_test();
    println!("=== v26.8-P v268p_pred016 — PRED-016 の 4 比を導出モデルで再採点 (spec §12.9) ===\n");
    println!("修正条項 1 (08a1321 で走行前凍結)。的: stagD → 2A_or / stagX → 4A_or /");
    println!("wilD → A_or / wilX → A_or/2 (A_or = −1/(160π²))。中心 = 尾部窓 (a ≤ 0.125)、");
    println!("系統 = |全域 − 尾部|。バー: |中心 − 1| ≤ 1% かつ spread ≤ 0.5% (4 比すべて)。\n");
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
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
    let al = alphas();
    let qset = [0.3f64, 0.6, 0.9, 1.2];
    let w_null = null_weights(&qset);
    let aor = a_oracle_1d();
    // 的 (spec §12.9): [stagD, wilD, stagX, wilX]
    let targets = [2.0 * aor, aor, 4.0 * aor, aor / 2.0];
    let names = ["stagD", "wilD", "stagX", "wilX"];

    let ladder = [0.35f64, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.032, 0.022];

    // チャネル計算: ch 0=stagD 1=wilD 2=stagX 3=wilX
    let run_ch = |ch: usize, a: f64, ngl: usize| -> f64 {
        let gl = gauss_legendre(ngl);
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let scale = a * qset[0];
        let edges = if ch % 2 == 0 {
            nest_edges(PI / 2.0, 0.0, PI, scale)
        } else {
            nest_edges(0.0, -PI, PI, scale)
        };
        let nodes = make_nodes(&edges, &gl);
        let raw = match ch {
            0 => bz_sum(&nodes, nthreads, |k| stag_integrand(k, &qs_lat, &w_null, false)),
            1 => bz_sum(&nodes, nthreads, |k| wil_integrand(&al, k, &qs_lat, &w_null, false)),
            2 => bz_sum(&nodes, nthreads, |k| stag_integrand(k, &qs_lat, &w_null, true)),
            _ => bz_sum(&nodes, nthreads, |k| wil_integrand(&al, k, &qs_lat, &w_null, true)),
        };
        raw / a.powi(4)
    };

    // ---- ladder 4 本 ----
    let mut lads: Vec<Vec<(f64, f64)>> = vec![Vec::new(); 4];
    println!("    [A(a)/的 表] a | stagD | wilD | stagX | wilX");
    for &a in &ladder {
        let mut row = String::new();
        for ch in 0..4 {
            let v = run_ch(ch, a, 14);
            lads[ch].push((a, v));
            row = format!("{} {:.4} |", row, v / targets[ch]);
        }
        println!("      a = {:.4}: {} ({} s)", a, row, t0.elapsed().as_secs());
    }

    // ---- [P0] 転記なし回帰 (a = 0.125 = index 3) ----
    {
        let mut worst = 0.0f64;
        let mut msg = String::new();
        for ch in 0..4 {
            let got = lads[ch][3].1 / targets[ch];
            let refv = REF_A0125[ch].1;
            let rel = (got / refv - 1.0).abs();
            worst = worst.max(rel);
            msg = format!("{} {}={:.4}/{:.4}", msg, names[ch], got, refv);
        }
        check(
            "[P0] 回帰: a = 0.125 rung が公表 JSON (v268b/c/x) と 2e-3 相対一致",
            worst < 2e-3,
            format!("max 相対差 = {:.1e} —{}", worst, msg),
        );
    }

    // ---- [P1] 求積自己整合 ----
    {
        let mut worst = 0.0f64;
        for ch in 0..4 {
            let fine = run_ch(ch, *ladder.last().unwrap(), 18);
            worst = worst.max((fine / lads[ch].last().unwrap().1 - 1.0).abs());
        }
        check(
            "[P1] 求積自己整合: a_min = 0.022 で GL 14 → 18 < 1e-3 (4 チャネル)",
            worst < 1e-3,
            format!("max 相対変化 = {:.1e} ({} s)", worst, t0.elapsed().as_secs()),
        );
    }

    // ---- fit: 導出モデル, 全域 + 尾部窓 ----
    let b_stag = |a: f64| vec![1.0, a * a * (1.0 / a).ln(), a * a];
    let b_wil = |a: f64| vec![1.0, a, a * a];
    let mut centrals = [0.0f64; 4];
    let mut spreads = [0.0f64; 4];
    let mut resids = [0.0f64; 4];
    println!("    [外挿表 (S₀/的)] ch | 全域 | 尾部 (a ≤ 0.125) | spread | 全域残差");
    for ch in 0..4 {
        let basis: &dyn Fn(f64) -> Vec<f64> = if ch % 2 == 0 { &b_stag } else { &b_wil };
        let (full, resid) = fit_a0_resid(&lads[ch], basis);
        let (tail, _) = fit_a0_resid(&lads[ch][3..], basis);
        let (rf, rt) = (full / targets[ch], tail / targets[ch]);
        centrals[ch] = rt;
        spreads[ch] = (rf - rt).abs();
        resids[ch] = resid;
        println!(
            "      {}: {:.4} | {:.4} | {:.4} | {:.2e}",
            names[ch], rf, rt, spreads[ch], resid
        );
    }

    // ---- [P2] モデル妥当性 ----
    {
        let worst = resids.iter().cloned().fold(0.0f64, f64::max);
        check(
            "[P2] モデル妥当性: 全域フィットの最大相対残差 < 0.3% (4 チャネル)",
            worst < 3e-3,
            format!("max 残差 = {:.2e}", worst),
        );
    }

    // ---- [P3] PRED-016 中心 ----
    {
        let worst = centrals
            .iter()
            .map(|&c| (c - 1.0).abs())
            .fold(0.0f64, f64::max);
        check(
            "[P3] PRED-016 中心: 4 比すべて |尾部窓 S₀/的 − 1| ≤ 1%",
            worst <= 0.01,
            format!("max |中心 − 1| = {:.4}", worst),
        );
    }

    // ---- [P4] PRED-016 系統 ----
    {
        let worst = spreads.iter().cloned().fold(0.0f64, f64::max);
        check(
            "[P4] PRED-016 系統: 4 比すべて |全域 − 尾部| ≤ 0.5%",
            worst <= 0.005,
            format!("max spread = {:.4}", worst),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.8-P".into())),
        ("kind".into(), Json::Str("pred016_rescore_derived_models".into())),
        ("pred".into(), Json::Str("PRED-016".into())),
        ("spec".into(), Json::Str("§12.9 (08a1321)".into())),
        (
            "ratios_central_tail".into(),
            Json::Arr(
                (0..4)
                    .map(|ch| {
                        Json::Obj(vec![
                            ("ch".into(), Json::Str(names[ch].into())),
                            ("central".into(), Json::Num(centrals[ch])),
                            ("spread".into(), Json::Num(spreads[ch])),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "ladders".into(),
            Json::Arr(
                (0..4)
                    .map(|ch| {
                        Json::Obj(vec![
                            ("ch".into(), Json::Str(names[ch].into())),
                            (
                                "rungs".into(),
                                Json::Arr(
                                    lads[ch]
                                        .iter()
                                        .map(|&(a, v)| {
                                            Json::Obj(vec![
                                                ("a".into(), Json::Num(a)),
                                                ("ratio".into(), Json::Num(v / targets[ch])),
                                            ])
                                        })
                                        .collect(),
                                ),
                            ),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v268p_pred016.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **PRED-016 scored-hit — 4 比が登録バー (1%・系統 0.5%) で成立** (operator/regulator universality の完成。gravitational response とは呼ばない)"
        } else {
            "FAIL あり — 分岐 (b) interim 維持 (§12.9 により再修正禁止) / (c) 器械。欄が一次ソース"
        }
    );
    println!(
        "\n総合判定: {} ({} s)",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" },
        t0.elapsed().as_secs()
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
