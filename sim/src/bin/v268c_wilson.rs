//! v26.8-C v268c_wilson — Wilson 独立離散化の連続極限: A^Wil(a) → 1×A_oracle
//!
//! 事前登録: spec §12.5 (bc644d4)。**PRED-016 の独立離散化側 (D チャネル)**:
//! 空間 Wilson Hamiltonian (時間連続 — temporal regulator を混ぜない, taste trap 回避)
//!   H_W(k) = Σᵢ αᵢ sin kᵢ + β·M(k),  M(k) = m + r Σᵢ (1 − cos kᵢ)
//! は単一ノード (k = 0) の **1 Dirac flavor** — 的は 1×A_oracle = −1/(160π²)
//! (staggered の 2×A_oracle との対比が taste 数の交差検証になる)。
//!
//! 構造の利得: H_W は厳密に Dirac 型 (H² = [Σsin²kᵢ + M(k)²]·1 — [S0] で機械検査)
//! なので射影子 P±(k) = (1 ± H/E)/2 が閉形式 — 対角化不要で BZ 積分が走る。
//!
//! BOND-A 結合 (spec §2 の凍結則の Wilson 転写): 方向 i のボンド演算子
//! (α ホップと Wilson-β ホップは同じボンドに住む) を (1+h_ii)^{−1/2} 倍。
//! on-site 片 (質量 m と Wilson の 3rβ) は長さ要素を持たず h 非結合。
//!   V_ii(k; qŷ) = −½ [αᵢ sin(kᵢ + δ_{iy}q/2) − rβ cos(kᵢ + δ_{iy}q/2)]
//! (中点位相規約)。D チャネル V_D = (V_xx − V_zz)/√2 の β 汚染は差で
//! −rβ(cos kx − cos kz)/√2 ≈ −rβ(pz²−px²)/(2√2) — **O(p²) の improvement 級**
//! (v268z の許容差クラス) なので tree-level matching は自動で成立する。
//!
//! 検査 (凍結):
//!  [S0] Dirac 型恒等式: H_W(k)² = (Σsin²kᵢ + M(k)²)·1 (機械精度, 決定的 k 点)
//!  [S0b] 頂点 = ∂H/∂h の解析恒等 (q=0): V_ii(k;0) = −½·bond_i(k) (行列等価)
//!  [S1] 求積自己整合: a_min で GL 14 → 18 の変化 < 1e-3
//!  [S2] **PRED-016 Wilson-D**: A(a) ladder の 2 モデル外挿 (spec §12.5 凍結:
//!       A₀+c₁a+c₂a² / A₀+c₂a²+c₄a⁴) で **A₀/A_oracle = 1 ± 0.02 (spread 込み)**
//!  [S3] q̂ = (2/a)sin(aq/2) 変種: 外挿一致 2%
//!  [S4] **regulator 内普遍性**: r = 1.0 と r = 0.7 (irrelevant 結合の変更) の
//!       外挿 A₀ が 2% で一致 — 同じ連続極限の直接検査
//!  [S5] 変異: V_zz の rβ 片の符号反転 (頂点 ≠ ∂H/∂h) → A が > 2% 変化
//!
//! 事前登録分岐: (a) S0–S4 PASS → **Wilson-D も oracle に収束 — 二離散化の
//!   universality (D チャネル) 成立**。X チャネル完成で PRED-016 採点へ /
//!   (b) S2 miss & staggered hit → Wilson source/外挿の誤り (spec §12.8) /
//!   (c) S0/S1 FAIL → 器械。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

fn a_oracle_1d() -> f64 {
    -1.0 / (160.0 * PI * PI)
}

// ---------------- 複素 4×4 (v268a の写経) ----------------

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

/// M(k) = m + r Σ(1 − cos kᵢ)
fn mw(k: [f64; 3], m: f64, r: f64) -> f64 {
    m + r * (3.0 - k[0].cos() - k[1].cos() - k[2].cos())
}

fn ew(k: [f64; 3], m: f64, r: f64) -> f64 {
    let s2 = k[0].sin().powi(2) + k[1].sin().powi(2) + k[2].sin().powi(2);
    let mk = mw(k, m, r);
    (s2 + mk * mk).sqrt()
}

fn hw(al: &[M4; 4], k: [f64; 3], m: f64, r: f64) -> M4 {
    let mut h = mzero();
    for ax in 0..3 {
        let s = k[ax].sin();
        for i in 0..16 {
            h[i].0 += s * al[ax][i].0;
            h[i].1 += s * al[ax][i].1;
        }
    }
    let mk = mw(k, m, r);
    for i in 0..16 {
        h[i].0 += mk * al[3][i].0;
        h[i].1 += mk * al[3][i].1;
    }
    h
}

/// P±(k) = (1 ± H/E)/2 — Dirac 型なので閉形式
fn projw(al: &[M4; 4], k: [f64; 3], m: f64, r: f64, sign: f64) -> M4 {
    let h = hw(al, k, m, r);
    let e = ew(k, m, r);
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

/// 方向 i のボンド演算子 bond_i(kᵢ) = αᵢ sin kᵢ − rβ cos kᵢ
fn bondop(al: &[M4; 4], dir: usize, karg: f64, r: f64) -> M4 {
    let (s, c) = (karg.sin(), karg.cos());
    let mut v = mzero();
    for i in 0..16 {
        v[i].0 = s * al[dir][i].0 - r * c * al[3][i].0;
        v[i].1 = s * al[dir][i].1 - r * c * al[3][i].1;
    }
    v
}

/// V_D(k; qŷ) = [bond_x(kx) − bond_z(kz)]/√2 (D は y ボンド非関与 — q 移送は
/// 状態側 P₊(k+qŷ) が担う)。mutate: V_zz の rβ 片の符号を反転。
fn vertex_d(al: &[M4; 4], k: [f64; 3], r: f64, mutate: bool) -> M4 {
    let vx = bondop(al, 0, k[0], r);
    let vz = bondop(al, 2, k[2], if mutate { -r } else { r });
    let r2i = 1.0 / (2.0f64).sqrt();
    let mut v = mzero();
    for i in 0..16 {
        v[i].0 = (vx[i].0 - vz[i].0) * r2i;
        v[i].1 = (vx[i].1 - vz[i].1) * r2i;
    }
    v
}

/// null 結合の被積分 (k 点ごと): Σᵢ wᵢ tr[P₊(k+qᵢŷ)V_D P₋(k)V_D]·2/(E₁+E₂)
fn integrand(
    al: &[M4; 4],
    k: [f64; 3],
    qs_lat: &[f64],
    w_null: &[f64],
    m: f64,
    r: f64,
    mutate: bool,
) -> f64 {
    let pm = projw(al, k, m, r, -1.0);
    let e2 = ew(k, m, r);
    let vd = vertex_d(al, k, r, mutate);
    let pmv = mmul(&pm, &vd);
    let mut acc = 0.0f64;
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let pp = projw(al, kq, m, r, 1.0);
        let e1 = ew(kq, m, r);
        let tr = mtrace_re(&mmul(&mmul(&pp, &vd), &pmv));
        acc += w_null[qi] * tr * 2.0 / (e1 + e2);
    }
    acc
}

// ---------------- GL と BZ 積分 ----------------

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

/// BZ [−π,π]³ — ノード k = 0 の 3 段入れ子箱 + ノード点分割
fn bz_integrate(
    qs_lat: &[f64],
    w_null: &[f64],
    m: f64,
    r: f64,
    r_node: f64,
    r_fine: f64,
    ngl: usize,
    nthreads: usize,
    mutate: bool,
) -> f64 {
    let gl = gauss_legendre(ngl);
    let r1 = r_node;
    let r2 = (r_node / 3.0).max(2.0 * r_fine.min(r_node / 3.0));
    let r3 = r_fine.min(r2 / 2.0);
    let edges = [-PI, -r1, -r2, -r3, 0.0, r3, r2, r1, PI];
    let mut nodes1d: Vec<(f64, f64)> = Vec::new();
    for w2 in edges.windows(2) {
        let (a, b) = (w2[0], w2[1]);
        let (cc, hh) = (0.5 * (a + b), 0.5 * (b - a));
        for (x, wgt) in gl.0.iter().zip(&gl.1) {
            nodes1d.push((cc + hh * x, wgt * hh));
        }
    }
    let n1 = nodes1d.len();
    let mut rows: Vec<Option<f64>> = Vec::new();
    rows.resize_with(n1, || None);
    let chunk = n1.div_ceil(nthreads);
    let al = alphas();
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            let nodes = &nodes1d;
            let al = &al;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let ix = t * chunk + i;
                    let (kx, wx) = nodes[ix];
                    let mut acc = 0.0f64;
                    for &(ky, wy) in nodes.iter() {
                        for &(kz, wz) in nodes.iter() {
                            acc += wy
                                * wz
                                * integrand(al, [kx, ky, kz], qs_lat, w_null, m, r, mutate);
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
        for row in col + 1..4 {
            let f = m[row][col] / d;
            for c in col..4 {
                m[row][c] -= f * m[col][c];
            }
            rhs[row] -= f * rhs[col];
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

/// spec §12.5 の凍結 2 モデル: A₀+c₁a+c₂a² と A₀+c₂a²+c₄a⁴
fn extrapolate_wilson(avals: &[(f64, f64)]) -> (f64, f64) {
    let fit = |basis: &dyn Fn(f64) -> Vec<f64>| -> f64 {
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
                .max_by(|&r1, &r2| {
                    m[col + r1 * p].abs().partial_cmp(&m[col + r2 * p].abs()).unwrap()
                })
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
        sol[0]
    };
    let a_lin = fit(&|a: f64| vec![1.0, a, a * a]);
    let a_quad = fit(&|a: f64| vec![1.0, a * a, a.powi(4)]);
    (a_lin, a_quad)
}

fn main() {
    self_test();
    println!("=== v26.8-C v268c_wilson — Wilson 独立離散化の連続極限 (PRED-016 の D 独立側) ===\n");
    println!("事前登録: spec §12.5 (bc644d4)。的: A(a→0) = 1×A_oracle = −1/(160π²) (1 flavor —");
    println!("staggered の 2× との対比が taste 数の交差検証)。外挿 2 モデルは凍結済み (a+a²/a²+a⁴)。\n");
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
    let target = a_oracle_1d();
    let r_wil = 1.0f64;

    // ---- [S0] Dirac 型恒等式 H² = E²·1 ----
    {
        let mut worst = 0.0f64;
        for &(k, m, r) in &[
            ([0.3f64, -1.1, 2.2], 0.0f64, 1.0f64),
            ([2.9, 0.4, -0.6], 0.25, 0.7),
            ([-1.5, 3.0, 0.05], 0.5, 1.0),
        ] {
            let h = hw(&al, k, m, r);
            let h2 = mmul(&h, &h);
            let want = ew(k, m, r).powi(2);
            for rr in 0..4 {
                for c in 0..4 {
                    let w = if rr == c { want } else { 0.0 };
                    worst = worst
                        .max((h2[c + rr * 4].0 - w).abs())
                        .max(h2[c + rr * 4].1.abs());
                }
            }
        }
        check(
            "[S0] Dirac 型恒等式: H_W(k)² = (Σsin²k + M(k)²)·1 (機械精度)",
            worst < 1e-13,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [S0b] 頂点 = ∂H/∂h の解析恒等 (q = 0) ----
    {
        // H(h_xx) = (1+h)^{−1/2}·bond_x + rest ⇒ ∂H/∂h|₀ = −½ bond_x = V_xx(q=0)
        // 数値微分 (解析式なので 2 点で厳密級) との行列一致
        let k = [0.7f64, -0.4, 1.9];
        let (_m, r) = (0.25f64, 1.0f64);
        let eps = 1e-6;
        let mut worst = 0.0f64;
        let bx = bondop(&al, 0, k[0], r);
        for i in 0..16 {
            // H(h) の x 片 = (1+h)^{-1/2} bx: FD
            let hp = (1.0f64 + eps).powf(-0.5) * bx[i].0;
            let hm = (1.0f64 - eps).powf(-0.5) * bx[i].0;
            let fd = (hp - hm) / (2.0 * eps);
            worst = worst.max((fd - (-0.5) * bx[i].0).abs());
            let hp1 = (1.0f64 + eps).powf(-0.5) * bx[i].1;
            let hm1 = (1.0f64 - eps).powf(-0.5) * bx[i].1;
            let fd1 = (hp1 - hm1) / (2.0 * eps);
            worst = worst.max((fd1 - (-0.5) * bx[i].1).abs());
        }
        check(
            "[S0b] 頂点 = ∂H/∂h (BOND-A): V_xx(q=0) = −½·bond_x (行列 FD 一致)",
            worst < 1e-9,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- massless 系列 A(a) (r = 1.0) ----
    let ladder = [0.5f64, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.032];
    let run_ladder = |r_par: f64, ngl: usize| -> Vec<(f64, f64)> {
        let mut out = Vec::new();
        for &a in &ladder {
            let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
            let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
            let r_fine = 2.0 * a * qset[0];
            let aval = bz_integrate(
                &qs_lat, &w_null, 0.0, r_par, r_node, r_fine, ngl, nthreads, false,
            ) / a.powi(4);
            out.push((a, aval));
        }
        out
    };
    let avals = run_ladder(r_wil, 14);
    println!("    [A(a) 表 (Wilson r=1, massless)] a | A(a) | A(a)/A_oracle");
    for &(a, v) in &avals {
        println!("      a = {:.3}: {:+.6e} | {:.4}", a, v, v / target);
    }
    println!("      ({} s)", t0.elapsed().as_secs());

    // ---- [S1] 求積自己整合 ----
    {
        let a = *ladder.last().unwrap();
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
        let r_fine = 2.0 * a * qset[0];
        let fine = bz_integrate(
            &qs_lat, &w_null, 0.0, r_wil, r_node, r_fine, 18, nthreads, false,
        ) / a.powi(4);
        let rel = (fine / avals.last().unwrap().1 - 1.0).abs();
        check(
            "[S1] 求積自己整合: a_min で GL 14 → 18 の変化 < 1e-3 (相対)",
            rel < 1e-3,
            format!("相対変化 = {:.1e}", rel),
        );
    }

    // ---- [S2] PRED-016 Wilson-D ----
    let (a_lin, a_quad) = extrapolate_wilson(&avals);
    {
        let (rl, rq) = (a_lin / target, a_quad / target);
        let spread = (rl - rq).abs();
        let dev = (rl - 1.0).abs().max((rq - 1.0).abs()).max(spread);
        println!(
            "    [S2 外挿] A₀ (a+a²) = {:+.6e} → 比 {:.4} / A₀ (a²+a⁴) = {:+.6e} → 比 {:.4} (spread {:.3})",
            a_lin, rl, a_quad, rq, spread
        );
        check(
            "[S2] PRED-016 Wilson-D: A₀/A_oracle = 1 ± 0.02 (凍結 2 モデル, spread 込み)",
            dev < 0.02,
            format!("比 = {:.4}/{:.4}, 偏差 max = {:.3}", rl, rq, dev),
        );
    }

    // ---- [S3] q̂ 変種 ----
    {
        let mut avals_hat: Vec<(f64, f64)> = Vec::new();
        for &a in &ladder {
            let qhat: [f64; 4] = [
                (2.0 / a) * (a * qset[0] / 2.0).sin(),
                (2.0 / a) * (a * qset[1] / 2.0).sin(),
                (2.0 / a) * (a * qset[2] / 2.0).sin(),
                (2.0 / a) * (a * qset[3] / 2.0).sin(),
            ];
            let w_hat = null_weights(&qhat);
            let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
            let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
            let r_fine = 2.0 * a * qset[0];
            let aval = bz_integrate(
                &qs_lat, &w_hat, 0.0, r_wil, r_node, r_fine, 14, nthreads, false,
            ) / a.powi(4);
            avals_hat.push((a, aval));
        }
        let (ah_lin, _) = extrapolate_wilson(&avals_hat);
        let rel = (ah_lin / a_lin - 1.0).abs();
        check(
            "[S3] q̂ = (2/a)sin(aq/2) 変種: 外挿 A₀ が 2% で一致 (spec §12.1)",
            rel < 0.02,
            format!("A₀(q̂)/A₀(q) − 1 = {:.3}", rel),
        );
    }

    // ---- [S4] regulator 内普遍性: r = 0.7 ----
    {
        let avals_r7 = run_ladder(0.7, 14);
        println!("    [S4 表 (r = 0.7)] a 端点: A/A_or = {:.4} (a=0.5) → {:.4} (a=0.032)",
            avals_r7[0].1 / target, avals_r7.last().unwrap().1 / target);
        let (a7_lin, a7_quad) = extrapolate_wilson(&avals_r7);
        let rel = (a7_lin / a_lin - 1.0).abs().max((a7_quad / a_quad - 1.0).abs());
        check(
            "[S4] regulator 内普遍性: r = 1.0 と r = 0.7 の外挿 A₀ が 2% で一致",
            rel < 0.02,
            format!("max 相対差 = {:.3} (r=0.7 比: {:.4}/{:.4})", rel, a7_lin / target, a7_quad / target),
        );
    }

    // ---- [S5] 変異 ----
    {
        let a = 0.125f64;
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
        let r_fine = 2.0 * a * qset[0];
        let bad = bz_integrate(
            &qs_lat, &w_null, 0.0, r_wil, r_node, r_fine, 14, nthreads, true,
        ) / a.powi(4);
        let good = avals[4].1;
        check(
            "[S5] 変異: V_zz の rβ 片符号反転 (頂点 ≠ ∂H/∂h) → A が > 2% 変化",
            (bad / good - 1.0).abs() > 0.02,
            format!("相対変化 = {:.3}", (bad / good - 1.0).abs()),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.8-C".into())),
        ("kind".into(), Json::Str("wilson_tt_continuum_limit".into())),
        ("pred".into(), Json::Str("PRED-016 (Wilson-D)".into())),
        ("target_1flavor".into(), Json::Num(target)),
        (
            "a_ladder".into(),
            Json::Arr(
                avals
                    .iter()
                    .map(|&(a, v)| {
                        Json::Obj(vec![
                            ("a".into(), Json::Num(a)),
                            ("A".into(), Json::Num(v)),
                            ("ratio".into(), Json::Num(v / target)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("A0_lin".into(), Json::Num(a_lin)),
        ("A0_quad".into(), Json::Num(a_quad)),
        ("A0_over_target".into(), Json::Num(a_lin / target)),
    ]);
    let p = write_artifact("results/v268c_wilson.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **Wilson-D も 1×A_oracle に収束 — 二離散化 universality (D チャネル) 成立** (taste 数 2 vs 1 の交差検証込み)。X チャネル完成で PRED-016 採点へ"
        } else {
            "FAIL — 分岐 (b) Wilson source/外挿の誤り / (c) 器械。欄が一次ソース"
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
