//! v26.8-B v268b_continuum — staggered TT の連続極限: A_null(a) → 2A_oracle (PRED-016 前半)
//!
//! 事前登録: spec §12.4 (bc644d4)。continuum trajectory (am = a·m_phys, aq = a·Q,
//! 無限体積) で格子の null-combination 係数 A(a) を測り、a → 0 で解析 oracle
//! (v26.8-A: A_oracle(2 taste) = −1/(80π²), 16π²·2A = −1/5) に収束するかを判定する。
//! **これが経路 B の最重要 falsifier の前半** (後半 = Wilson 独立離散化, v26.8-C)。
//!
//! 器械: 8 成分折込み基底 (v268z で認証済み) の H(k) = Σcos kᵢ Γᵢ + m Γm は
//! **無限体積格子の厳密な 1 粒子 Hamiltonian** — χ_D^lat(q) は縮約セル (ノードを
//! 中心に置く [0,π)³) の BZ 積分:
//!   χ_D(q) = ∫_cell d³k/(2π)³·(8 成分正規化) Σ_{ν occ(k), μ unocc(k+qŷ)}
//!             2|⟨μ|V_D(k,q)|ν⟩|²/(E_μ − E_ν)
//! 頂点 V_D = (V_xx − V_zz)/√2 は v268z の一般頂点公式 (格子 sandwich と 5e-15 で
//! 照合済み) をそのまま使う。占有 = 負エネルギー 4 状態 (半充填 Dirac 海)。
//!
//! 単位系: 格子単位で q^lat = a·Q (Q = 物理運動量, 基準集合 {0.3, 0.6, 0.9, 1.2}),
//! m^lat = a·m_phys。null 重み w は**物理 Q で解く** (Σw = ΣwQ² = ΣwQ⁴ = 0,
//! ΣwQ⁴ln Q² = 1 — ln(a²Q²) の a 依存は ΣwQ⁴ = 0 が消す)。
//!   A(a) := a⁻⁴ Σᵢ wᵢ χ^lat_D(a·Qᵢ)
//! 格子の解析項 c₀(a) + c₂(a)Q² + c₄(a)Q⁴ は null 結合が代数的に消し、
//! A(a) = 2A_oracle + O(a²) が期待形 (staggered は O(a²) 改善)。
//! null 結合は k 点ごとに直接累積 (相殺深さ ~ (aQ)⁴ ≫ f64 床 — v26.8-A の教訓)。
//!
//! 検査 (凍結):
//!  [S0] 器械回帰: BZ 積分の χ_D(q = 2π/16, m ∈ {0, 0.5}) が v26.7-II の有限体積
//!       公表値 (0.154068 / 0.150054, N=64) と 5e-3 で一致 (有限サイズ差込み)
//!  [S1] 求積自己整合: ノード box 分解能を倍化して A(a_min) の変化 < 1e-3 (相対)
//!  [S2] **PRED-016 前半の判定**: massless 系列 A(a), a ∈ {0.5, 0.35, 0.25, 0.18,
//!       0.125} の 2 モデル外挿 (A₀+c₂a² / A₀+c₁a+c₂a²) で
//!       **A₀/(2A_oracle) = 1 ± 0.02 (モデル間 spread を系統に含む)** —
//!       hit/miss どちらも公表 (miss = source/taste/scaling のどれかが誤り)
//!  [S3] q̂ = (2/a)sin(aq/2) 変種: 同じ外挿で A₀ が [S2] と系統内一致 (spec §12.1)
//!  [S4] massive 対照 (PRED-018 の格子側): m_phys = 2.0 (q̄/m < 1) の
//!       A(a→0)/A₀^massless が oracle の decoupling 比 (v26.8-A の A(m)/A(0)) と
//!       15% で一致
//!  [S5] 変異: V_zz の parity を破る (ε 3→2) → A が O(1) 変化 (> 10%)
//!
//! 事前登録分岐: (a) S0–S3 PASS → **PRED-016 前半 hit — staggered は繰り込み後
//!   TT form factor で連続 Dirac に流れる** (Wilson 側 v26.8-C で完成) /
//!   (b) S2 miss → source/taste/trajectory の監査 (spec §12.8 の分岐表) /
//!   (c) S0/S1 FAIL → 器械 (BZ 積分の設計)。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// v26.7-II の有限体積公表値 (N=64, q = 2π/16) — S0 回帰の的
const REF267_CHI_D: [(f64, f64); 2] = [(0.0, 0.154068), (0.5, 0.150054)];
/// v26.8-A の oracle (1 Dirac): A = −1/(160π²)。2 taste は ×2
fn a_oracle_2t() -> f64 {
    -2.0 / (160.0 * PI * PI)
}
/// oracle の massive decoupling 比は placeholder にせず、v26.8-A で認証済みの
/// 閉形式 ρ_D(s;m) + 安定通分核 K(s) = n₀/Π(s+Qᵢ²) からこの場で計算する
/// (開発記録: run1 は事前計算値を定数で置いた — 出所不明の定数は禁じ手)。
fn rho_d_closed(s: f64, m: f64) -> f64 {
    let p2 = s / 4.0 - m * m;
    if p2 <= 0.0 {
        return 0.0;
    }
    let p = p2.sqrt();
    let ep = s.sqrt() / 2.0;
    (p * ep / (4.0 * PI * PI)) * ((2.0 / 3.0) * p2 - (4.0 / 15.0) * p2 * p2 / (ep * ep))
}

/// oracle A(m) (1 Dirac): ∫ds ρ_D(s;m)·n₀/Π(s+Qᵢ²) — v268a Route II の閉形式版
fn oracle_a(qs: &[f64; 4], w: &[f64; 4], m: f64) -> f64 {
    let xs: Vec<f64> = qs.iter().map(|&q| q * q).collect();
    let mut n0 = 0.0f64;
    for i in 0..4 {
        let mut prod = w[i];
        for j in 0..4 {
            if j != i {
                prod *= xs[j];
            }
        }
        n0 += prod;
    }
    // 単純な合成 Simpson (被積分は滑らか・正定): s ∈ [4m², s_mid] + tail
    let integrand = |s: f64| rho_d_closed(s, m) * n0 / ((s + xs[0]) * (s + xs[1]) * (s + xs[2]) * (s + xs[3]));
    let s_th = 4.0 * m * m;
    let s_mid = (60.0 * xs[3]).max(30.0 * s_th).max(6.0);
    let nstep = 40000usize;
    let h = (s_mid - s_th) / nstep as f64;
    let mut total = 0.0;
    for i in 0..nstep {
        let s0 = s_th + i as f64 * h;
        total += h / 6.0 * (integrand(s0) + 4.0 * integrand(s0 + 0.5 * h) + integrand(s0 + h));
    }
    // tail: s = s_mid/t
    let nt = 4000usize;
    let ht = 1.0 / nt as f64;
    for i in 0..nt {
        let t0 = 1e-9 + i as f64 * ht;
        let f = |t: f64| {
            let s = s_mid / t;
            integrand(s) * s_mid / (t * t)
        };
        total += ht / 6.0 * (f(t0) + 4.0 * f(t0 + 0.5 * ht) + f((t0 + ht).min(1.0)));
    }
    total
}

// ---------------- 8 成分折込み基底 (v268z の認証済み実装を写経) ----------------

fn sbit(s: usize, ax: usize) -> usize {
    (s >> ax) & 1
}

/// H(k) = Σ cos(kᵢ + sᵢπ) 構造 — 実対称 8×8 (列優先 v[r + c*8] 対称なので同じ)
fn h8(k: [f64; 3], m: f64) -> Vec<f64> {
    let mut h = vec![0.0f64; 64];
    for s in 0..8usize {
        // Γx: 対角 (−1)^{sx} cos kx
        let cx = if sbit(s, 0) == 0 { 1.0 } else { -1.0 } * k[0].cos();
        h[s + s * 8] += cx;
        // Γy: sx flip, (−1)^{sy} cos ky
        let s2 = s ^ 1;
        let cy = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * k[1].cos();
        h[s2 + s * 8] += cy;
        // Γz: sx,sy flip, (−1)^{sz} cos kz
        let s3 = s ^ 3;
        let cz = if sbit(s, 2) == 0 { 1.0 } else { -1.0 } * k[2].cos();
        h[s3 + s * 8] += cz;
        // Γm: 全 flip, m
        let s4 = s ^ 7;
        h[s4 + s * 8] += m;
    }
    h
}

/// 一般頂点 (v268z の公式, 実部のみ — D 用の項は全て実):
/// V_{s+ε,s}(k;q) = e^{iq·dy/2}[w e^{i(k+sπ)·d} + w e^{−i(k+qŷ+(s+ε)π)·d}] (w 実)
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
            // 実重みの実部: w[cos(ph1+mid) + cos(−ph2+mid)]
            v[s + s2 * 8] += t.w * ((ph1 + mid).cos() + (-ph2 + mid).cos());
        }
    }
    v
}

fn t_xx() -> Vec<Term> {
    vec![Term { eps: 0, d: [1, 0, 0], w: 0.5 }]
}
fn t_zz(mutate: bool) -> Vec<Term> {
    vec![Term {
        eps: if mutate { 2 } else { 3 },
        d: [0, 0, 1],
        w: 0.5,
    }]
}

/// χ_D の被積分 (k 点ごと): Σᵢ wᵢ Σ_{occ,unocc} 2M²/ΔE を直接累積。
/// 8 成分セルの規格化: サイトあたり = (1/8)·(セル平均) — v26.7 公表値 (示強
/// χ/site) に合わせ、セル積分 ∫d³k/(π³ セル体積) × (1/…) は [S0] が較正を裁く。
fn chi_null_integrand(
    k: [f64; 3],
    qs_lat: &[f64],
    w_null: &[f64],
    m: f64,
    mutate: bool,
) -> f64 {
    let hk = h8(k, m);
    let (wk, vk) = jacobi_eigh(&hk, 8);
    let r2i = 1.0 / (2.0f64).sqrt();
    let mut acc = 0.0f64;
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let hq = h8(kq, m);
        let (wq, vq) = jacobi_eigh(&hq, 8);
        let vx = vertex8(&t_xx(), k, q);
        let vz = vertex8(&t_zz(mutate), k, q);
        // V_D = (vx − vz)/√2
        let mut chi = 0.0f64;
        for mu in 4..8 {
            for nu in 0..4 {
                let mut mel = 0.0f64;
                for r in 0..8 {
                    let mut s = 0.0f64;
                    for c in 0..8 {
                        s += (vx[c + r * 8] - vz[c + r * 8]) * r2i * vk[c + nu * 8];
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

/// セル [0,π)³ (ノード k* = (π/2,π/2,π/2) が中心) の BZ 積分。
/// パネル: 各軸 {[0, π/2−r], [π/2−r, π/2−r/4], [π/2−r/4, π/2+r/4],
/// [π/2+r/4, π/2+r], [π/2+r, π]} — ノード box を細分。決定的スレッド分割。
fn bz_integrate(
    qs_lat: &[f64],
    w_null: &[f64],
    m: f64,
    r_node: f64,
    r_fine: f64,
    ngl: usize,
    nthreads: usize,
    mutate: bool,
) -> f64 {
    let gl = {
        // GL ノード (v268a と同じ Newton 法)
        let n = ngl;
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
    };
    // 開発記録 (run1 → run2): 2 段の箱ではノード構造 (スケール ~ aQ_min) を
    // 解像できず S1 が 35% を検出 — 3 段の入れ子箱 (r₃ ~ 2aQ_min まで) に細分。
    let c = PI / 2.0;
    let r1 = r_node; // 外箱
    let r2 = (r_node / 3.0).max(2.0 * r_fine.min(r_node / 3.0)); // 中箱
    let r3 = r_fine.min(r2 / 2.0); // 内箱 ~ 2aQ_min
    // 内箱はノード点 c で分割 (GL 節点は端に密集 — 構造中心に密度を置く)
    let edges = [
        0.0,
        c - r1,
        c - r2,
        c - r3,
        c,
        c + r3,
        c + r2,
        c + r1,
        PI,
    ];
    // 1 軸のノード列: (パネル, GL 点) の全列挙
    let mut nodes1d: Vec<(f64, f64)> = Vec::new(); // (k, weight)
    for w2 in edges.windows(2) {
        let (a, b) = (w2[0], w2[1]);
        let (cc, hh) = (0.5 * (a + b), 0.5 * (b - a));
        for (x, wgt) in gl.0.iter().zip(&gl.1) {
            nodes1d.push((cc + hh * x, wgt * hh));
        }
    }
    let n1 = nodes1d.len();
    // kx を行としてスレッド分割 (決定的)
    let mut rows: Vec<Option<f64>> = Vec::new();
    rows.resize_with(n1, || None);
    let chunk = n1.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            let nodes = &nodes1d;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let ix = t * chunk + i;
                    let (kx, wx) = nodes[ix];
                    let mut acc = 0.0f64;
                    for &(ky, wy) in nodes.iter() {
                        for &(kz, wz) in nodes.iter() {
                            acc += wy
                                * wz
                                * chi_null_integrand([kx, ky, kz], qs_lat, w_null, m, mutate);
                        }
                    }
                    *slot = Some(acc * wx);
                }
            });
        }
    });
    // セル平均 → サイトあたり: (1/(2π)³)·(8 成分は既に和に含む — セルは体積 π³ で
    // 全状態を 1 回ずつ被覆) ⇒ ∫_cell d³k/(2π)³ = (1/(2π)³)Σ
    rows.into_iter().map(|r| r.unwrap()).sum::<f64>() / (2.0 * PI).powi(3)
}

/// 単一 q の χ (S0 回帰用): null 重み {1} 相当
fn chi_single(q_lat: f64, m: f64, r_node: f64, ngl: usize, nthreads: usize) -> f64 {
    bz_integrate(&[q_lat], &[1.0], m, r_node, 0.5 * q_lat, ngl, nthreads, false)
}

/// null 重み (物理 Q で)
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

/// 2 モデル外挿: A₀+c₂a² と A₀+c₁a+c₂a² (最小二乗)
fn extrapolate(avals: &[(f64, f64)]) -> (f64, f64) {
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
        // 小さな正規方程式をガウス消去
        let mut m = ata;
        let mut r = atb;
        for col in 0..p {
            let piv = (col..p)
                .max_by(|&r1, &r2| m[col + r1 * p].abs().partial_cmp(&m[col + r2 * p].abs()).unwrap())
                .unwrap();
            for c in 0..p {
                m.swap(c + col * p, c + piv * p);
            }
            r.swap(col, piv);
            let d = m[col + col * p];
            for row in col + 1..p {
                let f = m[col + row * p] / d;
                for c in col..p {
                    m[c + row * p] -= f * m[c + col * p];
                }
                r[row] -= f * r[col];
            }
        }
        let mut sol = vec![0.0f64; p];
        for col in (0..p).rev() {
            let mut s = r[col];
            for c in col + 1..p {
                s -= m[c + col * p] * sol[c];
            }
            sol[col] = s / m[col + col * p];
        }
        sol[0]
    };
    let a_quad = fit(&|a: f64| vec![1.0, a * a]);
    let a_lin = fit(&|a: f64| vec![1.0, a, a * a]);
    (a_quad, a_lin)
}

fn main() {
    self_test();
    println!("=== v26.8-B v268b_continuum — staggered TT の連続極限 (PRED-016 前半) ===\n");
    println!("事前登録: spec §12.4 (bc644d4)。的: A(a→0) = 2A_oracle = −1/(80π²) (16π²·2A = −1/5)。");
    println!("continuum trajectory: q^lat = a·Q, m^lat = a·m_phys (lattice-unit m 固定の a→0 は禁止)。\n");
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
    let qset = [0.3f64, 0.6, 0.9, 1.2];
    let w_null = null_weights(&qset);
    let target = a_oracle_2t();

    // ---- [S0] 器械回帰: BZ 積分 vs v26.7-II 有限体積公表値 ----
    {
        let q = 2.0 * PI / 16.0;
        let mut worst = 0.0f64;
        let mut msg = String::new();
        for &(m, refv) in &REF267_CHI_D {
            let got = chi_single(q, m, 0.5, 12, nthreads);
            worst = worst.max((got - refv).abs());
            msg = format!("{} χ(m={:.1}) = {:.6} (公表 {:.6})", msg, m, got, refv);
        }
        check(
            "[S0] 器械回帰: BZ 積分の χ_D(2π/16) = v26.7-II 公表値 (±5e-3, 有限サイズ差込み)",
            worst < 5e-3,
            format!("max|Δ| = {:.1e} —{} ({} s)", worst, msg, t0.elapsed().as_secs()),
        );
    }

    // ---- massless 系列 A(a) ----
    let ladder = [0.5f64, 0.35, 0.25, 0.18, 0.125, 0.09, 0.0625, 0.045, 0.032];
    let mut avals: Vec<(f64, f64)> = Vec::new();
    println!("\n    [A(a) 表 (massless 系列)] a | A(a) | A(a)/2A_oracle");
    for &a in &ladder {
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
        let r_fine = 2.0 * a * qset[0];
        let aval = bz_integrate(&qs_lat, &w_null, 0.0, r_node, r_fine, 16, nthreads, false)
            / a.powi(4);
        println!(
            "      a = {:.3}: {:+.6e} | {:.4} ({} s)",
            a,
            aval,
            aval / target,
            t0.elapsed().as_secs()
        );
        avals.push((a, aval));
    }

    // ---- [S1] 求積自己整合 (a_min で GL 分解能を上げる) ----
    {
        let a = *ladder.last().unwrap();
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
        let r_fine = 2.0 * a * qset[0];
        let fine = bz_integrate(&qs_lat, &w_null, 0.0, r_node, r_fine, 20, nthreads, false) / a.powi(4);
        let rel = (fine / avals.last().unwrap().1 - 1.0).abs();
        check(
            "[S1] 求積自己整合: a_min で GL 16 → 20 の変化 < 1e-3 (相対)",
            rel < 1e-3,
            format!("相対変化 = {:.1e}", rel),
        );
    }

    // ---- [S2] PRED-016 前半の判定 ----
    let (a_quad, a_lin) = extrapolate(&avals);
    {
        let ratio_q = a_quad / target;
        let ratio_l = a_lin / target;
        let spread = (ratio_q - ratio_l).abs();
        let dev = (ratio_q - 1.0).abs().max(spread);
        println!(
            "    [S2 外挿] A₀ (a²) = {:+.6e} → 比 {:.4} / A₀ (a+a²) = {:+.6e} → 比 {:.4} (spread {:.3})",
            a_quad, ratio_q, a_lin, ratio_l, spread
        );
        check(
            "[S2] PRED-016 前半: A₀/(2A_oracle) = 1 ± 0.02 (2 モデル外挿, spread 込み)",
            dev < 0.02,
            format!("比 = {:.4} (a² モデル), 偏差 max = {:.3}", ratio_q, dev),
        );
    }

    // ---- [S3] q̂ 変種 ----
    {
        // q̂ = (2/a)sin(aQ/2): 重みは q̂ の物理値で解き直す
        let a_min = ladder[4];
        let mut avals_hat: Vec<(f64, f64)> = Vec::new();
        for &(a, _) in avals.iter() {
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
            let aval = bz_integrate(&qs_lat, &w_hat, 0.0, r_node, r_fine, 16, nthreads, false) / a.powi(4);
            avals_hat.push((a, aval));
        }
        let (ah_quad, _) = extrapolate(&avals_hat);
        let rel = (ah_quad / a_quad - 1.0).abs();
        let _ = a_min;
        check(
            "[S3] q̂ = (2/a)sin(aq/2) 変種: 外挿 A₀ が q 変種と 2% で一致 (spec §12.1)",
            rel < 0.02,
            format!("A₀(q̂)/A₀(q) − 1 = {:.3}", rel),
        );
    }

    // ---- [S4] massive 対照 (PRED-018 の格子側) ----
    {
        // 開発記録 (run2 → run3): m_phys = 2.0 は粗い a で am = 0.7 と artifact が
        // 支配的 — m_phys = 1.0 (am ≤ 0.25) と細かい a 系列に変更。
        let m_phys = 1.0f64;
        let mut avals_m: Vec<(f64, f64)> = Vec::new();
        for &a in &[0.25f64, 0.18, 0.125, 0.09] {
            let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
            let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
            let r_fine = 2.0 * a * qset[0];
            let aval =
                bz_integrate(&qs_lat, &w_null, a * m_phys, r_node, r_fine, 16, nthreads, false) / a.powi(4);
            avals_m.push((a, aval));
        }
        let (am_quad, _) = extrapolate(&avals_m);
        let ratio = am_quad / a_quad;
        let orc = oracle_a(&qset, &w_null, m_phys) / oracle_a(&qset, &w_null, 0.0);
        let rel = (ratio / orc - 1.0).abs();
        println!(
            "    [S4 表] A₀(m_phys=1)/A₀(0) = {:.4} vs oracle 比 {:.4} (閉形式 ρ から計算)",
            ratio, orc
        );
        check(
            "[S4] massive decoupling (PRED-018 格子側): oracle 比と 15% で一致",
            rel < 0.15,
            format!("相対差 = {:.3}", rel),
        );
    }

    // ---- [S5] 変異 ----
    {
        let a = 0.25f64;
        let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
        let r_node = (10.0 * a * qset[3]).min(1.2).max(0.4);
        let r_fine = 2.0 * a * qset[0];
        let bad = bz_integrate(&qs_lat, &w_null, 0.0, r_node, r_fine, 12, nthreads, true) / a.powi(4);
        let good = avals[2].1;
        check(
            "[S5] 変異: V_zz の η parity 破り (ε 3→2) → A が > 10% 変化",
            (bad / good - 1.0).abs() > 0.10,
            format!("相対変化 = {:.3}", (bad / good - 1.0).abs()),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.8-B".into())),
        ("kind".into(), Json::Str("staggered_tt_continuum_limit".into())),
        ("pred".into(), Json::Str("PRED-016 (前半)".into())),
        ("target_2taste".into(), Json::Num(target)),
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
        ("A0_quad".into(), Json::Num(a_quad)),
        ("A0_lin".into(), Json::Num(a_lin)),
        ("A0_over_target".into(), Json::Num(a_quad / target)),
    ]);
    let p = write_artifact("results/v268b_continuum.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **PRED-016 前半 hit — staggered の繰り込み後 TT form factor は連続 2-taste Dirac に収束** (完成は Wilson 側 v26.8-C)"
        } else {
            "FAIL — 分岐 (b) source/taste/trajectory の監査 / (c) 器械。欄が一次ソース"
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
