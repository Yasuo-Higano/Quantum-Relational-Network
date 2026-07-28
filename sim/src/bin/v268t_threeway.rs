//! v26.8-T v268t_threeway — 三者一致の完結: 格子 spectral density と分散再構成
//!
//! 事前登録: spec §12.4 の「三者一致」残り経路 (doc v26.8 §1 の登録済み残作業)。
//! TT (D) チャネルの同一物理量 A (q⁴ln q² 係数) を三つの独立経路で:
//!   Route I  = 直接 null 射影 (v268b/v268p — 認証済み)
//!   Route II = **格子 spectral density σ^lat(E;q) の測定 + 分散再構成 (本ユニット)**
//!   Route III = 解析 oracle (v268a — ρ_D(s) = s²/(160π²)/Dirac, 三重検証)
//! Route II は単一数値 A より遥かに強い検査 — スペクトル密度の**形**と Lorentz
//! 回復 σ(E;q) = ρ(E²−q²) (oracle では厳密に証明済み・格子では a² 破れ) を測る。
//!
//! 器械 (認証済みコアの写経 + ヒストグラム):
//!   σ^lat(E;q) := ∫_cell d³k/(2π)³ Σ_pairs |⟨μ|V_D|ν⟩|² δ(E − ΔE)
//!   ⇒ χ^lat(q) = ∫dE 2σ^lat/E (v268a と同じ規約)。物理密度は
//!   σ_phys(E) = a⁻⁴ σ^lat(aE) → 2·ρ_D(E²−q²) (2 taste, a → 0)。
//!   ヒストグラムは物理 E 幅 δE で、スペクトル全域 (E^lat ≤ 4) を被覆 —
//!   分散再構成に解析 tail は不要。再構成は v268a の教訓どおり**安定核**
//!   A_rec = ∫ρ^lat(s)·K(s)ds, K = n₀/Π(s+Qᵢ²) (素朴 Σwᵢ∫ は桁落ち死)。
//!   ρ^lat(s) は q_ref = 1.2 の σ_phys から s = E²−q_ref² で転写。
//!
//! 点評価の器械 (run1 の教訓): σ の点評価は **level-set 殻積分** で行う —
//!   σ^lat(E;q) = (1/(2π)³) ∫dΩ r²(n̂)·Σ_pairs|M|² / |∂F/∂r|,
//!   F(p) = E(k*+p+qŷ) + E(k*+p) (massless 8 成分は全対が同一 ΔE = F)。
//! 射線ごとに二分法 (F は閾値殻の外で単一上向き交差)、勾配は閉形式
//! ∂E/∂kᵢ = −sin kᵢ cos kᵢ/E。**ヒストグラムの点評価は禁じ手** — 滑らかな
//! GL 求積格子の δ-ビニングは遠ノード域で節点間隔 (~0.09) ≫ 殻幅 (~aδE)
//! となり ~15% のサンプリングノイズを生む (run1 で実測。積分 T0/T3/T4 は
//! ノイズが平均化されるため健全 — 「滑らかな格子の δ-ビニングを点評価に
//! 使うな。殻は root-solve で切れ」)。
//!
//! 検査 (凍結):
//!  [T0] 測度閉環: 各 (a, q) で per-bin χ 寄与の総和 = 直接 χ (再配列恒等式,
//!       1e-12) — ヒストグラムの索引・規格化のバグ検出器
//!  [T1a] 殻積分の角度自己整合: nθ 48 → 72 で σ の変化 < 0.3%
//!  [T1] **Lorentz 回復**: s* ∈ {0.5, 1.0, 2.0} で σ_phys(√(s*+Q²); Q) が
//!       Q = 0.3 vs 1.2 で一致 (殻積分) — a = 0.045 で max 偏差 < 1% かつ
//!       dev(0.045) < 0.5·dev(0.09) (a² 破れの縮小)
//!  [T2] **点ごと密度**: a = 0.045, Q = 0.6, E ∈ [1, 3] で
//!       σ_phys/(2ρ_D(E²−Q²)) = 1 ± 2% (oracle の絶対規格込み) かつ
//!       max 偏差が a = 0.09 より小さい
//!  [T3] **分散再構成 (Route II = I)**: A_rec(a) が同一 a の直接 A と
//!       a = 0.09 で < 2%, a = 0.045 で < 1% 一致
//!  [T4] 変異: 転写不変量を s = E² − 2q² に置換 → A_rec が > 5% 逸脱
//!
//! 事前登録分岐: (a) T0–T3 PASS → **三者一致完結 — spectral 経路でも測定器が
//!   正しい** (regulator 論文の主装置が完成。QRN の証拠ではない) /
//!   (b) T1/T2 のみ破れ → Lorentz 回復が遅い (a² 係数の公表, 器械は健全) /
//!   (c) T0/T1a FAIL → 器械。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

fn rho_d_closed(s: f64) -> f64 {
    // massless: s²/(160π²) (v268a 認証・v268b 写経と同値)
    s * s / (160.0 * PI * PI)
}

// ============ staggered 8 成分 D チャネル (v268p の認証済み写経) ============

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

fn t_xx() -> Vec<Term> {
    vec![Term { eps: 0, d: [1, 0, 0], w: 0.5 }]
}
fn t_zz() -> Vec<Term> {
    vec![Term { eps: 3, d: [0, 0, 1], w: 0.5 }]
}

/// k 点 1 つの寄与: 4 つの q それぞれの (χ 直接値, ヒストグラム流し込み)。
/// hist は物理 E ビン (幅 de_phys, 開始 0) — hist[qi*nbins + bin] += |M|²·w_bz。
/// chi_bins には χ 寄与 (2|M|²/ΔE·w_bz) を同じビンに入れる (T0 の再配列恒等式)。
struct KPointOut {
    chi: [f64; 4],
}

#[allow(clippy::too_many_arguments)]
fn accumulate_k(
    k: [f64; 3],
    w_bz: f64,
    qs_lat: &[f64],
    a: f64,
    de_phys: f64,
    nbins: usize,
    hist: &mut [f64],
    chi_bins: &mut [f64],
) -> KPointOut {
    let hk = h8(k, 0.0);
    let (wk, vk) = jacobi_eigh(&hk, 8);
    let r2i = 1.0 / (2.0f64).sqrt();
    let mut out = KPointOut { chi: [0.0; 4] };
    for (qi, &q) in qs_lat.iter().enumerate() {
        let kq = [k[0], k[1] + q, k[2]];
        let hq = h8(kq, 0.0);
        let (wq, vq) = jacobi_eigh(&hq, 8);
        let vx = vertex8(&t_xx(), k, q);
        let vz = vertex8(&t_zz(), k, q);
        let vv: Vec<f64> = (0..64).map(|i| (vx[i] - vz[i]) * r2i).collect();
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
                let de = wq[mu] - wk[nu];
                let m2 = mel * mel * w_bz;
                out.chi[qi] += 2.0 * m2 / de;
                let e_phys = de / a;
                let bin = (e_phys / de_phys) as usize;
                if bin < nbins {
                    hist[qi * nbins + bin] += m2;
                    chi_bins[qi * nbins + bin] += 2.0 * m2 / de;
                }
            }
        }
    }
    out
}

// ============ GL / 入れ子 edges / 重み (v268p の写経) ============

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

/// 安定核の分子 n₀ = Σᵢ wᵢ Πⱼ≠ᵢ xⱼ (x = Q² or 2Q²) — n₃ = n₂ = n₁ = 0 は
/// null 制約から恒等 (v268a で証明済みの構造)
fn kernel_n0(w: &[f64; 4], xs: &[f64; 4]) -> f64 {
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
    n0
}

/// 1 つの a での測定一式: 4 q の (χ 直接, 物理 E ヒストグラム σ^lat 用, χ ビン)
struct Measured {
    chi: [f64; 4],
    hist: Vec<f64>,     // [qi*nbins + bin]: Σ|M|²·w_bz (σ^lat·δE^lat 相当)
    chi_bins: Vec<f64>, // 同ビンの χ 寄与 (T0)
}

fn measure(a: f64, qset: &[f64; 4], de_phys: f64, nbins: usize, ngl: usize, nthreads: usize) -> Measured {
    let gl = gauss_legendre(ngl);
    let qs_lat: Vec<f64> = qset.iter().map(|&q| a * q).collect();
    let edges = nest_edges(PI / 2.0, 0.0, PI, a * qset[0]);
    let nodes = make_nodes(&edges, &gl);
    let n1 = nodes.len();
    type Row = (Vec<f64>, Vec<f64>, [f64; 4]);
    let mut rows: Vec<Option<Row>> = Vec::new();
    rows.resize_with(n1, || None);
    let chunk = n1.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            let nodes = &nodes;
            let qs_lat = &qs_lat;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let ix = t * chunk + i;
                    let (kx, wx) = nodes[ix];
                    let mut hist = vec![0.0f64; 4 * nbins];
                    let mut chib = vec![0.0f64; 4 * nbins];
                    let mut chi = [0.0f64; 4];
                    for &(ky, wy) in nodes.iter() {
                        for &(kz, wz) in nodes.iter() {
                            let o = accumulate_k(
                                [kx, ky, kz],
                                wx * wy * wz,
                                qs_lat,
                                a,
                                de_phys,
                                nbins,
                                &mut hist,
                                &mut chib,
                            );
                            for qi in 0..4 {
                                chi[qi] += o.chi[qi];
                            }
                        }
                    }
                    *slot = Some((hist, chib, chi));
                }
            });
        }
    });
    // 決定的マージ (行順)
    let norm = 1.0 / (2.0 * PI).powi(3);
    let mut m = Measured {
        chi: [0.0; 4],
        hist: vec![0.0; 4 * nbins],
        chi_bins: vec![0.0; 4 * nbins],
    };
    for r in rows.into_iter() {
        let (h, cb, c) = r.unwrap();
        for i in 0..4 * nbins {
            m.hist[i] += h[i] * norm;
            m.chi_bins[i] += cb[i] * norm;
        }
        for qi in 0..4 {
            m.chi[qi] += c[qi] * norm;
        }
    }
    m
}

// ============ 殻積分 (点評価の器械 — run1 の教訓) ============

/// massless 8 成分の 1 粒子エネルギー E(k) = √(Σcos²kᵢ) と radial 勾配素材
fn e8m(k: [f64; 3]) -> f64 {
    (k[0].cos().powi(2) + k[1].cos().powi(2) + k[2].cos().powi(2)).sqrt()
}

/// 対エネルギー F(p) = E(k*+p+q_lat ŷ) + E(k*+p), k* = (π/2,π/2,π/2)
fn f_pair(p: [f64; 3], q_lat: f64) -> f64 {
    let c = PI / 2.0;
    let k = [c + p[0], c + p[1], c + p[2]];
    let kq = [k[0], k[1] + q_lat, k[2]];
    e8m(kq) + e8m(k)
}

/// n̂ 方向の radial 微分 ∂F/∂r (閉形式: ∂E/∂kᵢ = −sin kᵢ cos kᵢ/E)
fn df_dr(p: [f64; 3], n: [f64; 3], q_lat: f64) -> f64 {
    let c = PI / 2.0;
    let k = [c + p[0], c + p[1], c + p[2]];
    let kq = [k[0], k[1] + q_lat, k[2]];
    let (e1, e2) = (e8m(kq).max(1e-14), e8m(k).max(1e-14));
    let mut d = 0.0f64;
    for ax in 0..3 {
        d += n[ax] * (-kq[ax].sin() * kq[ax].cos() / e1 - k[ax].sin() * k[ax].cos() / e2);
    }
    d
}

/// Σ_pairs |⟨μ|V_D|ν⟩|² (massless — 全対が同一 ΔE なので単純和)
fn m2_sum(k: [f64; 3], q_lat: f64) -> f64 {
    let hk = h8(k, 0.0);
    let (_, vk) = jacobi_eigh(&hk, 8);
    let kq = [k[0], k[1] + q_lat, k[2]];
    let hq = h8(kq, 0.0);
    let (_, vq) = jacobi_eigh(&hq, 8);
    let vx = vertex8(&t_xx(), k, q_lat);
    let vz = vertex8(&t_zz(), k, q_lat);
    let r2i = 1.0 / (2.0f64).sqrt();
    let vv: Vec<f64> = (0..64).map(|i| (vx[i] - vz[i]) * r2i).collect();
    let mut t = 0.0f64;
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
            t += mel * mel;
        }
    }
    t
}

/// 殻積分 σ^lat(E^lat; q_lat) = (1/(2π)³)∫dΩ r²·Σ|M|²/|∂F/∂r|
fn sigma_shell(e_lat: f64, q_lat: f64, nth: usize, nph: usize) -> f64 {
    let gl = gauss_legendre(nth);
    let c = PI / 2.0;
    let mut acc = 0.0f64;
    for (ct, wt) in gl.0.iter().zip(&gl.1) {
        let st = (1.0 - ct * ct).sqrt();
        for j in 0..nph {
            let ph = (j as f64 + 0.5) * 2.0 * PI / nph as f64;
            let n = [st * ph.cos(), st * ph.sin(), *ct];
            // 二分法: F(0) = 閾値 < E_target, r_hi を倍化して括る
            let mut r_hi = e_lat;
            let mut guard = 0;
            while f_pair([r_hi * n[0], r_hi * n[1], r_hi * n[2]], q_lat) <= e_lat && guard < 40 {
                r_hi *= 1.5;
                guard += 1;
            }
            let mut lo = 0.0f64;
            let mut hi = r_hi;
            for _ in 0..80 {
                let mid = 0.5 * (lo + hi);
                if f_pair([mid * n[0], mid * n[1], mid * n[2]], q_lat) < e_lat {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            let r = 0.5 * (lo + hi);
            let p = [r * n[0], r * n[1], r * n[2]];
            let k = [c + p[0], c + p[1], c + p[2]];
            let m2 = m2_sum(k, q_lat);
            let dfr = df_dr(p, n, q_lat).abs().max(1e-12);
            acc += wt * (2.0 * PI / nph as f64) * r * r * m2 / dfr;
        }
    }
    acc / (2.0 * PI).powi(3)
}

/// 物理密度 (殻積分): σ_phys(E; Q) = a⁻⁴ σ^lat(aE; aQ)
fn sigma_phys_shell(a: f64, e: f64, qq: f64, nth: usize, nph: usize) -> f64 {
    sigma_shell(a * e, a * qq, nth, nph) / a.powi(4)
}

fn main() {
    self_test();
    println!("=== v26.8-T v268t_threeway — 三者一致の完結: 格子 spectral density と分散再構成 ===\n");
    println!("Route I = 直接 null 射影 (v268b/p 認証済み) / Route II = σ^lat 測定 + 安定核");
    println!("K(s) = n₀/Π(s+Qᵢ²) 再構成 (本ユニット) / Route III = 解析 ρ_D = s²/(160π²)。");
    println!("的: σ_phys → 2ρ_D (2 taste)・A_rec = A_direct。\n");
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
    let de_phys = 0.02f64;

    // ---- 測定 (a = 0.09, 0.045) ----
    let avals = [0.09f64, 0.045];
    let mut meas: Vec<(f64, usize, Measured)> = Vec::new();
    for &a in &avals {
        // スペクトル全域: E^lat ≤ 4 (max ΔE ~ 2√(3+m²) < 3.5)
        let nbins = (4.0 / a / de_phys).ceil() as usize;
        let m = measure(a, &qset, de_phys, nbins, 14, nthreads);
        println!(
            "    [測定] a = {:.3}: nbins = {} (E_phys ≤ {:.0}), χ(q₄)^lat = {:.6e} ({} s)",
            a,
            nbins,
            4.0 / a,
            m.chi[3],
            t0.elapsed().as_secs()
        );
        meas.push((a, nbins, m));
    }

    // ---- [T0] 測度閉環 (再配列恒等式) ----
    {
        let mut worst = 0.0f64;
        for (_, nbins, m) in &meas {
            for qi in 0..4 {
                let s: f64 = m.chi_bins[qi * nbins..(qi + 1) * nbins].iter().sum();
                worst = worst.max((s / m.chi[qi] - 1.0).abs());
            }
        }
        check(
            "[T0] 測度閉環: per-bin χ 寄与の総和 = 直接 χ (全 a, q — 再配列恒等式)",
            worst < 1e-12,
            format!("max 相対差 = {:.1e}", worst),
        );
    }

    // ---- [T1a] 殻積分の角度自己整合 ----
    {
        let s0 = sigma_phys_shell(0.045, 1.5, qset[1], 48, 96);
        let s1 = sigma_phys_shell(0.045, 1.5, qset[1], 72, 144);
        let rel = (s0 / s1 - 1.0).abs();
        check(
            "[T1a] 殻積分の角度自己整合: nθ×nφ = 48×96 → 72×144 で変化 < 0.3%",
            rel < 3e-3,
            format!("相対変化 = {:.1e} ({} s)", rel, t0.elapsed().as_secs()),
        );
    }

    // ---- [T1] Lorentz 回復 (殻積分) ----
    {
        let mut devs = [0.0f64; 2];
        for (ai, &a) in avals.iter().enumerate() {
            let mut worst = 0.0f64;
            for &sstar in &[0.5f64, 1.0, 2.0] {
                let e1 = (sstar + qset[0] * qset[0]).sqrt();
                let e2 = (sstar + qset[3] * qset[3]).sqrt();
                let s1 = sigma_phys_shell(a, e1, qset[0], 48, 96);
                let s2 = sigma_phys_shell(a, e2, qset[3], 48, 96);
                worst = worst.max((s1 / s2 - 1.0).abs());
            }
            devs[ai] = worst;
            println!(
                "    [T1 表] a = {:.3}: max |σ(√(s*+Q₁²);Q₁)/σ(√(s*+Q₄²);Q₄) − 1| = {:.5}",
                a, worst
            );
        }
        check(
            "[T1] Lorentz 回復 (殻): a = 0.045 で < 1% かつ dev(0.045) < 0.5·dev(0.09)",
            devs[1] < 0.01 && devs[1] < 0.5 * devs[0],
            format!("dev = {:.5} (a=0.09) → {:.5} (a=0.045)", devs[0], devs[1]),
        );
    }

    // ---- [T2] 点ごと密度 vs oracle (殻積分) ----
    {
        let mut worsts = [0.0f64; 2];
        for (ai, &a) in avals.iter().enumerate() {
            let mut worst = 0.0f64;
            let mut e = 1.0f64;
            while e <= 3.0 {
                let sp = sigma_phys_shell(a, e, qset[1], 48, 96);
                let orc = 2.0 * rho_d_closed(e * e - qset[1] * qset[1]);
                worst = worst.max((sp / orc - 1.0).abs());
                e += 0.1;
            }
            worsts[ai] = worst;
        }
        println!(
            "    [T2 表] max |σ_phys/(2ρ_D) − 1| (E ∈ [1,3], Q = 0.6): {:.5} (a=0.09) → {:.5} (a=0.045)",
            worsts[0], worsts[1]
        );
        check(
            "[T2] 点ごと密度 (殻): a = 0.045 で σ_phys/(2ρ_D) = 1 ± 2% かつ a=0.09 より改善",
            worsts[1] < 0.02 && worsts[1] < worsts[0],
            format!("max 偏差 = {:.5}", worsts[1]),
        );
    }

    // ---- [T3] 分散再構成 (Route II = Route I) ----
    {
        let mut msg = String::new();
        let mut devs = [0.0f64; 2];
        for (ai, (a, nbins, m)) in meas.iter().enumerate() {
            // Route I: 直接 null 射影 (同一測定の χ から)
            let a_direct: f64 = (0..4).map(|qi| w_null[qi] * m.chi[qi]).sum::<f64>() / a.powi(4);
            // Route II: ρ^lat(s) を q_ref = Q₄ から転写し安定核で再結合
            // A_rec = ∫ρ(s)K(s)ds = Σ_bins σ_phys(E)·K(E²−Q₄²)·2E·δE (s = E²−Q₄², ds = 2EdE)
            let xs = [
                qset[0] * qset[0],
                qset[1] * qset[1],
                qset[2] * qset[2],
                qset[3] * qset[3],
            ];
            let n0 = kernel_n0(&w_null, &xs);
            let mut a_rec = 0.0f64;
            for b in 0..*nbins {
                let e = (b as f64 + 0.5) * de_phys;
                let s = e * e - xs[3];
                if s <= 0.0 {
                    continue;
                }
                let rho = m.hist[3 * nbins + b] / (a.powi(4) * a * de_phys);
                let kker = n0 / ((s + xs[0]) * (s + xs[1]) * (s + xs[2]) * (s + xs[3]));
                a_rec += rho * kker * 2.0 * e * de_phys;
            }
            devs[ai] = (a_rec / a_direct - 1.0).abs();
            msg = format!(
                "{} a={:.3}: A_dir={:.4e} A_rec={:.4e} (Δ {:.4})",
                msg, a, a_direct, a_rec, devs[ai]
            );
            // T4 用に a = 0.045 の変異値を後で使う
            if ai == 1 {
                // 変異: s̃ = E² − 2q² 転写 (核も 2Q² で組む — 同じ安定構造)
                let xs2 = [2.0 * xs[0], 2.0 * xs[1], 2.0 * xs[2], 2.0 * xs[3]];
                let n02 = kernel_n0(&w_null, &xs2);
                let mut a_mut = 0.0f64;
                for b in 0..*nbins {
                    let e = (b as f64 + 0.5) * de_phys;
                    let st = e * e - xs2[3];
                    if st <= 0.0 {
                        continue;
                    }
                    let rho = m.hist[3 * nbins + b] / (a.powi(4) * a * de_phys);
                    let kker = n02 / ((st + xs2[0]) * (st + xs2[1]) * (st + xs2[2]) * (st + xs2[3]));
                    a_mut += rho * kker * 2.0 * e * de_phys;
                }
                let mdev = (a_mut / a_direct - 1.0).abs();
                check(
                    "[T4] 変異: 転写不変量 s = E² − 2q² → A_rec が > 5% 逸脱",
                    mdev > 0.05,
                    format!("逸脱 = {:.3}", mdev),
                );
            }
        }
        check(
            "[T3] 分散再構成: |A_rec/A_direct − 1| < 2% (a=0.09) / < 1% (a=0.045)",
            devs[0] < 0.02 && devs[1] < 0.01,
            format!("{}", msg),
        );
    }

    // ---- artifact ----
    let a_fine = avals[1];
    let mut spectral = Vec::new();
    let mut e = 0.8f64;
    while e <= 3.0 {
        spectral.push(Json::Obj(vec![
            ("E".into(), Json::Num(e)),
            (
                "sigma_over_2rho".into(),
                Json::Num(
                    sigma_phys_shell(a_fine, e, qset[1], 48, 96)
                        / (2.0 * rho_d_closed(e * e - qset[1] * qset[1])),
                ),
            ),
        ]));
        e += 0.2;
    }
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.8-T".into())),
        ("kind".into(), Json::Str("threeway_spectral_density".into())),
        ("a_fine".into(), Json::Num(a_fine)),
        ("spectral_q06".into(), Json::Arr(spectral)),
    ]);
    let p = write_artifact("results/v268t_threeway.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **三者一致完結 — 直接 null 射影・spectral density 分散再構成・解析 oracle が同一の A に一致し、Lorentz 回復と密度の形も点ごとに oracle と一致** (測定器の証明 — QRN・創発重力の証拠ではない)"
        } else {
            "FAIL あり — 分岐 (b) Lorentz 回復が遅い (a² 係数公表) / (c) 器械。欄が一次ソース"
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
