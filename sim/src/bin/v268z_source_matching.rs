//! v26.8-0b (Gate 1) v268z_source_matching — BOND-A 頂点の spin⊗taste matching (PRED-015)
//!
//! 事前登録: spec §12.2-3 (コミット bc644d4)。中心命題: **BOND-A の格子 source が
//! 連続極限で 2-taste Dirac の Belinfante stress tensor に流れるか**。失敗 (D) =
//! 経路 B は strain response — gravitational interpretation を終了する (Gate 1)。
//! X は §2 の未認証転写 — X のみの失敗は素朴転写の棄却 + 再設計 (非対称判定, 凍結済み)。
//!
//! 器械: 8 成分折込み基底 (シフト s ∈ {0,π}³, 縮約 BZ の k)。
//!   H(k) = cos kx·Γx + cos ky·Γy + cos kz·Γz + m·Γm,
//!   Γx = σ₃⊗I⊗I, Γy = σ₁⊗σ₃⊗I, Γz = σ₁⊗σ₁⊗σ₃, Γm = σ₁⊗σ₁⊗σ₁ (全て反交換)。
//!   Clifford 像 (単項 16 個の張る空間) = spin⊗I_taste の全体。taste 代数 = 可換子環
//!   ≅ M₂ (2 tastes)。taste-singlet 射影 = Clifford 単項への Frobenius 射影
//!   (I_taste⊗tr_taste/2 と等価)。ノード k* = (π/2,π/2,π/2)、α_i := −Γ_i, β := Γm。
//!
//! 一般頂点公式 (位置空間定義から厳密): 項 (parity ε, 変位 d, 重み w) の演算子
//!   O = Σ_x (−1)^{ε·x} e^{iq·(x+d/2)}[w c†(x)c(x+d) + w̄ c†(x+d)c(x)]
//! の (k → k+qŷ) 頂点行列は
//!   V_{s+ε,s}(k;q) = e^{iq·d_y/2}[w e^{i(k+sπ)·d} + w̄ e^{−i(k+q+(s+ε)π)·d}]
//!
//! 頂点:
//!   T_xx: (ε=0, d=x̂, w=1/2) → cos(kx+sxπ)·δ = cos kx·Γx (厳密 — Γx は Clifford 生成元
//!     そのもの ⇒ **対角頂点の taste-singlet 性は構成的定理**)
//!   T_zz: (ε=(1,1,0), d=ẑ, w=1/2) → cos kz·Γz / T_yy: (ε=(1,0,0), d=ŷ, w=1/2, 中点位相)
//!   X naive (BOND-A off-diagonal 則 = spec §2 の暫定転写):
//!     ½[(ε=0, d=ẑ, w=1/2) + (ε=(1,1,0), d=x̂, w=1/2)] — η_a 位相 × ĵ ホップ
//!   X split (本版で構成する修正転写): 4 隅 point-split
//!     piece1: ε=0, d = σx̂+2ρẑ, w = σρ/16 → sin(kx+sxπ)·sin 2kz = Γx sin kx sin 2kz
//!     piece2: ε=(1,1,0), d = σẑ+2ρx̂, w = −σρ/16 → Γz sin kz sin 2kx
//!     (符号交代 4 隅和が cos 項を消し sin·sin 積だけを残す — 厳密に taste-singlet)
//!
//! 検査 (凍結):
//!  [S0] Γ 代数 (反交換・Γ²=I) 厳密 + H(k)² = (Σcos²+m²)·I (数点 1e-14)
//!  [S1] Clifford 単項 16 個の Frobenius 直交性 (厳密) + taste 可換子環の次元 = 4 (M₂)
//!  [S2] 頂点公式の格子照合: N=8 周期格子の位置空間演算子を平面波 8 成分基底で挟み、
//!       解析 8×8 と一致 (1e-12) — T_xx/T_yy/T_zz/X_naive/X_split の全て
//!  [S3] **D matching**: taste 残差 r_D = 0 (< 1e-12, 全 sample — 構成的) かつ
//!       shape 残差 m_D(ε) が ladder {0.4,0.2,0.1,0.05} で単調減少 (O(ε²)) +
//!       最細 < 5e-3 (許容差 = improvement 級)
//!  [S4] **X naive**: r(ε) の記録 — 予想 O(1) (→ 0 しない) = 素朴転写の棄却 branch
//!  [S5] **X split**: r < 1e-12 (厳密 singlet) + m(ε) 単調減少 + 最細 < 5e-3
//!  [S6] Z_D(ε), Z_X(ε) の記録 (D と X に別々の Z — 立方既約表現の分裂, spec §12.2-4)
//!  [S7] PRED-015 採点: D 通過 ∧ 認証済み X 転写の存在 → hit (naive 棄却込み)
//!  [S8] 変異: X_split の corner 符号 1 個反転 → r が O(1) (> 0.1) に跳ねる
//!
//! 事前登録分岐 (spec §12.2-3): (a) S3 PASS → D は TreeLevelMatchedTTSource — Gate 1
//!   開通 (X は S5 の結果で naive棄却+修正認証を記録) / (b) S3 FAIL → **経路 B の
//!   gravitational interpretation 終了** (one-loop へ進まない) / (c) S0–S2 FAIL → 器械。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

// ---------------- 複素 8×8 ----------------

type M8 = Vec<(f64, f64)>; // 8×8, m[c + r*8]

fn mzero() -> M8 {
    vec![(0.0, 0.0); 64]
}

fn madd(a: &M8, b: &M8) -> M8 {
    a.iter().zip(b).map(|(x, y)| (x.0 + y.0, x.1 + y.1)).collect()
}

fn mscale(c: (f64, f64), a: &M8) -> M8 {
    a.iter().map(|x| (c.0 * x.0 - c.1 * x.1, c.0 * x.1 + c.1 * x.0)).collect()
}

fn mmul(a: &M8, b: &M8) -> M8 {
    let mut o = mzero();
    for r in 0..8 {
        for k in 0..8 {
            let av = a[k + r * 8];
            if av.0 == 0.0 && av.1 == 0.0 {
                continue;
            }
            for c in 0..8 {
                let bv = b[c + k * 8];
                o[c + r * 8].0 += av.0 * bv.0 - av.1 * bv.1;
                o[c + r * 8].1 += av.0 * bv.1 + av.1 * bv.0;
            }
        }
    }
    o
}

/// Frobenius 内積 ⟨A, B⟩ = Σ conj(A)·B
fn mdot(a: &M8, b: &M8) -> (f64, f64) {
    let (mut re, mut im) = (0.0, 0.0);
    for i in 0..64 {
        re += a[i].0 * b[i].0 + a[i].1 * b[i].1;
        im += a[i].0 * b[i].1 - a[i].1 * b[i].0;
    }
    (re, im)
}

fn mnorm(a: &M8) -> f64 {
    mdot(a, a).0.sqrt()
}

// ---------------- Γ 行列 (パウリのテンソル積, 実) ----------------

/// s = (sx, sy, sz) を 3 bit で: index = sx + 2sy + 4sz
fn sbit(s: usize, ax: usize) -> usize {
    (s >> ax) & 1
}

/// Γ を (flip mask, サイトごとの符号関数) で表現して構築
fn build_gamma(flip: usize, sign: impl Fn(usize) -> f64) -> M8 {
    let mut g = mzero();
    for s in 0..8 {
        let s2 = s ^ flip;
        g[s + s2 * 8] = (sign(s), 0.0); // 列 s → 行 s2
    }
    g
}

fn gammas() -> [M8; 4] {
    // Γx = σ₃(sx): 対角, 符号 (−1)^{sx}
    let gx = build_gamma(0, |s| if sbit(s, 0) == 0 { 1.0 } else { -1.0 });
    // Γy = σ₁(sx)⊗σ₃(sy): sx flip, 符号 (−1)^{sy}
    let gy = build_gamma(1, |s| if sbit(s, 1) == 0 { 1.0 } else { -1.0 });
    // Γz = σ₁(sx)⊗σ₁(sy)⊗σ₃(sz): sx,sy flip, 符号 (−1)^{sz}
    let gz = build_gamma(3, |s| if sbit(s, 2) == 0 { 1.0 } else { -1.0 });
    // Γm = σ₁⊗σ₁⊗σ₁: 全 flip, 符号 +1
    let gm = build_gamma(7, |_| 1.0);
    [gx, gy, gz, gm]
}

/// Clifford 単項 16 個 (I, Γi, ΓiΓj, ΓiΓjΓk, Γ5) — Frobenius 直交基底
fn clifford_monomials(g: &[M8; 4]) -> Vec<M8> {
    let mut out = Vec::new();
    for mask in 0..16u32 {
        let mut m = {
            let mut id = mzero();
            for i in 0..8 {
                id[i + i * 8] = (1.0, 0.0);
            }
            id
        };
        for i in 0..4 {
            if (mask >> i) & 1 == 1 {
                m = mmul(&m, &g[i]);
            }
        }
        out.push(m);
    }
    out
}

/// taste-singlet 射影: Clifford 単項への Frobenius 射影 (単項は直交, ‖C‖² = 8)
fn singlet_part(mono: &[M8], v: &M8) -> M8 {
    let mut o = mzero();
    for c in mono {
        let ip = mdot(c, v);
        o = madd(&o, &mscale((ip.0 / 8.0, ip.1 / 8.0), c));
    }
    o
}

// ---------------- 一般頂点 (解析 8×8) ----------------

/// 項: parity ε (bit mask), 変位 d, 複素重み w
struct Term {
    eps: usize,
    d: [i32; 3],
    w: (f64, f64),
}

/// V_{s+ε,s}(k; qŷ) = e^{iq·dy/2}[w e^{i(k+sπ)·d} + w̄ e^{−i(k+qŷ+(s+ε)π)·d}]
fn vertex_analytic(terms: &[Term], k: [f64; 3], q: f64) -> M8 {
    let mut v = mzero();
    for t in terms {
        for s in 0..8usize {
            let s2 = s ^ t.eps;
            let mut ph1 = 0.0f64; // (k+sπ)·d
            let mut ph2 = 0.0f64; // (k+qŷ+(s+ε)π)·d
            for ax in 0..3 {
                let ka = k[ax] + PI * sbit(s, ax) as f64;
                let ka2 = k[ax]
                    + if ax == 1 { q } else { 0.0 }
                    + PI * sbit(s2, ax) as f64;
                ph1 += ka * t.d[ax] as f64;
                ph2 += ka2 * t.d[ax] as f64;
            }
            let mid = 0.5 * q * t.d[1] as f64; // 中点位相 e^{iq·dy/2}
            // w e^{i(ph1+mid)} + w̄ e^{i(−ph2+mid)}
            let e1 = (ph1 + mid).cos();
            let f1 = (ph1 + mid).sin();
            let e2 = (-ph2 + mid).cos();
            let f2 = (-ph2 + mid).sin();
            let re = t.w.0 * e1 - t.w.1 * f1 + t.w.0 * e2 + t.w.1 * f2;
            let im = t.w.0 * f1 + t.w.1 * e1 + t.w.0 * f2 - t.w.1 * e2;
            v[s + s2 * 8].0 += re;
            v[s + s2 * 8].1 += im;
        }
    }
    v
}

fn t_xx() -> Vec<Term> {
    vec![Term { eps: 0, d: [1, 0, 0], w: (0.5, 0.0) }]
}
fn t_yy() -> Vec<Term> {
    vec![Term { eps: 1, d: [0, 1, 0], w: (0.5, 0.0) }]
}
fn t_zz() -> Vec<Term> {
    vec![Term { eps: 3, d: [0, 0, 1], w: (0.5, 0.0) }]
}
/// X naive: BOND-A off-diagonal 則 (spec §2) — η_x(≡1) 位相 z ホップ + η_z 位相 x ホップ
fn t_x_naive() -> Vec<Term> {
    vec![
        Term { eps: 0, d: [0, 0, 1], w: (0.25, 0.0) },
        Term { eps: 3, d: [1, 0, 0], w: (0.25, 0.0) },
    ]
}
/// X split: 4 隅 point-split (σρ 符号交代 — sin·sin 積のみ残り厳密 taste-singlet)。
/// 開発記録 (run1 → run2): piece1 の重みの相対符号を +σρ/16 と誤実装 — 2 片が
/// (−αx pz + αz px) の反対称結合 (回転生成子) になり、shape ゲート m ≈ 2.0 /
/// Z ≈ −0.25 が検出した (taste 残差 r は両片とも厳密 0 のため不感)。両片とも
/// w = −σρ/16 が対称 stress (αx pz + αz px ≈ 2×target, Z ≈ 2) を与える正しい符号。
fn t_x_split() -> Vec<Term> {
    let mut v = Vec::new();
    for sg in [1i32, -1] {
        for rh in [1i32, -1] {
            let c = (sg * rh) as f64 / 16.0;
            v.push(Term { eps: 0, d: [sg, 0, 2 * rh], w: (-c, 0.0) });
            v.push(Term { eps: 3, d: [2 * rh, 0, sg], w: (-c, 0.0) });
        }
    }
    v
}

/// D チャネル (xx − zz)/√2
fn vertex_d(k: [f64; 3], q: f64) -> M8 {
    let vx = vertex_analytic(&t_xx(), k, q);
    let vz = vertex_analytic(&t_zz(), k, q);
    let r2 = 1.0 / (2.0f64).sqrt();
    madd(&mscale((r2, 0.0), &vx), &mscale((-r2, 0.0), &vz))
}

// ---------------- 連続 Dirac 頂点 (折込み基底, α = −Γ) ----------------

/// T_ij^cont = ¼(α_i (2p+q)_j + α_j (2p+q)_i), q = qŷ
fn vertex_cont(g: &[M8; 4], ij: (usize, usize), p: [f64; 3], q: f64) -> M8 {
    let tp = |ax: usize| 2.0 * p[ax] + if ax == 1 { q } else { 0.0 };
    let ai = mscale((-0.25 * tp(ij.1), 0.0), &g[ij.0]);
    let aj = mscale((-0.25 * tp(ij.0), 0.0), &g[ij.1]);
    madd(&ai, &aj)
}

fn vertex_cont_d(g: &[M8; 4], p: [f64; 3], q: f64) -> M8 {
    let vx = vertex_cont(g, (0, 0), p, q);
    let vz = vertex_cont(g, (2, 2), p, q);
    let r2 = 1.0 / (2.0f64).sqrt();
    madd(&mscale((r2, 0.0), &vx), &mscale((-r2, 0.0), &vz))
}

// ---------------- N=8 周期格子の平面波照合 (S2) ----------------

/// 位置空間演算子を 8 成分平面波で挟む: ⟨k+qŷ, s'|O|k, s⟩
/// O = Σ_x (−1)^{ε·x} e^{iq(y+dy/2)} [w c†(x)c(x+d) + w̄ c†(x+d)c(x)]
fn sandwich_lattice(n: usize, terms: &[Term], kbase: [f64; 3], qy: f64) -> M8 {
    let ns = n * n * n;
    let mut v = mzero();
    // 平面波: ψ_{k,s}(x) = e^{i(k+sπ)·x}/√V
    let wavevec = |s: usize| -> [f64; 3] {
        [
            kbase[0] + PI * sbit(s, 0) as f64,
            kbase[1] + PI * sbit(s, 1) as f64,
            kbase[2] + PI * sbit(s, 2) as f64,
        ]
    };
    for t in terms {
        for s in 0..8usize {
            let s2 = s ^ t.eps;
            let kv = wavevec(s);
            let kv2 = {
                let mut a = wavevec(s2);
                a[1] += qy;
                a
            };
            let (mut re, mut im) = (0.0f64, 0.0f64);
            for x in 0..n {
                for y in 0..n {
                    for z in 0..n {
                        let xi = [x as f64, y as f64, z as f64];
                        let par = ((t.eps & 1) * x + ((t.eps >> 1) & 1) * y + ((t.eps >> 2) & 1) * z) % 2;
                        let pf = if par == 0 { 1.0 } else { -1.0 };
                        // 周期格子 + 格子可換な k なので wrap は平面波位相に自動吸収
                        let qph = qy * (xi[1] + 0.5 * t.d[1] as f64);
                        // 項 1: w·conj(ψ_{k+q,s'}(x))·ψ_{k,s}(x+d)
                        let ph1 = -(kv2[0] * xi[0] + kv2[1] * xi[1] + kv2[2] * xi[2])
                            + kv[0] * (xi[0] + t.d[0] as f64)
                            + kv[1] * (xi[1] + t.d[1] as f64)
                            + kv[2] * (xi[2] + t.d[2] as f64)
                            + qph;
                        let (c1, s1v) = (ph1.cos(), ph1.sin());
                        re += pf * (t.w.0 * c1 - t.w.1 * s1v);
                        im += pf * (t.w.0 * s1v + t.w.1 * c1);
                        // 項 2: w̄·conj(ψ_{k+q,s'}(x+d))·ψ_{k,s}(x)
                        let ph2 = -(kv2[0] * (xi[0] + t.d[0] as f64)
                            + kv2[1] * (xi[1] + t.d[1] as f64)
                            + kv2[2] * (xi[2] + t.d[2] as f64))
                            + kv[0] * xi[0]
                            + kv[1] * xi[1]
                            + kv[2] * xi[2]
                            + qph;
                        let (c2, s2v) = (ph2.cos(), ph2.sin());
                        re += pf * (t.w.0 * c2 + t.w.1 * s2v);
                        im += pf * (t.w.0 * s2v - t.w.1 * c2);
                    }
                }
            }
            v[s + s2 * 8].0 += re / ns as f64;
            v[s + s2 * 8].1 += im / ns as f64;
        }
    }
    v
}

fn main() {
    self_test();
    println!("=== v26.8-0b (Gate 1) — BOND-A 頂点の spin⊗taste matching (PRED-015) ===\n");
    println!("事前登録: spec §12.2-3 (bc644d4)。D 失敗 = 経路 B の gravitational 解釈終了 /");
    println!("X のみ失敗 = 素朴転写の棄却 + 再設計 (X は未認証 — 非対称判定は凍結済み)。\n");
    let t0 = std::time::Instant::now();
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
    let g = gammas();

    // ---- [S0] Γ 代数 + H(k)² ----
    {
        let mut worst = 0.0f64;
        for i in 0..4 {
            for j in 0..4 {
                let ab = madd(&mmul(&g[i], &g[j]), &mmul(&g[j], &g[i]));
                for r in 0..8 {
                    for c in 0..8 {
                        let want = if i == j && r == c { 2.0 } else { 0.0 };
                        worst = worst
                            .max((ab[c + r * 8].0 - want).abs())
                            .max(ab[c + r * 8].1.abs());
                    }
                }
            }
        }
        // H(k)² = (Σcos² + m²)I の数点検査
        let mut worst_h = 0.0f64;
        for (k, m) in [
            ([0.3f64, 1.1, -0.7], 0.5f64),
            ([1.2, 0.4, 2.0], 0.0),
            ([PI / 2.0, PI / 2.0, PI / 2.0], 0.25),
        ] {
            let mut h = mzero();
            for ax in 0..3 {
                h = madd(&h, &mscale((k[ax].cos(), 0.0), &g[ax]));
            }
            h = madd(&h, &mscale((m, 0.0), &g[3]));
            let h2 = mmul(&h, &h);
            let want = k[0].cos().powi(2) + k[1].cos().powi(2) + k[2].cos().powi(2) + m * m;
            for r in 0..8 {
                for c in 0..8 {
                    let w = if r == c { want } else { 0.0 };
                    worst_h = worst_h
                        .max((h2[c + r * 8].0 - w).abs())
                        .max(h2[c + r * 8].1.abs());
                }
            }
        }
        check(
            "[S0] Γ 代数 {Γa,Γb} = 2δ·I (厳密) + H(k)² = (Σcos²+m²)·I (1e-14)",
            worst < 1e-14 && worst_h < 1e-14,
            format!("max|Δ| = {:.1e} / {:.1e}", worst, worst_h),
        );
    }

    // ---- [S1] Clifford 単項の直交性 + taste 可換子環の次元 ----
    let mono = clifford_monomials(&g);
    {
        let mut worst = 0.0f64;
        for i in 0..16 {
            for j in 0..16 {
                let ip = mdot(&mono[i], &mono[j]);
                let want = if i == j { 8.0 } else { 0.0 };
                worst = worst.max((ip.0 - want).abs()).max(ip.1.abs());
            }
        }
        // taste 可換子環: twirl Q(M) = (1/16)Σ C M C⁻¹ の像の次元 (matrix unit で走査)
        // C は直交 (C² = ±I): C⁻¹ = C^T = ±C — Frobenius 射影で単項係数を測る方が単純:
        // M が可換子環 ⇔ 全単項と可換。次元 = 64 − rank(ad 写像) を数値で:
        // ここでは同値な検査として「Clifford 像 (dim 16) ⊗ 可換子環 (dim 4) = 64」を
        // twirl 像の Gram rank = 4 で確認する。
        let mut twirl_imgs: Vec<M8> = Vec::new();
        for u in 0..8 {
            for v in 0..8 {
                let mut e = mzero();
                e[v + u * 8] = (1.0, 0.0);
                let mut q = mzero();
                for c in &mono {
                    // C e C⁻¹, C⁻¹ = C^T (実直交, C² = ±I → C⁻¹ = ±C — 符号は共役で相殺)
                    let ct: M8 = {
                        let mut t = mzero();
                        for r in 0..8 {
                            for cc in 0..8 {
                                t[cc + r * 8] = c[r + cc * 8];
                            }
                        }
                        t
                    };
                    q = madd(&q, &mmul(c, &mmul(&e, &ct)));
                }
                // C C^T = I なので (1/16)ΣC e C^T が可換子環への射影 (Clifford 単項群の twirl)
                twirl_imgs.push(mscale((1.0 / 16.0, 0.0), &q));
            }
        }
        // Gram rank (閾値 1e-8)
        let nimg = twirl_imgs.len();
        let mut gram = vec![0.0f64; nimg * nimg];
        for i in 0..nimg {
            for j in 0..nimg {
                gram[j + i * nimg] = mdot(&twirl_imgs[i], &twirl_imgs[j]).0;
            }
        }
        // 簡易 rank: ピボット付きガウス消去
        let mut rank = 0usize;
        let mut gm = gram.clone();
        let mut used = vec![false; nimg];
        for _ in 0..nimg {
            let mut piv = None;
            let mut best = 1e-8;
            for r in 0..nimg {
                if !used[r] && gm[r + r * nimg].abs() > best {
                    best = gm[r + r * nimg].abs();
                    piv = Some(r);
                }
            }
            let Some(p) = piv else { break };
            used[p] = true;
            rank += 1;
            let d = gm[p + p * nimg];
            for r in 0..nimg {
                if used[r] {
                    continue;
                }
                let f = gm[p + r * nimg] / d;
                for c2 in 0..nimg {
                    gm[c2 + r * nimg] -= f * gm[c2 + p * nimg];
                }
            }
        }
        check(
            "[S1] Clifford 単項 16 個の直交性 (厳密) + taste 可換子環 dim = 4 (M₂ = 2 tastes)",
            worst < 1e-13 && rank == 4,
            format!("直交 max|Δ| = {:.1e} / twirl 像 rank = {}", worst, rank),
        );
    }

    // ---- [S2] 頂点公式の格子照合 (N=8 周期, 平面波 sandwich) ----
    {
        let n = 8usize;
        let kbase = [2.0 * PI * 1.0 / n as f64, 2.0 * PI * 2.0 / n as f64, -2.0 * PI * 1.0 / n as f64];
        let qy = 2.0 * PI * 1.0 / n as f64;
        let mut worst = 0.0f64;
        for (name, terms) in [
            ("T_xx", t_xx()),
            ("T_yy", t_yy()),
            ("T_zz", t_zz()),
            ("X_naive", t_x_naive()),
            ("X_split", t_x_split()),
        ] {
            let va = vertex_analytic(&terms, kbase, qy);
            let vl = sandwich_lattice(n, &terms, kbase, qy);
            let mut dev = 0.0f64;
            for i in 0..64 {
                dev = dev.max((va[i].0 - vl[i].0).abs()).max((va[i].1 - vl[i].1).abs());
            }
            worst = worst.max(dev);
            let _ = name;
        }
        check(
            "[S2] 解析 8×8 頂点 = N=8 周期格子の平面波 sandwich (全 5 頂点, 1e-12)",
            worst < 1e-12,
            format!("max|Δ| = {:.1e} ({} s)", worst, t0.elapsed().as_secs()),
        );
    }

    // ---- matching ladder ----
    let kstar = [PI / 2.0, PI / 2.0, PI / 2.0];
    let ladder = [0.4f64, 0.2, 0.1, 0.05];
    let pdirs: [[f64; 3]; 4] = [
        [0.577350, 0.577350, 0.577350],
        [1.0, 0.0, 0.0],
        [0.0, 0.0, 1.0],
        [0.324443, -0.757034, 0.567776],
    ];
    // (r, m, Z) を頂点ごとに測る
    let measure = |mk_lat: &dyn Fn([f64; 3], f64) -> M8,
                   mk_cont: &dyn Fn([f64; 3], f64) -> M8,
                   eps: f64|
     -> (f64, f64, f64) {
        let mut worst_r = 0.0f64;
        let mut worst_m = 0.0f64;
        let mut zsum = 0.0f64;
        let mut zn = 0usize;
        for pd in &pdirs {
            for qf in [0.0f64, 0.7] {
                let p = [eps * pd[0], eps * pd[1], eps * pd[2]];
                let q = eps * qf;
                if qf == 0.0 && pd[1] == 0.0 && pd[0] == 1.0 {
                    // p∥x̂, q=0 も有効 sample — 除外しない
                }
                let k = [kstar[0] + p[0], kstar[1] + p[1], kstar[2] + p[2]];
                let vl = mk_lat(k, q);
                let nl = mnorm(&vl);
                if nl < 1e-14 {
                    continue;
                }
                let vs = singlet_part(&mono, &vl);
                let resid = madd(&vl, &mscale((-1.0, 0.0), &vs));
                worst_r = worst_r.max(mnorm(&resid) / nl);
                let vc = mk_cont(p, q);
                let nc = mnorm(&vc);
                if nc < 1e-14 {
                    continue;
                }
                let z = mdot(&vc, &vs).0 / (nc * nc);
                let diff = madd(&vs, &mscale((-z, 0.0), &vc));
                worst_m = worst_m.max(mnorm(&diff) / nc);
                zsum += z;
                zn += 1;
            }
        }
        (worst_r, worst_m, zsum / zn as f64)
    };
    let lat_d = |k: [f64; 3], q: f64| vertex_d(k, q);
    let cont_d = |p: [f64; 3], q: f64| vertex_cont_d(&g, p, q);
    let lat_xn = |k: [f64; 3], q: f64| vertex_analytic(&t_x_naive(), k, q);
    let lat_xs = |k: [f64; 3], q: f64| vertex_analytic(&t_x_split(), k, q);
    let cont_x = |p: [f64; 3], q: f64| vertex_cont(&g, (0, 2), p, q);

    println!("\n    [matching ladder] ε | r_D | m_D | Z_D | r_Xnaive | r_Xsplit | m_Xsplit | Z_X");
    let mut rows = Vec::new();
    for &eps in &ladder {
        let (rd, md, zd) = measure(&lat_d, &cont_d, eps);
        let (rxn, _mxn, _zxn) = measure(&lat_xn, &cont_x, eps);
        let (rxs, mxs, zxs) = measure(&lat_xs, &cont_x, eps);
        println!(
            "      ε={:.2}: {:.2e} | {:.4e} | {:.4} | {:.3} | {:.2e} | {:.4e} | {:.4}",
            eps, rd, md, zd, rxn, rxs, mxs, zxs
        );
        rows.push((eps, rd, md, zd, rxn, rxs, mxs, zxs));
    }

    // ---- [S3] D matching (Gate 1 本体) ----
    {
        let worst_r = rows.iter().map(|r| r.1).fold(0.0f64, f64::max);
        let mono_dec = rows.windows(2).all(|w| w[0].2 > w[1].2);
        let finest_m = rows.last().unwrap().2;
        check(
            "[S3] D matching: r_D < 1e-12 (taste-singlet 構成的) ∧ m_D(ε) 単調減少 ∧ 最細 < 5e-3",
            worst_r < 1e-12 && mono_dec && finest_m < 5e-3,
            format!("max r_D = {:.1e}, m_D 最細 = {:.1e}", worst_r, finest_m),
        );
    }

    // ---- [S4] X naive の棄却 branch ----
    {
        let r_coarse = rows[0].4;
        let r_fine = rows.last().unwrap().4;
        println!(
            "    [S4 branch] r_Xnaive(ε): {:.3} (ε=0.4) → {:.3} (ε=0.05) — {} (素朴 off-diagonal 転写は taste-nonsinglet)",
            r_coarse,
            r_fine,
            if r_fine > 0.5 {
                "O(1) で残存 ⇒ **棄却**"
            } else {
                "予想外に減少 (要検討)"
            }
        );
        check(
            "[S4] X naive: r が O(1) で残存 (> 0.5 — 素朴転写の棄却が分解能をもって確定)",
            r_fine > 0.5,
            format!("r_Xnaive(最細) = {:.3}", r_fine),
        );
    }

    // ---- [S5] X split (修正転写) の認証 ----
    {
        let worst_r = rows.iter().map(|r| r.5).fold(0.0f64, f64::max);
        let mono_dec = rows.windows(2).all(|w| w[0].6 > w[1].6);
        let finest_m = rows.last().unwrap().6;
        check(
            "[S5] X split: r < 1e-12 (厳密 singlet) ∧ m(ε) 単調減少 ∧ 最細 < 5e-3",
            worst_r < 1e-12 && mono_dec && finest_m < 5e-3,
            format!("max r = {:.1e}, m 最細 = {:.1e}", worst_r, finest_m),
        );
    }

    // ---- [S6] Z の記録 ----
    {
        let zd = rows.last().unwrap().3;
        let zx = rows.last().unwrap().7;
        println!(
            "    [S6 記録] Z_D(最細) = {:.6} / Z_X(最細) = {:.6} — D と X は別の正規化 (立方分裂, spec §12.2-4)",
            zd, zx
        );
        check(
            "[S6] Z_D, Z_X が有限・非零 (別々の Z の記録)",
            zd.abs() > 0.1 && zx.abs() > 0.1,
            format!("Z_D = {:.4}, Z_X = {:.4}", zd, zx),
        );
    }

    // ---- [S7] PRED-015 採点 ----
    {
        let d_ok = rows.iter().map(|r| r.1).fold(0.0f64, f64::max) < 1e-12
            && rows.last().unwrap().2 < 5e-3;
        let x_ok = rows.iter().map(|r| r.5).fold(0.0f64, f64::max) < 1e-12
            && rows.last().unwrap().6 < 5e-3;
        let hit = d_ok && x_ok;
        println!(
            "    [S7 採点] PRED-015: D {} / X (修正転写) {} / naive 転写は棄却 ⇒ **{}**",
            if d_ok { "通過" } else { "失敗" },
            if x_ok { "通過" } else { "失敗" },
            if hit {
                "hit — BOND-A (対角) + point-split X は TreeLevelMatchedTTSource"
            } else {
                "miss — Gate 1 の分岐に従う"
            }
        );
        check(
            "[S7] PRED-015 採点確定 (D ∧ 認証済み X 転写の存在)",
            hit,
            format!("D = {}, X_split = {}", d_ok, x_ok),
        );
    }

    // ---- [S8] 変異: X_split の η parity 破り ----
    // 開発記録 (run1 → run2): 初版の変異 (corner 符号反転) は r に不感だった —
    // piece1 の全 Fourier 成分が Γx⊗I 構造 (1-link x 変位は常に (−1)^{sx} 重み) で
    // singlet のまま、壊れるのは shape (m) のみ。taste 構造を守るのは η parity なので、
    // 変異は piece2 の ε = (1,1,0) → (0,1,0) (η_z → 誤 parity) に変更。
    {
        let mut terms = t_x_split();
        for t in terms.iter_mut() {
            if t.eps == 3 {
                t.eps = 2; // (1,1,0) → (0,1,0)
            }
        }
        let eps = 0.1;
        let k = [kstar[0] + eps * 0.5, kstar[1] + eps * 0.5, kstar[2] + eps * 0.5];
        let v = vertex_analytic(&terms, k, eps * 0.7);
        let vs = singlet_part(&mono, &v);
        let resid = madd(&v, &mscale((-1.0, 0.0), &vs));
        let r = mnorm(&resid) / mnorm(&v);
        check(
            "[S8] 変異: X_split piece2 の η parity 破り (ε 3→2) → taste 残差 O(1) (> 0.1)",
            r > 0.1,
            format!("r = {:.3}", r),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.8-0b".into())),
        ("kind".into(), Json::Str("source_matching_certificate".into())),
        ("pred".into(), Json::Str("PRED-015".into())),
        (
            "ladder".into(),
            Json::Arr(
                rows.iter()
                    .map(|r| {
                        Json::Obj(vec![
                            ("eps".into(), Json::Num(r.0)),
                            ("r_D".into(), Json::Num(r.1)),
                            ("m_D".into(), Json::Num(r.2)),
                            ("Z_D".into(), Json::Num(r.3)),
                            ("r_X_naive".into(), Json::Num(r.4)),
                            ("r_X_split".into(), Json::Num(r.5)),
                            ("m_X_split".into(), Json::Num(r.6)),
                            ("Z_X".into(), Json::Num(r.7)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v268z_source_matching.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **Gate 1 開通 — D は TreeLevelMatchedTTSource、素朴 off-diagonal 転写は棄却、point-split X が認証** — v26.8-A (解析 oracle) へ進む資格が生じた"
        } else {
            "FAIL — 分岐 (b) D 失敗 = 経路 B の gravitational 解釈終了 / (c) 器械。欄が一次ソース"
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
