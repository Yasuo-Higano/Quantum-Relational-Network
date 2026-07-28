//! v26.9-C v269c_spin4d — 4D spin 分離と Gate 5 総括
//!
//! 事前登録: spec §12.6 の残り項目 — 10×10 symmetric-tensor kernel・spin-0/2
//! の 4D 完全分離。v26.9 arc の最終ユニット。
//!
//! 置換則 (v26.9-0/A で 3 独立導出 [V₀₀, V_yy, V₀y — d_y = 1, 2] と整合):
//!   **中点変調頂点 V_O(k;q) = V_O^{unmod}(k + (q/2)ŷ)** — 両向き和が y 変位
//!   d_y の中点位相 e^{iq·d_y/2} を一様に生む。これで v268z の X 構成 (x,z) を
//!   (x,y)/(y,z) に回した T_xy/T_yz point-split が直ちに得られる (tree Z = 2)。
//! 全 10 source: T₀₀ = h(k+q/2ŷ) / T₀ᵢ = −(1/2)sin(2kᵢ+qδ_iy)𝟙 /
//!   BOND-A 対角 3 種 / point-split 混合 3 種 (xz, xy, yz — 後 2 者は置換則)。
//!
//! 4D spin 分離 (本ユニットの核): 殻積分 (v268t の器械) で対スペクトル密度の
//! 10×10 行列 σ_{IJ}(E; Q) = ∫dΩ r² M_I M_J† / |∂F/∂r| を測り、Lorentz 4 元
//! q_L = (E, Q ŷ) (timelike: E > Q) の Barnes–Rivers 射影子 (η = diag(+,−,−,−),
//! ProjectorND.lean の d = 4 代数の数値実装) で分解する。連続予言:
//!   **σ = ρ₂(s)·P₂ + ρ₀(s)·P₀s** (吸収部は接触項を含まない → 厳密横断的:
//!   q_L^μ σ_{μν,ρσ} = 0)。格子では O(a²) 破れ → trajectory ladder で
//!   P₁/P₀w 重みが a² で消え、10×10 が 2 つのスペクトル関数に崩壊するか。
//!
//! 検査 (凍結):
//!  [C0] 置換則の独立照合: d_y = 1 (y ボンド) と d_y = 2 (2 ホップ) の両向き和
//!       を明示導出し、置換則の値と一致 (1e-15)
//!  [C1] **Belinfante 構造の分解 (run1 の発見)**: 保存フラックス V_Fx (正準
//!       T^{yx} = α_y p_x 型) は point-split の **piece2 (Γy 構造 × x 包絡) 単独**
//!       と O(ε²) 一致し、両 piece 平均 (Belinfante) とは回転流の分だけ O(1) で
//!       違う — 射影ブロック差: vs piece2 → 0 / vs Belinfante → 定数 (両方測る)
//!  [C2] full 4D Ward の 10 列拡張: 4 行 × 10 列 × q₀ ∈ {0.3, 0.9} = 80 恒等式
//!       (k 点ごと < 1e-10, 混合正規化)
//!  [C3] **4D spin-2 分離**: D = (T_xx−T_zz)/√2 と X = T_xz は timelike
//!       q_L = (E, Qŷ) でも厳密 P₂ (q·D = 0, tr_θ D = 0 が (t,y) 面の任意 q で
//!       成立 — ProjectorND の ŷ 定理の時間方向拡張)。殻積分で
//!       (i) 偏極縮退 σ_DD/σ_XX → 1 (O(a²), 縮小比 ~4)
//!       (ii) 直交性 |σ_DX|/√(σ_DD σ_XX) → 0
//!       (iii) 正準横断性 ΔE·M₀ν = q̂·M_Fν が殻上で厳密 (< 1e-12)
//!  [C4] 変異: T_xz の Z 補正 (÷2) を落とす → σ_XX = 4σ_DD (縮退比 4 倍逸脱)
//!
//! run1 の発見 (scheme 混合の定理): (正準 T⁰ᵢ, Belinfante T_ij) を混ぜた
//! 10×10 の横断性は**破れが a 非依存 (0.616)** — 横断性は一貫した scheme
//! でのみ成立し、格子の厳密横断性は**正準 (非対称・16 成分) テンソル**が担う
//! ([C3iii])。対称 (Belinfante) 10×10 の 2 スペクトル関数への崩壊には
//! T⁰ⁱ の対称化 (= エネルギー流との平均) が必要で、x/z エネルギー流は
//! 二重変調則 V₀₀(k; px̂+qŷ) = h(k+(p/2)x̂+(q/2)ŷ) で構成可能 —
//! **v26.9-D として登録** (本ユニットの範囲外)。
//!
//! Gate 5 総括 (欄の末尾に印字 — 確立/残りを凍結解釈で明記):
//!   確立 = h₀₀/h₀ᵢ source・q₀ ≠ 0・full 4D Ward (局所カレント + 接触項,
//!   80 恒等式)・正準テンソルの厳密横断性・spin-2 チャネルの 4D 分離
//!   (P₂ 縮退 + 直交)。残り = Belinfante 対称 10×10 の完全崩壊 (v26.9-D)・
//!   temporal h の二次変分 scheme (v27.0 設計項目)。**型名
//!   FullGravitationalVacuumPolarization は保留を維持** (1/Π 禁止も維持)。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

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

type C8 = Vec<(f64, f64)>;

fn c8_from_real(a: &[f64]) -> C8 {
    a.iter().map(|&x| (x, 0.0)).collect()
}

fn c8_mul(a: &C8, b: &C8) -> C8 {
    let mut o = vec![(0.0f64, 0.0f64); 64];
    for r in 0..8 {
        for kk in 0..8 {
            let av = a[kk + r * 8];
            if av.0 == 0.0 && av.1 == 0.0 {
                continue;
            }
            for c in 0..8 {
                let bv = b[c + kk * 8];
                o[c + r * 8].0 += av.0 * bv.0 - av.1 * bv.1;
                o[c + r * 8].1 += av.0 * bv.1 + av.1 * bv.0;
            }
        }
    }
    o
}

// ---------------- 非変調頂点 (置換則の素材) ----------------

struct Term {
    eps: usize,
    d: [i32; 3],
    w: f64,
}

/// 非変調 (q = 0) の折込み頂点: V = Σ_t w·2cos((k+sπ)·d) 構造 (両向き和)
fn vertex_unmod(terms: &[Term], k: [f64; 3]) -> Vec<f64> {
    let mut v = vec![0.0f64; 64];
    for t in terms {
        for s in 0..8usize {
            let s2 = s ^ t.eps;
            let mut ph = 0.0f64;
            for ax in 0..3 {
                ph += (k[ax] + PI * sbit(s, ax) as f64) * t.d[ax] as f64;
            }
            v[s + s2 * 8] += t.w * 2.0 * ph.cos();
        }
    }
    v
}

/// point-split 混合 stress (v268z の X 構成を軸対 (a1, a2) に一般化, tree Z = 2):
/// 4 隅 σρ 交代・両片 w = −σρ/16, eps は「a1 方向の 1 歩」と「a2 方向の 1 歩」の
/// 折込みフリップ (x:0 → 対角, ただし eps は 1 歩側の軸に対応する h8 フリップ)
fn t_split_terms(a1: usize, a2: usize) -> Vec<Term> {
    // h8 のフリップ: x → eps 0 (対角), y → 1, z → 3
    let flip = [0usize, 1, 3];
    let mut v = Vec::new();
    for sg in [1i32, -1] {
        for rh in [1i32, -1] {
            let c = (sg * rh) as f64 / 16.0;
            let mut d1 = [0i32; 3];
            d1[a1] = sg;
            d1[a2] = 2 * rh;
            v.push(Term { eps: flip[a1], d: d1, w: -c });
            let mut d2 = [0i32; 3];
            d2[a2] = sg;
            d2[a1] = 2 * rh;
            v.push(Term { eps: flip[a2], d: d2, w: -c });
        }
    }
    v
}

// ---------------- 10 source (置換則で一括変調) ----------------

/// source I = 0..10: 00, 0x, 0y, 0z, xx, yy, zz, xy, xz, yz
/// 返値: 変調頂点 (実) — 全て置換則 V(k;q) = V_unmod(k + q/2 ŷ)
/// zmut: [C4] 用 — xz の Z 補正 (÷2) を落とす
fn source(i: usize, k: [f64; 3], q: f64, m: f64, zmut: bool) -> C8 {
    let km = [k[0], k[1] + 0.5 * q, k[2]];
    let re: Vec<f64> = match i {
        0 => h8(km, m),
        1 | 2 | 3 => {
            let ax = i - 1;
            let val = -0.5 * (2.0 * km[ax]).sin();
            let mut v = vec![0.0f64; 64];
            for s in 0..8usize {
                v[s + s * 8] = val;
            }
            v
        }
        // 対角 BOND-A T_ii = h8 の該当片 (係数 1·cos(k+sπ)) そのもの —
        // vertex_unmod は 2cos 規約なので w = 0.5 (開発記録: run2 は w = 0.25 と
        // 半分にしており σ_XX = 2σ_DD の因子 4 異常で発覚 — 「規約系数は
        // 認証済みバイナリの q = 0 極限と突き合わせよ」)
        4 => vertex_unmod(&[Term { eps: 0, d: [1, 0, 0], w: 0.5 }], km),
        5 => vertex_unmod(&[Term { eps: 1, d: [0, 1, 0], w: 0.5 }], km),
        6 => vertex_unmod(&[Term { eps: 3, d: [0, 0, 1], w: 0.5 }], km),
        7 => {
            let v = vertex_unmod(&t_split_terms(0, 1), km);
            v.iter().map(|x| x / if zmut { 1.0 } else { 2.0 }).collect()
        }
        8 => {
            let v = vertex_unmod(&t_split_terms(0, 2), km);
            v.iter().map(|x| x / if zmut { 1.0 } else { 2.0 }).collect()
        }
        _ => {
            let v = vertex_unmod(&t_split_terms(1, 2), km);
            v.iter().map(|x| x / 2.0).collect()
        }
    };
    c8_from_real(&re)
}

/// 保存フラックス (v26.9-B と同じ一様構成)
fn flux(nu: usize, k: [f64; 3], q: f64, m: f64) -> C8 {
    let a = c8_from_real(&h8([k[0], k[1] + q, k[2]], m));
    let b = c8_from_real(&h8(k, m));
    let d = source(nu, k, q, m, false);
    let am = c8_mul(&a, &d);
    let mb = c8_mul(&d, &b);
    let dd = 2.0 * (0.5 * q).sin();
    (0..64)
        .map(|i| ((am[i].0 - mb[i].0) / dd, (am[i].1 - mb[i].1) / dd))
        .collect()
}

// ---------------- 殻積分 (v268t の器械 + 10 成分行列化) ----------------

fn e8m(k: [f64; 3]) -> f64 {
    (k[0].cos().powi(2) + k[1].cos().powi(2) + k[2].cos().powi(2)).sqrt()
}

fn f_pair(p: [f64; 3], q_lat: f64) -> f64 {
    let c = PI / 2.0;
    let k = [c + p[0], c + p[1], c + p[2]];
    e8m([k[0], k[1] + q_lat, k[2]]) + e8m(k)
}

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

/// 殻上 1 点: 対ごとの (M_D, M_X, M₀ν[4], M_Fν[4]) — massless
/// (D = (M₄ − M₆)/√2, X = M₈ [Z 補正済み; zmut で落とす])
#[allow(clippy::type_complexity)]
fn mels_dx(
    k: [f64; 3],
    q_lat: f64,
    zmut: bool,
) -> Vec<((f64, f64), (f64, f64), [(f64, f64); 4], [(f64, f64); 4])> {
    let hk = h8(k, 0.0);
    let (_, vk) = jacobi_eigh(&hk, 8);
    let kq = [k[0], k[1] + q_lat, k[2]];
    let hq = h8(kq, 0.0);
    let (_, vq) = jacobi_eigh(&hq, 8);
    let v4 = source(4, k, q_lat, 0.0, false);
    let v6 = source(6, k, q_lat, 0.0, false);
    let v8 = source(8, k, q_lat, 0.0, zmut);
    let dens: Vec<C8> = (0..4).map(|nu| source(nu, k, q_lat, 0.0, false)).collect();
    let flx: Vec<C8> = (0..4).map(|nu| flux(nu, k, q_lat, 0.0)).collect();
    let melc = |v: &C8, mu: usize, nu: usize| -> (f64, f64) {
        let (mut re, mut im) = (0.0f64, 0.0f64);
        for r in 0..8 {
            let (mut a, mut b) = (0.0f64, 0.0f64);
            for cc in 0..8 {
                a += v[cc + r * 8].0 * vk[cc + nu * 8];
                b += v[cc + r * 8].1 * vk[cc + nu * 8];
            }
            re += vq[r + mu * 8] * a;
            im += vq[r + mu * 8] * b;
        }
        (re, im)
    };
    let r2i = 1.0 / (2.0f64).sqrt();
    let mut out = Vec::with_capacity(16);
    for mu in 4..8 {
        for nu in 0..4 {
            let m4 = melc(&v4, mu, nu);
            let m6 = melc(&v6, mu, nu);
            let md = ((m4.0 - m6.0) * r2i, (m4.1 - m6.1) * r2i);
            let mx = melc(&v8, mu, nu);
            let mut m0 = [(0.0f64, 0.0f64); 4];
            let mut mf = [(0.0f64, 0.0f64); 4];
            for nn in 0..4 {
                m0[nn] = melc(&dens[nn], mu, nu);
                mf[nn] = melc(&flx[nn], mu, nu);
            }
            out.push((md, mx, m0, mf));
        }
    }
    out
}

/// 殻積分の D/X スペクトルと正準横断性: (σ_DD, σ_XX, σ_DX, 横断性 max 相対残差)
fn sigma_dx(a: f64, e_phys: f64, q_phys: f64, nth: usize, nph: usize, zmut: bool) -> (f64, f64, f64, f64) {
    let e_lat = a * e_phys;
    let q_lat = a * q_phys;
    let gl = {
        let n = nth;
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
    let c = PI / 2.0;
    let (mut sdd, mut sxx, mut sdx) = (0.0f64, 0.0f64, 0.0f64);
    let mut tviol = 0.0f64;
    let qhat_lat = 2.0 * (0.5 * q_lat).sin();
    for (ct, wt) in gl.0.iter().zip(&gl.1) {
        let st = (1.0 - ct * ct).sqrt();
        for j in 0..nph {
            let ph = (j as f64 + 0.5) * 2.0 * PI / nph as f64;
            let n = [st * ph.cos(), st * ph.sin(), *ct];
            let mut r_hi = e_lat;
            let mut guard = 0;
            while f_pair([r_hi * n[0], r_hi * n[1], r_hi * n[2]], q_lat) <= e_lat && guard < 40 {
                r_hi *= 1.5;
                guard += 1;
            }
            let (mut lo, mut hi) = (0.0f64, r_hi);
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
            let mels = mels_dx(k, q_lat, zmut);
            let dfr = df_dr(p, n, q_lat).abs().max(1e-12);
            let wgt = wt * (2.0 * PI / nph as f64) * r * r / dfr;
            for (md, mx, m0, mf) in &mels {
                sdd += wgt * (md.0 * md.0 + md.1 * md.1);
                sxx += wgt * (mx.0 * mx.0 + mx.1 * mx.1);
                sdx += wgt * (md.0 * mx.0 + md.1 * mx.1);
                for nn in 0..4 {
                    let vr = e_lat * m0[nn].0 - qhat_lat * mf[nn].0;
                    let vi = e_lat * m0[nn].1 - qhat_lat * mf[nn].1;
                    let den = (e_lat * e_lat * (m0[nn].0.powi(2) + m0[nn].1.powi(2))
                        + qhat_lat * qhat_lat * (mf[nn].0.powi(2) + mf[nn].1.powi(2)))
                    .sqrt()
                    .max(1e-10);
                    tviol = tviol.max((vr * vr + vi * vi).sqrt() / den);
                }
            }
        }
    }
    let norm = (2.0 * PI).powi(3);
    (sdd / norm, sxx / norm, sdx / norm, tviol)
}

fn main() {
    self_test();
    println!("=== v26.9-C v269c_spin4d — 4D spin 分離と Gate 5 総括 ===\n");
    println!("置換則で T_xy/T_yz point-split を構成 (全 10 source 完成)。殻積分の 10×10");
    println!("スペクトル行列を Lorentz BR 射影子で分解 — σ = ρ₂P₂ + ρ₀P₀s への崩壊を判定。\n");
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
    let m0 = 0.3f64;

    // ---- [C0] 置換則の独立照合 ----
    {
        // d_y = 1 (y ボンド): 明示両向き和 = (−1)^{s_y}cos(k_y + q/2) [v26.9-0 導出]
        // d_y = 2 (2 ホップ): 明示 = −(1/2)sin(2k_y + q) [v26.9-A 導出]
        let mut worst = 0.0f64;
        for &(k, q) in &[([0.7f64, -0.4, 1.9], 0.6f64), ([1.2, 2.1, 0.3], 0.15)] {
            // 置換則版
            let km = [k[0], k[1] + 0.5 * q, k[2]];
            let v1_rule = vertex_unmod(&[Term { eps: 1, d: [0, 1, 0], w: 0.25 }], km);
            let v2_rule = -0.5 * (2.0 * km[1]).sin();
            // 明示版
            for s in 0..8usize {
                let v1_exp = if sbit(s, 1) == 0 { 1.0 } else { -1.0 } * (k[1] + 0.5 * q).cos() * 0.5;
                worst = worst.max((v1_rule[s + (s ^ 1) * 8] - v1_exp).abs());
            }
            let v2_exp = -0.5 * (2.0 * k[1] + q).sin();
            worst = worst.max((v2_rule - v2_exp).abs());
        }
        check(
            "[C0] 置換則 V(k;q) = V_unmod(k+q/2ŷ): d_y = 1, 2 の明示両向き和と一致",
            worst < 1e-15,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [C1] Belinfante 構造の分解: V_Fx vs piece2 / vs Belinfante ----
    {
        let nvec = [0.55f64, 0.2, -0.81];
        let nn = (nvec[0] * nvec[0] + nvec[1] * nvec[1] + nvec[2] * nvec[2]).sqrt();
        let n = [nvec[0] / nn, nvec[1] / nn, nvec[2] / nn];
        let c = PI / 2.0;
        let mut rel_p2 = Vec::new();
        let mut rel_bel = Vec::new();
        for &eps in &[0.4f64, 0.2, 0.1, 0.05] {
            let mphys = 0.5 * eps;
            let k = [c + eps * n[0], c + eps * n[1], c + eps * n[2]];
            let qq = 0.8 * eps;
            let km = [k[0], k[1] + 0.5 * qq, k[2]];
            let vf = flux(1, k, qq, mphys);
            // piece2 のみ (eps = flip[y] = 1 の 4 項, 元の重み — Z 補正なし)
            let terms: Vec<Term> = t_split_terms(0, 1).into_iter().filter(|t| t.eps == 1).collect();
            let vp2 = c8_from_real(&vertex_unmod(&terms, km));
            let vbel = source(7, k, qq, mphys, false);
            let hk = h8(k, mphys);
            let (_, vk) = jacobi_eigh(&hk, 8);
            let kq = [k[0], k[1] + qq, k[2]];
            let hq = h8(kq, mphys);
            let (_, vq) = jacobi_eigh(&hq, 8);
            let block = |v: &C8| -> f64 {
                let mut sm = 0.0f64;
                for mu in 4..8 {
                    for nu in 0..4 {
                        let (mut re, mut im) = (0.0f64, 0.0f64);
                        for r in 0..8 {
                            let (mut a, mut b) = (0.0f64, 0.0f64);
                            for cc in 0..8 {
                                a += v[cc + r * 8].0 * vk[cc + nu * 8];
                                b += v[cc + r * 8].1 * vk[cc + nu * 8];
                            }
                            re += vq[r + mu * 8] * a;
                            im += vq[r + mu * 8] * b;
                        }
                        sm += re * re + im * im;
                    }
                }
                sm.sqrt()
            };
            let dp: C8 = (0..64).map(|i| (vf[i].0 - vp2[i].0, vf[i].1 - vp2[i].1)).collect();
            let db: C8 = (0..64).map(|i| (vf[i].0 - vbel[i].0, vf[i].1 - vbel[i].1)).collect();
            rel_p2.push((eps, block(&dp) / block(&vf)));
            rel_bel.push((eps, block(&db) / block(&vf)));
        }
        let mut msg = String::new();
        for i in 0..4 {
            msg = format!("{} [{}: p2 {:.2e} / Bel {:.2e}]", msg, rel_p2[i].0, rel_p2[i].1, rel_bel[i].1);
        }
        // run2 の発見: piece2 一致は O(ε²) どころか厳密恒等 (機械精度) —
        // [h, V₀x]/q̂ = piece2 が格子恒等式
        let ok = rel_p2.iter().all(|r| r.1 < 1e-12) && rel_bel[3].1 > 0.5;
        check(
            "[C1] 正準フラックス V_Fx = piece2 単独 (厳密恒等 < 1e-12) / Belinfante とは回転流の分 O(1)",
            ok,
            format!("{}", msg),
        );
    }

    // ---- [C2] Ward 4 行 × 10 列 (サンプル k) ----
    {
        // v26.9-B の ward_at_k 相当を 10 列で: 密度行 nu ∈ 0..4, 列 b ∈ 0..10
        let samples = [[0.5f64, 1.2, -0.8], [1.7, -0.4, 0.9]];
        let q = 0.4f64;
        let qhat = 2.0 * (0.5 * q).sin();
        let mut worst = 0.0f64;
        for &k in &samples {
            let mut recs = Vec::new();
            let mut cmax = 0.0f64;
            let (wk, vk) = jacobi_eigh(&h8(k, m0), 8);
            let kq = [k[0], k[1] + q, k[2]];
            let (wq, vq) = jacobi_eigh(&h8(kq, m0), 8);
            let kmq = [k[0], k[1] - q, k[2]];
            let (wm, vm) = jacobi_eigh(&h8(kmq, m0), 8);
            let melc = |v: &C8, va: &[f64], mu: usize, vb: &[f64], nu: usize| -> (f64, f64) {
                let (mut re, mut im) = (0.0f64, 0.0f64);
                for r in 0..8 {
                    let (mut a, mut b) = (0.0f64, 0.0f64);
                    for cc in 0..8 {
                        a += v[cc + r * 8].0 * vb[cc + nu * 8];
                        b += v[cc + r * 8].1 * vb[cc + nu * 8];
                    }
                    re += va[r + mu * 8] * a;
                    im += va[r + mu * 8] * b;
                }
                (re, im)
            };
            for nu_row in 0..4 {
                let a_v = source(nu_row, k, q, m0, false);
                let j_v = flux(nu_row, k, q, m0);
                let a_cr = source(nu_row, kmq, q, m0, false);
                let j_cr = flux(nu_row, kmq, q, m0);
                for b in 0..10 {
                    let b_rev = source(b, kq, -q, m0, false);
                    let b_cr = source(b, k, -q, m0, false);
                    for &q0 in &[0.3f64, 0.9] {
                        let (mut lhs, mut cp) = ((0.0f64, 0.0f64), (0.0f64, 0.0f64));
                        for mu in 4..8 {
                            for nuo in 0..4 {
                                let ma = melc(&a_v, &vq, mu, &vk, nuo);
                                let mj = melc(&j_v, &vq, mu, &vk, nuo);
                                let nb = melc(&b_rev, &vk, nuo, &vq, mu);
                                let de = wq[mu] - wk[nuo];
                                let den = de * de + q0 * q0;
                                let (gre, gim) = (de / den, q0 / den);
                                let pr = |x: (f64, f64), y: (f64, f64)| {
                                    (x.0 * y.0 - x.1 * y.1, x.0 * y.1 + x.1 * y.0)
                                };
                                let am = pr(nb, ma);
                                let jm = pr(nb, mj);
                                let cab = (am.0 * gre - am.1 * gim, am.0 * gim + am.1 * gre);
                                let cjb = (jm.0 * gre - jm.1 * gim, jm.0 * gim + jm.1 * gre);
                                lhs.0 += -q0 * cab.1 - qhat * cjb.0;
                                lhs.1 += q0 * cab.0 - qhat * cjb.1;
                                cp.0 += am.0;
                                cp.1 += am.1;
                            }
                        }
                        for mu in 4..8 {
                            for nuo in 0..4 {
                                let mb = melc(&b_cr, &vm, mu, &vk, nuo);
                                let na = melc(&a_cr, &vk, nuo, &vm, mu);
                                let nj = melc(&j_cr, &vk, nuo, &vm, mu);
                                let de = wm[mu] - wk[nuo];
                                let den = de * de + q0 * q0;
                                let (gre, gim) = (de / den, -q0 / den);
                                let pr = |x: (f64, f64), y: (f64, f64)| {
                                    (x.0 * y.0 - x.1 * y.1, x.0 * y.1 + x.1 * y.0)
                                };
                                let am = pr(na, mb);
                                let jm = pr(nj, mb);
                                let cab = (am.0 * gre - am.1 * gim, am.0 * gim + am.1 * gre);
                                let cjb = (jm.0 * gre - jm.1 * gim, jm.0 * gim + jm.1 * gre);
                                lhs.0 += -q0 * cab.1 - qhat * cjb.0;
                                lhs.1 += q0 * cab.0 - qhat * cjb.1;
                                cp.0 -= am.0;
                                cp.1 -= am.1;
                            }
                        }
                        cmax = cmax.max((cp.0 * cp.0 + cp.1 * cp.1).sqrt());
                        recs.push((lhs, cp));
                    }
                }
            }
            for (lhs, cp) in recs {
                let scale = (cp.0 * cp.0 + cp.1 * cp.1).sqrt() + 0.01 * cmax;
                worst = worst
                    .max(((lhs.0 + cp.0).powi(2) + (lhs.1 + cp.1).powi(2)).sqrt() / scale);
            }
        }
        check(
            "[C2] full 4D Ward の 10 列拡張: 4 行 × 10 列 × 2 周波数 = 80 恒等式 (2 サンプル k)",
            worst < 1e-10,
            format!("max 相対残差 = {:.1e}", worst),
        );
    }

    // ---- [C3] 4D spin-2 分離 (D/X は timelike q_L でも厳密 P₂) ----
    {
        let (e_phys, q_phys) = (1.5f64, 0.6);
        // P₂ 幾何: ⟨D|P₂|D⟩ = 1, ⟨X|P₂|X⟩ = 1/2 (θ_xx = θ_zz = −1, θ_xz = 0)
        // ⇒ 縮退関係は **σ_DD = 2σ_XX**。絶対アンカー: σ_DD^phys = 2ρ_D(s = E²−Q²)
        let s_inv = e_phys * e_phys - q_phys * q_phys;
        let rho_d = s_inv * s_inv / (160.0 * PI * PI);
        let mut rows = Vec::new();
        for &a in &[0.18f64, 0.09, 0.045] {
            let (sdd, sxx, sdx, tviol) = sigma_dx(a, e_phys, q_phys, 32, 64, false);
            let deg = (sdd / (2.0 * sxx) - 1.0).abs();
            let orth = sdx.abs() / (sdd * sxx).sqrt();
            let anchor = sdd / a.powi(4) / (2.0 * rho_d);
            println!(
                "    [C3 表] a = {:.3}: |σ_DD/2σ_XX − 1| = {:.5}, σ_DD/(2ρ_D) = {:.4}, |σ_DX|/√(σσ) = {:.2e}, 横断性 max = {:.1e} ({} s)",
                a, deg, anchor, orth, tviol, t0.elapsed().as_secs()
            );
            rows.push((a, deg, orth, tviol, anchor));
        }
        let r1 = rows[0].1 / rows[1].1;
        let r2 = rows[1].1 / rows[2].1;
        let ok_deg = rows[2].1 < 0.02 && (2.5..6.0).contains(&r1) && (2.5..6.0).contains(&r2);
        let ok_orth = rows[2].2 < 1e-6;
        let ok_tv = rows.iter().all(|r| r.3 < 1e-9);
        let ok_anchor = (rows[2].4 - 1.0).abs() < 0.02;
        check(
            "[C3] 4D spin-2 分離: σ_DD/2σ_XX → 1 (O(a²)) / σ_DX = 0 / 正準横断性 / oracle アンカー 2%",
            ok_deg && ok_orth && ok_tv && ok_anchor,
            format!(
                "縮退 {:.4} → {:.4} → {:.4} (比 {:.1}, {:.1}), 直交 {:.1e}, 横断 {:.1e}, アンカー {:.4}",
                rows[0].1, rows[1].1, rows[2].1, r1, r2, rows[2].2, rows[2].3, rows[2].4
            ),
        );
    }

    // ---- [C4] 変異 ----
    {
        let (sdd, sxx, _, _) = sigma_dx(0.045, 1.5, 0.6, 32, 64, true);
        let ratio = 2.0 * sxx / sdd;
        check(
            "[C4] 変異: T_xz の Z 補正 (÷2) を落とす → 2σ_XX/σ_DD = 4 (縮退の 4 倍逸脱)",
            (ratio - 4.0).abs() < 0.2,
            format!("2σ_XX/σ_DD = {:.3}", ratio),
        );
    }

    // ---- Gate 5 総括 ----
    println!("\n[Gate 5 総括 (spec §12.6 の項目別)]");
    println!("  確立: h₀₀/h₀ᵢ source / q₀ ≠ 0 (Matsubara Ward 64+80 恒等式) / 全 10 source");
    println!("        構成 (置換則) / 4D Ward = 局所カレント + 計算可能な接触項 / 正準テンソル");
    println!("        の厳密横断性 / spin-2 の 4D 分離 (P₂ 縮退 + 直交) / BOND-A = 保存 stress");
    println!("        (tree O(ε²)) / Belinfante 構造の格子分解 (正準 = piece 単独・対称化 = 平均)");
    println!("  残り: Belinfante 対称 10×10 の完全崩壊 (T⁰ⁱ の対称化 — 二重変調則で構成可,");
    println!("        v26.9-D) / temporal h の二次変分 scheme (v27.0 設計項目)。");
    println!("        **FullGravitationalVacuumPolarization 型は保留を維持・1/Π 禁止も維持**");
    println!("        (凍結解釈: 測定器の証明であり QRN・創発重力の証拠ではない)");

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.9-C".into())),
        ("kind".into(), Json::Str("spin4d_separation_gate5".into())),
        (
            "substitution_rule".into(),
            Json::Str("V(k;q) = V_unmod(k + q/2 ŷ) — d_y = 1,2 で独立照合".into()),
        ),
    ]);
    let p = write_artifact("results/v269c_spin4d.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **10 source 完成・full 4D Ward 80 恒等式・正準横断性厳密・spin-2 の 4D 分離成立 + Belinfante 構造の格子分解 — Gate 5 の Ward/分離部門確立** (対称 10×10 完全崩壊 = v26.9-D, temporal 二次変分 = v27.0 設計)"
        } else {
            "FAIL あり — 分岐 (b) 分離の破れ (公表) / (c) 器械。欄が一次ソース"
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
