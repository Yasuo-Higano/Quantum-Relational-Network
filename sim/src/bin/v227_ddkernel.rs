//! v22.7 DD 全経路核 — f64 κ 床の突破試験 (第二十三期, G 窓 named limit (i) への攻撃)
//!
//! v22.2/22.4 で確定した f64 κ 床: カナリア R₁ は d=2.0 で 2%・d=2.5 で 5.5% 狂い
//! (N 非依存 = 波束の絶対深さで噛む)。原因は κ = ln((1−c)/c) の f64 分解能
//! (c ~ e^{−κ} < 1e-16 が表現不能) — 従来の DD-Rayleigh (v22.2 v1) は f64 固有
//! ベクトルの縮退クラスタ混合で失敗した。本版は**経路全体を DD 化**する:
//!   [1] f64 jacobi で H の占有部分空間 → **DD 部分空間 Newton** 1 回で 1e-30 へ
//!       (X[e,o] = R̃[e,o]/(λ_o − λ_e); 閉殻ギャップ ≥ 0.18 が全分母を保証・二次収束)
//!   [2] C_A を DD で構成 (864²) → DD jacobi → 真の κ (床 ~70 まで)
//!   [3] Gauss 対 rank-2 回転・エントロピー・⟨K_A⟩ を全て DD で評価
//! 検定: R₁ = δS/δ⟨K_A⟩ (中心差分 ±α, ±2α + Richardson) を d ∈ {1.5, 2.5, 3.5} で測る。
//!
//! 開発記録 (run1, 保存): (i) 初版の [G1] はベクトル毎の残差を測っており、占有多重項
//!   内の f64 混合 (射影子には無害) まで拾って FAIL — 正しい対象は**部分空間残差**
//!   (I − VVᵀ)HV。(ii) run1 の実質は突破だった: DD C_A が c_min = 5.2e-18
//!   (κ_max = 39.8 > f64 飽和 32.2) を解像。(iii) d=2.5 の残差 1.7% は δ⟨K⟩ の
//!   α 非線形 (α=0.02→0.04 で ×2.43) が示す α² 汚染 — DD の分解能なら α を
//!   1 桁下げて漸近域に入れる (基準は 1±1% のまま = 精細化であって緩和ではない)。
//!
//! 装置ゲート v2: [G1] 部分空間残差 max|(I−VVᵀ)HV| < 1e-24,
//!   [G2] DD と f64 の C_A 固有値が非深部 (1e-8 < c < 1−1e-8) で一致 ±1e-10,
//!   [G3] d=1.5 (f64 生存域): R₁(DD) = 1 ± 1% (器械連続性)。
//! 事前登録 v2 (α = 0.005, 0.01 の漸近対 + d = 2.0 追加):
//!   (a) R₁(DD, d=2.5) = 1 ± 1% = **床の突破** — カナリア復活、深部が測定域になる
//!       (d=2.0 も 1±1% なら v22.4 の境界は「f64 固有」と書き換え) /
//!   (a′) ゲート PASS だが d=2.5 で外れ = 信号消失が真因 — named limit の書き換え /
//!   (b) ゲート故障 = 記録。d=3.5 と深部 R₂ (K_boost 分母, 変種 A) は記録。

use uft_sim::*;

// ---------------- double-double 演算 (Dekker/Knuth) ----------------
#[derive(Clone, Copy, Debug)]
struct Dd {
    hi: f64,
    lo: f64,
}

const DD0: Dd = Dd { hi: 0.0, lo: 0.0 };

fn qts(a: f64, b: f64) -> Dd {
    // quick_two_sum: |a| ≥ |b| 前提
    let s = a + b;
    Dd {
        hi: s,
        lo: b - (s - a),
    }
}
fn two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    let bb = s - a;
    Dd {
        hi: s,
        lo: (a - (s - bb)) + (b - bb),
    }
}
fn two_prod(a: f64, b: f64) -> Dd {
    let p = a * b;
    Dd {
        hi: p,
        lo: a.mul_add(b, -p),
    }
}
fn dd(a: f64) -> Dd {
    Dd { hi: a, lo: 0.0 }
}
fn dd_add(x: Dd, y: Dd) -> Dd {
    let s = two_sum(x.hi, y.hi);
    qts(s.hi, s.lo + x.lo + y.lo)
}
fn dd_neg(x: Dd) -> Dd {
    Dd {
        hi: -x.hi,
        lo: -x.lo,
    }
}
fn dd_sub(x: Dd, y: Dd) -> Dd {
    dd_add(x, dd_neg(y))
}
fn dd_mul(x: Dd, y: Dd) -> Dd {
    let p = two_prod(x.hi, y.hi);
    qts(p.hi, p.lo + x.hi * y.lo + x.lo * y.hi)
}
fn dd_mul_f(x: Dd, y: f64) -> Dd {
    let p = two_prod(x.hi, y);
    qts(p.hi, p.lo + x.lo * y)
}
fn dd_div(x: Dd, y: Dd) -> Dd {
    let q1 = x.hi / y.hi;
    let r1 = dd_sub(x, dd_mul_f(y, q1));
    let q2 = r1.hi / y.hi;
    let r2 = dd_sub(r1, dd_mul_f(y, q2));
    let q3 = r2.hi / y.hi;
    dd_add(qts(q1, q2), dd(q3))
}
fn dd_sqrt(a: Dd) -> Dd {
    if a.hi <= 0.0 {
        return DD0;
    }
    let x0 = a.hi.sqrt();
    // 1 回の DD Newton: x = (x + a/x)/2
    let x = dd(x0);
    dd_mul_f(dd_add(x, dd_div(a, x)), 0.5)
}
// DD exp: 範囲縮約 z = k·ln2 + f, |f| ≤ ln2/2, Taylor 26 項
fn dd_exp(z: Dd) -> Dd {
    let ln2 = Dd {
        hi: 0.6931471805599453,
        lo: 2.3190468138462996e-17,
    };
    if z.hi < -745.0 {
        return DD0;
    }
    if z.hi > 709.0 {
        return dd(f64::INFINITY);
    }
    let k = (z.hi / ln2.hi).round();
    let f = dd_sub(z, dd_mul_f(ln2, k));
    let mut term = dd(1.0);
    let mut sum = dd(1.0);
    for i in 1..27 {
        term = dd_mul_f(dd_mul(term, f), 1.0 / i as f64);
        sum = dd_add(sum, term);
    }
    // × 2^k (指数スケーリングは正確)
    let scale = (2.0f64).powi(k as i32);
    Dd {
        hi: sum.hi * scale,
        lo: sum.lo * scale,
    }
}
// DD ln: y₀ = f64 ln + DD Newton 1 回 (y ← y + a·e^{−y} − 1)
fn dd_ln(a: Dd) -> Dd {
    if a.hi <= 0.0 {
        return dd(-f64::INFINITY);
    }
    let y0 = a.hi.ln();
    let corr = dd_sub(dd_mul(a, dd_exp(dd(-y0))), dd(1.0));
    dd_add(dd(y0), corr)
}

// ---------------- DD cyclic Jacobi (対称, 固有値のみ + 直交系) ----------------
// 返り値: (固有値 昇順, 固有ベクトル列優先 v[i + k*n])
fn jacobi_dd(a_in: &[Dd], n: usize, sweeps_max: usize) -> (Vec<Dd>, Vec<Dd>) {
    let mut a = a_in.to_vec();
    let mut v = vec![DD0; n * n];
    for i in 0..n {
        v[i + i * n] = dd(1.0);
    }
    for _sw in 0..sweeps_max {
        let mut off = 0.0f64;
        for p in 0..n {
            for q in p + 1..n {
                off = off.max(a[p * n + q].hi.abs());
            }
        }
        if off < 1e-28 {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p * n + q];
                if apq.hi.abs() < 1e-30 {
                    continue;
                }
                let app = a[p * n + p];
                let aqq = a[q * n + q];
                // theta = (aqq − app)/(2 apq)
                let theta = dd_div(dd_sub(aqq, app), dd_mul_f(apq, 2.0));
                let t = {
                    let at = theta.hi.abs();
                    let den = dd_add(dd(at), dd_sqrt(dd_add(dd_mul(theta, theta), dd(1.0))));
                    let tt = dd_div(dd(1.0), den);
                    if theta.hi < 0.0 {
                        dd_neg(tt)
                    } else {
                        tt
                    }
                };
                let c = dd_div(dd(1.0), dd_sqrt(dd_add(dd_mul(t, t), dd(1.0))));
                let s = dd_mul(t, c);
                // 行列更新 (対称 Jacobi 回転)
                for i in 0..n {
                    let aip = a[i * n + p];
                    let aiq = a[i * n + q];
                    a[i * n + p] = dd_sub(dd_mul(c, aip), dd_mul(s, aiq));
                    a[i * n + q] = dd_add(dd_mul(s, aip), dd_mul(c, aiq));
                }
                for i in 0..n {
                    let api = a[p * n + i];
                    let aqi = a[q * n + i];
                    a[p * n + i] = dd_sub(dd_mul(c, api), dd_mul(s, aqi));
                    a[q * n + i] = dd_add(dd_mul(s, api), dd_mul(c, aqi));
                }
                for i in 0..n {
                    let vip = v[i + p * n];
                    let viq = v[i + q * n];
                    v[i + p * n] = dd_sub(dd_mul(c, vip), dd_mul(s, viq));
                    v[i + q * n] = dd_add(dd_mul(s, vip), dd_mul(c, viq));
                }
            }
        }
    }
    // 固有値と昇順ソート
    let mut idx: Vec<usize> = (0..n).collect();
    let evs: Vec<Dd> = (0..n).map(|i| a[i * n + i]).collect();
    idx.sort_by(|&i, &j| evs[i].hi.partial_cmp(&evs[j].hi).unwrap());
    let evs_s: Vec<Dd> = idx.iter().map(|&i| evs[i]).collect();
    let mut v_s = vec![DD0; n * n];
    for (k, &i) in idx.iter().enumerate() {
        for r in 0..n {
            v_s[r + k * n] = v[r + i * n];
        }
    }
    (evs_s, v_s)
}

// 3D staggered H (x 開放, y/z 周期) — v22.2 と同一 (成分は正確な半整数 → DD 誤差ゼロ)
fn build_h3d_periodic(n: usize) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let mut add = |i: usize, j: usize, t: f64| {
        h[j + i * ns] += t;
        h[i + j * ns] += t;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                if x + 1 < n {
                    add(i, idx(x + 1, y, z), 0.5);
                }
                let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                add(i, idx(x, (y + 1) % n, z), ey);
                let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                add(i, idx(x, y, (z + 1) % n), ez);
            }
        }
    }
    h
}

// S(C) を DD 固有値から (クランプ 1e-200)
fn entropy_dd(evs: &[Dd]) -> Dd {
    let mut s = DD0;
    for &c in evs {
        for &p in &[c, dd_sub(dd(1.0), c)] {
            if p.hi > 1e-200 && p.hi < 1.0 {
                s = dd_sub(s, dd_mul(p, dd_ln(p)));
            }
        }
    }
    s
}

fn main() {
    self_test();
    println!("=== v22.7 DD 全経路核 — f64 κ 床の突破試験 (第二十三期) ===\n");
    println!("事前登録: (a) R₁(DD, d=2.5) = 1 ± 1% = 床の突破 / (a′) ゲート PASS だが");
    println!("          外れ = 信号消失が真因 (named limit 書き換え) / (b) ゲート故障\n");
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
    let t0 = std::time::Instant::now();
    let n = 12usize;
    let ns = n * n * n;
    let nocc = ns / 2;

    // ---- [1] f64 対角化 + DD 部分空間 Newton ----
    let h = build_h3d_periodic(n);
    let (ev, vv) = jacobi_eigh(&h, ns);
    let gap = ev[nocc] - ev[nocc - 1];
    println!(
        "    f64 jacobi 完了: gap = {:.4} ({} s)",
        gap,
        t0.elapsed().as_secs()
    );
    // R̃[m, k] = Σ_ij V[i,m] H[i,j] V[j,k] (全モード × 占有) を DD で
    // (f64 固有基底では対角 ± 1e-15 — 非対角が部分空間誤差)
    // HV[i,k] (DD): k = 占有
    let mut hv = vec![DD0; ns * nocc];
    for i in 0..ns {
        for j in 0..ns {
            let hij = h[i * ns + j];
            if hij == 0.0 {
                continue;
            }
            for k in 0..nocc {
                hv[i * nocc + k] = dd_add(hv[i * nocc + k], dd_mul_f(dd(vv[j + k * ns]), hij));
            }
        }
    }
    // [診断] 補正前の部分空間残差 (f64 V): サンプル k で (I − VVᵀ)Hv_k
    {
        let mut res0 = 0.0f64;
        for &k in &[0usize, nocc / 2, nocc - 1] {
            let mut hvk: Vec<Dd> = (0..ns).map(|i| hv[i * nocc + k]).collect();
            for j in 0..nocc {
                let mut ip = DD0;
                for i in 0..ns {
                    ip = dd_add(ip, dd_mul_f(hvk[i], vv[i + j * ns]));
                }
                for i in 0..ns {
                    hvk[i] = dd_sub(hvk[i], dd_mul_f(ip, vv[i + j * ns]));
                }
            }
            for i in 0..ns {
                res0 = res0.max(hvk[i].hi.abs());
            }
        }
        println!("    [診断] 補正前残差 (f64 V) = {:.2e}", res0);
    }
    // X[m,k] = (⟨m|Hv_k⟩ − θ_k⟨m|v_k⟩)/(θ_k − λ_m) — 第 2 項が f64 基底の交差汚染を
    // 相殺する (run2/3 の教訓: これを欠くと補正が無効化し残差 1.8e-15 のまま)
    let nemp = ns - nocc;
    let mut xmat = vec![DD0; nemp * nocc];
    for m in 0..nemp {
        let vm = &vv[(nocc + m) * ns..(nocc + m + 1) * ns];
        for k in 0..nocc {
            let mut acc = DD0;
            let mut ovl = DD0;
            for i in 0..ns {
                acc = dd_add(acc, dd_mul_f(hv[i * nocc + k], vm[i]));
                ovl = dd_add(ovl, dd_mul_f(dd(vv[i + k * ns]), vm[i]));
            }
            acc = dd_sub(acc, dd_mul_f(ovl, ev[k]));
            let den = ev[k] - ev[nocc + m];
            xmat[m * nocc + k] = dd_mul_f(acc, 1.0 / den);
        }
    }
    {
        let xmax = xmat.iter().map(|x| x.hi.abs()).fold(0.0f64, f64::max);
        println!("    [診断] ‖X‖_max = {:.2e}", xmax);
    }
    // V' = V_occ + V_emp·X (DD, 列優先 ns × nocc)
    let mut vocc = vec![DD0; ns * nocc];
    for k in 0..nocc {
        for i in 0..ns {
            vocc[i + k * ns] = dd(vv[i + k * ns]);
        }
    }
    for m in 0..nemp {
        let vm = &vv[(nocc + m) * ns..(nocc + m + 1) * ns];
        for k in 0..nocc {
            let x = xmat[m * nocc + k];
            if x.hi == 0.0 && x.lo == 0.0 {
                continue;
            }
            for i in 0..ns {
                vocc[i + k * ns] = dd_add(vocc[i + k * ns], dd_mul_f(x, vm[i]));
            }
        }
    }
    // DD Gram-Schmidt (列直交化)
    for k in 0..nocc {
        for j in 0..k {
            let mut ip = DD0;
            for i in 0..ns {
                ip = dd_add(ip, dd_mul(vocc[i + j * ns], vocc[i + k * ns]));
            }
            for i in 0..ns {
                vocc[i + k * ns] = dd_sub(vocc[i + k * ns], dd_mul(ip, vocc[i + j * ns]));
            }
        }
        let mut nr = DD0;
        for i in 0..ns {
            nr = dd_add(nr, dd_mul(vocc[i + k * ns], vocc[i + k * ns]));
        }
        let inv = dd_div(dd(1.0), dd_sqrt(nr));
        for i in 0..ns {
            vocc[i + k * ns] = dd_mul(vocc[i + k * ns], inv);
        }
    }
    println!(
        "    DD 部分空間 Newton + 直交化 完了 ({} s)",
        t0.elapsed().as_secs()
    );
    // [G1] 部分空間残差: 占有多重項内の混合は射影子に無害なので (I − VVᵀ)HV のみ測る
    let mut res_max = 0.0f64;
    for &k in &[0usize, nocc / 4, nocc / 2, nocc - 1] {
        let mut hvk = vec![DD0; ns];
        for i in 0..ns {
            for j in 0..ns {
                let hij = h[i * ns + j];
                if hij == 0.0 {
                    continue;
                }
                hvk[i] = dd_add(hvk[i], dd_mul_f(vocc[j + k * ns], hij));
            }
        }
        for j in 0..nocc {
            let mut ip = DD0;
            for i in 0..ns {
                ip = dd_add(ip, dd_mul(vocc[i + j * ns], hvk[i]));
            }
            if ip.hi.abs() < 1e-40 {
                continue;
            }
            for i in 0..ns {
                hvk[i] = dd_sub(hvk[i], dd_mul(ip, vocc[i + j * ns]));
            }
        }
        for i in 0..ns {
            res_max = res_max.max(hvk[i].hi.abs());
        }
    }
    check(
        "[G1] DD 部分空間残差 (I−VVᵀ)HV < 1e-24",
        res_max < 1e-24,
        format!("max res = {:.1e} ({} s)", res_max, t0.elapsed().as_secs()),
    );
    if std::env::var("V227_G1ONLY").is_ok() {
        println!("(G1 デバッグ早期終了)");
        return;
    }

    // ---- [2] C_A (半空間) を DD で構成 → DD jacobi ----
    let idx3 = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut sel = Vec::new();
    for z in 0..n {
        for y in 0..n {
            for x in 0..n / 2 {
                sel.push(idx3(x, y, z));
            }
        }
    }
    let m = sel.len();
    let mut ca0 = vec![DD0; m * m];
    for k in 0..nocc {
        for (a, &ia) in sel.iter().enumerate() {
            let va = vocc[ia + k * ns];
            if va.hi.abs() < 1e-40 {
                continue;
            }
            for (b, &ib) in sel.iter().enumerate() {
                if b < a {
                    continue;
                }
                ca0[a * m + b] = dd_add(ca0[a * m + b], dd_mul(va, vocc[ib + k * ns]));
            }
        }
    }
    for a in 0..m {
        for b in 0..a {
            ca0[a * m + b] = ca0[b * m + a];
        }
    }
    println!("    DD C_A 構成完了 ({} s)", t0.elapsed().as_secs());
    let (cw_dd, cv_dd) = jacobi_dd(&ca0, m, 30);
    println!(
        "    DD jacobi C_A 完了: c_min = {:.3e}, κ_max = {:.1} ({} s)",
        cw_dd[0].hi,
        {
            let c = cw_dd[0];
            if c.hi > 0.0 {
                dd_ln(dd_div(dd_sub(dd(1.0), c), c)).hi
            } else {
                f64::INFINITY
            }
        },
        t0.elapsed().as_secs()
    );
    // [G2] 非深部固有値の f64 一致
    {
        let ca_f64: Vec<f64> = ca0.iter().map(|z| z.hi).collect();
        let (cw_f, _) = jacobi_eigh(&ca_f64, m);
        let mut dev = 0.0f64;
        for i in 0..m {
            let c = cw_dd[i].hi;
            if c > 1e-8 && c < 1.0 - 1e-8 {
                dev = dev.max((c - cw_f[i]).abs());
            }
        }
        check(
            "[G2] 非深部固有値の DD/f64 一致 ± 1e-10",
            dev < 1e-10,
            format!("max dev = {:.1e}", dev),
        );
    }
    // κ (DD)
    let kappa: Vec<Dd> = cw_dd
        .iter()
        .map(|&c| {
            let cc = if c.hi < 1e-200 {
                dd(1e-200)
            } else if c.hi > 1.0 - 1e-16 && dd_sub(dd(1.0), c).hi < 1e-200 {
                dd_sub(dd(1.0), dd(1e-200))
            } else {
                c
            };
            dd_ln(dd_div(dd_sub(dd(1.0), cc), cc))
        })
        .collect();

    // ---- [3] Gauss 対 rank-2 回転を DD で: d ∈ {1.5, 2.5, 3.5} ----
    // Cu, Cw (DD 全空間ベクトル) は C·u = Σ_occ v (v·u)
    let xplane = n as f64 / 2.0 - 0.5;
    let mid = n as f64 / 2.0;
    let wid = 1.0f64;
    let mk = |x0: f64| -> Vec<f64> {
        let mut p = vec![0.0f64; ns];
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let d2 = (x as f64 - x0).powi(2)
                        + (y as f64 - mid).powi(2)
                        + (z as f64 - mid).powi(2);
                    p[idx3(x, y, z)] = (-d2 / (2.0 * wid * wid)).exp();
                }
            }
        }
        let nr: f64 = p.iter().map(|a| a * a).sum::<f64>().sqrt();
        for a in p.iter_mut() {
            *a /= nr;
        }
        p
    };
    let mut r1_of_d: Vec<(f64, f64)> = Vec::new();
    let mut r2_of_d: Vec<(f64, f64)> = Vec::new();
    // K_boost ボンド (変種 A — 深部 R₂ の記録用)
    let mut bonds: Vec<(usize, usize, f64, f64)> = Vec::new();
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx3(x, y, z);
                if x + 1 < n {
                    let xi = xplane - (x as f64 + 0.5);
                    if xi > 1e-9 {
                        bonds.push((i, idx3(x + 1, y, z), 0.5, xi));
                    }
                }
                let xi_site = xplane - x as f64;
                if xi_site > 1e-9 {
                    let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                    bonds.push((i, idx3(x, (y + 1) % n, z), ey, xi_site));
                    let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                    bonds.push((i, idx3(x, y, (z + 1) % n), ez, xi_site));
                }
            }
        }
    }
    let two_pi = 2.0 * std::f64::consts::PI;
    for &d in &[1.5f64, 2.0, 2.5, 3.5] {
        // u, w を DD で直交化
        let uf = mk(xplane - d);
        let wf = mk(xplane - d - 1.0);
        let mut u: Vec<Dd> = uf.iter().map(|&z| dd(z)).collect();
        let mut w: Vec<Dd> = wf.iter().map(|&z| dd(z)).collect();
        // 正規化 u
        let mut nu = DD0;
        for i in 0..ns {
            nu = dd_add(nu, dd_mul(u[i], u[i]));
        }
        let inv = dd_div(dd(1.0), dd_sqrt(nu));
        for i in 0..ns {
            u[i] = dd_mul(u[i], inv);
        }
        // w ⊥ u + 正規化
        let mut ov = DD0;
        for i in 0..ns {
            ov = dd_add(ov, dd_mul(u[i], w[i]));
        }
        for i in 0..ns {
            w[i] = dd_sub(w[i], dd_mul(ov, u[i]));
        }
        let mut nw = DD0;
        for i in 0..ns {
            nw = dd_add(nw, dd_mul(w[i], w[i]));
        }
        let invw = dd_div(dd(1.0), dd_sqrt(nw));
        for i in 0..ns {
            w[i] = dd_mul(w[i], invw);
        }
        // Cu = Σ_occ v (v·u), Cw 同様 (DD)
        let mut cu = vec![DD0; ns];
        let mut cwv = vec![DD0; ns];
        for k in 0..nocc {
            let mut ipu = DD0;
            let mut ipw = DD0;
            for i in 0..ns {
                ipu = dd_add(ipu, dd_mul(vocc[i + k * ns], u[i]));
                ipw = dd_add(ipw, dd_mul(vocc[i + k * ns], w[i]));
            }
            for i in 0..ns {
                cu[i] = dd_add(cu[i], dd_mul(vocc[i + k * ns], ipu));
                cwv[i] = dd_add(cwv[i], dd_mul(vocc[i + k * ns], ipw));
            }
        }
        let mut uu = DD0;
        let mut uw = DD0;
        let mut ww = DD0;
        for i in 0..ns {
            uu = dd_add(uu, dd_mul(u[i], cu[i]));
            uw = dd_add(uw, dd_mul(u[i], cwv[i]));
            ww = dd_add(ww, dd_mul(w[i], cwv[i]));
        }
        // 中心差分 ±α, ±2α
        let a1 = 0.005f64; // DD 分解能を活かした漸近対 (run1: 0.02 は深部で α² 汚染)
        let mut rr = [0.0f64; 2];
        let mut r2r = [0.0f64; 2];
        for (ii, &al) in [a1, 2.0 * a1].iter().enumerate() {
            let mut svals = [DD0; 2];
            let mut kvals = [DD0; 2];
            for (jj, &sgn) in [1.0f64, -1.0].iter().enumerate() {
                let aa = (al * sgn).cos() - 1.0;
                let ss = (al * sgn).sin();
                // C'_A[a,b] = C_A + du_a Cu_b + dw_a Cw_b + du_b Cu_a + dw_b Cw_a
                //           + du_a (uu du_b + uw dw_b) + dw_a (uw du_b + ww dw_b)
                let mut ca = ca0.clone();
                let du: Vec<Dd> = sel
                    .iter()
                    .map(|&i| dd_add(dd_mul_f(u[i], aa), dd_mul_f(w[i], ss)))
                    .collect();
                let dw: Vec<Dd> = sel
                    .iter()
                    .map(|&i| dd_sub(dd_mul_f(w[i], aa), dd_mul_f(u[i], ss)))
                    .collect();
                let cu_a: Vec<Dd> = sel.iter().map(|&i| cu[i]).collect();
                let cw_a: Vec<Dd> = sel.iter().map(|&i| cwv[i]).collect();
                for a in 0..m {
                    for b in 0..m {
                        let mut t = ca[a * m + b];
                        t = dd_add(t, dd_mul(du[a], cu_a[b]));
                        t = dd_add(t, dd_mul(dw[a], cw_a[b]));
                        t = dd_add(t, dd_mul(du[b], cu_a[a]));
                        t = dd_add(t, dd_mul(dw[b], cw_a[a]));
                        let m1 = dd_add(dd_mul(uu, du[b]), dd_mul(uw, dw[b]));
                        let m2 = dd_add(dd_mul(uw, du[b]), dd_mul(ww, dw[b]));
                        t = dd_add(t, dd_add(dd_mul(du[a], m1), dd_mul(dw[a], m2)));
                        ca[a * m + b] = t;
                    }
                }
                // S (DD jacobi) と ⟨K⟩
                let (cw2, _) = jacobi_dd(&ca, m, 30);
                svals[jj] = entropy_dd(&cw2);
                // δ⟨K⟩ = Σ_k κ_k (v_kᵀ (C'−C) v_k)
                let mut kexp = DD0;
                for k in 0..m {
                    if kappa[k].hi.abs() < 1e-14 {
                        continue;
                    }
                    let mut acc = DD0;
                    for a in 0..m {
                        let va = cv_dd[a + k * m];
                        if va.hi.abs() < 1e-30 {
                            continue;
                        }
                        for b in 0..m {
                            let diff = dd_sub(ca[a * m + b], ca0[a * m + b]);
                            acc = dd_add(acc, dd_mul(dd_mul(va, diff), cv_dd[b + k * m]));
                        }
                    }
                    kexp = dd_add(kexp, dd_mul(kappa[k], acc));
                }
                kvals[jj] = kexp;
            }
            let ds = dd_mul_f(dd_sub(svals[0], svals[1]), 0.5);
            let dk = dd_mul_f(dd_sub(kvals[0], kvals[1]), 0.5);
            rr[ii] = dd_div(ds, dk).hi;
            // R₂ 記録: δ⟨K_boost⟩ を rank-2 更新式で全ボンド DD 評価
            let mut dkb = [DD0; 2];
            for (jj, &sgn) in [1.0f64, -1.0].iter().enumerate() {
                let aa = (al * sgn).cos() - 1.0;
                let ss = (al * sgn).sin();
                let mut acc = DD0;
                for &(i, j, t, xi) in &bonds {
                    let du_i = dd_add(dd_mul_f(u[i], aa), dd_mul_f(w[i], ss));
                    let dw_i = dd_sub(dd_mul_f(w[i], aa), dd_mul_f(u[i], ss));
                    let du_j = dd_add(dd_mul_f(u[j], aa), dd_mul_f(w[j], ss));
                    let dw_j = dd_sub(dd_mul_f(w[j], aa), dd_mul_f(u[j], ss));
                    let mut diff = dd_add(dd_mul(du_i, cu[j]), dd_mul(dw_i, cwv[j]));
                    diff = dd_add(diff, dd_add(dd_mul(du_j, cu[i]), dd_mul(dw_j, cwv[i])));
                    let m1 = dd_add(dd_mul(uu, du_j), dd_mul(uw, dw_j));
                    let m2 = dd_add(dd_mul(uw, du_j), dd_mul(ww, dw_j));
                    diff = dd_add(diff, dd_add(dd_mul(du_i, m1), dd_mul(dw_i, m2)));
                    acc = dd_add(acc, dd_mul_f(diff, two_pi * xi * t * 2.0));
                }
                dkb[jj] = acc;
            }
            let dkb_c = dd_mul_f(dd_sub(dkb[0], dkb[1]), 0.5);
            r2r[ii] = dd_div(ds, dkb_c).hi;
            println!(
                "      d={:.1} α={:.3}: δS = {:+.6e}, δ⟨K⟩ = {:+.6e}, 比 = {:.6}, R₂ 素 = {:.6} ({} s)",
                d,
                al,
                ds.hi,
                dk.hi,
                rr[ii],
                r2r[ii],
                t0.elapsed().as_secs()
            );
        }
        let r1 = (4.0 * rr[0] - rr[1]) / 3.0;
        let r2 = (4.0 * r2r[0] - r2r[1]) / 3.0;
        println!(
            "    N=12 d={:.1}: R₁(DD) = {:.5}, R₂(DD, 記録) = {:.5}",
            d, r1, r2
        );
        r1_of_d.push((d, r1));
        r2_of_d.push((d, r2));
    }

    // ---- 判定 ----
    let r15 = r1_of_d.iter().find(|r| r.0 == 1.5).unwrap().1;
    let r20 = r1_of_d.iter().find(|r| r.0 == 2.0).unwrap().1;
    let r25 = r1_of_d.iter().find(|r| r.0 == 2.5).unwrap().1;
    let r35 = r1_of_d.iter().find(|r| r.0 == 3.5).unwrap().1;
    check(
        "[G3] d=1.5 の器械連続性 R₁(DD) = 1 ± 1%",
        (r15 - 1.0).abs() < 0.01,
        format!("R₁ = {:.5}", r15),
    );
    let breakthrough = (r25 - 1.0).abs() < 0.01;
    println!(
        "\n[判定] {}",
        if nfail == 0 && breakthrough {
            "事前登録 (a): 床の突破 — DD 全経路核 + 漸近 α で d=2.5 のカナリアが復活。深部が測定域になった"
        } else if nfail == 0 {
            "事前登録 (a′): ゲートは健全だが d=2.5 で外れ — 床ではなく信号消失が真因 (named limit の書き換え)"
        } else {
            "事前登録 (b): ゲート故障 — 記録"
        }
    );
    println!(
        "    R₁(DD): d=1.5: {:.5} / d=2.0: {:.5} / d=2.5: {:.5} / d=3.5 (記録): {:.5}",
        r15, r20, r25, r35
    );
    if nfail == 0 && breakthrough && (r20 - 1.0).abs() < 0.01 {
        println!("    [境界書き換え] d=2.0 も 1±1% — v22.4 の境界は f64 固有の現象と確定");
    }
    let r2s: Vec<String> = r2_of_d
        .iter()
        .map(|(d, r)| format!("d={:.1}: {:.5}", d, r))
        .collect();
    println!("    [記録] R₂(DD 分子): {}", r2s.join(" / "));

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v22.7".into())),
        ("r1_d15".into(), Json::Num(r15)),
        ("r1_d20".into(), Json::Num(r20)),
        ("r1_d25".into(), Json::Num(r25)),
        ("r1_d35".into(), Json::Num(r35)),
        (
            "r2_records".into(),
            Json::Arr(
                r2_of_d
                    .iter()
                    .map(|&(d, r)| {
                        Json::Obj(vec![
                            ("d".into(), Json::Num(d)),
                            ("r2".into(), Json::Num(r)),
                        ])
                    })
                    .collect(),
            ),
        ),
        (
            "breakthrough".into(),
            Json::Bool(nfail == 0 && breakthrough),
        ),
    ]);
    let p = write_artifact("results/v227_ddkernel.json", &j.render());
    println!("\n[artifact] {}", p);
    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(a′)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
