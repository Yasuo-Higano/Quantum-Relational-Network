//! v23.3 バケット別 K_boost vs K_A — R₂ は悪条件比か (第二十四期)
//!
//! v23.2 の発見: δ⟨K_A⟩ は巨大なバケット間相殺の微小残差 (d=1.5 で個別項の 0.5%、
//! d=2.5 で 0.03%)。K_boost の NN 形は K_A と 1% 一致 (v23.1) なので、R₂ の O(1)
//! 逸脱は「sub-% の演算子差 × 相殺増幅」の可能性 — つまり R₂ は悪条件比で、
//! Bisognano–Wichmann は演算子レベルで成立している、という再枠組みが立つ。
//! 本版はバケット別に δ⟨K_boost⟩ と δ⟨K_A⟩ を比較して裁く:
//!   比較 1: 各 NN バケットの応答比 ρ_b = δK_boost[b]/δK_A[b] (x-NN, y-NN, z-NN)
//!   比較 2: 増幅率 A(d) = Σ_b |δK_A[b]| / |δ⟨K_A⟩| と R₂ 逸脱の整合
//! 事前登録: (a) 全 NN バケットで ρ_b = 1 ± 0.03 (d=2.5 まで) = 「K_boost は
//!   バケット単位で正しく、R₂ 逸脱は相殺増幅」— R₂ を悪条件比として退役、
//!   BW の演算子レベル成立を宣言 → 次版で G_lattice 読み出し (プロファイル正当化) /
//!   (a′) 一部バケットが 3–10% ずれ = その成分の改良形を同定 /
//!   (b) O(10%) 超のずれ = 相殺増幅説の反証。
//! 装置ゲート: [G1/G2/G4] は v23.2 と同一。d ∈ {0.5, 1.5, 2.5}。

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
        hi: std::f64::consts::LN_2,
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
    println!("=== v23.3 バケット別 K_boost vs K_A — R₂ は悪条件比か (第二十四期) ===\n");
    println!("事前登録: (a) 全 NN バケット比 1±3% = 相殺増幅の確定・BW 演算子成立 /");
    println!("          (a′) 3-10% = 成分改良へ / (b) >10% = 増幅説の反証\n");
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

    // ---- K_A 行列 (DD) の構成: K = Σ_k κ_k v_k v_kᵀ ----
    println!("    K_A 行列構成中 ({} s)", t0.elapsed().as_secs());
    let mut kmat = vec![DD0; m * m];
    for k in 0..m {
        if kappa[k].hi.abs() < 1e-13 {
            continue;
        }
        for a in 0..m {
            let va = dd_mul(kappa[k], cv_dd[a + k * m]);
            if va.hi.abs() < 1e-40 {
                continue;
            }
            for b in a..m {
                kmat[a * m + b] = dd_add(kmat[a * m + b], dd_mul(va, cv_dd[b + k * m]));
            }
        }
    }
    for a in 0..m {
        for b in 0..a {
            kmat[a * m + b] = kmat[b * m + a];
        }
    }
    println!("    K_A 行列完了 ({} s)", t0.elapsed().as_secs());
    // sel の逆引き: a → (x, y, z) — sel は z 外・y 中・x 内 (x ∈ 0..n/2)
    let half = n / 2;
    let coord = |a: usize| -> (usize, usize, usize) {
        let x = a % half;
        let y = (a / half) % n;
        let z = a / (half * n);
        (x, y, z)
    };
    let wrapd = |p: usize, q: usize, nn: usize| -> usize {
        let d = if p > q { p - q } else { q - p };
        d.min(nn - d)
    };

    // ---- 各 d: δC_A (rank-2, 中心差分) → バケット分解 ----
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
    let labels = [
        "対角 (a=b)",
        "x-NN",
        "y-NN",
        "z-NN",
        "x-NNN (Δx=2)",
        "横 NNN",
        "長距離 (それ以遠)",
    ];
    let mut budget_rows: Vec<(f64, [f64; 7])> = Vec::new();
    let mut ratio_rows: Vec<(f64, [f64; 3], f64, f64)> = Vec::new();
    for &d in &[0.5f64, 1.5, 2.5] {
        let uf = mk(xplane - d);
        let wf = mk(xplane - d - 1.0);
        let mut u: Vec<Dd> = uf.iter().map(|&z| dd(z)).collect();
        let mut w: Vec<Dd> = wf.iter().map(|&z| dd(z)).collect();
        let mut nu = DD0;
        for i in 0..ns {
            nu = dd_add(nu, dd_mul(u[i], u[i]));
        }
        let inv = dd_div(dd(1.0), dd_sqrt(nu));
        for i in 0..ns {
            u[i] = dd_mul(u[i], inv);
        }
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
        // δC_A の rank-2 中心差分 (α = ±0.005): 奇部のみ
        let al = 0.005f64;
        let mut buckets = [DD0; 7];
        let mut total = DD0;
        // δC[a,b] = (C(+α) − C(−α))/2 の閉形式: 奇部 = s·(u_a Cw_b + w?…) — 明示に両符号評価
        let sel_u: Vec<Dd> = sel.iter().map(|&i| u[i]).collect();
        let sel_w: Vec<Dd> = sel.iter().map(|&i| w[i]).collect();
        let sel_cu: Vec<Dd> = sel.iter().map(|&i| cu[i]).collect();
        let sel_cw: Vec<Dd> = sel.iter().map(|&i| cwv[i]).collect();
        let mk_d = |sgn: f64, a: usize, b: usize| -> Dd {
            let aa = (al * sgn).cos() - 1.0;
            let ss = (al * sgn).sin();
            let du_a = dd_add(dd_mul_f(sel_u[a], aa), dd_mul_f(sel_w[a], ss));
            let dw_a = dd_sub(dd_mul_f(sel_w[a], aa), dd_mul_f(sel_u[a], ss));
            let du_b = dd_add(dd_mul_f(sel_u[b], aa), dd_mul_f(sel_w[b], ss));
            let dw_b = dd_sub(dd_mul_f(sel_w[b], aa), dd_mul_f(sel_u[b], ss));
            let mut t = dd_add(dd_mul(du_a, sel_cu[b]), dd_mul(dw_a, sel_cw[b]));
            t = dd_add(t, dd_add(dd_mul(du_b, sel_cu[a]), dd_mul(dw_b, sel_cw[a])));
            let m1 = dd_add(dd_mul(uu, du_b), dd_mul(uw, dw_b));
            let m2 = dd_add(dd_mul(uw, du_b), dd_mul(ww, dw_b));
            dd_add(t, dd_add(dd_mul(du_a, m1), dd_mul(dw_a, m2)))
        };
        let mut kb_buckets = [DD0; 7]; // K_boost (線形 NN) の同時集計
        let two_pi = 2.0 * std::f64::consts::PI;
        for a in 0..m {
            let (xa, ya, za) = coord(a);
            for b in 0..m {
                let kv = kmat[a * m + b];
                let (xb, yb, zb) = coord(b);
                let dx = if xa > xb { xa - xb } else { xb - xa };
                let dy = wrapd(ya, yb, n);
                let dz = wrapd(za, zb, n);
                let bucket = if a == b {
                    0
                } else if dx == 1 && dy == 0 && dz == 0 {
                    1
                } else if dx == 0 && dy == 1 && dz == 0 {
                    2
                } else if dx == 0 && dy == 0 && dz == 1 {
                    3
                } else if dx == 2 && dy == 0 && dz == 0 {
                    4
                } else if dx == 0 && (dy + dz == 2) {
                    5
                } else {
                    6
                };
                let kb: f64 = if bucket == 1 {
                    let xmin = xa.min(xb) as f64;
                    let xi = xplane - (xmin + 0.5);
                    if xi > 1e-9 {
                        two_pi * 0.5 * xi
                    } else {
                        0.0
                    }
                } else if bucket == 2 || bucket == 3 {
                    let xi = xplane - xa as f64;
                    if xi > 1e-9 {
                        let t = if bucket == 2 {
                            if xa % 2 == 0 {
                                0.5
                            } else {
                                -0.5
                            }
                        } else if (xa + ya) % 2 == 0 {
                            0.5
                        } else {
                            -0.5
                        };
                        two_pi * t * xi
                    } else {
                        0.0
                    }
                } else {
                    0.0
                };
                if kv.hi.abs() < 1e-24 && kb == 0.0 {
                    continue;
                }
                let dc = dd_mul_f(dd_sub(mk_d(1.0, a, b), mk_d(-1.0, a, b)), 0.5);
                if dc.hi.abs() < 1e-30 && dc.lo.abs() < 1e-30 {
                    continue;
                }
                let contrib = dd_mul(kv, dc);
                total = dd_add(total, contrib);
                buckets[bucket] = dd_add(buckets[bucket], contrib);
                if kb != 0.0 {
                    kb_buckets[bucket] = dd_add(kb_buckets[bucket], dd_mul_f(dc, kb));
                }
            }
        }
        // [G4] 恒等式: バケット総和 = Σ_k κ_k (v_kᵀ δC v_k) (v22.7 方式・固有基底)
        {
            let mut kexp = DD0;
            // δC の rank-2 構造を使う: v_kᵀ δC v_k は (v·δu 型の内積) の組で O(m)
            for k in 0..m {
                if kappa[k].hi.abs() < 1e-13 {
                    continue;
                }
                let mut acc = DD0;
                for (jj, &sgn) in [1.0f64, -1.0].iter().enumerate() {
                    let aa = (al * sgn).cos() - 1.0;
                    let ss = (al * sgn).sin();
                    // vᵀ δC v = 2 (vᵀdu)(vᵀCu) + 2 (vᵀdw)(vᵀCw) + (vᵀdu)m1 + (vᵀdw)m2 の構成
                    let mut vdu = DD0;
                    let mut vdw = DD0;
                    let mut vcu = DD0;
                    let mut vcw = DD0;
                    for a in 0..m {
                        let vk = cv_dd[a + k * m];
                        if vk.hi.abs() < 1e-30 {
                            continue;
                        }
                        let du_a = dd_add(dd_mul_f(sel_u[a], aa), dd_mul_f(sel_w[a], ss));
                        let dw_a = dd_sub(dd_mul_f(sel_w[a], aa), dd_mul_f(sel_u[a], ss));
                        vdu = dd_add(vdu, dd_mul(vk, du_a));
                        vdw = dd_add(vdw, dd_mul(vk, dw_a));
                        vcu = dd_add(vcu, dd_mul(vk, sel_cu[a]));
                        vcw = dd_add(vcw, dd_mul(vk, sel_cw[a]));
                    }
                    let mut t = dd_mul_f(dd_add(dd_mul(vdu, vcu), dd_mul(vdw, vcw)), 2.0);
                    let m1 = dd_add(dd_mul(uu, vdu), dd_mul(uw, vdw));
                    let m2 = dd_add(dd_mul(uw, vdu), dd_mul(ww, vdw));
                    t = dd_add(t, dd_add(dd_mul(vdu, m1), dd_mul(vdw, m2)));
                    let sg = if jj == 0 { 0.5 } else { -0.5 };
                    acc = dd_add(acc, dd_mul_f(t, sg));
                }
                kexp = dd_add(kexp, dd_mul(kappa[k], acc));
            }
            check(
                &format!(
                    "[G4] d={:.1} バケット総和 = 固有基底 δ⟨K⟩ (恒等式 ± 1e-8 相対)",
                    d
                ),
                ((total.hi - kexp.hi) / kexp.hi).abs() < 1e-8,
                format!("site 基底 {:+.8e} vs 固有基底 {:+.8e}", total.hi, kexp.hi),
            );
        }
        println!(
            "    d={:.1}: δ⟨K_A⟩ = {:+.6e} ({} s)",
            d,
            total.hi,
            t0.elapsed().as_secs()
        );
        let mut fr = [0.0f64; 7];
        let mut absum = 0.0f64;
        for i in 0..7 {
            fr[i] = buckets[i].hi / total.hi;
            absum += buckets[i].hi.abs();
        }
        let amp = absum / total.hi.abs();
        let mut kb_total = DD0;
        for i in 0..7 {
            kb_total = dd_add(kb_total, kb_buckets[i]);
        }
        let mut ratios = [0.0f64; 3];
        for (idx, bi) in [1usize, 2, 3].iter().enumerate() {
            ratios[idx] = kb_buckets[*bi].hi / buckets[*bi].hi;
            println!(
                "      {}: K_A {:+.4e} | K_boost {:+.4e} | 比 ρ = {:.5}",
                labels[*bi], buckets[*bi].hi, kb_buckets[*bi].hi, ratios[idx]
            );
        }
        println!(
            "      増幅率 A = {:.1}, K_boost 総和 = {:+.4e} (1/R₂ 再構成 = {:.5})",
            amp,
            kb_total.hi,
            total.hi / kb_total.hi
        );
        budget_rows.push((d, fr));
        ratio_rows.push((d, ratios, amp, total.hi / kb_total.hi));
    }

    // ---- 判定 ----
    let mut max_rho_dev = 0.0f64;
    for &(d, ratios, amp, r2rec) in &ratio_rows {
        for r in ratios {
            max_rho_dev = max_rho_dev.max((r - 1.0).abs());
        }
        println!(
            "    d={:.1}: ρ(x/y/z-NN) = {:.4}/{:.4}/{:.4}, A = {:.0}, 1/R₂(再構成) = {:.4}",
            d, ratios[0], ratios[1], ratios[2], amp, r2rec
        );
    }
    let branch_a = nfail == 0 && max_rho_dev <= 0.03;
    let branch_ap = nfail == 0 && !branch_a && max_rho_dev <= 0.10;
    println!(
        "\n[判定] {}",
        if nfail > 0 {
            "装置ゲート故障 — 記録"
        } else if branch_a {
            "事前登録 (a): 全 NN バケット比 = 1 ± 3% — R₂ は相殺増幅の悪条件比と確定、BW は演算子レベルで成立。次版で G_lattice 読み出し"
        } else if branch_ap {
            "事前登録 (a′): バケット比 3–10% のずれ — 当該成分の改良形を同定へ"
        } else {
            "事前登録 (b): O(10%) 超のずれ — 相殺増幅説の反証"
        }
    );
    println!("    max|ρ−1| = {:.4}", max_rho_dev);

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v23.3".into())),
        (
            "rows".into(),
            Json::Arr(
                ratio_rows
                    .iter()
                    .map(|&(d, ratios, amp, r2rec)| {
                        Json::Obj(vec![
                            ("d".into(), Json::Num(d)),
                            ("rho_x".into(), Json::Num(ratios[0])),
                            ("rho_y".into(), Json::Num(ratios[1])),
                            ("rho_z".into(), Json::Num(ratios[2])),
                            ("amp".into(), Json::Num(amp)),
                            ("r2_reconstructed_inv".into(), Json::Num(r2rec)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("max_rho_dev".into(), Json::Num(max_rho_dev)),
        ("branch_a".into(), Json::Bool(branch_a)),
    ]);
    let p = write_artifact("results/v233_bucketmatch.json", &j.render());
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
