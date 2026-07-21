//! v23.1 ブースト生成子の格子形 — 二境界 (Cardy–Tonni) 重みの検定と K_A プロファイル直測 (第二十四期)
//!
//! v22.7 で G 窓の残高は「深部まで K_A と整合するブースト生成子の格子形」1 点になった。
//! 本版はその同定を two-pronged で行う:
//!  [P] **直測**: DD で構成した K_A = Σ κ v vᵀ の x ボンド要素 (横平均) のプロファイル
//!      T(ξ) を測り、線形 (2πξ) と放物 (2π·ξ(2ℓ−ξ)/(2ℓ)) の形状と比較する。
//!  [R] **検定**: R₂ = δS/δ⟨K_w⟩ を重み 3 種 — [A] 線形 ξ (v19.2 凍結) /
//!      [CT] 放物 ξ(2ℓ−ξ)/(2ℓ) (BCFT の二境界重み — A は切断面と開放端の両方を持つ) /
//!      [P] 実測プロファイル T(ξ) の NN 形 — で d ∈ {0.5, 1.0, 1.5, 2.0, 2.5} を測る。
//!      δS は DD (v22.7 の床突破により d=2.5 までカナリア有効域)。
//!
//! 理論動機: 深部破綻 (d ≳ 2) は波束が開放端に近づく領域と一致する。CT 重みなら
//! 端で重みが有限に留まり、線形 ξ は深部を 2ℓ/(ℓ+x) 倍 (端で 2 倍) 過重み付けする。
//! 浅域の予測: d=1.0 で −8%・d=1.5 で −12.5% — v22.5 実測 (−10%/−13.5%) と整合。
//!
//! 装置ゲート (v22.7 と同一): [G1] DD 部分空間残差 < 1e-24, [G2] 非深部固有値一致,
//! [G3] R₁(DD, d=1.5) = 1 ± 1% (以後、各 d の R₁ がカナリア)。
//! 事前登録: (a) いずれかの重み形で、カナリア有効な全 d ≤ 2.5 が |R₂ − 1| ≤ 0.05
//!   = **ブースト生成子の格子形を同定 — G 窓開通** (次版で条件付き G_lattice 読み出し) /
//!   (a′) max|R₂−1| が線形形の半分以下 ∧ 深部の符号反転解消 = 部分同定 /
//!   (b) 改善なし = 二境界仮説の反証 (残る容疑 = K_A の非近接構造)。
//!   [P] プロファイル形状比較は判定でなく測定 (記録)。

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
    println!("=== v23.1 ブースト生成子の格子形 — CT 重み検定 + K_A 直測 (第二十四期) ===\n");
    println!("事前登録: (a) ある重み形で全有効 d が |R₂−1| ≤ 0.05 = 格子形同定・G 窓開通 /");
    println!("          (a′) 逸脱半減+符号反転解消 / (b) 二境界仮説の反証\n");
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

    // ---- [P] K_A の x ボンドプロファイル直測 (横平均) ----
    // K_A = Σ_k κ_k v_k v_kᵀ の NN 要素 K[a(x,y,z), a(x+1,y,z)] を (y,z) 平均。
    // sel の並び: a = x + (n/2)·y + (n/2)·n·z (z 外・y 中・x 内)。
    {
        let half = n / 2;
        let mut prof = vec![0.0f64; half - 1];
        for xb in 0..half - 1 {
            let mut acc = 0.0f64;
            let mut cnt = 0usize;
            for z in 0..n {
                for y in 0..n {
                    let a = xb + half * (y + n * z);
                    let b = (xb + 1) + half * (y + n * z);
                    let mut kel = DD0;
                    for k in 0..m {
                        if kappa[k].hi.abs() < 1e-13 {
                            continue;
                        }
                        kel = dd_add(
                            kel,
                            dd_mul(kappa[k], dd_mul(cv_dd[a + k * m], cv_dd[b + k * m])),
                        );
                    }
                    acc += kel.hi.abs();
                    cnt += 1;
                }
            }
            prof[xb] = acc / cnt as f64;
        }
        // 形状比較: ボンド ξ = (n/2 − 1) − xb + 0.5 … xplane − (xb+0.5)
        let ell = half as f64;
        println!("    [P] K_A ボンドプロファイル (横平均, ξ = 切断面からの距離):");
        let refi = half - 2; // ξ = 1 の参照ボンド
        for xb in 0..half - 1 {
            let xi = (half as f64) - 0.5 - (xb as f64 + 0.5);
            let lin = xi / ((half as f64) - 1.0);
            let ct =
                xi * (2.0 * ell - xi) / (2.0 * ell) / ((1.0) * (2.0 * ell - 1.0) / (2.0 * ell));
            println!(
                "      ξ={:.1}: T = {:.6} (T/T_ref = {:.4}) | 線形形 {:.4} | CT 形 {:.4}",
                xi,
                prof[xb],
                prof[xb] / prof[refi],
                lin,
                ct
            );
        }
    }

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
    let mut r2_of_d: Vec<(f64, f64, f64, f64)> = Vec::new();
    // K_w ボンド 3 変種: [A] 線形 ξ / [CT] 放物 ξ(2ℓ−ξ)/(2ℓ) / [P] 実測プロファイル形
    // (実測形は上の [P] の横平均 NN 要素を重みに正規化 — 自己整合検定)
    let ell = (n / 2) as f64;
    let wgt_lin = |xi: f64| xi;
    let wgt_ct = move |xi: f64| xi * (2.0 * ell - xi) / (2.0 * ell);
    let mut bonds3: Vec<Vec<(usize, usize, f64, f64)>> = vec![Vec::new(), Vec::new(), Vec::new()];
    {
        // 実測プロファイルを重み関数化 (ボンド ξ → T(ξ)/T(1) — サイト重みは補間)
        let half = n / 2;
        let mut prof = vec![0.0f64; half - 1];
        for xb in 0..half - 1 {
            let mut acc = 0.0f64;
            let mut cnt = 0usize;
            for z in 0..n {
                for y in 0..n {
                    let a = xb + half * (y + n * z);
                    let b = (xb + 1) + half * (y + n * z);
                    let mut kel = DD0;
                    for k in 0..m {
                        if kappa[k].hi.abs() < 1e-13 {
                            continue;
                        }
                        kel = dd_add(
                            kel,
                            dd_mul(kappa[k], dd_mul(cv_dd[a + k * m], cv_dd[b + k * m])),
                        );
                    }
                    acc += kel.hi.abs();
                    cnt += 1;
                }
            }
            prof[xb] = acc / cnt as f64;
        }
        let tref = prof[half - 2];
        let wgt_meas = |xi: f64| -> f64 {
            // ボンド ξ ∈ {1..ℓ−1} の実測値を線形補間 (範囲外は端値)
            let pos = (half as f64) - 1.0 - xi; // xb (実数)
            let p0 = pos.floor().max(0.0) as usize;
            let p1 = (p0 + 1).min(half - 2);
            let t = (pos - p0 as f64).clamp(0.0, 1.0);
            let v = prof[p0.min(half - 2)] * (1.0 - t) + prof[p1] * t;
            v / tref
        };
        for x in 0..n {
            for y in 0..n {
                for z in 0..n {
                    let i = idx3(x, y, z);
                    if x + 1 < n {
                        let xi = xplane - (x as f64 + 0.5);
                        if xi > 1e-9 {
                            let j = idx3(x + 1, y, z);
                            bonds3[0].push((i, j, 0.5, wgt_lin(xi)));
                            bonds3[1].push((i, j, 0.5, wgt_ct(xi)));
                            bonds3[2].push((i, j, 0.5, wgt_meas(xi)));
                        }
                    }
                    let xi_site = xplane - x as f64;
                    if xi_site > 1e-9 {
                        let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                        let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                        let jy = idx3(x, (y + 1) % n, z);
                        let jz = idx3(x, y, (z + 1) % n);
                        bonds3[0].push((i, jy, ey, wgt_lin(xi_site)));
                        bonds3[0].push((i, jz, ez, wgt_lin(xi_site)));
                        bonds3[1].push((i, jy, ey, wgt_ct(xi_site)));
                        bonds3[1].push((i, jz, ez, wgt_ct(xi_site)));
                        bonds3[2].push((i, jy, ey, wgt_meas(xi_site)));
                        bonds3[2].push((i, jz, ez, wgt_meas(xi_site)));
                    }
                }
            }
        }
    }
    let two_pi = 2.0 * std::f64::consts::PI;
    for &d in &[0.5f64, 1.0, 1.5, 2.0, 2.5] {
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
        let mut r2r3 = [[0.0f64; 2]; 3];
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
            // R₂ 3 変種: δ⟨K_w⟩ を rank-2 更新式で全ボンド DD 評価
            for vi in 0..3 {
                let mut dkb = [DD0; 2];
                for (jj, &sgn) in [1.0f64, -1.0].iter().enumerate() {
                    let aa = (al * sgn).cos() - 1.0;
                    let ss = (al * sgn).sin();
                    let mut acc = DD0;
                    for &(i, j, t, xi) in &bonds3[vi] {
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
                r2r3[vi][ii] = dd_div(ds, dkb_c).hi;
            }
            println!(
                "      d={:.1} α={:.3}: δS = {:+.6e}, 比 = {:.6}, R₂[A/CT/P] = {:.5}/{:.5}/{:.5} ({} s)",
                d,
                al,
                ds.hi,
                rr[ii],
                r2r3[0][ii],
                r2r3[1][ii],
                r2r3[2][ii],
                t0.elapsed().as_secs()
            );
        }
        let r1 = (4.0 * rr[0] - rr[1]) / 3.0;
        let mut r2v = [0.0f64; 3];
        for vi in 0..3 {
            r2v[vi] = (4.0 * r2r3[vi][0] - r2r3[vi][1]) / 3.0;
        }
        println!(
            "    N=12 d={:.1}: R₁(DD) = {:.5}, R₂[A] = {:.5}, R₂[CT] = {:.5}, R₂[P] = {:.5}",
            d, r1, r2v[0], r2v[1], r2v[2]
        );
        r1_of_d.push((d, r1));
        r2_of_d.push((d, r2v[0], r2v[1], r2v[2]));
    }

    // ---- 判定 ----
    let r15 = r1_of_d.iter().find(|r| r.0 == 1.5).unwrap().1;
    check(
        "[G3] d=1.5 の器械連続性 R₁(DD) = 1 ± 1%",
        (r15 - 1.0).abs() < 0.01,
        format!("R₁ = {:.5}", r15),
    );
    // カナリア有効点 (|R₁−1| < 1%) 上での各重みの max|R₂−1| と符号反転の有無
    let mut maxdev = [0.0f64; 3];
    let mut signflip = [false; 3];
    let mut nvalid = 0usize;
    for (i, &(d, r1)) in r1_of_d.iter().enumerate() {
        if (r1 - 1.0).abs() >= 0.01 {
            println!("    [除外] d={:.1}: カナリア無効 (R₁ = {:.5})", d, r1);
            continue;
        }
        nvalid += 1;
        let row = r2_of_d[i];
        for (vi, r2) in [row.1, row.2, row.3].iter().enumerate() {
            maxdev[vi] = maxdev[vi].max((r2 - 1.0).abs());
            if *r2 < 0.0 {
                signflip[vi] = true;
            }
        }
    }
    for (vi, name) in ["A 線形", "CT 放物", "P 実測"].iter().enumerate() {
        println!(
            "    [{}] max|R₂−1| = {:.4}{}",
            name,
            maxdev[vi],
            if signflip[vi] {
                " (符号反転あり)"
            } else {
                ""
            }
        );
    }
    let best_new = maxdev[1].min(maxdev[2]);
    let branch_a = nvalid >= 4 && (maxdev[1] <= 0.05 || maxdev[2] <= 0.05);
    let branch_ap = !branch_a
        && best_new <= 0.5 * maxdev[0]
        && ((!signflip[1] && maxdev[1] <= maxdev[2]) || (!signflip[2] && maxdev[2] < maxdev[1]));
    println!(
        "\n[判定] {}",
        if nfail > 0 {
            "装置ゲート故障 — 記録"
        } else if branch_a {
            "事前登録 (a): ブースト生成子の格子形を同定 — 全有効 d で |R₂−1| ≤ 0.05。G 窓開通 (次版で条件付き G_lattice)"
        } else if branch_ap {
            "事前登録 (a′): 部分同定 — 逸脱半減・深部符号反転の解消"
        } else {
            "事前登録 (b): 二境界仮説の反証 — 残る容疑は K_A の非近接構造"
        }
    );
    let rows: Vec<String> = r2_of_d
        .iter()
        .zip(r1_of_d.iter())
        .map(|(&(d, a, ct, pm), &(_, r1))| {
            format!(
                "d={:.1}: R₁={:.5} A={:.4} CT={:.4} P={:.4}",
                d, r1, a, ct, pm
            )
        })
        .collect();
    for r in &rows {
        println!("    {}", r);
    }

    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v23.1".into())),
        (
            "rows".into(),
            Json::Arr(
                r2_of_d
                    .iter()
                    .zip(r1_of_d.iter())
                    .map(|(&(d, a, ct, pm), &(_, r1))| {
                        Json::Obj(vec![
                            ("d".into(), Json::Num(d)),
                            ("r1".into(), Json::Num(r1)),
                            ("r2_lin".into(), Json::Num(a)),
                            ("r2_ct".into(), Json::Num(ct)),
                            ("r2_meas".into(), Json::Num(pm)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("maxdev_lin".into(), Json::Num(maxdev[0])),
        ("maxdev_ct".into(), Json::Num(maxdev[1])),
        ("maxdev_meas".into(), Json::Num(maxdev[2])),
        ("branch_a".into(), Json::Bool(branch_a)),
    ]);
    let p = write_artifact("results/v231_boostform.json", &j.render());
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
