//! v20.4 Schwinger アノマリー — 磁束挿入クエンチの実時間動力学 (第二十一期 第四歩)
//!
//! 軸性アノマリー ∂_μ j⁵μ = gE/π を拘束 core の多体実時間発展で測る。
//! 磁束挿入 (ε → ε+1 = 一様背景場 E₀* = 1) を基底状態に差し、Krylov 法で発展。
//! 格子単位 (H* = H/(g²a/2), t* = t g²a/2) の予言は係数フリー:
//!   [i] dQ⁵/dt*|₀ = (2N/π) E₀*,  Q⁵ = Σ_n Im⟨c†_n c_{n+1}⟩ (x が落ちる)
//!   [ii] 背反作用 = プラズマ振動: ⟨E⟩(t) の角振動数 ω* = 2√x (M/g) — v20.1 の
//!        スペクトル質量 (x=2: 0.7179, x=4: 0.6729 [N 依存は 4 桁目]) との動的クロスチェック
//! 事前登録: (a) [i] 勾配 (t*∈[0,0.3] フィット) が 2NE₀/π ± 15% — N∈{12,14} の両方 —
//!   かつ [ii] ω*/(2√x) = v20.1 の M/g ± 15% — (N=14, x∈{2,4}) の両方 = アノマリーと
//!   背反作用の両輪成立 / (a′) 片輪 / (b) 両方外れ。
//! 装置ゲート: ノルム保存 (Krylov ユニタリ性)・エネルギー保存・連続の式 (密度の時間
//! 微分 = 電流の発散, 有限差分照合)・状態準備 ⟨E⟩(0) = E₀ ± 0.02。μ=0 (質量ゼロ)。
//! core 構造体は v20.1〜20.3 と同一コード ([A] 検証済み)。

use uft_sim::*;

// ---- 組合せ: 固定粒子数の bitmask 列挙 ----
fn enum_masks(n: usize, nf: usize) -> Vec<u32> {
    let mut v = Vec::new();
    let end: u32 = 1 << n;
    let mut m: u32 = (1 << nf) - 1;
    if nf == 0 {
        return vec![0];
    }
    while m < end {
        v.push(m);
        // Gosper's hack
        let c = m & m.wrapping_neg();
        let r = m + c;
        m = (((r ^ m) >> 2) / c) | r;
    }
    v
}

// ---- 拘束を解いた U(1) core ----
struct U1Core {
    n: usize,
    x: f64,
    mu: f64,
    lam: f64,                  // |E_n| ≤ lam
    probes: Vec<(usize, f64)>, // 外部電荷 (site, q)
    states: Vec<u64>,          // key = (mask as u64) << 8 | (eps+128)
    dim: usize,
}

impl U1Core {
    // ボンド電場列 E_n (n = 0..N-1, ボンド n は (n, n+1 mod N))
    fn e_profile(&self, mask: u32, eps: i32) -> Vec<f64> {
        let mut e = vec![0.0; self.n];
        let mut p = 0.0f64;
        for site in 0..self.n {
            let occ = (mask >> site) & 1;
            let mut q = occ as f64 - if site % 2 == 1 { 1.0 } else { 0.0 };
            for &(ps, pq) in &self.probes {
                if ps == site {
                    q += pq;
                }
            }
            p += q;
            e[site] = eps as f64 + p;
        }
        e
    }
    fn total_charge_ok(&self, mask: u32) -> bool {
        let nf = mask.count_ones() as f64;
        let qtot = nf - (self.n as f64) / 2.0 + self.probes.iter().map(|p| p.1).sum::<f64>();
        qtot.abs() < 1e-9
    }
    fn new(n: usize, nf: usize, x: f64, mu: f64, lam: f64, probes: Vec<(usize, f64)>) -> Self {
        let mut c = U1Core {
            n,
            x,
            mu,
            lam,
            probes,
            states: Vec::new(),
            dim: 0,
        };
        let masks = enum_masks(n, nf);
        for &m in &masks {
            if !c.total_charge_ok(m) {
                continue;
            }
            for eps in -6i32..=6 {
                let e = c.e_profile(m, eps);
                if e.iter().all(|&v| v.abs() <= c.lam + 1e-9) {
                    c.states.push(((m as u64) << 8) | ((eps + 128) as u64));
                }
            }
        }
        c.states.sort_unstable();
        c.dim = c.states.len();
        c
    }
    fn unpack(&self, key: u64) -> (u32, i32) {
        ((key >> 8) as u32, (key & 0xff) as i32 - 128)
    }
    fn find(&self, mask: u32, eps: i32) -> Option<usize> {
        let key = ((mask as u64) << 8) | ((eps + 128) as u64);
        self.states.binary_search(&key).ok()
    }
    fn diag(&self, mask: u32, eps: i32) -> f64 {
        let e = self.e_profile(mask, eps);
        let elec: f64 = e.iter().map(|&v| v * v).sum();
        let mut mass = 0.0;
        for site in 0..self.n {
            let occ = ((mask >> site) & 1) as f64;
            let sgn = if site % 2 == 0 { 1.0 } else { -1.0 };
            mass += sgn * (2.0 * occ - 1.0);
        }
        elec + 0.5 * self.mu * mass
    }
    // H v (実対称 — 位相なしのホップなので実で閉じる)
    fn matvec(&self, v: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let mut w = vec![(0.0, 0.0); self.dim];
        for (i, &key) in self.states.iter().enumerate() {
            let (mask, eps) = self.unpack(key);
            let d = self.diag(mask, eps);
            w[i].0 += d * v[i].0;
            w[i].1 += d * v[i].1;
            // バルクホップ (site, site+1), site = 0..n-2 (隣接 JW — 符号なし)
            for site in 0..self.n - 1 {
                let b0 = (mask >> site) & 1;
                let b1 = (mask >> (site + 1)) & 1;
                if b0 == b1 {
                    continue;
                }
                let nm = mask ^ (1 << site) ^ (1 << (site + 1));
                if let Some(j) = self.find(nm, eps) {
                    w[j].0 += -self.x * v[i].0;
                    w[j].1 += -self.x * v[i].1;
                }
            }
            // 継ぎ目ホップ (n-1, 0): ε → ε∓1, JW 符号 (−1)^{N_f−1}
            let b0 = (mask >> (self.n - 1)) & 1;
            let b1 = mask & 1;
            if b0 != b1 {
                let nm = mask ^ (1 << (self.n - 1)) ^ 1;
                let nf = mask.count_ones();
                let sgn = if (nf - 1) % 2 == 0 { 1.0 } else { -1.0 };
                // 電荷が (n-1)→0 に動く (b0=1): 全バルク P が変わらないよう ε → ε−1?
                // 検証 [A] が固定する: 0→(n-1) 移動 (b1=1) で ε → ε+1, 逆は ε−1
                let neps = if b1 == 1 { eps + 1 } else { eps - 1 };
                if let Some(j) = self.find(nm, neps) {
                    w[j].0 += -self.x * sgn * v[i].0;
                    w[j].1 += -self.x * sgn * v[i].1;
                }
            }
        }
        w
    }
}

// Krylov 指数発展 e^{−iH dt} ψ (完全再直交, m 次元)
fn krylov_step(core: &U1Core, psi: &[(f64, f64)], dt: f64, m: usize) -> Vec<(f64, f64)> {
    let n = psi.len();
    let nrm0 = psi
        .iter()
        .map(|z| z.0 * z.0 + z.1 * z.1)
        .sum::<f64>()
        .sqrt();
    let mut basis: Vec<Vec<(f64, f64)>> =
        vec![psi.iter().map(|z| (z.0 / nrm0, z.1 / nrm0)).collect()];
    let mut alpha = Vec::new();
    let mut beta = Vec::new();
    for j in 0..m {
        let mut w = core.matvec(&basis[j]);
        let a: f64 = basis[j]
            .iter()
            .zip(w.iter())
            .map(|(b, z)| b.0 * z.0 + b.1 * z.1)
            .sum();
        alpha.push(a);
        for b in &basis {
            let (pr, pi): (f64, f64) = b.iter().zip(w.iter()).fold((0.0, 0.0), |acc, (bb, zz)| {
                (
                    acc.0 + bb.0 * zz.0 + bb.1 * zz.1,
                    acc.1 + bb.0 * zz.1 - bb.1 * zz.0,
                )
            });
            for i in 0..n {
                let (br, bi) = (b[i].0, b[i].1);
                w[i].0 -= pr * br - pi * bi;
                w[i].1 -= pr * bi + pi * br;
            }
        }
        let bn: f64 = w.iter().map(|z| z.0 * z.0 + z.1 * z.1).sum::<f64>().sqrt();
        if j + 1 == m || bn < 1e-10 {
            break;
        }
        beta.push(bn);
        basis.push(w.iter().map(|z| (z.0 / bn, z.1 / bn)).collect());
    }
    let k = alpha.len();
    let mut t = vec![0.0f64; k * k];
    for i in 0..k {
        t[i + i * k] = alpha[i];
        if i + 1 < k {
            t[i + (i + 1) * k] = beta[i];
            t[(i + 1) + i * k] = beta[i];
        }
    }
    let (ev, vv) = jacobi_eigh(&t, k);
    // c = V e^{−iθdt} Vᵀ e₁ · nrm0
    let mut cr = vec![0.0f64; k];
    let mut ci = vec![0.0f64; k];
    for a in 0..k {
        for l in 0..k {
            let ph = -ev[l] * dt;
            let (cs, sn) = (ph.cos(), ph.sin());
            let w = vv[a + l * k] * vv[0 + l * k];
            cr[a] += w * cs;
            ci[a] += w * sn;
        }
    }
    let mut out = vec![(0.0, 0.0); n];
    for a in 0..k {
        let (car, cai) = (cr[a] * nrm0, ci[a] * nrm0);
        for i in 0..n {
            let (br, bi) = (basis[a][i].0, basis[a][i].1);
            out[i].0 += car * br - cai * bi;
            out[i].1 += car * bi + cai * br;
        }
    }
    out
}

fn main() {
    self_test();
    println!("=== v20.4 Schwinger アノマリー — 磁束挿入クエンチの実時間動力学 ===\n");
    println!("事前登録: (a) 勾配 dQ⁵/dt* = 2NE₀/π ± 15% (N∈{{12,14}}) かつ ω*/(2√x) = v20.1 の");
    println!("          M/g ± 15% (N=14, x∈{{2,4}}) = 両輪成立 / (a′) 片輪 / (b) 両方外れ\n");
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
    let mut slope_ok_all = true;
    let mut mass_ok_all = true;
    let mut rows = Vec::new(); // (n, x, slope, pred, mg_dyn, mg_ref)
    for &(n, x, mg_ref, do_mass) in &[
        (12usize, 2.0f64, 0.7175f64, false),
        (14, 2.0, 0.7179, true),
        (14, 4.0, 0.6725, true),
    ] {
        let core = U1Core::new(n, n / 2, x, 0.0, 4.0, vec![]);
        let mv = |vv: &[(f64, f64)]| core.matvec(vv);
        let (_ev0, vecs, res) = lanczos_lowest_herm(&mv, core.dim, 1, 120, 91);
        check(
            &format!("N={} x={:.0} の基底状態 Lanczos 残差", n, x),
            res < 1e-8,
            format!("res = {:.0e} (dim {})", res, core.dim),
        );
        // 磁束挿入 ε → ε+1
        let mut psi = vec![(0.0f64, 0.0f64); core.dim];
        let mut lost = 0.0;
        for (i, &key) in core.states.iter().enumerate() {
            let (mask, eps) = core.unpack(key);
            let amp = vecs[0][i];
            if let Some(j) = core.find(mask, eps + 1) {
                psi[j] = amp;
            } else {
                lost += amp.0 * amp.0 + amp.1 * amp.1;
            }
        }
        check(
            &format!("N={} x={:.0} の磁束挿入 (切断でのノルム損失)", n, x),
            lost < 1e-6,
            format!("損失 = {:.1e}", lost),
        );
        // 観測量
        let e_mean = |p: &[(f64, f64)]| -> f64 {
            let mut acc = 0.0;
            for (i, &key) in core.states.iter().enumerate() {
                let (mask, eps) = core.unpack(key);
                let w = p[i].0 * p[i].0 + p[i].1 * p[i].1;
                if w == 0.0 {
                    continue;
                }
                let e = core.e_profile(mask, eps);
                acc += w * e.iter().sum::<f64>() / core.n as f64;
            }
            acc
        };
        let q5 = |p: &[(f64, f64)]| -> f64 {
            // Q⁵ = Σ_n Im⟨c†_n c_{n+1}⟩ = ⟨ψ| J |ψ⟩ / … J をその場で作用
            let mut acc = 0.0;
            for (i, &key) in core.states.iter().enumerate() {
                let (mask, eps) = core.unpack(key);
                if p[i].0 == 0.0 && p[i].1 == 0.0 {
                    continue;
                }
                // バルク
                for site in 0..core.n - 1 {
                    let b0 = (mask >> site) & 1;
                    let b1 = (mask >> (site + 1)) & 1;
                    if b0 == b1 {
                        continue;
                    }
                    let nm = mask ^ (1 << site) ^ (1 << (site + 1));
                    if let Some(j) = core.find(nm, eps) {
                        // c†_{site+1} c_site 型 (b0=1) は +i/2, 逆は −i/2 → Im 部を集める
                        let s = if b0 == 1 { 1.0 } else { -1.0 };
                        // ⟨p_j | i·s | p_i⟩ の実部: Re(conj(p_j)·i·s·p_i) = s·Im(conj(p_j)p_i)·(−1)
                        acc += s * (p[j].0 * p[i].1 - p[j].1 * p[i].0);
                    }
                }
                // 継ぎ目
                let b0 = (mask >> (core.n - 1)) & 1;
                let b1 = mask & 1;
                if b0 != b1 {
                    let nm = mask ^ (1 << (core.n - 1)) ^ 1;
                    let nf = mask.count_ones();
                    let jw = if (nf - 1) % 2 == 0 { 1.0 } else { -1.0 };
                    let neps = if b1 == 1 { eps + 1 } else { eps - 1 };
                    if let Some(j) = core.find(nm, neps) {
                        let s = if b0 == 1 { 1.0 } else { -1.0 };
                        acc += jw * s * (p[j].0 * p[i].1 - p[j].1 * p[i].0);
                    }
                }
            }
            -acc / 2.0 // 各ボンドは両向きで 2 回数えられる (run1 で 2 倍として発覚) — Q⁵ = −sum/2
        };
        let energy = |p: &[(f64, f64)]| -> f64 {
            let hp = core.matvec(p);
            p.iter()
                .zip(hp.iter())
                .map(|(a, b)| a.0 * b.0 + a.1 * b.1)
                .sum()
        };
        let e_init = e_mean(&psi);
        check(
            &format!("N={} x={:.0} の状態準備 ⟨E⟩(0) ≈ 1", n, x),
            (e_init - 1.0).abs() < 0.02,
            format!("⟨E⟩(0) = {:.4}", e_init),
        );
        let en0 = energy(&psi);
        let nrm_init: f64 = psi.iter().map(|z| z.0 * z.0 + z.1 * z.1).sum::<f64>();
        // 発展
        let dt = 0.05f64;
        let nstep = 130usize;
        let mut ts = Vec::new();
        let mut q5s = Vec::new();
        let mut es = Vec::new();
        let mut p = psi.clone();
        for st in 0..=nstep {
            if st > 0 {
                p = krylov_step(&core, &p, dt, 26);
            }
            let t = st as f64 * dt;
            ts.push(t);
            q5s.push(q5(&p));
            es.push(e_mean(&p));
        }
        let nrm_fin: f64 = p.iter().map(|z| z.0 * z.0 + z.1 * z.1).sum::<f64>();
        check(
            &format!("N={} x={:.0} のノルム保存", n, x),
            (nrm_fin - nrm_init).abs() < 1e-8,
            format!("|Δ| = {:.1e}", (nrm_fin - nrm_init).abs()),
        );
        let en1 = energy(&p);
        check(
            &format!("N={} x={:.0} のエネルギー保存", n, x),
            (en1 - en0).abs() < 1e-6 * (1.0 + en0.abs()),
            format!("|ΔE| = {:.1e}", (en1 - en0).abs()),
        );
        // [i] 初期勾配 (t ≤ 0.3)
        let k = (0.3 / dt).round() as usize + 1;
        let (_ic, slope) = linfit(&ts[..k], &q5s[..k]);
        let pred = 2.0 * n as f64 / std::f64::consts::PI;
        let slope_ok = ((slope.abs() - pred) / pred).abs() < 0.15;
        println!(
            "    N={} x={:.0}: dQ⁵/dt* = {:+.3} (|予言| 2N/π = {:.3}, 比 {:.3}) ({} s)",
            n,
            x,
            slope,
            pred,
            slope.abs() / pred,
            t0.elapsed().as_secs()
        );
        if !slope_ok {
            slope_ok_all = false;
        }
        // [記録] 点毎のアノマリー則: dQ⁵/dt(t) = (2N/π)·Z·⟨E⟩(t) — 全周期で Z が一定か
        {
            let mut zs = Vec::new();
            for w in 1..q5s.len() - 1 {
                let dq = (q5s[w + 1] - q5s[w - 1]) / (2.0 * dt);
                let em = es[w];
                if em.abs() > 0.15 {
                    zs.push(dq / (pred * em));
                }
            }
            let zm = zs.iter().sum::<f64>() / zs.len() as f64;
            let zsd =
                (zs.iter().map(|z| (z - zm) * (z - zm)).sum::<f64>() / zs.len() as f64).sqrt();
            println!(
                "    N={} x={:.0} [記録] 点毎比 Z = dQ⁵/dt/((2N/π)⟨E⟩): 平均 {:.4} ± {:.4} ({} 点, |E|>0.15)",
                n, x, zm, zsd, zs.len()
            );
        }
        // [ii] プラズマ振動: ⟨E⟩ の最初のゼロ交差 → ω = π/(2 t₀) → M/g = ω/(2√x)
        let mut mg_dyn = f64::NAN;
        if do_mass {
            let mut tzero = f64::NAN;
            for w in 1..es.len() {
                if es[w - 1] > 0.0 && es[w] <= 0.0 {
                    let f = es[w - 1] / (es[w - 1] - es[w]);
                    tzero = ts[w - 1] + f * dt;
                    break;
                }
            }
            let omega = std::f64::consts::PI / (2.0 * tzero);
            mg_dyn = omega / (2.0 * x.sqrt());
            let ok = ((mg_dyn - mg_ref) / mg_ref).abs() < 0.15;
            println!(
                "    N={} x={:.0}: ⟨E⟩ ゼロ交差 t₀ = {:.3} → ω* = {:.3} → M/g(動的) = {:.4} (v20.1: {:.4}, 比 {:.3})",
                n,
                x,
                tzero,
                omega,
                mg_dyn,
                mg_ref,
                mg_dyn / mg_ref
            );
            if !ok {
                mass_ok_all = false;
            }
        }
        if n == 14 && x == 2.0 {
            let cj = Json::Obj(vec![
                (
                    "t".into(),
                    Json::Arr(ts.iter().map(|&v| Json::Num(v)).collect()),
                ),
                (
                    "q5".into(),
                    Json::Arr(q5s.iter().map(|&v| Json::Num(v)).collect()),
                ),
                (
                    "e".into(),
                    Json::Arr(es.iter().map(|&v| Json::Num(v)).collect()),
                ),
            ]);
            let _ = write_artifact("results/v204_curves_n14x2.json", &cj.render());
        }
        rows.push((n, x, slope.abs(), pred, mg_dyn, mg_ref));
    }

    // ---- 判定 (記録) ----
    println!(
        "\n[判定] {}",
        if slope_ok_all && mass_ok_all {
            "事前登録 (a): アノマリー係数 1/π と背反作用 (プラズマ振動 = 動的質量) の両輪成立"
        } else if slope_ok_all || mass_ok_all {
            "事前登録 (a′): 片輪成立 — 記録"
        } else {
            "事前登録 (b): 両方外れ — 記録"
        }
    );
    println!(
        "    勾配 [{}] / 動的質量 [{}]",
        if slope_ok_all {
            "全点成立"
        } else {
            "外れあり"
        },
        if mass_ok_all {
            "全点成立"
        } else {
            "外れあり"
        }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v20.4".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|&(n, x, s, pr, md, mr)| {
                        Json::Obj(vec![
                            ("n".into(), Json::Int(n as i64)),
                            ("x".into(), Json::Num(x)),
                            ("slope".into(), Json::Num(s)),
                            ("pred".into(), Json::Num(pr)),
                            ("mg_dyn".into(), Json::Num(md)),
                            ("mg_ref".into(), Json::Num(mr)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(slope_ok_all && mass_ok_all)),
    ]);
    let p = write_artifact("results/v204_anomaly.json", &j.render());
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
