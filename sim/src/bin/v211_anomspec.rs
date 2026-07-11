//! v21.1 アノマリー係数のスペクトル読み出し — 台帳の数が動力学から出る (第二十二期開幕)
//!
//! PROMPT/3 §3 の橋「anomaly-search の代数的一意性を dynamical core に接続する」の
//! 1+1D 誠実形。多フレーバー U(1) core (共有電場・フレーバー電荷 {q_f}) の
//! Schwinger ボソン質量は M² = (Σ_f q_f²) g²/π — アノマリー係数 Σq² が
//! (i) 台帳側では代数 (anomaly 方程式) の係数として、(ii) core 側では動力学
//! スペクトルの数として、同じ値で現れる。荷電集合 {1}, {1,1}, {2}, {1,2} で
//! 質量比 1 : √2 : 2 : √5 を検証する。
//!
//! 装置 (v2 = スペクトル関数法):
//!   - 基底 = (mask_f1, mask_f2, ε∈ℤ)。E_n = ε + Σ_f q_f P_n^f (per-bond 切断 |E|≤Λ)。
//!     継ぎ目ホップはフレーバー f の電荷分 ε → ε ∓ q_f (v20.1 の一般化)。
//!   - ボソンの同定と質量: シード (Ê_mean − ⟨E⟩)|0⟩ の 3 項 Lanczos (再直交なし) の
//!     Ritz 値と重み |s_j|² = Ê のスペクトル関数 — 「E が結合する状態の質量」を直接測る。
//!     ボソン極 = 重み ≥ 0.1 の最低 Ritz 値。フレーバー擬ゼロモードは E と結合しない。
//!   - 開発記録 (run1, 保存): 励起状態の素朴 Lanczos (k=8, m=170) は非収束 (res ~1e0)、
//!     {2} の E₁ はフラックス状態でボソンでない、cross-N 回帰帯 (±2%) は根拠曖昧で
//!     2.002% で落ちた — 全てゲートが捕捉。v2 で根治。
//! 装置ゲート: 基底状態残差 < 1e-7・ボソン極の重み ≥ 0.1・{1} の連続外挿が
//!   1/√π ± 10% (v20.1 と同型の自己完結アンカー)。
//! 事前登録: 連続外挿 (a₀ + a₁/√x, x∈{1,2,3,4}, N=10) の質量比が
//!   (a) √2, 2, √5 の全てに ±7% = アノマリー係数がスペクトルから読めた (橋の 1+1D 形成立) /
//!   (a′) 2 つ以上 = 部分成立 / (b) それ以外。

use uft_sim::*;

fn enum_masks(n: usize, nf: usize) -> Vec<u32> {
    let mut v = Vec::new();
    let end: u32 = 1 << n;
    if nf == 0 {
        return vec![0];
    }
    let mut m: u32 = (1 << nf) - 1;
    while m < end {
        v.push(m);
        let c = m & m.wrapping_neg();
        let r = m + c;
        m = (((r ^ m) >> 2) / c) | r;
    }
    v
}

struct MultiU1 {
    n: usize,
    x: f64,
    lam: f64,
    q: Vec<i32>,      // フレーバー電荷 (1 or 2 フレーバー)
    states: Vec<u64>, // key = m1 << (n+8) | m2 << 8 | (ε+128)
    dim: usize,
}

impl MultiU1 {
    fn e_profile(&self, m1: u32, m2: u32, eps: i32) -> Vec<f64> {
        let mut e = vec![0.0; self.n];
        let mut p = 0i32;
        for site in 0..self.n {
            let bg = if site % 2 == 1 { 1 } else { 0 };
            p += self.q[0] * (((m1 >> site) & 1) as i32 - bg);
            if self.q.len() > 1 {
                p += self.q[1] * (((m2 >> site) & 1) as i32 - bg);
            }
            e[site] = eps as f64 + p as f64;
        }
        e
    }
    fn new(n: usize, x: f64, lam: f64, q: Vec<i32>) -> Self {
        let mut c = MultiU1 {
            n,
            x,
            lam,
            q,
            states: Vec::new(),
            dim: 0,
        };
        let masks = enum_masks(n, n / 2);
        let one_flavor = c.q.len() == 1;
        let m2set: Vec<u32> = if one_flavor { vec![0] } else { masks.clone() };
        for &m1 in &masks {
            for &m2 in &m2set {
                for eps in -8i32..=8 {
                    let e = c.e_profile(m1, m2, eps);
                    if e.iter().all(|&v| v.abs() <= c.lam + 1e-9) {
                        c.states.push(
                            ((m1 as u64) << (c.n + 8)) | ((m2 as u64) << 8) | ((eps + 128) as u64),
                        );
                    }
                }
            }
        }
        c.states.sort_unstable();
        c.dim = c.states.len();
        c
    }
    fn unpack(&self, key: u64) -> (u32, u32, i32) {
        (
            (key >> (self.n + 8)) as u32,
            ((key >> 8) & ((1 << self.n) - 1)) as u32,
            (key & 0xff) as i32 - 128,
        )
    }
    fn find(&self, m1: u32, m2: u32, eps: i32) -> Option<usize> {
        let key = ((m1 as u64) << (self.n + 8)) | ((m2 as u64) << 8) | ((eps + 128) as u64);
        self.states.binary_search(&key).ok()
    }
    fn diag(&self, m1: u32, m2: u32, eps: i32) -> f64 {
        self.e_profile(m1, m2, eps).iter().map(|&v| v * v).sum()
    }
    fn e_mean_of(&self, m1: u32, m2: u32, eps: i32) -> f64 {
        self.e_profile(m1, m2, eps).iter().sum::<f64>() / self.n as f64
    }
    fn matvec(&self, v: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let mut w = vec![(0.0, 0.0); self.dim];
        let two_fl = self.q.len() > 1;
        for (i, &key) in self.states.iter().enumerate() {
            if v[i].0 == 0.0 && v[i].1 == 0.0 {
                continue;
            }
            let (m1, m2, eps) = self.unpack(key);
            let d = self.diag(m1, m2, eps);
            w[i].0 += d * v[i].0;
            w[i].1 += d * v[i].1;
            // フレーバーごとのホップ
            for fl in 0..self.q.len() {
                let mask = if fl == 0 { m1 } else { m2 };
                let qf = self.q[fl];
                // バルク
                for site in 0..self.n - 1 {
                    let b0 = (mask >> site) & 1;
                    let b1 = (mask >> (site + 1)) & 1;
                    if b0 == b1 {
                        continue;
                    }
                    let nm = mask ^ (1 << site) ^ (1 << (site + 1));
                    let (n1, n2) = if fl == 0 { (nm, m2) } else { (m1, nm) };
                    if let Some(j) = self.find(n1, n2, eps) {
                        w[j].0 += -self.x * v[i].0;
                        w[j].1 += -self.x * v[i].1;
                    }
                }
                // 継ぎ目: ε → ε ± q_f (JW 符号は自フレーバーの N_f)
                let b0 = (mask >> (self.n - 1)) & 1;
                let b1 = mask & 1;
                if b0 != b1 {
                    let nm = mask ^ (1 << (self.n - 1)) ^ 1;
                    let nf = mask.count_ones();
                    let sgn = if (nf - 1) % 2 == 0 { 1.0 } else { -1.0 };
                    let neps = if b1 == 1 { eps + qf } else { eps - qf };
                    let (n1, n2) = if fl == 0 { (nm, m2) } else { (m1, nm) };
                    if neps >= -8 && neps <= 8 {
                        if let Some(j) = self.find(n1, n2, neps) {
                            w[j].0 += -self.x * sgn * v[i].0;
                            w[j].1 += -self.x * sgn * v[i].1;
                        }
                    }
                }
            }
            let _ = two_fl;
        }
        w
    }
}

// 再開始 Lanczos (基底状態, v20.5 と同形)
fn ground(
    core: &MultiU1,
    m: usize,
    rounds: usize,
    tol: f64,
    seed: u64,
) -> (f64, Vec<(f64, f64)>, f64) {
    let n = core.dim;
    let mut rng = Rng::new(seed);
    let mut v: Vec<(f64, f64)> = (0..n).map(|_| (rng.gauss(), rng.gauss())).collect();
    let mut ev0 = 0.0;
    let mut res = f64::INFINITY;
    for _round in 0..rounds {
        let nrm = v.iter().map(|z| z.0 * z.0 + z.1 * z.1).sum::<f64>().sqrt();
        for z in v.iter_mut() {
            z.0 /= nrm;
            z.1 /= nrm;
        }
        let mut basis = vec![v.clone()];
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
            for _ in 0..2 {
                for b in &basis {
                    let (pr, pi): (f64, f64) =
                        b.iter().zip(w.iter()).fold((0.0, 0.0), |ac, (bb, zz)| {
                            (
                                ac.0 + bb.0 * zz.0 + bb.1 * zz.1,
                                ac.1 + bb.0 * zz.1 - bb.1 * zz.0,
                            )
                        });
                    for i in 0..n {
                        let (br, bi) = (b[i].0, b[i].1);
                        w[i].0 -= pr * br - pi * bi;
                        w[i].1 -= pr * bi + pi * br;
                    }
                }
            }
            let bn: f64 = w.iter().map(|z| z.0 * z.0 + z.1 * z.1).sum::<f64>().sqrt();
            if j + 1 == m || bn < 1e-12 {
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
        let (evs, vv) = jacobi_eigh(&t, k);
        ev0 = evs[0];
        let mut nv = vec![(0.0f64, 0.0f64); n];
        for a in 0..k {
            let cc = vv[a];
            for i in 0..n {
                nv[i].0 += cc * basis[a][i].0;
                nv[i].1 += cc * basis[a][i].1;
            }
        }
        let hv = core.matvec(&nv);
        res = hv
            .iter()
            .zip(nv.iter())
            .map(|(h, z)| {
                let dr = h.0 - ev0 * z.0;
                let di = h.1 - ev0 * z.1;
                dr * dr + di * di
            })
            .sum::<f64>()
            .sqrt();
        v = nv;
        if res < tol {
            break;
        }
    }
    (ev0, v, res)
}

// スペクトル関数: シード (Ê_mean − ⟨E⟩)|0⟩ の 3 項 Lanczos (再直交なし) →
// (ボソン極の ΔE, 極の重み)
fn boson_pole(core: &MultiU1, gs: &[(f64, f64)], e0: f64, m: usize) -> (f64, f64) {
    let n = core.dim;
    // シード
    let emean: f64 = (0..n)
        .map(|i| {
            let (m1, m2, eps) = core.unpack(core.states[i]);
            (gs[i].0 * gs[i].0 + gs[i].1 * gs[i].1) * core.e_mean_of(m1, m2, eps)
        })
        .sum();
    let mut v0: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let (m1, m2, eps) = core.unpack(core.states[i]);
            let w = core.e_mean_of(m1, m2, eps) - emean;
            (gs[i].0 * w, gs[i].1 * w)
        })
        .collect();
    let nrm = v0.iter().map(|z| z.0 * z.0 + z.1 * z.1).sum::<f64>().sqrt();
    for z in v0.iter_mut() {
        z.0 /= nrm;
        z.1 /= nrm;
    }
    // 3 項 Lanczos (基底 2 本のみ保持)
    let mut alpha = Vec::new();
    let mut beta = Vec::new();
    let mut vp: Vec<(f64, f64)> = vec![(0.0, 0.0); n];
    let mut vc = v0;
    for j in 0..m {
        let mut w = core.matvec(&vc);
        let a: f64 = vc
            .iter()
            .zip(w.iter())
            .map(|(b, z)| b.0 * z.0 + b.1 * z.1)
            .sum();
        alpha.push(a);
        for i in 0..n {
            w[i].0 -= a * vc[i].0;
            w[i].1 -= a * vc[i].1;
            if j > 0 {
                w[i].0 -= beta[j - 1] * vp[i].0;
                w[i].1 -= beta[j - 1] * vp[i].1;
            }
        }
        let bn: f64 = w.iter().map(|z| z.0 * z.0 + z.1 * z.1).sum::<f64>().sqrt();
        if j + 1 == m || bn < 1e-10 {
            break;
        }
        beta.push(bn);
        vp = vc;
        vc = w.iter().map(|z| (z.0 / bn, z.1 / bn)).collect();
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
    let (evs, vv) = jacobi_eigh(&t, k);
    // 重み |s_j|² = |U[0,j]|², 極 = 重み ≥ 0.1 の最低 Ritz
    for j in 0..k {
        let wj = vv[0 + j * k] * vv[0 + j * k];
        if wj >= 0.1 {
            return (evs[j] - e0, wj);
        }
    }
    (f64::NAN, 0.0)
}

fn main() {
    self_test();
    println!(
        "=== v21.1 アノマリー係数のスペクトル読み出し — M² = (Σq²) g²/π (第二十二期開幕) ===\n"
    );
    println!(
        "事前登録: 連続外挿の質量比 M({{1,1}})/M({{1}}), M({{2}})/M({{1}}), M({{1,2}})/M({{1}}) が"
    );
    println!(
        "          (a) √2, 2, √5 の全てに ±7% = 橋の 1+1D 形成立 / (a′) 2 つ = 部分 / (b) 外れ"
    );
    println!(
        "装置 v2: ボソン = Ê スペクトル関数の低端極 (重み ≥ 0.1)。run1 (保存) の教訓で再設計\n"
    );
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
    let n = 10usize;
    let xs = [1.0f64, 2.0, 3.0, 4.0];
    let sets: Vec<(Vec<i32>, &str, f64)> = vec![
        (vec![1], "{1}", 1.0),
        (vec![1, 1], "{1,1}", 2.0),
        (vec![2], "{2}", 4.0),
        (vec![1, 2], "{1,2}", 5.0),
    ];
    let mut mg_extrap = Vec::new();
    for (qs, label, sq2) in &sets {
        let lam = 3.0 + (*sq2 as f64).sqrt();
        let mut mgs = Vec::new();
        for &x in &xs {
            let core = MultiU1::new(n, x, lam, qs.clone());
            let (e0, gs, res) = ground(&core, 90, 10, 1e-8, 7);
            let (de, wpole) = boson_pole(&core, &gs, e0, 140);
            check(
                &format!("{} x={:.0} の基底残差・ボソン極", label, x),
                res < 1e-7 && wpole >= 0.1 && de.is_finite(),
                format!(
                    "res {:.0e}, 極重み {:.2}, ΔE = {:.4} (dim {}) ({} s)",
                    res,
                    wpole,
                    de,
                    core.dim,
                    t0.elapsed().as_secs()
                ),
            );
            mgs.push(de / (2.0 * x.sqrt()));
        }
        let lx: Vec<f64> = xs.iter().map(|x| 1.0 / x.sqrt()).collect();
        let (a0, _a1) = linfit(&lx, &mgs);
        println!(
            "    {}: M/g(x=1..4) = {:.4}, {:.4}, {:.4}, {:.4} → 外挿 {:.4} (予言 √Σq²/√π = {:.4})",
            label,
            mgs[0],
            mgs[1],
            mgs[2],
            mgs[3],
            a0,
            (sq2).sqrt() / std::f64::consts::PI.sqrt()
        );
        if *label == "{1}" {
            let target = 1.0 / std::f64::consts::PI.sqrt();
            check(
                "{1} の連続外挿 = 1/√π ± 10% (自己完結アンカー)",
                ((a0 - target) / target).abs() < 0.10,
                format!("a₀ = {:.4} vs {:.4}", a0, target),
            );
        }
        mg_extrap.push(a0);
    }

    // ---- 判定 ----
    let base = mg_extrap[0];
    let ratios = [
        mg_extrap[1] / base,
        mg_extrap[2] / base,
        mg_extrap[3] / base,
    ];
    let preds = [2.0f64.sqrt(), 2.0, 5.0f64.sqrt()];
    let hits: Vec<bool> = ratios
        .iter()
        .zip(preds.iter())
        .map(|(r, p)| ((r - p) / p).abs() < 0.07)
        .collect();
    let nhit = hits.iter().filter(|&&b| b).count();
    println!(
        "\n[判定] {}",
        if nhit == 3 {
            "事前登録 (a): アノマリー係数 Σq² がスペクトルから読めた — 台帳の代数と core の動力学が同じ数を出す (橋の 1+1D 形成立)"
        } else if nhit == 2 {
            "事前登録 (a′): 部分成立 — 記録"
        } else {
            "事前登録 (b): 外れ — 記録"
        }
    );
    println!(
        "    質量比: {:.4} (√2={:.4} {}), {:.4} (2 {}), {:.4} (√5={:.4} {})",
        ratios[0],
        preds[0],
        if hits[0] { "✓" } else { "×" },
        ratios[1],
        if hits[1] { "✓" } else { "×" },
        ratios[2],
        preds[2],
        if hits[2] { "✓" } else { "×" }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v21.1".into())),
        (
            "mg_extrap".into(),
            Json::Arr(mg_extrap.iter().map(|&v| Json::Num(v)).collect()),
        ),
        (
            "ratios".into(),
            Json::Arr(ratios.iter().map(|&v| Json::Num(v)).collect()),
        ),
        ("branch_a".into(), Json::Bool(nhit == 3)),
    ]);
    let p = write_artifact("results/v211_anomspec.json", &j.render());
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
