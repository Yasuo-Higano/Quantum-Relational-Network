//! v20.1 compact U(1) core — 格子 Schwinger 模型 (第二十一期開幕: PROMPT/3 v21 の第一歩)
//!
//! QrnCoreV2 (v15.3, Z₂) の U(1) 昇格。1+1D compact U(1) + staggered fermion (環) の
//! Gauss 拘束を基底の構造として厳密に解く: 基底 = (占有 bitmask, 巻き付き電場 ε ∈ ℤ)、
//! ボンド電場 E_n = ε + P_n (P_n = 累積電荷)。Z₂ の ε = ±1 が ε ∈ ℤ (切断 |E_n| ≤ Λ)
//! になるだけで、継ぎ目ホップが ε を ±1 する構造は同型。
//! H* = −x Σ (c†c + h.c.) + (μ/2) Σ (−1)^n (2n_c − 1) + Σ_n E_n²,
//! x = 1/(ga)², μ = 2m/(g²a), 物理エネルギー = (g²a/2) H*。M/g = ΔE*/(2√x)。
//!
//! 検証 (装置ゲート):
//!   [A] 拘束の正しさ: L=4, Λ=1 のリンク陽性表現 (3⁴×2⁴ = 1296 次元, U がボンド E を
//!       昇降) を稠密対角化し、(i) H が Gauss セクターを保つ、(ii) 物理 (中性・全 G=0)
//!       セクターのスペクトル = 拘束を解いた core のスペクトル (機械精度)、
//!       (iii) 外部電荷 ±1 セクターも一致。
//!   [B] w=0 (x=0) の対角厳密解と一致。Λ 収束 (3→4 で ΔM/g < 1e-4)。Lanczos 残差。
//! 事前登録 (分岐):
//!   (a) 連続外挿 M/g = a₀ + a₁/√x (x∈{1,2,3,4}, N=18) が 1/√π = 0.5642 ± 0.05 かつ
//!       遮蔽判別が成立 (m=0: ±½ 外部電荷の弦エネルギーが飽和 [遮蔽] /
//!       μ=1: 線形成長が持続 [分数電荷の閉じ込め] — Coleman の判別) = U(1) core 成立 /
//!   (a′) 片方のみ = 部分成立の記録 / (b) 両方外れ = 記録。
//!   カイラル凝縮 ⟨ψ̄ψ⟩/g (厳密値 −e^γ/(2π^{3/2}) = −0.1599) は記録 (外挿が遅いため)。

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
    fn lowest(&self, k: usize, m: usize, seed: u64) -> (Vec<f64>, f64) {
        let mv = |vv: &[(f64, f64)]| self.matvec(vv);
        let (ev, _vecs, res) = lanczos_lowest_herm(&mv, self.dim, k, m, seed);
        (ev, res)
    }
}

// ---- [A] リンク陽性表現 (L=4, Λ=1): fermion 2^4 ⊗ links 3^4 = 1296 ----
struct Explicit {
    n: usize,
    x: f64,
    mu: f64,
    dim: usize, // 16 * 81
}

impl Explicit {
    fn idx(&self, fmask: usize, links: &[i32]) -> usize {
        let mut li = 0usize;
        for &l in links {
            li = li * 3 + (l + 1) as usize;
        }
        fmask * 81 + li
    }
    fn links_of(&self, i: usize) -> (usize, Vec<i32>) {
        let fmask = i / 81;
        let mut li = i % 81;
        let mut links = vec![0i32; self.n];
        for k in (0..self.n).rev() {
            links[k] = (li % 3) as i32 - 1;
            li /= 3;
        }
        (fmask, links)
    }
    // Gauss 残差: G_n = E_n − E_{n−1} − q_n (中性真空基準)
    fn gauss(&self, fmask: usize, links: &[i32]) -> Vec<i32> {
        (0..self.n)
            .map(|site| {
                let occ = ((fmask >> site) & 1) as i32;
                let q = occ - if site % 2 == 1 { 1 } else { 0 };
                let e_here = links[site];
                let e_prev = links[(site + self.n - 1) % self.n];
                e_here - e_prev - q
            })
            .collect()
    }
    fn build_h(&self) -> Vec<f64> {
        let d = self.dim;
        let mut h = vec![0.0f64; d * d];
        for i in 0..d {
            let (fm, links) = self.links_of(i);
            // 対角: Σ E² + 質量
            let elec: f64 = links.iter().map(|&l| (l * l) as f64).sum();
            let mut mass = 0.0;
            for site in 0..self.n {
                let occ = ((fm >> site) & 1) as f64;
                let sgn = if site % 2 == 0 { 1.0 } else { -1.0 };
                mass += sgn * (2.0 * occ - 1.0);
            }
            h[i + i * d] += elec + 0.5 * self.mu * mass;
            // ホップ c†_{s+1} U†_s c_s + c†_s U_s c_{s+1} (E は「電荷が右へ動くと +1」)
            for site in 0..self.n {
                let s1 = (site + 1) % self.n;
                let b0 = (fm >> site) & 1;
                let b1 = (fm >> s1) & 1;
                if b0 == b1 {
                    continue;
                }
                let nfm = fm ^ (1 << site) ^ (1 << s1);
                // JW 符号: 間のサイト数 (環では継ぎ目のみ非自明: (−1)^{N_f−1})
                let sgn = if site == self.n - 1 {
                    let nf = (fm as u32).count_ones();
                    if (nf - 1) % 2 == 0 {
                        1.0
                    } else {
                        -1.0
                    }
                } else {
                    1.0
                };
                // 電荷移動の向きで E_site を ∓1 (Gauss 則 E_k − E_{k−1} = q_k と整合:
                // 電荷が右 (site → s1) へ動くと E_site は −1。∂E/∂t = −j)
                let de = if b0 == 1 { -1 } else { 1 };
                let mut nl = links.clone();
                nl[site] += de;
                if nl[site].abs() > 1 {
                    continue; // 切断
                }
                let j = self.idx(nfm, &nl);
                h[j + i * d] += -self.x * sgn;
            }
        }
        h
    }
}

fn main() {
    self_test();
    println!("=== v20.1 compact U(1) core — 格子 Schwinger 模型 (第二十一期開幕) ===\n");
    println!("事前登録: (a) M/g の連続外挿 (a₀ + a₁/√x, x∈{{1,2,3,4}}, N=18) = 1/√π ± 0.05 かつ");
    println!("          遮蔽判別 (m=0: ±½ 弦の飽和 / μ=1: 線形持続) の両立 = U(1) core 成立 /");
    println!("          (a′) 片方 = 部分成立 / (b) 両方外れ = 記録。凝縮は記録。");
    println!("装置ゲート: [A] リンク陽性表現との機械精度一致 (i)(ii)(iii) / [B] x=0 厳密解・");
    println!("          Λ 収束・Lanczos 残差\n");
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

    // ---- [A] 拘束の正しさ (L=4, Λ=1) ----
    {
        let n = 4usize;
        let (x, mu) = (0.7, 0.3);
        let ex = Explicit {
            n,
            x,
            mu,
            dim: 16 * 81,
        };
        let h = ex.build_h();
        // (i) H が Gauss セクターを保つ: H の非零要素 (i,j) で gauss(i) ≠ gauss(j) が無いこと
        let d = ex.dim;
        let mut mixed = 0usize;
        for i in 0..d {
            let gi = {
                let (fm, l) = ex.links_of(i);
                ex.gauss(fm, &l)
            };
            for j in 0..d {
                if h[j + i * d].abs() > 1e-12 {
                    let gj = {
                        let (fm, l) = ex.links_of(j);
                        ex.gauss(fm, &l)
                    };
                    if gi != gj {
                        mixed += 1;
                    }
                }
            }
        }
        check(
            "[A-i] H が Gauss セクターを保つ (リンク陽性表現)",
            mixed == 0,
            format!(
                "セクター混合要素 = {} ({} s)",
                mixed,
                t0.elapsed().as_secs()
            ),
        );
        // (ii) 物理セクター (全 G=0) のスペクトル = core のスペクトル
        let phys: Vec<usize> = (0..d)
            .filter(|&i| {
                let (fm, l) = ex.links_of(i);
                ex.gauss(fm, &l).iter().all(|&g| g == 0)
            })
            .collect();
        let m = phys.len();
        let mut hp = vec![0.0f64; m * m];
        for a in 0..m {
            for b in 0..m {
                hp[a + b * m] = h[phys[a] + phys[b] * d];
            }
        }
        let (evp, _) = jacobi_eigh(&hp, m);
        let core = U1Core::new(n, 2, x, mu, 1.0, vec![]);
        let mut hc = vec![0.0f64; core.dim * core.dim];
        for i in 0..core.dim {
            let mut unit = vec![(0.0, 0.0); core.dim];
            unit[i] = (1.0, 0.0);
            let col = core.matvec(&unit);
            for j in 0..core.dim {
                hc[j + i * core.dim] = col[j].0;
            }
        }
        let (evc, _) = jacobi_eigh(&hc, core.dim);
        let mdev = if m == core.dim {
            evp.iter()
                .zip(evc.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max)
        } else {
            9e9
        };
        check(
            "[A-ii] 物理セクター = core (次元・スペクトル機械精度)",
            m == core.dim && mdev < 1e-10,
            format!("dim {} vs {}, 最大偏差 {:.1e}", m, core.dim, mdev),
        );
        // (iii) 外部電荷 ±1 (サイト 0 に +1, サイト 2 に −1) セクターの一致
        let want: Vec<i32> = vec![1, 0, -1, 0];
        let sect: Vec<usize> = (0..d)
            .filter(|&i| {
                let (fm, l) = ex.links_of(i);
                ex.gauss(fm, &l) == want
            })
            .collect();
        let ms = sect.len();
        let mut hs = vec![0.0f64; ms * ms];
        for a in 0..ms {
            for b in 0..ms {
                hs[a + b * ms] = h[sect[a] + sect[b] * d];
            }
        }
        let (evs, _) = jacobi_eigh(&hs, ms);
        let core2 = U1Core::new(n, 2, x, mu, 1.0, vec![(0, 1.0), (2, -1.0)]);
        let mut hc2 = vec![0.0f64; core2.dim * core2.dim];
        for i in 0..core2.dim {
            let mut unit = vec![(0.0, 0.0); core2.dim];
            unit[i] = (1.0, 0.0);
            let col = core2.matvec(&unit);
            for j in 0..core2.dim {
                hc2[j + i * core2.dim] = col[j].0;
            }
        }
        let (evc2, _) = jacobi_eigh(&hc2, core2.dim);
        let mdev2 = if ms == core2.dim && ms > 0 {
            evs.iter()
                .zip(evc2.iter())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0f64, f64::max)
        } else {
            9e9
        };
        check(
            "[A-iii] 外部電荷 ±1 セクター = 外部電荷つき core",
            ms == core2.dim && ms > 0 && mdev2 < 1e-10,
            format!("dim {} vs {}, 最大偏差 {:.1e}", ms, core2.dim, mdev2),
        );
        // 電荷の超選択則: 中性でない裸 core は空でない (電荷 +1 は probe 無しでは構成不能)
        let odd = U1Core::new(n, 3, x, mu, 1.0, vec![]);
        check(
            "[A-iv] 電荷超選択則 (Q≠0 は probe 無しで構成不能)",
            odd.dim == 0,
            format!("dim = {}", odd.dim),
        );
    }

    // ---- [B] x=0 厳密解 + Λ 収束 ----
    {
        let core = U1Core::new(8, 4, 0.0, 0.6, 2.0, vec![]);
        // x=0: H 対角 — 最低値は基底列挙から直接
        let mut emin = f64::INFINITY;
        for &key in &core.states {
            let (mask, eps) = core.unpack(key);
            emin = emin.min(core.diag(mask, eps));
        }
        let (ev, res) = core.lowest(1, 60, 11);
        check(
            "[B-i] x=0 の対角厳密解と Lanczos の一致",
            (ev[0] - emin).abs() < 1e-9 && res < 1e-8,
            format!("E₀ = {:.9} vs {:.9} (res {:.1e})", ev[0], emin, res),
        );
        let c3 = U1Core::new(16, 8, 4.0, 0.0, 3.0, vec![]);
        let c4 = U1Core::new(16, 8, 4.0, 0.0, 4.0, vec![]);
        let (e3, r3) = c3.lowest(2, 140, 7);
        let (e4, r4) = c4.lowest(2, 140, 7);
        let dm = ((e3[1] - e3[0]) - (e4[1] - e4[0])).abs() / (2.0 * 4.0f64.sqrt());
        check(
            "[B-ii] Λ 収束 (x=4, N=16: Λ=3→4 で ΔM/g < 1e-4)",
            dm < 1e-4 && r3 < 1e-7 && r4 < 1e-7,
            format!(
                "|Δ(M/g)| = {:.2e} (res {:.0e}/{:.0e}, dim {}/{})",
                dm, r3, r4, c3.dim, c4.dim
            ),
        );
    }

    // ---- [C-1] ベクトル質量 M/g の連続外挿 (m=0) ----
    println!();
    let xs = [1.0f64, 2.0, 3.0, 4.0];
    let mut mg_n18 = Vec::new();
    for &x in &xs {
        let mut row = format!("    x={:.0}:", x);
        for &n in &[12usize, 16, 18] {
            let core = U1Core::new(n, n / 2, x, 0.0, 3.0, vec![]);
            let (ev, res) = core.lowest(2, 150, 21);
            let mg = (ev[1] - ev[0]) / (2.0 * x.sqrt());
            row += &format!("  N={}: M/g={:.4} (res {:.0e})", n, mg, res);
            if n == 18 {
                mg_n18.push(mg);
            }
        }
        println!("{} ({} s)", row, t0.elapsed().as_secs());
    }
    // a→0: M/g = a₀ + a₁/√x (最小二乗)
    let lx: Vec<f64> = xs.iter().map(|x| 1.0 / x.sqrt()).collect();
    let (a0, a1) = {
        let (ic, sl) = linfit(&lx, &mg_n18);
        (ic, sl)
    };
    let target = 1.0 / std::f64::consts::PI.sqrt();
    println!(
        "    連続外挿: M/g = {:.4} + {:.4}/√x → a₀ = {:.4} (厳密 1/√π = {:.4}, 差 {:+.4})",
        a0,
        a1,
        a0,
        target,
        a0 - target
    );
    let mass_ok = (a0 - target).abs() < 0.05;

    // ---- [C-2] 遮蔽判別 (N=16, x=2, ±½ 外部電荷) ----
    println!();
    let mut slopes = Vec::new(); // (μ, 初期勾配, 末端勾配)
    for &mu in &[0.0f64, 1.0] {
        let base = U1Core::new(16, 8, 2.0, mu, 3.5, vec![]);
        let (e0, _) = base.lowest(1, 120, 31);
        let mut es = Vec::new();
        for &r in &[2usize, 4, 6] {
            // r=8 は N=16 の対蹠点 (両回りの弦が縮退) なので避ける
            let core = U1Core::new(16, 8, 2.0, mu, 3.5, vec![(0, 0.5), (r, -0.5)]);
            let (ev, res) = core.lowest(1, 120, 31);
            es.push(ev[0] - e0[0]);
            println!(
                "    μ={:.0} r={}: E_string = {:.5} (res {:.0e}) ({} s)",
                mu,
                r,
                ev[0] - e0[0],
                res,
                t0.elapsed().as_secs()
            );
        }
        let s_init = (es[1] - es[0]) / 2.0;
        let s_tail = (es[2] - es[1]) / 2.0;
        println!(
            "    μ={:.0}: 勾配 初期 {:.4} → 末端 {:.4}",
            mu, s_init, s_tail
        );
        slopes.push((mu, s_init, s_tail));
    }
    let screen_ok = slopes[0].2.abs() < 0.10 * slopes[0].1.abs() // m=0: 飽和
        && slopes[1].2 > 0.5 * slopes[1].1; // μ=1: 線形持続
                                            // ---- [C-3] カイラル凝縮 (記録) ----
    {
        let n = 18usize;
        let x = 4.0;
        let core = U1Core::new(n, n / 2, x, 0.0, 3.0, vec![]);
        let mv = |vv: &[(f64, f64)]| core.matvec(vv);
        let (ev, vecs, _res) = lanczos_lowest_herm(&mv, core.dim, 1, 150, 41);
        let gs = &vecs[0];
        let mut cond = 0.0;
        for (i, &key) in core.states.iter().enumerate() {
            let (mask, _eps) = core.unpack(key);
            let w = gs[i].0 * gs[i].0 + gs[i].1 * gs[i].1;
            let mut s = 0.0;
            for site in 0..n {
                let occ = ((mask >> site) & 1) as f64;
                let sgn = if site % 2 == 0 { 1.0 } else { -1.0 };
                s += sgn * occ;
            }
            cond += w * s;
        }
        let psibar = cond * x.sqrt() / (n as f64);
        println!(
            "\n    [記録] カイラル凝縮 (N=18, x=4): ⟨ψ̄ψ⟩/g = {:.4} (厳密連続値 −0.1599, E₀* = {:.4})",
            psibar, ev[0]
        );
    }

    // ---- 判定 ----
    println!(
        "\n[判定] {}",
        if mass_ok && screen_ok {
            "事前登録 (a): M/g 外挿と遮蔽判別が両立 — compact U(1) core が Schwinger 物理を再現 (第二十一期の第一歩成立)"
        } else if mass_ok || screen_ok {
            "事前登録 (a′): 部分成立 — record"
        } else {
            "事前登録 (b): 両方外れ — 記録"
        }
    );
    println!(
        "    質量: a₀ = {:.4} vs 1/√π = {:.4} [{}] / 遮蔽: m=0 末端/初期 = {:.3}, μ=1 = {:.3} [{}]",
        a0,
        target,
        if mass_ok { "帯内" } else { "帯外" },
        slopes[0].2 / slopes[0].1,
        slopes[1].2 / slopes[1].1,
        if screen_ok { "成立" } else { "不成立" }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v20.1".into())),
        ("mg_extrap".into(), Json::Num(a0)),
        ("mg_slope".into(), Json::Num(a1)),
        (
            "mg_n18".into(),
            Json::Arr(mg_n18.iter().map(|&v| Json::Num(v)).collect()),
        ),
        (
            "screening".into(),
            Json::Arr(
                slopes
                    .iter()
                    .map(|&(mu, si, st)| {
                        Json::Obj(vec![
                            ("mu".into(), Json::Num(mu)),
                            ("s_init".into(), Json::Num(si)),
                            ("s_tail".into(), Json::Num(st)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("branch_a".into(), Json::Bool(mass_ok && screen_ok)),
    ]);
    let p = write_artifact("results/v201_u1core.json", &j.render());
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
