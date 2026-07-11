//! v20.2 遮蔽の雲と θ 階段 — U(1) 判別子の完成形 (第二十一期 第二歩)
//!
//! v20.1 は Schwinger 質量を 0.7% で的中させたが、遮蔽判別は弦の漸近形 (r ≳ 5ℓ ≈ 10)
//! が ED の環の外で (a′) に留まった。本版は判別子を局所観測量に書き換える:
//!   [雲] ±½ プローブ対 (r=6, N=16) の周りの誘導フェルミオン電荷
//!        I(w) = Σ_{dist≤w} δ⟨q_n⟩ — m=0 なら −½ へ補償 (遮蔽長 ~2 サイトで測定可能)、
//!        μ=1 なら補償は部分的 (だから弦が残る)。
//!        開発記録 (run1, results/v202_screen_run1.txt): 初走は r=8 (対蹠点) に置き、
//!        内外二本の弦 E=±½ の厳密縮退で基底状態が不定になり雲の対称性が壊れた
//!        (和則ゲートが正しく捕捉)。r=6 は縮退が 0.5 開き一意。
//!   [θ 階段] μ=1 の弦張力 σ(q), q ∈ {¼,½,¾,1}: θ 真空の周期性 ε(θ=2πq) から
//!        σ(1) ≈ 0 (整数電荷は対生成で遮蔽 = 弦切れ) かつ σ(¾) ≈ σ(¼) (折り返し) が
//!        compact U(1) の署名 (Coleman)。
//! 事前登録 (v3): 判別子 = 中点ボンドの電場 ⟨E_mid⟩ (弦が運ぶ磁束そのもの —
//! E = ε + 累積電荷なので staggered CDW 振動が累積で平滑化される)。
//!   (a) E_mid(μ=0)/E_mid(μ=1) < 0.5 かつ E_mid(μ=1) > 0.3 = 遮蔽判別完成 /
//!   (b) 外れ = 記録。θ 階段と誘導電荷雲 I(w) は記録 (階段の判定は v20.3 の新登録)。
//! 開発記録 3 走 (results/v202_screen_run{1,2,3}.txt — 各走の器械ゲートが実際の欠陥を捕捉):
//! run1: r=8 は対蹠点 — 内外二本の弦 E=±½ が厳密縮退し基底状態が不定 (和則ゲートが捕捉)。
//! run2/3: 鏡映-C 反対称はどの μ でも厳密でない (staggered の C は 1 サイト並進を伴い
//! 偶サイト固定プローブと非整合 — μ=0 で 3.1e-1, μ=1 で 8.2e-2 の破れは物理) → 記録に降格。
//! 窓積分 I(w) は CDW のノードに乗り閾値が脆い (I(1)=−0.26, I(2)=−0.16, I(3)=−0.33)。
//! E_mid 比の実測 0.521 は登録 0.5 を 4% 外れ — 閾値は動かさず (b) を記録する。
//! 装置ゲート: v20.1 実測への回帰 (M/g(x=2,N=16) = 0.7179 ± 1e-3 — [A] 検証済み器械と
//! 同一であることの錨)・Lanczos 残差。鏡映-C 破れは記録。
//! core 構造体は v20.1 と同一コード (SU(2) 一般化時に lib へ昇格予定)。

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

fn main() {
    self_test();
    println!("=== v20.2 遮蔽の雲と θ 階段 — U(1) 判別子の完成形 ===\n");
    println!(
        "事前登録 (v3): (a) E_mid(μ=0)/E_mid(μ=1) < 0.5 かつ E_mid(μ=1) > 0.3 = 遮蔽判別完成 /"
    );
    println!("          (b) 外れ = 記録。θ 階段・雲 I(w) は記録 (階段の判定は v20.3)\n");
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
    let (n, x) = (16usize, 2.0f64);

    // ---- 回帰ゲート: v20.1 の M/g(x=2, N=16) = 0.7179 ----
    {
        let core = U1Core::new(n, n / 2, x, 0.0, 3.0, vec![]);
        let (ev, res) = core.lowest(2, 140, 7);
        let mg = (ev[1] - ev[0]) / (2.0 * x.sqrt());
        check(
            "v20.1 回帰 (M/g = 0.7179 ± 1e-3)",
            (mg - 0.7179).abs() < 1e-3 && res < 1e-6,
            format!(
                "M/g = {:.4} (res {:.0e}) ({} s)",
                mg,
                res,
                t0.elapsed().as_secs()
            ),
        );
    }

    // ---- 密度・電場読み出しつき基底状態 ----
    let density = |core: &U1Core, seed: u64| -> (f64, Vec<f64>, Vec<f64>, f64) {
        let mv = |vv: &[(f64, f64)]| core.matvec(vv);
        let (ev, vecs, res) = lanczos_lowest_herm(&mv, core.dim, 1, 120, seed);
        let gs = &vecs[0];
        let mut dens = vec![0.0f64; core.n];
        let mut ef = vec![0.0f64; core.n];
        for (i, &key) in core.states.iter().enumerate() {
            let (mask, eps) = core.unpack(key);
            let w = gs[i].0 * gs[i].0 + gs[i].1 * gs[i].1;
            for site in 0..core.n {
                dens[site] += w * ((mask >> site) & 1) as f64;
            }
            let e = core.e_profile(mask, eps);
            for b in 0..core.n {
                ef[b] += w * e[b];
            }
        }
        (ev[0], dens, ef, res)
    };

    // ---- [雲] ±½ @ (0, 8), μ ∈ {0, 1} ----
    println!();
    let mut cloud = Vec::new(); // (μ, E_mid, I(2))
    for &mu in &[0.0f64, 1.0] {
        let vac = U1Core::new(n, n / 2, x, mu, 3.75, vec![]);
        let (_e0, d0, ef0, r0) = density(&vac, 51);
        let core = U1Core::new(n, n / 2, x, mu, 3.75, vec![(0, 0.5), (6, -0.5)]);
        let (_ep, dp, efp, rp) = density(&core, 52);
        let dq: Vec<f64> = (0..n).map(|s| dp[s] - d0[s]).collect();
        let de: Vec<f64> = (0..n).map(|b| efp[b] - ef0[b]).collect();
        let iw = |c: usize, w: usize| -> f64 {
            (0..n)
                .filter(|&s| {
                    let d = (s as i32 - c as i32).unsigned_abs() as usize;
                    d.min(n - d) <= w
                })
                .map(|s| dq[s])
                .sum()
        };
        // 鏡映-C はどの μ でも厳密対称でない (staggered の C は 1 サイト並進を伴い
        // 偶サイト固定プローブと非整合) — run2/3 の教訓。記録に降格
        let asym = (0..n)
            .map(|s| (dq[s] + dq[(6 + n - s) % n]).abs())
            .fold(0.0f64, f64::max);
        println!(
            "    [記録] μ={:.0} の鏡映-C 破れ = {:.1e} (staggered C の並進構造)",
            mu, asym
        );
        check(
            &format!("Lanczos 残差 (μ={:.0})", mu),
            r0 < 1e-6 && rp < 1e-6,
            format!("vac {:.0e} / probe {:.0e}", r0, rp),
        );
        let emid = de[2..4].iter().sum::<f64>() / 2.0; // 中点近傍ボンド 2,3 の平均
        println!(
            "    μ={:.0}: E_mid = {:+.4}, I(w=1,2,3) = {:+.4}, {:+.4}, {:+.4} ({} s)",
            mu,
            emid,
            iw(0, 1),
            iw(0, 2),
            iw(0, 3),
            t0.elapsed().as_secs()
        );
        cloud.push((mu, emid, iw(0, 2)));
    }
    let cloud_ok = (cloud[0].1 / cloud[1].1) < 0.5 && cloud[1].1 > 0.3;

    // ---- [θ 階段] μ=1: σ(q), q ∈ {¼,½,¾,1} ----
    println!();
    let mu = 1.0;
    let vac = U1Core::new(n, n / 2, x, mu, 3.75, vec![]);
    let (e0, r0) = vac.lowest(1, 120, 61);
    check(
        "θ 階段の真空 Lanczos 残差",
        r0 < 1e-6,
        format!("res = {:.0e}", r0),
    );
    let mut sigmas = Vec::new();
    for &q in &[0.25f64, 0.5, 0.75, 1.0] {
        let mut es = Vec::new();
        for &r in &[2usize, 4, 6] {
            let core = U1Core::new(n, n / 2, x, mu, 3.75, vec![(0, q), (r, -q)]);
            let (ev, _res) = core.lowest(1, 120, 71);
            es.push(ev[0] - e0[0]);
        }
        let sig = (es[2] - es[0]) / 4.0;
        println!(
            "    q={:.2}: E(2,4,6) = {:.4}, {:.4}, {:.4} → σ = {:.4} ({} s)",
            q,
            es[0],
            es[1],
            es[2],
            sig,
            t0.elapsed().as_secs()
        );
        sigmas.push(sig);
    }
    // θ 階段は本版では記録 (判定は v20.3 — r ≤ 6 は弦切れ/折り返しのクロスオーバー手前)

    // ---- 判定 ----
    println!(
        "\n[判定] {}",
        if cloud_ok {
            "事前登録 (a): E_mid コントラストが立った — 遮蔽判別完成 (m=0 は磁束が溶け、μ=1 は弦が磁束を運ぶ)"
        } else {
            "事前登録 (b): E_mid コントラスト外れ — 記録"
        }
    );
    println!(
        "    E_mid: {:+.4} (μ=0) vs {:+.4} (μ=1), 比 {:.3} [{}]",
        cloud[0].1,
        cloud[1].1,
        cloud[0].1 / cloud[1].1,
        if cloud_ok { "成立" } else { "不成立" }
    );
    println!(
        "    [記録] θ 階段 σ(¼,½,¾,1) = {:.4}, {:.4}, {:.4}, {:.4} (判定は v20.3)",
        sigmas[0], sigmas[1], sigmas[2], sigmas[3]
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v20.2".into())),
        ("emid_mu0".into(), Json::Num(cloud[0].1)),
        ("emid_mu1".into(), Json::Num(cloud[1].1)),
        ("iw2_mu0".into(), Json::Num(cloud[0].2)),
        ("iw2_mu1".into(), Json::Num(cloud[1].2)),
        (
            "sigmas".into(),
            Json::Arr(sigmas.iter().map(|&s| Json::Num(s)).collect()),
        ),
        ("branch_a".into(), Json::Bool(cloud_ok)),
    ]);
    let p = write_artifact("results/v202_screen.json", &j.render());
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
