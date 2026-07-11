//! v20.3 θ 階段の完成形 — 弦切れ・折り返し・½ 弦の持続 (第二十一期 第三歩)
//!
//! v20.2 は θ 階段 (σ(q), q∈{¼,½,¾,1}) を r ≤ 6 で測り、弦切れ (q=1) と折り返し (q=¾)
//! のクロスオーバー (算術: 切れ r* ≈ 2m_pair/σ ≈ 4-5・¾ のフェルミオン遮蔽 r* ≈ 5-6・
//! 巻き付きによる折り返しは r > 2N/3 で環では不可) の手前だった。本版は N=18・
//! r ∈ {2,4,6,8} に伸ばし、判定を閾値でなく**順序命題**で登録する (v20.2 の 3 走の教訓:
//! メソスコピック窓での漸近閾値の推定は脆い — 順序は頑健):
//!   ρ(q) = s_tail/s_init, s_init = [E(4)−E(2)]/2, s_tail = [E(8)−E(6)]/2。
//!   物理: q=1 は対生成で弦が切れ ρ→小 / q=¾ は整数遮蔽で ¼ に折れ始め ρ 中間 /
//!   q=½ は |±½| 縮退で整数遮蔽が効かず弦が持続 ρ≈1 — compact U(1) の署名。
//! 事前登録: (a) ρ(1) < ρ(¾) < ρ(½) (厳密順序) かつ ρ(1) < 0.5·ρ(½) かつ ρ(½) > 0.6
//!   = θ 構造成立 / (a′) ρ(½) > 0.6 かつ [順序 or 切れ] の片方 / (b) それ以外。
//! 装置ゲート: v20.1 回帰 (M/g(x=2, N=18) = 0.7179 ± 1e-3)・Lanczos 残差。
//! core 構造体は v20.1/v20.2 と同一コード ([A] 検証済み)。

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
    println!("=== v20.3 θ 階段の完成形 — 弦切れ・折り返し・½ 弦の持続 ===\n");
    println!("事前登録: (a) ρ(1) < ρ(¾) < ρ(½) かつ ρ(1) < 0.5ρ(½) かつ ρ(½) > 0.6 /");
    println!(
        "          (a′) ρ(½) > 0.6 + 片方 / (b) それ以外。ρ = 末端勾配/初期勾配 (N=18, r≤8)\n"
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
    let (n, x, mu) = (18usize, 2.0f64, 1.0f64);

    // ---- 回帰ゲート: v20.1 の M/g(x=2, N=18) = 0.7179 ----
    {
        let core = U1Core::new(n, n / 2, x, 0.0, 3.0, vec![]);
        let (ev, res) = core.lowest(2, 150, 7);
        let mg = (ev[1] - ev[0]) / (2.0 * x.sqrt());
        check(
            "v20.1 回帰 (M/g(N=18) = 0.7179 ± 1e-3)",
            (mg - 0.7179).abs() < 1e-3 && res < 1e-5,
            format!(
                "M/g = {:.4} (res {:.0e}, dim {}) ({} s)",
                mg,
                res,
                core.dim,
                t0.elapsed().as_secs()
            ),
        );
    }

    // ---- θ 階段 (μ=1) ----
    let vac = U1Core::new(n, n / 2, x, mu, 3.75, vec![]);
    let (e0, r0) = vac.lowest(1, 120, 61);
    check("真空 Lanczos 残差", r0 < 1e-6, format!("res = {:.0e}", r0));
    let mut rhos = Vec::new();
    let mut sig_init = Vec::new();
    for &q in &[0.25f64, 0.5, 0.75, 1.0] {
        let mut es = Vec::new();
        let mut resmax = 0.0f64;
        for &r in &[2usize, 4, 6, 8] {
            let core = U1Core::new(n, n / 2, x, mu, 3.75, vec![(0, q), (r, -q)]);
            let (ev, res) = core.lowest(1, 120, 71);
            resmax = resmax.max(res);
            es.push(ev[0] - e0[0]);
        }
        let s_i = (es[1] - es[0]) / 2.0;
        let s_t = (es[3] - es[2]) / 2.0;
        let rho = s_t / s_i;
        println!(
            "    q={:.2}: E(2,4,6,8) = {:.4}, {:.4}, {:.4}, {:.4} → 勾配 {:.4}→{:.4}, ρ = {:.3} (res≤{:.0e}) ({} s)",
            q, es[0], es[1], es[2], es[3], s_i, s_t, rho, resmax, t0.elapsed().as_secs()
        );
        check(
            &format!("q={:.2} の Lanczos 残差", q),
            resmax < 1e-5,
            format!("max res = {:.0e}", resmax),
        );
        rhos.push(rho);
        sig_init.push(s_i);
    }
    let (r14, r12, r34, r1) = (rhos[0], rhos[1], rhos[2], rhos[3]);
    let order_ok = r1 < r34 && r34 < r12;
    let break_ok = r1 < 0.5 * r12;
    let half_ok = r12 > 0.6;

    // ---- 判定 (記録) ----
    println!(
        "\n[判定] {}",
        if order_ok && break_ok && half_ok {
            "事前登録 (a): θ 構造成立 — 整数は切れ・¾ は折れ始め・½ 弦は持続 (compact U(1) の署名)"
        } else if half_ok && (order_ok || break_ok) {
            "事前登録 (a′): ½ 弦の持続 + 片方 — 部分成立の記録"
        } else {
            "事前登録 (b): 記録"
        }
    );
    println!(
        "    ρ(¼,½,¾,1) = {:.3}, {:.3}, {:.3}, {:.3} | 順序 [{}] 切れ [{}] ½ 持続 [{}]",
        r14,
        r12,
        r34,
        r1,
        if order_ok { "成立" } else { "不成立" },
        if break_ok { "成立" } else { "不成立" },
        if half_ok { "成立" } else { "不成立" }
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v20.3".into())),
        (
            "rhos".into(),
            Json::Arr(rhos.iter().map(|&v| Json::Num(v)).collect()),
        ),
        (
            "sig_init".into(),
            Json::Arr(sig_init.iter().map(|&v| Json::Num(v)).collect()),
        ),
        (
            "branch_a".into(),
            Json::Bool(order_ok && break_ok && half_ok),
        ),
    ]);
    let p = write_artifact("results/v203_staircase.json", &j.render());
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
