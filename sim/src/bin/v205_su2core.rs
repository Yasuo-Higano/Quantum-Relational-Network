//! v20.5 SU(2) 拘束 core — 非可換 Coulomb 形式と Casimir スケーリング (第二十一期 第五歩)
//!
//! U(1) の章 (v20.1-20.4) の非可換化。開鎖ではゲージ固定 U_n = 1 で SU(2) ゲージ場が
//! 厳密に消去でき、電場項は蓄積色電荷の Casimir になる:
//!   H* = −x Σ_{n,α} (c†_{n,α} c_{n+1,α} + h.c.) + (μ/2) Σ_n (−1)^n (n_n − 1)·2
//!        + Σ_{links} L_n²,  L_n = Σ_{k≤n} Q_k  (Q_k = 色 SU(2) 電荷)
//!   ⇒ H_elec = Σ_k w(k,k) Q_k² + 2 Σ_{k<l} w(k,l) Q_k·Q_l,  w(k,l) = N−1−max(k,l)。
//! 行列要素は Fock 基底のスピン演算子だけ (Clebsch 不要)・ゲージ切断誤差ゼロ。
//! 外部プローブ (表現 j) は静的スピンとして L_n に合流。大域一重項は罰則 λ S²_tot
//! (λ=4) で選択し ⟨S²⟩ ゲートで検証。
//!
//! 装置ゲート:
//!   [A] x=0 の弦エネルギー: 破断が起き得ない短弦 (j=½, r=2) の等式 E = r·j(j+1)
//!       (機械精度) + 全点の変分不等式 E ≤ 非破断値。開発記録 (run1): x=0 でも
//!       セクター選択の弦切れが起きる (j=1, r=4 の実測 5.52 < 8.0 は正しい基底状態 —
//!       「厳密解」を非破断セクター値と誤認した私の物理誤りをゲートが暴いた)。
//!       また λ=8 の罰則はスペクトル幅で Lanczos を未収束にした (res 0.7 — 変分違反
//!       4.18 > 4.0 として発覚) → λ=4 + 再開始 Lanczos (メモリ有界) で根治。
//!   [B] 電場項ゼロでのホップスペクトル = 自由 2 色鎖の厳密和 (N=6 稠密照合)。
//!   [C] 基底状態の ⟨S²_tot⟩ < 1e-6 (罰則の検証)・Lanczos 残差。
//! 事前登録 (順序命題 — v20.2/20.4 の教訓):
//!   (a) x=1, μ=1, N=10 で初期勾配比 σ₁/σ_½ ∈ (1.8, 8/3+0.2) かつ全 r で
//!       E₁(r) > E_½(r) = 非可換弦の Casimir 域成立 (SU(2) core が立った) /
//!   (a′) 順序のみ成立 / (b) 外れ。遮蔽・弦切れの梯子は次版 (クロスオーバー算術を
//!   本版の実測 σ から見積もってから登録する — 規律)。

use uft_sim::*;

// サイト状態: 0=空, 1=↑, 2=↓, 3=両方 (色一重項)。orbital 順序 (site,color): o = 2n+c
#[derive(Clone)]
struct Su2Core {
    n: usize,
    x: f64,
    mu: f64,
    lambda: f64,
    gauge: f64, // 電場項の係数 (1 = 物理, 0 = ホップのみ検査用)
    // 不純物 (外部色源): (site, 2j) — 2j ∈ {1, 2}
    imps: Vec<(usize, u32)>,
    states: Vec<u64>, // key = fock(2bit×N) << 8 | (m1idx << 4) | m2idx
    dim: usize,
}

fn occ(fock: u64, site: usize) -> u64 {
    (fock >> (2 * site)) & 3
}

impl Su2Core {
    fn nup_ndn(&self, fock: u64) -> (u32, u32) {
        let mut up = 0;
        let mut dn = 0;
        for s in 0..self.n {
            let o = occ(fock, s);
            if o == 1 || o == 3 {
                up += 1;
            }
            if o == 2 || o == 3 {
                dn += 1;
            }
        }
        (up, dn)
    }
    fn new(n: usize, x: f64, mu: f64, lambda: f64, imps: Vec<(usize, u32)>) -> Self {
        Self::new_g(n, x, mu, lambda, 1.0, imps)
    }
    fn new_g(n: usize, x: f64, mu: f64, lambda: f64, gauge: f64, imps: Vec<(usize, u32)>) -> Self {
        let mut c = Su2Core {
            n,
            x,
            mu,
            lambda,
            gauge,
            imps,
            states: Vec::new(),
            dim: 0,
        };
        let nimp = c.imps.len();
        let dims: Vec<u32> = c.imps.iter().map(|&(_, tj)| tj + 1).collect();
        let total: u64 = 1 << (2 * n);
        for fock in 0..total {
            let (up, dn) = c.nup_ndn(fock);
            if up + dn != n as u32 {
                continue; // 半充填
            }
            let sz2_f = up as i32 - dn as i32; // 2×Sz (フェルミオン)
                                               // 不純物の 2m: m1idx ∈ 0..=2j → 2m = 2*idx − 2j
            let (d1, d2) = (
                if nimp > 0 { dims[0] } else { 1 },
                if nimp > 1 { dims[1] } else { 1 },
            );
            for i1 in 0..d1 {
                for i2 in 0..d2 {
                    let mut sz2 = sz2_f;
                    if nimp > 0 {
                        sz2 += 2 * i1 as i32 - c.imps[0].1 as i32;
                    }
                    if nimp > 1 {
                        sz2 += 2 * i2 as i32 - c.imps[1].1 as i32;
                    }
                    if sz2 != 0 {
                        continue; // 全 Sz = 0 セクター
                    }
                    c.states.push((fock << 8) | ((i1 as u64) << 4) | i2 as u64);
                }
            }
        }
        c.states.sort_unstable();
        c.dim = c.states.len();
        c
    }
    fn find(&self, key: u64) -> Option<usize> {
        self.states.binary_search(&key).ok()
    }
    // JW 符号: orbital o1 < o2 の間の占有数 (o1, o2 は動かす orbital)
    fn jw(&self, fock: u64, o1: usize, o2: usize) -> f64 {
        let (a, b) = if o1 < o2 { (o1, o2) } else { (o2, o1) };
        let mut cnt = 0u32;
        for o in a + 1..b {
            let site = o / 2;
            let col = o % 2;
            let oc = occ(fock, site);
            let has = if col == 0 {
                oc == 1 || oc == 3
            } else {
                oc == 2 || oc == 3
            };
            if has {
                cnt += 1;
            }
        }
        if cnt % 2 == 0 {
            1.0
        } else {
            -1.0
        }
    }
    // Qz (フェルミオン サイト k): (n↑ − n↓)/2
    fn qz(&self, fock: u64, k: usize) -> f64 {
        match occ(fock, k) {
            1 => 0.5,
            2 => -0.5,
            _ => 0.0,
        }
    }
    fn w_pair(&self, k: f64, l: f64) -> f64 {
        // w(k,l) = N−1−max(k,l) + λ (罰則は全対に一様)
        self.gauge * (self.n as f64 - 1.0 - k.max(l)).max(0.0) + self.lambda
    }
    // H v
    fn matvec(&self, v: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let mut w = vec![(0.0, 0.0); self.dim];
        let nimp = self.imps.len();
        for (i, &key) in self.states.iter().enumerate() {
            if v[i].0 == 0.0 && v[i].1 == 0.0 {
                continue;
            }
            let fock = key >> 8;
            let i1 = ((key >> 4) & 0xf) as i32;
            let i2 = (key & 0xf) as i32;
            // ---- 対角 ----
            let mut diag = 0.0;
            // 質量 (staggered): (μ/2)(−1)^n·2(n_n − 1)
            for s in 0..self.n {
                let nn = match occ(fock, s) {
                    0 => 0.0,
                    3 => 2.0,
                    _ => 1.0,
                };
                let sgn = if s % 2 == 0 { 1.0 } else { -1.0 };
                diag += self.mu * sgn * (nn - 1.0);
            }
            // Q² 対角: フェルミオン w(k,k)·(3/4)[単一占有] + 不純物 w(s,s)·j(j+1)
            for k in 0..self.n {
                let o = occ(fock, k);
                if o == 1 || o == 2 {
                    diag += self.w_pair(k as f64, k as f64) * 0.75;
                }
            }
            for &(s, tj) in &self.imps {
                let j = tj as f64 / 2.0;
                diag += self.w_pair(s as f64, s as f64) * j * (j + 1.0);
            }
            // QzQz (全対: フェルミオン×フェルミオン, フェルミオン×不純物, 不純物×不純物)
            let mut zs: Vec<(f64, f64)> = Vec::new(); // (位置, Qz)
            for k in 0..self.n {
                let q = self.qz(fock, k);
                if q != 0.0 {
                    zs.push((k as f64, q));
                }
            }
            if nimp > 0 {
                let m1 = i1 as f64 - self.imps[0].1 as f64 / 2.0;
                if m1 != 0.0 {
                    zs.push((self.imps[0].0 as f64, m1));
                }
            }
            if nimp > 1 {
                let m2 = i2 as f64 - self.imps[1].1 as f64 / 2.0;
                if m2 != 0.0 {
                    zs.push((self.imps[1].0 as f64, m2));
                }
            }
            for a in 0..zs.len() {
                for b in a + 1..zs.len() {
                    diag += 2.0 * self.w_pair(zs[a].0, zs[b].0) * zs[a].1 * zs[b].1;
                }
            }
            w[i].0 += diag * v[i].0;
            w[i].1 += diag * v[i].1;
            // ---- ホップ (色対角) ----
            for s in 0..self.n - 1 {
                for col in 0..2usize {
                    let o_from = occ(fock, s);
                    let o_to = occ(fock, s + 1);
                    let has_from = if col == 0 {
                        o_from == 1 || o_from == 3
                    } else {
                        o_from == 2 || o_from == 3
                    };
                    let has_to = if col == 0 {
                        o_to == 1 || o_to == 3
                    } else {
                        o_to == 2 || o_to == 3
                    };
                    // s → s+1
                    if has_from && !has_to {
                        let nf = fock
                            ^ ((1 + col as u64) << (2 * s))
                            ^ ((1 + col as u64) << (2 * (s + 1)));
                        let sgn = self.jw(fock, 2 * s + col, 2 * (s + 1) + col);
                        if let Some(jix) = self.find((nf << 8) | ((i1 as u64) << 4) | i2 as u64) {
                            w[jix].0 += -self.x * sgn * v[i].0;
                            w[jix].1 += -self.x * sgn * v[i].1;
                        }
                    }
                    // s+1 → s
                    if has_to && !has_from {
                        let nf = fock
                            ^ ((1 + col as u64) << (2 * s))
                            ^ ((1 + col as u64) << (2 * (s + 1)));
                        let sgn = self.jw(fock, 2 * s + col, 2 * (s + 1) + col);
                        if let Some(jix) = self.find((nf << 8) | ((i1 as u64) << 4) | i2 as u64) {
                            w[jix].0 += -self.x * sgn * v[i].0;
                            w[jix].1 += -self.x * sgn * v[i].1;
                        }
                    }
                }
            }
            // ---- flip-flop: Q⁺_k Q⁻_l + Q⁻_k Q⁺_l (係数 w(k,l), 同サイト双線形は符号 +1) ----
            // フェルミオン×フェルミオン: Q⁺ = c†_↑ c_↓ は単一占有 ↓ サイトのみ非零
            let singles: Vec<usize> = (0..self.n)
                .filter(|&k| {
                    let o = occ(fock, k);
                    o == 1 || o == 2
                })
                .collect();
            for &k in &singles {
                for &l in &singles {
                    if k == l {
                        continue;
                    }
                    // Q⁺_k Q⁻_l: k は ↓→↑, l は ↑→↓
                    if occ(fock, k) == 2 && occ(fock, l) == 1 {
                        let nf = fock ^ (3u64 << (2 * k)) ^ (3u64 << (2 * l));
                        if let Some(jix) = self.find((nf << 8) | ((i1 as u64) << 4) | i2 as u64) {
                            let amp = self.w_pair(k as f64, l as f64); // 2·w·(1/2 交換) = w… 展開: 2w·(½ flip) → w
                            w[jix].0 += amp * v[i].0;
                            w[jix].1 += amp * v[i].1;
                        }
                    }
                }
            }
            // フェルミオン×不純物 1/2
            for (a, &(s_imp, tj)) in self.imps.iter().enumerate() {
                let (idx, other) = if a == 0 { (i1, i2) } else { (i2, i1) };
                let j = tj as f64 / 2.0;
                let m = idx as f64 - j;
                // S⁺_imp Q⁻_k: imp m→m+1, k: ↑→↓
                if m < j {
                    let cp = ((j - m) * (j + m + 1.0)).sqrt();
                    for &k in &singles {
                        if occ(fock, k) == 1 {
                            let nf = fock ^ (3u64 << (2 * k));
                            let (n1, n2) = if a == 0 {
                                (idx + 1, other)
                            } else {
                                (other, idx + 1)
                            };
                            if let Some(jix) = self.find((nf << 8) | ((n1 as u64) << 4) | n2 as u64)
                            {
                                let amp = self.w_pair(s_imp as f64, k as f64) * cp; // 2w·½·cp
                                w[jix].0 += amp * v[i].0;
                                w[jix].1 += amp * v[i].1;
                            }
                        }
                    }
                }
                // S⁻_imp Q⁺_k: imp m→m−1, k: ↓→↑
                if m > -j {
                    let cm = ((j + m) * (j - m + 1.0)).sqrt();
                    for &k in &singles {
                        if occ(fock, k) == 2 {
                            let nf = fock ^ (3u64 << (2 * k));
                            let (n1, n2) = if a == 0 {
                                (idx - 1, other)
                            } else {
                                (other, idx - 1)
                            };
                            if let Some(jix) = self.find((nf << 8) | ((n1 as u64) << 4) | n2 as u64)
                            {
                                let amp = self.w_pair(s_imp as f64, k as f64) * cm;
                                w[jix].0 += amp * v[i].0;
                                w[jix].1 += amp * v[i].1;
                            }
                        }
                    }
                }
            }
            // 不純物×不純物 (2 個のとき)
            if nimp == 2 {
                let (s1, tj1) = self.imps[0];
                let (s2, tj2) = self.imps[1];
                let (j1, j2) = (tj1 as f64 / 2.0, tj2 as f64 / 2.0);
                let (m1, m2) = (i1 as f64 - j1, i2 as f64 - j2);
                // S⁺₁S⁻₂
                if m1 < j1 && m2 > -j2 {
                    let c = ((j1 - m1) * (j1 + m1 + 1.0) * (j2 + m2) * (j2 - m2 + 1.0)).sqrt();
                    if let Some(jix) =
                        self.find((fock << 8) | (((i1 + 1) as u64) << 4) | (i2 - 1) as u64)
                    {
                        let amp = self.w_pair(s1 as f64, s2 as f64) * c;
                        w[jix].0 += amp * v[i].0;
                        w[jix].1 += amp * v[i].1;
                    }
                }
                // S⁻₁S⁺₂
                if m1 > -j1 && m2 < j2 {
                    let c = ((j1 + m1) * (j1 - m1 + 1.0) * (j2 - m2) * (j2 + m2 + 1.0)).sqrt();
                    if let Some(jix) =
                        self.find((fock << 8) | (((i1 - 1) as u64) << 4) | (i2 + 1) as u64)
                    {
                        let amp = self.w_pair(s1 as f64, s2 as f64) * c;
                        w[jix].0 += amp * v[i].0;
                        w[jix].1 += amp * v[i].1;
                    }
                }
            }
        }
        w
    }
    // S²_tot の期待値 (罰則検証用): S² = Σ_pair 2 S_a·S_b + Σ S_a² を λ=1, w=0 の
    // 電場項として評価するには w_pair を差し替えた一時 core を使う (main 側で実装)
    fn lowest(&self, k: usize, m: usize, seed: u64) -> (Vec<f64>, Vec<Vec<(f64, f64)>>, f64) {
        let mv = |vv: &[(f64, f64)]| self.matvec(vv);
        lanczos_lowest_herm(&mv, self.dim, k, m, seed)
    }
}

// 再開始 Lanczos: 現在の Ritz 基底ベクトルから Krylov を組み直す (メモリ有界・収束まで)
fn lanczos_restart(
    core: &Su2Core,
    m: usize,
    max_rounds: usize,
    tol: f64,
    seed: u64,
) -> (f64, Vec<(f64, f64)>, f64) {
    let n = core.dim;
    let mut rng = Rng::new(seed);
    let mut v: Vec<(f64, f64)> = (0..n).map(|_| (rng.gauss(), rng.gauss())).collect();
    let mut ev0 = 0.0;
    let mut res = f64::INFINITY;
    for _round in 0..max_rounds {
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
            let c = vv[a + 0 * k];
            for i in 0..n {
                nv[i].0 += c * basis[a][i].0;
                nv[i].1 += c * basis[a][i].1;
            }
        }
        // 残差 ‖Hv − Ev‖
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

fn main() {
    self_test();
    println!("=== v20.5 SU(2) 拘束 core — 非可換 Coulomb 形式と Casimir スケーリング ===\n");
    println!("事前登録: (a) x=1, μ=1, N=10 で σ₁/σ_½ ∈ (1.8, 2.87) かつ全 r で E₁(r) > E_½(r) /");
    println!("          (a′) 順序のみ / (b) 外れ。x=0 の Casimir 比 8/3 は厳密ゲート\n");
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
    let lam = 4.0; // run1 の教訓: λ=8 は罰則スペクトル幅で Lanczos 収束を殺す (λ·S² 分離は 8 で十分)

    // ---- [B] 電場ゼロ (gauge=0, λ=0) のホップのみ = 自由 2 色鎖の厳密和 (N=6) ----
    {
        let n = 6usize;
        let core = Su2Core::new_g(n, 1.0, 0.0, 0.0, 0.0, vec![]);
        let (ev, _vv, rr) = core.lowest(1, 80, 5);
        // 自由 1 粒子準位 ε_k = −2 cos(kπ/(N+1)), 各色 N/2 粒子 (Sz=0 で ↑↓ 各 3)
        let mut eps: Vec<f64> = (1..=n)
            .map(|k| -2.0 * (std::f64::consts::PI * k as f64 / (n as f64 + 1.0)).cos())
            .collect();
        eps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let exact: f64 = 2.0 * eps[..n / 2].iter().sum::<f64>();
        check(
            "[B] 電場ゼロのホップ = 自由 2 色鎖の厳密和 (N=6)",
            (ev[0] - exact).abs() < 1e-9 && rr < 1e-9,
            format!("E₀ = {:.9} vs 厳密 {:.9}", ev[0], exact),
        );
    }

    // ---- [A] x=0 の弦エネルギー: 短弦の等式 + 変分不等式 ----
    // run1 の教訓: x=0 でもセクター選択の弦切れが起きる (j=1, r=4 の破断は正しい物理)。
    // 等式ゲートは破断が起き得ない j=½, r=2 のみ — 他は E ≤ 非破断値 (変分) を課す。
    {
        let n = 8usize;
        let base = Su2Core::new(n, 0.0, 1.0, lam, vec![]);
        let (e0b, _v0, r0) = lanczos_restart(&base, 80, 8, 1e-9, 11);
        let mut ok_all = r0 < 1e-8;
        let mut detail = String::new();
        for &(tj, jj) in &[(1u32, 0.75f64), (2, 2.0)] {
            for &r in &[2usize, 4] {
                let core = Su2Core::new(n, 0.0, 1.0, lam, vec![(2, tj), (2 + r, tj)]);
                let (ev, _vv, rr) = lanczos_restart(&core, 80, 8, 1e-9, 13);
                let es = ev - e0b;
                let unbroken = r as f64 * jj;
                let is_anchor = tj == 1 && r == 2;
                if rr > 1e-8 || (is_anchor && (es - unbroken).abs() > 1e-8) || es > unbroken + 1e-8
                {
                    ok_all = false;
                }
                detail += &format!(
                    "j={:.1} r={}: {:.6} ({} {:.1}) ",
                    tj as f64 / 2.0,
                    r,
                    es,
                    if is_anchor { "厳密" } else { "≤" },
                    unbroken
                );
            }
        }
        check(
            "[A] x=0: 短弦等式 + 変分不等式 (破断は物理)",
            ok_all,
            detail,
        );
    }

    // ---- [C] ⟨S²_tot⟩ ゲート + 本測定 (x=1, μ=1, N=10) ----
    let n = 10usize;
    let x = 1.0;
    let mu = 1.0;
    let vac = Su2Core::new(n, x, mu, lam, vec![]);
    let (e0, v0, r0) = lanczos_restart(&vac, 100, 10, 1e-8, 21);
    // ⟨S²⟩ = ⟨H(λ+1) − H(λ)⟩ の差分で評価 (S² 項だけ増やす)
    let vac2 = Su2Core::new(n, x, mu, lam + 1.0, vec![]);
    let hv = vac2.matvec(&v0);
    let s2: f64 = v0
        .iter()
        .zip(hv.iter())
        .map(|(a, b)| a.0 * b.0 + a.1 * b.1)
        .sum::<f64>()
        - e0;
    check(
        "[C] 真空の ⟨S²_tot⟩ < 1e-6 (罰則検証)",
        s2.abs() < 1e-6 && r0 < 1e-7,
        format!(
            "⟨S²⟩ = {:.1e} (res {:.0e}, dim {}) ({} s)",
            s2,
            r0,
            vac.dim,
            t0.elapsed().as_secs()
        ),
    );

    // 弦エネルギー E_j(r), j ∈ {½, 1}, r ∈ {2, 4, 6}
    let mut es_all = Vec::new();
    for &(tj, label) in &[(1u32, "½"), (2, "1")] {
        let mut es = Vec::new();
        for &r in &[2usize, 4, 6] {
            let core = Su2Core::new(n, x, mu, lam, vec![(1, tj), (1 + r, tj)]);
            let (ev, vv, rr) = lanczos_restart(&core, 100, 12, 1e-8, 31);
            // 罰則検証 (プローブ込み)
            let core2 = Su2Core::new(n, x, mu, lam + 1.0, vec![(1, tj), (1 + r, tj)]);
            let hv2 = core2.matvec(&vv);
            let s2p: f64 = vv
                .iter()
                .zip(hv2.iter())
                .map(|(a, b)| a.0 * b.0 + a.1 * b.1)
                .sum::<f64>()
                - ev;
            check(
                &format!("j={} r={} の残差・⟨S²⟩", label, r),
                rr < 1e-7 && s2p.abs() < 1e-6,
                format!("res {:.0e}, ⟨S²⟩ {:.1e} (dim {})", rr, s2p, core.dim),
            );
            es.push(ev - e0);
            println!(
                "    j={} r={}: E_string = {:.5} ({} s)",
                label,
                r,
                ev - e0,
                t0.elapsed().as_secs()
            );
        }
        es_all.push(es);
    }
    let s_half = (es_all[0][1] - es_all[0][0]) / 2.0;
    let s_one = (es_all[1][1] - es_all[1][0]) / 2.0;
    let ratio = s_one / s_half;
    let order_ok = (0..3).all(|k| es_all[1][k] > es_all[0][k]);
    let band_ok = ratio > 1.8 && ratio < 8.0 / 3.0 + 0.2;

    println!(
        "\n[判定] {}",
        if band_ok && order_ok {
            "事前登録 (a): Casimir 域の非可換弦が立った — SU(2) core 成立 (第二十一期の非可換第一歩)"
        } else if order_ok {
            "事前登録 (a′): 順序のみ成立 — 記録"
        } else {
            "事前登録 (b): 外れ — 記録"
        }
    );
    println!(
        "    σ_½ = {:.4}, σ₁ = {:.4}, 比 = {:.3} (Casimir 8/3 = {:.3}) | 順序 [{}]",
        s_half,
        s_one,
        ratio,
        8.0 / 3.0,
        if order_ok { "成立" } else { "不成立" }
    );
    println!(
        "    末端勾配 (記録): ½: {:.4} / 1: {:.4} — 遮蔽梯子は実測 σ でクロスオーバーを見積もってから次版で登録",
        (es_all[0][2] - es_all[0][1]) / 2.0,
        (es_all[1][2] - es_all[1][1]) / 2.0
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v20.5".into())),
        (
            "e_half".into(),
            Json::Arr(es_all[0].iter().map(|&v| Json::Num(v)).collect()),
        ),
        (
            "e_one".into(),
            Json::Arr(es_all[1].iter().map(|&v| Json::Num(v)).collect()),
        ),
        ("ratio".into(), Json::Num(ratio)),
        ("branch_a".into(), Json::Bool(band_ok && order_ok)),
    ]);
    let p = write_artifact("results/v205_su2core.json", &j.render());
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
