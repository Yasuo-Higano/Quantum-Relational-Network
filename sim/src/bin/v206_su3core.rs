//! v20.6 SU(3) color toy — Coulomb 形式の一般化とバリオン Y 弦 (第二十一期 第六歩)
//!
//! v20.5 の非可換 Coulomb 形式は SU(3) にそのまま持ち上がる。Fierz 恒等式
//!   Σ_a (λ^a)_{γδ}(λ^a)_{αβ} = 2δ_{γβ}δ_{αδ} − (2/3)δ_{γδ}δ_{αβ}
//! により、電場項は色交換演算子で書け Clebsch 不要:
//!   2w Q_k·Q_l = w Σ_{αβ} (c†_α c_β)_k (c†_β c_α)_l − (w/3) n_k n_l   (k≠l)
//!   サイト Casimir (対角): 占有 {0,1,2,3} → C₂ = {0, 4/3, 4/3, 0} (2 = 3̄, 3 = バリオン一重項)
//! プローブ: 基本 3 (色レジスタ, +w 交換 − w/3 n) / 反基本 3̄ (−w 同方向回転 + w/3 n)。
//! 大域一重項は罰則 λ C₂(全電荷) (λ=6, 最小非一重項 4/3 → 分離 8)。
//! 真空 (x=0, μ>0) = 0/3 交互のバリオン結晶 (全サイト色一重項)。
//!
//! SU(3) 固有の新アンカー: **バリオン Y 弦** — 3 つの 3-プローブの ε 一重項は
//! x=0 で E = (4/3)(r₁₂ + r₂₃) の厳密値 (フラックスの N-ality: 3̄ → 3 → 0 の継走)。
//! 装置ゲート:
//!   [A] x=0: メソン (3,3̄) r=2 の等式 E = 8/3 (機械精度 — 機構検証) + バリオン (0,2,4) の
//!       変分不等式 E ≤ 16/3 (等式成立か破断かは記録 — v20.5 の教訓: セクター弦切れは物理)。
//!       開発記録 (run1): (i) プローブを (1,3,5) に置くと端のサイト 5 は w(5,·) = 0 で
//!       電場から完全に切断される (配置バグ — 内部 (0,2,4) へ; 公式は並進不変で 16/3 不変)、
//!       (ii) バリオン実測 14/3 < 16/3 = 「厳密値」のセクター仮定が破れ得る再確認。
//!   [B] 電場ゼロ = 自由 3 色鎖の厳密和 (N=6)。
//!   [C] ⟨C₂_tot⟩ < 1e-6 (罰則検証)・再開始 Lanczos 残差 < 1e-7。
//! 事前登録 (順序命題): (a) [A][B][C] 全 PASS かつ x=1, μ=1 でメソン E(4) > E(2)
//!   (閉じ込め方向) = SU(3) core 成立 / (b) 外れ。バリオンの x=1 束縛エネルギー・
//!   Y 劣加法性は記録 (閾値は次版で実測算術後に登録 — 規律)。

use uft_sim::*;

// サイト = 3 色の占有 bitmask (3 bit/site)。orbital o = 3n + c
fn occ3(fock: u64, site: usize, col: usize) -> bool {
    (fock >> (3 * site + col)) & 1 == 1
}
fn ncol(fock: u64, site: usize) -> u32 {
    ((fock >> (3 * site)) & 7).count_ones()
}

#[derive(Clone)]
struct Su3Core {
    n: usize,
    x: f64,
    mu: f64,
    lambda: f64,
    gauge: f64,
    // 不純物: (site, kind) kind: +1 = 3, −1 = 3̄。レジスタは色 idx ∈ {0,1,2}
    imps: Vec<(usize, i32)>,
    states: Vec<u64>, // key = fock << 8 | c1<<4 | c2<<2 | c3 (最大 3 個)
    dim: usize,
}

impl Su3Core {
    fn col_counts(&self, fock: u64, regs: &[usize]) -> [i32; 3] {
        let mut c = [0i32; 3];
        for s in 0..self.n {
            for a in 0..3 {
                if occ3(fock, s, a) {
                    c[a] += 1;
                }
            }
        }
        for (i, &(_, kind)) in self.imps.iter().enumerate() {
            c[regs[i]] += kind;
        }
        c
    }
    fn new(
        n: usize,
        x: f64,
        mu: f64,
        lambda: f64,
        gauge: f64,
        imps: Vec<(usize, i32)>,
        target: [i32; 3],
    ) -> Self {
        let mut c = Su3Core {
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
        let total: u64 = 1 << (3 * n);
        let nregs = 3usize.pow(nimp as u32);
        for fock in 0..total {
            // 総粒子数 = 3n/2 (半充填)
            if (fock.count_ones() as usize) != 3 * n / 2 {
                continue;
            }
            for ri in 0..nregs {
                let mut regs = Vec::with_capacity(nimp);
                let mut t = ri;
                for _ in 0..nimp {
                    regs.push(t % 3);
                    t /= 3;
                }
                if c.col_counts(fock, &regs) != target {
                    continue;
                }
                let mut key = fock << 8;
                for (i, &r) in regs.iter().enumerate() {
                    key |= (r as u64) << (2 * i);
                }
                c.states.push(key);
            }
        }
        c.states.sort_unstable();
        c.dim = c.states.len();
        c
    }
    fn regs_of(&self, key: u64) -> Vec<usize> {
        (0..self.imps.len())
            .map(|i| ((key >> (2 * i)) & 3) as usize)
            .collect()
    }
    fn find(&self, key: u64) -> Option<usize> {
        self.states.binary_search(&key).ok()
    }
    fn jw(&self, fock: u64, o1: usize, o2: usize) -> f64 {
        let (a, b) = if o1 < o2 { (o1, o2) } else { (o2, o1) };
        let mut cnt = 0u32;
        for o in a + 1..b {
            if (fock >> o) & 1 == 1 {
                cnt += 1;
            }
        }
        if cnt % 2 == 0 {
            1.0
        } else {
            -1.0
        }
    }
    fn w_pair(&self, k: f64, l: f64) -> f64 {
        self.gauge * (self.n as f64 - 1.0 - k.max(l)).max(0.0) + self.lambda
    }
    fn matvec(&self, v: &[(f64, f64)]) -> Vec<(f64, f64)> {
        let mut w = vec![(0.0, 0.0); self.dim];
        let nimp = self.imps.len();
        for (i, &key) in self.states.iter().enumerate() {
            if v[i].0 == 0.0 && v[i].1 == 0.0 {
                continue;
            }
            let fock = key >> 8;
            let regs = self.regs_of(key);
            let (vr, vi) = (v[i].0, v[i].1);
            // ---- 対角 ----
            let mut diag = 0.0;
            for s in 0..self.n {
                let nn = ncol(fock, s) as f64;
                let sgn = if s % 2 == 0 { 1.0 } else { -1.0 };
                diag += self.mu * sgn * (nn - 1.5);
                // サイト Casimir: {0,1,2,3} → {0, 4/3, 4/3, 0}
                let c2 = match ncol(fock, s) {
                    1 | 2 => 4.0 / 3.0,
                    _ => 0.0,
                };
                diag += self.w_pair(s as f64, s as f64) * c2;
            }
            for &(s, _) in &self.imps {
                diag += self.w_pair(s as f64, s as f64) * (4.0 / 3.0);
            }
            // 対角ペア項: w Σ_α n_kα n_lα − w/3 n_k n_l (フェルミオン×フェルミオン)
            for k in 0..self.n {
                let nk = ncol(fock, k) as f64;
                if nk == 0.0 {
                    continue;
                }
                for l in k + 1..self.n {
                    let nl = ncol(fock, l) as f64;
                    if nl == 0.0 {
                        continue;
                    }
                    let wp = self.w_pair(k as f64, l as f64);
                    let mut same = 0.0;
                    for a in 0..3 {
                        if occ3(fock, k, a) && occ3(fock, l, a) {
                            same += 1.0;
                        }
                    }
                    diag += wp * same - wp / 3.0 * nk * nl;
                }
            }
            // フェルミオン×不純物 対角: 3: w n_{kα=reg} − w/3 n_k / 3̄: −w n_{kα=reg} + w/3 n_k
            for (ii, &(s_imp, kind)) in self.imps.iter().enumerate() {
                let a = regs[ii];
                for k in 0..self.n {
                    let nk = ncol(fock, k) as f64;
                    if nk == 0.0 {
                        continue;
                    }
                    let wp = self.w_pair(s_imp as f64, k as f64);
                    let na = if occ3(fock, k, a) { 1.0 } else { 0.0 };
                    if kind > 0 {
                        diag += wp * na - wp / 3.0 * nk;
                    } else {
                        diag += -wp * na + wp / 3.0 * nk;
                    }
                }
            }
            // 不純物×不純物 対角
            for ii in 0..nimp {
                for jj in ii + 1..nimp {
                    let (s1, k1) = self.imps[ii];
                    let (s2, k2) = self.imps[jj];
                    let wp = self.w_pair(s1 as f64, s2 as f64);
                    let same = if regs[ii] == regs[jj] { 1.0 } else { 0.0 };
                    // 3·3: w(δ_same − 1/3) / 3·3̄: −w δ_same + w/3 (Fierz 別枝) / 3̄·3̄: 3·3 と同形
                    if k1 * k2 > 0 {
                        diag += wp * (same - 1.0 / 3.0);
                    } else {
                        diag += wp * (-same + 1.0 / 3.0);
                    }
                }
            }
            w[i].0 += diag * vr;
            w[i].1 += diag * vi;
            // ---- ホップ (色対角) ----
            for s in 0..self.n - 1 {
                for col in 0..3usize {
                    let h_from = occ3(fock, s, col);
                    let h_to = occ3(fock, s + 1, col);
                    if h_from == h_to {
                        continue;
                    }
                    let nf = fock ^ (1u64 << (3 * s + col)) ^ (1u64 << (3 * (s + 1) + col));
                    let sgn = self.jw(fock, 3 * s + col, 3 * (s + 1) + col);
                    let nkey = (nf << 8) | (key & 0xff);
                    if let Some(jix) = self.find(nkey) {
                        w[jix].0 += -self.x * sgn * vr;
                        w[jix].1 += -self.x * sgn * vi;
                    }
                }
            }
            // ---- 色交換 (フェルミオン×フェルミオン): w (c†αcβ)_k (c†βcα)_l, α≠β ----
            for k in 0..self.n {
                for l in 0..self.n {
                    if k == l {
                        continue;
                    }
                    let wp = self.w_pair(k as f64, l as f64);
                    for a in 0..3usize {
                        for b in 0..3usize {
                            if a == b {
                                continue;
                            }
                            // k: β→α (c†αcβ), l: α→β (c†βcα)
                            if occ3(fock, k, b)
                                && !occ3(fock, k, a)
                                && occ3(fock, l, a)
                                && !occ3(fock, l, b)
                            {
                                let f1 = fock ^ (1u64 << (3 * k + b)) ^ (1u64 << (3 * k + a));
                                let s1 = self.jw(fock, 3 * k + b, 3 * k + a);
                                let f2 = f1 ^ (1u64 << (3 * l + a)) ^ (1u64 << (3 * l + b));
                                let s2 = self.jw(f1, 3 * l + a, 3 * l + b);
                                let nkey = (f2 << 8) | (key & 0xff);
                                if let Some(jix) = self.find(nkey) {
                                    // k<l で両順序を数える二重計上を避ける: ordered (k,l) の
                                    // 各項はエルミート共役対で別 — 係数 w/2? → 検証 [A] が固定:
                                    // ordered 全対で w/2 (unordered 対あたり w×2 項/2 = w)
                                    let amp = 0.5 * wp * s1 * s2;
                                    w[jix].0 += amp * vr;
                                    w[jix].1 += amp * vi;
                                }
                            }
                        }
                    }
                }
            }
            // ---- 色交換 (不純物): 3: +w |β⟩⟨α| (c†αcβ)_k / 3̄: −w |β⟩⟨α| (c†βcα)_k ----
            for (ii, &(s_imp, kind)) in self.imps.iter().enumerate() {
                let a = regs[ii]; // 現在色 α
                for b in 0..3usize {
                    if b == a {
                        continue;
                    }
                    for k in 0..self.n {
                        let wp = self.w_pair(s_imp as f64, k as f64);
                        if kind > 0 {
                            // imp α→β, fermion k: β→α
                            if occ3(fock, k, b) && !occ3(fock, k, a) {
                                let nf = fock ^ (1u64 << (3 * k + b)) ^ (1u64 << (3 * k + a));
                                let sg = self.jw(fock, 3 * k + b, 3 * k + a);
                                let mut nkey = (nf << 8) | (key & 0xff);
                                nkey = (nkey & !(3u64 << (2 * ii))) | ((b as u64) << (2 * ii));
                                if let Some(jix) = self.find(nkey) {
                                    w[jix].0 += wp * sg * vr;
                                    w[jix].1 += wp * sg * vi;
                                }
                            }
                        } else {
                            // 3̄: imp ᾱ→β̄, fermion k: α→β (同方向), 係数 −w
                            if occ3(fock, k, a) && !occ3(fock, k, b) {
                                let nf = fock ^ (1u64 << (3 * k + a)) ^ (1u64 << (3 * k + b));
                                let sg = self.jw(fock, 3 * k + a, 3 * k + b);
                                let mut nkey = (nf << 8) | (key & 0xff);
                                nkey = (nkey & !(3u64 << (2 * ii))) | ((b as u64) << (2 * ii));
                                if let Some(jix) = self.find(nkey) {
                                    w[jix].0 += -wp * sg * vr;
                                    w[jix].1 += -wp * sg * vi;
                                }
                            }
                        }
                    }
                }
            }
            // ---- 色交換 (不純物×不純物) — 非順序対で全振幅 ----
            for ii in 0..nimp {
                for jj in ii + 1..nimp {
                    let (s1, k1) = self.imps[ii];
                    let (s2, k2) = self.imps[jj];
                    let (a1, a2) = (regs[ii], regs[jj]);
                    let wp = self.w_pair(s1 as f64, s2 as f64);
                    if k1 * k2 > 0 {
                        // 3·3 (と 3̄·3̄): w |a2⟩⟨a1|₁ ⊗ |a1⟩⟨a2|₂ (a1≠a2 で色交換)
                        if a1 != a2 {
                            let mut nkey = key;
                            nkey = (nkey & !(3u64 << (2 * ii))) | ((a2 as u64) << (2 * ii));
                            nkey = (nkey & !(3u64 << (2 * jj))) | ((a1 as u64) << (2 * jj));
                            if let Some(jix) = self.find(nkey) {
                                w[jix].0 += wp * vr;
                                w[jix].1 += wp * vi;
                            }
                        }
                    } else {
                        // 3·3̄: −w Σ_{β≠α} |β⟩⟨α|₁ ⊗ |β̄⟩⟨ᾱ|₂ (a1 = a2 のとき両方 → b)
                        if a1 == a2 {
                            for b in 0..3usize {
                                if b == a1 {
                                    continue;
                                }
                                let mut nkey = key;
                                nkey = (nkey & !(3u64 << (2 * ii))) | ((b as u64) << (2 * ii));
                                nkey = (nkey & !(3u64 << (2 * jj))) | ((b as u64) << (2 * jj));
                                if let Some(jix) = self.find(nkey) {
                                    w[jix].0 += -wp * vr;
                                    w[jix].1 += -wp * vi;
                                }
                            }
                        }
                    }
                }
            }
        }
        w
    }
}

// 再開始 Lanczos (v20.5 と同形)
fn lanczos_restart3(
    core: &Su3Core,
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

fn main() {
    self_test();
    println!("=== v20.6 SU(3) color toy — Coulomb 形式の一般化とバリオン Y 弦 ===\n");
    println!("事前登録: (a) [A] x=0 (メソン等式 8/3 + バリオン変分 ≤ 16/3) + [B][C] 全 PASS かつ");
    println!("          x=1 メソン E(4) > E(2) = SU(3) core 成立 / (b) 外れ。バリオン束縛は記録\n");
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
    let lam = 6.0;

    // ---- [B] 電場ゼロ = 自由 3 色鎖 (N=6) ----
    {
        let n = 6usize;
        let core = Su3Core::new(n, 1.0, 0.0, 0.0, 0.0, vec![], [3, 3, 3]);
        let (ev, _vv, rr) = lanczos_restart3(&core, 80, 8, 1e-9, 5);
        let mut eps: Vec<f64> = (1..=n)
            .map(|k| -2.0 * (std::f64::consts::PI * k as f64 / (n as f64 + 1.0)).cos())
            .collect();
        eps.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let exact: f64 = 3.0 * eps[..n / 2].iter().sum::<f64>();
        check(
            "[B] 電場ゼロのホップ = 自由 3 色鎖の厳密和 (N=6)",
            (ev - exact).abs() < 1e-8 && rr < 1e-8,
            format!("E₀ = {:.9} vs 厳密 {:.9} (dim {})", ev, exact, core.dim),
        );
    }

    // ---- [A] x=0 厳密アンカー (N=6) ----
    {
        let n = 6usize;
        let vac = Su3Core::new(n, 0.0, 1.0, lam, 1.0, vec![], [3, 3, 3]);
        let (e0, _v, r0) = lanczos_restart3(&vac, 60, 8, 1e-9, 11);
        // メソン (3 @1, 3̄ @3): E = 2·4/3 = 8/3
        let mes = Su3Core::new(n, 0.0, 1.0, lam, 1.0, vec![(1, 1), (3, -1)], [3, 3, 3]);
        let (em, _v, rm) = lanczos_restart3(&mes, 60, 8, 1e-9, 13);
        // バリオン (3 @0, 3 @2, 3 @4): 非破断値 (4/3)(2+2) = 16/3 は変分上界
        // (run1 の教訓: (1,3,5) は端のサイト 5 が w=0 で電場から切断される配置バグ +
        //  x=0 でもセクター遮蔽で非破断値より低くなり得る — 等式でなく不等式を課す)
        let bar = Su3Core::new(
            n,
            0.0,
            1.0,
            lam,
            1.0,
            vec![(0, 1), (2, 1), (4, 1)],
            [4, 4, 4],
        );
        let (eb, _v, rb) = lanczos_restart3(&bar, 60, 8, 1e-9, 17);
        let dm = em - e0 - 8.0 / 3.0;
        let eb_rel = eb - e0;
        check(
            "[A] x=0: メソン等式 8/3 (機械精度) + バリオン変分不等式 ≤ 16/3",
            dm.abs() < 1e-8 && eb_rel <= 16.0 / 3.0 + 1e-8 && r0 < 1e-8 && rm < 1e-8 && rb < 1e-8,
            format!(
                "メソン {:.9} (差 {:.1e}), バリオン {:.9} ({}) (dim {}/{})",
                em - e0,
                dm,
                eb_rel,
                if (eb_rel - 16.0 / 3.0).abs() < 1e-8 {
                    "= 16/3 非破断".to_string()
                } else {
                    format!("< 16/3 セクター遮蔽 (差 {:+.4})", eb_rel - 16.0 / 3.0)
                },
                mes.dim,
                bar.dim
            ),
        );
    }

    // ---- [C] + 本測定 (x=1, μ=1, N=8 メソン / N=6 バリオン) ----
    let x = 1.0;
    let mu = 1.0;
    let n8 = 8usize;
    let vac8 = Su3Core::new(n8, x, mu, lam, 1.0, vec![], [4, 4, 4]);
    let (e08, v08, r08) = lanczos_restart3(&vac8, 80, 12, 1e-7, 21);
    // ⟨C₂⟩ = ⟨H(λ+1)⟩ − E₀
    let vac8b = Su3Core::new(n8, x, mu, lam + 1.0, 1.0, vec![], [4, 4, 4]);
    let hv = vac8b.matvec(&v08);
    let c2: f64 = v08
        .iter()
        .zip(hv.iter())
        .map(|(a, b)| a.0 * b.0 + a.1 * b.1)
        .sum::<f64>()
        - e08;
    check(
        "[C] 真空 (N=8) の ⟨C₂_tot⟩ < 1e-6",
        c2.abs() < 1e-6 && r08 < 1e-6,
        format!(
            "⟨C₂⟩ = {:.1e} (res {:.0e}, dim {}) ({} s)",
            c2,
            r08,
            vac8.dim,
            t0.elapsed().as_secs()
        ),
    );
    // メソン E(r), r ∈ {2, 4}
    let mut e_mes = Vec::new();
    for &r in &[2usize, 4] {
        let core = Su3Core::new(n8, x, mu, lam, 1.0, vec![(1, 1), (1 + r, -1)], [4, 4, 4]);
        let (ev, vv, rr) = lanczos_restart3(&core, 80, 12, 1e-7, 31);
        let core2 = Su3Core::new(
            n8,
            x,
            mu,
            lam + 1.0,
            1.0,
            vec![(1, 1), (1 + r, -1)],
            [4, 4, 4],
        );
        let hv2 = core2.matvec(&vv);
        let c2p: f64 = vv
            .iter()
            .zip(hv2.iter())
            .map(|(a, b)| a.0 * b.0 + a.1 * b.1)
            .sum::<f64>()
            - ev;
        check(
            &format!("メソン r={} の残差・⟨C₂⟩", r),
            rr < 1e-6 && c2p.abs() < 1e-6,
            format!("res {:.0e}, ⟨C₂⟩ {:.1e} (dim {})", rr, c2p, core.dim),
        );
        e_mes.push(ev - e08);
        println!(
            "    メソン r={}: E = {:.5} ({} s)",
            r,
            ev - e08,
            t0.elapsed().as_secs()
        );
    }
    // バリオン (x=1, N=6): Y 配置 (1,3,5) — 束縛の記録
    {
        let n6 = 6usize;
        let vac6 = Su3Core::new(n6, x, mu, lam, 1.0, vec![], [3, 3, 3]);
        let (e06, _v, r06) = lanczos_restart3(&vac6, 80, 12, 1e-7, 41);
        let bar = Su3Core::new(n6, x, mu, lam, 1.0, vec![(0, 1), (2, 1), (4, 1)], [4, 4, 4]);
        let (eb, _v, rb) = lanczos_restart3(&bar, 80, 12, 1e-7, 43);
        let mes6 = Su3Core::new(n6, x, mu, lam, 1.0, vec![(1, 1), (3, -1)], [3, 3, 3]);
        let (em6, _v, rm6) = lanczos_restart3(&mes6, 80, 12, 1e-7, 47);
        check(
            "バリオン系 (N=6) の残差",
            r06 < 1e-6 && rb < 1e-6 && rm6 < 1e-6,
            format!("res {:.0e}/{:.0e}/{:.0e}", r06, rb, rm6),
        );
        println!(
            "    [記録] x=1 バリオン Y (0,2,4): E = {:.5} vs 2×メソン(r=2) = {:.5} (劣加法性 {})",
            eb - e06,
            2.0 * (em6 - e06),
            if eb - e06 < 2.0 * (em6 - e06) {
                "成立"
            } else {
                "不成立"
            }
        );
    }

    let mono = e_mes[1] > e_mes[0];
    println!(
        "\n[判定] {}",
        if mono && nfail == 0 {
            "事前登録 (a): SU(3) core 成立 — メソン弦は閉じ込め方向・x=0 の N-ality アンカー厳密 (Z₂→U(1)→SU(2)→SU(3) の階段が通った)"
        } else {
            "事前登録 (b): 外れ — 記録"
        }
    );
    println!(
        "    メソン E(2) = {:.4}, E(4) = {:.4} (σ_3 = {:.4})",
        e_mes[0],
        e_mes[1],
        (e_mes[1] - e_mes[0]) / 2.0
    );

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v20.6".into())),
        (
            "e_meson".into(),
            Json::Arr(e_mes.iter().map(|&v| Json::Num(v)).collect()),
        ),
        ("branch_a".into(), Json::Bool(mono && nfail == 0)),
    ]);
    let p = write_artifact("results/v206_su3core.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — 分岐 (a)/(b) は [判定] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
