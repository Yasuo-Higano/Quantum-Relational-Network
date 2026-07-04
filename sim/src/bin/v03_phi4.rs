//! v0.3 場と粒子の創発 — 2次元 φ⁴ 格子場理論
//! S = Σ_x [ -2κ φ_x Σ_μ φ_{x+μ̂} + φ_x² + λ(φ_x²-1)² ],  λ = 1.0
//!
//! v0.2 の「1次元弾性鎖」を2次元の膜にしたものが場の量子論。
//! κ (ホッピング≒逆格子間隔²) を変えると Z2 対称相 ↔ 自発的破れ相の相転移が起き、
//! 臨界点近傍で相関長 → ∞、そこで連続場理論 (質量ゼロに近いスカラー場) が創発する。
//!
//! 更新: Metropolis (|φ|) + 埋め込み Ising の Wolff クラスター (符号) — 臨界減速を回避

use uft_sim::*;

struct Phi4 {
    l: usize,
    kappa: f64,
    lambda: f64,
    phi: Vec<f64>,
    rng: Rng,
    step: f64,
}

impl Phi4 {
    fn new(l: usize, kappa: f64, lambda: f64, seed: u64) -> Self {
        Phi4 {
            l,
            kappa,
            lambda,
            phi: vec![1.0; l * l],
            rng: Rng::new(seed),
            step: 0.8,
        }
    }
    fn idx(&self, x: usize, y: usize) -> usize {
        (x % self.l) + (y % self.l) * self.l
    }
    fn neighbors(&self, i: usize) -> [usize; 4] {
        let (x, y) = (i % self.l, i / self.l);
        let l = self.l;
        [
            self.idx(x + 1, y),
            self.idx(x + l - 1, y),
            self.idx(x, y + 1),
            self.idx(x, y + l - 1),
        ]
    }
    fn metropolis_sweep(&mut self) -> f64 {
        let n = self.l * self.l;
        let mut acc = 0;
        for _ in 0..n {
            let i = self.rng.range(n);
            let old = self.phi[i];
            let new = old + self.step * (2.0 * self.rng.f64() - 1.0);
            let nb: f64 = self.neighbors(i).iter().map(|&j| self.phi[j]).sum();
            let s_old =
                -2.0 * self.kappa * old * nb + old * old + self.lambda * (old * old - 1.0).powi(2);
            let s_new =
                -2.0 * self.kappa * new * nb + new * new + self.lambda * (new * new - 1.0).powi(2);
            let ds = s_new - s_old;
            if ds <= 0.0 || self.rng.f64() < (-ds).exp() {
                self.phi[i] = new;
                acc += 1;
            }
        }
        acc as f64 / n as f64
    }
    /// 埋め込み Ising Wolff クラスター: 同符号ボンドを p=1-exp(-4κ|φφ'|) で結び符号反転
    fn wolff(&mut self) -> usize {
        let n = self.l * self.l;
        let seed = self.rng.range(n);
        let mut in_cluster = vec![false; n];
        let mut stack = vec![seed];
        in_cluster[seed] = true;
        let mut size = 0;
        while let Some(i) = stack.pop() {
            size += 1;
            for &j in self.neighbors(i).iter() {
                if !in_cluster[j] && self.phi[i] * self.phi[j] > 0.0 {
                    let p = 1.0 - (-4.0 * self.kappa * (self.phi[i] * self.phi[j]).abs()).exp();
                    if self.rng.f64() < p {
                        in_cluster[j] = true;
                        stack.push(j);
                    }
                }
            }
        }
        for i in 0..n {
            if in_cluster[i] {
                self.phi[i] = -self.phi[i];
            }
        }
        size
    }
    fn magnetization(&self) -> f64 {
        self.phi.iter().sum::<f64>() / (self.l * self.l) as f64
    }
}

fn run_point(l: usize, kappa: f64, seed: u64) -> (f64, f64, f64, f64) {
    let mut sys = Phi4::new(l, kappa, 1.0, seed);
    for _ in 0..3000 {
        sys.metropolis_sweep();
        sys.wolff();
        sys.wolff();
    }
    let nmeas = 40000;
    let (mut s_absm, mut s_m2, mut s_m4) = (0.0, 0.0, 0.0);
    for _ in 0..nmeas {
        sys.metropolis_sweep();
        sys.wolff();
        sys.wolff();
        let m = sys.magnetization();
        s_absm += m.abs();
        s_m2 += m * m;
        s_m4 += m * m * m * m;
    }
    let (absm, m2, m4) = (
        s_absm / nmeas as f64,
        s_m2 / nmeas as f64,
        s_m4 / nmeas as f64,
    );
    let binder = 1.0 - m4 / (3.0 * m2 * m2);
    let chi = (l * l) as f64 * (m2 - absm * absm);
    (absm, m2, binder, chi)
}

fn main() {
    println!("=== v0.3 場の量子論の創発: 2D φ⁴ 格子理論 (λ=1) ===\n");
    println!("[A] κ スキャン — Z2 対称相から自発的対称性の破れ (SSB) へ");
    let kappas = [0.20, 0.27, 0.30, 0.32, 0.33, 0.34, 0.35, 0.36];
    let mut binders: Vec<(f64, f64, f64)> = Vec::new(); // (κ, U16, U32)
    for (li, &l) in [16usize, 32].iter().enumerate() {
        println!("  L = {}", l);
        println!("  κ      ⟨|m|⟩    χ        Binder U");
        for (ki, &k) in kappas.iter().enumerate() {
            let (absm, _m2, u, chi) = run_point(l, k, 7000 + (li * 100 + ki) as u64);
            println!("  {:.3}  {:.4}  {:8.2}  {:.4}", k, absm, chi, u);
            if li == 0 {
                binders.push((k, u, 0.0));
            } else {
                binders[ki].2 = u;
            }
        }
    }
    println!("\n[B] Binder 累積量 U(L) の交差 = 臨界点 (熱力学極限の相転移の指紋)");
    println!("  κ      U(16)   U(32)   U16-U32");
    // 対称相では U≈0 (ノイズ)、破れ相では両者 2/3 に合流。
    // 交差 = U32-U16 が負から正へ転じる最後の点 (臨界点で大きい L ほど U が急峻に立ち上がる)
    let mut kc = f64::NAN;
    for w in binders.windows(2) {
        let (k0, u16_0, u32_0) = w[0];
        let (k1, u16_1, u32_1) = w[1];
        let d0 = u32_0 - u16_0;
        let d1 = u32_1 - u16_1;
        if d0 < 0.0 && d1 >= 0.0 && (u16_1 + u32_1) > 0.3 {
            kc = k0 + (k1 - k0) * (-d0) / (d1 - d0);
        }
    }
    for &(k, u16, u32) in &binders {
        println!("  {:.3}  {:.4}  {:.4}  {:+.4}", k, u16, u32, u16 - u32);
    }
    println!("  => 臨界点 κ_c ≈ {:.4} (線形補間)", kc);
    println!("     臨界点では相関長 ξ→∞: 格子の詳細が消え、連続なスカラー場理論が創発する。");
    println!("     (2D φ⁴ は Ising 普遍類: 臨界指数はミクロな作りに依らない)\n");

    println!("[C] 対称性の破れの実像 — 深い破れ相 κ=0.42, L=32 での秩序変数分布");
    {
        let mut sys = Phi4::new(32, 0.42, 1.0, 4242);
        for _ in 0..2000 {
            sys.metropolis_sweep();
            sys.wolff();
        }
        let mut hist = [0usize; 9];
        let nmeas = 20000;
        for _ in 0..nmeas {
            sys.metropolis_sweep();
            sys.wolff();
            let m = sys.magnetization();
            let bin = (((m + 1.8) / 0.4) as isize).clamp(0, 8) as usize;
            hist[bin] += 1;
        }
        println!("  m のヒストグラム (-1.8..1.8, 幅0.4):");
        for (b, &h) in hist.iter().enumerate() {
            let center = -1.6 + 0.4 * b as f64;
            let bar = "#".repeat((h * 60 / nmeas).min(60));
            println!("  {:+.1}: {:5} {}", center, h, bar);
        }
        println!("  => 二峰分布: 真空 (最低エネルギー状態) が2つある。系はどちらかを「選ぶ」。");
        println!("     これがヒッグス機構の原型 (SSB)。粒子の質量は真空の選択から生まれる。");
    }
    println!("\n結論: 「場」の量子論は経路積分格子の統計力学であり、粒子=場の励起、");
    println!("      真空=最低エネルギー状態、質量=対称性の破れ、連続時空=臨界点の普遍性として全て創発する。");
}
