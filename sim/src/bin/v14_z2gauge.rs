//! v1.4 ゲージ場は創発する — 3D Z2 格子ゲージ理論の閉じ込め/非閉じ込め転移
//!
//! 問い: ゲージ場 (光子・グルーオン) は基本的実在か、それとも量子ビットの海の集団現象か。
//! 3D Z2 格子ゲージ理論 (最も単純なゲージ理論) には 2 つの相がある:
//!   強結合 (β小): 閉じ込め相 — Wilson ループは面積則
//!   弱結合 (β大): 非閉じ込め相 (トーリック符号 = トポロジカル秩序) — 周長則
//! 転移点は 3D Ising との厳密な双対性で β_c = -½ ln tanh(0.221654) = 0.7613。
//!
//! 非閉じ込め相こそ「創発ゲージ場」の相: 微視的にはただのスピン系なのに、
//! 低エネルギーには変形不能なループ (電束) とゲージ構造が現れる (Wen の弦ネット凝縮)。
//! 双対性はもう一つの教訓: ゲージ理論とスピン系は同じ理論の二つの読み (群盲象の再演)。

use uft_sim::*;

const L: usize = 12;
const NS: usize = L * L * L;

struct Z2 {
    link: Vec<i8>, // link[dir*NS + site]
    beta: f64,
    rng: Rng,
}

fn site(x: usize, y: usize, z: usize) -> usize {
    (x % L) + (y % L) * L + (z % L) * L * L
}

impl Z2 {
    fn new(beta: f64, seed: u64) -> Self {
        Z2 {
            link: vec![1; 3 * NS],
            beta,
            rng: Rng::new(seed),
        }
    }
    fn coords(s: usize) -> (usize, usize, usize) {
        (s % L, (s / L) % L, s / (L * L))
    }
    fn shift(s: usize, dir: usize, by: usize) -> usize {
        let (x, y, z) = Self::coords(s);
        match dir {
            0 => site(x + by, y, z),
            1 => site(x, y + by, z),
            _ => site(x, y, z + by),
        }
    }
    /// プラケット (s; μ<ν 平面) の値
    fn plaq(&self, s: usize, mu: usize, nu: usize) -> i32 {
        let u1 = self.link[mu * NS + s];
        let u2 = self.link[nu * NS + Self::shift(s, mu, 1)];
        let u3 = self.link[mu * NS + Self::shift(s, nu, 1)];
        let u4 = self.link[nu * NS + s];
        (u1 * u2 * u3 * u4) as i32
    }
    /// リンク (dir, s) を含む 4 枚のプラケットの和
    fn staple_sum(&self, dir: usize, s: usize) -> i32 {
        let mut sum = 0;
        for nu in 0..3 {
            if nu == dir {
                continue;
            }
            let (mu2, nu2) = if dir < nu { (dir, nu) } else { (nu, dir) };
            sum += self.plaq(s, mu2, nu2);
            let sm = Self::shift(s, nu, L - 1); // s - ν
            sum += self.plaq(sm, mu2, nu2);
        }
        sum
    }
    fn sweep(&mut self) {
        for _ in 0..3 * NS {
            let r = self.rng.range(3 * NS);
            let (dir, s) = (r / NS, r % NS);
            let ds = 2.0 * self.beta * self.staple_sum(dir, s) as f64;
            if ds <= 0.0 || self.rng.f64() < (-ds).exp() {
                self.link[r] = -self.link[r];
            }
        }
    }
    fn mean_plaq(&self) -> f64 {
        let mut s = 0.0;
        for n in 0..NS {
            s += (self.plaq(n, 0, 1) + self.plaq(n, 0, 2) + self.plaq(n, 1, 2)) as f64;
        }
        s / (3 * NS) as f64
    }
    /// Wilson ループ W(R,T): 全位置・3 平面平均
    fn wilson(&self, r: usize, t: usize) -> f64 {
        let mut sum = 0.0;
        let planes = [(0usize, 1usize), (0, 2), (1, 2)];
        for &(mu, nu) in &planes {
            for s0 in 0..NS {
                let mut p = 1i32;
                let mut s = s0;
                for _ in 0..r {
                    p *= self.link[mu * NS + s] as i32;
                    s = Self::shift(s, mu, 1);
                }
                for _ in 0..t {
                    p *= self.link[nu * NS + s] as i32;
                    s = Self::shift(s, nu, 1);
                }
                for _ in 0..r {
                    s = Self::shift(s, mu, L - 1);
                    p *= self.link[mu * NS + s] as i32;
                }
                for _ in 0..t {
                    s = Self::shift(s, nu, L - 1);
                    p *= self.link[nu * NS + s] as i32;
                }
                sum += p as f64;
            }
        }
        sum / (3 * NS) as f64
    }
}

fn main() {
    println!(
        "=== v1.4 創発ゲージ場: 3D Z2 格子ゲージ理論 (L={}) ===\n",
        L
    );
    let beta_c = -0.5 * (0.221654f64.tanh()).ln();
    println!(
        "3D Ising 双対性による厳密な転移点: β_c = -½ ln tanh(0.221654) = {:.4}\n",
        beta_c
    );
    println!(
        "  β      ⟨P⟩      χ(2)      χ(3)     [χ(R)=Creutz比→弦張力σ]   強結合予言 -ln tanh β"
    );
    for (bi, &beta) in [0.60f64, 0.70, 0.74, 0.7613, 0.80, 0.90].iter().enumerate() {
        let mut sys = Z2::new(beta, 31337 + bi as u64);
        for _ in 0..1500 {
            sys.sweep();
        }
        let mut w = std::collections::HashMap::new();
        let sizes: [(usize, usize); 6] = [(1, 1), (2, 1), (2, 2), (3, 2), (3, 3), (4, 3)];
        let mut plq = Vec::new();
        let nmeas = 1500;
        for _ in 0..nmeas {
            for _ in 0..3 {
                sys.sweep();
            }
            plq.push(sys.mean_plaq());
            for &(r, t) in &sizes {
                *w.entry((r, t)).or_insert(0.0) += sys.wilson(r, t) / nmeas as f64;
            }
        }
        let (mp, _) = mean_err(&plq);
        let chi = |r: usize| -> f64 {
            let a = w[&(r, r)] * w[&(r - 1, r - 1)];
            let b = w[&(r, r - 1)] * w[&(r, r - 1)];
            if a > 0.0 && b > 0.0 {
                -(a / b).ln()
            } else {
                f64::NAN
            }
        };
        println!(
            "  {:.4} {:.4}   {:7.4}   {:7.4}                             {:.4}",
            beta,
            mp,
            chi(2),
            chi(3),
            -(beta.tanh().ln())
        );
    }
    println!("\n  読み方: χ(R) が R によらず正の一定値 → 面積則 (閉じ込め相)。");
    println!("          χ(3) ≪ χ(2) → 周長則へ (非閉じ込め相 = トポロジカル相 = 創発ゲージ場)。");
    println!(
        "          β_c={:.3} を境に χ(3) が急落することを確認せよ。",
        beta_c
    );
    println!("\n結論: ゲージ構造は基本法則に書き込まなくても、量子ビットの海の相として創発しうる");
    println!(
        "      (Wen: 光子は弦ネット凝縮の集団励起)。さらに双対性 (ゲージ理論 ⇔ Ising 磁石) は、"
    );
    println!(
        "      「同じ実在の異なる記述」が可能なことを示す — 基本変数は一意でない (群盲象の教訓)。"
    );
    println!("      QRN 公理 A4 の強化: ゲージ場・物質はネットワークの相と励起である。");
}
