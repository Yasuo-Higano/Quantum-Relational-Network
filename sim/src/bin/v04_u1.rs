//! v0.4 ゲージ原理 — 2次元 U(1) 格子ゲージ理論 (コンパクト QED)
//! リンク変数 U = e^{iθ} ∈ U(1)、Wilson 作用 S = β Σ_p (1 - cos θ_p)
//! θ_p = θx(n) + θy(n+x̂) - θx(n+ŷ) - θy(n)  (プラケット = 最小の曲率)
//!
//! 2D では厳密解が知られる:
//!   ⟨cos θ_p⟩ = I₁(β)/I₀(β),  Wilson ループ W(R,T) = [I₁(β)/I₀(β)]^{RT}
//! → 面積則 = 閉じ込め (クォークを引き離すコストが距離に比例)
//!
//! ゲージ場とは「内部空間の接続」であり、プラケットは「内部空間の曲率」。
//! GR の曲率 (v0.1) と同じ幾何学的構造が、電磁気・核力の正体である。

use uft_sim::*;

struct U1 {
    l: usize,
    beta: f64,
    th: Vec<f64>, // th[dir*L*L + site], dir 0=x, 1=y
    rng: Rng,
    step: f64,
}

impl U1 {
    fn new(l: usize, beta: f64, seed: u64) -> Self {
        U1 {
            l,
            beta,
            th: vec![0.0; 2 * l * l],
            rng: Rng::new(seed),
            step: 1.0 / (1.0 + beta).sqrt(),
        }
    }
    fn site(&self, x: usize, y: usize) -> usize {
        (x % self.l) + (y % self.l) * self.l
    }
    /// プラケット角 θ_p at (x,y)
    fn plaq(&self, x: usize, y: usize) -> f64 {
        let n = self.l * self.l;
        self.th[self.site(x, y)]
            + self.th[n + self.site(x + 1, y)]
            - self.th[self.site(x, y + 1)]
            - self.th[n + self.site(x, y)]
    }
    /// リンク (dir, x, y) を含む2枚のプラケットの cos 和
    fn staple_action(&self, dir: usize, x: usize, y: usize) -> f64 {
        let l = self.l;
        match dir {
            0 => self.plaq(x, y).cos() + self.plaq(x, y + l - 1).cos(),
            _ => self.plaq(x, y).cos() + self.plaq(x + l - 1, y).cos(),
        }
    }
    fn sweep(&mut self) -> f64 {
        let n = self.l * self.l;
        let mut acc = 0;
        for _ in 0..2 * n {
            let r = self.rng.range(2 * n);
            let (dir, s) = (r / n, r % n);
            let (x, y) = (s % self.l, s / self.l);
            let old = self.th[r];
            let s_old = self.staple_action(dir, x, y);
            self.th[r] = old + self.step * (2.0 * self.rng.f64() - 1.0);
            let s_new = self.staple_action(dir, x, y);
            let ds = -self.beta * (s_new - s_old);
            if ds <= 0.0 || self.rng.f64() < (-ds).exp() {
                acc += 1;
            } else {
                self.th[r] = old;
            }
        }
        acc as f64 / (2 * n) as f64
    }
    fn mean_plaq(&self) -> f64 {
        let mut s = 0.0;
        for y in 0..self.l {
            for x in 0..self.l {
                s += self.plaq(x, y).cos();
            }
        }
        s / (self.l * self.l) as f64
    }
    /// Wilson ループ W(R,T) の全平行移動平均
    fn wilson(&self, r: usize, t: usize) -> f64 {
        let n = self.l * self.l;
        let mut sum = 0.0;
        for y0 in 0..self.l {
            for x0 in 0..self.l {
                let mut a = 0.0;
                for i in 0..r {
                    a += self.th[self.site(x0 + i, y0)]; // 下辺 →
                    a -= self.th[self.site(x0 + i, y0 + t)]; // 上辺 ←
                }
                for j in 0..t {
                    a += self.th[n + self.site(x0 + r, y0 + j)]; // 右辺 ↑
                    a -= self.th[n + self.site(x0, y0 + j)]; // 左辺 ↓
                }
                sum += a.cos();
            }
        }
        sum / n as f64
    }
}

fn main() {
    println!("=== v0.4 ゲージ原理: 2D U(1) 格子ゲージ理論 ===\n");

    println!("[A] プラケット平均 (内部曲率の期待値) vs 厳密解 I₁(β)/I₀(β)");
    println!("  β     MC結果          厳密解");
    for (bi, &beta) in [0.5f64, 1.0, 2.0, 4.0].iter().enumerate() {
        let mut sys = U1::new(64, beta, 9000 + bi as u64);
        for _ in 0..1000 {
            sys.sweep();
        }
        let mut ps = Vec::new();
        for _ in 0..4000 {
            sys.sweep();
            ps.push(sys.mean_plaq());
        }
        let (m, e) = mean_err(&ps);
        let exact = bessel_i(1, beta) / bessel_i(0, beta);
        let ok = (m - exact).abs() < 4.0 * e.max(2e-4);
        println!("  {:.1}  {:.5}±{:.5}  {:.5}  {}", beta, m, e, exact, pass(ok));
    }

    println!("\n[B] Wilson ループと閉じ込め (β=2.0)");
    println!("    W(R,T) = ⟨e^{{i∮A·dl}}⟩ : 電荷対を距離 R で時間 T だけ保持する振幅");
    let beta = 2.0f64;
    let mut sys = U1::new(64, beta, 12321);
    for _ in 0..2000 {
        sys.sweep();
    }
    let loops: [(usize, usize); 7] = [(1, 1), (2, 1), (2, 2), (3, 2), (3, 3), (4, 3), (4, 4)];
    let mut sums = vec![Vec::new(); loops.len()];
    for _ in 0..6000 {
        sys.sweep();
        for (i, &(r, t)) in loops.iter().enumerate() {
            sums[i].push(sys.wilson(r, t));
        }
    }
    let ratio = bessel_i(1, beta) / bessel_i(0, beta);
    println!("  R×T   面積  W(MC)             W(厳密)=(I₁/I₀)^A");
    let mut areas = Vec::new();
    let mut lnw = Vec::new();
    for (i, &(r, t)) in loops.iter().enumerate() {
        let (m, e) = mean_err(&sums[i]);
        let exact = ratio.powi((r * t) as i32);
        let ok = (m - exact).abs() < 4.0 * e.max(5e-4);
        println!(
            "  {}×{}   {:2}   {:.5}±{:.5}  {:.5}  {}",
            r,
            t,
            r * t,
            m,
            e,
            exact,
            pass(ok)
        );
        if m > 0.0 {
            areas.push((r * t) as f64);
            lnw.push(m.ln());
        }
    }
    let (_, slope) = linfit(&areas, &lnw);
    println!(
        "  面積則フィット: -ln W = σ·A, 弦張力 σ(MC) = {:.4}  厳密 -ln(I₁/I₀) = {:.4}",
        -slope,
        -ratio.ln()
    );
    println!("  => W は「周長」でなく「面積」で減衰 = 電荷間のエネルギーが距離に比例して増える");
    println!("     = 閉じ込め。2D U(1) は QCD のクォーク閉じ込めと同じ機構を厳密に見せる。\n");

    println!("[C] ゲージ対称性の意味の確認 — ランダムゲージ変換で観測量は不変か");
    {
        let l = sys.l;
        let n = l * l;
        let w_before = sys.wilson(3, 3);
        let p_before = sys.mean_plaq();
        // ゲージ変換: θ_μ(x) → θ_μ(x) + α(x) - α(x+μ̂)
        let mut rng2 = Rng::new(555);
        let alpha: Vec<f64> = (0..n).map(|_| rng2.f64() * 6.283).collect();
        for y in 0..l {
            for x in 0..l {
                let s = sys.site(x, y);
                let sx = sys.site(x + 1, y);
                let sy = sys.site(x, y + 1);
                sys.th[s] += alpha[s] - alpha[sx];
                sys.th[n + s] += alpha[s] - alpha[sy];
            }
        }
        let w_after = sys.wilson(3, 3);
        let p_after = sys.mean_plaq();
        println!("  プラケット: {:+.6} → {:+.6}  (差 {:.1e})", p_before, p_after, (p_after - p_before).abs());
        println!("  W(3,3)   : {:+.6} → {:+.6}  (差 {:.1e})", w_before, w_after, (w_after - w_before).abs());
        println!("  => リンク変数(ポテンシャル A)は物理でない。閉ループ(ホロノミー=曲率)だけが物理。");
    }
    println!("\n結論: 力とは接続の曲率である。電磁気も核力も(そして v0.1 の重力も)同じ「幾何」の言語で書ける。");
    println!("      ゲージ原理は標準模型と GR を貫く唯一の設計原理 — 統一の最有力な手がかり。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
