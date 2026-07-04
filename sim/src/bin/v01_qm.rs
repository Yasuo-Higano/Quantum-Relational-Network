//! v0.1 量子力学の数値検証 — 1D シュレーディンガー方程式 (Crank–Nicolson 法)
//! ħ = m = 1 の単位系。CN 法は厳密にユニタリー(ノルム保存)。
//!
//! A: 自由波束の拡散     σ(t) = σ0 √(1+(t/(2σ0²))²) と比較
//! B: 矩形障壁のトンネル  透過率を平面波厳密式の運動量平均と比較
//! C: 調和振動子のコヒーレント状態  ⟨x⟩(t) = x0 cos t と比較

use uft_sim::*;

struct Cn {
    n: usize,
    dx: f64,
    dt: f64,
    x0: f64,
    v: Vec<f64>,
    psi: Vec<C64>,
}

impl Cn {
    fn new(x0: f64, x1: f64, n: usize, dt: f64, pot: impl Fn(f64) -> f64) -> Self {
        let dx = (x1 - x0) / n as f64;
        let v = (0..n).map(|i| pot(x0 + i as f64 * dx)).collect();
        Cn {
            n,
            dx,
            dt,
            x0,
            v,
            psi: vec![CZERO; n],
        }
    }
    fn x(&self, i: usize) -> f64 {
        self.x0 + i as f64 * self.dx
    }
    fn set_gaussian(&mut self, xc: f64, sigma: f64, k0: f64) {
        for i in 0..self.n {
            let x = self.x(i);
            let g = (-(x - xc) * (x - xc) / (4.0 * sigma * sigma)).exp();
            self.psi[i] = C64::expi(k0 * x).scale(g);
        }
        self.normalize();
    }
    fn normalize(&mut self) {
        let s: f64 = self.psi.iter().map(|p| p.norm2()).sum::<f64>() * self.dx;
        let a = 1.0 / s.sqrt();
        for p in self.psi.iter_mut() {
            *p = p.scale(a);
        }
    }
    fn norm(&self) -> f64 {
        self.psi.iter().map(|p| p.norm2()).sum::<f64>() * self.dx
    }
    /// (I + i dt H/2) ψ' = (I - i dt H/2) ψ を1ステップ解く
    fn step(&mut self) {
        let n = self.n;
        let lam = self.dt / (4.0 * self.dx * self.dx); // dt/2 * 1/(2dx²)
        let mut a = vec![CZERO; n];
        let mut b = vec![CZERO; n];
        let mut c = vec![CZERO; n];
        let mut d = vec![CZERO; n];
        for i in 0..n {
            // H ψ = -ψ''/2 + V ψ,  H_ii = 1/dx² + V_i, H_i,i±1 = -1/(2dx²)
            let hd = 1.0 / (self.dx * self.dx) + self.v[i];
            b[i] = CONE + CI.scale(self.dt * 0.5 * hd);
            if i > 0 {
                a[i] = -CI.scale(lam);
            }
            if i < n - 1 {
                c[i] = -CI.scale(lam);
            }
            let mut rhs = (CONE - CI.scale(self.dt * 0.5 * hd)) * self.psi[i];
            if i > 0 {
                rhs = rhs + CI.scale(lam) * self.psi[i - 1];
            }
            if i < n - 1 {
                rhs = rhs + CI.scale(lam) * self.psi[i + 1];
            }
            d[i] = rhs;
        }
        self.psi = solve_tridiag_c(&a, &b, &c, &d);
    }
    fn moments(&self) -> (f64, f64) {
        let mut m1 = 0.0;
        let mut m2 = 0.0;
        for i in 0..self.n {
            let p = self.psi[i].norm2() * self.dx;
            let x = self.x(i);
            m1 += x * p;
            m2 += x * x * p;
        }
        (m1, (m2 - m1 * m1).sqrt())
    }
    fn prob_right_of(&self, xb: f64) -> f64 {
        (0..self.n)
            .filter(|&i| self.x(i) > xb)
            .map(|i| self.psi[i].norm2() * self.dx)
            .sum()
    }
}

fn main() {
    self_test();
    println!("=== v0.1 量子力学: 1D Schrödinger (Crank–Nicolson) ===\n");

    // ---- A: 自由波束の拡散 ----
    println!("[A] 自由波束の拡散  σ(t) = σ0·sqrt(1+(t/(2σ0²))²)");
    let sigma0 = 2.0;
    let mut cn = Cn::new(-100.0, 100.0, 2000, 0.01, |_| 0.0);
    cn.set_gaussian(0.0, sigma0, 0.0);
    let mut max_rel = 0.0f64;
    println!("  t      σ(数値)   σ(理論)   相対誤差");
    for step in 0..=4000 {
        if step % 1000 == 0 {
            let t = step as f64 * cn.dt;
            let (_, s) = cn.moments();
            let s_th = sigma0 * (1.0 + (t / (2.0 * sigma0 * sigma0)).powi(2)).sqrt();
            let rel = ((s - s_th) / s_th).abs();
            max_rel = max_rel.max(rel);
            println!("  {:5.1}  {:8.5}  {:8.5}  {:.2e}", t, s, s_th, rel);
        }
        if step < 4000 {
            cn.step();
        }
    }
    println!("  ノルム保存: |1-norm| = {:.2e}", (cn.norm() - 1.0).abs());
    println!(
        "  => 最大相対誤差 {:.2e}  {}\n",
        max_rel,
        pass(max_rel < 5e-3)
    );

    // ---- B: トンネル効果 (Eckart 障壁 V = V0/cosh²x — 厳密透過率が既知) ----
    println!("[B] Eckart 障壁トンネル  V=1/cosh²(x), E≈0.5 (古典的には全反射)");
    let (k0, sig) = (1.0, 10.0);
    let mut cn = Cn::new(-300.0, 300.0, 6000, 0.02, |x: f64| 1.0 / x.cosh().powi(2));
    cn.set_gaussian(-60.0, sig, k0);
    for _ in 0..9000 {
        cn.step(); // t = 180: 透過波は x≈+120, 反射波は x≈-240
    }
    let t_num = cn.prob_right_of(10.0);
    // 厳密透過率 T(k) = sinh²(πk) / [sinh²(πk) + cosh²(π√7/2)] を波束の運動量分布で平均
    let sig_k = 1.0 / (2.0 * sig);
    let coshterm = (std::f64::consts::PI * (7.0f64).sqrt() / 2.0)
        .cosh()
        .powi(2);
    let mut t_th = 0.0;
    let mut wsum = 0.0;
    let nk = 4000;
    for j in 0..nk {
        let k: f64 = k0 - 5.0 * sig_k + 10.0 * sig_k * (j as f64 + 0.5) / nk as f64;
        if k <= 0.0 {
            continue;
        }
        let w = (-(k - k0) * (k - k0) / (2.0 * sig_k * sig_k)).exp();
        let sh = (std::f64::consts::PI * k).sinh().powi(2);
        t_th += w * sh / (sh + coshterm);
        wsum += w;
    }
    t_th /= wsum;
    let rel = ((t_num - t_th) / t_th).abs();
    println!(
        "  透過率: 数値 {:.4}  厳密解(平面波平均) {:.4}  相対差 {:.1}%",
        t_num,
        t_th,
        rel * 100.0
    );
    println!("  ノルム保存: |1-norm| = {:.2e}", (cn.norm() - 1.0).abs());
    println!(
        "  => 古典禁止領域を {:.1}% が通過(トンネル効果) {}\n",
        t_num * 100.0,
        pass(rel < 0.05)
    );

    // ---- C: 調和振動子コヒーレント状態 ----
    println!("[C] 調和振動子 V=x²/2 のコヒーレント状態  ⟨x⟩(t) = 3·cos t");
    let mut cn = Cn::new(-15.0, 15.0, 2000, 0.002, |x| 0.5 * x * x);
    cn.set_gaussian(3.0, 1.0 / (2.0f64).sqrt(), 0.0); // 基底状態幅 σ=1/√2 を x=3 へ変位
    let mut max_err = 0.0f64;
    let steps = (2.0 * std::f64::consts::PI / 0.002) as usize * 2; // 2 周期
    for step in 0..=steps {
        let t = step as f64 * 0.002;
        let (m1, _) = cn.moments();
        let err = (m1 - 3.0 * t.cos()).abs();
        max_err = max_err.max(err);
        if step % (steps / 8) == 0 {
            println!(
                "  t={:6.3}  ⟨x⟩={:8.4}  3cos(t)={:8.4}",
                t,
                m1,
                3.0 * t.cos()
            );
        }
        if step < steps {
            cn.step();
        }
    }
    println!(
        "  => 2周期での最大誤差 {:.2e}  {}",
        max_err,
        pass(max_err < 1e-2)
    );
    println!("\n結論: ユニタリー時間発展・重ね合わせ・トンネルという QM の核心を数値で確認。");
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}
