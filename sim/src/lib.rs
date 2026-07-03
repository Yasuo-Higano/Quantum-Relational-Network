//! 統一場理論プロジェクト共通数値ライブラリ(外部依存なし)

use std::ops::{Add, Div, Mul, Neg, Sub};

// ---------------- 乱数 (xorshift64*) ----------------
pub struct Rng {
    s: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng {
            s: seed.wrapping_mul(0x9E3779B97F4A7C15) | 1,
        }
    }
    pub fn u64(&mut self) -> u64 {
        let mut x = self.s;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.s = x;
        x.wrapping_mul(0x2545F4914F6CDD1D)
    }
    /// [0,1) の一様乱数
    pub fn f64(&mut self) -> f64 {
        (self.u64() >> 11) as f64 / (1u64 << 53) as f64
    }
    pub fn range(&mut self, n: usize) -> usize {
        (self.u64() % n as u64) as usize
    }
    /// 標準正規乱数 (Box–Muller)
    pub fn gauss(&mut self) -> f64 {
        let u1 = self.f64().max(1e-300);
        let u2 = self.f64();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    }
}

// ---------------- 複素数 ----------------
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct C64 {
    pub re: f64,
    pub im: f64,
}

pub const CZERO: C64 = C64 { re: 0.0, im: 0.0 };
pub const CONE: C64 = C64 { re: 1.0, im: 0.0 };
pub const CI: C64 = C64 { re: 0.0, im: 1.0 };

impl C64 {
    pub fn new(re: f64, im: f64) -> Self {
        C64 { re, im }
    }
    pub fn conj(self) -> Self {
        C64::new(self.re, -self.im)
    }
    pub fn norm2(self) -> f64 {
        self.re * self.re + self.im * self.im
    }
    pub fn abs(self) -> f64 {
        self.norm2().sqrt()
    }
    pub fn expi(theta: f64) -> Self {
        C64::new(theta.cos(), theta.sin())
    }
    pub fn scale(self, a: f64) -> Self {
        C64::new(self.re * a, self.im * a)
    }
}

impl Add for C64 {
    type Output = C64;
    fn add(self, o: C64) -> C64 {
        C64::new(self.re + o.re, self.im + o.im)
    }
}
impl Sub for C64 {
    type Output = C64;
    fn sub(self, o: C64) -> C64 {
        C64::new(self.re - o.re, self.im - o.im)
    }
}
impl Mul for C64 {
    type Output = C64;
    fn mul(self, o: C64) -> C64 {
        C64::new(
            self.re * o.re - self.im * o.im,
            self.re * o.im + self.im * o.re,
        )
    }
}
impl Div for C64 {
    type Output = C64;
    fn div(self, o: C64) -> C64 {
        let d = o.norm2();
        C64::new(
            (self.re * o.re + self.im * o.im) / d,
            (self.im * o.re - self.re * o.im) / d,
        )
    }
}
impl Neg for C64 {
    type Output = C64;
    fn neg(self) -> C64 {
        C64::new(-self.re, -self.im)
    }
}

/// 複素三重対角方程式 (Thomas 法): a[i]x[i-1] + b[i]x[i] + c[i]x[i+1] = d[i]
pub fn solve_tridiag_c(a: &[C64], b: &[C64], c: &[C64], d: &[C64]) -> Vec<C64> {
    let n = b.len();
    let mut cp = vec![CZERO; n];
    let mut dp = vec![CZERO; n];
    cp[0] = c[0] / b[0];
    dp[0] = d[0] / b[0];
    for i in 1..n {
        let m = b[i] - a[i] * cp[i - 1];
        cp[i] = c[i] / m;
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m;
    }
    let mut x = vec![CZERO; n];
    x[n - 1] = dp[n - 1];
    for i in (0..n - 1).rev() {
        x[i] = dp[i] - cp[i] * x[i + 1];
    }
    x
}

// ---------------- 実対称行列の固有値分解 (循環ヤコビ法) ----------------
/// a: n×n 実対称 (列優先 a[i + j*n])。戻り値: (固有値昇順, 固有ベクトル列優先 v[i + k*n] = k番目ベクトルの i 成分)
pub fn jacobi_eigh(a_in: &[f64], n: usize) -> (Vec<f64>, Vec<f64>) {
    let mut a = a_in.to_vec();
    let mut v = vec![0.0; n * n];
    for i in 0..n {
        v[i + i * n] = 1.0;
    }
    let norm0: f64 = a.iter().map(|x| x * x).sum::<f64>().max(1e-300);
    for _sweep in 0..200 {
        let mut off = 0.0;
        for p in 0..n {
            for q in (p + 1)..n {
                off += a[p + q * n] * a[p + q * n];
            }
        }
        if off < 1e-28 * norm0 {
            break;
        }
        for p in 0..n.saturating_sub(1) {
            for q in (p + 1)..n {
                let apq = a[p + q * n];
                if apq.abs() < 1e-300 {
                    continue;
                }
                let theta = (a[q + q * n] - a[p + p * n]) / (2.0 * apq);
                let t = theta.signum() / (theta.abs() + (theta * theta + 1.0).sqrt());
                let t = if theta == 0.0 { 1.0 } else { t };
                let c = 1.0 / (t * t + 1.0).sqrt();
                let s = t * c;
                // A <- G^T A G (G は (p,q) 面の回転)
                for k in 0..n {
                    let akp = a[k + p * n];
                    let akq = a[k + q * n];
                    a[k + p * n] = c * akp - s * akq;
                    a[k + q * n] = s * akp + c * akq;
                }
                for k in 0..n {
                    let apk = a[p + k * n];
                    let aqk = a[q + k * n];
                    a[p + k * n] = c * apk - s * aqk;
                    a[q + k * n] = s * apk + c * aqk;
                }
                for k in 0..n {
                    let vkp = v[k + p * n];
                    let vkq = v[k + q * n];
                    v[k + p * n] = c * vkp - s * vkq;
                    v[k + q * n] = s * vkp + c * vkq;
                }
            }
        }
    }
    let mut w: Vec<f64> = (0..n).map(|i| a[i + i * n]).collect();
    // 固有値で昇順ソート(固有ベクトル列も並べ替え)
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| w[i].partial_cmp(&w[j]).unwrap());
    let w_sorted: Vec<f64> = idx.iter().map(|&i| w[i]).collect();
    let mut v_sorted = vec![0.0; n * n];
    for (new_j, &old_j) in idx.iter().enumerate() {
        for i in 0..n {
            v_sorted[i + new_j * n] = v[i + old_j * n];
        }
    }
    w = w_sorted;
    (w, v_sorted)
}

/// 実対称行列関数 f(A) = V f(Λ) V^T
pub fn matfun_sym(a: &[f64], n: usize, f: impl Fn(f64) -> f64) -> Vec<f64> {
    let (w, v) = jacobi_eigh(a, n);
    let fw: Vec<f64> = w.iter().map(|&x| f(x)).collect();
    let mut out = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            let mut s = 0.0;
            for k in 0..n {
                s += v[i + k * n] * fw[k] * v[j + k * n];
            }
            out[i + j * n] = s;
        }
    }
    out
}

/// 行列積 (n×n, 列優先)
pub fn matmul(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    let mut c = vec![0.0; n * n];
    for j in 0..n {
        for k in 0..n {
            let bkj = b[k + j * n];
            if bkj == 0.0 {
                continue;
            }
            for i in 0..n {
                c[i + j * n] += a[i + k * n] * bkj;
            }
        }
    }
    c
}

// ---------------- 特殊関数 ----------------
/// 第1種変形ベッセル関数 I_nu(x) (整数次, 級数展開)
pub fn bessel_i(nu: u32, x: f64) -> f64 {
    let half = x / 2.0;
    let mut term = 1.0;
    for k in 1..=nu {
        term *= half / k as f64;
    }
    let mut sum = term;
    let mut k = 1.0;
    loop {
        term *= half * half / (k * (k + nu as f64));
        sum += term;
        if term < 1e-17 * sum || k > 500.0 {
            break;
        }
        k += 1.0;
    }
    sum
}

/// ln Γ(x) (Lanczos 近似, x > 0)
pub fn ln_gamma(x: f64) -> f64 {
    const G: [f64; 9] = [
        0.99999999999980993,
        676.5203681218851,
        -1259.1392167224028,
        771.32342877765313,
        -176.61502916214059,
        12.507343278686905,
        -0.13857109526572012,
        9.9843695780195716e-6,
        1.5056327351493116e-7,
    ];
    let x = x - 1.0;
    let mut a = G[0];
    let t = x + 7.5;
    for (i, &g) in G.iter().enumerate().skip(1) {
        a += g / (x + i as f64);
    }
    0.5 * (2.0 * std::f64::consts::PI).ln() + (x + 0.5) * t.ln() - t + a.ln()
}

// ---------------- 統計 ----------------
/// ビニング法による平均と標準誤差(自己相関を粗く考慮)
pub fn mean_err(xs: &[f64]) -> (f64, f64) {
    let n = xs.len();
    let mean = xs.iter().sum::<f64>() / n as f64;
    let nbin = 50.min(n);
    let bs = n / nbin;
    if bs == 0 {
        return (mean, 0.0);
    }
    let mut bm = Vec::with_capacity(nbin);
    for b in 0..nbin {
        let s: f64 = xs[b * bs..(b + 1) * bs].iter().sum();
        bm.push(s / bs as f64);
    }
    let var = bm.iter().map(|&x| (x - mean) * (x - mean)).sum::<f64>() / (nbin as f64 - 1.0);
    (mean, (var / nbin as f64).sqrt())
}

/// 最小二乗直線フィット y = a + b x。戻り値 (a, b)
pub fn linfit(x: &[f64], y: &[f64]) -> (f64, f64) {
    let n = x.len() as f64;
    let sx: f64 = x.iter().sum();
    let sy: f64 = y.iter().sum();
    let sxx: f64 = x.iter().map(|v| v * v).sum();
    let sxy: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    let b = (n * sxy - sx * sy) / (n * sxx - sx * sx);
    let a = (sy - b * sx) / n;
    (a, b)
}

// ---------------- 自己テスト ----------------
pub fn self_test() {
    // ヤコビ法: ランダム対称行列で A v = λ v を検証
    let n = 8;
    let mut rng = Rng::new(12345);
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..=i {
            let x = rng.gauss();
            a[i + j * n] = x;
            a[j + i * n] = x;
        }
    }
    let (w, v) = jacobi_eigh(&a, n);
    let mut max_res: f64 = 0.0;
    for k in 0..n {
        for i in 0..n {
            let mut av = 0.0;
            for j in 0..n {
                av += a[i + j * n] * v[j + k * n];
            }
            max_res = max_res.max((av - w[k] * v[i + k * n]).abs());
        }
    }
    assert!(max_res < 1e-9, "jacobi residual = {}", max_res);
    // ベッセル: I0(1)=1.2660658..., I1(2)=1.5906368...
    assert!((bessel_i(0, 1.0) - 1.2660658777520084).abs() < 1e-12);
    assert!((bessel_i(1, 2.0) - 1.5906368546373291).abs() < 1e-12);
    // ln Γ(5) = ln 24
    assert!((ln_gamma(5.0) - (24.0f64).ln()).abs() < 1e-10);
    eprintln!("[self_test] OK (jacobi residual {:.2e})", max_res);
}
