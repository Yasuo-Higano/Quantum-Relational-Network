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

// ---------------- SHA-256 (証明書用, FIPS 180-4) ----------------
/// SHA-256 ハッシュ。探索の解集合など「証明書」の同一性検証に使う。
pub fn sha256(data: &[u8]) -> [u8; 32] {
    const K: [u32; 64] = [
        0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5, 0x3956c25b, 0x59f111f1, 0x923f82a4,
        0xab1c5ed5, 0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3, 0x72be5d74, 0x80deb1fe,
        0x9bdc06a7, 0xc19bf174, 0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc, 0x2de92c6f,
        0x4a7484aa, 0x5cb0a9dc, 0x76f988da, 0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
        0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967, 0x27b70a85, 0x2e1b2138, 0x4d2c6dfc,
        0x53380d13, 0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85, 0xa2bfe8a1, 0xa81a664b,
        0xc24b8b70, 0xc76c51a3, 0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070, 0x19a4c116,
        0x1e376c08, 0x2748774c, 0x34b0bcb5, 0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
        0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208, 0x90befffa, 0xa4506ceb, 0xbef9a3f7,
        0xc67178f2,
    ];
    let mut h: [u32; 8] = [
        0x6a09e667, 0xbb67ae85, 0x3c6ef372, 0xa54ff53a, 0x510e527f, 0x9b05688c, 0x1f83d9ab,
        0x5be0cd19,
    ];
    let bitlen = (data.len() as u64).wrapping_mul(8);
    let mut msg = data.to_vec();
    msg.push(0x80);
    while msg.len() % 64 != 56 {
        msg.push(0);
    }
    msg.extend_from_slice(&bitlen.to_be_bytes());
    for chunk in msg.chunks(64) {
        let mut w = [0u32; 64];
        for i in 0..16 {
            w[i] = u32::from_be_bytes([
                chunk[4 * i],
                chunk[4 * i + 1],
                chunk[4 * i + 2],
                chunk[4 * i + 3],
            ]);
        }
        for i in 16..64 {
            let s0 = w[i - 15].rotate_right(7) ^ w[i - 15].rotate_right(18) ^ (w[i - 15] >> 3);
            let s1 = w[i - 2].rotate_right(17) ^ w[i - 2].rotate_right(19) ^ (w[i - 2] >> 10);
            w[i] = w[i - 16]
                .wrapping_add(s0)
                .wrapping_add(w[i - 7])
                .wrapping_add(s1);
        }
        let (mut a, mut b, mut c, mut d, mut e, mut f, mut g, mut hh) =
            (h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
        for i in 0..64 {
            let s1 = e.rotate_right(6) ^ e.rotate_right(11) ^ e.rotate_right(25);
            let ch = (e & f) ^ ((!e) & g);
            let t1 = hh
                .wrapping_add(s1)
                .wrapping_add(ch)
                .wrapping_add(K[i])
                .wrapping_add(w[i]);
            let s0 = a.rotate_right(2) ^ a.rotate_right(13) ^ a.rotate_right(22);
            let maj = (a & b) ^ (a & c) ^ (b & c);
            let t2 = s0.wrapping_add(maj);
            hh = g;
            g = f;
            f = e;
            e = d.wrapping_add(t1);
            d = c;
            c = b;
            b = a;
            a = t1.wrapping_add(t2);
        }
        h[0] = h[0].wrapping_add(a);
        h[1] = h[1].wrapping_add(b);
        h[2] = h[2].wrapping_add(c);
        h[3] = h[3].wrapping_add(d);
        h[4] = h[4].wrapping_add(e);
        h[5] = h[5].wrapping_add(f);
        h[6] = h[6].wrapping_add(g);
        h[7] = h[7].wrapping_add(hh);
    }
    let mut out = [0u8; 32];
    for i in 0..8 {
        out[4 * i..4 * i + 4].copy_from_slice(&h[i].to_be_bytes());
    }
    out
}

pub fn sha256_hex(data: &[u8]) -> String {
    sha256(data).iter().map(|b| format!("{:02x}", b)).collect()
}

// ---------------- JSON (機械可読な結果成果物用) ----------------
/// 最小の JSON 値。外部クレート禁止の規約下で、結果成果物 (results/*.json,
/// certificates/*.json) を機械可読にするための自作実装。
pub enum Json {
    Null,
    Bool(bool),
    Int(i64),
    Num(f64),
    Str(String),
    Arr(Vec<Json>),
    Obj(Vec<(String, Json)>),
}

pub fn json_escape(s: &str) -> String {
    let mut out = String::with_capacity(s.len() + 2);
    for c in s.chars() {
        match c {
            '"' => out.push_str("\\\""),
            '\\' => out.push_str("\\\\"),
            '\n' => out.push_str("\\n"),
            '\t' => out.push_str("\\t"),
            '\r' => out.push_str("\\r"),
            c if (c as u32) < 0x20 => out.push_str(&format!("\\u{:04x}", c as u32)),
            c => out.push(c),
        }
    }
    out
}

impl Json {
    pub fn render(&self) -> String {
        let mut s = String::new();
        self.write(&mut s, 0);
        s
    }
    fn write(&self, out: &mut String, indent: usize) {
        let pad = "  ".repeat(indent);
        match self {
            Json::Null => out.push_str("null"),
            Json::Bool(b) => out.push_str(if *b { "true" } else { "false" }),
            Json::Int(i) => out.push_str(&i.to_string()),
            Json::Num(x) => {
                if x.is_finite() {
                    out.push_str(&format!("{}", x));
                } else {
                    out.push_str("null");
                }
            }
            Json::Str(s) => {
                out.push('"');
                out.push_str(&json_escape(s));
                out.push('"');
            }
            Json::Arr(v) => {
                if v.is_empty() {
                    out.push_str("[]");
                    return;
                }
                out.push_str("[\n");
                for (i, e) in v.iter().enumerate() {
                    out.push_str(&pad);
                    out.push_str("  ");
                    e.write(out, indent + 1);
                    if i + 1 < v.len() {
                        out.push(',');
                    }
                    out.push('\n');
                }
                out.push_str(&pad);
                out.push(']');
            }
            Json::Obj(kv) => {
                if kv.is_empty() {
                    out.push_str("{}");
                    return;
                }
                out.push_str("{\n");
                for (i, (k, v)) in kv.iter().enumerate() {
                    out.push_str(&pad);
                    out.push_str("  \"");
                    out.push_str(&json_escape(k));
                    out.push_str("\": ");
                    v.write(out, indent + 1);
                    if i + 1 < kv.len() {
                        out.push(',');
                    }
                    out.push('\n');
                }
                out.push_str(&pad);
                out.push('}');
            }
        }
    }
}

/// リポジトリルート相対パスへ成果物を書き出す。バイナリは sim/ から実行される
/// 想定 (規約) だが、ルートから実行しても機能するようにルートを検出する。
/// 戻り値は実際に書いたパス。
pub fn write_artifact(rel: &str, content: &str) -> String {
    let root = if std::path::Path::new("../results").is_dir() {
        ".."
    } else {
        "."
    };
    let path = format!("{}/{}", root, rel);
    if let Some(dir) = std::path::Path::new(&path).parent() {
        let _ = std::fs::create_dir_all(dir);
    }
    std::fs::write(&path, content)
        .unwrap_or_else(|e| panic!("成果物 {} の書き出し失敗: {}", path, e));
    path
}

// ---------------- QRN core (v6.7) ----------------
// 統一理論としての説得力は「新しいシミュレーションを増やすこと」ではなく
// 「既存のシミュレーションを同じ core から出すこと」で上がる (改良方針 §7)。
// 共通の状態空間 (ガウスフェルミオン網 = 相関行列) と共通の読み出しをここに置き、
// 各 vXY バイナリは QrnModel の実装 + 読み出しの組合せとして書けるようにする。

/// QRN の共通状態: ガウスフェルミオン状態。全ての物理量が相関行列
/// C_ij = ⟨c†_i c_j⟩ (エルミート) から厳密に計算できる。
pub struct QrnState {
    pub n: usize,
    pub cre: Vec<f64>,
    pub cim: Vec<f64>,
}

/// QRN 模型: 状態空間の初期化と動力学。仮定と主張 (claims.yml の id) を明示する。
pub trait QrnModel {
    fn assumptions(&self) -> Vec<&'static str>;
    fn claims(&self) -> Vec<&'static str>;
    fn init(&self) -> QrnState;
    fn evolve(&self, s: &QrnState, t: f64) -> QrnState;
}

/// 二値エントロピー h(z) = −z ln z − (1−z) ln(1−z)
pub fn h2_entropy(z: f64) -> f64 {
    let z = z.clamp(1e-14, 1.0 - 1e-14);
    -z * z.ln() - (1.0 - z) * (1.0 - z).ln()
}

/// 実対称相関行列のエンタングルメントエントロピー
pub fn entropy_corr_real(c: &[f64], n: usize) -> f64 {
    let (w, _) = jacobi_eigh(c, n);
    w.iter().map(|&z| h2_entropy(z)).sum()
}

/// エルミート相関行列のエントロピー (実埋め込み: 固有値は 2 重に出る)
pub fn entropy_corr_herm(cre: &[f64], cim: &[f64], n: usize) -> f64 {
    let m = 2 * n;
    let mut a = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            a[i + j * m] = cre[i + j * n];
            a[i + (j + n) * m] = -cim[i + j * n];
            a[(i + n) + j * m] = cim[i + j * n];
            a[(i + n) + (j + n) * m] = cre[i + j * n];
        }
    }
    let (w, _) = jacobi_eigh(&a, m);
    0.5 * w.iter().map(|&z| h2_entropy(z)).sum::<f64>()
}

impl QrnState {
    /// 読み出し: 区間 [a, a+len) のエンタングルメントエントロピー
    pub fn readout_entropy(&self, a: usize, len: usize) -> f64 {
        let n = self.n;
        let mut cre = vec![0.0; len * len];
        let mut cim = vec![0.0; len * len];
        let mut has_im = false;
        for i in 0..len {
            for j in 0..len {
                cre[i + j * len] = self.cre[(a + i) % n + ((a + j) % n) * n];
                cim[i + j * len] = self.cim[(a + i) % n + ((a + j) % n) * n];
                if cim[i + j * len] != 0.0 {
                    has_im = true;
                }
            }
        }
        if has_im {
            entropy_corr_herm(&cre, &cim, len)
        } else {
            entropy_corr_real(&cre, len)
        }
    }

    /// 読み出し: 2 サイトブロック間の相互情報量行列 (nb×nb) と最大値
    pub fn readout_mi_blocks(&self) -> (Vec<f64>, usize, f64) {
        let n = self.n;
        let nb = n / 2;
        let sub = |sites: &[usize]| -> f64 {
            let k = sites.len();
            let mut cre = vec![0.0; k * k];
            let mut cim = vec![0.0; k * k];
            for (a, &sa) in sites.iter().enumerate() {
                for (b, &sb) in sites.iter().enumerate() {
                    cre[a + b * k] = self.cre[sa + sb * n];
                    cim[a + b * k] = self.cim[sa + sb * n];
                }
            }
            entropy_corr_herm(&cre, &cim, k)
        };
        let sblk: Vec<f64> = (0..nb).map(|b| sub(&[2 * b, 2 * b + 1])).collect();
        let mut mi = vec![0.0; nb * nb];
        let mut mi_max = 0.0f64;
        for i in 0..nb {
            for j in (i + 1)..nb {
                let m = (sblk[i] + sblk[j] - sub(&[2 * i, 2 * i + 1, 2 * j, 2 * j + 1])).max(0.0);
                mi[i + j * nb] = m;
                mi[j + i * nb] = m;
                mi_max = mi_max.max(m);
            }
        }
        (mi, nb, mi_max)
    }

    /// 読み出し: 密度 → フェルミ運動量 (ラッティンジャーの定理の顔 k_F = π⟨n⟩)
    pub fn readout_fermi_momentum(&self) -> f64 {
        let n = self.n;
        let dens: f64 = (0..n).map(|i| self.cre[i + i * n]).sum::<f64>() / n as f64;
        std::f64::consts::PI * dens
    }
}

/// 幾何読み出しの結果 (v6.4 の円環判定と同じ規準)
pub struct RingMetrics {
    pub adjacency: f64,
    pub rsd: f64,
    pub lam21: f64,
    pub mi_max: f64,
    pub ring: bool,
}

/// 読み出し: MI 行列 → 情報距離 −ln(MI/MI_max) → 古典的 MDS → 円環判定
pub fn readout_ring_geometry(mi: &[f64], nb: usize, mi_max: f64) -> RingMetrics {
    if mi_max < 1e-12 {
        return RingMetrics {
            adjacency: 0.0,
            rsd: f64::INFINITY,
            lam21: 0.0,
            mi_max,
            ring: false,
        };
    }
    let mut d2 = vec![0.0; nb * nb];
    for i in 0..nb {
        for j in 0..nb {
            if i != j {
                let m = (mi[i + j * nb] / mi_max).max(1e-300);
                let dd = -m.ln();
                d2[i + j * nb] = dd * dd;
            }
        }
    }
    let row_mean: Vec<f64> = (0..nb)
        .map(|i| (0..nb).map(|j| d2[i + j * nb]).sum::<f64>() / nb as f64)
        .collect();
    let tot: f64 = row_mean.iter().sum::<f64>() / nb as f64;
    let mut b = vec![0.0; nb * nb];
    for i in 0..nb {
        for j in 0..nb {
            b[i + j * nb] = -0.5 * (d2[i + j * nb] - row_mean[i] - row_mean[j] + tot);
        }
    }
    let (w, v) = jacobi_eigh(&b, nb);
    let (l1, l2) = (w[nb - 1], w[nb - 2]);
    let coords: Vec<(f64, f64)> = (0..nb)
        .map(|i| {
            (
                l1.max(0.0).sqrt() * v[i + (nb - 1) * nb],
                l2.max(0.0).sqrt() * v[i + (nb - 2) * nb],
            )
        })
        .collect();
    let mut order: Vec<usize> = (0..nb).collect();
    order.sort_by(|&a, &bq| {
        coords[a]
            .1
            .atan2(coords[a].0)
            .partial_cmp(&coords[bq].1.atan2(coords[bq].0))
            .unwrap()
    });
    let mut adjacent_ok = 0;
    for k in 0..nb {
        let a = order[k];
        let bq = order[(k + 1) % nb];
        let d = (a as isize - bq as isize).unsigned_abs();
        if d == 1 || d == nb - 1 {
            adjacent_ok += 1;
        }
    }
    let radii: Vec<f64> = coords
        .iter()
        .map(|&(x, y)| (x * x + y * y).sqrt())
        .collect();
    let rmean: f64 = radii.iter().sum::<f64>() / nb as f64;
    let rsd = (radii.iter().map(|r| (r - rmean).powi(2)).sum::<f64>() / nb as f64).sqrt()
        / rmean.max(1e-300);
    let adjacency = adjacent_ok as f64 / nb as f64;
    let lam21 = if l1 > 0.0 { l2 / l1 } else { 0.0 };
    let ring = adjacency >= 0.9 && rsd <= 0.10 && lam21 >= 0.9;
    RingMetrics {
        adjacency,
        rsd,
        lam21,
        mi_max,
        ring,
    }
}

/// 具体的な QrnModel: 円環自由フェルミオン鎖 (半充填基底状態 + 最近接ホッピング)。
/// v0.5/v0.7/v1.1/v4.1 などの土台になっている系を core として一箇所に定義する。
pub struct RingChain {
    pub n: usize, // N ≡ 2 (mod 4) で閉殻・実相関
}

impl QrnModel for RingChain {
    fn assumptions(&self) -> Vec<&'static str> {
        vec![
            "A0: 有限次元ヒルベルト空間とテンソル分解 (サイト)",
            "A1: ユニタリー発展 (二次ハミルトニアン H = -Σ c†_x c_{x+1} + h.c.)",
            "初期状態: 半充填の基底状態 (ガウス状態)",
            "幾何・因果・物質は公理に置かず、読み出しで創発させる",
        ]
    }
    fn claims(&self) -> Vec<&'static str> {
        vec![
            "QRN-CORE-001",
            "QRN-GEOM-003",
            "QRN-ENT-001",
            "QRN-CAUSAL-001",
        ]
    }
    fn init(&self) -> QrnState {
        let n = self.n;
        let nocc = n / 2;
        let c0 = |d: isize| -> f64 {
            let d = d.unsigned_abs();
            let d = d.min(n - d);
            if d == 0 {
                return nocc as f64 / n as f64;
            }
            (std::f64::consts::PI * d as f64 / 2.0).sin()
                / (n as f64 * (std::f64::consts::PI * d as f64 / n as f64).sin())
        };
        let mut cre = vec![0.0; n * n];
        for x in 0..n {
            for y in 0..n {
                cre[x + y * n] = c0(x as isize - y as isize);
            }
        }
        QrnState {
            n,
            cre,
            cim: vec![0.0; n * n],
        }
    }
    /// U(t) C U† を厳密に計算 (並進不変な自由発展)
    fn evolve(&self, s: &QrnState, t: f64) -> QrnState {
        let n = self.n;
        let two_pi = 2.0 * std::f64::consts::PI;
        let mut ud = vec![(0.0f64, 0.0f64); n];
        for (d, slot) in ud.iter_mut().enumerate() {
            let (mut a, mut b) = (0.0, 0.0);
            for j in 0..n {
                let k = two_pi * j as f64 / n as f64;
                let e = -2.0 * k.cos();
                let ph = k * d as f64 - e * t;
                a += ph.cos();
                b += ph.sin();
            }
            *slot = (a / n as f64, b / n as f64);
        }
        let mut ur = vec![0.0; n * n];
        let mut ui = vec![0.0; n * n];
        for x in 0..n {
            for y in 0..n {
                let d = (x + n - y) % n;
                ur[x + y * n] = ud[d].0;
                ui[x + y * n] = ud[d].1;
            }
        }
        let t1r = matmul(&ur, &s.cre, n);
        let t1r2 = matmul(&ui, &s.cim, n);
        let t1i = matmul(&ur, &s.cim, n);
        let t1i2 = matmul(&ui, &s.cre, n);
        let ar: Vec<f64> = t1r.iter().zip(&t1r2).map(|(a, b)| a - b).collect();
        let ai: Vec<f64> = t1i.iter().zip(&t1i2).map(|(a, b)| a + b).collect();
        let mut vr = vec![0.0; n * n];
        let mut vi = vec![0.0; n * n];
        for x in 0..n {
            for y in 0..n {
                vr[x + y * n] = ur[y + x * n];
                vi[x + y * n] = -ui[y + x * n];
            }
        }
        let b1 = matmul(&ar, &vr, n);
        let b2 = matmul(&ai, &vi, n);
        let b3 = matmul(&ar, &vi, n);
        let b4 = matmul(&ai, &vr, n);
        QrnState {
            n,
            cre: b1.iter().zip(&b2).map(|(a, b)| a - b).collect(),
            cim: b3.iter().zip(&b4).map(|(a, b)| a + b).collect(),
        }
    }
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
    // SHA-256: FIPS 180-4 の検定ベクトル
    assert_eq!(
        sha256_hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256_hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    eprintln!("[self_test] OK (jacobi residual {:.2e})", max_res);
}
