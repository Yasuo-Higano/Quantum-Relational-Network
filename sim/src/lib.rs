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

    /// 読み出し: 任意のサイト集合 A, B の相互情報量 I(A:B) = S(A)+S(B)−S(A∪B)
    pub fn readout_mi_regions(&self, ra: &[usize], rb: &[usize]) -> f64 {
        let sub = |sites: &[usize]| -> f64 {
            let k = sites.len();
            let n = self.n;
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
        let mut ab: Vec<usize> = ra.to_vec();
        ab.extend_from_slice(rb);
        sub(ra) + sub(rb) - sub(&ab)
    }

    /// 状態の健全性: 純粋なガウス状態なら C は射影子 (‖C²−C‖_max = 0)。
    /// v6.7 でバグを発見した検査 — 全ての core 模型がこれを通ること。
    pub fn purity_defect(&self) -> f64 {
        let n = self.n;
        let mut dmax: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                let (mut sre, mut sim) = (0.0, 0.0);
                for k in 0..n {
                    let (a, b) = (self.cre[i + k * n], self.cim[i + k * n]);
                    let (c, d) = (self.cre[k + j * n], self.cim[k + j * n]);
                    sre += a * c - b * d;
                    sim += a * d + b * c;
                }
                let dre = sre - self.cre[i + j * n];
                let dim = sim - self.cim[i + j * n];
                dmax = dmax.max((dre * dre + dim * dim).sqrt());
            }
        }
        dmax
    }

    /// 読み出し: ボンドエネルギー密度 h(x) = −2 Re C_{x,x+1} (生値; 基底状態値は引かない)
    pub fn readout_bond_energy(&self, x: usize) -> f64 {
        let n = self.n;
        -2.0 * self.cre[x % n + ((x + 1) % n) * n]
    }

    /// 読み出し: 次近接流 j(x) = 2 Im C_{x−1,x+1} — 半充填の 2 部格子対称性で密度が
    /// 凍結する系では、これが光円錐を運ぶ因果的観測量になる (v7.4 の教訓)
    pub fn readout_current(&self, x: usize) -> f64 {
        let n = self.n;
        2.0 * self.cim[(x + n - 1) % n + ((x + 1) % n) * n]
    }

    /// 読み出し: 光的エネルギー T₋₋ 型 (3 ボンド平均)。T_kk = 4T₋₋ 規格化
    /// (CLAUDE.md の落とし穴参照)。h_gs は基底状態のボンドエネルギー参照値。
    pub fn readout_null_energy(&self, x: usize, vf: f64, h_gs: f64) -> f64 {
        let (mut h, mut j) = (0.0, 0.0);
        for b in [x - 1, x, x + 1] {
            h += self.readout_bond_energy(b) - h_gs;
            j += self.readout_current(b);
        }
        h /= 3.0;
        j /= 3.0;
        2.0 * (h / vf + j / (vf * vf))
    }

    /// 読み出し: カイラル度 Σj/(v_F Σh) — 右向きの波束なら +1、定在波なら ≈0
    pub fn readout_chirality(&self, vf: f64, h_gs: f64) -> f64 {
        let n = self.n;
        let (mut sh, mut sj) = (0.0, 0.0);
        for x in 2..n - 2 {
            sh += self.readout_bond_energy(x) - h_gs;
            sj += self.readout_current(x);
        }
        sj / (vf * sh)
    }

    /// 読み出し: 自然なテンソル分解 — サイト間相互情報量 I(i:j) の貪欲最大マッチングで
    /// サイトを対に組む。「どのテンソル分解が自然か」は状態が決める (残高 7 の第一級化)。
    /// 戻り値: (対のリスト, 対に捕捉された Σ I の全ペア Σ I に対する割合)。
    pub fn readout_natural_partition(&self) -> (Vec<(usize, usize)>, f64) {
        let n = self.n;
        let s1: Vec<f64> = (0..n).map(|i| h2_entropy(self.cre[i + i * n])).collect();
        let mut mi = vec![0.0f64; n * n];
        let mut total = 0.0f64;
        for i in 0..n {
            for j in (i + 1)..n {
                // 2 サイトのエルミート相関行列のエントロピー
                let cre = [
                    self.cre[i + i * n],
                    self.cre[j + i * n],
                    self.cre[i + j * n],
                    self.cre[j + j * n],
                ];
                let cim = [
                    self.cim[i + i * n],
                    self.cim[j + i * n],
                    self.cim[i + j * n],
                    self.cim[j + j * n],
                ];
                let s2 = entropy_corr_herm(&cre, &cim, 2);
                let m = (s1[i] + s1[j] - s2).max(0.0);
                mi[i + j * n] = m;
                mi[j + i * n] = m;
                total += m;
            }
        }
        let mut used = vec![false; n];
        let mut pairs = Vec::with_capacity(n / 2);
        let mut captured = 0.0f64;
        for _ in 0..n / 2 {
            let (mut bi, mut bj, mut bm) = (0usize, 0usize, -1.0f64);
            for i in 0..n {
                if used[i] {
                    continue;
                }
                for j in (i + 1)..n {
                    if used[j] {
                        continue;
                    }
                    if mi[i + j * n] > bm {
                        bm = mi[i + j * n];
                        bi = i;
                        bj = j;
                    }
                }
            }
            used[bi] = true;
            used[bj] = true;
            captured += bm;
            pairs.push((bi, bj));
        }
        (pairs, if total > 0.0 { captured / total } else { 0.0 })
    }
}

/// 熱場二重 (TFD): 2 本の円環鎖 L/R をモードごとに縫った純粋状態。
/// C_LL(k)=f_k (フェルミ分布), C_RR(k)=1−f_k, C_LR(k)=√(f_k(1−f_k)) — 各 k の
/// 2×2 ブロックの固有値は {0,1} なので全体は厳密に純粋 (もつれで熱を装う)。
/// v1.2 (ER=EPR) の状態を core の語彙で再定義したもの。
pub struct TfdPair {
    pub n: usize, // 各鎖のサイト数
    pub beta: f64,
}

impl QrnModel for TfdPair {
    fn assumptions(&self) -> Vec<&'static str> {
        vec![
            "A0: 2 本の鎖 L/R のテンソル積 (接続性は公理に置かない)",
            "A1: ユニタリー発展 (H_L + H_R)",
            "初期状態: 熱場二重 (全系純粋、片側は厳密に熱的)",
        ]
    }
    fn claims(&self) -> Vec<&'static str> {
        vec!["QRN-CORE-002", "QRN-ER-001"]
    }
    fn init(&self) -> QrnState {
        let n = self.n;
        let two_pi = 2.0 * std::f64::consts::PI;
        // モード和で実空間ブロックを構成 (f_k = f_{-k} なので実対称)
        let mut ff = vec![0.0; n];
        let mut gg = vec![0.0; n];
        for d in 0..n {
            let (mut sf, mut sg) = (0.0, 0.0);
            for j in 0..n {
                let k = two_pi * j as f64 / n as f64;
                let e = -2.0 * k.cos();
                let f = 1.0 / (1.0 + (self.beta * e).exp());
                let c = (k * d as f64).cos();
                sf += f * c;
                sg += (f * (1.0 - f)).max(0.0).sqrt() * c;
            }
            ff[d] = sf / n as f64;
            gg[d] = sg / n as f64;
        }
        let m = 2 * n;
        let mut cre = vec![0.0; m * m];
        for x in 0..n {
            for y in 0..n {
                let d = (x + n - y) % n;
                cre[x + y * m] = ff[d]; // LL
                                        // RR = 1 − F (純粋性の条件)
                cre[(x + n) + (y + n) * m] = if d == 0 { 1.0 - ff[0] } else { -ff[d] };
                cre[x + (y + n) * m] = gg[d]; // LR
                cre[(x + n) + y * m] = gg[d];
            }
        }
        QrnState {
            n: m,
            cre,
            cim: vec![0.0; m * m],
        }
    }
    /// H_L + H_R の発展 (各鎖のリング発展をブロック対角に適用)
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
        let m = 2 * n;
        let mut ur = vec![0.0; m * m];
        let mut ui = vec![0.0; m * m];
        for x in 0..n {
            for y in 0..n {
                let d = (x + n - y) % n;
                ur[x + y * m] = ud[d].0;
                ui[x + y * m] = ud[d].1;
                ur[(x + n) + (y + n) * m] = ud[d].0;
                ui[(x + n) + (y + n) * m] = ud[d].1;
            }
        }
        let t1r = matmul(&ur, &s.cre, m);
        let t1r2 = matmul(&ui, &s.cim, m);
        let t1i = matmul(&ur, &s.cim, m);
        let t1i2 = matmul(&ui, &s.cre, m);
        let ar: Vec<f64> = t1r.iter().zip(&t1r2).map(|(a, b)| a - b).collect();
        let ai: Vec<f64> = t1i.iter().zip(&t1i2).map(|(a, b)| a + b).collect();
        let mut vr = vec![0.0; m * m];
        let mut vi = vec![0.0; m * m];
        for x in 0..m {
            for y in 0..m {
                vr[x + y * m] = ur[y + x * m];
                vi[x + y * m] = -ui[y + x * m];
            }
        }
        let b1 = matmul(&ar, &vr, m);
        let b2 = matmul(&ai, &vi, m);
        let b3 = matmul(&ar, &vi, m);
        let b4 = matmul(&ai, &vr, m);
        QrnState {
            n: m,
            cre: b1.iter().zip(&b2).map(|(a, b)| a - b).collect(),
            cim: b3.iter().zip(&b4).map(|(a, b)| a + b).collect(),
        }
    }
}

/// 成長する鎖 (v5.1 の機構を core の語彙で): 固定の n_max サイトのうち active サイト
/// だけがホッピングで結ばれ、新サイトは真空 (C=0) で到着する。
/// 成長はテンソル分解自体の変化なので QrnModel::evolve の外に専用メソッドを持つ
/// (「A0 のテンソル分解が動く」拡張は core の残高)。
pub struct GrowingChain {
    pub n_max: usize,
}

impl GrowingChain {
    /// 最初の n0 サイトの開鎖基底状態 (残りは真空)
    pub fn init(&self, n0: usize) -> QrnState {
        let n = self.n_max;
        let mut cre = vec![0.0; n * n];
        let c0 = Self::open_gs(n0);
        for i in 0..n0 {
            for j in 0..n0 {
                cre[i + j * n] = c0[i + j * n0];
            }
        }
        QrnState {
            n,
            cre,
            cim: vec![0.0; n * n],
        }
    }
    fn open_gs(n0: usize) -> Vec<f64> {
        let mut h = vec![0.0; n0 * n0];
        for x in 0..n0 - 1 {
            h[x + (x + 1) * n0] = -1.0;
            h[(x + 1) + x * n0] = -1.0;
        }
        let (_, v) = jacobi_eigh(&h, n0);
        let nocc = n0 / 2;
        let mut c = vec![0.0; n0 * n0];
        for m in 0..nocc {
            for i in 0..n0 {
                let vi = v[i + m * n0];
                for j in 0..n0 {
                    c[i + j * n0] += vi * v[j + m * n0];
                }
            }
        }
        c
    }
    /// active サイトの開鎖ハミルトニアンで時間 t 発展
    pub fn evolve_active(&self, s: &QrnState, active: usize, t: f64) -> QrnState {
        let n = self.n_max;
        let mut h = vec![0.0; active * active];
        for x in 0..active - 1 {
            h[x + (x + 1) * active] = -1.0;
            h[(x + 1) + x * active] = -1.0;
        }
        let (w, v) = jacobi_eigh(&h, active);
        let mut ur = vec![0.0; active * active];
        let mut ui = vec![0.0; active * active];
        for i in 0..active {
            for j in 0..active {
                let (mut sr, mut si) = (0.0, 0.0);
                for k in 0..active {
                    let ph = -w[k] * t;
                    sr += v[i + k * active] * ph.cos() * v[j + k * active];
                    si += v[i + k * active] * ph.sin() * v[j + k * active];
                }
                ur[i + j * active] = sr;
                ui[i + j * active] = si;
            }
        }
        let mut cre = vec![0.0; active * active];
        let mut cim = vec![0.0; active * active];
        for i in 0..active {
            for j in 0..active {
                cre[i + j * active] = s.cre[i + j * n];
                cim[i + j * active] = s.cim[i + j * n];
            }
        }
        let t1r = matmul(&ur, &cre, active);
        let t1r2 = matmul(&ui, &cim, active);
        let t1i = matmul(&ur, &cim, active);
        let t1i2 = matmul(&ui, &cre, active);
        let ar: Vec<f64> = t1r.iter().zip(&t1r2).map(|(a, b)| a - b).collect();
        let ai: Vec<f64> = t1i.iter().zip(&t1i2).map(|(a, b)| a + b).collect();
        let mut vr = vec![0.0; active * active];
        let mut vi = vec![0.0; active * active];
        for x in 0..active {
            for y in 0..active {
                vr[x + y * active] = ur[y + x * active];
                vi[x + y * active] = -ui[y + x * active];
            }
        }
        let b1 = matmul(&ar, &vr, active);
        let b2 = matmul(&ai, &vi, active);
        let b3 = matmul(&ar, &vi, active);
        let b4 = matmul(&ai, &vr, active);
        let mut out = QrnState {
            n,
            cre: s.cre.clone(),
            cim: s.cim.clone(),
        };
        for i in 0..active {
            for j in 0..active {
                out.cre[i + j * n] = b1[i + j * active] - b2[i + j * active];
                out.cim[i + j * n] = b3[i + j * active] + b4[i + j * active];
            }
        }
        out
    }
    /// 新 2 サイトを局所基底状態 (結合軌道に 1 粒子 = 純粋・半充填を保つ「真空」) で
    /// 到着させる: C_new = [[1/2,1/2],[1/2,1/2]] (射影子)。v5.1 と同一のプロトコル。
    pub fn arrive_pair_vacuum(&self, s: &QrnState, at: usize) -> QrnState {
        let n = self.n_max;
        let mut out = QrnState {
            n,
            cre: s.cre.clone(),
            cim: s.cim.clone(),
        };
        out.cre[at + at * n] = 0.5;
        out.cre[(at + 1) + (at + 1) * n] = 0.5;
        out.cre[at + (at + 1) * n] = 0.5;
        out.cre[(at + 1) + at * n] = 0.5;
        out
    }
    /// 対照シナリオ: 新 2 サイトを熱的 (最大混合 C=1/2·I) に到着させる
    pub fn arrive_pair_thermal(&self, s: &QrnState, at: usize) -> QrnState {
        let n = self.n_max;
        let mut out = QrnState {
            n,
            cre: s.cre.clone(),
            cim: s.cim.clone(),
        };
        out.cre[at + at * n] = 0.5;
        out.cre[(at + 1) + (at + 1) * n] = 0.5;
        out
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

impl RingChain {
    /// 基底状態のボンドエネルギー h_gs = −2 C⁰(1) (readout_null_energy の参照値)
    pub fn gs_bond_energy(&self) -> f64 {
        let n = self.n as f64;
        -2.0 / (n * (std::f64::consts::PI / n).sin())
    }
}

/// カイラル波束つき円環 — v4.1/v6.3 の QNEC 掃引状態を core の語彙で再定義した模型。
/// 円環 GS の占有窓 [j_F−20, j_F] と空窓 [j_F+1, j_F+22] をガウス重み (幅 sig) で
/// 束ねた波束対を角度 alpha で回転させた励起状態。standing=true は鏡像モード (k→−k)
/// を重ねた定在波 (非カイラル対照 — QNEC は成立し続けるが共動凍結が壊れる)。
pub struct PacketRing {
    pub ring: RingChain,
    pub xc: f64,    // 波束の実空間中心
    pub alpha: f64, // 回転角 (励起の強さ)
    pub sig: f64,   // モード窓のガウス幅
    pub standing: bool,
}

impl QrnModel for PacketRing {
    fn assumptions(&self) -> Vec<&'static str> {
        let mut a = self.ring.assumptions();
        a.push("励起: フェルミ面近傍のモード窓を波束に束ねた 1 粒子-1 空孔回転 (ガウス状態を保つ)");
        a
    }
    fn claims(&self) -> Vec<&'static str> {
        vec!["QRN-QNEC-001", "QRN-QNEC-002", "QRN-CORE-003"]
    }
    fn init(&self) -> QrnState {
        let n = self.ring.n;
        let two_pi = 2.0 * std::f64::consts::PI;
        let mut s0 = self.ring.init();
        let jq = n / 4;
        let sig = self.sig;
        let mut hre = vec![0.0; n];
        let mut him = vec![0.0; n];
        let mut pre = vec![0.0; n];
        let mut pim = vec![0.0; n];
        let (mut nh, mut np) = (0.0, 0.0);
        let xc = self.xc;
        let add_window =
            |lo: i64, hi: i64, ctr: f64, re: &mut Vec<f64>, im: &mut Vec<f64>, nrm: &mut f64| {
                for j in lo..=hi {
                    let wj = (-((j as f64 - ctr) * (j as f64 - ctr)) / (2.0 * sig * sig)).exp();
                    *nrm += wj * wj;
                    for x in 0..n {
                        let ph = two_pi * j as f64 * (x as f64 - xc) / n as f64;
                        re[x] += wj * ph.cos();
                        im[x] += wj * ph.sin();
                    }
                }
            };
        let (jh, jp) = (jq as i64 - 8, jq as i64 + 9);
        add_window(
            jq as i64 - 20,
            jq as i64,
            jh as f64,
            &mut hre,
            &mut him,
            &mut nh,
        );
        add_window(
            jq as i64 + 1,
            jq as i64 + 22,
            jp as f64,
            &mut pre,
            &mut pim,
            &mut np,
        );
        if self.standing {
            let ni = n as i64;
            add_window(
                ni - jq as i64,
                ni - jq as i64 + 20,
                (ni - jq as i64 + 8) as f64,
                &mut hre,
                &mut him,
                &mut nh,
            );
            add_window(
                ni - jq as i64 - 22,
                ni - jq as i64 - 1,
                (ni - jq as i64 - 9) as f64,
                &mut pre,
                &mut pim,
                &mut np,
            );
        }
        let (nh, np) = ((nh * n as f64).sqrt(), (np * n as f64).sqrt());
        for x in 0..n {
            hre[x] /= nh;
            him[x] /= nh;
            pre[x] /= np;
            pim[x] /= np;
        }
        let (s, c) = (self.alpha.sin(), self.alpha.cos());
        for x in 0..n {
            for y in 0..n {
                let hh_re = hre[x] * hre[y] + him[x] * him[y];
                let hh_im = him[x] * hre[y] - hre[x] * him[y];
                let pp_re = pre[x] * pre[y] + pim[x] * pim[y];
                let pp_im = pim[x] * pre[y] - pre[x] * pim[y];
                let hp_re = hre[x] * pre[y] + him[x] * pim[y];
                let hp_im = him[x] * pre[y] - hre[x] * pim[y];
                let ph_re = pre[x] * hre[y] + pim[x] * him[y];
                let ph_im = pim[x] * hre[y] - pre[x] * him[y];
                s0.cre[x + y * n] += -s * s * hh_re + s * s * pp_re + s * c * (hp_re + ph_re);
                s0.cim[x + y * n] += -s * s * hh_im + s * s * pp_im + s * c * (hp_im + ph_im);
            }
        }
        s0
    }
    fn evolve(&self, s: &QrnState, t: f64) -> QrnState {
        self.ring.evolve(s, t)
    }
}

// ---------------- 疎行列固有値 (Lanczos 法) ----------------
// v13.1 (残高 11): 稠密ヤコビ法は ~2600 次元 (実埋め込み) が実用上限だった。
// 大格子 (N=18 の T⁴ = 10 万サイト) の最低バンドには matvec だけで動く
// Lanczos 法が要る。ここでは複素エルミート作用素の最低 k 固有対を、
// **完全再直交化** (縮退帯で必須) と**残差検証**つきで求める。
// ベクトルは (re, im) 対の Vec — 各バイナリの C3v 表現と互換。

/// 複素ベクトルの内積 ⟨a|b⟩ = Σ conj(a) b
pub fn cdot(a: &[(f64, f64)], b: &[(f64, f64)]) -> (f64, f64) {
    let (mut re, mut im) = (0.0, 0.0);
    for i in 0..a.len() {
        let (ar, ai) = a[i];
        let (br, bi) = b[i];
        re += ar * br + ai * bi;
        im += ar * bi - ai * br;
    }
    (re, im)
}

/// Lanczos 法: matvec で与えられた複素エルミート作用素の最低 k 固有対。
/// 完全再直交化つき。戻り値: (固有値昇順 k 個, 固有ベクトル, 最大残差 ‖Hv−λv‖)。
/// m は Krylov 次元 (k の 10 倍以上を推奨)。決定論 (固定シード)。
pub fn lanczos_lowest_herm(
    matvec: &dyn Fn(&[(f64, f64)]) -> Vec<(f64, f64)>,
    n: usize,
    k: usize,
    m: usize,
    seed: u64,
) -> (Vec<f64>, Vec<Vec<(f64, f64)>>, f64) {
    let mut rng = Rng::new(seed);
    let m = m.min(n);
    // 初期ベクトル (正規化した乱数)
    let mut v: Vec<(f64, f64)> = (0..n).map(|_| (rng.gauss(), rng.gauss())).collect();
    let nrm = cdot(&v, &v).0.sqrt();
    for x in v.iter_mut() {
        x.0 /= nrm;
        x.1 /= nrm;
    }
    let mut basis: Vec<Vec<(f64, f64)>> = vec![v.clone()];
    let mut alpha: Vec<f64> = Vec::new();
    let mut beta: Vec<f64> = Vec::new();
    for j in 0..m {
        let mut w = matvec(&basis[j]);
        let a = cdot(&basis[j], &w).0; // エルミートなので実
        alpha.push(a);
        // w -= α v_j + β v_{j-1} は完全再直交化に吸収させる
        // 完全再直交化 (2 回 — "twice is enough")
        for _ in 0..2 {
            for b in &basis {
                let (pr, pi) = cdot(b, &w);
                for i in 0..n {
                    let (br, bi) = b[i];
                    w[i].0 -= pr * br - pi * bi;
                    w[i].1 -= pr * bi + pi * br;
                }
            }
        }
        let bnorm = cdot(&w, &w).0.sqrt();
        if j + 1 == m || bnorm < 1e-12 {
            if bnorm < 1e-12 && j + 1 < m {
                // 不変部分空間に到達 — 新しい乱数方向で再開
                let mut r: Vec<(f64, f64)> = (0..n).map(|_| (rng.gauss(), rng.gauss())).collect();
                for _ in 0..2 {
                    for b in &basis {
                        let (pr, pi) = cdot(b, &r);
                        for i in 0..n {
                            let (br, bi) = b[i];
                            r[i].0 -= pr * br - pi * bi;
                            r[i].1 -= pr * bi + pi * br;
                        }
                    }
                }
                let rn = cdot(&r, &r).0.sqrt();
                if rn < 1e-12 {
                    break; // 空間を張り切った
                }
                for x in r.iter_mut() {
                    x.0 /= rn;
                    x.1 /= rn;
                }
                beta.push(0.0);
                basis.push(r);
                continue;
            }
            break;
        }
        beta.push(bnorm);
        for x in w.iter_mut() {
            x.0 /= bnorm;
            x.1 /= bnorm;
        }
        basis.push(w);
    }
    let mm = alpha.len();
    // 三重対角 (再直交化・再開で正確には帯だが、Ritz には射影行列を使うのが安全):
    // T_ij = ⟨v_i|H|v_j⟩ を陽に作る (mm ≤ 数百なので matvec mm 回ぶんは既に払った —
    // ここでは α, β から三重対角を組む代わりに、基底での射影を再計算して頑健にする)
    let mut t = vec![0.0f64; mm * mm];
    for j in 0..mm {
        let w = matvec(&basis[j]);
        for i in 0..mm {
            t[i + j * mm] = cdot(&basis[i], &w).0;
        }
    }
    // 対称化 (数値ノイズ除去)
    for i in 0..mm {
        for j in (i + 1)..mm {
            let a = 0.5 * (t[i + j * mm] + t[j + i * mm]);
            t[i + j * mm] = a;
            t[j + i * mm] = a;
        }
    }
    let (w_t, v_t) = jacobi_eigh(&t, mm);
    let kk = k.min(mm);
    let mut evals = Vec::with_capacity(kk);
    let mut evecs = Vec::with_capacity(kk);
    let mut max_res = 0.0f64;
    for e in 0..kk {
        let mut vec: Vec<(f64, f64)> = vec![(0.0, 0.0); n];
        for j in 0..mm {
            let c = v_t[j + e * mm];
            for i in 0..n {
                vec[i].0 += c * basis[j][i].0;
                vec[i].1 += c * basis[j][i].1;
            }
        }
        let nrm = cdot(&vec, &vec).0.sqrt();
        for x in vec.iter_mut() {
            x.0 /= nrm;
            x.1 /= nrm;
        }
        // 残差 ‖Hv − λv‖
        let hv = matvec(&vec);
        let lam = w_t[e];
        let mut res = 0.0f64;
        for i in 0..n {
            let dr = hv[i].0 - lam * vec[i].0;
            let di = hv[i].1 - lam * vec[i].1;
            res += dr * dr + di * di;
        }
        max_res = max_res.max(res.sqrt());
        evals.push(lam);
        evecs.push(vec);
    }
    (evals, evecs, max_res)
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
