//! double-double (DD, 有効 ~31 桁) 演算の共有モジュール (第二十五期, v24.1 で昇格)
//!
//! v22.7 (κ 床突破) と v23.4 (λ(N) 系列) で実証済みの DD 経路 — Dekker/Knuth の
//! error-free transformation (two_sum / two_prod は FMA) — をバイナリ内包から
//! lib へ昇格したもの。追加要素:
//!   - 演算子オーバーロード (+, −, ×, /) と比較
//!   - 有理角三角関数 sinpi_frac/cospi_frac: sin(π·m/den) を整数 mod で象限縮約
//!     してから Taylor — 大角の 2π 縮約誤差を避け、開放鎖 sine モードの構成を
//!     DD 全精度で行う (v24.2 のブロック経路の基礎)
//!   - Real trait: f64 と Dd を同一コードで扱う (κ クランプ・収束閾値も型ごと)
//!   - 汎用 cyclic Jacobi (jacobi_real) — lib::jacobi_eigh / v234 jacobi_dd と
//!     同じ規約 (固有値昇順, 固有ベクトル列優先 v[i + k*n])
//!
//! 精度の要点: DD の相対分解能 ~1e-32。相関行列 C の固有値は c ≳ 1e-30 まで解像
//! でき、モジュラー核 κ = ln((1−c)/c) の床は κ_max = ln(1e30) ≈ 69.1
//! (f64 経路の κ_max ≈ 32.2 の約 2 倍)。κ 床が K_A ボンド観測量に及ぼす影響の
//! 較正 (ξ 信頼域) は v24.2 が行う。

#[derive(Clone, Copy, Debug)]
pub struct Dd {
    pub hi: f64,
    pub lo: f64,
}

pub const DD0: Dd = Dd { hi: 0.0, lo: 0.0 };
pub const DD1: Dd = Dd { hi: 1.0, lo: 0.0 };
/// π の DD 表現 (hi = f64 π, lo = π − hi)
pub const DD_PI: Dd = Dd {
    hi: std::f64::consts::PI,
    lo: 1.2246467991473532e-16,
};

#[inline]
fn qts(a: f64, b: f64) -> Dd {
    // quick_two_sum: |a| ≥ |b| 前提
    let s = a + b;
    Dd {
        hi: s,
        lo: b - (s - a),
    }
}
#[inline]
fn two_sum(a: f64, b: f64) -> Dd {
    let s = a + b;
    let bb = s - a;
    Dd {
        hi: s,
        lo: (a - (s - bb)) + (b - bb),
    }
}
#[inline]
fn two_prod(a: f64, b: f64) -> Dd {
    let p = a * b;
    Dd {
        hi: p,
        lo: a.mul_add(b, -p),
    }
}

#[inline]
pub fn dd(a: f64) -> Dd {
    Dd { hi: a, lo: 0.0 }
}

impl Dd {
    #[inline]
    pub fn abs(self) -> Dd {
        if self.hi < 0.0 || (self.hi == 0.0 && self.lo < 0.0) {
            -self
        } else {
            self
        }
    }
    pub fn sqrt(self) -> Dd {
        if self.hi <= 0.0 {
            return DD0;
        }
        let x0 = self.hi.sqrt();
        // 1 回の DD Newton: x = (x + a/x)/2
        let x = dd(x0);
        (x + self / x) * dd(0.5)
    }
    /// DD exp: 範囲縮約 z = k·ln2 + f, Taylor 26 項 (v23.4 実証)
    pub fn exp(self) -> Dd {
        let ln2 = Dd {
            hi: 0.6931471805599453,
            lo: 2.3190468138462996e-17,
        };
        if self.hi < -745.0 {
            return DD0;
        }
        if self.hi > 709.0 {
            return dd(f64::INFINITY);
        }
        let k = (self.hi / ln2.hi).round();
        let f = self - ln2.mul_f(k);
        let mut term = DD1;
        let mut sum = DD1;
        for i in 1..27 {
            // 1/i の f64 事前丸め (v234 の実装) は項ごとに ~1e-17 の相対誤差を
            // 持ち込み DD 全精度を壊す — 正確な整数 DD 除算で置換 (v24.1)
            term = (term * f) / dd(i as f64);
            sum = sum + term;
        }
        // × 2^k (指数スケーリングは正確)
        let scale = (2.0f64).powi(k as i32);
        Dd {
            hi: sum.hi * scale,
            lo: sum.lo * scale,
        }
    }
    /// DD ln: y₀ = f64 ln + DD Newton 1 回 (y ← y + a·e^{−y} − 1) (v23.4 実証)
    pub fn ln(self) -> Dd {
        if self.hi <= 0.0 {
            return dd(-f64::INFINITY);
        }
        let y0 = self.hi.ln();
        let corr = self * dd(-y0).exp() - DD1;
        dd(y0) + corr
    }
    #[inline]
    pub fn mul_f(self, y: f64) -> Dd {
        let p = two_prod(self.hi, y);
        qts(p.hi, p.lo + self.lo * y)
    }
}

impl std::ops::Add for Dd {
    type Output = Dd;
    #[inline]
    fn add(self, y: Dd) -> Dd {
        let s = two_sum(self.hi, y.hi);
        qts(s.hi, s.lo + self.lo + y.lo)
    }
}
impl std::ops::Neg for Dd {
    type Output = Dd;
    #[inline]
    fn neg(self) -> Dd {
        Dd {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
}
impl std::ops::Sub for Dd {
    type Output = Dd;
    #[inline]
    fn sub(self, y: Dd) -> Dd {
        self + (-y)
    }
}
impl std::ops::Mul for Dd {
    type Output = Dd;
    #[inline]
    fn mul(self, y: Dd) -> Dd {
        let p = two_prod(self.hi, y.hi);
        qts(p.hi, p.lo + self.hi * y.lo + self.lo * y.hi)
    }
}
impl std::ops::Div for Dd {
    type Output = Dd;
    #[inline]
    fn div(self, y: Dd) -> Dd {
        let q1 = self.hi / y.hi;
        let r1 = self - y.mul_f(q1);
        let q2 = r1.hi / y.hi;
        let r2 = r1 - y.mul_f(q2);
        let q3 = r2.hi / y.hi;
        qts(q1, q2) + dd(q3)
    }
}

/// sin(π·x) の Taylor (|x| ≤ 1/4 想定, 15 項で < 1e-33)
fn sin_pi_taylor(x: Dd) -> Dd {
    let z = x * DD_PI;
    let z2 = z * z;
    let mut term = z;
    let mut sum = z;
    for k in 1..16 {
        let d = (2 * k) * (2 * k + 1);
        term = -((term * z2) / dd(d as f64));
        sum = sum + term;
    }
    sum
}
/// cos(π·x) の Taylor (|x| ≤ 1/4 想定)
fn cos_pi_taylor(x: Dd) -> Dd {
    let z = x * DD_PI;
    let z2 = z * z;
    let mut term = DD1;
    let mut sum = DD1;
    for k in 1..17 {
        let d = (2 * k - 1) * (2 * k);
        term = -((term * z2) / dd(d as f64));
        sum = sum + term;
    }
    sum
}

/// sin(π·num/den) を DD 全精度で。整数 mod で象限縮約するので num が大きくても
/// 2π 縮約誤差が出ない (開放鎖 sine モード φ_n(x) の構成用)。den > 0。
pub fn dd_sinpi_frac(num: i64, den: i64) -> Dd {
    debug_assert!(den > 0);
    // m/den ∈ [0, 2)
    let mut m = num.rem_euclid(2 * den);
    let mut sign = 1.0;
    if m >= den {
        // sin(π(1+t)) = −sin(πt)
        m -= den;
        sign = -sign;
    }
    if 2 * m > den {
        // sin(π(1−t)) = sin(πt)
        m = den - m;
    }
    // ここで m/den ∈ [0, 1/2]
    let v = if 4 * m <= den {
        sin_pi_taylor(dd(m as f64) / dd(den as f64))
    } else {
        // sin(πm/den) = cos(π(den−2m)/(2den)), (den−2m)/(2den) ∈ [0, 1/4)
        cos_pi_taylor(dd((den - 2 * m) as f64) / dd((2 * den) as f64))
    };
    if sign < 0.0 {
        -v
    } else {
        v
    }
}
/// cos(π·num/den) = sin(π·(num + den/2)/den) — 整数で書けるよう 2 倍角規約で処理
pub fn dd_cospi_frac(num: i64, den: i64) -> Dd {
    // cos(πx) = sin(π(x + 1/2)) → num/den + 1/2 = (2num + den)/(2den)
    dd_sinpi_frac(2 * num + den, 2 * den)
}

// ---------------- Real trait: f64 と Dd の統一 ----------------

/// f64 / Dd を同一コードで扱う実スカラー。閾値 (Jacobi 収束・κ クランプ) も型ごと。
pub trait Real:
    Copy
    + std::ops::Add<Output = Self>
    + std::ops::Sub<Output = Self>
    + std::ops::Mul<Output = Self>
    + std::ops::Neg<Output = Self>
{
    const R0: Self;
    const R1: Self;
    fn from_f64(x: f64) -> Self;
    /// 近似値 (判定・出力用)
    fn hi(self) -> f64;
    fn divr(self, o: Self) -> Self;
    fn sqrtr(self) -> Self;
    fn lnr(self) -> Self;
    fn sinpi_frac(num: i64, den: i64) -> Self;
    fn cospi_frac(num: i64, den: i64) -> Self;
    /// Jacobi 収束閾値 (非対角の絶対値)
    fn jac_tol() -> f64;
    /// 相関固有値の物理クランプ床 (κ 床 = −ln(clamp))
    fn c_clamp() -> f64;
}

impl Real for f64 {
    const R0: f64 = 0.0;
    const R1: f64 = 1.0;
    #[inline]
    fn from_f64(x: f64) -> f64 {
        x
    }
    #[inline]
    fn hi(self) -> f64 {
        self
    }
    #[inline]
    fn divr(self, o: f64) -> f64 {
        self / o
    }
    #[inline]
    fn sqrtr(self) -> f64 {
        self.sqrt()
    }
    #[inline]
    fn lnr(self) -> f64 {
        self.ln()
    }
    fn sinpi_frac(num: i64, den: i64) -> f64 {
        dd_sinpi_frac(num, den).hi
    }
    fn cospi_frac(num: i64, den: i64) -> f64 {
        dd_cospi_frac(num, den).hi
    }
    fn jac_tol() -> f64 {
        1e-15
    }
    fn c_clamp() -> f64 {
        1e-14
    }
}

impl Real for Dd {
    const R0: Dd = DD0;
    const R1: Dd = DD1;
    #[inline]
    fn from_f64(x: f64) -> Dd {
        dd(x)
    }
    #[inline]
    fn hi(self) -> f64 {
        self.hi
    }
    #[inline]
    fn divr(self, o: Dd) -> Dd {
        self / o
    }
    #[inline]
    fn sqrtr(self) -> Dd {
        self.sqrt()
    }
    #[inline]
    fn lnr(self) -> Dd {
        self.ln()
    }
    fn sinpi_frac(num: i64, den: i64) -> Dd {
        dd_sinpi_frac(num, den)
    }
    fn cospi_frac(num: i64, den: i64) -> Dd {
        dd_cospi_frac(num, den)
    }
    fn jac_tol() -> f64 {
        1e-29
    }
    fn c_clamp() -> f64 {
        1e-30
    }
}

// ---------------- 汎用 cyclic Jacobi (実対称) ----------------

/// a: n×n 実対称 (列優先)。戻り値: (固有値昇順, 固有ベクトル列優先 v[i + k*n])。
/// jacobi_eigh (f64) / v234 jacobi_dd と同じ規約の統一実装。
/// 収束判定は Demmel–Veselić の相対基準 |a_pq| ≤ eps·√(|a_pp|·|a_qq|) —
/// 深部固有値 (c ~ clamp 床) を相対精度で解像するために必須 (絶対閾値だと
/// 小さい対角ブロックが回転されず κ = ln((1−c)/c) が壊れる)。
pub fn jacobi_real<T: Real>(a_in: &[T], n: usize, sweeps_max: usize) -> (Vec<T>, Vec<T>) {
    let mut a = a_in.to_vec();
    let mut v = vec![T::R0; n * n];
    for i in 0..n {
        v[i + i * n] = T::R1;
    }
    let eps = T::jac_tol();
    for _sw in 0..sweeps_max {
        let mut conv = true;
        'outer: for p in 0..n {
            for q in p + 1..n {
                let den = (a[p + p * n].hi().abs() * a[q + q * n].hi().abs())
                    .sqrt()
                    .max(1e-300);
                if a[p + q * n].hi().abs() > eps * den {
                    conv = false;
                    break 'outer;
                }
            }
        }
        if conv {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                let apq = a[p + q * n];
                let den = (a[p + p * n].hi().abs() * a[q + q * n].hi().abs())
                    .sqrt()
                    .max(1e-300);
                if apq.hi().abs() <= eps * den {
                    continue;
                }
                let app = a[p + p * n];
                let aqq = a[q + q * n];
                let theta = (aqq - app).divr(apq * T::from_f64(2.0));
                let t = {
                    let tt = T::R1.divr(
                        T::from_f64(theta.hi().abs())
                            + (theta * theta + T::R1).sqrtr(),
                    );
                    if theta.hi() < 0.0 {
                        -tt
                    } else {
                        tt
                    }
                };
                let c = T::R1.divr((t * t + T::R1).sqrtr());
                let s = t * c;
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
    let evs: Vec<T> = (0..n).map(|i| a[i + i * n]).collect();
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| evs[i].hi().partial_cmp(&evs[j].hi()).unwrap());
    let evs_s: Vec<T> = idx.iter().map(|&i| evs[i]).collect();
    let mut v_s = vec![T::R0; n * n];
    for (k, &i) in idx.iter().enumerate() {
        for r in 0..n {
            v_s[r + k * n] = v[r + i * n];
        }
    }
    (evs_s, v_s)
}

/// 相関行列 C (実対称, 固有値 ∈ [0,1]) からモジュラー核 K = ln((1−C)/C) を構成。
/// clamp は固有値の物理床 (型の既定は c_clamp(): f64 1e-14 → κ≤32.2 / Dd 1e-30 → κ≤69.1
/// — 床較正 (v24.2) では別値を渡して感度を測る)。
/// 戻り値: (K 行列 [列優先], 固有値ベクトル c 昇順)
pub fn modular_k<T: Real>(c_mat: &[T], n: usize, sweeps: usize, clamp: f64) -> (Vec<T>, Vec<T>) {
    let (cw, cv) = jacobi_real(c_mat, n, sweeps);
    let kappa: Vec<T> = cw
        .iter()
        .map(|&c| {
            let ch = c.hi();
            let cc = if ch < clamp {
                T::from_f64(clamp)
            } else if ch > 1.0 - clamp {
                T::R1 - T::from_f64(clamp)
            } else {
                c
            };
            ((T::R1 - cc).divr(cc)).lnr()
        })
        .collect();
    let mut k = vec![T::R0; n * n];
    for m in 0..n {
        let km = kappa[m];
        if km.hi().abs() < 1e-13 {
            continue;
        }
        for j in 0..n {
            let vj = cv[j + m * n] * km;
            if vj.hi() == 0.0 {
                continue;
            }
            for i in 0..n {
                k[i + j * n] = k[i + j * n] + cv[i + m * n] * vj;
            }
        }
    }
    (k, cw)
}

/// DD 演算の自己検証 (呼び出し側の [PASS] ゲートで使う)。全て通れば true。
pub fn dd_self_test() -> bool {
    let mut ok = true;
    // (1e16 + 1) − 1e16 = 1 (f64 では 0)
    let a = (dd(1e16) + DD1) - dd(1e16);
    ok &= (a.hi - 1.0).abs() < 1e-30;
    // sqrt(2)² = 2
    let s2 = dd(2.0).sqrt();
    let r = s2 * s2 - dd(2.0);
    ok &= r.hi.abs() < 1e-30;
    // ln(exp(x)) = x
    let x = dd(0.37);
    let d = x.exp().ln() - x;
    ok &= d.hi.abs() < 1e-30;
    // ln(1e-30·e) − ln(1e-30) = 1 (深部スケールの相対精度)
    let d2 = (dd(1e-30) * dd(1.0).exp()).ln() - dd(1e-30).ln() - DD1;
    ok &= d2.hi.abs() < 1e-28;
    // sinpi_frac: sin(π/6) = 1/2, sin²+cos² = 1 (大整数角)
    let s16 = dd_sinpi_frac(1, 6) - dd(0.5);
    ok &= s16.hi.abs() < 1e-31;
    let (s, c) = (dd_sinpi_frac(12345, 787), dd_cospi_frac(12345, 787));
    let one = s * s + c * c - DD1;
    ok &= one.hi.abs() < 1e-30;
    // 加法定理: sin(π(a+b)/d) = sin cos + cos sin
    let (aa, bb, ddn) = (341i64, 173i64, 1009i64);
    let lhs = dd_sinpi_frac(aa + bb, ddn);
    let rhs = dd_sinpi_frac(aa, ddn) * dd_cospi_frac(bb, ddn)
        + dd_cospi_frac(aa, ddn) * dd_sinpi_frac(bb, ddn);
    ok &= (lhs - rhs).hi.abs() < 1e-30;
    ok
}
