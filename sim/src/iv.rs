//! 区間演算 (v25.2) — 外向き丸めによる厳密包含
//!
//! 仮定 (ASM-IEEE754): binary64 の +,−,×,÷,√ は正しく丸められる (IEEE-754)。
//! 正しい丸めの結果は真値から高々 1/2 ulp — 丸め後 1 ulp の外向きナッジで
//! 真値の包含が保証される。**π と cos は libm を信頼しない**: π は
//! 「f64::consts::PI < π < next_up(PI)」の区間定数 (PI = 3.141592653589793116… は
//! 真値 3.141592653589793238… の切り捨て側にある — 既知の事実)、cos は剰余項
//! パディングつき Taylor 級数 (係数 1/(2k)! も区間で構成) で評価する。
//!
//! 用途: BZ moment の証明付き数値 (v252_bz_certificate)。設計上の注意:
//! 区間和は階層化する (行和 → 行の順序畳み込み) — 平坦な逐次和は部分和の
//! 桁上がりで外向きナッジが肥大する。

/// 次に大きい f64 (+∞ と NaN はそのまま)。ビット操作 — libm 非依存。
pub fn fnext_up(x: f64) -> f64 {
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x == 0.0 {
        return f64::from_bits(1); // 最小正 subnormal
    }
    let b = x.to_bits();
    if x > 0.0 {
        f64::from_bits(b + 1)
    } else {
        f64::from_bits(b - 1)
    }
}

/// 次に小さい f64
pub fn fnext_down(x: f64) -> f64 {
    -fnext_up(-x)
}

/// 閉区間 [lo, hi] — 不変条件 lo ≤ hi (構成時に検査)
#[derive(Clone, Copy, Debug)]
pub struct Iv {
    pub lo: f64,
    pub hi: f64,
}

pub const IV_ZERO: Iv = Iv { lo: 0.0, hi: 0.0 };
pub const IV_ONE: Iv = Iv { lo: 1.0, hi: 1.0 };

/// 点区間 (f64 で正確に表現される値)
pub fn iv(x: f64) -> Iv {
    assert!(x.is_finite());
    Iv { lo: x, hi: x }
}

/// 区間 π: f64::consts::PI は真の π の切り捨て — [PI, next_up(PI)] ∋ π (幅 1 ulp)
pub fn iv_pi() -> Iv {
    Iv {
        lo: std::f64::consts::PI,
        hi: fnext_up(std::f64::consts::PI),
    }
}

impl Iv {
    pub fn width(self) -> f64 {
        self.hi - self.lo
    }
    pub fn mid(self) -> f64 {
        0.5 * (self.lo + self.hi)
    }
    pub fn contains(self, x: f64) -> bool {
        self.lo <= x && x <= self.hi
    }
    pub fn contains_iv(self, o: Iv) -> bool {
        self.lo <= o.lo && o.hi <= self.hi
    }
    pub fn hull(self, o: Iv) -> Iv {
        Iv {
            lo: self.lo.min(o.lo),
            hi: self.hi.max(o.hi),
        }
    }
    pub fn add(self, o: Iv) -> Iv {
        Iv {
            lo: fnext_down(self.lo + o.lo),
            hi: fnext_up(self.hi + o.hi),
        }
    }
    pub fn sub(self, o: Iv) -> Iv {
        Iv {
            lo: fnext_down(self.lo - o.hi),
            hi: fnext_up(self.hi - o.lo),
        }
    }
    pub fn neg(self) -> Iv {
        Iv {
            lo: -self.hi,
            hi: -self.lo,
        }
    }
    pub fn mul(self, o: Iv) -> Iv {
        let c = [
            self.lo * o.lo,
            self.lo * o.hi,
            self.hi * o.lo,
            self.hi * o.hi,
        ];
        let lo = c.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Iv {
            lo: fnext_down(lo),
            hi: fnext_up(hi),
        }
    }
    /// 除算 — 0 ∉ o が前提 (破れたら abort = fail-closed)
    pub fn div(self, o: Iv) -> Iv {
        assert!(o.lo > 0.0 || o.hi < 0.0, "区間除算: 分母が 0 を含む");
        let c = [
            self.lo / o.lo,
            self.lo / o.hi,
            self.hi / o.lo,
            self.hi / o.hi,
        ];
        let lo = c.iter().cloned().fold(f64::INFINITY, f64::min);
        let hi = c.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        Iv {
            lo: fnext_down(lo),
            hi: fnext_up(hi),
        }
    }
    /// √ — lo ≥ 0 が前提。IEEE-754 の sqrt は正しい丸めなので 1 ulp ナッジで厳密
    pub fn sqrt(self) -> Iv {
        assert!(self.lo >= 0.0, "区間 sqrt: 負領域");
        Iv {
            lo: fnext_down(self.lo.sqrt()).max(0.0),
            hi: fnext_up(self.hi.sqrt()),
        }
    }
    /// ×0.5 (2 冪の乗算は正確 — ナッジ不要)
    pub fn half(self) -> Iv {
        Iv {
            lo: 0.5 * self.lo,
            hi: 0.5 * self.hi,
        }
    }
    /// ×2 (2 冪の乗算は正確 — ナッジ不要)
    pub fn dbl(self) -> Iv {
        Iv {
            lo: 2.0 * self.lo,
            hi: 2.0 * self.hi,
        }
    }
    /// ÷2^k (2 冪の除算は正確 — ナッジ不要。アンダーフロー域は対象外)
    pub fn div_pow2(self, k: u32) -> Iv {
        let f = (2.0f64).powi(-(k as i32));
        Iv {
            lo: self.lo * f,
            hi: self.hi * f,
        }
    }
    /// 平方 (x² ≥ 0 を使う — 0 跨ぎでも下限 0)
    pub fn sq(self) -> Iv {
        if self.lo >= 0.0 {
            Iv {
                lo: fnext_down(self.lo * self.lo),
                hi: fnext_up(self.hi * self.hi),
            }
        } else if self.hi <= 0.0 {
            Iv {
                lo: fnext_down(self.hi * self.hi),
                hi: fnext_up(self.lo * self.lo),
            }
        } else {
            let m = (-self.lo).max(self.hi);
            Iv {
                lo: 0.0,
                hi: fnext_up(m * m),
            }
        }
    }
}

/// cos の区間評価 (定義域 |x| ≤ 2.1 — BZ 格子 [0, π/2] 用)。
/// 交代級数 Σ(−1)^k x^{2k}/(2k)! を K=13 で打ち切り、剰余 |R| ≤ x²ᵏ⁺²/(2K+2)!
/// (x ≤ 2.1 で項は k≥1 から単調減少) を外側にパディングする。係数 1/(2k)! は
/// 整数の区間積から構成 (26! は f64 で正確に表せないため)。
pub fn cos_iv(x: Iv) -> Iv {
    assert!(x.lo >= -2.1 && x.hi <= 2.1, "cos_iv: 定義域外");
    const K: usize = 13;
    let x2 = x.sq();
    // c_k = 1/(2k)! を区間で
    let mut coef = Vec::with_capacity(K + 2);
    let mut fact = IV_ONE;
    let mut m = 0u64;
    for _ in 0..=(K + 1) {
        coef.push(IV_ONE.div(fact));
        fact = fact.mul(iv((m + 1) as f64)).mul(iv((m + 2) as f64));
        m += 2;
    }
    // Horner: P = c₀ − x²(c₁ − x²(c₂ − …))
    let mut acc = coef[K];
    for k in (0..K).rev() {
        acc = coef[k].sub(x2.mul(acc));
    }
    // 剰余パディング: |R| ≤ (x².hi)^{K+1}/(2K+2)!
    let mut p = iv(1.0);
    for _ in 0..=K {
        p = p.mul(Iv { lo: 0.0, hi: x2.hi });
    }
    let rem = p.mul(coef[K + 1]).hi;
    Iv {
        lo: fnext_down(acc.lo - rem),
        hi: fnext_up(acc.hi + rem),
    }
}

/// 算術幾何平均の区間版。任意ステップで min(a,b) ≤ AGM ≤ max(a,b) が成り立つので
/// hull(a,b) は常に厳密な包含 — 反復は幅を締めるだけで、幅が減らなくなったら返す。
pub fn agm_iv(a0: Iv, b0: Iv) -> Iv {
    assert!(a0.lo > 0.0 && b0.lo > 0.0, "agm_iv: 正領域のみ");
    let (mut a, mut b) = (a0, b0);
    let mut best = a.hull(b);
    for _ in 0..60 {
        let na = a.add(b).half();
        let nb = a.mul(b).sqrt();
        a = na;
        b = nb;
        let h = a.hull(b);
        if h.width() >= best.width() {
            return best;
        }
        best = h;
    }
    best
}

/// g(μ) の区間版 — 入力は μ² (BZ では μ² = cos²ky + cos²kz が一次量)。
/// g = 1/AGM(1, √(1+μ²))
pub fn g_iv_mu2(mu2: Iv) -> Iv {
    IV_ONE.div(agm_iv(IV_ONE, IV_ONE.add(mu2).sqrt()))
}

/// 自己検証 — 外部照合値 (数学定数) の包含と幅
pub fn iv_self_test() -> bool {
    let mut ok = true;
    // 外向き丸めの基本
    ok &= fnext_up(1.0) > 1.0 && fnext_down(1.0) < 1.0;
    ok &= fnext_up(0.0) > 0.0 && fnext_down(0.0) < 0.0;
    // 1/3 の包含
    let third = IV_ONE.div(iv(3.0));
    ok &= third.contains(1.0 / 3.0) && third.width() < 1e-15;
    // √2 の包含と平方の復元
    let r2 = iv(2.0).sqrt();
    ok &= r2.contains(std::f64::consts::SQRT_2);
    ok &= r2.sq().contains(2.0);
    // π: 幅 1 ulp
    let pi = iv_pi();
    ok &= pi.contains(std::f64::consts::PI) && pi.width() <= 5e-16;
    // cos の厳密値: cos(π/3) = 1/2, cos(π/2) = 0, cos(0) = 1
    ok &= cos_iv(pi.div(iv(3.0))).contains(0.5);
    ok &= cos_iv(pi.half()).contains(0.0);
    ok &= cos_iv(IV_ZERO).contains(1.0);
    ok &= cos_iv(pi.div(iv(3.0))).width() < 1e-14;
    // AGM(1,√2) = 1.198140234735592207… (ガウス定数の逆数, 既知の数学定数)
    let agm = agm_iv(IV_ONE, iv(2.0).sqrt());
    ok &= agm.contains(1.198_140_234_735_592_2) && agm.width() < 1e-14;
    // g(μ=1) = ガウス定数
    ok &= g_iv_mu2(IV_ONE).contains(0.834_626_841_674_073_2);
    ok
}
