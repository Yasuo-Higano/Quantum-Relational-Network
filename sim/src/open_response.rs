//! v35.2 開放系 signed-response 観測商 (PROMPT/16 §5–6)
//!
//! 有限モード・number-conserving・quasi-free (Markov) 開放系:
//!   ρ̇ = −i[H,ρ] + Σ_k D[L_k]ρ + Σ_l D[G_l]ρ,
//!   H = Σ h_{ab} c†_a c_b (h エルミート), L_k = Σ_a ℓ_{k,a} c_a (loss),
//!   G_l = Σ_a g_{l,a} c†_a (gain)
//! の normal covariance C_{ab} = ⟨c†_b c_a⟩ は閉じた affine ODE に従う:
//!   Ċ = X C + C X† + Y,  X = −ih − ½(Λᵀ + M),  Y = M,
//!   Λ = Σ_k ℓ_k ℓ_k† (loss は**転置** Λᵀ = conj(Λ) で入る),  M = Σ_l g_l g_l†
//! (反交換子の縮約を機械で追うと loss 側だけ転置が出る — この規約は v352 [G1] が
//!  dense Lindblad (2^N 次元 Jordan–Wigner) との一致で較正・凍結した。初版の
//!  Λ 非転置は [G1] が 5e-2 の不一致で検出 — 較正セルの存在理由)。
//!
//! 中心となる型の規律 (禁止変換 32/33):
//!   OpenSignedCovarianceProbe → **EffectiveDriftTopology** (曲率が読むのは X)
//!   EffectiveDriftTopology ↛ HamiltonianTopology  (禁止変換 32 — GQF-3 反例対:
//!     coherent hopping と collective loss が同一曲率。昇格の唯一の門は
//!     DissipatorLocalityCertificate = 検証済み「cross-node 散逸 drift = 0」)
//!   ChargeNonconservingResponse ↛ HamiltonianPairingWitness  (禁止変換 33 —
//!     GQF-5 反例対: loss と BdG pairing が同一の電荷応答値。門は
//!     DissipativeChargeConservationCertificate = 検証済み「散逸ゼロ/電荷保存散逸」)
//!
//! 観測商 (PROMPT/16 §6): 登録契約 C = (probe 族, 観測量族, 時刻族) の応答
//! Φ_C(L)_{ji}(t) が識別するのは生成子そのものではなく核で割った商 [L]_C。
//! quasi-free density 応答では少なくとも local phase gauge (X ↦ DXD†)・
//! global frequency (X ↦ X + iωI)・複素共役 (X ↦ conj X) が応答を保存する
//! (v352 [C7]/[C8] が機械実証)。よって reader は X でなく同値類を返す。
//!
//! 対応する Lean 定理: proofs/OpenQuotient.lean (GQF-1..5 の jet 恒等式・反例対 —
//! 16 定理)。本モジュールは解析側 (e^{Xt}・有限差分・有限 shot) を担う。

use crate::finite_data::cp_interval;
use crate::record_v2::{assess, CorrelationVerdict, RecordData};
use crate::{C64, CONE, CZERO};

pub type CMat = Vec<Vec<C64>>;

pub fn cs(re: f64, im: f64) -> C64 {
    C64 { re, im }
}

// ---------------------------------------------------------------- 複素行列の基本演算

pub fn mat_zero(n: usize) -> CMat {
    vec![vec![CZERO; n]; n]
}

pub fn mat_eye(n: usize) -> CMat {
    let mut m = mat_zero(n);
    for (i, row) in m.iter_mut().enumerate() {
        row[i] = CONE;
    }
    m
}

pub fn mat_add(a: &CMat, b: &CMat) -> CMat {
    a.iter()
        .zip(b)
        .map(|(ra, rb)| ra.iter().zip(rb).map(|(&x, &y)| x + y).collect())
        .collect()
}

pub fn mat_sub(a: &CMat, b: &CMat) -> CMat {
    a.iter()
        .zip(b)
        .map(|(ra, rb)| ra.iter().zip(rb).map(|(&x, &y)| x - y).collect())
        .collect()
}

pub fn mat_scale(s: C64, a: &CMat) -> CMat {
    a.iter()
        .map(|r| r.iter().map(|&x| s * x).collect())
        .collect()
}

pub fn mat_mul(a: &CMat, b: &CMat) -> CMat {
    let n = a.len();
    let m = b[0].len();
    let k = b.len();
    let mut c = vec![vec![CZERO; m]; n];
    for i in 0..n {
        for l in 0..k {
            let ail = a[i][l];
            if ail.re == 0.0 && ail.im == 0.0 {
                continue;
            }
            for j in 0..m {
                c[i][j] = c[i][j] + ail * b[l][j];
            }
        }
    }
    c
}

/// 共役転置 A†
pub fn mat_adj(a: &CMat) -> CMat {
    let n = a.len();
    let m = a[0].len();
    (0..m)
        .map(|j| (0..n).map(|i| a[i][j].conj()).collect())
        .collect()
}

/// 成分ごとの複素共役 (転置なし)
pub fn mat_conj(a: &CMat) -> CMat {
    a.iter()
        .map(|r| r.iter().map(|x| x.conj()).collect())
        .collect()
}

pub fn mat_trace(a: &CMat) -> C64 {
    let mut t = CZERO;
    for (i, row) in a.iter().enumerate() {
        t = t + row[i];
    }
    t
}

pub fn mat_max_abs(a: &CMat) -> f64 {
    a.iter()
        .flat_map(|r| r.iter())
        .map(|x| x.norm2().sqrt())
        .fold(0.0, f64::max)
}

/// スケーリング & スクエアリング + Taylor の複素行列指数 (小次元・決定的)
pub fn expm(a: &CMat) -> CMat {
    let n = a.len();
    let norm = a
        .iter()
        .map(|r| r.iter().map(|x| x.norm2().sqrt()).sum::<f64>())
        .fold(0.0, f64::max);
    let mut s = 0u32;
    let mut scale = 1.0;
    while norm * scale > 0.5 {
        s += 1;
        scale *= 0.5;
    }
    let b = mat_scale(cs(scale, 0.0), a);
    // Taylor: Σ B^k/k! — ‖B‖ ≤ 0.5 なので 20 項で ~1e-19
    let mut term = mat_eye(n);
    let mut sum = mat_eye(n);
    for k in 1..=20 {
        term = mat_mul(&term, &b);
        term = mat_scale(cs(1.0 / k as f64, 0.0), &term);
        sum = mat_add(&sum, &term);
    }
    for _ in 0..s {
        sum = mat_mul(&sum, &sum);
    }
    sum
}

// ---------------------------------------------------------------- quasi-free open model

/// number-conserving quasi-free 開放系の宣言。
/// 構成は `new` の門のみ — h の非エルミート・pairing 項は構成時拒否 (fail-closed)。
pub struct QuasiFreeOpenModel {
    pub n: usize,
    pub h: CMat,
    /// 各要素が jump ベクトル ℓ_k (L_k = Σ_a ℓ_{k,a} c_a)
    pub loss: Vec<Vec<C64>>,
    /// 各要素が gain ベクトル g_l (G_l = Σ_a g_{l,a} c†_a)
    pub gain: Vec<Vec<C64>>,
}

/// 構成時拒否の理由 (fail-closed — 強制回答しない)
#[derive(Debug)]
pub enum OpenLaneRefusal {
    /// h が非エルミート (宣言ミスまたは非物理)
    NonHermitianHamiltonian { max_dev: f64 },
    /// pairing (Nambu) 項の宣言 — 本 lane は number-conserving のみ (OutOfDomain)
    PairingOutOfDomain { max_delta: f64 },
}

impl QuasiFreeOpenModel {
    /// 唯一の構成門: h エルミート検査 + pairing 宣言の拒否。
    /// `delta` は宣言された pairing block (Nambu ならここに Δ が入る) — 非零なら
    /// 本 lane では OutOfDomain (禁止変換の運用形: 強制的に読まない)。
    pub fn new(
        h: CMat,
        loss: Vec<Vec<C64>>,
        gain: Vec<Vec<C64>>,
        delta: Option<&CMat>,
    ) -> Result<Self, OpenLaneRefusal> {
        let n = h.len();
        if let Some(d) = delta {
            let md = mat_max_abs(d);
            if md > 0.0 {
                return Err(OpenLaneRefusal::PairingOutOfDomain { max_delta: md });
            }
        }
        let dev = mat_max_abs(&mat_sub(&h, &mat_adj(&h)));
        if dev > 1e-12 {
            return Err(OpenLaneRefusal::NonHermitianHamiltonian { max_dev: dev });
        }
        Ok(QuasiFreeOpenModel { n, h, loss, gain })
    }

    /// Λᵀ = conj(Σ_k ℓ_k ℓ_k†) — covariance 規約 C_{ab} = ⟨c†_b c_a⟩ での
    /// loss 散逸 drift ([G1] で較正済みの転置規約)
    pub fn lambda_t(&self) -> CMat {
        let mut m = mat_zero(self.n);
        for l in &self.loss {
            for a in 0..self.n {
                for b in 0..self.n {
                    m[a][b] = m[a][b] + l[a].conj() * l[b];
                }
            }
        }
        m
    }

    /// M = Σ_l g_l g_l† (gain 注入行列 — [G1] で較正済み: gain 側は転置しない)
    pub fn gain_matrix(&self) -> CMat {
        let mut m = mat_zero(self.n);
        for g in &self.gain {
            for a in 0..self.n {
                for b in 0..self.n {
                    m[a][b] = m[a][b] + g[a] * g[b].conj();
                }
            }
        }
        m
    }

    /// 有効 drift X = −ih − ½(Λᵀ + M)
    pub fn effective_drift(&self) -> CMat {
        let ih = mat_scale(cs(0.0, -1.0), &self.h);
        let gamma = mat_add(&self.lambda_t(), &self.gain_matrix());
        mat_add(&ih, &mat_scale(cs(-0.5, 0.0), &gamma))
    }

    /// affine 注入 Y = M
    pub fn injection(&self) -> CMat {
        self.gain_matrix()
    }

    /// 総散逸 drift Γ = Λᵀ + M (dissipator locality 証明書の検査対象)
    pub fn dissipative_drift(&self) -> CMat {
        mat_add(&self.lambda_t(), &self.gain_matrix())
    }
}

/// C(t) = e^{Xt} C₀ e^{X†t} + ∫₀ᵗ e^{Xs} Y e^{X†s} ds — Van Loan block で厳密
/// (M = [[X, Y],[0, −X†]], e^{Mt} = [[F, G],[0, H]] → C(t) = F C₀ F† + G F†)
pub fn evolve_covariance(x: &CMat, y: &CMat, c0: &CMat, t: f64) -> CMat {
    let n = x.len();
    let mut big = mat_zero(2 * n);
    for i in 0..n {
        for j in 0..n {
            big[i][j] = x[i][j] * cs(t, 0.0);
            big[i][n + j] = y[i][j] * cs(t, 0.0);
            big[n + i][n + j] = -(x[j][i].conj()) * cs(t, 0.0);
        }
    }
    let e = expm(&big);
    let mut f = mat_zero(n);
    let mut g = mat_zero(n);
    for i in 0..n {
        for j in 0..n {
            f[i][j] = e[i][j];
            g[i][j] = e[i][n + j];
        }
    }
    let fd = mat_adj(&f);
    mat_add(&mat_mul(&mat_mul(&f, c0), &fd), &mat_mul(&g, &fd))
}

// ---------------------------------------------------------------- 曲率統計と型

/// 登録曲率統計 (厳密側): w_{ji} = ‖P_j X P_i‖²_F (P = 宣言 node 射影の族)。
/// 開放系でこれが読むのは **有効 drift の topology** であって Hamiltonian ではない。
pub fn curvature_exact(x: &CMat, nodes: &[Vec<usize>]) -> Vec<Vec<f64>> {
    let m = nodes.len();
    let mut w = vec![vec![0.0; m]; m];
    for (j, nj) in nodes.iter().enumerate() {
        for (i, ni) in nodes.iter().enumerate() {
            if i == j {
                continue;
            }
            let mut s = 0.0;
            for &a in nj {
                for &b in ni {
                    s += x[a][b].norm2();
                }
            }
            w[j][i] = s;
        }
    }
    w
}

/// OpenSignedCovarianceProbe の読み — 到達型は EffectiveDriftTopology のみ
#[derive(Debug, Clone)]
pub struct EffectiveDriftTopology {
    pub w: Vec<Vec<f64>>,
}

/// Hamiltonian topology — **唯一の構成子は promote_with_certificate** (禁止変換 32:
/// EffectiveDriftTopology からの直接変換 (From/Into) は存在しない)
#[derive(Debug, Clone)]
pub struct HamiltonianTopology {
    _gate: (),
    pub w: Vec<Vec<f64>>,
}

/// 検証済み dissipator locality 証明書 — 宣言でなく検査 (v33.2「門は較正」):
/// Γ = Λ + M の cross-node block の Frobenius ノルムが bar 以下であることを
/// 構成時に機械検査する。
#[derive(Debug)]
pub struct DissipatorLocalityCertificate {
    _gate: (),
    pub max_cross_block: f64,
    pub bar: f64,
}

#[derive(Debug)]
pub enum PromotionRefusal {
    /// cross-node 散逸 drift が bar 超 — GQF-3 反例対の型 (曲率は X しか読めない)
    OffDiagonalDissipator { max_cross_block: f64, bar: f64 },
}

pub fn certify_dissipator_locality(
    model: &QuasiFreeOpenModel,
    nodes: &[Vec<usize>],
    bar: f64,
) -> Result<DissipatorLocalityCertificate, PromotionRefusal> {
    let gamma = model.dissipative_drift();
    let mut max_cross = 0.0f64;
    for (j, nj) in nodes.iter().enumerate() {
        for (i, ni) in nodes.iter().enumerate() {
            if i == j {
                continue;
            }
            let mut s = 0.0;
            for &a in nj {
                for &b in ni {
                    s += gamma[a][b].norm2();
                }
            }
            max_cross = max_cross.max(s.sqrt());
        }
    }
    if max_cross <= bar {
        Ok(DissipatorLocalityCertificate {
            _gate: (),
            max_cross_block: max_cross,
            bar,
        })
    } else {
        Err(PromotionRefusal::OffDiagonalDissipator {
            max_cross_block: max_cross,
            bar,
        })
    }
}

/// GQF-4 の解錠: EffectiveDriftTopology + DissipatorLocalityCertificate →
/// HamiltonianTopology (‖P_j X P_i‖² = ‖P_j h P_i‖² が厳密に成立する領域)
pub fn promote_with_certificate(
    eff: &EffectiveDriftTopology,
    _cert: &DissipatorLocalityCertificate,
) -> HamiltonianTopology {
    HamiltonianTopology {
        _gate: (),
        w: eff.w.clone(),
    }
}

// ---------------------------------------------------------------- gauge (観測商)

/// local phase gauge: X ↦ D X D† (D = diag(e^{iφ_a}))
pub fn local_phase_gauge(x: &CMat, phis: &[f64]) -> CMat {
    let n = x.len();
    let mut y = mat_zero(n);
    for i in 0..n {
        for j in 0..n {
            let ph = cs((phis[i] - phis[j]).cos(), (phis[i] - phis[j]).sin());
            y[i][j] = ph * x[i][j];
        }
    }
    y
}

/// global frequency: X ↦ X + iωI
pub fn global_frequency_shift(x: &CMat, omega: f64) -> CMat {
    let n = x.len();
    let mut y = x.clone();
    for (i, row) in y.iter_mut().enumerate().take(n) {
        row[i] = row[i] + cs(0.0, omega);
    }
    y
}

/// 登録契約の全時刻応答表: Φ(t)_{ji} = Tr(P_j e^{Xt} P_i e^{X†t})
/// (1 次元 node なら |propagator 成分|²) — 観測商の同値判定に使う
pub fn response_table(x: &CMat, nodes: &[Vec<usize>], times: &[f64]) -> Vec<Vec<Vec<f64>>> {
    let n = x.len();
    times
        .iter()
        .map(|&t| {
            let f = expm(&mat_scale(cs(t, 0.0), x));
            let m = nodes.len();
            let mut tab = vec![vec![0.0; m]; m];
            for (j, nj) in nodes.iter().enumerate() {
                for (i, ni) in nodes.iter().enumerate() {
                    let mut s = 0.0;
                    for &a in nj {
                        for &b in ni {
                            let _ = n;
                            s += f[a][b].norm2();
                        }
                    }
                    tab[j][i] = s;
                }
            }
            tab
        })
        .collect()
}

// ---------------------------------------------------------------- 電荷応答 (GQF-5 の型)

/// 電荷応答統計 (number-conserving covariance lane): dQ/dt = Re Tr((X+X†)C + Y)
pub fn charge_response(x: &CMat, y: &CMat, c: &CMat) -> f64 {
    let xs = mat_add(x, &mat_adj(x));
    (mat_trace(&mat_mul(&xs, c)) + mat_trace(y)).re
}

/// 電荷非保存応答の marker — ここから pairing への変換は存在しない (禁止変換 33)
#[derive(Debug, Clone, Copy)]
pub struct ChargeNonconservingResponse(pub f64);

/// 検証済み「散逸が電荷を保存する」証明書 — 本 lane では Λ = M = 0 の機械検査
/// (線形 jump は必ず電荷を ±1 変えるため、quasi-free lane での電荷保存散逸は
/// 散逸ゼロに限る — 一般 GKLS の quadratic dephasing は本 lane 外)
#[derive(Debug)]
pub struct DissipativeChargeConservationCertificate {
    _gate: (),
    pub max_dissipation: f64,
}

pub fn certify_charge_conserving_dissipation(
    model: &QuasiFreeOpenModel,
) -> Result<DissipativeChargeConservationCertificate, f64> {
    let md = mat_max_abs(&model.dissipative_drift());
    if md == 0.0 {
        Ok(DissipativeChargeConservationCertificate {
            _gate: (),
            max_dissipation: md,
        })
    } else {
        Err(md)
    }
}

/// pairing witness — **唯一の構成子** (電荷非保存 + 散逸電荷保存証明書)
#[derive(Debug)]
pub struct HamiltonianPairingWitness {
    _gate: (),
    pub response: f64,
}

pub fn pairing_witness_with_certificate(
    resp: ChargeNonconservingResponse,
    _cert: &DissipativeChargeConservationCertificate,
) -> Option<HamiltonianPairingWitness> {
    if resp.0 != 0.0 {
        Some(HamiltonianPairingWitness {
            _gate: (),
            response: resp.0,
        })
    } else {
        None
    }
}

// ---------------------------------------------------------------- 有限 shot 読み (GQF-6)

/// 有限 shot の曲率読み (登録推定器 — Richardson 4 点):
///   K̂ = [8·Δ̂(δ) − Δ̂(2δ)] / (8εδ²),  Δ̂(t) = n̂⁺_j(t) − n̂⁻_j(t)
/// Richardson は δ³ 項を厳密に消し、残余は K̂ = w − 2c₄δ² —
/// |c₄| ≤ (2/3)R⁴ (t⁴ 係数の粗い上界) から登録バイアス B = (4/3)R⁴δ²。
/// 区間 = CP 区間の線形伝播 + B。裁定は RobustVerdict 意味論 (v34.3):
/// 区間全体で支持判定が一致するときだけ回答。
pub struct FiniteShotCurvatureReader {
    pub eps: f64,
    pub delta: f64,
    pub alpha: f64,
    /// 宣言された drift ノルム上界 (観測契約の一部 — バイアス上界に入る)
    pub x_norm_bound: f64,
    /// 支持判定の登録バー (w > tau = edge)
    pub tau: f64,
    pub min_shots: usize,
}

#[derive(Debug, Clone, PartialEq)]
pub enum CurvatureVerdict {
    RobustEdge,
    RobustNoEdge,
    Straddled,
    InsufficientObservation,
    /// iid 契約の破れ (相関 shot) — 読まないのが正答
    OutOfDomainCorrelated,
}

impl FiniteShotCurvatureReader {
    /// shot 群から w_{ji} の区間と裁定。probe ±/時刻 δ, 2δ の 4 系列の {0,1} shot。
    /// 系列はまず v35.1 の相関粒度ゲート (record_v2) を通す — IidConsistent で
    /// ないものは読まない (禁止変換 29/30/31 の運用)。
    pub fn read(
        &self,
        shots_plus_d: &[u8],
        shots_minus_d: &[u8],
        shots_plus_2d: &[u8],
        shots_minus_2d: &[u8],
    ) -> CurvatureVerdict {
        let series = [shots_plus_d, shots_minus_d, shots_plus_2d, shots_minus_2d];
        if series.iter().any(|s| s.len() < self.min_shots) {
            return CurvatureVerdict::InsufficientObservation;
        }
        for s in series {
            match assess(&RecordData::OrderedShots(s.to_vec()), self.alpha / 8.0) {
                CorrelationVerdict::IidConsistent(_) => {}
                _ => return CurvatureVerdict::OutOfDomainCorrelated,
            }
        }
        // 各系列の CP 区間 (同時 α: 4 系列 + 相関ゲートで α/8 ずつ)
        let iv = |s: &[u8]| {
            let k = s.iter().filter(|&&x| x == 1).count();
            cp_interval(k, s.len(), self.alpha / 8.0)
        };
        let (pd_lo, pd_hi) = iv(shots_plus_d);
        let (md_lo, md_hi) = iv(shots_minus_d);
        let (p2_lo, p2_hi) = iv(shots_plus_2d);
        let (m2_lo, m2_hi) = iv(shots_minus_2d);
        // Δ(δ) ∈ [pd_lo − md_hi, pd_hi − md_lo], Δ(2δ) 同様
        let dd = (pd_lo - md_hi, pd_hi - md_lo);
        let d2 = (p2_lo - m2_hi, p2_hi - m2_lo);
        // K̂ = [8Δ(δ) − Δ(2δ)]/(8εδ²) の区間演算
        let denom = 8.0 * self.eps * self.delta * self.delta;
        let k_lo = (8.0 * dd.0 - d2.1) / denom;
        let k_hi = (8.0 * dd.1 - d2.0) / denom;
        // 登録バイアス上界: |K̂ − w| ≤ (4/3) R⁴ δ² (Richardson 後の t⁴ 剰余)
        let bias = 4.0 / 3.0 * self.x_norm_bound.powi(4) * self.delta * self.delta;
        let (w_lo, w_hi) = (k_lo - bias, k_hi + bias);
        if w_lo > self.tau {
            CurvatureVerdict::RobustEdge
        } else if w_hi <= self.tau {
            CurvatureVerdict::RobustNoEdge
        } else {
            CurvatureVerdict::Straddled
        }
    }
}

// ---------------------------------------------------------------- 自己検証

pub fn open_response_self_test() -> Result<(), String> {
    // expm の既知値: exp([[0, 1],[−1, 0]]·θ) = 回転行列
    let th = 0.7f64;
    let a = vec![
        vec![CZERO, cs(th, 0.0)],
        vec![cs(-th, 0.0), CZERO],
    ];
    let e = expm(&a);
    if (e[0][0].re - th.cos()).abs() > 1e-14 || (e[0][1].re - th.sin()).abs() > 1e-14 {
        return Err(format!("expm 回転行列 dev {}", (e[0][0].re - th.cos()).abs()));
    }
    // Van Loan: X = 0, Y = I → C(t) = C0 + tI
    let n = 2;
    let x = mat_zero(n);
    let y = mat_eye(n);
    let c0 = mat_scale(cs(0.25, 0.0), &mat_eye(n));
    let ct = evolve_covariance(&x, &y, &c0, 0.5);
    if (ct[0][0].re - 0.75).abs() > 1e-13 {
        return Err(format!("Van Loan 積分 dev {}", (ct[0][0].re - 0.75).abs()));
    }
    // 禁止変換 32/33 の型: 直接構成子が非公開であること (コンパイル時性質 —
    // ここでは gate 関数の存在のみ確認)
    Ok(())
}
