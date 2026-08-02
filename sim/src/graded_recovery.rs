// graded_recovery — Majorana frame と Dirac pairing の識別境界 (v33.5, PROMPT/14)
//
// odd CAR ({γ_a, γ_b} = 2δ_ab) だけから通常の複素 fermion mode を一意に得ることは
// できない: 2N 本の Majorana generator の CAR は実直交変換 O(2N) に不変であり、
// Majorana を二本ずつ組にして annihilation/creation を定めるには、実 span 上の
// **直交複素構造 J (J² = −I)** が追加で必要である。
//
//   Majorana locality と Dirac-mode locality は同じ識別問題ではない —
//   後者は追加の U(1) charge / complex structure witness を必要とする。
//
// 出力階層 (凍結):
//   odd CAR only                    → MajoranaFrameOnly (O(2N) orbit —
//                                     Dirac モードは同値類。**禁止変換 20**:
//                                     MajoranaFrame → ComplexModeFactorization の
//                                     witness なし昇格は存在しない)
//   odd CAR + charge witness        → ComplexModeFactorization (U(N) gauge を
//                                     除いて回復 — J は一意・モード基底は自由)
//   縮退 witness (J² ≠ −I)          → Abstain(ComplexStructureUnresolved)
//   frame 上で線形でない witness    → Abstain(WitnessNotLinearOnFrame)
//
// witness からの複素構造の抽出: J_{ba} = ⟨γ_b, i[Q, γ_a]⟩ / ‖γ‖² — Q が
// 非縮退 U(1) charge なら J は実・反対称・J² = −I (直交複素構造)。複素構造が
// fermionic creation/annihilation 表現を定めるのは既存の数学的記述 (Araki 以来の
// 自己双対形式) と整合する — 本モジュールの寄与は識別可能性の fiber としての
// 型化・棄却・機械検証である。
//
// 注意: FermionicZ2Graded net の graded bracket graph は独立 Majorana 対で**空**
// (異なる γ は反可換 = graded-commuting) — 非可換グラフをいくら読んでも「どの
// 二本が一つの複素モードか」は決まらない (v335 [M5] が既存復元器の正直な Abstain
// を機械記録)。
//
// 一次ソース: docs/uft-v33.5.md / core.schema.yml (概念 + 禁止変換 20)。
// 整合は v335_graded_recovery が機械検査する。

use crate::operational_net::{anticommutator, commutator, hs_inner, hs_norm};
use crate::C64;

// ---------------------------------------------------------------- Majorana frame

/// 構成時拒否
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GradedRejection {
    /// γ_a がエルミートでない
    NotHermitian,
    /// {γ_a, γ_b} = 2δ_ab が破れる (CAR でない)
    CarViolation,
    /// 本数が奇数 (複素モードの候補にならない)
    OddCount,
}

impl GradedRejection {
    pub fn as_str(self) -> &'static str {
        match self {
            GradedRejection::NotHermitian => "not_hermitian",
            GradedRejection::CarViolation => "car_violation",
            GradedRejection::OddCount => "odd_count",
        }
    }
}

/// Majorana frame — CAR を構成時資格審査した 2N 本の族。O(2N) 回転で移り合う
/// frame は同じ CAR データを持つ (どの二本が一組かの情報は**ここに無い**)。
pub struct MajoranaFrame {
    gammas: Vec<Vec<C64>>,
    dim: usize,
}

impl MajoranaFrame {
    pub fn certify(gammas: Vec<Vec<C64>>, dim: usize) -> Result<Self, GradedRejection> {
        if gammas.len() % 2 != 0 {
            return Err(GradedRejection::OddCount);
        }
        for g in &gammas {
            for i in 0..dim {
                for j in 0..dim {
                    let a = g[i * dim + j];
                    let b = g[j * dim + i];
                    if (a.re - b.re).abs() > 1e-12 || (a.im + b.im).abs() > 1e-12 {
                        return Err(GradedRejection::NotHermitian);
                    }
                }
            }
        }
        for (a, ga) in gammas.iter().enumerate() {
            for (b, gb) in gammas.iter().enumerate() {
                let ac = anticommutator(ga, gb, dim);
                let want = if a == b { 2.0 } else { 0.0 };
                let mut dev = 0.0f64;
                for i in 0..dim {
                    for j in 0..dim {
                        let x = ac[i * dim + j];
                        let w = if i == j { want } else { 0.0 };
                        dev = dev.max((x.re - w).hypot(x.im));
                    }
                }
                if dev > 1e-9 {
                    return Err(GradedRejection::CarViolation);
                }
            }
        }
        Ok(MajoranaFrame { gammas, dim })
    }
    pub fn n_majorana(&self) -> usize {
        self.gammas.len()
    }
    pub fn dim(&self) -> usize {
        self.dim
    }
    pub fn gamma(&self, a: usize) -> &[C64] {
        &self.gammas[a]
    }
    /// 実係数線形結合 γ(v) = Σ v_a γ_a
    pub fn linear(&self, v: &[f64]) -> Vec<C64> {
        let n2 = self.dim * self.dim;
        let mut out = vec![C64::new(0.0, 0.0); n2];
        for (va, g) in v.iter().zip(self.gammas.iter()) {
            for (o, x) in out.iter_mut().zip(g.iter()) {
                *o = *o + x.scale(*va);
            }
        }
        out
    }
    /// 実直交回転 R (m×m, 行優先) を掛けた frame — CAR は resertify される
    pub fn rotated(&self, r: &[f64]) -> Result<MajoranaFrame, GradedRejection> {
        let m = self.gammas.len();
        let mut out = Vec::with_capacity(m);
        for i in 0..m {
            let row: Vec<f64> = (0..m).map(|j| r[i * m + j]).collect();
            out.push(self.linear(&row));
        }
        MajoranaFrame::certify(out, self.dim)
    }
}

// ---------------------------------------------------------------- 複素構造 witness

/// witness からの抽出の棄却理由
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GradedAbstainReason {
    /// witness の adjoint 作用が frame の実 span 上で閉じない (線形でない)
    WitnessNotLinearOnFrame,
    /// 抽出された J が直交複素構造でない (J² ≠ −I — 縮退/部分 witness)
    ComplexStructureUnresolved,
}

impl GradedAbstainReason {
    pub fn as_str(self) -> &'static str {
        match self {
            GradedAbstainReason::WitnessNotLinearOnFrame => "witness_not_linear_on_frame",
            GradedAbstainReason::ComplexStructureUnresolved => "complex_structure_unresolved",
        }
    }
}

/// 直交複素構造の証明書 — **唯一の構成は extract_complex_structure** (witness から
/// の抽出と J² = −I の機械検証)。MajoranaFrame 単独からの変換は存在しない
/// (禁止変換 20)。
pub struct ComplexStructureWitness {
    /// J (m×m 実, 行優先): i[Q, γ_a] = Σ_b J_{ba} γ_b
    j: Vec<f64>,
    m: usize,
    pub linearity_residual: f64,
    pub structure_residual: f64,
}

impl ComplexStructureWitness {
    pub fn j(&self) -> &[f64] {
        &self.j
    }
    pub fn m(&self) -> usize {
        self.m
    }
}

/// charge witness Q (偶演算子) から複素構造を抽出する:
///   M_a = i[Q, γ_a] を frame に展開 (実係数 J_{ba})。
///   - 展開残差 (span 外成分) が rel > 1e-9 → WitnessNotLinearOnFrame
///   - ‖J + Jᵀ‖ / ‖J² + I‖ が > 1e-9 → ComplexStructureUnresolved
pub fn extract_complex_structure(
    frame: &MajoranaFrame,
    q: &[C64],
) -> Result<ComplexStructureWitness, GradedAbstainReason> {
    let n = frame.dim();
    let m = frame.n_majorana();
    let gnorm2 = hs_norm(frame.gamma(0)).powi(2); // ‖γ‖² (= dim, CAR 資格済み)
    let mut j = vec![0.0f64; m * m];
    let mut lin_resid = 0.0f64;
    for a in 0..m {
        let c = commutator(q, frame.gamma(a), n);
        let ma: Vec<C64> = c.iter().map(|x| C64::new(-x.im, x.re)).collect(); // i·[Q, γ_a]
        let scale = hs_norm(&ma).max(1e-300);
        let mut coeffs = vec![0.0f64; m];
        for b in 0..m {
            let ip = hs_inner(frame.gamma(b), &ma);
            coeffs[b] = ip.re / gnorm2;
            // 虚部は線形性破れに数える
            lin_resid = lin_resid.max(ip.im.abs() / gnorm2.max(1e-300) / scale.max(1e-300));
        }
        let recon = frame.linear(&coeffs);
        let resid: f64 = recon
            .iter()
            .zip(ma.iter())
            .map(|(x, y)| (*x - *y).norm2())
            .sum::<f64>()
            .sqrt();
        lin_resid = lin_resid.max(resid / scale);
        for b in 0..m {
            j[b * m + a] = coeffs[b];
        }
    }
    if lin_resid > 1e-9 {
        return Err(GradedAbstainReason::WitnessNotLinearOnFrame);
    }
    // J の資格: 反対称・J² = −I
    let mut antisym = 0.0f64;
    for a in 0..m {
        for b in 0..m {
            antisym = antisym.max((j[a * m + b] + j[b * m + a]).abs());
        }
    }
    let mut j2 = vec![0.0f64; m * m];
    for a in 0..m {
        for b in 0..m {
            let mut s = 0.0;
            for k in 0..m {
                s += j[a * m + k] * j[k * m + b];
            }
            j2[a * m + b] = s;
        }
    }
    let mut structure = antisym;
    for a in 0..m {
        for b in 0..m {
            let want = if a == b { -1.0 } else { 0.0 };
            structure = structure.max((j2[a * m + b] - want).abs());
        }
    }
    if structure > 1e-9 {
        return Err(GradedAbstainReason::ComplexStructureUnresolved);
    }
    Ok(ComplexStructureWitness {
        j,
        m,
        linearity_residual: lin_resid,
        structure_residual: structure,
    })
}

// ---------------------------------------------------------------- graded recovery

/// graded lane の読み (凍結階層)
pub enum GradedRecoveryReading {
    /// odd CAR のみ — O(2N) orbit。Dirac モードは同値類 (昇格の門なし)
    MajoranaFrameOnly { n_majorana: usize },
    /// witness つき — U(N) gauge を除いて複素モードを回復。
    /// modes[i] = (a_i, a_i†) — CAR は呼び出し側で機械照合できる
    ComplexModeFactorization {
        n_modes: usize,
        modes: Vec<(Vec<C64>, Vec<C64>)>,
    },
    Abstain(GradedAbstainReason),
}

impl GradedRecoveryReading {
    pub fn as_str(&self) -> &'static str {
        match self {
            GradedRecoveryReading::MajoranaFrameOnly { .. } => "majorana_frame_only",
            GradedRecoveryReading::ComplexModeFactorization { .. } => {
                "complex_mode_factorization"
            }
            GradedRecoveryReading::Abstain(_) => "abstain",
        }
    }
}

/// graded recovery (凍結手順): witness なし → MajoranaFrameOnly (O(2N) orbit —
/// これ以上は読まない)。witness あり → J の不変平面 {v, Jv} を貪欲直交化で取り、
/// a_i = (γ(v_i) + i γ(J v_i))/2 を構成する (J は一意・平面の基底選択が U(N) gauge)。
pub fn recover_graded(
    frame: &MajoranaFrame,
    witness: Option<&ComplexStructureWitness>,
) -> GradedRecoveryReading {
    let m = frame.n_majorana();
    let w = match witness {
        None => {
            return GradedRecoveryReading::MajoranaFrameOnly { n_majorana: m };
        }
        Some(w) => w,
    };
    assert_eq!(w.m(), m, "witness と frame の本数が一致しない");
    let j = w.j();
    // 貪欲に J-不変平面を直交化して取る (J 直交 ⟹ 平面は自動的に直交)
    let mut planes: Vec<(Vec<f64>, Vec<f64>)> = Vec::new();
    let dot = |x: &[f64], y: &[f64]| -> f64 { x.iter().zip(y).map(|(a, b)| a * b).sum() };
    let jmul = |x: &[f64]| -> Vec<f64> {
        (0..m)
            .map(|a| (0..m).map(|b| j[a * m + b] * x[b]).sum())
            .collect()
    };
    for seed in 0..m {
        if planes.len() == m / 2 {
            break;
        }
        let mut v = vec![0.0f64; m];
        v[seed] = 1.0;
        for (p, q) in &planes {
            let cp = dot(&v, p);
            let cq = dot(&v, q);
            for t in 0..m {
                v[t] -= cp * p[t] + cq * q[t];
            }
        }
        let nv = dot(&v, &v).sqrt();
        if nv < 1e-9 {
            continue;
        }
        for t in v.iter_mut() {
            *t /= nv;
        }
        let jv = jmul(&v);
        planes.push((v, jv));
    }
    let mut modes = Vec::new();
    for (v, jv) in &planes {
        let gv = frame.linear(v);
        let gjv = frame.linear(jv);
        // a = (γ(v) + i γ(Jv))/2, a† = (γ(v) − i γ(Jv))/2
        let a: Vec<C64> = gv
            .iter()
            .zip(gjv.iter())
            .map(|(x, y)| C64::new((x.re - y.im) * 0.5, (x.im + y.re) * 0.5))
            .collect();
        let adag: Vec<C64> = gv
            .iter()
            .zip(gjv.iter())
            .map(|(x, y)| C64::new((x.re + y.im) * 0.5, (x.im - y.re) * 0.5))
            .collect();
        modes.push((a, adag));
    }
    GradedRecoveryReading::ComplexModeFactorization {
        n_modes: modes.len(),
        modes,
    }
}

// ---------------------------------------------------------------- 自己検査

/// graded_recovery の不変条件 (v335_graded_recovery が呼ぶ) — 1 モード (dim 2)
pub fn graded_recovery_self_test() -> Result<(), String> {
    let n = 2usize;
    let g1 = vec![
        C64::new(0.0, 0.0),
        C64::new(1.0, 0.0),
        C64::new(1.0, 0.0),
        C64::new(0.0, 0.0),
    ]; // X
    let g2 = vec![
        C64::new(0.0, 0.0),
        C64::new(0.0, -1.0),
        C64::new(0.0, 1.0),
        C64::new(0.0, 0.0),
    ]; // Y
    let frame = MajoranaFrame::certify(vec![g1.clone(), g2.clone()], n)
        .map_err(|e| e.as_str().to_string())?;
    // CAR 破れの拒否 (γ₂ を γ₁ に重ねる)
    if MajoranaFrame::certify(vec![g1.clone(), g1.clone()], n).is_ok() {
        return Err("CAR 破れが資格を通った".into());
    }
    // witness なし → MajoranaFrameOnly
    match recover_graded(&frame, None) {
        GradedRecoveryReading::MajoranaFrameOnly { n_majorana: 2 } => {}
        r => return Err(format!("witness なしの読みが {}", r.as_str())),
    }
    // Q = n = (1 + iγ₁γ₂)/2 = diag(0, 1)? — JW 規約で n = (I − Z)/2
    let q = vec![
        C64::new(0.0, 0.0),
        C64::new(0.0, 0.0),
        C64::new(0.0, 0.0),
        C64::new(1.0, 0.0),
    ];
    let w = extract_complex_structure(&frame, &q).map_err(|e| e.as_str().to_string())?;
    if w.structure_residual > 1e-9 {
        return Err(format!("J² = −I 残差 {:.1e}", w.structure_residual));
    }
    match recover_graded(&frame, Some(&w)) {
        GradedRecoveryReading::ComplexModeFactorization { n_modes: 1, modes } => {
            // CAR: {a, a†} = I・a² = 0
            let (a, ad) = &modes[0];
            let ac = anticommutator(a, ad, n);
            let mut dev = 0.0f64;
            for i in 0..n {
                for j2 in 0..n {
                    let x = ac[i * n + j2];
                    let want = if i == j2 { 1.0 } else { 0.0 };
                    dev = dev.max((x.re - want).hypot(x.im));
                }
            }
            let a2 = crate::operational_net::cmul(a, a, n);
            let na2 = hs_norm(&a2);
            if dev > 1e-9 || na2 > 1e-9 {
                return Err(format!("モード CAR 破れ dev {:.1e} / a² {:.1e}", dev, na2));
            }
        }
        r => return Err(format!("witness つきの読みが {}", r.as_str())),
    }
    // 縮退 witness (Q = 0) → ComplexStructureUnresolved
    let zero = vec![C64::new(0.0, 0.0); 4];
    match extract_complex_structure(&frame, &zero) {
        Err(GradedAbstainReason::ComplexStructureUnresolved) => {}
        _ => return Err("縮退 witness が棄却されない".into()),
    }
    Ok(())
}
