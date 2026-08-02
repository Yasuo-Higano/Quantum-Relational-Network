// laboratory_interface — Certified Laboratory Interface の型契約 (v33.2, PROMPT/14)
//
// 第三十三期テーゼ「可アクセス性は作用素単体の属性ではなく、系・制御器・測定器・
// 資源・誤差証明書の関係である」の正側の型実装。v33.1 が固定した境界 (primitive
// 選別の循環 — independently accessible な entangler 1 本で読みが併合する) を、
// 「accessible operations を入力として与えた」から「各 operation が、どの command・
// 校正・合成列・誤差・資源によって accessible なのかを証明した」へ進める。
//
// 設計の核:
//   1. **DeclaredOperation は資格なし** — AccessibleOperation への直接変換は存在
//      しない (**禁止変換 14**)。通過できる門は 3 つの出自証明書のみ:
//      較正 (CalibrationRecord → IndependentAddressabilityCertificate)・
//      合成 (SynthesisCertificate — 機械実行で検証される bracket/線形結合の列)・
//      トモグラフィ (TomographyCertificate — 情報完全状態族からの再構成残差)。
//      証明書は対象行列の sha256 に**結束**され、別の行列への流用は構成時拒否。
//      文字列 provenance は新 lane に存在しない (出自 = 型 + sha256)。
//   2. **数学的分解は独立 addressability を与えない** (**禁止変換 15**): 装置が
//      u(t)(X₁+X₂) しか操作できないとき、行列を X₁ と X₂ に分解して「二つの独立
//      primitive」とは扱えない — 同一 command を 2 標的に割り当てる較正は
//      Jacobian rank 1 < 2 で構成時拒否。ただし bracket 合成路 (例: Z₂ knob との
//      入れ子交換子) が実在すれば Synthesized として資格が立つ — **可アクセス性は
//      interface との関係で決まる** (同じ X₁ が interface A では不可・B では可)。
//   3. **context は役割ごとに別型** (v33.1 の予告の履行): ControlContext は独立
//      addressability・MeasurementContext は joint measurability (**可換性より
//      広い** — 非可換な unsharp 対も共同測定可能: 明示的 joint POVM の機械検証が
//      資格)・PreparationFamily は凸可達性・DriftRegime は安定性。いずれも
//      CertifiedCommutator (代数) からの変換は存在しない (**禁止変換 16** —
//      禁止変換 12 の役割別展開)。
//   4. **ResourceBudget は成分ごとの半順序** — 単一スカラーへの縮約 (Ord) を
//      実装しない (恣意的な重み付き和は新しい選択バイアス — v33.3 の
//      resource-indexed profile の土台)。
//
// 一次ソース: docs/uft-v33.2.md / core.schema.yml (概念 + 禁止変換 14/15/16)。
// 整合は v332_certified_interface が機械検査する。

use crate::operational_net::{
    commutator, hs_inner, hs_norm, push_ortho, CommutationGrading, ControlGenerator,
    DriftGenerator, MarkedRecoveryDetail, MeasurementEffect, OpId, OpKind, OperationalNet,
    OperatorParity, Preparation, PrimitiveOperation, RecoveryInputRejection,
};
use crate::{jacobi_eigh, sha256_hex, C64};

// ---------------------------------------------------------------- 共通素子

/// エルミート行列 (C64) の固有値 (2n×2n 実対称埋め込み — 各固有値 2 重)
fn herm_evals_c64(m: &[C64], n: usize) -> Vec<f64> {
    let d = 2 * n;
    let mut big = vec![0.0; d * d];
    for i in 0..n {
        for j in 0..n {
            big[i * d + j] = m[i * n + j].re;
            big[(i + n) * d + (j + n)] = m[i * n + j].re;
            big[i * d + (j + n)] = -m[i * n + j].im;
            big[(i + n) * d + j] = m[i * n + j].im;
        }
    }
    let (evals, _) = jacobi_eigh(&big, d);
    evals
}

/// 行列の正準直列化の sha256 (証明書の対象結束に使う — 桁は 1e-12 で固定)
pub fn matrix_hash(m: &[C64], n: usize) -> String {
    let mut s = format!("dim={};", n);
    for c in m {
        s.push_str(&format!("{:.12e},{:.12e};", c.re, c.im));
    }
    sha256_hex(s.as_bytes())
}

/// 誤差区間 [lo, hi] (0 ≤ lo ≤ hi) — 閾値との 3 値比較は跨ぎを Abstain 扱いにする
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BoundInterval {
    lo: f64,
    hi: f64,
}

impl BoundInterval {
    pub fn new(lo: f64, hi: f64) -> Result<Self, &'static str> {
        if !(lo.is_finite() && hi.is_finite()) || lo < 0.0 || hi < lo {
            return Err("区間が不正 (0 ≤ lo ≤ hi の有限値であること)");
        }
        Ok(BoundInterval { lo, hi })
    }
    pub fn lo(&self) -> f64 {
        self.lo
    }
    pub fn hi(&self) -> f64 {
        self.hi
    }
    /// バーとの 3 値: Ok(true) = 全区間がバー以下 / Ok(false) = 全区間が超過 / Err = 跨ぎ
    pub fn within(&self, bar: f64) -> Result<bool, &'static str> {
        if self.hi <= bar {
            Ok(true)
        } else if self.lo > bar {
            Ok(false)
        } else {
            Err("区間がバーを跨ぐ (Abstain — 強制判定の禁止)")
        }
    }
}

/// 構成時拒否 — 資格審査に通らない laboratory interface の申告
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InterfaceRejection {
    /// 標的族が退化 (一次従属) — addressability の対象にならない
    DegenerateTargets,
    /// command と標的の本数が一致しない (較正は 1:1 の割り当て)
    CommandTargetCountMismatch,
    /// command Jacobian の rank が標的数に足りない (数学的分解 ≠ 独立 addressability — 禁止変換 15)
    InsufficientCommandRank,
    /// σ_min がバー未満 (弱すぎる駆動は独立 addressability を名乗れない)
    SigmaBelowBar,
    /// cross-talk がバー超過
    CrosstalkExcess,
    /// cross-talk 区間がバーを跨ぐ (Abstain — 強制判定の禁止)
    CrosstalkMarginStraddled,
    /// 合成路の残差がバー超過 (合成は機械実行で検証される)
    SynthesisResidualExcess,
    /// Lie 閉包に標的が入らない — 合成路が存在しない
    NoSynthesisPath,
    /// 状態族が情報完全でない (トモグラフィの資格なし)
    NotInformationallyComplete,
    /// トモグラフィ残差がバー超過 (データと再構成の不整合)
    TomographyResidualExcess,
    /// 再構成が effect の資格 (0 ≤ E ≤ I) を通らない
    EffectQualificationFailed,
    /// 提案 joint POVM に負固有値 (この構成では共同測定の資格が立たない)
    JointCandidateNotPositive,
    /// 提案 joint POVM の総和が I でない
    JointSumNotIdentity,
    /// marginal が対象 effect を再現しない
    MarginalMismatch,
    /// 凸重みが不正 (負・和 ≠ 1)
    WeightsInvalid,
    /// 凸結合が標的状態を再現しない (凸可達性なし)
    MixtureMismatch,
    /// drift の時間変動がバー超過 (単一 regime を名乗れない)
    DriftRegimeUnstable,
    /// 証明書の対象 sha256 が行列と一致しない (証明書の流用禁止)
    CertificateTargetMismatch,
}

impl InterfaceRejection {
    pub fn as_str(self) -> &'static str {
        match self {
            InterfaceRejection::DegenerateTargets => "degenerate_targets",
            InterfaceRejection::CommandTargetCountMismatch => "command_target_count_mismatch",
            InterfaceRejection::InsufficientCommandRank => "insufficient_command_rank",
            InterfaceRejection::SigmaBelowBar => "sigma_below_bar",
            InterfaceRejection::CrosstalkExcess => "crosstalk_excess",
            InterfaceRejection::CrosstalkMarginStraddled => "crosstalk_margin_straddled",
            InterfaceRejection::SynthesisResidualExcess => "synthesis_residual_excess",
            InterfaceRejection::NoSynthesisPath => "no_synthesis_path",
            InterfaceRejection::NotInformationallyComplete => "not_informationally_complete",
            InterfaceRejection::TomographyResidualExcess => "tomography_residual_excess",
            InterfaceRejection::EffectQualificationFailed => "effect_qualification_failed",
            InterfaceRejection::JointCandidateNotPositive => "joint_candidate_not_positive",
            InterfaceRejection::JointSumNotIdentity => "joint_sum_not_identity",
            InterfaceRejection::MarginalMismatch => "marginal_mismatch",
            InterfaceRejection::WeightsInvalid => "weights_invalid",
            InterfaceRejection::MixtureMismatch => "mixture_mismatch",
            InterfaceRejection::DriftRegimeUnstable => "drift_regime_unstable",
            InterfaceRejection::CertificateTargetMismatch => "certificate_target_mismatch",
        }
    }
}

// ---------------------------------------------------------------- 宣言 (資格なし)

/// 意図する役割 — 型タグのみ (装置の実在を主張しない)
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RoleIntent {
    Control,
    Measurement,
    Preparation,
    Drift,
}

/// **宣言された操作 — 資格なし**。行列と意図だけを持ち、AccessibleOperation への
/// 直接変換は存在しない (禁止変換 14)。通過できる門は較正・合成・トモグラフィの
/// 3 証明書のみ。
#[derive(Clone, Debug)]
pub struct DeclaredOperation {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub dim: usize,
    pub intent: RoleIntent,
}

impl DeclaredOperation {
    pub fn matrix_c64(&self) -> Vec<C64> {
        (0..self.dim * self.dim)
            .map(|k| C64::new(self.re[k], self.im[k]))
            .collect()
    }
}

// ---------------------------------------------------------------- 較正: 独立 addressability

/// 証明つき rank — σ_r ≥ σ_bar と σ_{r+1} ≤ dust の**ギャップ**が立つときのみ構成
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CertifiedRank {
    pub rank: usize,
    pub sigma_r: f64,
    pub sigma_r_plus_1: f64,
}

/// 独立 addressability の証明書 (v33.2) — command Jacobian の rank・最小特異値・
/// cross-talk 区間・較正記録の sha256。対象標的族の sha256 に結束される。
#[derive(Clone, Debug)]
pub struct IndependentAddressabilityCertificate {
    pub command_jacobian_rank: CertifiedRank,
    pub smallest_singular_value: BoundInterval,
    pub crosstalk_bound: BoundInterval,
    pub calibration_hash: String,
    target_hashes: Vec<String>,
}

impl IndependentAddressabilityCertificate {
    pub fn covers_target(&self, m: &[C64], n: usize) -> bool {
        self.target_hashes.contains(&matrix_hash(m, n))
    }
}

/// 較正: commands[k] が targets[k] を駆動する 1:1 の割り当てを審査する。
///   - 標的族を HS 正規直交化 (退化は拒否)
///   - M̂_{ik} = ⟨T̂_i, Ĝ_k⟩ の特異値: 全 σ ≥ σ_bar で rank = 標的数 (ギャップ証明)
///   - cross-talk = max_k ‖Ĝ_k − ⟨T̂_k, Ĝ_k⟩T̂_k‖ (自分の標的以外への漏れ) の区間
/// **同一 command を複数標的に割り当てる「数学的分解」は rank 不足で拒否される**
/// (禁止変換 15 の執行点)。
pub fn certify_addressability(
    targets: &[Vec<C64>],
    commands: &[Vec<C64>],
    n: usize,
    sigma_bar: f64,
    crosstalk_bar: f64,
) -> Result<IndependentAddressabilityCertificate, InterfaceRejection> {
    let r = targets.len();
    if commands.len() != r {
        return Err(InterfaceRejection::CommandTargetCountMismatch);
    }
    // 標的の正規直交化 (退化検出)
    let mut tbasis: Vec<Vec<C64>> = Vec::new();
    for t in targets {
        if !push_ortho(&mut tbasis, t, 1e-9) {
            return Err(InterfaceRejection::DegenerateTargets);
        }
    }
    // 正規化 command と M̂ (r×r 複素)
    let ghat: Vec<Vec<C64>> = commands
        .iter()
        .map(|g| {
            let nn = hs_norm(g).max(1e-300);
            g.iter().map(|c| c.scale(1.0 / nn)).collect()
        })
        .collect();
    let that: Vec<Vec<C64>> = targets
        .iter()
        .map(|t| {
            let nn = hs_norm(t).max(1e-300);
            t.iter().map(|c| c.scale(1.0 / nn)).collect()
        })
        .collect();
    let mut m = vec![C64::new(0.0, 0.0); r * r];
    for i in 0..r {
        for k in 0..r {
            m[i * r + k] = hs_inner(&that[i], &ghat[k]);
        }
    }
    // 特異値 = sqrt(eig(M M†))
    let mut mmh = vec![C64::new(0.0, 0.0); r * r];
    for i in 0..r {
        for j in 0..r {
            let mut s = C64::new(0.0, 0.0);
            for k in 0..r {
                let a = m[i * r + k];
                let b = m[j * r + k]; // conj で内積
                s = s + C64::new(a.re * b.re + a.im * b.im, a.im * b.re - a.re * b.im);
            }
            mmh[i * r + j] = s;
        }
    }
    let evals = herm_evals_c64(&mmh, r);
    let mut sigmas: Vec<f64> = evals.iter().map(|&e| e.max(0.0).sqrt()).collect();
    sigmas.sort_by(|a, b| b.partial_cmp(a).unwrap());
    // 実埋め込みは各固有値をちょうど 2 重化する — 降順の偶数位置が真の多重集合
    // (dedup では真に縮退した特異値まで併合してしまう)
    let sigmas: Vec<f64> = sigmas.iter().step_by(2).cloned().collect();
    let sigma_min = sigmas.last().cloned().unwrap_or(0.0);
    let rank_above = sigmas.iter().filter(|&&s| s >= sigma_bar).count();
    if rank_above < r {
        // ギャップの下側 σ が dust (1e-9) を超えるなら rank 不足・以下なら退化寄り —
        // どちらも独立 addressability は立たない
        return Err(InterfaceRejection::InsufficientCommandRank);
    }
    if sigma_min < sigma_bar {
        return Err(InterfaceRejection::SigmaBelowBar);
    }
    // cross-talk: 自分の標的成分を除いた残り
    let mut xtalk = 0.0f64;
    for k in 0..r {
        let c = hs_inner(&that[k], &ghat[k]);
        let resid: Vec<C64> = ghat[k]
            .iter()
            .zip(that[k].iter())
            .map(|(g, t)| *g - c * *t)
            .collect();
        xtalk = xtalk.max(hs_norm(&resid));
    }
    let xt = BoundInterval::new((xtalk - 1e-12).max(0.0), xtalk + 1e-12).unwrap();
    match xt.within(crosstalk_bar) {
        Ok(true) => {}
        Ok(false) => return Err(InterfaceRejection::CrosstalkExcess),
        Err(_) => return Err(InterfaceRejection::CrosstalkMarginStraddled),
    }
    let sv = BoundInterval::new((sigma_min - 1e-12).max(0.0), sigma_min + 1e-12).unwrap();
    let mut record = format!("addressability;n={};r={};", n, r);
    for x in &m {
        record.push_str(&format!("{:.12e},{:.12e};", x.re, x.im));
    }
    Ok(IndependentAddressabilityCertificate {
        command_jacobian_rank: CertifiedRank {
            rank: rank_above,
            sigma_r: sigma_min,
            sigma_r_plus_1: 0.0,
        },
        smallest_singular_value: sv,
        crosstalk_bound: xt,
        calibration_hash: sha256_hex(record.as_bytes()),
        target_hashes: targets.iter().map(|t| matrix_hash(t, n)).collect(),
    })
}

// ---------------------------------------------------------------- 合成: bracket/線形結合の列

/// 合成の一手 — 機械実行で検証される (「レシピの主張」ではなく実行結果が資格)
#[derive(Clone, Debug)]
pub enum SynthStep {
    /// (1/i)[a, b] — 既得 op (workspace index) からエルミートを作る
    BracketOverI(usize, usize),
    /// Σ c_j · op_j (実係数 — エルミート保存)
    Linear(Vec<(f64, usize)>),
}

/// 合成証明書 — base (較正済み生成子) から steps を機械実行し、最終行列が標的に
/// 一致することの記録。対象行列の sha256 に結束。
#[derive(Clone, Debug)]
pub struct SynthesisCertificate {
    pub depth: usize,
    pub residual: BoundInterval,
    pub synthesis_hash: String,
    target_hash: String,
}

impl SynthesisCertificate {
    pub fn covers_target(&self, m: &[C64], n: usize) -> bool {
        self.target_hash == matrix_hash(m, n)
    }
}

pub fn certify_synthesis(
    base: &[Vec<C64>],
    steps: &[SynthStep],
    target: &[C64],
    n: usize,
    err_bar: f64,
) -> Result<SynthesisCertificate, InterfaceRejection> {
    let mut ws: Vec<Vec<C64>> = base.to_vec();
    for st in steps {
        let next = match st {
            SynthStep::BracketOverI(a, b) => {
                let c = commutator(&ws[*a], &ws[*b], n);
                // (1/i)·c = −i·c
                c.iter().map(|x| C64::new(x.im, -x.re)).collect::<Vec<C64>>()
            }
            SynthStep::Linear(terms) => {
                let mut m = vec![C64::new(0.0, 0.0); n * n];
                for (coef, idx) in terms {
                    for (mi, oi) in m.iter_mut().zip(ws[*idx].iter()) {
                        *mi = *mi + oi.scale(*coef);
                    }
                }
                m
            }
        };
        ws.push(next);
    }
    let last = ws.last().cloned().unwrap_or_else(|| vec![C64::new(0.0, 0.0); n * n]);
    let resid: f64 = last
        .iter()
        .zip(target.iter())
        .map(|(a, b)| (*a - *b).norm2())
        .sum::<f64>()
        .sqrt();
    if resid > err_bar {
        return Err(InterfaceRejection::SynthesisResidualExcess);
    }
    let mut record = format!("synthesis;n={};depth={};", n, steps.len());
    for st in steps {
        match st {
            SynthStep::BracketOverI(a, b) => record.push_str(&format!("B{},{};", a, b)),
            SynthStep::Linear(terms) => {
                record.push('L');
                for (c, i) in terms {
                    record.push_str(&format!("{:.12e}*{};", c, i));
                }
            }
        }
    }
    Ok(SynthesisCertificate {
        depth: steps.len(),
        residual: BoundInterval::new((resid - 1e-12).max(0.0), resid + 1e-12).unwrap(),
        synthesis_hash: sha256_hex(record.as_bytes()),
        target_hash: matrix_hash(target, n),
    })
}

/// (1/i)-bracket と線形結合で閉じる実 Lie 閉包の ONB (合成路の存在判定に使う)
pub fn lie_closure(gens: &[Vec<C64>], n: usize) -> Vec<Vec<C64>> {
    let mut basis: Vec<Vec<C64>> = Vec::new();
    for g in gens {
        push_ortho(&mut basis, g, 1e-9);
    }
    loop {
        let snapshot = basis.clone();
        let mut grew = false;
        for a in snapshot.iter() {
            for b in snapshot.iter() {
                let c = commutator(a, b, n);
                let h: Vec<C64> = c.iter().map(|x| C64::new(x.im, -x.re)).collect();
                if hs_norm(&h) > 1e-9 && push_ortho(&mut basis, &h, 1e-9) {
                    grew = true;
                }
            }
        }
        if !grew {
            break;
        }
    }
    basis
}

/// 標的が Lie 閉包 (合成可能域) に入るかの残差 — 入らなければ NoSynthesisPath
pub fn synthesis_path_residual(base: &[Vec<C64>], target: &[C64], n: usize) -> f64 {
    let cl = lie_closure(base, n);
    let mut v = target.to_vec();
    for b in &cl {
        let c = hs_inner(b, &v);
        for (vi, bi) in v.iter_mut().zip(b.iter()) {
            *vi = *vi - c * *bi;
        }
    }
    hs_norm(&v) / hs_norm(target).max(1e-300)
}

// ---------------------------------------------------------------- トモグラフィ (測定側の出自)

/// トモグラフィ証明書 — 情報完全状態族からの線形再構成の残差記録。
/// 再構成された effect の sha256 に結束。
#[derive(Clone, Debug)]
pub struct TomographyCertificate {
    pub residual: BoundInterval,
    pub frame_sigma_min: f64,
    pub data_hash: String,
    target_hash: String,
}

impl TomographyCertificate {
    pub fn covers_target(&self, m: &[C64], n: usize) -> bool {
        self.target_hash == matrix_hash(m, n)
    }
}

/// effect のトモグラフィ: 状態族 ρ_j と観測確率 p_j から最小二乗で Ê を再構成。
///   - 状態族が情報完全 (design 行列の σ_min > 1e-9) でなければ拒否
///   - LS 残差 (max_j |Tr(ρ_j Ê) − p_j|) がバーを超えれば拒否
///   - Ê が effect の資格 (0 ≤ Ê ≤ I) を通らなければ拒否
pub fn certify_effect_tomography(
    states: &[Vec<C64>],
    observed: &[f64],
    n: usize,
    residual_bar: f64,
) -> Result<(MeasurementEffect, TomographyCertificate), InterfaceRejection> {
    let m = states.len();
    assert_eq!(observed.len(), m, "状態数と観測数が一致しない");
    // エルミート ONB (n² 本): 対角 E_aa・(E_ab+E_ba)/√2・i(E_ab−E_ba)/√2
    let mut hbasis: Vec<Vec<C64>> = Vec::new();
    for a in 0..n {
        let mut e = vec![C64::new(0.0, 0.0); n * n];
        e[a * n + a] = C64::new(1.0, 0.0);
        hbasis.push(e);
    }
    let inv = 1.0 / (2.0f64).sqrt();
    for a in 0..n {
        for b in (a + 1)..n {
            let mut e = vec![C64::new(0.0, 0.0); n * n];
            e[a * n + b] = C64::new(inv, 0.0);
            e[b * n + a] = C64::new(inv, 0.0);
            hbasis.push(e);
            let mut f = vec![C64::new(0.0, 0.0); n * n];
            f[a * n + b] = C64::new(0.0, inv);
            f[b * n + a] = C64::new(0.0, -inv);
            hbasis.push(f);
        }
    }
    let d = hbasis.len(); // n²
    // design 行列 A_{jk} = Tr(ρ_j B_k) (実 — 両者エルミート)
    let mut a_mat = vec![0.0f64; m * d];
    for j in 0..m {
        for k in 0..d {
            a_mat[j * d + k] = hs_inner(&hbasis[k], &states[j]).re;
        }
    }
    // 正規方程式 G = AᵀA・b = Aᵀp
    let mut gram = vec![0.0f64; d * d];
    let mut rhs = vec![0.0f64; d];
    for j in 0..m {
        for k in 0..d {
            rhs[k] += a_mat[j * d + k] * observed[j];
            for l in 0..d {
                gram[k * d + l] += a_mat[j * d + k] * a_mat[j * d + l];
            }
        }
    }
    let (evals, vecs) = jacobi_eigh(&gram, d);
    let emax = evals.iter().cloned().fold(0.0f64, f64::max).max(1e-300);
    let emin = evals.iter().cloned().fold(f64::INFINITY, f64::min);
    if emin <= 1e-9 * emax {
        return Err(InterfaceRejection::NotInformationallyComplete);
    }
    // 解 x = G⁻¹ b (固有分解経由)
    let mut x = vec![0.0f64; d];
    for k in 0..d {
        let mut proj = 0.0;
        for t in 0..d {
            proj += vecs[t + k * d] * rhs[t];
        }
        let coef = proj / evals[k];
        for t in 0..d {
            x[t] += coef * vecs[t + k * d];
        }
    }
    let mut ehat = vec![C64::new(0.0, 0.0); n * n];
    for k in 0..d {
        for (ei, bi) in ehat.iter_mut().zip(hbasis[k].iter()) {
            *ei = *ei + bi.scale(x[k]);
        }
    }
    // 残差
    let mut resid = 0.0f64;
    for j in 0..m {
        let pred = hs_inner(&ehat, &states[j]).re;
        resid = resid.max((pred - observed[j]).abs());
    }
    if resid > residual_bar {
        return Err(InterfaceRejection::TomographyResidualExcess);
    }
    let eff = MeasurementEffect::certify(
        ehat.iter().map(|c| c.re).collect(),
        ehat.iter().map(|c| c.im).collect(),
        n,
    )
    .map_err(|_| InterfaceRejection::EffectQualificationFailed)?;
    let mut record = format!("tomography;n={};m={};", n, m);
    for p in observed {
        record.push_str(&format!("{:.12e};", p));
    }
    let cert = TomographyCertificate {
        residual: BoundInterval::new((resid - 1e-12).max(0.0), resid + 1e-12).unwrap(),
        frame_sigma_min: (emin / emax).sqrt(),
        data_hash: sha256_hex(record.as_bytes()),
        target_hash: matrix_hash(&ehat, n),
    };
    Ok((eff, cert))
}

// ---------------------------------------------------------------- 出自 (3 つの門)

/// 操作の出自 — **文字列 provenance の廃止**: 出自は型と sha256 結束の証明書で運ぶ
#[derive(Clone, Debug)]
pub enum OperationOrigin {
    DirectlyCalibrated(IndependentAddressabilityCertificate),
    Synthesized(SynthesisCertificate),
    TomographicallyInferred(TomographyCertificate),
}

impl OperationOrigin {
    pub fn as_str(&self) -> &'static str {
        match self {
            OperationOrigin::DirectlyCalibrated(_) => "directly_calibrated",
            OperationOrigin::Synthesized(_) => "synthesized",
            OperationOrigin::TomographicallyInferred(_) => "tomographically_inferred",
        }
    }
    fn covers_target(&self, m: &[C64], n: usize) -> bool {
        match self {
            OperationOrigin::DirectlyCalibrated(c) => c.covers_target(m, n),
            OperationOrigin::Synthesized(c) => c.covers_target(m, n),
            OperationOrigin::TomographicallyInferred(c) => c.covers_target(m, n),
        }
    }
}

// ---------------------------------------------------------------- 資源予算 (成分半順序)

/// 資源予算 — **成分ごとの半順序のみ** (Ord は実装しない: 恣意的な重み付き和での
/// 全順序化は新しい選択バイアスになる — PROMPT/14)。v33.3 の resource-indexed
/// factorization profile の土台。
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResourceBudget {
    pub max_duration: f64,
    pub max_amplitude: f64,
    pub max_bandwidth: f64,
    pub max_depth: f64,
    pub max_error: f64,
}

impl ResourceBudget {
    pub fn certify(
        max_duration: f64,
        max_amplitude: f64,
        max_bandwidth: f64,
        max_depth: f64,
        max_error: f64,
    ) -> Result<Self, &'static str> {
        let vals = [max_duration, max_amplitude, max_bandwidth, max_depth, max_error];
        if vals.iter().any(|v| !v.is_finite() || *v < 0.0) {
            return Err("資源予算は非負有限であること");
        }
        Ok(ResourceBudget {
            max_duration,
            max_amplitude,
            max_bandwidth,
            max_depth,
            max_error,
        })
    }
    /// 成分ごとの ≤ (半順序) — a ≤ b でも b ≤ a でもない比較不能対が存在する
    pub fn componentwise_le(&self, other: &ResourceBudget) -> bool {
        self.max_duration <= other.max_duration
            && self.max_amplitude <= other.max_amplitude
            && self.max_bandwidth <= other.max_bandwidth
            && self.max_depth <= other.max_depth
            && self.max_error <= other.max_error
    }
    /// 比較可能か (どちらか向きの ≤ が立つか)
    pub fn comparable(&self, other: &ResourceBudget) -> bool {
        self.componentwise_le(other) || other.componentwise_le(self)
    }
}

// ---------------------------------------------------------------- AccessibleOperation

/// **資格つき操作** — 役割 (数学的資格審査済みの OpKind)・出自 (3 門のいずれか)・
/// 独立 addressability・資源予算。構成は certify のみで、出自証明書の対象 sha256
/// が行列と一致しない場合 (証明書の流用) は構成時拒否。
#[derive(Clone, Debug)]
pub struct AccessibleOperation {
    kind: OpKind,
    parity: OperatorParity,
    origin: OperationOrigin,
    addressability: IndependentAddressabilityCertificate,
    budget: ResourceBudget,
}

impl AccessibleOperation {
    pub fn certify(
        kind: OpKind,
        parity: OperatorParity,
        origin: OperationOrigin,
        addressability: IndependentAddressabilityCertificate,
        budget: ResourceBudget,
    ) -> Result<Self, InterfaceRejection> {
        let (re, im, d) = kind.matrix();
        let m: Vec<C64> = (0..d * d).map(|k| C64::new(re[k], im[k])).collect();
        if !origin.covers_target(&m, d) {
            return Err(InterfaceRejection::CertificateTargetMismatch);
        }
        Ok(AccessibleOperation {
            kind,
            parity,
            origin,
            addressability,
            budget,
        })
    }
    pub fn origin(&self) -> &OperationOrigin {
        &self.origin
    }
    pub fn budget(&self) -> &ResourceBudget {
        &self.budget
    }
    pub fn addressability(&self) -> &IndependentAddressabilityCertificate {
        &self.addressability
    }
}

// ---------------------------------------------------------------- role-typed 文脈 4 型

/// joint measurability の証明書 — 明示的 joint POVM の機械検証 (**可換性より広い**:
/// 非可換な unsharp 対でも joint POVM が実在すれば資格が立つ)。CertifiedCommutator
/// からの変換は存在しない (禁止変換 16)。
#[derive(Clone, Debug)]
pub struct JointMeasurementCertificate {
    pub psd_min: f64,
    pub marginal_residual: f64,
    pub n_outcomes: usize,
}

/// 提案 joint POVM {G_j} と分割 marginals[i] ⊆ outcomes を検証する:
/// 各 G_j ⪰ 0・Σ G_j = I・Σ_{j∈marginals[i]} G_j = effects[i]。
pub fn certify_joint_measurement(
    effects: &[Vec<C64>],
    joint: &[Vec<C64>],
    marginals: &[Vec<usize>],
    n: usize,
) -> Result<JointMeasurementCertificate, InterfaceRejection> {
    assert_eq!(effects.len(), marginals.len(), "effect と marginal 分割の数が不一致");
    let mut psd_min = f64::INFINITY;
    for g in joint {
        let ev = herm_evals_c64(g, n);
        let mn = ev.iter().cloned().fold(f64::INFINITY, f64::min);
        psd_min = psd_min.min(mn);
    }
    if psd_min < -1e-12 {
        return Err(InterfaceRejection::JointCandidateNotPositive);
    }
    let mut sum = vec![C64::new(0.0, 0.0); n * n];
    for g in joint {
        for (s, x) in sum.iter_mut().zip(g.iter()) {
            *s = *s + *x;
        }
    }
    for i in 0..n {
        sum[i * n + i] = sum[i * n + i] - C64::new(1.0, 0.0);
    }
    if hs_norm(&sum) > 1e-12 {
        return Err(InterfaceRejection::JointSumNotIdentity);
    }
    let mut marg_resid = 0.0f64;
    for (e, idxs) in effects.iter().zip(marginals.iter()) {
        let mut acc = vec![C64::new(0.0, 0.0); n * n];
        for &j in idxs {
            for (a, x) in acc.iter_mut().zip(joint[j].iter()) {
                *a = *a + *x;
            }
        }
        let d: f64 = acc
            .iter()
            .zip(e.iter())
            .map(|(a, b)| (*a - *b).norm2())
            .sum::<f64>()
            .sqrt();
        marg_resid = marg_resid.max(d);
    }
    if marg_resid > 1e-12 {
        return Err(InterfaceRejection::MarginalMismatch);
    }
    Ok(JointMeasurementCertificate {
        psd_min,
        marginal_residual: marg_resid,
        n_outcomes: joint.len(),
    })
}

/// 凸可達性の証明書 — 校正済み準備族の凸結合が標的を再現することの機械検証
#[derive(Clone, Debug)]
pub struct ConvexReachabilityCertificate {
    pub weights: Vec<f64>,
    pub residual: f64,
}

pub fn certify_convex_reachability(
    preps: &[Preparation],
    weights: &[f64],
    target: &Preparation,
) -> Result<ConvexReachabilityCertificate, InterfaceRejection> {
    assert_eq!(preps.len(), weights.len(), "準備と重みの数が不一致");
    if weights.iter().any(|&w| !w.is_finite() || w < -1e-12)
        || (weights.iter().sum::<f64>() - 1.0).abs() > 1e-10
    {
        return Err(InterfaceRejection::WeightsInvalid);
    }
    let n = target.dim;
    let mut mix_re = vec![0.0f64; n * n];
    let mut mix_im = vec![0.0f64; n * n];
    for (p, &w) in preps.iter().zip(weights.iter()) {
        for k in 0..n * n {
            mix_re[k] += w * p.re[k];
            mix_im[k] += w * p.im[k];
        }
    }
    let mut resid = 0.0f64;
    for k in 0..n * n {
        resid += (mix_re[k] - target.re[k]).powi(2) + (mix_im[k] - target.im[k]).powi(2);
    }
    let resid = resid.sqrt();
    if resid > 1e-10 {
        return Err(InterfaceRejection::MixtureMismatch);
    }
    Ok(ConvexReachabilityCertificate {
        weights: weights.to_vec(),
        residual: resid,
    })
}

/// drift regime の安定性証明書 — 時間窓内のスナップショット間変動のバー検査
#[derive(Clone, Debug)]
pub struct StabilityCertificate {
    pub max_variation: f64,
    pub bar: f64,
}

pub fn certify_drift_stability(
    snapshots: &[DriftGenerator],
    variation_bar: f64,
) -> Result<StabilityCertificate, InterfaceRejection> {
    let mut var = 0.0f64;
    for (i, a) in snapshots.iter().enumerate() {
        for b in snapshots.iter().skip(i + 1) {
            let n = a.dim;
            let mut d = 0.0f64;
            for k in 0..n * n {
                d += (a.re[k] - b.re[k]).powi(2) + (a.im[k] - b.im[k]).powi(2);
            }
            var = var.max(d.sqrt());
        }
    }
    if var > variation_bar {
        return Err(InterfaceRejection::DriftRegimeUnstable);
    }
    Ok(StabilityCertificate {
        max_variation: var,
        bar: variation_bar,
    })
}

/// 制御文脈 — 独立 addressability 証明書を要求する (可換性だけでは構成できない)
pub struct ControlContext {
    pub members: Vec<OpId>,
    pub addressability: IndependentAddressabilityCertificate,
}

/// 測定文脈 — joint measurability 証明書を要求する (可換性より広い)
pub struct MeasurementContext {
    pub effects: Vec<MeasurementEffect>,
    pub joint: JointMeasurementCertificate,
}

/// 準備族 — 凸可達性証明書を要求する
pub struct PreparationFamily {
    pub preps: Vec<Preparation>,
    pub reach: ConvexReachabilityCertificate,
}

/// drift regime — 安定性証明書を要求する
pub struct DriftRegime {
    pub generator: DriftGenerator,
    pub stability: StabilityCertificate,
}

// ---------------------------------------------------------------- AccessibleOperationalNet

/// 資格つき操作だけで構成される operational net (v33.2)。
/// **DeclaredOperation を受け付ける口は存在しない** (禁止変換 14) — 入口は
/// admit(AccessibleOperation) のみ。制御文脈は独立 addressability 証明書つきで
/// 登録され、下層の OperationalNet (v32.2/v33.1 契約) に降ろされる。
/// 測定文脈・準備族・drift regime は型つきで保持される (因子分解への統合は
/// v33.4 の contextual factorization)。
pub struct AccessibleOperationalNet<G: CommutationGrading> {
    net: OperationalNet<G>,
    origins: Vec<OperationOrigin>,
    control_contexts: Vec<ControlContext>,
    measurement_contexts: Vec<MeasurementContext>,
    preparation_families: Vec<PreparationFamily>,
    drift_regimes: Vec<DriftRegime>,
}

impl<G: CommutationGrading> AccessibleOperationalNet<G> {
    pub fn new(dim: usize, threshold: f64) -> Self {
        AccessibleOperationalNet {
            net: OperationalNet::new(dim, threshold),
            origins: Vec::new(),
            control_contexts: Vec::new(),
            measurement_contexts: Vec::new(),
            preparation_families: Vec::new(),
            drift_regimes: Vec::new(),
        }
    }

    /// 唯一の操作入口 — AccessibleOperation (資格つき) のみ受け付ける
    pub fn admit(&mut self, op: AccessibleOperation) -> Result<OpId, &'static str> {
        let AccessibleOperation {
            kind,
            parity,
            origin,
            ..
        } = op;
        let id = self.net.add_primitive(PrimitiveOperation {
            kind,
            parity,
            // 旧 v32.2 契約の文字列フィールドは定数で中立化 — 新 lane の出自は
            // origins (型 + sha256) が運ぶ
            provenance: "certified_accessible_operation",
        })?;
        self.origins.push(origin);
        Ok(id)
    }

    pub fn set_commutator(
        &mut self,
        a: OpId,
        b: OpId,
        c: crate::operational_net::CertifiedCommutator,
    ) {
        self.net.set_commutator(a, b, c);
    }

    /// 制御文脈の登録 — 独立 addressability 証明書が member 全員を結束していること
    /// + 下層 net の可換子証明書要求 (v32.2) の両方を通す
    pub fn add_control_context(
        &mut self,
        members: &[OpId],
        addressability: IndependentAddressabilityCertificate,
    ) -> Result<usize, String> {
        for id in members {
            let (re, im, d) = self.net.primitive(*id).kind.matrix();
            let m: Vec<C64> = (0..d * d).map(|k| C64::new(re[k], im[k])).collect();
            if !addressability.covers_target(&m, d) {
                return Err(format!(
                    "証明書が member {:?} を結束していない (証明書の流用禁止)",
                    id
                ));
            }
        }
        self.net.add_context(members)?;
        self.control_contexts.push(ControlContext {
            members: members.to_vec(),
            addressability,
        });
        Ok(self.control_contexts.len() - 1)
    }

    pub fn attach_measurement_context(&mut self, ctx: MeasurementContext) -> usize {
        self.measurement_contexts.push(ctx);
        self.measurement_contexts.len() - 1
    }
    pub fn attach_preparation_family(&mut self, fam: PreparationFamily) -> usize {
        self.preparation_families.push(fam);
        self.preparation_families.len() - 1
    }
    pub fn attach_drift_regime(&mut self, reg: DriftRegime) -> usize {
        self.drift_regimes.push(reg);
        self.drift_regimes.len() - 1
    }

    pub fn origin(&self, id: OpId) -> &OperationOrigin {
        &self.origins[id.0 as usize]
    }
    pub fn n_control_contexts(&self) -> usize {
        self.control_contexts.len()
    }
    pub fn n_measurement_contexts(&self) -> usize {
        self.measurement_contexts.len()
    }

    /// v33.1 の修復入口 (contexts 必須・role 純度・被覆) を通した因子分解復元
    pub fn recover(&self) -> Result<MarkedRecoveryDetail, RecoveryInputRejection> {
        Ok(self.net.recovery_input()?.recover())
    }

    pub fn inner_net(&self) -> &OperationalNet<G> {
        &self.net
    }
}

// ---------------------------------------------------------------- 自己検査

/// laboratory_interface の不変条件 (v332_certified_interface が呼ぶ)
pub fn laboratory_interface_self_test() -> Result<(), String> {
    let n = 4usize;
    let px = [C64::new(0.0, 0.0), C64::new(1.0, 0.0), C64::new(1.0, 0.0), C64::new(0.0, 0.0)];
    let pz = [
        C64::new(1.0, 0.0),
        C64::new(0.0, 0.0),
        C64::new(0.0, 0.0),
        C64::new(-1.0, 0.0),
    ];
    let id2 = [C64::new(1.0, 0.0), C64::new(0.0, 0.0), C64::new(0.0, 0.0), C64::new(1.0, 0.0)];
    let kron2 = |a: &[C64], b: &[C64]| -> Vec<C64> {
        let mut out = vec![C64::new(0.0, 0.0); 16];
        for i1 in 0..2 {
            for j1 in 0..2 {
                for i2 in 0..2 {
                    for j2 in 0..2 {
                        out[(i1 * 2 + i2) * 4 + (j1 * 2 + j2)] = a[i1 * 2 + j1] * b[i2 * 2 + j2];
                    }
                }
            }
        }
        out
    };
    let x1 = kron2(&px, &id2);
    let x2 = kron2(&id2, &px);
    let z2 = kron2(&id2, &pz);
    let tied: Vec<C64> = x1.iter().zip(x2.iter()).map(|(a, b)| *a + *b).collect();
    // 1. 独立 knobs は資格・数学的分解 (同一 command × 2 標的) は rank 不足で拒否
    let ok = certify_addressability(
        &[x1.clone(), x2.clone()],
        &[x1.clone(), x2.clone()],
        n,
        0.5,
        0.1,
    );
    if ok.is_err() {
        return Err("独立 knobs の addressability が立たない".into());
    }
    let bad = certify_addressability(
        &[x1.clone(), x2.clone()],
        &[tied.clone(), tied.clone()],
        n,
        0.5,
        0.1,
    );
    match bad {
        Err(InterfaceRejection::InsufficientCommandRank) => {}
        Err(e) => return Err(format!("数学的分解の拒否理由が違う: {:?}", e)),
        Ok(_) => return Err("数学的分解が拒否されない".into()),
    }
    // 2. 合成: {X₁+X₂, Z₂} から X₁ (bracket 3 手 + 線形 1 手)・{X₁+X₂} 単独は路なし
    let base = vec![tied.clone(), z2.clone()];
    let steps = vec![
        SynthStep::BracketOverI(0, 1),          // (1/i)[X₁+X₂, Z₂] = −2Y₂ → idx 2
        SynthStep::BracketOverI(2, 1),          // (1/i)[−2Y₂, Z₂] = −4X₂ → idx 3
        SynthStep::Linear(vec![(1.0, 0), (0.25, 3)]), // (X₁+X₂) + 0.25·(−4X₂) = X₁
    ];
    let synth = certify_synthesis(&base, &steps, &x1, n, 1e-9);
    if synth.is_err() {
        return Err(format!("合成路の証明が立たない: {:?}", synth.err()));
    }
    if synthesis_path_residual(&[tied.clone()], &x1, n) < 0.1 {
        return Err("単独 tied command から X₁ への路が偽装できた".into());
    }
    // 3. 証明書の流用禁止 (X₁ の合成証明書を X₂ に付けられない)
    let cert = synth.unwrap();
    let mk_ctrl = |g: &[C64]| {
        OpKind::Control(
            ControlGenerator::certify(
                g.iter().map(|c| c.re).collect(),
                g.iter().map(|c| c.im).collect(),
                n,
            )
            .unwrap(),
        )
    };
    let addr = certify_addressability(&[x1.clone()], &[x1.clone()], n, 0.5, 0.1).unwrap();
    let budget = ResourceBudget::certify(1.0, 1.0, 1.0, 3.0, 1e-9).unwrap();
    match AccessibleOperation::certify(
        mk_ctrl(&x2),
        OperatorParity::Even,
        OperationOrigin::Synthesized(cert.clone()),
        addr.clone(),
        budget,
    ) {
        Err(InterfaceRejection::CertificateTargetMismatch) => {}
        _ => return Err("X₁ の合成証明書が X₂ に流用できた".into()),
    }
    if AccessibleOperation::certify(
        mk_ctrl(&x1),
        OperatorParity::Even,
        OperationOrigin::Synthesized(cert),
        addr,
        budget,
    )
    .is_err()
    {
        return Err("正当な合成 X₁ の受理が拒否された".into());
    }
    // 4. joint POVM: 非可換 unsharp 対 (η = 0.6, 0.6) は資格・(0.8, 0.8) は拒否
    let n2 = 2usize;
    let mk_joint = |eta: f64| -> (Vec<Vec<C64>>, Vec<Vec<C64>>) {
        let e_x: Vec<C64> = (0..4)
            .map(|k| id2[k].scale(0.5) + px[k].scale(0.5 * eta))
            .collect();
        let e_z: Vec<C64> = (0..4)
            .map(|k| id2[k].scale(0.5) + pz[k].scale(0.5 * eta))
            .collect();
        let mut joint = Vec::new();
        for s in [1.0f64, -1.0] {
            for t in [1.0f64, -1.0] {
                let g: Vec<C64> = (0..4)
                    .map(|k| {
                        id2[k].scale(0.25) + px[k].scale(0.25 * s * eta) + pz[k].scale(0.25 * t * eta)
                    })
                    .collect();
                joint.push(g);
            }
        }
        (vec![e_x, e_z], joint)
    };
    let (effects, joint) = mk_joint(0.6);
    let marg = vec![vec![0usize, 1], vec![0usize, 2]];
    if certify_joint_measurement(&effects, &joint, &marg, n2).is_err() {
        return Err("η = 0.6 の joint POVM が資格を通らない".into());
    }
    let (effects8, joint8) = mk_joint(0.8);
    match certify_joint_measurement(&effects8, &joint8, &marg, n2) {
        Err(InterfaceRejection::JointCandidateNotPositive) => {}
        r => return Err(format!("η = 0.8 の joint 候補が拒否されない: {:?}", r.err())),
    }
    // 5. 資源予算は半順序 — 比較不能対の実在
    let b1 = ResourceBudget::certify(1.0, 2.0, 1.0, 1.0, 1e-3).unwrap();
    let b2 = ResourceBudget::certify(2.0, 1.0, 1.0, 1.0, 1e-3).unwrap();
    if b1.comparable(&b2) {
        return Err("成分半順序に比較不能対が無い (全順序化されている)".into());
    }
    Ok(())
}
