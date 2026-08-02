// operational_net — 操作的文脈の型契約 (v32.2, PROMPT/13 §2)
//
// 第三十二期の中心テーゼ「局所性は状態に宿るのではなく、操作的文脈の可換分解と、
// その上の Liouvillian 応答の同値類としてのみ識別される」の型実装。v31.4 の E3
// no-go (状態単独では因子分解を選定できない — 因子分解は操作代数が運ぶ) を受けて、
// 旧 OperationalAlgebra (readout_contract — 文字列リストの雛形) を **marked family
// of subalgebras** に精密化する。
//
// 設計の核 (qrn_core / readout_contract と同型の封鎖):
//   1. **単一の閉じた大域代数は因子分解を運べない** (global-algebra erasure no-go,
//      v322 [N1]): site-local な完全生成族も mode-local な完全生成族も、*-閉包を
//      取れば同じ B(H) を生成し得る — どちらが site かの情報は閉包時に消える。
//      因子分解を運ぶのは「どの primitive が存在し、どれが独立に address 可能で、
//      どの部分集合が一つの操作文脈をなすか」という **marking** である。
//      型: GlobalClosure は代数の同型類 (次元) しか持たず、OperationalNet・
//      因子分解への変換は存在しない (**禁止変換 11** — v322 [N6] が impl From の
//      不在を監査)。
//   2. **準備・介入・測定は同じ数学的型ではない**: Preparation = 状態 (凸族,
//      tr = 1 / PSD を構成時に資格審査)、ControlGenerator = エルミート生成子
//      (Lie 閉包で合成)、MeasurementEffect = 作用素系の元 (0 ≤ E ≤ I — **積閉包を
//      要求しない**)、DriftGenerator = 制御不能な発展 (Control と別役割)。
//      相互の impl From は存在しない。
//   3. **可換性は証明書つき**: CertifiedCommutator は (graded) bracket ノルムの
//      区間 [lo, hi] で持ち、閾値を跨ぐ区間の裁定は Abstain — ノイズ下で辺を
//      強制しない (HOLD-7 棄却原則の継承)。
//   4. **grading は型**: OrdinaryCommutation の net は odd (fermionic) primitive を
//      構成時に拒否する — Jordan–Wigner 弦の見かけの非可換性を幾何と誤読する事故
//      (v322 [N5]: 独立モードの odd 対が ordinary 可換子では完全グラフに見える) を
//      型で遮断する。現行の資格は parity-even lane。
//
// 一次ソース: docs/uft-v32.2.md (no-go と型分離) / core.schema.yml (概念登録 +
// 禁止変換 11)。整合は v322_operational_net が機械検査する。

use crate::{jacobi_eigh, C64};
use std::collections::BTreeSet;
use std::marker::PhantomData;

mod sealed {
    pub trait Sealed {}
}

// ---------------------------------------------------------------- grading (2 タグ)

/// 可換性の等級付け。Ordinary は odd primitive を受け付けない (JW 誤読の遮断)。
pub trait CommutationGrading: sealed::Sealed {
    const NAME: &'static str;
    const ACCEPTS_ODD: bool;
}

/// 通常の可換性 ([A, B]) — parity-even lane 専用
pub enum OrdinaryCommutation {}
/// Z2 graded — odd×odd 対は反可換子 {A, B} で裁定する
pub enum FermionicZ2Graded {}

impl sealed::Sealed for OrdinaryCommutation {}
impl CommutationGrading for OrdinaryCommutation {
    const NAME: &'static str = "ordinary_commutation";
    const ACCEPTS_ODD: bool = false;
}
impl sealed::Sealed for FermionicZ2Graded {}
impl CommutationGrading for FermionicZ2Graded {
    const NAME: &'static str = "fermionic_z2_graded";
    const ACCEPTS_ODD: bool = true;
}

/// フェルミオンパリティ (Z2)。Even = 双線形など・Odd = 単一生成/消滅の線形結合
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OperatorParity {
    Even,
    Odd,
}

// ---------------------------------------------------------------- 役割 4 型 (相互 From なし)

fn herm_check(re: &[f64], im: &[f64], n: usize) -> Result<(), &'static str> {
    if re.len() != n * n || im.len() != n * n {
        return Err("行列の形が n×n でない");
    }
    for i in 0..n {
        if im[i * n + i].abs() > 1e-12 {
            return Err("エルミートなら対角虚部は 0");
        }
        for j in 0..n {
            if (re[i * n + j] - re[j * n + i]).abs() > 1e-12
                || (im[i * n + j] + im[j * n + i]).abs() > 1e-12
            {
                return Err("エルミートでない");
            }
        }
    }
    Ok(())
}

/// エルミート行列の固有値 (2n×2n 実対称埋め込み — 各固有値は 2 重化される)
fn herm_evals(re: &[f64], im: &[f64], n: usize) -> Vec<f64> {
    let m = 2 * n;
    let mut big = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            big[i * m + j] = re[i * n + j];
            big[(i + n) * m + (j + n)] = re[i * n + j];
            big[i * m + (j + n)] = -im[i * n + j];
            big[(i + n) * m + j] = im[i * n + j];
        }
    }
    let (evals, _) = jacobi_eigh(&big, m);
    evals
}

/// 準備 = 密度行列 (凸族の元)。構成は certify のみ — tr = 1・エルミート・PSD。
#[derive(Clone, Debug)]
pub struct Preparation {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub dim: usize,
}

impl Preparation {
    pub fn certify(re: Vec<f64>, im: Vec<f64>, dim: usize) -> Result<Self, &'static str> {
        herm_check(&re, &im, dim)?;
        let tr: f64 = (0..dim).map(|i| re[i * dim + i]).sum();
        if (tr - 1.0).abs() > 1e-10 {
            return Err("trace が 1 でない (状態でない)");
        }
        let evals = herm_evals(&re, &im, dim);
        if evals.iter().any(|&e| e < -1e-10) {
            return Err("PSD でない (状態でない)");
        }
        Ok(Preparation { re, im, dim })
    }
    /// 凸結合 — 準備の合成は混合 (積ではない)
    pub fn mix(a: &Preparation, b: &Preparation, p: f64) -> Result<Preparation, &'static str> {
        if !(0.0..=1.0).contains(&p) || a.dim != b.dim {
            return Err("凸係数または次元が不正");
        }
        let re: Vec<f64> = a
            .re
            .iter()
            .zip(&b.re)
            .map(|(x, y)| p * x + (1.0 - p) * y)
            .collect();
        let im: Vec<f64> = a
            .im
            .iter()
            .zip(&b.im)
            .map(|(x, y)| p * x + (1.0 - p) * y)
            .collect();
        Preparation::certify(re, im, a.dim)
    }
}

/// 介入 = エルミート生成子 (control)。合成は Lie bracket — 積閉包は要求しない。
#[derive(Clone, Debug)]
pub struct ControlGenerator {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub dim: usize,
}

impl ControlGenerator {
    pub fn certify(re: Vec<f64>, im: Vec<f64>, dim: usize) -> Result<Self, &'static str> {
        herm_check(&re, &im, dim)?;
        Ok(ControlGenerator { re, im, dim })
    }
}

/// 測定 = effect (0 ≤ E ≤ I)。作用素系の元 — **積に閉じる必要はない** (v322 [N3])。
#[derive(Clone, Debug)]
pub struct MeasurementEffect {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub dim: usize,
}

impl MeasurementEffect {
    pub fn certify(re: Vec<f64>, im: Vec<f64>, dim: usize) -> Result<Self, &'static str> {
        herm_check(&re, &im, dim)?;
        let evals = herm_evals(&re, &im, dim);
        if evals.iter().any(|&e| e < -1e-10 || e > 1.0 + 1e-10) {
            return Err("effect でない (0 ≤ E ≤ I が破れる)");
        }
        Ok(MeasurementEffect { re, im, dim })
    }
}

/// 制御不能な発展生成子 (drift) — Control と同じ数学だが役割が違う (別型)
#[derive(Clone, Debug)]
pub struct DriftGenerator {
    pub re: Vec<f64>,
    pub im: Vec<f64>,
    pub dim: usize,
}

impl DriftGenerator {
    pub fn certify(re: Vec<f64>, im: Vec<f64>, dim: usize) -> Result<Self, &'static str> {
        herm_check(&re, &im, dim)?;
        Ok(DriftGenerator { re, im, dim })
    }
}

/// primitive の役割 — 数学的型ごと分離 (文字列タグではない)
#[derive(Clone, Debug)]
pub enum OpKind {
    Prepare(Preparation),
    Control(ControlGenerator),
    Measure(MeasurementEffect),
    Drift(DriftGenerator),
}

impl OpKind {
    pub fn role_name(&self) -> &'static str {
        match self {
            OpKind::Prepare(_) => "preparation",
            OpKind::Control(_) => "control",
            OpKind::Measure(_) => "measurement",
            OpKind::Drift(_) => "drift",
        }
    }
    pub fn matrix(&self) -> (&[f64], &[f64], usize) {
        match self {
            OpKind::Prepare(p) => (&p.re, &p.im, p.dim),
            OpKind::Control(c) => (&c.re, &c.im, c.dim),
            OpKind::Measure(m) => (&m.re, &m.im, m.dim),
            OpKind::Drift(d) => (&d.re, &d.im, d.dim),
        }
    }
}

/// primitive operation — 役割・パリティ・出自を運ぶ (marking の最小単位)
#[derive(Clone, Debug)]
pub struct PrimitiveOperation {
    pub kind: OpKind,
    pub parity: OperatorParity,
    pub provenance: &'static str,
}

/// net 内の primitive の識別子
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct OpId(pub u32);

// ---------------------------------------------------------------- 可換子証明書

/// (graded) bracket ノルムの区間証明書。裁定は閾値との 3 値比較 —
/// 区間が閾値を跨ぐ場合は Abstain (辺の強制禁止)。
#[derive(Clone, Copy, Debug)]
pub struct CertifiedCommutator {
    lo: f64,
    hi: f64,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CommutatorVerdict {
    Commuting,
    NonCommuting,
    Abstain,
}

impl CommutatorVerdict {
    pub fn as_str(self) -> &'static str {
        match self {
            CommutatorVerdict::Commuting => "commuting",
            CommutatorVerdict::NonCommuting => "non_commuting",
            CommutatorVerdict::Abstain => "abstain",
        }
    }
}

impl CertifiedCommutator {
    pub fn new(lo: f64, hi: f64) -> Result<Self, &'static str> {
        if !(lo.is_finite() && hi.is_finite()) || lo < 0.0 || hi < lo {
            return Err("区間が不正 (0 ≤ lo ≤ hi の有限値であること)");
        }
        Ok(CertifiedCommutator { lo, hi })
    }
    pub fn lo(&self) -> f64 {
        self.lo
    }
    pub fn hi(&self) -> f64 {
        self.hi
    }
    pub fn verdict(&self, threshold: f64) -> CommutatorVerdict {
        if self.hi < threshold {
            CommutatorVerdict::Commuting
        } else if self.lo > threshold {
            CommutatorVerdict::NonCommuting
        } else {
            CommutatorVerdict::Abstain
        }
    }
}

// ---------------------------------------------------------------- OperationalNet

/// marked family of subalgebras — 「どの操作が存在し、どれが独立に address 可能で、
/// どの部分集合が一つの操作文脈をなすか」。閉包 (GlobalClosure) が消す情報の置き場。
pub struct OperationalNet<G: CommutationGrading> {
    primitives: Vec<PrimitiveOperation>,
    contexts: Vec<BTreeSet<u32>>,
    /// 上三角 (i < j) の平坦格納 — set_commutator で埋める
    commutators: Vec<Option<CertifiedCommutator>>,
    dim: usize,
    threshold: f64,
    _g: PhantomData<G>,
}

impl<G: CommutationGrading> OperationalNet<G> {
    pub fn new(dim: usize, threshold: f64) -> Self {
        assert!(threshold > 0.0, "閾値は正であること");
        OperationalNet {
            primitives: Vec::new(),
            contexts: Vec::new(),
            commutators: Vec::new(),
            dim,
            threshold,
            _g: PhantomData,
        }
    }

    pub fn grading() -> &'static str {
        G::NAME
    }
    pub fn n_primitives(&self) -> usize {
        self.primitives.len()
    }
    pub fn threshold(&self) -> f64 {
        self.threshold
    }
    pub fn primitive(&self, id: OpId) -> &PrimitiveOperation {
        &self.primitives[id.0 as usize]
    }
    pub fn contexts(&self) -> &[BTreeSet<u32>] {
        &self.contexts
    }

    /// primitive の登録 — **Ordinary lane は odd を構成時に拒否する** (JW 誤読の遮断)
    pub fn add_primitive(&mut self, p: PrimitiveOperation) -> Result<OpId, &'static str> {
        let (_, _, d) = p.kind.matrix();
        if d != self.dim {
            return Err("primitive の次元が net と一致しない");
        }
        if p.parity == OperatorParity::Odd && !G::ACCEPTS_ODD {
            return Err(
                "OrdinaryCommutation の net は odd (fermionic) primitive を受け付けない — FermionicZ2Graded lane を使うこと (JW 弦の幾何誤読の遮断)",
            );
        }
        self.primitives.push(p);
        let k = self.primitives.len();
        self.commutators.resize(k * (k - 1) / 2, None);
        Ok(OpId((k - 1) as u32))
    }

    fn pair_index(&self, a: u32, b: u32) -> usize {
        let (i, j) = if a < b { (a, b) } else { (b, a) };
        let (i, j) = (i as usize, j as usize);
        j * (j - 1) / 2 + i
    }

    /// (graded) bracket ノルムの証明書を登録する。odd×odd 対は G が Z2 graded の
    /// とき反可換子ノルムであること (測定側の契約 — provenance に残す)。
    pub fn set_commutator(&mut self, a: OpId, b: OpId, c: CertifiedCommutator) {
        assert!(a != b, "自分自身との対に証明書は置けない");
        let idx = self.pair_index(a.0, b.0);
        self.commutators[idx] = Some(c);
    }

    pub fn commutator_verdict(&self, a: OpId, b: OpId) -> Option<CommutatorVerdict> {
        self.commutators[self.pair_index(a.0, b.0)].map(|c| c.verdict(self.threshold))
    }

    /// 操作文脈 (可換分解の一区画) の登録 — 全対に Commuting 証明書があるときのみ
    pub fn add_context(&mut self, members: &[OpId]) -> Result<usize, String> {
        for (k, &a) in members.iter().enumerate() {
            for &b in members.iter().skip(k + 1) {
                match self.commutator_verdict(a, b) {
                    Some(CommutatorVerdict::Commuting) => {}
                    Some(CommutatorVerdict::NonCommuting) => {
                        return Err(format!("文脈に非可換対 ({:?}, {:?}) が含まれる", a, b))
                    }
                    Some(CommutatorVerdict::Abstain) => {
                        return Err(format!(
                            "対 ({:?}, {:?}) の可換子区間が閾値を跨ぐ (Abstain) — 文脈を構成できない",
                            a, b
                        ))
                    }
                    None => {
                        return Err(format!("対 ({:?}, {:?}) に可換子証明書がない", a, b))
                    }
                }
            }
        }
        self.contexts.push(members.iter().map(|id| id.0).collect());
        Ok(self.contexts.len() - 1)
    }

    /// 非可換グラフの連結成分 (v32.3 factorization recovery の入力)。
    /// 証明書のない対・Abstain 対が一つでもあれば棄却 (辺の推測禁止)。
    pub fn noncommutation_components(&self) -> Result<Vec<Vec<u32>>, FactorizationAbstainReason> {
        let k = self.primitives.len();
        if k == 0 {
            return Err(FactorizationAbstainReason::InsufficientOperationalGenerators);
        }
        let mut parent: Vec<u32> = (0..k as u32).collect();
        fn find(parent: &mut Vec<u32>, x: u32) -> u32 {
            let mut r = x;
            while parent[r as usize] != r {
                r = parent[r as usize];
            }
            let mut c = x;
            while parent[c as usize] != r {
                let nx = parent[c as usize];
                parent[c as usize] = r;
                c = nx;
            }
            r
        }
        for i in 0..k as u32 {
            for j in (i + 1)..k as u32 {
                match self.commutators[self.pair_index(i, j)].map(|c| c.verdict(self.threshold)) {
                    Some(CommutatorVerdict::NonCommuting) => {
                        let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                        if ri != rj {
                            parent[ri as usize] = rj;
                        }
                    }
                    Some(CommutatorVerdict::Commuting) => {}
                    Some(CommutatorVerdict::Abstain) | None => {
                        return Err(FactorizationAbstainReason::CommutatorMarginStraddled)
                    }
                }
            }
        }
        let mut groups: std::collections::BTreeMap<u32, Vec<u32>> = Default::default();
        for i in 0..k as u32 {
            let r = find(&mut parent, i);
            groups.entry(r).or_default().push(i);
        }
        Ok(groups.into_values().collect())
    }
}

// ---------------------------------------------------------------- GlobalClosure (禁止変換 11)

/// 生成族の *-閉包の**同型類だけ**を運ぶ型 — 意図的に基底・生成元の由来を持たない。
/// **OperationalNet / 因子分解への変換は存在しない** (禁止変換 11): 閉包を取ると
/// marking は消える (v322 [N1] — site 生成族も mode 生成族も同じ B(H) に閉じる)。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct GlobalClosure {
    pub dim_hilbert: usize,
    /// *-閉包の複素次元 (B(H) なら dim_hilbert²)
    pub dim_algebra: usize,
}

impl GlobalClosure {
    pub fn is_full(&self) -> bool {
        self.dim_algebra == self.dim_hilbert * self.dim_hilbert
    }
}

// ---------------------------------------------------------------- 因子分解の読み出し裁定 (v32.3 の出力型)

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FactorizationAbstainReason {
    /// 生成族が対象 sector の B(H) を生成しない/足りない (number operator のみ等)
    InsufficientOperationalGenerators,
    /// 可換子区間が閾値を跨ぎ、非可換グラフが確定しない
    CommutatorMarginStraddled,
    /// odd primitive を ordinary lane で読もうとした等の等級不整合
    GradingMismatch,
    /// 成分閉包が full matrix factor にならない (rank 判定不能を含む)
    ComponentNotFactor,
    /// v33.1: 成分間の共同 addressability の操作的証人 (共有文脈) が無い —
    /// 代数的可換性 (CertifiedCommutator) からの暗黙昇格の禁止 (禁止変換 12)。
    /// 全対に Commuting 証明書があっても、宣言済み文脈が成分対を跨がなければ
    /// 因子分解を Exact に昇格しない
    OperationalCompatibilityUnwitnessed,
}

impl FactorizationAbstainReason {
    pub fn as_str(self) -> &'static str {
        match self {
            FactorizationAbstainReason::InsufficientOperationalGenerators => {
                "insufficient_operational_generators"
            }
            FactorizationAbstainReason::CommutatorMarginStraddled => "commutator_margin_straddled",
            FactorizationAbstainReason::GradingMismatch => "grading_mismatch",
            FactorizationAbstainReason::ComponentNotFactor => "component_not_factor",
            FactorizationAbstainReason::OperationalCompatibilityUnwitnessed => {
                "operational_compatibility_unwitnessed"
            }
        }
    }
}

/// 因子分解読み出しの裁定 — 三値 + 中心非自明の direct-sum (PROMPT/13 §2)。
/// tensor product を強制しない: 中心が非自明なら SuperselectionSectors を返す。
#[derive(Clone, Debug, PartialEq)]
pub enum FactorizationReading {
    /// 局所 unitary × ノード置換を除いて一意回復
    ExactUpToLocalUnitaryAndPermutation { local_dims: Vec<usize> },
    /// H ≅ ⊕_α (C^{m_α} ⊗ C^{n_α}) — 中心非自明の正答 (tensor 強制の禁止)
    SuperselectionSectors { sectors: Vec<(usize, usize)> },
    /// 複数の非同値 gauge orbit
    EquivalenceClassOnly { class_desc: String },
    /// 棄却 (理由つき) — 非識別セルでの正しい読み出し結果
    Abstain(FactorizationAbstainReason),
}

impl FactorizationReading {
    pub fn as_str(&self) -> &'static str {
        match self {
            FactorizationReading::ExactUpToLocalUnitaryAndPermutation { .. } => {
                "exact_up_to_local_unitary_and_permutation"
            }
            FactorizationReading::SuperselectionSectors { .. } => "superselection_sectors",
            FactorizationReading::EquivalenceClassOnly { .. } => "equivalence_class_only",
            FactorizationReading::Abstain(_) => "abstain",
        }
    }
}

// ---------------------------------------------------------------- interaction hypergraph (v32.5 の出力型)

/// H = Σ_S H_S (局所 Hilbert–Schmidt 条件期待値による直交分解) の重み w_S = ‖H_S‖_F²。
/// 局所演算子基底の選択に依存せず block-local unitary で不変 (v32.5 で機械検査)。
/// |S| = 1: on-site / 2: graph edge / 3: correlated-hopping 級 hyperedge / ≥4: 高体。
#[derive(Clone, Debug, Default)]
pub struct InteractionHypergraph {
    pub n_nodes: usize,
    pub weights: std::collections::BTreeMap<Vec<u32>, f64>,
}

impl InteractionHypergraph {
    /// |S| = k の支持 (w_S > threshold)
    pub fn support_of_order(&self, k: usize, threshold: f64) -> Vec<Vec<u32>> {
        self.weights
            .iter()
            .filter(|(s, &w)| s.len() == k && w > threshold)
            .map(|(s, _)| s.clone())
            .collect()
    }
    pub fn total_weight(&self) -> f64 {
        self.weights.values().sum()
    }
}

// ---------------------------------------------------------------- 複素行列小道具 (8×8 級の共用素子)

pub fn cmul(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let mut out = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik.re == 0.0 && aik.im == 0.0 {
                continue;
            }
            for j in 0..n {
                out[i * n + j] = out[i * n + j] + aik * b[k * n + j];
            }
        }
    }
    out
}

pub fn cdag(a: &[C64], n: usize) -> Vec<C64> {
    let mut out = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        for j in 0..n {
            let x = a[j * n + i];
            out[i * n + j] = C64::new(x.re, -x.im);
        }
    }
    out
}

/// Hilbert–Schmidt 内積 tr(a† b)
pub fn hs_inner(a: &[C64], b: &[C64]) -> C64 {
    let mut s = C64::new(0.0, 0.0);
    for (x, y) in a.iter().zip(b.iter()) {
        // conj(x) * y
        s = s + C64::new(x.re * y.re + x.im * y.im, x.re * y.im - x.im * y.re);
    }
    s
}

pub fn hs_norm(a: &[C64]) -> f64 {
    a.iter().map(|x| x.re * x.re + x.im * x.im).sum::<f64>().sqrt()
}

/// 交換子 [a, b] = ab − ba
pub fn commutator(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let ab = cmul(a, b, n);
    let ba = cmul(b, a, n);
    ab.iter().zip(ba.iter()).map(|(x, y)| *x - *y).collect()
}

/// 反交換子 {a, b} = ab + ba
pub fn anticommutator(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let ab = cmul(a, b, n);
    let ba = cmul(b, a, n);
    ab.iter().zip(ba.iter()).map(|(x, y)| *x + *y).collect()
}

/// 正規直交基底 (HS) への直交化 — 残差が rel_tol を超えたら正規化して追加。
/// 返り値 = 追加されたか。二重直交化で数値安定化。
pub fn push_ortho(basis: &mut Vec<Vec<C64>>, cand: &[C64], rel_tol: f64) -> bool {
    let scale = hs_norm(cand).max(1e-300);
    let mut v = cand.to_vec();
    for _ in 0..2 {
        for b in basis.iter() {
            let c = hs_inner(b, &v);
            for (vk, bk) in v.iter_mut().zip(b.iter()) {
                *vk = *vk - c * *bk;
            }
        }
    }
    let r = hs_norm(&v);
    if r / scale <= rel_tol {
        return false;
    }
    let inv = 1.0 / r;
    for vk in v.iter_mut() {
        *vk = vk.scale(inv);
    }
    basis.push(v);
    true
}

/// 生成族の *-代数閉包の正規直交基底 (HS Gram–Schmidt 成長 — 決定的)。
/// seed = I + 生成元 + 随伴。成長は全対の積で、追加が止まるまで。
pub fn algebra_closure(gens: &[Vec<C64>], n: usize) -> Vec<Vec<C64>> {
    let rel_tol = 1e-9;
    let mut basis: Vec<Vec<C64>> = Vec::new();
    let mut ident = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        ident[i * n + i] = C64::new(1.0, 0.0);
    }
    push_ortho(&mut basis, &ident, rel_tol);
    for g in gens {
        push_ortho(&mut basis, g, rel_tol);
        push_ortho(&mut basis, &cdag(g, n), rel_tol);
    }
    loop {
        let snapshot = basis.clone();
        let mut grew = false;
        for a in snapshot.iter() {
            for b in snapshot.iter() {
                if basis.len() >= n * n {
                    break;
                }
                let p = cmul(a, b, n);
                if push_ortho(&mut basis, &p, rel_tol) {
                    grew = true;
                }
            }
        }
        if !grew || basis.len() >= n * n {
            break;
        }
    }
    basis
}

/// 閉包の同型類 (禁止変換 11 の左辺) — これ以上の情報は運ばない
pub fn closure_of(gens: &[Vec<C64>], n: usize) -> GlobalClosure {
    GlobalClosure {
        dim_hilbert: n,
        dim_algebra: algebra_closure(gens, n).len(),
    }
}

/// commutant {M : [A_k, M] = 0 ∀k} の複素次元 — 実 2n²×2n² Gram の零空間で数える
pub fn commutant_dim(gens: &[Vec<C64>], n: usize) -> usize {
    let n2 = n * n;
    let dim_r = 2 * n2;
    // G = Σ_g L_g† L_g (実対称)。L_g は M ↦ [A_g, M] の実表現。
    let mut gram = vec![0.0; dim_r * dim_r];
    // 実基底 e_t: t < n² → E_{kl} (実), t ≥ n² → iE_{kl}
    // 各生成子 A について L e_t = [A, e_t] を陽に計算し、G += L^T L を蓄積
    for a in gens {
        // 列ベクトル像を先に全部作る (dim_r 本)
        let mut cols: Vec<Vec<f64>> = Vec::with_capacity(dim_r);
        for t in 0..dim_r {
            let (kl, imag) = (t % n2, t >= n2);
            let mut m = vec![C64::new(0.0, 0.0); n2];
            m[kl] = if imag {
                C64::new(0.0, 1.0)
            } else {
                C64::new(1.0, 0.0)
            };
            let c = commutator(a, &m, n);
            let mut col = vec![0.0; dim_r];
            for (idx, x) in c.iter().enumerate() {
                col[idx] = x.re;
                col[n2 + idx] = x.im;
            }
            cols.push(col);
        }
        for s in 0..dim_r {
            for t in s..dim_r {
                let mut acc = 0.0;
                for r in 0..dim_r {
                    acc += cols[s][r] * cols[t][r];
                }
                gram[s * dim_r + t] += acc;
                if s != t {
                    gram[t * dim_r + s] += acc;
                }
            }
        }
    }
    let (evals, _) = jacobi_eigh(&gram, dim_r);
    let emax = evals.iter().cloned().fold(0.0f64, f64::max).max(1e-300);
    let nullity = evals.iter().filter(|&&e| e <= 1e-10 * emax).count();
    // 複素次元 = 実零空間次元 / 2
    nullity / 2
}

// ---------------------------------------------------------------- 自己検査

/// operational_net の不変条件 (v322_operational_net が呼ぶ)
pub fn operational_net_self_test() -> Result<(), String> {
    // 小道具: 2×2 Pauli
    let n = 2;
    let x = vec![
        C64::new(0.0, 0.0),
        C64::new(1.0, 0.0),
        C64::new(1.0, 0.0),
        C64::new(0.0, 0.0),
    ];
    let z = vec![
        C64::new(1.0, 0.0),
        C64::new(0.0, 0.0),
        C64::new(0.0, 0.0),
        C64::new(-1.0, 0.0),
    ];
    // 1. 閉包: {X} → span{I, X} = 2 / {X, Z} → M_2 = 4
    if algebra_closure(&[x.clone()], n).len() != 2 {
        return Err("閉包 {X} が次元 2 でない".into());
    }
    let cl = closure_of(&[x.clone(), z.clone()], n);
    if !(cl.dim_algebra == 4 && cl.is_full()) {
        return Err("閉包 {X, Z} が M_2 でない".into());
    }
    // 2. commutant: {X, Z} → スカラーのみ (複素次元 1) / {Z} → 対角 (複素次元 2)
    if commutant_dim(&[x.clone(), z.clone()], n) != 1 {
        return Err("commutant {X,Z} が 1 でない".into());
    }
    if commutant_dim(&[z.clone()], n) != 2 {
        return Err("commutant {Z} が 2 でない".into());
    }
    // 3. 役割の資格審査
    let rho_ok = Preparation::certify(vec![0.5, 0.0, 0.0, 0.5], vec![0.0; 4], 2);
    if rho_ok.is_err() {
        return Err("I/2 が準備の資格を通らない".into());
    }
    if Preparation::certify(vec![1.0, 0.0, 0.0, 1.0], vec![0.0; 4], 2).is_ok() {
        return Err("tr = 2 が準備を名乗れた".into());
    }
    if Preparation::certify(vec![1.5, 0.0, 0.0, -0.5], vec![0.0; 4], 2).is_ok() {
        return Err("非 PSD が準備を名乗れた".into());
    }
    if MeasurementEffect::certify(vec![2.0, 0.0, 0.0, 0.0], vec![0.0; 4], 2).is_ok() {
        return Err("E ≤ I 違反が effect を名乗れた".into());
    }
    // 4. 証明書の 3 値裁定
    let c_edge = CertifiedCommutator::new(1.0, 1.2).unwrap();
    let c_zero = CertifiedCommutator::new(0.0, 1e-9).unwrap();
    let c_str = CertifiedCommutator::new(1e-8, 1e-2).unwrap();
    let tau = 1e-6;
    if !(c_edge.verdict(tau) == CommutatorVerdict::NonCommuting
        && c_zero.verdict(tau) == CommutatorVerdict::Commuting
        && c_str.verdict(tau) == CommutatorVerdict::Abstain)
    {
        return Err("可換子証明書の 3 値裁定が誤り".into());
    }
    // 5. Ordinary net は odd を拒否・Z2 graded は受理
    let mk_odd = || PrimitiveOperation {
        kind: OpKind::Control(
            ControlGenerator::certify(vec![0.0, 1.0, 1.0, 0.0], vec![0.0; 4], 2).unwrap(),
        ),
        parity: OperatorParity::Odd,
        provenance: "self_test",
    };
    let mut net_o: OperationalNet<OrdinaryCommutation> = OperationalNet::new(2, tau);
    if net_o.add_primitive(mk_odd()).is_ok() {
        return Err("Ordinary net が odd primitive を受理した".into());
    }
    let mut net_g: OperationalNet<FermionicZ2Graded> = OperationalNet::new(2, tau);
    if net_g.add_primitive(mk_odd()).is_err() {
        return Err("Z2 graded net が odd primitive を拒否した".into());
    }
    // 6. 文脈は可換証明書がないと構成できない・非可換対を拒否
    let mk_even = |m: &[C64]| PrimitiveOperation {
        kind: OpKind::Control(
            ControlGenerator::certify(
                m.iter().map(|c| c.re).collect(),
                m.iter().map(|c| c.im).collect(),
                2,
            )
            .unwrap(),
        ),
        parity: OperatorParity::Even,
        provenance: "self_test",
    };
    let ix = net_o.add_primitive(mk_even(&x)).unwrap();
    let iz = net_o.add_primitive(mk_even(&z)).unwrap();
    if net_o.add_context(&[ix, iz]).is_ok() {
        return Err("証明書なしの文脈が構成できた".into());
    }
    let nxz = hs_norm(&commutator(&x, &z, 2));
    net_o.set_commutator(ix, iz, CertifiedCommutator::new(nxz * 0.99, nxz * 1.01).unwrap());
    if net_o.add_context(&[ix, iz]).is_ok() {
        return Err("非可換対を含む文脈が構成できた".into());
    }
    match net_o.noncommutation_components() {
        Ok(comps) if comps.len() == 1 && comps[0].len() == 2 => {}
        _ => return Err("非可換グラフの連結成分が誤り".into()),
    }
    // 7. 裁定型の名前一意性
    let names = [
        FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims: vec![2] }.as_str(),
        FactorizationReading::SuperselectionSectors { sectors: vec![(1, 2)] }.as_str(),
        FactorizationReading::EquivalenceClassOnly { class_desc: String::new() }.as_str(),
        FactorizationReading::Abstain(FactorizationAbstainReason::InsufficientOperationalGenerators)
            .as_str(),
    ];
    for (i, a) in names.iter().enumerate() {
        for b in names.iter().skip(i + 1) {
            if a == b {
                return Err("裁定名の重複".into());
            }
        }
    }
    Ok(())
}

// ================================================================ v33.1 境界監査と型スコープ修復 (PROMPT/14)
//
// 第三十三期テーゼ「可アクセス性は作用素単体の属性ではなく、系・制御器・測定器・
// 資源・誤差証明書の関係である」の型実装の開幕。v32.3 の復元器 (v323 — 凍結原本)
// は (net, gens) を並行に受け取り、contexts を一切参照しなかった。これは定理の
// 誤りではない — 定理の仮定「各ノードの局所生成子だけが選別された primitive
// family」が型に存在せず、net と gens の整合が呼び出し側の慣行任せだった、という
// 型スコープの空隙である (v331 [B1][B2] が反例つきで機械実証)。本節はその修復:
//
//   1. **net と別渡し gens の廃止**: 復元の生成子行列は net 自身の primitive から
//      のみ取る。復元器の唯一の型付き入口は MarkedRecoveryInput — 公開コンストラクタ
//      は OperationalNet::recovery_input のみで、行列を外から注入する経路はない。
//   2. **contexts は定理入力**: 文脈 0 (NoDeclaredContexts)・被覆不完全
//      (ContextCoverageIncomplete) は構成時拒否。成分間の共同 addressability は
//      宣言済み文脈の共有 (JointContextWitness) でのみ証人され、代数的可換性
//      (CertifiedCommutator) からの暗黙変換は存在しない (**禁止変換 12**) —
//      証人不在は Abstain(OperationalCompatibilityUnwitnessed)。
//   3. **role-mixed recovery の禁止**: Control 以外の役割を含む net の復元入力は
//      構成時拒否 (RoleMixedRecovery)。測定行列を数学的に可換グラフへ混ぜる経路を
//      型で塞ぐ。測定・準備・drift の文脈意味論 (joint measurability 等) の型化は
//      v33.2 (role-typed context atlas) の主題。
//   4. **Liouvillian lane の型分離**: v32.4 の応答恒等式 R⁽¹⁾/R⁽²⁾ の資格域は
//      HamiltonianCommutatorLiouvillian (L = −i[H,·]) であり、GklsLiouvillian への
//      昇格は存在しない (**禁止変換 13**)。lane の資格審査は導分 (Leibniz) 証明書 —
//      dissipator は Leibniz を破る (v331 [B5] が γ 比例の破れと γ→0 回復を実測)。
//
// 注意: primitive 選別の循環 (independently accessible な entangler を primitive に
// 加えるだけで因子分解の読みが併合する — v331 [B2]) は本修復では解けない。どの
// 操作が accessible かの証明 (DeclaredOperation → AccessibleOperation の資格) は
// v33.2 の Certified Laboratory Interface が担う。

/// 復元入力の構成時拒否 — Abstain (読みの裁定) とは別の型エラー。
/// 「不正な入力で走らない」ことと「走った上で棄却する」ことを混ぜない。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum RecoveryInputRejection {
    /// Control 以外の役割 (準備・測定・drift) が復元入力に混入 — role-mixed
    /// recovery の禁止。測定・準備の文脈意味論は v33.2 で別型になる
    RoleMixedRecovery,
    /// 宣言済み文脈が 0 — contexts は復元定理の入力である (v33.1 の修復点)
    NoDeclaredContexts,
    /// 文脈に属さない primitive がある — 文脈 atlas が primitive family を
    /// 被覆しない限り、復元は走らない
    ContextCoverageIncomplete,
}

impl RecoveryInputRejection {
    pub fn as_str(self) -> &'static str {
        match self {
            RecoveryInputRejection::RoleMixedRecovery => "role_mixed_recovery",
            RecoveryInputRejection::NoDeclaredContexts => "no_declared_contexts",
            RecoveryInputRejection::ContextCoverageIncomplete => "context_coverage_incomplete",
        }
    }
}

/// 成分対の共同 addressability の操作的証人 — **唯一の構成は宣言済み文脈の共有**
/// (OperationalNet::joint_context_witness)。CertifiedCommutator (代数的可換性) から
/// の変換は存在しない (禁止変換 12): 一般 POVM の joint measurability は可換性で
/// 特徴づけられず、可換な作用素対が同一実験文脈で共同実行可能とは限らない。
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct JointContextWitness {
    context_index: usize,
}

impl JointContextWitness {
    pub fn context_index(&self) -> usize {
        self.context_index
    }
}

/// v33.1 復元器の唯一の型付き入口 — 生成子行列は net の primitive からのみ取る
/// (gens 別渡しの廃止)。構成は OperationalNet::recovery_input のみ。
pub struct MarkedRecoveryInput<'a, G: CommutationGrading> {
    net: &'a OperationalNet<G>,
}

/// v33.1 復元の出力 — 読み (三値裁定) と gauge orbit 照合用の成分部分代数
pub struct MarkedRecoveryDetail {
    pub reading: FactorizationReading,
    /// Exact のとき: 成分ごとの traceless 部分代数 ONB (gauge orbit 照合用)
    pub component_subalgebras: Vec<Vec<Vec<C64>>>,
}

impl<G: CommutationGrading> OperationalNet<G> {
    /// 復元入力の資格審査 (v33.1) — 役割純度 (Control のみ)・文脈の存在・被覆を
    /// 構成時に検査する。通らなければ復元は**走らない** (Abstain ではなく拒否)。
    pub fn recovery_input(&self) -> Result<MarkedRecoveryInput<'_, G>, RecoveryInputRejection> {
        if self
            .primitives
            .iter()
            .any(|p| !matches!(p.kind, OpKind::Control(_)))
        {
            return Err(RecoveryInputRejection::RoleMixedRecovery);
        }
        if self.contexts.is_empty() {
            return Err(RecoveryInputRejection::NoDeclaredContexts);
        }
        for id in 0..self.primitives.len() as u32 {
            if !self.contexts.iter().any(|c| c.contains(&id)) {
                return Err(RecoveryInputRejection::ContextCoverageIncomplete);
            }
        }
        Ok(MarkedRecoveryInput { net: self })
    }

    /// 対 (a, b) の共同 addressability の操作的証人 — 両者を含む宣言済み文脈が
    /// あるときのみ Some。可換子証明書だけからは構成できない (禁止変換 12)。
    pub fn joint_context_witness(&self, a: OpId, b: OpId) -> Option<JointContextWitness> {
        self.contexts
            .iter()
            .position(|c| c.contains(&a.0) && c.contains(&b.0))
            .map(|context_index| JointContextWitness { context_index })
    }
}

/// エルミート行列 (C64 平坦格納) の固有値 — 2n×2n 実対称埋め込み (各固有値 2 重)
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

fn ident_c64(n: usize) -> Vec<C64> {
    let mut m = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

/// span(basis) 内で全生成子と可換な部分空間 (basis が閉包のとき = 中心) の
/// エルミート正規直交基底 (v32.3 kernel の lib 移植 — dust guard 込み [v33.0-A])
pub fn closure_center_basis(basis: &[Vec<C64>], gens: &[Vec<C64>], n: usize) -> Vec<Vec<C64>> {
    let d = basis.len();
    let dim_r = 2 * d;
    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(dim_r);
    for t in 0..dim_r {
        let m: Vec<C64> = if t < d {
            basis[t].clone()
        } else {
            basis[t - d].iter().map(|c| C64::new(-c.im, c.re)).collect()
        };
        let mut col = Vec::with_capacity(gens.len() * 2 * n * n);
        for g in gens {
            let c = commutator(&m, g, n);
            for x in &c {
                col.push(x.re);
                col.push(x.im);
            }
        }
        cols.push(col);
    }
    let mut gram = vec![0.0; dim_r * dim_r];
    for s in 0..dim_r {
        for t in s..dim_r {
            let mut acc = 0.0;
            for r in 0..cols[s].len() {
                acc += cols[s][r] * cols[t][r];
            }
            gram[s * dim_r + t] = acc;
            gram[t * dim_r + s] = acc;
        }
    }
    let (evals, vecs) = jacobi_eigh(&gram, dim_r);
    let emax = evals.iter().cloned().fold(0.0f64, f64::max).max(1e-300);
    let mut out: Vec<Vec<C64>> = Vec::new();
    for (k, &e) in evals.iter().enumerate() {
        if e > 1e-10 * emax {
            continue;
        }
        let mut m = vec![C64::new(0.0, 0.0); n * n];
        for t in 0..dim_r {
            let w = vecs[t + k * dim_r];
            if w.abs() < 1e-300 {
                continue;
            }
            let coeff = if t < d {
                C64::new(w, 0.0)
            } else {
                C64::new(0.0, w)
            };
            let b = &basis[t % d];
            for (mi, bi) in m.iter_mut().zip(b.iter()) {
                *mi = *mi + coeff * *bi;
            }
        }
        let mdag = cdag(&m, n);
        let h1: Vec<C64> = m
            .iter()
            .zip(mdag.iter())
            .map(|(a, b)| (*a + *b).scale(0.5))
            .collect();
        let h2: Vec<C64> = m
            .iter()
            .zip(mdag.iter())
            .map(|(a, b)| {
                let dd = *a - *b; // 反エルミート
                C64::new(dd.im * 0.5, -dd.re * 0.5) // /(2i)
            })
            .collect();
        // dust guard (v33.0-A 設計走行の発見・統一適用): 数値塵 (‖候補‖ ≈ 0) を
        // 正規化して基底に混入させない
        if hs_norm(&h1) > 1e-9 {
            push_ortho(&mut out, &h1, 1e-8);
        }
        if hs_norm(&h2) > 1e-9 {
            push_ortho(&mut out, &h2, 1e-8);
        }
    }
    out
}

/// 中心射影の族 (Lagrange 補間, v32.3 kernel の lib 移植): T = Σ √(k+2)·H_k の
/// 固有値クラスタごとに P_α = Π_{β≠α} (T − λ_β)/(λ_α − λ_β)
pub fn closure_central_projectors(center: &[Vec<C64>], n: usize) -> Option<Vec<Vec<C64>>> {
    let mut t = vec![C64::new(0.0, 0.0); n * n];
    for (k, h) in center.iter().enumerate() {
        let w = ((k + 2) as f64).sqrt();
        for (ti, hi) in t.iter_mut().zip(h.iter()) {
            *ti = *ti + hi.scale(w);
        }
    }
    let evals = herm_evals_c64(&t, n); // 各固有値 2 重 (実埋め込み)
    let scale = evals.iter().fold(0.0f64, |a, &b| a.max(b.abs())).max(1e-300);
    let mut distinct: Vec<f64> = Vec::new();
    for &e in &evals {
        if !distinct.iter().any(|&d| (d - e).abs() <= 1e-8 * scale) {
            distinct.push(e);
        }
    }
    distinct.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let mut projs = Vec::new();
    for (a, &la) in distinct.iter().enumerate() {
        let mut p = ident_c64(n);
        for (b, &lb) in distinct.iter().enumerate() {
            if a == b {
                continue;
            }
            let mut shifted = t.clone();
            for i in 0..n {
                shifted[i * n + i] = shifted[i * n + i] - C64::new(lb, 0.0);
            }
            p = cmul(&p, &shifted, n);
            let inv = 1.0 / (la - lb);
            for x in p.iter_mut() {
                *x = x.scale(inv);
            }
        }
        let p2 = cmul(&p, &p, n);
        let idem: f64 = p2
            .iter()
            .zip(p.iter())
            .map(|(a, b)| (*a - *b).norm2())
            .sum::<f64>()
            .sqrt();
        if idem > 1e-7 {
            return None;
        }
        projs.push(p);
    }
    let mut s = vec![C64::new(0.0, 0.0); n * n];
    for p in &projs {
        for (si, pi) in s.iter_mut().zip(p.iter()) {
            *si = *si + *pi;
        }
    }
    let idm = ident_c64(n);
    let dev: f64 = s
        .iter()
        .zip(idm.iter())
        .map(|(a, b)| (*a - *b).norm2())
        .sum::<f64>()
        .sqrt();
    if dev > 1e-7 {
        return None;
    }
    Some(projs)
}

impl<'a, G: CommutationGrading> MarkedRecoveryInput<'a, G> {
    /// 復元に使う生成子行列 — **net の primitive が唯一の出所** (外部注入なし)
    pub fn generator_matrices(&self) -> Vec<Vec<C64>> {
        self.net
            .primitives
            .iter()
            .map(|p| {
                let (re, im, d) = p.kind.matrix();
                (0..d * d).map(|k| C64::new(re[k], im[k])).collect()
            })
            .collect()
    }

    pub fn contexts(&self) -> &[BTreeSet<u32>] {
        &self.net.contexts
    }

    /// v33.1 復元器 — v32.3 の凍結決定手順 (成分 → joint 閉包 → 中心 → 三値裁定,
    /// v323 と同一・dust guard 込み) の前段に、成分間の共同 addressability の
    /// 操作的証人ゲート (宣言済み文脈の共有) を積む。証人不在は
    /// Abstain(OperationalCompatibilityUnwitnessed) — 代数的可換性からの暗黙昇格を
    /// しない (禁止変換 12 の運用形)。
    pub fn recover(&self) -> MarkedRecoveryDetail {
        let n = self.net.dim;
        let gens = self.generator_matrices();
        let abstain = |r: FactorizationAbstainReason| MarkedRecoveryDetail {
            reading: FactorizationReading::Abstain(r),
            component_subalgebras: Vec::new(),
        };
        // 1. 非可換グラフの連結成分 (Abstain 対があれば棄却 — v32.3 と同一)
        let comps = match self.net.noncommutation_components() {
            Ok(c) => c,
            Err(r) => return abstain(r),
        };
        // 2. v33.1 証人ゲート: 全成分対に共有文脈 (共同 addressability の操作的証人)
        for i in 0..comps.len() {
            for j in (i + 1)..comps.len() {
                let witnessed = self.net.contexts.iter().any(|ctx| {
                    comps[i].iter().any(|a| ctx.contains(a))
                        && comps[j].iter().any(|b| ctx.contains(b))
                });
                if !witnessed {
                    return abstain(FactorizationAbstainReason::OperationalCompatibilityUnwitnessed);
                }
            }
        }
        // 3. joint 閉包 (以降は v32.3 の凍結決定手順と同一)
        let joint = algebra_closure(&gens, n);
        // 4. joint が可換なら操作資源不足 (number operator のみ等)
        let mut commutative = true;
        'outer: for a in &gens {
            for b in &gens {
                if hs_norm(&commutator(a, b, n)) > 1e-9 {
                    commutative = false;
                    break 'outer;
                }
            }
        }
        if commutative {
            return abstain(FactorizationAbstainReason::InsufficientOperationalGenerators);
        }
        // 5. 中心
        let center = closure_center_basis(&joint, &gens, n);
        if center.is_empty() {
            return abstain(FactorizationAbstainReason::ComponentNotFactor);
        }
        if center.len() == 1 {
            // 中心自明: full ∧ 各成分 factor ∧ Π d_i = n → Exact
            if joint.len() != n * n {
                return abstain(FactorizationAbstainReason::InsufficientOperationalGenerators);
            }
            let mut dims = Vec::new();
            let mut subalgebras = Vec::new();
            for comp in &comps {
                let sub: Vec<Vec<C64>> = comp.iter().map(|&i| gens[i as usize].clone()).collect();
                let cl = algebra_closure(&sub, n);
                let d2 = cl.len();
                let d = (d2 as f64).sqrt().round() as usize;
                if d * d != d2 || d < 2 {
                    return abstain(FactorizationAbstainReason::ComponentNotFactor);
                }
                let comp_center = closure_center_basis(&cl, &sub, n);
                if comp_center.len() != 1 {
                    return abstain(FactorizationAbstainReason::ComponentNotFactor);
                }
                dims.push(d);
                let idn = ident_c64(n);
                let inorm = 1.0 / (n as f64).sqrt();
                let ihat: Vec<C64> = idn.iter().map(|c| c.scale(inorm)).collect();
                let mut traceless = Vec::new();
                for b in &cl {
                    let c = hs_inner(&ihat, b);
                    let t: Vec<C64> = b
                        .iter()
                        .zip(ihat.iter())
                        .map(|(bi, ii)| *bi - c * *ii)
                        .collect();
                    push_ortho(&mut traceless, &t, 1e-9);
                }
                subalgebras.push(traceless);
            }
            let prod: usize = dims.iter().product();
            if prod != n {
                return abstain(FactorizationAbstainReason::ComponentNotFactor);
            }
            let mut sorted_dims = dims.clone();
            sorted_dims.sort_unstable();
            return MarkedRecoveryDetail {
                reading: FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
                    local_dims: sorted_dims,
                },
                component_subalgebras: subalgebras,
            };
        }
        // 6. 中心非自明: sector 分割 (tensor 強制の禁止)
        let projs = match closure_central_projectors(&center, n) {
            Some(p) => p,
            None => return abstain(FactorizationAbstainReason::ComponentNotFactor),
        };
        let mut sectors = Vec::new();
        for p in &projs {
            let tr: f64 = (0..n).map(|i| p[i * n + i].re).sum();
            let b_dim = tr.round() as usize;
            if b_dim == 0 || (tr - b_dim as f64).abs() > 1e-7 {
                return abstain(FactorizationAbstainReason::ComponentNotFactor);
            }
            let mut restricted: Vec<Vec<C64>> = Vec::new();
            for b in &joint {
                let pbp = cmul(p, &cmul(b, p, n), n);
                // dust guard: 他 sector にしか台を持たない b の像 (≈ 0) を除外
                if hs_norm(&pbp) < 1e-9 {
                    continue;
                }
                push_ortho(&mut restricted, &pbp, 1e-8);
            }
            let m2 = restricted.len();
            let m = (m2 as f64).sqrt().round() as usize;
            if m * m != m2 || b_dim % m != 0 {
                return abstain(FactorizationAbstainReason::ComponentNotFactor);
            }
            sectors.push((m, b_dim / m));
        }
        sectors.sort_unstable();
        MarkedRecoveryDetail {
            reading: FactorizationReading::SuperselectionSectors { sectors },
            component_subalgebras: Vec::new(),
        }
    }
}

/// gauge orbit 照合 (v32.3 [F3] の判定器の lib 移植 — v33.3 の profile 安定性が使う):
/// 成分 traceless 部分代数の集合が置換で min-overlap ≈ 1 に matching できるか。
/// 返り値 = (同一 orbit か, 最良 min-overlap)。
pub fn same_gauge_orbit(a: &[Vec<Vec<C64>>], b: &[Vec<Vec<C64>>]) -> (bool, f64) {
    if a.len() != b.len() {
        return (false, 0.0);
    }
    let k = a.len();
    if k == 0 {
        return (true, 1.0);
    }
    let overlap = |u: &Vec<Vec<C64>>, w: &Vec<Vec<C64>>| -> f64 {
        if u.len() != w.len() {
            return 0.0;
        }
        let mut acc = 0.0;
        for x in w {
            for y in u {
                acc += hs_inner(y, x).norm2();
            }
        }
        acc / (u.len() as f64)
    };
    let mut perm: Vec<usize> = (0..k).collect();
    let mut best = 0.0f64;
    let mut found = false;
    loop {
        let mut minov = f64::INFINITY;
        for i in 0..k {
            minov = minov.min(overlap(&a[i], &b[perm[i]]));
        }
        best = best.max(minov);
        if minov >= 1.0 - 1e-9 {
            found = true;
            break;
        }
        let mut i = k as isize - 2;
        while i >= 0 && perm[i as usize] >= perm[(i + 1) as usize] {
            i -= 1;
        }
        if i < 0 {
            break;
        }
        let mut j = k - 1;
        while perm[j] <= perm[i as usize] {
            j -= 1;
        }
        perm.swap(i as usize, j);
        perm[(i as usize + 1)..].reverse();
    }
    (found, best)
}

// ---------------------------------------------------------------- Liouvillian lane の型分離 (v33.1)

fn herm_defect_c64(m: &[C64], n: usize) -> f64 {
    let mut d = 0.0f64;
    for i in 0..n {
        for j in 0..n {
            let a = m[i * n + j];
            let b = m[j * n + i];
            d = d.max((a.re - b.re).hypot(a.im + b.im));
        }
    }
    d
}

/// v32.4 応答恒等式 R⁽¹⁾ = −i Tr(B[H,A])・R⁽²⁾ = Tr([H,B][H,A]) の資格域
/// (Schrödinger 規約 ρ(t) = e^{−iHt} ρ e^{iHt} — v31.2/v32.4 と同一)。
/// **GklsLiouvillian への変換は存在しない (禁止変換 13)** — 可換子 lane の証明書を
/// 一般 GKLS 生成子へ昇格しない。lane の資格審査は classify_generator (導分証明書)。
pub struct HamiltonianCommutatorLiouvillian {
    h: Vec<C64>,
    dim: usize,
}

impl HamiltonianCommutatorLiouvillian {
    pub fn certify(h: Vec<C64>, dim: usize) -> Result<Self, &'static str> {
        if h.len() != dim * dim {
            return Err("行列の形が n×n でない");
        }
        if herm_defect_c64(&h, dim) > 1e-12 {
            return Err("H がエルミートでない");
        }
        Ok(HamiltonianCommutatorLiouvillian { h, dim })
    }
    pub fn hamiltonian(&self) -> &[C64] {
        &self.h
    }
    pub fn dim(&self) -> usize {
        self.dim
    }
    /// L(M) = −i[H, M] (Schrödinger 規約の生成子)
    pub fn apply(&self, m: &[C64]) -> Vec<C64> {
        let c = commutator(&self.h, m, self.dim);
        c.iter().map(|x| C64::new(x.im, -x.re)).collect()
    }
    /// R⁽¹⁾ = −i Tr(B[H,A]) — 可換子 lane 専用の型付き入口 (v32.4 [L0] の恒等式)
    pub fn r1(&self, b: &[C64], a: &[C64]) -> f64 {
        let c = commutator(&self.h, a, self.dim);
        let n = self.dim;
        let mut s = C64::new(0.0, 0.0);
        for i in 0..n {
            for k in 0..n {
                s = s + b[i * n + k] * c[k * n + i];
            }
        }
        s.im
    }
    /// R⁽²⁾ = Tr([H,B][H,A]) — 同上
    pub fn r2(&self, b: &[C64], a: &[C64]) -> f64 {
        let n = self.dim;
        let hb = commutator(&self.h, b, n);
        let ha = commutator(&self.h, a, n);
        let mut s = C64::new(0.0, 0.0);
        for i in 0..n {
            for k in 0..n {
                s = s + hb[i * n + k] * ha[k * n + i];
            }
        }
        s.re
    }
}

/// GKLS 生成子 L(ρ) = −i[H,ρ] + Σ_μ γ_μ (L_μ ρ L_μ† − ½{L_μ†L_μ, ρ}) の型。
/// v32.4 の応答恒等式の資格域**外** — HamiltonianCommutatorLiouvillian との相互
/// From は存在しない (禁止変換 13)。一般 GKLS の応答理論 (jump 表現の gauge を
/// 含む) は未構成で、次々期の独立テーマ (PROMPT/14)。
pub struct GklsLiouvillian {
    h: Vec<C64>,
    jumps: Vec<Vec<C64>>,
    rates: Vec<f64>,
    dim: usize,
}

impl GklsLiouvillian {
    pub fn certify(
        h: Vec<C64>,
        jumps: Vec<Vec<C64>>,
        rates: Vec<f64>,
        dim: usize,
    ) -> Result<Self, &'static str> {
        if h.len() != dim * dim || jumps.iter().any(|l| l.len() != dim * dim) {
            return Err("行列の形が n×n でない");
        }
        if herm_defect_c64(&h, dim) > 1e-12 {
            return Err("H がエルミートでない");
        }
        if rates.len() != jumps.len() {
            return Err("γ と jump の本数が一致しない");
        }
        if rates.iter().any(|&g| !g.is_finite() || g < 0.0) {
            return Err("GKLS 資格には γ_μ ≥ 0 が要る (Kossakowski 行列の PSD)");
        }
        Ok(GklsLiouvillian {
            h,
            jumps,
            rates,
            dim,
        })
    }
    pub fn dim(&self) -> usize {
        self.dim
    }
    /// Σ_μ γ_μ (dissipator の総強度 — 0 なら可換子 lane と一致する縮退点)
    pub fn dissipator_strength(&self) -> f64 {
        self.rates.iter().sum()
    }
    /// L(M) = −i[H,M] + Σ_μ γ_μ (L_μ M L_μ† − ½{L_μ†L_μ, M})
    pub fn apply(&self, m: &[C64]) -> Vec<C64> {
        let n = self.dim;
        let c = commutator(&self.h, m, n);
        let mut out: Vec<C64> = c.iter().map(|x| C64::new(x.im, -x.re)).collect();
        for (l, &g) in self.jumps.iter().zip(self.rates.iter()) {
            if g == 0.0 {
                continue;
            }
            let ldag = cdag(l, n);
            let lml = cmul(l, &cmul(m, &ldag, n), n);
            let ll = cmul(&ldag, l, n);
            let anti = anticommutator(&ll, m, n);
            for k in 0..n * n {
                out[k] = out[k] + (lml[k] - anti[k].scale(0.5)).scale(g);
            }
        }
        out
    }
}

/// 生成子の分類 (導分証明書) — 可換子 lane への資格審査の結果。
/// 裁定は相対欠陥 (÷ max‖L(E_ab)‖) ≤ 1e-9 で行い、報告値 leibniz_defect は
/// **生の (正規化しない) 最大欠陥** — γ 比例則の照合に使えるように (v331 [B5])。
pub enum GeneratorClassification {
    /// Leibniz・†-共変・unital が全て成立 — L = −i[Ĥ,·] を復元
    HamiltonianCommutator {
        h_hat: Vec<C64>,
        leibniz_defect: f64,
        /// max_ab ‖L(E_ab) − (−i)[Ĥ, E_ab]‖ / max‖L(E_ab)‖ (相対)
        reconstruction_residual: f64,
    },
    /// Leibniz が破れる — 可換子 lane に入る資格がない (GKLS dissipator 等)
    NonDerivation { leibniz_defect: f64 },
}

/// 線形写像 L: M_n → M_n の導分 (Leibniz) 証明書による分類。
/// L(AB) = L(A)B + A·L(B)・L(A†) = L(A)†・L(I) = 0 が行列単位の全対で相対
/// ≤ 1e-9 なら L = −i[Ĥ,·] (Ĥ はエルミート・中心を除いて一意 — traceless gauge)
/// を復元して返す。破れがあれば NonDerivation (可換子 lane への入口なし)。
pub fn classify_generator(l: &dyn Fn(&[C64]) -> Vec<C64>, n: usize) -> GeneratorClassification {
    let rel_tol = 1e-9;
    let unit = |a: usize, b: usize| -> Vec<C64> {
        let mut m = vec![C64::new(0.0, 0.0); n * n];
        m[a * n + b] = C64::new(1.0, 0.0);
        m
    };
    // images[a][b] = L(E_ab)
    let images: Vec<Vec<Vec<C64>>> = (0..n)
        .map(|a| (0..n).map(|b| l(&unit(a, b))).collect())
        .collect();
    let scale = images
        .iter()
        .flatten()
        .map(|m| hs_norm(m))
        .fold(0.0f64, f64::max)
        .max(1e-300);
    // Leibniz: E_ab·E_cd = δ_bc E_ad → L(E_ab E_cd) は precompute 済み
    let mut leib = 0.0f64;
    for a in 0..n {
        for b in 0..n {
            for c in 0..n {
                for d in 0..n {
                    // lhs = δ_bc L(E_ad)
                    // rhs = L(E_ab)·E_cd + E_ab·L(E_cd)
                    //   M·E_cd: 列 d に M の列 c / E_ab·M: 行 a に M の行 b
                    let mut defect = vec![C64::new(0.0, 0.0); n * n];
                    if b == c {
                        defect.copy_from_slice(&images[a][d]);
                    }
                    let m1 = &images[a][b];
                    for x in 0..n {
                        defect[x * n + d] = defect[x * n + d] - m1[x * n + c];
                    }
                    let m2 = &images[c][d];
                    for y in 0..n {
                        defect[a * n + y] = defect[a * n + y] - m2[b * n + y];
                    }
                    leib = leib.max(hs_norm(&defect));
                }
            }
        }
    }
    // †-共変: L(E_ba) = L(E_ab)† / unital: Σ_a L(E_aa) = 0
    let mut dag_defect = 0.0f64;
    for a in 0..n {
        for b in 0..n {
            let lhs = &images[b][a];
            let rhs = cdag(&images[a][b], n);
            let d: f64 = lhs
                .iter()
                .zip(rhs.iter())
                .map(|(x, y)| (*x - *y).norm2())
                .sum::<f64>()
                .sqrt();
            dag_defect = dag_defect.max(d);
        }
    }
    let mut uni = vec![C64::new(0.0, 0.0); n * n];
    for a in 0..n {
        for (u, x) in uni.iter_mut().zip(images[a][a].iter()) {
            *u = *u + *x;
        }
    }
    let uni_defect = hs_norm(&uni);
    let total = (leib.max(dag_defect).max(uni_defect)) / scale;
    if total > rel_tol {
        return GeneratorClassification::NonDerivation { leibniz_defect: leib };
    }
    // 復元: L(E_a0) = −i(H E_a0 − E_a0 H) → H_ca − δ_ca H_00 = i·(L(E_a0))_{c0}
    let mut h_hat = vec![C64::new(0.0, 0.0); n * n];
    for a in 0..n {
        for c in 0..n {
            let v = images[a][0][c * n]; // (L(E_a0))_{c,0}
            h_hat[c * n + a] = C64::new(-v.im, v.re); // i·v
        }
    }
    // traceless gauge (中心 H + cI は原理的に不可視 — v32.4 [L1]) + エルミート化
    let tr = (0..n).fold(C64::new(0.0, 0.0), |s, i| s + h_hat[i * n + i]);
    let sh = tr.scale(1.0 / n as f64);
    for i in 0..n {
        h_hat[i * n + i] = h_hat[i * n + i] - sh;
    }
    let hd = cdag(&h_hat, n);
    for (x, y) in h_hat.iter_mut().zip(hd.iter()) {
        *x = (*x + *y).scale(0.5);
    }
    // 残差: max_ab ‖L(E_ab) − (−i)[Ĥ, E_ab]‖ / scale
    let mut resid = 0.0f64;
    for a in 0..n {
        for b in 0..n {
            let c = commutator(&h_hat, &unit(a, b), n);
            let want: Vec<C64> = c.iter().map(|x| C64::new(x.im, -x.re)).collect();
            let d: f64 = images[a][b]
                .iter()
                .zip(want.iter())
                .map(|(x, y)| (*x - *y).norm2())
                .sum::<f64>()
                .sqrt();
            resid = resid.max(d);
        }
    }
    GeneratorClassification::HamiltonianCommutator {
        h_hat,
        leibniz_defect: leib,
        reconstruction_residual: resid / scale,
    }
}

// ---------------------------------------------------------------- v33.1 自己検査

/// v33.1 型スコープ修復の不変条件 (v331_scope_repair が呼ぶ)。既存の
/// operational_net_self_test には触れない — v32.2 契約の凍結出力を変えないため。
pub fn scope_repair_self_test() -> Result<(), String> {
    let n = 4usize;
    let tau = 1e-3;
    // 2 qubit の site 生成族 X₁, Z₁, X₂, Z₂ (dim 4)
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
    let gens = [
        kron2(&px, &id2),
        kron2(&pz, &id2),
        kron2(&id2, &px),
        kron2(&id2, &pz),
    ];
    let mk_ctrl = |g: &[C64]| PrimitiveOperation {
        kind: OpKind::Control(
            ControlGenerator::certify(
                g.iter().map(|c| c.re).collect(),
                g.iter().map(|c| c.im).collect(),
                n,
            )
            .unwrap(),
        ),
        parity: OperatorParity::Even,
        provenance: "scope_repair_self_test",
    };
    let build = |with_certs: bool| -> (OperationalNet<OrdinaryCommutation>, Vec<OpId>) {
        let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
        let ids: Vec<OpId> = gens.iter().map(|g| net.add_primitive(mk_ctrl(g)).unwrap()).collect();
        if with_certs {
            for a in 0..gens.len() {
                for b in (a + 1)..gens.len() {
                    let nu = hs_norm(&commutator(&gens[a], &gens[b], n));
                    net.set_commutator(
                        ids[a],
                        ids[b],
                        CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
                    );
                }
            }
        }
        (net, ids)
    };
    // 1. 文脈 0 → 構成時拒否 (contexts は定理入力)
    let (net0, _) = build(true);
    match net0.recovery_input() {
        Err(RecoveryInputRejection::NoDeclaredContexts) => {}
        _ => return Err("文脈 0 の net が復元入力を構成できた".into()),
    }
    // 2. 被覆不完全 → 構成時拒否
    let (mut net1, ids1) = build(true);
    net1.add_context(&[ids1[0], ids1[2]]).map_err(|e| e)?; // {X₁, X₂}
    match net1.recovery_input() {
        Err(RecoveryInputRejection::ContextCoverageIncomplete) => {}
        _ => return Err("被覆不完全の net が復元入力を構成できた".into()),
    }
    // 3. full atlas → Exact [2,2] (v32.3 手順の再現)
    let (mut net2, ids2) = build(true);
    net2.add_context(&[ids2[0], ids2[2]]).map_err(|e| e)?; // {X₁, X₂}
    net2.add_context(&[ids2[1], ids2[3]]).map_err(|e| e)?; // {Z₁, Z₂}
    let det = net2.recovery_input().map_err(|e| e.as_str())?.recover();
    if det.reading
        != (FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims: vec![2, 2] })
    {
        return Err(format!("full atlas の復元が {} (期待 Exact [2,2])", det.reading.as_str()));
    }
    // 4. singleton 文脈のみ → 証人不在で Abstain (代数的可換 ↛ 操作的両立)
    let (mut net3, ids3) = build(true);
    for id in &ids3 {
        net3.add_context(&[*id]).map_err(|e| e)?;
    }
    let det3 = net3.recovery_input().map_err(|e| e.as_str())?.recover();
    if det3.reading
        != FactorizationReading::Abstain(
            FactorizationAbstainReason::OperationalCompatibilityUnwitnessed,
        )
    {
        return Err(format!(
            "singleton 文脈の復元が {} (期待 Abstain(unwitnessed))",
            det3.reading.as_str()
        ));
    }
    // 5. JointContextWitness は共有文脈からのみ
    if net3.joint_context_witness(ids3[0], ids3[2]).is_some() {
        return Err("共有文脈なしに JointContextWitness が構成できた".into());
    }
    if net2.joint_context_witness(ids2[0], ids2[2]).is_none() {
        return Err("共有文脈があるのに JointContextWitness が構成できない".into());
    }
    // 6. role-mixed → 構成時拒否
    let (mut net4, _) = build(true);
    let n1: Vec<C64> = ident_c64(n)
        .iter()
        .zip(kron2(&pz, &id2).iter())
        .map(|(i, z)| (*i - *z).scale(0.5))
        .collect();
    net4.add_primitive(PrimitiveOperation {
        kind: OpKind::Measure(
            MeasurementEffect::certify(
                n1.iter().map(|c| c.re).collect(),
                n1.iter().map(|c| c.im).collect(),
                n,
            )
            .unwrap(),
        ),
        parity: OperatorParity::Even,
        provenance: "scope_repair_self_test",
    })
    .map_err(|e| e)?;
    match net4.recovery_input() {
        Err(RecoveryInputRejection::RoleMixedRecovery) => {}
        _ => return Err("role-mixed net が復元入力を構成できた".into()),
    }
    // 7. Liouvillian lane: 可換子は導分・GKLS は Leibniz 破れ
    let hx: Vec<C64> = kron2(&px, &id2);
    let lane = HamiltonianCommutatorLiouvillian::certify(hx.clone(), n)
        .map_err(|e| e.to_string())?;
    match classify_generator(&|m: &[C64]| lane.apply(m), n) {
        GeneratorClassification::HamiltonianCommutator {
            reconstruction_residual,
            ..
        } => {
            if reconstruction_residual > 1e-10 {
                return Err(format!("可換子 lane の復元残差 {:.1e}", reconstruction_residual));
            }
        }
        GeneratorClassification::NonDerivation { leibniz_defect } => {
            return Err(format!("可換子 lane が NonDerivation ({:.1e})", leibniz_defect));
        }
    }
    let sminus = vec![
        C64::new(0.0, 0.0),
        C64::new(1.0, 0.0),
        C64::new(0.0, 0.0),
        C64::new(0.0, 0.0),
    ];
    let gkls = GklsLiouvillian::certify(hx, vec![kron2(&sminus, &id2)], vec![0.5], n)
        .map_err(|e| e.to_string())?;
    match classify_generator(&|m: &[C64]| gkls.apply(m), n) {
        GeneratorClassification::NonDerivation { leibniz_defect } => {
            if leibniz_defect < 1e-3 {
                return Err(format!("GKLS の Leibniz 破れが小さすぎる ({:.1e})", leibniz_defect));
            }
        }
        GeneratorClassification::HamiltonianCommutator { .. } => {
            return Err("GKLS (γ > 0) が可換子 lane の資格を通った".into());
        }
    }
    if GklsLiouvillian::certify(ident_c64(n), vec![ident_c64(n)], vec![-1.0], n).is_ok() {
        return Err("負の γ が GKLS 資格を通った".into());
    }
    // 8. 拒否理由名の一意性
    let names = [
        RecoveryInputRejection::RoleMixedRecovery.as_str(),
        RecoveryInputRejection::NoDeclaredContexts.as_str(),
        RecoveryInputRejection::ContextCoverageIncomplete.as_str(),
        FactorizationAbstainReason::OperationalCompatibilityUnwitnessed.as_str(),
    ];
    for (i, a) in names.iter().enumerate() {
        for b in names.iter().skip(i + 1) {
            if a == b {
                return Err("拒否理由名の重複".into());
            }
        }
    }
    Ok(())
}
