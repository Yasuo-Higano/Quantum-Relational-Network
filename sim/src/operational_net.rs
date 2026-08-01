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
