// structured_backend — 構造化スケーリングの 3 lane (v33.6, PROMPT/14)
//
// 一般 dense 行列の *-閉包をそのまま高速化して「大型系へスケールした」と主張して
// はいけない — qubit 数に対して行列次元が指数増大する。lane を明示的に分ける:
//
//   GenericDenseBackend      — 小次元の完全一般系 (v33.1 の復元器)。dim バーを
//                              超える入力は **ScopeExceeded** (試行すらしない)
//   PauliSymplecticBackend   — Pauli/Clifford 型: 演算子 = GF(2)^{2n} ベクトル。
//                              可換性 = symplectic 形式・閉包次元 = 2^{dim V}・
//                              中心 = ω の radical — **行列を生成せず** 32–64 qubit
//                              の証明書を返す
//   MajoranaQuadraticBackend — quadratic fermion / Majorana 双線形: 実反対称行列
//                              A (2N×2N)。ブロック構造 = 支持グラフの成分・閉包 =
//                              Lie (so) 閉包。quadratic 閉包は full M_d を**与えない**
//                              (偶 Clifford = パリティ超選択) — dense 対応原理は
//                              「支持分割 = dense 成分・dense *-閉包 dim = 2^{2m−1}」
//
// 規律:
//   - 各 backend は自分の scope 外の入力を **ScopeExceeded** で拒否する (非 Pauli
//     和は PauliVector に構成不能・dense はサイズバー・quadratic は反対称行列のみ)。
//   - **generic dense の小系成果を大型一般系へ昇格しない (禁止変換 21)** — 大型の
//     主張は structured lane の scope 内でのみ立ち、scope 外は ScopeExceeded が正しい
//     裁定である。
//   - Pauli lane の「証明書」は宣言された Pauli 構造の厳密代数データであって測定
//     ではない (ノイズ lane は dense の区間証明書の領分 — v33.6 は交差させない)。
//
// 一次ソース: docs/uft-v33.6.md / core.schema.yml (概念 + 禁止変換 21)。
// 整合は v336_structured_backend が機械検査する。

use crate::operational_net::{FactorizationAbstainReason, FactorizationReading, RecoveryInputRejection};
use crate::C64;
use std::collections::BTreeSet;

// ---------------------------------------------------------------- scope 裁定

/// scope 外入力の正しい裁定 — 「できない」を型で言う
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StructuredScopeError {
    /// dense backend のサイズバー (dim ≤ 4096) を超えた
    DimensionTooLargeForDense,
    /// Pauli 文字列 (± 位相つき単項) でない (和・一般行列)
    NotPauliString,
    /// 反対称行列でない (quadratic lane の表現外)
    NotAntisymmetric,
}

impl StructuredScopeError {
    pub fn as_str(self) -> &'static str {
        match self {
            StructuredScopeError::DimensionTooLargeForDense => "dimension_too_large_for_dense",
            StructuredScopeError::NotPauliString => "not_pauli_string",
            StructuredScopeError::NotAntisymmetric => "not_antisymmetric",
        }
    }
}

/// dense backend の凍結サイズバー — これを超える一般行列入力は ScopeExceeded
pub const DENSE_DIM_BAR: usize = 4096;

pub fn dense_scope_guard(dim: usize) -> Result<(), StructuredScopeError> {
    if dim > DENSE_DIM_BAR {
        return Err(StructuredScopeError::DimensionTooLargeForDense);
    }
    Ok(())
}

// ---------------------------------------------------------------- GF(2) 線形代数

/// GF(2) ベクトル (ビット集合, 語長可変)
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gf2Vec {
    pub words: Vec<u64>,
    pub nbits: usize,
}

impl Gf2Vec {
    pub fn zeros(nbits: usize) -> Self {
        Gf2Vec {
            words: vec![0; (nbits + 63) / 64],
            nbits,
        }
    }
    pub fn get(&self, i: usize) -> bool {
        (self.words[i / 64] >> (i % 64)) & 1 == 1
    }
    pub fn set(&mut self, i: usize, v: bool) {
        if v {
            self.words[i / 64] |= 1u64 << (i % 64);
        } else {
            self.words[i / 64] &= !(1u64 << (i % 64));
        }
    }
    pub fn xor_assign(&mut self, other: &Gf2Vec) {
        for (a, b) in self.words.iter_mut().zip(other.words.iter()) {
            *a ^= b;
        }
    }
    pub fn dot(&self, other: &Gf2Vec) -> bool {
        let mut acc = 0u32;
        for (a, b) in self.words.iter().zip(other.words.iter()) {
            acc ^= (a & b).count_ones() & 1;
        }
        acc & 1 == 1
    }
    pub fn is_zero(&self) -> bool {
        self.words.iter().all(|&w| w == 0)
    }
}

/// GF(2) rank (ガウス消去 — 入力は複製して破壊)
pub fn gf2_rank(rows: &[Gf2Vec]) -> usize {
    let mut mat: Vec<Gf2Vec> = rows.to_vec();
    let nbits = rows.first().map(|r| r.nbits).unwrap_or(0);
    let mut rank = 0usize;
    for col in 0..nbits {
        let mut pivot = None;
        for r in rank..mat.len() {
            if mat[r].get(col) {
                pivot = Some(r);
                break;
            }
        }
        let Some(p) = pivot else { continue };
        mat.swap(rank, p);
        let pivot_row = mat[rank].clone();
        for (r, row) in mat.iter_mut().enumerate() {
            if r != rank && row.get(col) {
                row.xor_assign(&pivot_row);
            }
        }
        rank += 1;
        if rank == mat.len() {
            break;
        }
    }
    rank
}

// ---------------------------------------------------------------- Pauli symplectic backend

/// Pauli 演算子の GF(2) 表現 (x|z) — 位相は可換性・閉包次元に影響しないため持たない
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PauliVector {
    pub x: Gf2Vec,
    pub z: Gf2Vec,
    pub n_qubits: usize,
}

impl PauliVector {
    pub fn from_str(s: &str) -> Self {
        let n = s.len();
        let mut x = Gf2Vec::zeros(n);
        let mut z = Gf2Vec::zeros(n);
        for (i, c) in s.chars().enumerate() {
            match c {
                'I' => {}
                'X' => x.set(i, true),
                'Z' => z.set(i, true),
                'Y' => {
                    x.set(i, true);
                    z.set(i, true);
                }
                _ => panic!("未知の Pauli 文字"),
            }
        }
        PauliVector { x, z, n_qubits: n }
    }
    /// symplectic 形式 ω(P, Q) = x_P·z_Q + z_P·x_Q (true = 反可換)
    pub fn anticommutes(&self, other: &PauliVector) -> bool {
        self.x.dot(&other.z) ^ self.z.dot(&other.x)
    }
    /// (x|z) の結合ベクトル (rank 計算用)
    pub fn combined(&self) -> Gf2Vec {
        let n = self.n_qubits;
        let mut v = Gf2Vec::zeros(2 * n);
        for i in 0..n {
            if self.x.get(i) {
                v.set(i, true);
            }
            if self.z.get(i) {
                v.set(n + i, true);
            }
        }
        v
    }
    /// dense 行列 (2^n) が ± 単一 Pauli 文字列かの資格審査つき構成 —
    /// 和 (X₁+X₂ 等) は NotPauliString。小 n の照合用 (n ≤ 12)。
    pub fn from_dense(m: &[C64], dim: usize) -> Result<Self, StructuredScopeError> {
        let n = (dim as f64).log2().round() as usize;
        if 1usize << n != dim {
            return Err(StructuredScopeError::NotPauliString);
        }
        // 各 qubit の 1 文字を先頭行の非零パターンから同定するのは煩雑なため、
        // 全 4^n を舐めずに: 非零成分がちょうど 2^n 個・絶対値が全て等しいことを
        // 確認し、列位置のビットパターンから x を、符号/位相から z を読む。
        let mut nz = Vec::new();
        let mut amp = 0.0f64;
        for (k, c) in m.iter().enumerate() {
            let a = (c.re * c.re + c.im * c.im).sqrt();
            if a > 1e-12 {
                nz.push(k);
                if amp == 0.0 {
                    amp = a;
                } else if (a - amp).abs() > 1e-9 {
                    return Err(StructuredScopeError::NotPauliString);
                }
            }
        }
        if nz.len() != dim || (amp - 1.0).abs() > 1e-9 {
            return Err(StructuredScopeError::NotPauliString);
        }
        // 行 r の非零列は r ^ xmask で一定
        let r0 = nz[0] / dim;
        let c0 = nz[0] % dim;
        let xmask = r0 ^ c0;
        for &k in &nz {
            let (r, c) = (k / dim, k % dim);
            if r ^ c != xmask {
                return Err(StructuredScopeError::NotPauliString);
            }
        }
        // z は行ごとの符号パターン (±/±i) から: entry(r, r^xmask) = phase·(−1)^{z·r}
        let ent = |r: usize| m[r * dim + (r ^ xmask)];
        let base = ent(0);
        let mut zmask = 0usize;
        for b in 0..n {
            let e = ent(1 << b);
            // e/base = ±1 (z ビット)
            let ratio_re = (e.re * base.re + e.im * base.im) / (base.re * base.re + base.im * base.im);
            if (ratio_re - 1.0).abs() < 1e-9 {
                // z_b = 0
            } else if (ratio_re + 1.0).abs() < 1e-9 {
                zmask |= 1 << b;
            } else {
                return Err(StructuredScopeError::NotPauliString);
            }
        }
        // 全行の整合検査
        for r in 0..dim {
            let want_sign = if ((zmask & r).count_ones() & 1) == 1 {
                -1.0
            } else {
                1.0
            };
            let e = ent(r);
            let dev = ((e.re - want_sign * base.re).powi(2) + (e.im - want_sign * base.im).powi(2))
                .sqrt();
            if dev > 1e-9 {
                return Err(StructuredScopeError::NotPauliString);
            }
        }
        // ビット b ↔ qubit (n−1−b) (kron の並び — op3 と同一規約)
        let mut x = Gf2Vec::zeros(n);
        let mut z = Gf2Vec::zeros(n);
        for b in 0..n {
            let q = n - 1 - b;
            if (xmask >> b) & 1 == 1 {
                x.set(q, true);
            }
            if (zmask >> b) & 1 == 1 {
                z.set(q, true);
            }
        }
        Ok(PauliVector { x, z, n_qubits: n })
    }
}

/// Pauli net の宣言 (演算子 + 文脈) — 全て厳密代数データ
pub struct PauliNetSpec {
    pub n_qubits: usize,
    pub ops: Vec<PauliVector>,
    pub contexts: Vec<BTreeSet<usize>>,
}

/// PauliSymplecticBackend の復元 (v33.1 の凍結手順の GF(2) 実装):
///   資格 (文脈の存在・被覆) → 成分 (ω) → 証人ゲート → 可換なら Insufficient →
///   V = span・radical r → r = 0: full (dim V = 2n) ∧ 成分 factor (radical 0・偶) ∧
///   Π d = 2^n → Exact / r > 0: SuperselectionSectors [(m, n_α); 2^r]。
///   **2^n 次元の行列はどこにも現れない。**
pub fn recover_pauli_net(
    spec: &PauliNetSpec,
) -> Result<FactorizationReading, RecoveryInputRejection> {
    let k = spec.ops.len();
    let n = spec.n_qubits;
    if spec.contexts.is_empty() {
        return Err(RecoveryInputRejection::NoDeclaredContexts);
    }
    for i in 0..k {
        if !spec.contexts.iter().any(|c| c.contains(&i)) {
            return Err(RecoveryInputRejection::ContextCoverageIncomplete);
        }
    }
    // 成分 (ω による反可換グラフ)
    let mut parent: Vec<usize> = (0..k).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = x;
        while parent[c] != r {
            let nx = parent[c];
            parent[c] = r;
            c = nx;
        }
        r
    }
    let mut any_anticommute = false;
    for a in 0..k {
        for b in (a + 1)..k {
            if spec.ops[a].anticommutes(&spec.ops[b]) {
                any_anticommute = true;
                let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..k {
        groups.entry(find(&mut parent, i)).or_default().push(i);
    }
    let comps: Vec<Vec<usize>> = groups.into_values().collect();
    // 証人ゲート
    for i in 0..comps.len() {
        for j in (i + 1)..comps.len() {
            let witnessed = spec.contexts.iter().any(|ctx| {
                comps[i].iter().any(|a| ctx.contains(a)) && comps[j].iter().any(|b| ctx.contains(b))
            });
            if !witnessed {
                return Ok(FactorizationReading::Abstain(
                    FactorizationAbstainReason::OperationalCompatibilityUnwitnessed,
                ));
            }
        }
    }
    if !any_anticommute {
        return Ok(FactorizationReading::Abstain(
            FactorizationAbstainReason::InsufficientOperationalGenerators,
        ));
    }
    // V = span (基底に約す)・radical rank
    let vecs: Vec<Gf2Vec> = spec.ops.iter().map(|p| p.combined()).collect();
    let dim_v = gf2_rank(&vecs);
    // Gram_ω over ops (rank は基底の取り方に依らない — 生成集合の Gram の rank =
    // span 上の ω の rank)
    let mut gram_rows: Vec<Gf2Vec> = Vec::with_capacity(k);
    for a in 0..k {
        let mut row = Gf2Vec::zeros(k);
        for b in 0..k {
            if spec.ops[a].anticommutes(&spec.ops[b]) {
                row.set(b, true);
            }
        }
        gram_rows.push(row);
    }
    let rank_omega = gf2_rank(&gram_rows);
    let radical = dim_v - rank_omega;
    if radical == 0 {
        if dim_v != 2 * n {
            return Ok(FactorizationReading::Abstain(
                FactorizationAbstainReason::InsufficientOperationalGenerators,
            ));
        }
        // 各成分: dim V_c 偶・radical 0 → d_c = 2^{dim V_c / 2}
        let mut dims = Vec::new();
        let mut sum_dim = 0usize;
        for comp in &comps {
            let cv: Vec<Gf2Vec> = comp.iter().map(|&i| vecs[i].clone()).collect();
            let dv = gf2_rank(&cv);
            let mut grows: Vec<Gf2Vec> = Vec::with_capacity(comp.len());
            for (ia, &a) in comp.iter().enumerate() {
                let mut row = Gf2Vec::zeros(comp.len());
                for (ib, &b) in comp.iter().enumerate() {
                    let _ = ia;
                    if spec.ops[a].anticommutes(&spec.ops[b]) {
                        row.set(ib, true);
                    }
                }
                grows.push(row);
            }
            let r_omega = gf2_rank(&grows);
            if dv % 2 != 0 || dv != r_omega {
                return Ok(FactorizationReading::Abstain(
                    FactorizationAbstainReason::ComponentNotFactor,
                ));
            }
            dims.push(1usize << (dv / 2));
            sum_dim += dv;
        }
        if sum_dim != 2 * n {
            return Ok(FactorizationReading::Abstain(
                FactorizationAbstainReason::ComponentNotFactor,
            ));
        }
        dims.sort_unstable();
        return Ok(FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims: dims });
    }
    // 超選択: sectors = 2^r・m = 2^{(dim V − r)/2}・n_α = 2^n / (2^r · m)
    let m_pow2 = dim_v - radical;
    if m_pow2 % 2 != 0 {
        return Ok(FactorizationReading::Abstain(
            FactorizationAbstainReason::ComponentNotFactor,
        ));
    }
    let m = 1usize << (m_pow2 / 2);
    let n_sectors = 1usize << radical;
    let total = 1usize << n;
    if total % (n_sectors * m) != 0 {
        return Ok(FactorizationReading::Abstain(
            FactorizationAbstainReason::ComponentNotFactor,
        ));
    }
    let mult = total / (n_sectors * m);
    let sectors = vec![(m, mult); n_sectors];
    Ok(FactorizationReading::SuperselectionSectors { sectors })
}

// ---------------------------------------------------------------- Majorana quadratic backend

/// quadratic 生成子 — 実反対称 2N×2N (H = (i/4) Σ A_ab γ_a γ_b の係数)
#[derive(Clone, Debug)]
pub struct QuadraticGenerator {
    pub a: Vec<f64>,
    pub m: usize,
}

impl QuadraticGenerator {
    pub fn certify(a: Vec<f64>, m: usize) -> Result<Self, StructuredScopeError> {
        if a.len() != m * m {
            return Err(StructuredScopeError::NotAntisymmetric);
        }
        for i in 0..m {
            for j in 0..m {
                if (a[i * m + j] + a[j * m + i]).abs() > 1e-12 {
                    return Err(StructuredScopeError::NotAntisymmetric);
                }
            }
        }
        Ok(QuadraticGenerator { a, m })
    }
}

fn mat_comm(a: &[f64], b: &[f64], m: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; m * m];
    for i in 0..m {
        for k in 0..m {
            let aik = a[i * m + k];
            let bik = b[i * m + k];
            if aik == 0.0 && bik == 0.0 {
                continue;
            }
            for j in 0..m {
                out[i * m + j] += aik * b[k * m + j] - bik * a[k * m + j];
            }
        }
    }
    out
}

fn frob(a: &[f64]) -> f64 {
    a.iter().map(|x| x * x).sum::<f64>().sqrt()
}

/// 実行列族の線形 span への直交化 push (Frobenius) — 追加されたかを返す
fn push_ortho_real(basis: &mut Vec<Vec<f64>>, cand: &[f64], rel_tol: f64) -> bool {
    let scale = frob(cand).max(1e-300);
    let mut v = cand.to_vec();
    for _ in 0..2 {
        for b in basis.iter() {
            let c: f64 = b.iter().zip(v.iter()).map(|(x, y)| x * y).sum::<f64>()
                / b.iter().map(|x| x * x).sum::<f64>();
            for (vk, bk) in v.iter_mut().zip(b.iter()) {
                *vk -= c * bk;
            }
        }
    }
    let r = frob(&v);
    if r / scale <= rel_tol {
        return false;
    }
    basis.push(v);
    true
}

/// quadratic lane の読み — Majorana 支持ブロックと so 閉包
pub enum QuadraticBlockReading {
    Blocks {
        /// 各ブロックの Majorana 本数 (昇順)
        block_majoranas: Vec<usize>,
        /// 各ブロックの Lie 閉包が so(block) に到達したか
        lie_full: Vec<bool>,
    },
}

impl QuadraticBlockReading {
    /// dense 対応原理: ブロック (2m 本) の quadratic *-閉包の複素次元は 2^{2m−1}
    /// (偶 Clifford — full M_d ではない: パリティ超選択)
    pub fn predicted_dense_closure_dim(block_majoranas: usize) -> usize {
        1usize << (block_majoranas - 1)
    }
}

/// 支持グラフの成分 + ブロックごとの Lie (so) 閉包次元。**2^N は現れない。**
pub fn recover_quadratic_blocks(gens: &[QuadraticGenerator]) -> QuadraticBlockReading {
    assert!(!gens.is_empty(), "生成子が空");
    let m = gens[0].m;
    // 支持グラフ: a, b が結合 ⟺ ある生成子で A_ab ≠ 0
    let mut parent: Vec<usize> = (0..m).collect();
    fn find(parent: &mut Vec<usize>, x: usize) -> usize {
        let mut r = x;
        while parent[r] != r {
            r = parent[r];
        }
        let mut c = x;
        while parent[c] != r {
            let nx = parent[c];
            parent[c] = r;
            c = nx;
        }
        r
    }
    for g in gens {
        for i in 0..m {
            for j in (i + 1)..m {
                if g.a[i * m + j].abs() > 1e-12 {
                    let (ri, rj) = (find(&mut parent, i), find(&mut parent, j));
                    if ri != rj {
                        parent[ri] = rj;
                    }
                }
            }
        }
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..m {
        groups.entry(find(&mut parent, i)).or_default().push(i);
    }
    // 孤立 Majorana (どの生成子にも触られない) はブロックに数えない
    let touched: Vec<bool> = (0..m)
        .map(|i| {
            gens.iter()
                .any(|g| (0..m).any(|j| g.a[i * m + j].abs() > 1e-12))
        })
        .collect();
    let mut blocks: Vec<Vec<usize>> = groups
        .into_values()
        .filter(|blk| blk.iter().any(|&i| touched[i]))
        .collect();
    blocks.sort_by_key(|b| b.len());
    let mut block_majoranas = Vec::new();
    let mut lie_full = Vec::new();
    for blk in &blocks {
        let bm = blk.len();
        // ブロックに制限した生成子で Lie 閉包を成長
        let restrict = |g: &QuadraticGenerator| -> Vec<f64> {
            let mut out = vec![0.0f64; bm * bm];
            for (ii, &gi) in blk.iter().enumerate() {
                for (jj, &gj) in blk.iter().enumerate() {
                    out[ii * bm + jj] = g.a[gi * m + gj];
                }
            }
            out
        };
        let mut basis: Vec<Vec<f64>> = Vec::new();
        for g in gens {
            let r = restrict(g);
            if frob(&r) > 1e-12 {
                push_ortho_real(&mut basis, &r, 1e-9);
            }
        }
        loop {
            let snapshot = basis.clone();
            let mut grew = false;
            let target = bm * (bm - 1) / 2;
            for a in snapshot.iter() {
                for b in snapshot.iter() {
                    if basis.len() >= target {
                        break;
                    }
                    let c = mat_comm(a, b, bm);
                    if frob(&c) > 1e-9 && push_ortho_real(&mut basis, &c, 1e-9) {
                        grew = true;
                    }
                }
            }
            if !grew || basis.len() >= target {
                break;
            }
        }
        block_majoranas.push(bm);
        lie_full.push(basis.len() == bm * (bm - 1) / 2);
    }
    QuadraticBlockReading::Blocks {
        block_majoranas,
        lie_full,
    }
}

// ---------------------------------------------------------------- 自己検査

/// structured_backend の不変条件 (v336_structured_backend が呼ぶ)
pub fn structured_backend_self_test() -> Result<(), String> {
    // GF(2): rank と symplectic
    let x1 = PauliVector::from_str("XI");
    let z1 = PauliVector::from_str("ZI");
    let x2 = PauliVector::from_str("IX");
    if !x1.anticommutes(&z1) || x1.anticommutes(&x2) {
        return Err("symplectic 形式の値が誤り".into());
    }
    let r = gf2_rank(&[x1.combined(), z1.combined(), x2.combined()]);
    if r != 3 {
        return Err(format!("GF(2) rank = {} (期待 3)", r));
    }
    // Pauli 復元: 2 qubit site → Exact [2,2]
    let spec = PauliNetSpec {
        n_qubits: 2,
        ops: vec![
            PauliVector::from_str("XI"),
            PauliVector::from_str("ZI"),
            PauliVector::from_str("IX"),
            PauliVector::from_str("IZ"),
        ],
        contexts: vec![
            [0usize, 2].into_iter().collect(),
            [1usize, 3].into_iter().collect(),
        ],
    };
    match recover_pauli_net(&spec) {
        Ok(FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims })
            if local_dims == vec![2, 2] => {}
        r => {
            return Err(format!(
                "Pauli 復元が {:?}",
                r.map(|x| x.as_str().to_string()).map_err(|e| e.as_str())
            ))
        }
    }
    // scope: dense guard・非 Pauli 和
    if dense_scope_guard(1 << 20).is_ok() {
        return Err("dense サイズバーが機能しない".into());
    }
    let sum_x1x2: Vec<C64> = {
        let mk = |s: &str| -> Vec<C64> {
            // 2 qubit の文字列 → dense (検査用の最小実装)
            let p = |c: char| -> Vec<C64> {
                let (o, l) = (C64::new(0.0, 0.0), C64::new(1.0, 0.0));
                match c {
                    'I' => vec![l, o, o, l],
                    'X' => vec![o, l, l, o],
                    'Z' => vec![l, o, o, C64::new(-1.0, 0.0)],
                    _ => panic!(),
                }
            };
            let cs: Vec<char> = s.chars().collect();
            let (a, b) = (p(cs[0]), p(cs[1]));
            let mut out = vec![C64::new(0.0, 0.0); 16];
            for i1 in 0..2 {
                for j1 in 0..2 {
                    for i2 in 0..2 {
                        for j2 in 0..2 {
                            out[(i1 * 2 + i2) * 4 + (j1 * 2 + j2)] =
                                a[i1 * 2 + j1] * b[i2 * 2 + j2];
                        }
                    }
                }
            }
            out
        };
        let (a, b) = (mk("XI"), mk("IX"));
        a.iter().zip(b.iter()).map(|(x, y)| *x + *y).collect()
    };
    match PauliVector::from_dense(&sum_x1x2, 4) {
        Err(StructuredScopeError::NotPauliString) => {}
        _ => return Err("非 Pauli 和が PauliVector を名乗れた".into()),
    }
    // quadratic: 2 ブロック (4+2 Majorana)・so 閉包
    let m = 6usize;
    let mk_hop = |i: usize, j: usize| -> QuadraticGenerator {
        let mut a = vec![0.0f64; m * m];
        a[i * m + j] = 1.0;
        a[j * m + i] = -1.0;
        QuadraticGenerator::certify(a, m).unwrap()
    };
    let gens = vec![mk_hop(0, 1), mk_hop(1, 2), mk_hop(2, 3), mk_hop(4, 5)];
    let QuadraticBlockReading::Blocks {
        block_majoranas,
        lie_full,
    } = recover_quadratic_blocks(&gens);
    if block_majoranas != vec![2, 4] || lie_full != vec![true, true] {
        return Err(format!(
            "quadratic ブロックが {:?} / {:?}",
            block_majoranas, lie_full
        ));
    }
    if QuadraticBlockReading::predicted_dense_closure_dim(4) != 8 {
        return Err("dense 対応原理の予言が誤り".into());
    }
    // 非反対称の拒否
    if QuadraticGenerator::certify(vec![1.0; 4], 2).is_ok() {
        return Err("非反対称が quadratic 資格を通った".into());
    }
    Ok(())
}
