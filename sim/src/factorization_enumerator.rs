// factorization_enumerator — sector-aware complete finite factorization
// candidate enumerator (v34.4, PROMPT/15 §5)
//
// FollowUp の FAC-001 判定「代数不変量の診断は再現・任意の有限次元代数に対する
// complete tensor-factor candidate search は未提供」への応答。与えられた
// *-部分代数の族 (marked 成分の閉包など) から:
//   1. 中心 Z(A) の最小中心射影を列挙 (既存 closure_central_projectors)
//   2. 各 sector の Wedderburn 構造を**証明書つき**で分離
//      (dim A_α = n²・dim A′_α = m²・n·m = d_α・積 span = d_α²・二重可換子 A″= A)
//   3. multiplicity 空間と simple factor の分離 (n_α, m_α)
//   4. 族の可換部分集合から factorization candidate を列挙
//      (candidate = 互いに可換な simple 成分 + 残差可換子が simple のときの補因子)
//   5. 局所 unitary × 成分置換の witness (traceless 部分代数の部分空間 overlap
//      matching — OCS-1.0 §F3 と同じ意味論・バー 0.9)
//   6. 一意でなければ候補集合を返す (tie-break 禁止)
//
// 出力 6 型 (凍結):
//   UniqueFactorization        — 候補が gauge orbit を除いて 1 つ
//   FactorizationCandidateSet  — 非同値な候補が複数 (無制約選択の禁止)
//   SectorwiseFactorization    — 中心非自明: sector ごとの (n_α, m_α) が正答
//   IncompletePrimitiveSet     — 族に simple 因子候補が無い (abelian のみ等)
//   NontrivialCenterObstruction — 中心非自明なのに大域 tensor を要求した
//   ScopeExceeded              — dense バー超過 (試行しない — 正答)
//
// Wedderburn の同型 (A_α ≅ M_{n_α} ⊗ I_{m_α}) 自体は有限次元 *-代数の標準理論
// (C0) — 本モジュールが機械保証するのは上記証明書 (次元・積 span・二重可換子) と
// 候補列挙・orbit witness であり、同型の存在は証明書 + 標準理論から従う。

use crate::operational_net::{algebra_closure, closure_center_basis, closure_central_projectors};
use crate::C64;

// ---------------------------------------------------------------- 小道具 (HS 幾何)

fn hs_inner(a: &[C64], b: &[C64]) -> C64 {
    let mut s = C64::new(0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        s = s + x.conj() * *y;
    }
    s
}

fn hs_norm(a: &[C64]) -> f64 {
    a.iter().map(|x| x.norm2()).sum::<f64>().sqrt()
}

/// HS Gram–Schmidt: v を基底に直交化して加える (相対閾値) — 加われば true
fn push_ortho(basis: &mut Vec<Vec<C64>>, v: &[C64], rel_tol: f64) -> bool {
    let mut w: Vec<C64> = v.to_vec();
    let n0 = hs_norm(&w);
    if n0 < 1e-300 {
        return false;
    }
    for _ in 0..2 {
        for b in basis.iter() {
            let c = hs_inner(b, &w);
            for (wi, bi) in w.iter_mut().zip(b) {
                *wi = *wi - c * *bi;
            }
        }
    }
    let nrm = hs_norm(&w);
    if nrm > rel_tol * n0.max(1.0) {
        let inv = 1.0 / nrm;
        for wi in w.iter_mut() {
            *wi = wi.scale(inv);
        }
        basis.push(w);
        true
    } else {
        false
    }
}

fn matmul_c(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let mut c = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        for k in 0..n {
            let aik = a[i * n + k];
            if aik.norm2() == 0.0 {
                continue;
            }
            for j in 0..n {
                c[i * n + j] = c[i * n + j] + aik * b[k * n + j];
            }
        }
    }
    c
}

fn commutator_c(a: &[C64], b: &[C64], n: usize) -> Vec<C64> {
    let ab = matmul_c(a, b, n);
    let ba = matmul_c(b, a, n);
    ab.iter().zip(&ba).map(|(x, y)| *x - *y).collect()
}

/// 行列単位 E_ij の HS 正規直交族 (M_d の完全基底)
fn matrix_units(n: usize) -> Vec<Vec<C64>> {
    let mut out = Vec::with_capacity(n * n);
    for i in 0..n {
        for j in 0..n {
            let mut m = vec![C64::new(0.0, 0.0); n * n];
            m[i * n + j] = C64::new(1.0, 0.0);
            out.push(m);
        }
    }
    out
}

/// 全行列空間内の可換子零空間 = commutant のエルミート正規直交基底。
/// closure_center_basis の零空間しきい値は稠密 (無理数成分) の共役族で数値塵の
/// 方向を拾い得る — 返す基底の可換性を**機械再検証**して塵を落とす (証明書は
/// 再実行で検証する、の規律。v35.0-A の設計走行が発見した故障モード)。
pub fn commutant_basis(gens: &[Vec<C64>], n: usize) -> Vec<Vec<C64>> {
    let raw = closure_center_basis(&matrix_units(n), gens, n);
    let mut clean: Vec<Vec<C64>> = Vec::new();
    for m in raw {
        let mn = hs_norm(&m).max(1e-300);
        let mut worst = 0.0f64;
        for g in gens {
            let gn = hs_norm(g).max(1e-300);
            let c = commutator_c(&m, g, n);
            worst = worst.max(hs_norm(&c) / (mn * gn));
        }
        if worst < 1e-11 {
            push_ortho(&mut clean, &m, 1e-9);
        }
    }
    clean
}

/// (診断用) 基底元ごとの最悪可換子残差
pub fn commutant_residuals(gens: &[Vec<C64>], n: usize) -> Vec<f64> {
    let raw = closure_center_basis(&matrix_units(n), gens, n);
    raw.iter()
        .map(|m| {
            let mn = hs_norm(m).max(1e-300);
            gens.iter()
                .map(|g| hs_norm(&commutator_c(m, g, n)) / (mn * hs_norm(g).max(1e-300)))
                .fold(0.0, f64::max)
        })
        .collect()
}

/// 射影 p の像の正規直交列 → 等長 V (n×r, 列優先で r 本の n ベクトル)
fn range_isometry(p: &[C64], n: usize) -> Vec<Vec<C64>> {
    let mut cols: Vec<Vec<C64>> = Vec::new();
    for j in 0..n {
        let col: Vec<C64> = (0..n).map(|i| p[i * n + j]).collect();
        let nrm0 = col.iter().map(|x| x.norm2()).sum::<f64>().sqrt();
        if nrm0 < 1e-9 {
            continue;
        }
        let mut w = col;
        for _ in 0..2 {
            for b in cols.iter() {
                let mut c = C64::new(0.0, 0.0);
                for (x, y) in b.iter().zip(&w) {
                    c = c + x.conj() * *y;
                }
                for (wi, bi) in w.iter_mut().zip(b) {
                    *wi = *wi - c * *bi;
                }
            }
        }
        let nrm = w.iter().map(|x| x.norm2()).sum::<f64>().sqrt();
        if nrm > 1e-9 {
            let inv = 1.0 / nrm;
            for wi in w.iter_mut() {
                *wi = wi.scale(inv);
            }
            cols.push(w);
        }
    }
    cols
}

/// 圧縮 V† X V (V = r 本の列, 出力 r×r)
fn compress(x: &[C64], v: &[Vec<C64>], n: usize) -> Vec<C64> {
    let r = v.len();
    let mut out = vec![C64::new(0.0, 0.0); r * r];
    for a in 0..r {
        for b in 0..r {
            let mut s = C64::new(0.0, 0.0);
            for i in 0..n {
                for j in 0..n {
                    s = s + v[a][i].conj() * x[i * n + j] * v[b][j];
                }
            }
            out[a * r + b] = s;
        }
    }
    out
}

// ---------------------------------------------------------------- 出力型 (凍結 6 値)

/// sector ごとの Wedderburn 証明書
#[derive(Clone, Debug)]
pub struct SectorCertificate {
    pub sector_dim: usize,       // d_α
    pub simple_dim: usize,       // n_α (A_α ≅ M_{n_α} ⊗ I_{m_α})
    pub multiplicity: usize,     // m_α
    pub dims_ok: bool,           // dim A_α = n²・dim A′_α = m²・n·m = d_α
    pub product_span_full: bool, // span(A_α·A′_α) = d_α² (A ∨ A′ = B(H_α))
    pub double_commutant_ok: bool, // A″_α = A_α (次元一致)
}

impl SectorCertificate {
    pub fn certified(&self) -> bool {
        self.dims_ok && self.product_span_full && self.double_commutant_ok
    }
}

/// 列挙器の裁定 (凍結 6 値) — 単一分解を強制しない
#[derive(Debug)]
pub enum EnumeratorReading {
    UniqueFactorization {
        local_dims: Vec<usize>,
        /// 候補の成分 (traceless 部分の正規直交基底) — orbit 照合用
        components: Vec<Vec<Vec<C64>>>,
    },
    FactorizationCandidateSet {
        candidate_dims: Vec<Vec<usize>>,
        candidates: Vec<Vec<Vec<Vec<C64>>>>,
    },
    SectorwiseFactorization {
        sectors: Vec<SectorCertificate>,
    },
    IncompletePrimitiveSet,
    NontrivialCenterObstruction {
        n_sectors: usize,
    },
    ScopeExceeded,
}

impl EnumeratorReading {
    pub fn as_str(&self) -> &'static str {
        match self {
            EnumeratorReading::UniqueFactorization { .. } => "unique_factorization",
            EnumeratorReading::FactorizationCandidateSet { .. } => "factorization_candidate_set",
            EnumeratorReading::SectorwiseFactorization { .. } => "sectorwise_factorization",
            EnumeratorReading::IncompletePrimitiveSet => "incomplete_primitive_set",
            EnumeratorReading::NontrivialCenterObstruction { .. } => {
                "nontrivial_center_obstruction"
            }
            EnumeratorReading::ScopeExceeded => "scope_exceeded",
        }
    }
}

/// dense 列挙のスコープバー (これを超える dim は試行しない — ScopeExceeded が正答)
pub const ENUMERATOR_DENSE_BAR: usize = 64;

// ---------------------------------------------------------------- Wedderburn 解析

/// 閉包基底の traceless 部分の正規直交基底 (orbit 照合用の成分表現)
fn traceless_part(basis: &[Vec<C64>], n: usize) -> Vec<Vec<C64>> {
    let mut out: Vec<Vec<C64>> = Vec::new();
    let inv_n = 1.0 / n as f64;
    for b in basis {
        let mut tr = C64::new(0.0, 0.0);
        for i in 0..n {
            tr = tr + b[i * n + i];
        }
        let mut t = b.clone();
        for i in 0..n {
            t[i * n + i] = t[i * n + i] - tr.scale(inv_n);
        }
        push_ortho(&mut out, &t, 1e-9);
    }
    out
}

/// 単一 sector (中心自明が前提) の Wedderburn 証明書を組む
pub fn certify_sector(gens: &[Vec<C64>], d: usize) -> SectorCertificate {
    let a_basis = algebra_closure(gens, d);
    let ap_basis = commutant_basis(gens, d);
    let dim_a = a_basis.len();
    let dim_ap = ap_basis.len();
    let n_f = (dim_a as f64).sqrt().round() as usize;
    let m_f = (dim_ap as f64).sqrt().round() as usize;
    let dims_ok = n_f * n_f == dim_a && m_f * m_f == dim_ap && n_f * m_f == d;
    // 積 span: span{a·b} = d² (A ∨ A′ = B(H))
    let mut prod: Vec<Vec<C64>> = Vec::new();
    'outer: for a in &a_basis {
        for b in &ap_basis {
            push_ortho(&mut prod, &matmul_c(a, b, d), 1e-9);
            if prod.len() == d * d {
                break 'outer;
            }
        }
    }
    let product_span_full = prod.len() == d * d;
    // 二重可換子: dim A″ = dim A
    let app = commutant_basis(&ap_basis, d);
    let double_commutant_ok = app.len() == dim_a;
    SectorCertificate {
        sector_dim: d,
        simple_dim: n_f,
        multiplicity: m_f,
        dims_ok,
        product_span_full,
        double_commutant_ok,
    }
}

/// 族の joint 解析: 中心自明なら None、非自明なら sector 証明書列を返す
pub fn sectorwise_analysis(all_gens: &[Vec<C64>], d: usize) -> Option<Vec<SectorCertificate>> {
    let closure = algebra_closure(all_gens, d);
    let center = closure_center_basis(&closure, all_gens, d);
    if center.len() <= 1 {
        return None;
    }
    let projs = closure_central_projectors(&center, d)?;
    let mut sectors = Vec::new();
    for p in &projs {
        let v = range_isometry(p, d);
        let da = v.len();
        let gens_a: Vec<Vec<C64>> = all_gens.iter().map(|g| compress(g, &v, d)).collect();
        sectors.push(certify_sector(&gens_a, da));
    }
    sectors.sort_by_key(|s| (s.sector_dim, s.simple_dim));
    Some(sectors)
}

// ---------------------------------------------------------------- orbit 照合

/// 成分 (traceless 正規直交基底) 対の部分空間 overlap ∈ [0,1]
fn component_overlap(a: &[Vec<C64>], b: &[Vec<C64>]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return if a.len() == b.len() { 1.0 } else { 0.0 };
    }
    let mut s = 0.0;
    for x in a {
        for y in b {
            s += hs_inner(x, y).norm2();
        }
    }
    s / a.len().max(b.len()) as f64
}

/// 候補 (成分リスト) 同士の gauge orbit 照合 — 成分置換にわたる matching
/// (全成分の overlap ≥ bar を満たす置換が存在するか)。バーは OCS-1.0 の 0.9。
pub fn same_candidate_orbit(a: &[Vec<Vec<C64>>], b: &[Vec<Vec<C64>>], bar: f64) -> (bool, f64) {
    if a.len() != b.len() {
        return (false, 0.0);
    }
    let k = a.len();
    let mut perm: Vec<usize> = (0..k).collect();
    let mut best = 0.0f64;
    // k ≤ 5 想定 — 全置換
    fn permute(
        perm: &mut Vec<usize>,
        i: usize,
        a: &[Vec<Vec<C64>>],
        b: &[Vec<Vec<C64>>],
        best: &mut f64,
    ) {
        let k = perm.len();
        if i == k {
            let mut mn = f64::INFINITY;
            for (ia, &ib) in perm.iter().enumerate() {
                mn = mn.min(component_overlap(&a[ia], &b[ib]));
            }
            if mn > *best {
                *best = mn;
            }
            return;
        }
        for j in i..k {
            perm.swap(i, j);
            permute(perm, i + 1, a, b, best);
            perm.swap(i, j);
        }
    }
    permute(&mut perm, 0, a, b, &mut best);
    (best >= bar, best)
}

// ---------------------------------------------------------------- 候補列挙

/// 族の可換部分集合から factorization candidate を列挙する (証明書つき)。
/// demand_global_tensor = true のとき、中心非自明は NontrivialCenterObstruction。
pub fn enumerate_candidates(
    family: &[Vec<Vec<C64>>],
    d: usize,
    demand_global_tensor: bool,
) -> EnumeratorReading {
    if d > ENUMERATOR_DENSE_BAR {
        return EnumeratorReading::ScopeExceeded;
    }
    let all_gens: Vec<Vec<C64>> = family.iter().flat_map(|f| f.iter().cloned()).collect();
    // joint 中心の検査。abelian 族 (中心 = 閉包) は sector の皮を被せず
    // IncompletePrimitiveSet — rank-1 sector の列挙は因子候補として空虚 (v32.3 の
    // number-op-only = Insufficient の enumerator 版)
    {
        let closure = algebra_closure(&all_gens, d);
        let center = closure_center_basis(&closure, &all_gens, d);
        if center.len() > 1 && center.len() == closure.len() {
            return EnumeratorReading::IncompletePrimitiveSet;
        }
    }
    if let Some(sectors) = sectorwise_analysis(&all_gens, d) {
        if demand_global_tensor {
            return EnumeratorReading::NontrivialCenterObstruction {
                n_sectors: sectors.len(),
            };
        }
        return EnumeratorReading::SectorwiseFactorization { sectors };
    }
    // 各メンバーの閉包・simple 資格
    struct Member {
        basis: Vec<Vec<C64>>,
        simple_dim: usize,
        simple: bool,
    }
    let members: Vec<Member> = family
        .iter()
        .map(|gens| {
            let basis = algebra_closure(gens, d);
            let center = closure_center_basis(&basis, gens, d);
            let nf = (basis.len() as f64).sqrt().round() as usize;
            Member {
                simple: center.len() == 1 && nf * nf == basis.len() && nf > 1,
                simple_dim: nf,
                basis,
            }
        })
        .collect();
    // 可換性 (メンバー基底の全対)
    let k = members.len();
    let mut commute = vec![vec![false; k]; k];
    for i in 0..k {
        for j in 0..k {
            if i == j {
                commute[i][j] = true;
                continue;
            }
            let mut mx = 0.0f64;
            for a in &members[i].basis {
                for b in &members[j].basis {
                    mx = mx.max(hs_norm(&commutator_c(a, b, d)));
                }
            }
            commute[i][j] = mx < 1e-9;
        }
    }
    // 部分集合列挙: simple メンバーの**極大**可換集合 S — 補因子 = (∨S)′ が
    // simple なら候補。極大性を要求する理由: S の外に S と可換な simple メンバー
    // j があるとき、S 単独の補因子 (合成 commutant) は j の閉包を含む冗長な表現で、
    // 稠密 (無理数成分) の共役族では commutant の数値零空間が不安定になる
    // (v35.0-A の設計走行が発見)。明示のメンバーがあるならそれを使う。
    let mut cands: Vec<(Vec<usize>, Vec<Vec<Vec<C64>>>)> = Vec::new();
    'subset: for mask in 1u32..(1 << k) {
        let idx: Vec<usize> = (0..k).filter(|&i| (mask >> i) & 1 == 1).collect();
        if idx.iter().any(|&i| !members[i].simple) {
            continue;
        }
        for a in 0..idx.len() {
            for b in (a + 1)..idx.len() {
                if !commute[idx[a]][idx[b]] {
                    continue 'subset;
                }
            }
        }
        // 極大性: S の外の simple メンバー j が S 全体と可換なら S は非極大
        for j in 0..k {
            if !idx.contains(&j)
                && members[j].simple
                && idx.iter().all(|&i| commute[j][i])
            {
                continue 'subset;
            }
        }
        let prod_n: usize = idx.iter().map(|&i| members[i].simple_dim).product();
        if prod_n > d || d % prod_n != 0 {
            continue;
        }
        // 補因子: S の joint commutant
        let s_gens: Vec<Vec<C64>> = idx
            .iter()
            .flat_map(|&i| family[i].iter().cloned())
            .collect();
        let comm = commutant_basis(&s_gens, d);
        let m_res = d / prod_n;
        let mut comps: Vec<Vec<Vec<C64>>> = idx
            .iter()
            .map(|&i| traceless_part(&members[i].basis, d))
            .collect();
        let mut dims: Vec<usize> = idx.iter().map(|&i| members[i].simple_dim).collect();
        if m_res > 1 {
            // 補因子が simple (dim = m² で中心自明) のときのみ候補として資格
            if comm.len() != m_res * m_res {
                continue;
            }
            let comm_center = closure_center_basis(&comm, &comm, d);
            if comm_center.len() != 1 {
                continue;
            }
            comps.push(traceless_part(&comm, d));
            dims.push(m_res);
        } else if comm.len() != 1 {
            // 完全被覆なら補因子はスカラーのはず
            continue;
        }
        // 証明書: 全成分の積 span = d²
        let mut prod_basis: Vec<Vec<C64>> = vec![{
            let mut ident = vec![C64::new(0.0, 0.0); d * d];
            for i in 0..d {
                ident[i * d + i] = C64::new(1.0 / (d as f64).sqrt(), 0.0);
            }
            ident
        }];
        {
            let mut pool: Vec<Vec<C64>> = prod_basis.clone();
            let flat: Vec<&Vec<C64>> = comps.iter().flat_map(|c| c.iter()).collect();
            loop {
                let mut grew = false;
                let snapshot = pool.clone();
                for p in &snapshot {
                    for f in &flat {
                        let m = matmul_c(p, f, d);
                        if push_ortho(&mut pool, &m, 1e-9) {
                            grew = true;
                        }
                        if pool.len() == d * d {
                            break;
                        }
                    }
                    if pool.len() == d * d {
                        break;
                    }
                }
                if !grew || pool.len() == d * d {
                    break;
                }
            }
            prod_basis = pool;
        }
        if prod_basis.len() != d * d {
            continue;
        }
        // 既存候補との orbit 照合 (成分数・次元が同じで matching すれば同一候補)
        let mut sorted_dims = dims.clone();
        sorted_dims.sort();
        let mut dup = false;
        for (cd, cc) in &cands {
            let mut cds = cd.clone();
            cds.sort();
            if cds == sorted_dims && same_candidate_orbit(&comps, cc, 0.9).0 {
                dup = true;
                break;
            }
        }
        if !dup {
            cands.push((dims, comps));
        }
    }
    // 極大性: 成分数が最大の候補だけを残す (粗い候補は細い候補に吸収される —
    // 例: {A} 単独の [2,3] は {A,B} の [2,3] と同一 orbit で dedup 済み。
    // 残る非同値対は真の候補集合)
    if cands.is_empty() {
        return EnumeratorReading::IncompletePrimitiveSet;
    }
    let max_len = cands.iter().map(|(d0, _)| d0.len()).max().unwrap();
    let maximal: Vec<(Vec<usize>, Vec<Vec<Vec<C64>>>)> = cands
        .into_iter()
        .filter(|(d0, _)| d0.len() == max_len)
        .collect();
    if maximal.len() == 1 {
        let (dims, comps) = maximal.into_iter().next().unwrap();
        EnumeratorReading::UniqueFactorization {
            local_dims: dims,
            components: comps,
        }
    } else {
        EnumeratorReading::FactorizationCandidateSet {
            candidate_dims: maximal.iter().map(|(d0, _)| d0.clone()).collect(),
            candidates: maximal.into_iter().map(|(_, c)| c).collect(),
        }
    }
}

// ---------------------------------------------------------------- 自己検査

pub fn factorization_enumerator_self_test() -> Result<(), String> {
    // 2 qubit site 族 → Unique [2,2]
    let d = 4usize;
    let x = [0.0, 1.0, 1.0, 0.0];
    let z = [1.0, 0.0, 0.0, -1.0];
    let mk = |m2: &[f64; 4], site: usize| -> Vec<C64> {
        let mut out = vec![C64::new(0.0, 0.0); d * d];
        for i in 0..2 {
            for j in 0..2 {
                for k in 0..2 {
                    let (r, c) = if site == 0 {
                        (i * 2 + k, j * 2 + k)
                    } else {
                        (k * 2 + i, k * 2 + j)
                    };
                    out[r * d + c] = C64::new(m2[i * 2 + j], 0.0);
                }
            }
        }
        out
    };
    let fam = vec![
        vec![mk(&x, 0), mk(&z, 0)],
        vec![mk(&x, 1), mk(&z, 1)],
    ];
    match enumerate_candidates(&fam, d, false) {
        EnumeratorReading::UniqueFactorization { local_dims, .. } => {
            let mut s = local_dims.clone();
            s.sort();
            if s != vec![2, 2] {
                return Err(format!("self test dims {:?}", s));
            }
        }
        other => return Err(format!("self test verdict {}", other.as_str())),
    }
    Ok(())
}
