// contextual_factorization — 文脈 atlas 上の局所復元と overlap glue (v33.4, PROMPT/14)
//
// v33.1–v33.3 は「単一の大域 net」から読む器械だった。実際の実験室は**持ち場
// (chart)** の集まりである — 研究室 A は qubit 1,2 を・B は qubit 2,3 を制御し、
// 全系を一望する制御器はどこにもない。v33.4 はこの atlas 構造を一級市民にする:
//
//   1. **chart ごとの局所復元**: chart (実験者が宣言する primitive の部分集合) の
//      部分 net で、成分 → 閉包 → factor 資格の局所読みを行う。chart は自分の
//      持ち場しか語らない — 大域 fullness は chart の資格要件ではない。
//   2. **overlap 上の algebra matching**: 共有 primitive または交差非可換
//      (certified NonCommuting) で連結した成分対は、**同一の因子部分代数**
//      (traceless ONB の min-overlap ≈ 1) を指していなければならない。
//   3. **裁定 (凍結)**: 全 matching 整合 + 被覆完全 + 大域証人 → GluedExact
//      (v334 [A1] が直接大域復元との一致 = glue 定理を機械照合)。matching が
//      破れても「全被覆の整合 chart 群」が 2 つ以上あれば EquivalenceClassOnly
//      (site atlas vs DFT atlas — 無制約 tie-break の禁止)。それ以外は
//      Abstain(GlueInconsistent) — 宣言された chart を黙って捨てない。
//   4. **chart の局所 Exact は大域主張ではない** (**禁止変換 19**):
//      ChartLocalFactorization → GlobalFactorization の変換は存在しない。大域は
//      glue (matching + 被覆 + 証人) だけが与える。
//
// 一次ソース: docs/uft-v33.4.md / core.schema.yml (概念 + 禁止変換 19)。
// 整合は v334_contextual_factorization が機械検査する。

use crate::operational_net::{
    algebra_closure, closure_center_basis, hs_inner, hs_norm, push_ortho, CommutationGrading,
    CommutatorVerdict, OpId, OperationalNet, PrimitiveOperation,
};
use crate::C64;
use std::collections::BTreeSet;

// ---------------------------------------------------------------- chart

/// 実験者が宣言する持ち場 — net の primitive の部分集合
#[derive(Clone, Debug)]
pub struct ChartSpec {
    pub primitive_ids: Vec<OpId>,
}

/// chart 局所復元の失敗理由
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ChartFailure {
    /// 証明書欠落・跨ぎで成分が確定しない
    ComponentsUndetermined,
    /// chart 内の成分対に証人文脈がない (v33.1 の witness 規律の局所版)
    LocalWitnessMissing,
    /// ある成分の閉包が full matrix factor にならない
    ComponentNotFactor,
    /// chart 内に文脈が 1 つもない / 被覆しない
    ContextIncomplete,
}

impl ChartFailure {
    pub fn as_str(self) -> &'static str {
        match self {
            ChartFailure::ComponentsUndetermined => "components_undetermined",
            ChartFailure::LocalWitnessMissing => "local_witness_missing",
            ChartFailure::ComponentNotFactor => "component_not_factor",
            ChartFailure::ContextIncomplete => "context_incomplete",
        }
    }
}

/// chart 局所復元の結果 — **自分の持ち場の因子しか語らない** (大域主張ではない:
/// GlobalFactorization への変換は存在しない — 禁止変換 19)
pub struct ChartLocalFactorization {
    /// chart 内の因子: (次元 d, traceless 部分代数 ONB, 構成 primitive の大域 id)
    pub factors: Vec<(usize, Vec<Vec<C64>>, Vec<OpId>)>,
}

/// chart の局所復元 (凍結手順): chart 部分 net (証明書は大域 net から**継承** —
/// 再計算しない) の成分 → chart 内証人ゲート → 成分閉包の factor 資格。
pub fn recover_chart<G: CommutationGrading>(
    net: &OperationalNet<G>,
    chart: &ChartSpec,
    n: usize,
) -> Result<ChartLocalFactorization, ChartFailure> {
    let ids = &chart.primitive_ids;
    let k = ids.len();
    if k == 0 {
        return Err(ChartFailure::ContextIncomplete);
    }
    // chart 内文脈 = 大域文脈のうち chart に完全に含まれるもの
    let idset: BTreeSet<u32> = ids.iter().map(|i| i.0).collect();
    let charted_contexts: Vec<&BTreeSet<u32>> = net
        .contexts()
        .iter()
        .filter(|c| c.is_subset(&idset))
        .collect();
    if charted_contexts.is_empty()
        || ids
            .iter()
            .any(|i| !charted_contexts.iter().any(|c| c.contains(&i.0)))
    {
        return Err(ChartFailure::ContextIncomplete);
    }
    // 成分 (証明書は大域 net の記録 — Abstain/欠落は棄却)
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
    for a in 0..k {
        for b in (a + 1)..k {
            match net.commutator_verdict(ids[a], ids[b]) {
                Some(CommutatorVerdict::NonCommuting) => {
                    let (ra, rb) = (find(&mut parent, a), find(&mut parent, b));
                    if ra != rb {
                        parent[ra] = rb;
                    }
                }
                Some(CommutatorVerdict::Commuting) => {}
                Some(CommutatorVerdict::Abstain) | None => {
                    return Err(ChartFailure::ComponentsUndetermined)
                }
            }
        }
    }
    let mut groups: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
    for i in 0..k {
        let r = find(&mut parent, i);
        groups.entry(r).or_default().push(i);
    }
    let comps: Vec<Vec<usize>> = groups.into_values().collect();
    // chart 内証人ゲート (v33.1 の局所版): 成分対ごとに chart 内共有文脈
    for i in 0..comps.len() {
        for j in (i + 1)..comps.len() {
            let witnessed = charted_contexts.iter().any(|ctx| {
                comps[i].iter().any(|&a| ctx.contains(&ids[a].0))
                    && comps[j].iter().any(|&b| ctx.contains(&ids[b].0))
            });
            if !witnessed {
                return Err(ChartFailure::LocalWitnessMissing);
            }
        }
    }
    // 成分閉包 = full matrix factor M_d (d ≥ 2)・成分中心自明
    let matrix_of = |id: OpId| -> Vec<C64> {
        let (re, im, d) = net.primitive(id).kind.matrix();
        (0..d * d).map(|t| C64::new(re[t], im[t])).collect()
    };
    let mut factors = Vec::new();
    for comp in &comps {
        let gens: Vec<Vec<C64>> = comp.iter().map(|&i| matrix_of(ids[i])).collect();
        let cl = algebra_closure(&gens, n);
        let d2 = cl.len();
        let d = (d2 as f64).sqrt().round() as usize;
        if d * d != d2 || d < 2 {
            return Err(ChartFailure::ComponentNotFactor);
        }
        if closure_center_basis(&cl, &gens, n).len() != 1 {
            return Err(ChartFailure::ComponentNotFactor);
        }
        // traceless 部分代数 ONB (matching 用 — v32.3/v33.1 と同一の構成)
        let mut idn = vec![C64::new(0.0, 0.0); n * n];
        for i in 0..n {
            idn[i * n + i] = C64::new(1.0, 0.0);
        }
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
        factors.push((d, traceless, comp.iter().map(|&i| ids[i]).collect()));
    }
    Ok(ChartLocalFactorization { factors })
}

// ---------------------------------------------------------------- atlas の glue

/// atlas 裁定の棄却理由
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum AtlasAbstainReason {
    /// ある chart の局所復元が失敗
    ChartFailed(ChartFailure),
    /// matching が破れ、全被覆の整合 chart 群も存在しない
    GlueInconsistent,
    /// 整合だが因子次元の積が n に届かない (未 address 領域)
    CoverageIncomplete,
    /// 大域因子対の証人文脈がない (v33.1 の witness 規律の atlas 版)
    CompatibilityUnwitnessed,
}

impl AtlasAbstainReason {
    pub fn as_str(self) -> &'static str {
        match self {
            AtlasAbstainReason::ChartFailed(_) => "chart_failed",
            AtlasAbstainReason::GlueInconsistent => "glue_inconsistent",
            AtlasAbstainReason::CoverageIncomplete => "coverage_incomplete",
            AtlasAbstainReason::CompatibilityUnwitnessed => "compatibility_unwitnessed",
        }
    }
}

/// atlas 裁定 (凍結) — GluedExact / 複数の全被覆整合群 = EquivalenceClassOnly /
/// Abstain。単一 chart の局所 Exact をここへ昇格する経路は存在しない (禁止変換 19)。
pub enum AtlasReading {
    GluedExact {
        local_dims: Vec<usize>,
        /// 大域因子 (glue 後) の traceless 部分代数 — 直接大域復元との照合用
        factor_subalgebras: Vec<Vec<Vec<C64>>>,
    },
    EquivalenceClassOnly {
        n_consistent_atlases: usize,
    },
    Abstain(AtlasAbstainReason),
}

impl AtlasReading {
    pub fn as_str(&self) -> &'static str {
        match self {
            AtlasReading::GluedExact { .. } => "glued_exact",
            AtlasReading::EquivalenceClassOnly { .. } => "equivalence_class_only",
            AtlasReading::Abstain(_) => "abstain",
        }
    }
}

/// 部分代数 ONB 同士の一致度 (min 側でなく対称平均でもなく — v32.3 と同じ
/// overlap fraction: 次元不一致は 0)
fn subalgebra_overlap(u: &[Vec<C64>], w: &[Vec<C64>]) -> f64 {
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
}

/// atlas の glue (凍結手順):
///   1. 各 chart を局所復元 (失敗 = Abstain(ChartFailed))
///   2. 成分対の連結 = 共有 primitive ∨ 交差 NonCommuting。連結対の部分代数
///      matching (overlap ≥ 1 − 1e-9) を検査
///   3. matching が全て整合 → 因子クラスを union-find で束ね、次元一致・被覆
///      (Π d = n)・大域証人 (クラス対ごとに共有文脈) を検査 → GluedExact
///   4. matching が破れたら: 「内部整合な chart 群」に分割し、**全被覆** (Π d = n)
///      を達成する群を数える — 2 つ以上 → EquivalenceClassOnly / それ以外 →
///      Abstain(GlueInconsistent) (宣言 chart を黙って捨てない)
pub fn recover_atlas<G: CommutationGrading>(
    net: &OperationalNet<G>,
    charts: &[ChartSpec],
    n: usize,
) -> AtlasReading {
    // 1. 局所復元
    let mut locals: Vec<ChartLocalFactorization> = Vec::new();
    for ch in charts {
        match recover_chart(net, ch, n) {
            Ok(l) => locals.push(l),
            Err(f) => return AtlasReading::Abstain(AtlasAbstainReason::ChartFailed(f)),
        }
    }
    // 成分ノード (chart, factor) の平坦化
    struct Node {
        chart: usize,
        dim: usize,
        ids: BTreeSet<u32>,
    }
    let mut nodes: Vec<Node> = Vec::new();
    for (ci, l) in locals.iter().enumerate() {
        for (d, _, ids) in &l.factors {
            nodes.push(Node {
                chart: ci,
                dim: *d,
                ids: ids.iter().map(|i| i.0).collect(),
            });
        }
    }
    let sub = |t: usize| -> &Vec<Vec<C64>> {
        // nodes[t] に対応する部分代数 ONB
        let mut idx = t;
        for l in locals.iter() {
            if idx < l.factors.len() {
                return &l.factors[idx].1;
            }
            idx -= l.factors.len();
        }
        unreachable!()
    };
    // 2. 連結と matching
    let linked = |a: &Node, b: &Node| -> bool {
        if a.chart == b.chart {
            return false;
        }
        if a.ids.intersection(&b.ids).next().is_some() {
            return true;
        }
        for &x in &a.ids {
            for &y in &b.ids {
                if net.commutator_verdict(OpId(x), OpId(y))
                    == Some(CommutatorVerdict::NonCommuting)
                {
                    return true;
                }
            }
        }
        false
    };
    let m = nodes.len();
    let mut mismatch_charts: Vec<(usize, usize)> = Vec::new();
    let mut node_parent: Vec<usize> = (0..m).collect();
    fn nfind(parent: &mut Vec<usize>, x: usize) -> usize {
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
    for a in 0..m {
        for b in (a + 1)..m {
            if !linked(&nodes[a], &nodes[b]) {
                continue;
            }
            let ov = subalgebra_overlap(sub(a), sub(b));
            if ov >= 1.0 - 1e-9 {
                let (ra, rb) = (nfind(&mut node_parent, a), nfind(&mut node_parent, b));
                if ra != rb {
                    node_parent[ra] = rb;
                }
            } else {
                mismatch_charts.push((nodes[a].chart, nodes[b].chart));
            }
        }
    }
    // 補助: chart 部分集合の全被覆整合 assembly を試みる
    let try_assemble = |chart_set: &BTreeSet<usize>| -> Option<Vec<usize>> {
        // 対象 chart 群のノードで union-find (matching 必須)
        let sel: Vec<usize> = (0..m).filter(|&t| chart_set.contains(&nodes[t].chart)).collect();
        let mut parent: Vec<usize> = (0..sel.len()).collect();
        for (ia, &a) in sel.iter().enumerate() {
            for (ib, &b) in sel.iter().enumerate().skip(ia + 1) {
                if !linked(&nodes[a], &nodes[b]) {
                    continue;
                }
                let ov = subalgebra_overlap(sub(a), sub(b));
                if ov < 1.0 - 1e-9 {
                    return None; // 群の内部に不整合
                }
                let (ra, rb) = (nfind(&mut parent, ia), nfind(&mut parent, ib));
                if ra != rb {
                    parent[ra] = rb;
                }
            }
        }
        let mut classes: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
        for (i, &t) in sel.iter().enumerate() {
            let r = nfind(&mut parent, i);
            classes.entry(r).or_default().push(t);
        }
        let mut dims = Vec::new();
        for members in classes.values() {
            let d0 = nodes[members[0]].dim;
            if members.iter().any(|&t| nodes[t].dim != d0) {
                return None;
            }
            dims.push(d0);
        }
        let prod: usize = dims.iter().product();
        if prod != n {
            return None;
        }
        dims.sort_unstable();
        Some(dims)
    };
    if mismatch_charts.is_empty() {
        // 3. 整合 — クラス束ね・被覆・証人
        let mut classes: std::collections::BTreeMap<usize, Vec<usize>> = Default::default();
        for t in 0..m {
            let r = nfind(&mut node_parent, t);
            classes.entry(r).or_default().push(t);
        }
        let mut dims = Vec::new();
        let mut class_ids: Vec<BTreeSet<u32>> = Vec::new();
        let mut class_subs: Vec<Vec<Vec<C64>>> = Vec::new();
        for members in classes.values() {
            let d0 = nodes[members[0]].dim;
            if members.iter().any(|&t| nodes[t].dim != d0) {
                return AtlasReading::Abstain(AtlasAbstainReason::GlueInconsistent);
            }
            dims.push(d0);
            let mut ids = BTreeSet::new();
            for &t in members {
                ids.extend(nodes[t].ids.iter().cloned());
            }
            class_ids.push(ids);
            class_subs.push(sub(members[0]).clone());
        }
        let prod: usize = dims.iter().product();
        if prod != n {
            return AtlasReading::Abstain(AtlasAbstainReason::CoverageIncomplete);
        }
        // 大域証人: クラス対ごとに共有文脈 (v33.1 の atlas 版)
        for i in 0..class_ids.len() {
            for j in (i + 1)..class_ids.len() {
                let witnessed = net.contexts().iter().any(|ctx| {
                    class_ids[i].iter().any(|a| ctx.contains(a))
                        && class_ids[j].iter().any(|b| ctx.contains(b))
                });
                if !witnessed {
                    return AtlasReading::Abstain(AtlasAbstainReason::CompatibilityUnwitnessed);
                }
            }
        }
        let mut order: Vec<usize> = (0..dims.len()).collect();
        order.sort_by_key(|&i| dims[i]);
        let sorted_dims: Vec<usize> = order.iter().map(|&i| dims[i]).collect();
        let sorted_subs: Vec<Vec<Vec<C64>>> = order.into_iter().map(|i| class_subs[i].clone()).collect();
        return AtlasReading::GluedExact {
            local_dims: sorted_dims,
            factor_subalgebras: sorted_subs,
        };
    }
    // 4. 不整合 — 整合 chart 群 (matching 破れの無い連結群) への分割
    let nc = charts.len();
    let mut chart_parent: Vec<usize> = (0..nc).collect();
    // 「連結かつ整合」な chart 対を束ねる
    for a in 0..m {
        for b in (a + 1)..m {
            if nodes[a].chart == nodes[b].chart || !linked(&nodes[a], &nodes[b]) {
                continue;
            }
            if subalgebra_overlap(sub(a), sub(b)) >= 1.0 - 1e-9 {
                let (ra, rb) = (
                    nfind(&mut chart_parent, nodes[a].chart),
                    nfind(&mut chart_parent, nodes[b].chart),
                );
                if ra != rb {
                    chart_parent[ra] = rb;
                }
            }
        }
    }
    let mut chart_groups: std::collections::BTreeMap<usize, BTreeSet<usize>> = Default::default();
    for c in 0..nc {
        let r = nfind(&mut chart_parent, c);
        chart_groups.entry(r).or_default().insert(c);
    }
    let full_groups = chart_groups
        .values()
        .filter(|g| try_assemble(g).is_some())
        .count();
    if full_groups >= 2 {
        AtlasReading::EquivalenceClassOnly {
            n_consistent_atlases: full_groups,
        }
    } else {
        AtlasReading::Abstain(AtlasAbstainReason::GlueInconsistent)
    }
}

// ---------------------------------------------------------------- 自己検査

/// contextual_factorization の不変条件 (v334_contextual_factorization が呼ぶ)。
/// 2 qubit (dim 4) の最小 atlas で chart 復元・glue・禁止変換 19 の名目を検査。
pub fn contextual_factorization_self_test() -> Result<(), String> {
    use crate::operational_net::{
        CertifiedCommutator, ControlGenerator, OpKind, OperatorParity, OrdinaryCommutation,
    };
    use crate::operational_net::commutator;
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
    let gens = [
        kron2(&px, &id2),
        kron2(&pz, &id2),
        kron2(&id2, &px),
        kron2(&id2, &pz),
    ];
    let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, 1e-3);
    let mut ids = Vec::new();
    for g in &gens {
        ids.push(
            net.add_primitive(PrimitiveOperation {
                kind: OpKind::Control(
                    ControlGenerator::certify(
                        g.iter().map(|c| c.re).collect(),
                        g.iter().map(|c| c.im).collect(),
                        n,
                    )
                    .unwrap(),
                ),
                parity: OperatorParity::Even,
                provenance: "contextual_self_test",
            })
            .unwrap(),
        );
    }
    for a in 0..4 {
        for b in (a + 1)..4 {
            let nu = hs_norm(&commutator(&gens[a], &gens[b], n));
            net.set_commutator(
                ids[a],
                ids[b],
                CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
            );
        }
    }
    net.add_context(&[ids[0], ids[2]]).map_err(|e| e)?;
    net.add_context(&[ids[1], ids[3]]).map_err(|e| e)?;
    // chart A = {X₁, Z₁}: 文脈が chart 内に無い → ContextIncomplete
    let ch_a = ChartSpec {
        primitive_ids: vec![ids[0], ids[1]],
    };
    match recover_chart(&net, &ch_a, n) {
        Err(ChartFailure::ContextIncomplete) => {}
        r => {
            return Err(format!(
                "chart 内文脈なしが検出されない: {:?}",
                r.map(|_| "ok").err().map(|f| f.as_str())
            ))
        }
    }
    // 全体 chart = 全 primitive: 局所復元は因子 2 つ (M₂ × M₂)
    let ch_all = ChartSpec {
        primitive_ids: ids.clone(),
    };
    let l = recover_chart(&net, &ch_all, n).map_err(|f| f.as_str().to_string())?;
    let mut dims: Vec<usize> = l.factors.iter().map(|(d, _, _)| *d).collect();
    dims.sort_unstable();
    if dims != vec![2, 2] {
        return Err(format!("全体 chart の因子が {:?}", dims));
    }
    // atlas = {全体 chart} → GluedExact [2,2]
    match recover_atlas(&net, &[ch_all], n) {
        AtlasReading::GluedExact { local_dims, .. } if local_dims == vec![2, 2] => {}
        r => return Err(format!("最小 atlas の glue が {}", r.as_str())),
    }
    // 裁定名の一意性
    let names = [
        AtlasReading::GluedExact {
            local_dims: vec![2],
            factor_subalgebras: Vec::new(),
        }
        .as_str(),
        AtlasReading::EquivalenceClassOnly {
            n_consistent_atlases: 2,
        }
        .as_str(),
        AtlasReading::Abstain(AtlasAbstainReason::GlueInconsistent).as_str(),
    ];
    for (i, a) in names.iter().enumerate() {
        for b in names.iter().skip(i + 1) {
            if a == b {
                return Err("atlas 裁定名の重複".into());
            }
        }
    }
    Ok(())
}
