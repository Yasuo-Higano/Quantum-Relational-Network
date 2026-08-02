// resource_profile — Resource-Filtered OperationalNet と factorization profile (v33.3, PROMPT/14)
//
// 可アクセス性は resource budget に依存する。したがって v33.3 の中心対象は
// 「単一因子分解」ではなく **profile**:
//
//     OperationalFactorizationProfile: ResourceBudget ↦ FactorizationReading
//
// である。典型的には 資源不足 → Abstain / 局所制御可能 → [2,2,2] / entangler
// 可能 → [2,4] / 完全制御可能 → [8] (global) と読みが動く — **最終 budget だけを
// 見ると v32.2 の erasure no-go に戻る** (頂の読みは経路を消す)。
//
// 設計の核:
//   1. **通常の persistent homology と同一視しない**: 因子分解の読みは budget の
//      poset 上で分裂も併合も起こり得る (比較不能な budget 対が同じ dims で別の
//      gauge orbit を持つ — v333 [R2])。最初から「barcode」と呼ばず、有限 poset 上の
//      constructible profile として定義する。写像と安定性定理が成立した後にのみ
//      zigzag persistence 等へ進む (PROMPT/14)。
//   2. **昇格規則 (凍結)**: 昇格可能な局所性は「単一閾値で一瞬だけ出現した因子分解」
//      ではなく「budget perturbation に対して一定領域で同じ gauge orbit を保つ
//      因子分解」— 本版の凍結規則は **領域内に比較可能な対 (chain ≥ 2) が存在する
//      こと** (**禁止変換 17**: transient な読みを stable locality に昇格しない)。
//      grid の頂・縁で領域が切れる読みは「調査 grid 相対で transient」と正直に
//      記録する (grid の拡張だけが解除できる)。
//   3. **成分半順序をスカラーへ潰さない** (**禁止変換 18**): 恣意的な重み付き和に
//      よる全順序化は accessibility の集合を変え、読みを反転させ得る (v333 [R4] が
//      裁定反転の具体例を機械記録)。許されるのは成分ごとの狭義単調再パラメータ化
//      のみ — このとき profile は点ごとに不変 (command 再パラメータ化不変性)。
//   4. 予算は資格ではなく filter の座標: AccessibleOperation の資格 (出自・
//      addressability) は sha256 結束で不変のまま、予算だけを単位変換できる
//      (rebudgeted)。資格のない操作は budget をいくら積んでも accessible に
//      ならない (禁止変換 14 は profile 上でも維持)。
//
// 一次ソース: docs/uft-v33.3.md / core.schema.yml (概念 + 禁止変換 17/18)。
// 整合は v333_resource_profile が機械検査する。

use crate::laboratory_interface::{
    AccessibleOperation, AccessibleOperationalNet, IndependentAddressabilityCertificate,
    ResourceBudget,
};
use crate::operational_net::{
    commutator, hs_norm, same_gauge_orbit, CertifiedCommutator, CommutationGrading,
    FactorizationAbstainReason, FactorizationReading, OpId, RecoveryInputRejection,
};
use crate::C64;
use std::marker::PhantomData;

// ---------------------------------------------------------------- profile の点

/// budget 1 点での profile の値 — 読み (裁定) と、Exact のときの gauge orbit 照合材料
pub enum ProfilePoint {
    /// この budget で accessible な操作が 1 つも無い (資源不足の正直な記録)
    NoAccessibleOperations,
    /// 復元入力の構成時拒否 (被覆不全等) — Abstain とは区別して記録する
    InputRejected(RecoveryInputRejection),
    /// 復元の読み (三値裁定) + 成分 traceless 部分代数 (orbit 照合用)
    Reading {
        reading: FactorizationReading,
        component_subalgebras: Vec<Vec<Vec<C64>>>,
    },
}

impl ProfilePoint {
    pub fn as_str(&self) -> &'static str {
        match self {
            ProfilePoint::NoAccessibleOperations => "no_accessible_operations",
            ProfilePoint::InputRejected(_) => "input_rejected",
            ProfilePoint::Reading { reading, .. } => reading.as_str(),
        }
    }
    /// 同一クラスか — 読みの値が等しく、Exact なら gauge orbit も matching すること
    pub fn same_class(&self, other: &ProfilePoint) -> bool {
        match (self, other) {
            (ProfilePoint::NoAccessibleOperations, ProfilePoint::NoAccessibleOperations) => true,
            (ProfilePoint::InputRejected(a), ProfilePoint::InputRejected(b)) => a == b,
            (
                ProfilePoint::Reading {
                    reading: ra,
                    component_subalgebras: sa,
                },
                ProfilePoint::Reading {
                    reading: rb,
                    component_subalgebras: sb,
                },
            ) => {
                if ra != rb {
                    return false;
                }
                match ra {
                    FactorizationReading::ExactUpToLocalUnitaryAndPermutation { .. } => {
                        same_gauge_orbit(sa, sb).0
                    }
                    _ => true,
                }
            }
            _ => false,
        }
    }
}

// ---------------------------------------------------------------- Resource-Filtered interface

/// 資格つき操作 (コスト = AccessibleOperation の budget) と文脈レシピの族。
/// budget b での accessible 部分 interface = {op : cost(op) ≤ b (成分ごと)} から
/// AccessibleOperationalNet を組み、v33.1 修復入口で読む。
pub struct ResourceFilteredInterface<G: CommutationGrading> {
    dim: usize,
    threshold: f64,
    ops: Vec<AccessibleOperation>,
    matrices: Vec<Vec<C64>>,
    context_recipes: Vec<(Vec<usize>, IndependentAddressabilityCertificate)>,
    _g: PhantomData<G>,
}

impl<G: CommutationGrading> ResourceFilteredInterface<G> {
    pub fn new(dim: usize, threshold: f64) -> Self {
        ResourceFilteredInterface {
            dim,
            threshold,
            ops: Vec::new(),
            matrices: Vec::new(),
            context_recipes: Vec::new(),
            _g: PhantomData,
        }
    }

    /// 操作の登録 — AccessibleOperation (資格つき) のみ。コストは op.budget()。
    pub fn add_operation(&mut self, op: AccessibleOperation) -> usize {
        let (re, im, d) = op.kind().matrix();
        let m: Vec<C64> = (0..d * d).map(|k| C64::new(re[k], im[k])).collect();
        self.matrices.push(m);
        self.ops.push(op);
        self.ops.len() - 1
    }

    /// 文脈レシピの登録 — budget での制限は「レシピ ∩ accessible」で行う
    pub fn add_context_recipe(
        &mut self,
        members: Vec<usize>,
        cert: IndependentAddressabilityCertificate,
    ) -> Result<usize, String> {
        for &i in &members {
            if !cert.covers_target(&self.matrices[i], self.dim) {
                return Err(format!("証明書が member {} を結束していない", i));
            }
        }
        self.context_recipes.push((members, cert));
        Ok(self.context_recipes.len() - 1)
    }

    /// budget b で accessible な操作の添字 (成分ごとの ≤)
    pub fn accessible_indices(&self, b: &ResourceBudget) -> Vec<usize> {
        (0..self.ops.len())
            .filter(|&i| self.ops[i].budget().componentwise_le(b))
            .collect()
    }

    /// budget b での読み — accessible 部分 interface から net を組み v33.1 入口で復元
    pub fn reading_at(&self, b: &ResourceBudget) -> ProfilePoint {
        let idx = self.accessible_indices(b);
        if idx.is_empty() {
            return ProfilePoint::NoAccessibleOperations;
        }
        let mut net: AccessibleOperationalNet<G> =
            AccessibleOperationalNet::new(self.dim, self.threshold);
        let mut local: Vec<OpId> = Vec::new();
        for &i in &idx {
            local.push(
                net.admit(self.ops[i].clone())
                    .expect("資格つき操作の admit が失敗"),
            );
        }
        for (a, &ia) in idx.iter().enumerate() {
            for (bb, &ib) in idx.iter().enumerate().skip(a + 1) {
                let nu = hs_norm(&commutator(&self.matrices[ia], &self.matrices[ib], self.dim));
                net.set_commutator(
                    local[a],
                    local[bb],
                    CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
                );
            }
        }
        for (members, cert) in &self.context_recipes {
            let restricted: Vec<OpId> = members
                .iter()
                .filter_map(|m| idx.iter().position(|&i| i == *m).map(|p| local[p]))
                .collect();
            if restricted.is_empty() {
                continue;
            }
            net.add_control_context(&restricted, cert.clone())
                .expect("制限文脈の登録が失敗 (レシピは可換集合であること)");
        }
        match net.recover() {
            Ok(detail) => ProfilePoint::Reading {
                reading: detail.reading,
                component_subalgebras: detail.component_subalgebras,
            },
            Err(rej) => ProfilePoint::InputRejected(rej),
        }
    }

    /// 有限 grid 上の profile
    pub fn profile_over(&self, grid: &[ResourceBudget]) -> OperationalFactorizationProfile {
        OperationalFactorizationProfile {
            budgets: grid.to_vec(),
            points: grid.iter().map(|b| self.reading_at(b)).collect(),
        }
    }
}

// ---------------------------------------------------------------- profile と昇格規則

/// 有限 poset (grid) 上の constructible profile — barcode ではない
pub struct OperationalFactorizationProfile {
    pub budgets: Vec<ResourceBudget>,
    pub points: Vec<ProfilePoint>,
}

/// profile 上の同値クラス (読み + orbit) の領域と安定性
pub struct ProfileClass {
    /// 代表点の添字
    pub representative: usize,
    /// 領域 (このクラスに属する grid 点の添字)
    pub region: Vec<usize>,
    /// 凍結昇格規則: 領域内に比較可能な対 (chain ≥ 2) が存在する
    pub stable: bool,
}

impl OperationalFactorizationProfile {
    /// クラス分解 (読み + orbit の同値で grid 点を分類)
    pub fn classes(&self) -> Vec<ProfileClass> {
        let n = self.points.len();
        let mut assigned = vec![false; n];
        let mut out = Vec::new();
        for i in 0..n {
            if assigned[i] {
                continue;
            }
            let mut region = vec![i];
            assigned[i] = true;
            for j in (i + 1)..n {
                if !assigned[j] && self.points[i].same_class(&self.points[j]) {
                    region.push(j);
                    assigned[j] = true;
                }
            }
            let stable = region.iter().any(|&a| {
                region.iter().any(|&b| {
                    a != b
                        && (self.budgets[a].componentwise_le(&self.budgets[b])
                            || self.budgets[b].componentwise_le(&self.budgets[a]))
                })
            });
            out.push(ProfileClass {
                representative: i,
                region,
                stable,
            });
        }
        out
    }

    /// transient (昇格不能) な点の添字 — 禁止変換 17 の執行点
    pub fn transient_points(&self) -> Vec<usize> {
        self.classes()
            .into_iter()
            .filter(|c| !c.stable)
            .flat_map(|c| c.region)
            .collect()
    }
}

// ---------------------------------------------------------------- 自己検査

/// resource_profile の不変条件 (v333_resource_profile が呼ぶ)
pub fn resource_profile_self_test() -> Result<(), String> {
    use crate::laboratory_interface::{
        certify_addressability, OperationOrigin, ResourceBudget,
    };
    use crate::operational_net::{ControlGenerator, OpKind, OperatorParity, OrdinaryCommutation};
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
    let gv: Vec<Vec<C64>> = gens.to_vec();
    let cert = certify_addressability(&gv, &gv, n, 0.5, 0.1).map_err(|e| e.as_str().to_string())?;
    let cost = |d: f64| ResourceBudget::certify(1.0, 1.0, 1.0, d, 1e-9).unwrap();
    let mut iface: ResourceFilteredInterface<OrdinaryCommutation> =
        ResourceFilteredInterface::new(n, 1e-3);
    for g in &gv {
        let op = AccessibleOperation::certify(
            OpKind::Control(
                ControlGenerator::certify(
                    g.iter().map(|c| c.re).collect(),
                    g.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            OperatorParity::Even,
            OperationOrigin::DirectlyCalibrated(cert.clone()),
            cert.clone(),
            cost(1.0),
        )
        .map_err(|e| e.as_str().to_string())?;
        iface.add_operation(op);
    }
    iface
        .add_context_recipe(vec![0, 2], cert.clone())
        .map_err(|e| e)?;
    iface
        .add_context_recipe(vec![1, 3], cert.clone())
        .map_err(|e| e)?;
    // 資源不足 → NoAccessibleOperations / 十分 → Exact [2,2]
    let p0 = iface.reading_at(&cost(0.5));
    if !matches!(p0, ProfilePoint::NoAccessibleOperations) {
        return Err(format!("資源不足の点が {} (期待 no_accessible_operations)", p0.as_str()));
    }
    let p1 = iface.reading_at(&cost(1.0));
    match &p1 {
        ProfilePoint::Reading { reading, .. }
            if *reading
                == (FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
                    local_dims: vec![2, 2],
                }) => {}
        _ => return Err(format!("十分予算の点が {} (期待 Exact [2,2])", p1.as_str())),
    }
    // 同一クラス判定 (自明: 同じ点)
    if !p1.same_class(&iface.reading_at(&cost(2.0))) {
        return Err("同一読みの 2 点が同一クラスにならない".into());
    }
    // 昇格規則: chain {1.0, 2.0} で stable・単点 {0.5} は transient
    let prof = iface.profile_over(&[cost(0.5), cost(1.0), cost(2.0)]);
    let classes = prof.classes();
    let exact_stable = classes.iter().any(|c| {
        matches!(prof.points[c.representative], ProfilePoint::Reading { .. }) && c.stable
    });
    let empty_transient = classes.iter().any(|c| {
        matches!(
            prof.points[c.representative],
            ProfilePoint::NoAccessibleOperations
        ) && !c.stable
    });
    if !(exact_stable && empty_transient) {
        return Err("昇格規則 (chain ≥ 2) の判定が誤り".into());
    }
    // Abstain 理由名の整合 (profile 語彙の一意性)
    let names = [
        ProfilePoint::NoAccessibleOperations.as_str(),
        ProfilePoint::InputRejected(RecoveryInputRejection::NoDeclaredContexts).as_str(),
        FactorizationReading::Abstain(FactorizationAbstainReason::InsufficientOperationalGenerators)
            .as_str(),
    ];
    for (i, a) in names.iter().enumerate() {
        for b in names.iter().skip(i + 1) {
            if a == b {
                return Err("profile 点の語彙が重複".into());
            }
        }
    }
    Ok(())
}
