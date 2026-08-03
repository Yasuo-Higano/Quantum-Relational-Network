//! v34.4 sector-aware complete finite factorization enumerator (PROMPT/15 §5)
//!
//! 背景: FollowUp の FAC-001 判定は「代数不変量 (可換子・中心) の診断は再現、
//! 任意の有限次元代数に対する complete tensor-factor candidate search は未提供。
//! 非自明な center がある場合の全空間 tensor 昇格の禁止も確認」— 本版はその
//! enumerator を証明書つきで実装する (`sim/src/factorization_enumerator.rs`)。
//!
//! 手順 (与えられた *-部分代数の族から):
//!   1. joint 閉包と中心 → 最小中心射影の列挙 (Lagrange 補間 — v32.3 kernel 継承)
//!   2. 各 sector の Wedderburn 証明書: dim A_α = n²・dim A′_α = m²・n·m = d_α・
//!      積 span(A_α·A′_α) = d_α² (A ∨ A′ = B(H_α))・二重可換子 A″ = A
//!      (同型 A_α ≅ M_{n_α} ⊗ I_{m_α} の存在は証明書 + 標準構造定理 [C0])
//!   3. multiplicity m_α と simple factor n_α の分離
//!   4. 族の可換 simple 部分集合 + simple な補因子 (joint commutant) から
//!      候補列挙 — 積 span = d² の証明書つき
//!   5. 局所 unitary × 成分置換の witness (traceless 部分空間の overlap
//!      matching — OCS-1.0 §F3 と同一意味論・バー 0.9)
//!   6. 非同値候補が複数なら候補集合 (tie-break 禁止)
//!
//! 出力 6 型: UniqueFactorization / FactorizationCandidateSet /
//! SectorwiseFactorization / IncompletePrimitiveSet /
//! NontrivialCenterObstruction / ScopeExceeded — ScopeExceeded は正答。
//!
//! 検証: [W1] M₂⊗M₃ → Unique [2,3] (可換子 = 相方の閉包) / [W2] 3 qubit →
//! Unique [2,2,2] / [W3] multiplicity {a⊕a} → Unique [2,2] + Wedderburn 証明書
//! n = m = 2 / [W4] M₂⊕M₃ → Sectorwise [(2,1),(3,1)]・大域 tensor 要求は
//! Obstruction / [W4b] {a⊕a⊕b} → Sectorwise [(3,1),(2,2)] (sector 内多重度) /
//! [W5] site vs CNOT 共役の 2 bipartition → CandidateSet{2} (overlap 1/3) /
//! [W6] number op のみ (abelian) → IncompletePrimitiveSet / [W7] d = 128 →
//! ScopeExceeded / [W8] SWAP 置換 witness (同一 orbit) vs site×bell (非同値)。

use uft_sim::factorization_enumerator::{
    certify_sector, enumerate_candidates, factorization_enumerator_self_test,
    same_candidate_orbit, sectorwise_analysis, EnumeratorReading,
};
use uft_sim::self_test;
use uft_sim::C64;

// ---------------------------------------------------------------- 行列の構成

fn zeros(n: usize) -> Vec<C64> {
    vec![C64::new(0.0, 0.0); n * n]
}

/// kron (複素)
fn kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
    let n = na * nb;
    let mut c = zeros(n);
    for i in 0..na {
        for j in 0..na {
            for k in 0..nb {
                for l in 0..nb {
                    c[(i * nb + k) * n + (j * nb + l)] = a[i * na + j] * b[k * nb + l];
                }
            }
        }
    }
    c
}

fn eye(n: usize) -> Vec<C64> {
    let mut m = zeros(n);
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

fn re2(m: [f64; 4]) -> Vec<C64> {
    m.iter().map(|&x| C64::new(x, 0.0)).collect()
}

fn px() -> Vec<C64> {
    re2([0.0, 1.0, 1.0, 0.0])
}
fn pz() -> Vec<C64> {
    re2([1.0, 0.0, 0.0, -1.0])
}

/// qutrit の非可換生成子 (shift の実部/虚部 + clock 対角) — M₃ を生成
fn qutrit_gens() -> Vec<Vec<C64>> {
    let mut s = zeros(3); // shift |j⟩ → |j+1⟩
    s[1 * 3 + 0] = C64::new(1.0, 0.0);
    s[2 * 3 + 1] = C64::new(1.0, 0.0);
    s[0 * 3 + 2] = C64::new(1.0, 0.0);
    let mut a = zeros(3);
    let mut b = zeros(3);
    for i in 0..3 {
        for j in 0..3 {
            let v = s[i * 3 + j];
            let vt = s[j * 3 + i];
            a[i * 3 + j] = v + vt; // S + S†
            b[i * 3 + j] = C64::new(0.0, 1.0) * (v - vt); // i(S − S†)
        }
    }
    let mut d = zeros(3);
    d[0] = C64::new(2.0, 0.0);
    d[4] = C64::new(-1.0, 0.0);
    d[8] = C64::new(-1.0, 0.0);
    vec![a, b, d]
}

/// ブロック直和への埋め込み: blocks[i] = (開始 index, ブロック次元)、
/// ops[i] = そのブロックに置く行列 (None は 0)
fn block_embed(d: usize, blocks: &[(usize, usize)], ops: &[Option<&[C64]>]) -> Vec<C64> {
    let mut m = zeros(d);
    for (bi, &(off, bd)) in blocks.iter().enumerate() {
        if let Some(op) = ops[bi] {
            for i in 0..bd {
                for j in 0..bd {
                    m[(off + i) * d + (off + j)] = op[i * bd + j];
                }
            }
        }
    }
    m
}

/// site i の単一 qubit 演算子 (n_qubits 系)
fn site_op(op: &[C64], site: usize, n_qubits: usize) -> Vec<C64> {
    let mut m = vec![C64::new(1.0, 0.0)];
    let mut cur = 1usize;
    for s in 0..n_qubits {
        let next = if s == site {
            kron(&m, cur, op, 2)
        } else {
            kron(&m, cur, &eye(2), 2)
        };
        cur *= 2;
        m = next;
    }
    m
}

/// CNOT (制御 0, 標的 1) の 4×4
fn cnot() -> Vec<C64> {
    let mut m = zeros(4);
    m[0] = C64::new(1.0, 0.0);
    m[5] = C64::new(1.0, 0.0);
    m[2 * 4 + 3] = C64::new(1.0, 0.0);
    m[3 * 4 + 2] = C64::new(1.0, 0.0);
    m
}

/// SWAP の 4×4
fn swap2() -> Vec<C64> {
    let mut m = zeros(4);
    m[0] = C64::new(1.0, 0.0);
    m[1 * 4 + 2] = C64::new(1.0, 0.0);
    m[2 * 4 + 1] = C64::new(1.0, 0.0);
    m[15] = C64::new(1.0, 0.0);
    m
}

fn conj_by(u: &[C64], x: &[C64], n: usize) -> Vec<C64> {
    // U X U†
    let mut ux = zeros(n);
    for i in 0..n {
        for k in 0..n {
            let uik = u[i * n + k];
            if uik.norm2() == 0.0 {
                continue;
            }
            for j in 0..n {
                ux[i * n + j] = ux[i * n + j] + uik * x[k * n + j];
            }
        }
    }
    let mut out = zeros(n);
    for i in 0..n {
        for j in 0..n {
            let mut s = C64::new(0.0, 0.0);
            for k in 0..n {
                s = s + ux[i * n + k] * u[j * n + k].conj();
            }
            out[i * n + j] = s;
        }
    }
    out
}

fn main() {
    self_test();
    factorization_enumerator_self_test().expect("enumerator self test");
    println!("=== v34.4 sector-aware factorization candidate enumerator (PROMPT/15 §5) ===\n");
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!(
            "  [{}] {}  {}",
            if ok { "PASS" } else { "FAIL" },
            name,
            detail
        );
        if !ok {
            nfail += 1;
        }
    };

    // ---------------- [W1] M₂ ⊗ M₃ ----------------
    {
        let d = 6usize;
        let a_gens = vec![kron(&px(), 2, &eye(3), 3), kron(&pz(), 2, &eye(3), 3)];
        let b_gens: Vec<Vec<C64>> = qutrit_gens()
            .iter()
            .map(|g| kron(&eye(2), 2, g, 3))
            .collect();
        let fam = vec![a_gens, b_gens];
        let r = enumerate_candidates(&fam, d, false);
        let ok = match &r {
            EnumeratorReading::UniqueFactorization { local_dims, .. } => {
                let mut s = local_dims.clone();
                s.sort();
                s == vec![2, 3]
            }
            _ => false,
        };
        check(
            "[W1] M₂⊗M₃ (d=6): UniqueFactorization [2,3]",
            ok,
            format!("verdict = {}", r.as_str()),
        );
    }

    // ---------------- [W2]/[W8] 3 qubit と置換 witness ----------------
    let (w2_components, w8_swapped) = {
        let d = 8usize;
        let fam: Vec<Vec<Vec<C64>>> = (0..3)
            .map(|s| vec![site_op(&px(), s, 3), site_op(&pz(), s, 3)])
            .collect();
        let r = enumerate_candidates(&fam, d, false);
        let (ok, comps) = match r {
            EnumeratorReading::UniqueFactorization {
                ref local_dims,
                ref components,
            } => {
                let mut s = local_dims.clone();
                s.sort();
                (s == vec![2, 2, 2], components.clone())
            }
            _ => (false, Vec::new()),
        };
        check(
            "[W2] 3 qubit site 族 (d=8): UniqueFactorization [2,2,2]",
            ok,
            format!("成分 {} 個", comps.len()),
        );
        // SWAP₁₂ ⊗ I 共役の族 — 同じ候補の別表現 (orbit 一致 witness)
        let u = kron(&swap2(), 4, &eye(2), 2);
        let fam_sw: Vec<Vec<Vec<C64>>> = fam
            .iter()
            .map(|gens| gens.iter().map(|g| conj_by(&u, g, d)).collect())
            .collect();
        let r2 = enumerate_candidates(&fam_sw, d, false);
        let comps2 = match r2 {
            EnumeratorReading::UniqueFactorization { ref components, .. } => components.clone(),
            _ => Vec::new(),
        };
        (comps, comps2)
    };
    {
        let (same, ov) = same_candidate_orbit(&w2_components, &w8_swapped, 0.9);
        check(
            "[W8a] SWAP 置換 witness: 共役族の候補は同一 orbit (成分置換で matching)",
            same && ov > 1.0 - 1e-9,
            format!("min overlap (最良置換) = {:.12}", ov),
        );
    }

    // ---------------- [W3] multiplicity {a ⊕ a} ⊂ M₄ ----------------
    {
        let d = 4usize;
        let blocks = [(0usize, 2usize), (2, 2)];
        let a_gens = vec![
            block_embed(d, &blocks, &[Some(&px()), Some(&px())]),
            block_embed(d, &blocks, &[Some(&pz()), Some(&pz())]),
        ];
        let cert = certify_sector(&a_gens, d);
        let fam = vec![a_gens];
        let r = enumerate_candidates(&fam, d, false);
        let ok_r = match &r {
            EnumeratorReading::UniqueFactorization { local_dims, .. } => {
                let mut s = local_dims.clone();
                s.sort();
                s == vec![2, 2]
            }
            _ => false,
        };
        check(
            "[W3] multiplicity {a⊕a} ⊂ M₄: Wedderburn 証明書 n=2, m=2 + Unique [2,2]",
            cert.certified() && cert.simple_dim == 2 && cert.multiplicity == 2 && ok_r,
            format!(
                "dim A = n² ({}), dim A′ = m² ({}), 積 span = d² {}, A″=A {}, verdict = {}",
                cert.simple_dim * cert.simple_dim,
                cert.multiplicity * cert.multiplicity,
                cert.product_span_full,
                cert.double_commutant_ok,
                r.as_str()
            ),
        );
    }

    // ---------------- [W4] M₂ ⊕ M₃ (中心非自明) ----------------
    {
        let d = 5usize;
        let blocks = [(0usize, 2usize), (2, 3)];
        let q = qutrit_gens();
        let gens = vec![
            block_embed(d, &blocks, &[Some(&px()), None]),
            block_embed(d, &blocks, &[Some(&pz()), None]),
            block_embed(d, &blocks, &[None, Some(&q[0])]),
            block_embed(d, &blocks, &[None, Some(&q[1])]),
            block_embed(d, &blocks, &[None, Some(&q[2])]),
        ];
        let sectors = sectorwise_analysis(&gens, d);
        let fam = vec![gens];
        let r = enumerate_candidates(&fam, d, false);
        let r_demand = enumerate_candidates(&fam, d, true);
        let ok_sec = match &sectors {
            Some(s) => {
                s.len() == 2
                    && s.iter().all(|c| c.certified())
                    && s[0].sector_dim == 2
                    && s[0].simple_dim == 2
                    && s[0].multiplicity == 1
                    && s[1].sector_dim == 3
                    && s[1].simple_dim == 3
                    && s[1].multiplicity == 1
            }
            None => false,
        };
        check(
            "[W4] M₂⊕M₃ (d=5): Sectorwise [(2,1),(3,1)] 全証明書 PASS・tensor 要求は Obstruction",
            ok_sec
                && matches!(r, EnumeratorReading::SectorwiseFactorization { .. })
                && matches!(
                    r_demand,
                    EnumeratorReading::NontrivialCenterObstruction { n_sectors: 2 }
                ),
            format!("verdict = {} / demand → {}", r.as_str(), r_demand.as_str()),
        );
    }

    // ---------------- [W4b] {a ⊕ a ⊕ b} (sector 内 multiplicity) ----------------
    {
        let d = 7usize;
        let blocks = [(0usize, 2usize), (2, 2), (4, 3)];
        let q = qutrit_gens();
        let gens = vec![
            block_embed(d, &blocks, &[Some(&px()), Some(&px()), None]),
            block_embed(d, &blocks, &[Some(&pz()), Some(&pz()), None]),
            block_embed(d, &blocks, &[None, None, Some(&q[0])]),
            block_embed(d, &blocks, &[None, None, Some(&q[1])]),
            block_embed(d, &blocks, &[None, None, Some(&q[2])]),
        ];
        let sectors = sectorwise_analysis(&gens, d);
        let ok = match &sectors {
            Some(s) => {
                s.len() == 2
                    && s.iter().all(|c| c.certified())
                    && s[0].sector_dim == 3
                    && s[0].simple_dim == 3
                    && s[0].multiplicity == 1
                    && s[1].sector_dim == 4
                    && s[1].simple_dim == 2
                    && s[1].multiplicity == 2
            }
            None => false,
        };
        check(
            "[W4b] {a⊕a⊕b} (d=7): sector 内 multiplicity — [(3, m=1), (2, m=2)] 証明書 PASS",
            ok,
            match &sectors {
                Some(s) => format!(
                    "sectors = {:?}",
                    s.iter()
                        .map(|c| (c.sector_dim, c.simple_dim, c.multiplicity))
                        .collect::<Vec<_>>()
                ),
                None => "中心が見えていない".into(),
            },
        );
    }

    // ---------------- [W5] site vs CNOT 共役 — 候補集合 ----------------
    let (site_cand, bell_cand) = {
        let d = 4usize;
        let a_site = vec![site_op(&px(), 0, 2), site_op(&pz(), 0, 2)];
        let b_site = vec![site_op(&px(), 1, 2), site_op(&pz(), 1, 2)];
        let u = cnot();
        let a_bell: Vec<Vec<C64>> = a_site.iter().map(|g| conj_by(&u, g, d)).collect();
        let b_bell: Vec<Vec<C64>> = b_site.iter().map(|g| conj_by(&u, g, d)).collect();
        let fam = vec![a_site, b_site, a_bell, b_bell];
        let r = enumerate_candidates(&fam, d, false);
        let (ok, cands) = match r {
            EnumeratorReading::FactorizationCandidateSet {
                ref candidate_dims,
                ref candidates,
            } => (
                candidate_dims.len() == 2
                    && candidate_dims.iter().all(|c| {
                        let mut s = c.clone();
                        s.sort();
                        s == vec![2, 2]
                    }),
                candidates.clone(),
            ),
            _ => (false, Vec::new()),
        };
        check(
            "[W5] site vs CNOT 共役 (d=4): FactorizationCandidateSet{2} — どちらも [2,2]・tie-break なし",
            ok,
            format!("候補 {} 個", cands.len()),
        );
        if cands.len() == 2 {
            (cands[0].clone(), cands[1].clone())
        } else {
            (Vec::new(), Vec::new())
        }
    };
    {
        let (same, ov) = same_candidate_orbit(&site_cand, &bell_cand, 0.9);
        check(
            "[W8b] 非同値 witness: site 候補と bell 候補は orbit 不一致 (overlap ≪ 0.9)",
            !same && ov < 0.5,
            format!("min overlap (最良置換) = {:.4}", ov),
        );
    }

    // ---------------- [W6] abelian 族 → IncompletePrimitiveSet ----------------
    {
        let d = 4usize;
        let mut n1 = zeros(d);
        n1[1 * d + 1] = C64::new(1.0, 0.0);
        n1[3 * d + 3] = C64::new(1.0, 0.0);
        let mut n2 = zeros(d);
        n2[2 * d + 2] = C64::new(1.0, 0.0);
        n2[3 * d + 3] = C64::new(1.0, 0.0);
        let fam = vec![vec![n1, n2]];
        let r = enumerate_candidates(&fam, d, false);
        check(
            "[W6] number op のみ (abelian 閉包): IncompletePrimitiveSet — rank-1 sector の皮を被せない",
            matches!(r, EnumeratorReading::IncompletePrimitiveSet),
            format!("verdict = {}", r.as_str()),
        );
    }

    // ---------------- [W7] ScopeExceeded ----------------
    {
        let r = enumerate_candidates(&[], 128, false);
        check(
            "[W7] d = 128 > バー 64: ScopeExceeded (試行しない — 正答)",
            matches!(r, EnumeratorReading::ScopeExceeded),
            format!("verdict = {}", r.as_str()),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "complete finite factorization candidate search が証明書つきで立った —\n       中心・sector・multiplicity を分離し、非同値候補は集合のまま返し (tie-break\n       禁止)、abelian 族と scope 超過は正直に拒否する。FollowUp FAC-001 の\n       「未提供」への応答。"
        } else {
            "**enumerator の破れ** — 証明書・候補列挙・witness を修復せよ"
        }
    );
    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
