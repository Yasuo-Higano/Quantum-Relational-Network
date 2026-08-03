//! v33.6 構造化スケーリング — Pauli GF(2) / Majorana quadratic backend (PROMPT/14)
//!
//! 一般 dense の *-閉包をそのまま高速化して「大型系へスケールした」とは主張しない。
//! lane を分け、各 lane は自分の scope でのみ語り、scope 外は ScopeExceeded を返す。
//! 共有契約 = `sim/src/structured_backend.rs`。
//!
//!   [S0] 契約自己検査 — structured_backend_self_test (新) + 既存 5 契約の不変
//!   [S1] **n = 3 の全セルで dense と Pauli backend の裁定が完全一致**: site
//!        [2,2,2]・entangler [2,4]・number-only Insufficient・部分 net
//!        Insufficient・singleton 文脈 unwitnessed・部分 address 超選択
//!        [(2,2),(2,2)]・パリティ超選択 [(4,1),(4,1)] — GF(2) symplectic
//!        (成分 = ω グラフ・閉包次元 = 2^{dim V}・中心 = radical) が v33.1 凍結
//!        手順の裁定を行列なしで再現する
//!   [S2] **48 qubit の証明書 (2^48 次元行列はどこにも現れない)**: site 96 本 →
//!        Exact [2×48]・+entangler → [2×46, 4]・X₄₈ 欠落 → SuperselectionSectors
//!        [(2^47, 1)×2] — 全て GF(2) rank (96×96) の厳密線形代数 (計時は合否条件に
//!        しない — 壁時計は並列負荷依存で決定性規約に反する: v34.0-B 儀式の器械訂正)
//!   [S3] **Majorana quadratic backend と dense 対応原理**: 小 N = 3 で支持分割 =
//!        dense 非可換成分・dense *-閉包次元 = 2^{2m−1} (偶 Clifford — full M_d では
//!        ない: パリティ超選択) の予言が一致。大 N = 24 (48 Majorana) で 3 ブロック
//!        so(16) 閉包 (dim 120×3)・cross hop で [16,32] へ併合 (so(32) dim 496) —
//!        2^24 は現れない
//!   [S4] **scope 規律 (禁止変換 21)**: dense はサイズバー (dim > 4096) で
//!        ScopeExceeded・非 Pauli 和 (X₁+X₂) は PauliVector に構成不能・非反対称は
//!        quadratic 資格を通らない・Pauli from_dense は真の Pauli 文字列を往復復元
//!   [S5] 封鎖の schema/文書検査 — 概念登録 + 禁止変換 21・アンカー
//!
//! 実行: cargo run --release --bin v336_structured_backend

use std::fs;
use std::path::Path;
use uft_sim::operational_net::{
    algebra_closure, commutator, hs_norm, CertifiedCommutator, ControlGenerator,
    FactorizationReading, OpId, OpKind, OperationalNet, OperatorParity, OrdinaryCommutation,
    PrimitiveOperation,
};
use uft_sim::structured_backend::*;
use uft_sim::C64;

// ---------------------------------------------------------------- dense 素子 (n = 3 照合用)

fn pauli(which: char) -> Vec<C64> {
    let (o, l) = (C64::new(0.0, 0.0), C64::new(1.0, 0.0));
    match which {
        'I' => vec![l, o, o, l],
        'X' => vec![o, l, l, o],
        'Y' => vec![o, C64::new(0.0, -1.0), C64::new(0.0, 1.0), o],
        'Z' => vec![l, o, o, C64::new(-1.0, 0.0)],
        _ => panic!("未知の Pauli"),
    }
}

fn kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
    let n = na * nb;
    let mut out = vec![C64::new(0.0, 0.0); n * n];
    for i1 in 0..na {
        for j1 in 0..na {
            for i2 in 0..nb {
                for j2 in 0..nb {
                    out[(i1 * nb + i2) * n + (j1 * nb + j2)] = a[i1 * na + j1] * b[i2 * nb + j2];
                }
            }
        }
    }
    out
}

fn op3(s: &str) -> Vec<C64> {
    let cs: Vec<char> = s.chars().collect();
    let a = kron(&pauli(cs[0]), 2, &pauli(cs[1]), 2);
    kron(&a, 4, &pauli(cs[2]), 2)
}

fn build_dense_net(
    strs: &[&str],
    tau: f64,
    contexts: &[Vec<usize>],
) -> (OperationalNet<OrdinaryCommutation>, Vec<OpId>) {
    let n = 8usize;
    let gens: Vec<Vec<C64>> = strs.iter().map(|s| op3(s)).collect();
    let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
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
                provenance: "v336_control",
            })
            .unwrap(),
        );
    }
    for (a, ga) in gens.iter().enumerate() {
        for (b, gb) in gens.iter().enumerate().skip(a + 1) {
            let nu = hs_norm(&commutator(ga, gb, n));
            net.set_commutator(
                ids[a],
                ids[b],
                CertifiedCommutator::new((nu - 1e-12).max(0.0), nu + 1e-12).unwrap(),
            );
        }
    }
    for ctx in contexts {
        let members: Vec<OpId> = ctx.iter().map(|&i| ids[i]).collect();
        net.add_context(&members).unwrap();
    }
    (net, ids)
}

fn dense_reading(strs: &[&str], contexts: &[Vec<usize>]) -> String {
    let (net, _) = build_dense_net(strs, 1e-3, contexts);
    match net.recovery_input() {
        Ok(inp) => {
            let r = inp.recover().reading;
            format!("{:?}", r)
        }
        Err(e) => format!("rejected:{}", e.as_str()),
    }
}

fn pauli_reading(strs: &[&str], contexts: &[Vec<usize>]) -> String {
    let spec = PauliNetSpec {
        n_qubits: strs[0].len(),
        ops: strs.iter().map(|s| PauliVector::from_str(s)).collect(),
        contexts: contexts
            .iter()
            .map(|c| c.iter().cloned().collect())
            .collect(),
    };
    match recover_pauli_net(&spec) {
        Ok(r) => format!("{:?}", r),
        Err(e) => format!("rejected:{}", e.as_str()),
    }
}

fn main() {
    uft_sim::self_test();
    println!("=== v33.6 構造化スケーリング — Pauli GF(2) / Majorana quadratic backend (PROMPT/14) ===\n");
    let root = if Path::new("core.schema.yml").exists() {
        "."
    } else {
        ".."
    };
    let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p));
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---- [S0] 契約自己検査 ----
    {
        let mut bad = Vec::new();
        if let Err(e) = uft_sim::structured_backend::structured_backend_self_test() {
            bad.push(format!("structured_backend_self_test: {}", e));
        }
        if let Err(e) = uft_sim::graded_recovery::graded_recovery_self_test() {
            bad.push(format!("graded_recovery_self_test: {}", e));
        }
        if let Err(e) = uft_sim::contextual_factorization::contextual_factorization_self_test() {
            bad.push(format!("contextual_factorization_self_test: {}", e));
        }
        if let Err(e) = uft_sim::operational_net::operational_net_self_test() {
            bad.push(format!("operational_net_self_test: {}", e));
        }
        if let Err(e) = uft_sim::operational_net::scope_repair_self_test() {
            bad.push(format!("scope_repair_self_test: {}", e));
        }
        check(
            "[S0] 契約自己検査 — structured_backend (新) + 既存 4 契約の不変",
            bad.is_empty(),
            if bad.is_empty() {
                "lane は scope で語り、scope 外は ScopeExceeded".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S1] n = 3 全セルで dense = Pauli backend ----
    {
        let mut bad = Vec::new();
        let cells: Vec<(&str, Vec<&str>, Vec<Vec<usize>>)> = vec![
            (
                "site [2,2,2]",
                vec!["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"],
                vec![vec![0, 2, 4], vec![1, 3, 5]],
            ),
            (
                "entangler [2,4]",
                vec!["XII", "ZII", "IXI", "IZI", "IIX", "IIZ", "XXI"],
                vec![vec![0, 2, 4, 6], vec![1, 3, 5]],
            ),
            (
                "number-only Insufficient",
                vec!["ZII", "IZI", "IIZ"],
                vec![vec![0, 1, 2]],
            ),
            (
                "部分 net Insufficient",
                vec!["XII", "ZII"],
                vec![vec![0], vec![1]],
            ),
            (
                "singleton 文脈 unwitnessed",
                vec!["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"],
                vec![vec![0], vec![1], vec![2], vec![3], vec![4], vec![5]],
            ),
            (
                "部分 address 超選択 [(2,2),(2,2)]",
                vec!["XII", "ZII", "IZI"],
                vec![vec![0, 2], vec![1, 2]],
            ),
            (
                "パリティ超選択 [(4,1),(4,1)]",
                vec!["ZII", "XXI", "IZI", "IXX", "IIZ"],
                vec![vec![0, 2, 4], vec![1, 3]],
            ),
        ];
        let mut lines = Vec::new();
        for (name, strs, ctxs) in &cells {
            let d = dense_reading(strs, ctxs);
            let p = pauli_reading(strs, ctxs);
            if d != p {
                bad.push(format!("{}: dense = {} / pauli = {}", name, d, p));
            } else {
                lines.push(format!("{} ✓", name));
            }
        }
        check(
            "[S1] n = 3 全 7 セルで dense と Pauli backend の裁定が完全一致 (Exact/Insufficient/unwitnessed/超選択)",
            bad.is_empty(),
            if bad.is_empty() {
                lines.join(" / ")
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [S2] 48 qubit の証明書 (行列なし) ----
    {
        let mut bad = Vec::new();
        let nq = 48usize;
        let mk = |pos: usize, ch: char| -> String {
            let mut s = vec!['I'; nq];
            s[pos] = ch;
            s.iter().collect()
        };
        // site 96 本
        let mut strs: Vec<String> = Vec::new();
        for i in 0..nq {
            strs.push(mk(i, 'X'));
            strs.push(mk(i, 'Z'));
        }
        let ops: Vec<PauliVector> = strs.iter().map(|s| PauliVector::from_str(s)).collect();
        let ctx_x: std::collections::BTreeSet<usize> = (0..nq).map(|i| 2 * i).collect();
        let ctx_z: std::collections::BTreeSet<usize> = (0..nq).map(|i| 2 * i + 1).collect();
        let spec = PauliNetSpec {
            n_qubits: nq,
            ops: ops.clone(),
            contexts: vec![ctx_x.clone(), ctx_z.clone()],
        };
        let r1 = recover_pauli_net(&spec);
        let want48 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2; 48],
        };
        if r1.as_ref().map(|r| r == &want48) != Ok(true) {
            bad.push(format!("site 48q の読みが {:?}", r1.map(|r| r.as_str().to_string())));
        }
        // + entangler X₀X₁
        let mut strs2 = strs.clone();
        let mut e = vec!['I'; nq];
        e[0] = 'X';
        e[1] = 'X';
        strs2.push(e.iter().collect());
        let mut ctx_x2 = ctx_x.clone();
        ctx_x2.insert(96);
        let spec2 = PauliNetSpec {
            n_qubits: nq,
            ops: strs2.iter().map(|s| PauliVector::from_str(s)).collect(),
            contexts: vec![ctx_x2, ctx_z.clone()],
        };
        let r2 = recover_pauli_net(&spec2);
        let mut dims2 = vec![2usize; 46];
        dims2.push(4);
        dims2.sort_unstable();
        let want_ent = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: dims2,
        };
        if r2.as_ref().map(|r| r == &want_ent) != Ok(true) {
            bad.push(format!("+entangler の読みが {:?}", r2.map(|r| r.as_str().to_string())));
        }
        // X₄₈ 欠落 (qubit 47 は Z のみ) → 超選択 [(2^47, 1) × 2]
        let strs3: Vec<String> = strs[..95].to_vec(); // 最後の X₄₇ を落とす → ops: X0..Z46, X47? 並びは (X_i, Z_i) — 95 本 = X47 まで含み Z47 欠落
        // 並び注意: strs = [X0, Z0, X1, Z1, ..., X47, Z47] — 95 本にすると Z47 が欠落
        let ctx_x3: std::collections::BTreeSet<usize> = (0..nq).map(|i| 2 * i).collect();
        let ctx_z3: std::collections::BTreeSet<usize> = (0..(nq - 1)).map(|i| 2 * i + 1).collect();
        let spec3 = PauliNetSpec {
            n_qubits: nq,
            ops: strs3.iter().map(|s| PauliVector::from_str(s)).collect(),
            contexts: vec![ctx_x3, ctx_z3],
        };
        let r3 = recover_pauli_net(&spec3);
        // Z₄₇ 欠落 → X₄₇ は単独成分・radical に X₄₇ → sectors 2・m = 2^47・mult 1
        let want_ss = FactorizationReading::SuperselectionSectors {
            sectors: vec![(1usize << 47, 1); 2],
        };
        if r3.as_ref().map(|r| r == &want_ss) != Ok(true) {
            bad.push(format!("Z₄₇ 欠落の読みが {:?}", r3.map(|r| r.as_str().to_string())));
        }
        // 計時は合否条件にしない (壁時計は並列負荷に依存し「並列化で結果が変わらない」
        // 規約に反する — v34.0-B の儀式 [JOBS=12] が旧版の壁時計バー超過を検出した
        // 器械訂正。規模の主張は「行列を生成しない」構成そのもので担う)
        check(
            "[S2] 48 qubit の証明書 — Exact [2×48]・+entangler [2×46,4]・Z₄₇ 欠落 → [(2^47,1)×2] (行列なし)",
            bad.is_empty(),
            "GF(2) rank (96×96) の厳密線形代数のみ — 2^48 次元の行列はどこにも現れない (dense はこの入力を ScopeExceeded で拒否する [S4])".into(),
        );
    }

    // ---- [S3] Majorana quadratic backend と dense 対応原理 ----
    {
        let mut bad = Vec::new();
        // 小 N = 3 (6 Majorana): ブロック {γ1..γ4} (NN 3 本) + {γ5,γ6} (1 本)
        let m6 = 6usize;
        let mk_hop = |i: usize, j: usize, m: usize| -> QuadraticGenerator {
            let mut a = vec![0.0f64; m * m];
            a[i * m + j] = 1.0;
            a[j * m + i] = -1.0;
            QuadraticGenerator::certify(a, m).unwrap()
        };
        let gens6 = vec![
            mk_hop(0, 1, m6),
            mk_hop(1, 2, m6),
            mk_hop(2, 3, m6),
            mk_hop(4, 5, m6),
        ];
        let QuadraticBlockReading::Blocks {
            block_majoranas,
            lie_full,
        } = recover_quadratic_blocks(&gens6);
        if !(block_majoranas == vec![2, 4] && lie_full == vec![true, true]) {
            bad.push(format!("小 N ブロック {:?} / {:?}", block_majoranas, lie_full));
        }
        // dense 対応原理: 支持分割 = dense 非可換成分・dense *-閉包 dim = 2^{2m−1}
        let gam = ["XII", "YII", "ZXI", "ZYI", "ZZX", "ZZY"];
        let bil = |a: &str, b: &str| -> Vec<C64> {
            let (ga, gb) = (op3(a), op3(b));
            uft_sim::operational_net::cmul(&ga, &gb, 8)
                .iter()
                .map(|c| C64::new(-c.im, c.re))
                .collect()
        };
        let dense_hops = vec![
            bil(gam[0], gam[1]),
            bil(gam[1], gam[2]),
            bil(gam[2], gam[3]),
            bil(gam[4], gam[5]),
        ];
        // dense 非可換成分: hop 0-1, 1-2 は非可換・hop 3 は独立
        let mut linked01 = hs_norm(&commutator(&dense_hops[0], &dense_hops[1], 8)) > 1e-9;
        linked01 &= hs_norm(&commutator(&dense_hops[1], &dense_hops[2], 8)) > 1e-9;
        let sep3 = (0..3)
            .all(|i| hs_norm(&commutator(&dense_hops[i], &dense_hops[3], 8)) < 1e-12);
        let cl_block1 = algebra_closure(&dense_hops[..3].to_vec(), 8).len();
        let cl_block2 = algebra_closure(&[dense_hops[3].clone()].to_vec(), 8).len();
        let pred1 = QuadraticBlockReading::predicted_dense_closure_dim(4);
        let pred2 = QuadraticBlockReading::predicted_dense_closure_dim(2);
        if !(linked01 && sep3 && cl_block1 == pred1 && cl_block2 == pred2) {
            bad.push(format!(
                "対応原理: 成分 {}/{}・閉包 {} (予言 {}) / {} (予言 {})",
                linked01, sep3, cl_block1, pred1, cl_block2, pred2
            ));
        }
        // 大 N = 24 (48 Majorana): 3 ブロック × 16 → so(16) dim 120・cross hop で併合
        let m48 = 48usize;
        let mut gens48 = Vec::new();
        for blk in 0..3 {
            for k in 0..15 {
                gens48.push(mk_hop(16 * blk + k, 16 * blk + k + 1, m48));
            }
        }
        let QuadraticBlockReading::Blocks {
            block_majoranas: bm48,
            lie_full: lf48,
        } = recover_quadratic_blocks(&gens48);
        if !(bm48 == vec![16, 16, 16] && lf48 == vec![true, true, true]) {
            bad.push(format!("大 N ブロック {:?} / {:?}", bm48, lf48));
        }
        let mut gens48b = gens48.clone();
        gens48b.push(mk_hop(15, 16, m48));
        let QuadraticBlockReading::Blocks {
            block_majoranas: bm48b,
            lie_full: lf48b,
        } = recover_quadratic_blocks(&gens48b);
        if !(bm48b == vec![16, 32] && lf48b == vec![true, true]) {
            bad.push(format!("cross hop 後 {:?} / {:?}", bm48b, lf48b));
        }
        // 計時は合否条件にしない ([S2] と同じ器械訂正 — 儀式の並列負荷で検出)
        check(
            "[S3] Majorana quadratic backend — 小 N 対応原理 (支持分割 = dense 成分・閉包 dim = 2^{2m−1})・大 N = 48 Majorana so(16)³ → so(32) 併合",
            bad.is_empty(),
            format!(
                "N=3: [2,4] blocks・dense 閉包 {}/{} = 予言 {}/{} (偶 Clifford — full M_d ではない = パリティ超選択の quadratic 版) / N=24: [16,16,16] so 閉包 120×3・cross hop → [16,32] (so(32) dim 496) — 2^24 は現れない",
                cl_block1, cl_block2, pred1, pred2
            ),
        );
    }

    // ---- [S4] scope 規律 (禁止変換 21) ----
    {
        let mut bad = Vec::new();
        // dense のサイズバー
        if dense_scope_guard(1usize << 48) != Err(StructuredScopeError::DimensionTooLargeForDense)
        {
            bad.push("dense が 2^48 を拒否しない".into());
        }
        if dense_scope_guard(8).is_err() {
            bad.push("dense が小次元を拒否した".into());
        }
        // 非 Pauli 和は PauliVector に構成不能
        let x1 = op3("XII");
        let x2 = op3("IXI");
        let tied: Vec<C64> = x1.iter().zip(x2.iter()).map(|(a, b)| *a + *b).collect();
        match PauliVector::from_dense(&tied, 8) {
            Err(StructuredScopeError::NotPauliString) => {}
            r => bad.push(format!("非 Pauli 和の裁定が {:?}", r.err().map(|e| e.as_str()))),
        }
        // 真の Pauli 文字列は往復復元 (反可換パターンの一致)
        let rt = PauliVector::from_dense(&op3("XIZ"), 8);
        match rt {
            Ok(p) => {
                let want = PauliVector::from_str("XIZ");
                for s in ["ZII", "XII", "IXI", "IZI", "IIX", "IIZ"] {
                    let q = PauliVector::from_str(s);
                    if p.anticommutes(&q) != want.anticommutes(&q) {
                        bad.push(format!("往復復元の反可換パターン不一致 ({})", s));
                    }
                }
            }
            Err(e) => bad.push(format!("真の Pauli が拒否された ({})", e.as_str())),
        }
        // 非反対称は quadratic 資格外
        if QuadraticGenerator::certify(vec![0.0, 1.0, 1.0, 0.0], 2).is_ok() {
            bad.push("対称行列が quadratic 資格を通った".into());
        }
        check(
            "[S4] scope 規律 — dense は dim > 4096 で ScopeExceeded・非 Pauli 和は構成不能・真の Pauli は往復復元・非反対称は拒否",
            bad.is_empty(),
            "「できない」は型と裁定で言う — 小系 dense の成果を大型一般系へ昇格しない (禁止変換 21)".into(),
        );
    }

    // ---- [S5] 封鎖の schema/文書検査 ----
    {
        let mut bad = Vec::new();
        let src = rd("sim/src/structured_backend.rs").unwrap_or_default();
        for needle in [
            "禁止変換 21",
            "ScopeExceeded",
            "DENSE_DIM_BAR",
            "NotPauliString",
            "2^{2m−1}",
            "行列を生成せず",
        ] {
            if !src.contains(needle) {
                bad.push(format!("structured_backend.rs: 「{}」が無い", needle));
            }
        }
        let schema = rd("core.schema.yml").unwrap_or_default();
        for needle in [
            "- name: PauliSymplecticBackend",
            "- name: MajoranaQuadraticBackend",
            "- name: GenericDenseBackend",
            "- name: StructuredScopeError",
        ] {
            if !schema.contains(needle) {
                bad.push(format!("core.schema.yml: 「{}」が無い", needle));
            }
        }
        if !schema
            .contains("- from: GenericDenseSmallSystemResult\n  to: LargeSystemFactorizationClaim\n  reason:")
        {
            bad.push("禁止変換 21 が未登録".into());
        }
        let doc = rd("docs/uft-v33.6.md").unwrap_or_default();
        for needle in [
            "構造化スケーリング",
            "GF(2)",
            "symplectic",
            "禁止変換 21",
            "ScopeExceeded",
            "対応原理",
            "so(16)",
            "行列なし",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v33.6.md: 「{}」が無い", needle));
            }
        }
        check(
            "[S5] 封鎖の schema/文書 — 概念登録 + 禁止変換 21・アンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "3 lane の scope が schema/文書の三点で凍結された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "structured lane が大型系の証明書を返すようになった — dense と厳密一致する小系対応原理つきで、scope 外は ScopeExceeded。残るは HOLD-9 (v34.0-A/B) と Track X (D2-R campaign layer)"
        } else {
            "**backend 契約の破れ** — structured_backend と schema/文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
