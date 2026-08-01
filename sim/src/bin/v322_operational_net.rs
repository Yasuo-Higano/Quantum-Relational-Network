//! v32.2 OperationalNet — 操作的文脈の型分離と global-algebra erasure no-go (PROMPT/13 §2)
//!
//! 第三十二期の主研究の開幕。v31.4 の E3 no-go は「状態単独では因子分解を選定でき
//! ない — 因子分解は操作代数が運ぶ」を定理化したが、旧 OperationalAlgebra は文字列
//! リストの雛形で、積・随伴・閉包・可換子証明書を持たなかった。本版はこれを
//! **marked family of subalgebras (OperationalNet)** に精密化し、**第二の no-go**
//! を機械化する:
//!
//! > **目標定理 A (global-algebra erasure no-go)**: 異なる tensor factorization に
//! > 対応する二つの operational generator family が、同一の大域閉包代数 B(H) を
//! > 生成し得る。したがって大域代数の同型類だけから因子分解は識別不能であり、
//! > 因子分解を運ぶのは閉包ではなく marking (どの primitive が独立に address
//! > 可能で、どの部分集合が操作文脈をなすか) である。
//!
//!   [N0] operational_net_self_test — 役割 4 型の資格審査・証明書 3 値裁定・
//!        Ordinary lane の odd 拒否・文脈の可換子証明書要求・閉包/commutant 素子
//!   [N1] **erasure no-go の機械実証**: site 生成族 {X_i, Z_i} と DFT 共役族
//!        {V X_i V†, V Z_i V†} (V = DFT_8) が**ともに閉包 = M_8 (複素次元 64)**。
//!        一方 su(2) 部分代数の対応は存在しない (subspace overlap ≪ 1) —
//!        同一の閉包・異なる因子分解
//!   [N2] **marking は因子分解を運ぶ** (v32.3 recovery の最小 preview): 両 net の
//!        非可換グラフ連結成分 → 各成分閉包 = M_2 (dim 4)・成分間可換・joint = M_8
//!        → local_dims [2,2,2] — 同じ閉包から net ごとに別の因子分解が復元される
//!   [N3] **役割の型分離の意味論**: 測定 (作用素系) は積に閉じない (span 3 →
//!        積で 4)・準備は凸結合で閉じ積は状態でない (tr ≠ 1)・介入は Lie bracket
//!        で閉じる ([X,Z]/i ∈ su(2))
//!   [N4] **可換子証明書の 3 値裁定**: 明確な非可換 (‖[X₁,Z₁]‖ = 2√8)・明確な可換
//!        (ν = 0)・**閾値を跨ぐ区間 → Abstain** (辺の強制禁止 — 文脈構成も拒否)
//!   [N5] **Z2 grading の罠 (JW 弦の幾何誤読)**: 独立 3 モードの site-local odd
//!        (γ₁, γ₃, γ₅) は ordinary 可換子では完全グラフ K₃ に見える (‖[γ,γ']‖ =
//!        2√8 ≠ 0) が、graded (反可換子) では空グラフ (真の独立)。Ordinary net は
//!        odd を構成時に拒否 (型遮断)・parity-even 双線形は ordinary で安全
//!   [N6] 型レベル封鎖の source 検査 — 禁止変換 11 (GlobalClosure → OperationalNet/
//!        因子分解) の impl From 不在・役割型間の From 不在・odd 拒否ゲートの実在・
//!        qrn_core/readout_contract の封鎖不変
//!   [N7] schema/文書アンカー — core.schema.yml の概念登録 + 禁止変換 11・
//!        uft-v32.2.md の定理文
//!
//! 実行: cargo run --release --bin v322_operational_net

use std::fs;
use std::path::Path;
use uft_sim::operational_net::*;
use uft_sim::C64;

// ---------------------------------------------------------------- Pauli / kron / DFT 素子

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

/// 3 qubit (dim 8) の site 演算子: 文字列 "XIZ" 等
fn op3(s: &str) -> Vec<C64> {
    let cs: Vec<char> = s.chars().collect();
    let a = kron(&pauli(cs[0]), 2, &pauli(cs[1]), 2);
    kron(&a, 4, &pauli(cs[2]), 2)
}

fn site_pauli(p: char, site: usize) -> Vec<C64> {
    let mut s = ['I', 'I', 'I'];
    s[site] = p;
    op3(&s.iter().collect::<String>())
}

/// DFT_8: F[j][k] = e^{2πi jk/8}/√8 (決定的な非局所 unitary)
fn dft8() -> Vec<C64> {
    let n = 8;
    let inv = 1.0 / (n as f64).sqrt();
    let mut f = vec![C64::new(0.0, 0.0); n * n];
    for j in 0..n {
        for k in 0..n {
            f[j * n + k] =
                C64::expi(2.0 * std::f64::consts::PI * (j * k) as f64 / n as f64).scale(inv);
        }
    }
    f
}

fn conj_by(v: &[C64], a: &[C64], n: usize) -> Vec<C64> {
    cmul(&cmul(v, a, n), &cdag(v, n), n)
}

/// su(2) 部分空間 (正規直交 3 基底) 同士の overlap fraction —
/// 同一部分空間なら 1・直交なら 0
fn su2_overlap(a: &[Vec<C64>], b: &[Vec<C64>]) -> f64 {
    let mut acc = 0.0;
    for x in b {
        for y in a {
            acc += hs_inner(y, x).norm2();
        }
    }
    acc / 3.0
}

fn main() {
    uft_sim::self_test();
    println!(
        "=== v32.2 OperationalNet — 型分離と global-algebra erasure no-go (PROMPT/13 §2) ===\n"
    );
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
    let n = 8usize;
    let tau = 1e-6; // 可換子閾値 (凍結)

    // ---- [N0] 契約の自己検査 ----
    {
        let r = operational_net_self_test();
        check(
            "[N0] operational_net_self_test — 役割 4 型・証明書 3 値・odd 拒否・文脈要求・閉包/commutant",
            r.is_ok(),
            match &r {
                Ok(()) => "役割は型・可換性は証明書・grading は構成時ゲート".into(),
                Err(e) => e.clone(),
            },
        );
    }

    // site 族と DFT 共役族
    let v = dft8();
    let site_gens: Vec<Vec<C64>> = (0..3)
        .flat_map(|i| vec![site_pauli('X', i), site_pauli('Z', i)])
        .collect();
    let mode_gens: Vec<Vec<C64>> = site_gens.iter().map(|g| conj_by(&v, g, n)).collect();

    // ---- [N1] erasure no-go ----
    {
        let cl_site = closure_of(&site_gens, n);
        let cl_mode = closure_of(&mode_gens, n);
        // su(2) 部分代数の対応の不在 (因子分解の相違の証人)
        let norm = 1.0 / (n as f64).sqrt(); // ‖P‖_F = √8
        let su2_site: Vec<Vec<Vec<C64>>> = (0..3)
            .map(|i| {
                ['X', 'Y', 'Z']
                    .iter()
                    .map(|&p| {
                        site_pauli(p, i)
                            .iter()
                            .map(|c| c.scale(norm))
                            .collect::<Vec<C64>>()
                    })
                    .collect()
            })
            .collect();
        let su2_mode: Vec<Vec<Vec<C64>>> = (0..3)
            .map(|i| {
                ['X', 'Y', 'Z']
                    .iter()
                    .map(|&p| {
                        conj_by(&v, &site_pauli(p, i), n)
                            .iter()
                            .map(|c| c.scale(norm))
                            .collect::<Vec<C64>>()
                    })
                    .collect()
            })
            .collect();
        let mut max_ov = 0.0f64;
        for sb in &su2_mode {
            for sa in &su2_site {
                max_ov = max_ov.max(su2_overlap(sa, sb));
            }
        }
        // 自己整合: site 対 site の対応は対角で 1
        let self_ov = su2_overlap(&su2_site[0], &su2_site[0]);
        let ok = cl_site.is_full()
            && cl_mode.is_full()
            && cl_site.dim_algebra == 64
            && cl_mode.dim_algebra == 64
            && (self_ov - 1.0).abs() < 1e-12
            && max_ov < 0.9;
        check(
            "[N1] erasure no-go — site 族と DFT 共役族の閉包はともに M_8 (dim 64)・su(2) 対応は不在",
            ok,
            format!(
                "閉包 dim: site = {} / mode = {} (full = {})・su(2) subspace overlap: 対角自己 = {:.12}・site×mode 最大 = {:.6} (≪ 1 = 因子分解は別物)",
                cl_site.dim_algebra,
                cl_mode.dim_algebra,
                cl_site.is_full() && cl_mode.is_full(),
                self_ov,
                max_ov
            ),
        );
    }

    // ---- net の構築 (exact ノルムの証明書) ----
    let build_net = |gens: &[Vec<C64>]| -> OperationalNet<OrdinaryCommutation> {
        let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
        let mut ids = Vec::new();
        for g in gens {
            let re: Vec<f64> = g.iter().map(|c| c.re).collect();
            let im: Vec<f64> = g.iter().map(|c| c.im).collect();
            let p = PrimitiveOperation {
                kind: OpKind::Control(ControlGenerator::certify(re, im, n).unwrap()),
                parity: OperatorParity::Even,
                provenance: "v322_site_or_mode_control",
            };
            ids.push(net.add_primitive(p).unwrap());
        }
        for (a, ga) in gens.iter().enumerate() {
            for (b, gb) in gens.iter().enumerate().skip(a + 1) {
                let nu = hs_norm(&commutator(ga, gb, n));
                let cert = CertifiedCommutator::new(
                    (nu - 1e-12).max(0.0),
                    nu + 1e-12,
                )
                .unwrap();
                net.set_commutator(ids[a], ids[b], cert);
            }
        }
        net
    };
    let net_site = build_net(&site_gens);
    let net_mode = build_net(&mode_gens);

    // ---- [N2] marking が因子分解を運ぶ (recovery preview) ----
    {
        let mut bad = Vec::new();
        let mut dims_pair = Vec::new();
        for (name, net, gens) in [
            ("site", &net_site, &site_gens),
            ("mode", &net_mode, &mode_gens),
        ] {
            match net.noncommutation_components() {
                Ok(comps) => {
                    if comps.len() != 3 || comps.iter().any(|c| c.len() != 2) {
                        bad.push(format!("{}: 成分構造が {:?}", name, comps));
                        continue;
                    }
                    let mut dims = Vec::new();
                    for comp in &comps {
                        let sub: Vec<Vec<C64>> =
                            comp.iter().map(|&i| gens[i as usize].clone()).collect();
                        let d = algebra_closure(&sub, n).len();
                        if d != 4 {
                            bad.push(format!("{}: 成分閉包 dim = {} ≠ 4 (M_2)", name, d));
                        }
                        dims.push((d as f64).sqrt() as usize);
                    }
                    dims_pair.push(dims);
                }
                Err(e) => bad.push(format!("{}: 成分分解が棄却 {:?}", name, e)),
            }
        }
        let ok = bad.is_empty()
            && dims_pair.len() == 2
            && dims_pair.iter().all(|d| d == &vec![2, 2, 2]);
        check(
            "[N2] marking → 因子分解 (preview) — 両 net とも成分 3 × M_2・local_dims [2,2,2]・同一閉包から別の因子分解",
            ok,
            format!(
                "site/mode の成分閉包次元 → local_dims = {:?} (読みは ExactUpToLocalUnitaryAndPermutation の候補 — 本裁定は v32.3)",
                dims_pair
            ),
        );
    }

    // ---- [N3] 役割の型分離の意味論 ----
    {
        let mut bad = Vec::new();
        // (a) 測定: 作用素系は積に閉じない
        let id8: Vec<C64> = {
            let mut m = vec![C64::new(0.0, 0.0); n * n];
            for i in 0..n {
                m[i * n + i] = C64::new(1.0, 0.0);
            }
            m
        };
        let n1: Vec<C64> = id8
            .iter()
            .zip(site_pauli('Z', 0).iter())
            .map(|(i, z)| (*i - *z).scale(0.5))
            .collect();
        let n2: Vec<C64> = id8
            .iter()
            .zip(site_pauli('Z', 1).iter())
            .map(|(i, z)| (*i - *z).scale(0.5))
            .collect();
        let mut span: Vec<Vec<C64>> = Vec::new();
        for m in [&id8, &n1, &n2] {
            push_ortho(&mut span, m, 1e-9);
        }
        let d_span = span.len();
        let prod = cmul(&n1, &n2, n);
        let grew = push_ortho(&mut span, &prod, 1e-9);
        if !(d_span == 3 && grew && span.len() == 4) {
            bad.push(format!("測定 span {} → 積後 {} (期待 3 → 4)", d_span, span.len()));
        }
        // effect の資格
        let eff = MeasurementEffect::certify(
            n1.iter().map(|c| c.re).collect(),
            n1.iter().map(|c| c.im).collect(),
            n,
        );
        if eff.is_err() {
            bad.push("n₁ が effect 資格を通らない".into());
        }
        // (b) 準備: 凸結合は状態・積は状態でない
        let mut rho_a = vec![0.0; n * n];
        rho_a[0] = 1.0; // |000⟩⟨000|
        let rho_b = vec![
            {
                let mut d = vec![0.0; n * n];
                for i in 0..n {
                    d[i * n + i] = 1.0 / n as f64;
                }
                d
            },
            vec![0.0; n * n],
        ];
        let pa = Preparation::certify(rho_a.clone(), vec![0.0; n * n], n).unwrap();
        let pb = Preparation::certify(rho_b[0].clone(), rho_b[1].clone(), n).unwrap();
        if Preparation::mix(&pa, &pb, 0.5).is_err() {
            bad.push("凸結合が準備の資格を通らない".into());
        }
        // 行列積 ρ_a ρ_b は trace 1/8 — 準備を名乗れない
        let prod_re: Vec<f64> = {
            let a: Vec<C64> = rho_a.iter().map(|&r| C64::new(r, 0.0)).collect();
            let b: Vec<C64> = rho_b[0].iter().map(|&r| C64::new(r, 0.0)).collect();
            cmul(&a, &b, n).iter().map(|c| c.re).collect()
        };
        if Preparation::certify(prod_re, vec![0.0; n * n], n).is_ok() {
            bad.push("状態の積が準備を名乗れた (tr ≠ 1 のはず)".into());
        }
        // (c) 介入: Lie bracket で閉じる — [X₁, Z₁]/i = −2Y₁ ∈ su(2)₁
        let br = commutator(&site_pauli('X', 0), &site_pauli('Z', 0), n);
        let br_over_i: Vec<C64> = br.iter().map(|c| C64::new(c.im, -c.re)).collect(); // (1/i)·c = −i·c
        let h = ControlGenerator::certify(
            br_over_i.iter().map(|c| c.re).collect(),
            br_over_i.iter().map(|c| c.im).collect(),
            n,
        );
        if h.is_err() {
            bad.push("[X,Z]/i がエルミート資格を通らない".into());
        }
        let mut su2_basis: Vec<Vec<C64>> = Vec::new();
        for p in ['X', 'Y', 'Z'] {
            push_ortho(&mut su2_basis, &site_pauli(p, 0), 1e-9);
        }
        let grew2 = push_ortho(&mut su2_basis, &br_over_i, 1e-9);
        if grew2 {
            bad.push("[X₁,Z₁]/i が su(2)₁ の外に出た".into());
        }
        check(
            "[N3] 役割の型分離 — 測定は積に閉じない (3→4)・準備は凸で閉じ積は tr≠1・介入は Lie で閉じる",
            bad.is_empty(),
            if bad.is_empty() {
                "preparation / control / measurement は同じ数学的型ではない — 型と意味論が一致".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [N4] 可換子証明書の 3 値裁定 ----
    {
        let mut bad = Vec::new();
        let x1 = site_pauli('X', 0);
        let z1 = site_pauli('Z', 0);
        let x2 = site_pauli('X', 1);
        let nu_edge = hs_norm(&commutator(&x1, &z1, n));
        let nu_comm = hs_norm(&commutator(&x1, &x2, n));
        let c_edge = CertifiedCommutator::new(nu_edge - 1e-12, nu_edge + 1e-12).unwrap();
        let c_comm = CertifiedCommutator::new(0.0, nu_comm + 1e-12).unwrap();
        if c_edge.verdict(tau) != CommutatorVerdict::NonCommuting {
            bad.push(format!("明確な非可換 (ν = {:.6}) が NonCommuting にならない", nu_edge));
        }
        if c_comm.verdict(tau) != CommutatorVerdict::Commuting {
            bad.push(format!("明確な可換 (ν = {:.1e}) が Commuting にならない", nu_comm));
        }
        // 閾値を跨ぐ区間: 真値 0.9τ・誤差 ±0.6τ → [0.3τ, 1.5τ]
        let c_marginal = CertifiedCommutator::new(0.3 * tau, 1.5 * tau).unwrap();
        if c_marginal.verdict(tau) != CommutatorVerdict::Abstain {
            bad.push("跨ぎ区間が Abstain にならない".into());
        }
        // Abstain 対を含む文脈は構成できない・成分分解も棄却される
        let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
        let mk = |g: &[C64]| PrimitiveOperation {
            kind: OpKind::Control(
                ControlGenerator::certify(
                    g.iter().map(|c| c.re).collect(),
                    g.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            parity: OperatorParity::Even,
            provenance: "v322_n4",
        };
        let a = net.add_primitive(mk(&x1)).unwrap();
        let b = net.add_primitive(mk(&x2)).unwrap();
        net.set_commutator(a, b, c_marginal);
        let ctx_refused = net.add_context(&[a, b]).is_err();
        let comp_refused = matches!(
            net.noncommutation_components(),
            Err(FactorizationAbstainReason::CommutatorMarginStraddled)
        );
        if !(ctx_refused && comp_refused) {
            bad.push("Abstain 対で文脈/成分分解が拒否されない".into());
        }
        check(
            "[N4] 証明書 3 値 — 非可換 2√8・可換 0・跨ぎは Abstain (文脈構成・成分分解とも拒否 = 辺の強制禁止)",
            bad.is_empty(),
            if bad.is_empty() {
                format!(
                    "ν_edge = {:.6} (= 2√8 = {:.6})・ν_comm = {:.1e}・跨ぎ [0.3τ, 1.5τ] → abstain",
                    nu_edge,
                    2.0 * (8.0f64).sqrt(),
                    nu_comm
                )
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [N5] Z2 grading の罠 (JW 弦の幾何誤読) ----
    {
        let mut bad = Vec::new();
        // 3 モードの JW Majorana (odd): γ₁ = XII, γ₃ = ZXI, γ₅ = ZZX
        let g1 = op3("XII");
        let g3 = op3("ZXI");
        let g5 = op3("ZZX");
        let odd = [("mode1", &g1), ("mode2", &g3), ("mode3", &g5)];
        let mut ord_edges = 0;
        let mut graded_edges = 0;
        let mut min_ord = f64::INFINITY;
        let mut max_graded = 0.0f64;
        for (i, (_, a)) in odd.iter().enumerate() {
            for (_, b) in odd.iter().skip(i + 1) {
                let no = hs_norm(&commutator(a, b, n));
                let ng = hs_norm(&anticommutator(a, b, n));
                if no > tau {
                    ord_edges += 1;
                }
                if ng > tau {
                    graded_edges += 1;
                }
                min_ord = min_ord.min(no);
                max_graded = max_graded.max(ng);
            }
        }
        if !(ord_edges == 3 && graded_edges == 0) {
            bad.push(format!(
                "ordinary 辺 {} (期待 3 = K₃ 捏造)・graded 辺 {} (期待 0)",
                ord_edges, graded_edges
            ));
        }
        // 型遮断: Ordinary net は odd を受けない
        let mut net_o: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
        let odd_prim = PrimitiveOperation {
            kind: OpKind::Control(
                ControlGenerator::certify(
                    g1.iter().map(|c| c.re).collect(),
                    g1.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            parity: OperatorParity::Odd,
            provenance: "v322_jw_majorana",
        };
        if net_o.add_primitive(odd_prim.clone()).is_ok() {
            bad.push("Ordinary net が odd を受理した".into());
        }
        let mut net_g: OperationalNet<FermionicZ2Graded> = OperationalNet::new(n, tau);
        if net_g.add_primitive(odd_prim).is_err() {
            bad.push("Z2 graded net が odd を拒否した".into());
        }
        // parity-even 双線形 b_p = iγ_{2p-1}γ_{2p} は ordinary で可換 (安全)
        let g2 = op3("YII");
        let g6 = op3("ZZY");
        let b1: Vec<C64> = cmul(&g1, &g2, n).iter().map(|c| C64::new(-c.im, c.re)).collect();
        let b3: Vec<C64> = cmul(&g5, &g6, n).iter().map(|c| C64::new(-c.im, c.re)).collect();
        let even_comm = hs_norm(&commutator(&b1, &b3, n));
        if even_comm > 1e-12 {
            bad.push(format!("parity-even 双線形が ordinary で非可換 ({:.1e})", even_comm));
        }
        check(
            "[N5] JW 幾何誤読の遮断 — odd は ordinary で K₃ 捏造 (‖[γ,γ']‖ = 2√8)・graded で空・型は構成時拒否",
            bad.is_empty(),
            format!(
                "ordinary 辺 {}/3 (min ‖[·,·]‖ = {:.6})・graded 辺 {}/3 (max ‖{{·,·}}‖ = {:.1e})・even 双線形の ordinary 可換 = {:.1e}",
                ord_edges, min_ord, graded_edges, max_graded, even_comm
            ),
        );
    }

    // ---- [N6] 型レベル封鎖の source 検査 ----
    {
        let mut bad = Vec::new();
        let src = rd("sim/src/operational_net.rs").unwrap_or_default();
        if src.is_empty() {
            bad.push("operational_net.rs が読めない".to_string());
        }
        for forbidden in [
            "impl From<GlobalClosure",
            "impl From<Preparation",
            "impl From<ControlGenerator",
            "impl From<MeasurementEffect",
            "impl From<DriftGenerator",
        ] {
            if src.contains(forbidden) {
                bad.push(format!("禁止 impl {} が存在する", forbidden));
            }
        }
        for needle in [
            "ACCEPTS_ODD",
            "OrdinaryCommutation の net は odd (fermionic) primitive を受け付けない",
            "禁止変換 11",
        ] {
            if !src.contains(needle) {
                bad.push(format!("operational_net.rs: 「{}」が無い", needle));
            }
        }
        let rc = rd("sim/src/readout_contract.rs").unwrap_or_default();
        for needle in ["v32.2 で後継が確定", "operational_net::OperationalNet"] {
            if !rc.contains(needle) {
                bad.push(format!("readout_contract.rs: 「{}」が無い", needle));
            }
        }
        let lib = rd("sim/src/lib.rs").unwrap_or_default();
        if !lib.contains("pub mod operational_net;") {
            bad.push("lib.rs に operational_net が登録されていない".into());
        }
        if let Err(e) = uft_sim::qrn_core::qrn_core_self_test() {
            bad.push(format!("qrn_core_self_test: {}", e));
        }
        if let Err(e) = uft_sim::readout_contract::readout_contract_self_test() {
            bad.push(format!("readout_contract_self_test: {}", e));
        }
        check(
            "[N6] 封鎖の source 検査 — 禁止変換 11 の From 不在・役割型間 From 不在・odd ゲート実在・既存契約不変",
            bad.is_empty(),
            if bad.is_empty() {
                "閉包から marking へ戻る経路はコンパイル不能 — qrn_core/readout_contract の封鎖も不変".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [N7] schema/文書アンカー ----
    {
        let mut bad = Vec::new();
        let schema = rd("core.schema.yml").unwrap_or_default();
        for needle in [
            "- name: OperationalNet",
            "- name: GlobalClosure",
            "- name: PrimitiveOperation",
            "- name: CertifiedCommutator",
        ] {
            if !schema.contains(needle) {
                bad.push(format!("core.schema.yml: 「{}」が無い", needle));
            }
        }
        if !schema.contains("- from: GlobalClosure\n  to: OperationalNet\n  reason:") {
            bad.push("禁止変換 11 (GlobalClosure → OperationalNet) が未登録".into());
        }
        let doc = rd("docs/uft-v32.2.md").unwrap_or_default();
        for needle in [
            "global-algebra erasure no-go",
            "marked family of subalgebras",
            "禁止変換 11",
            "Jordan–Wigner",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v32.2.md: 「{}」が無い", needle));
            }
        }
        check(
            "[N7] schema/文書 — 概念登録 + 禁止変換 11 + 定理文のアンカー",
            bad.is_empty(),
            if bad.is_empty() {
                "no-go が型・schema・文書の三点で凍結された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "大域閉包は因子分解を消す — 因子分解を運ぶのは OperationalNet (marking) であり、閉包から戻る経路は型に存在しない"
        } else {
            "**契約の破れ** — operational_net と schema/文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
