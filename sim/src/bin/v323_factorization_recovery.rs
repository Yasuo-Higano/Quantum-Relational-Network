//! v32.3 Factorization recovery — marked operational recovery と三値裁定 (PROMPT/13 §2)
//!
//! v32.2 の erasure no-go (閉包は marking を消す) の**裏返しの正側**:
//!
//! > **目標定理 B (marked operational recovery)**: primitive generator の非可換
//! > グラフ (α∼β ⟺ ‖[O_α,O_β]‖ > 0, 証明書つき) について、(1) 異なる未知ノードの
//! > 局所操作は可換 (2) 各ノード内の非可換グラフは連結 (3) 各連結成分の *-閉包は
//! > full matrix factor M_{d_i} (4) 全成分の join が対象 sector の B(H) を生成
//! > (5) 非可換マージン > 観測誤差 — のとき、連結成分と閉包から H ≅ ⊗_i H_i を
//! > **局所 unitary × ノード置換まで**回復できる。中心が非自明なら tensor product
//! > を強制せず **superselection sector** H ≅ ⊕_α (C^{m_α} ⊗ C^{n_α}) を返す。
//!
//! 裁定は operational_net::FactorizationReading (v32.2 で先凍結):
//!   一つの許容 gauge orbit → ExactUpToLocalUnitaryAndPermutation /
//!   複数の非同値 orbit → EquivalenceClassOnly / 判定不能 → Abstain (理由 4 種)。
//!
//!   [F0] **exact oracle lane**: site net / DFT₈ 共役 net ((C²)⊗³) → Exact [2,2,2]・
//!        qutrit×qubit net (C⁶) → Exact [2,3] — 閉包 full・成分 = factor・Π d_i = n
//!   [F1] **不可能セルの正棄却**: number operator のみ ({Z₁,Z₂,Z₃} — 可換 joint) →
//!        Abstain(InsufficientOperationalGenerators) / 部分 net ({X₁,Z₁} — 中心自明
//!        だが full でない) → 同棄却 (未 address 自由度の存在)
//!   [F2] **中心非自明 → superselection (tensor 強制の禁止)**: 電荷つき 2 sector net
//!        (C² ⊕ C²⊗C²) → [(2,1),(2,2)] / 部分 address net {X₁, Z₁, Z₂} (qubit₁ 完全
//!        制御 + qubit₂ は古典ラベル Z₂ のみ + qubit₃ 未 address) → [(2,2),(2,2)] —
//!        **測定しかできない軸は超選択ラベルに・操作が届かない自由度は多重度 n_α に
//!        現れる**。なお {X₁, Z₁Z₂} の生成代数は次元 4 の因子 (中心自明・非局所
//!        符号化の 1 qubit) で、正しい読みは Insufficient — [F1] と同族
//!   [F3] **gauge orbit の裁定**: 局所 unitary × 置換で結ばれた 2 net → 同一 orbit
//!        (成分部分代数の overlap = 1 の完全 matching) / site net と DFT net →
//!        matching 不在 (max overlap 0.618) → **EquivalenceClassOnly** (無制約の
//!        tie-breaker で 1 つを選ばない — v31.4 疎性負制御の教訓)
//!   [F4] **noise abstention**: 可換子証明書を σ ノイズ区間で構成 — σ = 1e-6 は
//!        Exact 復元・σ = 5e-4 は可換対の区間が閾値を跨ぎ
//!        Abstain(CommutatorMarginStraddled) (辺の強制なし)
//!   [F5] **graded scope (フェルミオン)**: odd primitive は Ordinary lane が構成時
//!        拒否 (v32.2)。parity-even 完全 net (Majorana 双線形の path) の復元は
//!        tensor でなく **[(4,1),(4,1)] = パリティ超選択の機械発見** — 偶部分代数
//!        (dim 32) の中心 {I, Γ} が総パリティを運ぶ
//!   [F6] 文書アンカー — uft-v32.3.md の定理文・凍結決定手順
//!
//! 決定手順 (凍結): 成分 (Abstain 跨ぎで棄却) → joint 閉包 → 可換なら Insufficient →
//! 中心 → 自明: full ∧ 成分 factor ∧ Π d_i = n → Exact / full でない → Insufficient →
//! 非自明: 中心射影 (Lagrange 補間) で sector 分割 → 各 sector の制限代数 = M_{m_α}
//! (m_α² 次元) ∧ n_α = dim_α/m_α 整数 → SuperselectionSectors / 整数条件の破れ →
//! Abstain(ComponentNotFactor) (有限次元 *-代数は Wedderburn 分解を持つため、これは
//! 数値縮退の guard であって期待される出力ではない)。
//!
//! 実行: cargo run --release --bin v323_factorization_recovery

use std::fs;
use std::path::Path;
use uft_sim::operational_net::*;
use uft_sim::{Rng, C64};

// ---------------------------------------------------------------- Pauli / kron / DFT (v322 と同一素子)

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

fn ident(n: usize) -> Vec<C64> {
    let mut m = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

/// エルミート行列の固有値 (2n×2n 実対称埋め込み — 各固有値 2 重)
fn herm_evals(m: &[C64], n: usize) -> Vec<f64> {
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
    let (evals, _) = uft_sim::jacobi_eigh(&big, d);
    evals
}

// ---------------------------------------------------------------- 中心・sector 射影の素子

/// span(basis) 内で全生成子と可換な部分空間 (= 中心, basis が閉包のとき) の
/// **エルミート**正規直交基底を返す
fn center_basis(basis: &[Vec<C64>], gens: &[Vec<C64>], n: usize) -> Vec<Vec<C64>> {
    let d = basis.len();
    let dim_r = 2 * d;
    // 実座標 t (< d: basis[t], ≥ d: i·basis[t-d]) ごとに像 ⊕_g [M_t, g] を作り Gram
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
    let (evals, vecs) = uft_sim::jacobi_eigh(&gram, dim_r);
    let emax = evals.iter().cloned().fold(0.0f64, f64::max).max(1e-300);
    // 零空間ベクトル → 行列 → エルミート分解 → 正規直交化
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
                let d = *a - *b; // 反エルミート
                C64::new(d.im * 0.5, -d.re * 0.5) // /(2i)
            })
            .collect();
        // dust guard (v33.0-A 設計走行で発見・統一適用): 共役射影の数値塵
        // (‖候補‖ ≈ 0) を正規化して基底に混入させない
        if hs_norm(&h1) > 1e-9 {
            push_ortho(&mut out, &h1, 1e-8);
        }
        if hs_norm(&h2) > 1e-9 {
            push_ortho(&mut out, &h2, 1e-8);
        }
    }
    out
}

/// 中心射影の族 (Lagrange 補間): T = Σ √(k+2)·H_k の固有値クラスタごとに
/// P_α = Π_{β≠α} (T − λ_β)/(λ_α − λ_β)
fn central_projectors(center: &[Vec<C64>], n: usize) -> Option<Vec<Vec<C64>>> {
    let mut t = vec![C64::new(0.0, 0.0); n * n];
    for (k, h) in center.iter().enumerate() {
        let w = ((k + 2) as f64).sqrt();
        for (ti, hi) in t.iter_mut().zip(h.iter()) {
            *ti = *ti + hi.scale(w);
        }
    }
    let evals = herm_evals(&t, n); // 各固有値 2 重 (実埋め込み)
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
        let mut p = ident(n);
        for (b, &lb) in distinct.iter().enumerate() {
            if a == b {
                continue;
            }
            // p ← p · (T − λ_b)/(λ_a − λ_b)
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
        // 検証: 冪等・エルミート
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
    // ΣP = I
    let mut s = vec![C64::new(0.0, 0.0); n * n];
    for p in &projs {
        for (si, pi) in s.iter_mut().zip(p.iter()) {
            *si = *si + *pi;
        }
    }
    let idm = ident(n);
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

// ---------------------------------------------------------------- 復元器 (凍結決定手順)

struct RecoveryDetail {
    reading: FactorizationReading,
    /// Exact のとき: 成分ごとの traceless 部分代数 ONB (gauge orbit 照合用)
    component_subalgebras: Vec<Vec<Vec<C64>>>,
}

/// 目標定理 B の復元器 — 入力は marked net (可換子証明書つき) と生成子行列
fn recover_factorization<G: CommutationGrading>(
    net: &OperationalNet<G>,
    gens: &[Vec<C64>],
    n: usize,
) -> RecoveryDetail {
    let abstain = |r: FactorizationAbstainReason| RecoveryDetail {
        reading: FactorizationReading::Abstain(r),
        component_subalgebras: Vec::new(),
    };
    // 1. 非可換グラフの連結成分 (Abstain 対があれば棄却)
    let comps = match net.noncommutation_components() {
        Ok(c) => c,
        Err(r) => return abstain(r),
    };
    // 2. joint 閉包
    let joint = algebra_closure(gens, n);
    // 3. joint が可換なら操作資源不足 (number operator のみ等)
    let mut commutative = true;
    'outer: for a in gens {
        for b in gens {
            if hs_norm(&commutator(a, b, n)) > 1e-9 {
                commutative = false;
                break 'outer;
            }
        }
    }
    if commutative {
        return abstain(FactorizationAbstainReason::InsufficientOperationalGenerators);
    }
    // 4. 中心
    let center = center_basis(&joint, gens, n);
    if center.is_empty() {
        return abstain(FactorizationAbstainReason::ComponentNotFactor);
    }
    if center.len() == 1 {
        // 中心自明: full ∧ 各成分 factor ∧ Π d_i = n → Exact
        if joint.len() != n * n {
            // 未 address の自由度が残る (M_d ⊗ I の類)
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
            let comp_center = center_basis(&cl, &sub, n);
            if comp_center.len() != 1 {
                return abstain(FactorizationAbstainReason::ComponentNotFactor);
            }
            dims.push(d);
            // traceless 部分代数 ONB (gauge orbit 照合用)
            let idn = ident(n);
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
        return RecoveryDetail {
            reading: FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
                local_dims: sorted_dims,
            },
            component_subalgebras: subalgebras,
        };
    }
    // 5. 中心非自明: sector 分割 (tensor 強制の禁止)
    let projs = match central_projectors(&center, n) {
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
        // sector 制限代数 {P b P} の次元 = m²
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
    RecoveryDetail {
        reading: FactorizationReading::SuperselectionSectors { sectors },
        component_subalgebras: Vec::new(),
    }
}

/// gauge orbit 照合: 成分 traceless 部分代数の集合が置換で overlap ≈ 1 に
/// matching できるか (局所 unitary は部分代数を集合として不変に保つ)
fn same_gauge_orbit(a: &[Vec<Vec<C64>>], b: &[Vec<Vec<C64>>]) -> (bool, f64) {
    if a.len() != b.len() {
        return (false, 0.0);
    }
    let k = a.len();
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
    // 全置換 (k ≤ 3 想定)
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
        // 次の置換 (Heap 法の簡易 — 辞書順)
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

// ---------------------------------------------------------------- net 構築の小道具

/// exact ノルム (± noise) の証明書つき Ordinary net を作る。
/// sigma > 0 のとき ν̂ = |ν + σg|・区間 = [max(ν̂ − 6σ, 0), ν̂ + 6σ]。
fn build_net(
    gens: &[Vec<C64>],
    n: usize,
    tau: f64,
    sigma: f64,
    rng: &mut Rng,
) -> OperationalNet<OrdinaryCommutation> {
    let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, tau);
    let mut ids = Vec::new();
    for g in gens {
        let p = PrimitiveOperation {
            kind: OpKind::Control(
                ControlGenerator::certify(
                    g.iter().map(|c| c.re).collect(),
                    g.iter().map(|c| c.im).collect(),
                    n,
                )
                .unwrap(),
            ),
            parity: OperatorParity::Even,
            provenance: "v323_control",
        };
        ids.push(net.add_primitive(p).unwrap());
    }
    for (a, ga) in gens.iter().enumerate() {
        for (b, gb) in gens.iter().enumerate().skip(a + 1) {
            let nu = hs_norm(&commutator(ga, gb, n));
            let (lo, hi) = if sigma > 0.0 {
                let nu_hat = (nu + sigma * rng.gauss()).abs();
                ((nu_hat - 6.0 * sigma).max(0.0), nu_hat + 6.0 * sigma)
            } else {
                ((nu - 1e-12).max(0.0), nu + 1e-12)
            };
            net.set_commutator(ids[a], ids[b], CertifiedCommutator::new(lo, hi).unwrap());
        }
    }
    net
}

fn main() {
    uft_sim::self_test();
    println!("=== v32.3 Factorization recovery — 目標定理 B と三値裁定 (PROMPT/13 §2) ===\n");
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
    let tau = 1e-3; // 可換子閾値 (凍結 — dim ≤ 8 の Pauli 尺度で真の辺は ≥ 2√2)
    let mut rng = Rng::new(32301);

    let site_gens = |n_qubits: usize| -> Vec<Vec<C64>> {
        (0..n_qubits)
            .flat_map(|i| {
                let mut s = vec!['I'; 3];
                s[i] = 'X';
                let x = op3(&s.iter().collect::<String>());
                s[i] = 'Z';
                let z = op3(&s.iter().collect::<String>());
                vec![x, z]
            })
            .collect()
    };

    // ---- [F0] exact oracle lane ----
    let site = site_gens(3);
    let v = dft8();
    let mode: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&v, g, 8)).collect();
    let det_site = recover_factorization(&build_net(&site, 8, tau, 0.0, &mut rng), &site, 8);
    let det_mode = recover_factorization(&build_net(&mode, 8, tau, 0.0, &mut rng), &mode, 8);
    {
        // qutrit × qubit (C⁶): {path₃, diag₃} ⊗ I₂ と I₃ ⊗ {X, Z}
        let n6 = 6;
        let mut a3 = vec![C64::new(0.0, 0.0); 9];
        a3[1] = C64::new(1.0, 0.0);
        a3[3] = C64::new(1.0, 0.0);
        a3[5] = C64::new(1.0, 0.0);
        a3[7] = C64::new(1.0, 0.0);
        let mut b3 = vec![C64::new(0.0, 0.0); 9];
        b3[0] = C64::new(1.0, 0.0);
        b3[8] = C64::new(-1.0, 0.0);
        let i3 = ident(3);
        let i2 = ident(2);
        let gens6 = vec![
            kron(&a3, 3, &i2, 2),
            kron(&b3, 3, &i2, 2),
            kron(&i3, 3, &pauli('X'), 2),
            kron(&i3, 3, &pauli('Z'), 2),
        ];
        let det6 = recover_factorization(&build_net(&gens6, n6, tau, 0.0, &mut rng), &gens6, n6);
        let want222 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 2, 2],
        };
        let want23 = FactorizationReading::ExactUpToLocalUnitaryAndPermutation {
            local_dims: vec![2, 3],
        };
        let ok = det_site.reading == want222 && det_mode.reading == want222 && det6.reading == want23;
        check(
            "[F0] exact oracle lane — site/DFT 共役 net → Exact [2,2,2]・qutrit×qubit → Exact [2,3]",
            ok,
            format!(
                "site = {} {:?} / mode = {} / C⁶ = {} {:?}",
                det_site.reading.as_str(),
                if let FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims } =
                    &det_site.reading
                {
                    local_dims.clone()
                } else {
                    vec![]
                },
                det_mode.reading.as_str(),
                det6.reading.as_str(),
                if let FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims } =
                    &det6.reading
                {
                    local_dims.clone()
                } else {
                    vec![]
                }
            ),
        );
    }

    // ---- [F1] 不可能セルの正棄却 ----
    {
        let numbers = vec![op3("ZII"), op3("IZI"), op3("IIZ")];
        let det_num =
            recover_factorization(&build_net(&numbers, 8, tau, 0.0, &mut rng), &numbers, 8);
        let partial = vec![op3("XII"), op3("ZII")];
        let det_part =
            recover_factorization(&build_net(&partial, 8, tau, 0.0, &mut rng), &partial, 8);
        let want =
            FactorizationReading::Abstain(FactorizationAbstainReason::InsufficientOperationalGenerators);
        let ok = det_num.reading == want && det_part.reading == want;
        check(
            "[F1] 不可能セルの正棄却 — number operator のみ / 部分 net → Abstain(InsufficientOperationalGenerators)",
            ok,
            format!(
                "{{Z₁,Z₂,Z₃}} → {} / {{X₁,Z₁}} → {} (未 address 自由度で tensor を主張しない)",
                det_num.reading.as_str(),
                det_part.reading.as_str()
            ),
        );
    }

    // ---- [F2] 中心非自明 → superselection ----
    {
        // (a) 電荷つき 2 sector: C⁶ = C² ⊕ (C²⊗C²)
        let n6 = 6;
        let embed = |top: &[C64], bot: &[C64]| -> Vec<C64> {
            let mut m = vec![C64::new(0.0, 0.0); n6 * n6];
            for i in 0..2 {
                for j in 0..2 {
                    m[i * n6 + j] = top[i * 2 + j];
                }
            }
            for i in 0..4 {
                for j in 0..4 {
                    m[(2 + i) * n6 + (2 + j)] = bot[i * 4 + j];
                }
            }
            m
        };
        let x2 = pauli('X');
        let z2 = pauli('Z');
        let i2 = ident(2);
        let xi = kron(&x2, 2, &i2, 2);
        let zi = kron(&z2, 2, &i2, 2);
        let zero2 = vec![C64::new(0.0, 0.0); 4];
        let zero4 = vec![C64::new(0.0, 0.0); 16];
        let g1 = embed(&x2, &zero4);
        let g2 = embed(&z2, &zero4);
        let g3 = embed(&zero2, &xi);
        let g4 = embed(&zero2, &zi);
        // 電荷 (sector を分ける古典観測量)
        let mut charge = vec![C64::new(0.0, 0.0); n6 * n6];
        for i in 0..2 {
            charge[i * n6 + i] = C64::new(1.0, 0.0);
        }
        for i in 2..6 {
            charge[i * n6 + i] = C64::new(-1.0, 0.0);
        }
        let gens_a = vec![g1, g2, g3, g4, charge];
        let det_a = recover_factorization(&build_net(&gens_a, n6, tau, 0.0, &mut rng), &gens_a, n6);
        // (b) 部分 address net {X₁, Z₁, Z₂}: qubit₁ 完全制御・qubit₂ は古典ラベル
        //     Z₂ のみ (測定軸だけ届く)・qubit₃ 未 address → Z₂ = ±1 が超選択
        //     ラベル・qubit₃ が多重度 2
        let gens_b = vec![op3("XII"), op3("ZII"), op3("IZI")];
        let det_b = recover_factorization(&build_net(&gens_b, 8, tau, 0.0, &mut rng), &gens_b, 8);
        // 対照: {X₁, Z₁Z₂} は非局所符号化の 1 qubit (dim 4 の因子・中心自明) —
        //     正しい読みは Insufficient (超選択ではない)
        let gens_c = vec![op3("XII"), op3("ZZI")];
        let det_c = recover_factorization(&build_net(&gens_c, 8, tau, 0.0, &mut rng), &gens_c, 8);
        let want_a = FactorizationReading::SuperselectionSectors {
            sectors: vec![(2, 1), (2, 2)],
        };
        let want_b = FactorizationReading::SuperselectionSectors {
            sectors: vec![(2, 2), (2, 2)],
        };
        let want_c =
            FactorizationReading::Abstain(FactorizationAbstainReason::InsufficientOperationalGenerators);
        let ok = det_a.reading == want_a && det_b.reading == want_b && det_c.reading == want_c;
        check(
            "[F2] 中心非自明 → SuperselectionSectors (tensor 強制の禁止) — 電荷 net [(2,1),(2,2)]・部分 address [(2,2),(2,2)]",
            ok,
            format!(
                "電荷つき C²⊕C²⊗C² → {} {:?} / {{X₁,Z₁,Z₂}} → {} {:?} (測定だけの軸 = 超選択ラベル・未 address = 多重度) / 対照 {{X₁,Z₁Z₂}} → {} (非局所符号化 1 qubit — 中心自明 dim 4 の因子で Insufficient が正)",
                det_a.reading.as_str(),
                if let FactorizationReading::SuperselectionSectors { sectors } = &det_a.reading {
                    sectors.clone()
                } else {
                    vec![]
                },
                det_b.reading.as_str(),
                if let FactorizationReading::SuperselectionSectors { sectors } = &det_b.reading {
                    sectors.clone()
                } else {
                    vec![]
                },
                det_c.reading.as_str()
            ),
        );
    }

    // ---- [F3] gauge orbit の裁定 ----
    {
        // 局所 unitary × 置換: u_i = exp(iθ n·σ) を各 qubit に・SWAP₁₃ で置換
        let rot = |theta: f64, nx: f64, ny: f64, nz: f64| -> Vec<C64> {
            let (c, s) = (theta.cos(), theta.sin());
            // cosθ I + i sinθ (n·σ)
            vec![
                C64::new(c, s * nz),
                C64::new(s * ny, s * nx),
                C64::new(-s * ny, s * nx),
                C64::new(c, -s * nz),
            ]
        };
        let u1 = rot(0.3, 1.0, 0.0, 0.0);
        let u2 = rot(0.7, 0.0, 1.0, 0.0);
        let u3 = rot(1.1, 0.6, 0.0, 0.8);
        let u12 = kron(&u1, 2, &u2, 2);
        let u = kron(&u12, 4, &u3, 2);
        // SWAP₁₃: ビット反転置換行列
        let mut p13 = vec![C64::new(0.0, 0.0); 64];
        for b in 0..8usize {
            let rb = ((b & 1) << 2) | (b & 2) | ((b >> 2) & 1);
            p13[rb * 8 + b] = C64::new(1.0, 0.0);
        }
        let w = cmul(&p13, &u, 8);
        let rotated: Vec<Vec<C64>> = site.iter().map(|g| conj_by(&w, g, 8)).collect();
        let det_rot =
            recover_factorization(&build_net(&rotated, 8, tau, 0.0, &mut rng), &rotated, 8);
        let (same_rot, ov_rot) =
            same_gauge_orbit(&det_site.component_subalgebras, &det_rot.component_subalgebras);
        let (same_dft, ov_dft) =
            same_gauge_orbit(&det_site.component_subalgebras, &det_mode.component_subalgebras);
        // 裁定: matching あり → 一つの orbit (Exact のまま) / なし → EquivalenceClassOnly
        let verdict_dft = if same_dft {
            FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims: vec![2, 2, 2] }
        } else {
            FactorizationReading::EquivalenceClassOnly {
                class_desc: "site 因子分解と DFT 因子分解 — 完全だが非互換な 2 つの marking (無制約 tie-breaker で 1 つを選ばない)".into(),
            }
        };
        let ok = det_rot.reading.as_str() == "exact_up_to_local_unitary_and_permutation"
            && same_rot
            && ov_rot >= 1.0 - 1e-9
            && !same_dft
            && ov_dft < 0.9
            && verdict_dft.as_str() == "equivalence_class_only";
        check(
            "[F3] gauge orbit — 局所 unitary × SWAP は同一 orbit (overlap 1)・site vs DFT は EquivalenceClassOnly",
            ok,
            format!(
                "回転+置換 net: {} (matching overlap = {:.12}) / site×DFT: matching なし (best min-overlap = {:.6}) → {}",
                det_rot.reading.as_str(),
                ov_rot,
                ov_dft,
                verdict_dft.as_str()
            ),
        );
    }

    // ---- [F4] noise abstention ----
    {
        let det_clean =
            recover_factorization(&build_net(&site, 8, tau, 1e-6, &mut rng), &site, 8);
        let det_noisy =
            recover_factorization(&build_net(&site, 8, tau, 5e-4, &mut rng), &site, 8);
        let ok = det_clean.reading.as_str() == "exact_up_to_local_unitary_and_permutation"
            && det_noisy.reading
                == FactorizationReading::Abstain(
                    FactorizationAbstainReason::CommutatorMarginStraddled,
                );
        check(
            "[F4] noise abstention — σ=1e-6 は Exact・σ=5e-4 は可換対の区間が閾値を跨ぎ Abstain",
            ok,
            format!(
                "σ=1e-6 → {} / σ=5e-4 → {} (辺の強制なし — 假定 5 の機械化)",
                det_clean.reading.as_str(),
                det_noisy.reading.as_str()
            ),
        );
    }

    // ---- [F5] graded scope — parity-even net はパリティ超選択を発見する ----
    {
        // Majorana 双線形の path: h_k = iγ_k γ_{k+1} (k = 1..5) — 全て parity-even
        let gam = [op3("XII"), op3("YII"), op3("ZXI"), op3("ZYI"), op3("ZZX"), op3("ZZY")];
        let bil = |a: &[C64], b: &[C64]| -> Vec<C64> {
            cmul(a, b, 8).iter().map(|c| C64::new(-c.im, c.re)).collect()
        };
        let gens_f: Vec<Vec<C64>> = (0..5).map(|k| bil(&gam[k], &gam[k + 1])).collect();
        let det_f = recover_factorization(&build_net(&gens_f, 8, tau, 0.0, &mut rng), &gens_f, 8);
        let want = FactorizationReading::SuperselectionSectors {
            sectors: vec![(4, 1), (4, 1)],
        };
        // 偶部分代数の次元と中心の直接照合
        let cl = algebra_closure(&gens_f, 8);
        let cen = center_basis(&cl, &gens_f, 8);
        // 中心 2 次元の非自明成分がパリティ Γ = −γ₁γ₂γ₃γ₄γ₅γ₆·i³ (= ZZZ) に一致するか
        let gamma = op3("ZZZ");
        let idn = ident(8);
        let mut span_ok = false;
        if cen.len() == 2 {
            // span{I, ZZZ} との一致 (射影残差)
            let mut basis_ref: Vec<Vec<C64>> = Vec::new();
            push_ortho(&mut basis_ref, &idn, 1e-9);
            push_ortho(&mut basis_ref, &gamma, 1e-9);
            let mut resid = 0.0f64;
            for c in &cen {
                let mut v = c.clone();
                for b in &basis_ref {
                    let x = hs_inner(b, &v);
                    for (vi, bi) in v.iter_mut().zip(b.iter()) {
                        *vi = *vi - x * *bi;
                    }
                }
                resid = resid.max(hs_norm(&v));
            }
            span_ok = resid < 1e-8;
        }
        let ok = det_f.reading == want && cl.len() == 32 && cen.len() == 2 && span_ok;
        check(
            "[F5] parity-even フェルミオン net — 復元は tensor でなく [(4,1),(4,1)] = パリティ超選択の機械発見",
            ok,
            format!(
                "偶代数 dim = {} (= 32)・中心 dim = {} (span{{I, Z₁Z₂Z₃ = Γ}} 一致 = {})・裁定 = {} {:?}",
                cl.len(),
                cen.len(),
                span_ok,
                det_f.reading.as_str(),
                if let FactorizationReading::SuperselectionSectors { sectors } = &det_f.reading {
                    sectors.clone()
                } else {
                    vec![]
                }
            ),
        );
    }

    // ---- [F6] 文書アンカー ----
    {
        let mut bad = Vec::new();
        let doc = rd("docs/uft-v32.3.md").unwrap_or_default();
        for needle in [
            "marked operational recovery",
            "SuperselectionSectors",
            "EquivalenceClassOnly",
            "CommutatorMarginStraddled",
            "パリティ超選択",
            "tensor product を強制しない",
        ] {
            if !doc.contains(needle) {
                bad.push(format!("uft-v32.3.md: 「{}」が無い", needle));
            }
        }
        check(
            "[F6] 文書アンカー — 定理 B・三値裁定・凍結決定手順",
            bad.is_empty(),
            if bad.is_empty() {
                "復元は定理・裁定・文書の三点で凍結された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "marking からの因子分解復元が三値裁定つきで機械化された — 読めるときは gauge orbit まで・中心非自明は sector・不能は棄却"
        } else {
            "**復元器の破れ** — 決定手順と文書の整合を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
