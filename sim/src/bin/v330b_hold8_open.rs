//! v33.0-B HOLD-8 の開封 — SECRET 開示・holdout 初生成・本採点 (調整なし)
//!
//! SECRET を開示して holdout を初生成・本採点する。**調整なし**。
//!   [H0] SECRET がコミットメント fb05a4a0… と一致 + FROZEN-HOLD8 区間が
//!        v330a と逐語一致 (SHA-256) + train seed 不変
//!   [H1] holdout 20 セル (SECRET 導出シード) の本採点 — セル表と 3 計量
//!        (selective risk / impossibility recall / answerable recall) を記録
//!   [H2] 採点の自己整合 — 計量がセル表と一致し、強制回答が risk に数えられている
//!
//! 実行: cargo run --release --bin v330b_hold8_open

use uft_sim::operational_net::*;
use uft_sim::{sha256_hex, Rng, C64};

// ================================================================================
// FROZEN-HOLD8-BEGIN  (この区間は v330b と逐語一致 — [H0] が SHA-256 で照合する)
// ================================================================================

pub const HOLD8_COMMITMENT: &str =
    "fb05a4a07dd7feef78e0bcecdaeeff933a656d36a38ceb26c71047002621c582";
pub const HOLD8_TRAIN_SEED: u64 = 33001;

// ---- 凍結バー (開封後に変更しない) ----
pub const TAU_COMM: f64 = 1e-3; // 可換子閾値 (Pauli 尺度 — 真の辺は ≥ 2√2)
pub const NOISE_Z: f64 = 6.0; // ノイズ区間の z (証明書)
pub const GAP_GUARD: f64 = 1e-3; // gap 則スケールガード (v31.6 最終規則)
pub const BAR_W_REL: f64 = 1e-8; // 重み読みの相対バー (クリーン)
pub const BAR_HYPER_REL: f64 = 1e-6; // hyperedge 検出の相対バー (条件差/max)
pub const BAR_COHERENT: f64 = 0.05; // coherent 分離の最小差

// ---- 凍結素子 (v32.2–v32.5 の逐語コピー) ----

pub fn pauli(which: char) -> Vec<C64> {
    let (o, l) = (C64::new(0.0, 0.0), C64::new(1.0, 0.0));
    match which {
        'I' => vec![l, o, o, l],
        'X' => vec![o, l, l, o],
        'Y' => vec![o, C64::new(0.0, -1.0), C64::new(0.0, 1.0), o],
        'Z' => vec![l, o, o, C64::new(-1.0, 0.0)],
        _ => panic!("未知の Pauli"),
    }
}

pub fn kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
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

pub fn op3(s: &str) -> Vec<C64> {
    let cs: Vec<char> = s.chars().collect();
    let a = kron(&pauli(cs[0]), 2, &pauli(cs[1]), 2);
    kron(&a, 4, &pauli(cs[2]), 2)
}

pub fn ident(n: usize) -> Vec<C64> {
    let mut m = vec![C64::new(0.0, 0.0); n * n];
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}

pub fn add_scaled(a: &mut [C64], b: &[C64], s: f64) {
    for (x, y) in a.iter_mut().zip(b.iter()) {
        *x = *x + y.scale(s);
    }
}

pub fn dft8() -> Vec<C64> {
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

pub fn conj_by(v: &[C64], a: &[C64], n: usize) -> Vec<C64> {
    cmul(&cmul(v, a, n), &cdag(v, n), n)
}

pub fn rot2(theta: f64, nx: f64, ny: f64, nz: f64) -> Vec<C64> {
    let (c, s) = (theta.cos(), theta.sin());
    vec![
        C64::new(c, s * nz),
        C64::new(s * ny, s * nx),
        C64::new(-s * ny, s * nx),
        C64::new(c, -s * nz),
    ]
}

pub fn herm_evals(m: &[C64], n: usize) -> Vec<f64> {
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

pub fn jw_annihilators() -> Vec<Vec<C64>> {
    let sm = {
        let mut m = vec![C64::new(0.0, 0.0); 4];
        m[1] = C64::new(1.0, 0.0);
        m
    };
    let z = pauli('Z');
    let i2 = ident(2);
    let a1 = kron(&kron(&sm, 2, &i2, 2), 4, &i2, 2);
    let a2 = kron(&kron(&z, 2, &sm, 2), 4, &i2, 2);
    let a3 = kron(&kron(&z, 2, &z, 2), 4, &sm, 2);
    vec![a1, a2, a3]
}

pub fn number_op(k: usize) -> Vec<C64> {
    let ann = jw_annihilators();
    cmul(&cdag(&ann[k], 8), &ann[k], 8)
}

pub fn hop_op(i: usize, j: usize) -> Vec<C64> {
    let ann = jw_annihilators();
    let a = cmul(&cdag(&ann[i], 8), &ann[j], 8);
    let b = cmul(&cdag(&ann[j], 8), &ann[i], 8);
    a.iter().zip(b.iter()).map(|(x, y)| *x + *y).collect()
}

/// pair hopping Δ(c†_i c†_j + h.c.)
pub fn pair_op(i: usize, j: usize) -> Vec<C64> {
    let ann = jw_annihilators();
    let a = cmul(&cdag(&ann[i], 8), &cdag(&ann[j], 8), 8);
    let b = cmul(&ann[j], &ann[i], 8);
    a.iter().zip(b.iter()).map(|(x, y)| *x + *y).collect()
}

/// hopping H(θ) = Σ t_b (e^{iθ} a†_i a_j + h.c.)
pub fn ring_h(theta: f64, bonds: &[(usize, usize, f64)]) -> Vec<C64> {
    let ann = jw_annihilators();
    let n = 8;
    let mut h = vec![C64::new(0.0, 0.0); n * n];
    for &(i, j, t) in bonds {
        let hopm = cmul(&cdag(&ann[i], n), &ann[j], n);
        let phase = C64::expi(theta).scale(t);
        for (kx, hv) in hopm.iter().enumerate() {
            h[kx] = h[kx] + phase * *hv;
        }
        let hop2 = cmul(&cdag(&ann[j], n), &ann[i], n);
        let phase2 = C64::expi(-theta).scale(t);
        for (kx, hv) in hop2.iter().enumerate() {
            h[kx] = h[kx] + phase2 * *hv;
        }
    }
    h
}

pub fn r1_exact(h: &[C64], b: &[C64], a: &[C64], n: usize) -> f64 {
    let c = commutator(h, a, n);
    let mut s = C64::new(0.0, 0.0);
    for i in 0..n {
        for k in 0..n {
            s = s + b[i * n + k] * c[k * n + i];
        }
    }
    s.im
}

pub fn r2_exact(h: &[C64], b: &[C64], a: &[C64], n: usize) -> f64 {
    let hb = commutator(h, b, n);
    let ha = commutator(h, a, n);
    let mut s = C64::new(0.0, 0.0);
    for i in 0..n {
        for k in 0..n {
            s = s + hb[i * n + k] * ha[k * n + i];
        }
    }
    s.re
}

/// 積状態接ベクトル (n_i − 1/2)/4, 条件付き (n_i − 1/2)(I/2)P_v(k) — v32.5 と同一
pub fn tangent_uncond(i: usize) -> Vec<C64> {
    let n = 8;
    let mut a = vec![C64::new(0.0, 0.0); n * n];
    for idx in 0..n {
        let bi = (idx >> (2 - i)) & 1;
        a[idx * n + idx] = C64::new((bi as f64 - 0.5) * 0.25, 0.0);
    }
    a
}

pub fn tangent_cond(i: usize, k: usize, v: usize) -> Vec<C64> {
    let n = 8;
    let mut a = vec![C64::new(0.0, 0.0); n * n];
    for idx in 0..n {
        let bi = (idx >> (2 - i)) & 1;
        let bk = (idx >> (2 - k)) & 1;
        if bk != v {
            continue;
        }
        a[idx * n + idx] = C64::new((bi as f64 - 0.5) * 0.5, 0.0);
    }
    a
}

/// gap 支持 (v31.6 最終規則 — 逐語)
pub fn support_from_weights(w: &[f64], n: usize) -> Vec<Vec<bool>> {
    let mut vals: Vec<f64> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            vals.push(w[i * n + j].abs().max(1e-300));
        }
    }
    let mut sorted = vals.clone();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
    let guard = sorted[0] * GAP_GUARD;
    let mut cut: Option<usize> = None;
    let mut best_gap = 0.0;
    for k in 0..sorted.len() - 1 {
        if sorted[k + 1] < guard {
            break;
        }
        let gap = (sorted[k] / sorted[k + 1]).ln();
        if gap > best_gap {
            best_gap = gap;
            cut = Some(k);
        }
    }
    let thr = match cut {
        Some(k) if best_gap >= 3.0f64.ln() => (sorted[k] * sorted[k + 1]).sqrt(),
        _ => guard,
    };
    let mut adj = vec![vec![false; n]; n];
    for i in 0..n {
        for j in 0..n {
            if i != j && w[i * n + j].abs() > thr {
                adj[i][j] = true;
            }
        }
    }
    adj
}

// ---- 因子分解 recovery (v32.3 の逐語コピー) ----

pub fn center_basis(basis: &[Vec<C64>], gens: &[Vec<C64>], n: usize) -> Vec<Vec<C64>> {
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
    let (evals, vecs) = uft_sim::jacobi_eigh(&gram, dim_r);
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
                let d = *a - *b;
                C64::new(d.im * 0.5, -d.re * 0.5)
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

pub fn central_projectors(center: &[Vec<C64>], n: usize) -> Option<Vec<Vec<C64>>> {
    let mut t = vec![C64::new(0.0, 0.0); n * n];
    for (k, h) in center.iter().enumerate() {
        let w = ((k + 2) as f64).sqrt();
        for (ti, hi) in t.iter_mut().zip(h.iter()) {
            *ti = *ti + hi.scale(w);
        }
    }
    let evals = herm_evals(&t, n);
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

pub struct RecoveryDetail {
    pub reading: FactorizationReading,
    pub component_subalgebras: Vec<Vec<Vec<C64>>>,
}

pub fn recover_factorization<G: CommutationGrading>(
    net: &OperationalNet<G>,
    gens: &[Vec<C64>],
    n: usize,
) -> RecoveryDetail {
    let abstain = |r: FactorizationAbstainReason| RecoveryDetail {
        reading: FactorizationReading::Abstain(r),
        component_subalgebras: Vec::new(),
    };
    let comps = match net.noncommutation_components() {
        Ok(c) => c,
        Err(r) => return abstain(r),
    };
    let joint = algebra_closure(gens, n);
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
    let center = center_basis(&joint, gens, n);
    if center.is_empty() {
        return abstain(FactorizationAbstainReason::ComponentNotFactor);
    }
    if center.len() == 1 {
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
            let comp_center = center_basis(&cl, &sub, n);
            if comp_center.len() != 1 {
                return abstain(FactorizationAbstainReason::ComponentNotFactor);
            }
            dims.push(d);
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

pub fn same_gauge_orbit(a: &[Vec<Vec<C64>>], b: &[Vec<Vec<C64>>]) -> (bool, f64) {
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

/// exact ノルム (± comm ノイズ) の証明書つき Ordinary net (v32.3 と同一規約)
pub fn build_net(
    gens: &[Vec<C64>],
    n: usize,
    sigma: f64,
    rng: &mut Rng,
) -> OperationalNet<OrdinaryCommutation> {
    let mut net: OperationalNet<OrdinaryCommutation> = OperationalNet::new(n, TAU_COMM);
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
            provenance: "hold8_control",
        };
        ids.push(net.add_primitive(p).unwrap());
    }
    for (a, ga) in gens.iter().enumerate() {
        for (b, gb) in gens.iter().enumerate().skip(a + 1) {
            let nu = hs_norm(&commutator(ga, gb, n));
            let (lo, hi) = if sigma > 0.0 {
                let nu_hat = (nu + sigma * rng.gauss()).abs();
                ((nu_hat - NOISE_Z * sigma).max(0.0), nu_hat + NOISE_Z * sigma)
            } else {
                ((nu - 1e-12).max(0.0), nu + 1e-12)
            };
            net.set_commutator(ids[a], ids[b], CertifiedCommutator::new(lo, hi).unwrap());
        }
    }
    net
}

// ---- セル定義 (task = readout の入力 / expect = 採点器の真値) ----

pub enum Task {
    /// 因子分解: net + 生成子 (可換子ノイズ σ)
    Recover { gens: Vec<Vec<C64>>, sigma: f64 },
    /// odd 入力 (ordinary は構成時拒否) + graded 化した偶双線形
    RecoverGraded { odd: Vec<Vec<C64>>, even: Vec<Vec<C64>> },
    /// 非互換 2 net の照合
    Reconcile { gens_a: Vec<Vec<C64>>, gens_b: Vec<Vec<C64>> },
    /// 相互作用: 密度 (+条件付き) readout・coherent 可否・K ノイズ σ
    Interaction { h: Vec<C64>, coherent: bool, sigma_k: f64 },
    /// 対セル (H_a vs H_b): 密度 lane は同値裁定・coherent lane は分離
    Pair { h_a: Vec<C64>, h_b: Vec<C64>, coherent: bool },
}

pub enum Expect {
    FactorExact { dims: Vec<usize>, orbit_gens: Vec<Vec<C64>> },
    FactorSectors { sectors: Vec<(usize, usize)> },
    FactorAbstain(FactorizationAbstainReason),
    GradedSectors { sectors: Vec<(usize, usize)> },
    EquivalenceOnly,
    SameOrbit,
    Support {
        edges: Vec<(usize, usize)>,
        hyper: Vec<(usize, usize, usize)>,
        wexp: Vec<((usize, usize), f64)>,
        cond: Vec<((usize, usize, usize), f64, f64)>, // (i,j,k): K(0), K(1)
    },
    SupportAbstain,
    PairDensityEquivCoherentSplit { coherent: bool },
}

pub struct Cell {
    pub name: String,
    pub answerable: bool,
    pub task: Task,
    pub expect: Expect,
}

/// site 生成子 {X_i, Z_i}
pub fn site_gens() -> Vec<Vec<C64>> {
    (0..3)
        .flat_map(|i| {
            let mut s = ['I', 'I', 'I'];
            s[i] = 'X';
            let x = op3(&s.iter().collect::<String>());
            s[i] = 'Z';
            let z = op3(&s.iter().collect::<String>());
            vec![x, z]
        })
        .collect()
}

/// 隠しパラメータの局所 unitary × 置換 W = P_π (u₁⊗u₂⊗u₃)
pub fn hidden_w(rng: &mut Rng) -> Vec<C64> {
    let u1 = rot2(rng.f64() * 1.2 + 0.1, 1.0, 0.0, 0.0);
    let u2 = rot2(rng.f64() * 1.2 + 0.1, 0.0, 1.0, 0.0);
    let u3 = rot2(
        rng.f64() * 1.2 + 0.1,
        0.6,
        0.0,
        0.8,
    );
    let u = kron(&kron(&u1, 2, &u2, 2), 4, &u3, 2);
    // 置換: 3 ラベルの巡回シフト 0/1/2 回
    let shift = rng.range(3);
    let mut p = vec![C64::new(0.0, 0.0); 64];
    for b in 0..8usize {
        let bits = [(b >> 2) & 1, (b >> 1) & 1, b & 1];
        let mut nb = [0usize; 3];
        for q in 0..3 {
            nb[(q + shift) % 3] = bits[q];
        }
        let rb = (nb[0] << 2) | (nb[1] << 1) | nb[2];
        p[rb * 8 + b] = C64::new(1.0, 0.0);
    }
    cmul(&p, &u, 8)
}

/// セル一式の生成 (seed が train / holdout を分ける — 構造は凍結・数値は隠し)
pub fn make_cells(seed: u64) -> Vec<Cell> {
    let mut rng = Rng::new(seed);
    let mut cells = Vec::new();
    let site = site_gens();
    let conj_all =
        |w: &Vec<C64>, gens: &Vec<Vec<C64>>| -> Vec<Vec<C64>> {
            gens.iter().map(|g| conj_by(w, g, 8)).collect()
        };

    // F1: site net (隠し局所 U × 置換) → Exact [2,2,2]
    {
        let w = hidden_w(&mut rng);
        let gens = conj_all(&w, &site);
        cells.push(Cell {
            name: "F1-site".into(),
            answerable: true,
            expect: Expect::FactorExact { dims: vec![2, 2, 2], orbit_gens: gens.clone() },
            task: Task::Recover { gens, sigma: 0.0 },
        });
    }
    // F2: mode (DFT) net (隠し局所 U) → Exact [2,2,2] (site と別 orbit)
    {
        let w0 = hidden_w(&mut rng);
        let w = cmul(&w0, &dft8(), 8);
        let gens = conj_all(&w, &site);
        cells.push(Cell {
            name: "F2-mode".into(),
            answerable: true,
            expect: Expect::FactorExact { dims: vec![2, 2, 2], orbit_gens: gens.clone() },
            task: Task::Recover { gens, sigma: 0.0 },
        });
    }
    // F3: number operator のみ → Abstain(Insufficient)
    {
        let w = hidden_w(&mut rng);
        let nums: Vec<Vec<C64>> = vec![op3("ZII"), op3("IZI"), op3("IIZ")];
        let gens = conj_all(&w, &nums);
        cells.push(Cell {
            name: "F3-numberonly".into(),
            answerable: false,
            expect: Expect::FactorAbstain(
                FactorizationAbstainReason::InsufficientOperationalGenerators,
            ),
            task: Task::Recover { gens, sigma: 0.0 },
        });
    }
    // F4: 完全だが非互換な 2 net → EquivalenceClassOnly
    {
        let wa = hidden_w(&mut rng);
        let wb = cmul(&hidden_w(&mut rng), &dft8(), 8);
        cells.push(Cell {
            name: "F4-incompatible".into(),
            answerable: false,
            expect: Expect::EquivalenceOnly,
            task: Task::Reconcile {
                gens_a: conj_all(&wa, &site),
                gens_b: conj_all(&wb, &site),
            },
        });
    }
    // F5: 中心非自明 {X_a, Z_a, Z_b} → [(2,2),(2,2)]
    {
        let w = hidden_w(&mut rng);
        let raw = vec![op3("XII"), op3("ZII"), op3("IZI")];
        let gens = conj_all(&w, &raw);
        cells.push(Cell {
            name: "F5-superselect".into(),
            answerable: true,
            expect: Expect::FactorSectors { sectors: vec![(2, 2), (2, 2)] },
            task: Task::Recover { gens, sigma: 0.0 },
        });
    }
    // F6: odd 入力は ordinary 拒否・偶双線形 (隠し係数) は パリティ sector
    {
        let gam = [op3("XII"), op3("YII"), op3("ZXI"), op3("ZYI"), op3("ZZX"), op3("ZZY")];
        let odd: Vec<Vec<C64>> = vec![gam[0].clone(), gam[2].clone(), gam[4].clone()];
        let bil = |a: &[C64], b: &[C64]| -> Vec<C64> {
            cmul(a, b, 8).iter().map(|c| C64::new(-c.im, c.re)).collect()
        };
        let even: Vec<Vec<C64>> = (0..5)
            .map(|k| {
                let lam = 0.5 + rng.f64();
                bil(&gam[k], &gam[k + 1]).iter().map(|c| c.scale(lam)).collect()
            })
            .collect();
        cells.push(Cell {
            name: "F6-graded".into(),
            answerable: true,
            expect: Expect::GradedSectors { sectors: vec![(4, 1), (4, 1)] },
            task: Task::RecoverGraded { odd, even },
        });
    }
    // I1: quadratic (隠し重み・ring/chain の隠し選択)
    {
        let t1 = 0.6 + 0.8 * rng.f64();
        let t2 = 0.6 + 0.8 * rng.f64();
        let ring = rng.f64() > 0.5;
        let mut bonds = vec![(0usize, 1usize, t1), (1, 2, t2)];
        let mut edges = vec![(0usize, 1usize), (1, 2)];
        let mut wexp = vec![((0usize, 1usize), t1 * t1), ((1, 2), t2 * t2)];
        if ring {
            let t3 = 0.6 + 0.8 * rng.f64();
            bonds.push((0, 2, t3));
            edges.push((0, 2));
            wexp.push(((0, 2), t3 * t3));
        }
        edges.sort_unstable();
        wexp.sort_by_key(|(e, _)| *e);
        cells.push(Cell {
            name: format!("I1-quadratic{}", if ring { "-ring" } else { "-chain" }),
            answerable: true,
            expect: Expect::Support { edges, hyper: vec![], wexp, cond: vec![] },
            task: Task::Interaction { h: ring_h(0.0, &bonds), coherent: false, sigma_k: 0.0 },
        });
    }
    // I2: t-V (密度対角) — V は密度曲率に不可視 (ŵ = t²)
    {
        let t1 = 0.6 + 0.8 * rng.f64();
        let t2 = 0.6 + 0.8 * rng.f64();
        let v = 0.4 + 0.6 * rng.f64();
        let mut h = ring_h(0.0, &[(0, 1, t1), (1, 2, t2)]);
        let nn = cmul(&number_op(0), &number_op(1), 8);
        add_scaled(&mut h, &nn, v);
        cells.push(Cell {
            name: "I2-tv".into(),
            answerable: true,
            expect: Expect::Support {
                edges: vec![(0, 1), (1, 2)],
                hyper: vec![],
                wexp: vec![((0, 1), t1 * t1), ((1, 2), t2 * t2)],
                cond: vec![],
            },
            task: Task::Interaction { h, coherent: false, sigma_k: 0.0 },
        });
    }
    // I3: 相関 hopping V n₂ h₀₁ — 条件付きが |t + vV|² を分離
    {
        let t1 = 0.6 + 0.8 * rng.f64();
        let t2 = 0.6 + 0.8 * rng.f64();
        let v = 0.4 + 0.6 * rng.f64();
        let mut h = ring_h(0.0, &[(0, 1, t1), (1, 2, t2)]);
        let corr = cmul(&number_op(2), &hop_op(0, 1), 8);
        add_scaled(&mut h, &corr, v);
        cells.push(Cell {
            name: "I3-corrhop".into(),
            answerable: true,
            expect: Expect::Support {
                edges: vec![(0, 1), (1, 2)],
                hyper: vec![(0, 1, 2)],
                wexp: vec![],
                cond: vec![((0, 1, 2), t1 * t1, (t1 + v) * (t1 + v))],
            },
            task: Task::Interaction { h, coherent: false, sigma_k: 0.0 },
        });
    }
    // I4: pair hopping Δ(c†₀c†₁ + h.c.) + t h₁₂ — 支持は {0,1},{1,2}・hyper なし
    {
        let dl = 0.6 + 0.8 * rng.f64();
        let t2 = 0.6 + 0.8 * rng.f64();
        let mut h = ring_h(0.0, &[(1, 2, t2)]);
        add_scaled(&mut h, &pair_op(0, 1), dl);
        cells.push(Cell {
            name: "I4-pairhop".into(),
            answerable: true,
            expect: Expect::Support {
                edges: vec![(0, 1), (1, 2)],
                hyper: vec![],
                wexp: vec![],
                cond: vec![],
            },
            task: Task::Interaction { h, coherent: false, sigma_k: 0.0 },
        });
    }
    // I5a/I5b: 三体項 g Z₂h₀₁ の有無 (hyperedge 検出器の対)
    {
        let t1 = 0.6 + 0.8 * rng.f64();
        let g = 0.3 + 0.4 * rng.f64();
        for (name, gg) in [("I5a-threebody", g), ("I5b-threebody-null", 0.0)] {
            let mut h = ring_h(0.0, &[(0, 1, t1), (1, 2, 0.8)]);
            let z2h: Vec<C64> = {
                let z2 = op3("IIZ");
                cmul(&z2, &hop_op(0, 1), 8)
            };
            add_scaled(&mut h, &z2h, gg);
            let hyper = if gg > 0.0 { vec![(0, 1, 2)] } else { vec![] };
            let cond = vec![(
                (0usize, 1usize, 2usize),
                (t1 + gg) * (t1 + gg), // v = 0: Z₂ = +1
                (t1 - gg) * (t1 - gg), // v = 1: Z₂ = −1
            )];
            cells.push(Cell {
                name: name.into(),
                answerable: true,
                expect: Expect::Support {
                    edges: vec![(0, 1), (1, 2)],
                    hyper,
                    wexp: vec![],
                    cond,
                },
                task: Task::Interaction { h, coherent: false, sigma_k: 0.0 },
            });
        }
    }
    // I6: H ↔ −H 対 — 密度は同値裁定・coherent は分離
    {
        let t1 = 0.6 + 0.8 * rng.f64();
        let h = ring_h(0.0, &[(0, 1, t1), (1, 2, 0.9)]);
        let hneg: Vec<C64> = h.iter().map(|c| c.scale(-1.0)).collect();
        for (name, coh, ans) in [("I6d-signpair-density", false, false), ("I6c-signpair-coherent", true, true)]
        {
            cells.push(Cell {
                name: name.into(),
                answerable: ans,
                expect: Expect::PairDensityEquivCoherentSplit { coherent: coh },
                task: Task::Pair { h_a: h.clone(), h_b: hneg.clone(), coherent: coh },
            });
        }
    }
    // I7: 磁束対 (θ, −θ) — 密度は同値・coherent は分離
    {
        let th = 0.3 + 0.6 * rng.f64();
        let bonds = [(0usize, 1usize, 1.0), (1, 2, 1.0), (0, 2, 1.0)];
        let ha = ring_h(th, &bonds);
        let hb = ring_h(-th, &bonds);
        for (name, coh, ans) in [("I7d-fluxpair-density", false, false), ("I7c-fluxpair-coherent", true, true)]
        {
            cells.push(Cell {
                name: name.into(),
                answerable: ans,
                expect: Expect::PairDensityEquivCoherentSplit { coherent: coh },
                task: Task::Pair { h_a: ha.clone(), h_b: hb.clone(), coherent: coh },
            });
        }
    }
    // M1: 変成対 (局所 U × 置換) — 同一 gauge orbit
    {
        let w1 = hidden_w(&mut rng);
        let w2 = hidden_w(&mut rng);
        let base = conj_all(&w1, &site);
        let re = conj_all(&w2, &base);
        cells.push(Cell {
            name: "M1-metamorphic".into(),
            answerable: true,
            expect: Expect::SameOrbit,
            task: Task::Reconcile { gens_a: base, gens_b: re },
        });
    }
    // M3: 演算子基底の可逆再結合 — recovery 不変
    {
        let mut gens = Vec::new();
        for i in 0..3 {
            let mut s = ['I', 'I', 'I'];
            s[i] = 'X';
            let x = op3(&s.iter().collect::<String>());
            s[i] = 'Z';
            let z = op3(&s.iter().collect::<String>());
            // 可逆な実結合 2 本 (det ≠ 0 まで draw)
            loop {
                let (a, b) = (rng.f64() * 2.0 - 1.0, rng.f64() * 2.0 - 1.0);
                let (c, d) = (rng.f64() * 2.0 - 1.0, rng.f64() * 2.0 - 1.0);
                if (a * d - b * c).abs() < 0.3 {
                    continue;
                }
                let mut g1 = vec![C64::new(0.0, 0.0); 64];
                add_scaled(&mut g1, &x, a);
                add_scaled(&mut g1, &z, b);
                let mut g2 = vec![C64::new(0.0, 0.0); 64];
                add_scaled(&mut g2, &x, c);
                add_scaled(&mut g2, &z, d);
                gens.push(g1);
                gens.push(g2);
                break;
            }
        }
        cells.push(Cell {
            name: "M3-recombination".into(),
            answerable: true,
            expect: Expect::FactorExact { dims: vec![2, 2, 2], orbit_gens: site.clone() },
            task: Task::Recover { gens, sigma: 0.0 },
        });
    }
    // M4: 可換子 margin 以下のノイズ → Abstain
    {
        let w = hidden_w(&mut rng);
        let gens = conj_all(&w, &site);
        cells.push(Cell {
            name: "M4-commnoise".into(),
            answerable: false,
            expect: Expect::FactorAbstain(FactorizationAbstainReason::CommutatorMarginStraddled),
            task: Task::Recover { gens, sigma: 5e-4 },
        });
    }
    // M5: support margin 以下の弱辺 + K ノイズ → Abstain(InsufficientObservation)
    {
        let t1 = 0.9 + 0.4 * rng.f64();
        let tw = 0.015 + 0.01 * rng.f64(); // 弱辺
        let h = ring_h(0.0, &[(0, 1, t1), (1, 2, tw)]);
        cells.push(Cell {
            name: "M5-weakedge".into(),
            answerable: false,
            expect: Expect::SupportAbstain,
            task: Task::Interaction { h, coherent: false, sigma_k: 5e-4 },
        });
    }
    cells
}

// ---- readout (真値に触れない — task のみを見る) ----

pub enum Verdict {
    Factor(RecoveryDetail),
    GradedFactor { ordinary_refused: bool, reading: FactorizationReading },
    ReconcileSame { overlap: f64 },
    ReconcileClass,
    SupportRead {
        edges: Vec<(usize, usize)>,
        hyper: Vec<(usize, usize, usize)>,
        w: Vec<Vec<f64>>,
        cond: Vec<((usize, usize, usize), f64, f64)>,
    },
    SupportAbstained,
    PairRead { density_identical: bool, coherent_split: Option<bool> },
}

pub fn density_kernel(h: &[C64], sigma: f64, rng: &mut Rng) -> Vec<Vec<f64>> {
    let mut k = vec![vec![0.0; 3]; 3];
    for j in 0..3 {
        for i in 0..3 {
            if i == j {
                continue;
            }
            let val = r2_exact(h, &number_op(j), &tangent_uncond(i), 8);
            k[j][i] = val + if sigma > 0.0 { sigma * rng.gauss() } else { 0.0 };
        }
    }
    k
}

pub fn readout(task: &Task, seed: u64) -> Verdict {
    let mut rng = Rng::new(seed);
    match task {
        Task::Recover { gens, sigma } => {
            let net = build_net(gens, 8, *sigma, &mut rng);
            Verdict::Factor(recover_factorization(&net, gens, 8))
        }
        Task::RecoverGraded { odd, even } => {
            // ordinary lane は odd を構成時拒否する (型ゲート)
            let mut net_o: OperationalNet<OrdinaryCommutation> = OperationalNet::new(8, TAU_COMM);
            let refused = odd.iter().all(|g| {
                net_o
                    .add_primitive(PrimitiveOperation {
                        kind: OpKind::Control(
                            ControlGenerator::certify(
                                g.iter().map(|c| c.re).collect(),
                                g.iter().map(|c| c.im).collect(),
                                8,
                            )
                            .unwrap(),
                        ),
                        parity: OperatorParity::Odd,
                        provenance: "hold8_odd",
                    })
                    .is_err()
            });
            let net = build_net(even, 8, 0.0, &mut rng);
            let det = recover_factorization(&net, even, 8);
            Verdict::GradedFactor { ordinary_refused: refused, reading: det.reading }
        }
        Task::Reconcile { gens_a, gens_b } => {
            let mut r1 = Rng::new(seed ^ 0x9e37);
            let mut r2 = Rng::new(seed ^ 0x79b9);
            let da = recover_factorization(&build_net(gens_a, 8, 0.0, &mut r1), gens_a, 8);
            let db = recover_factorization(&build_net(gens_b, 8, 0.0, &mut r2), gens_b, 8);
            let (same, ov) =
                same_gauge_orbit(&da.component_subalgebras, &db.component_subalgebras);
            if same {
                Verdict::ReconcileSame { overlap: ov }
            } else {
                Verdict::ReconcileClass
            }
        }
        Task::Interaction { h, coherent: _, sigma_k } => {
            let k = density_kernel(h, *sigma_k, &mut rng);
            // 対称化した重み行列で gap 支持 (順序対は平均 — クリーンでは同値)
            let mut w = vec![0.0; 9];
            for i in 0..3 {
                for j in 0..3 {
                    if i != j {
                        w[i * 3 + j] = 0.5 * (k[i][j] + k[j][i]);
                    }
                }
            }
            // SupportNoiseCertificate (v32.1 凍結): σ·z が窓ガードを跨ぐなら棄却
            if *sigma_k > 0.0 {
                let maxw = w.iter().cloned().fold(0.0f64, f64::max);
                if NOISE_Z * *sigma_k > GAP_GUARD * maxw {
                    return Verdict::SupportAbstained;
                }
            }
            let adj = support_from_weights(&w, 3);
            let mut edges = Vec::new();
            for i in 0..3 {
                for j in (i + 1)..3 {
                    if adj[i][j] || adj[j][i] {
                        edges.push((i, j));
                    }
                }
            }
            // hyperedge 検出: 条件差 |K^{k,1} − K^{k,0}| > bar·max
            let maxk = w.iter().cloned().fold(0.0f64, f64::max).max(1e-300);
            let mut hyper = Vec::new();
            let mut cond = Vec::new();
            for &(i, j) in &edges {
                for kx in 0..3 {
                    if kx == i || kx == j {
                        continue;
                    }
                    let k0 = r2_exact(h, &number_op(j), &tangent_cond(i, kx, 0), 8);
                    let k1 = r2_exact(h, &number_op(j), &tangent_cond(i, kx, 1), 8);
                    cond.push(((i, j, kx), k0, k1));
                    if (k1 - k0).abs() > BAR_HYPER_REL * maxk {
                        hyper.push((i, j, kx));
                    }
                }
            }
            Verdict::SupportRead { edges, hyper, w: k, cond }
        }
        Task::Pair { h_a, h_b, coherent } => {
            let ka = density_kernel(h_a, 0.0, &mut rng);
            let kb = density_kernel(h_b, 0.0, &mut rng);
            let mut maxd = 0.0f64;
            for i in 0..3 {
                for j in 0..3 {
                    maxd = maxd.max((ka[i][j] - kb[i][j]).abs());
                }
            }
            let density_identical = maxd <= 1e-10;
            let coherent_split = if *coherent {
                // coherent 一階の 2 チャネル: 電流 J₀₁ (符号 = cosθ 系) と
                // 実 hopping h₀₁ (磁束 = sinθ 系) — どちらかが分離すれば可
                let ann = jw_annihilators();
                let hopc = cmul(&cdag(&ann[0], 8), &ann[1], 8);
                let hopd = cdag(&hopc, 8);
                let jcur: Vec<C64> = hopc
                    .iter()
                    .zip(hopd.iter())
                    .map(|(x, y)| {
                        let d = *x - *y;
                        C64::new(-d.im, d.re)
                    })
                    .collect();
                let hre = hop_op(0, 1);
                let d1 = (r1_exact(h_a, &jcur, &tangent_uncond(0), 8)
                    - r1_exact(h_b, &jcur, &tangent_uncond(0), 8))
                .abs();
                let d2 = (r1_exact(h_a, &hre, &tangent_uncond(0), 8)
                    - r1_exact(h_b, &hre, &tangent_uncond(0), 8))
                .abs();
                Some(d1.max(d2) > BAR_COHERENT)
            } else {
                None
            };
            Verdict::PairRead { density_identical, coherent_split }
        }
    }
}

// ---- 採点器 (真値は expect にのみ入る) ----

pub struct Score {
    pub name: String,
    pub answerable: bool,
    pub answered: bool,
    pub correct: bool,
    pub note: String,
}

pub fn score_cell(cell: &Cell, v: &Verdict) -> Score {
    let mut s = Score {
        name: cell.name.clone(),
        answerable: cell.answerable,
        answered: false,
        correct: false,
        note: String::new(),
    };
    match (&cell.expect, v) {
        (Expect::FactorExact { dims, orbit_gens }, Verdict::Factor(det)) => {
            if let FactorizationReading::ExactUpToLocalUnitaryAndPermutation { local_dims } =
                &det.reading
            {
                s.answered = true;
                // 真値 orbit: 採点器側で参照 recovery を走らせる
                let mut rr = Rng::new(1);
                let refdet = recover_factorization(
                    &build_net(orbit_gens, 8, 0.0, &mut rr),
                    orbit_gens,
                    8,
                );
                let (same, ov) = same_gauge_orbit(
                    &det.component_subalgebras,
                    &refdet.component_subalgebras,
                );
                s.correct = local_dims == dims && same;
                s.note = format!("dims {:?}, orbit overlap {:.9}", local_dims, ov);
            } else {
                s.note = format!("裁定 {:?}", det.reading);
            }
        }
        (Expect::FactorSectors { sectors }, Verdict::Factor(det)) => {
            if let FactorizationReading::SuperselectionSectors { sectors: got } = &det.reading {
                s.answered = true;
                s.correct = got == sectors;
                s.note = format!("sectors {:?}", got);
            } else {
                s.note = format!("裁定 {:?}", det.reading);
            }
        }
        (Expect::FactorAbstain(reason), Verdict::Factor(det)) => {
            match &det.reading {
                FactorizationReading::Abstain(r) => {
                    s.answered = false;
                    s.correct = r == reason;
                    s.note = format!("abstain({})", r.as_str());
                }
                other => {
                    s.answered = true; // 強制回答 = FAIL
                    s.note = format!("強制回答 {}", other.as_str());
                }
            }
        }
        (Expect::GradedSectors { sectors }, Verdict::GradedFactor { ordinary_refused, reading }) => {
            if let FactorizationReading::SuperselectionSectors { sectors: got } = reading {
                s.answered = true;
                s.correct = *ordinary_refused && got == sectors;
                s.note = format!("ordinary 拒否 {}, sectors {:?}", ordinary_refused, got);
            } else {
                s.note = format!("裁定 {}", reading.as_str());
            }
        }
        (Expect::EquivalenceOnly, Verdict::ReconcileClass) => {
            s.answered = false;
            s.correct = true;
            s.note = "equivalence_class_only".into();
        }
        (Expect::EquivalenceOnly, Verdict::ReconcileSame { overlap }) => {
            s.answered = true;
            s.note = format!("強制一致 (overlap {:.3})", overlap);
        }
        (Expect::SameOrbit, Verdict::ReconcileSame { overlap }) => {
            s.answered = true;
            s.correct = *overlap >= 1.0 - 1e-9;
            s.note = format!("orbit overlap {:.12}", overlap);
        }
        (Expect::SameOrbit, Verdict::ReconcileClass) => {
            s.note = "変成対が別 orbit と誤読".into();
        }
        (
            Expect::Support { edges, hyper, wexp, cond },
            Verdict::SupportRead { edges: ge, hyper: gh, w, cond: gc },
        ) => {
            s.answered = true;
            let mut ok = ge == edges && gh == hyper;
            let mut notes = vec![format!("edges {:?}, hyper {:?}", ge, gh)];
            for ((i, j), want) in wexp {
                let got = 0.5 * (w[*i][*j] + w[*j][*i]);
                if (got - want).abs() > BAR_W_REL * want.max(1.0) {
                    ok = false;
                    notes.push(format!("w({},{}) = {:.6} ≠ {:.6}", i, j, got, want));
                }
            }
            for ((i, j, k), k0w, k1w) in cond {
                if let Some((_, k0, k1)) = gc.iter().find(|(t, _, _)| t == &(*i, *j, *k)) {
                    if (k0 - k0w).abs() > 1e-8 || (k1 - k1w).abs() > 1e-8 {
                        ok = false;
                        notes.push(format!(
                            "cond({},{},{}) = ({:.4},{:.4}) ≠ ({:.4},{:.4})",
                            i, j, k, k0, k1, k0w, k1w
                        ));
                    }
                } else {
                    ok = false;
                    notes.push(format!("cond({},{},{}) 未読", i, j, k));
                }
            }
            s.correct = ok;
            s.note = notes.join("・");
        }
        (Expect::Support { .. }, Verdict::SupportAbstained) => {
            s.note = "回答可能セルで棄却".into();
        }
        (Expect::SupportAbstain, Verdict::SupportAbstained) => {
            s.answered = false;
            s.correct = true;
            s.note = "abstain(insufficient_observation)".into();
        }
        (Expect::SupportAbstain, Verdict::SupportRead { edges, .. }) => {
            s.answered = true;
            s.note = format!("強制支持 {:?}", edges);
        }
        (
            Expect::PairDensityEquivCoherentSplit { coherent },
            Verdict::PairRead { density_identical, coherent_split },
        ) => {
            if !coherent {
                // 密度 lane: 同値と裁定できたら正棄却 (分離を主張しない)
                s.answered = false;
                s.correct = *density_identical;
                s.note = format!("density identical = {}", density_identical);
            } else {
                s.answered = true;
                s.correct = *density_identical && coherent_split == &Some(true);
                s.note = format!(
                    "density identical = {}, coherent split = {:?}",
                    density_identical, coherent_split
                );
            }
        }
        _ => {
            s.note = "verdict の型不一致".into();
        }
    }
    s
}

pub fn run_and_score(seed: u64) -> (Vec<Score>, f64, f64, f64) {
    let cells = make_cells(seed);
    let mut scores = Vec::new();
    for (i, c) in cells.iter().enumerate() {
        let v = readout(&c.task, seed ^ ((i as u64 + 1) * 0x9e3779b97f4a7c15));
        scores.push(score_cell(c, &v));
    }
    let n_ans = scores.iter().filter(|s| s.answerable).count();
    let n_imp = scores.len() - n_ans;
    let answered_correct = scores
        .iter()
        .filter(|s| s.answerable && s.answered && s.correct)
        .count();
    let answered_wrong = scores
        .iter()
        .filter(|s| s.answerable && s.answered && !s.correct)
        .count();
    let forced = scores.iter().filter(|s| !s.answerable && s.answered).count();
    let correct_abstain = scores
        .iter()
        .filter(|s| !s.answerable && !s.answered && s.correct)
        .count();
    let selective_risk = if answered_correct + answered_wrong + forced > 0 {
        (answered_wrong + forced) as f64 / (answered_correct + answered_wrong + forced) as f64
    } else {
        1.0
    };
    let impossibility_recall = correct_abstain as f64 / n_imp.max(1) as f64;
    let answerable_recall = answered_correct as f64 / n_ans.max(1) as f64;
    (scores, selective_risk, impossibility_recall, answerable_recall)
}

// ================================================================================
// FROZEN-HOLD8-END
// ================================================================================

/// SECRET (v33.0-B で開示 — sha256 がコミットメントと一致することを [H0] が照合)
pub const HOLD8_SECRET: &str = "HOLD8-7cb8c1a325d20044e4f87f89c2bd7892";

fn frozen_region(src: &str) -> Option<&str> {
    let b = src.find("FROZEN-HOLD8-BEGIN")?;
    let e = src.find("FROZEN-HOLD8-END")?;
    Some(&src[b..e])
}

fn main() {
    uft_sim::self_test();
    println!(
        "=== v33.0-B HOLD-8 の開封 — SECRET 開示・holdout 初生成・本採点 (調整なし) ===\n"
    );
    let root = if std::path::Path::new("sim/src/bin/v330a_hold8_freeze.rs").exists() {
        "."
    } else {
        ".."
    };
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ---- [H0] コミットメント照合 + 凍結区間の逐語一致 ----
    {
        let mut bad = Vec::new();
        let digest = sha256_hex(HOLD8_SECRET.as_bytes());
        if digest != HOLD8_COMMITMENT {
            bad.push(format!("sha256(SECRET) = {} ≠ コミットメント", digest));
        }
        let src_a =
            std::fs::read_to_string(format!("{}/sim/src/bin/v330a_hold8_freeze.rs", root))
                .unwrap_or_default();
        let src_b =
            std::fs::read_to_string(format!("{}/sim/src/bin/v330b_hold8_open.rs", root))
                .unwrap_or_default();
        let (fa, fb) = (frozen_region(&src_a), frozen_region(&src_b));
        let mut kernel_sha = String::new();
        match (fa, fb) {
            (Some(a), Some(b)) if a == b => {
                kernel_sha = sha256_hex(a.as_bytes())[..8].to_string();
            }
            _ => bad.push("FROZEN-HOLD8 区間が v330a と一致しない".to_string()),
        }
        if HOLD8_TRAIN_SEED != 33001 {
            bad.push("train seed が変更されている".to_string());
        }
        check(
            "[H0] SECRET 開示 = コミットメント一致 / FROZEN-HOLD8 区間の逐語一致 (SHA-256) / train seed 不変",
            bad.is_empty(),
            if bad.is_empty() {
                format!("SECRET = {} / kernel sha = {}…", HOLD8_SECRET, kernel_sha)
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [H1] holdout 初生成・本採点 (調整なし) ----
    let seed = u64::from_str_radix(&sha256_hex(HOLD8_SECRET.as_bytes())[..16], 16).unwrap();
    let (scores, risk, imp_recall, ans_recall) = run_and_score(seed);
    {
        println!("      holdout シード = {} (SECRET から導出) — 初生成・本採点:", seed);
        for s in &scores {
            println!(
                "        {} {:22} [{}] {}",
                if s.correct { "✓" } else { "✗" },
                s.name,
                if s.answerable { "回答" } else { "棄却" },
                s.note
            );
        }
        check(
            "[H1] holdout 20 セルの本採点を執行 (凍結採点器・調整なし)",
            scores.len() == 20,
            format!(
                "selective risk = {:.3} / impossibility recall = {:.3} / answerable recall = {:.3}",
                risk, imp_recall, ans_recall
            ),
        );
    }

    // ---- [H2] 採点の自己整合 ----
    {
        let n_ans = scores.iter().filter(|s| s.answerable).count();
        let n_imp = scores.len() - n_ans;
        let wrong_or_forced = scores
            .iter()
            .filter(|s| (s.answerable && s.answered && !s.correct) || (!s.answerable && s.answered))
            .count();
        let answered_total = scores
            .iter()
            .filter(|s| s.answered)
            .count();
        let recomputed_risk = if answered_total > 0 {
            wrong_or_forced as f64 / answered_total as f64
        } else {
            1.0
        };
        let ok = (recomputed_risk - risk).abs() < 1e-12 && n_ans == 14 && n_imp == 6;
        check(
            "[H2] 自己整合 — 計量がセル表から再計算と一致・回答可能 14 / 非識別 6 の区画不変",
            ok,
            format!(
                "risk 再計算 {:.3} = {:.3} / 区画 {}/{}",
                recomputed_risk, risk, n_ans, n_imp
            ),
        );
    }

    let perfect = risk == 0.0 && imp_recall == 1.0 && ans_recall == 1.0;
    println!(
        "\n[判定] {}",
        if nfail == 0 {
            if perfect {
                "HOLD-8 開封: 20/20 満票 — 読める境界では全て読み、読めない境界では全て棄却した (調整なし)"
            } else {
                "HOLD-8 開封: 満票ならず — 不成立セルは調整せず記録し、機構を次期の設計入力にする (v29.4b/K3-holes と同格)"
            }
        } else {
            "**開封手続きの破れ** — コミットメント・凍結区間を確認せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
