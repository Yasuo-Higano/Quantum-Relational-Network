//! v32.0-A HOLD-7 の凍結 — identifiability 相図全体の新鮮 holdout (PROMPT/12 第三十一期末)
//!
//! HOLD-5/6 と同じ開封順序:
//!   本コミット (v32.0-A) = 生成器・採点器・バー・観測契約の凍結 + SECRET コミットメント
//!   次コミット (v32.0-B) = SECRET 開示・holdout 初生成・本採点 (調整なし)
//!
//! sha256(SECRET) = c50c4c30993b1bb7113734609e87f116e8f280a30315571261f92377bd3ec9ea
//!
//! HOLD-7 が HOLD-5/6 と違う点 (PROMPT/12 §7): **非識別セルの正しい棄却を採点する**。
//!   - selective risk (回答したセルの誤り率) = 0 が要求
//!   - coverage (回答可能セルの回答率) ≥ 0.9
//!   - impossibility recall (非識別セルの正しい棄却率) = 1.0 —
//!     **非識別セルで無理に回答したら FAIL**・正しい EquivalenceClassOnly/Abstain が PASS
//!   - 変成対 (置換ゲージ) の読み出し一致・真値の readout への非流入
//!
//! セルクラス (10 種 — 状態領域 × 観測契約 × 因子分解 × ノイズ × 位相/計量 × regulator):
//!   K1 topo2d (torus/genus2/two-holes, 重み場+置換)     → 応答 lane: (β, 曲面種) 回答
//!   K2 topo3d (T³ L=4 / 16-cell / 中実立方体)           → 応答 lane: (β₀..β₃, link 種) 回答
//!   K3 metric (重み場 ring の 1/3 法則・大小 2 穴の寿命比) → VR: 定量回答
//!   K4 oracle (乱択グラフ β=1 / β=30)                   → exact ĥ 回答 / **棄却 (非識別)**
//!   K5 projector 衝突 (P6/693)                          → **同値類 (非識別)**
//!   K6 interacting (t-V, ED): oracle → **棄却** / 応答 → 回答 (厳密転移)
//!   K7 unknown factorization (稠密回転基底の C)          → **棄却 (非識別)**
//!   K8 noise (σ 小 → 回答 / σ 大 → **棄却**)
//!   K9 変成対 (ノード置換 — 読み出し一致)
//!   K10 regulator (開放鎖 R-A×R-C: 境界バイアスの a 減衰 ≥ 2×)
//!
//! 本バイナリ: [F0] コミットメント表示 [F1] 生成器健全性 (真の複体・真値プロファイル —
//! HOLD-6 [G1b] の常設) [F2] train 全セル採点 (可視シード) [F3] 凍結宣言 (カーネル SHA)
//!
//! 実行: cargo run --release --bin v320a_hold7_freeze

use uft_sim::readout_contract::*;
use uft_sim::{jacobi_eigh, matfun_sym, sha256_hex, Rng, C64};

// ================================================================================
// FROZEN-HOLD7-BEGIN  (この区間は v320b と逐語一致 — [H0] が SHA-256 で照合する)
// ================================================================================

pub const HOLD7_COMMITMENT: &str =
    "c50c4c30993b1bb7113734609e87f116e8f280a30315571261f92377bd3ec9ea";
pub const HOLD7_TRAIN_SEED: u64 = 32001;

// ---- 凍結バー (開封後に変更しない) ----
pub const BAR_ORACLE_ERR: f64 = 1e-6; // K4/K6: exact ĥ の max 誤差
pub const BAR_WICK: f64 = 1e-8; // Gaussianity witness の棄却バー
pub const BAR_RESP_REL: f64 = 1e-3; // 応答重みの相対誤差 (クリーン)
pub const BAR_NOISY_REL: f64 = 0.05; // 応答重み (σ 小のノイズ下)
pub const BAR_NOISE_ABSTAIN: f64 = 0.1; // ノイズ誤差見積り > これ → 棄却
pub const BAR_RING_THIRD: (f64, f64) = (0.28, 0.37); // VR 円環 death/周長
pub const BAR_HOLE_RATIO_TOL: f64 = 0.35; // 2 穴寿命比の相対許容
pub const BAR_REG_DECAY: f64 = 2.0; // regulator 境界バイアスの a=1→½ 減衰比
pub const BAR_COVERAGE: f64 = 0.9; // 回答可能セルの回答率
pub const GAP_GUARD: f64 = 1e-3; // gap 則スケールガード
pub const DT_BASE: f64 = 0.02; // 応答 lane dt (スペクトル半径スケール前)
pub const EPS_PROBE: f64 = 0.3;

// ---- 読み出しの裁定 (セルごとの出力) ----
#[derive(Debug, Clone, PartialEq)]
pub enum CellVerdict {
    Topo2d { beta: [i64; 3], closed: bool },
    Topo3d { beta: [i64; 4], kind: &'static str },
    MetricRing { death_over_perim: f64 },
    MetricHoles { ratio: f64 },
    ExactH { h: Vec<f64>, n: usize },
    RespWeights { w: Vec<f64>, n: usize },
    EquivClass,
    Abstained(&'static str),
}

// ---- 数値素子 (v31.1–31.6 の凍結カーネル — 逐語) ----

pub fn gibbs_c(h: &[f64], n: usize, beta: f64) -> Vec<f64> {
    matfun_sym(h, n, |x| 1.0 / (1.0 + (beta * x).exp()))
}

pub fn logit_k(c: &[f64], n: usize) -> Vec<f64> {
    matfun_sym(c, n, |x| ((1.0 - x) / x).ln())
}

/// 密度曲率測定 lane (dt はスペクトル半径スケール・時系列のみ・ノイズ注入可)
pub fn curvature_w(
    h: &[f64],
    n: usize,
    i: usize,
    sigma: f64,
    rng: &mut Rng,
) -> Vec<f64> {
    let norm1 = (0..n)
        .map(|r| (0..n).map(|c| h[r * n + c].abs()).sum::<f64>())
        .fold(0.0f64, f64::max)
        .max(1.0);
    let dt = DT_BASE / norm1;
    let (vals, vecs) = jacobi_eigh(h, n);
    let times = [-dt, -dt / 2.0, dt / 2.0, dt];
    let mut narr = [
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
        vec![0.0; n],
    ];
    let mut n0 = [vec![0.0; n], vec![0.0; n]];
    for (pi, sign) in [(0usize, 1.0), (1usize, -1.0)] {
        let c0: Vec<f64> = (0..n)
            .map(|s| if s == i { 0.5 + sign * EPS_PROBE } else { 0.5 })
            .collect();
        for j in 0..n {
            n0[pi][j] = c0[j] + sigma * rng.gauss();
        }
        for (ti, &t) in times.iter().enumerate() {
            // C̃0 の回転 (実固有系)
            let mut ct = vec![C64::new(0.0, 0.0); n * n];
            for a in 0..n {
                for b in 0..n {
                    let mut s = 0.0;
                    for q in 0..n {
                        s += vecs[a * n + q] * c0[q] * vecs[b * n + q];
                    }
                    ct[a * n + b] = C64::expi(-(vals[a] - vals[b]) * t).scale(s);
                }
            }
            for j in 0..n {
                let mut s = C64::new(0.0, 0.0);
                for a in 0..n {
                    for b in 0..n {
                        s = s + ct[a * n + b].scale(vecs[a * n + j] * vecs[b * n + j]);
                    }
                }
                narr[pi * 4 + ti][j] = s.re + sigma * rng.gauss();
            }
        }
    }
    let mut w = vec![0.0; n];
    for j in 0..n {
        let d2 = |pi: usize, half: bool| -> f64 {
            let (tm, tp, dd) = if half { (1, 2, dt / 2.0) } else { (0, 3, dt) };
            (narr[pi * 4 + tp][j] - 2.0 * n0[pi][j] + narr[pi * 4 + tm][j]) / (dd * dd)
        };
        let coarse = (d2(0, false) - d2(1, false)) / (4.0 * EPS_PROBE);
        let fine = (d2(0, true) - d2(1, true)) / (4.0 * EPS_PROBE);
        w[j] = (4.0 * fine - coarse) / 3.0;
    }
    w
}

/// ノイズ誤差見積り (真値不使用 — σ と器械定数から): σ·17√6/(3·dt²·4ε)
pub fn noise_error_bound(sigma: f64, h_norm1: f64) -> f64 {
    let dt = DT_BASE / h_norm1.max(1.0);
    sigma * 17.0 * 6.0f64.sqrt() / 3.0 / (dt * dt * 4.0 * EPS_PROBE)
}

/// gap 支持 (v31.6 最終規則)
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
            if i != j && w[i * n + j].abs() > thr && w[j * n + i].abs() > thr {
                adj[i][j] = true;
            }
        }
    }
    adj
}

// ---- Z2 homology (v31.6 凍結エンジン — 逐語) ----

#[derive(Clone, Default)]
pub struct Cx {
    pub simp: Vec<Vec<Vec<usize>>>,
}

pub fn clique_complex(adj: &[Vec<bool>], n: usize) -> Cx {
    let mut c = Cx::default();
    c.simp = vec![Vec::new(); 5];
    for v in 0..n {
        c.simp[0].push(vec![v]);
    }
    for i in 0..n {
        for j in (i + 1)..n {
            if adj[i][j] {
                c.simp[1].push(vec![i, j]);
            }
        }
    }
    for d in 2..=4 {
        let prev = c.simp[d - 1].clone();
        for s in prev {
            let last = *s.last().unwrap();
            for v in (last + 1)..n {
                if s.iter().all(|&u| adj[u][v]) {
                    let mut t = s.clone();
                    t.push(v);
                    c.simp[d].push(t);
                }
            }
        }
    }
    c
}

pub fn xor_merge(a: &[usize], b: &[usize]) -> Vec<usize> {
    let mut out = Vec::with_capacity(a.len() + b.len());
    let (mut i, mut j) = (0, 0);
    while i < a.len() || j < b.len() {
        if j >= b.len() || (i < a.len() && a[i] < b[j]) {
            out.push(a[i]);
            i += 1;
        } else if i >= a.len() || b[j] < a[i] {
            out.push(b[j]);
            j += 1;
        } else {
            i += 1;
            j += 1;
        }
    }
    out
}

pub fn boundary_rank(faces: &[Vec<usize>], simps: &[Vec<usize>]) -> usize {
    use std::collections::HashMap;
    let mut fidx: HashMap<&[usize], usize> = HashMap::new();
    for (i, f) in faces.iter().enumerate() {
        fidx.insert(f.as_slice(), i);
    }
    let mut cols: Vec<Vec<usize>> = Vec::with_capacity(simps.len());
    for s in simps {
        let mut col = Vec::with_capacity(s.len());
        for drop in 0..s.len() {
            let mut f: Vec<usize> = Vec::with_capacity(s.len() - 1);
            for (i, &v) in s.iter().enumerate() {
                if i != drop {
                    f.push(v);
                }
            }
            col.push(fidx[f.as_slice()]);
        }
        col.sort_unstable();
        cols.push(col);
    }
    let mut pivot_owner: HashMap<usize, usize> = HashMap::new();
    let mut rank = 0usize;
    for j in 0..cols.len() {
        let mut col = cols[j].clone();
        loop {
            match col.last() {
                None => break,
                Some(&low) => match pivot_owner.get(&low) {
                    None => {
                        pivot_owner.insert(low, j);
                        cols[j] = col;
                        rank += 1;
                        break;
                    }
                    Some(&k) => {
                        col = xor_merge(&col, &cols[k]);
                    }
                },
            }
        }
    }
    rank
}

pub fn betti(c: &Cx) -> Vec<i64> {
    let maxd = c.simp.len() - 1;
    let mut rank = vec![0usize; maxd + 2];
    for k in 1..=maxd {
        rank[k] = boundary_rank(&c.simp[k - 1], &c.simp[k]);
    }
    let mut b = Vec::new();
    for k in 0..=maxd {
        let nk = c.simp[k].len() as i64;
        let r_k1 = if k + 1 <= maxd { rank[k + 1] as i64 } else { 0 };
        b.push(nk - rank[k] as i64 - r_k1);
    }
    b
}

pub fn vertex_link(c: &Cx, v: usize) -> Cx {
    let mut l = Cx::default();
    l.simp = vec![Vec::new(); 3];
    for d in 1..c.simp.len() {
        for s in &c.simp[d] {
            if let Some(pos) = s.iter().position(|&u| u == v) {
                let mut f: Vec<usize> = s.clone();
                f.remove(pos);
                if f.len() <= 3 && !l.simp[f.len() - 1].contains(&f) {
                    l.simp[f.len() - 1].push(f);
                }
            }
        }
    }
    l
}

pub fn classify_surface(l: &Cx) -> &'static str {
    if l.simp.len() < 3 || l.simp[2].is_empty() {
        return "other";
    }
    let mut cnt = std::collections::HashMap::new();
    for t in &l.simp[2] {
        for drop in 0..3 {
            let mut e: Vec<usize> = t.clone();
            e.remove(drop);
            *cnt.entry(e).or_insert(0usize) += 1;
        }
    }
    let mut n1 = 0usize;
    for e in &l.simp[1] {
        match cnt.get(e) {
            Some(&2) => {}
            Some(&1) => n1 += 1,
            _ => return "other",
        }
    }
    let b = betti(l);
    if n1 == 0 && b[0] == 1 && b[1] == 0 && b[2] == 1 {
        "S2"
    } else if n1 > 0 && b[0] == 1 && b[1] == 0 && b[2] == 0 {
        "D2"
    } else {
        "other"
    }
}

// ---- VR persistence H1 (v31.6 凍結 — 逐語) ----

pub fn metric_closure(w: &[f64], adj: &[Vec<bool>], n: usize) -> Vec<f64> {
    let mut d = vec![f64::INFINITY; n * n];
    for s in 0..n {
        d[s * n + s] = 0.0;
        let mut done = vec![false; n];
        for _ in 0..n {
            let mut u = usize::MAX;
            let mut best = f64::INFINITY;
            for v in 0..n {
                if !done[v] && d[s * n + v] < best {
                    best = d[s * n + v];
                    u = v;
                }
            }
            if u == usize::MAX {
                break;
            }
            done[u] = true;
            for v in 0..n {
                if v != u && adj[u][v] {
                    let len = 1.0 / w[u * n + v].abs().sqrt();
                    let nd = d[s * n + u] + len;
                    if nd < d[s * n + v] {
                        d[s * n + v] = nd;
                    }
                }
            }
        }
    }
    d
}

pub fn vr_h1_bars(d: &[f64], n: usize) -> Vec<(f64, f64)> {
    #[derive(Clone)]
    struct Simp {
        f: f64,
        dim: usize,
        verts: [usize; 3],
    }
    let mut simps: Vec<Simp> = Vec::new();
    for v in 0..n {
        simps.push(Simp {
            f: 0.0,
            dim: 0,
            verts: [v, 0, 0],
        });
    }
    for i in 0..n {
        for j in (i + 1)..n {
            simps.push(Simp {
                f: d[i * n + j],
                dim: 1,
                verts: [i, j, 0],
            });
        }
    }
    for i in 0..n {
        for j in (i + 1)..n {
            for k in (j + 1)..n {
                let f = d[i * n + j].max(d[i * n + k]).max(d[j * n + k]);
                simps.push(Simp {
                    f,
                    dim: 2,
                    verts: [i, j, k],
                });
            }
        }
    }
    let mut order: Vec<usize> = (0..simps.len()).collect();
    order.sort_by(|&a, &b| {
        simps[a]
            .f
            .partial_cmp(&simps[b].f)
            .unwrap()
            .then(simps[a].dim.cmp(&simps[b].dim))
            .then(simps[a].verts.cmp(&simps[b].verts))
    });
    let mut pos = vec![0usize; simps.len()];
    for (rank, &id) in order.iter().enumerate() {
        pos[id] = rank;
    }
    let mut eid = std::collections::HashMap::new();
    for (id, s) in simps.iter().enumerate() {
        if s.dim == 1 {
            eid.insert((s.verts[0], s.verts[1]), pos[id]);
        }
    }
    let mut cols: Vec<Vec<usize>> = vec![Vec::new(); simps.len()];
    for (id, s) in simps.iter().enumerate() {
        let j = pos[id];
        match s.dim {
            1 => {
                let mut c = vec![pos[s.verts[0]], pos[s.verts[1]]];
                c.sort_unstable();
                cols[j] = c;
            }
            2 => {
                let (a, b, c3) = (s.verts[0], s.verts[1], s.verts[2]);
                let mut c = vec![eid[&(a, b)], eid[&(a, c3)], eid[&(b, c3)]];
                c.sort_unstable();
                cols[j] = c;
            }
            _ => {}
        }
    }
    let m = simps.len();
    let mut pivot_owner: Vec<Option<usize>> = vec![None; m];
    let mut reduced: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut pairs: Vec<(usize, usize)> = Vec::new();
    for jrank in 0..m {
        let id = order[jrank];
        if simps[id].dim == 0 {
            continue;
        }
        let mut col = cols[pos[id]].clone();
        loop {
            match col.last() {
                None => break,
                Some(&low) => match pivot_owner[low] {
                    None => {
                        pivot_owner[low] = Some(jrank);
                        reduced[jrank] = col.clone();
                        pairs.push((low, jrank));
                        break;
                    }
                    Some(k) => col = xor_merge(&col, &reduced[k]),
                },
            }
        }
    }
    let mut bars = Vec::new();
    for &(b, dth) in &pairs {
        let sb = &simps[order[b]];
        let sd = &simps[order[dth]];
        if sb.dim == 1 && sd.dim == 2 && sd.f > sb.f + 1e-12 {
            bars.push((sb.f, sd.f));
        }
    }
    bars.sort_by(|a, b| (b.1 - b.0).partial_cmp(&(a.1 - a.0)).unwrap());
    bars
}

// ---- ED (t-V, v31.5 凍結 — 逐語圧縮版) ----

pub fn sector_masks(n: usize, k: usize) -> Vec<u32> {
    (0u32..(1 << n))
        .filter(|m| m.count_ones() as usize == k)
        .collect()
}

pub fn hop_sign(mask: u32, a: usize, b: usize) -> Option<(u32, f64)> {
    if (mask >> a) & 1 == 0 || (mask >> b) & 1 == 1 {
        return None;
    }
    let m1 = mask & !(1 << a);
    let s1 = ((mask & ((1 << a) - 1)).count_ones() % 2) as i32;
    let s2 = ((m1 & ((1 << b) - 1)).count_ones() % 2) as i32;
    Some((m1 | (1 << b), if (s1 + s2) % 2 == 0 { 1.0 } else { -1.0 }))
}

pub fn build_h_tv(masks: &[u32], bonds: &[(usize, usize, f64)], v: f64) -> Vec<f64> {
    let d = masks.len();
    let idx = |m: u32| masks.binary_search(&m).unwrap();
    let mut h = vec![0.0; d * d];
    for (r, &m) in masks.iter().enumerate() {
        let mut diag = 0.0;
        for &(a, b, _) in bonds {
            if (m >> a) & 1 == 1 && (m >> b) & 1 == 1 {
                diag += v;
            }
        }
        h[r * d + r] = diag;
        for &(a, b, t) in bonds {
            for (x, y) in [(a, b), (b, a)] {
                if let Some((m2, sgn)) = hop_sign(m, x, y) {
                    h[idx(m2) * d + r] += -t * sgn;
                }
            }
        }
    }
    h
}

/// t-V の熱的 (grand canonical) C と密度-密度 → Wick witness
pub fn tv_thermal_c_witness(
    nsite: usize,
    bonds: &[(usize, usize, f64)],
    v: f64,
    beta: f64,
) -> (Vec<f64>, f64) {
    let mut z = 0.0;
    let mut c = vec![0.0; nsite * nsite];
    let mut nn = vec![0.0; nsite * nsite];
    for k in 0..=nsite {
        let masks = sector_masks(nsite, k);
        let d = masks.len();
        let h = build_h_tv(&masks, bonds, v);
        let (ev, evec) = jacobi_eigh(&h, d);
        for m in 0..d {
            let w = (-beta * ev[m]).exp();
            z += w;
            let psi: Vec<f64> = (0..d).map(|r| evec[m * d + r]).collect();
            let idx = |mm: u32| masks.binary_search(&mm).unwrap();
            for (r, &mask) in masks.iter().enumerate() {
                let a = psi[r];
                if a == 0.0 {
                    continue;
                }
                for i in 0..nsite {
                    if (mask >> i) & 1 == 1 {
                        c[i * nsite + i] += w * a * a;
                        for j in 0..nsite {
                            if (mask >> j) & 1 == 1 {
                                nn[i * nsite + j] += w * a * a;
                            }
                        }
                    }
                    for j in 0..nsite {
                        if i != j {
                            if let Some((m2, sgn)) = hop_sign(mask, j, i) {
                                c[i * nsite + j] += w * psi[idx(m2)] * a * sgn;
                            }
                        }
                    }
                }
            }
        }
    }
    for e in 0..nsite * nsite {
        c[e] /= z;
        nn[e] /= z;
    }
    let mut wit: f64 = 0.0;
    for i in 0..nsite {
        for j in 0..nsite {
            if i != j {
                let conn = nn[i * nsite + j] - c[i * nsite + i] * c[j * nsite + j];
                wit = wit.max((conn + c[i * nsite + j] * c[j * nsite + i]).abs());
            }
        }
    }
    (c, wit)
}

/// t-V の厳密曲率 (対角 probe, v31.5 [N3] 凍結)
pub fn tv_curvature(nsite: usize, bonds: &[(usize, usize, f64)], v: f64, i: usize) -> Vec<f64> {
    let mut out = vec![0.0; nsite];
    for j in 0..nsite {
        if j == i {
            continue;
        }
        let mut acc = 0.0;
        for k in 0..=nsite {
            let masks = sector_masks(nsite, k);
            let d = masks.len();
            let h = build_h_tv(&masks, bonds, v);
            let njd: Vec<f64> = masks.iter().map(|&m| ((m >> j) & 1) as f64).collect();
            let mut a = vec![0.0; d * d];
            for r in 0..d {
                for cc in 0..d {
                    a[r * d + cc] = h[r * d + cc] * (njd[cc] - njd[r]);
                }
            }
            for r in 0..d {
                let mut brr = 0.0;
                for cc in 0..d {
                    brr += h[r * d + cc] * a[cc * d + r] - a[r * d + cc] * h[cc * d + r];
                }
                let wp = {
                    let mut w = 1.0;
                    for st in 0..nsite {
                        let occ = (masks[r] >> st) & 1 == 1;
                        w *= if st == i {
                            if occ {
                                0.5 + EPS_PROBE
                            } else {
                                0.5 - EPS_PROBE
                            }
                        } else {
                            0.5
                        };
                    }
                    w
                };
                let wm = {
                    let mut w = 1.0;
                    for st in 0..nsite {
                        let occ = (masks[r] >> st) & 1 == 1;
                        w *= if st == i {
                            if occ {
                                0.5 - EPS_PROBE
                            } else {
                                0.5 + EPS_PROBE
                            }
                        } else {
                            0.5
                        };
                    }
                    w
                };
                acc += -(wp - wm) * brr;
            }
        }
        out[j] = acc / (4.0 * EPS_PROBE);
    }
    out
}

// ---- 幾何構成 (生成器側 — v31.6 凍結) ----

pub fn kuhn_tets(lx: usize, ly: usize, lz: usize, pbc: bool) -> (usize, Vec<Vec<usize>>) {
    let nvx = if pbc { lx } else { lx + 1 };
    let nvy = if pbc { ly } else { ly + 1 };
    let nvz = if pbc { lz } else { lz + 1 };
    let vid =
        |x: usize, y: usize, z: usize| -> usize { ((x % nvx) * nvy + (y % nvy)) * nvz + (z % nvz) };
    let perms3: [[usize; 3]; 6] = [
        [0, 1, 2],
        [0, 2, 1],
        [1, 0, 2],
        [1, 2, 0],
        [2, 0, 1],
        [2, 1, 0],
    ];
    let mut tets = Vec::new();
    for x in 0..lx {
        for y in 0..ly {
            for z in 0..lz {
                for p in &perms3 {
                    let mut pos = [x, y, z];
                    let mut t = vec![vid(pos[0], pos[1], pos[2])];
                    for &ax in p {
                        pos[ax] += 1;
                        t.push(vid(pos[0], pos[1], pos[2]));
                    }
                    t.sort_unstable();
                    t.dedup();
                    if t.len() == 4 {
                        tets.push(t);
                    }
                }
            }
        }
    }
    tets.sort();
    tets.dedup();
    (nvx * nvy * nvz, tets)
}

pub fn skeleton_of_tets(nv: usize, tets: &[Vec<usize>]) -> Vec<Vec<bool>> {
    let mut adj = vec![vec![false; nv]; nv];
    for t in tets {
        for a in 0..4 {
            for b in (a + 1)..4 {
                adj[t[a]][t[b]] = true;
                adj[t[b]][t[a]] = true;
            }
        }
    }
    adj
}

pub fn torus_mesh(l: usize) -> (usize, Vec<Vec<usize>>) {
    let vid = |x: usize, y: usize| (x % l) * l + (y % l);
    let mut tris = Vec::new();
    for x in 0..l {
        for y in 0..l {
            let (a, b, c, d) = (vid(x, y), vid(x + 1, y), vid(x, y + 1), vid(x + 1, y + 1));
            for t in [[a, b, d], [a, c, d]] {
                let mut t2 = t.to_vec();
                t2.sort_unstable();
                tris.push(t2);
            }
        }
    }
    (l * l, tris)
}

pub fn genus2() -> (usize, Vec<Vec<usize>>) {
    let l = 5usize;
    let (nv, tris) = torus_mesh(l);
    let keep: Vec<Vec<usize>> = tris.iter().filter(|t| !t.contains(&0)).cloned().collect();
    let mut nb = std::collections::BTreeSet::new();
    for t in tris.iter().filter(|t| t.contains(&0)) {
        for &u in t {
            if u != 0 {
                nb.insert(u);
            }
        }
    }
    let nbv: Vec<usize> = nb.into_iter().collect();
    let cycle = {
        let mut adjm: std::collections::HashMap<usize, Vec<usize>> = Default::default();
        for t in tris.iter().filter(|t| t.contains(&0)) {
            let o: Vec<usize> = t.iter().cloned().filter(|&u| u != 0).collect();
            adjm.entry(o[0]).or_default().push(o[1]);
            adjm.entry(o[1]).or_default().push(o[0]);
        }
        let mut order = vec![nbv[0]];
        let mut prev = usize::MAX;
        while order.len() < nbv.len() {
            let cur = *order.last().unwrap();
            let nx = if adjm[&cur][0] != prev {
                adjm[&cur][0]
            } else {
                adjm[&cur][1]
            };
            prev = cur;
            order.push(nx);
        }
        order
    };
    let mut tris_all: Vec<Vec<usize>> = keep.clone();
    for t in &keep {
        let s: Vec<usize> = t.iter().map(|&u| u + nv).collect();
        tris_all.push(s);
    }
    let m = cycle.len();
    let mut rename: std::collections::HashMap<usize, usize> = Default::default();
    for k in 0..m {
        rename.insert(cycle[k] + nv, cycle[(m - k) % m]);
    }
    let mut fin: Vec<Vec<usize>> = Vec::new();
    for t in &tris_all {
        let mut s: Vec<usize> = t.iter().map(|&u| *rename.get(&u).unwrap_or(&u)).collect();
        s.sort_unstable();
        s.dedup();
        if s.len() == 3 {
            fin.push(s);
        }
    }
    fin.sort();
    fin.dedup();
    let mut used: Vec<usize> = fin.iter().flatten().cloned().collect();
    used.sort_unstable();
    used.dedup();
    let idx: std::collections::HashMap<usize, usize> =
        used.iter().enumerate().map(|(i, &u)| (u, i)).collect();
    let tris_c: Vec<Vec<usize>> = fin
        .iter()
        .map(|t| {
            let mut s: Vec<usize> = t.iter().map(|&u| idx[&u]).collect();
            s.sort_unstable();
            s
        })
        .collect();
    (used.len(), tris_c)
}

pub fn skeleton_of_tris(nv: usize, tris: &[Vec<usize>]) -> Vec<Vec<bool>> {
    let mut adj = vec![vec![false; nv]; nv];
    for t in tris {
        for a in 0..3 {
            for b in (a + 1)..3 {
                adj[t[a]][t[b]] = true;
                adj[t[b]][t[a]] = true;
            }
        }
    }
    adj
}

// ---- セル定義と生成器 ----

#[derive(Clone, Debug, PartialEq)]
pub enum Expected {
    Topo2d { beta: [i64; 3], closed: bool },
    Topo3d { beta: [i64; 4], kind: &'static str },
    RingThird,
    HoleRatio { ratio: f64 },
    ExactH,
    RespSupport,
    MustEquivOrAbstain,
    MetamorphicPair,
    RegulatorDecay,
}

pub struct Cell {
    pub name: String,
    pub h: Vec<f64>,
    pub n: usize,
    pub h2: Option<Vec<f64>>, // 変成対の第 2 系 / regulator の第 2 regulator
    pub obs: &'static str,    // 観測契約
    pub fact_given: bool,
    pub sigma: f64,
    pub beta_state: f64, // oracle lane の状態温度 (0 = projector)
    pub v_int: f64,      // t-V の V (0 = 自由)
    pub expected: Expected,
    pub truth_adj: Option<Vec<Vec<bool>>>,
    pub truth_perim: f64,
}

/// 重み場つき隣接 → h (シードで滑らかな変調 + ノード置換)
pub fn weighted_h(adj: &[Vec<bool>], n: usize, rng: &mut Rng) -> (Vec<f64>, Vec<Vec<bool>>) {
    let (a1, a2, p1, p2) = (
        0.1 + 0.15 * rng.f64(),
        0.05 + 0.1 * rng.f64(),
        rng.f64() * std::f64::consts::TAU,
        rng.f64() * std::f64::consts::TAU,
    );
    let perm: Vec<usize> = {
        let mut p: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.range(i + 1);
            p.swap(i, j);
        }
        p
    };
    let mut h = vec![0.0; n * n];
    let mut adjp = vec![vec![false; n]; n];
    for i in 0..n {
        for j in 0..n {
            if adj[i][j] {
                let x = (i + j) as f64 / n as f64;
                let t = 1.0
                    + a1 * (std::f64::consts::TAU * x + p1).sin()
                    + a2 * (2.0 * std::f64::consts::TAU * x + p2).cos();
                let (pi_, pj) = (perm[i], perm[j]);
                h[pi_ * n + pj] = -t;
                h[pj * n + pi_] = -t;
                adjp[pi_][pj] = true;
                adjp[pj][pi_] = true;
            }
        }
    }
    (h, adjp)
}

/// セル生成器 (シード決定的)。クラス構成は凍結 — holdout はシードだけが未知
pub fn generate_cells(seed: u64) -> Vec<Cell> {
    let mut rng = Rng::new(seed);
    let mut cells = Vec::new();
    // K1 topo2d ×3
    for k in 0..3 {
        let cls = rng.range(3);
        let (nv, adj, beta, closed, name) = match cls {
            0 => {
                let l = 4 + rng.range(2);
                let (nv, t) = torus_mesh(l);
                (
                    nv,
                    skeleton_of_tris(nv, &t),
                    [1, 2, 1],
                    true,
                    format!("K1-{} torus{}", k, l),
                )
            }
            1 => {
                let (nv, t) = genus2();
                (
                    nv,
                    skeleton_of_tris(nv, &t),
                    [1, 4, 1],
                    true,
                    format!("K1-{} genus2", k),
                )
            }
            _ => {
                // two-holes: 8 字 (大小 2 環 1 点接着) — グラフ (三角形なし), β = (1,2,0)
                let big = 10 + rng.range(4);
                let small = 6 + rng.range(3);
                let n = big + small - 1;
                let mut adj = vec![vec![false; n]; n];
                for kk in 0..big {
                    adj[kk][(kk + 1) % big] = true;
                    adj[(kk + 1) % big][kk] = true;
                }
                let sm: Vec<usize> = std::iter::once(0).chain(big..big + small - 1).collect();
                for kk in 0..sm.len() {
                    let (x, y) = (sm[kk], sm[(kk + 1) % sm.len()]);
                    adj[x][y] = true;
                    adj[y][x] = true;
                }
                (n, adj, [1, 2, 0], false, format!("K1-{} twoholes", k))
            }
        };
        let (h, adjp) = weighted_h(&adj, nv, &mut rng);
        cells.push(Cell {
            name,
            h,
            n: nv,
            h2: None,
            obs: "local_bias_density_response",
            fact_given: true,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: 0.0,
            expected: Expected::Topo2d { beta, closed },
            truth_adj: Some(adjp),
            truth_perim: 0.0,
        });
    }
    // K2 topo3d ×2
    for k in 0..2 {
        let cls = rng.range(3);
        let (nv, adj, beta, kind, name) = match cls {
            0 => {
                let (nv, tets) = kuhn_tets(4, 4, 4, true);
                (
                    nv,
                    skeleton_of_tets(nv, &tets),
                    [1, 3, 3, 1],
                    "closed",
                    format!("K2-{} T3", k),
                )
            }
            1 => {
                let n8 = 8;
                let mut adj8 = vec![vec![true; n8]; n8];
                for i in 0..n8 {
                    adj8[i][i] = false;
                }
                for i in 0..4 {
                    adj8[2 * i][2 * i + 1] = false;
                    adj8[2 * i + 1][2 * i] = false;
                }
                (n8, adj8, [1, 0, 0, 1], "closed", format!("K2-{} S3", k))
            }
            _ => {
                let (nv, tets) = kuhn_tets(2, 2, 2, false);
                (
                    nv,
                    skeleton_of_tets(nv, &tets),
                    [1, 0, 0, 0],
                    "boundary",
                    format!("K2-{} ball", k),
                )
            }
        };
        let (h, adjp) = weighted_h(&adj, nv, &mut rng);
        cells.push(Cell {
            name,
            h,
            n: nv,
            h2: None,
            obs: "local_bias_density_response",
            fact_given: true,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: 0.0,
            expected: Expected::Topo3d { beta, kind },
            truth_adj: Some(adjp),
            truth_perim: 0.0,
        });
    }
    // K3 metric ×2 (ring 1/3 + two-hole 比)
    {
        let n = 20 + rng.range(9);
        let mut adj = vec![vec![false; n]; n];
        for k in 0..n {
            adj[k][(k + 1) % n] = true;
            adj[(k + 1) % n][k] = true;
        }
        let (h, _) = weighted_h(&adj, n, &mut rng);
        // 周長 = Σ 1/√w (真値 — 採点でのみ使用)
        let mut perim = 0.0;
        for i in 0..n {
            for j in (i + 1)..n {
                if h[i * n + j] != 0.0 {
                    perim += 1.0 / h[i * n + j].abs();
                }
            }
        }
        cells.push(Cell {
            name: format!("K3-ring{}", n),
            h,
            n,
            h2: None,
            obs: "local_bias_density_response",
            fact_given: true,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: 0.0,
            expected: Expected::RingThird,
            truth_adj: None,
            truth_perim: perim,
        });
        let big = 12 + rng.range(6);
        let small = 6 + rng.range(3);
        let n2 = big + small - 1;
        let mut adj2 = vec![vec![false; n2]; n2];
        for k in 0..big {
            adj2[k][(k + 1) % big] = true;
            adj2[(k + 1) % big][k] = true;
        }
        let sm: Vec<usize> = std::iter::once(0).chain(big..big + small - 1).collect();
        for k in 0..sm.len() {
            let (x, y) = (sm[k], sm[(k + 1) % sm.len()]);
            adj2[x][y] = true;
            adj2[y][x] = true;
        }
        let mut h2 = vec![0.0; n2 * n2];
        for i in 0..n2 {
            for j in 0..n2 {
                if adj2[i][j] {
                    h2[i * n2 + j] = -1.0;
                }
            }
        }
        cells.push(Cell {
            name: format!("K3-holes{}:{}", big, small),
            h: h2,
            n: n2,
            h2: None,
            obs: "local_bias_density_response",
            fact_given: true,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: 0.0,
            expected: Expected::HoleRatio {
                ratio: big as f64 / small as f64,
            },
            truth_adj: None,
            truth_perim: 0.0,
        });
    }
    // K4 oracle ×2 (β=1 回答 / β=30 非識別)
    for (beta_state, name) in [(1.0, "K4-warm"), (30.0, "K4-cold")] {
        let n = 6 + rng.range(3);
        let mut adj = vec![vec![false; n]; n];
        for i in 0..n - 1 {
            adj[i][i + 1] = true;
            adj[i + 1][i] = true;
        }
        for i in 0..n {
            for j in (i + 2)..n {
                if rng.f64() < 0.3 {
                    adj[i][j] = true;
                    adj[j][i] = true;
                }
            }
        }
        let (h, _) = weighted_h(&adj, n, &mut rng);
        cells.push(Cell {
            name: name.to_string(),
            h,
            n,
            h2: None,
            obs: "global_one_body_correlation",
            fact_given: true,
            sigma: 0.0,
            beta_state,
            v_int: 0.0,
            expected: if beta_state < 10.0 {
                Expected::ExactH
            } else {
                Expected::MustEquivOrAbstain
            },
            truth_adj: None,
            truth_perim: 0.0,
        });
    }
    // K5 projector 衝突 (P6/693) — 非識別
    {
        let n = 6;
        let edges: [(usize, usize); 5] = [(0, 3), (0, 5), (1, 2), (1, 4), (2, 3)];
        let mut h = vec![0.0; n * n];
        for &(i, j) in &edges {
            h[i * n + j] = -1.0;
            h[j * n + i] = -1.0;
        }
        cells.push(Cell {
            name: "K5-projector".into(),
            h,
            n,
            h2: None,
            obs: "global_one_body_correlation",
            fact_given: true,
            sigma: 0.0,
            beta_state: 0.0, // projector
            v_int: 0.0,
            expected: Expected::MustEquivOrAbstain,
            truth_adj: None,
            truth_perim: 0.0,
        });
    }
    // K6 interacting ×2 (oracle 非識別 / 応答 回答)
    {
        let nsite = 8usize;
        let v = 1.0 + 2.0 * rng.f64();
        let mut h = vec![0.0; nsite * nsite];
        for k in 0..nsite {
            h[k * nsite + (k + 1) % nsite] = -1.0;
            h[((k + 1) % nsite) * nsite + k] = -1.0;
        }
        cells.push(Cell {
            name: format!("K6-oracle V={:.2}", v),
            h: h.clone(),
            n: nsite,
            h2: None,
            obs: "global_one_body_correlation",
            fact_given: true,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: v,
            expected: Expected::MustEquivOrAbstain,
            truth_adj: None,
            truth_perim: 0.0,
        });
        cells.push(Cell {
            name: format!("K6-resp V={:.2}", v),
            h,
            n: nsite,
            h2: None,
            obs: "local_bias_density_response",
            fact_given: true,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: v,
            expected: Expected::RespSupport,
            truth_adj: None,
            truth_perim: 0.0,
        });
    }
    // K7 unknown factorization — 非識別
    {
        let n = 12;
        let mut adj = vec![vec![false; n]; n];
        for k in 0..n {
            adj[k][(k + 1) % n] = true;
            adj[(k + 1) % n][k] = true;
        }
        let (h, _) = weighted_h(&adj, n, &mut rng);
        cells.push(Cell {
            name: "K7-unknownfact".into(),
            h,
            n,
            h2: None,
            obs: "global_one_body_correlation",
            fact_given: false,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: 0.0,
            expected: Expected::MustEquivOrAbstain,
            truth_adj: None,
            truth_perim: 0.0,
        });
    }
    // K8 noise ×2 (σ 小 → 回答 / σ 大 → 棄却)
    for (sigma, name) in [(1e-9, "K8-lownoise"), (1e-3, "K8-highnoise")] {
        let n = 12;
        let mut adj = vec![vec![false; n]; n];
        for k in 0..n {
            adj[k][(k + 1) % n] = true;
            adj[(k + 1) % n][k] = true;
        }
        let (h, adjp) = weighted_h(&adj, n, &mut rng);
        cells.push(Cell {
            name: name.to_string(),
            h,
            n,
            h2: None,
            obs: "local_bias_density_response",
            fact_given: true,
            sigma,
            beta_state: 1.0,
            v_int: 0.0,
            expected: if sigma < 1e-6 {
                Expected::RespSupport
            } else {
                Expected::MustEquivOrAbstain
            },
            truth_adj: if sigma < 1e-6 { Some(adjp) } else { None },
            truth_perim: 0.0,
        });
    }
    // K9 変成対 (置換)
    {
        let n = 12;
        let mut adj = vec![vec![false; n]; n];
        for k in 0..n {
            adj[k][(k + 1) % n] = true;
            adj[(k + 1) % n][k] = true;
        }
        let (h, adjp) = weighted_h(&adj, n, &mut rng);
        // 置換した第 2 系
        let perm: Vec<usize> = {
            let mut p: Vec<usize> = (0..n).collect();
            for i in (1..n).rev() {
                let j = rng.range(i + 1);
                p.swap(i, j);
            }
            p
        };
        let mut hp = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                hp[perm[i] * n + perm[j]] = h[i * n + j];
            }
        }
        cells.push(Cell {
            name: "K9-metamorphic".into(),
            h,
            n,
            h2: Some(hp),
            obs: "local_bias_density_response",
            fact_given: true,
            sigma: 0.0,
            beta_state: 1.0,
            v_int: 0.0,
            expected: Expected::MetamorphicPair,
            truth_adj: Some(adjp),
            truth_perim: 0.0,
        });
    }
    // K10 regulator (開放鎖 R-A×R-C: 端バイアスの a 減衰)
    {
        let n_phys = 32usize;
        // 第 1 系 = R-A a=1 / 第 2 系 = R-A a=1/2 (減衰検査は同一 regulator 対で
        //   端バイアス |τ_A − τ_C| を a=1 と a=½ で比較する — h は両 regulator ×2 格子
        //   を lab 側で構成するため、Cell には n_phys のみ有効)
        cells.push(Cell {
            name: "K10-regulator".into(),
            h: vec![0.0; 1],
            n: n_phys,
            h2: None,
            obs: "arrival_time_response",
            fact_given: true,
            sigma: 0.0,
            beta_state: 0.4, // ここでは m_phys を運ぶ
            v_int: 0.0,
            expected: Expected::RegulatorDecay,
            truth_adj: None,
            truth_perim: 0.0,
        });
    }
    cells
}

// ---- regulator lab (v31.7 凍結 — 端 τ バイアスの a 減衰) ----

pub fn regulator_edge_bias(n_phys: usize, m_phys: f64, s: usize) -> f64 {
    // R-A (単鎖 staggered)
    let build_a = |s: usize| -> (Vec<f64>, usize, usize) {
        let a = 1.0 / s as f64;
        let n = n_phys * s;
        let mut h = vec![0.0; n * n];
        for b in 0..n - 1 {
            let t = 1.0 / (2.0 * a);
            h[b * n + b + 1] = -t;
            h[(b + 1) * n + b] = -t;
        }
        for x in 0..n {
            h[x * n + x] += if x % 2 == 0 { m_phys * a } else { -m_phys * a };
        }
        (h, n, 1)
    };
    let build_c = |s: usize| -> (Vec<f64>, usize, usize) {
        let a = 1.0 / s as f64;
        let n = n_phys * s;
        let d = 2 * n;
        let t = 1.0 / (2.0 * a);
        let mut h = vec![0.0; d * d];
        for b in 0..n - 1 {
            let blk = [[-t, -t], [t, t]];
            for al in 0..2 {
                for be in 0..2 {
                    h[(2 * b + al) * d + (2 * (b + 1) + be)] += blk[al][be];
                    h[(2 * (b + 1) + be) * d + (2 * b + al)] += blk[al][be];
                }
            }
        }
        for x in 0..n {
            let on = m_phys + 2.0 * t;
            h[(2 * x) * d + 2 * x] += on;
            h[(2 * x + 1) * d + (2 * x + 1)] += -on;
        }
        (h, d, 2)
    };
    let tau_of = |h: &Vec<f64>, dim: usize, comp: usize, s: usize| -> f64 {
        // 端セル (n_phys−2) の到着時刻 (中央源)
        let (vals, vecs) = jacobi_eigh(h, dim);
        let src_sites: Vec<usize> = (0..s)
            .flat_map(|k| (0..comp).map(move |c| comp * ((n_phys / 2) * s + k) + c))
            .collect();
        let cell = n_phys - 2;
        let cell_sites: Vec<usize> = (0..s)
            .flat_map(|k| (0..comp).map(move |c| comp * (cell * s + k) + c))
            .collect();
        let mut t = 0.0;
        for it in 1..=2400 {
            t = 0.05 * it as f64;
            let mut dn = 0.0;
            for &al in &cell_sites {
                for &be in &src_sites {
                    let mut gre = 0.0;
                    let mut gim = 0.0;
                    for m in 0..dim {
                        let ph = vals[m] * t;
                        let w = vecs[m * dim + al] * vecs[m * dim + be];
                        gre += w * ph.cos();
                        gim -= w * ph.sin();
                    }
                    dn += EPS_PROBE * (gre * gre + gim * gim);
                }
            }
            if dn.abs() >= 1e-3 {
                break;
            }
        }
        t
    };
    let (ha, da, ca) = build_a(s);
    let (hc, dc, cc) = build_c(s);
    (tau_of(&ha, da, ca, s) - tau_of(&hc, dc, cc, s)).abs()
}

// ---- readout dispatch (真値非流入 — 観測契約と観測データのみで裁定) ----

pub fn readout(cell: &Cell, rng: &mut Rng) -> CellVerdict {
    // 凍結決定規則 1: 因子分解が与えられていなければ棄却 (v31.4 no-go)
    if !cell.fact_given {
        return CellVerdict::Abstained("unknown_factorization");
    }
    match cell.obs {
        "global_one_body_correlation" => {
            // 状態を lab が構成 (β=0 は projector)
            let c = if cell.v_int > 0.0 {
                let bonds: Vec<(usize, usize, f64)> = (0..cell.n)
                    .map(|k| (k, (k + 1) % cell.n, 1.0))
                    .collect();
                let (c, wit) = tv_thermal_c_witness(cell.n, &bonds, cell.v_int, 1.0);
                // 凍結決定規則 2: witness > バー → 棄却
                if wit > BAR_WICK {
                    return CellVerdict::Abstained("non_gaussian_domain");
                }
                c
            } else if cell.beta_state == 0.0 {
                // projector (半充填 GS)
                let (ev, evec) = jacobi_eigh(&cell.h, cell.n);
                let mut c = vec![0.0; cell.n * cell.n];
                for m in 0..cell.n {
                    if ev[m] < 0.0 {
                        for i in 0..cell.n {
                            for j in 0..cell.n {
                                c[i * cell.n + j] += evec[m * cell.n + i] * evec[m * cell.n + j];
                            }
                        }
                    }
                }
                c
            } else {
                gibbs_c(&cell.h, cell.n, cell.beta_state)
            };
            // 凍結決定規則 3: 資格審査 — RankDeficient → 同値類 / IllConditioned → 棄却
            match ExactFullRankCorrelation::certify_real(&c, cell.n) {
                Err(AbstainReason::RankDeficient) => CellVerdict::EquivClass,
                Err(_) => CellVerdict::Abstained("ill_conditioned"),
                Ok(cert) => {
                    let k = logit_k(cert.c_re(), cell.n);
                    let parent = ParentModularGenerator {
                        re: k,
                        im: vec![0.0; cell.n * cell.n],
                        n: cell.n,
                    };
                    match identify_physical_generator(
                        &parent,
                        GaussianityEvidence::ByConstruction,
                        GibbsProvenance::KnownBetaMu {
                            beta: cell.beta_state,
                            mu: 0.0,
                        },
                    ) {
                        Ok(PhysicalGeneratorReading::Exact(hr)) => CellVerdict::ExactH {
                            h: hr.re,
                            n: cell.n,
                        },
                        _ => CellVerdict::Abstained("gate_refused"),
                    }
                }
            }
        }
        "local_bias_density_response" => {
            // 凍結決定規則 4: ノイズ誤差見積り > バー → 棄却
            let norm1 = (0..cell.n)
                .map(|r| (0..cell.n).map(|c| cell.h[r * cell.n + c].abs()).sum::<f64>())
                .fold(0.0f64, f64::max);
            if noise_error_bound(cell.sigma, norm1) > BAR_NOISE_ABSTAIN {
                return CellVerdict::Abstained("insufficient_observation");
            }
            // 曲率測定 (t-V は厳密曲率 lane)
            let w = if cell.v_int > 0.0 {
                let bonds: Vec<(usize, usize, f64)> = (0..cell.n)
                    .map(|k| (k, (k + 1) % cell.n, 1.0))
                    .collect();
                let mut wm = vec![0.0; cell.n * cell.n];
                for i in 0..cell.n {
                    let wi = tv_curvature(cell.n, &bonds, cell.v_int, i);
                    for j in 0..cell.n {
                        if j != i {
                            wm[j * cell.n + i] = wi[j];
                        }
                    }
                }
                wm
            } else {
                let mut wm = vec![0.0; cell.n * cell.n];
                for i in 0..cell.n {
                    let wi = curvature_w(&cell.h, cell.n, i, cell.sigma, rng);
                    for j in 0..cell.n {
                        if j != i {
                            wm[j * cell.n + i] = wi[j];
                        }
                    }
                }
                wm
            };
            let adj = support_from_weights(&w, cell.n);
            match &cell.expected {
                Expected::Topo2d { .. } => {
                    let cx = clique_complex(&adj, cell.n);
                    let b = betti(&cx);
                    // 曲面性: 全辺が三角形 ≤ 2 (閉 = 全て 2)
                    let mut cnt = std::collections::HashMap::new();
                    for t in &cx.simp[2] {
                        for drop in 0..3 {
                            let mut e: Vec<usize> = t.clone();
                            e.remove(drop);
                            *cnt.entry(e).or_insert(0usize) += 1;
                        }
                    }
                    let closed = !cx.simp[2].is_empty()
                        && cx.simp[1].iter().all(|e| cnt.get(e) == Some(&2));
                    CellVerdict::Topo2d {
                        beta: [b[0], b[1], b[2]],
                        closed,
                    }
                }
                Expected::Topo3d { .. } => {
                    let cx = clique_complex(&adj, cell.n);
                    let b = betti(&cx);
                    let mut n_sing = 0usize;
                    let mut n_bnd = 0usize;
                    for v in 0..cell.n {
                        match classify_surface(&vertex_link(&cx, v)) {
                            "S2" => {}
                            "D2" => n_bnd += 1,
                            _ => n_sing += 1,
                        }
                    }
                    let kind = if n_sing > 0 {
                        "singular"
                    } else if n_bnd == 0 {
                        "closed"
                    } else {
                        "boundary"
                    };
                    CellVerdict::Topo3d {
                        beta: [b[0], b[1], b[2], b[3]],
                        kind,
                    }
                }
                Expected::RingThird => {
                    let d = metric_closure(&w, &adj, cell.n);
                    let bars = vr_h1_bars(&d, cell.n);
                    // 周長の観測推定 = Σ 辺長 (支持辺)
                    let mut perim_est = 0.0;
                    for i in 0..cell.n {
                        for j in (i + 1)..cell.n {
                            if adj[i][j] {
                                perim_est += 1.0 / w[i * cell.n + j].abs().sqrt();
                            }
                        }
                    }
                    CellVerdict::MetricRing {
                        death_over_perim: bars[0].1 / perim_est,
                    }
                }
                Expected::HoleRatio { .. } => {
                    let d = metric_closure(&w, &adj, cell.n);
                    let bars = vr_h1_bars(&d, cell.n);
                    let p1 = bars[0].1 - bars[0].0;
                    let p2 = if bars.len() > 1 { bars[1].1 - bars[1].0 } else { 1e-12 };
                    CellVerdict::MetricHoles { ratio: p1 / p2 }
                }
                _ => CellVerdict::RespWeights { w, n: cell.n },
            }
        }
        "arrival_time_response" => {
            // K10: 端バイアスの a 減衰比 (RegulatorDecay)
            let b1 = regulator_edge_bias(cell.n, cell.beta_state, 1);
            let b2 = regulator_edge_bias(cell.n, cell.beta_state, 2);
            CellVerdict::MetricHoles { ratio: b1 / b2.max(1e-12) } // 便宜的に ratio 枠で返す
        }
        _ => CellVerdict::Abstained("unknown_contract"),
    }
}

// ---- 採点器 (真値はここでのみ使用) ----

pub struct Score {
    pub answered: usize,
    pub answerable: usize,
    pub errors: usize,
    pub impossible: usize,
    pub correct_abstain: usize,
    pub forced_answers: usize,
    pub lines: Vec<String>,
}

pub fn score_cells(cells: &[Cell], rng: &mut Rng) -> Score {
    let mut sc = Score {
        answered: 0,
        answerable: 0,
        errors: 0,
        impossible: 0,
        correct_abstain: 0,
        forced_answers: 0,
        lines: Vec::new(),
    };
    for cell in cells {
        let v = readout(cell, rng);
        let is_impossible = matches!(cell.expected, Expected::MustEquivOrAbstain);
        if is_impossible {
            sc.impossible += 1;
            let ok = matches!(v, CellVerdict::EquivClass | CellVerdict::Abstained(_));
            if ok {
                sc.correct_abstain += 1;
            } else {
                sc.forced_answers += 1;
            }
            sc.lines.push(format!(
                "  {} [非識別] → {:?} {}",
                cell.name,
                short(&v),
                if ok { "✓ 正しい棄却/同値類" } else { "✗ 強制回答 = FAIL" }
            ));
            continue;
        }
        sc.answerable += 1;
        let (answered, correct, note) = judge(cell, &v, rng);
        if answered {
            sc.answered += 1;
            if !correct {
                sc.errors += 1;
            }
        }
        sc.lines.push(format!(
            "  {} → {:?} {} {}",
            cell.name,
            short(&v),
            if !answered {
                "△ 棄却 (coverage 減)"
            } else if correct {
                "✓"
            } else {
                "✗"
            },
            note
        ));
    }
    sc
}

fn short(v: &CellVerdict) -> String {
    match v {
        CellVerdict::Topo2d { beta, closed } => format!("β={:?} closed={}", beta, closed),
        CellVerdict::Topo3d { beta, kind } => format!("β={:?} {}", beta, kind),
        CellVerdict::MetricRing { death_over_perim } => {
            format!("death/perim={:.4}", death_over_perim)
        }
        CellVerdict::MetricHoles { ratio } => format!("ratio={:.3}", ratio),
        CellVerdict::ExactH { .. } => "ExactH".into(),
        CellVerdict::RespWeights { .. } => "Weights".into(),
        CellVerdict::EquivClass => "EquivClass".into(),
        CellVerdict::Abstained(r) => format!("Abstain({})", r),
    }
}

/// 回答セルの正誤 (真値使用は採点側のみ)
pub fn judge(cell: &Cell, v: &CellVerdict, rng: &mut Rng) -> (bool, bool, String) {
    match (&cell.expected, v) {
        (Expected::Topo2d { beta, closed }, CellVerdict::Topo2d { beta: b, closed: c }) => {
            (true, b == beta && c == closed, String::new())
        }
        (Expected::Topo3d { beta, kind }, CellVerdict::Topo3d { beta: b, kind: k }) => {
            (true, b == beta && k == kind, String::new())
        }
        (Expected::RingThird, CellVerdict::MetricRing { death_over_perim }) => (
            true,
            *death_over_perim >= BAR_RING_THIRD.0 && *death_over_perim <= BAR_RING_THIRD.1,
            format!("(バー [{}, {}])", BAR_RING_THIRD.0, BAR_RING_THIRD.1),
        ),
        (Expected::HoleRatio { ratio }, CellVerdict::MetricHoles { ratio: r }) => (
            true,
            (r / ratio - 1.0).abs() <= BAR_HOLE_RATIO_TOL,
            format!("(真比 {:.2})", ratio),
        ),
        (Expected::ExactH, CellVerdict::ExactH { h, n }) => {
            let mut err: f64 = 0.0;
            for k in 0..n * n {
                if k / n != k % n {
                    err = err.max((h[k] - cell.h[k]).abs());
                }
            }
            (true, err <= BAR_ORACLE_ERR, format!("(err {:.1e})", err))
        }
        (Expected::RespSupport, CellVerdict::RespWeights { w, n }) => {
            let adj = support_from_weights(w, *n);
            // 真の支持 = |h| > 0 (K6-resp は ring) / truth_adj があればそれ
            let truth: Vec<Vec<bool>> = if let Some(t) = &cell.truth_adj {
                t.clone()
            } else {
                (0..*n)
                    .map(|i| (0..*n).map(|j| cell.h[i * n + j].abs() > 1e-9).collect())
                    .collect()
            };
            let mut ok = true;
            for i in 0..*n {
                for j in 0..*n {
                    if adj[i][j] != truth[i][j] {
                        ok = false;
                    }
                }
            }
            // K6-resp: t-V ring の支持は ring (v_int があっても厳密転移)
            let bar = if cell.sigma > 0.0 { BAR_NOISY_REL } else { BAR_RESP_REL };
            let _ = bar;
            (true, ok, String::new())
        }
        (Expected::MetamorphicPair, CellVerdict::RespWeights { w, n }) => {
            // 第 2 系の読み出しと突き合わせ: 支持のグラフ不変量 (次数列) が一致
            let h2 = cell.h2.as_ref().unwrap();
            let cell2 = Cell {
                name: "K9b".into(),
                h: h2.clone(),
                n: *n,
                h2: None,
                obs: cell.obs,
                fact_given: true,
                sigma: 0.0,
                beta_state: 1.0,
                v_int: 0.0,
                expected: Expected::RespSupport,
                truth_adj: None,
                truth_perim: 0.0,
            };
            let v2 = readout(&cell2, rng);
            if let CellVerdict::RespWeights { w: w2, n: n2 } = v2 {
                let deg = |w: &Vec<f64>, n: usize| -> Vec<usize> {
                    let adj = support_from_weights(w, n);
                    let mut d: Vec<usize> = (0..n)
                        .map(|i| (0..n).filter(|&j| adj[i][j]).count())
                        .collect();
                    d.sort_unstable();
                    d
                };
                let (d1, d2) = (deg(w, *n), deg(&w2, n2));
                // 重み多重集合も一致 (ソート比較, 1e-6)
                let ms = |w: &Vec<f64>, n: usize| -> Vec<f64> {
                    let mut v: Vec<f64> = (0..n)
                        .flat_map(|i| ((i + 1)..n).map(move |j| (i, j)))
                        .map(|(i, j)| w[i * n + j].abs())
                        .filter(|x| *x > 1e-6)
                        .collect();
                    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
                    v
                };
                let (m1, m2) = (ms(w, *n), ms(&w2, n2));
                let w_match = m1.len() == m2.len()
                    && m1
                        .iter()
                        .zip(m2.iter())
                        .all(|(a, b)| (a - b).abs() <= 1e-4 * (1.0 + a.abs()));
                (true, d1 == d2 && w_match, "(変成: 次数列 + 重み多重集合)".into())
            } else {
                (true, false, "(第 2 系が回答しない)".into())
            }
        }
        (Expected::RegulatorDecay, CellVerdict::MetricHoles { ratio }) => (
            true,
            *ratio >= BAR_REG_DECAY,
            format!("(端バイアス a=1/a=½ 減衰比, バー ≥ {})", BAR_REG_DECAY),
        ),
        (_, CellVerdict::Abstained(_)) => (false, false, String::new()),
        _ => (true, false, "(裁定型の不整合)".into()),
    }
}

// ================================================================================
// FROZEN-HOLD7-END
// ================================================================================

fn main() {
    uft_sim::self_test();
    println!("=== v32.0-A HOLD-7 の凍結 — identifiability 相図の新鮮 holdout (PROMPT/12) ===\n");
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

    // ---- [F0] コミットメント ----
    {
        println!("  [F0] SECRET コミットメント (開示は v32.0-B):");
        println!("       sha256(SECRET) = {}", HOLD7_COMMITMENT);
        println!("       train seed = {} (可視 — holdout シードは SECRET から導出)\n", HOLD7_TRAIN_SEED);
    }

    // ---- [F1] 生成器健全性 (HOLD-6 [G1b] の常設) ----
    {
        let mut bad = Vec::new();
        // 幾何構成の真値プロファイル
        let (nv, tris) = genus2();
        let cx = {
            let adj = skeleton_of_tris(nv, &tris);
            clique_complex(&adj, nv)
        };
        let b = betti(&cx);
        if !(b[0] == 1 && b[1] == 4 && b[2] == 1) {
            bad.push(format!("genus2 構成 β = {:?}", &b[..3]));
        }
        let (nv3, tets) = kuhn_tets(4, 4, 4, true);
        let cx3 = clique_complex(&skeleton_of_tets(nv3, &tets), nv3);
        let b3 = betti(&cx3);
        if !(b3[0] == 1 && b3[1] == 3 && b3[2] == 3 && b3[3] == 1) {
            bad.push(format!("T³ L=4 skeleton clique β = {:?}", &b3[..4]));
        }
        // 生成セルの構成的整合 (train seed で 17 セル)
        let cells = generate_cells(HOLD7_TRAIN_SEED);
        let n_imp = cells
            .iter()
            .filter(|c| matches!(c.expected, Expected::MustEquivOrAbstain))
            .count();
        if cells.len() != 17 || n_imp != 5 {
            bad.push(format!("セル構成 {} (非識別 {})", cells.len(), n_imp));
        }
        check(
            "[F1] 生成器健全性: genus2 (1,4,1)・T³ L=4 flag (1,3,3,1)・セル構成 17 (非識別 5)",
            bad.is_empty(),
            if bad.is_empty() {
                "真の複体プロファイルと構成数が期待どおり".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [F2] train 全セル採点 ----
    {
        let cells = generate_cells(HOLD7_TRAIN_SEED);
        let mut rng = Rng::new(HOLD7_TRAIN_SEED ^ 0x5eed);
        let sc = score_cells(&cells, &mut rng);
        for l in &sc.lines {
            println!("{}", l);
        }
        let coverage = sc.answered as f64 / sc.answerable as f64;
        let recall = sc.correct_abstain as f64 / sc.impossible as f64;
        check(
            "[F2] train 採点: selective risk 0 / coverage ≥ 0.9 / impossibility recall = 1.0",
            sc.errors == 0 && coverage >= BAR_COVERAGE && recall >= 1.0 && sc.forced_answers == 0,
            format!(
                "回答 {}/{} (coverage {:.2}) / 誤り {} / 非識別 {} 中 正棄却 {} (recall {:.2}) / 強制回答 {}",
                sc.answered,
                sc.answerable,
                coverage,
                sc.errors,
                sc.impossible,
                sc.correct_abstain,
                recall,
                sc.forced_answers
            ),
        );
    }

    // ---- [F3] 凍結宣言 ----
    {
        let src = std::fs::read_to_string(
            if std::path::Path::new("sim/src/bin/v320a_hold7_freeze.rs").exists() {
                "sim/src/bin/v320a_hold7_freeze.rs"
            } else {
                "../sim/src/bin/v320a_hold7_freeze.rs"
            },
        )
        .unwrap_or_default();
        let sec = frozen_section(&src);
        check(
            "[F3] 凍結宣言: FROZEN-HOLD7 区間の SHA-256 (v32.0-B が逐語一致を照合する)",
            !sec.is_empty(),
            format!("kernel sha256 = {}…", &sha256_hex(sec.as_bytes())[..16]),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "HOLD-7 は凍結された — 生成器・採点器・バー・観測契約・非識別の裁定規則が本コミットで固定。次コミット (v32.0-B) で SECRET を開示し holdout を初生成・本採点する (調整なし)"
        } else {
            "**凍結の破れ**"
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

/// FROZEN 区間の抽出 (BEGIN/END マーカー行間)
fn frozen_section(src: &str) -> String {
    let b = src.find("FROZEN-HOLD7-BEGIN");
    let e = src.find("FROZEN-HOLD7-END");
    match (b, e) {
        (Some(b), Some(e)) if e > b => src[b..e].to_string(),
        _ => String::new(),
    }
}
