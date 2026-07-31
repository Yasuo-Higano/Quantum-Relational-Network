//! v31.6 計量 Vietoris–Rips persistence / 高 genus / 3D homology (PROMPT/12 第三十一期)
//!
//! **凍結済み核からのみ開始する** (v31.2 曲率則 = 状態非依存の重み ŵ = |t_ij|²)。
//! raw affinity を計量と呼ばない — 事前登録する計量は
//!   **辺長変換 ℓ_ij = 1/√ŵ_ij + 最短路閉包 (Dijkstra)**
//! であり、閉包が三角不等式を構成的に保証する (擬計量証明書 = 対称・非負・d(x,x)=0)。
//!
//! 主要部:
//!   A. **Z2 homology engine の 4-simplex 拡張**: β₃ = dim ker ∂₃ − rank ∂₄ —
//!      **im ∂₄ を必須とする** (K5 anchor: ∂₄ を落とすと β₃ = 1 に化ける)。
//!      anchors: K5 (1,0,0,0) / ∂Δ⁴ (1,0,0,1) / 16-cell (1,0,0,1) /
//!      Kuhn 三角化 T³ L=3 **(1,3,3,1)** / 中実 Kuhn 立方体 (1,0,0,0)
//!   B. **3-manifold link 分類**: 内部頂点 link = S² (閉曲面 β=(1,0,1))・境界頂点
//!      link = D² (円板)・特異 link は abstain/not-manifold (2 tet の 1 点接着で検査)
//!   C. **from-state 3D end-to-end**: 熱的 Gaussian 状態 + 曲率則測定 lane の重みから
//!      T³ (1,3,3,1)・S³ = 16-cell (1,0,0,1)・3-ball (1,0,0,0) を同定 (β₀..β₃ +
//!      link 多様体性 — 数え上げだけでは不十分)
//!   D. **計量 VR persistence**: full filtration + barcode (birth/death/persistence を
//!      採点 — 穴の個数だけではない):
//!      円環の 1/3 法則 (VR の H1 bar は死亡 ≈ 周長/3)・大小 2 穴の寿命分離
//!      (v29.6 の既知限界の解消)・genus-2 (β₁ = 4, 状態から)・narrow neck /
//!      close holes の敵対対照・**ノイズ安定性 bottleneck(D, D') ≤ sup|Δd|**
//!      (安定性定理の機械照合)
//!
//! 実行: cargo run --release --bin v316_metric_vr

use uft_sim::{jacobi_eigh, Rng, C64};

// ================================================================ 単体複体と Z2 rank

/// 次元別単体リスト (頂点 id 昇順タプル)
#[derive(Clone, Default)]
struct Complex {
    /// simp[d] = d-単体の頂点列 (長さ d+1, 昇順) の列
    simp: Vec<Vec<Vec<usize>>>,
}

impl Complex {
    /// グラフの clique complex (次元 ≤ 4 = 5-clique まで列挙)
    fn clique_complex(adj: &[Vec<bool>], n: usize) -> Complex {
        let mut c = Complex::default();
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
        // k-clique の拡張列挙 (最後の頂点より大きい共通隣接のみ追加 — 重複なし)
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

    fn betti(&self) -> Vec<i64> {
        // β_k = dim ker ∂_k − rank ∂_{k+1} = (n_k − rank ∂_k) − rank ∂_{k+1}
        let maxd = self.simp.len() - 1;
        let mut rank = vec![0usize; maxd + 2]; // rank[k] = rank ∂_k (∂_0 = 0)
        for k in 1..=maxd {
            rank[k] = boundary_rank(&self.simp[k - 1], &self.simp[k]);
        }
        let mut b = Vec::new();
        for k in 0..=maxd {
            let nk = self.simp[k].len() as i64;
            let r_k = rank[k] as i64;
            let r_k1 = if k + 1 <= maxd { rank[k + 1] as i64 } else { 0 };
            b.push(nk - r_k - r_k1);
        }
        b
    }
}

/// ∂_k の Z2 rank (rows = (k−1) 単体, cols = k 単体) — 列簡約
fn boundary_rank(faces: &[Vec<usize>], simps: &[Vec<usize>]) -> usize {
    use std::collections::HashMap;
    let mut fidx: HashMap<&[usize], usize> = HashMap::new();
    for (i, f) in faces.iter().enumerate() {
        fidx.insert(f.as_slice(), i);
    }
    // 列 = 境界面の行 id (昇順)
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
    // pivot 簡約
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

fn xor_merge(a: &[usize], b: &[usize]) -> Vec<usize> {
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

// ================================================================ link 分類 (3-complex)

/// 頂点 v の link (v を含む単体から v を除いた面の複体 — 2 次元まで)
fn vertex_link(c: &Complex, v: usize) -> Complex {
    let mut l = Complex::default();
    l.simp = vec![Vec::new(); 3];
    let mut push_unique = |d: usize, s: Vec<usize>, l: &mut Complex| {
        if !l.simp[d].contains(&s) {
            l.simp[d].push(s);
        }
    };
    for d in 1..c.simp.len() {
        for s in &c.simp[d] {
            if let Some(pos) = s.iter().position(|&u| u == v) {
                let mut f: Vec<usize> = s.clone();
                f.remove(pos);
                if f.len() <= 3 {
                    push_unique(f.len() - 1, f, &mut l);
                }
            }
        }
    }
    l
}

/// 2-complex の曲面分類: "S2" (閉球面) / "D2" (円板) / "other"
fn classify_surface(l: &Complex) -> &'static str {
    if l.simp.len() < 3 || l.simp[2].is_empty() {
        return "other";
    }
    // 各辺の三角形帰属数
    let mut cnt = std::collections::HashMap::new();
    for t in &l.simp[2] {
        for drop in 0..3 {
            let mut e: Vec<usize> = t.clone();
            e.remove(drop);
            *cnt.entry(e).or_insert(0usize) += 1;
        }
    }
    let mut n1 = 0usize;
    let mut bad = false;
    for e in &l.simp[1] {
        match cnt.get(e) {
            Some(&2) => {}
            Some(&1) => n1 += 1,
            _ => bad = true,
        }
    }
    if bad {
        return "other";
    }
    let b = l.betti();
    if n1 == 0 && b[0] == 1 && b[1] == 0 && b[2] == 1 {
        "S2"
    } else if n1 > 0 && b[0] == 1 && b[1] == 0 && b[2] == 0 {
        "D2"
    } else {
        "other"
    }
}

// ================================================================ 3D anchors (直接構成)

/// Kuhn 三角化: 各立方体 (x,y,z) を主対角に沿う 6 tet に割る。pbc = 周期境界
fn kuhn_tets(lx: usize, ly: usize, lz: usize, pbc: bool) -> (usize, Vec<Vec<usize>>) {
    let nvx = if pbc { lx } else { lx + 1 };
    let nvy = if pbc { ly } else { ly + 1 };
    let nvz = if pbc { lz } else { lz + 1 };
    let vid = |x: usize, y: usize, z: usize| -> usize {
        ((x % nvx) * nvy + (y % nvy)) * nvz + (z % nvz)
    };
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
                    // 単調路 0 → e_{p0} → e_{p0}+e_{p1} → (1,1,1)
                    let mut pos = [x, y, z];
                    let mut t = vec![vid(pos[0], pos[1], pos[2])];
                    for &ax in p {
                        pos[ax] += 1;
                        t.push(vid(pos[0], pos[1], pos[2]));
                    }
                    let mut t2 = t.clone();
                    t2.sort_unstable();
                    t2.dedup();
                    if t2.len() == 4 {
                        tets.push(t2);
                    }
                }
            }
        }
    }
    tets.sort();
    tets.dedup();
    (nvx * nvy * nvz, tets)
}

/// tet 集合から複体 (面閉包)
fn complex_from_tets(nv: usize, tets: &[Vec<usize>]) -> Complex {
    let mut c = Complex::default();
    c.simp = vec![Vec::new(); 4];
    for v in 0..nv {
        c.simp[0].push(vec![v]);
    }
    let mut edges = std::collections::BTreeSet::new();
    let mut tris = std::collections::BTreeSet::new();
    for t in tets {
        for a in 0..4 {
            for b in (a + 1)..4 {
                edges.insert(vec![t[a], t[b]]);
                for d in (b + 1)..4 {
                    tris.insert(vec![t[a], t[b], t[d]]);
                }
            }
        }
    }
    c.simp[1] = edges.into_iter().collect();
    c.simp[2] = tris.into_iter().collect();
    c.simp[3] = tets.to_vec();
    c
}

/// 1-skeleton 隣接
fn skeleton_adj(c: &Complex, n: usize) -> Vec<Vec<bool>> {
    let mut adj = vec![vec![false; n]; n];
    for e in &c.simp[1] {
        adj[e[0]][e[1]] = true;
        adj[e[1]][e[0]] = true;
    }
    adj
}

// ================================================================ 曲率則測定 lane (v31.2 — 実対称)

fn evolve_diag_c0(vals: &[f64], vecs: &[f64], c0: &[f64], n: usize, t: f64) -> Vec<f64> {
    let mut ct = vec![C64::new(0.0, 0.0); n * n];
    for a in 0..n {
        for b in 0..n {
            let mut s = 0.0;
            for i in 0..n {
                s += vecs[a * n + i] * c0[i] * vecs[b * n + i];
            }
            ct[a * n + b] = C64::expi(-(vals[a] - vals[b]) * t).scale(s);
        }
    }
    let mut diag = vec![0.0; n];
    for i in 0..n {
        let mut s = C64::new(0.0, 0.0);
        for a in 0..n {
            for b in 0..n {
                s = s + ct[a * n + b].scale(vecs[a * n + i] * vecs[b * n + i]);
            }
        }
        diag[i] = s.re;
    }
    diag
}

/// 密度曲率測定 lane: source i → 各 j の ŵ = |h_ji|² (Richardson, 時系列のみ)。
/// dt はスペクトル半径でスケールする (dt_eff = dt_base/‖h‖₁) — 高次数グラフ
/// (T³ 骨格は次数 14) では固定 dt の高階微分誤差が O(1) になる (1D 調整値の
/// 暗黙拡張は本期の禁止事項そのもの — 適応則を器械契約にする)
fn curvature_w(h: &[f64], n: usize, i: usize, eps: f64, dt_base: f64) -> Vec<f64> {
    let norm1 = (0..n)
        .map(|r| (0..n).map(|c| h[r * n + c].abs()).sum::<f64>())
        .fold(0.0f64, f64::max)
        .max(1.0);
    let dt = dt_base / norm1;
    let (vals, vecs) = jacobi_eigh(h, n);
    let times = [-dt, -dt / 2.0, dt / 2.0, dt];
    let mut nplus = [vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    let mut nminus = [vec![0.0; n], vec![0.0; n], vec![0.0; n], vec![0.0; n]];
    let mut n0 = [vec![0.0; n], vec![0.0; n]];
    for (pi, sign) in [(0usize, 1.0), (1usize, -1.0)] {
        let c0: Vec<f64> = (0..n)
            .map(|s| if s == i { 0.5 + sign * eps } else { 0.5 })
            .collect();
        n0[pi] = c0.clone();
        for (ti, &t) in times.iter().enumerate() {
            let d = evolve_diag_c0(&vals, &vecs, &c0, n, t);
            if pi == 0 {
                nplus[ti] = d;
            } else {
                nminus[ti] = d;
            }
        }
    }
    let mut w = vec![0.0; n];
    for j in 0..n {
        let d2 = |arr: &[Vec<f64>; 4], base: &Vec<f64>, half: bool| -> f64 {
            let (tm, tp, dd) = if half { (1, 2, dt / 2.0) } else { (0, 3, dt) };
            (arr[tp][j] - 2.0 * base[j] + arr[tm][j]) / (dd * dd)
        };
        let coarse = (d2(&nplus, &n0[0], false) - d2(&nminus, &n0[1], false)) / (4.0 * eps);
        let fine = (d2(&nplus, &n0[0], true) - d2(&nminus, &n0[1], true)) / (4.0 * eps);
        w[j] = (4.0 * fine - coarse) / 3.0;
    }
    w
}

/// 測定重み → 支持 (スケールガード付き gap 則, v31.3 と同一)
fn support_from_weights(w: &[f64], n: usize) -> Vec<Vec<bool>> {
    let mut vals: Vec<f64> = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            vals.push(w[i * n + j].abs().max(1e-300));
        }
    }
    let mut sorted = vals.clone();
    sorted.sort_by(|a, b| b.partial_cmp(a).unwrap());
    // 最終規則 (v31.6 で凍結): スケールガード窓 (max·1e-3) の**内側**の最大対数段差が
    // 有意 (≥ ln 3) ならそこで切る。有意な窓内段差がなければ窓内は単一クラスタ =
    // 全て辺として窓境界で切る。旧規則 2 種の故障を両方閉じる:
    //   (i) ガードなし: f64 尾部の発散段差を拾う (v31.3 で訂正済み)
    //   (ii) 跨ぎ段差を本命にする: 多段階系で物理段差を尾部跨ぎが上書き
    let guard = sorted[0] * 1e-3;
    let mut cut: Option<usize> = None;
    let mut best_gap = 0.0;
    for k in 0..sorted.len() - 1 {
        if sorted[k + 1] < guard {
            break; // 窓内段差のみ (両端 ≥ guard)
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
    // 対称化 (and)
    for i in 0..n {
        for j in 0..n {
            let s = adj[i][j] && adj[j][i];
            adj[i][j] = s;
        }
    }
    adj
}

// ================================================================ 計量と VR persistence

/// 事前登録計量: ℓ_ij = 1/√ŵ_ij (辺長変換) → 最短路閉包 (Dijkstra 全対)
fn metric_closure(w: &[f64], adj: &[Vec<bool>], n: usize) -> Vec<f64> {
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

/// VR persistence (H0/H1) — 標準列簡約。返り値: H1 の (birth, death) bars (有限のみ)
fn vr_h1_bars(d: &[f64], n: usize) -> Vec<(f64, f64)> {
    // 単体: 頂点 (f=0), 辺 (f=d_ij), 三角形 (f=max 対距離)
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
    // 辺 id 引き
    let mut eid = std::collections::HashMap::new();
    for (id, s) in simps.iter().enumerate() {
        if s.dim == 1 {
            eid.insert((s.verts[0], s.verts[1]), pos[id]);
        }
    }
    // 境界列 (整列順の行 id)
    let mut cols: Vec<Vec<usize>> = vec![Vec::new(); simps.len()];
    for (id, s) in simps.iter().enumerate() {
        let j = pos[id];
        match s.dim {
            1 => {
                // 頂点行 = 頂点 id の整列位置 (頂点は f=0 で先頭 n 個 — pos で引く)
                let mut c = vec![pos[s.verts[0]], pos[s.verts[1]]];
                c.sort_unstable();
                cols[j] = c;
            }
            2 => {
                let (a, b, c3) = (s.verts[0], s.verts[1], s.verts[2]);
                let mut c = vec![
                    eid[&(a, b)],
                    eid[&(a, c3)],
                    eid[&(b, c3)],
                ];
                c.sort_unstable();
                cols[j] = c;
            }
            _ => {}
        }
    }
    // 簡約
    let m = simps.len();
    let mut pivot_owner: Vec<Option<usize>> = vec![None; m];
    let mut reduced: Vec<Vec<usize>> = vec![Vec::new(); m];
    let mut pairs: Vec<(usize, usize)> = Vec::new(); // (birth_pos, death_pos)
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
    // H1 bars: birth = 辺 (dim1) が正 (どの対でも死んでいない) …対 (low=辺, killer=三角形)
    let mut bars = Vec::new();
    for &(b, dth) in &pairs {
        let sb = &simps[order[b]];
        let sd = &simps[order[dth]];
        if sb.dim == 1 && sd.dim == 2 {
            let (fb, fd) = (sb.f, sd.f);
            if fd > fb + 1e-12 {
                bars.push((fb, fd));
            }
        }
    }
    bars.sort_by(|a, b| (b.1 - b.0).partial_cmp(&(a.1 - a.0)).unwrap());
    bars
}

/// bottleneck 距離 (小さい図用の brute force — 有意 bar ≤ 6 に制限して呼ぶ)
fn bottleneck(a: &[(f64, f64)], b: &[(f64, f64)]) -> f64 {
    let ca: Vec<(f64, f64)> = a.to_vec();
    let cb: Vec<(f64, f64)> = b.to_vec();
    let na = ca.len();
    let nb = cb.len();
    assert!(na <= 6 && nb <= 6, "bottleneck brute force は小図専用");
    // 各 bar の対角距離 = persistence/2
    let diag = |p: &(f64, f64)| (p.1 - p.0) / 2.0;
    let dist = |p: &(f64, f64), q: &(f64, f64)| (p.0 - q.0).abs().max((p.1 - q.1).abs());
    // b 側の割当 (a の各要素 → b の要素 or 対角) を全探索
    let mut best = f64::INFINITY;
    let nb_opts = nb + 1;
    let total = (nb_opts as u64).pow(na as u32);
    for code in 0..total {
        let mut c = code;
        let mut used = vec![false; nb];
        let mut cost: f64 = 0.0;
        let mut ok = true;
        for i in 0..na {
            let choice = (c % nb_opts as u64) as usize;
            c /= nb_opts as u64;
            if choice == nb {
                cost = cost.max(diag(&ca[i]));
            } else {
                if used[choice] {
                    ok = false;
                    break;
                }
                used[choice] = true;
                cost = cost.max(dist(&ca[i], &cb[choice]));
            }
        }
        if !ok {
            continue;
        }
        for j in 0..nb {
            if !used[j] {
                cost = cost.max(diag(&cb[j]));
            }
        }
        best = best.min(cost);
    }
    best
}

// ================================================================ 2D メッシュ (genus-2 等)

/// 5×5 周期正方格子の Kuhn (右下がり対角) 三角化 T²
fn torus_mesh(l: usize) -> (usize, Vec<Vec<usize>>) {
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

/// 2-complex (三角形集合) → Complex
fn complex_from_tris(nv: usize, tris: &[Vec<usize>]) -> Complex {
    let mut c = Complex::default();
    c.simp = vec![Vec::new(); 3];
    for v in 0..nv {
        c.simp[0].push(vec![v]);
    }
    let mut edges = std::collections::BTreeSet::new();
    for t in tris {
        for a in 0..3 {
            for b in (a + 1)..3 {
                edges.insert(vec![t[a], t[b]]);
            }
        }
    }
    c.simp[1] = edges.into_iter().collect();
    c.simp[2] = tris.to_vec();
    c
}

/// genus-2: 2 つの T² (5×5) から頂点星を除き、境界六角形を同一視して接着
fn genus2_complex() -> (usize, Vec<Vec<usize>>) {
    let l = 5usize;
    let (nv, tris) = torus_mesh(l);
    // 除く頂点 = 0。その link 六角形を求める
    let star_of = |tris: &[Vec<usize>], v: usize| -> (Vec<Vec<usize>>, Vec<usize>) {
        let keep: Vec<Vec<usize>> = tris.iter().filter(|t| !t.contains(&v)).cloned().collect();
        // 境界頂点 = v と辺を張っていた頂点 (六角形)
        let mut nb = std::collections::BTreeSet::new();
        for t in tris.iter().filter(|t| t.contains(&v)) {
            for &u in t {
                if u != v {
                    nb.insert(u);
                }
            }
        }
        (keep, nb.into_iter().collect())
    };
    let (t1, b1) = star_of(&tris, 0);
    let (t2, b2) = star_of(&tris, 0);
    // 境界六角形を巡回順に並べる (keep 側の境界辺 = v の star にあった辺)
    let cycle_order = |tris_all: &[Vec<usize>], v: usize, nb: &[usize]| -> Vec<usize> {
        // v を含む三角形 {v,a,b} の (a,b) が境界辺
        let mut adj: std::collections::HashMap<usize, Vec<usize>> = Default::default();
        for t in tris_all.iter().filter(|t| t.contains(&v)) {
            let others: Vec<usize> = t.iter().cloned().filter(|&u| u != v).collect();
            adj.entry(others[0]).or_default().push(others[1]);
            adj.entry(others[1]).or_default().push(others[0]);
        }
        let mut order = vec![nb[0]];
        let mut prev = usize::MAX;
        while order.len() < nb.len() {
            let cur = *order.last().unwrap();
            let nexts = &adj[&cur];
            let nx = if nexts[0] != prev { nexts[0] } else { nexts[1] };
            prev = cur;
            order.push(nx);
        }
        order
    };
    let o1 = cycle_order(&tris, 0, &b1);
    let o2 = cycle_order(&tris, 0, &b2);
    // 頂点番号の付け替え: 複体 1 = そのまま (0 は除去済みだが番号は残す — 未使用)、
    // 複体 2 = +nv シフト。その後 o2[k] ↔ o1[k] を同一視 (逆向きで貼る)
    let mut map2 = |u: usize| -> usize { u + nv };
    let mut tris_all: Vec<Vec<usize>> = t1.clone();
    for t in &t2 {
        let mut s: Vec<usize> = t.iter().map(|&u| map2(u)).collect();
        s.sort_unstable();
        tris_all.push(s);
    }
    // 同一視: map2(o2[k]) → o1[rev k]
    let mut rename: std::collections::HashMap<usize, usize> = Default::default();
    let m = o1.len();
    for k in 0..m {
        rename.insert(o2[k] + nv, o1[(m - k) % m]);
    }
    let mut final_tris: Vec<Vec<usize>> = Vec::new();
    for t in &tris_all {
        let mut s: Vec<usize> = t
            .iter()
            .map(|&u| *rename.get(&u).unwrap_or(&u))
            .collect();
        s.sort_unstable();
        s.dedup();
        if s.len() == 3 {
            final_tris.push(s);
        }
    }
    final_tris.sort();
    final_tris.dedup();
    // 頂点を圧縮 (使われる id のみ)
    let mut used: Vec<usize> = final_tris.iter().flatten().cloned().collect();
    used.sort_unstable();
    used.dedup();
    let idx: std::collections::HashMap<usize, usize> =
        used.iter().enumerate().map(|(i, &u)| (u, i)).collect();
    let tris_c: Vec<Vec<usize>> = final_tris
        .iter()
        .map(|t| {
            let mut s: Vec<usize> = t.iter().map(|&u| idx[&u]).collect();
            s.sort_unstable();
            s
        })
        .collect();
    (used.len(), tris_c)
}

// ================================================================ main

fn main() {
    uft_sim::self_test();
    println!("=== v31.6 計量 VR persistence / 高 genus / 3D homology (PROMPT/12) ===");
    println!("(事前登録計量: ℓ = 1/√ŵ [曲率則の凍結重み] + 最短路閉包 — raw affinity を計量と呼ばない)\n");
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

    // ---- [V0] 3D homology engine anchors (∂₄ 必須) ----
    {
        // K5: clique complex = Δ⁴ (可縮) — ∂₄ を落とすと β₃ = 1 に化ける
        let n = 5;
        let adj = vec![vec![true; n]; n];
        let mut adj = adj;
        for i in 0..n {
            adj[i][i] = false;
        }
        let c = Complex::clique_complex(&adj, n);
        let b_k5 = c.betti();
        let n4 = c.simp[4].len();
        // ∂Δ⁴ = S³: K5 の 4-simplex を除いた複体
        let mut c_s3 = c.clone();
        c_s3.simp[4].clear();
        let b_s3 = c_s3.betti();
        // 16-cell (4 次元交差多面体の境界 = S³): 頂点 ±e_i, 反対対以外は全結線
        let n8 = 8;
        let mut adj8 = vec![vec![true; n8]; n8];
        for i in 0..n8 {
            adj8[i][i] = false;
        }
        for i in 0..4 {
            adj8[2 * i][2 * i + 1] = false;
            adj8[2 * i + 1][2 * i] = false;
        }
        let c16 = Complex::clique_complex(&adj8, n8);
        let b16 = c16.betti();
        // Kuhn T³ (L=3, pbc): tets 直接構成
        let (nv_t3, tets_t3) = kuhn_tets(3, 3, 3, true);
        let ct3 = complex_from_tets(nv_t3, &tets_t3);
        let bt3 = ct3.betti();
        // 中実 Kuhn 立方体 (L=2, open): 3-ball
        let (nv_bl, tets_bl) = kuhn_tets(2, 2, 2, false);
        let cbl = complex_from_tets(nv_bl, &tets_bl);
        let bbl = cbl.betti();
        let ok = b_k5[..4] == [1, 0, 0, 0]
            && n4 == 1
            && b_s3[..4] == [1, 0, 0, 1]
            && b16[..4] == [1, 0, 0, 1]
            && bt3[..4] == [1, 3, 3, 1]
            && bbl[..4] == [1, 0, 0, 0];
        check(
            "[V0] 3D engine anchors: K5 (1,0,0,0 — im ∂₄ 必須) / ∂Δ⁴ (1,0,0,1) / 16-cell (1,0,0,1) / Kuhn T³ (1,3,3,1) / 中実立方体 (1,0,0,0)",
            ok,
            format!(
                "K5 β={:?} (4-simplex {} 個 — 除くと β₃=1 = ∂Δ⁴ {:?}) / 16-cell {:?} / T³ {:?} ({} tets) / ball {:?}",
                &b_k5[..4],
                n4,
                &b_s3[..4],
                &b16[..4],
                &bt3[..4],
                tets_t3.len(),
                &bbl[..4]
            ),
        );
    }

    // ---- [V1] 3-manifold link 分類 ----
    {
        // T³: 全頂点の link = S² (閉 3 多様体)
        let (nv_t3, tets_t3) = kuhn_tets(3, 3, 3, true);
        let ct3 = complex_from_tets(nv_t3, &tets_t3);
        let mut all_s2 = true;
        for v in 0..nv_t3 {
            if classify_surface(&vertex_link(&ct3, v)) != "S2" {
                all_s2 = false;
            }
        }
        // 中実立方体 L=2 (27 頂点): 中心 (1,1,1) は内部 = S²・角は境界 = D²
        let (nv_bl, tets_bl) = kuhn_tets(2, 2, 2, false);
        let cbl = complex_from_tets(nv_bl, &tets_bl);
        let vid = |x: usize, y: usize, z: usize| (x * 3 + y) * 3 + z;
        let center = classify_surface(&vertex_link(&cbl, vid(1, 1, 1)));
        let corner = classify_surface(&vertex_link(&cbl, vid(0, 0, 0)));
        let mut n_int = 0usize;
        let mut n_bnd = 0usize;
        let mut n_sing = 0usize;
        for v in 0..nv_bl {
            match classify_surface(&vertex_link(&cbl, v)) {
                "S2" => n_int += 1,
                "D2" => n_bnd += 1,
                _ => n_sing += 1,
            }
        }
        // 特異対照: 2 つの tet が 1 頂点だけ共有 → link は非連結 → not-manifold
        let tets_sing = vec![vec![0, 1, 2, 3], vec![3, 4, 5, 6]];
        let csing = complex_from_tets(7, &tets_sing);
        let sing = classify_surface(&vertex_link(&csing, 3));
        check(
            "[V1] 3-manifold link 分類: T³ 全 27 頂点 link = S² / 中実立方体は内部 1 (S²)・境界 26 (D²)・特異 0 / 1 点接着 tet 対は not-manifold",
            all_s2 && center == "S2" && corner == "D2" && n_int == 1 && n_bnd == 26 && n_sing == 0
                && sing == "other",
            format!(
                "T³ link 全 S² = {} / ball: 内部 {} 境界 {} 特異 {} (中心 {} 角 {}) / 特異対照 = {}",
                all_s2, n_int, n_bnd, n_sing, center, corner, sing
            ),
        );
    }

    // ---- [V2] from-state 3D end-to-end (曲率則測定 lane → 支持 → β + link) ----
    {
        let eps = 0.3;
        let dt = 0.02;
        let mut results = Vec::new();
        let mut all_ok = true;
        // (a) T³ skeleton (27 頂点)
        // T³ from-state は L=4: L=3 では周期軸線 (0-1-2-0) が巻き付き 3-クリークを
        // 作り clique complex ≠ 三角化 (flag 性が破れる) — L≥4 で回復 (発見として記録)
        let cases: Vec<(&str, Complex, usize, [i64; 4], &str)> = {
            let (nv_t3, tets_t3) = kuhn_tets(4, 4, 4, true);
            let ct3 = complex_from_tets(nv_t3, &tets_t3);
            let (nv16, c16) = {
                let n8 = 8;
                let mut adj8 = vec![vec![true; n8]; n8];
                for i in 0..n8 {
                    adj8[i][i] = false;
                }
                for i in 0..4 {
                    adj8[2 * i][2 * i + 1] = false;
                    adj8[2 * i + 1][2 * i] = false;
                }
                (n8, Complex::clique_complex(&adj8, n8))
            };
            let (nv_bl, tets_bl) = kuhn_tets(2, 2, 2, false);
            let cbl = complex_from_tets(nv_bl, &tets_bl);
            vec![
                ("T3", ct3, nv_t3, [1, 3, 3, 1], "closed"),
                ("S3(16cell)", c16, nv16, [1, 0, 0, 1], "closed"),
                ("3-ball", cbl, nv_bl, [1, 0, 0, 0], "boundary"),
            ]
        };
        for (name, ctruth, nv, bexp, kind) in cases {
            // 隠れ h = −(1-skeleton 隣接)
            let adj_true = skeleton_adj(&ctruth, nv);
            let mut h = vec![0.0; nv * nv];
            for i in 0..nv {
                for j in 0..nv {
                    if adj_true[i][j] {
                        h[i * nv + j] = -1.0;
                    }
                }
            }
            // 測定 lane: 全 source の曲率重み (状態非依存 probe)
            let mut w = vec![0.0; nv * nv];
            for i in 0..nv {
                let wi = curvature_w(&h, nv, i, eps, dt);
                for j in 0..nv {
                    if j != i {
                        w[j * nv + i] = wi[j];
                    }
                }
            }
            let adj_est = support_from_weights(&w, nv);
            // 支持一致
            let mut miss = 0usize;
            let mut extra = 0usize;
            for i in 0..nv {
                for j in (i + 1)..nv {
                    match (adj_true[i][j], adj_est[i][j]) {
                        (true, false) => miss += 1,
                        (false, true) => extra += 1,
                        _ => {}
                    }
                }
            }
            // clique complex → β → link 分類
            let cest = Complex::clique_complex(&adj_est, nv);
            let best = cest.betti();
            let mut n_int = 0usize;
            let mut n_bnd = 0usize;
            let mut n_sing = 0usize;
            for v in 0..nv {
                match classify_surface(&vertex_link(&cest, v)) {
                    "S2" => n_int += 1,
                    "D2" => n_bnd += 1,
                    _ => n_sing += 1,
                }
            }
            let kind_est = if n_sing > 0 {
                "singular"
            } else if n_bnd == 0 {
                "closed"
            } else {
                "boundary"
            };
            let ok = miss == 0
                && extra == 0
                && best[..4] == bexp
                && kind_est == kind
                && n_sing == 0;
            all_ok &= ok;
            results.push(format!(
                "{}: 欠{}余{} β={:?} {} {}",
                name,
                miss,
                extra,
                &best[..4],
                kind_est,
                if ok { "✓" } else { "✗" }
            ));
        }
        check(
            "[V2] from-state 3D: 曲率則測定 lane → 支持 → clique → (β₀..β₃) + link 多様体性 — T³ (1,3,3,1)/S³ (1,0,0,1)/3-ball (1,0,0,0) を同定",
            all_ok,
            results.join(" / "),
        );
    }

    // ---- [V3] 円環の 1/3 法則 (計量 VR の定量アンカー) ----
    {
        let n = 24usize;
        let mut h = vec![0.0; n * n];
        let mut rng = Rng::new(316);
        let mut lens = Vec::new();
        for k in 0..n {
            let t = 0.8 + 0.4 * rng.f64();
            h[k * n + (k + 1) % n] = -t;
            h[((k + 1) % n) * n + k] = -t;
            lens.push(1.0 / t.sqrt());
        }
        // 測定重み → 計量 → VR
        let mut w = vec![0.0; n * n];
        for i in 0..n {
            let wi = curvature_w(&h, n, i, 0.3, 0.02);
            for j in 0..n {
                if j != i {
                    w[j * n + i] = wi[j];
                }
            }
        }
        let adj = support_from_weights(&w, n);
        let d = metric_closure(&w, &adj, n);
        let bars = vr_h1_bars(&d, n);
        let perim: f64 = lens.iter().sum();
        let main_bar = bars[0];
        let ratio = main_bar.1 / perim;
        let second = if bars.len() > 1 { bars[1].1 - bars[1].0 } else { 0.0 };
        check(
            "[V3] 円環の 1/3 法則: VR H1 主 bar の死亡/周長 ∈ [0.28, 0.37] (測地円の理論値 1/3)・第 2 bar は雑音床",
            ratio >= 0.28 && ratio <= 0.37 && second < 0.2 * (main_bar.1 - main_bar.0),
            format!(
                "death/周長 = {:.4} (bar [{:.3}, {:.3}], 周長 {:.3}) / 第 2 persistence {:.3}",
                ratio, main_bar.0, main_bar.1, perim, second
            ),
        );
    }

    // ---- [V4] 大小 2 穴の寿命分離 (v29.6 既知限界の解消) ----
    {
        // 8 字型: 大環 (16) と小環 (8) を 1 頂点で接着 — 穴スケール 2:1
        let big = 16usize;
        let small = 8usize;
        let n = big + small - 1;
        let mut h = vec![0.0; n * n];
        let mut bond = |a: usize, b: usize, h: &mut Vec<f64>| {
            h[a * n + b] = -1.0;
            h[b * n + a] = -1.0;
        };
        for k in 0..big {
            bond(k, (k + 1) % big, &mut h);
        }
        // 小環: 頂点 0 を共有し big..big+small-2 を巡回
        let sm: Vec<usize> = std::iter::once(0)
            .chain(big..big + small - 1)
            .collect();
        for k in 0..sm.len() {
            bond(sm[k], sm[(k + 1) % sm.len()], &mut h);
        }
        let mut w = vec![0.0; n * n];
        for i in 0..n {
            let wi = curvature_w(&h, n, i, 0.3, 0.02);
            for j in 0..n {
                if j != i {
                    w[j * n + i] = wi[j];
                }
            }
        }
        let adj = support_from_weights(&w, n);
        let d = metric_closure(&w, &adj, n);
        let bars = vr_h1_bars(&d, n);
        let p1 = bars[0].1 - bars[0].0;
        let p2 = if bars.len() > 1 { bars[1].1 - bars[1].0 } else { 0.0 };
        let p3 = if bars.len() > 2 { bars[2].1 - bars[2].0 } else { 0.0 };
        let ratio = p1 / p2.max(1e-12);
        check(
            "[V4] 大小 2 穴 (16:8): H1 bar 2 本が寿命で分離 (比 ∈ [1.5, 3])・第 3 bar は床 — 穴スケールの寿命分離 (v29.6 既知限界の解消)",
            bars.len() >= 2 && ratio >= 1.5 && ratio <= 3.0 && p3 < 0.3 * p2,
            format!(
                "persistence: {:.3} / {:.3} (比 {:.2}) / 第 3 {:.3}",
                p1, p2, ratio, p3
            ),
        );
    }

    // ---- [V5] genus-2 (状態から): β₁ = 4 + 閉曲面 ----
    {
        let (nv, tris) = genus2_complex();
        let ctruth = complex_from_tris(nv, &tris);
        let bt = ctruth.betti();
        // 構成の自己検査: (1,4,1)
        let construct_ok = bt[..3] == [1, 4, 1];
        // from state
        let adj_true = skeleton_adj(&ctruth, nv);
        let mut h = vec![0.0; nv * nv];
        for i in 0..nv {
            for j in 0..nv {
                if adj_true[i][j] {
                    h[i * nv + j] = -1.0;
                }
            }
        }
        let mut w = vec![0.0; nv * nv];
        for i in 0..nv {
            let wi = curvature_w(&h, nv, i, 0.3, 0.02);
            for j in 0..nv {
                if j != i {
                    w[j * nv + i] = wi[j];
                }
            }
        }
        let adj_est = support_from_weights(&w, nv);
        let mut miss = 0;
        let mut extra = 0;
        for i in 0..nv {
            for j in (i + 1)..nv {
                match (adj_true[i][j], adj_est[i][j]) {
                    (true, false) => miss += 1,
                    (false, true) => extra += 1,
                    _ => {}
                }
            }
        }
        let cest = Complex::clique_complex(&adj_est, nv);
        let be = cest.betti();
        // 曲面性 (2D link = 円): classify は 3D 用 — 2D は辺の三角形帰属 2 で閉
        let mut cnt = std::collections::HashMap::new();
        for t in &cest.simp[2] {
            for drop in 0..3 {
                let mut e: Vec<usize> = t.clone();
                e.remove(drop);
                *cnt.entry(e).or_insert(0usize) += 1;
            }
        }
        let closed_surface = cest.simp[1].iter().all(|e| cnt.get(e) == Some(&2));
        check(
            "[V5] genus-2 (2 T² の接着, 42 頂点): 構成 β = (1,4,1) + 状態から支持 欠0余0 → β₁ = 4・閉曲面 — 高 genus の初同定",
            construct_ok && miss == 0 && extra == 0 && be[..3] == [1, 4, 1] && closed_surface,
            format!(
                "構成 β = {:?} / 支持 欠{}余{} / 推定 β = {:?} / 閉曲面 = {} (頂点 {})",
                &bt[..3],
                miss,
                extra,
                &be[..3],
                closed_surface,
                nv
            ),
        );
    }

    // ---- [V6] 敵対対照: narrow neck / close holes ----
    {
        // narrow neck: 2 つの完全三角形化された小円板 (K4) を 2 辺の橋で接続 — β₁ = 0
        let n = 10usize;
        let mut adj = vec![vec![false; n]; n];
        let mut e = |a: usize, b: usize, adj: &mut Vec<Vec<bool>>| {
            adj[a][b] = true;
            adj[b][a] = true;
        };
        for i in 0..4 {
            for j in (i + 1)..4 {
                e(i, j, &mut adj);
            }
        }
        for i in 5..9 {
            for j in (i + 1)..9 {
                e(i, j, &mut adj);
            }
        }
        e(3, 4, &mut adj);
        e(4, 5, &mut adj);
        let mut h = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                if adj[i][j] {
                    h[i * n + j] = -1.0;
                }
            }
        }
        let mut w = vec![0.0; n * n];
        for i in 0..n {
            let wi = curvature_w(&h, n, i, 0.3, 0.02);
            for j in 0..n {
                if j != i {
                    w[j * n + i] = wi[j];
                }
            }
        }
        let adj_est = support_from_weights(&w, n);
        let d = metric_closure(&w, &adj_est, n);
        let bars = vr_h1_bars(&d, n);
        let max_p = bars.first().map(|b| b.1 - b.0).unwrap_or(0.0);
        // close holes: 同サイズ 2 穴が近接 (theta グラフ様 8 字, 8+8 で共有 1 頂点)
        let n2 = 15usize;
        let mut h2 = vec![0.0; n2 * n2];
        let mut bond2 = |a: usize, b: usize, h: &mut Vec<f64>| {
            h[a * n2 + b] = -1.0;
            h[b * n2 + a] = -1.0;
        };
        for k in 0..8 {
            bond2(k, (k + 1) % 8, &mut h2);
        }
        let sm: Vec<usize> = std::iter::once(0).chain(8..15).collect();
        for k in 0..sm.len() {
            bond2(sm[k], sm[(k + 1) % sm.len()], &mut h2);
        }
        let mut w2 = vec![0.0; n2 * n2];
        for i in 0..n2 {
            let wi = curvature_w(&h2, n2, i, 0.3, 0.02);
            for j in 0..n2 {
                if j != i {
                    w2[j * n2 + i] = wi[j];
                }
            }
        }
        let adj2 = support_from_weights(&w2, n2);
        let d2m = metric_closure(&w2, &adj2, n2);
        let bars2 = vr_h1_bars(&d2m, n2);
        let (q1, q2) = (
            bars2.first().map(|b| b.1 - b.0).unwrap_or(0.0),
            bars2.get(1).map(|b| b.1 - b.0).unwrap_or(0.0),
        );
        check(
            "[V6] 敵対対照: narrow neck (K4 対 + 橋) は H1 主 persistence < 0.5 (穴を捏造しない) / 近接同サイズ 2 穴は bar 2 本が同寿命で残る (合体しない)",
            max_p < 0.5 && bars2.len() >= 2 && q2 > 0.7 * q1,
            format!(
                "neck max persistence = {:.3} / 近接 2 穴 persistence = {:.3}, {:.3} (比 {:.2})",
                max_p,
                q1,
                q2,
                q2 / q1.max(1e-12)
            ),
        );
    }

    // ---- [V7] ノイズ安定性: bottleneck(D, D') ≤ sup|Δd| (安定性定理の機械照合) ----
    {
        // [V4] の 8 字系で辺長に乗法ノイズ → 計量摂動 → bottleneck
        let big = 16usize;
        let small = 8usize;
        let n = big + small - 1;
        let mut h = vec![0.0; n * n];
        let mut bond = |a: usize, b: usize, h: &mut Vec<f64>| {
            h[a * n + b] = -1.0;
            h[b * n + a] = -1.0;
        };
        for k in 0..big {
            bond(k, (k + 1) % big, &mut h);
        }
        let sm: Vec<usize> = std::iter::once(0).chain(big..big + small - 1).collect();
        for k in 0..sm.len() {
            bond(sm[k], sm[(k + 1) % sm.len()], &mut h);
        }
        let mut w = vec![0.0; n * n];
        for i in 0..n {
            let wi = curvature_w(&h, n, i, 0.3, 0.02);
            for j in 0..n {
                if j != i {
                    w[j * n + i] = wi[j];
                }
            }
        }
        let adj = support_from_weights(&w, n);
        let d0 = metric_closure(&w, &adj, n);
        let bars0: Vec<(f64, f64)> = vr_h1_bars(&d0, n).into_iter().take(4).collect();
        let mut rng = Rng::new(3167);
        let mut all_ok = true;
        let mut detail = String::new();
        for &sigma in &[0.02, 0.05] {
            let mut wn = w.clone();
            for i in 0..n {
                for j in (i + 1)..n {
                    let f = 1.0 + sigma * (2.0 * rng.f64() - 1.0);
                    wn[i * n + j] *= f;
                    wn[j * n + i] = wn[i * n + j];
                }
            }
            let dn = metric_closure(&wn, &adj, n);
            let mut supd: f64 = 0.0;
            for k in 0..n * n {
                supd = supd.max((dn[k] - d0[k]).abs());
            }
            let barsn: Vec<(f64, f64)> = vr_h1_bars(&dn, n).into_iter().take(4).collect();
            let bd = bottleneck(&bars0, &barsn);
            let ok = bd <= supd + 1e-12;
            all_ok &= ok;
            detail.push_str(&format!(
                "σ={}: bottleneck {:.4} ≤ sup|Δd| {:.4} {} ",
                sigma,
                bd,
                supd,
                if ok { "✓" } else { "✗" }
            ));
        }
        check(
            "[V7] ノイズ安定性: bottleneck(H1 図) ≤ sup|Δd| (persistence 安定性定理の機械照合, σ = 2%/5%)",
            all_ok,
            detail,
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "幾何読み出しが 3 次元と計量に到達した — β₃ (im ∂₄ 必須) と link 多様体性で T³/S³/3-ball を状態から同定し、事前登録計量の VR persistence が円環 1/3 法則・穴スケールの寿命分離・genus-2 を定量化し、安定性定理が bottleneck で機械照合された"
        } else {
            "**3D/計量読み出しの破れ**"
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
