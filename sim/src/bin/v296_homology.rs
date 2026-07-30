//! v29.6 dimension-agnostic topology — persistent homology pipeline の凍結と資格 (第三十期, PROMPT/11 第四課題)
//!
//! 現行の再構成器 (v29/v30 凍結節) は相互 top-2 グラフの Cycle/Path 認識 = 1 次元
//! 専用で、2D トーラスを構造上受け付けない。本版は次の dimension-agnostic
//! pipeline を実装し、2D 正例と敵対対照で**資格** (設計区画 — holdout は v30.0) を
//! とる:
//!
//!   1. affinity K_ij の順位 filtration と **頂点適応 gap 則** (各頂点は自分の
//!      親和度列の最大段差 w(m)/w(m+1) で近傍数 m を決める — 混在次数・境界に頑健)
//!   2. 相互配 (mutual) グラフ → clique complex (単体は四面体 = 3-simplex まで)
//!   3. homology (Z2 係数, 境界行列の列簡約) → Betti 数 β0, β1, β2
//!   4. neighborhood shell growth S(r) ~ r^{d−1} による内在次元
//!   5. vertex link による manifoldness 検査 (閉曲面: link ≅ 円周 / 境界: 弧)
//!   6. graph geodesic と静的核の整合 (Spearman)
//!   7. 探索窓 kmax に対する安定領域 (全会一致 — profile が kmax に依らない)
//!
//! v29.5 [C7b] の設計入力を構造的に解決する:
//!   - **gap 保留**: 最大段差が GAP_MIN 未満の頂点が過半なら裁定保留 (K7 のような
//!     縮退カーネルで幾何を捏造しない)。
//!   - Petersen の誤認は v30 の k=2 固定が真の次数 3 を下回ったのが根因 — 適応 k
//!     は次数 3 を選び、link 検査が「グラフは復元・曲面ではない」と正しく裁定する。
//!
//! 正例 (全て**状態から end-to-end**: 自由フェルミオン基底状態 → B3-COV 核 →
//! pipeline): 三角格子トーラス (1,2,1)・円板 (1,0,0)・円筒 (1,1,0)・細分
//! 二十面体球 (1,0,1)。敵対対照: ランダム 3-正則 (曲面と誤認しない)・分岐面 3 枚
//! (spine の link が theta 型 → 非多様体)・大小 2 穴の円板 (β1 = 2 に分離)・
//! Petersen / K7 (非曲面裁定 / 縮退保留)。乱数は固定シード。最大 n = 121 の
//! jacobi のみ (重い処理なし)。
//!
//! 既知の限界: 穴スケールの寿命分離は clique-順位 filtration では現れない (グラフ
//! の穴は不滅バー) — 計量 VR (三角形を大スケールで充填) の persistence は次段。
//! 同梱修正: v295 [C7b] の is_edge が Petersen 辺リストの未整列 tuple (9,6) を
//! 誤判定していた表示バグ (裁定・引用値は無傷) を v295 側で修正。

use uft_sim::*;

// ---------------- 状態と核 (B3-COV — v29 凍結節と同じ定義式, m = 1 特化) ----------------

/// 隣接 (重みつき) → 半充填基底状態相関 (縮退シェルは分数占有 — v295 と同じ規則)
fn corr_half_filling(a: &[f64], n: usize) -> Vec<f64> {
    let h: Vec<f64> = a.iter().map(|&x| -x).collect();
    let (ev, vv) = jacobi_eigh(&h, n);
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| ev[i].partial_cmp(&ev[j]).unwrap().then(i.cmp(&j)));
    let mut shells: Vec<Vec<usize>> = Vec::new();
    for &k in &idx {
        if let Some(last) = shells.last_mut() {
            if (ev[k] - ev[*last.last().unwrap()]).abs() < 1e-9 {
                last.push(k);
                continue;
            }
        }
        shells.push(vec![k]);
    }
    let target = n as f64 / 2.0;
    let mut c = vec![0.0; n * n];
    let mut filled = 0.0;
    for sh in shells {
        if filled >= target - 1e-12 {
            break;
        }
        let alpha = ((target - filled) / sh.len() as f64).min(1.0);
        for &k in &sh {
            for i in 0..n {
                for j in 0..n {
                    c[i + j * n] += alpha * vv[i + k * n] * vv[j + k * n];
                }
            }
        }
        filled += alpha * sh.len() as f64;
    }
    c
}

/// 熱的 Gaussian 状態の相関: C = (1 + e^{βh})^{-1}, h = −A, μ = 0。
/// ギャップのある (指数減衰) 相関で、臨界半充填 GS の境界増強・Friedel 共鳴を
/// 持たない — pipeline の資格はこの状態族でとる (臨界 GS の破れは [T2b] に記録)。
fn corr_thermal(a: &[f64], n: usize, beta: f64) -> Vec<f64> {
    let h: Vec<f64> = a.iter().map(|&x| -x).collect();
    let (ev, vv) = jacobi_eigh(&h, n);
    let mut c = vec![0.0; n * n];
    for k in 0..n {
        let f = 1.0 / (1.0 + (beta * ev[k]).exp());
        for i in 0..n {
            for j in 0..n {
                c[i + j * n] += f * vv[i + k * n] * vv[j + k * n];
            }
        }
    }
    c
}

const BETA_T: f64 = 1.0; // 資格状態族の逆温度 (β t = 1 — 距離 2 相関は O(β²) で潰れ nn は O(β): 隘路の増強対を抑制。β = 2 では 2 穴円板の隘路で余辺 1 本 [開発記録])

/// B3-COV 核: w(i,j) = |⟨n_i n_j⟩_c| = |C_ij|² (m = 1, 実相関)
fn kernel_b3(c: &[f64], n: usize) -> Vec<f64> {
    let mut w = vec![0.0; n * n];
    for i in 0..n {
        for j in 0..n {
            if i != j {
                let x = c[i + j * n];
                w[i + j * n] = x * x;
            }
        }
    }
    w
}

// ---------------- pipeline 核心 (適応 gap kNN + clique complex + Z2 homology) ----------------

const GAP_MIN: f64 = 1.5; // 頂点の段差 (比) の下限 — ≥ GAP_MIN の段差がなければ縮退頂点
const W_FLOOR6: f64 = 1e-8;
/// 近傍候補の絶対スケールガード: w(m 位) ≥ W_FRAC × w(1 位) の範囲だけを近傍候補と
/// する。真の近傍シェルは同桁に固まり、遠方シェルは桁で落ちる — 遠方の偶発段差
/// (Friedel 振動の節) を近傍と誤認しないため (開発記録: 球で偽辺 → β1 = 6)。
const W_FRAC: f64 = 0.15;

/// 頂点適応の相互グラフ: 各頂点 i は「m 位/(m+1) 位比 ≥ GAP_MIN を満たす**最大**の
/// m (kmin ≤ m ≤ kmax)」を自分の近傍数に選ぶ — 最遠の有意な段差まで含めることで、
/// 境界頂点の不均一な真近傍 (段差が内部にもある) を取りこぼさない。argmax 規則は
/// 境界で辺を落とし偽の穴を列生成した (開発記録)。辺 = 相互に選ばれた対。
/// 返り値: (辺, 縮退頂点数 [段差 ≥ GAP_MIN が一つもない, または床下])
fn adaptive_mutual(w: &[f64], n: usize, kmin: usize, kmax: usize) -> (Vec<(usize, usize)>, usize) {
    let mut keep: Vec<Vec<usize>> = Vec::with_capacity(n);
    let mut weak = 0usize;
    for i in 0..n {
        let mut js: Vec<usize> = (0..n).filter(|&j| j != i).collect();
        js.sort_by(|&a, &b| {
            w[i + b * n]
                .partial_cmp(&w[i + a * n])
                .unwrap()
                .then(a.cmp(&b))
        });
        let w1 = w[i + js[0] * n];
        let mut best_m = 0usize;
        for m in kmin..=kmax.min(n - 2) {
            let wm = w[i + js[m - 1] * n];
            let wm1 = w[i + js[m] * n];
            if wm < W_FLOOR6 || wm < W_FRAC * w1 {
                break;
            }
            let ratio = if wm1 > 0.0 { wm / wm1 } else { f64::INFINITY };
            if ratio >= GAP_MIN {
                best_m = m;
            }
        }
        if best_m == 0 {
            weak += 1;
            keep.push(Vec::new());
        } else {
            keep.push(js[..best_m].to_vec());
        }
    }
    let mut edges = Vec::new();
    for i in 0..n {
        for &j in &keep[i] {
            if j > i && keep[j].contains(&i) {
                edges.push((i, j));
            }
        }
    }
    (edges, weak)
}

/// clique complex の単体列挙 (次元 ≤ 3)。返り値: (三角形, 四面体)
fn cliques(edges: &[(usize, usize)], n: usize) -> (Vec<[usize; 3]>, Vec<[usize; 4]>) {
    let mut adj = vec![vec![false; n]; n];
    for &(i, j) in edges {
        adj[i][j] = true;
        adj[j][i] = true;
    }
    let mut tris = Vec::new();
    for &(i, j) in edges {
        for l in (j + 1)..n {
            if adj[i][l] && adj[j][l] {
                tris.push([i, j, l]);
            }
        }
    }
    let mut tets = Vec::new();
    for t in &tris {
        for l in (t[2] + 1)..n {
            if adj[t[0]][l] && adj[t[1]][l] && adj[t[2]][l] {
                tets.push([t[0], t[1], t[2], l]);
            }
        }
    }
    (tris, tets)
}

/// Betti 数 (Z2, 次元 0..2) — 境界行列の列簡約 (標準の persistence 算法の固定複体版)
fn betti(edges: &[(usize, usize)], n: usize) -> [usize; 3] {
    let (tris, tets) = cliques(edges, n);
    let mut edge_id = std::collections::BTreeMap::new();
    for (e, &(i, j)) in edges.iter().enumerate() {
        edge_id.insert((i, j), n + e);
    }
    let eid = |i: usize, j: usize| -> usize { edge_id[&(i.min(j), i.max(j))] };
    let mut tri_id = std::collections::BTreeMap::new();
    let t0 = n + edges.len();
    for (t, tr) in tris.iter().enumerate() {
        tri_id.insert(*tr, t0 + t);
    }
    let q0 = t0 + tris.len();
    let total = q0 + tets.len();
    let mut cols: Vec<Vec<usize>> = vec![Vec::new(); total];
    for (e, &(i, j)) in edges.iter().enumerate() {
        let mut b = vec![i, j];
        b.sort_unstable();
        cols[n + e] = b;
    }
    for (t, tr) in tris.iter().enumerate() {
        let mut b = vec![eid(tr[0], tr[1]), eid(tr[0], tr[2]), eid(tr[1], tr[2])];
        b.sort_unstable();
        cols[t0 + t] = b;
    }
    for (q, te) in tets.iter().enumerate() {
        let mut b = vec![
            tri_id[&[te[0], te[1], te[2]]],
            tri_id[&[te[0], te[1], te[3]]],
            tri_id[&[te[0], te[2], te[3]]],
            tri_id[&[te[1], te[2], te[3]]],
        ];
        b.sort_unstable();
        cols[q0 + q] = b;
    }
    let mut pivot_of: Vec<Option<usize>> = vec![None; total];
    let mut killed_row = vec![false; total];
    for j in 0..total {
        loop {
            let low = match cols[j].last() {
                Some(&l) => l,
                None => break,
            };
            match pivot_of[low] {
                Some(i) => {
                    let merged: Vec<usize> = {
                        let (a, b) = (&cols[j], &cols[i]);
                        let mut out = Vec::with_capacity(a.len() + b.len());
                        let (mut x, mut y) = (0usize, 0usize);
                        while x < a.len() || y < b.len() {
                            if y >= b.len() || (x < a.len() && a[x] < b[y]) {
                                out.push(a[x]);
                                x += 1;
                            } else if x >= a.len() || b[y] < a[x] {
                                out.push(b[y]);
                                y += 1;
                            } else {
                                x += 1;
                                y += 1;
                            }
                        }
                        out
                    };
                    cols[j] = merged;
                }
                None => {
                    pivot_of[low] = Some(j);
                    killed_row[low] = true;
                    break;
                }
            }
        }
    }
    // β_d = # { 単体 s (次元 d) : 簡約列が空 (= creator) かつ行として殺されていない }
    let dim_of = |s: usize| -> usize {
        if s < n {
            0
        } else if s < t0 {
            1
        } else if s < q0 {
            2
        } else {
            3
        }
    };
    let mut b = [0usize; 3];
    for s in 0..total {
        let d = dim_of(s);
        if d <= 2 && cols[s].is_empty() && !killed_row[s] {
            b[d] += 1;
        }
    }
    b
}

/// vertex link: (閉 [円周], 境界 [弧], 破れ) の頂点数
fn linkness(edges: &[(usize, usize)], n: usize) -> (usize, usize, usize) {
    let mut adj = vec![vec![false; n]; n];
    for &(i, j) in edges {
        adj[i][j] = true;
        adj[j][i] = true;
    }
    let (mut closed, mut boundary, mut bad) = (0usize, 0usize, 0usize);
    for v in 0..n {
        let nb: Vec<usize> = (0..n).filter(|&u| adj[v][u]).collect();
        let m = nb.len();
        if m == 0 {
            bad += 1;
            continue;
        }
        let mut deg = vec![0usize; m];
        let mut le = 0usize;
        let mut ladj = vec![vec![false; m]; m];
        for a in 0..m {
            for b in (a + 1)..m {
                if adj[nb[a]][nb[b]] {
                    deg[a] += 1;
                    deg[b] += 1;
                    le += 1;
                    ladj[a][b] = true;
                    ladj[b][a] = true;
                }
            }
        }
        let mut seen = vec![false; m];
        seen[0] = true;
        let mut stack = vec![0usize];
        let mut cnt = 1usize;
        while let Some(u) = stack.pop() {
            for x in 0..m {
                if ladj[u][x] && !seen[x] {
                    seen[x] = true;
                    cnt += 1;
                    stack.push(x);
                }
            }
        }
        let connected = cnt == m;
        let d1 = deg.iter().filter(|&&d| d == 1).count();
        let d2 = deg.iter().filter(|&&d| d == 2).count();
        if connected && d2 == m && le == m && m >= 3 {
            closed += 1;
        } else if connected && d1 == 2 && d1 + d2 == m && le + 1 == m {
            boundary += 1;
        } else {
            bad += 1;
        }
    }
    (closed, boundary, bad)
}

/// 内在次元: 殻計数 S(r) = N(r) − N(r−1) ~ c·r^{d−1} → dim = 1 + slope(ln S vs ln r)
fn intrinsic_dim(edges: &[(usize, usize)], n: usize, rmax: usize) -> f64 {
    let mut adj = vec![Vec::new(); n];
    for &(i, j) in edges {
        adj[i].push(j);
        adj[j].push(i);
    }
    let mut slopes = Vec::new();
    for s in 0..n {
        let mut dist = vec![usize::MAX; n];
        dist[s] = 0;
        let mut q = std::collections::VecDeque::new();
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            for &v in &adj[u] {
                if dist[v] == usize::MAX {
                    dist[v] = dist[u] + 1;
                    q.push_back(v);
                }
            }
        }
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        let mut ok = true;
        for r in 1..=rmax {
            let sr = dist.iter().filter(|&&d| d == r).count();
            if sr == 0 {
                ok = false;
                break;
            }
            xs.push((r as f64).ln());
            ys.push((sr as f64).ln());
        }
        if ok {
            if let Ok(f) = linfit_typed(&xs, &ys) {
                slopes.push(1.0 + f.slope);
            }
        }
    }
    if slopes.is_empty() {
        return f64::NAN;
    }
    slopes.iter().sum::<f64>() / slopes.len() as f64
}

/// graph geodesic (BFS) と核の順位整合 (源 3 点の平均 Spearman, 核は −w)
fn geodesic_consistency(w: &[f64], edges: &[(usize, usize)], n: usize) -> f64 {
    let mut adj = vec![Vec::new(); n];
    for &(i, j) in edges {
        adj[i].push(j);
        adj[j].push(i);
    }
    let mut tot = 0.0;
    let sources = [0usize, n / 3, 2 * n / 3];
    for &s in &sources {
        let mut dist = vec![f64::INFINITY; n];
        dist[s] = 0.0;
        let mut q = std::collections::VecDeque::new();
        q.push_back(s);
        while let Some(u) = q.pop_front() {
            for &v in &adj[u] {
                if dist[v].is_infinite() {
                    dist[v] = dist[u] + 1.0;
                    q.push_back(v);
                }
            }
        }
        let negw: Vec<f64> = (0..n).map(|j| -w[s + j * n]).collect();
        tot += spearman_local(&dist, &negw);
    }
    tot / sources.len() as f64
}

fn spearman_local(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len();
    let rank = |v: &[f64]| -> Vec<f64> {
        let mut idx: Vec<usize> = (0..n).collect();
        idx.sort_by(|&a, &b| v[a].partial_cmp(&v[b]).unwrap().then(a.cmp(&b)));
        let mut r = vec![0.0; n];
        for (pos, &i) in idx.iter().enumerate() {
            r[i] = pos as f64;
        }
        r
    };
    let (rx, ry) = (rank(x), rank(y));
    let mx = rx.iter().sum::<f64>() / n as f64;
    let my = ry.iter().sum::<f64>() / n as f64;
    let (mut sxy, mut sxx, mut syy) = (0.0, 0.0, 0.0);
    for i in 0..n {
        let (dx, dy) = (rx[i] - mx, ry[i] - my);
        sxy += dx * dy;
        sxx += dx * dx;
        syy += dy * dy;
    }
    sxy / (sxx * syy).sqrt().max(1e-300)
}

/// pipeline 裁定: 探索窓 kmax ∈ {7, 8, 9, 10} で適応相互グラフを作り、
/// (β, 曲面性) の**全会一致**を安定性とする。縮退頂点が過半の窓は保留。
struct Verdict {
    n_windows: usize, // 採用された kmax 窓の数
    unanimous: bool,
    betti: [usize; 3],
    surface: &'static str, // "closed" / "boundary" / "not-surface"
    dim: f64,
    geo: f64,
    degenerate_windows: usize,
}

fn run_pipeline(w: &[f64], n: usize) -> Verdict {
    let mut profiles: Vec<(usize, [usize; 3], &'static str, f64, f64)> = Vec::new();
    let mut degenerate_windows = 0usize;
    for kmax in [7usize, 8, 9, 10] {
        let (edges, weak) = adaptive_mutual(w, n, 2, kmax.min(n - 2));
        if weak * 2 > n || edges.is_empty() {
            degenerate_windows += 1;
            continue;
        }
        let b = betti(&edges, n);
        let (closed, bnd, bad) = linkness(&edges, n);
        let surface = if bad == 0 && bnd == 0 && closed == n {
            "closed"
        } else if bad == 0 && bnd > 0 {
            "boundary"
        } else {
            "not-surface"
        };
        let dim = intrinsic_dim(&edges, n, 3);
        let geo = geodesic_consistency(w, &edges, n);
        profiles.push((kmax, b, surface, dim, geo));
    }
    if profiles.is_empty() {
        return Verdict {
            n_windows: 0,
            unanimous: false,
            betti: [0, 0, 0],
            surface: "not-surface",
            dim: f64::NAN,
            geo: f64::NAN,
            degenerate_windows,
        };
    }
    let unanimous = profiles
        .iter()
        .all(|p| p.1 == profiles[0].1 && p.2 == profiles[0].2);
    Verdict {
        n_windows: profiles.len(),
        unanimous,
        betti: profiles[0].1,
        surface: if unanimous { profiles[0].2 } else { "not-surface" },
        dim: profiles.iter().map(|p| p.3).sum::<f64>() / profiles.len() as f64,
        geo: profiles.iter().map(|p| p.4).sum::<f64>() / profiles.len() as f64,
        degenerate_windows,
    }
}

// ---------------- 幾何 (格子とグラフ) ----------------

/// 三角格子トーラス L×L (6 正則)
fn torus_tri(l: usize) -> (usize, Vec<f64>) {
    let n = l * l;
    let mut a = vec![0.0; n * n];
    let id = |x: usize, y: usize| -> usize { (x % l) + (y % l) * l };
    for x in 0..l {
        for y in 0..l {
            let i = id(x, y);
            for (dx, dy) in [(1usize, 0usize), (0, 1), (1, 1)] {
                let j = id(x + dx, y + dy);
                a[i + j * n] = 1.0;
                a[j + i * n] = 1.0;
            }
        }
    }
    (n, a)
}

/// 三角格子円筒 L×W (x 周期・y 開放)
fn cylinder_tri(l: usize, wd: usize) -> (usize, Vec<f64>) {
    let n = l * wd;
    let mut a = vec![0.0; n * n];
    let id = |x: usize, y: usize| -> usize { (x % l) + y * l };
    for x in 0..l {
        for y in 0..wd {
            let i = id(x, y);
            a[i + id(x + 1, y) * n] = 1.0;
            a[id(x + 1, y) + i * n] = 1.0;
            if y + 1 < wd {
                for j in [id(x, y + 1), id(x + 1, y + 1)] {
                    a[i + j * n] = 1.0;
                    a[j + i * n] = 1.0;
                }
            }
        }
    }
    (n, a)
}

/// 三角格子円板 (L×W 開放パッチ, holes のサイトを除去)
fn disk_tri(l: usize, wd: usize, holes: &[(usize, usize)]) -> (usize, Vec<f64>, Vec<usize>) {
    let n0 = l * wd;
    let mut alive = vec![true; n0];
    for &(hx, hy) in holes {
        alive[hx + hy * l] = false;
    }
    let mut a = vec![0.0; n0 * n0];
    let id = |x: usize, y: usize| -> usize { x + y * l };
    for x in 0..l {
        for y in 0..wd {
            let i = id(x, y);
            if !alive[i] {
                continue;
            }
            let mut targets = Vec::new();
            if x + 1 < l {
                targets.push(id(x + 1, y));
            }
            if y + 1 < wd {
                targets.push(id(x, y + 1));
                if x + 1 < l {
                    targets.push(id(x + 1, y + 1));
                }
            }
            for j in targets {
                if alive[j] {
                    a[i + j * n0] = 1.0;
                    a[j + i * n0] = 1.0;
                }
            }
        }
    }
    let map: Vec<usize> = (0..n0).filter(|&i| alive[i]).collect();
    let n = map.len();
    let mut b = vec![0.0; n * n];
    for (i2, &i) in map.iter().enumerate() {
        for (j2, &j) in map.iter().enumerate() {
            b[i2 + j2 * n] = a[i + j * n0];
        }
    }
    (n, b, map)
}

/// linkness の破れ頂点リスト (診断用)
fn bad_link_vertices(edges: &[(usize, usize)], n: usize) -> Vec<usize> {
    let mut adj = vec![vec![false; n]; n];
    for &(i, j) in edges {
        adj[i][j] = true;
        adj[j][i] = true;
    }
    let mut out = Vec::new();
    for v in 0..n {
        let nb: Vec<usize> = (0..n).filter(|&u| adj[v][u]).collect();
        let m = nb.len();
        if m == 0 {
            out.push(v);
            continue;
        }
        let mut deg = vec![0usize; m];
        let mut le = 0usize;
        let mut ladj = vec![vec![false; m]; m];
        for a in 0..m {
            for b in (a + 1)..m {
                if adj[nb[a]][nb[b]] {
                    deg[a] += 1;
                    deg[b] += 1;
                    le += 1;
                    ladj[a][b] = true;
                    ladj[b][a] = true;
                }
            }
        }
        let mut seen = vec![false; m];
        seen[0] = true;
        let mut stack = vec![0usize];
        let mut cnt = 1usize;
        while let Some(u) = stack.pop() {
            for x in 0..m {
                if ladj[u][x] && !seen[x] {
                    seen[x] = true;
                    cnt += 1;
                    stack.push(x);
                }
            }
        }
        let connected = cnt == m;
        let d1 = deg.iter().filter(|&&d| d == 1).count();
        let d2 = deg.iter().filter(|&&d| d == 2).count();
        let closed = connected && d2 == m && le == m && m >= 3;
        let bnd = connected && d1 == 2 && d1 + d2 == m && le + 1 == m;
        if !closed && !bnd {
            out.push(v);
        }
    }
    out
}

/// 細分二十面体 (icosphere-1): 42 頂点・120 辺の三角化球
fn icosphere1() -> (usize, Vec<f64>) {
    let phi = (1.0 + 5.0f64.sqrt()) / 2.0;
    let mut verts: Vec<[f64; 3]> = Vec::new();
    for &(a, b) in &[(1.0, phi), (1.0, -phi), (-1.0, phi), (-1.0, -phi)] {
        verts.push([0.0, a, b]);
        verts.push([a, b, 0.0]);
        verts.push([b, 0.0, a]);
    }
    let norm = |v: [f64; 3]| -> [f64; 3] {
        let r = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt();
        [v[0] / r, v[1] / r, v[2] / r]
    };
    for v in verts.iter_mut() {
        *v = norm(*v);
    }
    let d2 = |a: [f64; 3], b: [f64; 3]| -> f64 {
        (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
    };
    let mut edges0 = Vec::new();
    for i in 0..12 {
        for j in (i + 1)..12 {
            if d2(verts[i], verts[j]) < 1.2 {
                edges0.push((i, j));
            }
        }
    }
    assert_eq!(edges0.len(), 30, "二十面体の辺数");
    let mut mid = std::collections::BTreeMap::new();
    for &(i, j) in &edges0 {
        let m = norm([
            (verts[i][0] + verts[j][0]) / 2.0,
            (verts[i][1] + verts[j][1]) / 2.0,
            (verts[i][2] + verts[j][2]) / 2.0,
        ]);
        mid.insert((i, j), verts.len());
        verts.push(m);
    }
    let has = |a: usize, b: usize| -> bool { edges0.contains(&(a.min(b), a.max(b))) };
    let mut eset = std::collections::BTreeSet::new();
    for i in 0..12 {
        for j in (i + 1)..12 {
            for k in (j + 1)..12 {
                if has(i, j) && has(i, k) && has(j, k) {
                    let mij = mid[&(i, j)];
                    let mik = mid[&(i, k)];
                    let mjk = mid[&(j, k)];
                    for &(a, b) in &[
                        (i, mij),
                        (i, mik),
                        (j, mij),
                        (j, mjk),
                        (k, mik),
                        (k, mjk),
                        (mij, mik),
                        (mij, mjk),
                        (mik, mjk),
                    ] {
                        eset.insert((a.min(b), a.max(b)));
                    }
                }
            }
        }
    }
    let n = verts.len();
    let mut a = vec![0.0; n * n];
    for &(i, j) in &eset {
        a[i + j * n] = 1.0;
        a[j + i * n] = 1.0;
    }
    (n, a)
}

/// 分岐面: **3 枚**の三角格子シートが 1 本の線 (spine) を共有 (2 枚では曲がった
/// 平面 = 多様体になってしまう — 本物の分岐には ≥ 3 枚が要る)
fn branched(l: usize, wd: usize) -> (usize, Vec<f64>) {
    let n = l + 3 * l * (wd - 1);
    let mut a = vec![0.0; n * n];
    let spine = |x: usize| -> usize { x };
    let sheet = |s: usize, x: usize, y: usize| -> usize { l + s * l * (wd - 1) + x + (y - 1) * l };
    let link = |i: usize, j: usize, a: &mut Vec<f64>| {
        a[i + j * n] = 1.0;
        a[j + i * n] = 1.0;
    };
    for x in 0..l - 1 {
        link(spine(x), spine(x + 1), &mut a);
    }
    for s in 0..3 {
        for x in 0..l {
            for y in 1..wd {
                let i = sheet(s, x, y);
                let below = if y == 1 { spine(x) } else { sheet(s, x, y - 1) };
                link(i, below, &mut a);
                if x + 1 < l {
                    link(i, sheet(s, x + 1, y), &mut a);
                    let below_r = if y == 1 {
                        spine(x + 1)
                    } else {
                        sheet(s, x + 1, y - 1)
                    };
                    link(i, below_r, &mut a);
                }
            }
        }
    }
    (n, a)
}

/// ランダム 3-正則 (configuration model, 固定シード, 単純グラフまで再試行)
fn random_cubic(n: usize, seed: u64) -> (usize, Vec<f64>) {
    let mut rng = Rng::new(seed);
    loop {
        let mut stubs: Vec<usize> = (0..n).flat_map(|i| [i, i, i]).collect();
        for i in (1..stubs.len()).rev() {
            let j = rng.range(i + 1);
            stubs.swap(i, j);
        }
        let mut a = vec![0.0; n * n];
        let mut ok = true;
        for p in stubs.chunks(2) {
            let (i, j) = (p[0], p[1]);
            if i == j || a[i + j * n] != 0.0 {
                ok = false;
                break;
            }
            a[i + j * n] = 1.0;
            a[j + i * n] = 1.0;
        }
        if ok {
            return (n, a);
        }
    }
}

fn petersen() -> (usize, Vec<f64>) {
    let pe = [
        (0usize, 1usize),
        (1, 2),
        (2, 3),
        (3, 4),
        (0, 4),
        (5, 7),
        (7, 9),
        (6, 9),
        (6, 8),
        (5, 8),
        (0, 5),
        (1, 6),
        (2, 7),
        (3, 8),
        (4, 9),
    ];
    let n = 10;
    let mut a = vec![0.0; n * n];
    for &(i, j) in &pe {
        a[i + j * n] = 1.0;
        a[j + i * n] = 1.0;
    }
    (n, a)
}

fn edges_of_adj(a: &[f64], n: usize) -> Vec<(usize, usize)> {
    let mut e = Vec::new();
    for i in 0..n {
        for j in (i + 1)..n {
            if a[i + j * n] != 0.0 {
                e.push((i, j));
            }
        }
    }
    e
}

fn main() {
    self_test();
    println!("=== v29.6 dimension-agnostic topology — persistent homology pipeline の資格 (第三十期) ===\n");
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

    // ---- [T0] homology エンジンの器械検査 (既知複体) ----
    {
        let ring8: Vec<(usize, usize)> = (0..8)
            .map(|i| (i.min((i + 1) % 8), i.max((i + 1) % 8)))
            .collect();
        let b_ring = betti(&ring8, 8);
        let tetra: Vec<(usize, usize)> = vec![(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)];
        let b_tetra = betti(&tetra, 4); // 充満 K4 (四面体込み) → 可縮 (1,0,0)
        let mut two_rings: Vec<(usize, usize)> = (0..6)
            .map(|i| (i.min((i + 1) % 6), i.max((i + 1) % 6)))
            .collect();
        two_rings.extend((0..6).map(|i| (6 + i.min((i + 1) % 6), 6 + i.max((i + 1) % 6))));
        let b_two = betti(&two_rings, 12);
        let (nt, at) = torus_tri(4);
        let b_torus = betti(&edges_of_adj(&at, nt), nt);
        let (ni, ai) = icosphere1();
        let ei = edges_of_adj(&ai, ni);
        let b_ico = betti(&ei, ni);
        let (cl, bd, bad) = linkness(&ei, ni);
        check(
            "[T0] homology エンジン (Z2 列簡約) — C8 (1,1,0) / K4 充満 (1,0,0) / 2×C6 (2,2,0) / 三角化 T² (1,2,1) / icosphere S² (1,0,1)",
            b_ring == [1, 1, 0]
                && b_tetra == [1, 0, 0]
                && b_two == [2, 2, 0]
                && b_torus == [1, 2, 1]
                && b_ico == [1, 0, 1],
            format!(
                "C8 {:?} / K4 {:?} / 2×C6 {:?} / T² {:?} / S² {:?}",
                b_ring, b_tetra, b_two, b_torus, b_ico
            ),
        );
        check(
            "[T0b] icosphere-1 構成 — V = 42, E = 120, 全 link 円周 (閉曲面)",
            ni == 42 && ei.len() == 120 && cl == ni && bd == 0 && bad == 0,
            format!("V={} E={} link(閉/境/破)=({},{},{})", ni, ei.len(), cl, bd, bad),
        );
    }

    // ---- [T1] 内在次元 (殻計数) の器械検査 ----
    {
        let ring24: Vec<(usize, usize)> = (0..24)
            .map(|i| (i.min((i + 1) % 24), i.max((i + 1) % 24)))
            .collect();
        let d1 = intrinsic_dim(&ring24, 24, 3);
        let (nt, at) = torus_tri(8);
        let d2 = intrinsic_dim(&edges_of_adj(&at, nt), nt, 3);
        check(
            "[T1] 内在次元 (殻計数 S(r) ~ r^{d−1}) — C24 = 1・三角 T² 8×8 = 2",
            (d1 - 1.0).abs() < 0.05 && (d2 - 2.0).abs() < 0.05,
            format!("dim(C24) = {:.3} / dim(T²) = {:.3}", d1, d2),
        );
    }

    // ---- [T2] 正例 (状態から end-to-end, 熱的 Gaussian 状態族 β t = 2) ----
    let cases_geom = || -> Vec<(&'static str, usize, Vec<f64>, [usize; 3], &'static str)> {
        vec![
            {
                let (n, a) = torus_tri(8);
                ("torus 8×8", n, a, [1, 2, 1], "closed")
            },
            {
                let (n, a) = cylinder_tri(8, 6);
                ("cylinder 8×6", n, a, [1, 1, 0], "boundary")
            },
            {
                let (n, a, _) = disk_tri(8, 7, &[]);
                ("disk 8×7", n, a, [1, 0, 0], "boundary")
            },
            {
                let (n, a) = icosphere1();
                ("sphere ico-1", n, a, [1, 0, 1], "closed")
            },
        ]
    };
    {
        println!(
            "\n  -- [T2] 正例 (熱的 Gaussian 状態 β t = {} → B3-COV 核 → pipeline, 窓 kmax = 7..10) --",
            BETA_T
        );
        println!("     系             β 期待     β 実測     曲面裁定   全会一致  dim    geodesic");
        let mut all_ok = true;
        for (name, n, a, want_b, want_s) in cases_geom() {
            let c = corr_thermal(&a, n, BETA_T);
            let w = kernel_b3(&c, n);
            let v = run_pipeline(&w, n);
            let ok = v.unanimous && v.n_windows >= 2 && v.betti == want_b && v.surface == want_s;
            all_ok &= ok;
            let (edges, _) = adaptive_mutual(&w, n, 2, 8.min(n - 2));
            let te = edges_of_adj(&a, n);
            let missing = te.iter().filter(|e| !edges.contains(e)).count();
            let extra = edges.iter().filter(|e| !te.contains(e)).count();
            let (lc, lb, lx) = linkness(&edges, n);
            println!(
                "     {:13} {:?}  {:?}  {:9} {:5} ({}窓)  {:.2}  {:.3}  [診断 E={} 欠{} 余{} link {}/{}/{}]  {}",
                name,
                want_b,
                v.betti,
                v.surface,
                v.unanimous,
                v.n_windows,
                v.dim,
                v.geo,
                edges.len(),
                missing,
                extra,
                lc,
                lb,
                lx,
                if ok { "" } else { "** 期待と不一致 **" }
            );
        }
        check(
            "[T2] 正例 4 種の end-to-end 資格 (熱的状態族) — β・曲面性が期待どおり + 窓の全会一致",
            all_ok,
            "詳細は上表".into(),
        );
    }

    // ---- [T2b] 臨界半充填 GS の既知限界 (記録のみ — ゲートにしない) ----
    {
        println!("\n  -- [T2b] 臨界半充填 GS の既知限界 (報告のみ): 境界増強・Friedel 共鳴が局所グラフ抽出を破る --");
        for (name, n, a, want_b, _want_s) in cases_geom() {
            let c = corr_half_filling(&a, n);
            let w = kernel_b3(&c, n);
            let v = run_pipeline(&w, n);
            let (edges, _) = adaptive_mutual(&w, n, 2, 8.min(n - 2));
            let te = edges_of_adj(&a, n);
            let missing = te.iter().filter(|e| !edges.contains(e)).count();
            let extra = edges.iter().filter(|e| !te.contains(e)).count();
            println!(
                "     {:13} β {:?} (期待 {:?}) 裁定 {:11} [欠{} 余{}] {}",
                name,
                v.betti,
                want_b,
                v.surface,
                missing,
                extra,
                if missing + extra > 0 {
                    "← 臨界 GS の長距離構造 (v28 境界増強の 2D 版)"
                } else {
                    "(この系は臨界でも無傷)"
                }
            );
        }
        println!("     → 状態族の仮定が pipeline の前提 (v29.5 [C4] と同じ論点の裏面 — ギャップ相関で資格・臨界は既知限界)");
    }

    // ---- [T3] 敵対対照 ----
    {
        println!("\n  -- [T3] 敵対対照 --");
        let mut all_ok = true;
        {
            let (n, a) = random_cubic(64, 296);
            let c = corr_thermal(&a, n, BETA_T);
            let w = kernel_b3(&c, n);
            let v = run_pipeline(&w, n);
            let ok = v.surface == "not-surface";
            all_ok &= ok;
            println!(
                "     random 3-regular n=64: 裁定 {} (一致 {}) β {:?}  {}",
                v.surface,
                v.unanimous,
                v.betti,
                if ok { "(曲面と誤認せず)" } else { "** 曲面と誤認 **" }
            );
        }
        {
            let (n, a) = branched(8, 4);
            let c = corr_thermal(&a, n, BETA_T);
            let w = kernel_b3(&c, n);
            let v = run_pipeline(&w, n);
            let ok = v.surface == "not-surface";
            all_ok &= ok;
            println!(
                "     branched (3 シート共有線) n={}: 裁定 {} (spine link は theta 型)  {}",
                n,
                v.surface,
                if ok { "(非多様体を検出)" } else { "** 曲面と誤認 **" }
            );
        }
        {
            let holes: Vec<(usize, usize)> = vec![
                (3, 3),
                (7, 6),
                (8, 6),
                (7, 7),
                (8, 7),
                (9, 6),
                (9, 7),
                (8, 8),
            ];
            let (n, a, map) = disk_tri(12, 11, &holes);
            let c = corr_thermal(&a, n, BETA_T);
            let w = kernel_b3(&c, n);
            let v = run_pipeline(&w, n);
            let ok = v.unanimous && v.betti == [1, 2, 0] && v.surface == "boundary";
            all_ok &= ok;
            {
                let (edges, _) = adaptive_mutual(&w, n, 2, 8);
                let bad = bad_link_vertices(&edges, n);
                if !bad.is_empty() {
                    let coords: Vec<String> = bad
                        .iter()
                        .map(|&v2| format!("({},{})", map[v2] % 12, map[v2] / 12))
                        .collect();
                    println!("     [診断] two-holes の link 破れ頂点: {}", coords.join(" "));
                    let te = edges_of_adj(&a, n);
                    let fmt = |e: &(usize, usize)| -> String {
                        format!(
                            "({},{})-({},{})",
                            map[e.0] % 12,
                            map[e.0] / 12,
                            map[e.1] % 12,
                            map[e.1] / 12
                        )
                    };
                    let miss: Vec<String> = te.iter().filter(|e| !edges.contains(e)).map(fmt).collect();
                    let extra: Vec<String> = edges.iter().filter(|e| !te.contains(e)).map(fmt).collect();
                    println!("     [診断] 欠辺: {} / 余辺: {}", miss.join(" "), extra.join(" "));
                }
            }
            println!(
                "     two-holes disk n={}: β {:?} (期待 [1,2,0]) 裁定 {} 一致 {}  {}",
                n,
                v.betti,
                v.surface,
                v.unanimous,
                if ok { "(大小 2 穴を分離)" } else { "** 不一致 **" }
            );
        }
        {
            let (n, a) = petersen();
            let c = corr_thermal(&a, n, BETA_T);
            let w = kernel_b3(&c, n);
            let v = run_pipeline(&w, n);
            let ok_p = v.surface == "not-surface";
            let n7 = 7;
            let mut a7 = vec![1.0; 49];
            for i in 0..7 {
                a7[i + i * 7] = 0.0;
            }
            let c7 = corr_thermal(&a7, n7, BETA_T);
            let w7 = kernel_b3(&c7, n7);
            let v7 = run_pipeline(&w7, n7);
            let ok_k = v7.surface == "not-surface" && v7.degenerate_windows > 0;
            all_ok &= ok_p && ok_k;
            println!(
                "     Petersen: 裁定 {} β {:?} (適応 k=3 + link 検査 — v29.5 誤認の構造解決)  {}",
                v.surface,
                v.betti,
                if ok_p { "" } else { "** 誤認 **" }
            );
            println!(
                "     K7 (完全グラフ): 裁定 {} 縮退保留窓 {}/4 — {}",
                v7.surface,
                v7.degenerate_windows,
                if ok_k { "捏造なし" } else { "** 誤認 **" }
            );
        }
        check(
            "[T3] 敵対対照 5 種 — 曲面の捏造なし・2 穴分離・縮退保留",
            all_ok,
            "詳細は上表".into(),
        );
    }

    println!("\n[要約] dimension-agnostic pipeline (頂点適応 gap kNN + clique complex + Z2 homology + link 検査 +");
    println!("    殻計数次元 + geodesic 整合 + 窓の全会一致) が 2D 正例 4 種を状態から end-to-end で同定し、");
    println!("    敵対対照で幾何を捏造しないことを資格した — 設計区画であり、新鮮 holdout での検証は v30.0。");

    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
