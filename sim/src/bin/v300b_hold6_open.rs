//! v30.0-B HOLD-6 の開封と本採点 — 第三十期の最終検証 (PROMPT/11)
//!
//! hold-0..7 (8 系) をここで**初めて生成・観測**する。凍結済み: 生成器・採点
//! pipeline・離散バー (FROZEN TOPO v32 = v30.0-A コミット 8024f17 — [H0] が節の
//! SHA-256 逐語一致を検査)。SECRET を本ファイルで開示 — sha256(SECRET) が
//! v30.0-A 公表のコミットメント fe3c9cbd… と一致し、train seed が凍結定数と
//! 一致することを [H0] が機械検査する (生成系列全体が第三者検証可能)。
//!
//! **バー外の結果は調整せず受理する** — 成績は [採点] 行、[FAIL] は器械破れ
//! (凍結節不一致・採点不能) のみ (v28.2 以来の慣例)。

use uft_sim::*;

// ================== FROZEN TOPO v32 (BEGIN) ==================
// 2D topology pipeline (v296 逐語コピー) + HOLD-6 生成器とバー。v30.0-A の
// コミットで凍結 — v30.0-B が本節の SHA-256 一致を検査する。節外での再定義禁止。

/// sha256(SECRET) — SECRET は v30.0-B で開示 (それまでリポジトリに存在しない)
pub const HOLD6_COMMITMENT: &str = "fe3c9cbd0c2d733f852422734ca4f212cd01e488fac763142ceef07d598d6b62";

/// train instance seed (= sha256(SECRET+":"+id) 先頭 8 バイト — 開示時に第三者検証)
pub const HOLD6_TRAIN_SEEDS: [(u64, &str); 5] = [
    (13734079978951196301, "train-torus"),
    (5520125636343895927, "train-cylinder"),
    (4702528085787679178, "train-disk"),
    (13593300980156879169, "train-holes"),
    (16551578033132347782, "train-sphere"),
];

// ---- 採点 pipeline (v296_homology.rs の逐語コピー) ----

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


// ---- HOLD-6 生成器 ----

/// HOLD-6 instance (seed から一意)。クラス: 0 torus / 1 cylinder / 2 disk /
/// 3 two-holes disk / 4 sphere (icosphere-1)。train はクラスを強制、holdout は
/// seed が決める。期待プロファイル (Betti, 曲面性) が離散バー。
pub struct Hold6Instance {
    pub class: usize,
    pub n: usize,
    pub adj: Vec<f64>,
    pub expect_betti: [usize; 3],
    pub expect_surface: &'static str,
}

/// ガウス束場 v(p) = 1 + Σ A_k exp(−d(p, c_k)²/w_k²)。d は per-axis min-image
/// (周期軸のみ)。全評価点で v ∈ [0.65, 1.45] になるまで棄却 (決定的)。
struct Field2 {
    amps: Vec<f64>,
    cents: Vec<[f64; 2]>,
    wids: Vec<f64>,
    period: [Option<f64>; 2],
}

impl Field2 {
    fn eval(&self, p: [f64; 2]) -> f64 {
        let mut v = 1.0;
        for i in 0..self.amps.len() {
            let mut d2 = 0.0;
            for ax in 0..2 {
                let mut dx = p[ax] - self.cents[i][ax];
                if let Some(l) = self.period[ax] {
                    dx -= l * (dx / l).round();
                }
                d2 += dx * dx;
            }
            v += self.amps[i] * (-d2 / (self.wids[i] * self.wids[i])).exp();
        }
        v
    }
}

fn sample_field2(rng: &mut Rng, lx: f64, ly: f64, period: [Option<f64>; 2]) -> Field2 {
    loop {
        let k = 2 + rng.range(3);
        let mut f = Field2 {
            amps: Vec::new(),
            cents: Vec::new(),
            wids: Vec::new(),
            period,
        };
        for _ in 0..k {
            f.amps.push(-0.2 + 0.45 * rng.f64());
            f.cents.push([lx * rng.f64(), ly * rng.f64()]);
            f.wids.push(2.0 + 3.0 * rng.f64());
        }
        let steps = 40usize;
        let ok = (0..=steps).all(|ix| {
            (0..=steps).all(|iy| {
                let v = f.eval([lx * ix as f64 / steps as f64, ly * iy as f64 / steps as f64]);
                (0.65..=1.45).contains(&v)
            })
        });
        if ok {
            return f;
        }
    }
}

/// 格子クラスの座標つき構成: alive な (x, y) の集合に近傍則
/// (+1,0), (0,+1), (+1,+1) [周期軸は wrap] で辺を張り、場の中点値を重みにする。
fn grid_weighted(
    l: usize,
    w: usize,
    px: bool,
    py: bool,
    holes: &[(usize, usize)],
    field: &Field2,
) -> (usize, Vec<f64>) {
    let n0 = l * w;
    let mut alive = vec![true; n0];
    for &(hx, hy) in holes {
        alive[hx + hy * l] = false;
    }
    let map: Vec<usize> = (0..n0).filter(|&i| alive[i]).collect();
    let inv: std::collections::BTreeMap<usize, usize> =
        map.iter().enumerate().map(|(k, &v)| (v, k)).collect();
    let n = map.len();
    let mut adj = vec![0.0; n * n];
    for &site in &map {
        let (x, y) = (site % l, site / l);
        for (dx, dy) in [(1usize, 0usize), (0, 1), (1, 1)] {
            let (nx, ny) = (x + dx, y + dy);
            let (tx, ty) = (
                if px { nx % l } else { nx },
                if py { ny % w } else { ny },
            );
            if tx >= l || ty >= w {
                continue;
            }
            let t = tx + ty * l;
            if !alive[t] {
                continue;
            }
            let (i, j) = (inv[&site], inv[&t]);
            // 中点 (周期軸は wrap を跨ぐ辺で min-image に一致するよう非 wrap 座標で評価)
            let v = field.eval([x as f64 + dx as f64 / 2.0, y as f64 + dy as f64 / 2.0]);
            adj[i + j * n] = v;
            adj[j + i * n] = v;
        }
    }
    (n, adj)
}

/// icosphere-1 の頂点座標 (v296 icosphere1 と同一の構成順)
pub fn icosphere1_coords() -> Vec<[f64; 3]> {
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
    for &(i, j) in &edges0 {
        let m = norm([
            (verts[i][0] + verts[j][0]) / 2.0,
            (verts[i][1] + verts[j][1]) / 2.0,
            (verts[i][2] + verts[j][2]) / 2.0,
        ]);
        verts.push(m);
    }
    verts
}

/// icosphere-1 の隣接 (v296 icosphere1 と同一のグラフ) — 座標から再構成
fn icosphere1_adj_from_coords(verts: &[[f64; 3]]) -> (usize, Vec<f64>) {
    let d2 = |a: [f64; 3], b: [f64; 3]| -> f64 {
        (a[0] - b[0]).powi(2) + (a[1] - b[1]).powi(2) + (a[2] - b[2]).powi(2)
    };
    // 細分後の最近接距離: 元辺の半分 (弦² ≈ 0.276) と面内の中点間 (弦² ≈ 0.305)
    // — 第 2 近接 (弦² ≈ 0.55) との間に明確なギャップがある閾値 0.4 で辺を張る
    let n = verts.len();
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            if d2(verts[i], verts[j]) < 0.4 {
                a[i + j * n] = 1.0;
                a[j + i * n] = 1.0;
            }
        }
    }
    (n, a)
}

pub fn hold6_instance(seed: u64, forced_class: Option<usize>) -> Hold6Instance {
    let mut rng = Rng::new(seed);
    let class = forced_class.unwrap_or_else(|| rng.range(5));
    let (n, mut adj, expect_betti, expect_surface): (usize, Vec<f64>, [usize; 3], &'static str) =
        match class {
            0 => {
                let l = 7 + rng.range(3);
                let f = sample_field2(&mut rng, l as f64, l as f64, [Some(l as f64), Some(l as f64)]);
                let (n, a) = grid_weighted(l, l, true, true, &[], &f);
                (n, a, [1, 2, 1], "closed")
            }
            1 => {
                let l = 7 + rng.range(3);
                let w = 5 + rng.range(3);
                let f = sample_field2(&mut rng, l as f64, w as f64, [Some(l as f64), None]);
                let (n, a) = grid_weighted(l, w, true, false, &[], &f);
                (n, a, [1, 1, 0], "boundary")
            }
            2 => {
                let l = 7 + rng.range(3);
                let w = 5 + rng.range(3);
                let f = sample_field2(&mut rng, l as f64, w as f64, [None, None]);
                let (n, a) = grid_weighted(l, w, false, false, &[], &f);
                (n, a, [1, 0, 0], "boundary")
            }
            3 => {
                // 大穴 = **凸ブロック** a×b (三角格子の隣接に対し格子直線との交差が
                // 連続 = 各境界頂点の失う近傍が連続 → link は単一の弧)。貪欲成長の
                // 凹クラスタは真の複体にピンチ点 (非多様体頂点) を作った (設計走行
                // の発見 — 開発記録)。
                let l = 11 + rng.range(3);
                let w = 10 + rng.range(3);
                let blocks: [(usize, usize); 5] = [(2, 2), (2, 3), (3, 2), (1, 4), (4, 1)];
                let (ba, bb) = blocks[rng.range(5)];
                loop {
                    let sx = 2 + rng.range(l - 3 - ba);
                    let sy = 2 + rng.range(w - 3 - bb);
                    let hx = 2 + rng.range(l - 4);
                    let hy = 2 + rng.range(w - 4);
                    // 小穴はブロックから Manhattan ≥ 4
                    let mut far = true;
                    let mut holes: Vec<(usize, usize)> = Vec::new();
                    for dx in 0..ba {
                        for dy in 0..bb {
                            let (cx, cy) = (sx + dx, sy + dy);
                            holes.push((cx, cy));
                            if (cx as isize - hx as isize).abs() + (cy as isize - hy as isize).abs() < 4 {
                                far = false;
                            }
                        }
                    }
                    if !far {
                        continue;
                    }
                    holes.push((hx, hy));
                    let f = sample_field2(&mut rng, l as f64, w as f64, [None, None]);
                    let (n, a) = grid_weighted(l, w, false, false, &holes, &f);
                    break (n, a, [1, 2, 0], "boundary");
                }
            }
            _ => {
                let coords = icosphere1_coords();
                let (n, mut a) = icosphere1_adj_from_coords(&coords);
                loop {
                    let k = 2 + rng.range(3);
                    let mut amps = Vec::new();
                    let mut cents: Vec<[f64; 3]> = Vec::new();
                    let mut wids = Vec::new();
                    for _ in 0..k {
                        amps.push(-0.2 + 0.45 * rng.f64());
                        let g = [rng.gauss(), rng.gauss(), rng.gauss()];
                        let r = (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]).sqrt().max(1e-9);
                        cents.push([g[0] / r, g[1] / r, g[2] / r]);
                        wids.push(0.6 + 0.6 * rng.f64());
                    }
                    let eval = |p: [f64; 3]| -> f64 {
                        let mut v = 1.0;
                        for i in 0..amps.len() {
                            let dot = (p[0] * cents[i][0] + p[1] * cents[i][1] + p[2] * cents[i][2])
                                .clamp(-1.0, 1.0);
                            let ang = dot.acos();
                            v += amps[i] * (-(ang * ang) / (wids[i] * wids[i])).exp();
                        }
                        v
                    };
                    if !coords.iter().all(|&p| (0.65..=1.45).contains(&eval(p))) {
                        continue;
                    }
                    for i in 0..n {
                        for j in (i + 1)..n {
                            if a[i + j * n] != 0.0 {
                                let m = [
                                    (coords[i][0] + coords[j][0]) / 2.0,
                                    (coords[i][1] + coords[j][1]) / 2.0,
                                    (coords[i][2] + coords[j][2]) / 2.0,
                                ];
                                let r = (m[0] * m[0] + m[1] * m[1] + m[2] * m[2]).sqrt().max(1e-9);
                                let v = eval([m[0] / r, m[1] / r, m[2] / r]);
                                a[i + j * n] = v;
                                a[j + i * n] = v;
                            }
                        }
                    }
                    break;
                }
                (n, a, [1, 0, 1], "closed")
            }
        };
    // ノード置換 (ラベル隠蔽 — pipeline は行列しか見ないが、tie-break の
    // 構成順依存を隠す)
    {
        let mut perm: Vec<usize> = (0..n).collect();
        for i in (1..n).rev() {
            let j = rng.range(i + 1);
            perm.swap(i, j);
        }
        let mut b = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                b[i + j * n] = adj[perm[i] + perm[j] * n];
            }
        }
        adj = b;
    }
    Hold6Instance {
        class,
        n,
        adj,
        expect_betti,
        expect_surface,
    }
}

/// 採点 (離散バー): 熱的状態 (βt = 1) の B3-COV 核 → pipeline →
/// (Betti, 曲面性) がクラス期待と厳密一致 + 全会一致 + 採用窓 ≥ 2
pub fn score_hold6(inst: &Hold6Instance) -> (bool, Verdict) {
    let c = corr_thermal(&inst.adj, inst.n, BETA_T);
    let w = kernel_b3(&c, inst.n);
    let v = run_pipeline(&w, inst.n);
    let ok = v.unanimous
        && v.n_windows >= 2
        && v.betti == inst.expect_betti
        && v.surface == inst.expect_surface;
    (ok, v)
}

// =================== FROZEN TOPO v32 (END) ===================

// ---------------- 実験側 (開封 — 凍結節の外) ----------------

/// SECRET (開示 — v30.0-A コミットメント fe3c9cbd… の原像)
const HOLD6_SECRET: &str = "HOLD6-b4b84a54679589ad2773200a3251c43d";

fn hold6_seed(id: &str) -> u64 {
    let s = format!("{}:{}", HOLD6_SECRET, id);
    let h = sha256(s.as_bytes());
    u64::from_be_bytes([h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]])
}

fn main() {
    self_test();
    println!("=== v30.0-B HOLD-6 の開封と本採点 (第三十期 最終検証) ===\n");
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

    // ---- [H0] 凍結の完全性 ----
    {
        let commit = sha256_hex(HOLD6_SECRET.as_bytes());
        let ok_commit = commit == HOLD6_COMMITMENT;
        let ok_train = HOLD6_TRAIN_SEEDS
            .iter()
            .all(|&(seed, id)| hold6_seed(id) == seed);
        let read_src = |p: &str| -> String {
            std::fs::read_to_string(p)
                .or_else(|_| std::fs::read_to_string(format!("../{}", p)))
                .unwrap_or_default()
        };
        let cut = |src: &str| -> Option<String> {
            let b = src.find("// ================== FROZEN TOPO v32 (BEGIN)")?;
            let e = src.find("// =================== FROZEN TOPO v32 (END)")?;
            Some(src[b..e].to_string())
        };
        let me = read_src("sim/src/bin/v300b_hold6_open.rs");
        let sa = read_src("sim/src/bin/v300a_hold6_freeze.rs");
        let (ok32, sha32) = match (cut(&sa), cut(&me)) {
            (Some(a), Some(m)) => (
                !a.is_empty() && sha256_hex(a.as_bytes()) == sha256_hex(m.as_bytes()),
                sha256_hex(a.as_bytes()),
            ),
            _ => (false, "?".into()),
        };
        check(
            "[H0] 凍結の完全性 — SECRET コミットメント・train seed・FROZEN TOPO v32 の逐語一致",
            ok_commit && ok_train && ok32,
            format!(
                "commitment {} / train {} / v32 {} ({}…)",
                ok_commit,
                ok_train,
                ok32,
                &sha32[..16.min(sha32.len())]
            ),
        );
    }

    // ---- 開封: hold-0..7 の seed を初導出 ----
    let hold_ids: [&str; 8] = [
        "hold-0", "hold-1", "hold-2", "hold-3", "hold-4", "hold-5", "hold-6", "hold-7",
    ];
    println!("\n  [開封] hold seed (本走行で初導出):");
    for id in hold_ids {
        println!("        {} = {}", id, hold6_seed(id));
    }

    // ---- [H1] 本採点 (離散バー — 成績は [採点], 器械破れのみ [FAIL]) ----
    let mut n_out = 0usize;
    let mut scoreable = true;
    {
        println!("\n  -- [H1] holdout 8 系 (クラスは seed が決定, 重み場 + 置換つき) --");
        println!("     instance  class      n    β 実測     曲面裁定   一致  窓   裁定");
        for id in hold_ids {
            let seed = hold6_seed(id);
            let inst = hold6_instance(seed, None);
            let (ok, v) = score_hold6(&inst);
            if v.n_windows == 0 {
                scoreable = false;
            }
            if !ok {
                n_out += 1;
            }
            println!(
                "     {:8} {:9} {:4}  {:?}  {:9} {:5}  {}   {}",
                id,
                ["torus", "cylinder", "disk", "two-holes", "sphere"][inst.class],
                inst.n,
                v.betti,
                v.surface,
                v.unanimous,
                v.n_windows,
                if ok { "バー内" } else { "**バー外**" }
            );
        }
        println!(
            "\n  [採点] [H1] holdout 8 系: バー外 {} 件 (離散バー = (Betti, 曲面性) 厳密一致 + 全会一致 + 窓 ≥ 2)",
            n_out
        );
    }
    check(
        "[H2] 採点可能性 — 全 8 系で pipeline が裁定を返した (全窓縮退保留なし)",
        scoreable,
        format!("採点不能 0 件想定 / バー外 {} 件は [採点] 行 (器械破れではない)", n_out),
    );

    println!(
        "\n[裁定] HOLD-6 holdout (事前凍結の離散バーでの本採点 — 調整なし): {}",
        if n_out == 0 {
            "全 8 系バー内 — dimension-agnostic pipeline が新鮮な重み場つき 2D holdout で生存"
        } else {
            "バー外あり — 該当系の棄却/限界を確定表に記録 (調整・再採点はしない)"
        }
    );

    {
        let j = Json::Obj(vec![
            ("version".into(), Json::Str("v30.0-B".into())),
            ("n_out".into(), Json::Num(n_out as f64)),
            ("all_in".into(), Json::Bool(n_out == 0)),
        ]);
        let p = write_artifact("results/v300b_hold6_open.json", &j.render());
        println!("\n[成果物] {}", p);
    }

    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
