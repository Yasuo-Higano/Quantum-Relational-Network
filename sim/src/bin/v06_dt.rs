//! v0.6b 2次元動的三角形分割 (DT) — 時空自体を量子的に揺らがせる
//!
//! 2D ユークリッド量子重力: 分配関数 Z = Σ_{S²の三角形分割} e^{-λN₂}。
//! 面積 (N₂) を固定すると Einstein-Hilbert 作用は位相不変量になるため、
//! **全ての三角形分割が等確率** = 純粋 2D 量子重力のアンサンブル。
//! エルゴード的な (2,2) フリップ移動でサンプリングする。
//!
//! 理論の予言 (KPZ / Liouville 重力):
//!   ハウスドルフ次元 d_H = 4  (2次元なのに! 球の体積が r⁴ で増える激しいフラクタル)
//!   スペクトル次元   d_s = 2  (拡散で測ると 2 次元)
//! 平坦な三角格子 (古典幾何) は d_H = d_s = 2。量子幾何はまるで違う。

use std::collections::HashSet;
use uft_sim::*;

struct Dt {
    tri: Vec<[u32; 3]>,  // 頂点
    nbr: Vec<[u32; 3]>,  // nbr[t][i] = 頂点 i の対辺の向こうの三角形
    deg: Vec<u32>,       // 頂点次数
    edges: HashSet<u64>, // 無向辺
    rng: Rng,
}

fn ekey(a: u32, b: u32) -> u64 {
    let (lo, hi) = if a < b { (a, b) } else { (b, a) };
    (hi as u64) << 32 | lo as u64
}

impl Dt {
    /// 正四面体 (S² の最小分割) から開始
    fn tetrahedron(seed: u64) -> Self {
        let tri = vec![[0, 1, 2], [0, 3, 1], [0, 2, 3], [1, 3, 2]];
        let mut dt = Dt {
            tri,
            nbr: vec![[0; 3]; 4],
            deg: vec![3; 4],
            edges: HashSet::new(),
            rng: Rng::new(seed),
        };
        dt.rebuild_neighbors();
        dt
    }
    fn rebuild_neighbors(&mut self) {
        use std::collections::HashMap;
        let mut map: HashMap<u64, Vec<(u32, u32)>> = HashMap::new();
        self.edges.clear();
        for (t, tv) in self.tri.iter().enumerate() {
            for i in 0..3 {
                let (a, b) = (tv[(i + 1) % 3], tv[(i + 2) % 3]);
                map.entry(ekey(a, b)).or_default().push((t as u32, i as u32));
                self.edges.insert(ekey(a, b));
            }
        }
        for (_, v) in map.iter() {
            assert!(v.len() == 2, "non-manifold edge");
            let (t0, i0) = v[0];
            let (t1, i1) = v[1];
            self.nbr[t0 as usize][i0 as usize] = t1;
            self.nbr[t1 as usize][i1 as usize] = t0;
        }
    }
    /// 1→3 移動: 三角形 t に新頂点 w を挿入
    fn insert(&mut self, t: usize) {
        let [v0, v1, v2] = self.tri[t];
        let [n0, n1, n2] = self.nbr[t];
        let w = self.deg.len() as u32;
        self.deg.push(3);
        let t1 = self.tri.len() as u32; // (v1,v2,w)
        let t2 = t1 + 1; // (v2,v0,w)
        let t0 = t as u32;
        self.tri[t] = [v0, v1, w];
        self.nbr[t] = [t1, t2, n2];
        self.tri.push([v1, v2, w]);
        self.nbr.push([t2, t0, n0]);
        self.tri.push([v2, v0, w]);
        self.nbr.push([t0, t1, n1]);
        // 外側の隣接を付け替え (n2 は t のまま)
        self.repoint(n0 as usize, t0, t1);
        self.repoint(n1 as usize, t0, t2);
        self.edges.insert(ekey(v0, w));
        self.edges.insert(ekey(v1, w));
        self.edges.insert(ekey(v2, w));
        self.deg[v0 as usize] += 1;
        self.deg[v1 as usize] += 1;
        self.deg[v2 as usize] += 1;
    }
    fn repoint(&mut self, t: usize, from: u32, to: u32) {
        if t == from as usize && from == to {
            return;
        }
        for i in 0..3 {
            if self.nbr[t][i] == from {
                self.nbr[t][i] = to;
                return;
            }
        }
        panic!("repoint: link not found");
    }
    /// (2,2) フリップ: 三角形 t の頂点スロット i の対辺をフリップ。合法なら true
    fn flip(&mut self, t: usize, i: usize) -> bool {
        let u = self.nbr[t][i] as usize;
        let c = self.tri[t][i];
        let p = self.tri[t][(i + 1) % 3];
        let q = self.tri[t][(i + 2) % 3];
        // u 内で t を指すスロット j (対頂点 d)
        let mut j = 4;
        for k in 0..3 {
            if self.nbr[u][k] == t as u32 {
                j = k;
                break;
            }
        }
        let d = self.tri[u][j];
        if c == d || self.edges.contains(&ekey(c, d)) {
            return false;
        }
        if self.deg[p as usize] <= 3 || self.deg[q as usize] <= 3 {
            return false;
        }
        // u = (d, q, p) の向き: slot j の次が q, その次が p であることを確認
        let (uq, up) = (self.tri[u][(j + 1) % 3], self.tri[u][(j + 2) % 3]);
        debug_assert!(uq == q && up == p, "orientation mismatch");
        // t の外側: A = nbr 対 p (辺 q-c), B = nbr 対 q (辺 c-p)
        let a_ = self.nbr[t][(i + 1) % 3];
        let b_ = self.nbr[t][(i + 2) % 3];
        // u の外側: C = nbr 対 q (辺 p-d), D = nbr 対 p (辺 d-q)
        let c_ = self.nbr[u][(j + 1) % 3];
        let d_ = self.nbr[u][(j + 2) % 3];
        // 新しい三角形: t' = (c,p,d) [nbr: 対c=C, 対p=u, 対d=B], u' = (c,d,q) [nbr: 対c=D, 対d=A, 対q=t]
        self.tri[t] = [c, p, d];
        self.nbr[t] = [c_, u as u32, b_];
        self.tri[u] = [c, d, q];
        self.nbr[u] = [d_, a_, t as u32];
        self.repoint(a_ as usize, t as u32, u as u32);
        self.repoint(c_ as usize, u as u32, t as u32);
        self.edges.remove(&ekey(p, q));
        self.edges.insert(ekey(c, d));
        self.deg[p as usize] -= 1;
        self.deg[q as usize] -= 1;
        self.deg[c as usize] += 1;
        self.deg[d as usize] += 1;
        true
    }
    fn sweep_flips(&mut self) -> f64 {
        let nt = self.tri.len();
        let mut acc = 0;
        for _ in 0..3 * nt {
            let t = self.rng.range(nt);
            let i = self.rng.range(3);
            if self.flip(t, i) {
                acc += 1;
            }
        }
        acc as f64 / (3 * nt) as f64
    }
    fn validate(&self) {
        for t in 0..self.tri.len() {
            for i in 0..3 {
                let u = self.nbr[t][i] as usize;
                let (a, b) = (self.tri[t][(i + 1) % 3], self.tri[t][(i + 2) % 3]);
                let mut found = false;
                for k in 0..3 {
                    if self.nbr[u][k] == t as u32 {
                        let (ua, ub) = (self.tri[u][(k + 1) % 3], self.tri[u][(k + 2) % 3]);
                        assert!(ua == b && ub == a, "edge mismatch");
                        found = true;
                    }
                }
                assert!(found, "mutual neighbor missing");
            }
        }
        // オイラー: V - E + F = 2
        let v = self.deg.len();
        let e = self.edges.len();
        let f = self.tri.len();
        assert!(v + f == e + 2, "Euler check failed");
    }
    /// 頂点隣接リスト (CSR)
    fn adjacency(&self) -> (Vec<u32>, Vec<u32>) {
        let nv = self.deg.len();
        let mut off = vec![0u32; nv + 1];
        for &e in self.edges.iter() {
            let (a, b) = ((e & 0xffffffff) as usize, (e >> 32) as usize);
            off[a + 1] += 1;
            off[b + 1] += 1;
        }
        for i in 0..nv {
            off[i + 1] += off[i];
        }
        let mut adj = vec![0u32; off[nv] as usize];
        let mut cur = off.clone();
        for &e in self.edges.iter() {
            let (a, b) = ((e & 0xffffffff) as u32, (e >> 32) as u32);
            adj[cur[a as usize] as usize] = b;
            cur[a as usize] += 1;
            adj[cur[b as usize] as usize] = a;
            cur[b as usize] += 1;
        }
        (off, adj)
    }
}

/// ランダムウォークの帰還確率 P(t)
fn return_prob(off: &[u32], adj: &[u32], walks: usize, tmax: usize, rng: &mut Rng) -> Vec<f64> {
    let nv = off.len() - 1;
    let mut p = vec![0f64; tmax + 1];
    for _ in 0..walks {
        let start = rng.range(nv) as u32;
        let mut pos = start;
        for t in 1..=tmax {
            let (s, e) = (off[pos as usize], off[pos as usize + 1]);
            pos = adj[(s + (rng.u64() % (e - s) as u64) as u32) as usize];
            if pos == start {
                p[t] += 1.0;
            }
        }
    }
    for t in 0..=tmax {
        p[t] /= walks as f64;
    }
    p
}

/// BFS 距離配列
fn bfs_dist(off: &[u32], adj: &[u32], src: u32) -> Vec<u32> {
    let nv = off.len() - 1;
    let mut dist = vec![u32::MAX; nv];
    let mut queue = std::collections::VecDeque::new();
    dist[src as usize] = 0;
    queue.push_back(src);
    while let Some(v) = queue.pop_front() {
        for k in off[v as usize]..off[v as usize + 1] {
            let w = adj[k as usize];
            if dist[w as usize] == u32::MAX {
                dist[w as usize] = dist[v as usize] + 1;
                queue.push_back(w);
            }
        }
    }
    dist
}

/// (球成長ピーク次元, 平均測地距離, スペクトル次元)
/// d_s は次数一様なウォークグラフ (woff,wadj) で測る (DT では 3-正則の双対グラフ)
fn measure_dims(
    off: &[u32],
    adj: &[u32],
    woff: &[u32],
    wadj: &[u32],
    rng: &mut Rng,
    label: &str,
    rmax: usize,
) -> (f64, f64, f64) {
    let nv = off.len() - 1;
    // d_H(局所): 球成長の実効次元のピーク + 平均測地距離
    let nsrc = 60;
    let mut mean_count = vec![0f64; rmax + 1];
    let mut rbar = 0.0f64;
    for _ in 0..nsrc {
        let dist = bfs_dist(off, adj, rng.range(nv) as u32);
        let mut c = vec![0f64; rmax + 1];
        let mut sumd = 0.0;
        for &d in dist.iter() {
            if (d as usize) <= rmax {
                c[d as usize] += 1.0;
            }
            sumd += d as f64;
        }
        rbar += sumd / nv as f64 / nsrc as f64;
        let mut cum = 0.0;
        for r in 0..=rmax {
            cum += c[r];
            mean_count[r] += cum / nsrc as f64;
        }
    }
    let mut dh_peak = 0.0f64;
    for r in 2..rmax {
        if mean_count[r + 1] < nv as f64 * 0.45 {
            let d_eff = (mean_count[r + 1].ln() - mean_count[r - 1].ln())
                / (((r + 1) as f64).ln() - ((r - 1) as f64).ln());
            dh_peak = dh_peak.max(d_eff);
        }
    }
    // d_s: 帰還確率 P(t) ~ t^{-d_s/2}, 有限体積フロア 1/V を差し引き、2つの時間窓で測る
    let vw = (woff.len() - 1) as f64;
    let p = return_prob(woff, wadj, 400_000, 240, rng);
    let ds_win = |t0: usize, t1: usize| -> f64 {
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        let mut t = t0;
        while t + 1 <= t1.min(p.len() - 2) {
            let avg = 0.5 * (p[t] + p[t + 1]) - 1.0 / vw;
            if avg > 0.0 {
                xs.push((t as f64 + 0.5).ln());
                ys.push(avg.ln());
            }
            t += (t0 / 4).max(2);
        }
        let (_, slope) = linfit(&xs, &ys);
        -2.0 * slope
    };
    let ds_early = ds_win(16, 48);
    let ds_late = ds_win(80, 240);
    println!(
        "  {} (V={}): d_H(局所ピーク)={:.2}, ⟨r⟩={:.2}, d_s[t=16-48]={:.2}, d_s[t=80-240]={:.2}",
        label, nv, dh_peak, rbar, ds_early, ds_late
    );
    (dh_peak, rbar, ds_late)
}

/// BFS 幾何測定のみ (d_H 局所ピークと平均距離)
fn measure_geom(off: &[u32], adj: &[u32], rng: &mut Rng, rmax: usize, nsrc: usize) -> (f64, f64) {
    let nv = off.len() - 1;
    let mut mean_count = vec![0f64; rmax + 1];
    let mut rbar = 0.0f64;
    for _ in 0..nsrc {
        let dist = bfs_dist(off, adj, rng.range(nv) as u32);
        let mut c = vec![0f64; rmax + 1];
        let mut sumd = 0.0;
        for &d in dist.iter() {
            if (d as usize) <= rmax {
                c[d as usize] += 1.0;
            }
            sumd += d as f64;
        }
        rbar += sumd / nv as f64 / nsrc as f64;
        let mut cum = 0.0;
        for r in 0..=rmax {
            cum += c[r];
            mean_count[r] += cum / nsrc as f64;
        }
    }
    let mut dh_peak = 0.0f64;
    for r in 2..rmax {
        if mean_count[r + 1] < nv as f64 * 0.45 {
            let d_eff = (mean_count[r + 1].ln() - mean_count[r - 1].ln())
                / (((r + 1) as f64).ln() - ((r - 1) as f64).ln());
            dh_peak = dh_peak.max(d_eff);
        }
    }
    (dh_peak, rbar)
}

fn main() {
    println!("=== v0.6b 動的三角形分割: 揺らぐ時空のフラクタル幾何 ===\n");
    let mut rng = Rng::new(777);

    // ---- 基準: 平坦な三角格子トーラス (古典幾何) ----
    println!("[A] 基準系: 平坦な三角格子トーラス 64×64 (古典的な 2 次元)");
    {
        let w = 64usize;
        let nv = w * w;
        let mut off = vec![0u32; nv + 1];
        let mut adjv: Vec<Vec<u32>> = vec![Vec::new(); nv];
        for y in 0..w {
            for x in 0..w {
                let i = x + y * w;
                let nb = [
                    ((x + 1) % w) + y * w,
                    ((x + w - 1) % w) + y * w,
                    x + ((y + 1) % w) * w,
                    x + ((y + w - 1) % w) * w,
                    ((x + 1) % w) + ((y + 1) % w) * w,
                    ((x + w - 1) % w) + ((y + w - 1) % w) * w,
                ];
                for &n in &nb {
                    adjv[i].push(n as u32);
                }
            }
        }
        let mut adj = Vec::new();
        for i in 0..nv {
            off[i + 1] = off[i] + adjv[i].len() as u32;
            adj.extend(&adjv[i]);
        }
        measure_dims(&off, &adj, &off, &adj, &mut rng, "平坦トーラス", 24);
        println!("    (理論値: d_H = d_s = 2 — 古典幾何では両者は一致する)\n");
    }

    // ---- DT アンサンブル ----
    println!("[B] 量子重力アンサンブル (S², 一様測度, フリップでサンプリング)");
    println!("    d_s は 3-正則の双対グラフで測定 (次数の乱れを排除)");
    let mut size_data: Vec<(f64, f64)> = Vec::new(); // (V, ⟨r⟩)
    for &n2 in &[8000usize, 32000, 128000] {
        let mut dt = Dt::tetrahedron(1234 + n2 as u64);
        // 成長: 挿入とフリップを混ぜて目標サイズへ
        while dt.tri.len() < n2 {
            let t = dt.rng.range(dt.tri.len());
            dt.insert(t);
            for _ in 0..6 {
                let t = dt.rng.range(dt.tri.len());
                let i = dt.rng.range(3);
                dt.flip(t, i);
            }
        }
        dt.validate();
        // 熱化
        let mut acc = 0.0;
        for _ in 0..250 {
            acc = dt.sweep_flips();
        }
        dt.validate();
        println!("  N₂={} (熱化後 flip 受理率 {:.2}) — オイラー検査/整合性 OK", n2, acc);
        // 測定 (配位を替えて平均): 幾何は 12 配位、d_s は最初の 2 配位
        let mut dhs = Vec::new();
        let mut dss = Vec::new();
        let mut rbars = Vec::new();
        let rmax = if n2 >= 100_000 { 48 } else { 32 };
        for m in 0..12 {
            for _ in 0..30 {
                dt.sweep_flips();
            }
            let (off, adj) = dt.adjacency();
            if m < 2 {
                // 双対グラフ (3-正則): 三角形の隣接そのもの
                let nt = dt.tri.len();
                let mut woff = vec![0u32; nt + 1];
                let mut wadj = vec![0u32; 3 * nt];
                for t in 0..nt {
                    woff[t + 1] = 3 * (t as u32 + 1);
                    for i in 0..3 {
                        wadj[3 * t + i] = dt.nbr[t][i];
                    }
                }
                let (dh, rbar, ds) = measure_dims(
                    &off,
                    &adj,
                    &woff,
                    &wadj,
                    &mut rng,
                    &format!("配位{}", m + 1),
                    rmax,
                );
                dhs.push(dh);
                dss.push(ds);
                rbars.push(rbar);
            } else {
                let (dh, rbar) = measure_geom(&off, &adj, &mut rng, rmax, 24);
                dhs.push(dh);
                rbars.push(rbar);
            }
        }
        let dh_m: f64 = dhs.iter().sum::<f64>() / dhs.len() as f64;
        let ds_m: f64 = dss.iter().sum::<f64>() / dss.len() as f64;
        let rb_m: f64 = rbars.iter().sum::<f64>() / rbars.len() as f64;
        let rb_sd: f64 = (rbars.iter().map(|x| (x - rb_m).powi(2)).sum::<f64>()
            / (rbars.len() - 1) as f64)
            .sqrt();
        size_data.push(((n2 / 2 + 2) as f64, rb_m));
        println!(
            "  == N₂={}: d_H(局所)={:.2}, ⟨r⟩={:.2}±{:.2} (配位間ゆらぎ!), d_s(遅)≈{:.2} (理論値 2) ==",
            n2, dh_m, rb_m, rb_sd, ds_m
        );
        if n2 >= 100_000 {
            // 長時間ランダムウォークで d_s の 2 への接近を見る
            let nt = dt.tri.len();
            let mut woff = vec![0u32; nt + 1];
            let mut wadj = vec![0u32; 3 * nt];
            for t in 0..nt {
                woff[t + 1] = 3 * (t as u32 + 1);
                for i in 0..3 {
                    wadj[3 * t + i] = dt.nbr[t][i];
                }
            }
            let p = return_prob(&woff, &wadj, 600_000, 1024, &mut rng);
            let fitw = |t0: usize, t1: usize| -> f64 {
                let (mut xs, mut ys) = (Vec::new(), Vec::new());
                let mut t = t0;
                while t + 1 <= t1 {
                    let avg = 0.5 * (p[t] + p[t + 1]) - 1.0 / nt as f64;
                    if avg > 0.0 {
                        xs.push((t as f64 + 0.5).ln());
                        ys.push(avg.ln());
                    }
                    t += t0 / 8;
                }
                let (_, s) = linfit(&xs, &ys);
                -2.0 * s
            };
            println!(
                "  [長時間拡散] d_s: t∈[32,128]:{:.2} → t∈[128,512]:{:.2} → t∈[512,1023]:{:.2}",
                fitw(32, 128),
                fitw(128, 512),
                fitw(512, 1023)
            );
            // Liouville 重力の厳密形 P(t) ~ ln(t)/t (d_s=2 + 対数補正) の検証
            // 対数補正があると局所傾きは d_s(t) = 2 - 2/ln t (t=1000 でも 1.71!)
            println!("  [対数補正の検証] P~ln(t)/t なら P·t/ln t は一定、P·t は ln t で成長:");
            print!("    t:          ");
            for &t in &[64usize, 128, 256, 512, 1000] {
                print!("{:8}", t);
            }
            print!("\n    P·t/ln t:   ");
            for &t in &[64usize, 128, 256, 512, 1000] {
                let pt = 0.5 * (p[t] + p[t + 1]) - 1.0 / nt as f64;
                print!("{:8.4}", pt * t as f64 / (t as f64).ln());
            }
            print!("\n    P·t (単純): ");
            for &t in &[64usize, 128, 256, 512, 1000] {
                let pt = 0.5 * (p[t] + p[t + 1]) - 1.0 / nt as f64;
                print!("{:8.4}", pt * t as f64);
            }
            println!("\n    d_s(t)=2-2/ln t の予言: t=64:{:.2}, t=256:{:.2}, t=1000:{:.2}",
                2.0 - 2.0 / (64f64).ln(), 2.0 - 2.0 / (256f64).ln(), 2.0 - 2.0 / (1000f64).ln());
        }
        println!();
    }
    // 大域 d_H: ⟨r⟩ ~ V^{1/d_H} のサイズスケーリング (3 点の素直なフィット + 隣接ペア)
    println!("[C] 大域ハウスドルフ次元: ⟨r⟩ ~ V^(1/d_H) のサイズスケーリング");
    println!(
        "  ⟨r⟩(V): {}",
        size_data
            .iter()
            .map(|d| format!("V={:.0}:⟨r⟩={:.2}", d.0, d.1))
            .collect::<Vec<_>>()
            .join("  ")
    );
    let xs: Vec<f64> = size_data.iter().map(|d| d.0.ln()).collect();
    let ys: Vec<f64> = size_data.iter().map(|d| d.1.ln()).collect();
    let (_, b) = linfit(&xs, &ys);
    println!("  3点フィット: d_H = {:.2}", 1.0 / b);
    for w in size_data.windows(2) {
        let dh = (w[1].0 / w[0].0).ln() / (w[1].1 / w[0].1).ln();
        println!("  ペア推定 V={:.0}→{:.0}: d_H = {:.2}", w[0].0, w[1].0, dh);
    }
    println!("  (厳密解は d_H=4。これらのサイズでは下から接近する — 既知の強い有限サイズ効果)");
    println!("\n結論: 量子的に揺らぐ 2D 幾何は、体積で測ると ~4 次元 (d_H→4)、拡散で測ると");
    println!("      2 次元 (d_s=2) のフラクタルである。「次元」は基本量ではなく測り方に依存する創発量。");
    println!("      (4D CDT では d_s が大スケール 4 → プランク域 2 へ走る — 次元の走行は量子重力の一般的予言)");
}
