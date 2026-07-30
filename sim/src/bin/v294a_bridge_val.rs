//! v29.4a 定量計量の器械確認とバー導出 — val 区画の 1 回使用 (第三十期, PROMPT/11)
//!
//! 目的: HOLD-5 の val-4..7 (v29.2 で「中間確認 1 回」と凍結した区画) を初めて
//! かつ**最後に** 1 回走らせ、(i) v29.3 凍結の S×C 資格採点器が新鮮な系でも健全で
//! あることを確認し、(ii) **定量計量の読み出し (局所速度場 v̂ の再構成と未使用源の
//! 到着時刻予言) のバーを val 実測から機械規則で導出して本コミットで凍結する**。
//! holdout (hold-8..15) はバー凍結後の v29.4b まで生成すらしない。
//!
//! SECRET の開示: 本ファイルに SECRET を埋め込む (= 開示)。sha256(SECRET) が
//! v29.2 公開のコミットメント cfeb1959… に一致すること、および導出した train seed
//! が v29.3 の公開定数と一致することを [Q0] が機械検査する — HOLD-5 生成系列の
//! 全体が第三者検証可能になる。開示時点で採点器 (v29.3) とバー (本コミット) は
//! 凍結済みのため、holdout に対する調整自由度は残らない。
//!
//! 定量読み出し (FROZEN QUANT v31 節 — 本コミットで凍結):
//!   較正源 2 点の到着時刻の差分から再構成順に沿う**局所 slowness ŝ** (v̂ = a/ŝ が
//!   速度場) を推定し、**未使用源 4 点**の到着時刻を τ_pred = Σ ŝ (経路和) で予言。
//!   採点 3 種: (a) Δ∞(v̂, v_true) = inf_α max |ln(v̂/(α v_true))| (真値照合 —
//!   採点側の権利)、(b) 未使用源の max |ln(t*/τ_pred)| (d ≥ 5, スケール自由度なし)、
//!   (c) regulator 間 Δ∞(v̂_R, v̂_R')。バー値は本走行の val 実測 × 1.5 (小数 2 桁
//!   切り上げ) — 規則自体を本節に凍結し、導出走行が提案値を印字する。
//!
//! バー導出の手続き (as-run 記録): (1) バー = NaN の導出走行で val 実測と提案値を
//! 印字 (results/v294a_bridge_val_derivation.txt) → (2) 提案値を定数へ転記 →
//! (3) 検算走行 (同一 val データ — 乱数固定・決定的) が [Q3] を検査 → コミット =
//! 凍結点。val の統計的使用は 1 回 (転記後の再走行は同一データの検算であり、
//! 設計の自由度はコード内の機械規則に固定されている)。

use std::collections::BTreeMap;
use uft_sim::*;

// ===================== FROZEN BRIDGE READOUT v29 (BEGIN) =====================
// この節は v281_bridge_train.rs と v282_bridge_holdout.rs で逐語同一であること
// (両ファイルのこの節の SHA-256 一致を v282 [H0] が検査する)。
// 節内の定数 = holdout 判定バー。節外での再定義・上書きは禁止。

pub const W_FLOOR: f64 = 1e-8; // これ未満の最大重みは「幾何なし」
pub const BRIDGE_DOM: f64 = 3.0; // 橋対の支配率 (双方の第 2 位の 3 倍以上)
pub const EIG_CLAMP: f64 = 1e-12; // モジュラー核の固有値クランプ
pub const FRONT_EPS: f64 = 1e-3; // 前線到着の閾値 |Δn| ≥ FRONT_EPS
pub const FRONT_DT: f64 = 0.5; // 前線走査の時間刻み (Gaussian)
                               // ---- holdout 判定バー (v28.2 の契約 — 変更禁止) ----
pub const BAR_ADJACENCY: f64 = 1.0; // 隣接復元率 (Cycle/Path の辺が全て真の隣接)
pub const BAR_MIRROR_RATE: f64 = 0.95; // TFD 鏡像橋の検出率
pub const BAR_MDS_PAIR: f64 = 0.02; // 円環 MDS の縮退対 |λ1/λ2 − 1|
pub const BAR_MDS_PAIR2: f64 = 0.05; // 第 2 対 |λ3/λ4 − 1|
pub const BAR_PATH_MDS: f64 = 1e-8; // 線分 MDS の λ2/λ1
pub const BAR_V_REL: f64 = 0.05; // 前線速度の真値からの相対偏差 (訓練資格)
pub const BAR_H4_R2: f64 = 0.98; // HOLD-4: t* vs d_B1 の線形 R²
pub const BAR_H4_VSPREAD: f64 = 0.03; // HOLD-4: 源間の v ばらつき (相対)
pub const BAR_H3_XI_RATIO: (f64, f64) = (0.8, 1.25); // HOLD-3: ξ_stag/ξ_wil の許容窓
pub const BAR_Z2_SPEARMAN: f64 = 0.85; // HOLD-2: B4 到着順 vs 再構成距離の順位相関
pub const BAR_Z2_UNUSED: f64 = 0.90; // HOLD-2: B1 幾何が B3 チャネルを予言 (d≤5)
pub const BAR_ROBUST_PHASE: f64 = 1e-10; // 局所位相不変性 (相対)

/// ノード分割つきガウス状態 — ノード i のモード a は行 i*m + a (列優先, 複素エルミート)。
/// 読み出しに渡る時点でノードラベルは置換済み (幾何情報を含まない)。
pub struct NodeState {
    pub nodes: usize,
    pub m: usize,
    pub cre: Vec<f64>,
    pub cim: Vec<f64>,
}

impl NodeState {
    pub fn dim(&self) -> usize {
        self.nodes * self.m
    }
    /// ノード集合のブロック (列優先 複素)
    pub fn block(&self, ns: &[usize]) -> (Vec<f64>, Vec<f64>, usize) {
        let k = ns.len() * self.m;
        let d = self.dim();
        let mut re = vec![0.0; k * k];
        let mut im = vec![0.0; k * k];
        for (bi, &ni) in ns.iter().enumerate() {
            for a in 0..self.m {
                let gi = ni * self.m + a;
                let li = bi * self.m + a;
                for (bj, &nj) in ns.iter().enumerate() {
                    for b in 0..self.m {
                        let gj = nj * self.m + b;
                        let lj = bj * self.m + b;
                        re[li + lj * k] = self.cre[gi + gj * d];
                        im[li + lj * k] = self.cim[gi + gj * d];
                    }
                }
            }
        }
        (re, im, k)
    }
}

/// B1 (Gaussian): MI(i,j) = S_i + S_j − S_ij
pub fn w_b1_gauss(st: &NodeState) -> Vec<f64> {
    let n = st.nodes;
    let mut s1 = vec![0.0; n];
    for (i, s) in s1.iter_mut().enumerate() {
        let (re, im, k) = st.block(&[i]);
        *s = entropy_corr_herm(&re, &im, k);
    }
    let mut w = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (re, im, k) = st.block(&[i, j]);
            let sij = entropy_corr_herm(&re, &im, k);
            let mi = (s1[i] + s1[j] - sij).max(0.0);
            w[i + j * n] = mi;
            w[j + i * n] = mi;
        }
    }
    w
}

/// 複素エルミート k×k の実埋め込み 2k×2k ([[Re, −Im], [Im, Re]])
pub fn embed_herm(re: &[f64], im: &[f64], k: usize) -> Vec<f64> {
    let m = 2 * k;
    let mut a = vec![0.0; m * m];
    for i in 0..k {
        for j in 0..k {
            a[i + j * m] = re[i + j * k];
            a[i + (j + k) * m] = -im[i + j * k];
            a[(i + k) + j * m] = im[i + j * k];
            a[(i + k) + (j + k) * m] = re[i + j * k];
        }
    }
    a
}

/// B2 (Gaussian): 2 ノードブロックのモジュラー核 k̂ = ln((1−C)/C) のノード間
/// off-diagonal Frobenius ノルム (実埋め込みで matfun — 埋め込みはノルム² を 2 倍にする)
pub fn w_b2_gauss(st: &NodeState) -> Vec<f64> {
    let n = st.nodes;
    let mut w = vec![0.0; n * n];
    let f = |x: f64| {
        let c = x.clamp(EIG_CLAMP, 1.0 - EIG_CLAMP);
        ((1.0 - c) / c).ln()
    };
    for i in 0..n {
        for j in (i + 1)..n {
            let (re, im, k) = st.block(&[i, j]);
            let a = embed_herm(&re, &im, k);
            let ka = matfun_sym(&a, 2 * k, f);
            // ノード間ブロック: 実埋め込みの (行 ∈ i 側, 列 ∈ j 側) と虚部側を全部拾い /2
            let mut nrm2 = 0.0;
            let m2 = 2 * k;
            for a_ in 0..st.m {
                for b_ in st.m..k {
                    // (Re ブロック) と (Im ブロック) の両側
                    for (ri, ci) in [(a_, b_), (a_, b_ + k), (a_ + k, b_), (a_ + k, b_ + k)] {
                        nrm2 += ka[ri + ci * m2] * ka[ri + ci * m2];
                    }
                }
            }
            let v = (nrm2 / 2.0).sqrt();
            w[i + j * n] = v;
            w[j + i * n] = v;
        }
    }
    w
}

/// B3 (Gaussian): 密度応答 |⟨n_i n_j⟩_c| = Σ_{a∈i, b∈j} |C_ab|²
pub fn w_b3_gauss(st: &NodeState) -> Vec<f64> {
    let n = st.nodes;
    let d = st.dim();
    let mut w = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let mut s = 0.0;
            for a in 0..st.m {
                for b in 0..st.m {
                    let gi = i * st.m + a;
                    let gj = j * st.m + b;
                    let (re, im) = (st.cre[gi + gj * d], st.cim[gi + gj * d]);
                    s += re * re + im * im;
                }
            }
            w[i + j * n] = s;
            w[j + i * n] = s;
        }
    }
    w
}

// ---- Z2 (多体) 核 — 決定的 rdm (lib の HashMap 版は反復順が per-process 乱数で
// ρ の浮動小数和の順序が走行ごとに変わる [v153 ドリフトの機構] — BTreeMap で凍結) ----

/// 決定的 rdm: lib.rs Z2CoreState::rdm と同一の符号規約・BTreeMap 群化
pub fn rdm_det(
    l: usize,
    masks: &[u32],
    psi: &[(f64, f64)],
    sites: &[usize],
) -> (Vec<f64>, Vec<f64>, usize) {
    let na = sites.len();
    let dima = 1usize << na;
    let ncomb = masks.len();
    let mut groups: BTreeMap<u64, Vec<(usize, (f64, f64))>> = BTreeMap::new();
    let in_a: Vec<Option<usize>> = (0..l).map(|s| sites.iter().position(|&t| t == s)).collect();
    for (mi_, &mask) in masks.iter().enumerate() {
        for ei in 0..2u64 {
            let idx = mi_ + ncomb * (ei as usize);
            let (ar, ai) = psi[idx];
            if ar == 0.0 && ai == 0.0 {
                continue;
            }
            let mut akey = 0usize;
            let mut bkey = 0u64;
            let mut bpos = 0u32;
            let mut sign = 1.0f64;
            let mut b_seen = 0u32;
            for site in 0..l {
                let occ = (mask >> site) & 1;
                match in_a[site] {
                    Some(k) => {
                        if occ == 1 {
                            akey |= 1 << k;
                            if b_seen % 2 == 1 {
                                sign = -sign;
                            }
                        }
                    }
                    None => {
                        if occ == 1 {
                            bkey |= 1 << bpos;
                            b_seen += 1;
                        }
                        bpos += 1;
                    }
                }
            }
            bkey |= ei << 63;
            groups
                .entry(bkey)
                .or_default()
                .push((akey, (sign * ar, sign * ai)));
        }
    }
    let mut re = vec![0.0; dima * dima];
    let mut im = vec![0.0; dima * dima];
    for (_, g) in groups {
        for &(a1, (x1, y1)) in &g {
            for &(a2, (x2, y2)) in &g {
                re[a1 + a2 * dima] += x1 * x2 + y1 * y2;
                im[a1 + a2 * dima] += y1 * x2 - x1 * y2;
            }
        }
    }
    (re, im, dima)
}

/// Z2 状態のノード置換ビュー: 読み出しは perm 後のノード id しか見ない。
/// sites 引数 (perm 後) を perm 前に写して rdm_det を呼ぶラッパ。
pub struct Z2View<'a> {
    pub l: usize,
    pub masks: &'a [u32],
    pub psi: &'a [(f64, f64)],
    pub to_orig: Vec<usize>, // perm 後 id → 元 id (状態表現の都合 — 読み出しには不可視)
}

impl<'a> Z2View<'a> {
    fn rdm(&self, ns: &[usize]) -> (Vec<f64>, Vec<f64>, usize) {
        let orig: Vec<usize> = ns.iter().map(|&i| self.to_orig[i]).collect();
        rdm_det(self.l, self.masks, self.psi, &orig)
    }
}

/// B1 (Z2): MI from 決定的 rdm
pub fn w_b1_z2(v: &Z2View) -> Vec<f64> {
    let n = v.l;
    let s1: Vec<f64> = (0..n)
        .map(|i| {
            let (re, im, d) = v.rdm(&[i]);
            entropy_rdm_c(&re, &im, d)
        })
        .collect();
    let mut w = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (re, im, d) = v.rdm(&[i, j]);
            let sij = entropy_rdm_c(&re, &im, d);
            let mi = (s1[i] + s1[j] - sij).max(0.0);
            w[i + j * n] = mi;
            w[j + i * n] = mi;
        }
    }
    w
}

/// B2 (Z2): w = ‖ln ρ_ij − ln(ρ_i ⊗ ρ_j)‖_F (4×4, 実埋め込み matfun, クランプ)
pub fn w_b2_z2(v: &Z2View) -> Vec<f64> {
    let n = v.l;
    let lnf = |x: f64| x.max(EIG_CLAMP).ln();
    let single: Vec<(Vec<f64>, Vec<f64>)> = (0..n)
        .map(|i| {
            let (re, im, _) = v.rdm(&[i]);
            (re, im)
        })
        .collect();
    let mut w = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (re, im, d) = v.rdm(&[i, j]); // d = 4, 基底 bit0 = ノード i, bit1 = ノード j
            let a = embed_herm(&re, &im, d);
            let ln_ij = matfun_sym(&a, 2 * d, lnf);
            // ρ_i ⊗ ρ_j (複素 4×4): kron して同様に ln
            let (ri, ii) = (&single[i].0, &single[i].1);
            let (rj, ij_) = (&single[j].0, &single[j].1);
            let mut kre = vec![0.0; d * d];
            let mut kim = vec![0.0; d * d];
            for a1 in 0..2 {
                for a2 in 0..2 {
                    for b1 in 0..2 {
                        for b2 in 0..2 {
                            // index = bit0 (i) + 2*bit1 (j)
                            let r = a1 + 2 * b1;
                            let c = a2 + 2 * b2;
                            let (xr, xi) = (ri[a1 + a2 * 2], ii[a1 + a2 * 2]);
                            let (yr, yi) = (rj[b1 + b2 * 2], ij_[b1 + b2 * 2]);
                            kre[r + c * d] += xr * yr - xi * yi;
                            kim[r + c * d] += xr * yi + xi * yr;
                        }
                    }
                }
            }
            let ap = embed_herm(&kre, &kim, d);
            let ln_p = matfun_sym(&ap, 2 * d, lnf);
            let mut nrm2 = 0.0;
            for t in 0..(2 * d) * (2 * d) {
                let diff = ln_ij[t] - ln_p[t];
                nrm2 += diff * diff;
            }
            let val = (nrm2 / 2.0).sqrt();
            w[i + j * n] = val;
            w[j + i * n] = val;
        }
    }
    w
}

/// B3 (Z2): |⟨n_i n_j⟩ − ⟨n_i⟩⟨n_j⟩| (rdm の対角から — 2×2 の占有確率 = re[1 + 1·2])
pub fn w_b3_z2(v: &Z2View) -> Vec<f64> {
    let n = v.l;
    let dens: Vec<f64> = (0..n)
        .map(|i| {
            let (re, _, _) = v.rdm(&[i]);
            re[3]
        })
        .collect();
    let mut w = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let (re, _, d) = v.rdm(&[i, j]);
            let nn = re[3 + 3 * d]; // |11⟩ 対角 (bit0 = i, bit1 = j)
            let val = (nn - dens[i] * dens[j]).abs();
            w[i + j * n] = val;
            w[j + i * n] = val;
        }
    }
    w
}

// ---- 幾何パイプライン (候補共通) ----

pub struct Comp {
    pub order: Vec<usize>, // 巡回/経路順のノード列
    pub topology: u8,      // 0 = Cycle, 1 = Path, 2 = Other
}

pub struct Recon {
    pub detected: bool,
    pub bridges: Vec<(usize, usize)>,
    pub comps: Vec<Comp>,
    pub edges: Vec<(usize, usize)>,
}

/// 行 i の (最大, 2 位) を (値, 添字) で返す (skip は除外)
fn row_top2(w: &[f64], n: usize, i: usize, skip: Option<usize>) -> ((f64, usize), (f64, usize)) {
    let (mut b1, mut b2) = ((-1.0, usize::MAX), (-1.0, usize::MAX));
    for j in 0..n {
        if j == i || Some(j) == skip {
            continue;
        }
        let v = w[i + j * n];
        if v > b1.0 {
            b2 = b1;
            b1 = (v, j);
        } else if v > b2.0 {
            b2 = (v, j);
        }
    }
    (b1, b2)
}

pub fn reconstruct(w: &[f64], n: usize) -> Recon {
    let wmax = w.iter().cloned().fold(0.0, f64::max);
    if wmax < W_FLOOR {
        return Recon {
            detected: false,
            bridges: vec![],
            comps: vec![],
            edges: vec![],
        };
    }
    // 橋対: 相互 top-1 かつ双方の第 2 位を BRIDGE_DOM 倍以上支配
    let mut bridge_of = vec![usize::MAX; n];
    let mut bridges = Vec::new();
    for i in 0..n {
        let (t1i, t2i) = row_top2(w, n, i, None);
        let j = t1i.1;
        if j == usize::MAX || j < i {
            continue;
        }
        let (t1j, t2j) = row_top2(w, n, j, None);
        if t1j.1 == i
            && t1i.0 >= BRIDGE_DOM * t2i.0.max(0.0)
            && t1j.0 >= BRIDGE_DOM * t2j.0.max(0.0)
        {
            bridge_of[i] = j;
            bridge_of[j] = i;
            bridges.push((i, j));
        }
    }
    // 空間グラフ: 橋のパートナーを除いた相互 top-2
    let mut top2 = vec![[usize::MAX; 2]; n];
    for i in 0..n {
        let skip = if bridge_of[i] != usize::MAX {
            Some(bridge_of[i])
        } else {
            None
        };
        let (t1, t2) = row_top2(w, n, i, skip);
        if t1.0 >= W_FLOOR {
            top2[i][0] = t1.1;
        }
        if t2.0 >= W_FLOOR {
            top2[i][1] = t2.1;
        }
    }
    let mut edges = Vec::new();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for &j in &top2[i] {
            if j == usize::MAX || j < i {
                continue;
            }
            if top2[j].contains(&i) {
                edges.push((i, j));
                adj[i].push(j);
                adj[j].push(i);
            }
        }
    }
    // 成分分解と分類
    let mut seen = vec![false; n];
    let mut comps = Vec::new();
    for s in 0..n {
        if seen[s] || adj[s].is_empty() {
            continue;
        }
        // BFS で成分収集
        let mut stack = vec![s];
        let mut nodes = Vec::new();
        seen[s] = true;
        while let Some(u) = stack.pop() {
            nodes.push(u);
            for &v in &adj[u] {
                if !seen[v] {
                    seen[v] = true;
                    stack.push(v);
                }
            }
        }
        let ne: usize = nodes.iter().map(|&u| adj[u].len()).sum::<usize>() / 2;
        let deg1 = nodes.iter().filter(|&&u| adj[u].len() == 1).count();
        let deg2 = nodes.iter().filter(|&&u| adj[u].len() == 2).count();
        let nn = nodes.len();
        let topology = if deg2 == nn && ne == nn && nn >= 3 {
            0 // Cycle
        } else if deg1 == 2 && deg1 + deg2 == nn && ne + 1 == nn && nn >= 2 {
            1 // Path
        } else {
            2 // Other
        };
        // 巡回/経路順に歩く
        let order = if topology <= 1 {
            let start = if topology == 1 {
                *nodes.iter().find(|&&u| adj[u].len() == 1).unwrap()
            } else {
                nodes[0]
            };
            let mut ord = vec![start];
            let mut prev = usize::MAX;
            let mut cur = start;
            loop {
                let mut nxt = usize::MAX;
                for &v in &adj[cur] {
                    if v != prev {
                        nxt = v;
                        break;
                    }
                }
                if nxt == usize::MAX || nxt == start {
                    break;
                }
                ord.push(nxt);
                prev = cur;
                cur = nxt;
                if ord.len() > nn {
                    break;
                }
            }
            ord
        } else {
            nodes.clone()
        };
        comps.push(Comp { order, topology });
    }
    // 幾何検出 = Cycle/Path (サイズ ≥ 3) がノードの 90% 以上を被覆
    let covered: usize = comps
        .iter()
        .filter(|c| c.topology <= 1 && c.order.len() >= 3)
        .map(|c| c.order.len())
        .sum();
    let detected = covered * 10 >= n * 9;
    Recon {
        detected,
        bridges,
        comps,
        edges,
    }
}

// ---- B4 (因果) — Gaussian: 局所擾乱の前線 ----

/// クエンチ: ノード q に粒子を注入 (行/列を切断し占有を 1 に — 実験側操作)。
/// 位相破壊 (対角 1/2) では二部格子+半充填の粒子正孔対称により**局所密度が厳密に
/// 不変**で前線が読めない (v6.7 の教訓「対称性で保護された観測量は光円錐を運ばない」
/// の再演 — 本監査の設計段階で再発見)。注入は対称性を破り密度前線を作る。
pub fn quench_node(st: &NodeState, q: usize) -> NodeState {
    let d = st.dim();
    let mut cre = st.cre.clone();
    let mut cim = st.cim.clone();
    for a in 0..st.m {
        let g = q * st.m + a;
        for t in 0..d {
            cre[g + t * d] = 0.0;
            cre[t + g * d] = 0.0;
            cim[g + t * d] = 0.0;
            cim[t + g * d] = 0.0;
        }
        cre[g + g * d] = 1.0;
    }
    NodeState {
        nodes: st.nodes,
        m: st.m,
        cre,
        cim,
    }
}

/// 実対称一体ハミルトニアン h の下で C(t) の対角 (ノード密度) を追い、
/// |Δn_j(t)| ≥ FRONT_EPS の最初の t を返す (到達しなければ +∞)。
/// 読み出しが見るのは Δn_j(t) 時系列のみ (h は状態族の生成 = 実験側)。
pub fn fronts_gauss(
    h: &[f64],
    st0: &NodeState,
    base_density: &[f64],
    sources: &[usize],
    nt: usize,
) -> Vec<Vec<f64>> {
    let d = st0.dim();
    let (ev, vv) = jacobi_eigh(h, d);
    let mut out = Vec::new();
    for &q in sources {
        let stq = quench_node(st0, q);
        // A = Vᵀ C' V (実部のみ — C' 実対称の系で使う)
        let mut tmp = vec![0.0; d * d];
        for a in 0..d {
            for t in 0..d {
                let mut s = 0.0;
                for r in 0..d {
                    s += vv[r + a * d] * stq.cre[r + t * d];
                }
                tmp[a + t * d] = s;
            }
        }
        let mut am = vec![0.0; d * d];
        for a in 0..d {
            for b in 0..d {
                let mut s = 0.0;
                for t in 0..d {
                    s += tmp[a + t * d] * vv[t + b * d];
                }
                am[a + b * d] = s;
            }
        }
        let mut tstar = vec![f64::INFINITY; st0.nodes];
        let mut remaining = st0.nodes;
        for it in 1..=nt {
            if remaining == 0 {
                break;
            }
            let t = FRONT_DT * it as f64;
            let (c, s): (Vec<f64>, Vec<f64>) =
                ev.iter().map(|&e| ((e * t).cos(), (e * t).sin())).unzip();
            // n_j(t) = Σ_ab V_ja V_jb A_ab cos((Ea−Eb)t)
            //        = (Vc)_j A (Vc)_jᵀ + (Vs)_j A (Vs)_jᵀ
            // Wc = V·diag(c), Ws = V·diag(s); Y = A Wcᵀ, Z = A Wsᵀ
            let mut yc = vec![0.0; d * d];
            let mut zs = vec![0.0; d * d];
            for a in 0..d {
                for j in 0..d {
                    let (mut sy, mut sz) = (0.0, 0.0);
                    for b in 0..d {
                        let vjb = vv[j + b * d];
                        sy += am[a + b * d] * vjb * c[b];
                        sz += am[a + b * d] * vjb * s[b];
                    }
                    yc[a + j * d] = sy;
                    zs[a + j * d] = sz;
                }
            }
            for node in 0..st0.nodes {
                if tstar[node].is_finite() {
                    continue;
                }
                let mut dev: f64 = 0.0;
                for a_ in 0..st0.m {
                    let j = node * st0.m + a_;
                    let mut nj = 0.0;
                    for a in 0..d {
                        let vja = vv[j + a * d];
                        nj += vja * (c[a] * yc[a + j * d] + s[a] * zs[a + j * d]);
                    }
                    dev = dev.max((nj - base_density[j]).abs());
                }
                if dev >= FRONT_EPS {
                    tstar[node] = t;
                    remaining -= 1;
                }
            }
        }
        out.push(tstar);
    }
    out
}

/// 順位相関 (Spearman, 同順位は添字順で決定的に処理)
pub fn spearman(x: &[f64], y: &[f64]) -> f64 {
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

// ====================== FROZEN BRIDGE READOUT v29 (END) ======================

// ================== FROZEN BRIDGE MATRIX v30 (BEGIN) ==================
// 新候補 (B5/B6)・HOLD-5 生成器・合成採点の定義とバー。v29.3 のコミットで凍結 —
// v29.4 (val/holdout) は本節の SHA-256 一致を検査する。節外での再定義禁止。

/// HOLD-5 train インスタンス seed (= sha256(SECRET+":"+id) 先頭 8 バイト, 公開可 —
/// SECRET 開示時に第三者が検証する)。val は v29.4 で開示、holdout は走行版まで非生成。
pub const HOLD5_TRAIN_SEEDS: [(u64, &str); 4] = [
    (16960242293642875863, "train-0"),
    (3691857198250009775, "train-1"),
    (2113394486573742281, "train-2"),
    (493408761443360655, "train-3"),
];

pub const B6_EPS: f64 = 1e-3; // COMMUTATOR/KINEMATIC の到着閾値
pub const FRONT_DT5: f64 = 0.25; // HOLD-5 の前線走査刻み (v(x) < 1.4 のため細かく)
pub const FRONT_NT5: usize = 400; // t_max = 100
/// S×C 合成セルの資格バー (順位相関 — train 設計後に凍結した値):
/// train 実測 (4 系 × 3 regulator × 2 質量 × 4 空間核 × 3 因果チャネル = 288 セル)
/// は全セル ≥ 0.9718 — 余裕をとって 0.90 で凍結
pub const BAR_SXC_SPEARMAN: f64 = 0.90;
/// 空間セルの資格バー: 隣接 100% (v29 節の BAR_ADJACENCY を流用)・位相は真値と一致

/// HOLD-5 の系仕様 (seed から一意に決定)
pub struct Hold5Spec {
    pub n_phys: usize, // 物理サイト数 (a = 1 での N)
    pub ring: bool,
    pub amps: Vec<f64>,
    pub centers: Vec<f64>,
    pub widths: Vec<f64>,
}

/// v(x) = 1 + Σ A_k exp(−((x−c_k)/w_k)²) — 滑らかな正値速度プロファイル
pub fn hold5_profile(seed: u64) -> Hold5Spec {
    let n_phys = 128usize;
    let mut rng = Rng::new(seed);
    let ring = rng.f64() < 0.5;
    loop {
        let k = 3 + rng.range(3); // K ∈ {3,4,5}
        let mut amps = Vec::new();
        let mut centers = Vec::new();
        let mut widths = Vec::new();
        for _ in 0..k {
            amps.push(-0.3 + 0.7 * rng.f64()); // [−0.3, 0.4]
            centers.push(n_phys as f64 * rng.f64());
            widths.push(8.0 + 12.0 * rng.f64()); // [8, 20]
        }
        let spec = Hold5Spec {
            n_phys,
            ring,
            amps,
            centers,
            widths,
        };
        // v(x) > 0.5 を保証 (棄却サンプリング — 決定的: 同じ rng 列を続ける)
        let ok = (0..(4 * n_phys)).all(|q| v_of_x(&spec, q as f64 / 4.0) > 0.5);
        if ok {
            return spec;
        }
    }
}

pub fn v_of_x(spec: &Hold5Spec, x: f64) -> f64 {
    let mut v = 1.0;
    let n = spec.n_phys as f64;
    for k in 0..spec.amps.len() {
        let mut dx = x - spec.centers[k];
        if spec.ring {
            // 周期距離 (リングでは v も周期的)
            dx = dx - n * (dx / n).round();
        }
        v += spec.amps[k] * (-(dx / spec.widths[k]).powi(2)).exp();
    }
    v
}

/// HOLD-5 の系 (実対称 H と基底状態 NodeState)。regulator:
///   0 = R-A 単鎖 (a = 1, N = n_phys, ボンド t = v(x+1/2)/2, 質量 ±m·(−1)^x)
///   1 = R-B 細分格子 (a = 1/2, N = 2·n_phys, t = v/(2a) = v, 質量 ±(m/2)·(−1)^x)
///   2 = R-C Wilson (a = 1, N = n_phys, 2 成分 — 実ゲージ i^x で実対称化。
///       hopping = v(x+1/2)/2·[σ₁ − r·σ₃]/…: sin 項 σ₁/2・r 項 −σ₃/2、
///       on-site σ₃·(m + r·v̄(x)) with v̄ = 隣接ボンドの平均)
/// 返り値: (H 実対称 (dim²), NodeState (基底状態), 真の隣接判定用のノード数)
pub fn hold5_system(spec: &Hold5Spec, regulator: u8, m_phys: f64) -> (Vec<f64>, NodeState) {
    match regulator {
        0 | 1 => {
            let scale = if regulator == 0 { 1usize } else { 2usize };
            let n = spec.n_phys * scale;
            let a = 1.0 / scale as f64;
            let mut h = vec![0.0; n * n];
            let nb = if spec.ring { n } else { n - 1 };
            for b in 0..nb {
                let x_phys = (b as f64 + 0.5) * a;
                let t = v_of_x(spec, x_phys) / (2.0 * a);
                let (i, j) = (b, (b + 1) % n);
                h[i + j * n] = -t;
                h[j + i * n] = -t;
            }
            let m_lat = m_phys * a;
            for x in 0..n {
                h[x + x * n] += if x % 2 == 0 { m_lat } else { -m_lat };
            }
            let c = gs_corr_real(&h, n, n / 2);
            (
                h,
                NodeState {
                    nodes: n,
                    m: 1,
                    cre: c,
                    cim: vec![0.0; n * n],
                },
            )
        }
        _ => {
            // R-C Wilson (実ゲージ): サイト x に 2 成分 (u, d)。
            // 一様極限 h(k) = σ₁ sin k + σ₃(m + r(1−cos k)), r = 1 の位置空間形を
            // i^x ゲージで実化: hopping ブロック = t·(σ₁/2 − σ₃/2)? 導出:
            //   sin 項: ψ†_x(−iσ₁/2)ψ_{x+1} + h.c. → (ゲージ後) ψ†_x(σ₁/2)ψ_{x+1} + h.c.
            //   r 項: −ψ†_x(σ₃/2)ψ_{x+1} + h.c. はゲージで −(i)·σ₃/2 → 虚?…
            // 実化は σ₁ hopping のみに作用させるため、r 項は素の実 hopping のまま
            // (−σ₃/2) とし、ゲージは sin 項の i を回す: i^x ゲージで
            //   −i σ₁/2 → σ₁/2 (x 依存符号は (i)^{x+1−x} = i, −i·i = 1)
            //   −σ₃/2 → −σ₃/2 · i^{1} … 虚になるため、代わりに r 項も同じゲージ下で
            //   実に保つには σ₃ hopping に i が付かないことが必要 — 付く。
            // 正しい実化: ゲージ U_x = (σ₃)^x を併用すると σ₁ ↔ −σ₁ が交代し
            //   i^x (σ₃)^x の複合ゲージで両項が実になる (下の分散照合 [M0] が
            //   数値的にこの構成の正しさを担保する)。実装は複合ゲージ後の値を直書き:
            //   hopping ブロック B_x = t_x · [[ −r/2, +1/2 ], [ −1/2, −r/2 ]]
            //   (非対称 2×2 — H 全体では h.c. と合わせて実対称)
            //   on-site ブロック = σ₃·(m + r·v̄_x)
            let n = spec.n_phys;
            let d = 2 * n;
            let r_w = 1.0;
            let mut h = vec![0.0; d * d];
            let nb = if spec.ring { n } else { n - 1 };
            let mut tvals = vec![0.0; n];
            for b in 0..nb {
                let t = v_of_x(spec, b as f64 + 0.5) / 2.0;
                tvals[b] = t;
                let (x, y) = (b, (b + 1) % n);
                // B = t · [[−r/2? — 実装は成分で]]: u†_x u_y: −t·r_w/2? …
                // 成分: (u_x, d_x) → ブロック [[a, b], [c, e]] を h[2x+α, 2y+β]
                // 実化 (site 非依存の基底回転 T = diag(1, i)): σ₁ → σ₂,
                // hopping ブロック B = −t·(iσ₂) − r·t·σ₃ = [[−rt, −t],[t, rt]], t = v/2。
                // 一様極限: B e^{ik} + Bᵀ e^{−ik} = v sin k·σ₂ − r v cos k·σ₃ ✓
                let blk = [[-r_w * t, -t], [t, r_w * t]];
                for al in 0..2 {
                    for be in 0..2 {
                        h[(2 * x + al) + (2 * y + be) * d] += blk[al][be];
                        h[(2 * y + be) + (2 * x + al) * d] += blk[al][be];
                    }
                }
            }
            for x in 0..n {
                let vbar = if spec.ring {
                    (tvals[x] + tvals[(x + n - 1) % n]) // = (v(x+½)+v(x−½))/2
                } else if x == 0 {
                    2.0 * tvals[0]
                } else if x == n - 1 {
                    2.0 * tvals[n - 2]
                } else {
                    tvals[x] + tvals[x - 1]
                };
                let on = m_phys + r_w * vbar;
                h[(2 * x) + (2 * x) * d] += on;
                h[(2 * x + 1) + (2 * x + 1) * d] += -on;
            }
            let c = gs_corr_real(&h, d, d / 2);
            (
                h,
                NodeState {
                    nodes: n,
                    m: 2,
                    cre: c,
                    cim: vec![0.0; d * d],
                },
            )
        }
    }
}

/// 実対称 H の基底状態相関行列 (最低 nocc 準位)
pub fn gs_corr_real(h: &[f64], n: usize, nocc: usize) -> Vec<f64> {
    let (ev, vv) = jacobi_eigh(h, n);
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| ev[a].partial_cmp(&ev[b]).unwrap());
    let mut c = vec![0.0; n * n];
    for &k in idx.iter().take(nocc) {
        for i in 0..n {
            for j in 0..n {
                c[i + j * n] += vv[i + k * n] * vv[j + k * n];
            }
        }
    }
    c
}

// ---- B5-QFI: BKM 核 (数保存 Gaussian の厳密 Fock 構成) ----

/// 実小行列 (k ≤ 4) の行列式 (Gauss 消去)
pub fn det_small(a: &[f64], k: usize) -> f64 {
    let mut m = a.to_vec();
    let mut det = 1.0;
    for col in 0..k {
        // pivot
        let mut p = col;
        for r in (col + 1)..k {
            if m[r + col * k].abs() > m[p + col * k].abs() {
                p = r;
            }
        }
        if m[p + col * k].abs() < 1e-300 {
            return 0.0;
        }
        if p != col {
            for cc in 0..k {
                m.swap(col + cc * k, p + cc * k);
            }
            det = -det;
        }
        det *= m[col + col * k];
        for r in (col + 1)..k {
            let f = m[r + col * k] / m[col + col * k];
            for cc in col..k {
                m[r + cc * k] -= f * m[col + cc * k];
            }
        }
    }
    det
}

/// 実小行列 (k ≤ 4) の逆行列 (Gauss-Jordan)。特異なら None
pub fn inv_small(a: &[f64], k: usize) -> Option<Vec<f64>> {
    let mut m = a.to_vec();
    let mut inv = vec![0.0; k * k];
    for i in 0..k {
        inv[i + i * k] = 1.0;
    }
    for col in 0..k {
        let mut p = col;
        for r in (col + 1)..k {
            if m[r + col * k].abs() > m[p + col * k].abs() {
                p = r;
            }
        }
        if m[p + col * k].abs() < 1e-12 {
            return None;
        }
        if p != col {
            for cc in 0..k {
                m.swap(col + cc * k, p + cc * k);
                inv.swap(col + cc * k, p + cc * k);
            }
        }
        let piv = m[col + col * k];
        for cc in 0..k {
            m[col + cc * k] /= piv;
            inv[col + cc * k] /= piv;
        }
        for r in 0..k {
            if r == col {
                continue;
            }
            let f = m[r + col * k];
            for cc in 0..k {
                m[r + cc * k] -= f * m[col + cc * k];
                inv[r + cc * k] -= f * inv[col + cc * k];
            }
        }
    }
    Some(inv)
}

/// 数保存 Gaussian の Fock 行列要素: ρ_{S'S} = det(1−G)·det(K[S'|S]),
/// K = G(1−G)^{-1} (|S'| = |S| のみ非零)。G = 実対称 k×k (固有値を [ε, 1−ε] に
/// クランプしてから構成 — fail-closed)。返り値 2^k × 2^k 実対称。
pub fn fock_rho(g: &[f64], k: usize) -> Vec<f64> {
    // クランプ: G = V diag(clamp λ) Vᵀ
    let (ev, vv) = jacobi_eigh(g, k);
    let mut gc = vec![0.0; k * k];
    for a in 0..k {
        let lam = ev[a].clamp(1e-10, 1.0 - 1e-10);
        for i in 0..k {
            for j in 0..k {
                gc[i + j * k] += lam * vv[i + a * k] * vv[j + a * k];
            }
        }
    }
    let mut one_minus = vec![0.0; k * k];
    for i in 0..k {
        for j in 0..k {
            one_minus[i + j * k] = (if i == j { 1.0 } else { 0.0 }) - gc[i + j * k];
        }
    }
    let det0 = det_small(&one_minus, k);
    let inv1mg = inv_small(&one_minus, k).expect("1−G 可逆 (クランプ済み)");
    // K = G (1−G)^{-1}
    let mut kk = vec![0.0; k * k];
    for i in 0..k {
        for j in 0..k {
            let mut s = 0.0;
            for t in 0..k {
                s += gc[i + t * k] * inv1mg[t + j * k];
            }
            kk[i + j * k] = s;
        }
    }
    let dim = 1usize << k;
    let mut rho = vec![0.0; dim * dim];
    // 占有部分集合の列挙 (ビット昇順) — 部分行列 K[S'|S]
    for sp in 0..dim {
        let np = (sp as u32).count_ones();
        for s in 0..dim {
            if (s as u32).count_ones() != np {
                continue;
            }
            let rows: Vec<usize> = (0..k).filter(|b| (sp >> b) & 1 == 1).collect();
            let cols: Vec<usize> = (0..k).filter(|b| (s >> b) & 1 == 1).collect();
            let kn = rows.len();
            let val = if kn == 0 {
                det0
            } else {
                let mut sub = vec![0.0; kn * kn];
                for (ri, &r) in rows.iter().enumerate() {
                    for (ci, &c) in cols.iter().enumerate() {
                        sub[ri + ci * kn] = kk[r + c * k];
                    }
                }
                det0 * det_small(&sub, kn)
            };
            rho[sp + s * dim] = val;
        }
    }
    rho
}

/// BKM 重み c(λa, λb) = (λa − λb)/(ln λa − ln λb) (等値では λ)
pub fn bkm_weight(la: f64, lb: f64) -> f64 {
    let (a, b) = (la.max(1e-14), lb.max(1e-14));
    if (a - b).abs() < 1e-12 * a.max(b) {
        a
    } else {
        (a - b) / (a.ln() - b.ln())
    }
}

/// B5-QFI 核: w(i,j) = ‖ρ_ij − ρ_i⊗ρ_j‖_BKM (BKM 計量は ρ_ij 基準)
pub fn w_b5_gauss(st: &NodeState) -> Vec<f64> {
    let n = st.nodes;
    let m = st.m;
    // 単ノード ρ (2^m 次)
    let singles: Vec<Vec<f64>> = (0..n)
        .map(|i| {
            let (re, _, k) = st.block(&[i]);
            fock_rho(&re, k)
        })
        .collect();
    let mut w = vec![0.0; n * n];
    let k2 = 2 * m;
    let dim = 1usize << k2;
    let dm = 1usize << m;
    for i in 0..n {
        for j in (i + 1)..n {
            let (re, _, _) = st.block(&[i, j]);
            let rho = fock_rho(&re, k2);
            // ρ_i ⊗ ρ_j (占有ビット: 下位 m ビット = ノード i, 上位 = ノード j —
            // fock_rho のビット規約 (block の行順 = i のモード, 次に j) と一致)
            let mut prod = vec![0.0; dim * dim];
            for a1 in 0..dm {
                for a2 in 0..dm {
                    for b1 in 0..dm {
                        for b2 in 0..dm {
                            prod[(a1 + (b1 << m)) + (a2 + (b2 << m)) * dim] +=
                                singles[i][a1 + a2 * dm] * singles[j][b1 + b2 * dm];
                        }
                    }
                }
            }
            // Δ = ρ − ρ_i⊗ρ_j を ρ の固有基底で BKM ノルム
            let (ev, vv) = jacobi_eigh(&rho, dim);
            let mut nrm2 = 0.0;
            // Δ̃_ab = v_aᵀ Δ v_b
            let mut delta = vec![0.0; dim * dim];
            for t in 0..dim * dim {
                delta[t] = rho[t] - prod[t];
            }
            for a in 0..dim {
                for b in 0..dim {
                    let mut s = 0.0;
                    for p in 0..dim {
                        let mut vp = 0.0;
                        for q in 0..dim {
                            vp += delta[p + q * dim] * vv[q + b * dim];
                        }
                        s += vv[p + a * dim] * vp;
                    }
                    nrm2 += s * s / bkm_weight(ev[a], ev[b]).max(1e-14);
                }
            }
            // BKM ノルム (計量の逆重みで測る情報計量型)。定義は v29.2 追補の
            // ‖·‖_BKM — 実装規約: ⟨Δ, Δ⟩ = Σ |Δ̃_ab|²/c(λa, λb) (KMB 内積の双対)
            let val = nrm2.sqrt();
            w[i + j * n] = val;
            w[j + i * n] = val;
        }
    }
    w
}

// ---- B6-COMMUTATOR: retarded [n_j(t), n_q] と運動学変種 ----

/// 実対称 H の下で源 q からの (a) |⟨[n_j(t), n_q]⟩| 前線 (状態依存) と
/// (b) |P(t)_{qj}| 前線 (運動学) の到着時刻を返す。
/// 導出: ⟨[n_j(t), n_q]⟩ = 2i·Im(P_{qj}(t)·w_j), w = P(t)† C e_q, P = e^{−iht}
pub fn fronts_commutator(
    h: &[f64],
    c0: &[f64],
    n: usize,
    q: usize,
    nt: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (ev, vv) = jacobi_eigh(h, n);
    let ce: Vec<f64> = (0..n)
        .map(|i| {
            let mut s = 0.0;
            for t in 0..n {
                s += c0[i + t * n] * if t == q { 1.0 } else { 0.0 };
            }
            s
        })
        .collect(); // C e_q
    let vq: Vec<f64> = (0..n).map(|a| vv[q + a * n]).collect();
    let vce: Vec<f64> = (0..n)
        .map(|a| {
            let mut s = 0.0;
            for i in 0..n {
                s += vv[i + a * n] * ce[i];
            }
            s
        })
        .collect(); // Vᵀ C e_q
    let mut t_comm = vec![f64::INFINITY; n];
    let mut t_kin = vec![f64::INFINITY; n];
    let mut remaining = 2 * n;
    for it in 1..=nt {
        if remaining == 0 {
            break;
        }
        let t = FRONT_DT5 * it as f64;
        let cs: Vec<f64> = ev.iter().map(|&e| (e * t).cos()).collect();
        let sn: Vec<f64> = ev.iter().map(|&e| (e * t).sin()).collect();
        for j in 0..n {
            if t_comm[j].is_finite() && t_kin[j].is_finite() {
                continue;
            }
            // P_{qj} = Σ_a v_{qa} e^{−iE_a t} v_{ja} / w_j = Σ_a v_{ja} e^{+iE_a t} (VᵀCe_q)_a
            let (mut pr, mut pi, mut wr, mut wi) = (0.0, 0.0, 0.0, 0.0);
            for a in 0..n {
                let vja = vv[j + a * n];
                pr += vq[a] * cs[a] * vja;
                pi -= vq[a] * sn[a] * vja;
                wr += vja * cs[a] * vce[a];
                wi += vja * sn[a] * vce[a];
            }
            if t_kin[j].is_infinite() && (pr * pr + pi * pi).sqrt() >= B6_EPS {
                t_kin[j] = t;
                remaining -= 1;
            }
            if t_comm[j].is_infinite() && 2.0 * (pr * wi + pi * wr).abs() >= B6_EPS {
                t_comm[j] = t;
                remaining -= 1;
            }
        }
    }
    (t_comm, t_kin)
}

/// HOLD-5 の再構成パイプライン (順序付き — Occam): まず**橋なし**の相互 top-2
/// グラフで再構成し、幾何が検出されればそれを採る。失敗した場合のみ v29 節の
/// 橋つき再構成 (TFD 型の橋構造用) にフォールバックする。
/// 設計根拠 (train 区画で確定): 不均一計量の強ボンド対が v29 の橋判定 (3× 支配)
/// を発火させ B1/B5 の位相を壊した — v28 の「B1 境界増強」のバルク版。橋の導入は
/// 必要になったときだけ行う。
pub fn reconstruct_v30(w: &[f64], n: usize) -> Recon {
    let wmax = w.iter().cloned().fold(0.0, f64::max);
    if wmax < W_FLOOR {
        return Recon {
            detected: false,
            bridges: vec![],
            comps: vec![],
            edges: vec![],
        };
    }
    // 橋なしの相互 top-2
    let mut top2 = vec![[usize::MAX; 2]; n];
    for i in 0..n {
        let (t1, t2) = row_top2(w, n, i, None);
        if t1.0 >= W_FLOOR {
            top2[i][0] = t1.1;
        }
        if t2.0 >= W_FLOOR {
            top2[i][1] = t2.1;
        }
    }
    let mut edges = Vec::new();
    let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
    for i in 0..n {
        for &j in &top2[i] {
            if j == usize::MAX || j < i {
                continue;
            }
            if top2[j].contains(&i) {
                edges.push((i, j));
                adj[i].push(j);
                adj[j].push(i);
            }
        }
    }
    let mut seen = vec![false; n];
    let mut comps = Vec::new();
    for s in 0..n {
        if seen[s] || adj[s].is_empty() {
            continue;
        }
        let mut stack = vec![s];
        let mut nodes = Vec::new();
        seen[s] = true;
        while let Some(u) = stack.pop() {
            nodes.push(u);
            for &v in &adj[u] {
                if !seen[v] {
                    seen[v] = true;
                    stack.push(v);
                }
            }
        }
        let ne: usize = nodes.iter().map(|u| adj[*u].len()).sum::<usize>() / 2;
        let deg1 = nodes.iter().filter(|&&u| adj[u].len() == 1).count();
        let deg2 = nodes.iter().filter(|&&u| adj[u].len() == 2).count();
        let nn = nodes.len();
        let topology = if deg2 == nn && ne == nn && nn >= 3 {
            0
        } else if deg1 == 2 && deg1 + deg2 == nn && ne + 1 == nn && nn >= 2 {
            1
        } else {
            2
        };
        let order = if topology <= 1 {
            let start = if topology == 1 {
                *nodes.iter().find(|&&u| adj[u].len() == 1).unwrap()
            } else {
                nodes[0]
            };
            let mut ord = vec![start];
            let mut prev = usize::MAX;
            let mut cur = start;
            loop {
                let mut nxt = usize::MAX;
                for &v in &adj[cur] {
                    if v != prev {
                        nxt = v;
                        break;
                    }
                }
                if nxt == usize::MAX || nxt == start {
                    break;
                }
                ord.push(nxt);
                prev = cur;
                cur = nxt;
                if ord.len() > nn {
                    break;
                }
            }
            ord
        } else {
            nodes.clone()
        };
        comps.push(Comp { order, topology });
    }
    let covered: usize = comps
        .iter()
        .filter(|c| c.topology <= 1 && c.order.len() >= 3)
        .map(|c| c.order.len())
        .sum();
    if covered * 10 >= n * 9 {
        return Recon {
            detected: true,
            bridges: vec![],
            comps,
            edges,
        };
    }
    // フォールバック: 橋つき (v29 凍結節の reconstruct)
    reconstruct(w, n)
}

// =================== FROZEN BRIDGE MATRIX v30 (END) ===================

// ================== FROZEN QUANT v31 (BEGIN) ==================
// 定量計量の読み出し・採点定義とバー。v29.4a のコミットで凍結 — v29.4b (holdout)
// は本節の SHA-256 一致を検査する。節外での再定義禁止。

/// SECRET (開示 — v29.2 コミットメント cfeb1959… の原像。開示と同時に採点器・
/// バーが凍結済みであることが git 履歴で検証可能)
pub const HOLD5_SECRET: &str = "HOLD5-d8572eb883fcec9985ffdbb10ea8e510";

/// instance seed = sha256(SECRET + ":" + id) の先頭 8 バイト (big-endian)
pub fn hold5_seed(id: &str) -> u64 {
    let s = format!("{}:{}", HOLD5_SECRET, id);
    let h = sha256(s.as_bytes());
    u64::from_be_bytes([h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]])
}

/// 定量系列の代表セル: 空間核 = B3-COV (第二十九期の生存 S 核。隣接 100% 資格の
/// 下では順序は S 核に依らず同一のため代表 1 核)、因果チャネル = DENSITY-FRONT。
/// 質量は m = 0 のみ (前線が鋭い系列)。資格 [Q2] は全セル (S×C×質量) で検査する。
pub const QUANT_S: &str = "B3-COV";
pub const QUANT_C: &str = "DENSITY-FRONT";

/// 較正源 2 点 (再構成順位置の 0.2, 0.8) と未使用の予言対象源 4 点
pub const CAL_FRACS: [f64; 2] = [0.2, 0.8];
pub const PRED_FRACS: [f64; 4] = [0.1, 0.4, 0.6, 0.9];
/// slowness 推定: 源からの hop 距離 ≥ D_NEAR_HOPS (近接場除外)、リングは
/// ≤ n·D_FAR_FRAC (対蹠点の両回り前線衝突を回避)。基線半幅 W_SLOW (bond)。
pub const D_NEAR_HOPS: usize = 4;
pub const D_FAR_FRAC: f64 = 0.4;
pub const W_SLOW: usize = 3;
/// τ 予言の採点対象: hop 距離 ≥ D_PRED_MIN かつ t* が走査窓内かつ経路の ŝ が全て定義
pub const D_PRED_MIN: usize = 5;
/// 定量系列の前線走査 (資格系列の FRONT_DT = 0.5 では Δt が量子化され per-bond
/// slowness が壊れるため細分化。算法は v29 節 fronts_gauss と同一)
pub const FRONT_DT_Q: f64 = 0.1;
pub const FRONT_NT_Q: usize = 2500; // t_max = 250

/// fronts_gauss (v29 凍結節) の刻み引数版 — 算法・閾値 (FRONT_EPS) は同一
pub fn fronts_gauss_q(
    h: &[f64],
    st0: &NodeState,
    base_density: &[f64],
    sources: &[usize],
    dt: f64,
    nt: usize,
) -> Vec<Vec<f64>> {
    let d = st0.dim();
    let (ev, vv) = jacobi_eigh(h, d);
    let mut out = Vec::new();
    for &q in sources {
        let stq = quench_node(st0, q);
        let mut tmp = vec![0.0; d * d];
        for a in 0..d {
            for t in 0..d {
                let mut s = 0.0;
                for r in 0..d {
                    s += vv[r + a * d] * stq.cre[r + t * d];
                }
                tmp[a + t * d] = s;
            }
        }
        let mut am = vec![0.0; d * d];
        for a in 0..d {
            for b in 0..d {
                let mut s = 0.0;
                for t in 0..d {
                    s += tmp[a + t * d] * vv[t + b * d];
                }
                am[a + b * d] = s;
            }
        }
        let mut tstar = vec![f64::INFINITY; st0.nodes];
        let mut remaining = st0.nodes;
        for it in 1..=nt {
            if remaining == 0 {
                break;
            }
            let t = dt * it as f64;
            let (c, s): (Vec<f64>, Vec<f64>) =
                ev.iter().map(|&e| ((e * t).cos(), (e * t).sin())).unzip();
            let mut yc = vec![0.0; d * d];
            let mut zs = vec![0.0; d * d];
            for a in 0..d {
                for j in 0..d {
                    let (mut sy, mut sz) = (0.0, 0.0);
                    for b in 0..d {
                        let vjb = vv[j + b * d];
                        sy += am[a + b * d] * vjb * c[b];
                        sz += am[a + b * d] * vjb * s[b];
                    }
                    yc[a + j * d] = sy;
                    zs[a + j * d] = sz;
                }
            }
            for node in 0..st0.nodes {
                if tstar[node].is_finite() {
                    continue;
                }
                let mut dev: f64 = 0.0;
                for a_ in 0..st0.m {
                    let j = node * st0.m + a_;
                    let mut nj = 0.0;
                    for a in 0..d {
                        let vja = vv[j + a * d];
                        nj += vja * (c[a] * yc[a + j * d] + s[a] * zs[a + j * d]);
                    }
                    dev = dev.max((nj - base_density[j]).abs());
                }
                if dev >= FRONT_EPS {
                    tstar[node] = t;
                    remaining -= 1;
                }
            }
        }
        out.push(tstar);
    }
    out
}

/// 再構成順に沿う局所 slowness ŝ_u (bond u = 順位置 u→u+1) を較正源の到着時刻から
/// 推定する。bond u ごとに、源から片側 (単調到着) の基線 [k−W, k+1+W] の時間差 /
/// bond 数を標本とし、全標本 (源 × 方向) を平均。標本なしの bond は NaN。
/// 読み出しが見るのは (順位置, t*) のみ — 真の座標・速度は使わない。
pub fn slowness_from_cal(cyclic: bool, cal: &[(usize, Vec<f64>)], n: usize) -> Vec<f64> {
    let nb = if cyclic { n } else { n - 1 };
    let kmax: usize = if cyclic {
        (n as f64 * D_FAR_FRAC) as usize
    } else {
        n - 1
    };
    let mut sum = vec![0.0; nb];
    let mut cnt = vec![0usize; nb];
    for &(qp, ref tv) in cal {
        for dir in [1isize, -1isize] {
            let pos = |k: usize| -> Option<usize> {
                let p = qp as isize + dir * k as isize;
                if cyclic {
                    Some(p.rem_euclid(n as isize) as usize)
                } else if p < 0 || p >= n as isize {
                    None
                } else {
                    Some(p as usize)
                }
            };
            for k in (D_NEAR_HOPS + W_SLOW)..kmax {
                let a = k - W_SLOW;
                let b = k + 1 + W_SLOW;
                if b > kmax {
                    continue;
                }
                let (pa, pb) = match (pos(a), pos(b)) {
                    (Some(x), Some(y)) => (x, y),
                    _ => continue,
                };
                let (ta, tb) = (tv[pa], tv[pb]);
                if !ta.is_finite() || !tb.is_finite() || tb <= ta {
                    continue;
                }
                let sbar = (tb - ta) / (b - a) as f64;
                let bond = match (pos(k), pos(k + 1)) {
                    (Some(p0), Some(p1)) => {
                        if cyclic {
                            if (p0 + 1) % n == p1 {
                                p0
                            } else {
                                p1
                            }
                        } else {
                            p0.min(p1)
                        }
                    }
                    _ => continue,
                };
                if bond < nb {
                    sum[bond] += sbar;
                    cnt[bond] += 1;
                }
            }
        }
    }
    (0..nb)
        .map(|b| {
            if cnt[b] > 0 {
                sum[b] / cnt[b] as f64
            } else {
                f64::NAN
            }
        })
        .collect()
}

/// 順位置 from→to を dir 方向に歩く経路の Σ ŝ (経路上に未定義 bond があれば None)
pub fn tau_arc(
    shat: &[f64],
    n: usize,
    cyclic: bool,
    from_pos: usize,
    to_pos: usize,
    dir: isize,
) -> Option<f64> {
    let nb = shat.len();
    let mut cur = from_pos as isize;
    let mut tau = 0.0;
    let mut steps = 0usize;
    while cur != to_pos as isize {
        let nxt = if cyclic {
            (cur + dir).rem_euclid(n as isize)
        } else {
            cur + dir
        };
        if !cyclic && (nxt < 0 || nxt >= n as isize) {
            return None;
        }
        let bond = if dir == 1 {
            cur.rem_euclid(n as isize) as usize
        } else {
            nxt.rem_euclid(n as isize) as usize
        };
        if bond >= nb || !shat[bond].is_finite() {
            return None;
        }
        tau += shat[bond];
        cur = nxt;
        steps += 1;
        if steps > n {
            return None;
        }
    }
    Some(tau)
}

/// Δ∞(u, w) = inf_{α>0} max |ln(u_b/(α w_b))| (両方有限・正の bond のみ。
/// 1D の Chebyshev 中心 α* = exp((max+min)/2) で厳密に最適化)
pub fn delta_inf(u: &[f64], w: &[f64]) -> f64 {
    let mut logs: Vec<f64> = u
        .iter()
        .zip(w.iter())
        .filter(|(a, b)| a.is_finite() && b.is_finite() && **a > 0.0 && **b > 0.0)
        .map(|(a, b)| (a / b).ln())
        .collect();
    if logs.is_empty() {
        return f64::INFINITY;
    }
    logs.sort_by(|x, y| x.partial_cmp(y).unwrap());
    let c = 0.5 * (logs[0] + logs[logs.len() - 1]);
    logs.iter().map(|l| (l - c).abs()).fold(0.0, f64::max)
}

/// 採点側: 順序 bond 上の ŝ を真の物理 bond へ写像して v̂ を得る (perm は採点側の
/// 権利)。資格 (隣接 100%) を通ったセルでのみ呼ぶ。R-B は scale = 2 (物理 bond
/// あたり格子 bond 2 本 — セル内平均)、v̂ = a_lat/ŝ, a_lat = 1/scale。
pub fn vhat_on_true_bonds(
    order: &[usize],
    perm: &[usize],
    shat: &[f64],
    scale: usize,
    n_phys: usize,
    cyclic: bool,
) -> Vec<f64> {
    let n = order.len();
    let nb = shat.len();
    let nbp = if cyclic { n_phys } else { n_phys - 1 };
    let mut sum = vec![0.0; nbp];
    let mut cnt = vec![0usize; nbp];
    for u in 0..nb {
        if !shat[u].is_finite() || shat[u] <= 0.0 {
            continue;
        }
        let a = perm[order[u]];
        let b = perm[order[(u + 1) % n]];
        let dd = (a as isize - b as isize).unsigned_abs();
        let lat = if dd == 1 {
            a.min(b)
        } else if cyclic && dd == n - 1 {
            n - 1
        } else {
            continue;
        };
        let pb = lat / scale;
        if pb < nbp {
            sum[pb] += (1.0 / scale as f64) / shat[u];
            cnt[pb] += 1;
        }
    }
    (0..nbp)
        .map(|b| {
            if cnt[b] > 0 {
                sum[b] / cnt[b] as f64
            } else {
                f64::NAN
            }
        })
        .collect()
}

// ---- バー (凍結値) ----
// 導出規則 (本節に凍結): バー = val-4..7 実測の最大値 × 1.5 を小数 2 桁へ切り上げ。
// 導出走行 (バー = NaN) が実測と提案値を印字し、提案値を転記した検算走行を経て
// コミット = 凍結。バーは合格閾値であり、報告値は常に実測 — バーの余裕は主張の
// 強さを下げる方向にしか働かない (v29.3 の BAR_SXC_SPEARMAN と同じ手続き)。
/// Δ∞(v̂, v_true) の許容上限 (val 実測 max 0.2377 [val-7 R-B, 開放鎖境界] × 1.5 → 0.36)
pub const BAR_DINF_TRUE: f64 = 0.36;
/// 未使用源 τ 予言の max|ln(t*/τ_pred)| の許容上限 (val 実測 max 0.4998 × 1.5 → 0.75。
/// max 統計は d = D_PRED_MIN 近傍の近接場過渡が支配し弱い — 中央値 0.03–0.14 が
/// 分布の実態。統計量は導出走行前に宣言済みのため変更しない [val の二次使用禁止])
pub const BAR_TAU_PRED: f64 = 0.75;
/// regulator 間 Δ∞(v̂_R, v̂_R') の許容上限 (val 実測 max 0.1521 × 1.5 → 0.23)
pub const BAR_DINF_XREG: f64 = 0.23;

// =================== FROZEN QUANT v31 (END) ===================

// ---------------- 実験側 (val の 1 回走行 — 凍結節の外) ----------------

fn make_perm(n: usize, seed: u64) -> Vec<usize> {
    let mut rng = Rng::new(seed);
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = rng.range(i + 1);
        p.swap(i, j);
    }
    p
}

fn shuffle_state(st: &NodeState, perm: &[usize]) -> NodeState {
    let n = st.nodes;
    let d = st.dim();
    let mut cre = vec![0.0; d * d];
    let mut cim = vec![0.0; d * d];
    for i in 0..n {
        for a in 0..st.m {
            let gi_new = i * st.m + a;
            let gi_old = perm[i] * st.m + a;
            for j in 0..n {
                for b in 0..st.m {
                    let gj_new = j * st.m + b;
                    let gj_old = perm[j] * st.m + b;
                    cre[gi_new + gj_new * d] = st.cre[gi_old + gj_old * d];
                    cim[gi_new + gj_new * d] = st.cim[gi_old + gj_old * d];
                }
            }
        }
    }
    NodeState {
        nodes: n,
        m: st.m,
        cre,
        cim,
    }
}

fn adjacency_accuracy(edges: &[(usize, usize)], perm: &[usize], n: usize, ring: bool) -> f64 {
    if edges.is_empty() {
        return 0.0;
    }
    let truth = |a: usize, b: usize| {
        let d = (a as isize - b as isize).unsigned_abs();
        if ring {
            d.min(n - d) == 1
        } else {
            d == 1
        }
    };
    edges
        .iter()
        .filter(|&&(i, j)| truth(perm[i], perm[j]))
        .count() as f64
        / edges.len() as f64
}

/// 1 ジョブ = (instance, regulator)。資格 (2 質量 × S4 × C3) + 定量 (m=0) を計算。
struct JobOut {
    id: String,
    rname: &'static str,
    qual_fails: Vec<String>,
    min_sp: f64,
    quant_ok: bool,
    dinf_true: f64,
    tau_worst: f64,
    tau_med: f64,
    n_pairs: usize,
    coverage: f64,
    vhat: Vec<f64>, // 物理 bond 上の v̂ (xreg 用)
}

fn run_job(id: &str, seed: u64, rname: &'static str, reg: u8, scale: usize) -> JobOut {
    let spec = hold5_profile(seed);
    let mut qual_fails = Vec::new();
    let mut min_sp = 1.0f64;
    let mut quant_in: Option<(Vec<usize>, Vec<usize>, NodeState, Vec<f64>, Vec<f64>)> = None;
    for m_phys in [0.0, 0.05] {
        let (h, st) = hold5_system(&spec, reg, m_phys);
        let n = st.nodes;
        let d = st.dim();
        let perm = make_perm(n, seed ^ 0x5a5a ^ (reg as u64) ^ ((m_phys * 100.0) as u64));
        let shuf = shuffle_state(&st, &perm);
        let mut hp = vec![0.0; d * d];
        for i in 0..n {
            for a in 0..st.m {
                for j in 0..n {
                    for b in 0..st.m {
                        hp[(i * st.m + a) + (j * st.m + b) * d] =
                            h[(perm[i] * st.m + a) + (perm[j] * st.m + b) * d];
                    }
                }
            }
        }
        let spatial: [(&str, Vec<f64>); 4] = [
            ("B1", w_b1_gauss(&shuf)),
            ("B2", w_b2_gauss(&shuf)),
            ("B3-COV", w_b3_gauss(&shuf)),
            ("B5-QFI", w_b5_gauss(&shuf)),
        ];
        let sources = [0usize, n / 3, 2 * n / 3];
        let base: Vec<f64> = (0..d).map(|j| shuf.cre[j + j * d]).collect();
        let ts_den = fronts_gauss(&hp, &shuf, &base, &sources, FRONT_NT5);
        let mut ts_com: Vec<Vec<f64>> = Vec::new();
        let mut ts_kin: Vec<Vec<f64>> = Vec::new();
        for &q in &sources {
            let (mut tc_n, mut tk_n) = (vec![f64::INFINITY; n], vec![f64::INFINITY; n]);
            for a in 0..st.m {
                let (tc, tk) = fronts_commutator(&hp, &shuf.cre, d, q * st.m + a, FRONT_NT5);
                for node in 0..n {
                    for b in 0..st.m {
                        tc_n[node] = tc_n[node].min(tc[node * st.m + b]);
                        tk_n[node] = tk_n[node].min(tk[node * st.m + b]);
                    }
                }
            }
            ts_com.push(tc_n);
            ts_kin.push(tk_n);
        }
        for (sname, w) in &spatial {
            let r = reconstruct_v30(w, n);
            let topo_want = if spec.ring { 0u8 } else { 1u8 };
            let topo_ok = r.comps.len() == 1
                && r.comps[0].topology == topo_want
                && r.comps[0].order.len() == n
                && r.bridges.is_empty();
            let adj = adjacency_accuracy(&r.edges, &perm, n, spec.ring);
            if !(topo_ok && adj >= BAR_ADJACENCY) {
                qual_fails.push(format!(
                    "{} {} m={} {}: 空間資格 FAIL (topo {} adj {:.3})",
                    id, rname, m_phys, sname, topo_ok, adj
                ));
                continue;
            }
            let mut pos = vec![0usize; n];
            for (t, &nd) in r.comps[0].order.iter().enumerate() {
                pos[nd] = t;
            }
            let dist = |i: usize, j: usize| -> f64 {
                let dd = (pos[i] as isize - pos[j] as isize).unsigned_abs();
                if spec.ring {
                    dd.min(n - dd) as f64
                } else {
                    dd as f64
                }
            };
            for (cname, ts) in [
                ("DENSITY-FRONT", &ts_den),
                ("COMMUTATOR", &ts_com),
                ("KINEMATIC", &ts_kin),
            ] {
                let mut sps = Vec::new();
                for (si, &q) in sources.iter().enumerate() {
                    let dv: Vec<f64> = (0..n).map(|j| dist(q, j)).collect();
                    sps.push(spearman(&dv, &ts[si]));
                }
                let sp = sps.iter().sum::<f64>() / sps.len() as f64;
                min_sp = min_sp.min(sp);
                if sp < BAR_SXC_SPEARMAN {
                    qual_fails.push(format!(
                        "{} {} m={} {}×{}: Spearman {:.4} < {}",
                        id, rname, m_phys, sname, cname, sp, BAR_SXC_SPEARMAN
                    ));
                }
            }
            if m_phys == 0.0 && *sname == QUANT_S {
                quant_in = Some((
                    r.comps[0].order.clone(),
                    perm.clone(),
                    NodeState {
                        nodes: shuf.nodes,
                        m: shuf.m,
                        cre: shuf.cre.clone(),
                        cim: shuf.cim.clone(),
                    },
                    hp.clone(),
                    base.clone(),
                ));
            }
        }
    }
    // ---- 定量 (m = 0, S = B3-COV, C = DENSITY-FRONT) ----
    let (mut quant_ok, mut dinf_true, mut tau_worst, mut tau_med, mut n_pairs, mut coverage) =
        (false, f64::NAN, f64::NAN, f64::NAN, 0usize, 0.0f64);
    let mut vhat: Vec<f64> = Vec::new();
    if let Some((order, perm, shuf, hp, base)) = quant_in {
        let n = shuf.nodes;
        let mut posmap = vec![0usize; n];
        for (t, &nd) in order.iter().enumerate() {
            posmap[nd] = t;
        }
        let src_of_frac = |f: f64| -> usize { order[(f * n as f64) as usize % n] };
        let cal_nodes: Vec<usize> = CAL_FRACS.iter().map(|&f| src_of_frac(f)).collect();
        let pred_nodes: Vec<usize> = PRED_FRACS.iter().map(|&f| src_of_frac(f)).collect();
        let ts_cal = fronts_gauss_q(&hp, &shuf, &base, &cal_nodes, FRONT_DT_Q, FRONT_NT_Q);
        let ts_pred = fronts_gauss_q(&hp, &shuf, &base, &pred_nodes, FRONT_DT_Q, FRONT_NT_Q);
        let to_pos = |ts: &Vec<f64>| -> Vec<f64> {
            let mut out = vec![f64::INFINITY; n];
            for node in 0..n {
                out[posmap[node]] = ts[node];
            }
            out
        };
        let cal_data: Vec<(usize, Vec<f64>)> = cal_nodes
            .iter()
            .zip(ts_cal.iter())
            .map(|(&q, ts)| (posmap[q], to_pos(ts)))
            .collect();
        let shat = slowness_from_cal(spec.ring, &cal_data, n);
        coverage = shat.iter().filter(|s| s.is_finite()).count() as f64 / shat.len() as f64;
        // (a) 真値照合 (採点側)
        vhat = vhat_on_true_bonds(&order, &perm, &shat, scale, spec.n_phys, spec.ring);
        let nbp = vhat.len();
        let v_true: Vec<f64> = (0..nbp).map(|b| v_of_x(&spec, b as f64 + 0.5)).collect();
        dinf_true = delta_inf(&vhat, &v_true);
        // (b) 未使用源の τ 予言 (真値不使用 — 読み出し量のみ)
        let mut lnabs: Vec<f64> = Vec::new();
        for (si, &q) in pred_nodes.iter().enumerate() {
            let qp = posmap[q];
            let tv = to_pos(&ts_pred[si]);
            for p in 0..n {
                let dd = (p as isize - qp as isize).unsigned_abs();
                let d_hop = if spec.ring { dd.min(n - dd) } else { dd };
                if d_hop < D_PRED_MIN || !tv[p].is_finite() {
                    continue;
                }
                let tau = if spec.ring {
                    let t1 = tau_arc(&shat, n, true, qp, p, 1);
                    let t2 = tau_arc(&shat, n, true, qp, p, -1);
                    match (t1, t2) {
                        (Some(x), Some(y)) => Some(x.min(y)),
                        (Some(x), None) => Some(x),
                        (None, Some(y)) => Some(y),
                        _ => None,
                    }
                } else {
                    tau_arc(&shat, n, false, qp, p, if p > qp { 1 } else { -1 })
                };
                if let Some(tau) = tau {
                    if tau > 0.0 {
                        lnabs.push((tv[p] / tau).ln().abs());
                    }
                }
            }
        }
        n_pairs = lnabs.len();
        if n_pairs > 0 {
            lnabs.sort_by(|a, b| a.partial_cmp(b).unwrap());
            tau_worst = lnabs[n_pairs - 1];
            tau_med = lnabs[n_pairs / 2];
            quant_ok = true;
        }
    }
    JobOut {
        id: id.to_string(),
        rname,
        qual_fails,
        min_sp,
        quant_ok,
        dinf_true,
        tau_worst,
        tau_med,
        n_pairs,
        coverage,
        vhat,
    }
}

fn main() {
    self_test();
    println!("=== v29.4a 定量計量の器械確認とバー導出 — val 区画の 1 回使用 (第三十期) ===\n");
    let derivation_mode = BAR_DINF_TRUE.is_nan();
    if derivation_mode {
        println!("[モード] 導出走行 (バー未凍結 — val 実測から規則適用の提案値を印字)\n");
    }
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

    // ---- [Q0] SECRET 開示の検証 ----
    {
        let commit = sha256_hex(HOLD5_SECRET.as_bytes());
        let ok_commit =
            commit == "cfeb1959f3ba63b17c8ea7d1261f0f24570d33f7f8a2b3956a3df761bd301ab2";
        let ok_train = HOLD5_TRAIN_SEEDS
            .iter()
            .all(|&(seed, id)| hold5_seed(id) == seed);
        check(
            "[Q0] SECRET 開示 — sha256(SECRET) = v29.2 コミットメント / train seed = v29.3 公開定数",
            ok_commit && ok_train,
            format!("commitment {} / train 一致 {}", ok_commit, ok_train),
        );
        println!(
            "        val seed (本版で開示): {}",
            ["val-4", "val-5", "val-6", "val-7"]
                .iter()
                .map(|id| format!("{}={}", id, hold5_seed(id)))
                .collect::<Vec<_>>()
                .join(" ")
        );
        println!("        hold-8..15 の seed は v29.4b まで印字も生成もしない (凍結順序の維持)");
    }

    // ---- [Q1] 凍結節の照合 (v29 = v281 / v30 = v293) ----
    {
        let read_src = |p: &str| -> String {
            std::fs::read_to_string(p)
                .or_else(|_| std::fs::read_to_string(format!("../{}", p)))
                .unwrap_or_default()
        };
        let cut = |src: &str, beg: &str, end: &str| -> Option<String> {
            let b = src.find(beg)?;
            let e = src.find(end)?;
            Some(src[b..e].to_string())
        };
        let me = read_src("sim/src/bin/v294a_bridge_val.rs");
        let s281 = read_src("sim/src/bin/v281_bridge_train.rs");
        let s293 = read_src("sim/src/bin/v293_bridge_matrix.rs");
        const B29: &str = "// ===================== FROZEN BRIDGE READOUT v29 (BEGIN)";
        const E29: &str = "// ====================== FROZEN BRIDGE READOUT v29 (END)";
        const B30: &str = "// ================== FROZEN BRIDGE MATRIX v30 (BEGIN)";
        const E30: &str = "// =================== FROZEN BRIDGE MATRIX v30 (END)";
        let ok29 = match (cut(&s281, B29, E29), cut(&me, B29, E29)) {
            (Some(a), Some(b)) => {
                !a.is_empty() && sha256_hex(a.as_bytes()) == sha256_hex(b.as_bytes())
            }
            _ => false,
        };
        let (ok30, sha30) = match (cut(&s293, B30, E30), cut(&me, B30, E30)) {
            (Some(a), Some(b)) => (
                !a.is_empty() && sha256_hex(a.as_bytes()) == sha256_hex(b.as_bytes()),
                sha256_hex(a.as_bytes()),
            ),
            _ => (false, "?".into()),
        };
        check(
            "[Q1] 凍結節の逐語一致 (v29 = v281 / v30 = v293, SHA-256)",
            ok29 && ok30
                && sha30 == "93e1aa5cc0aed19ff8488fb8aa98f07ffd28ebda10cb61fd99b85f1a6e59fd8c",
            format!(
                "v29 {} / v30 {} (v30 = {}…)",
                ok29,
                ok30,
                &sha30[..16.min(sha30.len())]
            ),
        );
    }

    // ---- val 12 ジョブ (4 instance × 3 regulator) — 独立スレッド (決定性: 結果は
    //      添字順に回収・印字。ジョブ間に共有可変状態なし) ----
    let val_ids = ["val-4", "val-5", "val-6", "val-7"];
    let regs: [(&'static str, u8, usize); 3] = [("R-A", 0, 1), ("R-B", 1, 2), ("R-C", 2, 1)];
    let mut outs: Vec<Vec<JobOut>> = Vec::new();
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for id in val_ids {
            let seed = hold5_seed(id);
            let mut row = Vec::new();
            for &(rname, reg, scale) in &regs {
                row.push(s.spawn(move || run_job(id, seed, rname, reg, scale)));
            }
            handles.push(row);
        }
        for row in handles {
            outs.push(row.into_iter().map(|h| h.join().unwrap()).collect());
        }
    });

    // ---- [Q2] 資格 (v29.3 凍結採点器の健全性 — 新鮮な val 4 系) ----
    {
        let mut fails: Vec<String> = Vec::new();
        let mut min_sp = 1.0f64;
        for row in &outs {
            for o in row {
                fails.extend(o.qual_fails.iter().cloned());
                min_sp = min_sp.min(o.min_sp);
            }
        }
        check(
            "[Q2] val 4 系 × 3 regulator × 2 質量 × S(4) × C(3) = 288 セルの資格 (凍結バー)",
            fails.is_empty(),
            if fails.is_empty() {
                format!("全セル バー内 — 最小 Spearman = {:.4}", min_sp)
            } else {
                format!("バー外 {} 件: {:?}", fails.len(), &fails[..fails.len().min(6)])
            },
        );
    }

    // ---- [Q3] 定量実測とバー ----
    let (mut max_dinf, mut max_tau, mut max_xreg) = (0.0f64, 0.0f64, 0.0f64);
    let mut quant_all_ok = true;
    {
        println!("\n  -- 定量実測 (S = {}, C = {}, m = 0) --", QUANT_S, QUANT_C);
        println!("     instance reg    Δ∞(v̂,v_true)  τ予言max|ln| (中央値)   対数  被覆率");
        for row in &outs {
            for o in row {
                if !o.quant_ok {
                    quant_all_ok = false;
                    println!(
                        "     {} {}: 定量スキップ (資格 FAIL または対なし)",
                        o.id, o.rname
                    );
                    continue;
                }
                max_dinf = max_dinf.max(o.dinf_true);
                max_tau = max_tau.max(o.tau_worst);
                println!(
                    "     {:6} {:4}   {:.4}         {:.4}      ({:.4})   {:4}  {:.3}",
                    o.id, o.rname, o.dinf_true, o.tau_worst, o.tau_med, o.n_pairs, o.coverage
                );
            }
        }
        println!("\n     regulator 間 Δ∞(v̂, v̂′) (物理 bond 上, 共有 bond のみ):");
        for row in &outs {
            let mut parts = Vec::new();
            for i in 0..row.len() {
                for j in (i + 1)..row.len() {
                    if row[i].quant_ok && row[j].quant_ok {
                        let x = delta_inf(&row[i].vhat, &row[j].vhat);
                        max_xreg = max_xreg.max(x);
                        parts.push(format!("{}-{} {:.4}", row[i].rname, row[j].rname, x));
                    }
                }
            }
            println!("     {:6}  {}", row[0].id, parts.join("  "));
        }
        println!(
            "\n  [実測] max Δ∞(v̂,v_true) = {:.4} / max τ予言 = {:.4} / max regulator 間 = {:.4}",
            max_dinf, max_tau, max_xreg
        );
        if derivation_mode {
            let ceil2 = |x: f64| (x * 1.5 * 100.0).ceil() / 100.0;
            println!(
                "  [導出] 規則 (×1.5, 小数 2 桁切り上げ) の提案バー: BAR_DINF_TRUE = {:.2} / BAR_TAU_PRED = {:.2} / BAR_DINF_XREG = {:.2}",
                ceil2(max_dinf),
                ceil2(max_tau),
                ceil2(max_xreg)
            );
            println!("  [導出] 提案値を FROZEN QUANT v31 の定数へ転記し、検算走行の後にコミット (= 凍結) する");
        } else {
            check(
                "[Q3] val 定量実測が凍結バー内 (バー = 実測 × 1.5 規則 — 検算)",
                quant_all_ok
                    && max_dinf <= BAR_DINF_TRUE
                    && max_tau <= BAR_TAU_PRED
                    && max_xreg <= BAR_DINF_XREG,
                format!(
                    "バー: Δ∞true ≤ {} / τ ≤ {} / xreg ≤ {}",
                    BAR_DINF_TRUE, BAR_TAU_PRED, BAR_DINF_XREG
                ),
            );
        }
    }

    // ---- 成果物 ----
    {
        let frozen_sha = {
            let src = std::fs::read_to_string("sim/src/bin/v294a_bridge_val.rs")
                .or_else(|_| std::fs::read_to_string("../sim/src/bin/v294a_bridge_val.rs"))
                .unwrap_or_default();
            let b = src.find("// ================== FROZEN QUANT v31 (BEGIN)");
            let e = src.find("// =================== FROZEN QUANT v31 (END)");
            match (b, e) {
                (Some(b), Some(e)) => sha256_hex(src[b..e].as_bytes()),
                _ => "?".into(),
            }
        };
        println!(
            "\n[凍結] FROZEN QUANT v31 節 SHA-256 = {} (v29.4b が照合)",
            frozen_sha
        );
        if !derivation_mode {
            let mut cells = Vec::new();
            for row in &outs {
                for o in row {
                    cells.push(Json::Obj(vec![
                        ("id".into(), Json::Str(o.id.clone())),
                        ("reg".into(), Json::Str(o.rname.into())),
                        ("dinf_true".into(), Json::Num(o.dinf_true)),
                        ("tau_worst".into(), Json::Num(o.tau_worst)),
                        ("tau_med".into(), Json::Num(o.tau_med)),
                        ("n_pairs".into(), Json::Num(o.n_pairs as f64)),
                        ("coverage".into(), Json::Num(o.coverage)),
                    ]));
                }
            }
            let j = Json::Obj(vec![
                ("version".into(), Json::Str("v29.4a".into())),
                ("frozen_quant_sha256".into(), Json::Str(frozen_sha)),
                ("val_max_dinf_true".into(), Json::Num(max_dinf)),
                ("val_max_tau_pred".into(), Json::Num(max_tau)),
                ("val_max_xreg".into(), Json::Num(max_xreg)),
                ("bar_dinf_true".into(), Json::Num(BAR_DINF_TRUE)),
                ("bar_tau_pred".into(), Json::Num(BAR_TAU_PRED)),
                ("bar_dinf_xreg".into(), Json::Num(BAR_DINF_XREG)),
                ("cells".into(), Json::Arr(cells)),
            ]);
            let p = write_artifact("results/v294a_bridge_val.json", &j.render());
            println!("[成果物] {}", p);
        }
    }

    println!(
        "\n[判定] {}",
        if derivation_mode {
            "導出走行完了 — 提案バーを転記し検算走行へ (val の統計的使用はこの 1 回)"
        } else if nfail == 0 {
            "val は 1 回で閉じた — 定量バー凍結。holdout (hold-8..15) の初開封は v29.4b"
        } else {
            "val で破れ — 修正条項の手続き (holdout に進まない)"
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
