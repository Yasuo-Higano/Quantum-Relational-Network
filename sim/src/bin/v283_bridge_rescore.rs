//! v28.3 bridge holdout の修正採点 (第二十九期 III 前半) — 修正条項 1/2 の執行
//!
//! 事前凍結 = docs/uft-v28.2.md §3 (コミット 656d4cc — 本バイナリの実装・走行より先):
//!   修正条項 1: ξ 推定 = 距離 d ∈ {3,…,20} ごとの上包絡 ln(max w) の線形フィット
//!               (副格子盲目化 — staggered の偶距離 C≡0 チャネル混入を両系対称に排除)
//!   修正条項 2: HOLD-4 開放鎖の距離供給 = HOLD-1 通過の空間候補で番号最小 (B1→B2→B3)
//!   不変条項: 凍結節 (核・パイプライン・全バー) は v28.1 のまま。候補の調整なし。
//!
//! 検査: [R0] 凍結節 SHA-256 = v281 (読み出し規則の無変更) + TRAIN-1 アンカー
//!       [R1] 修正と無関係な全セル (H1×3, H2×5, H4 リング) の裁定が第一採点
//!            (v28.2 as-run) と一致 — 修正が波及していないことの回帰
//!       [R2] HOLD-3 の修正採点 (修正条項 1) — バーは不変 (ξ 比 ∈ [0.8, 1.25])
//!       [R3] HOLD-4 開放鎖の修正採点 (修正条項 2) — バー不変 (R² ≥ 0.98)
//! 最終成績表 = 第二十九期の bridge 比較の確定表 (v28.2 の第一採点表は保存のまま)。
//!
//! [FAIL] は器械破れ (R0/R1) のみ。候補の裁定は [採点] 行。

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


// ---------------- 実験側 (v282 と同一の系構成 + 修正条項の採点) ----------------

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

fn adjacency_accuracy(
    edges: &[(usize, usize)],
    perm: &[usize],
    truth: &dyn Fn(usize, usize) -> bool,
) -> f64 {
    if edges.is_empty() {
        return 0.0;
    }
    let ok = edges
        .iter()
        .filter(|&&(i, j)| truth(perm[i], perm[j]))
        .count();
    ok as f64 / edges.len() as f64
}

fn gs_correlation(h: &[f64], n: usize, nocc: usize) -> Vec<f64> {
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

fn hold1_system(n: usize) -> (Vec<f64>, NodeState) {
    let mut h = vec![0.0; n * n];
    for x in 0..n - 1 {
        let t = 1.0 + 0.15 * (-((x as f64 - 50.0) / 12.0).powi(2)).exp();
        h[x + (x + 1) * n] = -t;
        h[(x + 1) + x * n] = -t;
    }
    let c = gs_correlation(&h, n, n / 2);
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

fn hold3_staggered(n: usize, m: f64) -> NodeState {
    let mut h = vec![0.0; n * n];
    for x in 0..n {
        let y = (x + 1) % n;
        h[x + y * n] = -0.5;
        h[y + x * n] = -0.5;
        h[x + x * n] = if x % 2 == 0 { m } else { -m };
    }
    let c = gs_correlation(&h, n, n / 2);
    NodeState {
        nodes: n,
        m: 1,
        cre: c,
        cim: vec![0.0; n * n],
    }
}

fn hold3_wilson(n: usize, m: f64) -> NodeState {
    let d = 2 * n;
    let mut cre = vec![0.0; d * d];
    let mut cim = vec![0.0; d * d];
    let two_pi = 2.0 * std::f64::consts::PI;
    for kk in 0..n {
        let k = two_pi * kk as f64 / n as f64;
        let hx = k.sin();
        let hz = m + (1.0 - k.cos());
        let e = (hx * hx + hz * hz).sqrt().max(1e-300);
        let p = [
            [(1.0 - hz / e) * 0.5, -hx / e * 0.5],
            [-hx / e * 0.5, (1.0 + hz / e) * 0.5],
        ];
        for x in 0..n {
            for y in 0..n {
                let ph = k * (x as f64 - y as f64);
                let (cph, sph) = (ph.cos() / n as f64, ph.sin() / n as f64);
                for a in 0..2 {
                    for b in 0..2 {
                        let gi = x * 2 + a;
                        let gj = y * 2 + b;
                        cre[gi + gj * d] += p[a][b] * cph;
                        cim[gi + gj * d] += p[a][b] * sph;
                    }
                }
            }
        }
    }
    NodeState {
        nodes: n,
        m: 2,
        cre,
        cim,
    }
}

/// 修正条項 1: ξ = 距離 d ∈ {3,…,20} ごとの上包絡 ln(max_{対 at d} w) の線形フィット
fn xi_envelope(w: &[f64], n: usize, order: &[usize], cyclic: bool) -> f64 {
    let mut pos = vec![0usize; n];
    for (t, &node) in order.iter().enumerate() {
        pos[node] = t;
    }
    let dist = |i: usize, j: usize| -> usize {
        let d = (pos[i] as isize - pos[j] as isize).unsigned_abs();
        if cyclic {
            d.min(n - d)
        } else {
            d
        }
    };
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for dd in 3..=20usize {
        let mut mx = f64::NEG_INFINITY;
        for i in 0..n {
            for j in (i + 1)..n {
                if dist(i, j) == dd {
                    mx = mx.max(w[i + j * n]);
                }
            }
        }
        if mx > 0.0 && mx.is_finite() {
            xs.push(dd as f64);
            ys.push(mx.ln());
        }
    }
    let (slope, _) = linfit(&xs, &ys);
    -1.0 / slope
}

fn front_fit(pairs: &[(f64, f64)]) -> (f64, f64) {
    let (mut sxy, mut sxx) = (0.0, 0.0);
    for &(d, t) in pairs {
        sxy += d * t;
        sxx += t * t;
    }
    let v = sxy / sxx.max(1e-300);
    let mean_d = pairs.iter().map(|p| p.0).sum::<f64>() / pairs.len().max(1) as f64;
    let (mut ss_res, mut ss_tot) = (0.0, 0.0);
    for &(d, t) in pairs {
        ss_res += (d - v * t) * (d - v * t);
        ss_tot += (d - mean_d) * (d - mean_d);
    }
    (v, 1.0 - ss_res / ss_tot.max(1e-300))
}

fn main() {
    self_test();
    println!("=== v28.3 bridge holdout の修正採点 — 修正条項 1/2 の執行 (第二十九期 III 前半) ===\n");
    let mut n_instr_fail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!(
            "  [{}] {}  {}",
            if ok { "PASS" } else { "FAIL" },
            name,
            detail
        );
        if !ok {
            n_instr_fail += 1;
        }
    };
    let mut score: BTreeMap<&str, Vec<(String, bool)>> = BTreeMap::new();
    let mut judge = |score: &mut BTreeMap<&str, Vec<(String, bool)>>,
                     cand: &'static str,
                     sys: &str,
                     bar_ok: bool,
                     detail: String| {
        println!(
            "  [採点] {} × {} → {}  {}",
            cand,
            sys,
            if bar_ok { "バー内" } else { "バー外 → 棄却" },
            detail
        );
        score.entry(cand).or_default().push((sys.to_string(), bar_ok));
    };

    // ---- [R0] 凍結節 = v281 (無変更) ----
    {
        let read_src = |p: &str| -> String {
            std::fs::read_to_string(p)
                .or_else(|_| std::fs::read_to_string(format!("../{}", p)))
                .unwrap_or_default()
        };
        let cut = |src: &str| -> Option<String> {
            let b = src.find("// ===================== FROZEN BRIDGE READOUT v29 (BEGIN)")?;
            let e = src.find("// ====================== FROZEN BRIDGE READOUT v29 (END)")?;
            Some(src[b..e].to_string())
        };
        let s281 = read_src("sim/src/bin/v281_bridge_train.rs");
        let s283 = read_src("sim/src/bin/v283_bridge_rescore.rs");
        let sha_ok = match (cut(&s281), cut(&s283)) {
            (Some(a), Some(b)) => sha256_hex(a.as_bytes()) == sha256_hex(b.as_bytes()),
            _ => false,
        };
        check(
            "[R0] 凍結節 SHA-256 = v281 (核・パイプライン・バーの無変更)",
            sha_ok,
            format!("一致 {}", sha_ok),
        );
    }

    // ---- 無関係セルの再計算 (H1, H2, H4 リング — 裁定は v28.2 as-run と一致すること) ----
    // v28.2 の第一採点 (results/v282_bridge_holdout.txt): B1-H1 外 / B2-H1 内 / B3-H1 内 /
    // H2: B1 内・B2 内・B3 内・未使用 内 (0.9705)・B4 内 (0.8857) / H4 リング 内 (R² 0.99876)
    let expected_run1: [(&str, bool); 9] = [
        ("B1-H1", false),
        ("B2-H1", true),
        ("B3-H1", true),
        ("B1-H2", true),
        ("B2-H2", true),
        ("B3-H2", true),
        ("UNUSED-H2", true),
        ("B4-H2", true),
        ("B4-H4RING", true),
    ];
    let mut run2: Vec<(&str, bool)> = Vec::new();

    // H1
    let n1 = 101;
    let (h1_h, h1_st) = hold1_system(n1);
    let h1_perm = make_perm(n1, 31415);
    let h1_shuf = shuffle_state(&h1_st, &h1_perm);
    let h1_truth = |a: usize, b: usize| (a as isize - b as isize).unsigned_abs() == 1;
    let mut h1_orders: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    let mut h1_pass: BTreeMap<&str, bool> = BTreeMap::new();
    for (name, w) in [
        ("B1", w_b1_gauss(&h1_shuf)),
        ("B2", w_b2_gauss(&h1_shuf)),
        ("B3", w_b3_gauss(&h1_shuf)),
    ] {
        let r = reconstruct(&w, n1);
        let path_ok =
            r.comps.len() == 1 && r.comps[0].topology == 1 && r.comps[0].order.len() == n1;
        let adj = adjacency_accuracy(&r.edges, &h1_perm, &h1_truth);
        let ok = path_ok && adj >= BAR_ADJACENCY && r.bridges.is_empty();
        run2.push((
            match name {
                "B1" => "B1-H1",
                "B2" => "B2-H1",
                _ => "B3-H1",
            },
            ok,
        ));
        h1_pass.insert(name, ok);
        if path_ok {
            h1_orders.insert(name, r.comps[0].order.clone());
        }
    }

    // H2
    let l = 14;
    let ring2 = Z2GaugeRing::try_new(l, 7, 1.0, 0.6, 0.2, vec![]).expect("Z2 l=14");
    let (_e0, gs2, res2) = ring2.ground_state(11);
    check(
        "[R0b] Z2 基底状態の器械 (Lanczos 収束)",
        res2 < 1e-6,
        format!("残差 {:.1e}", res2),
    );
    let perm2 = make_perm(l, 2718);
    let view2 = Z2View {
        l,
        masks: &gs2.masks,
        psi: &gs2.psi,
        to_orig: perm2.clone(),
    };
    let truth2 = |a: usize, b: usize| {
        let d = (a as isize - b as isize).unsigned_abs();
        d.min(l - d) == 1
    };
    let mut b1_order2: Vec<usize> = Vec::new();
    for (name, w) in [
        ("B1", w_b1_z2(&view2)),
        ("B2", w_b2_z2(&view2)),
        ("B3", w_b3_z2(&view2)),
    ] {
        let r = reconstruct(&w, l);
        let cyc = r.comps.len() == 1 && r.comps[0].topology == 0 && r.comps[0].order.len() == l;
        let adj = adjacency_accuracy(&r.edges, &perm2, &truth2);
        let ok = cyc && adj >= BAR_ADJACENCY && r.bridges.is_empty();
        run2.push((
            match name {
                "B1" => "B1-H2",
                "B2" => "B2-H2",
                _ => "B3-H2",
            },
            ok,
        ));
        if name == "B1" && cyc {
            b1_order2 = r.comps[0].order.clone();
        }
    }
    {
        let mut pos = vec![0usize; l];
        for (t, &nd) in b1_order2.iter().enumerate() {
            pos[nd] = t;
        }
        let wb3 = w_b3_z2(&view2);
        let mut xs = Vec::new();
        let mut ys = Vec::new();
        for i in 0..l {
            for j in (i + 1)..l {
                let dd = (pos[i] as isize - pos[j] as isize).unsigned_abs();
                let d = dd.min(l - dd);
                if (1..=5).contains(&d) {
                    xs.push(-(d as f64));
                    ys.push(wb3[i + j * l]);
                }
            }
        }
        let sp = spearman(&xs, &ys);
        run2.push(("UNUSED-H2", !b1_order2.is_empty() && sp >= BAR_Z2_UNUSED));
        // B4 (H2)
        let base: Vec<f64> = (0..l).map(|s| gs2.density(s)).collect();
        let mut s = ring2.apply_bond_op(&gs2, 0);
        let nrm = cdot(&s.psi, &s.psi).0.sqrt();
        for x in s.psi.iter_mut() {
            x.0 /= nrm;
            x.1 /= nrm;
        }
        let mut tstar = vec![f64::INFINITY; l];
        for it in 1..=100 {
            s = ring2.step(&s, EvolutionParameter(0.1));
            let t = 0.1 * it as f64;
            for site in 0..l {
                if tstar[site].is_finite() {
                    continue;
                }
                if (s.density(site) - base[site]).abs() >= FRONT_EPS {
                    tstar[site] = t;
                }
            }
        }
        let src0 = (0..l).find(|&i| perm2[i] == 0).unwrap();
        let src1 = (0..l).find(|&i| perm2[i] == 1).unwrap();
        let d_by_new: Vec<f64> = (0..l)
            .map(|i| {
                let d0 = {
                    let dd = (pos[i] as isize - pos[src0] as isize).unsigned_abs();
                    dd.min(l - dd)
                };
                let d1 = {
                    let dd = (pos[i] as isize - pos[src1] as isize).unsigned_abs();
                    dd.min(l - dd)
                };
                d0.min(d1) as f64
            })
            .collect();
        let t_by_new: Vec<f64> = (0..l).map(|i| tstar[perm2[i]]).collect();
        let sp4 = spearman(&d_by_new, &t_by_new);
        run2.push(("B4-H2", !b1_order2.is_empty() && sp4 >= BAR_Z2_SPEARMAN));
    }

    // H4 リング
    {
        let n = 202;
        let ring = RingChain { n };
        let g = ring.init();
        let st = NodeState {
            nodes: n,
            m: 1,
            cre: g.cre.clone(),
            cim: g.cim.clone(),
        };
        let perm = make_perm(n, 20260729);
        let shuf = shuffle_state(&st, &perm);
        let w = w_b1_gauss(&shuf);
        let r = reconstruct(&w, n);
        let cyc = r.comps.len() == 1 && r.comps[0].topology == 0;
        let mut pos = vec![0usize; n];
        if cyc {
            for (t, &nd) in r.comps[0].order.iter().enumerate() {
                pos[nd] = t;
            }
        }
        let mut h = vec![0.0; n * n];
        for x in 0..n {
            let y = (x + 1) % n;
            h[x + y * n] = -1.0;
            h[y + x * n] = -1.0;
        }
        let mut hp = vec![0.0; n * n];
        for i in 0..n {
            for j in 0..n {
                hp[i + j * n] = h[perm[i] + perm[j] * n];
            }
        }
        let sources: Vec<usize> = (0..8).map(|k| k * n / 8).collect();
        let base: Vec<f64> = (0..n).map(|j| shuf.cre[j + j * n]).collect();
        let ts = fronts_gauss(&hp, &shuf, &base, &sources, 130);
        let mut all_pairs: Vec<(f64, f64)> = Vec::new();
        let mut vs = Vec::new();
        for (si, &q) in sources.iter().enumerate() {
            let mut pairs = Vec::new();
            for j in 0..n {
                let dd = (pos[j] as isize - pos[q] as isize).unsigned_abs();
                let d = dd.min(n - dd) as f64;
                let t = ts[si][j];
                if t.is_finite() && (3.0..=80.0).contains(&d) {
                    pairs.push((d, t));
                }
            }
            let (v, _) = front_fit(&pairs);
            vs.push(v);
            all_pairs.extend(pairs);
        }
        let (_v_all, r2) = front_fit(&all_pairs);
        let vmean = vs.iter().sum::<f64>() / vs.len() as f64;
        let vspread = vs
            .iter()
            .map(|v| (v / vmean - 1.0).abs())
            .fold(0.0, f64::max);
        run2.push((
            "B4-H4RING",
            cyc && r2 >= BAR_H4_R2 && vspread <= BAR_H4_VSPREAD,
        ));
    }

    // ---- [R1] 無関係セルの裁定一致 ----
    {
        let mut bad = Vec::new();
        for (name, want) in expected_run1 {
            match run2.iter().find(|(n, _)| *n == name) {
                Some((_, got)) if *got == want => {}
                Some((_, got)) => bad.push(format!("{}: run1={} run2={}", name, want, got)),
                None => bad.push(format!("{}: 再計算なし", name)),
            }
        }
        check(
            "[R1] 修正と無関係な 9 セルの裁定が第一採点 (v28.2 as-run) と一致",
            bad.is_empty(),
            if bad.is_empty() {
                "修正条項は波及していない".into()
            } else {
                format!("{:?}", bad)
            },
        );
        // 成績表へ転記 (無関係セルは第一採点の裁定のまま)
        judge(&mut score, "B1", "HOLD-1 開放端鎖", h1_pass["B1"], "第一採点どおり (橋誤検出 — 候補固有)".into());
        judge(&mut score, "B2", "HOLD-1 開放端鎖", h1_pass["B2"], "第一採点どおり".into());
        judge(&mut score, "B3", "HOLD-1 開放端鎖", h1_pass["B3"], "第一採点どおり".into());
        for (cand, key) in [("B1", "B1-H2"), ("B2", "B2-H2"), ("B3", "B3-H2")] {
            let ok = run2.iter().find(|(n, _)| *n == key).unwrap().1;
            judge(&mut score, match cand { "B1" => "B1", "B2" => "B2", _ => "B3" }, "HOLD-2 Z2 (相互作用)", ok, "第一採点どおり".into());
        }
        let un = run2.iter().find(|(n, _)| *n == "UNUSED-H2").unwrap().1;
        judge(&mut score, "B1", "HOLD-2 未使用チャネル", un, "第一採点どおり".into());
        let b4h2 = run2.iter().find(|(n, _)| *n == "B4-H2").unwrap().1;
        judge(&mut score, "B4", "HOLD-2 因果 (相互作用系)", b4h2, "第一採点どおり".into());
        let b4h4 = run2.iter().find(|(n, _)| *n == "B4-H4RING").unwrap().1;
        judge(&mut score, "B4", "HOLD-4 リング (Lorentzian)", b4h4, "第一採点どおり".into());
    }

    // ---- [R2] HOLD-3 の修正採点 (修正条項 1 — 上包絡 ξ) ----
    {
        let n = 402;
        let m_phys = 0.05;
        let st_s = hold3_staggered(n, m_phys);
        let st_w = hold3_wilson(n, m_phys);
        let perm_s = make_perm(n, 161803);
        let perm_w = make_perm(n, 141421);
        let shuf_s = shuffle_state(&st_s, &perm_s);
        let shuf_w = shuffle_state(&st_w, &perm_w);
        let truth = |a: usize, b: usize| {
            let d = (a as isize - b as isize).unsigned_abs();
            d.min(n - d) == 1
        };
        println!(
            "  (導出参考: κ_s = arsinh(ma) = {:.6} / κ_w ≈ ln(1+ma) = {:.6} — 期待比 ~{:.3})",
            (0.05f64).asinh(),
            (1.05f64).ln(),
            (1.05f64).ln() / (0.05f64).asinh()
        );
        for (name, ws, ww) in [
            ("B1", w_b1_gauss(&shuf_s), w_b1_gauss(&shuf_w)),
            ("B2", w_b2_gauss(&shuf_s), w_b2_gauss(&shuf_w)),
            ("B3", w_b3_gauss(&shuf_s), w_b3_gauss(&shuf_w)),
        ] {
            let mut oks = [false; 2];
            let mut xis = [0.0f64; 2];
            for (t, (w, perm)) in [(&ws, &perm_s), (&ww, &perm_w)].into_iter().enumerate() {
                let r = reconstruct(w, n);
                let cyc =
                    r.comps.len() == 1 && r.comps[0].topology == 0 && r.comps[0].order.len() == n;
                let adj = adjacency_accuracy(&r.edges, perm, &truth);
                oks[t] = cyc && adj >= BAR_ADJACENCY && r.bridges.is_empty();
                if cyc {
                    xis[t] = xi_envelope(w, n, &r.comps[0].order, true);
                }
            }
            let ratio = if xis[1] != 0.0 { xis[0] / xis[1] } else { 0.0 };
            let ratio_ok = ratio >= BAR_H3_XI_RATIO.0 && ratio <= BAR_H3_XI_RATIO.1;
            judge(
                &mut score,
                match name {
                    "B1" => "B1",
                    "B2" => "B2",
                    _ => "B3",
                },
                "HOLD-3 二離散化 (修正条項 1)",
                oks[0] && oks[1] && ratio_ok,
                format!(
                    "stag ξ={:.2} / Wil ξ={:.2} / ξ比 {:.3} (窓 [{}, {}] — バー不変)",
                    xis[0], xis[1], ratio, BAR_H3_XI_RATIO.0, BAR_H3_XI_RATIO.1
                ),
            );
        }
    }

    // ---- [R3] HOLD-4 開放鎖の修正採点 (修正条項 2 — 距離供給 B1→B2→B3) ----
    {
        let supplier = ["B1", "B2", "B3"]
            .iter()
            .find(|c| *h1_pass.get(*c).unwrap_or(&false))
            .copied();
        let (ok, detail) = if let Some(sup) = supplier {
            let order = &h1_orders[sup];
            let mut pos1 = vec![0usize; n1];
            for (t, &nd) in order.iter().enumerate() {
                pos1[nd] = t;
            }
            let mut hp1 = vec![0.0; n1 * n1];
            for i in 0..n1 {
                for j in 0..n1 {
                    hp1[i + j * n1] = h1_h[h1_perm[i] + h1_perm[j] * n1];
                }
            }
            let src1: Vec<usize> = [10usize, 50, 90]
                .iter()
                .map(|&o| (0..n1).find(|&i| h1_perm[i] == o).unwrap())
                .collect();
            let base1: Vec<f64> = (0..n1).map(|j| h1_shuf.cre[j + j * n1]).collect();
            let ts1 = fronts_gauss(&hp1, &h1_shuf, &base1, &src1, 130);
            let mut pairs1: Vec<(f64, f64)> = Vec::new();
            for (si, &q) in src1.iter().enumerate() {
                for j in 0..n1 {
                    let d = (pos1[j] as isize - pos1[q] as isize).unsigned_abs() as f64;
                    let t = ts1[si][j];
                    if t.is_finite() && d >= 3.0 {
                        pairs1.push((d, t));
                    }
                }
            }
            let (v, r2) = front_fit(&pairs1);
            (
                r2 >= BAR_H4_R2,
                format!(
                    "距離供給 = {} / v = {:.4} R² = {:.5} ({} 点, バー R² ≥ {})",
                    sup,
                    v,
                    r2,
                    pairs1.len(),
                    BAR_H4_R2
                ),
            )
        } else {
            (false, "供給候補なし (H1 全滅)".into())
        };
        judge(&mut score, "B4", "HOLD-4 開放鎖 (修正条項 2)", ok, detail);
    }

    // ---- 確定成績表 ----
    println!("\n---- 確定成績 (第一採点 + 修正条項 — 候補・核・バーは v28.1 凍結のまま) ----");
    let mut survived: Vec<String> = Vec::new();
    let mut rejected: Vec<String> = Vec::new();
    for (cand, rows) in &score {
        let all_ok = rows.iter().all(|(_, ok)| *ok);
        let summary: Vec<String> = rows
            .iter()
            .map(|(s, ok)| format!("{}:{}", s, if *ok { "内" } else { "外" }))
            .collect();
        println!(
            "  {} — {} → **{}**",
            cand,
            summary.join(" / "),
            if all_ok { "生存" } else { "棄却" }
        );
        if all_ok {
            survived.push(cand.to_string());
        } else {
            rejected.push(cand.to_string());
        }
    }

    {
        let j = Json::Obj(vec![
            ("version".into(), Json::Str("v28.3".into())),
            (
                "survived".into(),
                Json::Arr(survived.iter().map(|s| Json::Str(s.clone())).collect()),
            ),
            (
                "rejected".into(),
                Json::Arr(rejected.iter().map(|s| Json::Str(s.clone())).collect()),
            ),
        ]);
        let p = write_artifact("results/v283_bridge_rescore.json", &j.render());
        println!("\n[成果物] {}", p);
    }

    println!(
        "\n[判定] 確定: 生存 {:?} / 棄却 {:?}。生存候補が示すのは「同じ網から時空 (的構造) が出る」toy 機構 (C3) — bridge law 登録簿は空のまま (登録には 2 模型 × 2 regulator を超える一意性・Lorentzian 完全性・未使用チャネルの事前予言の常設化が要る)。",
        survived, rejected
    );
    println!(
        "\n総合判定: {}",
        if n_instr_fail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if n_instr_fail > 0 {
        std::process::exit(1);
    }
}
