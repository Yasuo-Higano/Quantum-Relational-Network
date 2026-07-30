//! v29.3 合成 bridge — S×C 全組合せの訓練資格と新凍結節 (第三十期, PROMPT/11)
//!
//! 事前凍結 = bridge_candidates.yml 追補 v29.2 (コミット 5e48875 — 本実装より先):
//! B5-QFI / B6-COMMUTATOR の定義・S×C 独立セル採点 (フォールバック禁止)・
//! HOLD-5 生成器・seed コミットメント cfeb1959…。
//!
//! 本バイナリは **HOLD-5 の train インスタンス (train-0..3) だけ**を走らせる —
//! 設計・較正が許される区画。val-4..7 は v29.4 で 1 回、hold-8..15 は走行版まで
//! 生成すらしない (seed も本ファイルに置かない)。
//!
//! 構成:
//!   [凍結節 v29] B1/B2/B3-COV 核とパイプライン — v281 と逐語同一 (SHA-256 検査)
//!   [凍結節 v30] 新規: HOLD-5 生成器・B5-QFI (BKM 核)・B6-COMMUTATOR
//!     (retarded [n(t),n] + 運動学変種)・合成採点の定義とバー — 本コミットが凍結点
//!     (v29.4 は本節の SHA-256 一致を検査してから val/holdout を走らせる)
//!   [M0] 器械検査: Fock 構成 ρ_{S'S} = det(1−G)·det(K_{S'S}) の恒等検査
//!     (トレース 1・冪等 G→純粋・1 モード解析形)・BKM 重みの対称性・
//!     Wilson 実ゲージ (i^x) の分散照合
//!   [M1] train-0..3 × 3 regulator × 2 mass の S×C 資格表 (設計区画 — 表を見て
//!     バーを凍結し、本コミット後は変更しない)
//!
//! 正名: S ∈ {B1, B2, B3-COV, B5-QFI} (空間核) × C ∈ {DENSITY-FRONT, COMMUTATOR,
//! KINEMATIC} (因果チャネル)。各セルは自セルの S の再構成距離のみを使う。

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

// ---------------- 実験側 (train 資格の採点 — 凍結節の外) ----------------

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
    n: usize,
    ring: bool,
) -> f64 {
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

fn main() {
    self_test();
    println!("=== v29.3 合成 bridge — S×C 全組合せの訓練資格 (HOLD-5 train, 第三十期) ===\n");
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

    // ---- [M0] 器械検査 ----
    {
        // (a) Fock 構成の恒等: 1 モード G = [λ] → ρ = diag(1−λ, λ)
        let rho1 = fock_rho(&[0.3], 1);
        let ok_a = (rho1[0] - 0.7).abs() < 1e-12 && (rho1[3] - 0.3).abs() < 1e-12;
        // (b) 2 モードのトレース 1 と既知対角: G = diag(0.2, 0.6)
        let g2 = vec![0.2, 0.0, 0.0, 0.6];
        let r2 = fock_rho(&g2, 2);
        let tr: f64 = (0..4).map(|i| r2[i + i * 4]).sum();
        let ok_b = (tr - 1.0).abs() < 1e-12
            && (r2[0] - 0.8 * 0.4).abs() < 1e-12
            && (r2[15] - 0.2 * 0.6).abs() < 1e-12;
        // (c) 非対角 G の一体縮約が G に戻る: G ランダム実対称 (固有値内部)
        let g3 = vec![0.4, 0.15, 0.15, 0.55];
        let r3 = fock_rho(&g3, 2);
        // ⟨c†_0 c_0⟩ = ρ(01|01)+ρ(11|11) 型の対角 + 非対角 ⟨c†_0 c_1⟩ = ρ の 1 粒子ブロック
        let g00 = r3[1 + 1 * 4] + r3[3 + 3 * 4];
        let g01 = r3[1 + 2 * 4]; // ⟨c†_1 c_0⟩ 型 — 符号規約込みで |·| 照合
        let ok_c = (g00 - 0.4).abs() < 1e-10 && (g01.abs() - 0.15).abs() < 1e-10;
        // (d) BKM 重みの対称性・退化極限
        let ok_d = (bkm_weight(0.3, 0.3) - 0.3).abs() < 1e-12
            && (bkm_weight(0.2, 0.4) - bkm_weight(0.4, 0.2)).abs() < 1e-15;
        check(
            "[M0a] Fock 構成 ρ = det(1−G)·det(K[S'|S]) の恒等検査 (1/2 モード・縮約・BKM 重み)",
            ok_a && ok_b && ok_c && ok_d,
            format!("1モード {} / trace&対角 {} / 縮約 {} / BKM {}", ok_a, ok_b, ok_c, ok_d),
        );
        // (e) Wilson 実ゲージの分散照合: 一様 v = 1, m = 0.05 の R-C スペクトルが
        //     解析分散 ±√(sin²k + (m + 1 − cos k)²) と一致
        let spec_u = Hold5Spec {
            n_phys: 32,
            ring: true,
            amps: vec![],
            centers: vec![],
            widths: vec![],
        };
        let (hw, _) = hold5_system(&spec_u, 2, 0.05);
        let (mut evw, _) = jacobi_eigh(&hw, 64);
        evw.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut want = Vec::new();
        for kk in 0..32 {
            let k = 2.0 * std::f64::consts::PI * kk as f64 / 32.0;
            let e = ((k.sin()).powi(2) + (0.05 + 1.0 - k.cos()).powi(2)).sqrt();
            want.push(e);
            want.push(-e);
        }
        want.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let maxdev = evw
            .iter()
            .zip(want.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        check(
            "[M0b] Wilson 実ゲージ構成の分散照合 (一様極限 = 解析 ±E(k))",
            maxdev < 1e-10,
            format!("max|ΔE| = {:.2e}", maxdev),
        );
    }

    // ---- [M1] train-0..3 × R-A/B/C × m ∈ {0, 0.05} の S×C 資格表 ----
    let mut all_ok = true;
    let mut min_spearman = 1.0f64;
    let mut table: Vec<String> = Vec::new();
    for &(seed, id) in HOLD5_TRAIN_SEEDS.iter() {
        let spec = hold5_profile(seed);
        for (rname, reg) in [("R-A", 0u8), ("R-B", 1), ("R-C", 2)] {
            for m_phys in [0.0, 0.05] {
                let (h, st) = hold5_system(&spec, reg, m_phys);
                let n = st.nodes;
                let perm = make_perm(n, seed ^ 0x5a5a ^ (reg as u64) ^ ((m_phys * 100.0) as u64));
                let shuf = shuffle_state(&st, &perm);
                // H も同じラベルで (状態族の生成は実験側)
                let d = st.dim();
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
                // 空間核 4 種
                let spatial: [(&str, Vec<f64>); 4] = [
                    ("B1", w_b1_gauss(&shuf)),
                    ("B2", w_b2_gauss(&shuf)),
                    ("B3-COV", w_b3_gauss(&shuf)),
                    ("B5-QFI", w_b5_gauss(&shuf)),
                ];
                // 因果チャネル: 源 = 置換後 id {0, n/3, 2n/3}
                let sources = [0usize, n / 3, 2 * n / 3];
                let base: Vec<f64> = (0..st.dim()).map(|j| shuf.cre[j + j * st.dim()]).collect();
                // DENSITY-FRONT は既存 fronts_gauss (ノード単位) — ここではノード=モード群
                let ts_den = fronts_gauss(&hp, &shuf, &base, &sources, FRONT_NT5);
                // COMMUTATOR/KINEMATIC はモード単位 → ノード集約 (ノード内最小 t*)
                let mut ts_com: Vec<Vec<f64>> = Vec::new();
                let mut ts_kin: Vec<Vec<f64>> = Vec::new();
                for &q in &sources {
                    let (mut tc_n, mut tk_n) = (vec![f64::INFINITY; n], vec![f64::INFINITY; n]);
                    for a in 0..st.m {
                        let (tc, tk) =
                            fronts_commutator(&hp, &shuf.cre, d, q * st.m + a, FRONT_NT5);
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
                    let s_ok = topo_ok && adj >= BAR_ADJACENCY;
                    if !s_ok {
                        all_ok = false;
                        table.push(format!(
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
                        if spec.ring { dd.min(n - dd) as f64 } else { dd as f64 }
                    };
                    for (cname, ts) in [
                        ("DENSITY-FRONT", &ts_den),
                        ("COMMUTATOR", &ts_com),
                        ("KINEMATIC", &ts_kin),
                    ] {
                        let mut sps = Vec::new();
                        for (si, &q) in sources.iter().enumerate() {
                            let dv: Vec<f64> = (0..n).map(|j| dist(q, j)).collect();
                            let tv: Vec<f64> = ts[si].clone();
                            sps.push(spearman(&dv, &tv));
                        }
                        let sp = sps.iter().sum::<f64>() / sps.len() as f64;
                        min_spearman = min_spearman.min(sp);
                        if sp < BAR_SXC_SPEARMAN {
                            all_ok = false;
                            table.push(format!(
                                "{} {} m={} {}×{}: Spearman {:.4} < {}",
                                id, rname, m_phys, sname, cname, sp, BAR_SXC_SPEARMAN
                            ));
                        }
                    }
                }
            }
        }
    }
    check(
        "[M1] train 4 系 × 3 regulator × 2 質量 × S(4) × C(3) = 288 セルの資格 (隣接 100%・Spearman ≥ 0.90)",
        all_ok,
        if all_ok {
            format!("全セル バー内 — 最小 Spearman = {:.4}", min_spearman)
        } else {
            format!("バー外 {} 件: {:?}", table.len(), &table[..table.len().min(6)])
        },
    );

    // ---- 成果物 ----
    {
        let frozen_sha = {
            let src = std::fs::read_to_string("sim/src/bin/v293_bridge_matrix.rs")
                .or_else(|_| std::fs::read_to_string("../sim/src/bin/v293_bridge_matrix.rs"))
                .unwrap_or_default();
            let b = src.find("// ================== FROZEN BRIDGE MATRIX v30 (BEGIN)");
            let e = src.find("// =================== FROZEN BRIDGE MATRIX v30 (END)");
            match (b, e) {
                (Some(b), Some(e)) => sha256_hex(src[b..e].as_bytes()),
                _ => "?".into(),
            }
        };
        println!("\n[凍結] FROZEN BRIDGE MATRIX v30 節 SHA-256 = {} (v29.4 が照合)", frozen_sha);
        let j = Json::Obj(vec![
            ("version".into(), Json::Str("v29.3".into())),
            ("frozen_matrix_sha256".into(), Json::Str(frozen_sha)),
            ("min_spearman_train".into(), Json::Num(min_spearman)),
            ("bar_sxc_spearman".into(), Json::Num(BAR_SXC_SPEARMAN)),
        ]);
        let p = write_artifact("results/v293_bridge_matrix.json", &j.render());
        println!("[成果物] {}", p);
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "S×C 全組合せ (B5-QFI/B6-COMMUTATOR 込み) が train で資格 — 本コミットが凍結点、val (1 回) と holdout は v29.4"
        } else {
            "資格 FAIL — 凍結前の設計区画なので原因を particularに記録して再設計する"
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
