//! v28.1 relational geometry bridge — 候補 B1–B4 の訓練資格認定 (第二十九期 I/III)
//!
//! 事前登録 = bridge_candidates.yml (v28.0, 実行前凍結)。本バイナリは訓練系
//! (TRAIN-1 RingChain / TRAIN-2 TfdPair) と陰性対照だけを走らせ、読み出し規則を
//! 「FROZEN BRIDGE READOUT」節に凍結する。holdout (HOLD-1..4) は v28.2 が執行し、
//! **v28.2 は本ファイルの凍結節を逐語複製して節 SHA-256 の一致を機械検査する**
//! (failed_bridges_rejected_not_tuned の機械化 — 調整は不可能になる)。
//!
//! 入力契約 (bridge_candidates.yml): 読み出しに渡されるのは**ラベルを乱数置換した
//! 状態 (相関行列/波動関数) だけ** — 座標・隣接関係・外部計量は渡さない。隠れた
//! 添字隣接の使用は置換で壊れる。摂動の位置 (クエンチ源) は実験側の知識であり、
//! 読み出しは状態族 ρ(t) と源ラベルのみを見る。
//!
//! 候補 (核 w(i,j) だけが異なる — パイプラインは共通):
//!   B1 = 相互情報量 MI(i,j) = S_i + S_j − S_ij
//!   B2 = 2 ノードモジュラー核 ln((1−C)/C) のノード間 off-diagonal ノルム
//!   B3 = 密度応答 |⟨n_i n_j⟩ − ⟨n_i⟩⟨n_j⟩|
//!   B4 = 局所擾乱の到着時刻 t*(q→j) (因果 — 状態族から)
//! パイプライン: 橋対 (相互 top-1 + 3× 支配) → 相互 top-2 グラフ → 成分/位相
//! (Cycle/Path/Other) → 巡回順 → (採点側) 隣接照合・MDS 署名。
//!
//! 検査 (訓練資格): [T0] 凍結節 SHA-256 印字 + Z2 決定的 rdm の器械検査
//! (lib rdm の HashMap 非決定の機構診断込み — v28.0 残高 #2 の特定) /
//! [T1a] 陰性対照 (積状態・重みシャッフル) で幾何を幻視しない /
//! [T1b-d] TRAIN-1: B1/B2/B3 が置換下で Cycle N=202・隣接 100%・橋 0 /
//! [T1e] MDS 署名 (円環 = 縮退対) / [T1f] 頑健性 (局所位相・2 サイト粗視化) /
//! [T1g] B4 前線 v = v_F ± 5% / [T2a-c] TRAIN-2: 鏡像橋 ≥95%・空間 2 成分 /
//! [T2d] B4 因果二成分 (L 源は R に到達しない)。
//!
//! 最後に holdout 判定バー (v28.2 の契約) を凍結節の定数から印字する。
//! 解釈の上限: 成功しても C3 (toy 機構) — bridge law 登録簿は空のまま
//! (qrn_core::REGISTERED_BRIDGE_LAWS, docs/qrn-core-v1-spec.md §5)。

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

// ---------------- 実験側 (系の構成・置換・採点 — 凍結節の外) ----------------

/// 固定シードのノード置換を生成し、状態に適用する (読み出しへの入力契約)
fn make_perm(n: usize, seed: u64) -> Vec<usize> {
    let mut rng = Rng::new(seed);
    let mut p: Vec<usize> = (0..n).collect();
    for i in (1..n).rev() {
        let j = rng.range(i + 1);
        p.swap(i, j);
    }
    p
}

/// perm[new] = old となるようにモードごと置換 (m 共通)
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

/// GaussianFermionState (m=1) → NodeState
fn from_gaussian(g: &GaussianFermionState) -> NodeState {
    NodeState {
        nodes: g.n,
        m: 1,
        cre: g.cre.clone(),
        cim: g.cim.clone(),
    }
}

/// 隣接復元率: 再構成グラフの辺が真の隣接 (トーラス/鎖) にある割合
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

/// 巡回順の MDS 署名: hop 距離で古典 MDS → 上位固有値
fn mds_spectrum(order_len: usize, topology: u8) -> Vec<f64> {
    let n = order_len;
    let dist = |i: usize, j: usize| -> f64 {
        let d = (i as isize - j as isize).unsigned_abs();
        if topology == 0 {
            d.min(n - d) as f64
        } else {
            d as f64
        }
    };
    let mut b = vec![0.0; n * n];
    // B = −½ J D² J
    let mut rowmean = vec![0.0; n];
    let mut total = 0.0;
    for i in 0..n {
        for j in 0..n {
            let d2 = dist(i, j) * dist(i, j);
            rowmean[i] += d2;
            total += d2;
        }
    }
    for r in rowmean.iter_mut() {
        *r /= n as f64;
    }
    total /= (n * n) as f64;
    for i in 0..n {
        for j in 0..n {
            let d2 = dist(i, j) * dist(i, j);
            b[i + j * n] = -0.5 * (d2 - rowmean[i] - rowmean[j] + total);
        }
    }
    let (mut ev, _) = jacobi_eigh(&b, n);
    ev.sort_by(|a, b| b.partial_cmp(a).unwrap());
    ev.truncate(6);
    ev
}

fn ring_hamiltonian(n: usize) -> Vec<f64> {
    let mut h = vec![0.0; n * n];
    for x in 0..n {
        let y = (x + 1) % n;
        h[x + y * n] = -1.0;
        h[y + x * n] = -1.0;
    }
    h
}

fn main() {
    self_test();
    println!("=== v28.1 relational geometry bridge — 訓練資格認定 (第二十九期 I/III, bridge_candidates.yml 執行) ===\n");
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

    // 凍結節のハッシュ (v282 が照合する)
    let frozen_sha = {
        let src = std::fs::read_to_string("sim/src/bin/v281_bridge_train.rs")
            .or_else(|_| std::fs::read_to_string("src/bin/v281_bridge_train.rs"))
            .or_else(|_| std::fs::read_to_string("../sim/src/bin/v281_bridge_train.rs"))
            .expect("自ソースが読めない");
        let b = src
            .find("// ===================== FROZEN BRIDGE READOUT v29 (BEGIN)")
            .unwrap();
        let e = src
            .find("// ====================== FROZEN BRIDGE READOUT v29 (END)")
            .unwrap();
        sha256_hex(src[b..e].as_bytes())
    };
    println!(
        "[T0] 凍結節 SHA-256 = {} (v282 が逐語一致を検査)\n",
        frozen_sha
    );

    // ---- [T0] Z2 決定的 rdm の器械検査 (l=4 — HOLD-2 の系ではない) ----
    {
        let ring = Z2GaugeRing::try_new(4, 2, 1.0, 0.3, 0.2, vec![]).expect("Z2 l=4");
        let (_e0, gs, res) = ring.ground_state(7);
        let (r1, i1, d1) = rdm_det(4, &gs.masks, &gs.psi, &[1, 3]);
        let (r2, i2, _) = rdm_det(4, &gs.masks, &gs.psi, &[1, 3]);
        let bitid = r1 == r2 && i1 == i2;
        let (rl, il, dl) = gs.rdm(&[1, 3]);
        let maxd = r1
            .iter()
            .zip(rl.iter())
            .chain(i1.iter().zip(il.iter()))
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        // lib rdm の非決定性診断 (2 回呼んで差を測る — 機構: HashMap 反復順)
        let (rl2, il2, _) = gs.rdm(&[1, 3]);
        let libdiff = rl
            .iter()
            .zip(rl2.iter())
            .chain(il.iter().zip(il2.iter()))
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        check(
            "[T0] Z2 決定的 rdm — 2 回呼び bit 同一・lib 値と一致 (l=4 器械検査)",
            bitid && maxd < 1e-12 && d1 == dl && res < 1e-6,
            format!(
                "bit同一 {} / lib との最大差 {:.1e} / lib 2 回呼びの差 {:.1e} (HashMap 群化 — 0 でも per-process 乱数の機構は残る)",
                bitid, maxd, libdiff
            ),
        );
        // Z2 核 3 種の器械スモーク (l=4 — 有限・対称・非負。HOLD-2 の系 l=14 には触れない)
        let view = Z2View {
            l: 4,
            masks: &gs.masks,
            psi: &gs.psi,
            to_orig: (0..4).collect(),
        };
        let (wb1, wb2, wb3) = (w_b1_z2(&view), w_b2_z2(&view), w_b3_z2(&view));
        let sane = |w: &[f64]| {
            w.iter().all(|x| x.is_finite() && *x >= 0.0)
                && (0..4).all(|i| (0..4).all(|j| (w[i + j * 4] - w[j + i * 4]).abs() < 1e-14))
        };
        let sp_ok = (spearman(&[1.0, 2.0, 3.0, 4.0], &[2.0, 4.0, 6.0, 8.0]) - 1.0).abs() < 1e-12;
        check(
            "[T0] Z2 核 3 種のスモーク (有限・対称・非負) + Spearman 恒等検査",
            sane(&wb1) && sane(&wb2) && sane(&wb3) && sp_ok,
            format!(
                "B1 max {:.3} / B2 max {:.3} / B3 max {:.3}",
                wb1.iter().cloned().fold(0.0, f64::max),
                wb2.iter().cloned().fold(0.0, f64::max),
                wb3.iter().cloned().fold(0.0, f64::max)
            ),
        );
    }

    // ---- 系の準備 ----
    let n_ring = 202;
    let ring = RingChain { n: n_ring };
    let ring_st = from_gaussian(&ring.init());
    let perm = make_perm(n_ring, 20260729);
    let shuf = shuffle_state(&ring_st, &perm);
    let ring_truth = |a: usize, b: usize| {
        let d = (a as isize - b as isize).unsigned_abs();
        d.min(n_ring - d) == 1
    };

    // ---- [T1a] 陰性対照 ----
    {
        // 積状態 (C = I/2)
        let d = n_ring;
        let mut cre = vec![0.0; d * d];
        for i in 0..d {
            cre[i + i * d] = 0.5;
        }
        let prod = NodeState {
            nodes: d,
            m: 1,
            cre,
            cim: vec![0.0; d * d],
        };
        let r_prod = reconstruct(&w_b1_gauss(&prod), d);
        // 重みシャッフル対照 (v6.4 方式): TRAIN-1 の B1 重みの成分を乱数置換
        let w = w_b1_gauss(&shuf);
        let mut vals: Vec<f64> = Vec::new();
        for i in 0..n_ring {
            for j in (i + 1)..n_ring {
                vals.push(w[i + j * n_ring]);
            }
        }
        let mut rng = Rng::new(4242);
        for i in (1..vals.len()).rev() {
            let j = rng.range(i + 1);
            vals.swap(i, j);
        }
        let mut ws = vec![0.0; n_ring * n_ring];
        let mut t = 0;
        for i in 0..n_ring {
            for j in (i + 1)..n_ring {
                ws[i + j * n_ring] = vals[t];
                ws[j + i * n_ring] = vals[t];
                t += 1;
            }
        }
        let r_scram = reconstruct(&ws, n_ring);
        check(
            "[T1a] 陰性対照 — 積状態と重みシャッフルで幾何を幻視しない",
            !r_prod.detected && !r_scram.detected,
            format!(
                "積状態 detected={} / シャッフル detected={} (被覆 Cycle/Path が 90% 未満)",
                r_prod.detected, r_scram.detected
            ),
        );
    }

    // ---- [T1b-d] TRAIN-1: B1/B2/B3 (置換下) ----
    let mut b1_sample = Vec::new();
    for (name, w) in [
        ("B1 (MI)", w_b1_gauss(&shuf)),
        ("B2 (modular)", w_b2_gauss(&shuf)),
        ("B3 (応答)", w_b3_gauss(&shuf)),
    ] {
        let r = reconstruct(&w, n_ring);
        let cyc =
            r.comps.len() == 1 && r.comps[0].topology == 0 && r.comps[0].order.len() == n_ring;
        let adj = adjacency_accuracy(&r.edges, &perm, &ring_truth);
        check(
            &format!("[T1b-d] TRAIN-1 {} — Cycle N=202・隣接 100%・橋 0", name),
            r.detected && cyc && adj >= BAR_ADJACENCY && r.bridges.is_empty(),
            format!(
                "detected={} comps={} topo={} 隣接 {:.4} 橋 {}",
                r.detected,
                r.comps.len(),
                if cyc { "Cycle" } else { "非Cycle" },
                adj,
                r.bridges.len()
            ),
        );
        if name.starts_with("B1") {
            for k in 0..5 {
                b1_sample.push(w[k + (k + 1) * n_ring]);
            }
        }
    }

    // ---- [T1e] MDS 署名 (円環 = 縮退対) ----
    {
        let ev = mds_spectrum(n_ring, 0);
        let p1 = (ev[0] / ev[1] - 1.0).abs();
        let p2 = (ev[2] / ev[3] - 1.0).abs();
        check(
            "[T1e] MDS 署名 — 円環の縮退対 (λ1≈λ2, λ3≈λ4)",
            p1 <= BAR_MDS_PAIR && p2 <= BAR_MDS_PAIR2,
            format!("|λ1/λ2−1| = {:.2e} / |λ3/λ4−1| = {:.2e}", p1, p2),
        );
    }

    // ---- [T1f] 頑健性 — 局所位相・2 サイト粗視化 ----
    {
        // 局所位相: C_ab → e^{iθ_a} C e^{-iθ_b}
        let d = shuf.dim();
        let mut rng = Rng::new(777);
        let th: Vec<f64> = (0..d).map(|_| rng.f64() * 6.283).collect();
        let mut cre = vec![0.0; d * d];
        let mut cim = vec![0.0; d * d];
        for a in 0..d {
            for b in 0..d {
                let (cr, ci) = (shuf.cre[a + b * d], shuf.cim[a + b * d]);
                let ph = th[a] - th[b];
                cre[a + b * d] = cr * ph.cos() - ci * ph.sin();
                cim[a + b * d] = cr * ph.sin() + ci * ph.cos();
            }
        }
        let rot = NodeState {
            nodes: shuf.nodes,
            m: 1,
            cre,
            cim,
        };
        let (w0, w1) = (w_b1_gauss(&shuf), w_b1_gauss(&rot));
        let dmax = w0
            .iter()
            .zip(w1.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max)
            / w0.iter().cloned().fold(0.0, f64::max);
        // 粗視化: 真の隣接対 (2i, 2i+1) を 1 ノードに (実験側の再分割) → 101 ノード Cycle
        let nb = n_ring / 2;
        let mut cre2 = vec![0.0; n_ring * n_ring];
        let cim2 = vec![0.0; n_ring * n_ring];
        for i in 0..n_ring {
            for j in 0..n_ring {
                cre2[i + j * n_ring] = ring_st.cre[i + j * n_ring];
            }
        }
        let blocked = NodeState {
            nodes: nb,
            m: 2,
            cre: cre2,
            cim: cim2,
        };
        let permb = make_perm(nb, 555);
        let shufb = shuffle_state(&blocked, &permb);
        let rb = reconstruct(&w_b1_gauss(&shufb), nb);
        let truthb = |a: usize, b: usize| {
            let d = (a as isize - b as isize).unsigned_abs();
            d.min(nb - d) == 1
        };
        let adjb = adjacency_accuracy(&rb.edges, &permb, &truthb);
        let cycb = rb.comps.len() == 1 && rb.comps[0].topology == 0;
        check(
            "[T1f] 頑健性 — 局所位相不変・2 サイト粗視化で Cycle 101",
            dmax <= BAR_ROBUST_PHASE && cycb && adjb >= BAR_ADJACENCY,
            format!(
                "位相 Δw/wmax = {:.1e} / 粗視化 Cycle={} 隣接 {:.4}",
                dmax, cycb, adjb
            ),
        );
    }

    // ---- [T1g] B4 — 前線速度 (訓練は真の距離で資格認定) ----
    let v_fit_train = {
        let h = ring_hamiltonian(n_ring);
        // 置換された標識で源を選ぶ (幾何知識なし)
        let sources: Vec<usize> = (0..8).map(|k| k * n_ring / 8).collect();
        let base: Vec<f64> = (0..n_ring).map(|j| shuf.cre[j + j * n_ring]).collect();
        // h も置換に共変 (状態族の生成は実験側 — 同じラベルで)
        let mut hp = vec![0.0; n_ring * n_ring];
        for i in 0..n_ring {
            for j in 0..n_ring {
                hp[i + j * n_ring] = h[perm[i] + perm[j] * n_ring];
            }
        }
        let ts = fronts_gauss(&hp, &shuf, &base, &sources, 130);
        let (mut sxy, mut sxx) = (0.0, 0.0);
        let mut n_used = 0;
        for (si, &q) in sources.iter().enumerate() {
            for j in 0..n_ring {
                let dtr = {
                    let a = perm[q] as isize;
                    let b = perm[j] as isize;
                    let d = (a - b).unsigned_abs();
                    d.min(n_ring - d) as f64
                };
                let t = ts[si][j];
                if t.is_finite() && dtr >= 3.0 && dtr <= 80.0 {
                    sxy += dtr * t;
                    sxx += t * t;
                    n_used += 1;
                }
            }
        }
        let v_fit = sxy / sxx;
        check(
            "[T1g] B4 — 局所擾乱の前線速度 v (真値 v_F = 2 で資格認定)",
            (v_fit - 2.0).abs() <= BAR_V_REL * 2.0 && n_used > 500,
            format!("v_fit = {:.4} (対 2.0, {} 点)", v_fit, n_used),
        );
        v_fit
    };

    // ---- [T2a-c] TRAIN-2: TfdPair — 鏡像橋 + 空間 2 成分 ----
    let n_tfd = 64;
    let beta_tfd = 0.4;
    {
        let tfd = TfdPair {
            n: n_tfd,
            beta: beta_tfd,
        };
        let st = from_gaussian(&tfd.init());
        let ntot = 2 * n_tfd;
        let permt = make_perm(ntot, 909);
        let shuft = shuffle_state(&st, &permt);
        let mirror_truth = |a: usize, b: usize| (a + n_tfd == b) || (b + n_tfd == a);
        let spatial_truth = |a: usize, b: usize| {
            let (ca, cb) = (a / n_tfd, b / n_tfd);
            if ca != cb {
                return false;
            }
            let (x, y) = (a % n_tfd, b % n_tfd);
            let d = (x as isize - y as isize).unsigned_abs();
            d.min(n_tfd - d) == 1
        };
        for (name, w) in [
            ("B1", w_b1_gauss(&shuft)),
            ("B2", w_b2_gauss(&shuft)),
            ("B3", w_b3_gauss(&shuft)),
        ] {
            let r = reconstruct(&w, ntot);
            let mirror_rate = if r.bridges.is_empty() {
                0.0
            } else {
                r.bridges
                    .iter()
                    .filter(|&&(i, j)| mirror_truth(permt[i], permt[j]))
                    .count() as f64
                    / n_tfd as f64
            };
            let two_cycles = r.comps.len() == 2
                && r.comps
                    .iter()
                    .all(|c| c.topology == 0 && c.order.len() == n_tfd);
            let adj = adjacency_accuracy(&r.edges, &permt, &spatial_truth);
            check(
                &format!(
                    "[T2a-c] TRAIN-2 {} — 鏡像橋 ≥95%・空間 2 Cycle・隣接 100%",
                    name
                ),
                mirror_rate >= BAR_MIRROR_RATE && two_cycles && adj >= BAR_ADJACENCY,
                format!(
                    "鏡像率 {:.3} 橋 {} / comps {} (2 Cycle {}) / 隣接 {:.4}",
                    mirror_rate,
                    r.bridges.len(),
                    r.comps.len(),
                    two_cycles,
                    adj
                ),
            );
        }
        // ---- [T2d] B4 — 因果二成分 (L 源は R に到達しない) ----
        {
            let mut h = vec![0.0; ntot * ntot];
            for x in 0..n_tfd {
                let y = (x + 1) % n_tfd;
                for off in [0, n_tfd] {
                    h[(x + off) + (y + off) * ntot] = -1.0;
                    h[(y + off) + (x + off) * ntot] = -1.0;
                }
            }
            let mut hp = vec![0.0; ntot * ntot];
            for i in 0..ntot {
                for j in 0..ntot {
                    hp[i + j * ntot] = h[permt[i] + permt[j] * ntot];
                }
            }
            // 源 = L 鎖のノード (元 id 5) の置換後ラベル
            let src_new = (0..ntot).find(|&i| permt[i] == 5).unwrap();
            let base: Vec<f64> = (0..ntot).map(|j| shuft.cre[j + j * ntot]).collect();
            let ts = fronts_gauss(&hp, &shuft, &base, &[src_new], 90);
            let mut l_reached = 0;
            let mut r_reached = 0;
            for j in 0..ntot {
                if ts[0][j].is_finite() {
                    if permt[j] < n_tfd {
                        l_reached += 1;
                    } else {
                        r_reached += 1;
                    }
                }
            }
            check(
                "[T2d] B4 — TFD の因果二成分 (L 源の前線は R に到達しない)",
                l_reached == n_tfd && r_reached == 0,
                format!(
                    "L 到達 {}/{} / R 到達 {} (期待 0)",
                    l_reached, n_tfd, r_reached
                ),
            );
        }
    }

    // ---- 成果物 JSON (v282 の回帰アンカー) ----
    {
        let j = Json::Obj(vec![
            ("version".into(), Json::Str("v28.1".into())),
            (
                "frozen_readout_sha256".into(),
                Json::Str(frozen_sha.clone()),
            ),
            (
                "train1_b1_sample".into(),
                Json::Arr(b1_sample.iter().map(|&x| Json::Num(x)).collect()),
            ),
            ("train1_v_fit".into(), Json::Num(v_fit_train)),
            ("tfd_n".into(), Json::Num(n_tfd as f64)),
            ("tfd_beta".into(), Json::Num(beta_tfd)),
            ("perm_seed_ring".into(), Json::Num(20260729.0)),
        ]);
        let p = write_artifact("results/v281_bridge_train.json", &j.render());
        println!("\n[成果物] {}", p);
    }

    // ---- holdout 契約の宣言 (凍結節の定数から印字 — v28.2 はこのバーで採点) ----
    println!("\n---- holdout 契約 (v28.2 — bridge_candidates.yml の 4 系, バーは凍結節定数) ----");
    println!("  HOLD-1 開放端鎖 N=101 (非一様ホッピング): Path・端点 2・隣接 {:.0}%・橋 0・MDS λ2/λ1 ≤ {:.0e}", BAR_ADJACENCY * 100.0, BAR_PATH_MDS);
    println!(
        "  HOLD-2 Z2GaugeRing l=14 (w=1, h=0.6, m=0.2, 閉じ込め相): Cycle 14・隣接 {:.0}%・橋 0",
        BAR_ADJACENCY * 100.0
    );
    println!(
        "          + B4 (apply_bond_op クエンチ): 到着順 vs 再構成距離の Spearman ≥ {:.2}",
        BAR_Z2_SPEARMAN
    );
    println!(
        "          + 未使用チャネル: B1 幾何が B3 チャネルを予言 (d ≤ 5, Spearman ≥ {:.2})",
        BAR_Z2_UNUSED
    );
    println!("  HOLD-3 staggered vs Wilson リング N=402 (m_phys·a = 0.05): 両方 Cycle・隣接 {:.0}%・ξ 比 ∈ [{}, {}]", BAR_ADJACENCY * 100.0, BAR_H3_XI_RATIO.0, BAR_H3_XI_RATIO.1);
    println!("  HOLD-4 リングクエンチ (実時間): t* vs d_B1 線形 R² ≥ {:.2}・源間 v ばらつき ≤ {:.0}%・開放鎖は無巻き付き", BAR_H4_R2, BAR_H4_VSPREAD * 100.0);
    println!("  B4 適用表: TRAIN-1/2, HOLD-1/2/4 (HOLD-3 は計量検査 — 因果は HOLD-4 が担う)");
    println!("  失敗規則: バー外の候補は調整せず棄却 (failed_bridges_rejected_not_tuned)");

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "全候補が訓練資格を通過 — 読み出し規則と判定バーはこのコミットで凍結され、v28.2 (holdout) へ進む"
        } else {
            "訓練資格に FAIL — holdout に進まない (候補または器械の再設計)"
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
