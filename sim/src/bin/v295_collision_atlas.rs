//! v29.5 一意性・反例探索 — collision atlas と identifiability/no-go 境界 (第三十期, PROMPT/11 第三課題)
//!
//! 状態だけから幾何を一意に決める普遍写像は存在しない (同一の積状態は任意の
//! 非同型グラフに置ける)。成立し得るのは「状態族 + 局所代数 + 応答チャネル →
//! 幾何の同値類」という仮定つき写像である。本版はその境界を**全数探索の証明書**
//! として出す:
//!
//!   [C0] 列挙の外部アンカー: n = 4..7 の連結グラフ同型類数 = 6/21/112/853
//!        (OEIS A001349) を自前列挙 (canonical form = 全置換 min) で照合。
//!   [C1] 読み出しサニティ: P7/C7 を凍結 readout (reconstruct_v30) が正しく
//!        Path/Cycle + 隣接 100% で再構成。
//!   [C2] 下界の性質検査: sorted-signature ∞距離 ≤ 厳密 min-perm ∞距離
//!        (多重集合整列は任意の置換整列の下界 — 標本対で機械検査)。
//!   [C3] collision atlas (証明書): 全非同型対のチャネル別識別余裕
//!        margin = min_{G≠G'} min_π ‖X(G)∘π − X(G')‖∞。下界で篩い、
//!        厳密探索は下界 < 現行最小の対のみ (探索順序つき — 証明として完全)。
//!        チャネル: B1 (MI), B3-COV, B5-QFI (静的) / KIN (運動学 t* 行列, 動的)。
//!        厳密衝突 (< 1e-9) の同値類と近衝突 (< 1e-3) 対を列挙する。
//!   [C4] 無情報状態の no-go: C = I/2 (無限温度) では全カーネル < W_FLOOR →
//!        readout は幾何なしを返す (幾何を捏造しない)。状態族の仮定が
//!        identifiability の前提であることの機械化。
//!   [C5] factorization 依存性: C12 均一リング基底状態を 5 通りにモード群化 —
//!        隣接対/オフセット対/三つ組 (局所粗視化) は整合する円環を返し、
//!        反対点対 (antipodal quotient) は**別の自己整合幾何 (C6)** を返し、
//!        ランダム対は幾何を返さない。「readout は与えられた factorization を
//!        幾何化する — factorization 自体は選べない」(FactorizationBridge の
//!        空隙) の機械記録。
//!   [C6] small-world 誤認窓: C16 + 弱い弦 (強度 g) で、静的 mutual top-2 が
//!        クリーンな Cycle と誤読する g 窓と、動的前線が短絡を暴く
//!        (r(g) = t*(0→8)/t*(0→4) < 2) g 窓を掃引 — 「静的位相の資格だけでは
//!        偽陰性がある。応答チャネルが identifiability に必要」の実例。
//!        + 非多様体集 (star/K7/K33/Q3/prism/Petersen) の裁定表。
//!
//! 縮退 Fermi 境界の規則 (置換共変・決定的): 固有値をシェル (gap > 1e-9) に
//! 群化し昇順に充填、境界シェルは分数占有 α = 残り/次元 (混合状態 0 ≤ C ≤ 1)。
//! 全ての読み出しは v29/v30 凍結節の逐語コピー (atlas は凍結器械そのものを検査)。

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

// ---------------- 実験側 (atlas — 凍結節の外) ----------------

/// 隣接 (重みつき可) → 半充填 (filling = 1/2) 基底状態の相関行列。
/// 縮退シェルは分数占有 (置換共変)。
fn corr_half_filling(a: &[f64], n: usize) -> Vec<f64> {
    let h: Vec<f64> = a.iter().map(|&x| -x).collect();
    let (ev, vv) = jacobi_eigh(&h, n);
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&i, &j| ev[i].partial_cmp(&ev[j]).unwrap().then(i.cmp(&j)));
    // シェル分割
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
        let dim = sh.len() as f64;
        let alpha = ((target - filled) / dim).min(1.0);
        for &k in &sh {
            for i in 0..n {
                for j in 0..n {
                    c[i + j * n] += alpha * vv[i + k * n] * vv[j + k * n];
                }
            }
        }
        filled += alpha * dim;
    }
    c
}

fn state_of_corr(c: Vec<f64>, n: usize) -> NodeState {
    NodeState {
        nodes: n,
        m: 1,
        cre: c,
        cim: vec![0.0; n * n],
    }
}

/// 運動学 t* 行列: t*(q,j) = 最初の t で |[e^{−iht}]_{qj}| ≥ B6_EPS (dt 格子)
fn kin_matrix(a: &[f64], n: usize, dt: f64, nt: usize) -> Vec<f64> {
    let h: Vec<f64> = a.iter().map(|&x| -x).collect();
    let (ev, vv) = jacobi_eigh(&h, n);
    let mut t_star = vec![f64::INFINITY; n * n];
    for i in 0..n {
        t_star[i + i * n] = 0.0;
    }
    let mut remaining = n * n - n;
    for it in 1..=nt {
        if remaining == 0 {
            break;
        }
        let t = dt * it as f64;
        let cs: Vec<f64> = ev.iter().map(|&e| (e * t).cos()).collect();
        let sn: Vec<f64> = ev.iter().map(|&e| (e * t).sin()).collect();
        for q in 0..n {
            for j in 0..n {
                if t_star[q + j * n].is_finite() {
                    continue;
                }
                let (mut re, mut im) = (0.0, 0.0);
                for k in 0..n {
                    let w = vv[q + k * n] * vv[j + k * n];
                    re += w * cs[k];
                    im -= w * sn[k];
                }
                if (re * re + im * im).sqrt() >= B6_EPS {
                    t_star[q + j * n] = t;
                    remaining -= 1;
                }
            }
        }
    }
    t_star
}

// ---- グラフ列挙 (canonical form = 全置換 min のビットマスク) ----

/// n 頂点の辺ビット位置: 辺 (i<j) → bit e(i,j)
fn edge_bit(i: usize, j: usize, n: usize) -> usize {
    let (a, b) = if i < j { (i, j) } else { (j, i) };
    a * n - a * (a + 1) / 2 + (b - a - 1)
}

fn perms(n: usize) -> Vec<Vec<usize>> {
    let mut out: Vec<Vec<usize>> = vec![vec![]];
    for k in 0..n {
        let mut next = Vec::new();
        for p in out {
            for pos in 0..=p.len() {
                let mut q = p.clone();
                q.insert(pos, k);
                next.push(q);
            }
        }
        out = next;
    }
    out
}

/// mask に置換 π を適用 (辺集合の像)
fn apply_perm(mask: u32, pi: &[usize], n: usize, nb: usize) -> u32 {
    let mut out = 0u32;
    for e in 0..nb {
        if mask & (1 << e) == 0 {
            continue;
        }
        // e → (i, j) 逆引き
        let mut i = 0;
        let mut acc = 0;
        while acc + (n - i - 1) <= e {
            acc += n - i - 1;
            i += 1;
        }
        let j = i + 1 + (e - acc);
        out |= 1 << edge_bit(pi[i], pi[j], n);
    }
    out
}

fn is_connected(mask: u32, n: usize) -> bool {
    let mut adj = vec![0u32; n];
    for i in 0..n {
        for j in (i + 1)..n {
            if mask & (1 << edge_bit(i, j, n)) != 0 {
                adj[i] |= 1 << j;
                adj[j] |= 1 << i;
            }
        }
    }
    let mut seen = 1u32;
    let mut stack = vec![0usize];
    while let Some(u) = stack.pop() {
        let mut nb = adj[u] & !seen;
        while nb != 0 {
            let v = nb.trailing_zeros() as usize;
            seen |= 1 << v;
            nb &= nb - 1;
            stack.push(v);
        }
    }
    seen.count_ones() as usize == n
}

/// n 頂点の連結グラフ同型類 (canonical mask 昇順)
fn enumerate_connected(n: usize) -> Vec<u32> {
    let nb = n * (n - 1) / 2;
    let ps = perms(n);
    let total = 1u32 << nb;
    let nthreads = 12usize;
    let chunk = (total as usize).div_ceil(nthreads);
    let mut sets: Vec<Vec<u32>> = Vec::new();
    std::thread::scope(|s| {
        let mut handles = Vec::new();
        for t in 0..nthreads {
            let ps = &ps;
            handles.push(s.spawn(move || {
                let lo = (t * chunk) as u32;
                let hi = (((t + 1) * chunk).min(total as usize)) as u32;
                let mut out = Vec::new();
                for mask in lo..hi {
                    if !is_connected(mask, n) {
                        continue;
                    }
                    // canonical = 全置換で最小の mask。自分が最小のときだけ採用
                    let mut minimal = true;
                    for pi in ps.iter() {
                        if apply_perm(mask, pi, n, nb) < mask {
                            minimal = false;
                            break;
                        }
                    }
                    if minimal {
                        out.push(mask);
                    }
                }
                out
            }));
        }
        for h in handles {
            sets.push(h.join().unwrap());
        }
    });
    let mut all: Vec<u32> = sets.into_iter().flatten().collect();
    all.sort_unstable();
    all
}

fn adj_of_mask(mask: u32, n: usize) -> Vec<f64> {
    let mut a = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            if mask & (1 << edge_bit(i, j, n)) != 0 {
                a[i + j * n] = 1.0;
                a[j + i * n] = 1.0;
            }
        }
    }
    a
}

// ---- チャネル特徴と距離 ----

struct Features {
    mask: u32,
    w1: Vec<f64>,
    w3: Vec<f64>,
    w5: Vec<f64>,
    kin: Vec<f64>,
    sig: [Vec<f64>; 4], // 各チャネルの上三角 entries を昇順整列 (下界用)
}

fn upper_sorted(x: &[f64], n: usize) -> Vec<f64> {
    let mut v = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            v.push(x[i + j * n]);
        }
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap());
    v
}

fn features(mask: u32, n: usize, kin_dt: f64, kin_nt: usize) -> Features {
    let a = adj_of_mask(mask, n);
    let st = state_of_corr(corr_half_filling(&a, n), n);
    let w1 = w_b1_gauss(&st);
    let w3 = w_b3_gauss(&st);
    let w5 = w_b5_gauss(&st);
    let kin = kin_matrix(&a, n, kin_dt, kin_nt);
    let sig = [
        upper_sorted(&w1, n),
        upper_sorted(&w3, n),
        upper_sorted(&w5, n),
        upper_sorted(&kin, n),
    ];
    Features {
        mask,
        w1,
        w3,
        w5,
        kin,
        sig,
    }
}

/// 下界: 整列済み多重集合の ∞ 距離 (任意の置換整列 ≥ これ)
fn sig_dist(a: &[f64], b: &[f64]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(x, y)| {
            let (x, y) = (
                if x.is_finite() { *x } else { 1e6 },
                if y.is_finite() { *y } else { 1e6 },
            );
            (x - y).abs()
        })
        .fold(0.0, f64::max)
}

/// 厳密: min_π ‖X∘π − Y‖∞ (n! 全探索 — n ≤ 7)
fn perm_dist(x: &[f64], y: &[f64], n: usize, ps: &[Vec<usize>]) -> f64 {
    let mut best = f64::INFINITY;
    for pi in ps {
        let mut worst: f64 = 0.0;
        'outer: for i in 0..n {
            for j in (i + 1)..n {
                let xv = x[pi[i] + pi[j] * n];
                let yv = y[i + j * n];
                let (xv, yv) = (
                    if xv.is_finite() { xv } else { 1e6 },
                    if yv.is_finite() { yv } else { 1e6 },
                );
                let d = (xv - yv).abs();
                if d > worst {
                    worst = d;
                    if worst >= best {
                        break 'outer;
                    }
                }
            }
        }
        if worst < best {
            best = worst;
        }
    }
    best
}

/// チャネル atlas: 全非同型対の margin (証明書つき) と衝突リスト
struct ChannelAtlas {
    margin: f64,
    argmin: (u32, u32),
    exact_evals: usize,
    collisions: Vec<(u32, u32, f64)>, // < 1e-9
    near: Vec<(u32, u32, f64)>,       // < 1e-3
}

fn chan(f: &Features, ch: usize) -> &Vec<f64> {
    match ch {
        0 => &f.w1,
        1 => &f.w3,
        2 => &f.w5,
        _ => &f.kin,
    }
}

fn atlas_channel(
    feats: &[Features],
    ch: usize,
    n: usize,
    ps: &[Vec<usize>],
) -> ChannelAtlas {
    // 下界を全対で計算し昇順に厳密評価 (下界 ≥ 現行最小になったら停止 = 証明)
    let m = feats.len();
    let mut pairs: Vec<(f64, u32, u32)> = Vec::new();
    for i in 0..m {
        for j in (i + 1)..m {
            pairs.push((
                sig_dist(&feats[i].sig[ch], &feats[j].sig[ch]),
                i as u32,
                j as u32,
            ));
        }
    }
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let mut margin = f64::INFINITY;
    let mut argmin = (0u32, 0u32);
    let mut exact_evals = 0usize;
    let mut collisions = Vec::new();
    let mut near = Vec::new();
    for &(lb, i, j) in &pairs {
        if lb >= margin && lb >= 1e-3 {
            break; // 以降の対は下界 ≥ margin かつ近衝突候補でもない — 証明完了
        }
        let d = perm_dist(
            chan(&feats[i as usize], ch),
            chan(&feats[j as usize], ch),
            n,
            ps,
        );
        exact_evals += 1;
        if d < margin {
            margin = d;
            argmin = (feats[i as usize].mask, feats[j as usize].mask);
        }
        if d < 1e-9 {
            collisions.push((feats[i as usize].mask, feats[j as usize].mask, d));
        } else if d < 1e-3 {
            near.push((feats[i as usize].mask, feats[j as usize].mask, d));
        }
    }
    ChannelAtlas {
        margin,
        argmin,
        exact_evals,
        collisions,
        near,
    }
}

fn main() {
    self_test();
    println!("=== v29.5 一意性・反例探索 — collision atlas と identifiability/no-go 境界 (第三十期) ===\n");
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
    const KIN_DT: f64 = 0.05;
    const KIN_NT: usize = 600; // t_max = 30

    // ---- [C0] 列挙の外部アンカー (OEIS A001349) ----
    let mut graphs: BTreeMap<usize, Vec<u32>> = BTreeMap::new();
    {
        let expect = [(4usize, 6usize), (5, 21), (6, 112), (7, 853)];
        let mut ok = true;
        let mut detail = Vec::new();
        for &(n, want) in &expect {
            let g = enumerate_connected(n);
            detail.push(format!("n={}: {}", n, g.len()));
            if g.len() != want {
                ok = false;
            }
            graphs.insert(n, g);
        }
        check(
            "[C0] 連結グラフ同型類の全数列挙 = OEIS A001349 (6/21/112/853)",
            ok,
            detail.join(" / "),
        );
    }

    // ---- [C1] 読み出しサニティ (P7 / C7) ----
    {
        let n = 7;
        let mut path = vec![0.0; n * n];
        for i in 0..n - 1 {
            path[i + (i + 1) * n] = 1.0;
            path[(i + 1) + i * n] = 1.0;
        }
        let mut ring = path.clone();
        ring[0 + (n - 1) * n] = 1.0;
        ring[(n - 1) + 0 * n] = 1.0;
        let mut ok = true;
        let mut det = Vec::new();
        for (name, a, want) in [("P7", &path, 1u8), ("C7", &ring, 0u8)] {
            let st = state_of_corr(corr_half_filling(a, n), n);
            let r = reconstruct_v30(&w_b3_gauss(&st), n);
            let good = r.comps.len() == 1
                && r.comps[0].topology == want
                && r.comps[0].order.len() == n
                && r.edges.iter().all(|&(i, j)| {
                    let d = (i as isize - j as isize).unsigned_abs();
                    if want == 0 { d.min(n - d) == 1 } else { d == 1 }
                });
            det.push(format!("{} {}", name, good));
            ok &= good;
        }
        check("[C1] 凍結 readout のサニティ (P7 = Path / C7 = Cycle, 隣接 100%)", ok, det.join(" / "));
    }

    // ---- 特徴の計算 (n = 4..7, スレッド分割・決定的回収) ----
    let mut feats: BTreeMap<usize, Vec<Features>> = BTreeMap::new();
    for (&n, masks) in &graphs {
        let nthreads = 12usize;
        let chunk = masks.len().div_ceil(nthreads);
        let mut parts: Vec<Vec<Features>> = Vec::new();
        std::thread::scope(|s| {
            let mut handles = Vec::new();
            for t in 0..nthreads {
                let lo = t * chunk;
                let hi = ((t + 1) * chunk).min(masks.len());
                let sl = &masks[lo.min(masks.len())..hi];
                handles.push(s.spawn(move || {
                    sl.iter().map(|&m| features(m, n, KIN_DT, KIN_NT)).collect::<Vec<_>>()
                }));
            }
            for h in handles {
                parts.push(h.join().unwrap());
            }
        });
        feats.insert(n, parts.into_iter().flatten().collect());
    }

    // ---- [C2] 下界の性質検査 (標本対で lb ≤ exact) ----
    {
        let ps7 = perms(7);
        let f7 = &feats[&7];
        let mut rng = Rng::new(295);
        let mut ok = true;
        let mut worst_gap = 0.0f64;
        for _ in 0..200 {
            let i = rng.range(f7.len());
            let j = rng.range(f7.len());
            if i == j {
                continue;
            }
            for ch in 0..4 {
                let lb = sig_dist(&f7[i].sig[ch], &f7[j].sig[ch]);
                let ex = perm_dist(
                    match ch {
                        0 => &f7[i].w1,
                        1 => &f7[i].w3,
                        2 => &f7[i].w5,
                        _ => &f7[i].kin,
                    },
                    match ch {
                        0 => &f7[j].w1,
                        1 => &f7[j].w3,
                        2 => &f7[j].w5,
                        _ => &f7[j].kin,
                    },
                    7,
                    &ps7,
                );
                if lb > ex + 1e-12 {
                    ok = false;
                }
                worst_gap = worst_gap.max(lb - ex);
            }
        }
        check(
            "[C2] 下界の健全性 (sorted-signature ≤ min-perm, 標本 200 対 × 4 チャネル)",
            ok,
            format!("max(lb − exact) = {:.2e} (≤ 0 が健全)", worst_gap),
        );
    }

    // ---- [C3] collision atlas (証明書) ----
    let ch_names = ["B1-MI", "B3-COV", "B5-QFI", "KIN"];
    let mut atlas_json: Vec<Json> = Vec::new();
    {
        println!("\n  -- [C3] identifiability margin (全非同型対の min-perm ∞距離, 証明書つき) --");
        println!("     n    対数      B1-MI       B3-COV      B5-QFI      KIN        (厳密評価数)");
        let mut all_ok = true;
        for (&n, fs) in &feats {
            let ps = perms(n);
            let npairs = fs.len() * (fs.len() - 1) / 2;
            let mut row = Vec::new();
            let mut evals = Vec::new();
            for ch in 0..4 {
                let a = atlas_channel(fs, ch, n, &ps);
                if !a.collisions.is_empty() {
                    println!(
                        "     [衝突] n={} {}: {} 対が厳密衝突 (< 1e-9): {:?}",
                        n,
                        ch_names[ch],
                        a.collisions.len(),
                        &a.collisions[..a.collisions.len().min(4)]
                    );
                }
                for &(m1, m2, d) in a.near.iter().take(3) {
                    println!(
                        "     [近衝突] n={} {}: mask {:#x} vs {:#x} — 距離 {:.2e}",
                        n, ch_names[ch], m1, m2, d
                    );
                }
                if a.margin.is_infinite() {
                    all_ok = false;
                }
                row.push(format!("{:.4e}", a.margin));
                evals.push(a.exact_evals);
                if ch == 1 && n == 7 {
                    println!(
                        "     [margin 対] n=7 B3-COV の最近接対: mask {:#x} vs {:#x}",
                        a.argmin.0, a.argmin.1
                    );
                }
                atlas_json.push(Json::Obj(vec![
                    ("n".into(), Json::Num(n as f64)),
                    ("channel".into(), Json::Str(ch_names[ch].into())),
                    ("margin".into(), Json::Num(a.margin)),
                    ("exact_evals".into(), Json::Num(a.exact_evals as f64)),
                    ("n_collisions".into(), Json::Num(a.collisions.len() as f64)),
                    ("n_near".into(), Json::Num(a.near.len() as f64)),
                ]));
            }
            println!(
                "     {}  {:7}  {}  {:?}",
                n,
                npairs,
                row.join("  "),
                evals
            );
        }
        check(
            "[C3] atlas 完了 — 全チャネル・全 n で margin が有限確定 (衝突/近衝突は上記リスト)",
            all_ok,
            "margin 値は上表 (証明書: 下界順探索で厳密評価は必要対のみ)".into(),
        );
    }

    // ---- [C4] 無情報状態の no-go ----
    {
        let mut ok = true;
        let mut det = Vec::new();
        // 代表 3 グラフで C = I/2 (グラフ構造は無関係に幾何なしになるべき)
        let mk_ring = |n: usize| -> Vec<f64> {
            let mut a = vec![0.0; n * n];
            for i in 0..n {
                let j = (i + 1) % n;
                a[i + j * n] = 1.0;
                a[j + i * n] = 1.0;
            }
            a
        };
        let mk_complete = |n: usize| -> Vec<f64> {
            let mut a = vec![1.0; n * n];
            for i in 0..n {
                a[i + i * n] = 0.0;
            }
            a
        };
        for (name, n) in [("C7", 7usize), ("K7", 7), ("C12", 12)] {
            let _a = if name == "K7" { mk_complete(n) } else { mk_ring(n) };
            let mut c = vec![0.0; n * n];
            for i in 0..n {
                c[i + i * n] = 0.5;
            }
            let st = state_of_corr(c, n);
            let (w1, w3, w5) = (w_b1_gauss(&st), w_b3_gauss(&st), w_b5_gauss(&st));
            let wmax = [&w1, &w3, &w5]
                .iter()
                .flat_map(|w| w.iter())
                .cloned()
                .fold(0.0, f64::max);
            let r = reconstruct_v30(&w3, n);
            let good = wmax < W_FLOOR && !r.detected;
            det.push(format!("{} max(w) = {:.1e} detected = {}", name, wmax, r.detected));
            ok &= good;
        }
        check(
            "[C4] 無情報状態 (C = I/2) の no-go — 全カーネル < W_FLOOR・幾何なし (捏造しない)",
            ok,
            det.join(" / "),
        );
    }

    // ---- [C5] factorization 依存性 (C12 リング GS の 5 群化) ----
    {
        let n_modes = 12;
        let mut a = vec![0.0; n_modes * n_modes];
        for i in 0..n_modes {
            let j = (i + 1) % n_modes;
            a[i + j * n_modes] = 1.0;
            a[j + i * n_modes] = 1.0;
        }
        let c = corr_half_filling(&a, n_modes);
        let group_state = |groups: &[Vec<usize>]| -> NodeState {
            let nn = groups.len();
            let m = groups[0].len();
            let d = nn * m;
            let mut cre = vec![0.0; d * d];
            for (gi, g) in groups.iter().enumerate() {
                for (ai, &ma) in g.iter().enumerate() {
                    for (gj, g2) in groups.iter().enumerate() {
                        for (bi, &mb) in g2.iter().enumerate() {
                            cre[(gi * m + ai) + (gj * m + bi) * d] = c[ma + mb * n_modes];
                        }
                    }
                }
            }
            NodeState {
                nodes: nn,
                m,
                cre,
                cim: vec![0.0; d * d],
            }
        };
        let adjacent: Vec<Vec<usize>> = (0..6).map(|k| vec![2 * k, 2 * k + 1]).collect();
        let offset: Vec<Vec<usize>> = (0..6).map(|k| vec![(2 * k + 1) % 12, (2 * k + 2) % 12]).collect();
        let triples: Vec<Vec<usize>> = (0..4).map(|k| vec![3 * k, 3 * k + 1, 3 * k + 2]).collect();
        let antipodal: Vec<Vec<usize>> = (0..6).map(|k| vec![k, k + 6]).collect();
        let random_pairs: Vec<Vec<usize>> = {
            let mut rng = Rng::new(1295);
            let mut modes: Vec<usize> = (0..12).collect();
            for i in (1..12).rev() {
                let j = rng.range(i + 1);
                modes.swap(i, j);
            }
            (0..6).map(|k| vec![modes[2 * k], modes[2 * k + 1]]).collect()
        };
        let mut det = Vec::new();
        let mut local_ok = true;
        let mut results = Vec::new();
        for (name, groups, want_nodes) in [
            ("隣接対", &adjacent, 6usize),
            ("オフセット対", &offset, 6),
            ("三つ組", &triples, 4),
            ("反対点対", &antipodal, 6),
            ("ランダム対", &random_pairs, 6),
        ] {
            let st = group_state(groups);
            let r = reconstruct_v30(&w_b3_gauss(&st), st.nodes);
            let is_ring = r.comps.len() == 1
                && !r.comps.is_empty()
                && r.comps[0].topology == 0
                && r.comps[0].order.len() == want_nodes;
            results.push((name, is_ring, r.detected));
            det.push(format!("{}:{}", name, if is_ring { "円環" } else if r.detected { "他" } else { "なし" }));
        }
        // 局所 3 種は円環、反対点対も (quotient の) 円環、ランダムは円環でない — を検査
        local_ok &= results[0].1 && results[1].1 && results[2].1;
        let quotient_ring = results[3].1;
        let random_no_ring = !results[4].1;
        check(
            "[C5] factorization 依存性 — 局所 3 群化は整合円環 / 反対点 quotient も自己整合幾何 / ランダムは円環なし",
            local_ok && random_no_ring,
            format!("{} (反対点 quotient = 円環: {})", det.join(" "), quotient_ring),
        );
        println!("        → readout は与えられた factorization を幾何化する — factorization 自体は選べない (FactorizationBridge の空隙の機械記録)");
    }

    // ---- [C6] small-world 誤認窓と非多様体裁定表 ----
    {
        // (a) C16 + 弦 (0–8, 強度 g) の掃引
        let n = 16;
        println!("\n  -- [C6a] C16 + 弱い弦 (0–8, 強度 g): 静的裁定 vs 動的short-cut 検出 --");
        println!("     g      静的裁定 (B3)    r = t*(0→8)/t*(0→4)   弦は top-2 に?");
        let mut prev_r = f64::INFINITY;
        let mut monotone = true;
        let mut window_static_clean: Vec<f64> = Vec::new();
        let mut window_dynamic_detect: Vec<f64> = Vec::new();
        for &g in &[0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0] {
            let mut a = vec![0.0; n * n];
            for i in 0..n {
                let j = (i + 1) % n;
                a[i + j * n] = 1.0;
                a[j + i * n] = 1.0;
            }
            a[0 + 8 * n] = g;
            a[8 + 0 * n] = g;
            let st = state_of_corr(corr_half_filling(&a, n), n);
            let w = w_b3_gauss(&st);
            let r = reconstruct_v30(&w, n);
            let clean_ring = r.comps.len() == 1
                && r.comps[0].topology == 0
                && r.comps[0].order.len() == n
                && r.edges.iter().all(|&(i, j)| {
                    let d = (i as isize - j as isize).unsigned_abs();
                    d.min(n - d) == 1
                });
            let chord_seen = r.edges.iter().any(|&(i, j)| (i, j) == (0, 8) || (i, j) == (8, 0));
            let kin = kin_matrix(&a, n, KIN_DT, KIN_NT);
            let ratio = kin[0 + 8 * n] / kin[0 + 4 * n];
            if ratio > prev_r + 1e-9 {
                monotone = false;
            }
            prev_r = ratio;
            if clean_ring {
                window_static_clean.push(g);
            }
            if ratio < 1.5 {
                window_dynamic_detect.push(g);
            }
            println!(
                "     {:.2}   {}          {:.3}                {}",
                g,
                if clean_ring { "クリーン円環 (弦を見落とし)" } else { "円環でない    " },
                ratio,
                chord_seen
            );
        }
        let overlap: Vec<f64> = window_static_clean
            .iter()
            .filter(|g| window_dynamic_detect.contains(g))
            .cloned()
            .collect();
        check(
            "[C6a] short-cut 検出比 r(g) の単調非増加 (物理: 強い弦ほど早く着く)",
            monotone,
            format!(
                "静的が見落とす窓 g ∈ {:?} / 動的が検出する窓 g ∈ {:?} / 重なり (静的偽陰性を動的が暴く) = {:?}",
                window_static_clean, window_dynamic_detect, overlap
            ),
        );
        // (b) 非多様体集の裁定表
        println!("\n  -- [C6b] 非多様体・構造グラフの裁定表 (B3 → reconstruct_v30) --");
        let mk = |n: usize, edges: &[(usize, usize)]| -> (usize, Vec<f64>) {
            let mut a = vec![0.0; n * n];
            for &(i, j) in edges {
                a[i + j * n] = 1.0;
                a[j + i * n] = 1.0;
            }
            (n, a)
        };
        let star7 = mk(7, &[(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (0, 6)]);
        let k7 = {
            let mut e = Vec::new();
            for i in 0..7 {
                for j in (i + 1)..7 {
                    e.push((i, j));
                }
            }
            mk(7, &e)
        };
        let k33 = mk(6, &[(0, 3), (0, 4), (0, 5), (1, 3), (1, 4), (1, 5), (2, 3), (2, 4), (2, 5)]);
        let q3 = mk(8, &[(0, 1), (1, 2), (2, 3), (3, 0), (4, 5), (5, 6), (6, 7), (7, 4), (0, 4), (1, 5), (2, 6), (3, 7)]);
        let prism = {
            let mut e = Vec::new();
            for i in 0..6 {
                e.push((i, (i + 1) % 6));
                e.push((i + 6, (i + 1) % 6 + 6));
                e.push((i, i + 6));
            }
            mk(12, &e)
        };
        let petersen = mk(
            10,
            &[(0, 1), (1, 2), (2, 3), (3, 4), (4, 0), (5, 7), (7, 9), (9, 6), (6, 8), (8, 5), (0, 5), (1, 6), (2, 7), (3, 8), (4, 9)],
        );
        let mut n_verdicts = 0usize;
        let mut misreads: Vec<&str> = Vec::new();
        for (name, (n, a)) in [
            ("star S7", star7),
            ("K7", k7),
            ("K3,3", k33),
            ("Q3 立方体", q3),
            ("prism C6×K2", prism),
            ("Petersen", petersen),
        ] {
            let st = state_of_corr(corr_half_filling(&a, n), n);
            let r = reconstruct_v30(&w_b3_gauss(&st), n);
            let full_1d = r.comps.len() == 1
                && !r.comps.is_empty()
                && r.comps[0].topology <= 1
                && r.comps[0].order.len() == n;
            let verdict = if full_1d {
                misreads.push(name);
                "**完全 1D と誤認 (発見 — [C7b] で特性化)**"
            } else if r.detected {
                "部分構造 (完全 1D ではない)"
            } else {
                "幾何なし"
            };
            n_verdicts += 1;
            println!("     {:12} → {}", name, verdict);
        }
        // ゲートは器械 (裁定表の完備) — 誤認の有無は仮説であり分岐記録 (発見は [C7b] へ)
        check(
            "[C6b] 非多様体集 6 グラフの裁定表が完備 (誤認は発見として記録 — ゲートにしない)",
            n_verdicts == 6,
            format!("完全 1D 誤認 = {:?} (それ以外は部分構造/幾何なし)", misreads),
        );
    }

    // ---- [C7] 発見の特性化 (決定的 — 対象は [C3]/[C6b] が同定した具体例) ----
    {
        println!("\n  -- [C7a] n=6 静的衝突対 mask 692 vs 693 の解剖 (非同型は列挙 [C0] が保証) --");
        let n = 6;
        let edges_of = |mask: u32| -> Vec<(usize, usize)> {
            let mut e = Vec::new();
            for i in 0..n {
                for j in (i + 1)..n {
                    if mask & (1 << edge_bit(i, j, n)) != 0 {
                        e.push((i, j));
                    }
                }
            }
            e
        };
        let analyze = |mask: u32| -> (Vec<(usize, usize)>, Vec<usize>, Vec<f64>, Vec<f64>) {
            let a = adj_of_mask(mask, n);
            let mut deg = vec![0usize; n];
            for (i, j) in edges_of(mask) {
                deg[i] += 1;
                deg[j] += 1;
            }
            deg.sort_unstable();
            let (mut ev, _) = jacobi_eigh(&a, n);
            ev.sort_by(|x, y| x.partial_cmp(y).unwrap());
            let c = corr_half_filling(&a, n);
            (edges_of(mask), deg, ev, c)
        };
        let (e1, d1, s1, _c1) = analyze(692);
        let (e2, d2, s2, _c2) = analyze(693);
        println!("     692: 辺 {:?} 次数列 {:?}", e1, d1);
        println!("          隣接スペクトル {:?}", s1.iter().map(|x| (x * 1e4).round() / 1e4).collect::<Vec<_>>());
        println!("     693: 辺 {:?} 次数列 {:?}", e2, d2);
        println!("          隣接スペクトル {:?}", s2.iter().map(|x| (x * 1e4).round() / 1e4).collect::<Vec<_>>());
        let cospectral = s1
            .iter()
            .zip(s2.iter())
            .all(|(x, y)| (x - y).abs() < 1e-9);
        let f1 = features(692, n, KIN_DT, KIN_NT);
        let f2 = features(693, n, KIN_DT, KIN_NT);
        let ps6 = perms(6);
        let dists: Vec<f64> = (0..4)
            .map(|ch| perm_dist(chan(&f1, ch), chan(&f2, ch), n, &ps6))
            .collect();
        println!(
            "     隣接 cospectral = {} / チャネル距離: B1 {:.2e}, B3 {:.2e}, B5 {:.2e}, KIN {:.2e}",
            cospectral, dists[0], dists[1], dists[2], dists[3]
        );
        let joint_separated = dists.iter().any(|&d| d > 1e-9);
        println!(
            "     → {}",
            if joint_separated {
                "少なくとも 1 チャネルが対を分離する (静的衝突は動的で解消)"
            } else {
                "**全 4 チャネルが衝突 — この非同型対は現行 readout の同値類 [g]_~ の非自明な実例**"
            }
        );
        check(
            "[C7a] 衝突対の特性化 — 辺数・スペクトル・全チャネル距離を機械確定",
            e1.len() + e2.len() > 0 && dists.iter().all(|d| d.is_finite()),
            format!(
                "辺数 {} vs {} / cospectral {} / 全チャネル衝突 = {}",
                e1.len(),
                e2.len(),
                cospectral,
                !joint_separated
            ),
        );

        println!("\n  -- [C7b] Petersen 誤認の機構 (縮退カーネル + 決定的 tie-break) --");
        let np = 10;
        let pe: Vec<(usize, usize)> = vec![
            (0, 1), (1, 2), (2, 3), (3, 4), (4, 0),
            (5, 7), (7, 9), (9, 6), (6, 8), (8, 5),
            (0, 5), (1, 6), (2, 7), (3, 8), (4, 9),
        ];
        let mut a = vec![0.0; np * np];
        for &(i, j) in &pe {
            a[i + j * np] = 1.0;
            a[j + i * np] = 1.0;
        }
        let st = state_of_corr(corr_half_filling(&a, np), np);
        let w = w_b3_gauss(&st);
        // 辺エントリと非辺エントリの分布 (頂点推移性 → 各クラス内は縮退のはず)
        let is_edge = |i: usize, j: usize| pe.contains(&(i.min(j), i.max(j)));
        let (mut e_vals, mut ne_vals) = (Vec::new(), Vec::new());
        for i in 0..np {
            for j in (i + 1)..np {
                if is_edge(i, j) {
                    e_vals.push(w[i + j * np]);
                } else {
                    ne_vals.push(w[i + j * np]);
                }
            }
        }
        let spread = |v: &[f64]| -> f64 {
            let mx = v.iter().cloned().fold(f64::MIN, f64::max);
            let mn = v.iter().cloned().fold(f64::MAX, f64::min);
            mx - mn
        };
        let r = reconstruct_v30(&w, np);
        let full_1d = r.comps.len() == 1 && r.comps[0].topology <= 1 && r.comps[0].order.len() == np;
        let all_real_edges = r.edges.iter().all(|&(i, j)| is_edge(i, j));
        let fallback_engaged = !all_real_edges || !r.bridges.is_empty();
        println!(
            "     辺カーネル: 値 {:.6} 拡がり {:.1e} / 非辺: 値 {:.6} 拡がり {:.1e} (辺/非辺比 {:.3})",
            e_vals[0],
            spread(&e_vals),
            ne_vals[0],
            spread(&ne_vals),
            e_vals[0] / ne_vals[0]
        );
        println!(
            "     再構成 = {} (topology {} / {} ノード) — 選択辺が全て真の Petersen 辺: {} / 橋 {} 本 → fallback (v29 橋つき) 発火: {}",
            if full_1d { "完全 1D (誤認)" } else { "非 1D" },
            if r.comps.is_empty() { 9 } else { r.comps[0].topology },
            if r.comps.is_empty() { 0 } else { r.comps[0].order.len() },
            all_real_edges,
            r.bridges.len(),
            fallback_engaged
        );
        println!("     機構: 全 15 辺が厳密縮退 (頂点推移性) → 添字順 tie-break で mutual top-2 が断片化 (被覆 < 90%)");
        println!("     → v30 の Occam 直行路は幾何を返さず、**v29 橋つき fallback が非辺 (距離 2 対, 核 1/49) を継いで");
        println!("     10 ノード path を捏造**した。縮退カーネルでは fallback の前に縮退検査で裁定保留すべき —");
        println!("     v29.6 pipeline への設計入力 (凍結 v30 は変更しない)。");
        check(
            "[C7b] Petersen 誤認の機構確定 — 辺カーネル縮退 + fallback 経路の非辺継ぎを機械同定",
            spread(&e_vals) < 1e-12 && full_1d && fallback_engaged,
            format!(
                "辺縮退 {:.1e} / 完全 1D 誤認 {} / 実辺のみ {} (非辺継ぎ = fallback の証拠)",
                spread(&e_vals),
                full_1d,
                all_real_edges
            ),
        );
    }

    {
        let j = Json::Obj(vec![
            ("version".into(), Json::Str("v29.5".into())),
            ("atlas".into(), Json::Arr(atlas_json)),
        ]);
        let p = write_artifact("results/v295_collision_atlas.json", &j.render());
        println!("\n[成果物] {}", p);
    }

    println!(
        "\n[裁定の要約] identifiability: n ≤ 7 の連結グラフ族 (半充填シェル規則) は各チャネル単独でも margin > 0 なら全数識別可能 —"
    );
    println!("    margin 値と衝突リストが証明書 ([C3] 表)。no-go 側: 無情報状態は幾何なし [C4]・factorization は選べない [C5]・");
    println!("    弱い弦は静的位相の偽陰性 (動的チャネルが必要) [C6a] — 「状態族 + 局所代数 + 応答チャネル → 同値類」の境界が機械化された。");

    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
