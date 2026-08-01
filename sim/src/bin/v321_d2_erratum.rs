//! v32.1 Unit D2 プロトコルの型修復と反例登録 (第三十二期 開始, PROMPT/13 §1)
//!
//! 外部再現 Unit D2 (v31.7 凍結) は「任意の連結グラフ」に対して静的 B3 gap 支持の
//! 厳密再現を要求していた。しかし熱的相関 C = (I+e^{−A})⁻¹ は A の解析関数であり、
//! 連結グラフでは一般に非辺にも歩道由来の非零相関が乗る — 凍結 gap 則 (スケール
//! ガード窓・窓内最大対数段差 ≥ ln 3) が全グラフで直接辺を分離する定理は存在しない。
//! 本器械は、**外部走行 0 件の時点で見つかったプロトコル反例を登録し** (黙って修正
//! しない — superseded_before_external_run)、後継プロトコル D2-S (margin certificate
//! scoped) / D2-R (応答 end-to-end) の資格を機械化する:
//!
//!   [E0] graph6 反例 `F}oXO` の独立デコード: 7 頂点 11 辺 連結 (PROMPT/13 の
//!        NetworkX 走査と別実装・別言語での一致)
//!   [E1] 凍結 gap 則 (v31.6 最終規則 = v320a FROZEN-HOLD7 の逐語) の故障再現:
//!        境界段差 1.9790 ≥ ln 3 なのに非辺尾部段差 2.5130 が窓内最大 →
//!        非辺 (2,5),(3,6) を余剰報告 (11 → 13 辺・欠落 0)。B3 参照値の照合 ≤ 1e-8
//!   [E2] n = 4..7 連結同型類全数 (6/21/112/853 = OEIS A001349) のスキャン:
//!        凍結則の故障は**余剰のみ・欠落 0**。n=7 は 22/853 (PROMPT/13 の独立走査
//!        と一致)・F}oXO は故障集合の元
//!   [E3] B3SupportMarginCertificate ⟺ 凍結則 exact の全数照合 (n=4..7 全 992,
//!        例外 0) — D2-S の scope「事前証明マージン族」が凍結則の成功域と過不足
//!        なく一致し、証明書は真値グラフから走行前に計算できる
//!   [E4] D2-R (曲率則応答 lane, v31.2 恒等式) は全静的故障を修復: 故障グラフ全数で
//!        ŵ = |h_ij|² 支持が欠0余0・ResponseSupportMargin (辺/非辺比) を記録
//!   [E5] ノイズ裁定と SupportNoiseCertificate: σ = 1e-3 → 誤差見積り > バー 0.1 →
//!        Abstain (凍結決定規則 4)。**発見 = 重みバー (0.1) を通る σ = 1e-9 でも、
//!        gap 則の窓ガード (max·1e-3) をノイズ最大値が跨いで余剰辺を作り得る**
//!        (F}oXO: 見積り 4.6e-4 vs ガード 1e-3 は ~1.6σ — 実測で余剰発生)。
//!        D2-R に凍結する支持段の証明書: noise_error_bound·√(2 ln(m/10⁻⁶)) ≤
//!        GAP_GUARD·max ŵ (m = ペア数) — F}oXO σ=1e-9 は不成立 (正しく棄却)・
//!        ring12 σ=1e-9 と F}oXO σ=1e-12 は成立 + 支持一致 (グラフ依存が正しく出る)
//!   [E6] 報告契約の型修復: unit-d-report.schema.json は実 JSON Schema (draft
//!        2020-12)。自作最小 validator の負制御 — pass 報告 適合・failed 報告 適合
//!        (正直な失敗は適合)・必須欠落 不適合・無効な能力昇格 (D1 の能力主張 /
//!        D2-R の語彙外能力) 不適合
//!   [E7] 版分離の整合: protocols/v27.4 = 凍結原本と byte 一致 (sha256)・
//!        protocol-index の supersession 記録・tolerances の参照値 = 本器械の計算値・
//!        replications.yml の 6 条件と external_replications = 0 の不変・旧文面の
//!        逐語保存 (protocols/v31.7/d2-v1-superseded.md)
//!
//! 教訓 (v29.4a/v32.0-B K3-holes と同族): 「プロトコルの主張域は、器械が証明できる
//! 域に事前に切ること」。反例は外部再現の失敗ではなく、実行前の設計入力である。
//!
//! 実行: cargo run --release --bin v321_d2_erratum

use std::fs;
use std::path::Path;
use uft_sim::{jacobi_eigh, matfun_sym, sha256_hex, C64, Rng};

// ---------------------------------------------------------------- 凍結器械 (v320a FROZEN-HOLD7 の逐語コピー)

const GAP_GUARD: f64 = 1e-3; // gap 則スケールガード
const DT_BASE: f64 = 0.02; // 応答 lane dt (スペクトル半径スケール前)
const EPS_PROBE: f64 = 0.3;
const BAR_NOISE_ABSTAIN: f64 = 0.1; // ノイズ誤差見積り > これ → 棄却

fn gibbs_c(h: &[f64], n: usize, beta: f64) -> Vec<f64> {
    matfun_sym(h, n, |x| 1.0 / (1.0 + (beta * x).exp()))
}

/// gap 支持 (v31.6 最終規則 — v320a と逐語一致)
fn support_from_weights(w: &[f64], n: usize) -> Vec<Vec<bool>> {
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
            if i != j && w[i * n + j].abs() > thr {
                adj[i][j] = true;
            }
        }
    }
    adj
}

/// 曲率測定 (v31.2 恒等式の測定 lane — v320a と逐語一致)
fn curvature_w(h: &[f64], n: usize, i: usize, sigma: f64, rng: &mut Rng) -> Vec<f64> {
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

/// ノイズ誤差見積り (真値不使用 — v320a と逐語一致)
fn noise_error_bound(sigma: f64, h_norm1: f64) -> f64 {
    let dt = DT_BASE / h_norm1.max(1.0);
    sigma * 17.0 * 6.0f64.sqrt() / 3.0 / (dt * dt * 4.0 * EPS_PROBE)
}

/// SupportNoiseCertificate (v32.1 で D2-R に凍結する支持段の裁定 — 全量 operational):
/// 重み誤差見積りを 1 読みのノイズスケールとみなし、m 読みの最大値の Gauss 尾
/// (超過確率 δ = 1e-6) が gap 則の窓ガード (GAP_GUARD·max ŵ) を跨がないことを
/// 要求する。凍結 gap 則は**順序対ごと**に閾値判定する (w[i][j] と w[j][i] は独立な
/// 測定) ため m = n(n−1)。跨ぐ場合、重みは読めても**支持**は保証できない — 棄却が正。
fn support_noise_certificate(
    sigma: f64,
    h_norm1: f64,
    max_w: f64,
    n: usize,
) -> (f64, f64, bool) {
    let m = (n * (n - 1)) as f64;
    let z = (2.0 * (m / 1e-6).ln()).sqrt();
    let lhs = noise_error_bound(sigma, h_norm1) * z;
    let rhs = GAP_GUARD * max_w;
    (lhs, rhs, lhs <= rhs)
}

// ---------------------------------------------------------------- グラフ列挙 (v29.5/v31.1 と同一手法)

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

fn apply_perm(mask: u32, pi: &[usize], n: usize, nb: usize) -> u32 {
    let mut out = 0u32;
    for e in 0..nb {
        if mask & (1 << e) == 0 {
            continue;
        }
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

/// n 頂点の連結グラフ同型類 (canonical mask = 全置換で最小)。スレッド分割は
/// 決定的 (結果は分割に依存しない — 最後にソート)。
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
                a[i * n + j] = 1.0;
                a[j * n + i] = 1.0;
            }
        }
    }
    a
}

// ---------------------------------------------------------------- graph6 デコーダ (独立実装)

/// graph6 (McKay 標準形式, n ≤ 62): 先頭バイト = n + 63、以後 6 bit/文字で
/// 上三角を**列優先** (j = 1..n−1, i = 0..j−1) に詰める。
fn g6_decode(s: &str) -> (usize, Vec<f64>) {
    let b = s.as_bytes();
    let n = (b[0] - 63) as usize;
    let mut bits: Vec<bool> = Vec::new();
    for &ch in &b[1..] {
        let v = ch - 63;
        for k in (0..6).rev() {
            bits.push((v >> k) & 1 == 1);
        }
    }
    let mut a = vec![0.0; n * n];
    let mut idx = 0;
    for j in 1..n {
        for i in 0..j {
            if bits[idx] {
                a[i * n + j] = 1.0;
                a[j * n + i] = 1.0;
            }
            idx += 1;
        }
    }
    (n, a)
}

// ---------------------------------------------------------------- 凍結則の適用と margin certificate

/// B3 核: B_ij = |C_ij|² (C = (I+e^{−A})⁻¹, β = 1, h = −A — UNIT_D の D2 凍結設定)
fn b3_kernel(a: &[f64], n: usize) -> Vec<f64> {
    let h: Vec<f64> = a.iter().map(|x| -x).collect();
    let c = gibbs_c(&h, n, 1.0);
    let mut b = vec![0.0; n * n];
    for i in 0..n * n {
        b[i] = c[i] * c[i];
    }
    b
}

struct RuleOutcome {
    missing: usize,
    surplus: usize,
    surplus_pairs: Vec<(usize, usize)>,
    reported_edges: usize,
}

fn run_frozen_rule(b: &[f64], a_true: &[f64], n: usize) -> RuleOutcome {
    let adj = support_from_weights(b, n);
    let mut missing = 0;
    let mut surplus = 0;
    let mut surplus_pairs = Vec::new();
    let mut reported = 0;
    for i in 0..n {
        for j in (i + 1)..n {
            let t = a_true[i * n + j] > 0.5;
            let r = adj[i][j];
            if r {
                reported += 1;
            }
            if t && !r {
                missing += 1;
            }
            if !t && r {
                surplus += 1;
                surplus_pairs.push((i, j));
            }
        }
    }
    RuleOutcome {
        missing,
        surplus,
        surplus_pairs,
        reported_edges: reported,
    }
}

/// B3SupportMarginCertificate (D2-S の scope 定義 — 真値グラフから走行前に計算可能):
///   Case A (窓内に非辺が残る場合):
///     (i) 全真辺が窓内 (min_true ≥ guard) ∧ (ii) 分離 (min_true > max_non) ∧
///     (iii) 境界段差 ln(min_true/max_non) ≥ ln 3 ∧ (iv) 境界段差が窓内段差の
///     厳密最大 — このとき凍結則は境界で切る (unique admissible gap)。
///   Case B (非辺が全て窓外 max_non < guard):
///     (i) 全真辺が窓内 (min_true > guard) ∧ (ii) 窓内 (辺のみ) に段差 ≥ ln 3 が
///     ない — このとき凍結則は窓境界で切る (単一クラスタ)。
/// 定理: certificate 成立 ⟹ 凍結則は支持を欠0余0で返す ([E3] が全数照合)。
struct MarginCert {
    holds: bool,
    case_b: bool,
    min_true_edge: f64,
    max_non_edge: f64,
    boundary_gap: f64,
}

fn margin_certificate(b: &[f64], a_true: &[f64], n: usize) -> MarginCert {
    let mut vals: Vec<(f64, bool)> = Vec::new(); // (B, 真辺か)
    for i in 0..n {
        for j in (i + 1)..n {
            vals.push((b[i * n + j].abs().max(1e-300), a_true[i * n + j] > 0.5));
        }
    }
    vals.sort_by(|x, y| y.0.partial_cmp(&x.0).unwrap());
    let guard = vals[0].0 * GAP_GUARD;
    let min_e = vals
        .iter()
        .filter(|v| v.1)
        .map(|v| v.0)
        .fold(f64::INFINITY, f64::min);
    let max_ne = vals
        .iter()
        .filter(|v| !v.1)
        .map(|v| v.0)
        .fold(0.0f64, f64::max);
    let boundary_gap = (min_e / max_ne.max(1e-300)).ln();
    let win: Vec<(f64, bool)> = vals.iter().filter(|v| v.0 >= guard).cloned().collect();
    let case_b = max_ne < guard;
    let holds = if !case_b {
        let all_edges_in = min_e >= guard;
        let separated = min_e > max_ne;
        let mut strict_max = separated;
        if separated {
            for k in 0..win.len() - 1 {
                let is_boundary = win[k].1 && !win[k + 1].1;
                if is_boundary {
                    continue;
                }
                if (win[k].0 / win[k + 1].0).ln() >= boundary_gap {
                    strict_max = false;
                    break;
                }
            }
        }
        all_edges_in && separated && boundary_gap >= 3.0f64.ln() && strict_max
    } else {
        let all_edges_in = min_e > guard;
        let mut no_sig = true;
        for k in 0..win.len().saturating_sub(1) {
            if (win[k].0 / win[k + 1].0).ln() >= 3.0f64.ln() {
                no_sig = false;
                break;
            }
        }
        all_edges_in && no_sig
    };
    MarginCert {
        holds,
        case_b,
        min_true_edge: min_e,
        max_non_edge: max_ne,
        boundary_gap,
    }
}

/// D2-R 応答 lane: 全ノード probe の曲率重み行列 (v320a の読み出しと同一の詰め方)
fn response_weights(h: &[f64], n: usize, sigma: f64, rng: &mut Rng) -> Vec<f64> {
    let mut wm = vec![0.0; n * n];
    for i in 0..n {
        let wi = curvature_w(h, n, i, sigma, rng);
        for j in 0..n {
            if j != i {
                wm[j * n + i] = wi[j];
            }
        }
    }
    wm
}

// ---------------------------------------------------------------- 最小 JSON パーサ (値木)

#[derive(Clone, Debug, PartialEq)]
enum Jv {
    Null,
    Bool(bool),
    Num(f64),
    Str(String),
    Arr(Vec<Jv>),
    Obj(Vec<(String, Jv)>),
}

struct Jp<'a> {
    b: &'a [u8],
    i: usize,
}

impl<'a> Jp<'a> {
    fn ws(&mut self) {
        while self.i < self.b.len() && (self.b[self.i] as char).is_whitespace() {
            self.i += 1;
        }
    }
    fn peek(&mut self) -> Option<u8> {
        self.ws();
        self.b.get(self.i).copied()
    }
    fn eat(&mut self, c: u8) -> Result<(), String> {
        self.ws();
        if self.b.get(self.i) == Some(&c) {
            self.i += 1;
            Ok(())
        } else {
            Err(format!("位置 {}: '{}' が必要", self.i, c as char))
        }
    }
    fn string(&mut self) -> Result<String, String> {
        self.eat(b'"')?;
        let mut out = String::new();
        while let Some(&c) = self.b.get(self.i) {
            self.i += 1;
            match c {
                b'"' => return Ok(out),
                b'\\' => {
                    let e = *self.b.get(self.i).ok_or("エスケープが切れた")?;
                    self.i += 1;
                    match e {
                        b'"' => out.push('"'),
                        b'\\' => out.push('\\'),
                        b'/' => out.push('/'),
                        b'n' => out.push('\n'),
                        b't' => out.push('\t'),
                        b'r' => out.push('\r'),
                        b'b' | b'f' => out.push(' '),
                        b'u' => {
                            let hex = std::str::from_utf8(&self.b[self.i..self.i + 4])
                                .map_err(|_| "\\u の hex が不正")?;
                            let cp =
                                u32::from_str_radix(hex, 16).map_err(|_| "\\u の hex が不正")?;
                            out.push(char::from_u32(cp).unwrap_or('?'));
                            self.i += 4;
                        }
                        _ => return Err(format!("不明なエスケープ \\{}", e as char)),
                    }
                }
                _ => out.push(c as char),
            }
        }
        Err("文字列が閉じない".into())
    }
    fn value(&mut self) -> Result<Jv, String> {
        match self.peek() {
            Some(b'{') => {
                self.eat(b'{')?;
                let mut m = Vec::new();
                if self.peek() == Some(b'}') {
                    self.eat(b'}')?;
                    return Ok(Jv::Obj(m));
                }
                loop {
                    let k = self.string()?;
                    self.eat(b':')?;
                    let v = self.value()?;
                    m.push((k, v));
                    match self.peek() {
                        Some(b',') => {
                            self.eat(b',')?;
                        }
                        _ => {
                            self.eat(b'}')?;
                            return Ok(Jv::Obj(m));
                        }
                    }
                }
            }
            Some(b'[') => {
                self.eat(b'[')?;
                let mut a = Vec::new();
                if self.peek() == Some(b']') {
                    self.eat(b']')?;
                    return Ok(Jv::Arr(a));
                }
                loop {
                    a.push(self.value()?);
                    match self.peek() {
                        Some(b',') => {
                            self.eat(b',')?;
                        }
                        _ => {
                            self.eat(b']')?;
                            return Ok(Jv::Arr(a));
                        }
                    }
                }
            }
            Some(b'"') => Ok(Jv::Str(self.string()?)),
            Some(_) => {
                self.ws();
                let start = self.i;
                while let Some(&c) = self.b.get(self.i) {
                    if (c as char).is_whitespace() || c == b',' || c == b'}' || c == b']' {
                        break;
                    }
                    self.i += 1;
                }
                let tok = std::str::from_utf8(&self.b[start..self.i]).unwrap_or("");
                match tok {
                    "true" => Ok(Jv::Bool(true)),
                    "false" => Ok(Jv::Bool(false)),
                    "null" => Ok(Jv::Null),
                    _ => tok
                        .parse::<f64>()
                        .map(Jv::Num)
                        .map_err(|_| format!("位置 {}: 不正なトークン '{}'", start, tok)),
                }
            }
            None => Err("空の JSON".into()),
        }
    }
}

fn json_parse(text: &str) -> Result<Jv, String> {
    let mut p = Jp {
        b: text.as_bytes(),
        i: 0,
    };
    let v = p.value()?;
    p.ws();
    if p.i != p.b.len() {
        return Err(format!("位置 {}: 末尾に余分な内容", p.i));
    }
    Ok(v)
}

// ---------------------------------------------------------------- 最小 JSON Schema validator
// 対応語彙 (本スキーマが使う部分集合): type / const / enum / required / properties /
// additionalProperties(false) / items / minItems / maxItems / minLength / minimum /
// maximum / minProperties / pattern (最小 regex) / $ref → #/$defs/*

fn jget<'a>(o: &'a Jv, k: &str) -> Option<&'a Jv> {
    if let Jv::Obj(m) = o {
        m.iter().find(|(kk, _)| kk == k).map(|(_, v)| v)
    } else {
        None
    }
}

/// 最小 regex: ^ と $ で全体固定・原子 = リテラル / \d / [文字クラス (範囲可)]、
/// 量指定子 = {m} / {m,n} (原子直後のみ)。unit-d schema の 2 パターン
/// (^[0-9a-f]{7,40}$ / ^[0-9]{4}-[0-9]{2}-[0-9]{2}$) を被覆する。
#[derive(Clone)]
enum RxAtom {
    Lit(char),
    Digit,
    Class(Vec<(char, char)>),
}

fn rx_compile(p: &str) -> Option<Vec<(RxAtom, usize, usize)>> {
    let b: Vec<char> = p.chars().collect();
    if b.first() != Some(&'^') || b.last() != Some(&'$') {
        return None;
    }
    let mut i = 1;
    let mut out = Vec::new();
    while i < b.len() - 1 {
        let atom = match b[i] {
            '\\' => {
                i += 1;
                if b.get(i) == Some(&'d') {
                    i += 1;
                    RxAtom::Digit
                } else {
                    return None;
                }
            }
            '[' => {
                i += 1;
                let mut cls = Vec::new();
                while i < b.len() - 1 && b[i] != ']' {
                    if i + 2 < b.len() && b[i + 1] == '-' && b[i + 2] != ']' {
                        cls.push((b[i], b[i + 2]));
                        i += 3;
                    } else {
                        cls.push((b[i], b[i]));
                        i += 1;
                    }
                }
                if b.get(i) != Some(&']') {
                    return None;
                }
                i += 1;
                RxAtom::Class(cls)
            }
            c => {
                i += 1;
                RxAtom::Lit(c)
            }
        };
        let (mut lo, mut hi) = (1usize, 1usize);
        if b.get(i) == Some(&'{') {
            let close = (i + 1..b.len()).find(|&k| b[k] == '}')?;
            let inner: String = b[i + 1..close].iter().collect();
            if let Some((a, z)) = inner.split_once(',') {
                lo = a.parse().ok()?;
                hi = z.parse().ok()?;
            } else {
                lo = inner.parse().ok()?;
                hi = lo;
            }
            i = close + 1;
        }
        out.push((atom, lo, hi));
    }
    Some(out)
}

fn rx_atom_match(a: &RxAtom, c: char) -> bool {
    match a {
        RxAtom::Lit(l) => *l == c,
        RxAtom::Digit => c.is_ascii_digit(),
        RxAtom::Class(cls) => cls.iter().any(|(lo, hi)| c >= *lo && c <= *hi),
    }
}

fn rx_match_from(atoms: &[(RxAtom, usize, usize)], s: &[char], si: usize) -> bool {
    if atoms.is_empty() {
        return si == s.len();
    }
    let (a, lo, hi) = &atoms[0];
    // 貪欲 + バックトラック
    let mut k = 0;
    while k < *hi && si + k < s.len() && rx_atom_match(a, s[si + k]) {
        k += 1;
    }
    while k + 1 > *lo {
        if rx_match_from(&atoms[1..], s, si + k) {
            return true;
        }
        if k == 0 {
            break;
        }
        k -= 1;
    }
    if *lo == 0 {
        return rx_match_from(&atoms[1..], s, si);
    }
    false
}

fn rx_match(pattern: &str, s: &str) -> bool {
    match rx_compile(pattern) {
        Some(atoms) => {
            let cs: Vec<char> = s.chars().collect();
            rx_match_from(&atoms, &cs, 0)
        }
        None => true, // 未対応 pattern は素通し (本スキーマでは発生しない)
    }
}

fn jv_eq(a: &Jv, b: &Jv) -> bool {
    match (a, b) {
        (Jv::Num(x), Jv::Num(y)) => x == y,
        (Jv::Str(x), Jv::Str(y)) => x == y,
        (Jv::Bool(x), Jv::Bool(y)) => x == y,
        (Jv::Null, Jv::Null) => true,
        _ => false,
    }
}

fn resolve<'a>(root: &'a Jv, schema: &'a Jv) -> &'a Jv {
    if let Some(Jv::Str(r)) = jget(schema, "$ref") {
        if let Some(name) = r.strip_prefix("#/$defs/") {
            if let Some(defs) = jget(root, "$defs") {
                if let Some(s) = jget(defs, name) {
                    return s;
                }
            }
        }
    }
    schema
}

fn validate(root: &Jv, schema0: &Jv, val: &Jv, path: &str, errs: &mut Vec<String>) {
    let schema = resolve(root, schema0);
    if let Some(c) = jget(schema, "const") {
        if !jv_eq(c, val) {
            errs.push(format!("{}: const 不一致", path));
        }
    }
    if let Some(Jv::Arr(opts)) = jget(schema, "enum") {
        if !opts.iter().any(|o| jv_eq(o, val)) {
            errs.push(format!("{}: enum 外の値", path));
        }
    }
    if let Some(Jv::Str(t)) = jget(schema, "type") {
        let ok = match t.as_str() {
            "object" => matches!(val, Jv::Obj(_)),
            "array" => matches!(val, Jv::Arr(_)),
            "string" => matches!(val, Jv::Str(_)),
            "number" => matches!(val, Jv::Num(_)),
            "integer" => matches!(val, Jv::Num(x) if x.fract() == 0.0),
            "boolean" => matches!(val, Jv::Bool(_)),
            "null" => matches!(val, Jv::Null),
            _ => false,
        };
        if !ok {
            errs.push(format!("{}: type {} でない", path, t));
            return;
        }
    }
    match val {
        Jv::Obj(m) => {
            if let Some(Jv::Arr(req)) = jget(schema, "required") {
                for r in req {
                    if let Jv::Str(k) = r {
                        if !m.iter().any(|(kk, _)| kk == k) {
                            errs.push(format!("{}: 必須フィールド {} が欠落", path, k));
                        }
                    }
                }
            }
            let props = jget(schema, "properties");
            let addl = jget(schema, "additionalProperties");
            for (k, v) in m {
                match props.and_then(|p| jget(p, k)) {
                    Some(s) => validate(root, s, v, &format!("{}/{}", path, k), errs),
                    None => {
                        if let Some(Jv::Bool(false)) = addl {
                            errs.push(format!("{}: 追加プロパティ {} は禁止", path, k));
                        }
                    }
                }
            }
            if let Some(Jv::Num(mp)) = jget(schema, "minProperties") {
                if (m.len() as f64) < *mp {
                    errs.push(format!("{}: minProperties {} 未満", path, mp));
                }
            }
        }
        Jv::Arr(a) => {
            if let Some(items) = jget(schema, "items") {
                for (i, v) in a.iter().enumerate() {
                    validate(root, items, v, &format!("{}/{}", path, i), errs);
                }
            }
            if let Some(Jv::Num(mi)) = jget(schema, "minItems") {
                if (a.len() as f64) < *mi {
                    errs.push(format!("{}: minItems {} 未満", path, mi));
                }
            }
            if let Some(Jv::Num(ma)) = jget(schema, "maxItems") {
                if (a.len() as f64) > *ma {
                    errs.push(format!("{}: maxItems {} 超過", path, ma));
                }
            }
        }
        Jv::Str(s) => {
            if let Some(Jv::Num(ml)) = jget(schema, "minLength") {
                if (s.chars().count() as f64) < *ml {
                    errs.push(format!("{}: minLength {} 未満", path, ml));
                }
            }
            if let Some(Jv::Str(p)) = jget(schema, "pattern") {
                if !rx_match(p, s) {
                    errs.push(format!("{}: pattern {} 不一致", path, p));
                }
            }
        }
        Jv::Num(x) => {
            if let Some(Jv::Num(mn)) = jget(schema, "minimum") {
                if x < mn {
                    errs.push(format!("{}: minimum {} 未満", path, mn));
                }
            }
            if let Some(Jv::Num(mx)) = jget(schema, "maximum") {
                if x > mx {
                    errs.push(format!("{}: maximum {} 超過", path, mx));
                }
            }
        }
        _ => {}
    }
}

// ---------------------------------------------------------------- 報告 fixture (負制御用)

/// 適合すべき pass 報告の例 (全 5 単位)
const REPORT_PASS: &str = r#"{
  "schema_version": "v32.1",
  "replication": {
    "author": "External Reproducer",
    "repository": "https://example.org/qrn-replication",
    "commit": "0123456789abcdef0123456789abcdef01234567",
    "date": "2026-09-01",
    "language": "Julia",
    "protocol_frozen_commit": "aaaabbbbccccddddeeeeffff0000111122223333",
    "shared_numerical_kernel": false
  },
  "units": {
    "D1": { "status": "pass", "c_max_abs_err": 3.1e-9, "support_missing": 0, "support_surplus": 0, "logit_offdiag_max_err": 4.0e-11 },
    "D2S": {
      "status": "pass",
      "graphs": [
        { "graph6": "EhEG", "n": 6, "weighted": false, "certificate_holds": true,
          "min_true_edge_b3": 0.021, "max_non_edge_b3": 0.0002, "boundary_gap_ln": 4.65,
          "verdict": "support_exact", "support_missing": 0, "support_surplus": 0 }
      ],
      "mandatory_negative_control": {
        "graph6": "F}oXO", "frozen_rule_reported_edges": 13,
        "support_missing": 0, "support_surplus": 2
      }
    },
    "D2R": {
      "status": "pass",
      "claimed_capabilities": ["spatial_topology_given_factorization"],
      "cells": [
        { "cell": "kuhn_t3_L4plus", "verdict": "support_exact_topology_match",
          "betti": [1, 3, 3, 1], "links": "closed",
          "support_missing": 0, "support_surplus": 0,
          "response_margin_min_edge": 0.94, "curvature_error_bound": 1.2e-6 },
        { "cell": "cell16_s3", "verdict": "support_exact_topology_match",
          "betti": [1, 0, 0, 1], "links": "closed",
          "support_missing": 0, "support_surplus": 0 },
        { "cell": "solid_kuhn_ball", "verdict": "support_exact_topology_match",
          "betti": [1, 0, 0, 0], "links": "boundary",
          "support_missing": 0, "support_surplus": 0 },
        { "cell": "t3_L3_flag_break_negative", "verdict": "negative_control_confirmed" },
        { "cell": "independent_sparse_weighted", "verdict": "support_exact_topology_match",
          "support_missing": 0, "support_surplus": 0,
          "response_margin_min_edge": 0.31, "curvature_error_bound": 2.0e-6 },
        { "cell": "high_noise_abstain", "verdict": "abstain", "abstain_reason": "insufficient_observation" }
      ]
    },
    "D3": { "status": "pass", "p6_u693_gauge_distance": 1.4e-13,
            "projector_exact_qualification_rejected": true,
            "factorization_nonuniqueness_confirmed": true },
    "D4": { "status": "pass", "uniform_ring_rel_err": 2.1e-5, "eps_pair_max_rel_diff": 3.0e-6,
            "gauge_invariance_ok": true, "tv_transfer_v_independent": true }
  },
  "failures": []
}"#;

/// 適合すべき failed 報告の例 — 正直な失敗はスキーマ適合であること
const REPORT_FAILED_OK: &str = r#"{
  "schema_version": "v32.1",
  "replication": {
    "author": "External Reproducer",
    "repository": "https://example.org/qrn-replication",
    "commit": "0123456789abcdef0123456789abcdef01234567",
    "date": "2026-09-01",
    "language": "C++",
    "protocol_frozen_commit": "aaaabbbbccccddddeeeeffff0000111122223333",
    "shared_numerical_kernel": false
  },
  "units": {
    "D2R": {
      "status": "fail",
      "claimed_capabilities": [],
      "cells": [
        { "cell": "kuhn_t3_L4plus", "verdict": "mismatch",
          "betti": [1, 2, 3, 1], "links": "closed",
          "support_missing": 0, "support_surplus": 1 },
        { "cell": "cell16_s3", "verdict": "support_exact_topology_match",
          "betti": [1, 0, 0, 1], "links": "closed",
          "support_missing": 0, "support_surplus": 0 },
        { "cell": "solid_kuhn_ball", "verdict": "support_exact_topology_match",
          "betti": [1, 0, 0, 0], "links": "boundary",
          "support_missing": 0, "support_surplus": 0 },
        { "cell": "t3_L3_flag_break_negative", "verdict": "negative_control_confirmed" },
        { "cell": "independent_sparse_weighted", "verdict": "support_exact_topology_match",
          "support_missing": 0, "support_surplus": 0 },
        { "cell": "high_noise_abstain", "verdict": "abstain", "abstain_reason": "insufficient_observation" }
      ]
    }
  },
  "failures": [
    { "unit": "D2R", "description": "kuhn_t3_L4plus で支持に余剰 1 辺 — 実装差か本リポジトリ主張の反証かの切り分けを issue で協議" }
  ]
}"#;

/// 不適合であるべき報告 1: 必須フィールド protocol_frozen_commit の欠落
const REPORT_MISSING_REQ: &str = r#"{
  "schema_version": "v32.1",
  "replication": {
    "author": "External Reproducer",
    "repository": "https://example.org/qrn-replication",
    "commit": "0123456789abcdef0123456789abcdef01234567",
    "date": "2026-09-01",
    "language": "Julia",
    "shared_numerical_kernel": false
  },
  "units": {
    "D1": { "status": "pass", "c_max_abs_err": 3.1e-9, "support_missing": 0, "support_surplus": 0, "logit_offdiag_max_err": 4.0e-11 }
  },
  "failures": []
}"#;

/// 不適合であるべき報告 2: D1 (数値再現のみ) が能力を主張 — 無効な能力昇格
const REPORT_BAD_CAP_D1: &str = r#"{
  "schema_version": "v32.1",
  "replication": {
    "author": "External Reproducer",
    "repository": "https://example.org/qrn-replication",
    "commit": "0123456789abcdef0123456789abcdef01234567",
    "date": "2026-09-01",
    "language": "Julia",
    "protocol_frozen_commit": "aaaabbbbccccddddeeeeffff0000111122223333",
    "shared_numerical_kernel": false
  },
  "units": {
    "D1": { "status": "pass", "c_max_abs_err": 3.1e-9, "support_missing": 0, "support_surplus": 0,
            "logit_offdiag_max_err": 4.0e-11,
            "claimed_capabilities": ["spatial_topology_given_factorization"] }
  },
  "failures": []
}"#;

/// 不適合であるべき報告 3: D2-R が語彙外の能力 (clock_calibration) を主張
const REPORT_BAD_CAP_D2R: &str = r#"{
  "schema_version": "v32.1",
  "replication": {
    "author": "External Reproducer",
    "repository": "https://example.org/qrn-replication",
    "commit": "0123456789abcdef0123456789abcdef01234567",
    "date": "2026-09-01",
    "language": "Julia",
    "protocol_frozen_commit": "aaaabbbbccccddddeeeeffff0000111122223333",
    "shared_numerical_kernel": false
  },
  "units": {
    "D2R": {
      "status": "pass",
      "claimed_capabilities": ["clock_calibration"],
      "cells": [
        { "cell": "kuhn_t3_L4plus", "verdict": "support_exact_topology_match" },
        { "cell": "cell16_s3", "verdict": "support_exact_topology_match" },
        { "cell": "solid_kuhn_ball", "verdict": "support_exact_topology_match" },
        { "cell": "t3_L3_flag_break_negative", "verdict": "negative_control_confirmed" },
        { "cell": "independent_sparse_weighted", "verdict": "support_exact_topology_match" },
        { "cell": "high_noise_abstain", "verdict": "abstain", "abstain_reason": "insufficient_observation" }
      ]
    }
  },
  "failures": []
}"#;

// ---------------------------------------------------------------- 補助

/// TOLERANCES 系 (平坦 key: value) から値を拾う (v274 と同型)
fn tol_value(txt: &str, key: &str) -> Option<String> {
    for line in txt.lines() {
        if let Some(rest) = line.strip_prefix(&format!("{}:", key)) {
            let v = rest.trim().trim_matches('"');
            return Some(v.to_string());
        }
    }
    None
}

fn main() {
    uft_sim::self_test();
    println!("=== v32.1 Unit D2 プロトコルの型修復と反例登録 (第三十二期, PROMPT/13 §1) ===\n");
    let root = if Path::new("replications.yml").exists() {
        "."
    } else {
        ".."
    };
    let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p));
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    const G6_COUNTEREXAMPLE: &str = "F}oXO";
    // PROMPT/13 §1.2 の記録値 (NetworkX + numpy の独立走査):
    const REF_MIN_TRUE_EDGE: f64 = 0.01807046;
    const REF_MAX_NON_EDGE: f64 = 0.00249751;
    const REF_NEXT_TIER: f64 = 0.00020236;
    const REF_MAX_B: f64 = 0.03296722;

    // ---- [E0] graph6 反例の独立デコード ----
    let (n7, a_cex) = g6_decode(G6_COUNTEREXAMPLE);
    {
        let mut edges = Vec::new();
        for i in 0..n7 {
            for j in (i + 1)..n7 {
                if a_cex[i * n7 + j] > 0.5 {
                    edges.push((i, j));
                }
            }
        }
        // 連結性: mask に詰め替えて既存の is_connected で検査
        let mut mask = 0u32;
        for &(i, j) in &edges {
            mask |= 1 << edge_bit(i, j, n7);
        }
        let expect_edges = [
            (0, 1),
            (0, 2),
            (0, 3),
            (0, 4),
            (1, 2),
            (1, 3),
            (1, 4),
            (2, 6),
            (3, 5),
            (4, 5),
            (4, 6),
        ];
        let ok = n7 == 7 && edges.len() == 11 && is_connected(mask, 7) && edges == expect_edges;
        check(
            "[E0] graph6 `F}oXO` の独立デコード — 7 頂点 11 辺 連結・辺表一致",
            ok,
            format!("n = {}, |E| = {}, 辺表 = {:?}", n7, edges.len(), edges),
        );
    }

    // ---- [E1] 凍結 gap 則の故障再現 (反例の登録) ----
    let b_cex = b3_kernel(&a_cex, n7);
    let cert_cex = margin_certificate(&b_cex, &a_cex, n7);
    {
        let out = run_frozen_rule(&b_cex, &a_cex, n7);
        // 窓内最大段差 (凍結則が選ぶ切断) を再計算
        let mut vals: Vec<f64> = Vec::new();
        for i in 0..n7 {
            for j in (i + 1)..n7 {
                vals.push(b_cex[i * n7 + j]);
            }
        }
        vals.sort_by(|a, b| b.partial_cmp(a).unwrap());
        let guard = vals[0] * GAP_GUARD;
        let mut best_gap = 0.0f64;
        for k in 0..vals.len() - 1 {
            if vals[k + 1] < guard {
                break;
            }
            best_gap = best_gap.max((vals[k] / vals[k + 1]).ln());
        }
        let refs_ok = (cert_cex.min_true_edge - REF_MIN_TRUE_EDGE).abs() <= 1e-8
            && (cert_cex.max_non_edge - REF_MAX_NON_EDGE).abs() <= 1e-8
            && (vals[0] - REF_MAX_B).abs() <= 1e-8
            && ((cert_cex.max_non_edge / REF_NEXT_TIER).ln() - best_gap).abs() <= 1e-4;
        let ok = out.missing == 0
            && out.surplus == 2
            && out.reported_edges == 13
            && out.surplus_pairs == vec![(2, 5), (3, 6)]
            && (cert_cex.boundary_gap - 1.9790).abs() <= 5e-5
            && (best_gap - 2.5130).abs() <= 5e-5
            && cert_cex.boundary_gap >= 3.0f64.ln()
            && best_gap > cert_cex.boundary_gap
            && refs_ok;
        check(
            "[E1] 凍結 gap 則の故障再現 — 境界段差 1.9790 ≥ ln 3 なのに尾部段差 2.5130 が最大 → 余剰 (2,5),(3,6)",
            ok,
            format!(
                "報告 {} 辺 (真 11)・欠 {} 余 {} {:?}・min_edge B = {:.8}・max_nonedge B = {:.8}・境界段差 {:.4}・窓内最大段差 {:.4}",
                out.reported_edges,
                out.missing,
                out.surplus,
                out.surplus_pairs,
                cert_cex.min_true_edge,
                cert_cex.max_non_edge,
                cert_cex.boundary_gap,
                best_gap
            ),
        );
    }

    // ---- [E2] n = 4..7 連結同型類全数のスキャン ----
    // (故障グラフの (n, mask) を [E4] の修復対象として収集)
    let mut fail_graphs: Vec<(usize, u32)> = Vec::new();
    let mut scan_summary: Vec<(usize, usize, usize, usize)> = Vec::new(); // (n, 類数, 故障, 欠落型故障)
    {
        let expect_counts = [(4usize, 6usize), (5, 21), (6, 112), (7, 853)];
        let mut counts_ok = true;
        let mut missing_total = 0usize;
        let mut cex_found = false;
        // 反例の canonical mask (同型類スキャンで同定するため)
        let nb7 = 21;
        let ps7 = perms(7);
        let mut cex_mask = 0u32;
        for i in 0..7 {
            for j in (i + 1)..7 {
                if a_cex[i * 7 + j] > 0.5 {
                    cex_mask |= 1 << edge_bit(i, j, 7);
                }
            }
        }
        let cex_canon = ps7
            .iter()
            .map(|pi| apply_perm(cex_mask, pi, 7, nb7))
            .min()
            .unwrap();
        for &(n, want) in &expect_counts {
            let classes = enumerate_connected(n);
            counts_ok &= classes.len() == want;
            let mut nfail_n = 0usize;
            let mut nmiss_n = 0usize;
            for &mask in &classes {
                let a = adj_of_mask(mask, n);
                let b = b3_kernel(&a, n);
                let out = run_frozen_rule(&b, &a, n);
                if out.missing > 0 {
                    nmiss_n += 1;
                }
                if out.missing > 0 || out.surplus > 0 {
                    nfail_n += 1;
                    fail_graphs.push((n, mask));
                    if n == 7 && mask == cex_canon {
                        cex_found = true;
                    }
                }
            }
            missing_total += nmiss_n;
            scan_summary.push((n, classes.len(), nfail_n, nmiss_n));
        }
        let n7fail = scan_summary.last().map(|x| x.2).unwrap_or(0);
        let ok = counts_ok && missing_total == 0 && n7fail == 22 && cex_found;
        check(
            "[E2] n=4..7 連結同型類全数 (OEIS A001349) — 凍結則の故障は余剰のみ (欠落 0)・n=7 で 22/853・F}oXO は故障集合の元",
            ok,
            format!(
                "類数/故障 (欠落型): {:?} — PROMPT/13 の独立走査 (NetworkX atlas, 22/853) と一致 = {}",
                scan_summary,
                n7fail == 22
            ),
        );
    }

    // ---- [E3] B3SupportMarginCertificate ⟺ 凍結則 exact の全数照合 ----
    {
        let mut cert_but_not_exact = 0usize;
        let mut exact_but_not_cert = 0usize;
        let mut n_cert = 0usize;
        let mut n_cert_case_b = 0usize;
        let mut n_exact = 0usize;
        let mut n_all = 0usize;
        for &n in &[4usize, 5, 6, 7] {
            let classes = enumerate_connected(n);
            for &mask in &classes {
                let a = adj_of_mask(mask, n);
                let b = b3_kernel(&a, n);
                let out = run_frozen_rule(&b, &a, n);
                let cert = margin_certificate(&b, &a, n);
                let exact = out.missing == 0 && out.surplus == 0;
                n_all += 1;
                if cert.holds {
                    n_cert += 1;
                    if cert.case_b {
                        n_cert_case_b += 1;
                    }
                }
                if exact {
                    n_exact += 1;
                }
                if cert.holds && !exact {
                    cert_but_not_exact += 1;
                }
                if exact && !cert.holds {
                    exact_but_not_cert += 1;
                }
            }
        }
        let ok = cert_but_not_exact == 0 && exact_but_not_cert == 0 && !cert_cex.holds;
        check(
            "[E3] B3SupportMarginCertificate ⟺ 凍結則 exact (n=4..7 全数・例外 0) — F}oXO は certificate 不成立",
            ok,
            format!(
                "全 {} グラフ: cert {} (うち Case B {}) = exact {}・cert∧¬exact = {}・exact∧¬cert = {}・F}}oXO cert = {}",
                n_all,
                n_cert,
                n_cert_case_b,
                n_exact,
                cert_but_not_exact,
                exact_but_not_cert,
                cert_cex.holds
            ),
        );
    }

    // ---- [E4] D2-R 応答 lane が全静的故障を修復 ----
    {
        let mut rng = Rng::new(32101);
        let mut all_exact = true;
        let mut min_margin = f64::INFINITY;
        let mut max_edge_relerr = 0.0f64;
        for &(n, mask) in &fail_graphs {
            let a = adj_of_mask(mask, n);
            let h: Vec<f64> = a.iter().map(|x| -x).collect();
            let wm = response_weights(&h, n, 0.0, &mut rng);
            let adj = support_from_weights(&wm, n);
            let mut min_e = f64::INFINITY;
            let mut max_ne = 0.0f64;
            for i in 0..n {
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let t = a[i * n + j] > 0.5;
                    if t != adj[i][j] {
                        all_exact = false;
                    }
                    if t {
                        min_e = min_e.min(wm[i * n + j].abs());
                        max_edge_relerr = max_edge_relerr.max((wm[i * n + j].abs() - 1.0).abs());
                    } else {
                        max_ne = max_ne.max(wm[i * n + j].abs());
                    }
                }
            }
            min_margin = min_margin.min(min_e / max_ne.max(1e-300));
        }
        // F}oXO 単体も明示的に修復 (fail_graphs は canonical 形なので原表示でも確認)
        let h_cex: Vec<f64> = a_cex.iter().map(|x| -x).collect();
        let wm = response_weights(&h_cex, n7, 0.0, &mut rng);
        let adj = support_from_weights(&wm, n7);
        let mut cex_exact = true;
        for i in 0..n7 {
            for j in 0..n7 {
                if i != j && (a_cex[i * n7 + j] > 0.5) != adj[i][j] {
                    cex_exact = false;
                }
            }
        }
        let ok = all_exact && cex_exact && min_margin >= 1e6 && max_edge_relerr <= 1e-6;
        check(
            "[E4] D2-R 応答 lane (曲率則) は全静的故障を修復 — ŵ = |h_ij|² 支持 欠0余0",
            ok,
            format!(
                "故障 {} グラフ + F}}oXO: 全て欠0余0 = {}・ResponseSupportMargin (辺/非辺比) ≥ {:.1e}・辺重み |ŵ−1| ≤ {:.1e}",
                fail_graphs.len(),
                all_exact && cex_exact,
                min_margin,
                max_edge_relerr
            ),
        );
    }

    // ---- [E5] ノイズ裁定と SupportNoiseCertificate ----
    {
        let mut bad = Vec::new();
        let h_cex: Vec<f64> = a_cex.iter().map(|x| -x).collect();
        let norm1 = |h: &[f64], n: usize| {
            (0..n)
                .map(|r| (0..n).map(|c| h[r * n + c].abs()).sum::<f64>())
                .fold(0.0f64, f64::max)
        };
        let n1_cex = norm1(&h_cex, n7);
        // (a) σ = 1e-3: 凍結決定規則 4 が棄却
        let bound_hi = noise_error_bound(1e-3, n1_cex);
        if bound_hi <= BAR_NOISE_ABSTAIN {
            bad.push(format!("σ=1e-3 が棄却されない (見積り {:.3})", bound_hi));
        }
        // 支持の測定 + certificate 評価 (max ŵ は測定値から — 真値不使用)
        let mut eval = |h: &[f64], n: usize, a_true: &[f64], sigma: f64, seed: u64| {
            let mut rng = Rng::new(seed);
            let wm = response_weights(h, n, sigma, &mut rng);
            let max_w = (0..n * n)
                .filter(|k| k / n != k % n)
                .map(|k| wm[k].abs())
                .fold(0.0f64, f64::max);
            let (lhs, rhs, cert) = support_noise_certificate(sigma, norm1(h, n), max_w, n);
            let adj = support_from_weights(&wm, n);
            // 支持一致は順序対全体で採点 (凍結則の閾値判定は順序対ごと —
            // ノイズ下では片三角のみの余剰が起こり得る)
            let mut missing = 0usize;
            let mut surplus = 0usize;
            for i in 0..n {
                for j in 0..n {
                    if i == j {
                        continue;
                    }
                    let t = a_true[i * n + j] > 0.5;
                    if t && !adj[i][j] {
                        missing += 1;
                    }
                    if !t && adj[i][j] {
                        surplus += 1;
                    }
                }
            }
            (lhs, rhs, cert, missing, surplus)
        };
        // (b) F}oXO σ = 1e-9: 重みバー (0.1) は回答を許すが、certificate は不成立 —
        //     実測でも余剰辺が出る (欠落は出ない)。重みバー単独では支持を守れない。
        let bound_lo = noise_error_bound(1e-9, n1_cex);
        if bound_lo > BAR_NOISE_ABSTAIN {
            bad.push("σ=1e-9 が重みバーで棄却された (回答可のはず)".into());
        }
        let (l_b, r_b, cert_b, miss_b, sur_b) = eval(&h_cex, n7, &a_cex, 1e-9, 32102);
        if cert_b {
            bad.push(format!(
                "F}}oXO σ=1e-9 で certificate が成立してしまった ({:.2e} ≤ {:.2e})",
                l_b, r_b
            ));
        }
        if miss_b != 0 || sur_b == 0 {
            bad.push(format!(
                "F}}oXO σ=1e-9 の実測が想定外 (欠 {} 余 {} — ガード跨ぎの余剰のみが出るはず)",
                miss_b, sur_b
            ));
        }
        // (c) ring12 (単位重み) σ = 1e-9: certificate 成立 + 支持一致 — グラフ依存
        //     (次数 → dt → 増幅) が正しく裁定に出る
        let n12 = 12;
        let mut a_ring = vec![0.0; n12 * n12];
        for k in 0..n12 {
            a_ring[k * n12 + (k + 1) % n12] = 1.0;
            a_ring[((k + 1) % n12) * n12 + k] = 1.0;
        }
        let h_ring: Vec<f64> = a_ring.iter().map(|x| -x).collect();
        let (l_c, r_c, cert_c, miss_c, sur_c) = eval(&h_ring, n12, &a_ring, 1e-9, 32103);
        if !(cert_c && miss_c == 0 && sur_c == 0) {
            bad.push(format!(
                "ring12 σ=1e-9: cert = {} ({:.2e} vs {:.2e})・欠 {} 余 {}",
                cert_c, l_c, r_c, miss_c, sur_c
            ));
        }
        // (d) F}oXO σ = 1e-12: certificate 成立 + 支持一致
        let (l_d, r_d, cert_d, miss_d, sur_d) = eval(&h_cex, n7, &a_cex, 1e-12, 32104);
        if !(cert_d && miss_d == 0 && sur_d == 0) {
            bad.push(format!(
                "F}}oXO σ=1e-12: cert = {} ({:.2e} vs {:.2e})・欠 {} 余 {}",
                cert_d, l_d, r_d, miss_d, sur_d
            ));
        }
        check(
            "[E5] ノイズ裁定 — 規則 4 の棄却 (σ=1e-3)・SupportNoiseCertificate (ガード比) が支持段を守る",
            bad.is_empty(),
            format!(
                "σ=1e-3 見積り {:.1} > {} 棄却 / F}}oXO σ=1e-9: 重みバー通過でも cert 不成立 ({:.2e} > {:.2e}) で実測 余{} 欠{} / ring12 σ=1e-9: cert 成立 支持一致 / F}}oXO σ=1e-12: cert 成立 支持一致{}",
                bound_hi,
                BAR_NOISE_ABSTAIN,
                l_b,
                r_b,
                sur_b,
                miss_b,
                if bad.is_empty() {
                    String::new()
                } else {
                    format!(" — {:?}", bad)
                }
            ),
        );
    }

    // ---- [E6] 報告契約の型修復 — 実 JSON Schema + validator 負制御 ----
    {
        let mut bad = Vec::new();
        match rd("reproducer/protocols/v32.1/unit-d-report.schema.json") {
            Err(_) => bad.push("スキーマファイルが読めない".to_string()),
            Ok(txt) => match json_parse(&txt) {
                Err(e) => bad.push(format!("スキーマが JSON として不正: {}", e)),
                Ok(schema) => {
                    if jget(&schema, "$schema")
                        != Some(&Jv::Str(
                            "https://json-schema.org/draft/2020-12/schema".to_string(),
                        ))
                    {
                        bad.push("$schema が draft 2020-12 でない".to_string());
                    }
                    let cases: [(&str, &str, bool); 5] = [
                        ("pass 報告", REPORT_PASS, true),
                        ("failed 報告 (正直な失敗)", REPORT_FAILED_OK, true),
                        ("必須欠落 (protocol_frozen_commit)", REPORT_MISSING_REQ, false),
                        ("無効な能力昇格 (D1)", REPORT_BAD_CAP_D1, false),
                        ("無効な能力昇格 (D2-R 語彙外)", REPORT_BAD_CAP_D2R, false),
                    ];
                    for (name, txt, want_valid) in cases {
                        match json_parse(txt) {
                            Err(e) => bad.push(format!("{}: fixture が不正 JSON: {}", name, e)),
                            Ok(v) => {
                                let mut errs = Vec::new();
                                validate(&schema, &schema, &v, "#", &mut errs);
                                let valid = errs.is_empty();
                                if valid != want_valid {
                                    bad.push(format!(
                                        "{}: 期待 {} だが {} ({:?})",
                                        name,
                                        if want_valid { "適合" } else { "不適合" },
                                        if valid { "適合" } else { "不適合" },
                                        errs.first()
                                    ));
                                }
                            }
                        }
                    }
                }
            },
        }
        check(
            "[E6] unit-d-report.schema.json (draft 2020-12) — pass/failed 適合・必須欠落/能力昇格 不適合",
            bad.is_empty(),
            if bad.is_empty() {
                "正直な失敗は適合・能力の水増しは不適合 — 契約が型で言えた".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    // ---- [E7] 版分離の整合 ----
    {
        let mut bad = Vec::new();
        // (i) protocols/v27.4 = 凍結原本と byte 一致
        for (copy, orig) in [
            (
                "reproducer/protocols/v27.4/abc-report.schema.json",
                "reproducer/EXPECTED_SCHEMA.json",
            ),
            (
                "reproducer/protocols/v27.4/abc-tolerances.yml",
                "reproducer/TOLERANCES.yml",
            ),
        ] {
            match (rd(copy), rd(orig)) {
                (Ok(c), Ok(o)) => {
                    if sha256_hex(c.as_bytes()) != sha256_hex(o.as_bytes()) {
                        bad.push(format!("{} が原本 {} と byte 不一致", copy, orig));
                    }
                }
                _ => bad.push(format!("{} または {} が読めない", copy, orig)),
            }
        }
        // (ii) protocol-index の supersession 記録
        let idx = rd("reproducer/protocols/v32.1/protocol-index.yml").unwrap_or_default();
        for needle in [
            "superseded_before_external_run",
            "F}oXO",
            "D2-S",
            "D2-R",
            "v321_d2_erratum",
        ] {
            if !idx.contains(needle) {
                bad.push(format!("protocol-index.yml: 「{}」が無い", needle));
            }
        }
        // (iii) 旧文面の逐語保存
        let arch = rd("reproducer/protocols/v31.7/d2-v1-superseded.md").unwrap_or_default();
        for needle in [
            "連結グラフ G (重みつき可) を独立に選ぶ",
            "superseded_before_external_run",
        ] {
            if !arch.contains(needle) {
                bad.push(format!("d2-v1-superseded.md: 「{}」が無い", needle));
            }
        }
        // (iv) tolerances の参照値 = 本器械の計算値
        let tol = rd("reproducer/protocols/v32.1/unit-d-tolerances.yml").unwrap_or_default();
        let tget = |k: &str| tol_value(&tol, k).and_then(|v| v.parse::<f64>().ok());
        let pairs: [(&str, f64, f64); 6] = [
            ("d2s_reference_min_true_edge_b3", cert_cex.min_true_edge, 1e-8),
            ("d2s_reference_max_non_edge_b3", cert_cex.max_non_edge, 1e-8),
            ("d2s_reference_boundary_gap_ln", cert_cex.boundary_gap, 5e-5),
            ("d2s_negative_control_reported_edges", 13.0, 0.0),
            ("d2s_scan_n7_classes", 853.0, 0.0),
            ("d2s_scan_n7_rule_failures", 22.0, 0.0),
        ];
        for (k, want, tolr) in pairs {
            match tget(k) {
                Some(v) if (v - want).abs() <= tolr => {}
                v => bad.push(format!("unit-d-tolerances.yml: {} = {:?} ≠ {:.8}", k, v, want)),
            }
        }
        if tol_value(&tol, "d2s_negative_control_graph6").as_deref() != Some("F}oXO") {
            bad.push("unit-d-tolerances.yml: 負制御 graph6 が F}oXO でない".into());
        }
        // (v) UNIT_D.md が版分離を指す
        let unit_d = rd("reproducer/UNIT_D.md").unwrap_or_default();
        for needle in ["superseded_before_external_run", "D2-S", "D2-R", "protocols/v32.1/"] {
            if !unit_d.contains(needle) {
                bad.push(format!("UNIT_D.md: 「{}」が無い", needle));
            }
        }
        // (vi) replications.yml — 6 条件と external 0 の不変 + 正誤表
        let rep = rd("replications.yml").unwrap_or_default();
        for c in [
            "different_author",
            "independent_repository",
            "no_shared_numerical_kernel",
            "protocol_frozen_before_run",
            "commit_hash_recorded",
            "result_including_failures_public",
        ] {
            if !rep.contains(c) {
                bad.push(format!("replications.yml: 条件 {} が壊れた", c));
            }
        }
        if !rep.contains("external_replications: 0") {
            bad.push("external_replications = 0 が壊れた".into());
        }
        for needle in ["ERR-D2-V1", "superseded_before_external_run"] {
            if !rep.contains(needle) {
                bad.push(format!("replications.yml: 「{}」が無い", needle));
            }
        }
        // (vii) reproducer/CLAIMS.md の D2-R 行
        let claims_md = rd("reproducer/CLAIMS.md").unwrap_or_default();
        for needle in ["D2-R", "QRN-BRIDGE-013"] {
            if !claims_md.contains(needle) {
                bad.push(format!("reproducer/CLAIMS.md: 「{}」が無い", needle));
            }
        }
        check(
            "[E7] 版分離の整合 — v27.4 byte 一致・supersession 記録・tolerances = 計算値・台帳不変・旧文面保存",
            bad.is_empty(),
            if bad.is_empty() {
                "凍結原本は不変のまま、新契約が版として分離された".into()
            } else {
                format!("{:?}", bad)
            },
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "D2-v1 は superseded_before_external_run として登録され、D2-S (margin certificate scoped) / D2-R (応答 end-to-end) が凍結された — 反例は外部再現失敗ではなく実行前の設計入力である"
        } else {
            "**プロトコル修復の破れ** — reproducer/protocols と台帳を修正せよ"
        }
    );
    println!("\n総合判定: {}", if nfail == 0 { "[PASS]" } else { "[FAIL]" });
    if nfail > 0 {
        std::process::exit(1);
    }
}
