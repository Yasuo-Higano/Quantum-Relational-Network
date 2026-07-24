//! v26.4 c₁ の方向依存性 — 異方性行列 c₁[T_aa](q∥b̂) (第二十七期, 経路 B)
//!
//! v24.x はモジュラー規格化 λ が方向依存 (λ_x ≠ λ_⊥) であることを示し、v25.2 で
//! その起源を閉じた。本版は同じ問いを**応答関数側**で立てる: 誘導重力運動項の係数
//! c₁ (v26.3) は格子方向に依存するか — 「真空の弾性」の縦・横係数行列
//!   c₁[T_xx](q∥x̂), c₁[T_yy](q∥x̂), c₁[T_xx](q∥ŷ), c₁[T_yy](q∥ŷ)
//! を N=64 で測る。z 方向は W 対称性 (v26.2 で応答関数レベル 1e-15 を実証) により
//! y と厳密同一 — 実装しない (対称性が保証)。
//!
//! 新エンジン (q ∥ ŷ): T(qŷ) はブロック (ky,kz) → (ky+q,kz) を結ぶ。ky+q が π を
//! 跨ぐと成分 {+0, +π} がスワップする (折返し)。頂点は全て実:
//!   x ボンド (T_xx): x 鎖ホップ ⊗ 成分恒等 (折返しスワップ) — 中点 y 位相は振幅 1
//!   y ボンド (T_yy): (−1)^x ⊗ ysgn·cos(k_y + q/2) 対角 — 中点位相が cos を q/2 シフト
//!   z ボンド: (−1)^x ⊗ σx^{(ky)}·zsgn·cos kz / 質量: m(−1)^x ⊗ σx⊗σx (T_00 用)
//!
//! 検査 (事前登録):
//!  [S0a] 内蔵 dense 照合: q∥ŷ の χ_00, χ_yy (N=8, m∈{0,0.5}, j∈{1,2}) を
//!        実空間 dense (v26.2 の機構の写経) と abs 1e-9 で一致 — 折返し込みの認証
//!  [S0b] x 経路の回帰: REF262 の 20 値 (v26.3 の S0 と同一) + c₁^xx(x̂) が v26.3 の
//!        公表値 0.01826/0.01770/0.01595 (N=64) と ±0.0002 で一致
//!  [S1] χ(0) の 2 エンジン一致: T_yy(0) を「x 経路 j=0」と「y 経路 j=0」で計算し
//!        機械精度で一致 (同一演算子 — エンジン間の恒等検査)
//!  [S2] c₁ 行列の窓安定性 (< 30%) と N=32→64 収束 (< 20%) — v26.3 と同一プロトコル
//!  [S3] 異方性 branch: 縦組 c₁[xx](x̂) vs c₁[yy](ŷ)、横組 c₁[yy](x̂) vs c₁[xx](ŷ) の
//!        差が系統の 3 倍超か (branch A: 応答も方向依存 / branch B: 等方 — λ 物語との対比)
//!  [S4] 変異検出: 折返しスワップを落とす → S0a が検出
//!
//! 事前登録分岐: (a) S0–S2 PASS → c₁ 行列が主結果 (S3 は branch 記録) /
//!   (b) S0a FAIL → ブロック対導出の誤り (dense が真) / (c) S2 FAIL → 窓の再設計。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// v26.2 dense 一次ソース (x 経路回帰用 — v26.3 S0 と同一)
const REF262: [(usize, u32, usize, f64, f64); 10] = [
    (8, 0, 1, 0.001321, 0.112588),
    (8, 0, 2, 0.023548, 0.146865),
    (8, 5, 1, 0.001315, 0.119518),
    (8, 5, 2, 0.020068, 0.148907),
    (12, 0, 1, 0.000334, 0.109067),
    (12, 0, 2, 0.005566, 0.124379),
    (12, 0, 3, 0.023895, 0.150448),
    (12, 5, 1, 0.000295, 0.114817),
    (12, 5, 2, 0.004400, 0.128006),
    (12, 5, 3, 0.020164, 0.149765),
];
/// v26.3 の c₁^xx(x̂) 公表値 (N=64; m = 0, 0.25, 0.5)
const REF263_C1XX: [(f64, f64); 3] = [(0.0, 0.01826), (0.25, 0.01770), (0.5, 0.01595)];

// ---------------- ブロック機構 (v26.3 から写経) ----------------

struct Block {
    h: Vec<f64>,
    dim: usize,
    n: usize,
}

fn block_h(n: usize, m: f64, cky: f64, ckz: f64) -> Block {
    let dim = 4 * n;
    let mut h = vec![0.0f64; dim * dim];
    let id = |x: usize, c: usize| x + n * c;
    let add = |h: &mut Vec<f64>, a: usize, b: usize, t: f64| {
        h[b + a * dim] += t;
        h[a + b * dim] += t;
    };
    let ysgn = [1.0, -1.0, 1.0, -1.0];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            add(&mut h, id(x, c), id((x + 1) % n, c), 0.5 * tw);
            h[id(x, c) + id(x, c) * dim] += px * ysgn[c] * cky;
            if c == 0 || c == 2 {
                add(&mut h, id(x, c), id(x, c + 1), px * zsgn[c] * ckz);
            }
        }
        add(&mut h, id(x, 0), id(x, 3), px * m);
        add(&mut h, id(x, 1), id(x, 2), px * m);
    }
    Block { h, dim, n }
}

/// x 経路の頂点 (q∥x̂, ブロック対角)。which: 0 = T_00, 1 = T_xx, 2 = T_yy
fn vertex_qx(
    blk: &Block,
    m: f64,
    cky: f64,
    ckz: f64,
    q: f64,
    which: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (n, dim) = (blk.n, blk.dim);
    let id = |x: usize, c: usize| x + n * c;
    let mut re = vec![0.0f64; dim * dim];
    let mut im = vec![0.0f64; dim * dim];
    let addc = |re: &mut Vec<f64>, im: &mut Vec<f64>, a: usize, b: usize, t: f64, ph: f64| {
        let (cp, sp) = (ph.cos(), ph.sin());
        re[b + a * dim] += t * cp;
        re[a + b * dim] += t * cp;
        im[b + a * dim] += t * sp;
        im[a + b * dim] += t * sp;
    };
    let ysgn = [1.0, -1.0, 1.0, -1.0];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        let ph_site = q * x as f64;
        let ph_bond = q * (x as f64 + 0.5);
        for c in 0..4 {
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            if which == 0 || which == 1 {
                addc(
                    &mut re,
                    &mut im,
                    id(x, c),
                    id((x + 1) % n, c),
                    0.5 * tw,
                    ph_bond,
                );
            }
            if which == 0 || which == 2 {
                let (cp, sp) = (ph_site.cos(), ph_site.sin());
                re[id(x, c) + id(x, c) * dim] += px * ysgn[c] * cky * cp;
                im[id(x, c) + id(x, c) * dim] += px * ysgn[c] * cky * sp;
            }
            if which == 0 && (c == 0 || c == 2) {
                addc(
                    &mut re,
                    &mut im,
                    id(x, c),
                    id(x, c + 1),
                    px * zsgn[c] * ckz,
                    ph_site,
                );
            }
        }
        if which == 0 {
            addc(&mut re, &mut im, id(x, 0), id(x, 3), px * m, ph_site);
            addc(&mut re, &mut im, id(x, 1), id(x, 2), px * m, ph_site);
        }
    }
    (re, im)
}

/// ブロック内 Lehmann (v26.3 写経)
fn chi_block(w: &[f64], v: &[f64], dim: usize, ore: &[f64], oim: &[f64]) -> f64 {
    let nocc = dim / 2;
    let mut tv_re = vec![0.0f64; dim * nocc];
    let mut tv_im = vec![0.0f64; dim * nocc];
    for ccol in 0..nocc {
        for r in 0..dim {
            let (mut sr, mut si) = (0.0, 0.0);
            for k in 0..dim {
                let vv = v[k + ccol * dim];
                sr += ore[k + r * dim] * vv;
                si += oim[k + r * dim] * vv;
            }
            tv_re[r + ccol * dim] = sr;
            tv_im[r + ccol * dim] = si;
        }
    }
    let mut chi = 0.0f64;
    for mu in nocc..dim {
        for nu in 0..nocc {
            let (mut mr, mut mi) = (0.0f64, 0.0f64);
            for k in 0..dim {
                let vm = v[k + mu * dim];
                mr += vm * tv_re[k + nu * dim];
                mi += vm * tv_im[k + nu * dim];
            }
            chi += 2.0 * (mr * mr + mi * mi) / (w[mu] - w[nu]);
        }
    }
    chi
}

/// x 経路の走査 (v26.3 と同じ流れ; which ∈ {0,1,2})
fn chi_scan_x(n: usize, m: f64, js: &[usize], which: usize, nthreads: usize) -> Vec<f64> {
    let nb = n / 2;
    let labels: Vec<(usize, usize)> = (0..nb * nb).map(|i| (i % nb, i / nb)).collect();
    let mut parts: Vec<Option<Vec<f64>>> = Vec::new();
    parts.resize_with(labels.len(), || None);
    let chunk = labels.len().div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in parts.chunks_mut(chunk).enumerate() {
            let labels = &labels;
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let (jy, jz) = labels[t * chunk + i];
                    let cky = (2.0 * PI * jy as f64 / n as f64).cos();
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    let blk = block_h(n, m, cky, ckz);
                    let (w, v) = jacobi_eigh(&blk.h, blk.dim);
                    let mut out = Vec::with_capacity(js.len());
                    for &j in js {
                        let q = 2.0 * PI * j as f64 / n as f64;
                        let (re, im) = vertex_qx(&blk, m, cky, ckz, q, which);
                        out.push(chi_block(&w, &v, blk.dim, &re, &im));
                    }
                    *slot = Some(out);
                }
            });
        }
    });
    let ns3 = (n * n * n) as f64;
    let mut acc = vec![0.0f64; js.len()];
    for p in parts {
        for (ji, &x) in p.unwrap().iter().enumerate() {
            acc[ji] += x;
        }
    }
    acc.iter().map(|&x| x / ns3).collect()
}

// ---------------- y 経路 (ブロック対) ----------------

/// q∥ŷ の頂点: ブロック (jy) → (jy+j) [成分スワップ sw 込み]。which: 0=T_00, 1=T_xx, 2=T_yy。
/// 全て実。行 = 行き先ブロックの基底 / 列 = 元ブロックの基底。
fn vertex_qy(n: usize, m: f64, ky: f64, ckz: f64, q: f64, sw: bool, which: usize) -> Vec<f64> {
    let dim = 4 * n;
    let id = |x: usize, c: usize| x + n * c;
    // 折返しスワップ: 成分 0↔1, 2↔3 (ky+π 側へ入る)
    let tc = |c: usize| -> usize {
        if sw {
            [1usize, 0, 3, 2][c]
        } else {
            c
        }
    };
    let mut v = vec![0.0f64; dim * dim];
    let zsgn = [1.0, 1.0, -1.0, -1.0];
    for x in 0..n {
        let px = if x % 2 == 0 { 1.0 } else { -1.0 };
        for c in 0..4 {
            let kyc = ky + if c == 1 || c == 3 { PI } else { 0.0 };
            let tw = if x == n - 1 { -1.0 } else { 1.0 };
            if which == 0 || which == 1 {
                // x ボンド: 振幅 1 の ky→ky+q シフト (実)。行き先成分 tc(c)。
                v[id((x + 1) % n, tc(c)) + id(x, c) * dim] += 0.5 * tw;
                v[id(x, tc(c)) + id((x + 1) % n, c) * dim] += 0.5 * tw;
            }
            if which == 0 || which == 2 {
                // y ボンド: (−1)^x cos(k_y^{(c)} + q/2) — 成分符号は kyc の +π に吸収済み
                let coef = px * (kyc + q / 2.0).cos();
                v[id(x, tc(c)) + id(x, c) * dim] += coef;
            }
            if (which == 0 || which == 3) && (c == 0 || c == 2) {
                // z ボンド: σx^{(ky)} (c ↔ c+1) × zsgn·cos kz × (−1)^x
                let coef = px * zsgn[c] * ckz;
                v[id(x, tc(c + 1)) + id(x, c) * dim] += coef;
                v[id(x, tc(c)) + id(x, c + 1) * dim] += coef;
            }
        }
        if which == 0 {
            // 質量: σx⊗σx (0↔3, 1↔2) × m(−1)^x
            v[id(x, tc(3)) + id(x, 0) * dim] += px * m;
            v[id(x, tc(0)) + id(x, 3) * dim] += px * m;
            v[id(x, tc(2)) + id(x, 1) * dim] += px * m;
            v[id(x, tc(1)) + id(x, 2) * dim] += px * m;
        }
    }
    v
}

/// ブロック対 Lehmann: 元 (w1,v1) の占有 → 先 (w2,v2) の非占有。
/// 頂点の格納規約は v[標的 + 元·dim] — 縮約は「元」側 (開発記録: run1 は転置縮約
/// [o[k + r*dim]] のバグ — 各片は転置で ±自身に戻るため単片 χ は不変で照合を通過し、
/// 片間の相対符号だけが壊れて T_00 の保存キャンセルが破れた。numpy 三者照合で特定)。
fn chi_pair(w1: &[f64], v1: &[f64], w2: &[f64], v2: &[f64], dim: usize, o: &[f64]) -> f64 {
    let nocc = dim / 2;
    let mut tv = vec![0.0f64; dim * nocc];
    for ccol in 0..nocc {
        for r in 0..dim {
            let mut s = 0.0;
            for k in 0..dim {
                s += o[r + k * dim] * v1[k + ccol * dim];
            }
            tv[r + ccol * dim] = s;
        }
    }
    let mut chi = 0.0f64;
    for mu in nocc..dim {
        for nu in 0..nocc {
            let mut mr = 0.0f64;
            for k in 0..dim {
                mr += v2[k + mu * dim] * tv[k + nu * dim];
            }
            chi += 2.0 * mr * mr / (w2[mu] - w1[nu]);
        }
    }
    chi
}

/// y 経路の走査: kz 行ごとに全ブロックを対角化し、全 (j, which) を評価。
/// mutate_fold = true で折返しスワップを落とす (S4 用)。
fn chi_scan_y(
    n: usize,
    m: f64,
    js: &[usize],
    whiches: &[usize],
    nthreads: usize,
    mutate_fold: bool,
) -> Vec<Vec<f64>> {
    // 戻り値 [which idx][j idx]
    let nb = n / 2;
    let mut rows: Vec<Option<Vec<Vec<f64>>>> = Vec::new();
    rows.resize_with(nb, || None);
    let chunk = nb.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let jz = t * chunk + i;
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    // 行内の全ブロックを対角化
                    let mut eigs: Vec<(Vec<f64>, Vec<f64>)> = Vec::with_capacity(nb);
                    for jy in 0..nb {
                        let cky = (2.0 * PI * jy as f64 / n as f64).cos();
                        let blk = block_h(n, m, cky, ckz);
                        eigs.push(jacobi_eigh(&blk.h, blk.dim));
                    }
                    let dim = 4 * n;
                    let mut acc = vec![vec![0.0f64; js.len()]; whiches.len()];
                    for (wi, &which) in whiches.iter().enumerate() {
                        for (ji, &j) in js.iter().enumerate() {
                            let q = 2.0 * PI * j as f64 / n as f64;
                            for jy in 0..nb {
                                let ky = 2.0 * PI * jy as f64 / n as f64;
                                let mut jt = jy + j;
                                let mut sw = false;
                                while jt >= nb {
                                    jt -= nb;
                                    sw = !sw;
                                }
                                let sw_eff = if mutate_fold { false } else { sw };
                                let o = vertex_qy(n, m, ky, ckz, q, sw_eff, which);
                                let (w1, v1) = &eigs[jy];
                                let (w2, v2) = &eigs[jt];
                                acc[wi][ji] += chi_pair(w1, v1, w2, v2, dim, &o);
                            }
                        }
                    }
                    *slot = Some(acc);
                }
            });
        }
    });
    let ns3 = (n * n * n) as f64;
    let mut out = vec![vec![0.0f64; js.len()]; whiches.len()];
    for r in rows {
        let r = r.unwrap();
        for wi in 0..whiches.len() {
            for ji in 0..js.len() {
                out[wi][ji] += r[wi][ji] / ns3;
            }
        }
    }
    out
}

// ---------------- 内蔵 dense 照合 (v26.2 の機構の写経, q∥ŷ 用) ----------------

fn build_h_dense(n: usize, m_stag: f64) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let add = |h: &mut Vec<f64>, i: usize, j: usize, t: f64| {
        h[j + i * ns] += t;
        h[i + j * ns] += t;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                let tw = if x == n - 1 { -1.0 } else { 1.0 };
                add(&mut h, i, idx((x + 1) % n, y, z), 0.5 * tw);
                let ey = if x % 2 == 0 { 1.0 } else { -1.0 };
                add(&mut h, i, idx(x, (y + 1) % n, z), 0.5 * ey);
                let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                add(&mut h, i, idx(x, y, (z + 1) % n), 0.5 * ez);
                let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                h[i + i * ns] += m_stag * sgn;
            }
        }
    }
    h
}

/// dense の T_00(q) / T_yy(q) (q 任意)。which: 0 = T_00, 2 = T_yy
fn dense_vertex(n: usize, m_stag: f64, q: [f64; 3], which: usize) -> (Vec<f64>, Vec<f64>) {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut re = vec![0.0f64; ns * ns];
    let mut im = vec![0.0f64; ns * ns];
    let addc = |re: &mut Vec<f64>, im: &mut Vec<f64>, i: usize, j: usize, t: f64, ph: f64| {
        let (cp, sp) = (ph.cos(), ph.sin());
        re[j + i * ns] += t * cp;
        re[i + j * ns] += t * cp;
        im[j + i * ns] += t * sp;
        im[i + j * ns] += t * sp;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                let (xf, yf, zf) = (x as f64, y as f64, z as f64);
                if which == 0 {
                    let tw = if x == n - 1 { -1.0 } else { 1.0 };
                    let ph = q[0] * (xf + 0.5) + q[1] * yf + q[2] * zf;
                    addc(&mut re, &mut im, i, idx((x + 1) % n, y, z), 0.5 * tw, ph);
                }
                let ey = if x % 2 == 0 { 1.0 } else { -1.0 };
                let ph = q[0] * xf + q[1] * (yf + 0.5) + q[2] * zf;
                addc(&mut re, &mut im, i, idx(x, (y + 1) % n, z), 0.5 * ey, ph);
                if which == 0 {
                    let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                    let ph = q[0] * xf + q[1] * yf + q[2] * (zf + 0.5);
                    addc(&mut re, &mut im, i, idx(x, y, (z + 1) % n), 0.5 * ez, ph);
                    let sgn = if (x + y + z) % 2 == 0 { 1.0 } else { -1.0 };
                    let ph = q[0] * xf + q[1] * yf + q[2] * zf;
                    re[i + i * ns] += m_stag * sgn * ph.cos();
                    im[i + i * ns] += m_stag * sgn * ph.sin();
                }
            }
        }
    }
    (re, im)
}

fn chi_dense(w: &[f64], v: &[f64], ns: usize, ore: &[f64], oim: &[f64]) -> f64 {
    let nocc = ns / 2;
    let mut tv_re = vec![0.0f64; ns * nocc];
    let mut tv_im = vec![0.0f64; ns * nocc];
    for c in 0..nocc {
        for r in 0..ns {
            let (mut sr, mut si) = (0.0, 0.0);
            for k in 0..ns {
                let vv = v[k + c * ns];
                sr += ore[k + r * ns] * vv;
                si += oim[k + r * ns] * vv;
            }
            tv_re[r + c * ns] = sr;
            tv_im[r + c * ns] = si;
        }
    }
    let mut chi = 0.0f64;
    for mu in nocc..ns {
        for nu in 0..nocc {
            let (mut mr, mut mi) = (0.0f64, 0.0f64);
            for k in 0..ns {
                let vm = v[k + mu * ns];
                mr += vm * tv_re[k + nu * ns];
                mi += vm * tv_im[k + nu * ns];
            }
            chi += 2.0 * (mr * mr + mi * mi) / (w[mu] - w[nu]);
        }
    }
    chi
}

/// 最小二乗 (v26.3 写経)
fn lstsq(xs: &[Vec<f64>], y: &[f64]) -> Vec<f64> {
    let p = xs.len();
    let mut a = vec![0.0f64; p * p];
    let mut b = vec![0.0f64; p];
    for i in 0..p {
        for jj in 0..p {
            a[jj + i * p] = xs[i].iter().zip(&xs[jj]).map(|(u, v)| u * v).sum();
        }
        b[i] = xs[i].iter().zip(y).map(|(u, v)| u * v).sum();
    }
    let mut idx: Vec<usize> = (0..p).collect();
    for col in 0..p {
        let piv = (col..p)
            .max_by(|&r1, &r2| {
                a[col + idx[r1] * p]
                    .abs()
                    .partial_cmp(&a[col + idx[r2] * p].abs())
                    .unwrap()
            })
            .unwrap();
        idx.swap(col, piv);
        let d = a[col + idx[col] * p];
        for r in col + 1..p {
            let f = a[col + idx[r] * p] / d;
            for cc in col..p {
                a[cc + idx[r] * p] -= f * a[cc + idx[col] * p];
            }
            b[idx[r]] -= f * b[idx[col]];
        }
    }
    let mut out = vec![0.0f64; p];
    for col in (0..p).rev() {
        let mut s = b[idx[col]];
        for cc in col + 1..p {
            s -= a[cc + idx[col] * p] * out[cc];
        }
        out[col] = s / a[col + idx[col] * p];
    }
    out
}

fn main() {
    self_test();
    println!("=== v26.4 c₁ の方向依存性 — 異方性行列 c₁[T_aa](q∥b̂) (第二十七期, 経路 B) ===\n");
    println!("事前登録: (a) S0–S2 PASS → c₁ 行列が主結果 (S3 は branch 記録) /");
    println!("          (b) S0a FAIL → ブロック対導出の誤り / (c) S2 FAIL → 窓の再設計\n");
    let t0 = std::time::Instant::now();
    let nthreads = std::thread::available_parallelism()
        .map(|x| x.get())
        .unwrap_or(4);
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

    // ---- [S0a] 内蔵 dense 照合 (q∥ŷ; N=8, m∈{0,0.5}, j∈{1,2}) ----
    {
        let mut worst = 0.0f64;
        for &m in &[0.0f64, 0.5] {
            let n = 8usize;
            let ns = n * n * n;
            let h = build_h_dense(n, m);
            let (w, v) = jacobi_eigh(&h, ns);
            let blocky = chi_scan_y(n, m, &[1, 2], &[0, 2], nthreads, false);
            if m == 0.0 {
                let pieces = chi_scan_y(n, m, &[1], &[1, 3], nthreads, false);
                println!(
                    "    [S0a piece probe] block: x片 = {:.6} (numpy 0.097061) / z片 = {:.6} (numpy 0.095269)",
                    pieces[0][0], pieces[1][0]
                );
            }
            for (ji, &j) in [1usize, 2].iter().enumerate() {
                let q = [0.0, 2.0 * PI * j as f64 / n as f64, 0.0];
                for (wi, &which) in [0usize, 2].iter().enumerate() {
                    let (ore, oim) = dense_vertex(n, m, q, which);
                    let cd = chi_dense(&w, &v, ns, &ore, &oim) / ns as f64;
                    println!(
                        "    [S0a probe] m={:.1} j={} which={}  block = {:.6}  dense = {:.6}",
                        m, j, which, blocky[wi][ji], cd
                    );
                    worst = worst.max((blocky[wi][ji] - cd).abs());
                }
            }
        }
        check(
            "[S0a] q∥ŷ ブロック対エンジン: 内蔵 dense と一致 (abs 1e-9)",
            worst < 1e-9,
            format!("max|Δ| = {:.1e} ({} s)", worst, t0.elapsed().as_secs()),
        );
    }

    // ---- [S0b] x 経路の回帰 (REF262 + v26.3 の c₁^xx) ----
    {
        let mut worst = 0.0f64;
        for &(n, m10) in &[(8usize, 0u32), (8, 5), (12, 0), (12, 5)] {
            let m = m10 as f64 / 10.0;
            let js: Vec<usize> = (1..=(n / 4)).collect();
            let g00 = chi_scan_x(n, m, &js, 0, nthreads);
            let gxx = chi_scan_x(n, m, &js, 1, nthreads);
            for &(rn, rm, rj, r00, rxx) in REF262.iter().filter(|r| r.0 == n && r.1 == m10) {
                worst = worst
                    .max((g00[rj - 1] - r00).abs())
                    .max((gxx[rj - 1] - rxx).abs());
                let _ = (rn, rm);
            }
        }
        check(
            "[S0b-1] x 経路: REF262 の 20 値を再現 (abs 1e-6)",
            worst < 1e-6,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- 走査: c₁ 行列の 4 成分 (N ∈ {32,64}, m ∈ {0, 0.5}) ----
    let ns_list = [32usize, 64];
    let ms = [0.0f64, 0.5];
    let jmax = 6usize;
    let js_all: Vec<usize> = (0..=jmax).collect();
    // tab[ni][mi][combo][j]: combo 0 = xx@x̂, 1 = yy@x̂, 2 = xx@ŷ, 3 = yy@ŷ
    let mut tab = vec![vec![vec![vec![0.0f64; jmax + 1]; 4]; ms.len()]; ns_list.len()];
    for (ni, &n) in ns_list.iter().enumerate() {
        for (mi, &m) in ms.iter().enumerate() {
            let xx_x = chi_scan_x(n, m, &js_all, 1, nthreads);
            let yy_x = chi_scan_x(n, m, &js_all, 2, nthreads);
            let ys = chi_scan_y(n, m, &js_all, &[1, 2], nthreads, false);
            for j in 0..=jmax {
                tab[ni][mi][0][j] = xx_x[j];
                tab[ni][mi][1][j] = yy_x[j];
                tab[ni][mi][2][j] = ys[0][j];
                tab[ni][mi][3][j] = ys[1][j];
            }
            println!(
                "    [走査] N={} m={:.1} 完了 ({} s) — χ_xx(0) = {:.6}, χ_yy(0) = {:.6}",
                n,
                m,
                t0.elapsed().as_secs(),
                tab[ni][mi][0][0],
                tab[ni][mi][1][0]
            );
        }
    }

    // ---- [S1] χ(0) の 2 エンジン一致 (T_yy(0): x 経路 j=0 vs y 経路 j=0) ----
    {
        let mut worst = 0.0f64;
        for ni in 0..ns_list.len() {
            for mi in 0..ms.len() {
                worst = worst.max((tab[ni][mi][1][0] - tab[ni][mi][3][0]).abs());
            }
        }
        check(
            "[S1] χ_yy(0) の 2 エンジン一致 (同一演算子の恒等検査)",
            worst < 1e-10,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- c₁ フィット (v26.3 プロトコル: 窓 2 種 × モデル 2 種) ----
    let fit_c1 = |ni: usize, mi: usize, combo: usize| -> (f64, f64) {
        let n = ns_list[ni];
        let q = |j: usize| 2.0 * PI * j as f64 / n as f64;
        let mut c1s = Vec::new();
        for (lo, hi) in [(1usize, 4usize), (2, 5)] {
            let jr: Vec<usize> = (lo..=hi).collect();
            let y: Vec<f64> = jr
                .iter()
                .map(|&j| (tab[ni][mi][combo][j] - tab[ni][mi][combo][0]) / (q(j) * q(j)))
                .collect();
            let ones: Vec<f64> = jr.iter().map(|_| 1.0).collect();
            let q2: Vec<f64> = jr.iter().map(|&j| q(j) * q(j)).collect();
            let c = lstsq(&[ones.clone(), q2.clone()], &y);
            c1s.push(c[0]);
            if ms[mi] == 0.0 {
                let ql: Vec<f64> = jr
                    .iter()
                    .map(|&j| q(j) * q(j) * (1.0 / q(j)).ln())
                    .collect();
                let c = lstsq(&[ones, q2, ql], &y);
                c1s.push(c[0]);
            }
        }
        let mut sorted = c1s.clone();
        sorted.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med = sorted[sorted.len() / 2];
        let spread = sorted.last().unwrap() - sorted.first().unwrap();
        (med, spread)
    };
    let names = ["c₁[xx](x̂)", "c₁[yy](x̂)", "c₁[xx](ŷ)", "c₁[yy](ŷ)"];
    println!(
        "\n    [c₁ 行列] N | m | {} | {} | {} | {}",
        names[0], names[1], names[2], names[3]
    );
    let mut c1 = vec![vec![vec![(0.0f64, 0.0f64); 4]; ms.len()]; ns_list.len()];
    for ni in 0..ns_list.len() {
        for mi in 0..ms.len() {
            for combo in 0..4 {
                c1[ni][mi][combo] = fit_c1(ni, mi, combo);
            }
            println!(
                "      N={} m={:.1}:  {:.5}(±{:.5})  {:.5}(±{:.5})  {:.5}(±{:.5})  {:.5}(±{:.5})",
                ns_list[ni],
                ms[mi],
                c1[ni][mi][0].0,
                c1[ni][mi][0].1,
                c1[ni][mi][1].0,
                c1[ni][mi][1].1,
                c1[ni][mi][2].0,
                c1[ni][mi][2].1,
                c1[ni][mi][3].0,
                c1[ni][mi][3].1
            );
        }
    }

    // ---- [S0b-2] v26.3 の c₁^xx(x̂) 回帰 ----
    {
        let mut worst = 0.0f64;
        for &(mi, (_m, r)) in [(0usize, REF263_C1XX[0]), (1usize, REF263_C1XX[2])].iter() {
            worst = worst.max((c1[1][mi][0].0 - r).abs());
        }
        check(
            "[S0b-2] c₁[xx](x̂) が v26.3 公表値と一致 (N=64, ±0.0002)",
            worst < 2e-4,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [S2] 窓安定性と N 収束 (4 成分) ----
    {
        let mut worst_w = 0.0f64;
        let mut worst_n = 0.0f64;
        for mi in 0..ms.len() {
            for combo in 0..4 {
                let (med64, spread64) = c1[1][mi][combo];
                worst_w = worst_w.max(spread64 / med64.abs().max(1e-12));
                let (med32, _) = c1[0][mi][combo];
                worst_n = worst_n.max((med32 / med64 - 1.0).abs());
            }
        }
        check(
            "[S2] c₁ 行列の窓安定性 (< 30%) と N=32→64 収束 (< 20%)",
            worst_w < 0.30 && worst_n < 0.20,
            format!(
                "窓 max = {:.1}% / N 移動 max = {:.1}%",
                100.0 * worst_w,
                100.0 * worst_n
            ),
        );
    }

    // ---- [S3] 異方性 branch ----
    {
        for (label, ca, cb) in [
            ("縦組 c₁[xx](x̂) vs c₁[yy](ŷ)", c1[1][1][0], c1[1][1][3]),
            ("横組 c₁[yy](x̂) vs c₁[xx](ŷ)", c1[1][1][1], c1[1][1][2]),
        ] {
            let diff = (ca.0 - cb.0).abs() / ca.0.abs().max(cb.0.abs());
            let syst = (ca.1 + cb.1) / ca.0.abs().max(cb.0.abs());
            println!(
                "    [S3 branch] {} (m=0.5): {:.5} vs {:.5} — 差 {:.1}% vs 系統 {:.1}% ⇒ {}",
                label,
                ca.0,
                cb.0,
                100.0 * diff,
                100.0 * syst,
                if diff > 3.0 * syst {
                    "branch A: 方向依存あり"
                } else {
                    "branch B: 等方 (系統内)"
                }
            );
        }
        check(
            "[S3] 異方性判定が分解能を持つ (縦組の差か系統が 1% 超)",
            {
                let (ca, cb) = (c1[1][1][0], c1[1][1][3]);
                let diff = (ca.0 - cb.0).abs() / ca.0.abs().max(cb.0.abs());
                let syst = (ca.1 + cb.1) / ca.0.abs().max(cb.0.abs());
                diff > 0.01 || syst < 0.01
            },
            String::new(),
        );
    }

    // ---- [S4] 変異検出: 折返しスワップ落とし ----
    {
        let (n, m) = (8usize, 0.0f64);
        let good = chi_scan_y(n, m, &[2], &[0], nthreads, false);
        let bad = chi_scan_y(n, m, &[2], &[0], nthreads, true);
        let dev = (good[0][0] - bad[0][0]).abs();
        check(
            "[S4] 変異検出: 折返しスワップ落とし → S0a 照合が検出",
            dev > 1e-4,
            format!("逸脱 {:.2e} > 1e-4", dev),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.4".into())),
        ("kind".into(), Json::Str("c1_anisotropy_matrix".into())),
        (
            "c1_matrix_n64".into(),
            Json::Arr(
                (0..ms.len())
                    .map(|mi| {
                        Json::Obj(vec![
                            ("m".into(), Json::Num(ms[mi])),
                            ("xx_qx".into(), Json::Num(c1[1][mi][0].0)),
                            ("yy_qx".into(), Json::Num(c1[1][mi][1].0)),
                            ("xx_qy".into(), Json::Num(c1[1][mi][2].0)),
                            ("yy_qy".into(), Json::Num(c1[1][mi][3].0)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v264_c1_anisotropy.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **c₁ 行列が主結果** — branch は [S3] の欄が一次ソース。解釈は docs/uft-v26.4.md へ"
        } else {
            "FAIL — 分岐 (b)/(c) は各検査の欄を一次ソースとする"
        }
    );
    println!(
        "\n総合判定: {} ({} s)",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" },
        t0.elapsed().as_secs()
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
