//! v26.5 spin-2/spin-0 分解 — graviton plus チャネルの c₁ (第二十七期, 経路 B)
//!
//! 誘導重力 10 要件の (7) tensor structure へ: 静的応答を横断偏極で分解する。
//! q ∥ ŷ の横方向指数 {x, z} を使うと、**plus 偏極の graviton チャネル**
//!   D = (T_xx − T_zz)/√2  (spin-2 の + 偏極),  S = (T_xx + T_zz)/√2  (横 spin-0)
//! が **v26.4 で認証済みの演算子だけ**で組める (T_xy 型の奇運動量頂点 — staggered の
//! taste 構造に絡む point-split 問題 — を回避)。新規の器械は交差相関
//! χ_AB = Σ 2·M_A M_B/ΔE のみ (頂点は実)。χ_D = (χ_xx + χ_zz)/2 − χ_xz、
//! χ_S = (χ_xx + χ_zz)/2 + χ_xz。
//!
//! 検査 (事前登録):
//!  [S0a] 交差相関の内蔵 dense 照合: (χ_xx, χ_zz, χ_xz)(qŷ) (N=8, m∈{0,0.5}, j∈{1,2})
//!        を実空間 dense と abs 1e-9 で一致
//!  [S0b] x↔z 対称回帰: c₁[zz](ŷ) と c₁[xx](ŷ) がともに v26.4 の横値
//!        (−0.01556 [m=0] / −0.01264 [m=0.5], ±0.0003) — 横 2 演算子の等価性
//!  [S1] 正値性: χ_D(q) ≥ 0 かつ χ_S(q) ≥ 0 (交差項の Cauchy–Schwarz 整合)
//!  [S2] c₁ の窓安定性 (< 30%) と N=32→64 収束 (< 20%)
//!  [S3] branch: c₁^(2) (graviton plus) と c₁^(0) の符号・大きさ・質量走行の記録
//!  [S4] 変異検出: T_zz 頂点の折返しスワップ落とし → S0a 照合が検出
//!
//! 事前登録分岐: (a) S0–S2 PASS → c₁^(2), c₁^(0) が主結果 / (b) S0 FAIL → 交差器械の
//!   誤り / (c) S2 FAIL → 窓の再設計。

use uft_sim::*;

const PI: f64 = std::f64::consts::PI;

/// v26.4 の横チャネル c₁ 公表値 (N=64; m = 0, 0.5) — x↔z 対称回帰の的
const REF264_TRANS: [(f64, f64); 2] = [(0.0, -0.01556), (0.5, -0.01264)];

// ---------------- ブロック機構 (v26.4 の認証済み実装を写経) ----------------

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

/// q∥ŷ の頂点 (ブロック対, 実行列, 格納 v[標的 + 元·dim])。
/// which: 1 = T_xx (x ボンド片), 3 = T_zz (z ボンド片)。sw = 折返しスワップ。
fn vertex_qy(n: usize, ky: f64, ckz: f64, q: f64, sw: bool, which: usize) -> Vec<f64> {
    let dim = 4 * n;
    let id = |x: usize, c: usize| x + n * c;
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
            if which == 1 {
                v[id((x + 1) % n, tc(c)) + id(x, c) * dim] += 0.5 * tw;
                v[id(x, tc(c)) + id((x + 1) % n, c) * dim] += 0.5 * tw;
            }
            if which == 2 {
                v[id(x, tc(c)) + id(x, c) * dim] += px * (kyc + q / 2.0).cos();
            }
            if which == 3 && (c == 0 || c == 2) {
                let coef = px * zsgn[c] * ckz;
                v[id(x, tc(c + 1)) + id(x, c) * dim] += coef;
                v[id(x, tc(c)) + id(x, c + 1) * dim] += coef;
            }
        }
    }
    v
}

/// ブロック対の交差 Lehmann: (χ_AA, χ_BB, χ_AB)。縮約は「元」側 (v26.4 の教訓:
/// 転置縮約は単片で不可視 — 格納規約 v[標的 + 元·dim] に対し o[r + k·dim] で読む)。
fn chi_pair_cross(
    w1: &[f64],
    v1: &[f64],
    w2: &[f64],
    v2: &[f64],
    dim: usize,
    oa: &[f64],
    ob: &[f64],
) -> (f64, f64, f64) {
    let nocc = dim / 2;
    let tv = |o: &[f64]| -> Vec<f64> {
        let mut t = vec![0.0f64; dim * nocc];
        for ccol in 0..nocc {
            for r in 0..dim {
                let mut s = 0.0;
                for k in 0..dim {
                    s += o[r + k * dim] * v1[k + ccol * dim];
                }
                t[r + ccol * dim] = s;
            }
        }
        t
    };
    let ta = tv(oa);
    let tb = tv(ob);
    let (mut caa, mut cbb, mut cab) = (0.0f64, 0.0f64, 0.0f64);
    for mu in nocc..dim {
        for nu in 0..nocc {
            let (mut ma, mut mb) = (0.0f64, 0.0f64);
            for k in 0..dim {
                let vm = v2[k + mu * dim];
                ma += vm * ta[k + nu * dim];
                mb += vm * tb[k + nu * dim];
            }
            let de = w2[mu] - w1[nu];
            caa += 2.0 * ma * ma / de;
            cbb += 2.0 * mb * mb / de;
            cab += 2.0 * ma * mb / de;
        }
    }
    (caa, cbb, cab)
}

/// kz 行ごとに対角化し (χ_xx, χ_zz, χ_xz)/V を全 j で。
/// mutate = true: T_zz 頂点の折返しスワップを落とす (S4 用)。
fn chi_scan_cross(
    n: usize,
    m: f64,
    js: &[usize],
    nthreads: usize,
    mutate: bool,
) -> Vec<(f64, f64, f64)> {
    let nb = n / 2;
    let mut rows: Vec<Option<Vec<(f64, f64, f64)>>> = Vec::new();
    rows.resize_with(nb, || None);
    let chunk = nb.div_ceil(nthreads);
    std::thread::scope(|sc| {
        for (t, sl) in rows.chunks_mut(chunk).enumerate() {
            sc.spawn(move || {
                for (i, slot) in sl.iter_mut().enumerate() {
                    let jz = t * chunk + i;
                    let ckz = (2.0 * PI * jz as f64 / n as f64).cos();
                    let mut eigs: Vec<(Vec<f64>, Vec<f64>)> = Vec::with_capacity(nb);
                    for jy in 0..nb {
                        let cky = (2.0 * PI * jy as f64 / n as f64).cos();
                        let blk = block_h(n, m, cky, ckz);
                        eigs.push(jacobi_eigh(&blk.h, blk.dim));
                    }
                    let dim = 4 * n;
                    let mut acc = vec![(0.0f64, 0.0f64, 0.0f64); js.len()];
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
                            let ox = vertex_qy(n, ky, ckz, q, sw, 1);
                            let oz = vertex_qy(n, ky, ckz, q, if mutate { false } else { sw }, 3);
                            let (w1, v1) = &eigs[jy];
                            let (w2, v2) = &eigs[jt];
                            let (a, b, c) = chi_pair_cross(w1, v1, w2, v2, dim, &ox, &oz);
                            acc[ji].0 += a;
                            acc[ji].1 += b;
                            acc[ji].2 += c;
                        }
                    }
                    *slot = Some(acc);
                }
            });
        }
    });
    let ns3 = (n * n * n) as f64;
    let mut out = vec![(0.0f64, 0.0f64, 0.0f64); js.len()];
    for r in rows {
        for (ji, &(a, b, c)) in r.unwrap().iter().enumerate() {
            out[ji].0 += a / ns3;
            out[ji].1 += b / ns3;
            out[ji].2 += c / ns3;
        }
    }
    out
}

// ---------------- 内蔵 dense 照合 ----------------

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

/// dense の T_xx(qŷ) / T_zz(qŷ) (which: 1 = x ボンド片, 3 = z ボンド片)
fn dense_piece(n: usize, qy: f64, which: usize) -> (Vec<f64>, Vec<f64>) {
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
                let ph = qy * y as f64;
                if which == 1 {
                    let tw = if x == n - 1 { -1.0 } else { 1.0 };
                    addc(&mut re, &mut im, i, idx((x + 1) % n, y, z), 0.5 * tw, ph);
                }
                if which == 3 {
                    let ez = if (x + y) % 2 == 0 { 1.0 } else { -1.0 };
                    addc(&mut re, &mut im, i, idx(x, y, (z + 1) % n), 0.5 * ez, ph);
                }
            }
        }
    }
    (re, im)
}

/// dense の交差 Lehmann (複素頂点対): (χ_AA, χ_BB, χ_AB) — χ_AB = Σ 2Re[M_A* M_B]/ΔE
fn chi_dense_cross(
    w: &[f64],
    v: &[f64],
    ns: usize,
    a: &(Vec<f64>, Vec<f64>),
    b: &(Vec<f64>, Vec<f64>),
) -> (f64, f64, f64) {
    let nocc = ns / 2;
    let tv = |o: &(Vec<f64>, Vec<f64>)| -> (Vec<f64>, Vec<f64>) {
        let mut tre = vec![0.0f64; ns * nocc];
        let mut tim = vec![0.0f64; ns * nocc];
        for c in 0..nocc {
            for r in 0..ns {
                let (mut sr, mut si) = (0.0, 0.0);
                for k in 0..ns {
                    let vv = v[k + c * ns];
                    sr += o.0[k + r * ns] * vv;
                    si += o.1[k + r * ns] * vv;
                }
                tre[r + c * ns] = sr;
                tim[r + c * ns] = si;
            }
        }
        (tre, tim)
    };
    let (are, aim) = tv(a);
    let (bre, bim) = tv(b);
    let (mut caa, mut cbb, mut cab) = (0.0f64, 0.0f64, 0.0f64);
    for mu in nocc..ns {
        for nu in 0..nocc {
            let (mut mar, mut mai, mut mbr, mut mbi) = (0.0f64, 0.0f64, 0.0f64, 0.0f64);
            for k in 0..ns {
                let vm = v[k + mu * ns];
                mar += vm * are[k + nu * ns];
                mai += vm * aim[k + nu * ns];
                mbr += vm * bre[k + nu * ns];
                mbi += vm * bim[k + nu * ns];
            }
            let de = w[mu] - w[nu];
            caa += 2.0 * (mar * mar + mai * mai) / de;
            cbb += 2.0 * (mbr * mbr + mbi * mbi) / de;
            cab += 2.0 * (mar * mbr + mai * mbi) / de;
        }
    }
    (caa, cbb, cab)
}

/// 最小二乗 (v26.3/26.4 写経)
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
    println!(
        "=== v26.5 spin-2/spin-0 分解 — graviton plus チャネルの c₁ (第二十七期, 経路 B) ===\n"
    );
    println!("事前登録: (a) S0–S2 PASS → c₁^(2), c₁^(0) が主結果 (S3 は branch 記録) /");
    println!("          (b) S0 FAIL → 交差器械の誤り / (c) S2 FAIL → 窓の再設計\n");
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

    // ---- [S0a] 交差相関の内蔵 dense 照合 (N=8) ----
    {
        let mut worst = 0.0f64;
        for &m in &[0.0f64, 0.5] {
            let n = 8usize;
            let ns = n * n * n;
            let h = build_h_dense(n, m);
            let (w, v) = jacobi_eigh(&h, ns);
            let blocky = chi_scan_cross(n, m, &[1, 2], nthreads, false);
            for (ji, &j) in [1usize, 2].iter().enumerate() {
                let qy = 2.0 * PI * j as f64 / n as f64;
                let a = dense_piece(n, qy, 1);
                let b = dense_piece(n, qy, 3);
                let (dxx, dzz, dxz) = chi_dense_cross(&w, &v, ns, &a, &b);
                let (bx, bz, bc) = blocky[ji];
                worst = worst
                    .max((bx - dxx / ns as f64).abs())
                    .max((bz - dzz / ns as f64).abs())
                    .max((bc - dxz / ns as f64).abs());
            }
        }
        check(
            "[S0a] 交差相関 (χ_xx, χ_zz, χ_xz)(qŷ): 内蔵 dense と一致 (abs 1e-9)",
            worst < 1e-9,
            format!("max|Δ| = {:.1e} ({} s)", worst, t0.elapsed().as_secs()),
        );
    }

    // ---- 走査 (N ∈ {32, 64}, m ∈ {0, 0.5}) ----
    let ns_list = [32usize, 64];
    let ms = [0.0f64, 0.5];
    let jmax = 6usize;
    let js_all: Vec<usize> = (0..=jmax).collect();
    let mut tab = vec![vec![vec![(0.0f64, 0.0f64, 0.0f64); jmax + 1]; ms.len()]; ns_list.len()];
    for (ni, &n) in ns_list.iter().enumerate() {
        for (mi, &m) in ms.iter().enumerate() {
            let got = chi_scan_cross(n, m, &js_all, nthreads, false);
            for j in 0..=jmax {
                tab[ni][mi][j] = got[j];
            }
            println!(
                "    [走査] N={} m={:.1} 完了 ({} s) — χ_xx(0) = {:.6}, χ_xz(0) = {:.6}",
                n,
                m,
                t0.elapsed().as_secs(),
                tab[ni][mi][0].0,
                tab[ni][mi][0].2
            );
        }
    }

    // ---- [S1] 正値性 (χ_D, χ_S ≥ 0) ----
    {
        let mut worst = f64::INFINITY;
        for ni in 0..ns_list.len() {
            for mi in 0..ms.len() {
                for j in 0..=jmax {
                    let (xx, zz, xz) = tab[ni][mi][j];
                    let d = 0.5 * (xx + zz) - xz;
                    let s = 0.5 * (xx + zz) + xz;
                    worst = worst.min(d).min(s);
                }
            }
        }
        check(
            "[S1] 正値性: χ_D(q) ≥ 0 かつ χ_S(q) ≥ 0 (交差項の整合)",
            worst > -1e-10,
            format!("min = {:.1e}", worst),
        );
    }

    // ---- c₁ フィット (チャネル: xx, zz, D = spin-2 plus, S = 横 spin-0) ----
    let chan = |t: (f64, f64, f64), ch: usize| -> f64 {
        match ch {
            0 => t.0,
            1 => t.1,
            2 => 0.5 * (t.0 + t.1) - t.2,
            _ => 0.5 * (t.0 + t.1) + t.2,
        }
    };
    let fit_c1 = |ni: usize, mi: usize, ch: usize| -> (f64, f64) {
        let n = ns_list[ni];
        let q = |j: usize| 2.0 * PI * j as f64 / n as f64;
        let mut c1s = Vec::new();
        for (lo, hi) in [(1usize, 4usize), (2, 5)] {
            let jr: Vec<usize> = (lo..=hi).collect();
            let y: Vec<f64> = jr
                .iter()
                .map(|&j| (chan(tab[ni][mi][j], ch) - chan(tab[ni][mi][0], ch)) / (q(j) * q(j)))
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
    println!(
        "\n    [c₁ 分解表] N | m | c₁[xx](ŷ) | c₁[zz](ŷ) | c₁^(2) (plus) | c₁^(0) (横トレース)"
    );
    let mut c1 = vec![vec![[(0.0f64, 0.0f64); 4]; ms.len()]; ns_list.len()];
    for ni in 0..ns_list.len() {
        for mi in 0..ms.len() {
            for ch in 0..4 {
                c1[ni][mi][ch] = fit_c1(ni, mi, ch);
            }
            println!(
                "      N={} m={:.1}:  {:+.5}(±{:.5})  {:+.5}(±{:.5})  {:+.5}(±{:.5})  {:+.5}(±{:.5})",
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

    // ---- [S0b] x↔z 対称回帰 (v26.4 の横値) ----
    {
        let mut worst = 0.0f64;
        for (mi, &(_m, r)) in REF264_TRANS.iter().enumerate() {
            worst = worst
                .max((c1[1][mi][0].0 - r).abs())
                .max((c1[1][mi][1].0 - r).abs());
        }
        check(
            "[S0b] x↔z 対称回帰: c₁[xx](ŷ), c₁[zz](ŷ) = v26.4 の横値 (±0.0003)",
            worst < 3e-4,
            format!("max|Δ| = {:.1e}", worst),
        );
    }

    // ---- [S2] 窓安定性と N 収束 ----
    {
        let mut worst_w = 0.0f64;
        let mut worst_n = 0.0f64;
        for mi in 0..ms.len() {
            for ch in 0..4 {
                let (med64, spread64) = c1[1][mi][ch];
                worst_w = worst_w.max(spread64 / med64.abs().max(1e-12));
                let (med32, _) = c1[0][mi][ch];
                worst_n = worst_n.max((med32 / med64 - 1.0).abs());
            }
        }
        check(
            "[S2] c₁ 分解の窓安定性 (< 30%) と N=32→64 収束 (< 20%)",
            worst_w < 0.30 && worst_n < 0.20,
            format!(
                "窓 max = {:.1}% / N 移動 max = {:.1}%",
                100.0 * worst_w,
                100.0 * worst_n
            ),
        );
    }

    // ---- [S3] branch 記録 ----
    {
        let (d0, sd0) = c1[1][0][2];
        let (d5, sd5) = c1[1][1][2];
        let (s5, _ss5) = c1[1][1][3];
        println!(
            "    [S3 branch] c₁^(2): {:+.5} (m=0) → {:+.5} (m=0.5) — 質量走行 {:.1}% (系統 {:.1}%)",
            d0,
            d5,
            100.0 * (d5 / d0 - 1.0).abs(),
            100.0 * (sd0.abs() / d0.abs()).max(sd5.abs() / d5.abs())
        );
        println!(
            "    [S3 branch] m=0.5: c₁^(2) = {:+.5} / c₁^(0) = {:+.5} — graviton チャネルの符号は {}",
            d5,
            s5,
            if d5 > 0.0 {
                "正 (非 ghost 的)"
            } else {
                "負 (要解釈)"
            }
        );
        check(
            "[S3] spin-2 チャネルの判定が分解能を持つ (|c₁^(2)| が系統の 3 倍超)",
            d5.abs() > 3.0 * sd5.abs(),
            format!("|c₁^(2)| = {:.5}, 系統 = {:.5}", d5.abs(), sd5.abs()),
        );
    }

    // ---- [S4] 変異検出 ----
    {
        let (n, m) = (8usize, 0.0f64);
        let ns = n * n * n;
        let h = build_h_dense(n, m);
        let (w, v) = jacobi_eigh(&h, ns);
        let qy = 2.0 * PI * 2.0 / n as f64;
        let a = dense_piece(n, qy, 1);
        let b = dense_piece(n, qy, 3);
        let (_dxx, _dzz, dxz) = chi_dense_cross(&w, &v, ns, &a, &b);
        let bad = chi_scan_cross(n, m, &[2], nthreads, true);
        let dev = (bad[0].2 - dxz / ns as f64).abs();
        check(
            "[S4] 変異検出: T_zz 頂点の折返しスワップ落とし → S0a 照合が検出",
            dev > 1e-4,
            format!("逸脱 {:.2e} > 1e-4", dev),
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v26.5".into())),
        ("kind".into(), Json::Str("spin2_spin0_decomposition".into())),
        (
            "c1_n64".into(),
            Json::Arr(
                (0..ms.len())
                    .map(|mi| {
                        Json::Obj(vec![
                            ("m".into(), Json::Num(ms[mi])),
                            ("xx_qy".into(), Json::Num(c1[1][mi][0].0)),
                            ("zz_qy".into(), Json::Num(c1[1][mi][1].0)),
                            ("spin2_plus".into(), Json::Num(c1[1][mi][2].0)),
                            ("spin0_trans".into(), Json::Num(c1[1][mi][3].0)),
                        ])
                    })
                    .collect(),
            ),
        ),
    ]);
    let p = write_artifact("results/v265_spin2.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "事前登録 (a): **c₁^(2) (graviton plus), c₁^(0) が主結果** — branch は [S3] の欄が一次ソース。解釈は docs/uft-v26.5.md へ"
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
