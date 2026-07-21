//! 3+1D staggered 半空間の厳密ブロック理論 (第二十五期, v24.1 で導入)
//!
//! 対象は v22.2/v23.x 系の格子 (x 開放・y/z 周期・N 偶数, ホップ ±1/2):
//!   H = D_x + D_y + D_z,
//!   D_x = (T_x+T_x†)/2, D_y = (−1)^x (T_y+T_y†)/2, D_z = (−1)^{x+y} (T_z+T_z†)/2.
//!
//! 鍵となる厳密構造 (v24.1 が dense jacobi との多重集合照合で検証):
//!   (0) H = D_x + (−1)^x M,  M = (T_y+T_y†)/2 + (−1)^y (T_z+T_z†)/2 は x 非依存。
//!   (1) 横ブロック: 基底 {e^{i·ky·y}, e^{i(ky+π)y}} ⊗ e^{i·kz·z} で M は 2×2 実対称
//!       [[cos ky, cos kz],[cos kz, −cos ky]] → μ± = ±√(cos²ky + cos²kz)。
//!       (−1)^y が ky ↔ ky+π を結ぶため、ky は半 BZ (q = 0..N/2−1) を走る。
//!   (2) x 鎖 (横固有値 μ ごと): 開放鎖 sine モード φ_n(x) = √(2/(N+1))·sin(πn(x+1)/(N+1))
//!       に対する厳密恒等式 (−1)^x φ_n(x) = φ_{N+1−n}(x) により、ペア {n, N+1−n}
//!       (n = 1..N/2) の 2×2 [[cos k_n, μ],[μ, −cos k_n]] に閉じる →
//!       E±(n, μ) = ±√(cos²k_n + μ²) — **開放境界でも厳密**。
//!   (3) よってバルク分散は E(k) = ±√(cos²kx + cos²ky + cos²kz):
//!       零集合は cos kx = cos ky = cos kz = 0 の **8 点ノード (余次元 3)** =
//!       2 taste × 4 成分 Dirac。Fermi 面 (余次元 1) ではない。
//!
//! 帰結: 半空間 A = {x < N/2} の相関行列 C_A は (q, p) ラベルの N²/2 個の
//! N×N 実対称ブロックに厳密分解され、dense N³ 対角化 (v23.4 の限界 N=16) を
//! 経ずに N = 128 級へ到達できる。ブロック行列式は全て解析形 (sine と 2×2) —
//! DD 経路 (dd.rs) なら C の固有値を c ≳ 1e-30 まで解像する。
//!
//! ブロック基底の添字: (x, a) → 2x + a。a = 0 は ky 成分, a = 1 は ky+π 成分。
//! 実空間ボンドの復元 (v24.2 で dense と照合):
//!   K_x(x → x+1; y) = A_x(x) + (−1)^y B_x(x),
//!     A_x = (1/N²) Σ_{q,p} [K00 + K11](x, x+1),  B_x = (1/N²) Σ [K01 + K10](x, x+1)
//!   K_y(x; y→y+1)   = (1/N²) Σ cos(ky)·[K00 − K11](x, x)   (y 非依存 — 厳密)
//!   K_z(x; z→z+1; y)= (1/N²) Σ cos(kz)·{[K00+K11](x,x) + (−1)^y·2·K01(x,x)}
//! BW 予言 (2π x T⁰⁰ の格子転写): A_x(ξ) → π·ξ, B_x → 0, K_y → π·ξ·(−1)^x,
//! z は (−1)^y 成分 (η_z = (−1)^{x+y}) が π·ξ·(−1)^x, y 一様成分 → 0。

use crate::dd::Real;

/// dense 3D staggered H (x 開放, y/z 周期) — v22.2/v23.4/v23.6 と同一規約 (照合用)
pub fn build_h3d_open_x(n: usize) -> Vec<f64> {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut h = vec![0.0f64; ns * ns];
    let mut add = |i: usize, j: usize, t: f64| {
        h[j + i * ns] += t;
        h[i + j * ns] += t;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                if x + 1 < n {
                    add(i, idx(x + 1, y, z), 0.5);
                }
                let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                add(i, idx(x, (y + 1) % n, z), ey);
                let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                add(i, idx(x, y, (z + 1) % n), ez);
            }
        }
    }
    h
}

/// dense D_x, D_y, D_z (全周期 or x 開放) — 反交換・H² 検証用 (v24.1)
pub fn build_d_matrices(n: usize, open_x: bool) -> [Vec<f64>; 3] {
    let ns = n * n * n;
    let idx = |x: usize, y: usize, z: usize| x + n * (y + n * z);
    let mut dx = vec![0.0f64; ns * ns];
    let mut dy = vec![0.0f64; ns * ns];
    let mut dz = vec![0.0f64; ns * ns];
    let add = |m: &mut Vec<f64>, i: usize, j: usize, t: f64| {
        m[j + i * ns] += t;
        m[i + j * ns] += t;
    };
    for x in 0..n {
        for y in 0..n {
            for z in 0..n {
                let i = idx(x, y, z);
                if x + 1 < n {
                    add(&mut dx, i, idx(x + 1, y, z), 0.5);
                } else if !open_x {
                    add(&mut dx, i, idx(0, y, z), 0.5);
                }
                let ey = if x % 2 == 0 { 0.5 } else { -0.5 };
                add(&mut dy, i, idx(x, (y + 1) % n, z), ey);
                let ez = if (x + y) % 2 == 0 { 0.5 } else { -0.5 };
                add(&mut dz, i, idx(x, y, (z + 1) % n), ez);
            }
        }
    }
    [dx, dy, dz]
}

/// x 開放格子の解析スペクトル (多重集合, 昇順):
/// {± √(cos²(πn/(N+1)) + cos²(2πqy/N) + cos²(2πqz/N))}, n = 1..N/2, qy,qz ∈ 0..N−1
pub fn spectrum_analytic_open(n: usize) -> Vec<f64> {
    let mut ev = Vec::with_capacity(n * n * n);
    for np in 1..=n / 2 {
        let c1 = (std::f64::consts::PI * np as f64 / (n as f64 + 1.0)).cos();
        for qy in 0..n {
            let cy = (2.0 * std::f64::consts::PI * qy as f64 / n as f64).cos();
            for qz in 0..n {
                let cz = (2.0 * std::f64::consts::PI * qz as f64 / n as f64).cos();
                let e = (c1 * c1 + cy * cy + cz * cz).sqrt();
                ev.push(e);
                ev.push(-e);
            }
        }
    }
    ev.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ev
}

/// 全周期格子の H² の解析スペクトル (多重集合, 昇順): {Σᵢ cos²(2πqᵢ/N)}
pub fn h2_spectrum_periodic(n: usize) -> Vec<f64> {
    let mut ev = Vec::with_capacity(n * n * n);
    let c: Vec<f64> = (0..n)
        .map(|q| (2.0 * std::f64::consts::PI * q as f64 / n as f64).cos())
        .collect();
    for qx in 0..n {
        for qy in 0..n {
            for qz in 0..n {
                ev.push(c[qx] * c[qx] + c[qy] * c[qy] + c[qz] * c[qz]);
            }
        }
    }
    ev.sort_by(|a, b| a.partial_cmp(b).unwrap());
    ev
}

/// x 開放格子の解析ギャップ (HOMO–LUMO = 2·min E) — シェル構造 (N mod 4) の起源
pub fn gap_analytic(n: usize) -> f64 {
    let sp = spectrum_analytic_open(n);
    let mut emin = f64::MAX;
    for &e in &sp {
        if e > 0.0 && e < emin {
            emin = e;
        }
    }
    2.0 * emin
}

/// ブロックラベル (q, p): ky = 2πq/N (q = 0..N/2−1, {ky, ky+π} ペア代表), kz = 2πp/N
pub fn blocks(n: usize) -> Vec<(usize, usize)> {
    let mut v = Vec::with_capacity(n * n / 2);
    for q in 0..n / 2 {
        for p in 0..n {
            v.push((q, p));
        }
    }
    v
}

/// ブロック (q,p) の占有モード行列 F (2N 行 × N 列, 列優先 f[row + col·2N])。
/// 行 = 2x + a (a=0: ky 成分, a=1: ky+π 成分)。列 = 占有モード (E<0)。
/// 全て解析形 (sine + 2×2 固有問題) — 型 T の全精度で構成。
pub fn block_f<T: Real>(n: usize, q: usize, p: usize) -> Vec<T> {
    let ni = n as i64;
    let cy = T::cospi_frac(2 * q as i64, ni);
    let cz = T::cospi_frac(2 * p as i64, ni);
    let rho = (cy * cy + cz * cz).sqrtr();
    let norm_x = (T::from_f64(2.0).divr(T::from_f64(n as f64 + 1.0))).sqrtr();
    let mut f = vec![T::R0; 2 * n * n];
    let mut col = 0usize;
    for band in 0..2usize {
        // 横固有対 (μ, (u,w))
        let (mu, u, w) = if rho.hi() == 0.0 {
            // 縮退規約: (1,0) と (0,1)
            if band == 0 {
                (T::R0, T::R1, T::R0)
            } else {
                (T::R0, T::R0, T::R1)
            }
        } else {
            let s_plus = band == 0;
            let mu = if s_plus { rho } else { -rho };
            let (u0, w0) = if cy.hi() >= 0.0 {
                if s_plus {
                    (rho + cy, cz)
                } else {
                    (cz, -(rho + cy))
                }
            } else {
                if s_plus {
                    (cz, rho - cy)
                } else {
                    (rho - cy, -cz)
                }
            };
            let nrm = (u0 * u0 + w0 * w0).sqrtr();
            (mu, u0.divr(nrm), w0.divr(nrm))
        };
        for np in 1..=n / 2 {
            // x ペア {np, N+1−np}: [[c1, μ],[μ, −c1]] の E₋ = −√(c1²+μ²)
            let c1 = T::cospi_frac(np as i64, ni + 1);
            let r = (c1 * c1 + mu * mu).sqrtr();
            let a0 = mu;
            let b0 = -(r + c1);
            let nrm = (a0 * a0 + b0 * b0).sqrtr();
            let alpha = a0.divr(nrm);
            let beta = b0.divr(nrm);
            for x in 0..n {
                let m1 = (np as i64) * (x as i64 + 1);
                let m2 = (ni + 1 - np as i64) * (x as i64 + 1);
                let fx = (alpha * T::sinpi_frac(m1, ni + 1)
                    + beta * T::sinpi_frac(m2, ni + 1))
                    * norm_x;
                f[2 * x + col * 2 * n] = fx * u;
                f[2 * x + 1 + col * 2 * n] = fx * w;
            }
            col += 1;
        }
    }
    debug_assert_eq!(col, n);
    f
}

/// F から x 部分領域 xsel 上の相関行列 C = F_sel F_selᵀ (次元 2·|xsel|, 列優先)。
/// 添字 (i, a) → 2i + a (i は xsel 内の順位)。
pub fn c_from_f<T: Real>(f: &[T], n: usize, xsel: &[usize]) -> Vec<T> {
    let dim = 2 * xsel.len();
    let nfill = n;
    let rows: Vec<usize> = xsel
        .iter()
        .flat_map(|&x| [2 * x, 2 * x + 1])
        .collect();
    let mut c = vec![T::R0; dim * dim];
    for k in 0..nfill {
        let fk = &f[k * 2 * n..(k + 1) * 2 * n];
        for a in 0..dim {
            let va = fk[rows[a]];
            if va.hi() == 0.0 {
                continue;
            }
            for b in a..dim {
                c[a + b * dim] = c[a + b * dim] + va * fk[rows[b]];
            }
        }
    }
    for a in 0..dim {
        for b in 0..a {
            c[a + b * dim] = c[b + a * dim];
        }
    }
    c
}

/// ブロック K (xsel = 連続 x 範囲, 次元 2·nx) からボンド観測量のブロック寄与を抽出。
/// 実空間復元の重み (1/N² と cos 因子) は呼び出し側の全ブロック和で付ける。
pub struct BlockBondSums {
    /// x-NN (i→i+1) の a 対角和 [K00+K11] — taste 一様成分 A_x の素材
    pub x_diag: Vec<f64>,
    /// x-NN の a 反対角和 [K01+K10] — taste 交替成分 B_x の素材
    pub x_off: Vec<f64>,
    /// x-NNN (i→i+2) の a 対角和 (演算子整合用)
    pub x2_diag: Vec<f64>,
    /// x-NNN の a 反対角和
    pub x2_off: Vec<f64>,
    /// サイト対角 (i) の [K00 − K11] — y ボンド復元は cos(ky)·これ
    pub on_sig: Vec<f64>,
    /// サイト対角 (i) の [K00 + K11] — z ボンドの y 一様成分は cos(kz)·これ
    pub on_uni: Vec<f64>,
    /// サイト (i) の 2·K01 — z ボンドの (−1)^y 成分は cos(kz)·これ
    pub on_alt: Vec<f64>,
}

pub fn block_bond_sums<T: Real>(k: &[T], nx: usize) -> BlockBondSums {
    let dim = 2 * nx;
    let el = |i: usize, a: usize, j: usize, b: usize| k[(2 * i + a) + (2 * j + b) * dim].hi();
    let mut s = BlockBondSums {
        x_diag: vec![0.0; nx.saturating_sub(1)],
        x_off: vec![0.0; nx.saturating_sub(1)],
        x2_diag: vec![0.0; nx.saturating_sub(2)],
        x2_off: vec![0.0; nx.saturating_sub(2)],
        on_sig: vec![0.0; nx],
        on_uni: vec![0.0; nx],
        on_alt: vec![0.0; nx],
    };
    for i in 0..nx {
        s.on_sig[i] = el(i, 0, i, 0) - el(i, 1, i, 1);
        s.on_uni[i] = el(i, 0, i, 0) + el(i, 1, i, 1);
        s.on_alt[i] = 2.0 * el(i, 0, i, 1);
        if i + 1 < nx {
            s.x_diag[i] = el(i, 0, i + 1, 0) + el(i, 1, i + 1, 1);
            s.x_off[i] = el(i, 0, i + 1, 1) + el(i, 1, i + 1, 0);
        }
        if i + 2 < nx {
            s.x2_diag[i] = el(i, 0, i + 2, 0) + el(i, 1, i + 2, 1);
            s.x2_off[i] = el(i, 0, i + 2, 1) + el(i, 1, i + 2, 0);
        }
    }
    s
}

/// ブロック (q,p) の相関行列 → モジュラー核 K (次元 2·|xsel|) を一括構成
pub fn block_k<T: Real>(n: usize, q: usize, p: usize, xsel: &[usize], clamp: f64) -> Vec<T> {
    let f = block_f::<T>(n, q, p);
    let c = c_from_f(&f, n, xsel);
    let (k, _) = crate::dd::modular_k(&c, 2 * xsel.len(), 60, clamp);
    k
}

/// 半空間走査の集計 (実空間復元済み, 全ブロック和 × 1/N²)
pub struct HalfScan {
    /// エンタングルメントエントロピー S(C_A)
    pub s_total: f64,
    /// x-NN ボンドの taste 一様成分 A_x(i), i→i+1 (BW 予言: π·ξ_b)
    pub ax: Vec<f64>,
    /// x-NN の taste 交替成分 B_x(i) ((−1)^y 係数; BW 予言: 0)
    pub bx: Vec<f64>,
    /// x-NNN (i→i+2) の一様/交替成分 (演算子整合用)
    pub ax2: Vec<f64>,
    pub bx2: Vec<f64>,
    /// y ボンド K_y(i) (y 非依存; BW 予言: π·ξ_s·(−1)^i·η 形)
    pub ky: Vec<f64>,
    /// z ボンドの y 一様成分 (BW 予言: 0) と (−1)^y 成分 (BW 本体)
    pub az: Vec<f64>,
    pub bz: Vec<f64>,
    /// サイト対角の一様/(K00−K11)/交替成分 (診断・整合用)
    pub on_uni: Vec<f64>,
    pub on_sig: Vec<f64>,
    pub on_alt: Vec<f64>,
    /// クランプに達した固有値の本数と、非クランプ域の max |κ|
    pub n_clamped: usize,
    pub kappa_max: f64,
}

/// 半空間 A = {x < N/2} の全ブロック走査。ブロックごとの結果を保存してから
/// ブロック順に縮約するので、**スレッド数に依らず結果は決定的** (PROMPT/4)。
pub fn half_space_scan<T: Real>(n: usize, clamp: f64, nthreads: usize) -> HalfScan {
    let bl = blocks(n);
    let nx = n / 2;
    let xsel: Vec<usize> = (0..nx).collect();
    struct BlockOut {
        s: f64,
        sums: BlockBondSums,
        cy: f64,
        cz: f64,
        n_clamped: usize,
        kappa_max: f64,
    }
    let mut per: Vec<Option<BlockOut>> = Vec::new();
    per.resize_with(bl.len(), || None);
    let nt = nthreads.max(1);
    let chunk = (bl.len() + nt - 1) / nt;
    std::thread::scope(|sc| {
        for (t, sl) in per.chunks_mut(chunk).enumerate() {
            let bl = &bl;
            let xsel = &xsel;
            sc.spawn(move || {
                for (j, slot) in sl.iter_mut().enumerate() {
                    let bi = t * chunk + j;
                    let (q, p) = bl[bi];
                    let f = block_f::<T>(n, q, p);
                    let c = c_from_f(&f, n, xsel);
                    let (k, cw) = crate::dd::modular_k(&c, 2 * nx, 60, clamp);
                    let mut s_blk = 0.0f64;
                    let mut ncl = 0usize;
                    let mut kmax = 0.0f64;
                    for &cev in &cw {
                        let ch = cev.hi().clamp(0.0, 1.0);
                        s_blk += crate::h2_entropy(ch);
                        // クランプ判定と κ は T 精度で行う — f64 では 1.0 − 1e-30 が
                        // 1.0 に丸まり、1 側の深部モードを見落とす (v24.2 の開発記録)
                        let dist0 = cev.hi();
                        let dist1 = (T::R1 - cev).hi();
                        if dist0 < clamp || dist1 < clamp {
                            ncl += 1;
                        } else {
                            let kap = ((T::R1 - cev).divr(cev)).lnr().hi().abs();
                            kmax = kmax.max(kap);
                        }
                    }
                    let sums = block_bond_sums(&k, nx);
                    let cy = crate::dd::dd_cospi_frac(2 * q as i64, n as i64).hi;
                    let cz = crate::dd::dd_cospi_frac(2 * p as i64, n as i64).hi;
                    *slot = Some(BlockOut {
                        s: s_blk,
                        sums,
                        cy,
                        cz,
                        n_clamped: ncl,
                        kappa_max: kmax,
                    });
                }
            });
        }
    });
    let w = 1.0 / (n * n) as f64;
    let mut out = HalfScan {
        s_total: 0.0,
        ax: vec![0.0; nx - 1],
        bx: vec![0.0; nx - 1],
        ax2: vec![0.0; nx.saturating_sub(2)],
        bx2: vec![0.0; nx.saturating_sub(2)],
        ky: vec![0.0; nx],
        az: vec![0.0; nx],
        bz: vec![0.0; nx],
        on_uni: vec![0.0; nx],
        on_sig: vec![0.0; nx],
        on_alt: vec![0.0; nx],
        n_clamped: 0,
        kappa_max: 0.0,
    };
    for b in per.into_iter().map(|o| o.expect("block 未計算")) {
        out.s_total += b.s;
        out.n_clamped += b.n_clamped;
        out.kappa_max = out.kappa_max.max(b.kappa_max);
        for i in 0..nx - 1 {
            out.ax[i] += w * b.sums.x_diag[i];
            out.bx[i] += w * b.sums.x_off[i];
        }
        for i in 0..nx.saturating_sub(2) {
            out.ax2[i] += w * b.sums.x2_diag[i];
            out.bx2[i] += w * b.sums.x2_off[i];
        }
        for i in 0..nx {
            out.ky[i] += w * b.cy * b.sums.on_sig[i];
            out.az[i] += w * b.cz * b.sums.on_uni[i];
            out.bz[i] += w * b.cz * b.sums.on_alt[i];
            out.on_uni[i] += w * b.sums.on_uni[i];
            out.on_sig[i] += w * b.sums.on_sig[i];
            out.on_alt[i] += w * b.sums.on_alt[i];
        }
    }
    out
}

/// stag ブロック理論の内部自己検証 (軽量): ペア恒等式と充填数。
pub fn stag_self_test() -> bool {
    let mut ok = true;
    // (−1)^x sin(πn(x+1)/(N+1)) = sin(π(N+1−n)(x+1)/(N+1)) — DD で < 1e-30
    let n = 8i64;
    for nn in 1..=n {
        for x in 0..n {
            let lhs = crate::dd::dd_sinpi_frac(nn * (x + 1), n + 1)
                .mul_f(if x % 2 == 0 { 1.0 } else { -1.0 });
            let rhs = crate::dd::dd_sinpi_frac((n + 1 - nn) * (x + 1), n + 1);
            ok &= (lhs - rhs).hi.abs() < 1e-30;
        }
    }
    // ブロック占有モードの直交規格性 (N=6, 全ブロック, DD): |FᵀF − 1| < 1e-29
    let nb = 6usize;
    for (q, p) in blocks(nb) {
        let f = block_f::<crate::dd::Dd>(nb, q, p);
        for k1 in 0..nb {
            for k2 in k1..nb {
                let mut ip = crate::dd::DD0;
                for r in 0..2 * nb {
                    ip = ip + f[r + k1 * 2 * nb] * f[r + k2 * 2 * nb];
                }
                let target = if k1 == k2 { 1.0 } else { 0.0 };
                ok &= (ip.hi - target).abs() < 1e-29;
            }
        }
    }
    ok
}
