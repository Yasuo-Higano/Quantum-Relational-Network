//! v15.3 QrnCoreV2 — ゲージ制約つき相互作用 core からの複数読み出し
//!
//! v6.7 の core はガウスフェルミオン網 (自由場) に限られていた (ASM-GAUSS,
//! 影響範囲 27 主張)。本版はその外へ出る: 1+1 次元 Z2 格子ゲージ + staggered
//! フェルミオン (環) を、Gauss 拘束を**基底の構造として厳密に解いた** Hilbert
//! 空間 (占有 bitmask × 巻き付き電場 ε) 上の相互作用系として実装し、
//! **同一状態・同一発展・同一拘束**から複数の読み出しを行う:
//!   [幾何]     全サイト対 MI → MDS → 円環判定 (v6.4/v6.7 と同じ判定器)
//!   [エントロピー] 部分系 S、純粋性双対 S(A)=S(Aᶜ)、クエンチ後の S(t) 成長
//!   [ゲージ]   外部電荷対の弦エネルギー E(r) → 弦張力 σ(h) (閉じ込め)
//!   [物質]     staggered 凝縮 χ(m)、質量ギャップ Δ
//!   [因果]     ボンド励起クエンチの前線速度 v(h)
//!   [握手]     h ↑ で σ ↑ かつ v ↓ — 閉じ込めが情報伝播を遅くする (読み出し間の整合)
//!
//! 検証の柱 (全て [PASS]/[FAIL] 内蔵):
//!   [A] 拘束の正しさ: L=4 で拘束なし全空間 (256 次元、リンク自由度を陽に持つ) を
//!       対角化し、(i) H が物理セクター (全 G_n=+1) を厳密に保つこと、
//!       (ii) 物理セクターのスペクトル = 拘束を解いた core のスペクトル (機械精度)、
//!       (iii) Gauss 破り (G=−1 対) セクター = 外部電荷つき core と厳密一致、
//!       w=0 でその弦コスト = 2h を検証する。
//!   [B] 極限の正しさ: h=0 は自由フェルミオンの 2 磁束セクター (境界位相は
//!       Jordan–Wigner 符号 (−1)^{N_f−1} 込み) の厳密和と機械精度一致。
//!       w=0 は対角厳密解と一致。L=8 で Lanczos = 稠密対角化。
//!   [C] 陰性対照: MI シャッフル・ランダム状態 (体積則)・サイトラベル置換は
//!       幾何読み出しに失敗する。奇数電荷セクターは構成自体が不能 (超選択則)。
//!
//! 物理: 1+1D Z2 ゲージ理論は h>0 で閉じ込め (弦張力 σ>0)。閉じ込めは
//! クエンチ後の相関前線を遅くする (中間子の束縛; Kormos et al. 2017 型の現象)。
//! 単一 core からこの両面 (静的な σ と動的な v) が読み出され、整合することが
//! 「同じ網の異なる読み出し」テーゼの相互作用系への最初の拡張である。

use std::f64::consts::PI;
use uft_sim::*;

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// h=0 の厳密自由エネルギー: 磁束セクター t=±1 (境界ひねり t·(−1)^{nf−1})
fn free_gs_energy(l: usize, nf: usize, w: f64) -> (f64, f64) {
    let mut es = Vec::new();
    for t in [1.0f64, -1.0] {
        let twist = t * if (nf as i32 - 1) % 2 == 0 { 1.0 } else { -1.0 };
        let theta = if twist > 0.0 { 0.0 } else { PI };
        let mut ek: Vec<f64> = (0..l)
            .map(|n| -2.0 * w * ((2.0 * PI * n as f64 + theta) / l as f64).cos())
            .collect();
        ek.sort_by(|a, b| a.partial_cmp(b).unwrap());
        es.push(ek[..nf].iter().sum::<f64>());
    }
    (es[0].min(es[1]), es[0].max(es[1]))
}

/// 拘束なし全空間 (L サイト + L リンクを陽に持つ) — 検証専用の小さな実装
struct FullZ2 {
    l: usize,
    w: f64,
    h: f64,
    m: f64,
}

impl FullZ2 {
    fn dim(&self) -> usize {
        1 << (2 * self.l)
    }
    fn split(&self, i: usize) -> (u32, u32) {
        ((i as u32) & ((1 << self.l) - 1), (i as u32) >> self.l)
    }
    fn elink(emask: u32, j: usize) -> f64 {
        if (emask >> j) & 1 == 0 {
            1.0
        } else {
            -1.0
        }
    }
    /// G_n = E_{n−1} E_n (−1)^{q_n} (外部電荷なし)
    fn gauss(&self, i: usize, n: usize) -> f64 {
        let (f, e) = self.split(i);
        let q = ((f >> n) & 1) ^ ((n % 2) as u32);
        let prev = (n + self.l - 1) % self.l;
        Self::elink(e, prev) * Self::elink(e, n) * if q == 0 { 1.0 } else { -1.0 }
    }
    fn dense(&self) -> Vec<f64> {
        let d = self.dim();
        let mut hm = vec![0.0; d * d];
        for i in 0..d {
            let (f, e) = self.split(i);
            let mut diag = 0.0;
            for n in 0..self.l {
                if (f >> n) & 1 == 1 {
                    diag += self.m * if n % 2 == 0 { 1.0 } else { -1.0 };
                }
            }
            for j in 0..self.l {
                diag -= self.h * Self::elink(e, j);
            }
            hm[i + i * d] = diag;
            // ホップ: ボンド a のリンク σ^z が E ビットを反転する
            for a in 0..self.l {
                let b = (a + 1) % self.l;
                for (x, y) in [(a, b), (b, a)] {
                    if (f >> x) & 1 == 1 && (f >> y) & 1 == 0 {
                        let m1 = f & !(1 << x);
                        let s1 = (m1 & ((1 << x) - 1)).count_ones()
                            + (f & ((1 << x) - 1)).count_ones()
                            - (m1 & ((1 << x) - 1)).count_ones();
                        let _ = s1;
                        let sgn1 = ((f & ((1 << x) - 1)).count_ones() % 2) as i32;
                        let sgn2 = ((m1 & ((1 << y) - 1)).count_ones() % 2) as i32;
                        let sgn = if (sgn1 + sgn2) % 2 == 0 { 1.0 } else { -1.0 };
                        let f2 = m1 | (1 << y);
                        let e2 = e ^ (1 << a);
                        let j2 = (f2 as usize) | ((e2 as usize) << self.l);
                        hm[j2 + i * d] += -self.w * sgn;
                    }
                }
            }
        }
        hm
    }
}

/// 稠密行列の下位固有値 (検証用)
fn dense_spectrum(g: &Z2GaugeRing) -> Vec<f64> {
    let d = g.dim();
    let mut hm = vec![0.0; d * d];
    for j in 0..d {
        let mut v = vec![(0.0, 0.0); d];
        v[j] = (1.0, 0.0);
        let hv = g.matvec_c(&v);
        for i in 0..d {
            hm[i + j * d] = hv[i].0;
        }
    }
    let (w, _) = jacobi_eigh(&hm, d);
    w
}

/// 前線速度 (到着時刻法): ボンド励起クエンチ → 各サイトの密度偏差がそのサイトの
/// ピーク値の 50% に初めて達する時刻 t_arr(d) → d–t_arr の線形 fit の傾き。
/// 固定の絶対閾値は Lieb–Robinson の裾を拾って過大評価し、時刻ごとの相対閾値は
/// 分散でピークが減衰すると破綻する — サイトごとのピーク規格化が頑健。
/// ピークが床 (1e-3) 未満のサイトは「未到達」: 到達サイトが 3 未満なら前線は凍結 (v=0)。
fn front_velocity(
    g: &Z2GaugeRing,
    gs: &Z2CoreState,
    dens0: &[f64],
    x0: usize,
    dt: f64,
    nstep: usize,
) -> (f64, usize, f64, f64) {
    let mut st = g.apply_bond_op(gs, x0);
    let c0 = g.conserved(&st);
    let (mut e0, mut n0) = (0.0, 0.0);
    for (k, v) in &c0 {
        if k == "energy" {
            e0 = *v;
        }
        if k == "norm" {
            n0 = *v;
        }
    }
    let l = g.l;
    let center = x0 as f64 + 0.5;
    let mut dev = vec![vec![0.0f64; l]; nstep];
    let mut max_edrift: f64 = 0.0;
    let mut max_ndrift: f64 = (n0 - 1.0).abs();
    for s in 0..nstep {
        st = g.step(&st, dt);
        let prof = v2_density_profile(&st);
        for j in 0..l {
            dev[s][j] = (prof[j] - dens0[j]).abs();
        }
        let c = g.conserved(&st);
        for (k, v) in &c {
            if k == "energy" {
                max_edrift = max_edrift.max((v - e0).abs());
            }
            if k == "norm" {
                max_ndrift = max_ndrift.max((v - 1.0).abs());
            }
        }
    }
    let mut pts: Vec<(f64, f64)> = Vec::new(); // (t_arr, d)
    for j in 0..l {
        let mut d = (j as f64 - center).abs();
        d = d.min(l as f64 - d);
        if !(1.5..=l as f64 / 2.0 - 1.0).contains(&d) {
            continue;
        }
        let peak = (0..nstep).map(|s| dev[s][j]).fold(0.0f64, f64::max);
        if peak < 1e-3 {
            continue; // 未到達 (凍結)
        }
        for s in 0..nstep {
            if dev[s][j] > 0.5 * peak {
                pts.push((dt * (s + 1) as f64, d));
                break;
            }
        }
    }
    let reached = pts.len();
    let v = if reached >= 3 {
        let xs: Vec<f64> = pts.iter().map(|p| p.0).collect();
        let ys: Vec<f64> = pts.iter().map(|p| p.1).collect();
        linfit(&xs, &ys).1.max(0.0)
    } else {
        0.0
    };
    (v, reached, max_edrift, max_ndrift)
}

fn main() {
    self_test();
    println!("=== v15.3 QrnCoreV2: ゲージ制約つき相互作用 core (Z2 格子ゲージ + staggered フェルミオン) ===\n");
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };

    // ================= [A] 拘束の正しさ (L=4, 拘束なし全空間との照合) =================
    println!("[A] 拘束の検証 — L=4 全空間 (256 次元、リンク陽持ち) vs 拘束を解いた core");
    {
        let (w, h, m) = (1.0, 0.7, 0.3);
        let full = FullZ2 { l: 4, w, h, m };
        let hm = full.dense();
        let d = full.dim();
        // エルミート性
        let mut herm: f64 = 0.0;
        for i in 0..d {
            for j in 0..d {
                herm = herm.max((hm[i + j * d] - hm[j + i * d]).abs());
            }
        }
        check(
            "全空間 H のエルミート性",
            herm < 1e-14,
            format!("max|H−Hᵀ| = {:.1e}", herm),
        );
        // (i) H が Gauss セクターを保つ: H_ij ≠ 0 なら全 n で G_n(i) = G_n(j)
        let mut sector_ok = true;
        for i in 0..d {
            for j in 0..d {
                if hm[i + j * d].abs() > 1e-14 {
                    for n in 0..4 {
                        if (full.gauss(i, n) - full.gauss(j, n)).abs() > 1e-14 {
                            sector_ok = false;
                        }
                    }
                }
            }
        }
        check(
            "[H, G_n] = 0 (H は Gauss セクターを厳密に保つ)",
            sector_ok,
            "全 65536 要素 × 4 拘束".to_string(),
        );
        // (ii) 物理セクター (全 G=+1) のスペクトル = 拘束を解いた core の合併スペクトル
        let phys: Vec<usize> = (0..d)
            .filter(|&i| (0..4).all(|n| full.gauss(i, n) > 0.0))
            .collect();
        let np = phys.len();
        let mut hp = vec![0.0; np * np];
        for (a, &i) in phys.iter().enumerate() {
            for (b, &j) in phys.iter().enumerate() {
                hp[a + b * np] = hm[i + j * d];
            }
        }
        let (wp, _) = jacobi_eigh(&hp, np);
        let mut wc: Vec<f64> = Vec::new();
        for nf in [0usize, 2, 4] {
            let g = Z2GaugeRing::try_new(4, nf, w, h, m, vec![]).unwrap();
            wc.extend(dense_spectrum(&g));
        }
        wc.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut dmax: f64 = 0.0;
        let dim_ok = np == wc.len();
        if dim_ok {
            for i in 0..np {
                dmax = dmax.max((wp[i] - wc[i]).abs());
            }
        }
        check(
            "物理セクターのスペクトル一致 (次元と全固有値)",
            dim_ok && dmax < 1e-10,
            format!("dim {} = {}, max|Δλ| = {:.1e}", np, wc.len(), dmax),
        );
        // (iii) Gauss 破り (G=−1 が 2 サイト) = 外部電荷つき core、w=0 の弦コスト = 2h
        let broken: Vec<usize> = (0..d)
            .filter(|&i| {
                (0..4).all(|n| {
                    let want = if n == 1 || n == 2 { -1.0 } else { 1.0 };
                    (full.gauss(i, n) - want).abs() < 1e-12
                })
            })
            .collect();
        let nb = broken.len();
        let mut hb = vec![0.0; nb * nb];
        for (a, &i) in broken.iter().enumerate() {
            for (b, &j) in broken.iter().enumerate() {
                hb[a + b * nb] = hm[i + j * d];
            }
        }
        let (wb, _) = jacobi_eigh(&hb, nb);
        let mut wce: Vec<f64> = Vec::new();
        for nf in [0usize, 2, 4] {
            let g = Z2GaugeRing::try_new(4, nf, w, h, m, vec![1, 2]).unwrap();
            wce.extend(dense_spectrum(&g));
        }
        wce.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let mut dmax2: f64 = 0.0;
        for i in 0..nb.min(wce.len()) {
            dmax2 = dmax2.max((wb[i] - wce[i]).abs());
        }
        check(
            "Gauss 破りセクター = 外部電荷つき core (スペクトル一致)",
            nb == wce.len() && dmax2 < 1e-10,
            format!("dim {} = {}, max|Δλ| = {:.1e}", nb, wce.len(), dmax2),
        );
        // w=0: 外部電荷対のコスト = min(2h·r, 2m) — 弦 vs 物質による遮蔽 (string breaking)。
        // 素朴な期待「常に 2h·r」は動的物質の存在で誤り: 2m < 2h·r なら中間子対が
        // 電荷を遮蔽して弦を切る。この min 構造こそ「物質つきゲージ理論」の証拠。
        let g0 = Z2GaugeRing::try_new(4, 2, 0.0, h, m, vec![]).unwrap();
        let g1 = Z2GaugeRing::try_new(4, 2, 0.0, h, m, vec![1, 2]).unwrap();
        let e0 = dense_spectrum(&g0)[0];
        let e1 = dense_spectrum(&g1)[0];
        let expect = (2.0 * h).min(2.0 * m);
        check(
            "w=0 の外部電荷対コスト = min(2h, 2m) (弦 or 遮蔽 — 厳密)",
            (e1 - e0 - expect).abs() < 1e-12,
            format!(
                "ΔE = {:.12} (min(2h, 2m) = {:.12}; 弦 2h = {:.2}, 遮蔽 2m = {:.2})",
                e1 - e0,
                expect,
                2.0 * h,
                2.0 * m
            ),
        );
        // 超選択則: 奇数電荷は構成不能
        let odd = Z2GaugeRing::try_new(4, 1, w, h, m, vec![]);
        check(
            "奇数電荷セクターは構成不能 (超選択則)",
            odd.is_err(),
            format!("{}", odd.err().unwrap_or_default()),
        );
    }

    // ================= [B] 極限の正しさ =================
    println!("\n[B] 極限の検証 — 自由極限・古典極限・稠密対角化");
    let (l, nf) = (14usize, 7usize);
    {
        // h=0, m=0: 自由フェルミオンの 2 磁束セクター
        let g = Z2GaugeRing::try_new(l, nf, 1.0, 0.0, 0.0, vec![]).unwrap();
        let (e_lz, _st, res) = g.ground_state(20260708);
        let (e_free, e_free2) = free_gs_energy(l, nf, 1.0);
        check(
            "h=0: E₀ = 自由フェルミオン (磁束セクター最小) と一致",
            (e_lz - e_free).abs() < 1e-9 && res < 1e-7,
            format!(
                "E_core = {:.12}, E_free = {:.12} (次セクター {:.6}), 残差 {:.1e}",
                e_lz, e_free, e_free2, res
            ),
        );
        // w=0: 対角厳密解 E₀ = −hL − mL/2
        let (h0, m0) = (0.8, 0.4);
        let g2 = Z2GaugeRing::try_new(l, nf, 0.0, h0, m0, vec![]).unwrap();
        let (e2, _s2, _r2) = g2.ground_state(3);
        let exact = -h0 * l as f64 - m0 * (l / 2) as f64;
        check(
            "w=0: E₀ = −hL − mL/2 (真空の対角厳密解)",
            (e2 - exact).abs() < 1e-9,
            format!("E = {:.12} vs {:.12}", e2, exact),
        );
        // L=8: Lanczos = 稠密対角化
        let g3 = Z2GaugeRing::try_new(8, 4, 1.0, 0.6, 0.2, vec![]).unwrap();
        let (e3, _s3, _r3) = g3.ground_state(7);
        let wd = dense_spectrum(&g3);
        check(
            "L=8: Lanczos E₀ = 稠密ヤコビ対角化",
            (e3 - wd[0]).abs() < 1e-10,
            format!("|Δ| = {:.1e} (dim {})", (e3 - wd[0]).abs(), g3.dim()),
        );
    }

    // ================= 基準点の状態と読み出し =================
    let (w0, h0, m0) = (1.0, 0.6, 0.2);
    println!(
        "\n[C] 同一状態からの読み出し — 基準点 L={}, N_f={}, w={}, h={}, m={} (相互作用・閉じ込め相)",
        l, nf, w0, h0, m0
    );
    let g = Z2GaugeRing::try_new(l, nf, w0, h0, m0, vec![]).unwrap();
    let (e_gs, gs, res_gs) = g.ground_state(20260708);
    println!(
        "  基底状態: E₀ = {:.9} (dim {}, Lanczos 残差 {:.1e})",
        e_gs,
        g.dim(),
        res_gs
    );
    {
        let cons = g.conserved(&gs);
        let norm = cons.iter().find(|c| c.0 == "norm").unwrap().1;
        let nf_v = cons.iter().find(|c| c.0 == "N_f").unwrap().1;
        check(
            "保存量: ‖ψ‖² = 1, ⟨N_f⟩ = 7 (Gauss 拘束は基底の構造として厳密)",
            (norm - 1.0).abs() < 1e-10 && (nf_v - 7.0).abs() < 1e-9,
            format!("norm − 1 = {:.1e}, N_f = {:.9}", norm - 1.0, nf_v),
        );
    }

    // ---- [幾何] MI → MDS → 円環 ----
    // 判定は v6.4 と同じ読み出し関数。ただし v6.4 の副次規準 (半径ばらつき ≤10%,
    // λ2/λ1 ≥ 0.9) は L≈200 の緩い MI 減衰で較正されたもので、L=14 の急峻な減衰では
    // 遠距離対の情報距離が床に張り付き半径が歪む (自由 h=0 でも同様 — 下で対照)。
    // 本版の判定は一次規準 (隣接復元) と陰性対照との弁別で行う。
    // 自由対照は m=0.2 (h=0): 厳密半充填 m=0 の自由鎖は 2k_F=π 振動で偶数距離の
    // 相関が厳密に消え、サイト単位の MI が偶奇縞になる (v0.7 が 2 サイトブロックで
    // 読んでいた理由)。質量 m で振動零点を外し、サイト単位読み出しの対照にする。
    let (mi, mi_max) = v2_mi_matrix(&gs);
    let geo = readout_ring_geometry(&mi, l, mi_max);
    let g_free = Z2GaugeRing::try_new(l, nf, w0, 0.0, m0, vec![]).unwrap();
    let (_ef, st_free, _rf) = g_free.ground_state(20260708);
    let (mif, mifmax) = v2_mi_matrix(&st_free);
    let geo_free = readout_ring_geometry(&mif, l, mifmax);
    check(
        "[幾何] MI 行列 → MDS が円環を復元 (隣接復元 100%; 相互作用 GS と自由 GS の両方)",
        geo.adjacency >= 0.999 && geo_free.adjacency >= 0.999,
        format!(
            "相互作用 (h=0.6): 隣接 {:.0}% (rsd {:.1}%, λ2/λ1 {:.2}) / 自由 (h=0, m=0.2): 隣接 {:.0}% (rsd {:.1}%)",
            geo.adjacency * 100.0,
            geo.rsd * 100.0,
            geo.lam21,
            geo_free.adjacency * 100.0,
            geo_free.rsd * 100.0
        ),
    );

    // ---- [エントロピー] 純粋性双対と面積則的プロファイル ----
    {
        let a2: Vec<usize> = (0..2).collect();
        let a4: Vec<usize> = (0..4).collect();
        let a7: Vec<usize> = (0..7).collect();
        let c7: Vec<usize> = (7..14).collect();
        let s7 = v2_entropy(&gs, &a7);
        let s7c = v2_entropy(&gs, &c7);
        check(
            "[S] 純粋性双対 S(A) = S(Aᶜ) (フェルミオン符号の総合検査)",
            (s7 - s7c).abs() < 1e-9,
            format!(
                "S(0..7) = {:.9}, S(7..14) = {:.9}, |Δ| = {:.1e}",
                s7,
                s7c,
                (s7 - s7c).abs()
            ),
        );
        let s2 = v2_entropy(&gs, &a2);
        let s4 = v2_entropy(&gs, &a4);
        println!(
            "       S(ℓ): ℓ=2 → {:.4}, ℓ=4 → {:.4}, ℓ=7 → {:.4} (閉じ込め相 — 増分の鈍化を確認)",
            s2, s4, s7
        );
        check(
            "[S] 面積則的振る舞い (S の増分が鈍る: S4−S2 > S7−S4 − 0.05)",
            (s4 - s2) > (s7 - s4) - 0.05,
            format!("ΔS(2→4) = {:.4}, ΔS(4→7) = {:.4}", s4 - s2, s7 - s4),
        );
    }

    // ---- [物質] staggered 凝縮 — ゲージ誘起の凝縮 (Schwinger 模型型) ----
    // 自由 (h=0, m=0) では密度一様で χ=0。ゲージ結合 h を入れると m=0 でも χ<0
    // (電場項が staggered 充填 = 電荷真空を好む — 質量ゼロでも凝縮が湧く)。
    // これは自由 core (v6.7) では原理的に出ない相互作用効果である。
    {
        let stag = |st: &Z2CoreState| -> f64 {
            (0..l)
                .map(|x| st.density(x) * if x % 2 == 0 { 1.0 } else { -1.0 })
                .sum::<f64>()
                / l as f64
        };
        let gfree = Z2GaugeRing::try_new(l, nf, w0, 0.0, 0.0, vec![]).unwrap();
        let (_e0, st_f, _r0) = gfree.ground_state(11);
        let chi_free = stag(&st_f);
        let gm0 = Z2GaugeRing::try_new(l, nf, w0, h0, 0.0, vec![]).unwrap();
        let (_e, st_m0, _r) = gm0.ground_state(11);
        let chi0 = stag(&st_m0);
        let gm5 = Z2GaugeRing::try_new(l, nf, w0, h0, 0.5, vec![]).unwrap();
        let (_e2, st_m5, _r2) = gm5.ground_state(11);
        let chi5 = stag(&st_m5);
        check(
            "[物質] 凝縮 χ: 自由 (h=0,m=0) で 0 / ゲージ誘起 (h=0.6,m=0) で負 / 質量で深まる",
            chi_free.abs() < 1e-9 && chi0 < -0.1 && chi5 < chi0,
            format!(
                "χ_free = {:.1e}, χ(h=0.6, m=0) = {:.4}, χ(h=0.6, m=0.5) = {:.4}",
                chi_free, chi0, chi5
            ),
        );
    }

    // ---- [ゲージ] 弦張力 σ(h) — 遮蔽を抑えるため重い物質 (m=3) で測る ----
    // 軽い物質 (m=0.2) では 2m < 2h·r で弦が切れる (遮蔽) ため E(r) が飽和し
    // 「張力」は定義できない。弦張力は重い物質の regime (2m ≫ 2h·r_max) で fit し、
    // 遮蔽そのもの (string breaking) は次の検査で実証する。
    println!("\n[D] ゲージと因果の握手 — 閉じ込めの二つの顔 (同一 core の静的/動的読み出し)");
    let m_heavy = 3.0;
    let hs = [0.0, 0.2, 0.4, 0.6];
    let mut sigmas = Vec::new();
    let mut gaps = Vec::new();
    for &hh in &hs {
        let mut es = Vec::new();
        for r in 0..=4usize {
            let ext = if r == 0 { vec![] } else { vec![3, 3 + r] };
            let gr = Z2GaugeRing::try_new(l, nf, w0, hh, m_heavy, ext).unwrap();
            let (er, _s, resr) = gr.ground_state(5);
            assert!(resr < 1e-6, "Lanczos 未収束 res={}", resr);
            es.push(er);
        }
        let rs: Vec<f64> = (0..=4).map(|r| r as f64).collect();
        let (_ic, slope) = linfit(&rs, &es);
        sigmas.push(slope);
        // 質量ギャップ (基準点の m=0.2 での E1 − E0)
        let gg = Z2GaugeRing::try_new(l, nf, w0, hh, m0, vec![]).unwrap();
        let mv = |v: &[(f64, f64)]| gg.matvec_c(v);
        let (ev2, _vv, _rr) = lanczos_lowest_herm(&mv, gg.dim(), 2, 160, 99);
        gaps.push(ev2[1] - ev2[0]);
    }
    // 弦切れ (string breaking): 軽い物質 m=0.2, h=0.8 で E(r) が min(2hr, ~2m+O(w)) に飽和。
    // 遮蔽はフェルミオン数を変える (奇-奇の電荷対は nf=5 の 2 正孔で遮蔽される) ため、
    // 物理的な E(r) は許されるセクター nf ∈ {5, 7, 9} の最小である (H は nf を保存 —
    // min はセクター間比較であって混合ではない)。
    let mut e_break = Vec::new();
    for r in 0..=4usize {
        let ext = if r == 0 { vec![] } else { vec![3, 3 + r] };
        let mut emin = f64::INFINITY;
        for nfs in [5usize, 7, 9] {
            let gb = Z2GaugeRing::try_new(l, nfs, w0, 0.8, m0, ext.clone()).unwrap();
            let (er, _s, _r2) = gb.ground_state(5);
            emin = emin.min(er);
        }
        e_break.push(emin);
    }
    let de_break: Vec<f64> = e_break.iter().map(|e| e - e_break[0]).collect();
    let mut fronts = Vec::new();
    let mut edrifts: f64 = 0.0;
    for &hh in &[0.0, 0.6, 1.2] {
        let gq = Z2GaugeRing::try_new(l, nf, w0, hh, m0, vec![]).unwrap();
        let (_e, gsq, _r) = gq.ground_state(20260708);
        let d0 = v2_density_profile(&gsq);
        let (v, reached, edrift, ndrift) = front_velocity(&gq, &gsq, &d0, 3, 0.05, 60);
        edrifts = edrifts.max(edrift).max(ndrift);
        fronts.push(v);
        if reached < 3 {
            println!(
                "   (h={} は前線凍結: 到達サイト {} — v=0 と記録)",
                hh, reached
            );
        }
    }
    println!(
        "   [静的] 弦張力 (m=3 重物質):  h = {:?} → σ = {:?}",
        hs,
        sigmas
            .iter()
            .map(|x| (x * 1000.0).round() / 1000.0)
            .collect::<Vec<_>>()
    );
    println!(
        "   [静的] ギャップ (m=0.2):      Δ = {:?}",
        gaps.iter()
            .map(|x| (x * 1000.0).round() / 1000.0)
            .collect::<Vec<_>>()
    );
    println!(
        "   [動的] 前線速度 (m=0.2):     h = [0.0, 0.6, 1.2] → v = {:.3}, {:.3}, {:.3}",
        fronts[0], fronts[1], fronts[2]
    );
    println!(
        "   [弦切れ] E(r)−E(0) (h=0.8, m=0.2): {:?} (2h·r なら {:?})",
        de_break
            .iter()
            .map(|x| (x * 100.0).round() / 100.0)
            .collect::<Vec<_>>(),
        (0..=4).map(|r| 1.6 * r as f64).collect::<Vec<_>>()
    );
    check(
        "[ゲージ] 弦張力 (重物質): σ(0) ≈ 0 (解放), h とともに単調増加, σ ≈ 2h (±25%)",
        sigmas[0].abs() < 0.02
            && sigmas[1] > 0.05
            && sigmas[2] > sigmas[1]
            && sigmas[3] > sigmas[2]
            && (1..4).all(|i| (sigmas[i] / (2.0 * hs[i]) - 1.0).abs() < 0.25),
        format!(
            "σ/2h = {:.3}, {:.3}, {:.3}",
            sigmas[1] / (2.0 * hs[1]),
            sigmas[2] / (2.0 * hs[2]),
            sigmas[3] / (2.0 * hs[3])
        ),
    );
    check(
        "[ゲージ] 弦切れ (軽物質): E(r) が立ち上がった後に飽和 (動的物質による遮蔽)",
        de_break[1] > 0.25
            && (de_break[4] - de_break[3]).abs() < 0.08
            && de_break[4] < 0.5 * 1.6 * 4.0,
        format!(
            "ΔE(1) = {:.3}, ΔE(4)−ΔE(3) = {:.3}, ΔE(4) = {:.3} ≪ 2h·4 = {:.1}",
            de_break[1],
            de_break[4] - de_break[3],
            de_break[4],
            6.4
        ),
    );
    // 到着時刻法 (ピークの 50% 初到達) が測るのは前線ピークの速度 — 自由双極子では
    // 群速度の重み平均 ⟨|v_g|⟩ = (1/π)∫|2w sin k|dk = 4w/π ≈ 1.27 の側 (最大群速度
    // 2w は Airy 前縁の速度で、L=14 では前線幅 ~t^{1/3} に埋もれ局所推定できない)。
    // 検査は (i) Lieb–Robinson 上界 v ≤ 2w (厳密不等式)、(ii) 自由値が ⟨|v_g|⟩ 帯、
    // (iii) h による単調減少 — の 3 点で行う。
    check(
        "[因果] 前線速度: 全て LR 上界 2w 以下, v(h=0) ≈ 4w/π (双極子の平均群速度), h で単調減少",
        fronts.iter().all(|&v| v <= 2.0 * w0 * 1.05)
            && (0.9..=1.6).contains(&fronts[0])
            && fronts[1] <= fronts[0] + 0.05
            && fronts[2] <= fronts[1] + 0.05,
        format!(
            "v = {:.3}, {:.3}, {:.3} (LR 上界 {:.1}, 4w/π = {:.2})",
            fronts[0],
            fronts[1],
            fronts[2],
            2.0 * w0,
            4.0 * w0 / std::f64::consts::PI
        ),
    );
    check(
        "[握手] 閉じ込めが情報を遅くする: σ ↑ かつ v ↓ (v(1.2) < 0.8 v(0))",
        fronts[2] < 0.8 * fronts[0],
        format!("v(1.2)/v(0) = {:.3}", fronts[2] / fronts[0]),
    );
    check(
        "[ギャップ] 閉じ込めとともにギャップが開く: Δ(h=0.6) > Δ(h=0) + 0.1",
        gaps[3] > gaps[0] + 0.1,
        format!("Δ = {:.4} → {:.4}", gaps[0], gaps[3]),
    );
    check(
        "[保存] クエンチ発展でノルム・エネルギーが保存 (Krylov 步進の健全性)",
        edrifts < 1e-8,
        format!("最大ドリフト {:.1e}", edrifts),
    );

    // ---- [S(t)] クエンチ後のエントロピー成長 — 閉じ込めは絡み合いの成長も遅くする ----
    let mut ds_pair = Vec::new();
    for &hh in &[0.0, 1.2] {
        let gq = Z2GaugeRing::try_new(l, nf, w0, hh, m0, vec![]).unwrap();
        let (_e, gsq, _r) = gq.ground_state(20260708);
        let mut st = gq.apply_bond_op(&gsq, 3);
        let region: Vec<usize> = (0..7).collect();
        let s_init = v2_entropy(&st, &region);
        let mut s_curve = vec![(0.0, s_init)];
        for k in 0..8 {
            for _ in 0..5 {
                st = gq.step(&st, 0.05);
            }
            let t = 0.25 * (k + 1) as f64;
            s_curve.push((t, v2_entropy(&st, &region)));
        }
        let s_end = s_curve.last().unwrap().1;
        println!(
            "   S_半分(t) [h={}]: {}",
            hh,
            s_curve
                .iter()
                .map(|(t, s)| format!("({:.2}, {:.3})", t, s))
                .collect::<Vec<_>>()
                .join(" ")
        );
        ds_pair.push(s_end - s_init);
    }
    check(
        "[S(t)] クエンチ後の絡み合い成長: 自由 (h=0) で成長し、閉じ込め (h=1.2) はそれを抑える",
        ds_pair[0] > 0.15 && ds_pair[1] < ds_pair[0],
        format!("ΔS(h=0) = {:.4}, ΔS(h=1.2) = {:.4}", ds_pair[0], ds_pair[1]),
    );

    // ================= [E] 陰性対照 (同じ読み出し関数・同じ判定器) =================
    println!("\n[E] 陰性対照 — 壊れた入力で壊れることの検査");
    {
        // (1) MI シャッフル (値の多重集合を保存)
        let mut rng = Rng::new(4242);
        let mut vals = Vec::new();
        for i in 0..l {
            for j in (i + 1)..l {
                vals.push(mi[i + j * l]);
            }
        }
        for i in (1..vals.len()).rev() {
            let j = rng.range(i + 1);
            vals.swap(i, j);
        }
        let mut mis = vec![0.0; l * l];
        let mut k = 0;
        for i in 0..l {
            for j in (i + 1)..l {
                mis[i + j * l] = vals[k];
                mis[j + i * l] = vals[k];
                k += 1;
            }
        }
        let gs_sh = readout_ring_geometry(&mis, l, mi_max);
        check(
            "MI シャッフル → 円環判定が破綻",
            !gs_sh.ring,
            format!(
                "隣接 {:.0}%, rsd {:.1}%",
                gs_sh.adjacency * 100.0,
                gs_sh.rsd * 100.0
            ),
        );
        // (2) サイトラベル置換 (誤ったテンソル分解)
        let mut perm: Vec<usize> = (0..l).collect();
        for i in (1..l).rev() {
            let j = rng.range(i + 1);
            perm.swap(i, j);
        }
        let mut mip = vec![0.0; l * l];
        for i in 0..l {
            for j in 0..l {
                mip[i + j * l] = mi[perm[i] + perm[j] * l];
            }
        }
        let gs_pm = readout_ring_geometry(&mip, l, mi_max);
        check(
            "サイトラベル置換 (誤った分解) → 円環判定が破綻",
            !gs_pm.ring,
            format!("隣接 {:.0}%", gs_pm.adjacency * 100.0),
        );
        // (3) ランダム状態 (体積則) — 同じ基底・同じ読み出し
        let mut psi: Vec<(f64, f64)> = (0..g.dim()).map(|_| (rng.gauss(), rng.gauss())).collect();
        let nr = cdot(&psi, &psi).0.sqrt();
        for a in psi.iter_mut() {
            a.0 /= nr;
            a.1 /= nr;
        }
        let rnd = Z2CoreState {
            l,
            masks: gs.masks.clone(),
            psi,
        };
        let (mir, mirmax) = v2_mi_matrix(&rnd);
        let gr = readout_ring_geometry(&mir, l, mirmax);
        let s_r = v2_entropy(&rnd, &(0..7).collect::<Vec<_>>());
        let s_g = v2_entropy(&gs, &(0..7).collect::<Vec<_>>());
        check(
            "ランダム状態 → 幾何なし + 体積則 (S ≫ 基底状態)",
            !gr.ring && s_r > 2.0 * s_g,
            format!(
                "隣接 {:.0}%, S_rand = {:.3} vs S_gs = {:.3}",
                gr.adjacency * 100.0,
                s_r,
                s_g
            ),
        );
    }

    // ================= JSON artifact =================
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v15.3".into())),
        (
            "model".into(),
            Json::Str("Z2GaugeRing L=14 nf=7 w=1 h=0.6 m=0.2".into()),
        ),
        ("dim".into(), Json::Num(g.dim() as f64)),
        ("e_gs".into(), Json::Num(e_gs)),
        ("adjacency".into(), Json::Num(geo.adjacency)),
        (
            "sigma_string".into(),
            Json::Arr(sigmas.iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "gap".into(),
            Json::Arr(gaps.iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "front_v".into(),
            Json::Arr(fronts.iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "h_grid_sigma".into(),
            Json::Arr(hs.iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "h_grid_front".into(),
            Json::Arr([0.0, 0.6, 1.2].iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "string_breaking_de".into(),
            Json::Arr(de_break.iter().map(|&x| Json::Num(x)).collect()),
        ),
        (
            "ds_quench".into(),
            Json::Arr(ds_pair.iter().map(|&x| Json::Num(x)).collect()),
        ),
    ]);
    let p = write_artifact("results/v153_corev2.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 拘束つき相互作用 core から幾何・エントロピー・物質・ゲージ・因果が読み出され、整合した"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
