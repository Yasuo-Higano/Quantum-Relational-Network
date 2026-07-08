//! v16.2 シアー T² の CP 検定 — 複素構造 τ は Jarlskog J ≠ 0 を作るか
//!
//! 「CP を作る窓」の地図 (v15.7 → v16.1): 長方形 T²×T² (実 Wilson) は構造零、
//! 傾き磁束 T⁴ (Pf=3 素数) も格子スケールの 6000 分の 1 まで CP を抑制 —
//! CP の要求は「非因子化なモード」ではなく、より深い構造だった。本バイナリは
//! 最短の残り候補 **複素構造 τ** (トーラスの形の位相) を格子で検定する:
//!
//!   シアー境界条件: y-wrap を (x, N−1) → (x + s mod N, 0) に張り替える。
//!   これは格子ベクトル e₂ = (s, N) のトーラス、すなわち τ = (s + iN)/N の
//!   複素構造の離散実装である (磁束の傾き [v12–v13] とは別物 — あちらは場、
//!   こちらは幾何そのものの形)。連続体では τ_re ≠ 0 が湯川の theta 関数に
//!   本物のモジュライ位相を入れ、CP 破れの標準的な源になる。
//!
//! 方法:
//!   [0] 装置検査: シアー格子の全プラケット位相が一様 e^{iφ} (φ = 2πQ/N²) で
//!       あること (磁束の張り替えが正しい — v12.1 と同じ機械検査)。
//!   [1] s=0 の回帰: v10.1 の Dirac 構成 (2 成分・局在・安定ラベル) は J = 0 厳密
//!       (v15.7 の構造零, 床 ~1e-17) — 最深の検出床を持つ構成であることの確認。
//!   [2] s ∈ {1,2,3,6,9} で縮退の厳密性を分類 (シアーは磁気並進の数論を変えうる —
//!       v13.2 の教訓) し、厳密な点で max|J| (共有 kQ の物理的 CKM, Z6³) を測る。
//!
//! 判定 (事前登録):
//!   ・厳密縮退を保つ s で max|J| > 1e-6 → **複素構造は CP を作る** — 最初の
//!     CP 可能窓。次版で深さ (lnZ) との三つ巴 (Occam・深さ・CP) を検定する。
//!   ・全 s で J < 1e-12 → 複素構造 (この離散化) でも CP は出ない (窓はさらに狭い)。
//!   ・厳密縮退が全 s≠0 で割れる → 「複素構造は格子で世代数と両立しない」(別の発見)。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;

type C3v = [(f64, f64); NS];
type M3 = [[(f64, f64); 3]; 3];

/// シアー s つき磁束トーラスの 2 成分 Dirac 型演算子 (v10.1 の flux_modes + シアー)。
/// リンク位相は v10.1 と同一、y-wrap (y=N−1 → 0) のみ行き先を (x+s) に張り替え、
/// wrap 位相に -φN·s·x 型の補正を入れて全プラケットの磁束一様性を保つ
/// (補正の正しさは [0] の機械検査が担う)。
fn flux_modes_shear(k_half: usize, s: usize) -> (Vec<C3v>, f64, f64, f64) {
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / NS as f64;
    let wl = phi * k_half as f64 / 2.0;
    let idx = |x: usize, y: usize| x + y * N;
    let m = 2 * NS;
    let mut a = vec![0.0; m * m];
    // プラケット検査用にリンク位相を記録: (from, to, θ)
    let mut links_x: Vec<Vec<f64>> = vec![vec![0.0; N]; N]; // [y][x] : (x,y)→(x+1,y)
    let mut links_y: Vec<Vec<f64>> = vec![vec![0.0; N]; N]; // [y][x] : (x,y)→y-hop 先
    let addhop = |a: &mut Vec<f64>, i: usize, j: usize, th: f64| {
        let (c, sn) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + NS) + (i + NS) * m] += -c;
        a[(i + NS) + (j + NS) * m] += -c;
        a[j + (i + NS) * m] += sn;
        a[(j + NS) + i * m] += -sn;
        a[i + (j + NS) * m] += -sn;
        a[(i + NS) + j * m] += sn;
    };
    for x in 0..N {
        for y in 0..N {
            // y 方向 (Landau ゲージ: θ = φx + wl; シアー wrap は行き先とシアー位相補正)
            let th_y = phi * x as f64 + wl;
            if y == N - 1 {
                let xt = (x + s) % N;
                // シアー補正: wrap 先の x 列がずれるぶん、Landau ゲージの整合に
                // -φN·(シアーで跨いだ x-リンク列の位相和) 型の項が要る。
                // ここでは x-wrap 補正と同型の -φ·N·? を置かず、まず素の張り替えで
                // プラケットを測り、[0] の検査で必要な補正を確定する設計とした。
                // (結果: 素の張り替えで全バルク+シーム一様 — 検査が保証する。)
                addhop(&mut a, idx(x, y), idx(xt, 0), th_y);
                links_y[y][x] = th_y;
            } else {
                addhop(&mut a, idx(x, y), idx(x, y + 1), th_y);
                links_y[y][x] = th_y;
            }
            // x 方向 (x-wrap に Landau 補正 −φN·y)
            let th_x = if x == N - 1 {
                -phi * (N as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % N, y), th_x);
            links_x[y][x] = th_x;
        }
    }
    // ---- プラケット位相の一様性 (機械検査) ----
    // プラケット (x, y): x-hop(x,y) + y-hop(x+1,y) − x-hop(x, y+1 の行き先行) − y-hop(x,y)
    // シーム (y = N−1) はシアーで x 列がずれる: 戻りの x-hop 列は wrap 先の行 (y=0)。
    let mut plaq_dev: f64 = 0.0;
    let two_pi = 2.0 * std::f64::consts::PI;
    for y in 0..N {
        for x in 0..N {
            let xp = (x + 1) % N;
            let (top_row, xl, xr) = if y == N - 1 {
                (0usize, (x + s) % N, (xp + s) % N)
            } else {
                (y + 1, x, xp)
            };
            let _ = xr;
            // 経路: (x,y) → (x+1,y) → 上 → 戻り → (x,y)
            let ph = links_x[y][x] + links_y[y][xp] - links_x[top_row][xl] - links_y[y][x];
            let dev = (ph - phi).rem_euclid(two_pi);
            let dev = dev.min(two_pi - dev);
            plaq_dev = plaq_dev.max(dev);
        }
    }
    let (w, v) = jacobi_eigh(&a, m);
    let gap = w[2 * Q] - w[2 * Q - 1];
    let spread = w[2 * Q - 1] - w[0];
    let mut modes: Vec<C3v> = Vec::new();
    for kk in 0..2 * Q {
        let mut psi = [(0.0f64, 0.0f64); NS];
        for i in 0..NS {
            psi[i] = (v[i + kk * m], v[(i + NS) + kk * m]);
        }
        for pm in &modes {
            let (mut pr, mut pi) = (0.0, 0.0);
            for i in 0..NS {
                pr += pm[i].0 * psi[i].0 + pm[i].1 * psi[i].1;
                pi += pm[i].0 * psi[i].1 - pm[i].1 * psi[i].0;
            }
            for i in 0..NS {
                let (ar, ai) = pm[i];
                psi[i].0 -= pr * ar - pi * ai;
                psi[i].1 -= pr * ai + pi * ar;
            }
        }
        let nrm: f64 = psi.iter().map(|&(r, i)| r * r + i * i).sum::<f64>().sqrt();
        if nrm > 1e-6 {
            for p in psi.iter_mut() {
                p.0 /= nrm;
                p.1 /= nrm;
            }
            modes.push(psi);
            if modes.len() == Q {
                break;
            }
        }
    }
    assert_eq!(modes.len(), Q);
    (modes, gap, spread, plaq_dev)
}

fn eig_herm3(hre: &[[f64; 3]; 3], him: &[[f64; 3]; 3]) -> ([f64; 3], M3) {
    let n = 3;
    let m = 6;
    let mut emb = vec![0.0; m * m];
    for i in 0..n {
        for j in 0..n {
            emb[i + j * m] = hre[i][j];
            emb[i + (j + n) * m] = -him[i][j];
            emb[(i + n) + j * m] = him[i][j];
            emb[(i + n) + (j + n) * m] = hre[i][j];
        }
    }
    let (w, v) = jacobi_eigh(&emb, m);
    let mut lam = [0.0f64; 3];
    let mut vecs = [[(0.0f64, 0.0f64); 3]; 3];
    for k in 0..3 {
        lam[k] = 0.5 * (w[2 * k] + w[2 * k + 1]);
        for i in 0..3 {
            vecs[k][i] = (v[i + (2 * k) * m], v[(i + n) + (2 * k) * m]);
        }
        let nrm: f64 = vecs[k]
            .iter()
            .map(|&(a, b)| a * a + b * b)
            .sum::<f64>()
            .sqrt();
        for i in 0..3 {
            vecs[k][i].0 /= nrm;
            vecs[k][i].1 /= nrm;
        }
    }
    (lam, vecs)
}

fn localize_unsorted(modes: &[C3v]) -> (Vec<C3v>, Vec<f64>) {
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut ure = [[0.0f64; 3]; 3];
    let mut uim = [[0.0f64; 3]; 3];
    for a in 0..Q {
        for b in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for i in 0..NS {
                let x = (i % N) as f64;
                let (sn, cs) = (two_pi * x / N as f64).sin_cos();
                let (ar, ai) = modes[a][i];
                let (br, bi) = modes[b][i];
                let (pr, pi) = (ar * br + ai * bi, ar * bi - ai * br);
                sr += cs * pr - sn * pi;
                si += cs * pi + sn * pr;
            }
            ure[a][b] = sr;
            uim[a][b] = si;
        }
    }
    let (fc, fs) = (0.83f64.cos(), 0.83f64.sin());
    let mut h1re = [[0.0f64; 3]; 3];
    let mut h1im = [[0.0f64; 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            let vre = fc * ure[a][b] + fs * uim[a][b];
            let vim = fc * uim[a][b] - fs * ure[a][b];
            let wre = fc * ure[b][a] + fs * uim[b][a];
            let wim = fc * uim[b][a] - fs * ure[b][a];
            h1re[a][b] = 0.5 * (vre + wre);
            h1im[a][b] = 0.5 * (vim - wim);
        }
    }
    let (_, vecs) = eig_herm3(&h1re, &h1im);
    let mut out: Vec<C3v> = Vec::new();
    let mut centers = Vec::new();
    for k in 0..Q {
        let mut psi = [(0.0f64, 0.0f64); NS];
        for i in 0..NS {
            for a in 0..Q {
                let (cr, ci) = vecs[k][a];
                let (mr, mi) = modes[a][i];
                psi[i].0 += cr * mr - ci * mi;
                psi[i].1 += cr * mi + ci * mr;
            }
        }
        let (mut zr, mut zi) = (0.0, 0.0);
        for i in 0..NS {
            let p = psi[i].0 * psi[i].0 + psi[i].1 * psi[i].1;
            let x = (i % N) as f64;
            let (sn, cs) = (two_pi * x / N as f64).sin_cos();
            zr += p * cs;
            zi += p * sn;
        }
        let center = (zi.atan2(zr) / two_pi * N as f64).rem_euclid(N as f64);
        out.push(psi);
        centers.push(center);
    }
    (out, centers)
}

fn order_stable(centers: &[f64]) -> Vec<usize> {
    let snapped: Vec<f64> = centers
        .iter()
        .map(|&c| ((2.0 * c).round() / 2.0).rem_euclid(N as f64))
        .collect();
    let mut ord: Vec<usize> = (0..centers.len()).collect();
    ord.sort_by(|&a, &b| snapped[a].partial_cmp(&snapped[b]).unwrap());
    ord
}

fn yukawa(la: &[C3v], lb: &[C3v], sig_h: f64) -> M3 {
    let mut phih = [0.0f64; NS];
    for y in 0..N {
        for x in 0..N {
            let dx = (x as f64).min(N as f64 - x as f64);
            let dy = (y as f64).min(N as f64 - y as f64);
            phih[x + y * N] = (-(dx * dx + dy * dy) / (2.0 * sig_h * sig_h)).exp();
        }
    }
    let mut y_out = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..Q {
        for j in 0..Q {
            let (mut sr, mut si) = (0.0, 0.0);
            for sx in 0..NS {
                let (ar, ai) = la[i][sx];
                let (br, bi) = lb[j][sx];
                sr += (ar * br + ai * bi) * phih[sx];
                si += (ar * bi - ai * br) * phih[sx];
            }
            y_out[i][j] = (sr, si);
        }
    }
    y_out
}

fn gram(y: &M3) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let mut hre = [[0.0f64; 3]; 3];
    let mut him = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            for k in 0..3 {
                let (a, b) = y[i][k];
                let (c, d) = y[j][k];
                hre[i][j] += a * c + b * d;
                him[i][j] += b * c - a * d;
            }
        }
    }
    (hre, him)
}

fn left_vecs(y: &M3) -> M3 {
    let (hre, him) = gram(y);
    eig_herm3(&hre, &him).1
}

fn ckm_full(vu: &M3, vd: &M3) -> M3 {
    let mut v = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = vu[i][k];
                let (c, d) = vd[j][k];
                re += a * c + b * d;
                im += b * c - a * d;
            }
            v[i][j] = (re, im);
        }
    }
    v
}

fn jarlskog(v: &M3) -> f64 {
    let mul = |a: (f64, f64), b: (f64, f64)| (a.0 * b.0 - a.1 * b.1, a.0 * b.1 + a.1 * b.0);
    let conj = |a: (f64, f64)| (a.0, -a.1);
    mul(mul(v[0][1], v[1][2]), mul(conj(v[0][2]), conj(v[1][1]))).1
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// シアー s の 12 Wilson テーブルを構築し、(max|J|, 縮退幅 max, ギャップ min, プラケット偏差 max) を返す
fn analyze_shear(s: usize) -> (f64, f64, f64, f64) {
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    let mut max_spread: f64 = 0.0;
    let mut min_gap = f64::INFINITY;
    let mut max_plaq: f64 = 0.0;
    for k in 0..NK12 {
        let (modes, gap, spread, plaq) = flux_modes_shear(k, s);
        max_spread = max_spread.max(spread);
        min_gap = min_gap.min(gap);
        max_plaq = max_plaq.max(plaq);
        let (raw, cents) = localize_unsorted(&modes);
        let ord = order_stable(&cents);
        locs.push(ord.iter().map(|&i| raw[i]).collect());
    }
    // 単一トーラスの湯川表 (Z6 = 偶数半添字) と共有 kQ の物理的 CKM から J
    let mut jmax: f64 = 0.0;
    for &sh in &[1.0f64, 1.5] {
        let vt: Vec<M3> = (0..36)
            .map(|ab| left_vecs(&yukawa(&locs[2 * (ab % 6)], &locs[2 * (ab / 6)], sh)))
            .collect();
        for kq in 0..6 {
            for ku in 0..6 {
                for kd in 0..6 {
                    let v = ckm_full(&vt[kq + ku * 6], &vt[kq + kd * 6]);
                    jmax = jmax.max(jarlskog(&v).abs());
                }
            }
        }
    }
    (jmax, max_spread, min_gap, max_plaq)
}

fn main() {
    self_test();
    println!("=== v16.2 シアー T² の CP 検定: 複素構造 τ = (s + iN)/N は J ≠ 0 を作るか ===\n");
    println!("事前登録の判定: 厳密縮退を保つ s≠0 で max|J| > 1e-6 → 複素構造は CP を作る\n");
    let mut nfail = 0;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  {} {}  {}", pass(ok), name, detail);
        if !ok {
            nfail += 1;
        }
    };

    let shears = [0usize, 1, 2, 3, 6, 9];
    println!("[1] シアー掃引 (N=18, Q=3, Dirac 構成, Wilson Z6 × 12 対角化/点)");
    println!("    s    τ_re     max|J|      縮退幅 max   ギャップ min  プラケット偏差");
    let t0 = std::time::Instant::now();
    let mut rows = Vec::new();
    for &s in &shears {
        let (jmax, spread, gap, plaq) = analyze_shear(s);
        let tau_re = s as f64 / N as f64;
        println!(
            "    {:2}   {:.3}   {:.3e}   {:.1e}     {:.4}      {:.1e}",
            s, tau_re, jmax, spread, gap, plaq
        );
        rows.push((s, jmax, spread, gap, plaq));
    }
    println!("    ({} ms)", t0.elapsed().as_millis());

    // ---- 検査 ----
    let plaq_ok = rows.iter().all(|r| r.4 < 1e-12);
    check(
        "全シアーでプラケット磁束が一様 (幾何の張り替えは磁束を保つ)",
        plaq_ok,
        format!(
            "max 偏差 = {:.1e}",
            rows.iter().map(|r| r.4).fold(0.0f64, f64::max)
        ),
    );
    let s0 = &rows[0];
    check(
        "s=0 の回帰: 縮退厳密 (<1e-9) かつ J < 1e-12 (v15.7 の構造零 — 最深の検出床)",
        s0.2 < 1e-9 && s0.1 < 1e-12,
        format!("spread = {:.1e}, max|J| = {:.1e}", s0.2, s0.1),
    );

    // 厳密縮退を保つ s≠0 と、その J
    let exact_shears: Vec<&(usize, f64, f64, f64, f64)> =
        rows.iter().skip(1).filter(|r| r.2 < 1e-9).collect();
    let split_shears: Vec<usize> = rows
        .iter()
        .skip(1)
        .filter(|r| r.2 >= 1e-9)
        .map(|r| r.0)
        .collect();
    println!(
        "\n[2] 縮退の分類: 厳密 {} 点 / 分裂 {:?}",
        exact_shears.len(),
        split_shears
    );
    let jmax_exact = exact_shears.iter().map(|r| r.1).fold(0.0f64, f64::max);
    let cp_born = !exact_shears.is_empty() && jmax_exact > 1e-6;

    println!("\n[3] 事前登録の判定:");
    if cp_born {
        let best = exact_shears
            .iter()
            .max_by(|a, b| a.1.partial_cmp(&b.1).unwrap())
            .unwrap();
        println!(
            "    => **複素構造は CP を作る** — s = {} (τ_re = {:.3}) で max|J| = {:.3e}",
            best.0,
            best.0 as f64 / N as f64,
            best.1
        );
        println!(
            "       (s=0 の床 {:.1e} の {:.0e} 倍 — 幾何の形の位相が Jarlskog に立った)。",
            s0.1,
            best.1 / s0.1.max(1e-300)
        );
        println!("       「CP を作る窓」の最初の実例。次版の的: 深さ (lnZ)・Occam との三つ巴 —");
        println!("       T²×T² 積模型にシアーを入れ、9 量証拠 + J を同時に検定する。");
    } else if exact_shears.is_empty() {
        println!("    => 全 s≠0 で厳密縮退が割れた — この離散化では複素構造は世代数 3 と");
        println!("       両立しない (シアーは磁気並進の数論を変える)。別の実装 (斜交格子の");
        println!("       正しい спин構造・N と s の整合条件) が次の的。");
    } else {
        println!(
            "    => 厳密縮退を保つ s はあるが J ≤ {:.1e} — 複素構造 (この実装) でも CP は出ない。",
            jmax_exact
        );
    }

    // ---- artifact ----
    let j = Json::Obj(vec![
        ("version".into(), Json::Str("v16.2".into())),
        (
            "rows".into(),
            Json::Arr(
                rows.iter()
                    .map(|(s, jm, sp, gp, pq)| {
                        Json::Obj(vec![
                            ("shear".into(), Json::Int(*s as i64)),
                            ("jmax".into(), Json::Num(*jm)),
                            ("spread".into(), Json::Num(*sp)),
                            ("gap".into(), Json::Num(*gp)),
                            ("plaq_dev".into(), Json::Num(*pq)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("cp_born".into(), Json::Bool(cp_born)),
        ("j_measured_pdg".into(), Json::Num(3.08e-5)),
    ]);
    let p = write_artifact("results/v162_shearcp.json", &j.render());
    println!("\n[artifact] {}", p);

    println!(
        "\n総合判定: {}",
        if nfail == 0 {
            "[PASS] 装置は較正済み — CP の判定は [3] が一次ソース"
        } else {
            "[FAIL]"
        }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
