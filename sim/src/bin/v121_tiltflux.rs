//! v12.1 傾き磁束 T⁴ — 指数 = Q₁Q₂ + ts の第一原理 3 世代 (構造編)
//!
//! これまでの積模型 (v7.2〜v11.3) は「T² のゼロモード Q 個 × Q 個から 3 世代を対で
//! 選ぶ」という**射影のアンザッツ**を必要とし、その対 σ は既約な自由度だった
//! (機構のはしご: S₃ 対 −19.86 に向き模型でも 1 nat 届かない)。
//!
//! 本バイナリは問いの土台を替える: **傾いた磁束の T⁴ では射影が要らない**。
//! 磁束 2-形式 F (平面 (x₁y₁)=Q₁, (x₂y₂)=Q₂, (x₁y₂)=t, (y₁x₂)=s) の Dirac 指数は
//! Pf(F) = Q₁Q₂ + t·s。(Q₁,Q₂,t,s) = (2,2,1,−1) なら指数 = 4 − 1 = **3** —
//! 3 世代が「9 (この磁束では 4) 個から選ぶ」のでなく**指数から直接**出る。
//!
//! 検証 (N=6 の T⁴ 格子, site 1296, 対角化 2592×2592 ×2):
//!  [1] 磁束の実装: 全 6 平面のプラケット位相が設計値 — 傾きも対照も機械精度
//!  [2] 指数 (傾きあり): 最低バンド縮退 = 3 (幅 ≪ ギャップ)
//!  [3] 指数 (対照 t=s=0): 縮退 = 4 (= Q₁Q₂) — 傾きが 4 → 3 に変えることの直接実証
//!  [4] 非因子化: バンド平均の同時分布 P(x₁,x₂) のトーラス間相互情報量が
//!      対照では 0 (機械精度)、傾きありでは正 — ゼロモードは 2 枚のトーラスに
//!      絡んで住む。「対・向き・整列」は傾き幾何の関数に置き換わる
//!
//! 将来課題 (正直に): 湯川の証拠比較には各セクターの Wilson 線ごとの対角化が要る
//! (傾きがあると並進は Wilson 線自体をシフトするため、単純な並進では代用できない —
//!  本バイナリの開発中にこの綱渡りを検出した)。36 配置 × 10 分 ≈ 6 時間の
//! バッチが次版の的である。

use uft_sim::*;

const N: usize = 6;
const NS4: usize = N * N * N * N;

type Flux = [i64; 4]; // (Q1, Q2, t, s)

fn idx(x1: usize, y1: usize, x2: usize, y2: usize) -> usize {
    x1 + N * (y1 + N * (x2 + N * y2))
}

/// リンク位相 θ_dir(site)。dir: 0=x1, 1=y1, 2=x2, 3=y2。境界で多価性を打ち消す。
fn link_phase(f: &Flux, x1: usize, y1: usize, x2: usize, y2: usize, dir: usize) -> f64 {
    let nn = (N * N) as f64;
    let two_pi = 2.0 * std::f64::consts::PI;
    let (p1, p2, pt, ps) = (
        two_pi * f[0] as f64 / nn,
        two_pi * f[1] as f64 / nn,
        two_pi * f[2] as f64 / nn,
        two_pi * f[3] as f64 / nn,
    );
    let nf = N as f64;
    match dir {
        0 => {
            if x1 == N - 1 {
                -(p1 * nf * y1 as f64 + pt * nf * y2 as f64)
            } else {
                0.0
            }
        }
        1 => {
            let mut th = p1 * x1 as f64;
            if y1 == N - 1 {
                th += -(ps * nf * x2 as f64);
            }
            th
        }
        2 => {
            let mut th = ps * y1 as f64;
            if x2 == N - 1 {
                th += -(p2 * nf * y2 as f64);
            }
            th
        }
        3 => p2 * x2 as f64 + pt * x1 as f64,
        _ => unreachable!(),
    }
}

fn plaquette(f: &Flux, x: [usize; 4], da: usize, db: usize) -> f64 {
    let step = |x: [usize; 4], d: usize| -> [usize; 4] {
        let mut y = x;
        y[d] = (y[d] + 1) % N;
        y
    };
    let lp = |x: [usize; 4], d: usize| link_phase(f, x[0], x[1], x[2], x[3], d);
    let xa = step(x, da);
    let xb = step(x, db);
    let th = lp(x, da) + lp(xa, db) - lp(xb, da) - lp(x, db);
    let two_pi = 2.0 * std::f64::consts::PI;
    let mut t = th % two_pi;
    if t > std::f64::consts::PI {
        t -= two_pi;
    }
    if t <= -std::f64::consts::PI {
        t += two_pi;
    }
    t
}

fn check_plaquettes(f: &Flux) -> f64 {
    let two_pi = 2.0 * std::f64::consts::PI;
    let nn = (N * N) as f64;
    let want = [
        (0usize, 1usize, two_pi * f[0] as f64 / nn),
        (2, 3, two_pi * f[1] as f64 / nn),
        (0, 3, two_pi * f[2] as f64 / nn),
        (1, 2, two_pi * f[3] as f64 / nn),
        (0, 2, 0.0),
        (1, 3, 0.0),
    ];
    let mut max_dev: f64 = 0.0;
    for x1 in 0..N {
        for y1 in 0..N {
            for x2 in 0..N {
                for y2 in 0..N {
                    for &(da, db, w) in &want {
                        let p = plaquette(f, [x1, y1, x2, y2], da, db);
                        max_dev = max_dev.max((p - w).abs());
                    }
                }
            }
        }
    }
    max_dev
}

/// 最低 nb 複素モードと実埋め込みスペクトルを返す (対角化 2·NS4)
fn t4_modes(f: &Flux, nb: usize) -> (Vec<Vec<(f64, f64)>>, Vec<f64>) {
    let m = 2 * NS4;
    let mut a = vec![0.0f64; m * m];
    let mut addhop = |i: usize, j: usize, th: f64| {
        let (c, s) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + NS4) + (i + NS4) * m] += -c;
        a[(i + NS4) + (j + NS4) * m] += -c;
        a[j + (i + NS4) * m] += s;
        a[(j + NS4) + i * m] += -s;
        a[i + (j + NS4) * m] += -s;
        a[(i + NS4) + j * m] += s;
    };
    for x1 in 0..N {
        for y1 in 0..N {
            for x2 in 0..N {
                for y2 in 0..N {
                    let i = idx(x1, y1, x2, y2);
                    addhop(
                        i,
                        idx((x1 + 1) % N, y1, x2, y2),
                        link_phase(f, x1, y1, x2, y2, 0),
                    );
                    addhop(
                        i,
                        idx(x1, (y1 + 1) % N, x2, y2),
                        link_phase(f, x1, y1, x2, y2, 1),
                    );
                    addhop(
                        i,
                        idx(x1, y1, (x2 + 1) % N, y2),
                        link_phase(f, x1, y1, x2, y2, 2),
                    );
                    addhop(
                        i,
                        idx(x1, y1, x2, (y2 + 1) % N),
                        link_phase(f, x1, y1, x2, y2, 3),
                    );
                }
            }
        }
    }
    let (w, v) = jacobi_eigh(&a, m);
    let mut modes: Vec<Vec<(f64, f64)>> = Vec::new();
    for kk in 0..2 * nb {
        let mut psi: Vec<(f64, f64)> = (0..NS4)
            .map(|i| (v[i + kk * m], v[(i + NS4) + kk * m]))
            .collect();
        for pm in &modes {
            let (mut pr, mut pi) = (0.0, 0.0);
            for i in 0..NS4 {
                pr += pm[i].0 * psi[i].0 + pm[i].1 * psi[i].1;
                pi += pm[i].0 * psi[i].1 - pm[i].1 * psi[i].0;
            }
            for i in 0..NS4 {
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
            if modes.len() == nb {
                break;
            }
        }
    }
    (modes, w)
}

/// バンド平均の同時分布 P(x₁,x₂) のトーラス間相互情報量 (基底回転に不変)
fn band_mi(modes: &[Vec<(f64, f64)>]) -> f64 {
    let mut p = vec![0.0f64; N * N];
    for m in modes {
        for x1 in 0..N {
            for y1 in 0..N {
                for x2 in 0..N {
                    for y2 in 0..N {
                        let (r, i) = m[idx(x1, y1, x2, y2)];
                        p[x1 + x2 * N] += r * r + i * i;
                    }
                }
            }
        }
    }
    let tot: f64 = p.iter().sum();
    for v in p.iter_mut() {
        *v /= tot;
    }
    let mut p1 = vec![0.0f64; N];
    let mut p2 = vec![0.0f64; N];
    for x1 in 0..N {
        for x2 in 0..N {
            p1[x1] += p[x1 + x2 * N];
            p2[x2] += p[x1 + x2 * N];
        }
    }
    let mut mi = 0.0;
    for x1 in 0..N {
        for x2 in 0..N {
            let q = p[x1 + x2 * N];
            if q > 1e-300 {
                mi += q * (q / (p1[x1] * p2[x2])).ln();
            }
        }
    }
    mi
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

fn main() {
    self_test();
    println!("=== v12.1 傾き磁束 T⁴: 指数 = Q₁Q₂ + ts の第一原理 3 世代 (構造編) ===\n");
    let f_tilt: Flux = [2, 2, 1, -1];
    let f_ctrl: Flux = [2, 2, 0, 0];
    println!(
        "格子 N={} (site {}), 傾き (2,2,1,−1) vs 対照 (2,2,0,0)",
        N, NS4
    );

    // ---- [1] 磁束の実装検査 ----
    println!("\n[1] プラケット検査 (全 site × 全 6 平面)");
    let d_tilt = check_plaquettes(&f_tilt);
    let d_ctrl = check_plaquettes(&f_ctrl);
    let ok_flux = d_tilt < 1e-12 && d_ctrl < 1e-12;
    println!(
        "    max|Δ位相|: 傾き {:.2e} / 対照 {:.2e} (< 1e-12)  {}",
        d_tilt,
        d_ctrl,
        pass(ok_flux)
    );

    // ---- [2] 指数 (傾き): 縮退 3 ----
    println!("\n[2] 指数 (傾き 2,2,1,−1): 最低バンド縮退 = Pf(F) = 3 か (対角化 ~10 分)");
    let t0 = std::time::Instant::now();
    let (modes_t, w_t) = t4_modes(&f_tilt, 3);
    let spread_t = w_t[5] - w_t[0];
    let gap_t = w_t[6] - w_t[5];
    let ok_idx3 = spread_t < 1e-8 && gap_t > 0.02;
    println!(
        "    幅 {:.2e} ≪ ギャップ {:.4} → 縮退 3  {}  ({} ms)",
        spread_t,
        gap_t,
        pass(ok_idx3),
        t0.elapsed().as_millis()
    );

    // ---- [3] 指数 (対照): 縮退 4 ----
    println!("\n[3] 指数 (対照 2,2,0,0): 最低バンド縮退 = Q₁Q₂ = 4 か (対角化 ~10 分)");
    let t1 = std::time::Instant::now();
    let (modes_c, w_c) = t4_modes(&f_ctrl, 4);
    let spread_c = w_c[7] - w_c[0];
    let gap_c = w_c[8] - w_c[7];
    let ok_idx4 = spread_c < 1e-8 && gap_c > 0.02;
    println!(
        "    幅 {:.2e} ≪ ギャップ {:.4} → 縮退 4  {}  ({} ms)",
        spread_c,
        gap_c,
        pass(ok_idx4),
        t1.elapsed().as_millis()
    );
    println!("    => 傾き (t,s)=(1,−1) が指数を 4 → 3 に変えた: Pf(F) = Q₁Q₂ + ts の直接実証");

    // ---- [4] 非因子化: トーラス間相互情報量 ----
    println!("\n[4] ゼロモードの非因子化 (バンド平均 P(x₁,x₂) のトーラス間 MI)");
    let mi_t = band_mi(&modes_t);
    let mi_c = band_mi(&modes_c);
    let ok_mi = mi_c < 1e-10 && mi_t > 100.0 * mi_c.max(1e-12) && mi_t > 1e-3;
    println!(
        "    MI(傾き) = {:.4} nats / MI(対照) = {:.2e} (対照は積 = 0)  {}",
        mi_t,
        mi_c,
        pass(ok_mi)
    );
    println!("    => 傾きのゼロモードは 2 枚のトーラスに絡んで住む — 対のアンザッツが");
    println!("       果たしていた役割 (どのモードとどのモードが同じ 4D 場か) は、");
    println!("       傾き幾何そのものが果たす。");

    let all_ok = ok_flux && ok_idx3 && ok_idx4 && ok_mi;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-012".into())),
        ("lattice".into(), Json::Int(N as i64)),
        ("plaquette_dev".into(), Json::Num(d_tilt.max(d_ctrl))),
        ("tilt_spread".into(), Json::Num(spread_t)),
        ("tilt_gap".into(), Json::Num(gap_t)),
        ("ctrl_spread".into(), Json::Num(spread_c)),
        ("ctrl_gap".into(), Json::Num(gap_c)),
        ("index_tilt".into(), Json::Int(3)),
        ("index_ctrl".into(), Json::Int(4)),
        ("band_mi_tilt".into(), Json::Num(mi_t)),
        ("band_mi_ctrl".into(), Json::Num(mi_c)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v121_tiltflux.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {}", pass(all_ok));
    println!("\n結論: 傾き磁束の T⁴ は、射影のアンザッツなしに指数 = Q₁Q₂ + ts = 3 の");
    println!("      3 世代を実現し、そのゼロモードは非因子化 (トーラス間 MI > 0)。");
    println!("      「対 σ とは何か」への構成的な答えの候補が立った。湯川の証拠比較は");
    println!("      Wilson 線ごとの対角化 (36 × ~10 分) を要する次版の的である。");
    if !all_ok {
        std::process::exit(1);
    }
}
