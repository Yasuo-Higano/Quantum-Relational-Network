//! v6.4 幾何読み出しの陰性対照 — 「何でも円環に見える装置」ではないことの実証
//!
//! v0.7 は円環格子の基底状態の相互情報量 (MI) だけから MDS で円環を 100% 復元した
//! (claims.yml: QRN-GEOM-003, C3)。だが「幾何を持たない状態を入れても幾何が出てしまう」
//! なら、この読み出しは何も言っていない。本バイナリは同一パイプライン
//! (ブロック MI → 情報距離 −ln(MI/MI_max) → 古典的 MDS → 円環判定) に
//!   [P]  正例: 円環格子の基底状態 (v0.7 と同じ)
//!   [N1] MI 行列の対称シャッフル (値の多重集合は同じ・局所構造だけ破壊)
//!   [N2] 非局所ハミルトニアン (GOE ランダム対称行列) の基底状態
//!   [N3] volume-law ランダム状態 (ランダム基底のランダム半分充填射影)
//!   [N4] 古典的・無相関状態 (ランダム占有の積状態 — MI ≡ 0)
//! を通し、正例だけが円環判定を通ることを検査する。
//! 円環判定: 隣接復元率 ≥ 90% かつ 半径ばらつき ≤ 10% かつ MDS 上位 2 固有値の
//! 縮退度 λ2/λ1 ≥ 0.9 (円 = 2 次元等方埋め込み)。正例はさらに MI 減衰べき −2±0.3。
//! 乱数は全て固定シード。

use uft_sim::*;

fn h2(z: f64) -> f64 {
    let z = z.clamp(1e-14, 1.0 - 1e-14);
    -z * z.ln() - (1.0 - z) * (1.0 - z).ln()
}
fn entropy_real(c: &[f64], n: usize) -> f64 {
    let (w, _) = jacobi_eigh(c, n);
    w.iter().map(|&z| h2(z)).sum()
}

const N: usize = 202;
const NB: usize = N / 2;

/// 2 サイトブロック間の MI 行列 (NB×NB) を一般の実相関行列から計算
fn block_mi(c: &[f64], mi: &mut [f64]) -> f64 {
    let mut sblk = vec![0.0; NB];
    let mut c2 = [0.0; 4];
    for b in 0..NB {
        let s = [2 * b, 2 * b + 1];
        for i in 0..2 {
            for j in 0..2 {
                c2[i + j * 2] = c[s[i] + s[j] * N];
            }
        }
        sblk[b] = entropy_real(&c2, 2);
    }
    let mut c4 = [0.0; 16];
    let mut mi_max = 0.0f64;
    for i in 0..NB {
        for j in (i + 1)..NB {
            let s = [2 * i, 2 * i + 1, 2 * j, 2 * j + 1];
            for a in 0..4 {
                for b in 0..4 {
                    c4[a + b * 4] = c[s[a] + s[b] * N];
                }
            }
            let m = (sblk[i] + sblk[j] - entropy_real(&c4, 4)).max(0.0);
            mi[i + j * NB] = m;
            mi[j + i * NB] = m;
            mi_max = mi_max.max(m);
        }
    }
    mi_max
}

struct Metrics {
    adjacency: f64,
    rsd: f64,
    lam21: f64,
    mi_max: f64,
    ring: bool,
}

/// MI 行列 → 情報距離 → MDS → 円環判定
fn ring_metrics(mi: &[f64], mi_max: f64) -> Metrics {
    if mi_max < 1e-12 {
        // 相関が皆無 — 幾何は定義できない (それ自体が正しい検出)
        return Metrics {
            adjacency: 0.0,
            rsd: f64::INFINITY,
            lam21: 0.0,
            mi_max,
            ring: false,
        };
    }
    let mut d2 = vec![0.0; NB * NB];
    for i in 0..NB {
        for j in 0..NB {
            if i != j {
                let m = (mi[i + j * NB] / mi_max).max(1e-300);
                let dd = -m.ln();
                d2[i + j * NB] = dd * dd;
            }
        }
    }
    let row_mean: Vec<f64> = (0..NB)
        .map(|i| (0..NB).map(|j| d2[i + j * NB]).sum::<f64>() / NB as f64)
        .collect();
    let tot: f64 = row_mean.iter().sum::<f64>() / NB as f64;
    let mut b = vec![0.0; NB * NB];
    for i in 0..NB {
        for j in 0..NB {
            b[i + j * NB] = -0.5 * (d2[i + j * NB] - row_mean[i] - row_mean[j] + tot);
        }
    }
    let (w, v) = jacobi_eigh(&b, NB);
    let (l1, l2) = (w[NB - 1], w[NB - 2]);
    let coords: Vec<(f64, f64)> = (0..NB)
        .map(|i| {
            (
                l1.max(0.0).sqrt() * v[i + (NB - 1) * NB],
                l2.max(0.0).sqrt() * v[i + (NB - 2) * NB],
            )
        })
        .collect();
    let mut order: Vec<usize> = (0..NB).collect();
    order.sort_by(|&a, &bq| {
        coords[a]
            .1
            .atan2(coords[a].0)
            .partial_cmp(&coords[bq].1.atan2(coords[bq].0))
            .unwrap()
    });
    let mut adjacent_ok = 0;
    for k in 0..NB {
        let a = order[k];
        let bq = order[(k + 1) % NB];
        let d = (a as isize - bq as isize).unsigned_abs();
        if d == 1 || d == NB - 1 {
            adjacent_ok += 1;
        }
    }
    let radii: Vec<f64> = coords
        .iter()
        .map(|&(x, y)| (x * x + y * y).sqrt())
        .collect();
    let rmean: f64 = radii.iter().sum::<f64>() / NB as f64;
    let rsd = (radii.iter().map(|r| (r - rmean).powi(2)).sum::<f64>() / NB as f64).sqrt()
        / rmean.max(1e-300);
    let adjacency = adjacent_ok as f64 / NB as f64;
    let lam21 = if l1 > 0.0 { l2 / l1 } else { 0.0 };
    let ring = adjacency >= 0.9 && rsd <= 0.10 && lam21 >= 0.9;
    Metrics {
        adjacency,
        rsd,
        lam21,
        mi_max,
        ring,
    }
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
    println!("=== v6.4 幾何読み出しの陰性対照: 正例だけが円環になるか ===\n");

    // ---- [P] 正例: 円環格子の基底状態 (v0.7 と同一) ----
    let c1 = |d: usize| -> f64 {
        if d == 0 {
            0.5
        } else {
            (std::f64::consts::PI * d as f64 / 2.0).sin()
                / (N as f64 * (std::f64::consts::PI * d as f64 / N as f64).sin())
        }
    };
    let mut c_true = vec![0.0; N * N];
    for x in 0..N {
        for y in 0..N {
            let mut d = (x as isize - y as isize).unsigned_abs();
            d = d.min(N - d);
            c_true[x + y * N] = c1(d);
        }
    }
    let mut mi_true = vec![0.0; NB * NB];
    let mm_true = block_mi(&c_true, &mut mi_true);
    let m_p = ring_metrics(&mi_true, mm_true);
    // MI 減衰べき (正例のみ意味を持つ)
    let slope = {
        let (mut xs, mut ys) = (Vec::new(), Vec::new());
        for d in 2..=40usize {
            let chord = (NB as f64 / std::f64::consts::PI)
                * (std::f64::consts::PI * d as f64 / NB as f64).sin();
            xs.push(chord.ln());
            ys.push(mi_true[d * NB].ln());
        }
        linfit(&xs, &ys).1
    };

    // ---- [N1] MI 行列の対称シャッフル ----
    let mut rng = Rng::new(6404);
    let mut vals: Vec<f64> = Vec::new();
    for i in 0..NB {
        for j in (i + 1)..NB {
            vals.push(mi_true[i + j * NB]);
        }
    }
    for k in (1..vals.len()).rev() {
        let r = rng.range(k + 1);
        vals.swap(k, r);
    }
    let mut mi_shuf = vec![0.0; NB * NB];
    let mut idx = 0;
    for i in 0..NB {
        for j in (i + 1)..NB {
            mi_shuf[i + j * NB] = vals[idx];
            mi_shuf[j + i * NB] = vals[idx];
            idx += 1;
        }
    }
    let m_n1 = ring_metrics(&mi_shuf, mm_true);

    // ---- [N2] 非局所ハミルトニアン (GOE) の基底状態 ----
    let goe_state = |rng: &mut Rng, lowest: bool| -> Vec<f64> {
        let mut h = vec![0.0; N * N];
        for i in 0..N {
            for j in 0..=i {
                let x = rng.gauss();
                h[i + j * N] = x;
                h[j + i * N] = x;
            }
        }
        let (_, v) = jacobi_eigh(&h, N);
        // 占有モード: 基底状態なら最低 N/2 個、volume-law ならランダムな N/2 個
        let mut occ: Vec<usize> = (0..N).collect();
        if !lowest {
            for k in (1..N).rev() {
                let r = rng.range(k + 1);
                occ.swap(k, r);
            }
        }
        let occ = &occ[..N / 2];
        let mut c = vec![0.0; N * N];
        for &m in occ {
            for i in 0..N {
                let vi = v[i + m * N];
                if vi == 0.0 {
                    continue;
                }
                for j in 0..N {
                    c[i + j * N] += vi * v[j + m * N];
                }
            }
        }
        c
    };
    let c_goe = goe_state(&mut rng, true);
    let mut mi_goe = vec![0.0; NB * NB];
    let mm_goe = block_mi(&c_goe, &mut mi_goe);
    let m_n2 = ring_metrics(&mi_goe, mm_goe);

    // ---- [N3] volume-law ランダム状態 (ランダム基底のランダム充填) ----
    let c_vol = goe_state(&mut rng, false);
    let mut mi_vol = vec![0.0; NB * NB];
    let mm_vol = block_mi(&c_vol, &mut mi_vol);
    let m_n3 = ring_metrics(&mi_vol, mm_vol);

    // ---- [N4] 古典的・無相関状態 (ランダム占有の積状態) ----
    let mut c_cls = vec![0.0; N * N];
    for i in 0..N {
        c_cls[i + i * N] = rng.f64();
    }
    let mut mi_cls = vec![0.0; NB * NB];
    let mm_cls = block_mi(&c_cls, &mut mi_cls);
    let m_n4 = ring_metrics(&mi_cls, mm_cls);

    // ---- 表とその判定 ----
    println!(
        "  状態                          隣接復元率   半径ばらつき  λ2/λ1   MI_max     円環判定"
    );
    let rows: [(&str, &Metrics, bool); 5] = [
        ("P  円環格子の基底状態", &m_p, true),
        ("N1 MI シャッフル", &m_n1, false),
        ("N2 非局所 H (GOE) 基底状態", &m_n2, false),
        ("N3 volume-law ランダム状態", &m_n3, false),
        ("N4 古典的・無相関 (積状態)", &m_n4, false),
    ];
    let mut all_ok = true;
    for (name, m, expect_ring) in &rows {
        println!(
            "  {:28} {:6.1}%     {:9.1}%   {:5.2}   {:.2e}   {}",
            name,
            100.0 * m.adjacency,
            100.0 * m.rsd,
            m.lam21,
            m.mi_max,
            if m.ring { "円環" } else { "なし" }
        );
        if m.ring != *expect_ring {
            all_ok = false;
        }
    }
    println!(
        "\n  MI 減衰べき (正例): {:.2} (自由フェルミオンの理論値 -2)",
        slope
    );
    let ok_slope = (slope + 2.0).abs() < 0.3;
    println!(
        "\n  => 正例のみ円環と判定される (対照 4 種は全て失敗する)  {}",
        pass(all_ok)
    );
    println!("  => 正例の MI 減衰べき -2±0.3  {}", pass(ok_slope));

    // ---- JSON artifact ----
    let met_json = |m: &Metrics| {
        Json::Obj(vec![
            ("adjacency".into(), Json::Num(m.adjacency)),
            (
                "radius_scatter".into(),
                Json::Num(if m.rsd.is_finite() { m.rsd } else { -1.0 }),
            ),
            ("lambda2_over_lambda1".into(), Json::Num(m.lam21)),
            ("mi_max".into(), Json::Num(m.mi_max)),
            ("ring_detected".into(), Json::Bool(m.ring)),
        ])
    };
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-GEOM-004".into())),
        ("seed".into(), Json::Int(6404)),
        ("lattice_size".into(), Json::Int(N as i64)),
        (
            "criterion".into(),
            Json::Str("adjacency>=0.9 && radius_scatter<=0.10 && lam2/lam1>=0.9".into()),
        ),
        ("positive".into(), met_json(&m_p)),
        (
            "controls".into(),
            Json::Arr(vec![
                Json::Obj(vec![
                    ("name".into(), Json::Str("shuffled_mi".into())),
                    ("expected".into(), Json::Str("reconstruction fails".into())),
                    ("metrics".into(), met_json(&m_n1)),
                ]),
                Json::Obj(vec![
                    ("name".into(), Json::Str("nonlocal_goe_ground_state".into())),
                    ("expected".into(), Json::Str("reconstruction fails".into())),
                    ("metrics".into(), met_json(&m_n2)),
                ]),
                Json::Obj(vec![
                    ("name".into(), Json::Str("volume_law_random_state".into())),
                    ("expected".into(), Json::Str("reconstruction fails".into())),
                    ("metrics".into(), met_json(&m_n3)),
                ]),
                Json::Obj(vec![
                    (
                        "name".into(),
                        Json::Str("classical_uncorrelated_product".into()),
                    ),
                    (
                        "expected".into(),
                        Json::Str("no geometry definable (MI=0)".into()),
                    ),
                    ("metrics".into(), met_json(&m_n4)),
                ]),
            ]),
        ),
        ("mi_decay_exponent_positive".into(), Json::Num(slope)),
        ("pass".into(), Json::Bool(all_ok && ok_slope)),
    ]);
    let p = write_artifact("results/v64_controls.json", &j.render());
    println!("  機械可読な結果: {}", p);

    println!("\n総合判定: {}", pass(all_ok && ok_slope));
    println!("\n結論: MI→MDS の幾何読み出しは、局所構造を持つ状態 (正例) だけを円環と判定し、");
    println!(
        "      同じ MI 値の集合でも配置を壊すと (N1)、もつれがあっても局所性がないと (N2,N3)、"
    );
    println!("      相関そのものがないと (N4)、それぞれ異なる仕方で失敗する。");
    println!("      「距離=相関の減衰率」という辞書 (A2) は、減衰の *局所的パターン* が");
    println!("      あって初めて幾何を生む — 都合よく何でも幾何に見えているわけではない。");
    if !(all_ok && ok_slope) {
        std::process::exit(1);
    }
}
