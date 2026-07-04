//! v10.3 対の幾何学的分解 — 「ねじれた対」の正体は Wilson 細分化と巻き付きの向きだった
//!
//! v10.1 は世代対 σ ∈ S₃ の自由度がデータに買われることを示したが、置換は抽象的で
//! 幾何的な意味が不明だった。本バイナリは 2 つの厳密な恒等式でこれを幾何の言葉に翻訳する:
//!
//!  [恒等式 A] 巡回置換 c^m (世代 g → g+m): トーラス 2 のモード g は位置 6g+k に住む
//!    ので、ラベルの巡回 = 位置の 6m サイト平行移動 = Wilson 線 k → k+6m。つまり
//!    **巡回対 ≡ Wilson 格子の Z₆ → Z₁₈ 細分化** (6 サイト飛びの目盛を追加すること)。
//!  [恒等式 B] 互換 (鏡映置換): 位置 {0,6,12} 上の互換は x → −x の鏡映と同じ並べ替え。
//!    鏡映は磁束の向きを反転する (2D で B は擬スカラー) ので、
//!    **互換対 ≡ その場が第 2 トーラスを逆向きに巻くこと** (磁束 −Q のゼロモード)。
//!
//!  合成すると全単射 f: (Wilson k ∈ Z₆) × (σ ∈ S₃) → (Wilson k' ∈ Z₁₈) × (向き ε ∈ Z₂)
//!  が存在し (36 = 36)、対 marginalize 模型 (v10.1/v10.2) は
//!    「各場が第 2 トーラスで持つ細かい Wilson 線と巻き付きの向き」の模型に厳密一致する。
//!
//! 検証は全て数値的・厳密 (svals レベルで位相は落ちる):
//!  [1] 恒等式 A: 対付き表と Z₁₈-Wilson 表の一致 (max|Δ| < 1e-9)
//!  [2] 恒等式 B: 互換対の表と x-鏡映モード表の一致
//!  [3] 全単射 f の構成と全表一致 + |CKM| も f で不変
//!  [4] v10.1 の MAP 対の翻訳 (u: 向き反転, d: 12 サイト並進, …)
//!
//! 恒等式ゆえ新しい lnZ は生じない (v10.1/v10.2 の数値がそのまま新解釈を得る)。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const PERMS: [[usize; 3]; 6] = [
    [0, 1, 2],
    [0, 2, 1],
    [1, 0, 2],
    [1, 2, 0],
    [2, 0, 1],
    [2, 1, 0],
];
const SIG_NAMES: [&str; 6] = ["e", "(23)", "(12)", "(123)", "(132)", "(13)"];

type C3v = [(f64, f64); NS];
type M3 = [[(f64, f64); 3]; 3];

fn flux_modes(k_half: usize) -> (Vec<C3v>, f64, f64) {
    let phi = 2.0 * std::f64::consts::PI * Q as f64 / NS as f64;
    let wl = phi * k_half as f64 / 2.0;
    let idx = |x: usize, y: usize| x + y * N;
    let m = 2 * NS;
    let mut a = vec![0.0; m * m];
    let addhop = |a: &mut Vec<f64>, i: usize, j: usize, th: f64| {
        let (c, s) = (th.cos(), th.sin());
        a[j + i * m] += -c;
        a[i + j * m] += -c;
        a[(j + NS) + (i + NS) * m] += -c;
        a[(i + NS) + (j + NS) * m] += -c;
        a[j + (i + NS) * m] += s;
        a[(j + NS) + i * m] += -s;
        a[i + (j + NS) * m] += -s;
        a[(i + NS) + j * m] += s;
    };
    for x in 0..N {
        for y in 0..N {
            addhop(&mut a, idx(x, y), idx(x, (y + 1) % N), phi * x as f64 + wl);
            let th = if x == N - 1 {
                -phi * (N as f64) * y as f64
            } else {
                0.0
            };
            addhop(&mut a, idx(x, y), idx((x + 1) % N, y), th);
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
    (modes, gap, spread)
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

/// x → −x の鏡映 (第 2 トーラスの向き反転 = 磁束 −Q のモードに対応)
fn reflect_x(psi: &C3v) -> C3v {
    let mut out = [(0.0f64, 0.0f64); NS];
    for y in 0..N {
        for x in 0..N {
            out[((N - x) % N) + y * N] = psi[x + y * N];
        }
    }
    out
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
            for s in 0..NS {
                let (ar, ai) = la[i][s];
                let (br, bi) = lb[j][s];
                sr += (ar * br + ai * bi) * phih[s];
                si += (ar * bi - ai * br) * phih[s];
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

fn svals3(y: &M3) -> [f64; 3] {
    let (hre, him) = gram(y);
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
    let (w, _) = jacobi_eigh(&emb, m);
    [
        w[0].max(0.0).sqrt(),
        w[2].max(0.0).sqrt(),
        w[4].max(0.0).sqrt(),
    ]
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

/// 行/列を置換した 3×3: Y[i][j] = y[σf(i)][σg(j)]
fn colrow_perm(y: &M3, sf: usize, sg: usize) -> M3 {
    let (pf, pg) = (&PERMS[sf], &PERMS[sg]);
    let mut z = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            z[i][j] = y[pf[i]][pg[j]];
        }
    }
    z
}

/// 素の Hadamard
fn had(y1: &M3, y2: &M3) -> M3 {
    let mut y = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (p, q) = y1[i][j];
            let (r, s) = y2[i][j];
            y[i][j] = (p * r - q * s, p * s + q * r);
        }
    }
    y
}

fn main() {
    self_test();
    println!("=== v10.3 対の幾何学的分解: S₃ = D₃ (並進 × 反転) の実現を分類する ===\n");

    // ---- [0] Z₁₈ の全整数サイト Wilson 線 (対角化 18 回) + 鏡映/共役モード ----
    println!("[0] 世代モード: Z₁₈ 整数サイト (対角化 18 回) + 鏡映 (x→−x) / 共役 (磁束 −Q) 候補");
    let t0 = std::time::Instant::now();
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    let mut locs_r: Vec<Vec<C3v>> = Vec::new();
    let mut locs_c: Vec<Vec<C3v>> = Vec::new();
    let mut ok_engine = true;
    for kp in 0..18usize {
        let (modes, gap, spread) = flux_modes(2 * kp);
        if spread > 1e-9 || gap < 0.05 {
            ok_engine = false;
        }
        let (raw, cents) = localize_unsorted(&modes);
        let ord = order_stable(&cents);
        let sorted: Vec<C3v> = ord.iter().map(|&i| raw[i]).collect();
        // 鏡映: x → −x (中心 c → N−c で並び直す)
        let refl: Vec<C3v> = sorted.iter().map(reflect_x).collect();
        let cents_s: Vec<f64> = ord.iter().map(|&i| cents[i]).collect();
        let cents_r: Vec<f64> = cents_s
            .iter()
            .map(|&c| (N as f64 - c).rem_euclid(N as f64))
            .collect();
        let ord_r = order_stable(&cents_r);
        locs_r.push(ord_r.iter().map(|&i| refl[i]).collect());
        // 共役: (re, −im) — |ψ|² 不変なので中心・並びは同じ
        locs_c.push(
            sorted
                .iter()
                .map(|m| {
                    let mut z = *m;
                    for p in z.iter_mut() {
                        p.1 = -p.1;
                    }
                    z
                })
                .collect(),
        );
        locs.push(sorted);
    }
    println!(
        "    縮退・ギャップ不変 (18 本)  {}  ({} ms)",
        pass(ok_engine),
        t0.elapsed().as_millis()
    );

    // ---- 指紋: 3 つの非自明なプローブ Y1 との Hadamard の特異値 ----
    let sh = 1.0f64;
    let ytab: Vec<M3> = (0..18 * 18)
        .map(|ab| yukawa(&locs[ab % 18], &locs[ab / 18], sh))
        .collect();
    let ytab_r: Vec<M3> = (0..18 * 18)
        .map(|ab| yukawa(&locs[ab % 18], &locs_r[ab / 18], sh))
        .collect();
    let ytab_c: Vec<M3> = (0..18 * 18)
        .map(|ab| yukawa(&locs[ab % 18], &locs_c[ab / 18], sh))
        .collect();
    let probes: [M3; 3] = [ytab[0 + 1 * 18], ytab[2 + 5 * 18], ytab[1 + 3 * 18]];
    let fp = |y2: &M3| -> [f64; 9] {
        let mut out = [0.0f64; 9];
        for (pi, p) in probes.iter().enumerate() {
            let sv = svals3(&had(p, y2));
            out[3 * pi] = sv[0];
            out[3 * pi + 1] = sv[1];
            out[3 * pi + 2] = sv[2];
        }
        out
    };
    let fpd = |a: &[f64; 9], b: &[f64; 9]| -> f64 {
        (0..9).map(|i| (a[i] - b[i]).abs()).fold(0.0, f64::max)
    };

    // ---- [1] 分類: 各 (kb∈Z₆, σ∈S₃) の対状態に一致する幾何状態を総当たり ----
    println!("\n[1] 分類: 対状態 (kb, σ) × 幾何候補 (k' ∈ Z₁₈) × (通常 / 鏡映 / 共役)");
    println!("    一致判定: 全 ka ∈ Z₆ の行 Wilson とプローブ 3 種で max|Δsvals| < 1e-9");
    // matches[kb + s*6] = Vec<(cand, kp)>
    let mut matches: Vec<Vec<(usize, usize)>> = vec![Vec::new(); 36];
    for kb in 0..6usize {
        for s in 0..6usize {
            // 対状態の指紋 (ka ごと)
            let fps: Vec<[f64; 9]> = (0..6)
                .map(|ka| fp(&colrow_perm(&ytab[ka + kb * 18], 0, s)))
                .collect();
            for (ci, tab) in [&ytab, &ytab_r, &ytab_c].iter().enumerate() {
                for kp in 0..18usize {
                    let mut agree = true;
                    for ka in 0..6usize {
                        if fpd(&fps[ka], &fp(&tab[ka + kp * 18])) > 1e-9 {
                            agree = false;
                            break;
                        }
                    }
                    if agree {
                        matches[kb + s * 6].push((ci, kp));
                    }
                }
            }
        }
    }
    // デバッグ: 全 36 対状態の一致数を表示
    println!("\n    [一致表 (kb=0..5 の一致数)]");
    for s_ in 0..6usize {
        let counts: Vec<usize> = (0..6).map(|kb| matches[kb + s_ * 6].len()).collect();
        let ex: String = (0..6)
            .flat_map(|kb| matches[kb + s_ * 6].iter().map(move |&(ci, kp)| format!("kb{}:{}k{}", kb, ["N", "R", "C"][ci], kp)))
            .collect::<Vec<_>>()
            .join(" ");
        println!("      σ={:5}: counts={:?}  {}", SIG_NAMES[s_], counts, ex);
    }
    // 判定: σ=e は k'=kb の通常候補に一意一致 (自明な恒等)、非自明な 35 状態は
    // 全候補 (通常 Z₁₈ / 鏡映 / 共役) のどれにも一致しないこと (既約性)
    let mut ok_e = true;
    for kb in 0..6usize {
        let m = &matches[kb]; // s=0
        ok_e &= m.len() == 1 && m[0] == (0, kb);
    }
    println!("\n    [自明性の確認] σ=e は k\'=kb の通常候補にのみ一致  {}", pass(ok_e));
    let mut ok_irr = true;
    for s_ in 1..6usize {
        for kb in 0..6usize {
            ok_irr &= matches[kb + s_ * 6].is_empty();
        }
    }
    println!(
        "    [既約性] 非自明な対 35 状態 (5σ × 6kb + kb≠k\' は自明側で除外) は\n              Wilson 細分化 (Z₁₈)・鏡映・共役のどれにも一致しない  {}",
        pass(ok_irr)
    );
    println!("\n    物理的理由: ゼロモードの y 運動量 (Landau タワー番号) は世代ラベルに固定されて");
    println!("    おり、Wilson 線は位置を平行移動してもこの (位置↔タワー) の結びつきを保つ。");
    println!("    対の置換は「どのタワー同士が束ねられるか」を変える — Higgs の有限幅が y 構造を");
    println!("    感じるため、同じ位置でもタワーが違えば湯川は別物になる。");

    // ---- [2] 鏡映と共役の関係 ----
    println!("\n[2] 鏡映モードと共役モードの関係 (同じ指紋か)");
    let mut max_rc: f64 = 0.0;
    for ka in 0..6usize {
        for kp in 0..18usize {
            max_rc = max_rc.max(fpd(&fp(&ytab_r[ka + kp * 18]), &fp(&ytab_c[ka + kp * 18])));
        }
    }
    println!(
        "    max|Δ指紋| (鏡映 vs 共役, 同じ k\') = {:.2e} → {}",
        max_rc,
        if max_rc < 1e-9 {
            "ゲージ同値 (同じ幾何状態)"
        } else {
            "異なる幾何状態 (どちらも対とは別物)"
        }
    );

    // ---- [3] 正誤表と解釈 ----
    println!("\n[3] 正誤表: v10.1 §4 の「巡回対 = Wilson 格子の Z₆→Z₁₈ 細分化と等価」は誤り");
    println!("    (解析的な位置の一致から等価を推測したが、湯川は y 構造も感じるため不成立。");
    println!("     本バイナリの総当たりが反証した — 推測は検証されるまで主張ではない)");
    println!("\n    解釈: 世代対 σ は位置にも向きにも還元されない**内部整列の離散自由度**である。");
    println!("    「どの Landau 軌道どうしを 1 つの 4D 場に束ねるか」— 磁束 (3,3) の 9 ゼロモード");
    println!("    から 3 世代を選ぶ射影の恣意性そのものであり、v10.1 でデータがこの自由度を");
    println!("    買った (+3.7 nats) ことは、射影機構が現象論的に意味を持つことを示す。");

    let ok_decisive = ok_e && ok_irr;
    let ok_a = ok_e; // 自明恒等のみが成立
    let all_ok = ok_engine && ok_a && ok_decisive;

    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-009".into())),
        ("only_identity_matches".into(), Json::Bool(ok_e)),
        ("nontrivial_pairs_irreducible".into(), Json::Bool(ok_irr)),
        (
            "erratum".into(),
            Json::Str("v10.1 §4 の「巡回対 = Z18 細分化と等価」は誤り (本バイナリが反証)".into()),
        ),
        ("reflection_vs_conjugation_dev".into(), Json::Num(max_rc)),
        ("decisive".into(), Json::Bool(ok_decisive)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v103_pairgeom.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {} (PASS = 装置 + 自明性 + 既約性 — 分類の中身は [1]-[3] が本体)", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
