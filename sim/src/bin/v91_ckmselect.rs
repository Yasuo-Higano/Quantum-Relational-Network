//! v9.1 CKM 込みの幾何選択 — 「質量は次元を数え、混合は傾きを測る」の緊張を解く (v9.0 残高 2/9)
//!
//! v8.1 の発見: 質量比 6 つの証拠は T³ を選ぶが、その MAP は CKM を失う。原因は
//! 尤度に混合が入っていないこと。本バイナリは**全 9 量 (質量 6 + CKM 3)** を尤度に
//! 入れた幾何選択を行う。CKM は u/d セクターを K_Q で結合するため尤度の完全な
//! 因子化は壊れるが、
//!   Z₉ = (1/|K|⁵|σ|) Σ_σ [ Σ_{K_Q} Σ_{K_u,K_d} L_u·L_d·L_CKM ] × [ Σ_{K_L,K_e} L_e ]
//! と e セクターだけ因子化を保てば、クォーク部は (K_Q,K_u,K_d) の三重和で厳密に
//! 計算できる (ペアごとの湯川・特異値・固有ベクトルをキャッシュ; T³×Z₆ で 216³ ≈ 10⁷)。
//!
//! 設計判断 (正直に): 全 9 量を学習に使うので out-of-sample 検証は消える。本バイナリは
//! 「どの幾何が全データを説明するか」という選択の問いに答えるもので、予測の検証は
//! v7.2/v8.1 (質量のみ学習・CKM 予測) が引き続き担う。
//!
//! 退化検査: CKM 項の重みを 0 にすると、三重和経路の lnZ が v8.1 の因子化経路の値
//! (T¹Z₆ −53.77 / T²Z₆ −20.41 / T²Z₁₂ −20.42 / T³Z₆ −17.41) と厳密一致すること。
//! T³×Z₁₂ は三重和が 1728³ ≈ 5×10⁹ で重いため対象外 (明示的な打ち切り)。

use uft_sim::*;

const N: usize = 18;
const NS: usize = N * N;
const Q: usize = 3;
const NK12: usize = 12;
const EPS_OBS: [f64; 9] = [
    1.3e-5, 3.7e-3, 1.1e-3, 2.2e-2, 2.9e-4, 5.9e-2, 0.225, 0.041, 0.0037,
];

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

fn localize(modes: &[C3v]) -> Vec<C3v> {
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
    let mut ord: Vec<usize> = (0..Q).collect();
    ord.sort_by(|&a, &b| centers[a].partial_cmp(&centers[b]).unwrap());
    ord.iter().map(|&i| out[i]).collect()
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

fn had_prod(a: &M3, b: &M3) -> M3 {
    let mut y = [[(0.0f64, 0.0f64); 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (p, q) = a[i][j];
            let (r, s) = b[i][j];
            y[i][j] = (p * r - q * s, p * s + q * r);
        }
    }
    y
}

/// Y から (質量比 2 つの対数, 左固有ベクトル) を返す
fn mass_and_vecs(y: &M3) -> ([f64; 2], M3) {
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
    let (lam, vecs) = eig_herm3(&hre, &him);
    let sv = [
        lam[0].max(0.0).sqrt(),
        lam[1].max(0.0).sqrt(),
        lam[2].max(0.0).sqrt(),
    ];
    (
        [
            (sv[0].max(1e-300) / sv[2].max(1e-300)).ln(),
            (sv[1].max(1e-300) / sv[2].max(1e-300)).ln(),
        ],
        vecs,
    )
}

/// CKM の 3 要素 |V_us|, |V_cb|, |V_ub|
fn ckm3(vu: &M3, vd: &M3) -> [f64; 3] {
    let mut ckm = [[0.0f64; 3]; 3];
    for i in 0..3 {
        for j in 0..3 {
            let (mut re, mut im) = (0.0, 0.0);
            for k in 0..3 {
                let (a, b) = vu[i][k];
                let (c, d) = vd[j][k];
                re += a * c + b * d;
                im += a * d - b * c;
            }
            ckm[i][j] = (re * re + im * im).sqrt();
        }
    }
    [ckm[0][1], ckm[1][2], ckm[0][2]]
}

fn lse(v: &[f64]) -> f64 {
    let m = v.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    m + v.iter().map(|&x| (x - m).exp()).sum::<f64>().ln()
}

fn pass(ok: bool) -> &'static str {
    if ok {
        "[PASS]"
    } else {
        "[FAIL]"
    }
}

struct Out9 {
    lnz: f64,
    map_ll: f64,
    preds: [f64; 9],
}

/// 全 9 量 (use_ckm=true) または質量のみ (false; 退化検査用) の証拠と MAP。
fn eval9(
    nt: usize,
    ks: &[usize],
    locs: &[Vec<C3v>],
    sig_grid: &[f64],
    sigma: f64,
    use_ckm: bool,
) -> Out9 {
    let nk = ks.len();
    let nc = nk.pow(nt as u32);
    let decoded: Vec<Vec<usize>> = (0..nc)
        .map(|mut c| {
            let mut v = Vec::with_capacity(nt);
            for _ in 0..nt {
                v.push(ks[c % nk]);
                c /= nk;
            }
            v
        })
        .collect();
    let norm1 = -(sigma * (2.0 * std::f64::consts::PI).sqrt()).ln(); // 観測 1 つ分
    let tgt: Vec<f64> = EPS_OBS.iter().map(|x| x.ln()).collect();
    let mut lnz_terms = Vec::new();
    let mut best = (f64::NEG_INFINITY, 0.0f64, [0usize; 5]);
    for &sh in sig_grid {
        let ytab: Vec<M3> = (0..NK12 * NK12)
            .map(|ab| yukawa(&locs[ab % NK12], &locs[ab / NK12], sh))
            .collect();
        // ペアキャッシュ: (A,B) → (質量対数比 2, 左固有ベクトル)
        let pair: Vec<([f64; 2], M3)> = (0..nc * nc)
            .map(|ab| {
                let (a, b) = (ab % nc, ab / nc);
                let mut y = ytab[decoded[a][0] + decoded[b][0] * NK12];
                for t in 1..nt {
                    y = had_prod(&y, &ytab[decoded[a][t] + decoded[b][t] * NK12]);
                }
                mass_and_vecs(&y)
            })
            .collect();
        // e セクター (因子化)
        let le: Vec<f64> = (0..nc * nc)
            .map(|ab| {
                let r = &pair[ab].0;
                -((r[0] - tgt[4]).powi(2) + (r[1] - tgt[5]).powi(2)) / (2.0 * sigma * sigma)
                    + 2.0 * norm1
            })
            .collect();
        // e セクターの最良 (MAP 用に 1 回だけ計算)
        let mut le_best = f64::NEG_INFINITY;
        let mut le_besti = 0usize;
        for (i, &v) in le.iter().enumerate() {
            if v > le_best {
                le_best = v;
                le_besti = i;
            }
        }
        // クォーク部: KQ ごとに (Ku, Kd) の二重和
        let mut per_q = Vec::with_capacity(nc);
        for kq in 0..nc {
            let mut acc_m = f64::NEG_INFINITY;
            let mut acc_s = 0.0f64;
            for ku in 0..nc {
                let (ru, vu) = &pair[kq + ku * nc];
                let llu = -((ru[0] - tgt[0]).powi(2) + (ru[1] - tgt[1]).powi(2))
                    / (2.0 * sigma * sigma)
                    + 2.0 * norm1;
                for kd in 0..nc {
                    let (rd, vd) = &pair[kq + kd * nc];
                    let lld = -((rd[0] - tgt[2]).powi(2) + (rd[1] - tgt[3]).powi(2))
                        / (2.0 * sigma * sigma)
                        + 2.0 * norm1;
                    let mut ll = llu + lld;
                    if use_ckm {
                        let c = ckm3(vu, vd);
                        for m in 0..3 {
                            let d = c[m].max(1e-300).ln() - tgt[6 + m];
                            ll += -d * d / (2.0 * sigma * sigma) + norm1;
                        }
                    }
                    // ストリーミング logsumexp
                    if ll > acc_m {
                        acc_s = acc_s * (acc_m - ll).exp() + 1.0;
                        acc_m = ll;
                    } else {
                        acc_s += (ll - acc_m).exp();
                    }
                    // MAP (e の最良と合成)
                    let tot = ll + le_best;
                    if tot > best.0 {
                        best = (tot, sh, [kq, ku, kd, le_besti % nc, le_besti / nc]);
                    }
                }
            }
            per_q.push(acc_m + acc_s.ln());
        }
        lnz_terms.push(lse(&per_q) + lse(&le));
    }
    let nkf = nk as f64;
    let lnz = lse(&lnz_terms) - (5.0 * (nt as f64) * nkf.ln() + (sig_grid.len() as f64).ln());
    // MAP の 9 量
    let sh = best.1;
    let build = |a: usize, b: usize| -> M3 {
        let mut y = yukawa(&locs[decoded[a][0]], &locs[decoded[b][0]], sh);
        for t in 1..nt {
            y = had_prod(&y, &yukawa(&locs[decoded[a][t]], &locs[decoded[b][t]], sh));
        }
        y
    };
    let yu = build(best.2[0], best.2[1]);
    let yd = build(best.2[0], best.2[2]);
    let ye = build(best.2[3], best.2[4]);
    let (ru, vu) = mass_and_vecs(&yu);
    let (rd, vd) = mass_and_vecs(&yd);
    let (re_, _) = mass_and_vecs(&ye);
    let c = ckm3(&vu, &vd);
    let preds = [
        ru[0].exp(),
        ru[1].exp(),
        rd[0].exp(),
        rd[1].exp(),
        re_[0].exp(),
        re_[1].exp(),
        c[0],
        c[1],
        c[2],
    ];
    Out9 {
        lnz,
        map_ll: best.0,
        preds,
    }
}

fn main() {
    self_test();
    println!("=== v9.1 CKM 込みの幾何選択: 質量と混合の緊張を解く ===\n");
    let sigma = (2.0f64).ln();
    let sig_grid = [1.0f64, 1.5, 2.0, 2.5];
    let names = [
        "m_u/m_t", "m_c/m_t", "m_d/m_b", "m_s/m_b", "m_e/m_τ", "m_μ/m_τ", "|V_us|", "|V_cb|",
        "|V_ub|",
    ];
    println!("[1] 世代モード (Z₁₂, 対角化 12 回)");
    let t0 = std::time::Instant::now();
    let mut locs: Vec<Vec<C3v>> = Vec::new();
    let mut ok_engine = true;
    for k in 0..NK12 {
        let (modes, gap, spread) = flux_modes(k);
        if spread > 1e-9 || gap < 0.05 {
            ok_engine = false;
        }
        locs.push(localize(&modes));
    }
    println!("    縮退・ギャップ不変  {}  ({} ms)", pass(ok_engine), t0.elapsed().as_millis());

    let z6: Vec<usize> = (0..NK12).step_by(2).collect();
    let z12: Vec<usize> = (0..NK12).collect();

    // ---- [2] 退化検査: CKM を外すと v8.1 の値に一致 ----
    println!("\n[2] 退化検査 (CKM 項の重み 0 → v8.1 の質量のみ証拠と一致するはず)");
    let v81_vals = [
        ("T¹×Z₆", 1usize, &z6, -53.77f64),
        ("T²×Z₆", 2, &z6, -20.41),
        ("T³×Z₆", 3, &z6, -17.41),
    ];
    let mut ok_degen = true;
    for (label, nt, ks, expect) in v81_vals {
        let o = eval9(nt, ks, &locs, &sig_grid, sigma, false);
        let ok = (o.lnz - expect).abs() < 0.05;
        if !ok {
            ok_degen = false;
        }
        println!("    {}: {:7.2} vs v8.1 {:7.2}  {}", label, o.lnz, expect, pass(ok));
    }

    // ---- [3] 全 9 量の証拠 ----
    println!("\n[3] 全 9 量 (質量 6 + CKM 3) の証拠 — T³×Z₁₂ は三重和 5×10⁹ のため対象外 (明示)");
    let mut rows: Vec<(String, Out9)> = Vec::new();
    for (label, nt, ks) in [
        ("T¹ × Z₆ ", 1usize, &z6),
        ("T² × Z₆ ", 2, &z6),
        ("T² × Z₁₂", 2, &z12),
        ("T³ × Z₆ ", 3, &z6),
    ] {
        let t1 = std::time::Instant::now();
        let o = eval9(nt, ks, &locs, &sig_grid, sigma, true);
        println!(
            "    {}: lnZ₉ = {:7.2}  (MAP lnL {:7.2}, {} ms)",
            label,
            o.lnz,
            o.map_ll,
            t1.elapsed().as_millis()
        );
        rows.push((label.trim().to_string(), o));
    }
    let winner = rows
        .iter()
        .enumerate()
        .max_by(|a, b| a.1 .1.lnz.partial_cmp(&b.1 .1.lnz).unwrap())
        .unwrap()
        .0;
    println!("    => 全 9 量の勝者: {} (lnZ₉ = {:.2})", rows[winner].0, rows[winner].1.lnz);

    // ---- [4] 勝者の MAP ----
    println!("\n[4] 勝者 ({}) の MAP 予測 (全 9 量学習 — out-of-sample ではない点に注意)", rows[winner].0);
    let mut ok9 = 0;
    for k in 0..9 {
        let ratio = rows[winner].1.preds[k] / EPS_OBS[k];
        let within = ratio > 0.2 && ratio < 5.0;
        if within {
            ok9 += 1;
        }
        println!(
            "    {:8}  予測 {:9.2e}   実測 {:8.2e}   比 {:6.2} {}",
            names[k],
            rows[winner].1.preds[k],
            EPS_OBS[k],
            ratio,
            if within { "✓" } else { " " }
        );
    }
    println!("    => 9 量中 {} が 5 倍以内", ok9);

    // ---- JSON / 判定 ----
    let all_ok = ok_engine && ok_degen;
    let j = Json::Obj(vec![
        ("claim_id".into(), Json::Str("QRN-YUK-005".into())),
        (
            "models".into(),
            Json::Arr(
                rows.iter()
                    .map(|(l, o)| {
                        Json::Obj(vec![
                            ("label".into(), Json::Str(l.clone())),
                            ("lnZ9".into(), Json::Num(o.lnz)),
                            ("map_lnL".into(), Json::Num(o.map_ll)),
                        ])
                    })
                    .collect(),
            ),
        ),
        ("winner".into(), Json::Str(rows[winner].0.clone())),
        ("winner_within_factor5".into(), Json::Int(ok9)),
        ("pass".into(), Json::Bool(all_ok)),
    ]);
    let p = write_artifact("results/v91_ckmselect.json", &j.render());
    println!("\n  機械可読な結果: {}", p);
    println!("\n総合判定: {} (幾何選択の結果は [3] の表が本体)", pass(all_ok));
    if !all_ok {
        std::process::exit(1);
    }
}
