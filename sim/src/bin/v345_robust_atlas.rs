//! v34.5 robust atlas — exact reader の同時信頼集合への持ち上げ (PROMPT/15 §5)
//!
//! 背景: 第三十三期の証明書 (addressability σ_min・cross-talk・glue overlap・
//! charge witness) は exact 値の申告だった。本版は v34.3 の同時信頼集合意味論を
//! v33 の各 reader に持ち上げる — 統計量ごとの区間 (有限ショット CP または宣言
//! 誤差の Bonferroni 同時割当) の**全域で裁定が一致するときだけ**回答し、バーを
//! 跨げば Straddled を返す。
//!
//! 持ち上げの器械 (各統計量 → 区間 → 裁定):
//!   [R0] 合成 reader の被覆定理の実例 — 2 統計量 (σ 系 + xtalk 系, Bonferroni
//!        α/2) の全結果空間 31² を厳密列挙し P(誤証明書) ≤ α を機械確認
//!   [R1] addressability — σ_min の同時下界 σ_lo = σ_min(M̂) − ‖ΔM‖_F (Weyl) と
//!        cross-talk の worst-case 上界 (平均ではない — 禁止変換 26)。資格 /
//!        tied 拒否 / Straddled の 3 セル + box 角の全数で下界の厳密性を照合
//!   [R2] glue — 候補 orbit 間距離 (overlap) の区間裁定: 同一 orbit (o = 1)・
//!        非同値 (o = 1/3, v34.4 の site×CNOT 値)・境界跨ぎ → Straddled
//!   [R3] charge witness — skew 作用素 K の spectral-gap 証明書: 区間全域で
//!        σ_min(K) > 0 (zero crossing なし) のときのみ J = K(−K²)^{−1/2} を構成
//!        (J² = −I 機械検証)。縮退 witness (部分 charge) は区間が 0 を含み拒否
//!   [R4] resource cost — interval cost と budget 半順序: 確実採用 (hi ≤ b)・
//!        確実排除 (lo > b)・跨ぎ (点 Straddled) の 3 値で set-valued profile。
//!        中点潰し (点推定) は跨ぎ点で読みを反転する (禁止変換 18/22 の合流)
//!   [R5] structured 一致 — 同じ区間意味論を dense / Pauli GF(2) の両 lane に
//!        載せ、edge/non-edge/Straddled・成分・閉包/中心の対応 (2^{dim V},
//!        2^{radical}) がセルごとに一致 (v33.6 の対応原理の confidence 版)
//!
//! 統計モデル (登録契約): 各統計量 s ∈ [0,1] は N ショットの Bernoulli(s) で
//! 観測され (synthetic lane — 決定的な代表カウント k = round(N·s) を記録として
//! 使う)、区間は Clopper–Pearson (v34.3)。宣言誤差の box は Bonferroni 同時割当。
//! これは synthetic であり実測ノイズではない (v34.6 の real-data lane と区別)。

use uft_sim::factorization_enumerator::same_candidate_orbit;
use uft_sim::finite_data::{binom_pmf, cp_interval};
use uft_sim::graded_recovery::MajoranaFrame;
use uft_sim::operational_net::{algebra_closure, closure_center_basis, commutator};
use uft_sim::structured_backend::PauliVector;
use uft_sim::{jacobi_eigh, self_test, C64};

// ---------------------------------------------------------------- 小道具

fn zeros(n: usize) -> Vec<C64> {
    vec![C64::new(0.0, 0.0); n * n]
}
fn eye(n: usize) -> Vec<C64> {
    let mut m = zeros(n);
    for i in 0..n {
        m[i * n + i] = C64::new(1.0, 0.0);
    }
    m
}
fn kron(a: &[C64], na: usize, b: &[C64], nb: usize) -> Vec<C64> {
    let n = na * nb;
    let mut c = zeros(n);
    for i in 0..na {
        for j in 0..na {
            for k in 0..nb {
                for l in 0..nb {
                    c[(i * nb + k) * n + (j * nb + l)] = a[i * na + j] * b[k * nb + l];
                }
            }
        }
    }
    c
}
fn pauli2(which: char) -> Vec<C64> {
    match which {
        'X' => vec![
            C64::new(0.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
        ],
        'Y' => vec![
            C64::new(0.0, 0.0),
            C64::new(0.0, -1.0),
            C64::new(0.0, 1.0),
            C64::new(0.0, 0.0),
        ],
        'Z' => vec![
            C64::new(1.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(0.0, 0.0),
            C64::new(-1.0, 0.0),
        ],
        _ => eye(2),
    }
}
/// Pauli 文字列 → dense 行列
fn pauli_dense(s: &str) -> Vec<C64> {
    let mut m = vec![C64::new(1.0, 0.0)];
    let mut cur = 1usize;
    for ch in s.chars() {
        let p = pauli2(ch);
        let next = kron(&m, cur, &p, 2);
        cur *= 2;
        m = next;
    }
    m
}
fn hs_norm(a: &[C64]) -> f64 {
    a.iter().map(|x| x.norm2()).sum::<f64>().sqrt()
}
/// 2×2 実行列の最小特異値 (閉形式: MᵀM の固有値)
fn sigma_min_2x2(m: &[f64; 4]) -> f64 {
    let g = [
        m[0] * m[0] + m[2] * m[2],
        m[0] * m[1] + m[2] * m[3],
        m[1] * m[0] + m[3] * m[2],
        m[1] * m[1] + m[3] * m[3],
    ];
    let tr = g[0] + g[3];
    let disc = ((g[0] - g[3]) * (g[0] - g[3]) + 4.0 * g[1] * g[2]).sqrt();
    ((tr - disc) / 2.0).max(0.0).sqrt()
}
/// 実対称行列の固有値 (昇順)
fn eig_sym(a: &[f64], n: usize) -> Vec<f64> {
    let (mut ev, _) = jacobi_eigh(a, n);
    ev.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ev
}

/// 3 値の区間裁定
#[derive(PartialEq, Debug, Clone, Copy)]
enum IntervalVerdict {
    Above,
    Below,
    Straddled,
}
fn interval_vs_bar(lo: f64, hi: f64, bar: f64) -> IntervalVerdict {
    if lo > bar {
        IntervalVerdict::Above
    } else if hi <= bar {
        IntervalVerdict::Below
    } else {
        IntervalVerdict::Straddled
    }
}

fn main() {
    self_test();
    println!("=== v34.5 robust atlas — exact reader の同時信頼集合への持ち上げ (PROMPT/15 §5) ===\n");
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

    // ---------------- [R0] 合成 reader の被覆 (厳密列挙) ----------------
    println!("[R0] 合成 robust reader の被覆定理の実例 (全結果空間の厳密列挙)");
    {
        // 2 統計量: s1 (強度系, バー 0.5 を上回ることが資格) と s2 (漏れ系,
        // バー 0.1 を下回ることが資格)。各 N = 30 ショット・Bonferroni α/2。
        let (n, alpha, bar1, bar2) = (30usize, 0.05f64, 0.5f64, 0.1f64);
        let certify = |k1: usize, k2: usize| -> bool {
            let (lo1, _) = cp_interval(k1, n, alpha / 2.0);
            let (_, hi2) = cp_interval(k2, n, alpha / 2.0);
            lo1 > bar1 && hi2 <= bar2
        };
        // 真値 (資格が偽になる側の境界近傍) で誤証明書確率を厳密和
        let mut worst = 0.0f64;
        for &(s1, s2) in &[
            (0.49, 0.05),  // 強度が実は不足
            (0.7, 0.101),  // 漏れが実は超過
            (0.5, 0.1),    // 両方境界上
            (0.49, 0.101), // 両方不足
        ] {
            let qualified_truth = s1 > bar1 && s2 <= bar2;
            if qualified_truth {
                continue;
            }
            let mut p_wrong = 0.0;
            for k1 in 0..=n {
                for k2 in 0..=n {
                    if certify(k1, k2) {
                        p_wrong += binom_pmf(n, k1, s1) * binom_pmf(n, k2, s2);
                    }
                }
            }
            worst = worst.max(p_wrong);
        }
        check(
            "[R0] 資格が偽の全境界セルで P(誤証明書) ≤ α (31² 結果の厳密和・Bonferroni α/2)",
            worst <= 0.05 + 1e-12,
            format!("max P(誤証明書) = {:.5}", worst),
        );
    }

    // ---------------- [R1] robust addressability ----------------
    println!("\n[R1] addressability — σ_min 同時下界 (Weyl) と worst-case cross-talk");
    {
        let (sigma_bar, xtalk_bar) = (0.5f64, 0.1f64);
        // 応答行列 M̂ (2 command × 2 標的) と成分ごとの同時半幅 w → ‖ΔM‖_F ≤ 2w
        let read = |m: &[f64; 4], w: f64| -> (IntervalVerdict, IntervalVerdict, f64, f64) {
            let s = sigma_min_2x2(m);
            let rf = 2.0 * w; // sqrt(4 成分 × w²)
            let s_verdict = interval_vs_bar(s - rf, s + rf, sigma_bar);
            let xt = m[1].abs().max(m[2].abs());
            let xt_verdict = match interval_vs_bar(xt - w, xt + w, xtalk_bar) {
                IntervalVerdict::Above => IntervalVerdict::Above, // 漏れ超過 (拒否)
                IntervalVerdict::Below => IntervalVerdict::Below, // 漏れ資格
                IntervalVerdict::Straddled => IntervalVerdict::Straddled,
            };
            (s_verdict, xt_verdict, s - rf, xt + w)
        };
        // セル 1: 明確な資格
        let (sv1, xv1, slo1, xhi1) = read(&[1.0, 0.03, 0.03, 1.0], 0.02);
        // セル 2: tied (rank 1) — σ 区間が完全にバー未満で確実拒否
        let (sv2, _, _, _) = read(&[0.7071, 0.7071, 0.7071, 0.7071], 0.02);
        // セル 3: 境界跨ぎ → Straddled (強制判定なし)
        let (sv3, _, _, _) = read(&[0.55, 0.02, 0.02, 0.55], 0.04);
        check(
            "[R1a] 3 セル: 資格 (σ_lo > 0.5 ∧ xtalk_hi ≤ 0.1) / tied 確実拒否 / 跨ぎ Straddled",
            sv1 == IntervalVerdict::Above
                && xv1 == IntervalVerdict::Below
                && sv2 == IntervalVerdict::Below
                && sv3 == IntervalVerdict::Straddled,
            format!(
                "資格セル: σ_lo = {:.3}, xtalk_hi = {:.3} / tied: σ 区間 ⊂ [0, 0.5) / 跨ぎ: Straddled",
                slo1, xhi1
            ),
        );
        // 下界の厳密性: box の全 16 角で σ_min ≥ σ_lo
        let m0 = [1.0f64, 0.03, 0.03, 1.0];
        let w = 0.02;
        let rf = 2.0 * w;
        let s_lo = sigma_min_2x2(&m0) - rf;
        let mut min_corner = f64::INFINITY;
        for mask in 0..16u32 {
            let mut mc = m0;
            for (b, e) in mc.iter_mut().enumerate() {
                *e += if (mask >> b) & 1 == 1 { w } else { -w };
            }
            min_corner = min_corner.min(sigma_min_2x2(&mc));
        }
        check(
            "[R1b] 同時下界の厳密性: box 全 16 角の σ_min ≥ σ_lo (Weyl ‖ΔM‖_F 束縛)",
            min_corner >= s_lo - 1e-12,
            format!("min corner = {:.4} ≥ σ_lo = {:.4}", min_corner, s_lo),
        );
        // worst-case vs 平均 (禁止変換 26 の robust 版): 平均が通る box でも
        // worst-case 上界は拒否する
        let leak = [1.0f64, 0.30, 0.01, 1.0];
        let mean_offdiag = (leak[1].abs() + leak[2].abs()) / 2.0;
        let (_, xv, _, xhi) = read(&leak, 0.02);
        check(
            "[R1c] worst-case 上界: 平均 0.155…でなく max+w = 0.32 で判定 → 拒否 (平均証明書は禁止)",
            mean_offdiag < 0.2 && xv == IntervalVerdict::Above && xhi > xtalk_bar,
            format!("mean = {:.3} / worst-case 上界 = {:.3} > {}", mean_offdiag, xhi, xtalk_bar),
        );
    }

    // ---------------- [R2] robust glue (orbit 距離の区間) ----------------
    println!("\n[R2] glue — 候補 orbit 間 overlap の区間裁定");
    {
        // 統計量: v34.4 の overlap (site vs site = 1.0 / site vs CNOT 共役 = 1/3)
        let mk_comp = |strs: [&str; 3]| -> Vec<Vec<C64>> {
            strs.iter()
                .map(|s| {
                    let m = pauli_dense(s);
                    let inv = 1.0 / hs_norm(&m);
                    m.iter().map(|x| x.scale(inv)).collect()
                })
                .collect()
        };
        let site = vec![mk_comp(["XI", "YI", "ZI"]), mk_comp(["IX", "IY", "IZ"])];
        // CNOT 共役: X⊗I → X⊗X, Z⊗I → Z⊗I / I⊗X → I⊗X, I⊗Z → Z⊗Z
        let bell = vec![mk_comp(["XX", "YX", "ZI"]), mk_comp(["IX", "ZY", "ZZ"])];
        let (_, o_same) = same_candidate_orbit(&site, &site, 0.9);
        let (_, o_diff) = same_candidate_orbit(&site, &bell, 0.9);
        let orbit_bar = 0.9;
        let w = 0.03; // 宣言半幅 (同時)
        let v_same = interval_vs_bar(o_same - w, o_same + w, orbit_bar);
        let v_diff = interval_vs_bar(o_diff - w, o_diff + w, orbit_bar);
        let v_borderline = interval_vs_bar(0.88 - 0.05, 0.88 + 0.05, orbit_bar);
        check(
            "[R2] 3 セル: matching (o=1)・非同値 (o=1/3)・境界 0.88±0.05 → Straddled",
            v_same == IntervalVerdict::Above
                && v_diff == IntervalVerdict::Below
                && v_borderline == IntervalVerdict::Straddled
                && (o_diff - 1.0 / 3.0).abs() < 1e-9,
            format!(
                "o(same) = {:.6}, o(site,bell) = {:.6} (= 1/3 厳密), 境界セル = Straddled",
                o_same, o_diff
            ),
        );
    }

    // ---------------- [R3] charge witness の spectral gap と J ----------------
    println!("\n[R3] charge witness — skew K の spectral-gap 証明書と J の構成");
    {
        // 3 モード (6 Majorana, dim 8): JW frame
        let d = 8usize;
        let gammas: Vec<Vec<C64>> = vec![
            pauli_dense("XII"),
            pauli_dense("YII"),
            pauli_dense("ZXI"),
            pauli_dense("ZYI"),
            pauli_dense("ZZX"),
            pauli_dense("ZZY"),
        ];
        let frame = MajoranaFrame::certify(gammas, d).expect("CAR frame");
        let m = frame.n_majorana();
        let gnorm2 = d as f64;
        // K の抽出 (witness Q から): K_{ba} = ⟨γ_b, i[Q, γ_a]⟩ / ‖γ‖²
        let extract_k = |q: &[C64]| -> Vec<f64> {
            let mut kmat = vec![0.0f64; m * m];
            for a in 0..m {
                let c = commutator(q, frame.gamma(a), d);
                let ma: Vec<C64> = c.iter().map(|x| C64::new(-x.im, x.re)).collect();
                for b in 0..m {
                    let mut ip = C64::new(0.0, 0.0);
                    for (x, y) in frame.gamma(b).iter().zip(&ma) {
                        ip = ip + x.conj() * *y;
                    }
                    kmat[b * m + a] = ip.re / gnorm2;
                }
            }
            kmat
        };
        // 完全 charge Q = Σ n_i = (3·I − Z₁ − Z₂ − Z₃)/2
        let mut q_full = zeros(d);
        for (i, s) in ["ZII", "IZI", "IIZ"].iter().enumerate() {
            let z = pauli_dense(s);
            for (qq, zz) in q_full.iter_mut().zip(&z) {
                *qq = *qq - zz.scale(0.5);
            }
            let _ = i;
        }
        for i in 0..d {
            q_full[i * d + i] = q_full[i * d + i] + C64::new(1.5, 0.0);
        }
        // 部分 charge (モード 1 のみ) — K が縮退 (4 本の零特異値)
        let mut q_part = zeros(d);
        {
            let z = pauli_dense("ZII");
            for (qq, zz) in q_part.iter_mut().zip(&z) {
                *qq = *qq - zz.scale(0.5);
            }
            for i in 0..d {
                q_part[i * d + i] = q_part[i * d + i] + C64::new(0.5, 0.0);
            }
        }
        let sigma_of = |kmat: &[f64]| -> (f64, f64) {
            // KᵀK の固有値 → (σ_min, σ_max)
            let mut g = vec![0.0f64; m * m];
            for i in 0..m {
                for j in 0..m {
                    let mut s = 0.0;
                    for l in 0..m {
                        s += kmat[l * m + i] * kmat[l * m + j];
                    }
                    g[i * m + j] = s;
                }
            }
            let ev = eig_sym(&g, m);
            (ev[0].max(0.0).sqrt(), ev[m - 1].max(0.0).sqrt())
        };
        let k_full = extract_k(&q_full);
        let k_part = extract_k(&q_part);
        let (s_full, _) = sigma_of(&k_full);
        let (s_part, s_part_max) = sigma_of(&k_part);
        // 区間: 成分ごと半幅 w ⇒ ‖ΔK‖_F ≤ w·m (m² 成分の box)
        let w = 0.02;
        let rf = w * m as f64;
        let gap_full = interval_vs_bar(s_full - rf, s_full + rf, 0.0);
        let gap_part = interval_vs_bar(s_part - rf, s_part + rf, 0.0);
        check(
            "[R3a] spectral gap: 完全 charge は σ_min − ‖ΔK‖_F > 0 (zero crossing なし)・部分 charge は 0 を跨ぐ",
            gap_full == IntervalVerdict::Above
                && gap_part == IntervalVerdict::Straddled
                && s_part < 1e-9
                && s_part_max > 0.9,
            format!(
                "σ_min(full) = {:.4} (下界 {:.4} > 0) / σ(part) ∈ [{:.1e}, {:.3}] — 縮退で構成拒否",
                s_full,
                s_full - rf,
                s_part,
                s_part_max
            ),
        );
        // gap が立つときだけ J = K(−K²)^{−1/2} を構成し J² = −I を機械検証
        let j = {
            let mut g = vec![0.0f64; m * m];
            for i in 0..m {
                for j2 in 0..m {
                    let mut s = 0.0;
                    for l in 0..m {
                        s += k_full[l * m + i] * k_full[l * m + j2];
                    }
                    g[i * m + j2] = s;
                }
            }
            let (ev, vecs) = jacobi_eigh(&g, m);
            // S^{-1/2} = V diag(1/√λ) Vᵀ
            let mut sinv = vec![0.0f64; m * m];
            for i in 0..m {
                for j2 in 0..m {
                    let mut s = 0.0;
                    for l in 0..m {
                        s += vecs[i * m + l] * vecs[j2 * m + l] / ev[l].max(1e-300).sqrt();
                    }
                    sinv[i * m + j2] = s;
                }
            }
            let mut jm = vec![0.0f64; m * m];
            for i in 0..m {
                for j2 in 0..m {
                    let mut s = 0.0;
                    for l in 0..m {
                        s += k_full[i * m + l] * sinv[l * m + j2];
                    }
                    jm[i * m + j2] = s;
                }
            }
            jm
        };
        let mut j2_dev = 0.0f64;
        let mut antisym_dev = 0.0f64;
        for i in 0..m {
            for j2 in 0..m {
                let mut s = 0.0;
                for l in 0..m {
                    s += j[i * m + l] * j[l * m + j2];
                }
                let want = if i == j2 { -1.0 } else { 0.0 };
                j2_dev = j2_dev.max((s - want).abs());
                antisym_dev = antisym_dev.max((j[i * m + j2] + j[j2 * m + i]).abs());
            }
        }
        check(
            "[R3b] gap 資格下の構成: J = K(−K²)^{−1/2} — 実反対称・J² = −I ≤ 1e-12",
            j2_dev < 1e-12 && antisym_dev < 1e-12,
            format!("‖J²+I‖∞ = {:.1e}, ‖J+Jᵀ‖∞ = {:.1e}", j2_dev, antisym_dev),
        );
    }

    // ---------------- [R4] interval cost と set-valued profile ----------------
    println!("\n[R4] resource cost — interval cost の 3 値採用と set-valued profile");
    {
        // 2 qubit: 単発操作 cost 1.0 ± 0.1・entangler cost 2.0 ± 0.3 (amp 軸のみ)
        #[derive(PartialEq, Clone, Copy, Debug)]
        enum Adm {
            Certain,
            Excluded,
            Uncertain,
        }
        let admit = |cost: (f64, f64), b: f64| -> Adm {
            if cost.1 <= b {
                Adm::Certain
            } else if cost.0 > b {
                Adm::Excluded
            } else {
                Adm::Uncertain
            }
        };
        let single = (0.9f64, 1.1f64);
        let ent = (1.7f64, 2.3f64);
        let grid = [0.5f64, 1.2, 1.5, 1.9, 2.4, 3.0];
        // 読み: 単発のみ → [2,2] / 単発 + entangler → [4] / 何もなし → 資源不足
        let mut readings: Vec<&str> = Vec::new();
        for &b in &grid {
            let s = admit(single, b);
            let e = admit(ent, b);
            readings.push(match (s, e) {
                (Adm::Excluded, Adm::Excluded) => "no_accessible",
                (Adm::Certain, Adm::Excluded) => "[2,2]",
                (Adm::Certain, Adm::Certain) => "[8→4]", // 併合読み
                _ => "straddled",
            });
        }
        let want = [
            "no_accessible",
            "[2,2]",
            "[2,2]",
            "straddled",
            "[8→4]",
            "[8→4]",
        ];
        check(
            "[R4a] grid 6 点の 3 値採用: 資源不足 → [2,2] (chain 2 = stable) → 跨ぎ Straddled → 併合 (chain 2 = stable)",
            readings == want,
            format!("readings = {:?}", readings),
        );
        // 中点潰し (点推定) の負制御: b = 1.9 で entangler 中点 2.0 > 1.9 → 排除と
        // 誤読 → [2,2] を「確定」— robust は Straddled が正答 (禁止変換 18/22 合流)
        let midpoint_reading = if (ent.0 + ent.1) / 2.0 <= 1.9 { "[8→4]" } else { "[2,2]" };
        check(
            "[R4b] 中点潰しの負制御: 点推定は b = 1.9 で確定読みを返す (robust は Straddled)",
            midpoint_reading == "[2,2]" && readings[3] == "straddled",
            format!("点推定 = {} / robust = {}", midpoint_reading, readings[3]),
        );
        // 単調 (lax) 性: 確実読みは budget 鎖に沿って粗くなる ([2,2] → [8→4])
        let ok_mono = readings[1] == "[2,2]" && readings[4] == "[8→4]";
        check(
            "[R4c] set-valued profile の単調性 (定義の器械化 — 関手性定理は次期)",
            ok_mono,
            "確実読みの鎖: 不足 ≤ [2,2] ≤ 併合 (跨ぎ点は class から除外して記録)".into(),
        );
    }

    // ---------------- [R5] structured 一致 (同じ区間意味論) ----------------
    println!("\n[R5] dense / Pauli GF(2) — 同じ confidence-region 意味論で裁定一致");
    {
        let n_qubits = 3usize;
        let d = 8usize;
        let ops_str = ["XII", "ZII", "IXI", "IZI", "IIX", "IIZ"];
        let dense_ops: Vec<Vec<C64>> = ops_str.iter().map(|s| pauli_dense(s)).collect();
        let pauli_ops: Vec<PauliVector> = ops_str.iter().map(|s| PauliVector::from_str(s)).collect();
        // 真の統計量 s ∈ {0,1}: dense = ‖[a,b]‖/‖2ab‖ 正規化, Pauli = 反可換 bit
        let k_pairs: Vec<(usize, usize)> = (0..6)
            .flat_map(|i| ((i + 1)..6).map(move |j| (i, j)))
            .collect();
        let mut agree = true;
        let mut verdicts_dense: Vec<IntervalVerdict> = Vec::new();
        // N = 5000 では Bonferroni α/15 の CP 上限 (≈1.3e-3) が τ = 1e-3 を跨ぐ —
        // 非辺の certify には 8000 ショットが要る (分解能は観測量の関数, v343 [F2b])
        let (n_shots, alpha, tau) = (8000usize, 0.05f64, 1e-3f64);
        let m_stats = k_pairs.len();
        for &(i, j) in &k_pairs {
            let c = commutator(&dense_ops[i], &dense_ops[j], d);
            let s_dense: f64 = if hs_norm(&c) > 1.0 { 1.0 } else { 0.0 };
            let s_pauli: f64 = if pauli_ops[i].anticommutes(&pauli_ops[j]) {
                1.0
            } else {
                0.0
            };
            if (s_dense - s_pauli).abs() > 1e-12 {
                agree = false;
            }
            // 代表カウント (決定的 synthetic 記録) → CP 区間 (Bonferroni α/m)
            let k_obs = (n_shots as f64 * s_dense).round() as usize;
            let (lo, hi) = cp_interval(k_obs, n_shots, alpha / m_stats as f64);
            verdicts_dense.push(interval_vs_bar(lo, hi, tau));
        }
        check(
            "[R5a] 統計量の lane 一致: dense ‖[·,·]‖ 正規化 = Pauli 反可換 bit (全 15 対)",
            agree,
            "同じ真値 → 同じ記録 → 同じ区間 → 同じ graph 裁定 (構成から)".into(),
        );
        // site 族: 全 site 内対が Above・site 間対が Below (N = 5000 で判定可能)
        let mut ok_site = true;
        for (idx, &(i, j)) in k_pairs.iter().enumerate() {
            let same_site = i / 2 == j / 2;
            let want = if same_site {
                IntervalVerdict::Above
            } else {
                IntervalVerdict::Below
            };
            if verdicts_dense[idx] != want {
                ok_site = false;
            }
        }
        check(
            "[R5b] site 族 (N = 8000, α/15): site 内 = edge・site 間 = non-edge が全対で確定",
            ok_site,
            format!("15 対の graph 裁定が確定 (τ = {} — N = 5000 では上限 1.3e-3 が跨ぐ)", tau),
        );
        // 不足ショット: N = 100 では k = 0 でも hi > τ — 両 lane とも Straddled
        let (_, hi100) = cp_interval(0, 100, alpha / m_stats as f64);
        check(
            "[R5c] N = 100 では非可換 0 でも certify 不能 (hi > τ) — 両 lane とも Straddled が正答",
            interval_vs_bar(0.0, hi100, tau) == IntervalVerdict::Straddled,
            format!("hi(k=0, N=100) = {:.4} > τ = {}", hi100, tau),
        );
        // 対応原理 (確定 graph 上): dense 閉包/中心 vs 2^{dim V}/2^{radical}
        // (a) 1 site 成分 {X₁, Z₁}: dim V = 2 → 閉包 4・radical 0 → 中心 1
        let comp_gens = vec![dense_ops[0].clone(), dense_ops[1].clone()];
        let closure = algebra_closure(&comp_gens, d);
        let center = closure_center_basis(&closure, &comp_gens, d);
        // (b) 部分 address {X₁, Z₁, Z₂}: radical {Z₂} → 中心 2 (dense: span{I, Z₂})
        let pa_gens = vec![
            dense_ops[0].clone(),
            dense_ops[1].clone(),
            dense_ops[3].clone(),
        ];
        let pa_closure = algebra_closure(&pa_gens, d);
        let pa_center = closure_center_basis(&pa_closure, &pa_gens, d);
        // GF(2) 側: rank と radical
        use uft_sim::structured_backend::{gf2_rank, Gf2Vec};
        let gf2_of = |idxs: &[usize]| -> (usize, usize) {
            let rows: Vec<Gf2Vec> = idxs.iter().map(|&i| pauli_ops[i].combined()).collect();
            let dim_v = gf2_rank(&rows);
            // radical: V 内で全行と symplectic 直交 — 全 2^|idxs| 組合せの XOR を走査
            // (小規模 — 正確): span の全元を列挙して ω 直交をカウント
            let mut span: Vec<Gf2Vec> = vec![Gf2Vec::zeros(2 * n_qubits)];
            for r in &rows {
                let snapshot = span.clone();
                for s in snapshot {
                    let mut x = s.clone();
                    x.xor_assign(r);
                    if !span.iter().any(|t| {
                        let mut y = t.clone();
                        y.xor_assign(&x);
                        y.is_zero()
                    }) {
                        span.push(x);
                    }
                }
            }
            let mut rad_elems = 0usize;
            for v in &span {
                if v.is_zero() {
                    continue;
                }
                // symplectic form: PauliVector::anticommutes 相当を Gf2 で
                let ortho = rows.iter().all(|r| {
                    // ω(v, r) = x_v·z_r + z_v·x_r — Gf2Vec::dot は combined 上の
                    // symplectic 実装を backend に合わせる: ここでは anticommute bit を
                    // ビット直積で再現するために PauliVector 経由は使えないので
                    // combined の前半/後半を突き合わせる
                    let nb = n_qubits;
                    let mut acc = false;
                    for b in 0..nb {
                        let (xv, zv) = (v.get(b), v.get(nb + b));
                        let (xr, zr) = (r.get(b), r.get(nb + b));
                        if xv && zr {
                            acc = !acc;
                        }
                        if zv && xr {
                            acc = !acc;
                        }
                    }
                    !acc
                });
                if ortho {
                    rad_elems += 1;
                }
            }
            // radical の GF(2) 次元: 元の個数 = 2^r − 1 (0 を除く)
            let r = (rad_elems + 1).trailing_zeros() as usize;
            (dim_v, r)
        };
        let (v1, r1) = gf2_of(&[0, 1]);
        let (v2, r2) = gf2_of(&[0, 1, 3]);
        let ok_corr = closure.len() == (1usize << v1)
            && center.len() == (1usize << r1)
            && pa_closure.len() == (1usize << v2)
            && pa_center.len() == (1usize << r2);
        check(
            "[R5d] 対応原理: dense 閉包/中心 = 2^{dim V}/2^{radical} (site 成分と部分 address)",
            ok_corr,
            format!(
                "site: 閉包 {} = 2^{}, 中心 {} = 2^{} / 部分 address: 閉包 {} = 2^{}, 中心 {} = 2^{}",
                closure.len(),
                v1,
                center.len(),
                r1,
                pa_closure.len(),
                v2,
                pa_center.len(),
                r2
            ),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "exact reader が同時信頼集合の上に持ち上がった — σ_min の同時下界・worst-case\n       cross-talk・orbit 距離の区間・spectral-gap 証明書つき J・interval cost の\n       3 値採用・dense/structured の同一区間意味論。裁定は集合の全域で一致する\n       ときだけ返り、跨ぎは Straddled が正答である。synthetic lane であり実測\n       ノイズではない (v34.6 の real-data lane と区別)。"
        } else {
            "**robust atlas の破れ** — 下界・区間・lane 一致を修復せよ"
        }
    );
    println!(
        "\n総合判定: {}",
        if nfail == 0 { "[PASS]" } else { "[FAIL]" }
    );
    if nfail > 0 {
        std::process::exit(1);
    }
}
