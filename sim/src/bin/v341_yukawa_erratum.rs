//! v34.1 Yukawa erratum — 素数 Pfaffian 推論の反証の正式受理 (PROMPT/15 §1)
//!
//! 背景: v13.2 (QRN-YUK-015 発見 2) と論文 geometric-yukawa は「Pf(F) = 3 が素数
//! であることが 2 つの磁気固有面 (skew singular value f±) を不均等にし、一方を
//! 平坦化する」という普遍推論を含んでいた。FollowUp (cross-model clean-room,
//! reports/final_report.md §10 / docs/derivations.md §12) は整数反対称行列
//!   F = [[0,1,1,1],[-1,0,1,-1],[-1,-1,0,1],[-1,1,-1,0]]
//! が Pf(F) = 3 かつ FᵀF = 3I₄ を満たすことを示した — 両 skew singular value は
//! ともに √3 で等しく、素数性は積 f₊f₋ = |Pf| を固定するが不等性を強制しない。
//!
//! 数学的事実 (4×4 実反対称 F): f₊² + f₋² = Σ_{i<j} F_ij² =: Σ、f₊f₋ = |Pf(F)|。
//! よって (f₊ − f₋)² = Σ − 2|Pf| で、|Pf| = 3 のとき等スケール ⟺ Σ = 6。
//! Pfaffian は unimodular 基底変換 (det S = 1) で不変な代数・位相量、skew scale は
//! 計量 (どの基底で FᵀF を組むか) に依存する幾何量 — 前者から後者への昇格には
//! metric/lattice compatibility bridge が必要である。
//!
//! 検証 (全て [PASS]/[FAIL] 内蔵):
//!   [E1] 反例の厳密受理 — Pf = 3・FᵀF = 3I₄・det F = 9 = Pf² (全て i64 厳密) +
//!        jacobi_eigh で固有値 4 本とも 3 (数値照合)。素朴推論の予測 ({3,1} 型
//!        スケール ⇒ FᵀF 固有値 {9,9,1,1}) との定量的乖離。
//!   [E2] 位相 ↛ 計量の両方向 — unimodular S (det 1) で F' = SᵀFS は Pf 不変 3 の
//!        まま Σ' = 8 → f₊ − f₋ = √2 (不等)。同じ Pf = 3 で等/不等の両方が実現。
//!   [E3] 族内定理 (窓なし) — v13.2 の走査族 (座標 2-平面 + 傾き対, Pf = Q₁Q₂+ts,
//!        f13 = f24 = 0) では a²+b²+c²+d² = 6 ∧ ab+cd = ±3 の整数解が存在しない
//!        (Σ=6 ⇒ |成分| ≤ 2 の有限全数で完全)。Σ ≥ 2|Pf| = 6 と合わせ、族内では
//!        Pf = ±3 ⇒ (f₊−f₋)² = Σ−6 ≥ 1 — 最小ギャップ (2,1,1,1) 型は f₊−f₋ = 1
//!        厳密。v13.2 の 7 走査点の gap² も列挙 (全て ≥ 4)。反例は 6 成分全て
//!        非零で族外 — 族内観測は正しく、その普遍化だけが誤りだった。
//!   [E4] 台帳・論文の正誤表アンカー — claims.yml の分割 (QRN-YUK-015 正誤表 +
//!        QRN-YUK-034 refuted_as_stated)・META-013 Erratum・論文 2 本の撤回文の
//!        blockquote 化 (裸の主張として残存しない)・replications.yml の REP-001
//!        (cross_model_clean_room / external_replications = 0 維持) を機械検査。

use std::fs;
use std::path::Path;
use uft_sim::{jacobi_eigh, self_test};

type M4 = [[i64; 4]; 4];

/// FollowUp docs/derivations.md §12 の反例行列 (逐語転記)
const F_CEX: M4 = [[0, 1, 1, 1], [-1, 0, 1, -1], [-1, -1, 0, 1], [-1, 1, -1, 0]];

/// v13.2 の走査 7 点 (sim/src/bin/v132_deeptilt.rs の family 逐語転記)
const V132_FAMILY: [[i64; 4]; 7] = [
    [3, 1, 0, 0],
    [2, 2, 1, -1],
    [4, 1, 1, -1],
    [3, 3, 1, -6],
    [3, 3, 2, -3],
    [3, 3, 3, -2],
    [3, 3, 6, -1],
];

/// 4×4 反対称行列の Pfaffian (上三角 6 成分の閉形式, i64 厳密)
fn pf4(f: &M4) -> i64 {
    f[0][1] * f[2][3] - f[0][2] * f[1][3] + f[0][3] * f[1][2]
}

/// 反対称性の検査 (対角 0・F = −Fᵀ)
fn is_antisym(f: &M4) -> bool {
    (0..4).all(|i| (0..4).all(|j| f[i][j] == -f[j][i]))
}

/// 行列積 (i64 厳密)
fn matmul_i(a: &M4, b: &M4) -> M4 {
    let mut c = [[0i64; 4]; 4];
    for (i, row) in c.iter_mut().enumerate() {
        for (j, cij) in row.iter_mut().enumerate() {
            *cij = (0..4).map(|k| a[i][k] * b[k][j]).sum();
        }
    }
    c
}

/// 転置
fn transpose(a: &M4) -> M4 {
    let mut t = [[0i64; 4]; 4];
    for (i, row) in a.iter().enumerate() {
        for (j, &v) in row.iter().enumerate() {
            t[j][i] = v;
        }
    }
    t
}

/// det (i64, 第 1 行の余因子展開)
fn det4(m: &M4) -> i64 {
    fn det3(m: &[[i64; 3]; 3]) -> i64 {
        m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
            - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
            + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
    }
    let mut d = 0i64;
    for j in 0..4 {
        let mut sub = [[0i64; 3]; 3];
        for (si, i) in (1..4).enumerate() {
            let mut sj = 0;
            for jj in 0..4 {
                if jj == j {
                    continue;
                }
                sub[si][sj] = m[i][jj];
                sj += 1;
            }
        }
        let sgn = if j % 2 == 0 { 1 } else { -1 };
        d += sgn * m[0][j] * det3(&sub);
    }
    d
}

/// 上三角成分の二乗和 Σ_{i<j} F_ij² = f₊² + f₋²
fn sum_sq_upper(f: &M4) -> i64 {
    let mut s = 0;
    for i in 0..4 {
        for j in (i + 1)..4 {
            s += f[i][j] * f[i][j];
        }
    }
    s
}

/// 走査族 (Q₁,Q₂,t,s) ↦ 磁束行列 — f12 = Q₁ (12-平面), f34 = Q₂ (34-平面),
/// f14 = t, f23 = s (傾き対), f13 = f24 = 0。Pf = Q₁Q₂ + ts (v13.2 の式)。
fn family_matrix(q: &[i64; 4]) -> M4 {
    [
        [0, q[0], 0, q[2]],
        [-q[0], 0, q[3], 0],
        [0, -q[3], 0, q[1]],
        [-q[2], 0, -q[1], 0],
    ]
}

/// FᵀF の固有値 (jacobi_eigh, 昇順ソート済み)
fn ftf_eigs(f: &M4) -> Vec<f64> {
    let g = matmul_i(&transpose(f), f);
    let mut a = vec![0.0f64; 16];
    for i in 0..4 {
        for j in 0..4 {
            a[i * 4 + j] = g[i][j] as f64;
        }
    }
    let (mut ev, _) = jacobi_eigh(&a, 4);
    ev.sort_by(|x, y| x.partial_cmp(y).unwrap());
    ev
}

fn main() {
    self_test();
    println!("=== v34.1 Yukawa erratum — 素数 Pfaffian 推論の反証の正式受理 (PROMPT/15 §1) ===\n");
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

    // ---------------- [E1] 反例の厳密受理 ----------------
    println!("[E1] FollowUp 反例の整数再検算 (受理の数学的決着)");
    {
        let f = &F_CEX;
        let anti = is_antisym(f);
        let pf = pf4(f);
        let g = matmul_i(&transpose(f), f);
        let mut is_3i = true;
        for (i, row) in g.iter().enumerate() {
            for (j, &v) in row.iter().enumerate() {
                let want = if i == j { 3 } else { 0 };
                if v != want {
                    is_3i = false;
                }
            }
        }
        let det = det4(f);
        let sig = sum_sq_upper(f);
        check(
            "[E1a] 反対称・Pf(F) = 3 (i64 厳密)",
            anti && pf == 3,
            format!("Pf = {} (af−be+cd = 1·1 − 1·(−1) + 1·1)", pf),
        );
        check(
            "[E1b] FᵀF = 3I₄ (i64 厳密) ⇒ 両 skew singular value = √3 (等スケール)",
            is_3i,
            format!("Σ_{{i<j}}F² = {} = 2|Pf| (等号 ⟺ f₊ = f₋)", sig),
        );
        check(
            "[E1c] det F = 9 = Pf² (整合)",
            det == 9 && det == pf * pf,
            format!("det = {}", det),
        );
        let ev = ftf_eigs(f);
        let dev = ev.iter().map(|x| (x - 3.0).abs()).fold(0.0, f64::max);
        check(
            "[E1d] jacobi_eigh: FᵀF の固有値 4 本とも 3 (数値照合)",
            dev < 1e-12,
            format!("max|λ−3| = {:.2e}", dev),
        );
        // 素朴推論 (「素数 ⇒ 整数分解 {3,1} 型で不均等」) の予測との乖離:
        // {f₊,f₋} = {3,1} なら FᵀF 固有値は {9,9,1,1}。実際は {3,3,3,3}。
        let naive_dev = (ev[3] - 9.0).abs().min((ev[0] - 1.0).abs());
        check(
            "[E1e] 素朴推論の予測 ({9,9,1,1}) と厳密値 ({3,3,3,3}) の乖離",
            naive_dev > 1.9,
            format!(
                "予測との最小距離 {:.3} (等スケール √3 は整数スケールのどの分解でもない)",
                naive_dev
            ),
        );
    }

    // ---------------- [E2] 位相 ↛ 計量 (両方向) ----------------
    println!("\n[E2] Pfaffian (位相・代数) と skew scale (計量) の分離 — unimodular 変換");
    {
        // S = I + E₁₂ (剪断, det = 1): 同じ整数交代形式の別の格子基底表示
        let s: M4 = [[1, 1, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]];
        let det_s = det4(&s);
        let fp = matmul_i(&matmul_i(&transpose(&s), &F_CEX), &s);
        let pf_p = pf4(&fp);
        let sig_p = sum_sq_upper(&fp);
        check(
            "[E2a] det S = 1 (unimodular) で Pf(SᵀFS) = det(S)·Pf(F) = 3 (不変)",
            det_s == 1 && pf_p == 3 && is_antisym(&fp),
            format!("Pf(F') = {}", pf_p),
        );
        let gap2 = sig_p - 2 * pf_p.abs();
        check(
            "[E2b] Σ' = 8 ⇒ (f₊−f₋)² = 2 — 同じ Pf = 3 で不等スケールに反転",
            sig_p == 8 && gap2 == 2,
            format!("Σ' = {}, gap² = {} (F は 6/0 だった)", sig_p, gap2),
        );
        let ev = ftf_eigs(&fp);
        // f±² = 4 ± √7
        let want_hi = 4.0 + 7.0f64.sqrt();
        let want_lo = 4.0 - 7.0f64.sqrt();
        let dev = (ev[3] - want_hi).abs().max((ev[0] - want_lo).abs());
        check(
            "[E2c] jacobi_eigh: F'ᵀF' 固有値 = 4 ± √7 (各 2 重) — 数値照合",
            dev < 1e-12,
            format!(
                "f₊ = {:.6}, f₋ = {:.6}, 積 = {:.12} (= |Pf| = 3)",
                ev[3].sqrt(),
                ev[0].sqrt(),
                (ev[3] * ev[0]).sqrt()
            ),
        );
        // 結論: Pf = 3 は等スケール (F) とも不等スケール (F') とも両立 —
        // Pfaffian 単独では skew 異方性をどちらの向きにも決められない。
    }

    // ---------------- [E3] 族内定理 (窓なし) ----------------
    println!("\n[E3] v13.2 走査族の内側では不等が整数論的に強制される (窓なしの定理)");
    {
        // (a) Σ = a²+b²+c²+d² = 6 ∧ ab+cd = ±3 の整数解は存在しない。
        //     Σ = 6 ⇒ |成分| ≤ 2 なので範囲 [-2,2] の全数で完全 (窓仮定なし)。
        let mut n_sum6 = 0usize;
        let mut n_equal_scale = 0usize;
        for a in -2i64..=2 {
            for b in -2i64..=2 {
                for c in -2i64..=2 {
                    for d in -2i64..=2 {
                        if a * a + b * b + c * c + d * d == 6 {
                            n_sum6 += 1;
                            if (a * b + c * d).abs() == 3 {
                                n_equal_scale += 1;
                            }
                        }
                    }
                }
            }
        }
        check(
            "[E3a] 族内の等スケール解 (Σ=6 ∧ |Pf|=3) は存在しない — 全数列挙",
            n_sum6 > 0 && n_equal_scale == 0,
            format!("Σ=6 の整数点 {} 個を走査、|ab+cd| = 3 は 0 個", n_sum6),
        );
        // ⇒ 族内で |Pf| = 3 なら Σ ≥ 2|Pf| = 6 (AM-GM) かつ Σ ≠ 6 (上の全数) かつ
        //    Σ ∈ ℤ なので Σ ≥ 7、つまり (f₊−f₋)² = Σ − 6 ≥ 1。範囲の上限は不要。

        // (b) 最小ギャップの実現点: |成分| ≤ 6 で Pf = 3 の族点を列挙し最小 Σ を確認
        let mut min_sig = i64::MAX;
        let mut n_pf3 = 0usize;
        let mut min_pts: Vec<[i64; 4]> = Vec::new();
        for a in -6i64..=6 {
            for b in -6i64..=6 {
                for c in -6i64..=6 {
                    for d in -6i64..=6 {
                        if a * b + c * d == 3 {
                            n_pf3 += 1;
                            let sig = a * a + b * b + c * c + d * d;
                            if sig < min_sig {
                                min_sig = sig;
                                min_pts.clear();
                            }
                            if sig == min_sig {
                                min_pts.push([a, b, c, d]);
                            }
                        }
                    }
                }
            }
        }
        check(
            "[E3b] 族内 (|成分| ≤ 6) の最小 Σ = 7 ⇒ 最小ギャップ f₊−f₋ = 1 厳密",
            min_sig == 7,
            format!(
                "Pf=3 の族点 {} 個中、最小 Σ = {} の実現 {} 個 (例 {:?} — f₊f₋=3, f₊−f₋=1 ⇒ f± = (√13±1)/2)",
                n_pf3,
                min_sig,
                min_pts.len(),
                min_pts.first().unwrap()
            ),
        );
        // 数値照合: (2,1,1,1) の f₊ − f₋ = 1
        let fm = family_matrix(&[2, 1, 1, 1]);
        let ev = ftf_eigs(&fm);
        let gap = ev[3].sqrt() - ev[0].sqrt();
        check(
            "[E3c] (2,1,1,1) の jacobi 照合: f₊ − f₋ = 1.000000000000",
            pf4(&fm) == 3 && (gap - 1.0).abs() < 1e-12,
            format!("gap = {:.12}", gap),
        );
        // (c) v13.2 の 7 走査点: 全て族内 (f13 = f24 = 0)・Pf = 3・gap² ≥ 4
        let mut all_ok = true;
        let mut gaps = Vec::new();
        for q in &V132_FAMILY {
            let fm = family_matrix(q);
            let pf = pf4(&fm);
            let sig = sum_sq_upper(&fm);
            let gap2 = sig - 6;
            if pf != 3 || gap2 < 4 {
                all_ok = false;
            }
            gaps.push(gap2);
        }
        check(
            "[E3d] v13.2 の 7 走査点は全て Pf = 3・gap² = Σ−6 ∈ {4,13,25,49}",
            all_ok,
            format!("gap² = {:?}", gaps),
        );
        // (d) 反例は族外 (f13 = f24 = 0 が破れる) — 「族内観測は正しく普遍化が誤り」
        let outside = F_CEX[0][2] != 0 || F_CEX[1][3] != 0;
        check(
            "[E3e] 反例 F は走査族の外 (f13 = f24 = 0 が破れ、6 成分全て非零)",
            outside,
            format!("f13 = {}, f24 = {}", F_CEX[0][2], F_CEX[1][3]),
        );
    }

    // ---------------- [E4] 台帳・論文の正誤表アンカー ----------------
    println!("\n[E4] 台帳・論文の正誤表アンカー (分割・撤回・受理の機械検査)");
    {
        let root = if Path::new("claims.yml").exists() {
            "."
        } else {
            ".."
        };
        let rd = |p: &str| fs::read_to_string(format!("{}/{}", root, p)).unwrap_or_default();
        let claims = rd("claims.yml");
        // QRN-YUK-034 (refuted_as_stated の正式記録) と QRN-YUK-015 の分割
        let yuk034 = claims.contains("- id: QRN-YUK-034") && claims.contains("refuted_as_stated");
        let yuk015 =
            claims.contains("正誤表 (v34.1)") && claims.contains("走査した admissible family");
        let meta013 = claims.contains("Erratum (v34.1)");
        check(
            "[E4a] claims.yml — QRN-YUK-034 登録・QRN-YUK-015 の走査族限定への分割・META-013 の Erratum",
            yuk034 && yuk015 && meta013,
            format!("034 {} / 015 正誤表 {} / META-013 {}", yuk034, yuk015, meta013),
        );
        // 論文 2 本: 撤回文が裸の主張として残存しない (引用は blockquote のみ許す)
        let mut paper_ok = true;
        let mut detail = String::new();
        for p in [
            "paper/geometric-yukawa-full.md",
            "paper/geometric-yukawa.md",
        ] {
            let text = rd(p);
            let bare: Vec<&str> = text
                .lines()
                .filter(|l| l.contains("being prime, forces") && !l.trim_start().starts_with('>'))
                .collect();
            let has_erratum = text.contains("Erratum (v34.1)")
                && text.contains("does not imply unequal magnetic skew scales");
            if !bare.is_empty() || !has_erratum {
                paper_ok = false;
            }
            detail.push_str(&format!(
                "{}: 裸の残存 {} / erratum {}; ",
                p,
                bare.len(),
                has_erratum
            ));
        }
        check(
            "[E4b] 論文 2 本 — 撤回文は blockquote 引用のみ・Erratum (v34.1) 節あり",
            paper_ok,
            detail,
        );
        // replications.yml — FollowUp の受理 (external には数えない)
        let rep = rd("replications.yml");
        let rep_ok = rep.contains("- id: REP-001")
            && rep.contains("replication_kind: cross_model_clean_room")
            && rep.contains("verdict: partially_replicated")
            && rep.contains("counts_as_external_replication: false")
            && rep.contains("different_author: false")
            && rep.lines().any(|l| l.trim() == "external_replications: 0");
        check(
            "[E4c] replications.yml — REP-001 (cross_model_clean_room) 登録・external_replications = 0 維持",
            rep_ok,
            "多次元独立性 profile での登録 (一次元の external カウントは動かさない)".into(),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "反例は整数演算で受理された — 素数性は積を固定するが不等性を強制しない。\n       族内不等 (f₊−f₋)² ≥ 1 は生き残り、その普遍化だけが撤回された。\n       位相・代数的不変量から計量依存量への昇格には bridge が必要である。"
        } else {
            "**erratum が不完全** — 反例の受理・台帳の分割・論文の撤回を完了せよ"
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
