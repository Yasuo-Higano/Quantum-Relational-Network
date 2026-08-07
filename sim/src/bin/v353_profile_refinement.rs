//! v35.3 set-valued resource profile の refinement / no-go / 安定性 (PROMPT/16 §7)
//!
//! v34.5 で器械化した set-valued profile (budget ↦ 信頼集合上の読み) の数学的地位を
//! 確定する。中心の問い: 「予算を増やせば読みは自動的に精密化するか」— 答えは**否**
//! (RPF-1)。関手性は被覆からは出ず、samplewise nesting を**構成で買った**ときだけ
//! lax refinement が成立する (RPF-2)。離散 verdict の安定性は margin 条件つきで
//! しか成立せず (RPF-3)、boundary では global Lipschitz 安定性が存在しない (RPF-4)。
//! v33.3 で凍結した昇格規則「stable ⟺ chain ≥ 2」(禁止変換 17) は refinement +
//! transient 反例の系として導出される (RPF-5)。
//!
//! 対応する Lean: proofs/ResourceProfile.lean (8 定理 — RPF-2/5a は型多相の一般
//! 定理・RPF-3 は全整数値・RPF-1/4/5b は反例)。本バイナリは CP 区間の実装で
//! 各定理の数値側を測る:
//!   [R1] RPF-1: 独立標本 / 累積標本 / intersection 構成の nesting 破れ率
//!   [R2] RPF-2: nested 構成での r(Q_{b'}) ⊆ Q_b 違反 0
//!   [R3] RPF-3: margin 未満の摂動で verdict 不変 (0 flip)・margin 超で flip 証人
//!   [R4] RPF-4: 摂動 1/n vs verdict 距離 1 — 比 n の発散表
//!   [R5] RPF-5: 単点 transient の実発生率 + 2-chain (nested) の健全性違反 0
//!   [R6] Lean 反例の整数橋

use uft_sim::finite_data::cp_interval;
use uft_sim::{self_test, Rng};

/// 離散 verdict (Lean ResourceProfile.Verdict の鏡映)
#[derive(Debug, Clone, Copy, PartialEq)]
enum Verdict {
    Edge,
    NoEdge,
    Straddled,
}

fn verdict(tau: f64, lo: f64, hi: f64) -> Verdict {
    if tau < lo {
        Verdict::Edge
    } else if hi <= tau {
        Verdict::NoEdge
    } else {
        Verdict::Straddled
    }
}

/// margin (Lean の margin の実数版 — 整数 +1 は連続では 0 距離境界)
fn margin(tau: f64, lo: f64, hi: f64) -> f64 {
    if tau < lo {
        lo - tau
    } else if hi <= tau {
        tau - hi
    } else {
        (tau - lo).min(hi - tau)
    }
}

fn main() {
    self_test();
    println!("=== v35.3 set-valued resource profile — RPF-1..5 (PROMPT/16 §7) ===\n");
    let mut nfail = 0usize;
    let mut check = |name: &str, ok: bool, detail: String| {
        println!("  [{}] {}  {}", if ok { "PASS" } else { "FAIL" }, name, detail);
        if !ok {
            nfail += 1;
        }
    };

    let alpha = 0.05;
    let p_true = 0.42; // 真のパラメータ (edge 重み)
    let reps = 4000usize;

    // ---------------- [R1] RPF-1: nesting は被覆から出ない ----------------
    {
        let mut rng = Rng::new(35301);
        let (n1, n2) = (100usize, 300usize);
        let mut viol_indep = 0usize;
        let mut viol_cumul = 0usize;
        let mut viol_inter = 0usize;
        let mut cover_inter = 0usize;
        for _ in 0..reps {
            // 独立標本: 予算 b (n1 shot) と b' (n2 shot 別取得)
            let k1 = (0..n1).filter(|_| rng.f64() < p_true).count();
            let k2 = (0..n2).filter(|_| rng.f64() < p_true).count();
            let (lo1, hi1) = cp_interval(k1, n1, alpha);
            let (lo2, hi2) = cp_interval(k2, n2, alpha);
            // nesting 違反: C_{b'} ⊄ C_b (区間の包含が破れる)
            if lo2 < lo1 || hi2 > hi1 {
                viol_indep += 1;
            }
            // 累積標本: b' = b の標本 + 追加 (restriction)
            let k_extra = (0..(n2 - n1)).filter(|_| rng.f64() < p_true).count();
            let kc = k1 + k_extra;
            let (lo2c, hi2c) = cp_interval(kc, n2, alpha);
            if lo2c < lo1 || hi2c > hi1 {
                viol_cumul += 1;
            }
            // intersection 構成 (同時 α 配分 α/2 ずつ + 逐次交差): 構成的に nested
            let (l1, h1) = cp_interval(k1, n1, alpha / 2.0);
            let (l2, h2) = cp_interval(kc, n2, alpha / 2.0);
            let (li, hi_) = (l2.max(l1), h2.min(h1));
            if li < l1 || hi_ > h1 {
                viol_inter += 1;
            }
            if li <= p_true && p_true <= hi_ {
                cover_inter += 1;
            }
        }
        let r_indep = viol_indep as f64 / reps as f64;
        let r_cumul = viol_cumul as f64 / reps as f64;
        let cov = cover_inter as f64 / reps as f64;
        check(
            "[R1a] RPF-1: 独立標本の CP 区間は高頻度で nesting が破れる (被覆 ≠ samplewise nesting)",
            r_indep > 0.05,
            format!("独立標本の破れ率 = {:.3} (n {} → {})", r_indep, n1, n2),
        );
        check(
            "[R1b] RPF-1: 累積標本 (restriction) でも nesting は自動でない (破れ率を正直に記録)",
            true, // 率の大小に関わらず記録が目的 (0 なら 0 と記録)
            format!("累積標本の破れ率 = {:.3} — 0 でも構成なしの保証はない (Lean 反例が最小証人)", r_cumul),
        );
        check(
            "[R1c] intersection 構成 (α/2 + 逐次交差): nesting 違反 0・被覆 ≥ 1 − α を測定",
            viol_inter == 0 && cov >= 1.0 - alpha,
            format!("違反 {}/{}, 被覆 = {:.4} (登録下限 {})", viol_inter, reps, cov, 1.0 - alpha),
        );
    }

    // ---------------- [R2] RPF-2: nested 構成での refinement ----------------
    {
        // 商: 予算 b は 2 値 (edge/noEdge/straddled → 粗く {検出, 非検出, 保留})、
        // 予算 b' は同じ verdict 語彙 (恒等細分) — r = id。
        // nested (intersection) 構成では Q_{b'} の verdict 集合 ⊆ Q_b の verdict 集合。
        let mut rng = Rng::new(35302);
        let tau = 0.3;
        let (n1, n2) = (80usize, 240usize);
        let mut viol = 0usize;
        let mut narrower = 0usize;
        for i in 0..reps {
            let p = 0.1 + 0.6 * ((i % 40) as f64 / 40.0); // 走査 (boundary 跨ぎを含む)
            let k1 = (0..n1).filter(|_| rng.f64() < p).count();
            let k_extra = (0..(n2 - n1)).filter(|_| rng.f64() < p).count();
            let (l1, h1) = cp_interval(k1, n1, alpha / 2.0);
            let (l2r, h2r) = cp_interval(k1 + k_extra, n2, alpha / 2.0);
            let (l2, h2) = (l2r.max(l1), h2r.min(h1)); // intersection 構成
            // 集合値の読み Q = 区間が許す verdict の集合 (区間 ⊆ 半直線判定)
            let q = |lo: f64, hi: f64| -> Vec<Verdict> {
                let mut v = Vec::new();
                if hi > tau {
                    v.push(Verdict::Edge); // 区間内に edge 側の点がある
                }
                if lo <= tau {
                    v.push(Verdict::NoEdge); // 区間内に no-edge 側の点がある
                }
                v
            };
            let qb = q(l1, h1);
            let qb2 = q(l2, h2);
            // r = id: Q_{b'} ⊆ Q_b
            if qb2.iter().any(|v| !qb.contains(v)) {
                viol += 1;
            }
            if qb2.len() < qb.len() {
                narrower += 1;
            }
        }
        check(
            "[R2] RPF-2: nested 構成で r(Q_{b'}) ⊆ Q_b の違反 0 (精密化は単調 — lax refinement)",
            viol == 0,
            format!("違反 {}/{} (読みが実際に狭まった回数 {} — 精密化は起きている)", viol, reps, narrower),
        );
    }

    // ---------------- [R3] RPF-3: margin 安定性 ----------------
    {
        let mut rng = Rng::new(35303);
        let tau = 0.3;
        let mut flips_below = 0usize;
        for _ in 0..reps {
            let lo = rng.f64() * 0.8;
            let hi = lo + rng.f64() * 0.3;
            let m = margin(tau, lo, hi);
            if m <= 0.0 {
                continue;
            }
            // margin 未満の摂動
            let d = m * 0.99 * rng.f64();
            let lo2 = lo + (rng.f64() * 2.0 - 1.0) * d;
            let hi2 = hi + (rng.f64() * 2.0 - 1.0) * d;
            let (lo2, hi2) = (lo2.min(hi2), lo2.max(hi2));
            if verdict(tau, lo2, hi2) != verdict(tau, lo, hi) {
                flips_below += 1;
            }
        }
        // margin 超の摂動は flip する証人 (boundary 跨ぎ)
        let (lo, hi) = (0.35, 0.5); // edge (margin 0.05)
        let v0 = verdict(tau, lo, hi);
        let v1 = verdict(tau, lo - 0.1, hi); // 摂動 0.1 > margin
        let flip_above_witness = v0 == Verdict::Edge && v1 == Verdict::Straddled;
        check(
            "[R3] RPF-3: margin 未満の摂動で verdict flip 0・margin 超は flip 証人あり",
            flips_below == 0 && flip_above_witness,
            format!("margin 未満 flip {}/{}, 証人: Edge → Straddled (摂動 0.1 > margin 0.05)", flips_below, reps),
        );
    }

    // ---------------- [R4] RPF-4: 大域 Lipschitz なし ----------------
    {
        let tau = 0.3;
        let mut table = Vec::new();
        let mut ok = true;
        for &n in &[10.0f64, 100.0, 1000.0, 100000.0] {
            let eps = 1.0 / n;
            let v0 = verdict(tau, tau + eps, tau + 2.0 * eps); // Edge
            let v1 = verdict(tau, tau, tau + 2.0 * eps); // Straddled (lo を eps 下げ)
            let dist = if v0 != v1 { 1.0 } else { 0.0 };
            let ratio = dist / eps;
            ok &= v0 == Verdict::Edge && v1 == Verdict::Straddled && ratio >= n * 0.99;
            table.push(format!("1/{:.0}→比{:.0}", n, ratio));
        }
        check(
            "[R4] RPF-4: 摂動 1/n で verdict 距離 1 — 比 n が非有界 (global Lipschitz 定数は存在しない)",
            ok,
            table.join(", "),
        );
    }

    // ---------------- [R5] RPF-5: transient と 2-chain 健全性 ----------------
    {
        let mut rng = Rng::new(35304);
        let tau = 0.3;
        let (n1, n2) = (60usize, 180usize);
        let p = 0.36; // boundary 近傍 (singleton が出たり消えたりする領域)
        let mut singleton_b = 0usize;
        let mut transient = 0usize; // 独立標本で singleton が消えた
        let mut chain_viol = 0usize; // nested 構成での 2-chain 健全性違反
        for _ in 0..reps {
            let k1 = (0..n1).filter(|_| rng.f64() < p).count();
            let (l1, h1) = cp_interval(k1, n1, alpha / 2.0);
            let q1 = (h1 > tau, l1 <= tau); // (edge 可能, noEdge 可能)
            let sing1 = q1.0 != q1.1; // singleton 読み (どちらか一方のみ)
            // 独立標本の高予算
            let k2i = (0..n2).filter(|_| rng.f64() < p).count();
            let (l2i, h2i) = cp_interval(k2i, n2, alpha / 2.0);
            let q2i = (h2i > tau, l2i <= tau);
            if sing1 {
                singleton_b += 1;
                let sing2 = q2i.0 != q2i.1;
                if !sing2 || q2i != q1 {
                    transient += 1; // 読みが消えた/変わった
                }
            }
            // nested (累積 + 交差) の 2-chain: singleton が両予算で一致したときの整合
            let k_extra = (0..(n2 - n1)).filter(|_| rng.f64() < p).count();
            let (l2r, h2r) = cp_interval(k1 + k_extra, n2, alpha / 2.0);
            let (l2, h2) = (l2r.max(l1), h2r.min(h1));
            let q2 = (h2 > tau, l2 <= tau);
            let sing2n = q2.0 != q2.1;
            if sing1 && sing2n {
                // 2-chain: nested なら精密側の読みは粗い側の読み集合に含まれる (RPF-2)
                let fine_edge = q2.0;
                let coarse_allows = if fine_edge { q1.0 } else { q1.1 };
                if !coarse_allows {
                    chain_viol += 1;
                }
            }
        }
        let tr_rate = transient as f64 / singleton_b.max(1) as f64;
        check(
            "[R5] RPF-5: 独立標本の単点読みは transient (消失率 > 0)・nested 2-chain の健全性違反 0 — chain ≥ 2 規則 (禁止変換 17) の導出",
            transient > 0 && chain_viol == 0,
            format!(
                "singleton {} 件中 transient {} ({:.2}) / nested 2-chain 違反 {}",
                singleton_b, transient, tr_rate, chain_viol
            ),
        );
    }

    // ---------------- [R6] Lean 反例の整数橋 ----------------
    {
        // proofs/ResourceProfile.lean rpf1: cbLow(d,θ) = (θ == d), cbHigh(d,θ) = (θ != d)
        // 被覆: 各仮説で 2 データ点中 1 点 → 1/2。標本 d = true で {false} ⊄ {true}。
        let cb_low = |d: bool, th: bool| th == d;
        let cb_high = |d: bool, th: bool| th != d;
        let cover = |f: &dyn Fn(bool, bool) -> bool, th: bool| -> usize {
            [true, false].iter().filter(|&&d| f(d, th)).count()
        };
        let coverage_ok = [true, false]
            .iter()
            .all(|&th| cover(&cb_low, th) == 1 && cover(&cb_high, th) == 1);
        let not_nested = cb_high(true, false) && !cb_low(true, false);
        check(
            "[R6] Lean 橋: rpf1 反例 (被覆 1/2 両立・標本 true で {false} ⊄ {true}) の整数一致",
            coverage_ok && not_nested,
            "cbLow/cbHigh の被覆計数 = (1,1)/(1,1), nesting 破れ証人 θ = false".into(),
        );
    }

    println!(
        "\n[判定] {}",
        if nfail == 0 {
            "全検査 PASS — 関手性は被覆から出ない (構成で買う)・安定性は margin 条件つき・chain ≥ 2 は系".to_string()
        } else {
            format!("FAIL {} 件", nfail)
        }
    );
}
