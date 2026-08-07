//! v35.1 RECORD v2 — 記録の相関粒度型 (PROMPT/16 §3.2)
//!
//! v34.6 の recorded lane (RECORD.schema.json v1) は各チャネルに shot 配列を必須と
//! したが、lane 文書は「shot 配列または counts + 取得時刻」を許しており、順序を
//! 失った aggregate counts では HOLD-10 が必須化した遷移数ゲート (serial
//! correlation 検査 — 同一周辺分布の持続的 Markov 鎖は split-half を通過する,
//! v343 [F6-29]) が実行できない。データ取得後に schema を直すことは「結果を見た
//! 後の変更」に見えるため、**取得前の今**、相関検査の可能な粒度を型で分離する:
//!
//!   OrderedShots        — split-half + 遷移数ゲートまで実行可能 (iid 資格に到達可能)
//!   TimestampedBatches  — batch 間 drift・overdispersion のみ検査可能
//!                         (batch 内の serial correlation は観測不能 — iid 資格に到達不能)
//!   AggregateCounts     — 相関検査不能 — 常に CorrelationUnassessed
//!
//! 禁止変換 30: AggregateCounts ↛ IidCertificate — 順序情報を失った counts に iid
//!   証明書を発行しない (禁止変換 29 [モデル条件付き証明書] の粒度版)。
//! 禁止変換 31: TimestampedBatches ↛ IidCertificate — batch 間検査を全て通過しても
//!   到達点は CorrelationUnresolved (「検出できる範囲に破れがない」であって
//!   「iid と整合」ではない)。
//!
//! `IidCertificate` は private フィールドを持ち、本モジュールの `assess()` が
//! OrderedShots で全ゲートを通過したときのみ構成する (門は較正 — v33.2 の規律)。
//! ゲートの登録意味論は HOLD-10 (v350a) の凍結形をそのまま採用する:
//!   split-half     — 前半/後半の Clopper–Pearson 区間 (各 α/2) の disjoint 判定
//!   遷移数         — t = #(s_i ≠ s_{i+1}) の CP 区間 vs q = 2p(1−p) の CP 伝播区間
//!   batch drift    — batch 対ごとの CP 区間 (Bonferroni α/(2·対数)) の disjoint 判定
//!   overdispersion — Pearson T = Σ (k_b − n_b p̂)²/(n_b p̂(1−p̂)) vs χ²_{m−1} 分位
//!                    (χ² は近似であることを登録 — 境界較正は v351 [P5] が固定シードで機械検査)

use crate::finite_data::cp_interval;

// ---------------------------------------------------------------- 記録データの相関粒度 3 型

/// TimestampedBatches の 1 バッチ (統計に使うのは計数のみ — 時刻は schema 側のメタデータ)
#[derive(Clone, Debug)]
pub struct BatchRecord {
    pub n_shots: usize,
    pub n_ones: usize,
}

/// 記録データの粒度 — 相関検査の可能性は粒度の属性であり、後から昇格できない
#[derive(Clone, Debug)]
pub enum RecordData {
    /// 取得順の 0/1 列 (遷移数ゲートまで実行可能)
    OrderedShots(Vec<u8>),
    /// タイムスタンプつきバッチ計数 (batch 間検査のみ可能)
    TimestampedBatches(Vec<BatchRecord>),
    /// 集計計数のみ (相関検査不能)
    AggregateCounts { n_shots: usize, n_ones: usize },
}

/// iid 契約検査の通過証 — **OrderedShots からのみ構成可能** (private フィールドが門)。
/// 「登録ゲート (split-half + 遷移数, 各 α/2) で破れが検出されなかった」の記録であり、
/// iid の証明ではない (検出力は n の関数 — v34.5 [R5] の分解能命題と同じ)。
#[derive(Clone, Debug)]
pub struct IidCertificate {
    _gate: (), // private — 本モジュールの assess() だけが構成できる
    pub alpha: f64,
    pub n_shots: usize,
    pub transitions: usize,
}

/// 相関裁定 — 粒度ごとの到達可能集合が異なる:
///   OrderedShots       → {IidConsistent, DriftDetected, SerialCorrelationDetected, CorrelationUnresolved}
///   TimestampedBatches → {BatchDriftDetected, OverdispersionDetected, CorrelationUnresolved}
///   AggregateCounts    → {CorrelationUnassessed}
#[derive(Clone, Debug)]
pub enum CorrelationVerdict {
    /// 登録ゲート全通過 (OrderedShots のみ) — iid 証明書を運ぶ
    IidConsistent(IidCertificate),
    /// split-half CP 区間が disjoint (前半/後半の率が両立しない)
    DriftDetected { k_first: usize, k_second: usize },
    /// 遷移数が iid 期待 q = 2p(1−p) の伝播区間と disjoint (Markov 相関の検出)
    SerialCorrelationDetected { transitions: usize, n: usize },
    /// batch 対の CP 区間が disjoint
    BatchDriftDetected { pair: (usize, usize) },
    /// Pearson 過分散統計量が χ²_{m−1} 分位を超過
    OverdispersionDetected { t_obs: f64, t_crit: f64 },
    /// 検出なし・ただし iid 資格には到達しない (TimestampedBatches の最大到達点 /
    /// OrderedShots の長さ不足)
    CorrelationUnresolved,
    /// 検査不能 (AggregateCounts は常にこれ — 禁止変換 30)
    CorrelationUnassessed,
}

impl CorrelationVerdict {
    pub fn as_str(&self) -> &'static str {
        match self {
            CorrelationVerdict::IidConsistent(_) => "iid_consistent",
            CorrelationVerdict::DriftDetected { .. } => "drift_detected",
            CorrelationVerdict::SerialCorrelationDetected { .. } => "serial_correlation_detected",
            CorrelationVerdict::BatchDriftDetected { .. } => "batch_drift_detected",
            CorrelationVerdict::OverdispersionDetected { .. } => "overdispersion_detected",
            CorrelationVerdict::CorrelationUnresolved => "correlation_unresolved",
            CorrelationVerdict::CorrelationUnassessed => "correlation_unassessed",
        }
    }
}

/// iid 証明書の取り出し — IidConsistent 以外からは何も出ない (禁止変換 30/31 の器械形)
pub fn iid_certificate(v: &CorrelationVerdict) -> Option<&IidCertificate> {
    match v {
        CorrelationVerdict::IidConsistent(c) => Some(c),
        _ => None,
    }
}

// ---------------------------------------------------------------- 登録ゲート定数

/// OrderedShots の登録最小長 (split-half の各半分に 20 shot — これ未満は
/// CorrelationUnresolved: 「検査できないのに資格を出さない」)
pub const MIN_ORDERED_SHOTS: usize = 40;

// ---------------------------------------------------------------- 裁定 (登録意味論)

/// 相関裁定の唯一の入口。α は粒度内の検査全体の同時予算
/// (Bonferroni: OrderedShots = split-half α/2 + 遷移数 α/2 /
///  TimestampedBatches = batch drift α/2 + overdispersion α/2)。
pub fn assess(data: &RecordData, alpha: f64) -> CorrelationVerdict {
    match data {
        RecordData::OrderedShots(shots) => assess_ordered(shots, alpha),
        RecordData::TimestampedBatches(batches) => assess_batches(batches, alpha),
        // 禁止変換 30: 順序を失った counts は検査不能 — 無条件で Unassessed
        RecordData::AggregateCounts { .. } => CorrelationVerdict::CorrelationUnassessed,
    }
}

fn assess_ordered(shots: &[u8], alpha: f64) -> CorrelationVerdict {
    let n = shots.len();
    if n < MIN_ORDERED_SHOTS {
        return CorrelationVerdict::CorrelationUnresolved;
    }
    let a = alpha / 2.0; // 2 ゲートへの Bonferroni 割当
    // (1) split-half (v346 drift gate / v350a h10_gate_splithalf と同形)
    let half = n / 2;
    let k1 = shots[..half].iter().filter(|&&s| s == 1).count();
    let k2 = shots[half..].iter().filter(|&&s| s == 1).count();
    let (lo1, hi1) = cp_interval(k1, half, a / 2.0);
    let (lo2, hi2) = cp_interval(k2, n - half, a / 2.0);
    if hi1 < lo2 || hi2 < lo1 {
        return CorrelationVerdict::DriftDetected {
            k_first: k1,
            k_second: k2,
        };
    }
    // (2) 遷移数ゲート (v350a h10_gate_correlation と同形)
    let k = shots.iter().filter(|&&s| s == 1).count();
    let t = shots.windows(2).filter(|w| w[0] != w[1]).count();
    let (pl, pu) = cp_interval(k, n, a / 2.0);
    let q_of = |p: f64| 2.0 * p * (1.0 - p);
    // q は p ∈ [pl, pu] 上で単調でない (頂 0.5) — 区間像を正確に取る
    let (mut ql, mut qu) = (q_of(pl).min(q_of(pu)), q_of(pl).max(q_of(pu)));
    if pl <= 0.5 && 0.5 <= pu {
        qu = 0.5;
    }
    ql = ql.min(qu);
    let (tl, tu) = cp_interval(t, n - 1, a / 2.0);
    if tu < ql || qu < tl {
        return CorrelationVerdict::SerialCorrelationDetected { transitions: t, n };
    }
    CorrelationVerdict::IidConsistent(IidCertificate {
        _gate: (),
        alpha,
        n_shots: n,
        transitions: t,
    })
}

fn assess_batches(batches: &[BatchRecord], alpha: f64) -> CorrelationVerdict {
    let m = batches.len();
    if m < 2 || batches.iter().any(|b| b.n_shots == 0 || b.n_ones > b.n_shots) {
        return CorrelationVerdict::CorrelationUnresolved;
    }
    let a = alpha / 2.0; // 2 ゲートへの Bonferroni 割当
    // (1) batch 間 drift: 対ごとの CP 区間 disjoint (α/2 を対数で Bonferroni)
    let npairs = m * (m - 1) / 2;
    let a_pair = a / npairs as f64;
    let ivs: Vec<(f64, f64)> = batches
        .iter()
        .map(|b| cp_interval(b.n_ones, b.n_shots, a_pair))
        .collect();
    for i in 0..m {
        for j in (i + 1)..m {
            let (l1, h1) = ivs[i];
            let (l2, h2) = ivs[j];
            if h1 < l2 || h2 < l1 {
                return CorrelationVerdict::BatchDriftDetected { pair: (i, j) };
            }
        }
    }
    // (2) overdispersion: Pearson T vs χ²_{m−1} の (1 − α/2) 分位
    let ntot: usize = batches.iter().map(|b| b.n_shots).sum();
    let ktot: usize = batches.iter().map(|b| b.n_ones).sum();
    let phat = ktot as f64 / ntot as f64;
    if phat > 0.0 && phat < 1.0 {
        let t_obs: f64 = batches
            .iter()
            .map(|b| {
                let e = b.n_shots as f64 * phat;
                let v = b.n_shots as f64 * phat * (1.0 - phat);
                (b.n_ones as f64 - e) * (b.n_ones as f64 - e) / v
            })
            .sum();
        let t_crit = chi2_quantile((m - 1) as f64, 1.0 - a);
        if t_obs > t_crit {
            return CorrelationVerdict::OverdispersionDetected { t_obs, t_crit };
        }
    }
    // 禁止変換 31: 全通過でも到達点は Unresolved (batch 内相関は観測不能)
    CorrelationVerdict::CorrelationUnresolved
}

// ---------------------------------------------------------------- χ² 分位 (正則化不完全ガンマ)

/// 正則化下側不完全ガンマ P(a, x) — 級数 (x < a+1) / 連分数 (x ≥ a+1)。
/// 決定的・std のみ (Numerical Recipes の標準法)。
pub fn gamma_p(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 0.0;
    }
    let ln_ga = crate::ln_gamma(a);
    if x < a + 1.0 {
        // 級数展開
        let mut ap = a;
        let mut sum = 1.0 / a;
        let mut del = sum;
        for _ in 0..500 {
            ap += 1.0;
            del *= x / ap;
            sum += del;
            if del.abs() < sum.abs() * 1e-15 {
                break;
            }
        }
        sum * (-x + a * x.ln() - ln_ga).exp()
    } else {
        // 連分数 (Lentz 法)
        let tiny = 1e-300;
        let mut b = x + 1.0 - a;
        let mut c = 1.0 / tiny;
        let mut d = 1.0 / b;
        let mut h = d;
        for i in 1..500 {
            let an = -(i as f64) * (i as f64 - a);
            b += 2.0;
            d = an * d + b;
            if d.abs() < tiny {
                d = tiny;
            }
            c = b + an / c;
            if c.abs() < tiny {
                c = tiny;
            }
            d = 1.0 / d;
            let del = d * c;
            h *= del;
            if (del - 1.0).abs() < 1e-15 {
                break;
            }
        }
        1.0 - h * (-x + a * x.ln() - ln_ga).exp()
    }
}

/// χ²_ν の p 分位 — P(ν/2, x/2) = p を二分法で解く (決定的)
pub fn chi2_quantile(nu: f64, p: f64) -> f64 {
    let f = |x: f64| gamma_p(nu / 2.0, x / 2.0);
    let (mut lo, mut hi) = (0.0f64, 1.0f64);
    while f(hi) < p {
        hi *= 2.0;
        if hi > 1e8 {
            break;
        }
    }
    for _ in 0..200 {
        let mid = 0.5 * (lo + hi);
        if f(mid) < p {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

// ---------------------------------------------------------------- 自己検証

/// モジュール自己検証 — χ² 分位の教科書値・粒度ごとの到達可能集合
pub fn record_v2_self_test() -> Result<(), String> {
    // χ² 分位の教科書値 (Abramowitz–Stegun): χ²₁(0.95) = 3.8415, χ²₄(0.99) = 13.2767
    let q1 = chi2_quantile(1.0, 0.95);
    if (q1 - 3.841459).abs() > 1e-4 {
        return Err(format!("chi2_quantile(1, .95) = {} ≠ 3.8415", q1));
    }
    let q4 = chi2_quantile(4.0, 0.99);
    if (q4 - 13.2767).abs() > 1e-3 {
        return Err(format!("chi2_quantile(4, .99) = {} ≠ 13.2767", q4));
    }
    // AggregateCounts は常に Unassessed (禁止変換 30)
    let agg = RecordData::AggregateCounts {
        n_shots: 100000,
        n_ones: 50000,
    };
    match assess(&agg, 0.01) {
        CorrelationVerdict::CorrelationUnassessed => {}
        v => return Err(format!("aggregate が {} を返した", v.as_str())),
    }
    // 短すぎる OrderedShots は Unresolved (資格を出さない)
    match assess(&RecordData::OrderedShots(vec![0, 1, 0, 1]), 0.01) {
        CorrelationVerdict::CorrelationUnresolved => {}
        v => return Err(format!("短列が {} を返した", v.as_str())),
    }
    Ok(())
}
