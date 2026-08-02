# QRN v34.0 — HOLD-9: 出自 × 文脈整合 × 資源依存局所性の holdout と第三十三期統合

**Version**: v34.0-A (2026-08-02 — 凍結半)
**Sim**: `sim/src/bin/v340a_hold9_freeze.rs` → `results/v340a_hold9_freeze.txt`
(3 検査 PASS + train 20 セル満票)
**位置づけ**: PROMPT/14 の最終版 — 第三十三期の全器械 (v33.1 修復入口 / v33.2
laboratory interface / v33.3 resource profile / v33.4 atlas glue / v33.5 graded
recovery / v33.6 structured backend) を凍結し、**操作の出自 × 文脈整合 ×
資源依存局所性**の識別可能性境界を新鮮 holdout で採点する。

---

## 1. HOLD-9 が採点するもの (HOLD-8 の後継)

従来の 3 計量:

```text
selective risk        = 0   (回答セルの誤り 0 — 強制回答も誤りに数える)
impossibility recall  = 1   (非識別セルの正しい棄却/同値類/構成時拒否)
answerable recall     = 1   (回答可能セルは全て答える)
```

に、第三十三期の 5 計量 (PROMPT/14) を加える:

```text
origin_certificate_coverage          = 1.0  (admit された全操作に出自証明書)
context_witness_coverage             = 1.0  (全 Exact 読みが証人ゲート通過)
raw_operation_promotions             = 0    (証明書なし/流用の昇格 0)
scope_violations                     = 0    (scope 外の裁定強行 0)
transient_factorization_promotions   = 0    (単点読みの局所性昇格 0)
```

## 2. 開封順序 (HOLD-5..8 と同一)

```text
v34.0-A (本コミット) = 生成器・採点器・バー・lib pin の凍結
                      + SECRET コミットメント公表 + train 採点 (可視シード 34001)
  → v34.0-B = SECRET 開示・holdout 初生成・本採点 (調整なし) + 期末完全儀式
sha256(SECRET) = ef6a8cd97b5d7693f4a4ffdb11ccdf42dbf5a971b0c79ba5f0f72f7b15739fcd
```

凍結カーネル (FROZEN-HOLD9 区間 — v340b と逐語一致を [H0] が照合) は第三十三期の
器械を**呼ぶ**駆動部であり、器械本体は lib 6 モジュール。**lib pin**: 凍結時の
sha256-16 (operational_net `7898d244…` / laboratory_interface `a11f188f…` /
resource_profile `f2fe4b96…` / contextual_factorization `e540eea6…` /
graded_recovery `c1ebb02c…` / structured_backend `f57d9bdf…`) を [A0]/[H0] が
照合する — 凍結から開封の間に器械が変わらないことの機械保証 (HOLD-8 の逐語コピー
方式の lib 版)。凍結バー: τ = 1e-3・σ_bar = 0.5・xtalk_bar = 0.1・合成/CAR/J バー
1e-9・orbit 非同値 < 0.9・transient 規則 chain ≥ 2 — 全て第三十三期の各版で凍結済み
の値の再固定。

## 3. セル 5 群 20 セル (隠しパラメータ = 変成 U×置換・対の抽選・γ・コスト・O(6) 回転・汚染振幅・qubit 数)

| 群 | セル | 要求 |
|---|---|---|
| 入力完全性 | IC1 raw + 流用証明書 | **CertificateTargetMismatch 拒否** (raw_operation_promotions = 0) |
| | IC2 net-owned gens | 復元行列 = net primitive **byte 同一** + Exact [2,2,2] |
| | IC3 role-mixed | **RoleMixedRecovery 構成時拒否** |
| | IC4 GKLS → 可換子 lane | **NonDerivation 拒否** (Leibniz 破れ) |
| accessibility | AC1 独立 knobs (変成) | Exact [2,2,2] (出自証明書つき) |
| | AC2 tied (隠し対) | 分解は rank 拒否 + 正直 net は **Abstain** |
| | AC3 同一閉包・別 interface | 両方 Exact・**orbit matching 不在** (非同値) |
| | AC4 合成証明書の流用 | **CertificateTargetMismatch 拒否** |
| context・resource | CR1 budget profile | 資源不足 → [2,2,2] → [2,4] (隠しコスト) |
| | CR2 一意 glue | GluedExact [2,2,2] (frame 回転不変) |
| | CR3 複数 glue | **EquivalenceClassOnly{2}** (tie-break なし) |
| | CR4 cross-talk 跨ぎ | **Straddled 棄却** (強制判定なし) |
| graded | GR1 odd CAR only | **MajoranaFrameOnly** (O(2N) orbit で止まる) |
| | GR2 charge witness | complex modes (**Σâ†â = Q ≤ 1e-9**) |
| | GR3 汚染 witness | **WitnessNotLinearOnFrame 棄却** |
| | GR4 ordinary odd | **構成時拒否** (JW 誤読の遮断) |
| scale・変成 | SC1 dense = Pauli | 隠しセル + 隠し置換で**裁定完全一致** |
| | SC2 大型 structured | Pauli Exact (隠し 32–60 qubit・行列なし)・**dense は ScopeExceeded** |
| | SC3 大域共役の共変 | dims 不変・orbit は W-共役で厳密対応 |
| | SC4 transient 非昇格 | stable = [2,2,2] のみ・単点 [2,4]/[8] は transient (**昇格 0**) |

**train (シード 34001, 可視) は 20 セル満票**: selective risk 0.000・impossibility
recall 1.000 (12/12)・answerable recall 1.000 (8/8)・強制回答 0・出自被覆 31/31・
証人被覆 7/7・昇格/違反 0。設計走行で 12 の独立シードでも満票 (生成器の
シード頑健性 — 隠しパラメータの縮退による器械故障がないことの設計時確認)。

**注**: 外部再現 (Track X, D2-R) の成否は HOLD-9 のセルに**含まれない** — 内部
SECRET から作ったセルは外部独立性ではない (PROMPT/14 非交渉)。

(v34.0-B — SECRET 開示・holdout 本採点・期末完全儀式 — は開封時に追記する)
