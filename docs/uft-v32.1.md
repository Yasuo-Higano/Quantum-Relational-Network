# QRN v32.1 — Unit D2 プロトコルの型修復と反例登録 (第三十二期 開始)

**Version**: v32.1 (2026-08-01)
**Sim**: `sim/src/bin/v321_d2_erratum.rs` → `results/v321_d2_erratum.txt` (8 検査 PASS)
**位置づけ**: PROMPT/13 第一版 — 第三十二期の中心テーゼ
**「局所性は状態に宿るのではなく、操作的文脈の可換分解と、その上の Liouvillian
応答の同値類としてのみ識別される」**の開幕として、全主張の独立性 blocker である
外部再現 Unit D2 のプロトコル欠陥を修復・凍結し、公募を開始可能にする
(Track X — 優先順位 0。主研究 OperationalNet [v32.2–] とは並列)。

---

## 1. 何が壊れていたか — 2 つの型不整合

### 1.1 提出契約の型不整合

`SPEC.md` は全報告が `EXPECTED_SCHEMA.json` に従うと規定するが、同ファイルは
A–C 用の「フィールド説明を値に書いた例示 JSON」であり、Unit D (D1–D4) の報告を
機械検証できる契約が存在しなかった — **第三者が「機械的に妥当な D2 報告」を作れない**。

### 1.2 D2-v1 の「任意の連結グラフ」は強すぎる (実行前反例)

熱的相関 C = (I+e^{−A})⁻¹ は A の解析関数であり、連結グラフでは一般に非辺にも
歩道由来の非零相関が乗る。B3 の凍結 gap 則 (スケールガード窓 max·10⁻³・窓内最大
対数段差 ≥ ln 3 で切断 — v31.6 最終規則) が全グラフで直接辺を分離する定理はない。

**反例 graph6 `F}oXO`** (7 頂点 11 辺 — PROMPT/13 の NetworkX+numpy 独立走査を、
本器械が独立デコード + 逐語凍結則で機械確定 [E0][E1]):

```text
max B               = 0.03296722
min 真辺 B          = 0.01807046
max 非辺 B          = 0.00249751   (対 (2,5), (3,6) — 自己同型で縮退)
次層 (非辺)         = 0.00020236
境界段差 ln(0.01807046/0.00249751) = 1.9790 ≥ ln 3   ← 正しい切断はここ
尾部段差 ln(0.00249751/0.00020236) = 2.5130          ← 「最大段差」則はこちらを選ぶ
→ 凍結則の報告 = 13 辺 (余剰 (2,5),(3,6)・欠落 0)
```

**全数スキャン [E2]**: n = 4..7 連結同型類 (6/21/112/853 = OEIS A001349 照合) で
故障は **n=7 の 22/853 のみ・全て余剰のみ (欠落 0)**・n ≤ 6 は 0。
これは**外部再現の失敗ではなく、外部走行 0 件の時点で発見されたプロトコル反例**である
— D2-v1 を黙って修正せず `superseded_before_external_run` として登録した
(replications.yml `ERR-D2-V1`・旧文面の逐語保存 = `reproducer/protocols/v31.7/d2-v1-superseded.md`)。

## 2. 修復 — 版分離と 2 lane 化

```text
reproducer/protocols/
  v27.4/  abc-report.schema.json / abc-tolerances.yml   (凍結原本の版付き複製 — byte 一致を常設監査)
  v31.7/  d2-v1-superseded.md                           (旧 D2 の逐語保存 + 反例の要旨)
  v32.1/  unit-d-report.schema.json / unit-d-tolerances.yml
          d2-static.md / d2-response.md / protocol-index.yml   (凍結 — sha256 認証 v310 [R9])
```

### 2.1 D2-S — 静的 B3 lane (StableEstimate scoped; blocker は解除しない)

主張を「**事前証明された分離マージンを持つグラフ族**に対する支持と下流位相の再現」
に縮めた。回答資格 **B3SupportMarginCertificate** (真値グラフから走行前に計算可能):

- **Case A** (窓内に非辺が残る): 全真辺が窓内 ∧ 分離 ∧ 境界段差 ≥ ln 3 ∧
  境界段差が窓内段差の**厳密最大** (unique admissible gap)。
- **Case B** (非辺が全て窓外): 全真辺が窓内 (strict) ∧ 辺内部に段差 ≥ ln 3 なし。

**[E3] certificate ⟺ 凍結則 exact を n=4..7 全 992 グラフで機械照合 — cert 970
(うち Case B 14) = exact 970・例外 0**。F}oXO は certificate 不成立 (条件「厳密最大」
違反) の**必須負制御** — certificate 不成立グラフの正答は `rejected_no_certificate`
であり強制回答は fail (v31.0 の採点原則「棄却は一級市民」の外部単位への適用)。

### 2.2 D2-R — 応答 lane (end_to_end; **geometry blocker 解除の本命**)

曲率則 (n̈_j⁺−n̈_j⁻)(0)/(4ε) = ‖P_j h P_i‖_F² は任意の一体生成子に対する**有限行列
恒等式** (v31.2) であり、静的 B3 の長距離相関混入を受けない。D4 (恒等式確認) を
トポロジーまで接続した end-to-end lane に昇格:

```text
独立 h 生成 → 独立 probe 準備 (C± = I/2 ± εP_i) → 密度時系列のみ
→ 曲率 → ŵ = |h_ij|² → gap 支持 → clique complex → Z2 homology (β₃ は ∂₄ 必須)
→ vertex link → 必須 6 セル (T³ L≥4 / 16-cell / 3-ball / T³ L=3 負制御 /
   独立 sparse weighted / 高ノイズ Abstain)
```

**[E4] 応答 lane は全 22 反例 + F}oXO を修復**: 支持 欠0余0・辺/非辺比 ≥ 8.6e8・
|ŵ − 1| ≤ 9.5e-10。成功で解除されるのは **spatial_topology_given_factorization の
独立性 blocker のみ** — E3 因子分解・PRED-019・bridge 登録は別問題のまま。

### 2.3 報告契約の型修復 [E6]

`unit-d-report.schema.json` は**実 JSON Schema (draft 2020-12)**。自作最小 validator
(JSON パーサ + schema 部分語彙 + 最小 regex — 外部クレートなし) の負制御:

| fixture | 期待 | 機構 |
|---|---|---|
| pass 報告 | 適合 | — |
| **failed 報告** | **適合** | 正直な失敗は一級市民 (失敗の隠蔽が反証可能性を壊す) |
| 必須欠落 (protocol_frozen_commit) | 不適合 | required |
| D1 が capability 主張 | 不適合 | additionalProperties: false |
| D2-R が語彙外能力 (clock_calibration) | 不適合 | items enum (許可は spatial_topology_given_factorization のみ) |

## 3. 発見 — SupportNoiseCertificate (重みバーは支持を守らない) [E5]

凍結重みバー (誤差見積り ≤ 0.1 で回答) を**通る** σ = 1e-9 でも、凍結則は順序対ごと
に閾値判定するため、ノイズ最大値が窓ガード (max·10⁻³) を跨いで余剰辺を作る:
F}oXO は見積り 4.6e-4 がガード 1e-3 の ~1.6σ しかなく、**実測で余剰 1 (欠落 0)**。
支持段の回答資格を独立に凍結した:

```text
noise_error_bound(σ, ‖h‖₁) · √(2·ln(m/10⁻⁶)) ≤ 10⁻³ · max ŵ    (m = n(n−1))
```

F}oXO σ=1e-9 は不成立 (2.7e-3 > 1.0e-3 — 正しく棄却)・ring12 σ=1e-9 と F}oXO
σ=1e-12 は成立 + 支持一致。次数 → dt = 0.02/‖h‖₁ → 1/dt² 増幅のグラフ依存が裁定に
正しく現れる (HOLD-7 K8-lownoise [n=12 重み環] の合格とも整合)。教訓は v32.0-B
K3-holes と同じ —「**定量バーは観測量の応答則から導出する**」。非辺の裁定は点推定の
閾値比較ではなく、ノイズ最大値の信頼上界がガードを跨がないことで行う。

## 4. 恒久化

- `v310_readout_semantics` に **[R9]** を追加 (ALWAYS_RUN): 凍結一式 6 ファイルの
  sha256-16 認証・v27.4 版付き複製 = 原本の byte 一致・supersession 台帳の実在。
  凍結ファイルの変更は認証値の意識的更新 (= 版分離) を要する。
- `v321_d2_erratum` [E7]: tolerances の凍結参照値 = 本器械の計算値 (転記ミスの機械
  検出 — v274 [P2]/[P3] と同型)・replications.yml の 6 条件と external_replications
  = 0 の不変。

## 5. 正直な残高

- **external_replications = 0 のまま** — 本版は受け皿の修復であって独立外部再現では
  ない。公募 (Track X) は本版の凍結直後から主研究と並列で開始する。
- D2-S は blocker を解除しない (StableEstimate scoped)。解除候補は D2-R のみで、
  その対象能力も **GivenNodeFactorization を入力に取る** E2 水準 — E3 (因子分解の
  操作的選定) は第三十二期の主研究 (v32.2 OperationalNet) の課題。
- certificate の全数照合は β=1・0/1 隣接・n ≤ 7 — 重みつき・別 β の certificate
  成立分布は未走査 (定義自体は任意重み・任意 β で走行前計算可能)。

## 6. 第三十二期の残り (PROMPT/13 の版計画)

| 版 | 主題 |
|---|---|
| v32.2 | OperationalNet semantics — 型分離 + global-algebra erasure no-go |
| v32.3 | Factorization recovery — commutant/center 証明書・三値裁定・graded scope |
| v32.4 | Liouvillian response hierarchy — 一階/二階恒等式・H ↔ −H no-go |
| v32.5 | Interaction hypergraph — H_S 直交分解・相関 hopping・conditional probes |
| v32.6 | VR exactness — 離散円環 bar 定理・端点規約・H2 persistence |
| v33.0 | HOLD-8 — factorization × interaction order × observation contract |
