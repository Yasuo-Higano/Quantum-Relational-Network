# D2-S — 静的 B3 lane (protocol v32.1, StableEstimate scoped)

**この lane は geometry 能力の blocker を解除しない** (解除候補は
[d2-response.md](d2-response.md) の D2-R のみ)。D2-v1 (「任意の連結グラフ」 —
[../v31.7/d2-v1-superseded.md](../v31.7/d2-v1-superseded.md)) の実行前反例を受けて、
静的 lane の主張は次に縮む:

> B3 核が直接辺と非辺の間に**事前証明された分離マージン**を持つグラフ族に対し、
> 凍結 gap 則が支持 (と下流位相) を再現する。

## 凍結読み出し (D2-v1 と同一 — 規則は変更しない)

- 状態: C = (I + e^{βh})⁻¹, β = 1, h = −(重み行列)。
- 核: B_ij = |C_ij|²。
- 凍結 gap 則 (言語中立):
  1. 非対角値 (i < j) を降順に並べ、窓 = { B ≥ max(B)·10⁻³ } (スケールガード)。
  2. 窓内の隣接対数段差 ln(B_k/B_{k+1}) の最大が ln 3 以上なら、その位置で切断
     (閾値 = 段差両端の幾何平均)。なければ閾値 = ガード (窓全体を辺とする)。
  3. 各順序対 (i,j), i ≠ j を独立に |B_ij| > 閾値 で判定する。

## B3SupportMarginCertificate (回答資格 — 走行前に真値グラフから計算できる)

**Case A** (窓内に非辺が残る場合) — 次の全てが成立すること:
1. 全真辺が窓内: min_{(i,j)∈E} B_ij ≥ guard
2. 分離: min_edge B > max_nonedge B
3. 境界段差 ln(min_edge/max_nonedge) ≥ ln 3
4. 境界段差が窓内段差の**厳密最大** (unique admissible gap)

**Case B** (非辺が全て窓外) — 全真辺が窓内 (strict) かつ窓内 (辺のみ) に
段差 ≥ ln 3 がないこと。

**定理** (`sim/src/bin/v321_d2_erratum.rs` [E3] が n=4..7 全 992 グラフで機械照合・
例外 0): certificate 成立 ⟺ 凍結則の支持が欠0余0。certificate 不成立グラフの正答は
`rejected_no_certificate` — **強制回答は fail と採点する** (非識別セルの正しい棄却は
一級市民 — v31.0 採点原則)。

## 必須負制御 — graph6 `F}oXO`

7 頂点 11 辺 (辺表: (0,1)(0,2)(0,3)(0,4)(1,2)(1,3)(1,4)(2,6)(3,5)(4,5)(4,6))。
凍結参照値 (β = 1 厳密, 8 桁丸め — [unit-d-tolerances.yml](unit-d-tolerances.yml)):

```text
max B               = 0.03296722
min 真辺 B          = 0.01807046
max 非辺 B          = 0.00249751   (対 (2,5), (3,6) — 自己同型で縮退)
次層 (非辺)         = 0.00020236
境界段差            = ln(0.01807046/0.00249751) = 1.9790 ≥ ln 3
窓内最大段差 (尾部) = ln(0.00249751/0.00020236) = 2.5130
```

境界段差が有意 (≥ ln 3) にもかかわらず非辺尾部の段差が窓内最大となり、凍結則は
**13 辺 (余剰 (2,5),(3,6)・欠落 0)** を返す。再現者はこの故障を再現し、certificate
不成立 (条件 4 違反) を確認して `mandatory_negative_control` に報告すること
(schema 必須 — [unit-d-report.schema.json](unit-d-report.schema.json))。

## 参考 (本リポジトリの全数走査 — v321_d2_erratum)

- n = 7 連結同型類 853 (OEIS A001349 照合): 凍結則の故障 22 — **全て余剰のみ・欠落 0**。
- n = 4..6 (6/21/112): 故障 0。
- certificate 成立は 992 中 970 (うち Case B 14) — 成立 ⟺ exact に例外 0。

機構: 熱的相関 C = (I+e^{−A})⁻¹ は A の解析関数なので、連結グラフでは一般に非辺にも
歩道由来の非零相関が乗る。「窓内最大段差」は辺/非辺境界でなく非辺の内部層構造を
拾い得る — B3 の gap 選択が全グラフで直接辺を分離する定理は存在しない。

## 報告

[unit-d-report.schema.json](unit-d-report.schema.json) の `units.D2S` に従う。
`graphs` には certificate 成立グラフ 1 つ以上 (独立に選ぶこと) + 検証結果、
`mandatory_negative_control` には F}oXO の故障再現を記す。
失敗・不一致も同じ形式で提出する。
