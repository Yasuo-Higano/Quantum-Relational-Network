# QRN v27.1 — 核意味論監査: 層分類・多軸台帳・昇格禁止 (第二十八期 開幕)

**Version**: v27.1
**Date**: 2026-07-28
**Sim**: `sim/src/bin/v271_core_audit.rs` → `results/v271_core_audit.txt` (13 検査 PASS)
**新規文書**: [qrn-core-v1-spec.md](qrn-core-v1-spec.md) / [qrn-terminology.md](qrn-terminology.md) / `core.schema.yml`
**位置づけ**: PROMPT/10 §2 の semantic audit。**物理 run なし** — 第二十八期の目的は
新しい物理モジュールを増やすことではなく、理論本体・測定器・橋渡し仮説を
型レベルで分離すること。

---

## 0. 中心判定 — §D.3 の成果物は「core」ではなく adapter + instrument

v27.0 §D.3 が「QRN-Core v1 の公理系」と呼んだもの (格子 = regulator・BOND-A・
中点変調・Belinfante λ = −1/8) を層分類にかけた結果、これらは QRN の存在論的
公理ではなく **Matter-on-Background Adapter の格子契約** (layer: adapter) である
(claims.yml QRN-META-030, spec §1)。QRN の中心命題「量子相関網から時空・重力・
物質・因果が読み出される」の核は:

```yaml
qrn_core:
  kinematics: defined            # 定義のみ — 実装はガウス toy 1 族 (C3)
  dynamics: model_family_only    # 固有原理は未定義
  geometry_bridge: conjectural   # toy 実演のみ — bridge law は未確立
  gravity_bridge: unsupported    # 支持する結果 0 件 (v27.0: graviton pole なし)
  empirical_prediction: none
```

この状態表示は `v271_core_audit` [S0]/[S11] が spec と schema の両方で機械照合する。

## 1. 成果物

1. **[qrn-terminology.md](qrn-terminology.md)** — 7 層 (core/dynamics/adapter/
   instrument/bridge/phenomenology/meta) の定義、型境界 (ノード 3 型・計量 2 型・
   時間 4 型・証拠 4 型)、禁止変換 7 種、旧名の限定 (QrnState →
   GaussianFermionState / QrnModel → GaussianToyModel — コード改名は v27.2)。
2. **[qrn-core-v1-spec.md](qrn-core-v1-spec.md)** — QRN-Core v1 の 4 層分割仕様、
   composite graviton 封印条件の拡張 (CG0–CG8)、Core v1 完了条件。
3. **core.schema.yml** — 全 67 概念の層分類 (機械可読)・6 軸語彙・禁止変換・
   PRED 分類表。**未分類 0 件** (v27.2 への進行ゲート [S1])。
4. **claims.yml の多軸化** — 全 214 主張に 6 軸 (layer / evidence_kind /
   independence / universality / data_relation / physical_scope) を付与。
   書式検査は `v61_ledger` [7]、意味論 (昇格禁止 R1–R7) は `v271_core_audit`。

## 2. 過大解釈の修正 (本監査で最も重要な既存主張の修正)

### 2.1 QRN-C0-001 —「第一法則 = 線形化 Einstein 方程式の同値性」の限定

旧: 「δS=δ⟨K⟩ と線形化 Einstein 方程式の同値性 (Jacobson 1995 / Faulkner et al.
2013) は確立された研究成果である」— **無条件の同値として引用していた**。

新: 前提と結論を区別する:
- **Faulkner et al. 2013** (arXiv:1312.7856): ホログラフィック CFT の全球領域への
  第一法則 ⟺ 双対バルクの線形化 Einstein 方程式 — **半古典ホログラフィック双対・
  RT/Wald 辞書が前提**。
- **Jacobson 1995** (gr-qc/9504004): 局所 Rindler 地平線への δQ=TdS と面積
  エントロピー仮定からの熱力学的導出 — **entanglement equilibrium の論文ではない**。
- **Jacobson 2015** (arXiv:1505.04753): entanglement equilibrium (小測地球・固定
  体積・真空停留性・CFT 的仮定) からの導出。

### 2.2 QRN-GRAV-001 (v0.7) — instrument benchmark への降格

v0.7 が実際に検証したのは **1+1 次元ガウス系におけるエンタングルメント第一法則
(QI の一般恒等式) とモジュラー核の数値再現**であり、Einstein 方程式ではない。
本系 (平坦 1+1D 自由鎖) は FGHMV の前提 (ホログラフィック双対) を満たさない。
A3 の直接証拠から **C1 instrument benchmark / modular-response benchmark** へ降格
(layer: instrument)。uft-v0.7.md §2b と uft-v1.0.md A3 に監査注記 (原文は保存)。

**この修正は研究を弱めない — QRN 固有の bridge law がまだ存在しないことを正確に
露出させる。**

## 3. 多軸台帳の分布 (2026-07-28, 全 214 主張)

| layer | 件数 | | evidence_kind | 件数 |
|---|---|---|---|---|
| phenomenology | 78 | | reproduction | 107 |
| instrument | 62 | | theorem | 38 |
| bridge | 31 | | interpretation | 29 |
| meta | 30 | | calibration | 24 |
| dynamics | 7 | | mechanism_demo | 12 |
| adapter | 6 | | internal_holdout | 4 |
| core | 0 | | natural_observation / external_replication | **0 / 0** |

**layer: core の主張は 0 件** — 「QRN 本体」を支える主張はまだ一つも無い、が
台帳の正直な姿である (定義は書けたが主張はない)。independent_author は C0
(外部の確立結果) の 7 件のみ = **独立外部再現 0** の機械化 [S5]。
scored-hit の γ_UT (PRED-007)・θ13 (PRED-011) は**公知測定値の holdout**
(data_relation: preregistered_holdout) であり「自然の観測量の的中」に数えない —
的中 0 は [S10] が natural×future_observation×scored-hit = 0 として機械保証。

## 4. 昇格禁止規則 (恒久 — スイート常時実行層)

```text
R1/R2  natural_observation・external_replication = 0 件           [S4]
R3     independent_author は C0 のみ                              [S5]
R4     core/dynamics/bridge (C0 以外) は toy/effective_model のみ  [S6]
R5     future_observation を根拠にする主張は存在しない             [S7]
R6     continuum_universal は adapter/instrument (または C0) 専有  [S8]
R7     internal_holdout ⇒ preregistered_holdout                   [S9]
+ 型レベル禁止変換 7 種 (impl From の不在を [S12] が全走査)
+ 文書アンカー (状態表示・監査注記・限定・降格・README 残高) [S11]
```

## 5. 開発記録 (教訓)

- [S12] 初走 FAIL: 禁止 impl パターンの検査が**監査自身のパターン定義文字列を
  自己検出**した — v61_ledger [6] と同型の exempt (監査基盤自身は対象外) で解決。
  「文字列走査の監査は自分自身を除外してから走らせよ」。
- 多軸の自動割り当てでは、期統合テーゼ (META) が本文中の偶発キーワード
  ("Lean"・"python") で same_author_clean_room に誤判定された — META は
  same_implementation に固定。キーワード規則は族単位のレビューが必須。
- v26.6 (QRN-GRAV-039) は「連続極限」の字句で continuum_universal に誤判定 —
  主結果は「bare c₁ は非共変 regulator 汚染に支配」なので scheme_dependent が正。

## 6. スイート (2026-07-28)

増分スイート (`make suite`): **実行 8 本 (監査層 + v271) PASS 64 / 引用 155 本
PASS 1075 (sha256 不変) / 総計 PASS 1139 / FAIL 0** — results/suite_incremental.txt。
開発記録: 素の `cargo fmt` が旧期バイナリ 21 本を再整形しスイートが「ソース変更」
と誤認 → 整形と上書きされた結果を revert し、以後は新規ファイルのみ整形する
(「フォーマッタはリポジトリ全体に無差別にかけない」)。v25.2 凍結対象は無傷。

## 7. 残り (第二十八期)

- v27.2: 型付き core contract (Rust 型状態 + Lean 禁止昇格定理・旧名改名 —
  sim/src 共有部の変更を一括し、全スイート再検証の儀式は v28.0 で実施)。
- v27.3: instruments.yml (adapter/metrology の凍結)。
- v27.4: reproducer/ + replications.yml (外部再現単位)。
- v28.0: Core v1 完了条件の判定と期統合。
