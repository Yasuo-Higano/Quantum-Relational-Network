# QRN v27.3 — 器械台帳の凍結: Adapter v1 / Metrology Suite v1

**Version**: v27.3
**Date**: 2026-07-29
**Sim**: `sim/src/bin/v273_instruments.rs` → `results/v273_instruments.txt` (8 検査 PASS)
**新規台帳**: [instruments.yml](../instruments.yml) (器械 22 + 常設回帰 5)
**位置づけ**: PROMPT/10 §5。uft-v27.0.md §D.3 の器械群を正名で凍結する —
**QRN-Matter-on-Background Adapter v1** (外部計量と物質の結合規約, 8 器械) と
**QRN-Metrology Suite v1** (較正済み測定器, 14 器械)。物理 run なし。

---

## 1. 台帳の構造 — 各器械 16 フィールド

全 22 器械に登録: id / family (adapter|metrology) / concept (core.schema.yml との
層一致を機械検査) / name / input_type / output_type / normalization / regulator /
continuum_contract / calibration_source / negative_controls / known_failure_modes /
allowed_claims (claims.yml の id — 実在検査) / forbidden_interpretations /
certificate_code / **certificate_sha256_16**。

**認証ハッシュの意味**: 器械 (認証バイナリのソース) の無断変更は v273 の FAIL に
なる。器械を改良するときは再認証 (較正の再実行 + ハッシュの意識的更新) を要する —
v27.1 の cargo fmt 事故 (旧期バイナリ 21 本の無断整形) のような事象を、以後は
監査が捕まえる。

## 2. adapter 8 器械 / metrology 14 器械

- **adapter**: BOND-A 結合則・中点変調・Belinfante 改良 (λ = −1/8)・接触完備化
  (独立 2 実装)・counterterm Λ√g (要件 0: 背景停留)・staggered (2-taste)・
  Wilson (1-flavor 対照)・continuum trajectory。
- **metrology**: null 結合 ladder・殻積分 (root-solve)・Matsubara Ward (144 恒等式)・
  導出モデル外挿・spectral measure (不変量 M² 採点)・解析 oracle (三重導出)・
  f-sum rule・BR 射影辞書 (Lean 証明済み)・厳密ブロック分解・dd 経路・区間演算
  (v25.2 凍結対象)・クランプ梯子・universality 監査 (PRED-016 常設化)・fork 監査
  (1/Π の一回限り使用記録)。

各器械の known_failure_modes には開発記録の教訓を凍結した (「殻は root-solve で
切れ」「外挿モデルは観測量ごとに導出」「√2 停留は位相不一致」「不定計量の縮約は
重みまで検査」「和則は密度と測度の組で一つ」等) — 較正の失敗史も台帳の一部である。

## 3. 常設回帰 (§D.3 要件 4 の執行)

suite の再走行 (バイナリは suite 対象のまま) に加えて、v273 が**記録照合**を毎走行
実施する二重化:

| 回帰 | 記録 | 凍結判定文 |
|---|---|---|
| REG-UNIVERSALITY-4RATIO | v268p | PRED-016 scored-hit — 4 比が登録バーで成立 |
| REG-SCALAR-SUMRULE | v268s | PRED-017 scored-hit — 和則が taste 数どおり収束 |
| REG-WARD-64 | v269w | 64 恒等式が格子上で厳密に閉じる |
| REG-WARD-KERNEL | v270a | full 4D kernel の厳密 Ward 充足 |
| REG-FORK-EXTERNAL | v270c | 分岐 (b) external metric 確定 — graviton を作らず |

記録から凍結判定文が消える・FAIL が混入する・認証ハッシュが動く、のいずれでも
suite が赤くなる。

## 4. 検査 (v273_instruments, 8 検査 PASS)

[I0] 解析 (器械 22 + 回帰 5) / [I1] 必須 16 フィールド完備・family 語彙 /
[I2] **concept ↔ core.schema.yml の層一致** (adapter → adapter, metrology →
instrument — 台帳間の分類が単一の層辞書に従う) / [I3] 較正記録・認証ファイルの
実在 + SHA-256 一致 / [I4] allowed_claims の実在 / [I5] 禁止解釈の全器械明示 /
[I6] 常設回帰 5 本の記録照合 / [I7] PROMPT/10 §5 指定 9 器械 + Wilson 対照の被覆。

恒久解釈 (spec §12.8/§13.3 の凍結を台帳の 1 行目に転記): **全器械の成果は
「測定器が正しい」ことの証明であり、QRN・創発重力の証拠ではない。**

## 5. 残り (第二十八期)

- v27.4: reproducer/ + replications.yml (外部再現単位)。
- v28.0: 完全再走の儀式 + Core v1 完了条件の判定 + 期統合。
