# QRN v34.2 — standalone operational core spec (OCS-1.0): 論文だけで閉じる

**Version**: v34.2 (2026-08-03)
**Sim**: `sim/src/bin/v342_spec_manifest.rs` → `results/v342_spec_manifest.txt` (8 検査 PASS)
**Paper**: `paper/operational-core-spec.md` (**OCS-1.0** — sha256 = bb3b4d4bff722225… で凍結) +
`paper/operational-core-closure.yml` (publication closure manifest, 20 エントリ)
**位置づけ**: PROMPT/15 §3。FollowUp の最終判定「supplied paper set does not contain a
complete core specification」への応答 — **多くの paper claim が repository replay で
しか再現できない**という指摘を、独立実装の入力を spec ひとつに閉じることで解消する。

---

## 1. 何を作ったか

**Operational Locality from Certified Quantum Interfaces: Identifiability,
Contextual Descent, and Finite-Data Guarantees — Core Specification (OCS-1.0)**
(作業名「証明付き量子インターフェースからの操作的局所性」)。QRN の統一描像から
切り離した量子情報・量子制御・system identification の独立成果として、第三十三期の
operational core を**ソースにも出力値にもアクセスせず実装できる形**で定義する:

| spec 節 | 内容 (凍結) |
|---|---|
| §0 | 再現性 4 等級 (repository replay / paper-closed / clean-room / organizationally external)・spec 版バンプの規約 (結果に合わせない) |
| §1 | 有限次元系・役割 4 型 (準備/制御/測定/drift — 暗黙変換なし)・grading 2 レーン |
| §2 | 出自 3 証明書 (addressability: rank + σ_min ≥ 0.5 + cross-talk ≤ 0.1 区間 / 合成: bracket 語の機械再実行 ≤ 1e-9 / トモグラフィ) + sha256 結束 (流用拒否) |
| §3 | role-typed 文脈 4 型・joint measurability は joint POVM 証人のみ (可換性から推論しない) |
| §4 | ResourceBudget 5 成分半順序・スカラー化禁止・昇格規則 chain ≥ 2 |
| §5 | marked recovery 凍結手順 (τ = 1e-3 グラフ → 成分 → witness ゲート → 閉包 → 中心 → 裁定)・gauge orbit matching (バー 0.9)・裁定語彙 |
| §6 | chart 局所復元・glue 条件・glue 定理・cocycle 不整合 Abstain |
| §7 | Majorana frame 資格・O(2N) vs U(N)・charge witness → J → モード回復 |
| §8 | **probe 型分離** (下記 §2)・「局所摂動」の語の単独使用禁止 |
| §9 | Pauli GF(2) lane (ω・radical = 中心)・quadratic lane (2^{2m−1})・ScopeExceeded |
| §10 | 規範インスタンス N1–N5 (構成のみ・出力値なし) |
| §11 | no-go / 正定理 T1–T10 (各 falsifier・禁止解釈つき) |
| §12 | 裁定の全順序 (fail-closed 7 段: OutOfDomain > 構成時拒否 > Insufficient > Straddled > EquivClass > Sectors > Exact) |
| §13 | 凍結バー全表 | 
| §14 | 閉じない領域の正直な宣言 (holdout harness・campaign・有限データ意味論) |

## 2. FollowUp の probe 型分離の受理 — 主リポジトリ側の機械実証 [M4]

FollowUp が clean-room で発見・実装した 2 つの分離を、規範インスタンスとして spec に
凍結し、**本リポジトリでも初めて機械実証**した:

```text
[M4a] N1 quench null (T9): h = σ_x・Γ₀ = I/2 で
      SignedInitialCovarianceProbe: (n̈⁺−n̈⁻)/(4ε) = ‖P₂hP₁‖² = 1 (厳密)
      HamiltonianQuench (h ± εP₁):  応答 = 0 (厳密 — 状態が全生成子と可換)
      → 同一の h で probe 型が答えを変える: 型は観測契約の必須フィールド
        (SignedInitialCovarianceProbe ↛ HamiltonianQuench)
[M4b] N2 pairing (T10): 2 モード JW・積状態 p₁ = 1/2 ± ε:
      H = c₁†c₂ + h.c. + Δ(c₁†c₂† + h.c.) → 応答 = 1 − Δ² 厳密
      (Δ = 0/0.3/0.7 → 1/0.91/0.51・対角 V n₁n₂ は不変)
      → 数保存応答則の証明書は charge witness を運ぶ
        (NumberConservingResponse ↛ BCS/PairingResponse)
[M4c] N3 Busch 対: R_ab = (I + aησ_x + bησ_z)/4 の joint POVM 資格 ⟺ 2η² ≤ 1
      (η = 0.6: min eig +0.0379 資格 / η = 0.8: −0.0328 拒否・‖[E^x,E^z]‖ = 0.255)
      → joint measurability は可換性より広く、証人は joint POVM のみ
```

導出注記 (N2): 応答の 3 部分解 — hopping [T,[T,n₂]] = −2(n₁−n₂) が +1・
pairing [P,[P,n₂]] = 2Δ²(n₁+n₂−1) が −Δ²・交差項は厳密 0。ε に線形なので
任意の ε で厳密 (数値も 1e-12 以下で一致)。

## 3. publication closure manifest — 「どの等級で閉じているか」の機械化

`paper/operational-core-closure.yml`: 20 エントリ (OCS-K1/A1–A3/X1/R1/F1/T1/T3/T4/
G1/M1/P1–P3/B1/B2/D1 = paper_closed 18・OCS-H1 [HOLD harness]/OCS-Z1 [campaign] =
repository_replay_only 2)。各エントリは claim_ids (claims.yml へ)・assumption_ids・
observation_contract・equivalence_relation・allowed_conclusion・falsifier・
`expected_outputs_in_spec: false` (全行) を持ち、v342 [M2] が参照実在と
「正直な境界 ≥ 1 件」(閉包の過大主張の禁止) を常設検査する。

**再現性の 4 等級を区別する** (spec §0): 現状の core claim 群は等級 1 (repository
replay) + 本版で等級 2 (paper-closed) の主張が立った。等級 3 (clean-room) は
FollowUp 型の実走行のみが与え、等級 4 (organizationally external) は 0 のまま。

## 4. 凍結の規律

- spec の sha256 (`bb3b4d4bff722225…`) を v342 [M0] に pin — 変更は「意図的な二重
  更新 + OCS 版バンプ」のみ (instruments.yml と同じ規律)。
- **実行順序** (spec §0): 実装は互いの結果を見ずに spec から書く。曖昧性が見つかれば
  **結果に合わせず spec の版を上げる**。「論文更新で FollowUp が走る」経路の入力が
  これで一意になった。
- 有限データ (shot noise) 意味論は **OCS-1.0 のスコープ外と明記** (§14) — v34.3 の
  第四 no-go と Robust Promotion Theorem が証明されてから spec 化する (定理より先に
  語彙を凍結しない)。

## 5. 残高

bridge law 登録簿は全能力で空・PRED-019 未登録・自然の的中 0・external_replications
= 0 — 全て不変。paper_closed は「独立実装可能」の主張であって「独立実装が完了した」
事実ではない (QRN-META-044 limitations)。

## 6. 次 (v34.3)

有限データ昇格不能定理 (第四の no-go — Le Cam 二点下限の有限版 Lean 証明) と
Robust Promotion Theorem (同時信頼集合 C_α(D) 上の裁定 RobustExact /
EquivalenceClassOnly / Straddled / InsufficientObservation / OutOfDomain・
P(wrong promotion) ≤ α と selective risk ≤ α/P(answer) の区別)・禁止変換 22–29。
