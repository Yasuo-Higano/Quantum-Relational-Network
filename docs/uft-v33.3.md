# QRN v33.3 — Resource-Filtered OperationalNet: poset 上の factorization profile

**Version**: v33.3 (2026-08-02)
**Sim**: `sim/src/bin/v333_resource_profile.rs` → `results/v333_resource_profile.txt`
(7 検査 PASS) / 共有契約 = `sim/src/resource_profile.rs` (新設)
**位置づけ**: PROMPT/14 第三版。可アクセス性は resource budget に依存する — ゆえに
v33.3 の中心対象は「単一因子分解」ではなく **profile**:

> **OperationalFactorizationProfile: ResourceBudget ↦ FactorizationReading**

である。budget は v33.2 の 5 成分半順序 (時間・振幅・帯域・深さ・誤差) のまま扱い、
有限 poset (grid) 上の constructible profile として読む。

---

## 1. budget chain — 読みは budget の関数 [R1]

同一 interface (site 6 本 [深さ 1]・X₁X₂ [深さ 2]・X₂X₃ [深さ 3]) の深さ鎖:

```text
depth 0.5 → NoAccessibleOperations   (資源不足の正直な記録 — Abstain と区別)
depth 1   → Exact [2,2,2]            (局所制御)
depth 2   → Exact [2,4]              (entangler が sites 1,2 を併合)
depth 3   → Exact [8]                (完全制御 = global, 局所性なし)
```

単一の「正解因子分解」は存在しない — 資格 (出自証明書) のない操作は budget を
いくら積んでも現れない (禁止変換 14 は profile 上でも維持)。

## 2. poset は barcode ではない [R2]

比較不能な budget 対 (amp 2, depth 1) と (amp 1, depth 2) が、**同じ dims [2,4] で
別の gauge orbit** を持つ:

- (2,1): X₁X₂ が accessible → orbit {12|3}
- (1,2): X₂X₃ が accessible → orbit {1|23}
- matching 不在 (最良 min-overlap 0.2500)・join (2,2) → [8]

分裂・併合は一次元の出生死滅 (birth–death) で書けない — 最初から「barcode」と
呼ばず、**有限 poset 上の constructible profile** として定義する。zigzag persistence
等への昇格は、写像と安定性定理が成立した後にのみ行う (PROMPT/14)。

## 3. 昇格規則 — transient を局所性に昇格しない (禁止変換 17) [R3]

> 昇格可能な局所性は、単一閾値で一瞬だけ出現した因子分解ではなく、誤差・budget
> perturbation に対して一定領域で同じ gauge orbit を保つ因子分解とする。

凍結規則: **stable ⟺ (読み + gauge orbit) の同値クラスの領域に比較可能な対
(chain ≥ 2) が存在する**。深さ鎖 {1, 1.5, 2, 3} の機械裁定:

| クラス | 領域 | 裁定 |
|---|---|---|
| Exact [2,2,2] (orbit 一致) | {1, 1.5} — chain | **stable** |
| Exact [2,4] | {2} 単点 | transient (昇格しない) |
| Exact [8] | {3} 単点 (grid の頂) | transient — **調査 grid 相対の正直な記録** |

transient_factorization_promotions = 0 (HOLD-9 の採点語彙)。誤差軸の perturbation
は可換子証明書の区間が担う (跨ぎは Straddled — v32.2 以来の契約)。

## 4. command 再パラメータ化不変性とスカラー潰しの禁止 (禁止変換 18) [R4]

- **成分ごとの狭義単調再パラメータ化** φ = (amp → 3·amp, depth → depth²) を
  コストと grid の両方に適用すると profile は**点ごとに不変** (accessibility の
  成分比較は単調写像と可換)。単位の取り方は読みを変えない。
- **恣意的な重み付き和での全順序化は不変でない**: スカラー amp + depth に潰すと
  budget (2,1) [スカラー 3] で X₂X₃ [コスト (1,2), スカラー 3] が accessible に
  なり、読みが [2,4] から [8] へ**反転**する。ResourceBudget → ScalarResourceCost
  の変換を禁止 (**禁止変換 18**) — 恣意的な重みは新しい選択バイアスである。

## 5. 頂は経路を消す — erasure 対照 [R5]

コスト構造だけが異なる 2 つの interface:

- P (site-first): 低予算 → Exact [2,2,2] → 頂 [8]
- Q (entangler-first): 低予算 → Abstain (entangler 2 本は可換 joint) → 頂 [8]

**top budget の読みは合流し ([8])、経路の情報 (どの操作から局所性が組み上がったか)
を消す** — 最終 budget だけを見ると v32.2 の erasure no-go に戻る。profile が
その情報を運ぶ。v32.2「閉包は marking を消す」の resource 版。

## 6. 型契約の登録 [R6]

- `sim/src/resource_profile.rs` — `ResourceFilteredInterface` (資格つき操作 +
  文脈レシピ・budget filter)・`ProfilePoint` 3 種 (NoAccessibleOperations /
  InputRejected / Reading — 資源不足・構成時拒否・裁定を混ぜない)・
  `OperationalFactorizationProfile` (classes / transient_points)・昇格規則
  (chain ≥ 2)。`same_gauge_orbit` を operational_net.rs へ lib 移植 (v32.3 [F3]
  判定器 — v323/v330 の凍結局所コピーは不変)。
- `core.schema.yml` に概念 4 種 + **禁止変換 17/18** を登録。

## 7. 正直な残高

- grid は有限・調査者が宣言する — stable/transient は **grid 相対**の裁定であり、
  grid の拡張だけが頂の transient を解除できる (連続 budget 空間の安定性定理は
  未構成)。
- 誤差軸の perturbation 安定性は可換子証明書の区間 (Straddled) に委ねたまま —
  budget と誤差の合成摂動の走査は HOLD-9 の変成セル。
- profile 間の写像 (interface の粗視化・restriction に対する関手性) と安定性定理は
  未構成 — それが立つまで zigzag/multiparameter persistence の語彙を使わない。
- 文脈は control lane のみ (v33.2 の測定・準備・drift 文脈は保持まで) — overlap
  整合と glue は v33.4。
- bridge law 登録簿は全能力で空・PRED-019 未登録・自然の的中 0・
  `external_replications = 0` のまま — Unit D2-R の公募が引き続き最優先。
