# QRN v34.4 — sector-aware complete finite factorization enumerator

**Version**: v34.4 (2026-08-03)
**Sim**: `sim/src/bin/v344_factorization_enumerator.rs` →
`results/v344_factorization_enumerator.txt` (10 検査 PASS)・共有契約
`sim/src/factorization_enumerator.rs` 新設
**位置づけ**: PROMPT/15 §5。FollowUp FAC-001 の判定「代数不変量の診断は再現・
**complete tensor-factor candidate search は未提供**・非自明 center の tensor 昇格
禁止も確認」への応答。

---

## 1. なにを作ったか — 証明書つき候補列挙器

与えられた *-部分代数の族 (marked 成分の閉包など) 𝒜, 𝒜′, Z(𝒜) から:

```text
1. joint 中心 → 最小中心射影の列挙 (v32.3 の Lagrange kernel を継承)
2. 各 sector α の Wedderburn 証明書:
     dim A_α = n_α² ・ dim A′_α = m_α² ・ n_α·m_α = d_α
     積 span(A_α · A′_α) = d_α²  (A ∨ A′ = B(H_α))
     二重可換子 A″_α = A_α
   — 同型 A_α ≅ M_{n_α} ⊗ I_{m_α} の存在は証明書 + 標準構造定理 (C0)
3. multiplicity 空間 m_α と simple factor n_α の分離
4. 族の可換 simple 部分集合 S + simple 補因子 ((∨S)′ が simple のとき) から
   candidate を列挙 — 全成分の積 span = d² の証明書つき
5. 局所 unitary × 成分置換の witness (traceless 部分空間の overlap matching —
   OCS-1.0 §F3 と同一意味論・バー 0.9)
6. 非同値な候補は集合のまま返す (tie-break 禁止)
```

出力 6 型 (凍結): `UniqueFactorization` / `FactorizationCandidateSet` /
`SectorwiseFactorization` / `IncompletePrimitiveSet` /
`NontrivialCenterObstruction` / `ScopeExceeded`。

## 2. 検証セル

| セル | 構成 | 裁定 |
|---|---|---|
| W1 | M₂⊗M₃ (d=6, qubit×qutrit) | Unique [2,3] |
| W2 | 3 qubit site 族 (d=8) | Unique [2,2,2] |
| W3 | multiplicity {a⊕a} ⊂ M₄ | 証明書 n=2, m=2・Unique [2,2] |
| W4 | M₂⊕M₃ (d=5, 中心非自明) | Sectorwise [(2,1),(3,1)]・**tensor 要求 → NontrivialCenterObstruction** |
| W4b | {a⊕a⊕b} (d=7) | Sectorwise [(3,1),(2,2)] — sector 内 multiplicity |
| W5 | site vs CNOT 共役の 2 bipartition (d=4) | **FactorizationCandidateSet{2}** (両方 [2,2]・overlap 1/3) |
| W6 | number op のみ (abelian 閉包) | **IncompletePrimitiveSet** (rank-1 sector の皮を被せない) |
| W7 | d = 128 > バー 64 | ScopeExceeded (試行しない — 正答) |
| W8 | SWAP 置換 / site×bell | witness: 同一 orbit 1.000000000000 / 非同値 0.3333 |

W6 の裁定は設計判断: abelian 族は形式的には「rank-1 sector ×d の直和」だが、
それを SectorwiseFactorization と呼ぶのは因子候補として空虚 — v32.3 の
number-op-only = Insufficient の enumerator 版として IncompletePrimitiveSet を
返す (裁定の皮を被せない)。

## 3. 位置づけと接続

- **E3 no-go との整合**: 候補列挙は「宣言された族の可換部分集合が生成する分解」
  に限る — 無制約の全 bipartition 走査はしない (因子分解は marking が運ぶ。
  無制約列挙は d=4 でも連続無限個の [2,2] を持ち、選択バイアスの再導入になる)。
- **v32.3/v33.x recovery との分担**: recovery は marked net (可換子グラフ +
  文脈証人) から読む。enumerator は代数入力に対する構造解析 + 候補列挙 —
  recovery の出力 (成分閉包) を入力にでき、超選択 (中心非自明) の構造を
  multiplicity 込みで証明書化する。
- **v34.5 への入力**: 本 enumerator の裁定 (Unique/CandidateSet/…) を有限データの
  同時信頼集合上に持ち上げる (集合の全要素で候補集合が一致するときのみ回答)。

## 4. 残高

bridge law 空・PRED-019 未登録・自然の的中 0・external 0 — 不変。

## 5. 次 (v34.5)

robust atlas — addressability σ_min の同時下界・cross-talk の worst-case 上界・
glue overlap 区間・charge witness の spectral gap・J の zero-crossing なし条件・
interval cost — exact reader (v33 器械) を信頼集合意味論 (v34.3) の上に持ち上げ、
dense/structured の裁定一致を検査する。
