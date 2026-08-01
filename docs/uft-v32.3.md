# QRN v32.3 — Factorization recovery: marked operational recovery と三値裁定

**Version**: v32.3 (2026-08-01)
**Sim**: `sim/src/bin/v323_factorization_recovery.rs` →
`results/v323_factorization_recovery.txt` (7 検査 PASS)
**位置づけ**: PROMPT/13 第三版 — v32.2 の erasure no-go (閉包は marking を消す) の
**裏返しの正側**。E3 (因子分解の操作的選定) が v31.4 の no-go 以来はじめて
**構成的な復元器 + 三値裁定**を持つ。

---

## 1. 目標定理 B — marked operational recovery

> primitive generator の非可換グラフ (α∼β ⟺ ‖[O_α,O_β]‖ > 0, 証明書つき) について、
> (1) 異なる未知ノードの局所操作は可換 (2) 各ノード内の非可換グラフは連結
> (3) 各連結成分の *-閉包は full matrix factor M_{d_i} (4) 全成分の join が対象
> sector の B(H) を生成 (5) 非可換マージン > 観測誤差 — のとき、連結成分と閉包から
> **H ≅ ⊗_i H_i を局所 unitary × ノード置換まで**回復できる。中心が非自明なら
> **tensor product を強制しない** — superselection sector
> H ≅ ⊕_α (C^{m_α} ⊗ C^{n_α}) を返す。

### 凍結決定手順 (recover_factorization)

```text
1. 非可換グラフの連結成分 — Abstain 対 (区間が閾値を跨ぐ) があれば
   Abstain(CommutatorMarginStraddled)
2. joint *-閉包 (HS Gram–Schmidt 成長)
3. joint が可換 → Abstain(InsufficientOperationalGenerators)
4. 中心 = {M ∈ 閉包 : [M, gens] = 0} (実 Gram 零空間・エルミート ONB)
5. 中心自明 (dim 1):
     full (dim n²) ∧ 各成分 factor (dim d_i², 中心 1) ∧ Π d_i = n
       → ExactUpToLocalUnitaryAndPermutation { local_dims }
     full でない → Abstain(InsufficientOperationalGenerators)  (未 address 自由度)
6. 中心非自明 (dim c ≥ 2): 中心射影 (√(k+2) 重みの生成元 T の Lagrange 補間 —
   冪等・ΣP = I を検証) で sector 分割 → 各 sector の制限代数 dim = m_α² ∧
   n_α = dim_α/m_α 整数 → SuperselectionSectors { (m_α, n_α) }
     整数条件の破れ → Abstain(ComponentNotFactor)
     (有限次元 *-代数は Wedderburn 分解を持つため、これは数値縮退の guard)
```

## 2. 結果 (7 検査 PASS — 全て operational_net::FactorizationReading で裁定)

| 検査 | セル | 裁定 |
|---|---|---|
| [F0] | site net / DFT₈ 共役 net ((C²)⊗³) | **Exact [2,2,2]** (両方 — 同一閉包・別 marking の各々が正しく復元) |
| [F0] | qutrit×qubit net (C⁶) | **Exact [2,3]** (異次元因子・Π d_i = n) |
| [F1] | number operator のみ {Z₁,Z₂,Z₃} | **Abstain(InsufficientOperationalGenerators)** (可換 joint) |
| [F1] | 部分 net {X₁,Z₁} | 同上 (中心自明だが full でない = 未 address 自由度で tensor を主張しない) |
| [F2] | 電荷つき C² ⊕ C²⊗C² | **SuperselectionSectors [(2,1),(2,2)]** |
| [F2] | 部分 address {X₁,Z₁,Z₂} | **[(2,2),(2,2)]** — 測定しかできない軸 (Z₂) は超選択ラベルに・未 address (qubit₃) は多重度 n_α に |
| [F2] | 対照 {X₁, Z₁Z₂} | Abstain(Insufficient) — 生成代数は dim 4 の因子 (非局所符号化の 1 qubit・中心自明)。超選択と誤読しない |
| [F3] | 局所 unitary × SWAP₁₃ 共役 net | Exact + **同一 gauge orbit** (成分部分代数 matching overlap = 1.000000000000) |
| [F3] | site net × DFT net | matching 不在 (best min-overlap 0.5625) → **EquivalenceClassOnly** |
| [F4] | 可換子証明書 σ = 1e-6 / 5e-4 | Exact / **Abstain(CommutatorMarginStraddled)** (辺の強制なし — 仮定 5 の機械化) |
| [F5] | parity-even フェルミオン net (Majorana 双線形 path) | **SuperselectionSectors [(4,1),(4,1)] = パリティ超選択の機械発見** |

### [F3] gauge orbit の裁定 — fiber の三値

因子分解の同値 (局所 unitary × 置換) は「成分部分代数 (traceless ONB) の集合の一致」
と厳密に等価 — matching する置換があれば**一つの許容 gauge orbit** (Exact)、
なければ **EquivalenceClassOnly** (完全だが非互換な 2 つの marking から無制約の
tie-breaker で 1 つを選ばない — v31.4 疎性負制御の教訓の実装)。事前登録された動的
functional による fiber 内 tie-break は将来課題として**意図的に不実装**。

### [F5] パリティ超選択の機械発見

Majorana 双線形 path {iγ_kγ_{k+1}} (全て parity-even — Ordinary lane が受理する唯一の
フェルミオン入口) の生成代数 = 偶部分代数 (dim 32)・中心 = span{I, Γ = Z₁Z₂Z₃}
(照合残差 < 1e-8)。復元は Fock 空間の tensor 因子ではなく **[(4,1),(4,1)]** —
**フェルミオンのパリティ超選択則が、操作的復元器の出力として自動的に現れた**。
odd primitive は v32.2 の型ゲートが Ordinary lane から排除済み (JW 誤読の遮断)。

## 3. 正直な残高

- toy (dim ≤ 8)・可換子は exact ノルム ± ノイズ模型の区間。ノイズは可換子ノルム
  読みへの加法 Gauss (σ, 6σ 区間) — 実測定の系統誤差モデルは HOLD-8 の変成セルで。
- [F2] 対照セルが示すとおり、「部分 address」には**中心が立つ場合 (超選択)** と
  **立たない場合 (非局所符号化の因子 — Insufficient)** の両方がある — 復元器は
  これを中心次元で機械的に区別する (推測しない)。
- graded lane の recovery (odd 演算子からのモード構造復元) は未実装 — 現行資格は
  parity-even lane。奇演算子の情報は偶双線形経由でのみ入る。
- 復元は「与えられた primitive 集合」に対する読み出し — **どの operations が
  physically accessible かは依然入力** (E3 の完全解決ではなく、E3 の
  operational fiber の機械化)。bridge law 登録は行わない。
