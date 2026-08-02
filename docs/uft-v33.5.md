# QRN v33.5 — Graded recovery の正しい境界: Majorana locality ≠ Dirac locality

**Version**: v33.5 (2026-08-02)
**Sim**: `sim/src/bin/v335_graded_recovery.rs` → `results/v335_graded_recovery.txt`
(7 検査 PASS) / 共有契約 = `sim/src/graded_recovery.rs` (新設)
**位置づけ**: PROMPT/14 第五版。graded lane の課題を「odd 演算子からモード構造を
復元」と置くのは誤り — **odd CAR だけから通常の複素 fermion mode を一意に得る
ことはできない**。まず no-go を立て、witness の門で正側を立てる:

> **Majorana locality と Dirac-mode locality は同じ識別問題ではない。後者は追加の
> U(1) charge / complex structure witness を必要とする。**

---

## 1. odd CAR は O(2N) 不変 — pairing の情報は無い [M1]

2N 本の Majorana ({γ_a, γ_b} = 2δ_ab) の CAR は実直交変換 O(2N) に不変である。
機械実証 (JW 6 本, dim 8): CAR 資格は原 frame・Givens 回転 (γ₂,γ₃ 混合) frame で
完全一致。graded bracket ノルムは全対**厳密 0** (空グラフ) — ordinary 可換子では
‖[γ,γ']‖ = 2√8 の K₆ に見える (v32.2 [N5] の罠の再確認)。**graded graph をいくら
読んでも「どの二本が一組か」は決まらない。**

## 2. Dirac pairing no-go (禁止変換 20) [M2]

witness なしの読みは **MajoranaFrameOnly (O(2N) orbit)** で止まる。標準 pairing
(12)(34)(56) と Givens 回転 pairing はともに**完全な mode-CAR** (偏差 ≤ 2.2e-16)
を満たしながら異なるモード (n̂₁ の差 0.97) を与える — CAR データの関数では選べない。
`MajoranaFrame → ComplexModeFactorization` の witness なし変換は存在しない
(**禁止変換 20**)。

## 3. 正側 — charge witness から複素構造を抽出する [M3]

Majorana の実 span 上の**直交複素構造 J (J² = −I)** が pairing を定める (複素構造が
fermionic creation/annihilation 表現を定めるのは自己双対 CAR 形式の標準的な数学 —
本版の寄与は識別可能性 fiber としての型化・棄却・機械検証)。抽出 (凍結):

```text
J_{ba} = ⟨γ_b, i[Q, γ_a]⟩ / ‖γ‖²   (charge witness Q の adjoint 作用の frame 展開)
資格: frame 外漏れ (線形性) ≤ 1e-9・‖J + Jᵀ‖・‖J² + I‖ ≤ 1e-9
```

Q = Σ n_i で J は実・反対称・**J² = −I 残差 4.4e-16**。J の不変平面 {v, Jv} から
a_i = (γ(v) + iγ(Jv))/2 を構成 — **3 モードの mode-CAR 6.3e-16・Σ â†â = Q 残差
1.2e-15** (回転 frame 座標でも 2.2e-15 — physical content は frame 選択に依らず、
残る自由度は U(N) gauge)。**witness が O(2N) orbit から pairing を選ぶ。**

## 4. 縮退・非線形 witness → Abstain [M4]

- **部分 charge n₁**: J'² = diag(−1,−1,0,…) ≠ −I — 残り 4 本の pairing が決まらず
  `Abstain(ComplexStructureUnresolved)`。
- **quartic 汚染 Q + 0.3 γ₁γ₂γ₃γ₄**: adjoint 作用が frame の実 span の外 (cubic
  Majorana) へ漏れ `Abstain(WitnessNotLinearOnFrame)`。
- **微小汚染 (1e-12)**: バー内で資格 (区間規律 — 拒否は閾値超過のみ)。

## 5. 既存復元器は捏造しない [M5]

- ordinary net は odd を**構成時拒否** (v32.2 ゲートの継承 — JW 誤読の遮断)。
- graded net (反可換子証明書 — 全対 0) の marked recovery (v33.1 入口) は成分が
  全て単本・成分閉包 dim 2 (= factor でない) → **Abstain(ComponentNotFactor)**。
  モード構造は graded graph から出ない — witness 経路 (§3) が唯一の門。

## 6. 型契約の登録 [M6]

- `sim/src/graded_recovery.rs` — `MajoranaFrame` (CAR 構成時資格 + O(2N) 回転) /
  `ComplexStructureWitness` (唯一の構成 = extract_complex_structure) /
  `GradedRecoveryReading` 凍結階層 (MajoranaFrameOnly / ComplexModeFactorization /
  Abstain 2 理由)。
- `core.schema.yml` に概念 3 種 + **禁止変換 20** を登録。

## 7. 正直な残高

- witness は偶・二次 (quadratic) 域の charge を想定 — 一般の偶演算子からの J 推定
  (最良近似複素構造 + 誤差区間) は未構成。ノイズ下の J 資格 (Straddled) は HOLD-9
  変成セルの主題。
- pairing 後の**モード局所性** (どのモードがどの空間サイトか) は別問題 — それは
  v33.1–v33.4 の control/context 機構が担う (graded 版の chart glue は未実装)。
- BCS 型 (数非保存) の witness — J を定めるのは U(1) charge に限らない (ギャップ
  つき quadratic Hamiltonian も J を定める) — は次期候補。
- bridge law 登録簿は全能力で空・PRED-019 未登録・自然の的中 0・
  `external_replications = 0` のまま — Unit D2-R の公募が引き続き最優先。
