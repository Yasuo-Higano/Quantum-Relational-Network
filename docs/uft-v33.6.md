# QRN v33.6 — 構造化スケーリング: Pauli GF(2) / Majorana quadratic backend

**Version**: v33.6 (2026-08-02)
**Sim**: `sim/src/bin/v336_structured_backend.rs` → `results/v336_structured_backend.txt`
(6 検査 PASS) / 共有契約 = `sim/src/structured_backend.rs` (新設)
**位置づけ**: PROMPT/14 第六版。一般 dense 行列の *-閉包をそのまま高速化して
「大型系へスケールした」とは主張しない — qubit 数に対して行列次元が指数増大する。
lane を明示的に分け、各 lane は自分の scope でのみ語る。

---

## 1. 3 lane の分離

| lane | scope | 表現 | 大型入力 |
|---|---|---|---|
| GenericDenseBackend | 小次元の完全一般系 | 行列 (v33.1 復元器) | dim > 4096 は **ScopeExceeded** |
| PauliSymplecticBackend | Pauli/Clifford 型 | GF(2)^{2n} (x\|z) + symplectic ω | **行列なし**で 48 qubit の証明書 |
| MajoranaQuadraticBackend | quadratic fermion | 実反対称 2N×2N | **行列なし**で 48 Majorana |

## 2. Pauli backend — v33.1 凍結手順の GF(2) 実装 [S1][S2]

可換性 = symplectic 形式 ω(P,Q) = x_P·z_Q + z_P·x_Q・成分 = ω グラフ・閉包次元 =
2^{dim V} (V = GF(2) span)・**中心 = ω の radical**。資格 (文脈の存在・被覆)・証人
ゲートも v33.1 と同一に実装する。

- **n = 3 の全 7 セルで dense と裁定が完全一致**: site [2,2,2]・entangler [2,4]・
  number-only Insufficient・部分 net Insufficient・singleton 文脈 unwitnessed・
  部分 address 超選択 [(2,2),(2,2)] (radical = {Z₂})・パリティ超選択
  [(4,1),(4,1)] (radical = {Z₁Z₂Z₃}) — superselection の sector 構造まで GF(2)
  radical が正しく数える。
- **48 qubit の証明書 (壁時計バー ≤ 5 s)**: site 96 本 → Exact [2×48]・
  +entangler → [2×46, 4]・Z₄₈ 欠落 → SuperselectionSectors [(2⁴⁷, 1)×2] — 全て
  96×96 の GF(2) rank で決まり、**2^48 次元の行列はどこにも現れない**。
- 正直な限定: Pauli lane の「証明書」は宣言された Pauli 構造の**厳密代数データ**で
  あって測定ではない — ノイズ・区間証明書の lane は dense の領分 (交差させない)。

## 3. Majorana quadratic backend と対応原理 [S3]

表現 = 実反対称 A (H = (i/4)Σ A_ab γ_a γ_b)。ブロック = 支持グラフの成分・閉包 =
Lie (so) 閉包。**quadratic 閉包は full M_d を与えない** (偶 Clifford = パリティ
超選択の quadratic 版) — dense との対応原理は:

```text
支持分割 = dense 非可換成分・dense *-閉包の複素次元 = 2^{2m−1} (2m 本のブロック)
```

- 小 N = 3 で機械照合: blocks [2,4]・dense 閉包 8/2 = 予言 8/2 ✓。
- 大 N = 24 (48 Majorana): 3 ブロック × NN hopping 15 本 → **so(16) 閉包 (dim 120)
  × 3**・cross hop 1 本で [16,32] へ併合 (so(32) dim 496) — 2^24 は現れない。

## 4. scope 規律 — ScopeExceeded は正答 (禁止変換 21) [S4]

- dense: dim > 4096 (凍結バー) は ScopeExceeded — **試行すらしない**。
- Pauli: 非 Pauli 和 (X₁+X₂) は `PauliVector` に**構成不能** (from_dense が
  NotPauliString で拒否・真の Pauli 文字列は反可換パターンまで往復復元)。
- quadratic: 非反対称は資格外 (quartic は表現の型に存在しない)。
- **禁止変換 21**: GenericDenseSmallSystemResult → LargeSystemFactorizationClaim。
  大型の主張は structured lane の scope 内でのみ立ち、その資格は小系での dense
  との裁定完全一致 [S1] と対応原理 [S3] が与える。

## 5. 型契約の登録 [S5]

- `sim/src/structured_backend.rs` — `Gf2Vec`/`gf2_rank` (GF(2) 線形代数)・
  `PauliVector` (from_str / 資格つき from_dense)・`PauliNetSpec`/`recover_pauli_net`
  (v33.1 手順の GF(2) 実装 — 同じ `FactorizationReading` を返し dense と直接比較
  可能)・`QuadraticGenerator` (反対称資格)・`recover_quadratic_blocks`・
  `StructuredScopeError` 3 種・`DENSE_DIM_BAR = 4096`。
- `core.schema.yml` に概念 4 種 + **禁止変換 21** を登録。

## 6. 正直な残高

- Pauli lane は宣言された厳密 Pauli 構造の lane — 測定ノイズ下の資格 (区間・
  Straddled) は dense lane の領分で、両 lane の跨ぎ (ノイズつき大型) は未構成。
- quadratic lane の読みはブロック分割 + so 閉包まで — graded witness (v33.5 の J)
  との統合 (大型 quadratic 系の complex mode 回復) は次期候補。
- 48 qubit / 48 Majorana は「行列を生成しない」ことの実演であって物理主張ではない
  (physical_scope: toy のまま)。resource profile / atlas glue の structured lane 版
  は未実装。
- bridge law 登録簿は全能力で空・PRED-019 未登録・自然の的中 0・
  `external_replications = 0` のまま — Unit D2-R の公募が引き続き最優先。
