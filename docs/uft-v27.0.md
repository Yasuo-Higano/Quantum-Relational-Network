# QRN v27.0 — full 4D kernel と dynamic metric fork (spec §13)

**位置づけ**: 第二十七期の最終 arc。事前登録 = spec §13 (b46bce3, 改訂 3 —
実装前凍結): (13.1) lapse-shift 型 temporal 結合 (N·ε + Nⁱ·πᵢ_Bel,
λ = −1/8)、(13.2) full 4D kernel の厳密 Ward と連続 universality、
(13.3) **dynamic metric fork の判断プロトコル — 1/Π 解禁 3 条件・三分岐・
Weinberg–Witten 破れ仮定の事前明記義務**、(13.4) ユニット計画 (0/A/B/C/D —
D = 期統合)。前提: v26.9 arc 完了 (動的 Ward 144 恒等式・Belinfante 完成・
massless σ = ρ₂P₂ 完全崩壊)。

## A. v270a_kernel4d — full 4D kernel の厳密 Ward (5 検査 PASS)

**構成** (v26.9-B の恒等式を核の定義に転化):

  k̂^{0ν,B}(iq₀,q) := C_{T⁰ν_Bel, B}(iq₀)
  k̂^{yν,B}(iq₀,q) := C_{J_ν, B}(iq₀) − ⟨[T⁰ν_Bel(q), B(−q)]⟩ / q̂

(T⁰ν_Bel = λ = −1/8 改良込み Belinfante 密度・J_ν = 厳密流束・第 2 項 =
接触完備化)。恒等式 iq₀C − q̂C_J = −⟨[A,B]⟩ により **iq₀k̂^{0ν} − q̂k̂^{yν} = 0
がカットオフ有限のまま構成的に厳密** — v26.6 静的核 (17 検査) の 4D 拡張。

| 検査 | 結果 |
|---|---|
| A0 接触完備化 X_ν = −⟨[T⁰ν,B]⟩/q̂ の q → 0 正則性 | 成長比 1.000 (非発散) |
| **A1 厳密 4D Ward (4 行 × 10 列 × 2 周波数, 12³)** | **6.8e-14** |
| A2 Onsager 対称性 (反対称改良項の微小破れ込み) | 7.5e-9 |
| **A3 流束行 → Belinfante stress 行 (全 4 行)** | O(ε²): 例 ν=x 1.5e-2 → 9.5e-4 |
| A4 変異 (λ → 0) | ν=x 行 1.30 停留 (正版の 1370 倍) |

A3 が核の**対称テンソル整合性**: Ward-厳密な流束行が連続極限で Belinfante
stress 行に一致 → 核は連続極限で「対称・保存・Ward-厳密」の三条件を同時に
満たす。開発記録: (i) A0 run1 は「定数収束」を誤要求 — 交換子が q³ で消える
対 (X ∝ q² → 0) も正則 (正しい判定 = 非発散)。(ii) Onsager は反対称改良頂点
により厳密には成立しない (破れ ~1e-8 — 高次微小)。

**型名**: spec §13.2-A の Ward 要件は充足。ただし
`FullGravitationalVacuumPolarization` の発行は §13.3 の全条件 (v27.0-B の
連続 universality + fork 予言の事前登録) 通過まで保留。**1/Π 禁止維持**。

## B. v270b_universality — 温度セクターの連続 universality (5 検査 PASS)

P₂ 幾何 + Lean 証明済みチャネル辞書から、静的 (q_L = (0, Qŷ)) kernel の
q⁴ln q² 係数は massless で P₂ form factor A₂ 一つに決まる: θ₀₀ = 1・
θ_yy = θ₀y = 0 より **A_00/A_D = P₂_{00,00} = 2/3 (厳密)・A_0y = A_yy = 0**。
null 結合 ladder (§12.9 凍結プロトコル・導出モデル {1, a²ln(1/a), a²}) を
channels {D, 00, 0y, yy} に適用:

| 検査 | 結果 |
|---|---|
| B0 A_D 回帰 (a = 0.125 vs v268p 公表) | 1.022148 vs 1.022197 (5e-5) |
| **B1 A_00/(⅔·A_D)** | **0.9974** (全域 0.9975, spread 0.0000) |
| B2 \|A_0y/A_D\| | **0.0021** (時間縦 = 純ゲージ) |
| B3 \|A_yy/A_D\| | **0.0024** (L = yy = 純ゲージ — Lean 定理の one-loop 版) |
| B4 変異 (V₀₀ ×1.02) | 1.0377 (ゲート外 ✓) |

**温度セクター (h₀₀/h₀y source) の one-loop 応答が空間セクターと同一の
繰り込み後 form factor に流れる** — spec §13.3 の 1/Π 解禁条件 (ii) 充足。
縦・時間チャネルの消滅は「純ゲージチャネルに独立な物理自由度なし」の
one-loop 実証 (v26.6 の bare c₁ 汚染 [L チャネル同桁] と対照的 — 汚染は
繰り込みで消える scheme 量だったことの最終確認)。

## C. v270c_fork — dynamic metric fork 判定の執行 (5 検査 PASS)

§13.3 の凍結プロトコル執行 — 解禁条件 (i)(ii) 充足を受け、**1/Π をこの監査に
限って計算**。繰り込み後 Π₂^ren(q²) := χ_D(q) − c₀ − c₂q² − c₄q⁴ (A₂ =
2A_oracle 凍結、a = 0.045/0.09, Q ∈ [0.15, 0.9]):

| 検査 | 結果 |
|---|---|
| F0 A₂ 凍結フィット残差 | **1.7e-12** (凍結 A₂ がデータを厳密記述) |
| F1 形状 Π₂^ren/(A₂Q⁴lnQ²) | ∈ [0.986, 1.014] (零交差なし) |
| **F2 no-pole** \|Π₂^ren/Q²\| → 0 | 縮小比 0.235 vs Q²lnQ² 導出値 **0.232** |
| **F3 Sakharov 走行** c₂(0.045)/c₂(0.09) | **4.000** (= a⁻² ちょうど) |
| F4 変異 (A₂ ×2) | 残差 380 倍 |

**fork 執行 (確定)**: **分岐 (b) external metric** —
1. 自由場 matter loop の普遍部分は A₂q⁴ln q² のみ — **1/Π₂ に有限留数の
   massless pole なし (graviton は生成されない)**。
2. pole を作りうる q² 項 (Sakharov 誘導 Newton 定数) は bare c₂ で **a⁻² で
   走る regulator 量** — 本プログラム自身の universality 基準で非普遍。
   v26.3 以来の「bare c₁/c₂ は scheme 量」の主張の最終的・定量的実証 (比 4.000)。
3. composite graviton 路線 (分岐 a) は「相互作用による普遍 q² 項の生成 +
   Weinberg–Witten 破れ仮定の明示」を要件とする**将来の別プログラムに封印**。
4. **1/Π の常用は解禁しない** (本監査限り)。metric は外部 regulator のまま —
   経路 B は「機械証明・source matching・二離散化・4D Ward を備えた Dirac
   stress-tensor vacuum-polarization audit」として完結し、QRN-Core v1
   (matter-on-background) へ接続する。

開発記録: F2 run1 のバー 0.2 は log 増強を落とした誤較正 — 導出値 0.232 との
5% 照合に較正 (測定 0.235 は 0.9% 一致 — 判定がむしろ強くなった)。

## D. 期統合 — 第二十七期の総括と QRN-Core v1 の定義要件

### D.1 第二十七期の総括 (v26.1 → v27.0-C)

**期のテーゼ: 「予言の台帳は赤字でも、測定器は黒字にできる」**。

- **v26.1 (flavor prequential)**: 凍結規則の台帳で S_program − S_FN = −55.9
  — 予測プログラムとしての flavor 統一は**赤字のまま公表** (D_J 支配)。
- **v26.2–26.5 (静的応答監査)**: BOND-A の 10 要件・c₁ 行列 (縦 +0.016/
  横 −0.013/spin-2 −0.014)・bare c₁ は非共変 regulator 汚染に支配 (縦 =
  spin-2 同桁)。
- **v26.6–26.7 (核・スペクトル)**: 完全静的核 (17 検査)・Lean BR 射影 13+14
  定理・**PRED-013 的中 (プログラム初の凍結新規予言: 相互作用が q⁴ 保護を
  破る)**・自由場 no-pole (PRED-014)。
- **v26.8 (連続極限 universality — 最重要 falsifier)**: 解析 oracle
  −1/(160π²) 三重導出・**4 比 (D/X × stag/Wil) = 1 ± 0.14% (PRED-016)・
  スカラー和則 (PRED-017)・三者一致 (spectral 経路, a² Lorentz 回復)** —
  全て scored-hit。
- **v26.9 (4D covariance closure)**: 動的 Ward 144 恒等式 (機械精度)・
  **Belinfante 構造の格子分解** (正準 = piece 単独・改良 λ = −1/8 は su(2)
  代数の予言どおり)・massless σ = ρ₂P₂ 完全崩壊・負の結果 1 件凍結。
- **v27.0 (fork)**: 厳密 4D kernel (Ward 6.8e-14)・温度セクター universality
  (A_00/A_D = 2/3 的中・純ゲージ消滅)・**分岐 (b) external metric 確定** —
  graviton pole なし・誘導 Newton 定数候補は a⁻² 走行の regulator 量。

**正直な残高 (不変の分類)**: 的中 = 計算実験内の機構予言 (PRED-013–018)。
**自然の観測量の的中 0・独立外部再現 0**。経路 B の成果は「測定器が正しい」
ことの証明であり、QRN・創発重力の証拠ではない (spec §12.8/§13.3 凍結解釈)。

### D.2 未解決問題の残高 (第二十八期への持ち越し)

1. **QRN-Core v1 の定義** (モラトリアム解除条件 — D.3 の要件で着手)。
2. 相互作用系の q² 項生成の有無 (v26.7 PRED-013 の続き — composite 路線の
   封印解除条件でもある)。
3. flavor 台帳の赤字 (前方予言 τ = 1/12 + i/2 は登録済み — 待ち)。
4. massive ρ₀ (trace 復活) の 4D 分解・Wilson 温度チャネル。
5. 外部: anomaly-search / modular-BW の投稿・第三者再実装の公募 (最優先
   マイルストーン — PROMPT/9)。

### D.3 QRN-Core v1 の定義要件 (登録 — 定義完了までモラトリアム維持)

QRN-Core v1 = 「この期までに**認証された器械と公理だけ**で書かれた、
matter-on-background の最小核」。定義完了の要件:

1. **公理系**: 格子 = regulator (存在論)・BOND-A 結合則・置換則 (中点変調)・
   Belinfante 改良 (λ = −1/8) — 全て機械検査可能な形で 1 文書に凍結。
2. **器械スイート**: 認証済み測定器の一覧 (null 結合・殻積分・Matsubara
   Ward・接触項 2 実装・導出モデル外挿) と各々の較正記録 (開発記録の教訓集)。
3. **主張の等級**: C0–C5 台帳 (claims.yml) の継続 + 「自然の観測量 0」の
   明示を維持。
4. **falsifier**: 4 比 universality・和則・Ward 恒等式群を常設回帰化
   (suite の監査層)。

### D.4 スイート (2026-07-28)

増分スイート (`make suite`, 台帳 = results/suite_manifest.tsv):
**実行 33 本 PASS 268 / FAIL 0・引用 129 本 PASS 856 / FAIL 0 (ソース不変を
sha256 確認)・総計 PASS 1124 / FAIL 0** — results/v270_suite.txt。
rustc 1.94.0。v252_manifest (v25.2 凍結台帳) を含む監査層は全て健全。

## 1. 残り

- なし (v27.0 arc 完結)。第二十八期 = QRN-Core v1 の定義 (D.3 の要件) から。

## 2. 成果物

A: `sim/src/bin/v270a_kernel4d.rs` / `results/v270a_kernel4d.txt`
(5 検査 PASS) / `results/v270a_kernel4d.json`。claims: QRN-GRAV-057。
B: `sim/src/bin/v270b_universality.rs` / `results/v270b_universality.txt`
(5 検査 PASS) / `results/v270b_universality.json`。claims: QRN-GRAV-058。
C: `sim/src/bin/v270c_fork.rs` / `results/v270c_fork.txt`
(5 検査 PASS) / `results/v270c_fork.json`。claims: QRN-GRAV-059。
D: `results/v270_suite.txt` (総計 PASS 1124 / FAIL 0)・本文書 §D
(期統合)・CLAUDE.md 到達点更新。claims: QRN-META-028。
