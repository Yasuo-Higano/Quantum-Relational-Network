# 量子情報網理論 v26.6 — 重力真空偏極の監査 I: 接触項込みの完全核と「縦チャネルは純ゲージ」

**第二十七期第 6 版。判定 (a) — 17 検査 PASS (454 s)。**
本版から誘導重力経路 B は PROMPT/8 の凍結仕様
[paper/grav-vacuum-polarization-spec.md](../paper/grav-vacuum-polarization-spec.md)
(コミット 157ca53 — **実装より先にコミット**) に従う。仕様が凍結したもの:
**格子の存在論 = regulator** (bare c₁ は scheme 量、比較対象は繰り込み後 form factor)、
結合則 scheme BOND-A、検査 S0–S9 の数値ゲート、**禁止暗黙変換 6 種**
(bond modulation ≠ vierbein / plus ≠ spin-2 / 質量依存 ≠ regulator 依存 /
1/χ ≠ graviton propagator / q² 係数 ≠ Newton 定数 / T00 保存 ≠ 保存 Tμν)、
要件 0 (背景停留) の新設と要件 10 の書き換え。v26.5 までの正名も凍結した:

> v26.5 は「誘導重力を発見した段階」ではなく、自由フェルミオン格子に対する外部計量
> 応答の測定器を較正し、素朴な Newton 定数解釈を棄却した段階である。

## 1. 器械 — 完全静的核 (接触項・counterterm 込み)

結合則 (spec §2): 方向 i のボンド振幅 t → t·(1 + h_ii(中点))^{−1/2} (質量項は正準規格化
で h 非結合)。一次頂点 V_A = −½T_A(q)、**二次接触項 (seagull) S_AB = (3/4)δ_dir Σ c_A c_B
t(c†c+h.c.)**、counterterm Γ_ren = E − Λ Σ_x √g (Λ = −e_y)。完全静的核:

> k̂_AB(q) = (3/4) e_A δ_AB − ¼ χ_AB(q) − (Λ/4)(J − 2I)_AB (A,B ∈ {xx,yy,zz}, q ∥ ŷ)

**射影代数は Lean で証明してから使う** ([proofs/Projector.lean](../proofs/Projector.lean),
定理 13 本 — 一発通過): Barnes–Rivers 射影 (整数化 m = 2q⁴P) の完全性・冪等・直交・
trace=rank・Ward 収縮・ゲージ消去/再現を格子 {0..8}³ + 次数勘定で機械検証し、
ŷ 方向のチャネル辞書 **D = (T_xx−T_zz)/√2 と X = T_xz が spin-2、S = (T_xx+T_zz)/√2 が
P0s、L = T_yy が P0w = 純ゲージ (2E_yy = q⊗ŷ + ŷ⊗q)** を定理化した。

器械認証 (全て機械精度):
- [S0] **FD Hessian 完全性**: exact ボンド因子の H[ε] を再対角化した d²E/dh² が
  ⟨S⟩ − χ_V と一致 — 対角 5 + 交差 3 モード × m 2 種で worst 0.005×許容
  (相対 ~5e-8)。**接触項 + connected = 完全な二次変分**の端から端の証明。
- [S1] block 3×3 χ 行列 = dense **5.1e-14** / [S1b] cos モード写像 2χ^TT_cos = χ^complex
  **8.0e-16** / [S0c] cos/sin 等値 (y 並進の端から端)。
- [S2] tadpole: block/dense = 解析 k 和 **6.5e-14**、e_y = e_z **1.6e-15**。
- [S3a] **背景停留 (要件 0)**: 単一の Λ = −e_y で |dΓ_ren/dε| = **7.4e-10** (uni-yy/zz)。
  [S3b] x 残差 Δ_x(N) = |e_x/e_y−1| は m=0.5 で 4.3e-3 → 1.5e-5 → 1.6e-9 → **1.8e-15**
  (N = 8→64; m=0 は 1.3e-2 → 2.3e-6) — twist 圧力異方性は有限サイズで、gapped は
  指数的に消える。
- [S4] y 並進不変性 |⟨T_i(2qŷ)⟩|/V = **1.8e-16** — **接触項が q 非依存 ((3/8)e_i V)
  であることの厳密な根拠**。ゆえにこの scheme では接触項は c₁ (q² 係数) を動かさない —
  v26.3–26.5 の c₁ 結論は接触完備化後もそのまま成立する ([S5] 回帰 max|Δ| = 4.1e-6)。

## 2. 主結果 — 縦チャネルの c₁ は純ゲージ汚染で、spin-2 と同桁

q ∥ ŷ の対角ソース空間 (xx, yy, zz) は Lean 辞書で D (spin-2) ⊕ S (P0s) ⊕ L (P0w) に
分解される。**L = yy はゲージモードそのもの**なので、diffeo 不変な連続理論では
(停留背景のまわりで) 完全核の縦列 **K_LL, K_SL, K_DL は恒等的に零**でなければ
ならない。格子の実測 (N=64, χ 単位; k̂ の係数は ×(−¼)):

| チャネル (Lean 辞書) | c₁^χ (m=0) | c₁^χ (m=0.5) | 質量走行 | 連続極限の要求 |
|---|---|---|---|---|
| D = spin-2 plus | −0.01935 (±0.00014) | −0.01444 (±0.00001) | 25.4% | 物理候補 (要繰り込み) |
| S = P0s 横トレース | −0.01177 (±0.00004) | −0.01083 (±0.00000) | 8.0% | 物理候補 (要繰り込み) |
| **L = yy (縦 = 純ゲージ)** | **+0.01833 (±0.00007)** | **+0.01595 (±0.00000)** | 13.0% | **恒等 0** |
| SL (P0s×P0w 混合) | +0.00187 (±0.00004) | +0.00234 (±0.00000) | 25.2% | **恒等 0** |

- **[S6 branch a] R(m) = |c₁^χ[L]| / |c₁^χ[D]| = 0.947 (m=0) / 1.105 (m=0.5)** —
  縦 (純ゲージ) チャネルの q² 係数は、重力を作るはずの spin-2 チャネルと**同じ桁で
  非零かつ質量とともに走る**。この量は連続極限で消えるべき純粋な非共変 regulator
  汚染だから、**同桁である事実は「bare c₁ の値には連続的意味がない」ことの機構的
  な定量化**である (v26.3–26.5 の「c₁ は bare 量」に、汚染の物差しがついた)。
- **q⁰ 側の汚染**: k̂_LL(0) = ½e_y − ¼χ_yy(0) = **−0.1256** (m=0, N=64;
  e_y = −0.19897, χ_yy(0) = +0.10460) ≠ 0 — 停留 (一点) を Λ で消しても、
  二点核の縦成分は q⁰ から汚染される。
- **帰結 (v26.8 の前提)**: 連続比較には少なくとも **3 種の繰り込み条件**
  (Λ [一点] / q⁰ 縦核 / q² 縦核) を課した後の量だけが意味を持つ。spec §6 の
  要件 10 書き換え (「同一繰り込み条件後の form factor と spectral density の一致。
  bare c₁ の一致は要求しない」) はこの実測に裏づけられた。
- [S7] スカラー混合 (新測定): c₁^χ[SL] = +0.00187 → +0.00234 (走行 25.2%)、
  χ_SL(0) = −0.0740 → −0.0590 — P0s×P0w 遷移核も縦列の一部として連続極限で
  消えるべき量。**対して K_DS, K_DL は m=0.5 で機械零 (~1e-16)、m=0 でも ~1e-7** —
  spin-2 とスカラー系の混合はバルクで消えており (x↔z 対称の核レベル実証)、
  汚染は縦列に局在する。
- [S8] uniform 連続性: block 一様変形の FD = 組立 (相対 4.5e-8)、
  |k̂(q₁)−k̂(0)| = 1.5e-4 ≪ 5% scale。

## 3. 開発記録 (器械の教訓)

1. **解析 tadpole の band-pairing 二重計上**: S2a の参照実装 (解析 k 和) が N³ 個の
   k 全てに負エネルギー枝を数えていた — band pairing (k ↔ k+π 系) で占有状態
   (N³/2) をちょうど 2 重計上する ×2 のバグ。**走行前の独立 python 照合**
   (dense e_y = −0.198 vs 解析 −0.396) が検出し、ゲート発火前に修正した。
   器械側 (dense/block) は無変更 — 参照側のバグ。「参照実装も器械と同格に疑え」。
2. run1 は S0 の逐次 FD 走行中に打ち切り (検査出力なし)、dense FD を 118 点の
   並列バッチ (14 threads, 決定的分割) に再構成した run2 (151 s) が公表 run。
   物理値の変更はない (run1 は S0 以降に到達していない)。

## 4. 次 (登録課題 — spec §7–9 で凍結済み)

- **v26.7**: 動的 spectral measure ρ_AB^(N)(s) (pole 判定 = residue の体積 scaling) +
  自由場 benchmark 5 項目。**副実験 v267b_q4break** (相互作用による χ_00 の q⁴→q²、
  [predictions.yml PRED-013](../predictions.yml) に登録済み — モデル交換禁止)。
  cross 偏極 X = T_xz の point-split 器械と taste 認証。
- **v26.8**: 明示的 a の連続極限 (固定 m_phys, q_phys, L_phys)・2-taste Dirac の
  解析 form factor 照合・**q⁴ln q² 係数**・Wilson fermion 独立離散化 —
  **経路 B の最重要 falsifier**。
- **v26.9**: 動的 metric の三分岐 (外部応答 / 誘導 EFT / composite)。composite には
  Weinberg–Witten の破る仮定の事前明記が必要 (spec §6)。

## 5. 成果物

`sim/src/bin/v266_vacuum_pol.rs` / `results/v266_vacuum_pol.txt` (17 検査 PASS) /
`results/v266_vacuum_pol.json` / `proofs/Projector.lean` (定理 13 本) /
`paper/grav-vacuum-polarization-spec.md` (凍結仕様) / predictions.yml PRED-013。
claims: QRN-GRAV-039。
