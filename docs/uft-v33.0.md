# QRN v33.0 — HOLD-8: 識別可能性境界の holdout と第三十二期統合

**Version**: v33.0-A (2026-08-01 — 凍結半)
**Sim**: `sim/src/bin/v330a_hold8_freeze.rs` → `results/v330a_hold8_freeze.txt`
(3 検査 PASS + train 20 セル満票)
**位置づけ**: PROMPT/13 の最終版 — 第三十二期の全器械 (v32.1 SupportNoiseCertificate /
v32.2 OperationalNet / v32.3 recovery / v32.4 応答階層 / v32.5 hypergraph) を凍結し、
**factorization × interaction order × observation contract** の識別可能性境界を
新鮮 holdout で採点する。

---

## 1. HOLD-8 が採点するもの (HOLD-7 の後継)

「正解グラフを何本当てたか」ではなく**識別可能性の境界**:

```text
selective risk        = 0   (回答セルの誤り 0 — 強制回答も誤りに数える)
impossibility recall  = 1   (非識別セルの正しい棄却/同値類)
answerable recall     = 1   (回答可能セルは全て答える — coverage 単独は最大化しない)
```

## 2. 開封順序 (HOLD-5/6/7 と同一)

```text
v33.0-A (本コミット) = 生成器・採点器・バー・棄却規則の凍結
                      + SECRET コミットメント公表 + train 採点 (可視シード 33001)
  → v33.0-B = SECRET 開示・holdout 初生成・本採点 (調整なし) + 期末完全儀式
sha256(SECRET) = fb05a4a07dd7feef78e0bcecdaeeff933a656d36a38ceb26c71047002621c582
```

凍結カーネル (FROZEN-HOLD8 区間 — v330b と逐語一致を [H0] が照合) は第三十二期の
器械の逐語コピー: OperationalNet 構築 (可換子の区間証明書, τ = 1e-3)・recovery
決定手順 (成分 → 閉包 → 中心 → 三値裁定)・R⁽¹⁾/R⁽²⁾ 応答核・条件付き接ベクトル・
最終 gap 則・SupportNoiseCertificate (σ·z 対 窓ガード)・hyperedge 検出バー
(条件差 > 1e-6·max)・coherent 分離バー (0.05, 電流 + 実 hopping の 2 チャネル)。

**器械訂正 (設計走行で発見・v32.3 kernel に統一適用)**: 共役射影の数値塵
(‖PbP‖ ≈ 0 の候補が正規化されて基底に混入し restricted 次元を膨らませる) を
dust guard (閾値 1e-9) で除外 — v323 の既存 8 検査は全て不変 (exact 射影では
no-op)。v31.6 の gap 抽出器訂正と同じ「器械契約への昇格」型の修正。

## 3. セルクラス (8 種 20 セル — 隠しパラメータは重み・角度・置換・V・θ・σ)

| クラス | 内容 | 要求 |
|---|---|---|
| F1/F2 ×2 | site net / mode (DFT) net (隠し局所 U × 置換) | **Exact [2,2,2]** + 正しい gauge orbit |
| F3 ×1 | number operator のみ | **Abstain(InsufficientOperationalGenerators)** |
| F4 ×1 | 完全だが非互換な 2 net | **EquivalenceClassOnly** (強制一致 = FAIL) |
| F5 ×1 | 中心非自明 {X_a, Z_a, Z_b} | **SuperselectionSectors [(2,2),(2,2)]** (tensor 強制禁止) |
| F6 ×1 | odd 入力 + 偶双線形 (隠し係数) | ordinary **構成時拒否** + graded 化で [(4,1),(4,1)] |
| I1–I5 ×6 | quadratic (ring/chain 隠し) / t-V / 相関 hopping / pair hopping / 三体 ± null | 支持 欠0余0・**ŵ = t² (V 不可視)**・hyperedge 検出 (null 対で不発)・条件付き K(v) = \|t+vV\|² |
| I6/I7 ×4 | H↔−H 対・磁束 ±θ 対 | **密度 lane = 同値裁定 (分離主張 = FAIL)**・coherent lane = 分離 |
| M1/M3 ×2 | 変成対 (局所 U×置換)・基底の可逆再結合 | 同一 orbit (overlap 1)・recovery 不変 |
| M4/M5 ×2 | 可換子 margin 以下のノイズ・support margin 以下の弱辺 | **Abstain** (Straddled / InsufficientObservation) |

**train (シード 33001, 可視) は 20 セル満票**: selective risk 0.000・
impossibility recall 1.000・answerable recall 1.000。回答 14/14 (F1/F2/F5/F6・
I1–I5・I6c/I7c・M1/M3)・正棄却 6/6 (F3/F4・I6d/I7d・M4/M5)。

*(v33.0-B の開封結果・期末儀式・期統合は開封後に本文書へ追記する — 調整なし)*
