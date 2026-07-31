# QRN v31.4 — 不変 operator response atlas / factorization no-go

**Version**: v31.4 (2026-07-31)
**Sim**: `sim/src/bin/v314_operator_atlas.rs` → `results/v314_operator_atlas.txt`
(8 検査 PASS)
**位置づけ**: PROMPT/12 第五版 — B3/B4 を統一する不変量 atlas と、E3 (因子分解発見)
の no-go 定理化。

---

## 1. Part A — whitened operator atlas

各ノードの parity-even 局所演算子基底 O_i^a = c†M^a c (2 軌道ノード: n₁, n₂, Re hop —
実対称 quadratic 空間は 3 次元で完備) に対し:

- 静的: X_ij^{ab} = ⟨O_i^a O_j^b⟩_c = **Tr(M^a (I−C) M^b C)** (Wick — 4 重ループ vs
  行列積の二重実装 5.6e-17 [A0])
- 動的: R_ij^{ab}(t) = i⟨[O_i^a(t), O_j^b]⟩ = i Tr(C [M^a(t), M^b]) (短時間展開との
  二重経路 1.9e-11 [A2])
- 局所 Gram G_i = X_ii → whitening X̂_ij = G_i^{−1/2} X_ij G_j^{−1/2}

**不変性 [A1]**: 特異スペクトル・作用素/Frobenius/核ノルムは局所演算子基底の
**任意の可逆再結合 L_i** (unitary に限らない) で不変 (7.8e-15)。時間積分ノルム
∫‖R̂‖dt も同様 (3.3e-16 [A2])。

**B3/B4 の位置づけ [A3]**: B3-COV は atlas の (n,n) 成分そのもの (|X^{nn}| = C² 厳密
1.4e-17)、B4 は動的 (n,n) 成分の到着時刻圧縮 (Spearman 0.991)。**「B3 と B4 の統一」が
atlas の成分/圧縮として実現した** (PROMPT/11 第二課題の残りの解)。

**rank 欠損 [A4]**: 従属演算子を含む基底では `ObservableSupportCertificate
{rank 3, nullspace 1}` を返し、支持制限 whitening の不変量が清浄基底と厳密一致
(0.0e0) — 無条件擬似逆行列は不使用。

## 2. Part B — factorization no-go (E3 の定理化)

**定理 (機械実例つき)**: Gaussian 状態の静的 one-body 共分散だけからは spatial
factorization を一意選択できない。

- 同一の大域状態 (ring12 熱的 Gaussian) が、因子分解の選び方で 3 つの異なる幾何を
  返す [A6]: **site 基底 → ring 12 辺 / eigenmode 基底 → 幾何なし (モード間相関
  厳密 0 = 4.3e-16 [A5]) / pair 回転基底 → 別の自己整合幾何 (24 辺)** — どれも
  「同じ状態の正しい読み」。
- **負制御 [A5]**: 自然な state-only 選択基準 (親生成子 K の疎性) は mode 基底
  (nnz 12 < site 24) を選ぶ — **最も疎な基底は幾何を自明化する**。state-only で
  factorization が「読めた」場合は hidden basis convention の流入を疑うべきことの
  機械実例。
- **OperationalAlgebra が幾何を選ぶ [A7]**: site-local probe (v31.2 曲率則) の応答は
  ring を返し (2.2e-16)、mode-local probe の応答は厳密 0 (1.3e-30) — **因子分解は
  状態ではなく操作代数 (準備・介入・測定・両立性) が運ぶ**。
  `OperationalAlgebra {preparations, interventions, measurements, compatibility}` を
  contract 型として実演。

v29.5 [C5] の「factorization は選定不能」という空隙が、対角化論法 + 機械実例 +
負制御で**定理化**された。FactorizationGivenObservables 能力の登録には状態以外の
操作的入力が必須 — RelationalDecompositionGoal (v31.0 schema) の成立条件が確定。

## 3. 正直な残高

- atlas の実対称枠: Im hop 演算子 (複素 M) は実 Wick 枠の外 — 複素側の gauge 共変性は
  v31.2 [L3] が担保。一般の複素 atlas は未実装。
- no-go は「静的 one-body 共分散から一意選定不可」の言明 — 高次相関・非 Gaussian
  資源・環境結合による選定は未走査 (open)。
- pair 回転基底の幾何 (24 辺) は「読める」が site 幾何と等価ではない — どの
  factorization が「物理的」かは E4 (自然との橋) の問い。
- 開発記録: 初版の 4 本目演算子 n₁−n₂ は span{n₁, n₂} 内 (Gram 特異で assert 発火)。
  2 サイトの実対称 quadratic は 3 次元 — 基底の数え間違いも Gram の正定 assert が
  捕捉した (支持証明書経路の存在意義)。
