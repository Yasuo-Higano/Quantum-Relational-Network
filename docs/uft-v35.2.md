# QRN v35.2 — 開放系 signed-response 観測商: GQF 定理群 (quasi-free lane)

**位置づけ**: PROMPT/16 (第三十五期「外部化と観測商」) §5–6・優先度 1-B。実データの
科学的解釈の前提となる open-system scope guard。期テーゼ「登録観測契約が識別する
のは Hamiltonian や jump operators そのものではなく、応答写像の核で割った観測商で
ある」の第一の定理化。

**一次ソース**: `results/v352_open_quotient.txt` (18 検査 PASS)・
`results/v352_lean.txt` (`proofs/OpenQuotient.lean` 16 定理・終了コード 0)。

---

## 1. lane の定義と規約の較正

有限モード・number-conserving・quasi-free (Markov) 開放系:

```text
ρ̇ = −i[H,ρ] + Σ_k D[L_k]ρ + Σ_l D[G_l]ρ
H = Σ h_ab c†_a c_b (h エルミート),  L_k = Σ_a ℓ_{k,a} c_a,  G_l = Σ_a g_{l,a} c†_a
```

の normal covariance C_ab = ⟨c†_b c_a⟩ は閉じた affine ODE に従う:

```text
Ċ = X C + C X† + Y,   X = −ih − ½(Λᵀ + M),   Y = M
Λ = Σ_k ℓ_k ℓ_k†,   M = Σ_l g_l g_l†
```

**開発記録 (規約較正 — [G1] の存在理由)**: 初版実装は loss を Λ (非転置) で入れて
いた。[G1] (dense Lindblad = 2^N 次元 Jordan–Wigner の RK4 積分との比較) が複素
振幅の loss で **5e-2 の不一致**を検出し、反交換子の縮約を機械で追い直して
**loss 側だけ転置 Λᵀ = conj(Λ) が出る**規約を確定した (gain 側 M は非転置)。
修正後は N=2 (loss+gain 複素) / N=3 (loss 複素) とも max|Δ| ≤ 3e-15。
covariance 閉包そのものを第一原理 (2^N 量子発展) から機械較正したことになる。

解は Van Loan block で厳密に評価する (`sim/src/open_response.rs`):
e^{[[X,Y],[0,−X†]]t} = [[F,G],[0,H]] → C(t) = F C₀ F† + G F†。

## 2. GQF 定理群 (Lean 16 定理 + 機械実証)

`proofs/OpenQuotient.lean` は行列恒等式を**全ての整数値で成立する一般定理**として
証明する (simp の積分配・可換正規化 + omega の原子化線形算術)。全整数点で成立する
多項式恒等式は係数恒等式であり (ℤⁿ の Zariski 稠密性)、ℝ/ℂ への拡張は標準論法 —
LocalBias.lean (v31.2) の格子 4 点法より強い形。反例対は整数インスタンスの decide。

| 定理 | 内容 | Lean | Rust |
|---|---|---|---|
| **GQF-1** | affine 消去: 符号つき差分から Y が恒等的に消え Ḋ = XD + DX† | `gqf1_affine_cancellation`・`gqf1_probe_difference` | [G2a] Y ≠ 0 で Δ(t) = 2εe^{Xt}Pe^{X†t} (4e-17)・[G2b] 生軌道は Y 依存 (負制御) |
| **GQF-2** | 曲率 (n̈⁺−n̈⁻)(0)/(4ε) = **‖P_j X P_i‖²_F** (有効 drift を読む) | `gqf2_block_frobenius`・`gqf2_first_derivative_zero`・`gqf2_curvature_2x2`・`gqf2_curvature_block` (2 軌道 node) | [G3a] jet 恒等式 = frob (厳密 0)・[G3b] Richardson 有限差分 1e-7 |
| **GQF-3** | **Hamiltonian 昇格 no-go**: coherent hopping (h≠0, ℓ=0) と collective loss (h=0, Λ=2(1,−1)(1,−1)ᵀ) が同一曲率 1 — 曲率統計から Hamiltonian support への写像は存在しない | `gqf3_same_curvature`・`gqf3_different_hamiltonian`・`gqf3_no_map` (∀f) | [C1] coherent w=1 昇格可・[C2] loss w=1 昇格拒否 |
| **GQF-4** | 正側の門: dissipator が node-block 対角なら cross block は −ih に厳密一致 — 曲率は ‖P_j h P_i‖² へ還元 (昇格の解錠) | `gqf4_block_diagonal_reduction`・`gqf4_neg_i_isometry`・`gqf4_scale_two`・`gqf4_scale_norm` | [C3] 対角 loss で閉鎖系値と差 0.0 |
| **GQF-5** | **charge attribution no-go**: 電荷応答値 −8 が loss (pairing なし) と BdG pairing (散逸なし・pure BCS 状態) で同一 — 値から pairing は選べない。正側: Δ = 0 の閉鎖 Nambu は電荷統計恒等 0 → 散逸ゼロ証明書下で ≠ 0 ⇒ Δ ≠ 0 | `gqf5_same_charge_response`・`gqf5_no_map`・`gqf5_closed_charge_conservation` | [G4a] 整数橋・[G4b] dense dN/dt×4 = 宣言統計・[C10] 証明書の門 |
| **GQF-6** | Robust open promotion: 有限 shot は同時信頼集合の全域で支持判定が一致するときだけ回答 | (v34.3 FiniteDataNoGo の Robust Promotion を継承) | [C8] RobustEdge/RobustNoEdge/Straddled/Insufficient の 4 裁定 |

**GQF-6 の登録推定器**: Richardson 4 点 K̂ = [8Δ̂(δ) − Δ̂(2δ)]/(8εδ²) は δ³ 項を厳密に
消し、残余は K̂ = w − 2c₄δ² (|c₄| ≤ (2/3)R⁴, R = 宣言 drift ノルム上界) — 区間 =
CP 伝播 + 登録バイアス (4/3)R⁴δ²。shot 系列はまず v35.1 の相関粒度ゲートを通す
([C9]: 持続 Markov shot は OutOfDomainCorrelated — 読まない)。

## 3. 型の規律 (禁止変換 32/33)

```text
OpenSignedCovarianceProbe → EffectiveDriftTopology        (曲率が読むのは X)
EffectiveDriftTopology   ↛ HamiltonianTopology            (禁止変換 32)
EffectiveDriftTopology + DissipatorLocalityCertificate → HamiltonianTopology
ChargeNonconservingResponse ↛ HamiltonianPairingWitness   (禁止変換 33)
ChargeNonconservingResponse + DissipativeChargeConservationCertificate
                                          → HamiltonianPairingWitness
```

- `HamiltonianTopology`/`HamiltonianPairingWitness` は private 門フィールドを持ち、
  証明書つき昇格関数だけが構成できる (From/Into 不在 — v27.2/v33.2 の規律)。
- 証明書は宣言でなく**機械検査**: `certify_dissipator_locality` は Γ = Λᵀ + M の
  cross-node block ノルムを構成時に検査 ([C2] で collective loss は拒否)。
  `certify_charge_conserving_dissipation` は本 lane では散逸ゼロの検査 (線形 jump
  は必ず電荷を ±1 変えるため — quadratic dephasing は lane 外)。
- pairing (Δ ≠ 0) の宣言は構成時 OutOfDomain ([C4] — 強制回答しない)。
- 生成子が probe に依存する場合 (back-action) は二 ε signed-linearity ゲートが
  検出して OutOfDomain ([C5]: 正直 lane は比 1 厳密・κεP back-action は 1.7e-2 逸脱)。

## 4. 観測商 [X]_C と gauge 族

登録契約 C = (probe 射影族, 観測量族, 時刻族) の応答 Φ_C が識別するのは X の
同値類のみ。少なくとも次が全時刻応答を保存する ([C6]/[C7] 機械実証):

- **local phase gauge**: X ↦ D X D† (D = diag(e^{iφ_a})) — 曲率・応答表とも保存
- **global frequency**: X ↦ X + iωI
- **複素共役**: X ↦ conj(X) — 応答表が厳密同一で X ≠ conj(X) (‖差‖ = 0.65) —
  契約は X を商までしか識別しない (EquivalenceClassOnly)

**未解決 (正直な残高)**: これらが generic case の観測同値の全てか (商の完全分類)
は開いた数学問題 — v35.4 の Packet A で FollowUp に反例探索を依頼する。

## 5. no-go の観測契約 scope (誇張しないための明示)

- GQF-3 の契約は **CurvatureOnlyOpenResponse** (曲率統計のみ)。反例対の全時刻
  応答は異なる (coherent は振動・loss は減衰) — 全時刻契約での Hamiltonian
  support 識別可能性は開いた問題であり、本 no-go はそれを主張しない。
  既知の gauge 族 (§4) は Hamiltonian support を保存するため、全時刻契約に
  対する support 反例は現時点で存在しない。
- GQF-5 の契約は **単一状態・単一時点の電荷応答値**。電荷軌道全体は loss
  (単調減衰) と pairing (振動) を区別し得る。
- どちらも「広い観測契約族を排除する新しい厳密 no-go」の**候補**である —
  正式に完成条件 3 に数えるのは、観測契約の固定 (本文書 + OpenQuotient.lean) に
  加えて FollowUp の反例探索を通してから (PROMPT/16 §10・方針 §2.5)。

## 6. 検査一覧 (v352_open_quotient — 18 PASS)

| 検査 | 内容 |
|---|---|
| [G1a/b] | covariance 閉包 vs dense Lindblad (N=2 loss+gain / N=3 複素 loss) ≤ 3e-15 |
| [G2a/b] | GQF-1 affine 消去 4e-17・生軌道の Y 依存 (負制御) |
| [G3a/b] | GQF-2 jet 恒等式 (厳密 0)・Richardson 1e-7 |
| [G4a/b] | Lean 整数橋 (曲率 1,1・電荷 −8,−8)・宣言スケールの dense 較正 (dN/dt = −2) |
| [C1–C3] | coherent 昇格可 / collective loss 昇格拒否 (禁止変換 32) / 対角 loss 還元 (GQF-4) |
| [C4–C5] | pairing 構成時 OutOfDomain / 二 ε ゲート (back-action 検出) |
| [C6–C7] | local phase + 周波数 gauge / 複素共役 → EquivalenceClassOnly |
| [C8–C9] | 有限 shot 4 裁定 (GQF-6) / 相関 shot 拒否 (v35.1 ゲート継承) |
| [C10] | GQF-5 証明書の門 (loss は witness 不能・散逸ゼロ証明書のみ解錠) |

## 7. 限界と非主張

- lane は**有限モード・number-conserving・quasi-free・Markov**。一般 GKLS
  (Kossakowski 同値類・jump gauge)・quadratic dephasing・non-Markov は lane 外
  (構成時拒否か OutOfDomain — 捏造しない)。
- Lean の一般定理は ℤ[i] 値の恒等式 (→ 多項式恒等式)。ℤ→ℂ の稠密性論法自体は
  未形式化。ODE の解析 (存在・収束) も未形式化 — jet 恒等式と数値較正で挟む。
- GQF-5 の pairing 側は宣言スケール規約 (R2 = 2R, QN = 2Q) — [G4b] が dense
  量子発展の dN/dt との整合 (×4) を機械確認する。
- これは synthetic lane の定理と器械であり、実データ・外部再現・bridge law・
  自然の的中はいずれも増えていない (正直な残高は不変)。
