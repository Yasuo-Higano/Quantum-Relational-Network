# QRN-Core v1 仕様 — 意味論的閉包と型付き境界 (第二十八期)

**Version**: v27.1 起草 (2026-07-28) / v27.2 で型実装 / v27.3 で器械凍結
**位置づけ**: PROMPT/10 §1–§2 の成果物。uft-v27.0.md §D.3 の「QRN-Core v1 の
定義要件」を執行する — ただし §D.3 の成果物を**無条件に QRN-Core v1 と呼ばない**。
それが定義しているのが core / adapter / instrument のどれかを判定し (判定結果:
adapter + instrument — §2)、理論本体・測定器・橋渡し仮説を型レベルで分離する。
用語法と型境界は [qrn-terminology.md](qrn-terminology.md)、機械可読の分類は
[core.schema.yml](../core.schema.yml) (監査: `v271_core_audit`)。

---

## 0. 凍結された到達点 (前提 — 変更禁止)

1. v27.0 の dynamic metric fork は**分岐 (b) external metric で確定済み**。
2. 自由場 matter loop には有限留数を持つ graviton pole は**ない**。
3. bare c₂ は a⁻² 走行する regulator 量であり、「誘導 Newton 定数」とは呼ばない。
4. BOND-A・4D Ward・Belinfante λ=−1/8・staggered/Wilson universality は
   **認証済み測定器**であり、QRN または創発重力の証拠ではない。
5. **自然の観測量の的中 0・独立外部再現 0** を維持する (機械監査 R1–R3)。
6. composite graviton・1/Π の常用・dynamic metric の再開は**禁止**
   (封印解除条件は spec §13.3 と本仕様 §6)。

## 1. 判定 — §D.3 は何を定義していたか

v27.0 §D.3 の要件 4 項を層分類 (qrn-terminology.md §1) にかけた結果:

| §D.3 の要件 | 実体 | 層 |
|---|---|---|
| 1. 「公理系」(格子 = regulator・BOND-A・中点変調・Belinfante λ=−1/8) | 格子離散化と source 結合の契約 | **adapter** (QRN の公理ではない) |
| 2. 器械スイート (null 結合・殻積分・Matsubara Ward・接触項・導出モデル外挿) | 較正済み測定器と較正記録 | **instrument** |
| 3. 主張の等級 (C0–C5 + 自然の的中 0 明示) | 台帳 | **meta** |
| 4. falsifier (4 比 universality・和則・Ward の常設回帰) | 常設監査 | **instrument (meta)** |

**結論**: §D.3 が定義するのは「外部計量上の Dirac 物質を格子正則化で正しく測定
するための認証済み計量応答カーネル」= **QRN-Matter-on-Background Adapter v1 +
QRN-Metrology Suite v1** である。これは価値があるが、QRN の中心命題

> 量子相関網から、時空・重力・物質・因果が読み出される

の核 (QRN-Kinematics / QRN-Dynamics / QRN-Bridge) そのものではない。
よって第二十八期の成果物は次の 4 つの名で呼ぶ:

```text
QRN-Kinematics v1                     (§2 — 定義済み・実装は toy のみ)
QRN-Matter-on-Background Adapter v1   (§3 — 認証済み, v27.3 で凍結)
QRN-Metrology Suite v1                (§3 — 認証済み, v27.3 で凍結)
QRN-Bridge Hypotheses v0              (§4 — 仮説段階)
```

「QRN-Core v1」という名称は、上記 4 層 + QRN-Dynamics の**状態表示つき複合**
としてのみ使う (§5 の status block が機械検査アンカー)。

## 2. QRN-Kinematics v1 (layer: core) — 状態: defined

QRN の存在論的骨格。**定義**であり、現時点の実装はガウスフェルミオン toy
1 族 (C3) しかないことを明示する。

- **StateSpace**: 有限次元 Hilbert 空間の族と状態。現行実装 =
  `GaussianFermionState` (旧 QrnState — 相関行列 C_ij = ⟨c†_i c_j⟩)。
  一般の (非ガウス・非フェルミオン) 状態空間は**未実装**。
- **ObservableAlgebra**: 局所演算子代数とその部分代数束。現行実装 =
  相関行列の部分行列 (ガウス系では十分)。一般代数は未実装。
- **RelationalDecomposition**: テンソル分解は入力ではなく読み出しである
  (v11.4 の動的テンソル分解 — MI 貪欲マッチング)。関係要素の型は
  `RelationalNodeId` — 格子点 `RegulatorSiteId` と別型 (禁止変換 §2.1)。
- **AllowedEquivalences**: 基底変換・分解の組み替えのうち物理を変えない同値類。
  現状は toy ごとの個別実装 — 一般規則は未定義 (Unknown)。

## 3. QRN-Dynamics (layer: dynamics) — 状態: model_family_only

**QRN 固有の動力学原理は存在しない。** あるのは模型族だけである:

- ガウス toy 族: RingChain / TfdPair / GrowingChain / PacketRing
  (`GaussianToyModel` = 旧 QrnModel)。`evolve(s, t)` の `t` は
  **EvolutionParameter (外部発展パラメータ)** であり A1 の創発時間ではない。
- 拘束模型族: SM 直積 core (v22.3)・SU core 階段 (v20–21) — 「拘束を基底で解く」
  語彙の模型群。
- InitialConditionRule: **未定義** (時間の矢の初期条件問題は v5.1 の機構どまり)。
- ConstraintAlgebra: 模型ごとの個別実装 — 一般代数は未定義。

**空の dynamics を「Core 完成」と判定してはならない** (完了条件 §7)。

## 4. Matter-on-Background Adapter v1 / Metrology Suite v1 (layer: adapter / instrument)

経路 B (v26.2–v27.0) で認証済み。各器械の登録簿 (input/output type,
normalization, continuum contract, calibration source, negative controls,
known failure modes, allowed claims, forbidden interpretations, certificate)
は **v27.3 の instruments.yml に凍結**する。対象:

- adapter: BOND-A 結合則・中点変調 (置換則)・Belinfante 改良 (λ = −1/8)・
  接触項 (2 実装照合)・counterterm Λ√g・staggered/Wilson 二離散化・
  continuum trajectory (a→0, 固定 m_phys)。
- instrument: null 結合 ladder・殻積分 (root-solve)・Matsubara Ward・
  導出モデル外挿 ({1, a²ln(1/a), a²} 型)・spectral measure・f-sum rule・
  解析 oracle (三重導出)・厳密ブロック分解・dd/iv 経路・監査台帳群。

常設回帰 (suite の監査層): 4 比 universality (PRED-016)・スカラー和則
(PRED-017)・4D Ward 恒等式群・v252_manifest (凍結不変性)・v271_core_audit
(本仕様の意味論)。

**許される主張**: 「測定器が正しい」(operator/regulator universality)。
**禁止される解釈**: QRN の証拠・創発重力の証拠・誘導 Newton 定数
(spec §12.8/§13.3 の凍結解釈)。

## 5. QRN-Bridge Hypotheses v0 (layer: bridge) — 状態: conjectural

状態・相関から計量・因果・時計・物質を読み出す規則。**全て仮説段階** —
toy 実演 (C3) はあるが、bridge law (regulator に依らない一意な読み出し規則)
は一つも確立していない。

- GeometryReadout: B1 mutual-information distance (v0.7 円環 100% — toy) /
  B2 modular-flow geometry (v19–25 の BW 系 — 器械は較正済み、bridge は未) /
  B3 relative-entropy / quantum-Fisher response (v2.1, v15.6 — toy) /
  B4 commutator・Lieb–Robinson causal geometry (v1.1, v6.7 前線 — toy)。
- CausalReadout: 光円錐の創発 (v1.1) — toy。Lorentzian 因果構造の再構成は未。
- ClockReadout: **未構成** (ModularParameter / OperationalClockReading /
  ProperTime の 3 型を v27.2 で分離定義 — 中身は空)。
- MatterReadout: k_F・密度 (v6.7) — toy。
- **GravityBridge: unsupported** — δS=δ⟨K⟩ の toy 検証 (v0.7/v15.6/v19 系) は
  エンタングルメント第一法則 (QI の一般恒等式) の数値再現であり、Einstein
  方程式の検証ではない (uft-v0.7.md / uft-v1.0.md A3 の監査注記、
  QRN-C0-001 の限定条項)。EmergentMetricCandidate を支持する結果は 0 件。

### bridge の成功条件 (次期以降の物理課題 — 実行前に bridge_candidates.yml へ凍結)

1. 座標・隣接関係・外部計量を入力しない。
2. 二つ以上の微視的に異なる模型で同じ巨視的幾何へ収束する。
3. 基底変換や tensor factorization の軽微な変更に頑健。
4. 空間距離だけでなく Lorentzian 因果構造を再構成する。
5. 未使用の応答チャネルを事前予言する。
6. 読み出された計量が、source として入力した計量のコピーになっていない。
7. 失敗した bridge は調整せず棄却する。

PRED-019 は、QRN 固有の数値または明確な二者択一が解析的に導出できるまで
**登録しない**。

## 6. composite graviton 路線の封印条件 (拡張)

v27.0 の封印解除条件 (相互作用による普遍 q² 項 + Weinberg–Witten 破れ仮定の
明示) は**必要条件であって十分条件ではない**。再開には少なくとも:

```text
CG0  regulator-independent な q² 項          CG5  soft-graviton limit の整合
CG1  正の有限留数を持つ孤立 spin-2 pole      CG6  Weinberg–Witten 回避機構
CG2  ghost・余分な spin-0 mode なし           CG7  kinematic nonlocality または同等の構造
CG3  3-point での diffeomorphism Ward         CG8  外部計量を入力せずに同じ metric sector
CG4  全物質と重力自己エネルギーへの普遍結合        を得る
```

(CG7 の根拠: 局所的な微視的運動学からの創発重力は kinematic nonlocality を
要するという一般制約 — arXiv:1409.2509。) この路線は第二十八期では再開しない。

## 7. QRN-Core v1 の状態表示 (機械検査アンカー — v271_core_audit R10)

```yaml
qrn_core:
  kinematics: defined
  dynamics: model_family_only
  geometry_bridge: conjectural
  gravity_bridge: unsupported
  empirical_prediction: none
```

**Core v1 完了条件** (PROMPT/10 §7 — 全て満たすまで「Core v1 完了」と言わない):

- ontology / regulator / background / instrument / bridge が分離済み (v27.1)
- 状態空間・観測量・発展則の状態が明示済み (v27.1 §2–§3)
- 外部時間と創発時間が別型 (v27.2)
- 外部計量と創発計量が別型 (v27.2)
- 旧 Gaussian core の適用範囲が限定済み (v27.1 用語法 §3, 改名は v27.2)
- v0.7 / A3 の過大解釈が修正済み (v27.1)
- 自然観測 0・独立外部再現 0 が保持される (常設監査 R1–R3)
- 全監査 PASS・make suite PASS
- 負の結果と未定義部分が文書から消されていない

**完了しても得られるのは「QRN 固有の理論核を構築できるかを初めて失敗可能な形で
問える状態」であって、統一理論ではない。** 現時点で完成しているのは強い計算物理・
理論監査プラットフォームである — この一文を弱める改訂を禁止する。
