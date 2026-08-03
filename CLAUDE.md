# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクトの性質

通常のソフトウェア開発ではなく**物理学の研究プロジェクト**。現代物理の各理論を Rust による第一原理シミュレーションで再現・検証しながら、統一的描像「量子情報網理論 (QRN)」をバージョン付き文書として積み上げている。原点の指示は `PROMPT/0.md`。

現在の到達点は `docs/uft-v34.2.md`(**第三十四期 (PROMPT/15)「申告証明書から観測証明書へ」進行中 — v34.1 = FollowUp 受理 (REP-001 cross_model_clean_room / Partially Replicated / external_replications = 0 維持) + Yukawa erratum (QRN-YUK-034 refuted_as_stated・QRN-YUK-015 走査族限定に分割・論文 2 本に Erratum 節)・v34.2 = standalone operational core spec (paper/operational-core-spec.md OCS-1.0 sha256 pin bb3b4d4b…・closure manifest 20 エントリ [paper_closed 18/replay_only 2]・probe 型分離の機械実証 [quench null / pairing 1−Δ² / Busch 2η²≤1] — v342_spec_manifest が常設監査)・v34.3 = 有限データ昇格不能定理 (第四の no-go — proofs/FiniteDataNoGo.lean 10 定理: Le Cam 下限 (1−TV)/2 + Robust Promotion 誤昇格 ≤ α) + `sim/src/finite_data.rs` (RobustVerdict 5 値・Clopper–Pearson・禁止変換 22–29 の型) + v343 厳密反例 17 検査 (selective risk ≠ coverage [P(wrong|ans) 0.318 > α]・0 誤りの上限 0/9→28.3%・1% には 299 セル = HOLD-10C 設計根拠)・v34.4 = sector-aware factorization enumerator (`sim/src/factorization_enumerator.rs` — Wedderburn 証明書 [n²/m²/nm=d/積 span/A″=A]・multiplicity 分離・候補列挙 + orbit witness・出力 6 型 [Unique/CandidateSet/Sectorwise/IncompletePrimitiveSet/NontrivialCenterObstruction/ScopeExceeded] — FollowUp FAC-001 への応答)・v34.5 = robust atlas (v345 — σ_min 同時下界 [Weyl]・worst-case cross-talk・glue overlap 区間・spectral-gap 証明書つき J 構成・interval cost 3 値 + set-valued profile・dense/Pauli 同一区間意味論 — synthetic lane であり実測ノイズではない)・v34.6 = real-data lane (v346 — DataProvenance 型 [synthetic ↛ experimental]・reproducer/real_data/ 受け皿 [schema/事前登録 sha256/drift gate/台帳 recorded_runs = 0]・D2R_PACKET 配布パケット・v35.0 完成条件凍結 [外部報告/実データ予測/新 no-go のいずれか])**。Rust 209 本 + Lean 定理 87 本 [11 ファイル]。README の到達点行と Rust 本数・Lean 定理数は v151_audit が機械検査するため、版を進めるたびにバンプする。**第三十四期の版計画 (PROMPT/15 §8)**: v34.2 = standalone operational core paper/spec (結果値なしで独立実装可能に — FollowUp の specification-limited への応答・SignedInitialCovarianceProbe ↛ HamiltonianQuench / NumberConservingResponse ↛ BCSResponse の型分離)・v34.3 = 有限データ no-go (Le Cam 二点下限の Lean 化) + Robust Promotion Theorem (同時信頼集合 C_α(D) → RobustExact/EquivalenceClassOnly/Straddled/InsufficientObservation/OutOfDomain) + 禁止変換 22–29 (PointEstimate ↛ Certificate / MarginalIntervals ↛ JointConfidenceRegion ほか)・v34.4 = sector-aware factorization enumerator・v34.5 = robust atlas/glue/J/profile (信頼集合意味論)・v34.6 = real-data lane + D2-R 実募集・v35.0-A/B = HOLD-10S (30 semantic) + HOLD-10C (600+ coverage, population-risk 上限)。非交渉: HOLD-9 記録不変・external 0 維持・点推定/marginal/GoF/単一時点較正の証明書への暗黙変換禁止・synthetic を experimental と呼ばない・一般 GKLS/gravity/PRED-019 は中心にしない)。**期テーゼ「可アクセス性は作用素単体の属性ではなく、系・制御器・測定器・資源・誤差証明書の関係である。局所性は、証明付き laboratory interface が生成する role-typed context atlas 上で整合する因子分解の、資源スケールにわたる安定な同値類としてのみ識別される」— 型 (禁止変換 12–21)・定理 (第三の no-go)・負制御・holdout の全てで閉じた**。**v33.1 = 境界監査と型スコープ修復** (`operational_net.rs` v33.1 節): contexts 盲目性の機械実証 (v32.3 参照手順は contexts ∅/atlas で読み完全一致 — 入力でなかった)・entangler 負制御 (X₁X₂ 1 本で [2,2,2]→[2,4]・閉包同一 M₈ = **primitive 選別の循環**)・**MarkedRecoveryInput** (gens 別渡し廃止 — 生成子行列は net の primitive のみ・RoleMixedRecovery/NoDeclaredContexts/ContextCoverageIncomplete は Abstain でなく構成時型エラー)・禁止変換 12 (CertifiedCommutator ↛ JointContextWitness — singleton 文脈は全対 definite 証明書でも Abstain(OperationalCompatibilityUnwitnessed))・13 (可換子 lane ↛ GKLS — Leibniz 導分証明書: lane は厳密 0 + Ĥ 復元 残差 0・GKLS は γ 厳密比例で破れ・R⁽¹⁾ 公式は GKLS 測定と乖離 1.2)。**v33.2 = Certified Laboratory Interface** (`laboratory_interface.rs`): 宣言 ≠ 資格 (禁止変換 14) — 門は較正 (IndependentAddressabilityCertificate: 証明つき rank・σ_min・cross-talk 区間)・合成 (SynthesisCertificate: bracket/線形列の機械実行検証)・トモグラフィの 3 出自証明書のみで **sha256 結束** (流用は構成時拒否)・文字列 provenance 廃止。tied control no-go (禁止変換 15: u(t)(X₁+X₂) の数学的分解は rank 1<2 で拒否・正直 net は Abstain)・**可アクセス性は interface との関係** (同じ X₁ が {X₁+X₂} では路なし 0.707・+Z₂ で Synthesized depth 3)・**controller-free decomposition no-go (E3-A, 第三の no-go)**: 同一 (H, drift, ρ) で 4 interface が非同値 ([2,2,2]α/[2,2,2]β 0.5625/[2,4]/Abstain) — 状態単独 ✗ (v31.4)・閉包 ✗ (v32.2)・controller-free ✗ の三段完成。role-typed 文脈 4 型 (禁止変換 16 — **joint measurability は可換性より広い**: 非可換 unsharp η=0.6 が joint POVM で資格・0.8 は正値性破れ = Busch iff の器械化)・ResourceBudget は 5 成分半順序 (Ord なし)。**v33.3 = Resource-Filtered OperationalNet** (`resource_profile.rs`): 中心対象は **profile: budget ↦ 読み** (資源不足 → [2,2,2] → [2,4] → [8])・**poset は barcode でない** (比較不能対が同 dims 別 orbit)・昇格規則 (禁止変換 17): stable ⟺ chain ≥ 2 (単点 transient 昇格 0・grid 頂も grid 相対 transient)・スカラー潰し (禁止変換 18) は accessibility を変え読みを [2,4]→[8] に反転・**頂は経路を消す** (erasure no-go の resource 版)。**v33.4 = Contextual factorization** (`contextual_factorization.rs`): chart (持ち場) 局所復元 (証明書は大域 net から継承・chart 内証人ゲート・factor 資格) + overlap glue — **glue 定理 (= 直接大域復元と読み・orbit 一致)**・変成不変 (u M₂ u† = M₂)・cocycle 不整合 (CZ₂₃ 捻り overlap 1/3) → Abstain(GlueInconsistent) = 禁止変換 19 (chart 局所 Exact ↛ 大域)・複数 glue (site vs DFT atlas) → EquivalenceClassOnly{2}・witness 境界は両 lane 一致。**v33.5 = Graded recovery の境界** (`graded_recovery.rs`): **Majorana locality ≠ Dirac locality** — odd CAR は O(2N) 不変で graded graph は全対厳密 0 (pairing 情報なし)・witness なしは MajoranaFrameOnly (禁止変換 20)・charge witness Q から複素構造 J = ⟨γ_b, i[Q,γ_a]⟩/‖γ‖² (実・反対称・J²=−I 4.4e-16) → 3 モード回復 (mode-CAR 6.3e-16・Σâ†â=Q 1.2e-15・U(N) を除く)・部分 charge は ComplexStructureUnresolved・quartic 汚染は WitnessNotLinearOnFrame・graded net の marked recovery は Abstain(ComponentNotFactor) (捏造しない)。**v33.6 = 構造化スケーリング** (`structured_backend.rs`): lane 分離 — Pauli GF(2) symplectic (可換性 = ω・閉包 2^{dim V}・**中心 = radical**: n=3 全 7 セルで dense と裁定完全一致 [超選択 sector まで]・48 qubit 証明書を行列なしで [2×48]/[2×46,4]/[(2⁴⁷,1)×2])・Majorana quadratic (対応原理: 支持分割 = dense 成分・dense 閉包 = 2^{2m−1} [偶 Clifford — full M_d を与えない]・48 Majorana so(16)³ → so(32))・**ScopeExceeded は正答** (禁止変換 21: dense dim > 4096 は試行しない・非 Pauli 和は構成不能)。**Track X = D2-R campaign layer** (`reproducer/campaigns/d2r-v1/` — v310 [R10] 常設監査): 数値 kernel なしの公募受け皿 — 事前登録 schema (凍結 sha 結束)・validator 較正 fixture (正1負2)・AMBIGUITIES (凍結文は追記明確化のみ)・OUTREACH/FUNNEL 台帳 (**数は実記録のみ・捏造禁止・protocol_viewed は not_instrumented**)・凍結された約束 (**有効報告は一件で足りる — 二件要求へ変更しない・外部再現の成否を内部 holdout に入れない**)。**v34.0 = HOLD-9 期完結**: 凍結 (sha256(SECRET) = ef6a8cd9… 公表・FROZEN-HOLD9 逐語一致・**lib pin 6 モジュール sha256-16**・train 20 セル満票・設計 12 シード頑健) → 開封 **20/20 満票 — selective risk 0.000・impossibility recall 1.000 (11/11)・answerable recall 1.000 (9/9)・強制回答 0 + 新 5 計量 (出自被覆 31/31・証人被覆 7/7・raw 昇格 0・scope 違反 0・transient 昇格 0)** (調整なし)。セル 5 群 (入力完全性 4・accessibility 4・context/resource 4・graded 4・scale/変成 4)。器械訂正 2 件 (監査注記つき・凍結区間/採点/裁定に触れない scaffolding): [H2] 区画定数 8/12→9/11 (初版 [H2] 自身が正しく FAIL して検出)・v336 壁時計バー (儀式の JOBS=12 並列負荷で超過 — 「並列化で結果が変わらない」規約 [PROMPT/4] は check 条件にも適用される。計時を合否から除外)。**儀式 = 全 203 本の完全再計算 (~18.2h) — 既存 183 非監査 + 監査 12 の PASS/FAIL は期前後で完全一致 (ドリフト 0 件)・台帳確定後の再検証で総計 PASS 1429/FAIL 0 (儀式中の v336 FAIL 2 は壁時計器械バグの正検出 — v33.0/v28.0 の儀式 + 再検証と同型)**。**正直な残高: bridge law 登録簿は全能力で空・PRED-019 未登録・自然の観測量の的中 0・external_replications = 0 のまま — blocker は D2-R の外部報告のみ (campaign layer で受け皿拡充済み・funnel 全段 0 の正直な記録)**。E3 は A (controller-free no-go) と B (certified interface 下の復元) で閉じたが、**証明書の内容 (較正・合成・資源上限) は依然として実験者の申告**であり実測ノイズ下の区間較正は未走査。第三十四期候補は uft-v34.0.md §6 (D2-R 公募 [最優先]・実測ノイズ下の証明書 [addressability/J/glue の区間資格と Straddled]・structured lane の統合 [Pauli/quadratic 上の profile・glue・graded witness]・一般 GKLS 応答 [jump gauge・Kossakowski 同値類]・BCS 型 witness・profile の関手性と安定性定理)。禁止継続: BridgeLawCertificate 登録・PRED-019 登録・ProperTime/full Lorentzian の表現・開封後のバー変更・DeclaredOperation → AccessibleOperation の暗黙変換・budget のスカラー化・transient の昇格・1/Π 常用。**v25 系列は凍結済み** (`v252_manifest` が常設監査。v25.3 は開始しない)。**経路 B (誘導重力監査) は完結** (v27.0: fork = 分岐 (b) external metric — graviton pole なし・c₂ は regulator 量。composite 路線は CG0–CG8 の封印)。**第二十八期の意味論の閉包は継続有効**: 層分類 (`core.schema.yml`)・claims 6 軸 + 昇格禁止 R1–R7・型付き契約 (`qrn_core.rs`/`readout_contract.rs`/`operational_net.rs`)・器械台帳 (`instruments.yml` 22 器械)・外部再現単位 (`reproducer/` + protocols/ 版台帳 + `replications.yml` — 同一 AI は独立でない)。監査層 = v61/v151/v271–v274/v310 ほか (suite の ALWAYS_RUN, 後段二相)。旧世代の地図は docs/uft-vX.0.md 系列。

**文書・コードコメント・コミットメッセージは全て日本語。**

## コマンド

```bash
cd sim
cargo build --release                      # 全バイナリをビルド (外部依存なし、std のみ)
./target/release/v01_qm                    # 単一シミュレーションの実行
cargo run --release --bin v34_unruh        # 同上 (ビルド込み)
./target/release/v34_unruh > ../results/v34_unruh.txt   # 結果の保存 (stdout をリダイレクト)

# 全スイートはルートの Makefile から (PROMPT/5)。台帳 results/suite_manifest.tsv が
# ソース不変のバイナリを判定し、前回結果を「引用」する。
cd ..
make suite-status                          # 実行/引用の判定だけ表示
make suite                                 # 増分 (変更分 + 監査層のみ実行。統合時は OUT= で保存先指定)
make suite-full OUT=results/vXX0_full_suite.txt JOBS=8  # 完全再計算 (数期に一度の儀式)
```

`cargo test` は無い。検証は各バイナリに内蔵された厳密解・観測値との `[PASS]`/`[FAIL]` 比較で行う。`lib.rs` の `self_test()`(ヤコビ法・ベッセル等の自己検証)を主要バイナリが起動時に呼ぶ。

## 構成

- `docs/uft-vX.Y.md` — バージョン付き理論文書。vX.0 は期の統合文書(公理系の改訂と未解決問題の残高)。
- `sim/src/lib.rs` — 唯一の共有コード。自作数値ライブラリ: xorshift64* 乱数 `Rng`、複素数 `C64`、Thomas 法 `solve_tridiag_c`、循環ヤコビ固有値分解 `jacobi_eigh`、行列関数 `matfun_sym`、Bessel `bessel_i`、`ln_gamma`、ビニング統計 `mean_err`、`linfit`。
- `sim/src/bin/vXY_topic.rs` — 1 実験 = 1 バイナリ。名前を文書バージョンに対応させる (例: v3.4 ↔ `v34_unruh.rs`)。冒頭の `//!` コメントに物理的背景・方法・検証内容を書く。
- `results/vXY_topic.txt` — 各バイナリの stdout を保存したもの。**文書中の全数値の一次ソース**。文書に数値を書く前に必ずここへ保存する。
- `README.md` — バージョン一覧表(文書・バイナリ・主結果の対応)。

## 作業規約

1. **外部クレート禁止**。必要な数値計算は `lib.rs` に自作して追加する。
2. **全シミュレーションに `[PASS]`/`[FAIL]` 判定を内蔵する** — 厳密解・観測値・双対性・別法いずれかとの比較。乱数は固定シード (`Rng::new(seed)`) で再現可能にする。
3. **バージョンごとに 1 コミット**。メッセージは `vX.Y: 日本語タイトル` 形式(`git log --oneline` 参照)。
4. 新バージョンを作ったら README のバージョン一覧表に行を追加し、**`claims.yml` に主張を追記する**(等級 C0–C5・証拠・限界。`v61_ledger` が常に PASS であること)。
5. **知的誠実性**: 文書では「教科書的に確立された事実」「確立された研究成果」「本プロジェクトの仮説」を明確に区別する。失敗した検証や設計ミスも削除せず記録する(例: v0.7「1 量子は無限小でない」)。

## 数値計算の既知の落とし穴

- モジュラー核 ln((1−C)/C) は区間長 ℓ≳40 で f64 の限界を超える (ζ~e⁻³⁰)。
- エンタングルメント第一法則 δS=δ⟨K⟩ の検証は「無限小のハミルトニアン摂動」で行う。有限の 1 量子励起では不等式 S_rel≥0 になってしまい一致しない。
- 理論空間の全数探索では群論係数に注意(過去のバグ: 弱二重項の次元の二重計上)。ヌル方向の規格化は T_kk = 4T₋₋。
- 固定 kη での評価はスケール不変スペクトルと自己相似の区別がつかない(v2.4 の教訓)。
- モジュラー核 κ = ln((1−c)/c) は固有値 c の 1e-14 クランプで κ ≤ 32.2 に飽和する (f64 の分解能床)。真の κ ~ 100 の深部モードに重みを持つ摂動 (拡がった固有モードの回転など) は δ⟨K⟩ が第一次で系統的に過小になり、第一法則の比に α 非依存の定数バイアスが乗る (v19.1: 半空間 N=12 で +5.4%)。壁局在の波束対など、深部に指数的にしか触れない摂動を使うこと。
- 厳密縮退帯の局在モードを中心の昇順で並べると、中心が wrap 境界 (格子座標 0≡N) に厳密に乗るモードで ±ε の綱渡りになる。単一トーラスの特異値は不変だが積模型の世代対が変わり lnZ が 1〜4 nats 動く(v9.2 の発見)。中心を 0.5 サイト格子にスナップしてからソートすること。

## 追加規則
- シミュレーションは原則としてrustで行うが、定理証明器が得意な分野についてはcoqやlean4などを自由に導入して構わない。
- 図表のような人間向け資料の作成にはpythonなどの便利な環境やフリーソフトウェア・ライブラリーを利用して構わない。
- **重い計算は必ず Rust で行う** (PROMPT/4, 2026-07-15)。python は集計・可視化のみで、重い処理をさせない (v22.1 の numpy DMRG が 20 時間超 → Rust 書き直しで ~20 倍の教訓)。Rust のマルチスレッド化は CPU 数まで可 — ただし**並列化で結果が変わらないこと** (決定性が守れないなら無理にマルチスレッド化しない。独立ジョブのスレッド分割は常に安全)。
- 長い計算の結果を待たずに進められる独立研究は、単スレッドのプログラムを複数作ってマルチプロセスで同時実行する (PROMPT/4)。

## DAG/依存グラフ
Prolog(swi-prolog)を使って処理できるようにする。
```
JSON
 ↓
Pythonで読み込み・正規化
 ↓
facts.plに変換
 ↓
Prologで推論
 ↓
DOT/Mermaid/Markdown/JSONレポート出力
```

# 計算の高速化
pythonはnumpyを使ったとしても計算が遅い。pythonは基本的に計算結果の可視化などに使い、重計算は行わないようにする。
重い計算はrustで行う。CPUのコア数を取得して、マルチスレッドで計算しても良い。ただし、計算結果が並列動作によって変わらないようにすること。
