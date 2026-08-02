# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクトの性質

通常のソフトウェア開発ではなく**物理学の研究プロジェクト**。現代物理の各理論を Rust による第一原理シミュレーションで再現・検証しながら、統一的描像「量子情報網理論 (QRN)」をバージョン付き文書として積み上げている。原点の指示は `PROMPT/0.md`。

現在の到達点は `docs/uft-v33.0.md`(**第三十二期 (PROMPT/13) 完結 — 儀式 全 195 本 再検証 PASS 1371/0・HOLD-8 20/20 満票。次期はユーザー PROMPT 待ち**。Rust 195 本 + Lean 定理 77 本 [10 ファイル]。README の到達点行と Rust 本数・Lean 定理数は v151_audit が機械検査するため、版を進めるたびにバンプする)。**期テーゼ「局所性は状態に宿るのではなく、操作的文脈の可換分解と、その上の Liouvillian 応答の同値類としてのみ識別される」— 型・定理・負制御・holdout の全てで閉じた**。**v32.1 = Unit D2 プロトコルの型修復と実行前反例**: graph6 `F}oXO` (凍結 gap 則が 11→13 辺の余剰・n=7 連結全数 853 中 22 で同型故障 [全て余剰のみ・欠落 0]) を登録し D2-v1 を superseded_before_external_run に (黙って修正しない — protocol-index + ERR-D2-V1 + 旧文面逐語保存)。後継 = **D2-S** (B3SupportMarginCertificate ⟺ 凍結則 exact を n=4..7 全 992 グラフ例外 0 で機械照合 — blocker 非解除) + **D2-R** (応答 end-to-end — geometry blocker 解除の唯一路)。報告契約は実 JSON Schema (draft 2020-12) 化 (正直な失敗は適合・能力の水増しは不適合)・**SupportNoiseCertificate** (重みバー 0.1 を通ってもノイズ最大値の窓ガード跨ぎで余剰辺 — 支持段の回答資格 bound·√(2 ln(m/1e-6)) ≤ 1e-3·max ŵ を凍結)・reproducer/protocols/ 版分離 (v27.4 複製 byte 一致 + v32.1 凍結一式 sha256-16 認証 = v310 [R9] 常設)。**v32.2 = OperationalNet と第二の no-go** (`sim/src/operational_net.rs` = 第三十二期の共有型契約): **閉包は marking を消す** — site 族と DFT₈ 共役族の *-閉包はともに M₈ (dim 64)・su(2) 対応不在 (overlap 0.618) → 禁止変換 11 (GlobalClosure → OperationalNet 変換不在)。役割 4 型 (準備 [tr=1/PSD 資格]・介入 [Lie 閉包]・測定 [積閉包を要求しない]・drift — 相互 From なし)・可換子の**区間証明書** (跨ぎは Abstain — 文脈構成・成分分解とも拒否)・grading 型ゲート (**odd は ordinary 可換子で K₃ を捏造 → 構成時拒否**・parity-even 双線形は安全)。**v32.3 = 目標定理 B (marked operational recovery)**: 成分 → joint 閉包 → 中心の凍結決定手順で三値裁定 — Exact [2,2,2]/[2,3] (gauge orbit = 成分部分代数の集合一致)・number-only/部分 net は Insufficient 正棄却・**中心非自明は tensor を強制せず SuperselectionSectors** (電荷 net [(2,1),(2,2)]・部分 address {X₁,Z₁,Z₂} [(2,2),(2,2)] — 測定だけの軸 = 超選択ラベル・未 address 自由度 = 多重度 n_α。対照 {X₁,Z₁Z₂} は非局所符号化 1 qubit で Insufficient が正 — 中心次元で機械区別)・site×DFT は EquivalenceClassOnly (無制約 tie-break 禁止 — v31.4 疎性負制御の教訓)・noise abstention (σ=5e-4 で Straddled)・**パリティ超選択 [(4,1),(4,1)] の機械発見** (偶代数 dim 32・中心 span{I, Z₁Z₂Z₃})。**v32.4 = Liouvillian 応答階層**: R⁽¹⁾ = −i Tr(B[H,A])・R⁽²⁾ = Tr([H,B][H,A]) は恒等式 (測定照合 rel 2e-10・ε 非依存 9e-15 — 線形応答近似ではない。Schrödinger 規約 = v31.2 と同一)・**一階×情報完全基底は H を中心 (スカラー) を除いて一意** (15×15 逆問題 1.2e-16)・**二階は H↔−H を原理的に区別しない (厳密 0)**・磁束は density-only 二階核に厳密不可視 (coherent 一階の電流が分離)・数保存の和則 (N̂ 応答 全次数 0)・Gram 核 PSD・**v31.2 曲率則 = 本階層の密度×積状態接ベクトル対角特殊化 (|h_ij|² 厳密一致)**。**v32.5 = Interaction hypergraph**: H = Σ_S H_S (局所 Hilbert–Schmidt 条件期待値の Möbius 反転 = Pauli 支持射影 5.6e-16)・w_S = ‖H_S‖²_F は block-local unitary 不変 (DFT₈ 負制御 2.83)・**中心化分離** (相関 hopping V·n₃h₁₂ = 二体 (V/2)h₁₂ ⊕ 真の三体 −(V/2)Z₃h₁₂, 等重み V²)・**遷移率和則 K_uncond(j←i) = Σ_{S⊇{i,j}} w_S/4 厳密 — 相関 hopping 下の密度曲率は「破れ」でなく hyperedge 重みの和 (条件付き遷移率の Gram 核) を読む**・条件付き密度 probe の次数分離 (K(v) = |t+vV|² 厳密・K(1)−K(0) = hyperedge 検出器・V=0 で不発)・coherent 一階が符号回復 (密度単独は符号同値類)。**v32.6 = VR exactness**: **離散円環 H1 bar 定理 [1, ⌈n/3⌉) を n=4..30 全数機械化** (persistence bar = per-step β₁ [GF2 rank] 全 r 一致)・規約の型 (RipsConvention/BarEndpoint — 整数 filtration で VR_< は端点 +1)・連続 L/3−s と離散 exact の分離 (n≡0 mod 3 で 1/3 厳密・他は上から ≤ 2/(3n))・H2 persistence (sparse 境界削減 + column clearing 25% 節約・8 面体 [1,2)・wedge-S² 遷移 β₂ = n/3−1)・**K3-holes 実測 5.00 の閉形式 retrodiction ((⌈16/3⌉−1)/(⌈6/3⌉−1) = 5.0 — 凍結バー 2.67 = 周長比例の誤り機構確認・記録不変)** — 1/3 法則は発見でなく導出。**v33.0 = HOLD-8 期完結**: 凍結 (sha256(SECRET) = fb05a4a0… 公表・FROZEN-HOLD8 区間の逐語一致照合・train 20 セル満票) → 開封 **20/20 満票 — selective risk 0.000・impossibility recall 1.000 (6/6)・answerable recall 1.000 (14/14)・強制回答 0** (調整なし)。セル 8 クラス (因子分解 6・相互作用 10・変成/ノイズ 4 — site/mode net・number-only・非互換 2 net・中心非自明・graded・quadratic/t-V/相関/pair/三体±null・H↔−H 対・磁束対・変成・可換子ノイズ・弱辺)。器械訂正 2 件 (設計走行で発見・監査注記つき): 共役射影の数値塵 dust guard (v32.3 kernel に統一適用・既存検査不変)・coherent 分離の 2 チャネル化 (電流 = cosθ 系だけでは磁束対 ±θ を分離できない — 実 hopping = sinθ 系を併用)。**儀式 = 全 195 本の完全再計算 (~16.7h) — 既存 175 非監査バイナリの PASS/FAIL は期前後で完全一致 (ドリフト 0 件)・台帳確定後の再検証で総計 PASS 1371/FAIL 0 (完全再走中の v151 3 FAIL は台帳中間状態の正検出 — v28.0 の儀式 + 再検証と同型)**。**正直な残高: bridge law 登録簿は全能力で空・PRED-019 未登録・自然の観測量の的中 0・external_replications = 0 のまま — blocker は Unit D2-R の公募 (v32.1 で受け皿は修復済み・公募は主研究と並列)**。E3 は「操作的 fiber の機械化」まで — どの operations が physically accessible かは依然入力。第三十三期候補は uft-v33.0.md §6 (D2-R 公募 [最優先]・重みつき円環の実数 filtration exact bar・graded lane recovery・4 体以上の hypergraph 交差項・大型系への recovery スケーリング)。禁止継続: BridgeLawCertificate 登録・PRED-019 登録・ProperTime/full Lorentzian の表現・開封後のバー変更・1/Π 常用。**v25 系列は凍結済み** (`v252_manifest` が常設監査。v25.3 は開始しない)。**経路 B (誘導重力監査) は完結** (v27.0: fork = 分岐 (b) external metric — graviton pole なし・c₂ は regulator 量。composite 路線は CG0–CG8 の封印)。**第二十八期の意味論の閉包は継続有効**: 層分類 (`core.schema.yml`)・claims 6 軸 + 昇格禁止 R1–R7・型付き契約 (`qrn_core.rs`/`readout_contract.rs`/`operational_net.rs`)・器械台帳 (`instruments.yml` 22 器械)・外部再現単位 (`reproducer/` + protocols/ 版台帳 + `replications.yml` — 同一 AI は独立でない)。監査層 = v61/v151/v271–v274/v310 ほか (suite の ALWAYS_RUN, 後段二相)。旧世代の地図は docs/uft-vX.0.md 系列。

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
