# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## プロジェクトの性質

通常のソフトウェア開発ではなく**物理学の研究プロジェクト**。現代物理の各理論を Rust による第一原理シミュレーションで再現・検証しながら、統一的描像「量子情報網理論 (QRN)」をバージョン付き文書として積み上げている。原点の指示は `PROMPT/0.md`。

現在の到達点は `docs/uft-v29.6.md`(**第三十期 (PROMPT/11) 進行中 — v29.6 dimension-agnostic topology まで完了**。Rust 175 本 + Lean 定理 62 本 [8 ファイル]。README の到達点行と Rust 本数・Lean 定理数は v151_audit が機械検査するため、版を進めるたびにバンプする)。**v29.1 の訂正確定**: HOLD-3 採点器の linfit 取り違え (ξ = −1/切片) を再現・訂正 — **B2・B3-COV・B4-DENSITY-FRONT 生存 / B1 のみ棄却** (恒久対策: `LinearFit` 型 + 変成テスト常設・Z2 rdm BTreeMap 決定化)。**v29.2–v29.4b = HOLD-5 の完全執行**: seed コミットメント (v29.2) → S×C 採点器凍結 (v29.3, train 288 セル満票) → 定量バー凍結 + SECRET 開示 (v29.4a, val 1 回) → **holdout 初開封・本採点 (v29.4b): 資格 576 セル満票・速度場 v̂ 真値照合 24/24 (円環 Δ∞ ≤ 7%)・未使用源 τ 予言 24/24・regulator 間 23/24 — hold-12 R-A×R-C の 1 対バー外を不成立として確定 (調整なし。教訓 = バー導出標本は位相層化)**。**v29.5 = 一意性の境界を全数機械化**: n ≤ 7 連結グラフ全数の collision atlas — 静的核の厳密衝突対 (P6 vs 単環 — 半充填静的核は sign(A) まで) を発見し動的チャネルが分離 (「静的単独不可・応答併用可」の最小実例)、Petersen 誤認 (縮退→fallback 捏造)・factorization 選定不能・弱い弦の静的偽陰性を証明書化。**v29.6 = dimension-agnostic pipeline 資格**: 適応 gap 則 + clique complex + Z2 homology 自作 + link 多様体性で、熱的 Gaussian 状態の核だけから torus (1,2,1)/cylinder/disk/sphere (1,0,1) を end-to-end 同定 (隣接 欠0余0・敵対対照 5 種で捏造なし)。発見 = 臨界半充填 GS は境界増強・Friedel 共鳴で局所抽出を破る (v28 境界増強の 2D 版)。**残り**: v30.0 scoped capability certificate (qrn_core 変更 = 儀式・新鮮 holdout・能力別 certificate)。禁止継続: BridgeLawCertificate 登録・PRED-019・ProperTime/full Lorentzian の表現・バー変更。**v25 系列は凍結済み** (`v252_manifest` が常設監査。v25.3 は開始しない)。**経路 B (誘導重力監査) は完結** (v27.0: fork = 分岐 (b) external metric — graviton pole なし・c₂ は regulator 量。composite 路線は CG0–CG8 [qrn-core-v1-spec.md §6] の封印・1/Π 常用禁止)。**第二十八期 (PROMPT/10) の成果 = 意味論の閉包**: 層分類 (`core.schema.yml` — §D.3 の正名 = Matter-on-Background Adapter v1 + Metrology Suite v1)・claims 6 軸 + 昇格禁止 R1–R7 (自然の的中 0・独立外部再現 0 の機械化)・型付き契約 (`sim/src/qrn_core.rs` — 別型・空の bridge 登録簿・居住不能型。旧 QrnState → GaussianFermionState)・器械台帳 (`instruments.yml` — 認証 SHA-256 つき 22 器械)・外部再現単位 (`reproducer/` + `replications.yml` — 同一 AI は独立でない)。監査層 = v61/v151/v271–v274 ほか (suite の ALWAYS_RUN, 後段二相)。**QRN-Core v1 は「完了」= 失敗可能な問いが立つ状態であって統一理論ではない — layer: core の主張 0 件・bridge law 0 本が正直な残高**。**第二十九期 = relational geometry bridge の比較** (`bridge_candidates.yml` に実行前凍結 — 外部計量を入力せず 2 微視的模型 × 2 regulator の無調整一致が bar。PRED-019 は解析的導出まで登録しない)。旧世代の地図は docs/uft-vX.0.md 系列。

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
