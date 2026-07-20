# 再現性ガイド (REPRODUCIBILITY)

本リポジトリの全数値は再現可能である。文書中の数値の一次ソースは `results/` に
保存された各バイナリの stdout であり、以下の手順で再生成できる。

## 環境

- Rust (stable)。開発時: `rustc 1.94.0`。**外部クレート依存なし** (std のみ) なので
  ネットワーク不要でビルドできる。数値ライブラリ (乱数・固有値分解・Bessel・SHA-256 等)
  は `sim/src/lib.rs` に自作し、起動時の `self_test()` で厳密値と照合される。
- OS 依存性なし (f64 演算のみ。乱数は自作 xorshift64* で全プラットフォーム同一)。

## ビルドと実行

```bash
cd sim
cargo build --release                    # 全バイナリ (外部依存なし)
./target/release/v62_atlas               # 単一実験の実行
cd ..
make suite                               # 全スイート増分 (ソース不変のバイナリは前回結果を引用)
make suite-full OUT=results/vXX0_full_suite.txt  # 完全再計算 (数期に一度の決定性・ドリフト検査)
```

増分実行の判定は台帳 `results/suite_manifest.tsv` (bin ソースと lib.rs / Cargo.toml の
SHA-256・前回の PASS/FAIL・実行日・結果ファイル) で行う。lib.rs が変わると全数再実行、
リポジトリ状態 (claims.yml, docs/, results/*.json) を入力に読む監査層は常時再実行される。
詳細は `tools/suite.sh` 冒頭のコメントを参照。
(旧来の `for b in target/release/v*; do $b; done` は .d 依存ファイルを拾って exit=126 を
混入させるため廃止 — v23.0 §5 の記録参照。)

## 決定性

- 乱数は全て固定シード (`Rng::new(seed)`)。シードは各バイナリのソース冒頭付近に明記。
- モンテカルロ系は再実行で bit 単位で同一の出力になる。
- 例外: 実行時間表示 (`ms`) のみ環境依存。

## 成果物の種類

| 種類 | 場所 | 内容 |
|---|---|---|
| stdout スナップショット | `results/*.txt` | 文書中の全数値の一次ソース |
| 機械可読な結果 | `results/*.json` | claim_id・誤差予算・対照の期待/観測 (v6.3 以降) |
| 探索の証明書 | `certificates/*.json`, `certificates/v62_sha256.txt` | 探索領域の定義・全解の正準形・解集合の SHA-256 |
| 主張台帳 | `claims.yml` | 全主張の等級 (C0–C5)・証拠・限界 |

## 検証手順

1. **主張台帳**: `./target/release/v61_ledger` — 証拠ファイルの実在・全バイナリの
   被覆などを機械検査 (総合 [PASS] で終了コード 0)。
2. **各実験**: 各バイナリは厳密解・観測値・対照との比較を `[PASS]`/`[FAIL]` として
   内蔵する。監査系 (v6.2+) は FAIL 時に非ゼロで終了する。
3. **探索の証明書**: `v62_atlas` を再実行し、出力される SHA-256 が
   `certificates/v62_sha256.txt` と一致することを確認する (解集合の同一性)。

## 誤差の考え方

- 許容誤差は経験値ではなく予算 (微分打ち切り + 有限サイズ + 丸め) として分解する
  (v6.3, `docs/uft-v6.3.md`)。
- 探索系は浮動小数を使わない (i64 厳密整数、v6.2)。
- 裾の積分など素朴 MC が破綻する場合は点推定を捨てて保守的な上界で置き換え、
  その事実を結果に記録する (v6.5)。

## 実行時間の目安 (Apple Silicon, release)

ほとんどのバイナリは数秒以内。重いもの: `v63_qnec_budget` ~45 秒 (N=402 の厳密対角化)。
