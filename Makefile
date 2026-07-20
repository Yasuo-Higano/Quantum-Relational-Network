# 全スイートの実行 (PROMPT/5「計算結果引用」)
#
#   make suite                              増分: 変更されたバイナリ + 監査層のみ実行、
#                                           ソース不変のバイナリは前回結果を引用
#   make suite-full OUT=results/v250_full_suite.txt
#                                           完全再計算 — 数期に一度の儀式 (固定シード
#                                           決定性・lib.rs 波及・末桁ドリフトの検査)
#   make suite-status                       実行/引用の判定だけ表示 (何も走らせない)
#   make seed SUITE_FILE=results/v240_full_suite.txt
#                                           完全再走の直後 (ソース不変) に台帳を初期化
#
#   JOBS=8 で独立バイナリを並列実行 (各バイナリは独立プロセスなので結果は並列度に
#   依らない。sim/cache が冷えた初回のみ JOBS=1 を推奨 — tools/suite.sh 冒頭を参照)
#
# 台帳: results/suite_manifest.tsv (bin ごとの src sha256・PASS/FAIL・実行日・結果ファイル)

SHELL := /bin/sh
JOBS ?= 1
OUT ?=
SUITE_FILE ?= results/v240_full_suite.txt

.PHONY: build suite suite-full suite-status seed

build:
	cd sim && cargo build --release

suite: build
	JOBS=$(JOBS) sh tools/suite.sh incremental $(OUT)

suite-full: build
	JOBS=$(JOBS) sh tools/suite.sh full $(OUT)

suite-status:
	sh tools/suite.sh status

seed:
	sh tools/suite.sh seed $(SUITE_FILE)
