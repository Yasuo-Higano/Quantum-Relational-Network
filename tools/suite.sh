#!/bin/sh
# tools/suite.sh — 全スイートの増分実行ランナー (PROMPT/5「計算結果引用」2026-07-20)
#
# 各バイナリのソース (sim/src/bin/<bin>.rs) と共有部 (sim/src/lib.rs, sim/Cargo.toml) の
# SHA-256 を台帳 results/suite_manifest.tsv に記録し、不変なバイナリの再計算を省いて
# 前回実行の結果を「引用」する。README の「実行 113 + 引用 12」の区分をハッシュ台帳で
# 汎用化したもの。完全再走 (固定シード決定性・末桁ドリフトの儀式 — v153/v154 で実績)
# は数期に一度 `full` モードで行う。
#
# 使い方 (リポジトリのルートから。通常は Makefile 経由):
#   sh tools/suite.sh status              判定のみ表示 (何も実行しない)
#   sh tools/suite.sh incremental [OUT]   増分実行   (既定 OUT=results/suite_incremental.txt)
#   sh tools/suite.sh full [OUT]          完全再計算 (既定 OUT=results/suite_full.txt)
#   sh tools/suite.sh seed SUITE_FILE     既存の完全スイート出力から台帳を初期化
#                                         (完全再走の直後・ソース不変のときのみ正当)
#   環境変数 JOBS=N — 並列実行数 (既定 1)。各バイナリは独立プロセスなので結果は
#   並列度に依らない (PROMPT/4)。ただし sim/cache が冷えた初回は v16xx/v17xx 系が
#   キャッシュ生成で競合し得るため、クリーン環境での初回フル再走は JOBS=1 を推奨。
#
# 引用の判定規則 (incremental):
#   [1] lib.rs / Cargo.toml のハッシュ変化 → 全数再実行 (共有部の波及は個別判定不能)
#   [2] 台帳に無い / bin ソースのハッシュ変化 → 再実行
#   [3] 前回 FAIL または exit≠0 → 再実行 (FAIL の引用は無意味)
#   [4] リポジトリ状態 (claims.yml, docs/, results/*.json, explore/*.json) を入力に
#       読む監査・照合層 (ALWAYS_RUN) → ソース不変でも常に再実行 (計 ~16 秒)
#   [5] それ以外 → 引用 (前回結果の PASS 数を集計に算入)
#   rustc のバージョン変化は警告のみ (数値ドリフト検査は full の儀式の役割)。

set -eu

# リポジトリ状態を入力に読むため、ソース不変でも引用できない監査・照合層
ALWAYS_RUN="v61_ledger v151_audit v213_dmrgaudit v214_bridgeaudit v217_fiveconditions v221_dmrgex v235_lambdainf"

MANIFEST=results/suite_manifest.tsv
BINDIR=sim/src/bin
REL=sim/target/release

mode=${1:-status}

[ -f sim/Cargo.toml ] || { echo "エラー: リポジトリのルートから実行してください" >&2; exit 2; }

sha() {
    if command -v sha256sum >/dev/null 2>&1; then sha256sum "$1" | awk '{print $1}'
    else shasum -a 256 "$1" | awk '{print $1}'; fi
}

count() { grep -o "$1" "$2" 2>/dev/null | wc -l | tr -d ' '; }

mtime_date() {
    if stat -f %m "$1" >/dev/null 2>&1; then date -r "$(stat -f %m "$1")" +%Y-%m-%d
    else date -d "@$(stat -c %Y "$1")" +%Y-%m-%d; fi
}

in_always_run() {
    case " $ALWAYS_RUN " in *" $1 "*) return 0 ;; *) return 1 ;; esac
}

manifest_header() {
    printf '# 全スイート実行台帳 — tools/suite.sh が生成・更新 (PROMPT/5)\n'
    printf '# GLOBAL\tlib.rs sha256\tCargo.toml sha256\trustc\n'
    printf '# <bin>\tsrc sha256\tPASS\tFAIL\texit\t秒\t実行日\t結果ファイル\n'
}

# ---- run-one: 1 バイナリを実行して結果とメタを書く (親から xargs/ループで呼ばれる) ----
if [ "$mode" = run-one ]; then
    bin=${2:?run-one にはバイナリ名が要る}
    : "${SUITE_TMP:?run-one は suite.sh 本体から呼ばれる想定}"
    t0=$(date +%s)
    # stderr も保存する — lib.rs の self_test は eprintln! で、既存の結果ファイル
    # (v240 スイート含む) は stderr 込みの慣行
    if "$REL/$bin" > "results/$bin.txt" 2>&1; then ec=0; else ec=$?; fi
    t1=$(date +%s)
    pass=$(count '\[PASS\]' "results/$bin.txt")
    fail=$(count '\[FAIL\]' "results/$bin.txt")
    printf '%s\t%s\t%s\t%s\t%s\t%s\n' "$bin" "$pass" "$fail" "$ec" "$((t1 - t0))" "$(date +%Y-%m-%d)" > "$SUITE_TMP/meta_$bin"
    printf '[実行] %-24s PASS %-4s FAIL %-2s exit=%s (%ss)\n' "$bin" "$pass" "$fail" "$ec" "$((t1 - t0))"
    exit 0
fi

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# ---- 現在のハッシュ ----
lib_sha=$(sha sim/src/lib.rs)
cargo_sha=$(sha sim/Cargo.toml)
rustc_ver=$(rustc -V 2>/dev/null || echo "rustc 不明")

bins=""
for f in "$BINDIR"/v*.rs; do
    b=$(basename "$f" .rs)
    bins="$bins $b"
    printf '%s\t%s\n' "$b" "$(sha "$f")" >> "$tmp/hashes.tsv"
done

cur_sha() { awk -F'\t' -v b="$1" '$1 == b { print $2; exit }' "$tmp/hashes.tsv"; }
old_row() {
    if [ -f "$MANIFEST" ]; then awk -F'\t' -v b="$1" '$1 == b { print; exit }' "$MANIFEST"; fi
}

# ---- seed: 既存の完全スイート出力から台帳を初期化 ----
if [ "$mode" = seed ]; then
    suite_file=${2:?使い方: sh tools/suite.sh seed <full_suite.txt>}
    [ -f "$suite_file" ] || { echo "エラー: $suite_file が無い" >&2; exit 2; }
    echo "前提: $suite_file が現在のソース状態での完全再走であること (直後の seed のみ正当)"
    # 区画ごとの PASS/FAIL/exit を数える。[引用] 区画 (増分集約) は実行記録でないので捨てる。
    awk '
        /^===== sim\/target\/release\// {
            if (index($0, "[引用]") > 0) { cur = ""; next }
            b = $2; sub(/^sim\/target\/release\//, "", b)
            cur = b; p[cur] += 0; f[cur] += 0; e[cur] = "?"; next
        }
        /^----- exit=/ { x = $0; sub(/^----- exit=/, "", x); sub(/ -----.*/, "", x); e[cur] = x; cur = ""; next }
        /^----- cite -----/ { cur = ""; next }
        cur != "" { p[cur] += gsub(/\[PASS\]/, "&"); f[cur] += gsub(/\[FAIL\]/, "&") }
        END { for (b in p) printf "%s\t%d\t%d\t%s\n", b, p[b], f[b], e[b] }
    ' "$suite_file" | sort > "$tmp/sections.tsv"
    [ -s "$tmp/sections.tsv" ] || { echo "エラー: $suite_file に実行区画が無い" >&2; exit 2; }
    d=$(mtime_date "$suite_file")
    {
        manifest_header
        printf 'GLOBAL\t%s\t%s\t%s\n' "$lib_sha" "$cargo_sha" "$rustc_ver"
        for b in $bins; do
            row=$(awk -F'\t' -v b="$b" '$1 == b { print; exit }' "$tmp/sections.tsv")
            if [ -n "$row" ]; then
                pass=$(printf '%s\n' "$row" | cut -f2)
                fail=$(printf '%s\n' "$row" | cut -f3)
                ec=$(printf '%s\n' "$row" | cut -f4)
                printf '%s\t%s\t%s\t%s\t%s\t-\t%s\t%s\n' "$b" "$(cur_sha "$b")" "$pass" "$fail" "$ec" "$d" "$suite_file"
            else
                echo "警告: $suite_file に $b の区画が無い — 初回の増分実行で実行される" >&2
            fi
        done
    } > "$MANIFEST"
    awk -F'\t' '!/^#/ && $1 != "GLOBAL" { n++; p += $3; f += $4 } END {
        printf "台帳を初期化: %d 本 (PASS %d / FAIL %d) → %s\n", n, p, f, "'"$MANIFEST"'" }' "$MANIFEST"
    exit 0
fi

# ---- 実行計画 (status / incremental / full 共通) ----
global_changed=""
rustc_note=""
if [ -f "$MANIFEST" ]; then
    g=$(awk -F'\t' '$1 == "GLOBAL" { print; exit }' "$MANIFEST")
    old_lib=$(printf '%s\n' "$g" | cut -f2)
    old_cargo=$(printf '%s\n' "$g" | cut -f3)
    old_rustc=$(printf '%s\n' "$g" | cut -f4-)
    [ "$old_lib" = "$lib_sha" ] || global_changed="lib.rs 変更"
    [ "$old_cargo" = "$cargo_sha" ] || global_changed="${global_changed:+$global_changed, }Cargo.toml 変更"
    if [ "$old_rustc" != "$rustc_ver" ]; then
        rustc_note="警告: rustc が台帳と違う (台帳: $old_rustc / 現在: $rustc_ver) — 末桁ドリフト検査には make suite-full を推奨"
    fi
else
    global_changed="台帳なし (初回)"
fi

for b in $bins; do
    act=cite
    reason=""
    if [ "$mode" = full ]; then
        act=run; reason="完全再計算"
    elif [ -n "$global_changed" ]; then
        act=run; reason="$global_changed"
    elif in_always_run "$b"; then
        act=run; reason="監査層 (リポジトリ状態を読むため常時実行)"
    else
        row=$(old_row "$b")
        if [ -z "$row" ]; then
            act=run; reason="台帳に無い"
        else
            osha=$(printf '%s\n' "$row" | cut -f2)
            ofail=$(printf '%s\n' "$row" | cut -f4)
            oexit=$(printf '%s\n' "$row" | cut -f5)
            if [ "$osha" != "$(cur_sha "$b")" ]; then
                act=run; reason="ソース変更"
            elif [ "$ofail" != 0 ] || [ "$oexit" != 0 ]; then
                act=run; reason="前回 FAIL/exit≠0"
            fi
        fi
    fi
    printf '%s\t%s\t%s\n' "$b" "$act" "$reason" >> "$tmp/plan.tsv"
done

n_run=$(awk -F'\t' '$2 == "run" { n++ } END { print n + 0 }' "$tmp/plan.tsv")
n_cite=$(awk -F'\t' '$2 == "cite" { n++ } END { print n + 0 }' "$tmp/plan.tsv")

if [ "$mode" = status ]; then
    [ -n "$rustc_note" ] && echo "$rustc_note"
    awk -F'\t' '{ printf "%s %-24s %s\n", ($2 == "run" ? "[実行]" : "[引用]"), $1, $3 }' "$tmp/plan.tsv"
    echo "---- 実行 $n_run 本 / 引用 $n_cite 本 (台帳: $MANIFEST)"
    exit 0
fi

[ "$mode" = incremental ] || [ "$mode" = full ] || { echo "エラー: 不明なモード $mode" >&2; exit 2; }

out=${2:-}
if [ -z "$out" ]; then
    if [ "$mode" = full ]; then out=results/suite_full.txt; else out=results/suite_incremental.txt; fi
fi

# ---- 実行フェーズ ----
awk -F'\t' '$2 == "run" { print $1 }' "$tmp/plan.tsv" > "$tmp/runlist"
while IFS= read -r b; do
    [ -x "$REL/$b" ] || { echo "エラー: $REL/$b が無い — 先に make build" >&2; exit 2; }
done < "$tmp/runlist"

[ -n "$rustc_note" ] && echo "$rustc_note"
JOBS=${JOBS:-1}
echo "実行 $n_run 本 / 引用 $n_cite 本 (JOBS=$JOBS) → $out"
SUITE_TMP=$tmp
export SUITE_TMP
if [ "$JOBS" -gt 1 ]; then
    xargs -P "$JOBS" -n1 sh tools/suite.sh run-one < "$tmp/runlist"
else
    while IFS= read -r b; do sh tools/suite.sh run-one "$b"; done < "$tmp/runlist"
fi

# ---- 集約と台帳の更新 ----
newman=$tmp/manifest.new
{
    manifest_header
    printf 'GLOBAL\t%s\t%s\t%s\n' "$lib_sha" "$cargo_sha" "$rustc_ver"
} > "$newman"

: > "$out.tmp"
run_pass=0; run_fail=0; cite_pass=0; cite_fail=0; bad=""
for b in $bins; do
    act=$(awk -F'\t' -v b="$b" '$1 == b { print $2; exit }' "$tmp/plan.tsv")
    if [ "$act" = run ]; then
        [ -f "$tmp/meta_$b" ] || { echo "エラー: $b の実行メタが無い (実行フェーズの異常)" >&2; exit 2; }
        meta=$(cat "$tmp/meta_$b")
        pass=$(printf '%s\n' "$meta" | cut -f2)
        fail=$(printf '%s\n' "$meta" | cut -f3)
        ec=$(printf '%s\n' "$meta" | cut -f4)
        secs=$(printf '%s\n' "$meta" | cut -f5)
        d=$(printf '%s\n' "$meta" | cut -f6)
        printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\tresults/%s.txt\n' \
            "$b" "$(cur_sha "$b")" "$pass" "$fail" "$ec" "$secs" "$d" "$b" >> "$newman"
        {
            printf '===== %s/%s =====\n' "$REL" "$b"
            cat "results/$b.txt"
            printf -- '----- exit=%s -----\n' "$ec"
        } >> "$out.tmp"
        run_pass=$((run_pass + pass)); run_fail=$((run_fail + fail))
        if [ "$ec" != 0 ] || [ "$fail" != 0 ]; then bad="$bad $b"; fi
    else
        row=$(old_row "$b")
        printf '%s\n' "$row" >> "$newman"
        osha=$(printf '%s\n' "$row" | cut -f2)
        pass=$(printf '%s\n' "$row" | cut -f3)
        fail=$(printf '%s\n' "$row" | cut -f4)
        d=$(printf '%s\n' "$row" | cut -f7)
        ref=$(printf '%s\n' "$row" | cut -f8)
        short=$(printf '%s\n' "$osha" | cut -c1-12)
        {
            printf '===== %s/%s ===== [引用]\n' "$REL" "$b"
            printf '[引用] ソース不変 (src sha256 %s…) — %s の実行結果 %s を引用 (PASS %s / FAIL %s)\n' \
                "$short" "$d" "$ref" "$pass" "$fail"
            printf -- '----- cite -----\n'
        } >> "$out.tmp"
        cite_pass=$((cite_pass + pass)); cite_fail=$((cite_fail + fail))
    fi
done

{
    printf '\n===== 集計 =====\n'
    printf '実行 %s 本: PASS %s / FAIL %s\n' "$n_run" "$run_pass" "$run_fail"
    printf '引用 %s 本: PASS %s / FAIL %s (ソース不変を sha256 で確認 — 台帳 %s)\n' "$n_cite" "$cite_pass" "$cite_fail" "$MANIFEST"
    printf '総計: PASS %s / FAIL %s\n' "$((run_pass + cite_pass))" "$((run_fail + cite_fail))"
    printf 'lib.rs %s / Cargo.toml %s / %s\n' "$lib_sha" "$cargo_sha" "$rustc_ver"
    if [ -n "$bad" ]; then printf 'FAIL または exit≠0:%s\n' "$bad"; fi
} >> "$out.tmp"

mv "$out.tmp" "$out"
mv "$newman" "$MANIFEST"

echo "---- 実行 $n_run 本 (PASS $run_pass) + 引用 $n_cite 本 (PASS $cite_pass) = 総計 PASS $((run_pass + cite_pass)) / FAIL $((run_fail + cite_fail))"
echo "---- 集約: $out / 台帳: $MANIFEST"
if [ -n "$bad" ]; then
    echo "FAIL または exit≠0:$bad" >&2
    exit 1
fi
