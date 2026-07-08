% ============================================================
% report.pl — QRN 依存グラフの Prolog 推論とレポート出力 (v15.2)
% ============================================================
% 実行: sh dag/run.sh  (リポジトリ直下で
%        python3 dag/json_to_facts.py && swipl -q dag/report.pl)
%
% 行うこと:
%   [1] 構造検査: 非循環・等級単調・孤児なし (Prolog 独立推論)
%   [2] Rust 監査 (v151_audit) の導出値との全数照合:
%       深さ 96+ / 被支持閉包 96+ / 仮定の影響範囲 / 反証条件の射程
%   [3] レポート出力: report.md / report.json / qrn_dag.dot / qrn_dag.mmd
% 全て決定論 (facts.pl の順序に従い、タイムスタンプなし)。
% 照合が 1 件でも崩れたら exit 1 ([FAIL])。

:- prolog_load_context(directory, D),
   atom_concat(D, '/facts.pl', F), consult(F),
   atom_concat(D, '/rules.pl', R), consult(R).

:- initialization(main, main).

% ---------- 表示ヘルパ ----------
passfail(true, '[PASS]').
passfail(false, '[FAIL]').

check(Name, Bad, Ok) :-
    ( Bad == [] -> Ok = true ; Ok = false ),
    passfail(Ok, T),
    ( Bad == []
    -> format("  ~w ~w~n", [T, Name])
    ;  format("  ~w ~w: ~w~n", [T, Name, Bad])
    ).

% ---------- 照合 ----------
mismatch_depth(X-P-R) :- claim(X, _, _), depth_of(X, P), rust_depth(X, R), P =\= R.
mismatch_closure(X-P-R) :- claim(X, _, _), closure_count(X, P), rust_closure(X, R), P =\= R.
mismatch_basm(A-P-R) :- assumption(A, _, _, _), blast_asm(A, P), rust_blast_asm(A, R), P =\= R.
mismatch_bfal(F-P-R) :- falsifier(F, _), blast_fal(F, P), rust_blast_fal(F, R), P =\= R.

% ---------- DOT 出力 ----------
level_color(c0, '"#bdbdbd"').
level_color(c1, '"#7cb342"').
level_color(c2, '"#1e88e5"').
level_color(c3, '"#fb8c00"').
level_color(c4, '"#fdd835"').
level_color(c5, '"#e53935"').

write_dot(S) :-
    format(S, "// 自動生成 (dag/report.pl) — QRN 主張依存グラフ~n", []),
    format(S, "digraph qrn {~n  rankdir=BT;~n  node [shape=box, style=filled, fontname=\"Helvetica\"];~n", []),
    forall(claim(X, V, L),
           ( level_color(L, Col),
             upcase_atom(L, LU),
             format(S, "  \"~w\" [fillcolor=~w, label=\"~w\\n~w ~w\"];~n", [X, Col, X, LU, V]) )),
    forall(dep(X, Y), format(S, "  \"~w\" -> \"~w\";~n", [X, Y])),
    format(S, "}~n", []).

% ---------- Mermaid 出力 ----------
mmd_id(X, M) :-
    atom_chars(X, Cs),
    maplist([C, D]>>( C == '-' -> D = '_' ; D = C ), Cs, Ds),
    atom_chars(M, Ds).

write_mmd(S) :-
    format(S, "%% 自動生成 (dag/report.pl) — QRN 主張依存グラフ~n", []),
    format(S, "graph BT~n", []),
    forall(claim(X, _, L),
           ( mmd_id(X, M), format(S, "  ~w[\"~w (~w)\"]:::~w~n", [M, X, L, L]) )),
    forall(dep(X, Y),
           ( mmd_id(X, MX), mmd_id(Y, MY), format(S, "  ~w --> ~w~n", [MX, MY]) )),
    format(S, "  classDef c0 fill:#bdbdbd~n", []),
    format(S, "  classDef c1 fill:#7cb342~n", []),
    format(S, "  classDef c2 fill:#1e88e5~n", []),
    format(S, "  classDef c3 fill:#fb8c00~n", []),
    format(S, "  classDef c4 fill:#fdd835~n", []),
    format(S, "  classDef c5 fill:#e53935~n", []).

% ---------- Markdown レポート ----------
write_md(S) :-
    aggregate_all(count, claim(_, _, _), NC),
    aggregate_all(count, dep(_, _), NE),
    aggregate_all(count, assumption(_, _, _, _), NA),
    aggregate_all(count, falsifier(_, _), NF),
    aggregate_all(max(D), (claim(X0, _, _), depth_of(X0, D)), MaxD),
    format(S, "# QRN 依存グラフ — Prolog 推論レポート~n~n", []),
    format(S, "**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**~n", []),
    format(S, "Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。~n~n", []),
    format(S, "主張 ~w / 依存辺 ~w / 仮定 ~w / 反証条件 ~w / 最大深さ ~w~n~n", [NC, NE, NA, NF, MaxD]),
    format(S, "## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)~n~n", []),
    format(S, "| 仮定 | type | 閉包 |~n|---|---|---|~n", []),
    findall(N-A, blast_asm(A, N), LA),
    sort(0, @>=, LA, SA),
    forall(member(N-A, SA),
           ( assumption(A, T, _, _), format(S, "| ~w | ~w | ~w |~n", [A, T, N]) )),
    format(S, "~n## 反証条件の射程 (発火すると落ちる主張の閉包 — 降順)~n~n", []),
    format(S, "| 反証条件 | status | 閉包 |~n|---|---|---|~n", []),
    findall(N-F, blast_fal(F, N), LF),
    sort(0, @>=, LF, SF),
    forall(member(N-F, SF),
           ( falsifier(F, St), format(S, "| ~w | ~w | ~w |~n", [F, St, N]) )),
    format(S, "~n## 深さ別の主張数~n~n| 深さ | 主張数 |~n|---|---|~n", []),
    forall(between(0, MaxD, D2),
           ( aggregate_all(count, (claim(X2, _, _), depth_of(X2, D2)), ND),
             format(S, "| ~w | ~w |~n", [D2, ND]) )).

% ---------- JSON レポート ----------
write_json(S) :-
    aggregate_all(count, claim(_, _, _), NC),
    aggregate_all(count, dep(_, _), NE),
    format(S, "{~n  \"generated_by\": \"dag/report.pl (swipl)\",~n", []),
    format(S, "  \"n_claims\": ~w,~n  \"n_edges\": ~w,~n", [NC, NE]),
    findall(A-N, blast_asm(A, N), LA),
    format(S, "  \"blast_asm\": {", []),
    write_pairs(S, LA),
    format(S, "},~n", []),
    findall(F-N, blast_fal(F, N), LF),
    format(S, "  \"blast_fal\": {", []),
    write_pairs(S, LF),
    format(S, "},~n", []),
    findall(X-D, (claim(X, _, _), depth_of(X, D)), LD),
    format(S, "  \"depth\": {", []),
    write_pairs(S, LD),
    format(S, "}~n}~n", []).

write_pairs(_, []).
write_pairs(S, [K-V]) :- format(S, "\"~w\": ~w", [K, V]).
write_pairs(S, [K-V, H | T]) :- format(S, "\"~w\": ~w, ", [K, V]), write_pairs(S, [H | T]).

% ---------- main ----------
main :-
    format("=== v15.2 依存グラフの Prolog 推論 (dag/report.pl) ===~n~n", []),
    % [1] 構造検査
    findall(X, cycle_node(X), Cyc),
    check('非循環性 (depends_tc に自己到達なし)', Cyc, Ok1),
    findall(V, mono_violation(V), Mono),
    check('等級単調性 rank(dep) =< rank(claim)', Mono, Ok2),
    findall(A, orphan_asm(A), OA),
    findall(F, orphan_fal(F), OF),
    append(OA, OF, Orph),
    check('孤児なし (局所仮定・反証条件は全て参照される)', Orph, Ok3),
    % [2] Rust 監査との全数照合
    findall(M, mismatch_depth(M), MD),
    check('深さ: Prolog = Rust (全主張)', MD, Ok4),
    findall(M, mismatch_closure(M), MC),
    check('被支持閉包: Prolog = Rust (全主張)', MC, Ok5),
    findall(M, mismatch_basm(M), MA),
    check('仮定の影響範囲: Prolog = Rust (全仮定)', MA, Ok6),
    findall(M, mismatch_bfal(M), MF),
    check('反証条件の射程: Prolog = Rust (全反証条件)', MF, Ok7),
    % [3] レポート出力
    setup_call_cleanup(open('dag/qrn_dag.dot', write, S1), write_dot(S1), close(S1)),
    setup_call_cleanup(open('dag/qrn_dag.mmd', write, S2), write_mmd(S2), close(S2)),
    setup_call_cleanup(open('dag/report.md', write, S3), write_md(S3), close(S3)),
    setup_call_cleanup(open('dag/report.json', write, S4), write_json(S4), close(S4)),
    format("~n  出力: dag/qrn_dag.dot, dag/qrn_dag.mmd, dag/report.md, dag/report.json~n", []),
    ( maplist(==(true), [Ok1, Ok2, Ok3, Ok4, Ok5, Ok6, Ok7])
    -> format("~n総合判定: [PASS] Prolog 独立推論は Rust 監査と全数一致~n", []), halt(0)
    ;  format("~n総合判定: [FAIL]~n", []), halt(1)
    ).
