% ============================================================
% rules.pl — QRN 主張依存グラフの推論規則 (v15.2)
% ============================================================
% facts.pl (claims.graph.json から生成) の上で、依存の推移閉包・深さ・
% 影響範囲・等級単調性・非循環性を導出する。これらは v151_audit (Rust) と
% **独立の実装**であり、report.pl が全数照合する (独立実装の相互検証)。

:- table depends_tc/2.
:- table depth_of/2.

% 依存の推移閉包: depends_tc(X, Y) — Y が落ちれば X も落ちる
depends_tc(X, Y) :- dep(X, Y).
depends_tc(X, Z) :- dep(X, Y), depends_tc(Y, Z).

% 深さ: 根 (依存なし) は 0、それ以外は 1 + max(依存先の深さ)
depth_of(X, 0) :- claim(X, _, _), \+ dep(X, _).
depth_of(X, D) :-
    claim(X, _, _), dep(X, _),
    aggregate_all(max(DY), (dep(X, Y), depth_of(Y, DY)), DM),
    D is DM + 1.

% 主張の被支持閉包: X が落ちると (推移的に) 落ちる主張の数 (X 自身を除く)
closure_count(X, N) :-
    claim(X, _, _),
    ( setof(C, depends_tc(C, X), L) -> length(L, N) ; N = 0 ).

% 仮定 A を抜くと落ちる主張: A を直接使う主張と、その依存閉包
falls_by_asm(A, C) :- asm_of(C, A).
falls_by_asm(A, C) :- asm_of(D, A), depends_tc(C, D).
blast_asm(A, N) :-
    assumption(A, _, _, _),
    ( setof(C, falls_by_asm(A, C), L) -> length(L, N) ; N = 0 ).

% 反証条件 F が発火すると落ちる主張
falls_by_fal(F, C) :- fal_of(C, F).
falls_by_fal(F, C) :- fal_of(D, F), depends_tc(C, D).
blast_fal(F, N) :-
    falsifier(F, _),
    ( setof(C, falls_by_fal(F, C), L) -> length(L, N) ; N = 0 ).

% 等級の強さ順位 (C3 機構と C4 現象論は横並び — ASM-EDGE-SEMANTICS)
rank(c0, 0).
rank(c1, 1).
rank(c2, 2).
rank(c3, 3).
rank(c4, 3).
rank(c5, 4).

% 等級単調性違反: 依存先の順位が自分より高い辺
mono_violation(X-Y) :-
    dep(X, Y),
    claim(X, _, LX), claim(Y, _, LY),
    rank(LX, RX), rank(LY, RY),
    RY > RX.

% 循環に関与する主張
cycle_node(X) :- claim(X, _, _), depends_tc(X, X).

% 孤児: どの主張からも参照されない局所仮定・反証条件
orphan_asm(A) :- assumption(A, _, local, _), \+ asm_of(_, A).
orphan_fal(F) :- falsifier(F, _), \+ fal_of(_, F).
