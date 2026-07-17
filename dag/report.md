# QRN 依存グラフ — Prolog 推論レポート

**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**
Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。

主張 166 / 依存辺 258 / 仮定 38 / 反証条件 15 / 最大深さ 31

## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)

| 仮定 | type | 閉包 |
|---|---|---|
| ASM-LATTICE | framework | 118 |
| ASM-SEED | design | 73 |
| ASM-PDG | data | 66 |
| ASM-TORUS | model | 61 |
| ASM-PRIOR | design | 61 |
| ASM-WILSON-GRID | design | 58 |
| ASM-OVERLAP | model | 58 |
| ASM-LOWDIM | framework | 58 |
| ASM-STABLE-LABEL | convention | 57 |
| ASM-DIAGPAIR | model | 56 |
| ASM-GAUSS | framework | 55 |
| ASM-SIGMA-DATA | model | 50 |
| ASM-MODK | definition | 38 |
| ASM-KTM | model | 38 |
| ASM-GAUGE-GROUP | model | 34 |
| ASM-ANOMALY-COEFS | data | 33 |
| ASM-CHIRALITY | model | 32 |
| ASM-ALL-CHARGED | definition | 32 |
| ASM-WINDOW-V31 | window | 27 |
| ASM-EFT-VALIDITY | model | 27 |
| ASM-WINDOW-U1SQ | window | 14 |
| ASM-SMCONTENT | data | 13 |
| ASM-Z2-MINIMAL | model | 11 |
| ASM-EDGE-SEMANTICS | design | 11 |
| ASM-INIT | model | 9 |
| ASM-WINDOW-EXT | window | 8 |
| ASM-LEAN-TRUST | trust | 8 |
| ASM-WINDOW-V43 | window | 7 |
| ASM-OBS-FRACTIONAL | observational | 7 |
| ASM-NET-REAL | ontology | 7 |
| ASM-DOF-GROWTH | model | 7 |
| ASM-WINDOW-PAIR | window | 5 |
| ASM-WINDOW-U1CUBE | window | 4 |
| ASM-ORBIFOLD | model | 3 |
| ASM-WINDOW-EXC | window | 2 |
| ASM-WICK | framework | 1 |
| ASM-ENVARIANCE | framework | 1 |
| ASM-QM | framework | 0 |

## 反証条件の射程 (発火すると落ちる主張の閉包 — 降順)

| 反証条件 | status | 閉包 |
|---|---|---|
| FAL-SUITE | open | 159 |
| FAL-CONTINUUM | open | 76 |
| FAL-CKM-OOS | open | 59 |
| FAL-BMV | open | 44 |
| FAL-CEX-WINDOW | open | 38 |
| FAL-AREALAW | open | 32 |
| FAL-EXOTIC-CHIRAL | open | 27 |
| FAL-QNEC | open | 14 |
| FAL-NEUTRINO | open | 12 |
| FAL-PAGE | open | 5 |
| FAL-LORENTZ | open | 5 |
| FAL-DIMFLOW | open | 4 |
| FAL-GLOBALSYM | open | 3 |
| FAL-COSMO | open | 2 |
| FAL-SUSY | open | 1 |

## 深さ別の主張数

| 深さ | 主張数 |
|---|---|
| 0 | 39 |
| 1 | 28 |
| 2 | 16 |
| 3 | 12 |
| 4 | 9 |
| 5 | 4 |
| 6 | 2 |
| 7 | 3 |
| 8 | 6 |
| 9 | 6 |
| 10 | 5 |
| 11 | 4 |
| 12 | 4 |
| 13 | 1 |
| 14 | 1 |
| 15 | 3 |
| 16 | 4 |
| 17 | 1 |
| 18 | 2 |
| 19 | 1 |
| 20 | 1 |
| 21 | 3 |
| 22 | 2 |
| 23 | 1 |
| 24 | 1 |
| 25 | 1 |
| 26 | 1 |
| 27 | 1 |
| 28 | 1 |
| 29 | 1 |
| 30 | 1 |
| 31 | 1 |
