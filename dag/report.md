# QRN 依存グラフ — Prolog 推論レポート

**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**
Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。

主張 132 / 依存辺 191 / 仮定 38 / 反証条件 15 / 最大深さ 25

## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)

| 仮定 | type | 閉包 |
|---|---|---|
| ASM-LATTICE | framework | 86 |
| ASM-SEED | design | 67 |
| ASM-PDG | data | 59 |
| ASM-TORUS | model | 54 |
| ASM-PRIOR | design | 54 |
| ASM-WILSON-GRID | design | 51 |
| ASM-OVERLAP | model | 51 |
| ASM-STABLE-LABEL | convention | 50 |
| ASM-DIAGPAIR | model | 50 |
| ASM-SIGMA-DATA | model | 43 |
| ASM-LOWDIM | framework | 36 |
| ASM-GAUSS | framework | 33 |
| ASM-KTM | model | 32 |
| ASM-GAUGE-GROUP | model | 28 |
| ASM-ANOMALY-COEFS | data | 27 |
| ASM-CHIRALITY | model | 26 |
| ASM-ALL-CHARGED | definition | 26 |
| ASM-WINDOW-V31 | window | 21 |
| ASM-EFT-VALIDITY | model | 21 |
| ASM-MODK | definition | 16 |
| ASM-INIT | model | 9 |
| ASM-WINDOW-U1SQ | window | 8 |
| ASM-WINDOW-EXT | window | 8 |
| ASM-WINDOW-V43 | window | 7 |
| ASM-SMCONTENT | data | 7 |
| ASM-OBS-FRACTIONAL | observational | 7 |
| ASM-NET-REAL | ontology | 7 |
| ASM-LEAN-TRUST | trust | 7 |
| ASM-DOF-GROWTH | model | 7 |
| ASM-Z2-MINIMAL | model | 5 |
| ASM-WINDOW-PAIR | window | 5 |
| ASM-EDGE-SEMANTICS | design | 5 |
| ASM-WINDOW-U1CUBE | window | 4 |
| ASM-ORBIFOLD | model | 3 |
| ASM-WINDOW-EXC | window | 2 |
| ASM-WICK | framework | 1 |
| ASM-ENVARIANCE | framework | 1 |
| ASM-QM | framework | 0 |

## 反証条件の射程 (発火すると落ちる主張の閉包 — 降順)

| 反証条件 | status | 閉包 |
|---|---|---|
| FAL-SUITE | open | 126 |
| FAL-CKM-OOS | open | 53 |
| FAL-CONTINUUM | open | 44 |
| FAL-BMV | open | 37 |
| FAL-CEX-WINDOW | open | 32 |
| FAL-EXOTIC-CHIRAL | open | 21 |
| FAL-AREALAW | open | 10 |
| FAL-QNEC | open | 8 |
| FAL-NEUTRINO | open | 6 |
| FAL-PAGE | open | 5 |
| FAL-LORENTZ | open | 4 |
| FAL-DIMFLOW | open | 4 |
| FAL-GLOBALSYM | open | 3 |
| FAL-COSMO | open | 2 |
| FAL-SUSY | open | 1 |

## 深さ別の主張数

| 深さ | 主張数 |
|---|---|
| 0 | 35 |
| 1 | 22 |
| 2 | 13 |
| 3 | 7 |
| 4 | 7 |
| 5 | 3 |
| 6 | 1 |
| 7 | 2 |
| 8 | 5 |
| 9 | 5 |
| 10 | 4 |
| 11 | 3 |
| 12 | 3 |
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
