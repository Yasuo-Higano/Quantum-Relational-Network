# QRN 依存グラフ — Prolog 推論レポート

**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**
Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。

主張 184 / 依存辺 299 / 仮定 39 / 反証条件 15 / 最大深さ 33

## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)

| 仮定 | type | 閉包 |
|---|---|---|
| ASM-LATTICE | framework | 132 |
| ASM-SEED | design | 75 |
| ASM-PDG | data | 68 |
| ASM-LOWDIM | framework | 64 |
| ASM-TORUS | model | 63 |
| ASM-PRIOR | design | 63 |
| ASM-GAUSS | framework | 61 |
| ASM-WILSON-GRID | design | 60 |
| ASM-OVERLAP | model | 60 |
| ASM-STABLE-LABEL | convention | 59 |
| ASM-DIAGPAIR | model | 58 |
| ASM-SIGMA-DATA | model | 52 |
| ASM-MODK | definition | 48 |
| ASM-KTM | model | 40 |
| ASM-GAUGE-GROUP | model | 36 |
| ASM-ANOMALY-COEFS | data | 35 |
| ASM-CHIRALITY | model | 34 |
| ASM-ALL-CHARGED | definition | 34 |
| ASM-WINDOW-V31 | window | 29 |
| ASM-EFT-VALIDITY | model | 29 |
| ASM-WINDOW-U1SQ | window | 16 |
| ASM-SMCONTENT | data | 15 |
| ASM-Z2-MINIMAL | model | 13 |
| ASM-EDGE-SEMANTICS | design | 13 |
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
| ASM-IEEE754 | trust | 2 |
| ASM-WICK | framework | 1 |
| ASM-ENVARIANCE | framework | 1 |
| ASM-QM | framework | 0 |

## 反証条件の射程 (発火すると落ちる主張の閉包 — 降順)

| 反証条件 | status | 閉包 |
|---|---|---|
| FAL-SUITE | open | 176 |
| FAL-CONTINUUM | open | 87 |
| FAL-CKM-OOS | open | 61 |
| FAL-BMV | open | 46 |
| FAL-CEX-WINDOW | open | 40 |
| FAL-AREALAW | open | 38 |
| FAL-EXOTIC-CHIRAL | open | 29 |
| FAL-QNEC | open | 16 |
| FAL-NEUTRINO | open | 14 |
| FAL-PAGE | open | 5 |
| FAL-LORENTZ | open | 5 |
| FAL-DIMFLOW | open | 4 |
| FAL-GLOBALSYM | open | 3 |
| FAL-COSMO | open | 2 |
| FAL-SUSY | open | 1 |

## 深さ別の主張数

| 深さ | 主張数 |
|---|---|
| 0 | 41 |
| 1 | 30 |
| 2 | 20 |
| 3 | 13 |
| 4 | 11 |
| 5 | 5 |
| 6 | 2 |
| 7 | 3 |
| 8 | 6 |
| 9 | 6 |
| 10 | 5 |
| 11 | 4 |
| 12 | 4 |
| 13 | 2 |
| 14 | 2 |
| 15 | 4 |
| 16 | 5 |
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
| 32 | 1 |
| 33 | 1 |
