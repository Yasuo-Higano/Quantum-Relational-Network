# QRN 依存グラフ — Prolog 推論レポート

**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**
Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。

主張 245 / 依存辺 416 / 仮定 40 / 反証条件 15 / 最大深さ 34

## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)

| 仮定 | type | 閉包 |
|---|---|---|
| ASM-LATTICE | framework | 176 |
| ASM-SEED | design | 86 |
| ASM-GAUSS | framework | 86 |
| ASM-LOWDIM | framework | 84 |
| ASM-PDG | data | 70 |
| ASM-TORUS | model | 65 |
| ASM-PRIOR | design | 65 |
| ASM-WILSON-GRID | design | 62 |
| ASM-OVERLAP | model | 62 |
| ASM-STABLE-LABEL | convention | 61 |
| ASM-DIAGPAIR | model | 60 |
| ASM-SIGMA-DATA | model | 54 |
| ASM-MODK | definition | 49 |
| ASM-KTM | model | 42 |
| ASM-GAUGE-GROUP | model | 37 |
| ASM-EDGE-SEMANTICS | design | 36 |
| ASM-ANOMALY-COEFS | data | 36 |
| ASM-IEEE754 | trust | 35 |
| ASM-CHIRALITY | model | 35 |
| ASM-ALL-CHARGED | definition | 35 |
| ASM-WINDOW-V31 | window | 30 |
| ASM-EFT-VALIDITY | model | 30 |
| ASM-LEAN-TRUST | trust | 28 |
| ASM-Z2-MINIMAL | model | 27 |
| ASM-INIT | model | 23 |
| ASM-LAYER-SEMANTICS | convention | 22 |
| ASM-DOF-GROWTH | model | 21 |
| ASM-WINDOW-U1SQ | window | 17 |
| ASM-SMCONTENT | data | 16 |
| ASM-WINDOW-EXT | window | 8 |
| ASM-WINDOW-V43 | window | 7 |
| ASM-OBS-FRACTIONAL | observational | 7 |
| ASM-NET-REAL | ontology | 7 |
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
| FAL-SUITE | open | 237 |
| FAL-CONTINUUM | open | 102 |
| FAL-CKM-OOS | open | 63 |
| FAL-BMV | open | 48 |
| FAL-CEX-WINDOW | open | 41 |
| FAL-AREALAW | open | 39 |
| FAL-EXOTIC-CHIRAL | open | 30 |
| FAL-QNEC | open | 17 |
| FAL-NEUTRINO | open | 15 |
| FAL-PAGE | open | 5 |
| FAL-LORENTZ | open | 5 |
| FAL-DIMFLOW | open | 4 |
| FAL-GLOBALSYM | open | 3 |
| FAL-COSMO | open | 2 |
| FAL-SUSY | open | 1 |

## 深さ別の主張数

| 深さ | 主張数 |
|---|---|
| 0 | 42 |
| 1 | 33 |
| 2 | 25 |
| 3 | 18 |
| 4 | 17 |
| 5 | 13 |
| 6 | 9 |
| 7 | 6 |
| 8 | 10 |
| 9 | 9 |
| 10 | 8 |
| 11 | 6 |
| 12 | 7 |
| 13 | 5 |
| 14 | 4 |
| 15 | 6 |
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
| 34 | 1 |
