# QRN 依存グラフ — Prolog 推論レポート

**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**
Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。

主張 155 / 依存辺 231 / 仮定 38 / 反証条件 15 / 最大深さ 29

## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)

| 仮定 | type | 閉包 |
|---|---|---|
| ASM-LATTICE | framework | 107 |
| ASM-SEED | design | 71 |
| ASM-PDG | data | 64 |
| ASM-TORUS | model | 59 |
| ASM-PRIOR | design | 59 |
| ASM-WILSON-GRID | design | 56 |
| ASM-OVERLAP | model | 56 |
| ASM-STABLE-LABEL | convention | 55 |
| ASM-DIAGPAIR | model | 54 |
| ASM-SIGMA-DATA | model | 48 |
| ASM-LOWDIM | framework | 48 |
| ASM-GAUSS | framework | 45 |
| ASM-KTM | model | 36 |
| ASM-GAUGE-GROUP | model | 32 |
| ASM-ANOMALY-COEFS | data | 31 |
| ASM-CHIRALITY | model | 30 |
| ASM-ALL-CHARGED | definition | 30 |
| ASM-MODK | definition | 28 |
| ASM-WINDOW-V31 | window | 25 |
| ASM-EFT-VALIDITY | model | 25 |
| ASM-WINDOW-U1SQ | window | 12 |
| ASM-SMCONTENT | data | 11 |
| ASM-Z2-MINIMAL | model | 9 |
| ASM-INIT | model | 9 |
| ASM-EDGE-SEMANTICS | design | 9 |
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
| FAL-SUITE | open | 148 |
| FAL-CONTINUUM | open | 65 |
| FAL-CKM-OOS | open | 57 |
| FAL-BMV | open | 42 |
| FAL-CEX-WINDOW | open | 36 |
| FAL-EXOTIC-CHIRAL | open | 25 |
| FAL-AREALAW | open | 22 |
| FAL-QNEC | open | 12 |
| FAL-NEUTRINO | open | 10 |
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
| 2 | 15 |
| 3 | 11 |
| 4 | 9 |
| 5 | 4 |
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
| 26 | 1 |
| 27 | 1 |
| 28 | 1 |
| 29 | 1 |
