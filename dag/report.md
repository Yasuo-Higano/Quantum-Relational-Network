# QRN 依存グラフ — Prolog 推論レポート

**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**
Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。

主張 171 / 依存辺 272 / 仮定 38 / 反証条件 15 / 最大深さ 32

## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)

| 仮定 | type | 閉包 |
|---|---|---|
| ASM-LATTICE | framework | 123 |
| ASM-SEED | design | 74 |
| ASM-PDG | data | 67 |
| ASM-LOWDIM | framework | 63 |
| ASM-TORUS | model | 62 |
| ASM-PRIOR | design | 62 |
| ASM-GAUSS | framework | 60 |
| ASM-WILSON-GRID | design | 59 |
| ASM-OVERLAP | model | 59 |
| ASM-STABLE-LABEL | convention | 58 |
| ASM-DIAGPAIR | model | 57 |
| ASM-SIGMA-DATA | model | 51 |
| ASM-MODK | definition | 43 |
| ASM-KTM | model | 39 |
| ASM-GAUGE-GROUP | model | 35 |
| ASM-ANOMALY-COEFS | data | 34 |
| ASM-CHIRALITY | model | 33 |
| ASM-ALL-CHARGED | definition | 33 |
| ASM-WINDOW-V31 | window | 28 |
| ASM-EFT-VALIDITY | model | 28 |
| ASM-WINDOW-U1SQ | window | 15 |
| ASM-SMCONTENT | data | 14 |
| ASM-Z2-MINIMAL | model | 12 |
| ASM-EDGE-SEMANTICS | design | 12 |
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
| FAL-SUITE | open | 164 |
| FAL-CONTINUUM | open | 81 |
| FAL-CKM-OOS | open | 60 |
| FAL-BMV | open | 45 |
| FAL-CEX-WINDOW | open | 39 |
| FAL-AREALAW | open | 37 |
| FAL-EXOTIC-CHIRAL | open | 28 |
| FAL-QNEC | open | 15 |
| FAL-NEUTRINO | open | 13 |
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
