# QRN 依存グラフ — Prolog 推論レポート

**このファイルは `sh dag/run.sh` が生成する。手で編集しない。**
Prolog (swipl) による独立推論であり、Rust 監査 `v151_audit` の導出値と全数照合済み。

主張 128 / 依存辺 187 / 仮定 38 / 反証条件 15 / 最大深さ 22

## 仮定の影響範囲 (抜くと落ちる主張の閉包 — 降順)

| 仮定 | type | 閉包 |
|---|---|---|
| ASM-LATTICE | framework | 82 |
| ASM-SEED | design | 63 |
| ASM-PDG | data | 55 |
| ASM-TORUS | model | 50 |
| ASM-PRIOR | design | 50 |
| ASM-WILSON-GRID | design | 47 |
| ASM-OVERLAP | model | 47 |
| ASM-STABLE-LABEL | convention | 46 |
| ASM-DIAGPAIR | model | 46 |
| ASM-SIGMA-DATA | model | 39 |
| ASM-LOWDIM | framework | 36 |
| ASM-GAUSS | framework | 33 |
| ASM-KTM | model | 28 |
| ASM-GAUGE-GROUP | model | 24 |
| ASM-ANOMALY-COEFS | data | 23 |
| ASM-CHIRALITY | model | 22 |
| ASM-ALL-CHARGED | definition | 22 |
| ASM-WINDOW-V31 | window | 17 |
| ASM-EFT-VALIDITY | model | 17 |
| ASM-MODK | definition | 16 |
| ASM-INIT | model | 9 |
| ASM-WINDOW-EXT | window | 8 |
| ASM-WINDOW-V43 | window | 7 |
| ASM-OBS-FRACTIONAL | observational | 7 |
| ASM-NET-REAL | ontology | 7 |
| ASM-LEAN-TRUST | trust | 7 |
| ASM-DOF-GROWTH | model | 7 |
| ASM-Z2-MINIMAL | model | 5 |
| ASM-WINDOW-PAIR | window | 5 |
| ASM-EDGE-SEMANTICS | design | 5 |
| ASM-WINDOW-U1SQ | window | 4 |
| ASM-WINDOW-U1CUBE | window | 4 |
| ASM-SMCONTENT | data | 3 |
| ASM-ORBIFOLD | model | 3 |
| ASM-WINDOW-EXC | window | 2 |
| ASM-WICK | framework | 1 |
| ASM-ENVARIANCE | framework | 1 |
| ASM-QM | framework | 0 |

## 反証条件の射程 (発火すると落ちる主張の閉包 — 降順)

| 反証条件 | status | 閉包 |
|---|---|---|
| FAL-SUITE | open | 122 |
| FAL-CKM-OOS | open | 49 |
| FAL-CONTINUUM | open | 40 |
| FAL-BMV | open | 33 |
| FAL-CEX-WINDOW | open | 28 |
| FAL-EXOTIC-CHIRAL | open | 17 |
| FAL-AREALAW | open | 10 |
| FAL-QNEC | open | 8 |
| FAL-PAGE | open | 5 |
| FAL-LORENTZ | open | 4 |
| FAL-DIMFLOW | open | 4 |
| FAL-GLOBALSYM | open | 3 |
| FAL-NEUTRINO | open | 2 |
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
| 22 | 1 |
