# QRN 証拠マトリクス (evidence matrix)

**このファイルは `v151_audit --write` が生成する。手で編集しない。**
機械可読版は [claims.graph.json](claims.graph.json)、辺の定義は [claims.graph.yml](claims.graph.yml)。

主張 111 件 / 依存辺 136 本 / 仮定 38 件 / 反証条件 15 件。
等級順位 C0 < C1 < C2 < {C3,C4} < C5 の単調性・非循環性は CI で機械検証される。

## 主張 × 証拠・依存

凡例: 証拠 C=code R=result D=doc L=Lean。「被支持」= この主張が落ちると (依存閉包で) 落ちる主張の数。

| ID | 等級 | 版 | 証拠 | deps | asm | fal | 被支持 |
|---|---|---|---|---|---|---|---|
| QRN-C0-001 | C0 | v0.7 | D | 0 | 0 | 0 | 6 |
| QRN-C0-002 | C0 | v4.1 | D | 0 | 0 | 0 | 8 |
| QRN-C0-003 | C0 | v0.8 | D | 0 | 0 | 0 | 1 |
| QRN-C0-004 | C0 | v3.4 | D | 0 | 0 | 0 | 3 |
| QRN-C0-005 | C0 | v2.3 | D | 0 | 0 | 0 | 49 |
| QRN-C0-006 | C0 | v3.2 | D | 0 | 0 | 0 | 33 |
| QRN-QM-001 | C1 | v0.1 | CR | 0 | 2 | 1 | 0 |
| QRN-GR-001 | C1 | v0.1 | CR | 0 | 0 | 1 | 0 |
| QRN-STAT-001 | C1 | v0.2 | CR | 0 | 3 | 1 | 0 |
| QRN-FIELD-001 | C1 | v0.3 | CR | 0 | 2 | 1 | 0 |
| QRN-GAUGE-001 | C1 | v0.4 | CR | 0 | 2 | 1 | 0 |
| QRN-ENT-001 | C1 | v0.5 | CR | 0 | 2 | 2 | 5 |
| QRN-GEOM-001 | C1 | v0.6 | CR | 0 | 2 | 2 | 0 |
| QRN-GEOM-002 | C1 | v0.6 | CR | 0 | 2 | 2 | 0 |
| QRN-GEOM-003 | C3 | v0.7 | CR | 0 | 3 | 2 | 12 |
| QRN-GRAV-001 | C1 | v0.7 | CR | 1 | 3 | 2 | 3 |
| QRN-GRAV-002 | C1 | v0.7 | CR | 0 | 2 | 1 | 0 |
| QRN-BH-001 | C1 | v0.8 | CR | 1 | 1 | 2 | 0 |
| QRN-CAUSAL-001 | C1 | v1.1 | CR | 0 | 3 | 2 | 11 |
| QRN-CAUSAL-002 | C3 | v1.1 | CR | 1 | 3 | 1 | 8 |
| QRN-ER-001 | C3 | v1.2 | CR | 0 | 2 | 1 | 4 |
| QRN-QEC-001 | C1 | v1.3 | CR | 0 | 0 | 1 | 4 |
| QRN-GAUGE-002 | C1 | v1.4 | CR | 0 | 2 | 1 | 0 |
| QRN-COSMO-001 | C1 | v1.5 | CR | 0 | 1 | 1 | 2 |
| QRN-COSMO-002 | C4 | v1.5 | CR | 1 | 1 | 2 | 0 |
| QRN-COSMO-003 | C3 | v1.5 | CR | 1 | 1 | 2 | 0 |
| QRN-BORN-001 | C1 | v1.6 | CR | 0 | 3 | 1 | 3 |
| QRN-BORN-002 | C1 | v1.6 | CR | 0 | 1 | 1 | 0 |
| QRN-GRAV-003 | C1 | v2.1 | CR | 0 | 3 | 1 | 3 |
| QRN-KK-001 | C1 | v2.2 | CR | 0 | 2 | 1 | 0 |
| QRN-MATTER-001 | C1 | v2.3 | CR | 1 | 2 | 1 | 31 |
| QRN-MATTER-002 | C2 | v2.3 | CR | 1 | 3 | 2 | 0 |
| QRN-COSMO-004 | C1 | v2.4 | CR | 0 | 0 | 1 | 0 |
| QRN-EXP-001 | C4 | v2.5 | CR | 0 | 0 | 2 | 12 |
| QRN-GAUGE-003 | C2 | v3.1 | CR | 1 | 6 | 3 | 15 |
| QRN-GAUGE-004 | C1 | v3.1 | CR | 0 | 1 | 1 | 0 |
| QRN-GAUGE-005 | C4 | v3.1 | CR | 0 | 1 | 2 | 0 |
| QRN-YUK-001 | C4 | v3.2 | CR | 1 | 3 | 2 | 31 |
| QRN-RG-001 | C1 | v3.3 | CR | 0 | 3 | 1 | 3 |
| QRN-UNRUH-001 | C1 | v3.4 | CR | 1 | 1 | 1 | 0 |
| QRN-QNEC-001 | C1 | v4.1 | CR | 1 | 4 | 2 | 7 |
| QRN-NU-001 | C4 | v4.2 | CR | 1 | 3 | 2 | 0 |
| QRN-GAUGE-006 | C2 | v4.3 | CR | 1 | 7 | 3 | 6 |
| QRN-DS-001 | C1 | v4.4 | CR | 1 | 1 | 1 | 1 |
| QRN-DS-002 | C4 | v4.4 | CR | 1 | 1 | 1 | 0 |
| QRN-ARROW-001 | C3 | v5.1 | CR | 1 | 3 | 1 | 6 |
| QRN-GAUGE-007 | C2 | v5.2 | CR | 1 | 5 | 2 | 7 |
| QRN-EFT-001 | C1 | v5.3 | CR | 0 | 3 | 1 | 0 |
| QRN-GAUGE-008 | C2 | v6.2 | CR | 3 | 7 | 2 | 4 |
| QRN-GAUGE-009 | C2 | v6.2 | CR | 1 | 4 | 2 | 0 |
| QRN-GAUGE-010 | C2 | v6.2 | CR | 2 | 5 | 2 | 0 |
| QRN-QNEC-002 | C1 | v6.3 | CR | 2 | 4 | 3 | 3 |
| QRN-GEOM-004 | C3 | v6.4 | CR | 1 | 3 | 1 | 1 |
| QRN-YUK-002 | C4 | v6.5 | CR | 2 | 3 | 2 | 30 |
| QRN-EXP-002 | C5 | v6.6 | D | 1 | 0 | 1 | 0 |
| QRN-CORE-001 | C3 | v6.7 | CR | 2 | 3 | 2 | 4 |
| QRN-GAUGE-011 | C2 | v6.8 | LR | 1 | 6 | 1 | 1 |
| QRN-GAUGE-012 | C2 | v7.1 | CR | 0 | 5 | 2 | 1 |
| QRN-YUK-003 | C4 | v7.2 | CR | 3 | 7 | 2 | 28 |
| QRN-EXP-003 | C4 | v7.3 | CR | 1 | 1 | 2 | 10 |
| QRN-CORE-002 | C3 | v7.4 | CR | 3 | 3 | 1 | 3 |
| QRN-GAUGE-013 | C2 | v7.5 | LR | 1 | 6 | 1 | 1 |
| QRN-YUK-004 | C4 | v8.1 | CR | 1 | 7 | 2 | 26 |
| QRN-GAUGE-014 | C2 | v8.2 | CR | 0 | 5 | 2 | 3 |
| QRN-META-001 | C5 | v1.0 | D | 5 | 1 | 6 | 2 |
| QRN-META-002 | C5 | v2.0 | D | 3 | 1 | 1 | 2 |
| QRN-META-003 | C5 | v4.0 | D | 4 | 2 | 2 | 2 |
| QRN-META-004 | C5 | v5.0 | D | 4 | 1 | 3 | 1 |
| QRN-META-005 | C5 | v6.0 | D | 1 | 1 | 3 | 0 |
| QRN-META-006 | C5 | v7.0 | D | 5 | 0 | 1 | 0 |
| QRN-META-007 | C5 | v8.0 | D | 5 | 1 | 1 | 0 |
| QRN-META-008 | C5 | v9.0 | D | 2 | 1 | 1 | 0 |
| QRN-YUK-005 | C4 | v9.1 | CR | 1 | 7 | 2 | 24 |
| QRN-YUK-006 | C4 | v9.2 | CR | 3 | 6 | 1 | 23 |
| QRN-CORE-003 | C1 | v9.4 | CR | 1 | 4 | 1 | 1 |
| QRN-META-009 | C5 | v10.0 | DR | 3 | 0 | 1 | 0 |
| QRN-YUK-007 | C4 | v10.1 | CR | 1 | 7 | 2 | 21 |
| QRN-YUK-008 | C4 | v10.2 | CR | 1 | 7 | 2 | 1 |
| QRN-YUK-009 | C2 | v10.3 | CR | 1 | 4 | 2 | 4 |
| QRN-META-010 | C5 | v11.0 | DR | 3 | 0 | 1 | 0 |
| QRN-GAUGE-015 | C2 | v11.1 | CR | 2 | 2 | 2 | 1 |
| QRN-YUK-010 | C3 | v11.2 | CR | 1 | 7 | 1 | 2 |
| QRN-YUK-011 | C4 | v11.3 | CR | 1 | 8 | 2 | 1 |
| QRN-CORE-004 | C3 | v11.4 | CR | 2 | 2 | 2 | 1 |
| QRN-META-011 | C5 | v12.0 | DR | 4 | 0 | 1 | 0 |
| QRN-YUK-012 | C3 | v12.1 | CR | 1 | 2 | 1 | 6 |
| QRN-YUK-013 | C4 | v12.2 | CR | 2 | 5 | 1 | 5 |
| QRN-YUK-014 | C4 | v12.3 | CR | 1 | 5 | 1 | 4 |
| QRN-META-012 | C5 | v13.0 | DR | 3 | 0 | 1 | 0 |
| QRN-TOOL-001 | C1 | v13.1 | CR | 0 | 1 | 1 | 3 |
| QRN-YUK-015 | C4 | v13.2 | CR | 3 | 5 | 1 | 2 |
| QRN-META-013 | C5 | v14.0 | DR | 3 | 0 | 1 | 0 |
| QRN-GAUGE-016 | C2 | v14.4 | LR | 1 | 6 | 1 | 1 |
| QRN-GAUGE-017 | C2 | v14.5 | LR | 1 | 6 | 1 | 1 |
| QRN-META-014 | C5 | v15.0 | DR | 2 | 0 | 1 | 0 |
| QRN-META-015 | C2 | v15.1 | CDR | 0 | 1 | 1 | 2 |
| QRN-TOOL-002 | C2 | v15.2 | CDR | 1 | 1 | 1 | 1 |
| QRN-CORE-005 | C3 | v15.3 | CDR | 1 | 3 | 2 | 1 |
| QRN-CONT-001 | C1 | v15.4 | CDR | 1 | 4 | 2 | 1 |
| QRN-SEL-001 | C4 | v15.5 | CDR | 1 | 7 | 1 | 1 |
| QRN-GRAV-004 | C1 | v15.6 | CDR | 1 | 4 | 2 | 1 |
| QRN-PRED-001 | C4 | v15.7 | CDR | 2 | 8 | 3 | 8 |
| QRN-META-016 | C5 | v16.0 | DR | 7 | 0 | 1 | 0 |
| QRN-YUK-016 | C4 | v16.1 | CDR | 1 | 5 | 1 | 0 |
| QRN-YUK-017 | C3 | v16.2 | CDR | 2 | 4 | 1 | 6 |
| QRN-YUK-018 | C4 | v16.3 | CDR | 2 | 7 | 2 | 5 |
| QRN-YUK-019 | C4 | v16.4 | CDR | 1 | 7 | 2 | 4 |
| QRN-YUK-020 | C4 | v16.5 | CDR | 1 | 7 | 2 | 0 |
| QRN-YUK-021 | C4 | v16.6 | CDR | 1 | 7 | 1 | 0 |
| QRN-YUK-022 | C4 | v16.7 | CDR | 1 | 8 | 2 | 1 |
| QRN-YUK-023 | C4 | v16.8 | CDR | 2 | 8 | 2 | 0 |

## 仮定の影響範囲 — これを抜くと何が落ちるか

「直接」= この仮定を明示的に使う主張の数。「閉包」= 依存を遡って落ちる主張の総数。

| 仮定 | type | status | 直接 | 閉包 | 閉包に含まれる主張 (抜粋) |
|---|---|---|---|---|---|
| ASM-LATTICE | framework | active | 26 | 65 | QRN-ARROW-001, QRN-BORN-001, QRN-CAUSAL-001, QRN-CAUSAL-002, … |
| ASM-SEED | design | active | 13 | 46 | QRN-BH-001, QRN-BORN-001, QRN-COSMO-003, QRN-FIELD-001, … |
| ASM-PDG | data | active | 27 | 38 | QRN-COSMO-001, QRN-COSMO-002, QRN-COSMO-003, QRN-DS-002, … |
| ASM-LOWDIM | framework | active | 23 | 34 | QRN-ARROW-001, QRN-BORN-001, QRN-CAUSAL-001, QRN-CAUSAL-002, … |
| ASM-TORUS | model | active | 25 | 33 | QRN-KK-001, QRN-MATTER-001, QRN-META-007, QRN-META-008, … |
| ASM-PRIOR | design | active | 23 | 33 | QRN-META-006, QRN-META-007, QRN-META-008, QRN-META-009, … |
| ASM-GAUSS | framework | active | 20 | 31 | QRN-ARROW-001, QRN-CAUSAL-001, QRN-CAUSAL-002, QRN-CONT-001, … |
| ASM-OVERLAP | model | active | 22 | 30 | QRN-META-007, QRN-META-008, QRN-META-009, QRN-META-010, … |
| ASM-WILSON-GRID | design | active | 22 | 30 | QRN-META-007, QRN-META-008, QRN-META-009, QRN-META-010, … |
| ASM-DIAGPAIR | model | falsified | 4 | 29 | QRN-META-007, QRN-META-008, QRN-META-009, QRN-META-010, … |
| ASM-STABLE-LABEL | convention | active | 16 | 29 | QRN-META-007, QRN-META-008, QRN-META-009, QRN-META-010, … |
| ASM-GAUGE-GROUP | model | active | 14 | 23 | QRN-GAUGE-003, QRN-GAUGE-004, QRN-GAUGE-006, QRN-GAUGE-007, … |
| ASM-SIGMA-DATA | model | active | 11 | 22 | QRN-META-010, QRN-META-011, QRN-META-012, QRN-META-013, … |
| ASM-ANOMALY-COEFS | data | active | 14 | 22 | QRN-GAUGE-003, QRN-GAUGE-006, QRN-GAUGE-007, QRN-GAUGE-008, … |
| ASM-CHIRALITY | model | active | 11 | 21 | QRN-GAUGE-003, QRN-GAUGE-006, QRN-GAUGE-007, QRN-GAUGE-008, … |
| ASM-ALL-CHARGED | definition | active | 11 | 21 | QRN-GAUGE-003, QRN-GAUGE-006, QRN-GAUGE-007, QRN-GAUGE-008, … |
| ASM-EFT-VALIDITY | model | active | 3 | 16 | QRN-GAUGE-003, QRN-GAUGE-006, QRN-GAUGE-007, QRN-GAUGE-008, … |
| ASM-WINDOW-V31 | window | active | 3 | 16 | QRN-GAUGE-003, QRN-GAUGE-006, QRN-GAUGE-007, QRN-GAUGE-008, … |
| ASM-MODK | definition | active | 7 | 14 | QRN-CORE-003, QRN-GRAV-001, QRN-GRAV-002, QRN-GRAV-003, … |
| ASM-KTM | model | active | 2 | 11 | QRN-EXP-003, QRN-META-007, QRN-META-016, QRN-PRED-001, … |
| ASM-INIT | model | active | 1 | 9 | QRN-ARROW-001, QRN-CAUSAL-002, QRN-CORE-002, QRN-CORE-004, … |
| ASM-WINDOW-EXT | window | active | 4 | 8 | QRN-GAUGE-007, QRN-GAUGE-008, QRN-GAUGE-010, QRN-GAUGE-013, … |
| ASM-DOF-GROWTH | model | active | 2 | 7 | QRN-ARROW-001, QRN-CORE-002, QRN-CORE-004, QRN-META-004, … |
| ASM-WINDOW-V43 | window | active | 4 | 7 | QRN-GAUGE-006, QRN-GAUGE-008, QRN-GAUGE-010, QRN-GAUGE-016, … |
| ASM-OBS-FRACTIONAL | observational | active | 1 | 7 | QRN-GAUGE-006, QRN-GAUGE-008, QRN-GAUGE-010, QRN-GAUGE-016, … |
| ASM-LEAN-TRUST | trust | active | 4 | 7 | QRN-GAUGE-011, QRN-GAUGE-013, QRN-GAUGE-016, QRN-GAUGE-017, … |
| ASM-NET-REAL | ontology | active | 7 | 7 | QRN-META-001, QRN-META-002, QRN-META-003, QRN-META-004, … |
| ASM-WINDOW-PAIR | window | active | 1 | 5 | QRN-META-010, QRN-META-013, QRN-YUK-009, QRN-YUK-015, … |
| ASM-WINDOW-U1CUBE | window | active | 1 | 4 | QRN-GAUGE-014, QRN-GAUGE-015, QRN-META-008, QRN-META-011 |
| ASM-ORBIFOLD | model | active | 2 | 3 | QRN-META-011, QRN-YUK-010, QRN-YUK-011 |
| ASM-WINDOW-U1SQ | window | active | 2 | 3 | QRN-GAUGE-009, QRN-GAUGE-012, QRN-META-007 |
| ASM-Z2-MINIMAL | model | active | 2 | 3 | QRN-CONT-001, QRN-CORE-005, QRN-META-016 |
| ASM-EDGE-SEMANTICS | design | active | 2 | 3 | QRN-META-015, QRN-META-016, QRN-TOOL-002 |
| ASM-SMCONTENT | data | active | 2 | 2 | QRN-GAUGE-009, QRN-MATTER-002 |
| ASM-WINDOW-EXC | window | active | 1 | 2 | QRN-GAUGE-015, QRN-META-011 |
| ASM-WICK | framework | active | 1 | 1 | QRN-STAT-001 |
| ASM-ENVARIANCE | framework | active | 1 | 1 | QRN-BORN-002 |
| ASM-QM | framework | active | 0 | 0 |  |

## 反証条件の射程 — これが発火すると何が落ちるか

| 反証条件 | status | 直接 | 閉包 | 条件 (要約) |
|---|---|---|---|---|
| FAL-SUITE | open | 95 | 105 | 再現スイートの回帰 FAIL — コードの再実行が主張の数値を再現しない (乱数は固定シード)。 |
| FAL-CKM-OOS | open | 12 | 32 | out-of-sample の CKM/PMNS 予測が系統的に失敗する (holdout 量が帯を外れる)。 |
| FAL-CEX-WINDOW | open | 15 | 27 | 明示された探索窓の内側で反例 (条件を満たす非 SM 解、または SM が条件を満たさないこと) が提示される。 |
| FAL-CONTINUUM | open | 9 | 23 | readout の残差 residual(N)=A·N^(−p)+B の fit で B≠0 が確立する (格子を細かくしても消えない系統残差)。 |
| FAL-BMV | open | 7 | 16 | BMV 型実験で、QRN の予測位相 Δφ=Gm²τΔx²/(ħd³) にて C≡0 かつ (C,V) が古典包絡 (C=0, V≤e^(−Δφ/2)) の内側に留まる。 |
| FAL-EXOTIC-CHIRAL | open | 4 | 16 | SM 1 世代の外のカイラル物質 (第 4 世代・分数電荷ハドロン・エキゾチック表現) が発見される。 |
| FAL-AREALAW | open | 4 | 8 | アナログ系 (冷却原子等) で面積則・エンタングルメント第一法則の系統的破れが観測される。 |
| FAL-QNEC | open | 5 | 8 | QNEC 型不等式の真の破れ (誤差予算を超える負ギャップが N を上げても縮まない) が数値または理論で確立する。 |
| FAL-PAGE | open | 4 | 5 | ブラックホール蒸発で情報喪失が確立する (Page 曲線からの系統的逸脱)。 |
| FAL-LORENTZ | open | 2 | 4 | プランクスケールのローレンツ不変性の破れが主要次数で確立する (GRB 光子のエネルギー依存遅延など)。 |
| FAL-DIMFLOW | open | 2 | 4 | 超高エネルギーで実効スペクトル次元が 2 へ低下しないことが確立する。 |
| FAL-GLOBALSYM | open | 1 | 3 | 厳密な大域的対称性の存在が確立する (例: 陽子の厳密な安定性)。 |
| FAL-COSMO | open | 2 | 2 | w(z) の精密観測が常在 Λ 型の揺らぎを排除する。 |
| FAL-SUSY | open | 1 | 1 | 超対称粒子の排除限界が MSSM 統一に必要なスケールを完全に排除する。 |
| FAL-NEUTRINO | open | 1 | 1 | ニュートリノが Dirac と確定する (0νββ の排除)、または質量順序・δ_CP がシーソー+FN の帯を外れる。 |

