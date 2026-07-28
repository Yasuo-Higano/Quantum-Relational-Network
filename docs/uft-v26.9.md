# QRN v26.9 — 4D covariance closure (Gate 5)

**位置づけ**: 第二十七期、誘導重力 経路 B の最終監査段階。事前登録は spec
§12.6 (bc644d4 + 08a1321): v26.9 = **4D covariance closure** — h₀₀・h₀ᵢ・
q₀ ≠ 0・10×10 symmetric-tensor kernel・4D Ward q_μΠ^{μν,ρσ} = 0・
contact/tadpole 込み second variation・spin-0/2 完全分離。**Gate 5 を通る
まで 1/Π・graviton propagator・dynamic metric へ進まない** (dynamic metric
fork は v27.0)。発行可能な型名は Gate 5 まで従来 4 種のみ。

前提 (v26.8 で完結済み): 空間セクターの 4 比 universality (PRED-016
scored-hit, max 0.14%)・スカラー和則 (PRED-017 scored-hit)・三者一致
(spectral density の a² Lorentz 回復)。**全て「測定器が正しい」ことの証明で
あり、QRN・創発重力の証拠ではない** (spec §12.8 の凍結解釈)。

## 0. v269e_energy_ward — エネルギーセクターの 4D Ward (厳密演算子恒等式)

h₀₀ セクターの定義閉包と、プロジェクト初の **q₀ ≠ 0 (動的) Ward 恒等式**。
7 検査 PASS。

**導出 (本ユニットの核)**: 中点変調規約 (e^{iq(y+1/2)} on y-bond / e^{iqy}
on-site・横ボンド) の格子エネルギー密度 T₀₀ は折込み 8 成分基底で

  **V₀₀(k;q) = h(k + (q/2)ŷ) — 厳密**

(y ボンドの両向き和が (−1)^{s_y}cos(k_y + q/2)、x/z ボンド・m 項は q 非依存。
η_z の (−1)^y は s_y flip として eps に吸収)。**副発見**: v268z の一般頂点
公式の (s+ε) 位相規約は **q 方向のボンドには適用できない** (両片が相殺する
誤った 0 を返す — v26.8 の D/X は全て q ⊥ ボンドだったため露見しなかった)。

**連続の式 (J_E の構成的存在証明)**: C(k;q) := h(k+qŷ)V₀₀ − V₀₀h(k) に対し
V_J(k;q) := i·C/(2 sin(q/2)) が **Fourier 台 |n| ≤ 2 の三角多項式** (残差
8.0e-14) — 到達距離 ≤ 2 サイトの局所エネルギー流が厳密に存在する。q → 0
極限存在 (O(q) 収束 2.4e-4)。

**動的 Ward (f-和則型, Euclidean q₀)**: χ_A(iq₀,q) = Σ 2|M_A|²ΔE/(ΔE²+q₀²)
に対し、対ごとの厳密恒等式 (E_μ−E_ν)M₀₀ = −i·q̂·M_J (1.6e-14) の帰結として

  **q₀²·χ₀₀(iq₀,q) + q̂²·χ_JJ(iq₀,q) = M₁(q)**   (q̂ = 2sin(q/2))

が BZ 積分で **4.4e-16** — カットオフ有限のまま周波数方向の Ward が厳密に
成立。M₁ = Σ2|M₀₀|²ΔE は等時刻交換子 (f-和則) 値。

| 検査 | 結果 |
|---|---|
| E0 V₀₀(k;0) = h(k)・転送ブロック整合 | 0.0 (厳密) |
| E1 taste-singlet (指標純度) | 1.8e-17 |
| E2 tree matching ‖V₀₀ − h_lin‖/‖h_lin‖ = O(ε²) | 縮小比 3.99, 4.00 |
| E3 J_E 局所性 (Fourier 台 ≤ 2) | 8.0e-14 |
| E4 対ごと恒等式 ΔE·M₀₀ = −i q̂ M_J | 1.6e-14 |
| E5 変異 (m 重み ×1.01 → 可分性破れ) | 発散比 24.5 (正版 0.81) |
| E6 動的 Ward (q = 0.4, q₀ ∈ {0.3, 0.9}) | ≤ 4.4e-16 |

開発記録: (i) E2 は run1-2 で固有基底の単一行列要素比を使い、±E 各 4 重縮退
の部分空間任意回転で非単調 — Frobenius (基底不変量) に置換。教訓:
「**縮退があるときは基底不変量で測れ**」。(ii) E3 の q → 0 収束判定は
V_J の O(q) 依存を O(q²) と誤想定 — 倍化 q 差分に修正 (物理は無変更)。

## A. v269m_momentum_ward — 運動量セクターの 4D Ward と Belinfante 対称性

h₀ᵢ の source T₀ᵢ の定義閉包 (6 検査 PASS)。**構成**: taste 安全な 2 サイト
分割 (1 サイトシフトは taste 混合・2 サイトは η² = 1 で自明 singlet):
**V₀y(k;q) = −(1/2)sin(2k_y + q)·𝟙** — node 極限 = 中点運動量 p_y + q/2、
q = 0 で [h, V₀y] = 0 厳密 (全運動量保存)。

**格子連続の式の符号解析**: ∂_tρ(q) = i[H,ρ(q)] = +2i sin(q/2)·J(q)
(中点 Fourier 規約) ⇒ **V_S = C/(2sin(q/2)) は実行列** (i は付かない)。

| 検査 | 結果 |
|---|---|
| M0 恒等式 ([h,V₀y(0)] = 0 含む) | 5.6e-17 |
| M1 V_S の局所性 (Fourier 台 ≤ 3) | 1.2e-15 |
| **M2 保存 stress = BOND-A T_yy (Z なし)** | rel: 6.0e-2 → 9.0e-4 (縮小比 4.1, 4.0, 4.0 = O(ε²)) |
| **M3 Belinfante T₀y = J_E^y (on-shell)** | rel: 1.3e-1 → 1.9e-3 (縮小比 4.1, 4.0, 4.0 = O(ε²)) |
| M4 動的 Ward q₀²χ₀y + q̂²χ_SS = M₁ | ≤ 1.7e-15 |
| M5 変異 (taste-nonsinglet 密度) | 発散比 19.4 (正版 0.66) |

**Gate 5 の主要成分が 2 つ確立**: (i) **BOND-A (metric 変分) の stress =
保存流由来の stress が tree レベル・自由パラメータなし (Z = 1) で一致** —
「BOND-A source は保存 T_μν に流れる」の momentum-行版。(ii) **Belinfante
対称性 T₀y = T_y0 (= J_E) が on-shell O(ε²) で成立** — 4D 対称テンソル
kernel の対称性の前提が正当化された。

開発記録: (i) run1 は V_S := iC/q̂ と置き、実の BOND-A 頂点との射影ブロック差
が **rel = √2 で停留 — 位相直交の指紋** (i は規約の取り違え。符号解析で
J = C/(2sin(q/2)) 実行列と確定)。教訓: 「**√2 の停留は大きさでなく位相の
不一致を疑え**」。(ii) run1 の変異 (h.c. 片 ×1.02) は複素スカラー×𝟙 のまま
h と可換で不発 — s_y 依存重み (taste-nonsinglet) に再設計。

## B. v269w_ward4d — full 4D Ward: 64 恒等式の一括機械検査

energy/momentum 行の器械を全 Ward 行に一般化 (6 検査 PASS)。**一様構成**:
密度 4 種 D_ν (V₀₀ = h(k+q/2ŷ), V₀ᵢ = −(1/2)sin(2kᵢ+qδ_iy)𝟙) に対し保存
フラックス **V_Fν := [h(k+qŷ)D_ν − D_νh(k)]/(2sin(q/2))** が自動生成され、
全 4 本が局所 (Fourier 台 ≤ 3, 6.3e-14)。列 = 密度 4 + BOND-A stress 4
(V_xx, V_yy [中点変調], V_zz, V_xz [point-split])。

**Matsubara Ward**: iq₀·C_{AB}(iq₀,q) − q̂·C_{J_A B}(iq₀,q) = −⟨[A(q),B(−q)]⟩
が対ごとの分数分解で厳密 — **4 行 × 8 列 × q₀ ∈ {0.3, 0.9} = 64 恒等式が
k 点ごと 1.5e-12・12³ BZ 積分 4.8e-14 で閉じる**。接触項は独立 2 実装
(対和 = −⟨[A,B]⟩ / 占有トレース = +⟨[A,B]⟩) が 1.3e-13 で照合。

**Ward の物理の所在 (本ユニットの概念的成果)**: J_A := [h,A]/q̂ と定義した
時点で Matsubara 恒等式は任意の双線形に対する厳密再配列 — **物理的内容は
(i) フラックスの局所性 (保存が破れると q → 0 で 1/q̂ 発散 = W4 変異が実演,
比 13.6) と (ii) 接触項の構造に宿る** (連続の Ward = 保存則 + カレント局所性
と同じ分業)。副定理: 接触項の対和とトレース公式は k 点ごとには一致しない —
差の occ-occ 片は k → k+q ラベル替えで BZ 和のみ相殺 (格子上で厳密にするには
**q を格子と可約に取る** — q = π/6 = 12³ 格子の 2 刻み)。

さらに [W3]: ‖V_Fy − V_yy^A‖_F (off-shell **全行列** Frobenius) が実効指数
**ε^2.0** — BOND-A と保存 stress の差は off-shell でも O(ε²) (v26.9-A の
on-shell 判定より強い)。

開発記録: run1 は (i) 接触項 2 実装を k 点ごとに比較 (occ-occ 片で不一致 —
上の副定理で解決)、(ii) 変異ゲートを Ward 残差に置いた (恒等式は任意 A で
成立するため不発 — 局所性破れ = 発散比に再設計)、(iii) 対和とトレースの
符号 (差/和 = 2.0 の指紋で確定)。

## C. v269c_spin4d — 4D spin 分離と Gate 5 総括

v26.9 arc の最終ユニット (6 検査 PASS)。**置換則** (v26.9-0/A の 3 独立導出と
整合、d_y = 1, 2 で明示照合 1.7e-16): 中点変調頂点 V(k;q) = V_unmod(k+q/2ŷ)。
これで T_xy/T_yz point-split が v268z の X 構成の軸回転で得られ、**全 10
source が完成**。full 4D Ward は 4 行 × 10 列 × 2 周波数 = **80 恒等式
(4.3e-13)**。

**Belinfante 構造の格子分解 (発見)**:
- **厳密恒等式: 正準フラックス V_Fx = point-split の piece2 (Γy 構造 × x 包絡)
  単独** — 機械精度 3.8e-16 の格子演算子恒等。[h, V₀x]/q̂ は「片脚」。
- Belinfante (両 piece 平均) とは回転流 (spin current) の分だけ O(1) (0.81) で
  違う — 差は保存する改良項で、これが格子上で目に見える形で分離した。
- scheme 混合の定理 (run1): (正準 T⁰ᵢ, Belinfante T_ij) を混ぜた 10×10 の
  横断性破れは **a 非依存 (0.616)** — 横断性は一貫 scheme でのみ。格子の厳密
  横断性は**正準 (非対称 16 成分) テンソル**が担う (ΔE·M₀ν = q̂·M_Fν, 殻上
  1e-10)。

**4D spin-2 分離**: D = (T_xx−T_zz)/√2 と X = T_xz は timelike q_L = (E, Qŷ)
でも**厳密 P₂** (q·D = 0, tr_θD = 0 が (t,y) 面の任意 q で成立 — ProjectorND
の ŷ 定理の時間方向拡張)。殻積分 (E, Q) = (1.5, 0.6), a ∈ {0.18, 0.09, 0.045}:

| 検査 | 値 |
|---|---|
| 偏極縮退 \|σ_DD/2σ_XX − 1\| (⟨D\|P₂\|D⟩ = 1, ⟨X\|P₂\|X⟩ = ½) | 0.0118 → 0.0029 → **0.0007** (縮小比 4.0, 4.0 = O(a²)) |
| 直交性 \|σ_DX\|/√(σσ) | **1.8e-16** |
| oracle 絶対アンカー σ_DD/(2ρ_D(E²−Q²)) | **1.0007** |
| 変異 (X の Z 補正落とし) | 2σ_XX/σ_DD = 3.997 (= 4 予言どおり) |

開発記録: (i) run1 の混合 scheme 10×10 (上の定理に転化)。(ii) run2 の対角
BOND-A 重み半分ミス — σ_XX = 2σ_DD の因子 4 異常で発覚 (「規約係数は認証済み
バイナリの q = 0 極限と突き合わせよ」)。(iii) 殻積分の (2π)³ 正規化落ち —
アンカー 248.2 ≈ (2π)³ = 248.05 の指紋で即特定。

### Gate 5 総括 (凍結解釈)

**確立**: h₀₀/h₀ᵢ source (厳密構成)・q₀ ≠ 0 (Matsubara Ward 64+80 恒等式)・
全 10 source (置換則)・full 4D Ward = 局所カレント + 計算可能な接触項・
正準テンソルの厳密横断性・spin-2 の 4D 分離 (P₂ 縮退 + 直交 + oracle
アンカー)・BOND-A = 保存 stress (tree, O(ε²))・Belinfante 構造の格子分解。

**残り**: (1) Belinfante 対称 10×10 の完全崩壊 — T⁰ⁱ の対称化 (= エネルギー
流との平均) が必要で、x/z エネルギー流は**二重変調則
V₀₀(k; px̂+qŷ) = h(k+(p/2)x̂+(q/2)ŷ)** で構成可能 → **v26.9-D として登録**。
(2) temporal h の二次変分 scheme (kinetic term 部門) → v27.0 設計項目。
**型名 FullGravitationalVacuumPolarization は保留を維持・1/Π 禁止も維持**
(凍結解釈: すべて測定器の証明であり QRN・創発重力の証拠ではない)。

## 1. 残り (v26.9 arc)

- **v26.9-D**: Belinfante 対称 10×10 の完全崩壊 (二重変調則で T⁰ⁱ を対称化 —
  σ = ρ₂P₂ + ρ₀P₀s への完全崩壊と ρ₀ の oracle 照合)。
- その後 Gate 5 の最終判定 → v27.0 dynamic metric fork の分岐判断。

## 2. 成果物

0: `sim/src/bin/v269e_energy_ward.rs` / `results/v269e_energy_ward.txt`
(7 検査 PASS) / `results/v269e_energy_ward.json`。claims: QRN-GRAV-051。
A: `sim/src/bin/v269m_momentum_ward.rs` / `results/v269m_momentum_ward.txt`
(6 検査 PASS) / `results/v269m_momentum_ward.json`。claims: QRN-GRAV-052。
B: `sim/src/bin/v269w_ward4d.rs` / `results/v269w_ward4d.txt`
(6 検査 PASS) / `results/v269w_ward4d.json`。claims: QRN-GRAV-053。
C: `sim/src/bin/v269c_spin4d.rs` / `results/v269c_spin4d.txt`
(6 検査 PASS) / `results/v269c_spin4d.json`。claims: QRN-GRAV-054。
