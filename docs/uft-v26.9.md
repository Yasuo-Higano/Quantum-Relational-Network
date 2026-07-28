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

## 1. 残り (v26.9 arc)

- **運動量セクター**: h₀ᵢ ↔ T₀ᵢ の格子構成と、T₀ᵢ = Tᵢ₀ = J_E (Belinfante
  対称性) の格子検証。運動量保存は離散並進 — Ward の形が energy 行と異なる。
- 10×10 symmetric-tensor kernel (q₀ ≠ 0) の組み立てと 4D 射影子分解
  (ProjectorND.lean の d = 4 代数を使う)。
- 4D Ward q_μΠ^{μν,ρσ} = contact 項 (tadpole/seagull 第二変分) の一括検証
  — **Gate 5 の判定**。
- spin-0/spin-2 form factor の 4D 分離と連続極限。

## 2. 成果物

0: `sim/src/bin/v269e_energy_ward.rs` / `results/v269e_energy_ward.txt`
(7 検査 PASS) / `results/v269e_energy_ward.json`。claims: QRN-GRAV-051。
