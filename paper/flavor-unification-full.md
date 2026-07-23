# Falsification-Driven Development of a Magnetized-Torus Flavor Model: A Prequential Audit of Adaptive Model Selection

**Draft v1 (v26.1).** 対象誌: JHEP / PRD (方法論成分が強い場合は SciPost)。本稿は三部作 [geometric-yukawa-full.md](geometric-yukawa-full.md) / [cp-complex-structure-full.md](cp-complex-structure-full.md) / [measure-dissolution-full.md](measure-dissolution-full.md) を **companion に再位置づけ**し、プログラム横断の予測力の監査を主結果とする。規則・台帳・数値の一次ソース: [flavor-unification-plan.md](flavor-unification-plan.md) (凍結 aa66a3e) / [prequential_ledger.yml](prequential_ledger.yml) (凍結 c1121a3) / `results/v261_prequential.txt`。**三部作の個別投稿は本稿の完成をもって解除されるが、本稿と切り離した「三つの独立した成功」としての提示は行わない。**

---

## Abstract

Three companion papers developed a magnetized-torus flavor model through explicit falsification: a rectangular T²×T² construction with no order-one coefficients (masses and CKM magnitudes), its refutation by a held-out observable (the Jarlskog invariant vanishes structurally), repair by complex structure (shear), and the dissolution of the residual measure problem into geometry (τ = 1/12 + i/2). Each paper paid its own within-stage holdout costs; none paid the *program-level* cost of adaptive model selection — each stage was chosen using the failures of the previous one. Here we pay that cost in full. We freeze a prequential protocol and a machine-readable chronology ledger (which observables first influenced selection at which commit), then score S = Σ_k ln p(D_k | M_{k−1}) using the era evidence machinery verbatim (eight instrument regressions reproduce the archived lnZ values, MAP configurations, and the frozen |V_td| prediction band bit-for-bit at print precision). The result is negative and we report it as the main result: against the program's own 2013-style Froggatt–Nielsen baseline (literature charges, complex O(1) coefficients marginalized), the geometric program's fresh-data score is S_prog − S_FN = −55.9 (bracketing the strict value −∞ set by the structural CP zero); against anarchy, −41.8. The deficit is dominated by the CP block — for twelve versions the program explored a geometry family that could not make CP while both baselines could (FN nearly centers |J|: ln-density +8.53 of a maximum +9.83) — and even excluding CP the program has not yet out-predicted FN (−1.96 over |V_td|, |V_ts|, β, γ). What survives is delimited precisely: the within-stage evidence contests remain valid as model comparison on contemporaneous data; the random-coefficient-free construction remains qualitatively distinct (FN marginalizes 18 complex O(1) coefficients); the converged point τ = 1/12 + i/2 has not yet faced fresh data and is the program's registered forward prediction; and the prediction-registration machinery that made this audit decidable (a frozen |V_td| band scored as a miss, a frozen J = 0 scored as window-refuted, both *before* the repair) is, we argue, the transferable result.

## 1. Introduction: three papers, one adaptive process

The companion papers established, in order: (i) fermion mass hierarchies and CKM magnitudes from magnetized tori with no random coefficients, with the discrete data (fluxes, Wilson lines, generation pairing) shown irreducible; (ii) that CP violation requires complex structure — the rectangular geometry has a *structural* Jarlskog zero, discovered when a pre-registered holdout prediction was scored; (iii) that the remaining "measure problem" dissolves into geometry at an interior complex-structure point where all twelve observables sit within a factor of five. Read separately, each paper reports a success with its own holdouts. Read as a sequence, each stage's model class was selected *using the observed failures of the previous stage*: the shear family exists because J refuted the rectangular family; the τ_im scan exists because measure corrections kept compensating a geometric lever. Within-stage holdouts cannot price this adaptivity. The prequential (predictive-sequential) score does: every datum is scored exactly once, against the model frozen *before* that datum influenced any choice.

Two features make this audit decidable rather than rhetorical. First, the program's own prediction ledger (v15.7) had frozen and honestly scored two predictions before the repair: |V_td| (band [0.0002, 0.0070] at 95%; measurement 0.0086 — a recorded **miss**) and J (exactly 0; measurement 3.08×10⁻⁵ — recorded **window-refuted**). These fix the "first seen" boundary. Second, the era evidence machinery is deterministic and archived, so the frozen models' posterior predictive densities can be recomputed today and *regressed against the archived numbers* before any new quantity is trusted.

## 2. The chronology ledger

The ledger (machine-readable, frozen before computation) reconstructs from the repository history: (a) for each observable, the version at which it first influenced selection; (b) for each model element, the observables used to justify its introduction; (c) the discrete alternatives actually tried at each stage (the search cost). Key entries: the CKM triple used from the start is {|V_us|, |V_cb|, |V_ub|} (seen at v3.2, in the likelihood from v9.1); |V_td| and J were first *scored* at v15.7 (the frozen predictions above), so the shear model's later |V_td| "holdout hit" is a reproduction of already-seen data and receives no prequential credit — the conservative rule applied at its most consequential point. The fresh blocks after the ledger's cut are |V_ts| and the unitarity-triangle angles (β, γ), first reported at v16.5, crediting the shear (1,1) model frozen at v16.4. A registered caveat: the nine-quantity evidence winner (three tori) was never extended to the shear family — untried alternatives are not charged (undercounting rule), but the discontinuity is disclosed.

## 3. Protocol

S = Σ_k ln p(D_k | M_{k−1}, D_{<k}) with p the posterior predictive density of the era-frozen model: era posterior weights (the archived likelihood machinery, verbatim) mixed with the era observation kernel (lognormal, σ = ln 2, the program's uniform convention; kernel-width sensitivity ln 1.5 / ln 3 reported throughout — all *differences* below are sign-stable across it). The structural CP zero is scored two ways, both published: strictly (ln p = −∞, the honest reading registered already at v16.4) and as an upper bracket (the best member of the 23-element tried set — tilted T⁴, window max |J| ≤ 2.6×10⁻⁸ — under a uniform mixture). Baselines: anarchy (complex log-uniform O(1) entries) and the v3.2 Froggatt–Nielsen model (literature charges, ε = 0.22), both by deterministic Monte Carlo (10⁶ draws; half-sample agreement ≤ 0.05 nats). The capital block D1 (mass ratios and the original CKM triple) is common training data by convention, Δ ≡ 0.

## 4. Instrument regressions (all pass)

The recomputed engines reproduce: rect lnZ₉ = −19.863 (archived −19.86), the archived MAP (lnL = −5.874, σ_H = 1, permutation pair (23)/(132)), the machine-zero MAP |J| = 8.4×10⁻¹⁷, and the frozen |V_td| band at all five archived quantiles; shear (1,1) lnZ₁₀ = −24.293 (archived), the archived MAP configuration (σ_H = 1, Wilson (7,32,9), σ = (5,1)) and its |V_td|, |V_ts|, β, γ values to print precision. A development record is kept: the first run's gate on the 2.5% quantile mis-fired against a one-significant-figure archived anchor (tolerance recalibrated to the anchor's rounding granularity), and the Monte Carlo half-sample gate at 2×10⁵ draws missed by 0.011 nats (draws increased to 10⁶); the physics numbers were identical across both runs.

## 5. Results

Posterior predictive ln-densities (kernel ln 2):

| block | rect (M1) | shear (1,1) | anarchy | FN (v3.2) |
|---|---|---|---|---|
| \|J\| | **−∞** (floor bracket −57.3) | +9.29 | +3.44 | **+8.53** |
| \|V_td\| | +2.73 | +3.92 | −2.40 | +3.42 |
| \|V_ts\| | +2.18 | +2.41 | −0.88 | +2.13 |
| β | −10.02 | −4.48 | −4.56 | −4.49 |
| γ | −11.15 | −6.85 | −5.38 | −5.29 |

Credited prequential increments (program = the ledger's credit targets):

| block | Δ(prog − anarchy) | Δ(prog − FN) |
|---|---|---|
| D_J (bracket: mixture −45.4) | −48.87 | −53.95 |
| D_VTD (rect) | +5.14 | −0.68 |
| D_VTS (shear) | +3.29 | +0.28 |
| D_β (shear) | +0.08 | +0.00 |
| D_γ (shear) | −1.47 | −1.56 |
| **total** | **−41.83 / strict −∞** | **−55.91 / strict −∞** |

Three readings. (1) The deficit is dominated by CP: for twelve versions the program searched a family that could not make CP while its own FN baseline nearly centered |J| — Wolfenstein scaling comes free with complex O(1) coefficients and hierarchical charges, and real-Wilson-line geometry started structurally behind. (2) Even excluding CP, the program has not yet out-predicted FN on fresh data (−1.96 over four blocks): it narrowly wins |V_ts| (+0.28), ties β, loses γ (−1.56) and |V_td| (−0.68). Against anarchy the non-CP blocks total +7.04 — the ability to make hierarchies at all is real predictive content. (3) The within-stage results of the companions (the +306-nat reversal once J is admitted to the likelihood; the 12/12 scorecard at the converged point) are *not* invalidated: they are model comparison on contemporaneous data. What this audit prices is prediction, and on that ledger the program is in deficit — the total of the costs the individual holdouts did not pay.

## 6. What survives, precisely

(i) The converged point τ = 1/12 + i/2 has faced no fresh data; it stands as the program's registered forward prediction, to be scored on the next independent observable updates. (ii) The construction remains coefficient-free: FN's density over |J| is bought by marginalizing 18 complex O(1) coefficients, the geometry's by discrete data alone — a qualitative difference the prequential score does not capture and which we do not convert into nats. (iii) The registration machinery itself: the audit was decidable only because a frozen band had been scored as a miss and a frozen zero as a refutation *at the time*. We propose this as the transferable practice: model-development programs should maintain exactly such a ledger, and report their prequential totals with their headline fits.

## 7. Limitations

The kernel convention (lognormal σ = ln 2 on all observables, including angles in degrees) is the program's uniform-rule extension of its era likelihood; absolute densities move with the kernel, differences are sign-stable across ln 1.5–ln 3. The D_J bracket uses the best tried member under a uniform mixture, not member-wise posterior integrals (the strict value is −∞). The FN baseline's literature charges embed that literature's own knowledge of the data — FN's hidden selection cost is unpriced here, so the anarchy column bounds the program's deficit from the generous side. The conservative first-seen rule can be relaxed (crediting the shear model's |V_td| at +3.92); the totals' signs do not change.

## 8. Reproducibility

Everything runs from the dependency-free Rust repository: `v261_prequential` (era engines transcribed verbatim, eight regressions gating validity), `paper/prequential_ledger.yml` (frozen chronology with evidence pointers), `paper/flavor-unification-plan.md` (frozen protocol), `results/v261_prequential.{txt,json}`. Claims and limitations are machine-tracked (QRN-YUK-033).

## References

1–3. The companion papers (this repository): geometric-yukawa-full.md; cp-complex-structure-full.md; measure-dissolution-full.md.
4. A. P. Dawid, "Present position and potential developments: some personal views. Statistical theory: the prequential approach," J. R. Statist. Soc. A 147, 278 (1984).
5. C. D. Froggatt and H. B. Nielsen, Nucl. Phys. B 147, 277 (1979).
6. Repository: Quantum-Relational-Network (claims.yml; predictions.yml PRED-002/003; results/).
