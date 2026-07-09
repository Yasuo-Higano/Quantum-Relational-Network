# The Measure Problem That Dissolved: A Pre-Registered Contest of Data-Blind Priors over Magnetized-Torus Yukawa Models

**Draft v1 (v17.12).** 対象誌: JHEP / PRD (letter 級)。companion: [geometric-yukawa-full.md](geometric-yukawa-full.md) [1], [cp-complex-structure-full.md](cp-complex-structure-full.md) [2]。全数値の一次ソースは `results/`。書誌は照合済み (外部参照は [6] を v17.12 で、他は v14.3/v17.2 照合を引き継ぐ)。

---

## Abstract

If the discrete data of a compactification — fluxes, Wilson lines, generation pairing, complex structure — are irreducible, the residual scientific question is the *measure*: which data-blind prior should weigh the model space? We report a complete, pre-registered arc on this question, run over the magnetized-torus Yukawa program of the companion papers. A hierarchy-tilted prior, catastrophically rejected on the rectangular geometry, becomes the program's first surviving selection principle on the CP-capable shear family (+2.3 nats); its anatomy (survival window β ∈ [0.25, 1.7], collapse of hard selection at −8386 nats) and identity (bin-free and binned density-flattening priors reproduce its peak within 0.05 nats) reduce it to a *measure correction*: the configuration population over-represents shallow mass ratios. We then registered six candidate measures with three success criteria fixed in advance — data-blind construction, convention invariance, and non-degradation of never-fitted holdouts — together with a continuum marginalization of the Higgs width (removing the last prior-grid convention) and a moment-based analytic form (which succeeds exactly on the geometry family whose configuration density is moment-compatible, and fails on the heavy-tailed diagonal family). The judgment kills the entire correction class: the only candidate to pass the evidence and stability criteria degrades the unitarity-triangle holdout and is rejected by the third criterion, whose first live use this is. The cause is a lever the corrections had been compensating for: scanning the *imaginary* part of the complex structure finds an interior evidence valley at τ = 1/12 + i/2, where the ten-observable evidence reaches the program's best value, the two chronic tensions (a two-fold Cabibbo excess and a low Jarlskog invariant) resolve simultaneously at MAP (factors 1.01 and 1.24) and at the posterior-band level, all twelve observables sit within a factor of five with six at the ten-percent level — **and every measure correction is now strictly harmful**. The measure problem, posed properly, dissolved into geometry: the uniform counting measure is correct once τ is right, and what remains is a thermodynamic prior's stubborn, sub-threshold positive signal (+0.3 to +0.6 nats, sign-stable across four generations of windows) — the sharpest surviving hint of a dynamical selection principle.

---

## 1. Introduction

The companion papers [1, 2] closed with a deferral: the discrete data of the compactification "await a common selection principle." The first attempt (nine observables, rectangular geometry) rejected or left undecided every candidate — minimum description length, self-stability, hierarchy maximization (at −494 nats), thermodynamic condensation (+0.32, the lone positive sign). This paper reports what happened when the question was asked again, properly: on the CP-capable geometry family of [2], with the Jarlskog invariant in the likelihood, with candidates and success criteria registered before the runs, and with every instrument convention (prior grids, histogram bins) either eliminated or promoted to a tested hypothesis.

The arc has a beginning (a principle survives for the first time), a middle (its identity is found, its class is registered and judged), and an end we did not anticipate: the measure problem dissolves. The discipline that carried the companion papers — pre-registered branch verdicts, regression anchors as gates, negative results reported with the same machinery as positive ones — is what made the dissolution *visible*: had we tuned candidates window-by-window, the compensation structure described in §6 would have read as success.

## 2. A principle survives, and is dissected (QRN-SEL-002/003)

On the 21-geometry shear family (both tori sheared by s₁, s₂ ∈ {0..5} at lattice size N = 36), with the ten-observable evidence of [2], the four principles of the first contest were re-run as two-level data-blind priors. Three verdicts repeated (MDL now *inherits the CP penalty* through its preference for the rectangular geometry, −2.16; self-stability rejected again; thermodynamic again positive, +0.41). One reversed spectacularly: the hierarchy-tilted prior e^{β·depth}, worst of the first contest, **survives at +2.31 nats** — the program's first surviving selection principle. Its per-geometry anatomy explains the reversal: on symmetric (equal-shear) geometries its deep tail overshoots the observed hierarchy exactly as on the rectangle (−514.9, reproducing the first contest's −514.4), while on asymmetric geometries the tilt concentrates prior mass where the likelihood lives.

Dissection followed. The survival curve over β peaks at β = 1.0 (Δ = +3.34), collapses below the survival threshold at β* = 2.0, and the analytic β → ∞ endpoint — hard selection of each geometry's deepest configuration — is a catastrophe at −8386 nats. Asymmetric tails are two orders of magnitude shallower than symmetric ones, but both overshoot the observation. "The universe does not select the deepest" survives the change of family; **what survives is a tilt of order one.**

## 3. The survivor is a measure correction (QRN-SEL-004)

Since e^{β·depth} = Π rᵢ^{−β}, the surviving tilt reads as a conversion toward a scale-invariant (log-uniform, Jeffreys-type [6]) measure on mass ratios. The direct test builds data-blind histogram-flattening priors from the configuration population itself (1D in depth; 2D in the two log-ratios; three bin widths): flattening reproduces the peak within **0.05 nats**, all variants landing in a +2.9…+3.3 band, while the naive mechanism is refuted — the configuration density is *not* exponential in depth (measured slope −0.19, not −1). The operative content is a **measure correction**: the configuration measure of the shear family over-represents shallow mass ratios, and removing this over-representation — by tilt, by flattening, in one or two dimensions — is worth +3 nats, with no observation consulted.

## 4. Registration and instruments (measures.yml; QRN-MSR-001/002)

At this point the program registered the question properly. Six candidate measures (uniform baseline; Jeffreys/Fisher; MDL; thermodynamic; depth tilt, β marginalized; histogram flattening) were entered in a ledger with three success criteria fixed before any judgment run:

- **S1 (data-blind)**: constructible without observations, *with every convention — bins, grids — counted as part of the construction*;
- **S2 (convention invariance)**: candidate rankings survive changes of bin width, prior grid, lattice size, and geometry window;
- **S3 (holdout non-degradation)**: a surviving measure must not worsen never-fitted observables, registered in the prediction ledger before the run (the unitarity-triangle angles β and γ here).

Two instruments were built to make the criteria testable. First, the Higgs-width prior grid — the one convention shown earlier to flip sub-nat rankings — was replaced by a continuum marginalization (trapezoid over σ ∈ [1.2, 8.0]; the convergence gate failed at the first grid density and was passed by doubling it, with the failure recorded). The marginalization confirmed that two earlier conclusions were not grid artifacts, and re-confirmed the program's most robust continuous structure: the evidence concentrates at physical Higgs width σ/√(area) = 1/18 across every geometry, lattice, and grid tested. Second, a *bin-free* analytic form for the flattening measure was constructed from sample moments (a quadratic tilt w = ½(x−μ)ᵀΣ⁻¹(x−μ), zero fitting freedom). It succeeds exactly where moments exist: on the asymmetric family the depth marginal is near-Gaussian (quantile-quantile RMS 0.25) and the quadratic tilt reproduces binned flattening to 0.10 nats — S1 fully satisfied, bins eliminated. On the diagonal family the sample variance is dominated by the anomalous deep tail (QQ ≈ 52) and the quadratic tilt destroys itself (−350 nats) — the same tail that killed hard selection, now killing moment estimation. *Whether the measure correction has an analytic form is itself a property of the geometry family.*

## 5. The judgment (QRN-MSR-003)

All candidates were then judged simultaneously, on three windows (asymmetric, diagonal, and the best-evidence rectangular-lattice geometry then known), under the continuum σ marginalization, against the pre-registered criteria. **Survivors: none.**

- MDL: rejected a third time (−1.64, sign-stable).
- Depth: rejected — the diagonal windows' tails collapse it (S2 sign-unstable: +2.84 vs −494).
- Flattening (1D and 2D): rejected — **the best-geometry window reverses its advantage** (−2.15 and −0.38), breaking S2.
- The bin-free quadratic tilt: passes evidence (+1.93) and S2 on its qualified window, and is **rejected by S3** — it degrades the unitarity-triangle holdout from 1.96 to 4.53 rad. This is S3's first live kill, and the sharpest lesson of the arc: *a measure that raises the evidence and a measure that improves predictions are different objects.*
- Thermodynamic: undecided at +0.49 — but sign-positive in *every* window, as it has been in every era (+0.32, +0.41, +0.49).

## 6. The dissolution: geometry was the missing variable (QRN-YUK-027/028/029)

The window that broke the flattening class was the clue. The aspect ratio of the torus — the *imaginary* part of the complex structure, fixed at 1 in every previous version — turned out to be a live lever: scanning it found monotone improvement toward wider tori and then an interior valley bottom at **τ = 1/12 + i/2** (lattice 36×18, diagonal shear (3,3), physical Higgs width 1/18). At this point the ten-observable evidence reaches the program's best value (−18.43, beating the rectangular nine-observable champion by 1.4 nats on the harder scale); the two chronic tensions resolve *simultaneously* — the Cabibbo element, once overpredicted by 2.5, sits at factor 1.01, and the Jarlskog invariant at 1.24, with both 68% posterior bands covering the measurements for the first time; the full scorecard holds all twelve observables within a factor of five, six of them at the ten-percent level.

And at this point, **every measure correction is strictly harmful** (1D flattening −3.73; 2D −1.23; the quadratic tilt disqualified and self-destructive; the depth tilt −1004). Re-examined across the τ_im scan, the pattern is monotone: the better the geometry, the less the correction helps, until it hurts. The +3 nats that the correction class earned in §3 was *compensation for a mis-set modulus* — the configuration measure of the τ_im = 1 section over-represents shallow ratios precisely because that section's geometry cannot reach the observed hierarchy pattern efficiently. Fix τ, and the uniform counting measure is not merely adequate; it is preferred.

We record the arc's conclusion in the form the ledger demands: **the measure problem, as posed, dissolved into geometry.** The surviving question is no longer "which prior on configurations" but "which prior on τ" — one level up, where the program's evidence now concentrates on a specific rational point (τ_re = 1/12 from one scan, τ_im = 1/2 from another) whose arithmetic (n_x = 2n_y; exact three-fold degeneracy at every aspect ratio tested, to 10⁻¹³) belongs to the same open lattice number theory as the companion's index conditions.

## 7. What survives: the thermodynamic whisper

One signal refused to die at any stage: the thermodynamic prior (weight ∝ e^{γ ln‖Y‖_F}, a minimal proxy for condensation energy) is positive in every window of every era — +0.32 (rectangular, nine observables), +0.41 (shear family, ten), +0.44/+0.40/+0.63 (the three judgment windows), +0.51 (the τ valley bottom) — and never reaches the +1 survival line. Four generations of windows, one sign, no verdict. It is the only candidate whose content is *dynamical* rather than statistical (it rewards strong coupling, not deep hierarchy), and its registered upgrade — replacing the Frobenius proxy by a one-loop vacuum energy with the Kaluza-Klein tower — is now the sharpest open item on the selection-principle ledger.

## 8. Limitations

- All judgments are conditioned on the likelihood conventions of the companion program (lognormal, σ = ln 2; ten observables) and on the shear-family model space; "dissolved" means *within this program's model space and observables*, not as a theorem.
- The τ valley was located by two one-dimensional scans (τ_re, then τ_im); the two-dimensional landscape and the arithmetic of its rational points are unexplored.
- The residual weak spots at the valley bottom are m_c/m_t (factor 2.04 — the companion's historic weak spot, untouched by the τ lever), a mild high bias of the |V_cb|/|V_ts| posterior bands, and a global sign: all CP-odd angles come out with magnitudes near measurement but negative sign, suggesting an orientation (τ → −τ) test that the |J|-based likelihood cannot see. [The orientation test was subsequently performed: the anti-unitary branch pair is exact to machine precision, the measured sign(J) > 0 fixes the branch, and the never-fitted γ = +66.8° then falls within the experimental uncertainty of 65.9° ± 3.5°.]
- The thermodynamic candidate's proxy is not yet a one-loop energy; its persistence could still be an artifact of the proxy's correlation with overlap magnitudes.

## 9. Reproducibility

Deterministic and dependency-free as in the companions (Rust std-only; built-in PASS/FAIL gates; machine-validated claim ledger; 469 PASS / 0 FAIL at the period-17 integration preceding this arc). Practices specific to this arc: (i) a *measure ledger* (measures.yml) registering candidates, conventions, and success criteria before judgment; (ii) *holdout-based rejection* (S3) wired into the same run that freezes the holdout baseline before scoring — the ordering is machine-enforced in the program's stdout; (iii) a *mode-table disk cache* with explicitly non-primary status, which cut the marginal cost of a full judgment (fifteen prior variants × three windows × 69 Higgs widths × 1.7M-configuration sums) to about one minute; (iv) convergence gates on numerical integration, with a recorded instance of the gate failing and being passed by refinement rather than relaxation.

## References

*(内部参照 [1]–[5] は本計画。外部 [6][7] — [6] は v17.12 で Web 照合、[7] は v14.3 照合を引き継ぐ)*

[1] This program: "Yukawa Hierarchies from Magnetized Tori without Order-One Coefficients," paper/geometric-yukawa-full.md.
[2] This program: "CP Violation Requires Complex Structure," paper/cp-complex-structure-full.md.
[3] This program, claim ledger and primary sources: claims.yml, measures.yml, predictions.yml, results/ (ids QRN-SEL-002…004, QRN-MSR-001…003, QRN-YUK-027…029 cited in text).
[4] This program: "Anomaly-Search" companion, paper/anomaly-search-full.md.
[5] 本計画のリポジトリ: Quantum Relational Network — 全数値の一次ソース。
[6] H. Jeffreys, "An invariant form for the prior probability in estimation problems," Proc. R. Soc. Lond. A 186 (1946) 453–461.
[7] R. Trotta, "Bayes in the sky: Bayesian inference and model selection in cosmology," Contemp. Phys. 49 (2008) 71–104.
