# Yukawa Hierarchies from Magnetized Tori without Order-One Coefficients: Bayesian Selection of the Number of Extra Dimensions

**Draft v1 (v14.1) + Erratum (v34.1) 適用済み。** 対象誌: JHEP / PRD。骨子は [geometric-yukawa.md](geometric-yukawa.md)、全数値の一次ソースは `results/`、主張の等級と限界は [claims.yml](../claims.yml) (id を本文に併記)。**v34.1 の正誤表 (§Erratum): §5c の素数 Pfaffian 論法は refuted as stated — 撤回文と訂正は文末の Erratum 節を参照。**

---

## Abstract

Froggatt–Nielsen fits of fermion mass hierarchies involve chosen integer charges and marginalized order-one coefficients. We replace both by geometry: generations are the exact zero modes of a lattice Dirac operator on a magnetized torus (degeneracy = flux = 3), sectors are distinguished by discrete Wilson lines, and Yukawa matrices are computed as overlap integrals with no random coefficients. A single torus is shown to be intrinsically too shallow — its hierarchy depth is set by the flux alone — failing the up-quark hierarchy by two orders of magnitude. The factorizable product T²×T² squares the suppression and reproduces five of six mass ratios and two of three out-of-sample CKM elements within a factor of five, with Bayesian evidence exceeding the anarchic bound by 15 nats. Treating the geometry itself as a hypothesis space, evidence from mass ratios alone selects three tori over two, while coarse Wilson-line lattices are always preferred; including the CKM elements in the evidence — computed exactly by a triple sum over the shared Wilson line — resolves the apparent mass–mixing tension in favor of three tori on all nine observables. The diagonal generation-pairing ansatz is itself testable: marginalizing over all inter-torus pairings beats it by 3.7 nats even after the Occam penalty, and the pairing proves irreducible to any Wilson-line refinement, mirror, or flux-orientation flip — it aligns Landau towers, not positions. Within the explicitly scanned finite-lattice family, all exactly threefold-degenerate tilted four-torus candidates were too shallow: inside that four-parameter flux family, integer arithmetic forces unequal magnetic skew scales ((f₊−f₋)² ≥ 1 whenever Pf = 3, a window-free theorem), so suppression only ever acts once. Primality of the Pfaffian alone, however, does not imply unequal magnetic skew scales — an exact integer counterexample with Pf(F) = 3 and FᵀF = 3I₄ has both skew singular values equal to √3 (Erratum v34.1) — so a general exclusion of unpaired index-three constructions remains open pending explicit lattice-arithmetic and metric-compatibility conditions. The pairing is thus established as irreducible discrete data within every construction scanned — translation, reflection, conjugation, orbifold projection, orientation, and the scanned tilted family — while its universal irreducibility is downgraded from established to open. All results are reproducible from a dependency-free Rust repository with a machine-validated claim ledger.

---

## 1. Introduction

The fermion mass hierarchy spans five orders of magnitude in the mass ratios and a further structured pattern in the CKM matrix. The Froggatt–Nielsen (FN) mechanism [1] organizes this by integer charges and powers of a small parameter ε, but every FN fit carries two kinds of hidden freedom: the *choice* of integer charges, and the *marginalization* over order-one coefficients. A Bayesian model comparison performed within this program (claim QRN-YUK-002) found that a free-charge FN model beats structureless anarchy by ln B ≥ 23 even after paying the Occam penalty for ~7.5×10⁵ charge assignments — but the charges themselves remain a selection, and the O(1) coefficients remain a prior.

This paper asks the sharper question: **can both freedoms be replaced by geometry?** The setting is the standard magnetized-compactification picture [2,3]: chiral generations arise as degenerate zero modes of the Dirac operator on a torus threaded by Q units of magnetic flux; their number is fixed by the index theorem, their wavefunction profiles by the geometry, and their Yukawa couplings by overlap integrals — with *no random coefficients anywhere*. Our contribution is to make this picture fully computational and fully Bayesian on a lattice, so that:

1. every Yukawa matrix is a deterministic function of a small set of *discrete* geometric data (fluxes, Wilson lines, generation pairing);
2. the geometric hypotheses themselves — how many tori, how fine a Wilson-line lattice, how the generations pair across tori — are placed in the model space and weighed by evidence;
3. negative results (constructions that exist but cannot fit the data) are computed and reported with the same machinery as positive ones.

Point 3 turned out to be the scientific core of the paper. Over the course of the analysis, three "innocent" assumptions were successively promoted to hypotheses and tested: the sorting convention that labels generations (§2.4), the diagonal pairing of generations across tori (§5b), and the very idea that pairing is a shadow of the factorized approximation, removable by an honest index-3 geometry (§5c). The first was stabilized after an independent reimplementation exposed a machine-epsilon knife edge in the published numbers; the second was *rejected by the data* at 3.7 nats; the third was excluded within the scanned finite-lattice family — the structural primality argument originally offered for it was itself later refuted as stated by an exact integer counterexample (Erratum v34.1), a fourth instance of the same lesson. What survives is a scoped statement: **within every construction scanned, the generation pairing is irreducible discrete data of the compactification**; its universal irreducibility is open.

All numbers in this paper are generated by a dependency-free Rust code base (single shared numerical library, fixed seeds, built-in PASS/FAIL verification against exact results), and every claim is tracked in a machine-validated ledger with explicit evidence files and limitations. The full suite currently reports 314 PASS / 0 FAIL across 59 programs.

## 2. Construction

### 2.1 Zero modes on the lattice

We work on an N×N lattice torus (N = 18 throughout the two-dimensional analysis) with a U(1) hopping Hamiltonian in the Landau gauge, A_y = φx with φ = 2πQ/N². The Q = 3 lowest states form an exactly degenerate band (spread ~10⁻¹³) separated by a finite gap (0.115 at N = 18) — the lattice avatar of the index theorem's Q zero modes (claim QRN-MATTER-001). These three states are the three generations.

### 2.2 Localization and Wilson lines

Within the degenerate band we diagonalize the position operator X̂ = e^{2πix/N}; a generic phase (φ₀ = 0.83) avoids the accidental cos-degeneracy of the symmetric choice. The resulting generation wavefunctions are Gaussians of width ≈ 2.9 sites centered 6 sites apart. A discrete Wilson line shifts all centers rigidly — one lattice site per unit — giving each matter sector (Q, u, d, L, e) an *address*: its Wilson-line integer k ∈ Z₆ (or half-integers for the Z₁₂ refinement).

*Development record.* Our first implementation attempted lattice magnetic translations to realize sector shifts; these close only when N | 2Q, which fails at (N, Q) = (18, 3). Wilson lines implement the same physics exactly and are the construction we retain. We report this dead end because the obstruction (a lattice-specific number-theoretic condition) recurs in §5c in a stronger form.

### 2.3 Yukawa couplings

The Higgs is a periodic Gaussian of width σ_H centered at the origin; Yukawa matrices are the overlap integrals Y_ij = Σ_x ψ̄_i^{(a)} ψ_j^{(b)} φ_H. Given the fluxes, the Wilson lines, and σ_H, *every entry of every Yukawa matrix is determined*. The likelihood on the six mass ratios (m_u/m_t, m_c/m_t, m_d/m_b, m_s/m_b, m_e/m_τ, m_μ/m_τ) and, where stated, the three CKM magnitudes (|V_us|, |V_cb|, |V_ub|), is lognormal with σ = ln 2 — identical to the FN baseline of QRN-YUK-002, so evidence values are directly comparable across the paper.

### 2.4 A reproducibility lesson: label stability

The localized modes are sorted by center to align generation labels across sectors. One mode sits *exactly* on the wrap boundary of the coordinate (computed center ~10⁻¹² from 0 ≡ N), so the sort is a machine-epsilon knife edge: the published product-model evidences are conditioned on which side of the rounding the implementation happens to fall. An independent numpy reimplementation (built to validate the figures) reproduced the single-torus numbers exactly and *disagreed* on the product model by 1.65 nats, which is how the knife edge was found (claim QRN-YUK-006). The fix — snap centers to the half-site lattice before sorting — makes labels convention-stable; the geometry-selection winner and complete ranking are unchanged on both sides of the edge. We adopt the stable convention throughout and quote published-convention values where they are the historical reference. *Fixed seeds do not imply convention-independent reproducibility; degenerate-band sort orders can become physics downstream.*

## 3. The single-torus no-go

The depth of the single-torus Yukawa hierarchy is set by the flux alone: the attainable minimum of σ₁/σ₃ at σ_H = 1 is ≈ 3×10⁻³ (floor at 2.98×10⁻³, lattice-size invariant), two orders of magnitude short of m_u/m_t = 1.3×10⁻⁵. Its full evidence is lnZ = −53.8, *below* the rigorous anarchic upper bound of −35.4 — a structureless model with 18 random coefficients explains the data better than a single magnetized torus with zero random coefficients (claim QRN-YUK-003). This is a principled negative result: geometry with too little structure is worse than no structure.

## 4. The factorizable product T²×T²

Two tori square the suppression: with per-sector Wilson lines on each factor (6¹⁰ configurations) and the diagonal generation pairing (each generation is the product of same-label modes), the attainable set reaches the up-quark ratio. The maximum-a-posteriori (MAP) geometry reproduces five of six mass ratios within a factor of five (the exception is m_c/m_t, ratio 14 — a persistent weak spot we report rather than repair) and, *out of sample*, two of three CKM magnitudes (|V_cb| within 3%). The evidence is lnZ = −20.4 (published convention; −18.8 stable): **15 nats above the anarchic bound with zero random coefficients** (QRN-YUK-003). Against the free-charge FN model (lnZ = −12.2) the residual is −8.2 nats, which we attribute to the coarseness of the Wilson lattice and the pairing ansatz — quantified next.

| model | free parameters | lnZ (6 ratios) |
|---|---|---|
| M0 anarchy | 18 O(1) coefficients | ≤ −35.4 (bound) |
| M1 FN free charges | 10 integers + 18 O(1) | −12.2 |
| M2 FN literature charges | 18 O(1) | −7.6 |
| M2geo single T² | 5 Wilson + σ_H | −53.8 |
| **M2geo² T²×T²** | 10 Wilson + σ_H | **−20.4** |

## 5. Selecting the geometry itself

Treating {number of tori} × {Wilson lattice} as a hypothesis space with uniform priors: evidence from the six mass ratios selects **T³ over T² by +3.0 nats**, and every lattice refinement (Z₆ → Z₁₂) *loses* — the data pays for more dimensions but never for finer Wilson lattices (QRN-YUK-004). The T³ mass-only MAP, however, loses the CKM structure badly, while a 9-observable point evaluation prefers T²: an apparent mass–mixing tension.

The tension is an artifact of point estimation. Computing the *evidence* on all nine observables — exact despite the broken factorization, via a triple sum over the shared left-handed Wilson line with the lepton sector factorized — T³ wins again (−25.56 vs −27.13 published; −21.84 vs −23.61 stable; QRN-YUK-005). The joint MAP places eight of nine observables within a factor of five. Evidence integrates over the attainable set; a bad point does not condemn a good model. This is the third instance in this program (after the anarchic-evidence and calibration episodes) where the *measurement method*, not the model, generated the apparent physics.

## 5b. The generation pairing is a physical degree of freedom

The product construction silently assumed *diagonal pairing*: generation i on torus 1 pairs with generation i on torus 2. Since the (3,3)-flux T⁴ actually has nine zero modes, choosing three diagonal products is a projection ansatz. We promote the pairing to a hypothesis: σ_F ∈ S₃ per matter field (gauge-fixing the global relabeling), with uniform prior — an Occam penalty of 4 ln 6 = 7.17 nats.

**The data rejects the diagonal ansatz.** Marginalizing the pairings *gains* 3.7 nats net of the penalty, on masses alone and on all nine observables (QRN-YUK-007). The MAP pairing is non-diagonal. The geometry selection survives: with pairings marginalized on both sides — the T³ case requiring a certified-truncation evidence sum over ~10¹³ terms, with a rigorous remainder bound of width ~10⁻¹⁰ nats — three tori still beat two (−17.89 vs −19.86; QRN-YUK-008).

What *is* the pairing, geometrically? An exhaustive classification against 54 candidate realizations (all Z₁₈ Wilson refinements, mirror reflection, complex conjugation = flux-orientation flip) matches none of the 30 non-trivial pairing states (QRN-YUK-009): the pairing permutes *which Landau towers are bundled into one four-dimensional field*, an internal alignment invisible to position or orientation. A uniform Z₃ shift-orbifold does derive the 9→3 projection — the projector's eigenspaces are exactly anti-diagonal pairing families — but its gauge-invariant prediction (all sectors in the same family) loses to the mixed-parity data by 2.0 nats (QRN-YUK-010); a field-dependent magnetization-orientation model narrows the gap but stops 1.0 nat short of the abstract S₃ pairing at identical prior volume (QRN-YUK-011). The mechanism ladder — diagonal (−23.6) < uniform orbifold (−21.8) < orientation (−20.9) < S₃ pairing (−19.9) — measures in nats how far each physical dressing falls short.

## 5c. The honest-index alternative: finite-scan exclusion (and an erratum on the prime-3 argument)

One alternative would dissolve the pairing entirely: tilt the four-torus flux so that the Dirac index itself equals three. With flux data (Q₁, Q₂, t, s) on the coordinate 2-planes, the index is Pf(F) = Q₁Q₂ + ts; (2,2,1,−1) gives 4 − 1 = 3, realized on the lattice as an exactly three-fold degenerate lowest band (spread 10⁻¹³, verified against a t = 0 control with degeneracy 4). Its zero modes are non-factorized across the tori (inter-torus mutual information > 0): the geometry itself bundles the towers. If this construction fit the data, σ would be exposed as a shadow of the factorized approximation.

It does not fit, for two instructive reasons.

**Lattice index arithmetic.** Exact lattice degeneracy is *not* implied by the continuum index: of seven index-3 flux matrices scanned, only three retain an exact triple at N = 18, and the flux (2,2,1,−1) that is exact at N = 6 *splits* at N = 18 by 4×10⁻². Which (flux, N) pairs protect the degeneracy is a number-theoretic question we leave open (QRN-YUK-015, finding 1).

**The family-level flattening (erratum v34.1 — the original "prime-3" argument is refuted as stated).** The flux two-form decomposes into two magnetic eigenplanes with skew-eigenvalues f₊, f₋ satisfying f₊² + f₋² = Σ_{i<j}F_ij² and f₊f₋ = |Pf(F)|. The original version of this section argued that the index constraint f₊f₋ = 3 with 3 prime forces extreme asymmetry; that inference is **false in general** — an exact integer counterexample with Pf(F) = 3 and FᵀF = 3I₄, hence f₊ = f₋ = √3, is exhibited in the Erratum below. What *is* true, and what the scan actually probed, is a family-level statement: within the scanned four-parameter family (two coordinate-plane fluxes plus one tilt pair, Pf = Q₁Q₂ + ts), integer arithmetic forces (f₊−f₋)² = Σ − 6 ≥ 1 — the equal-scale point Σ = 6 has no integer solution in this family (a window-free enumeration, `v341_yukawa_erratum` [E3]) — so suppression only ever acts once *inside the family*. Every exactly-degenerate member of the family sits at the *single-tower* depth floor (ln r₁ ∈ [−5.8, −5.2], versus the −11.3 required), and the best evidence achieved — with a purpose-built sparse Lanczos eigensolver scaling the computation to 2×10⁵ dimensions at 6,700× the dense-Jacobi speed — is lnZ₉ = −34.9, fifteen nats short of the S₃-pairing model (QRN-YUK-012–015). The T²×T² ansatz achieves its depth by *multiplying* two tower suppressions, which the scanned family cannot do; whether some integer flux outside the family (the counterexample class has all six components non-zero) admits exact lattice degeneracy *and* two-tower depth is open, and is precisely the lattice-arithmetic question of finding 1 joined to a metric-compatibility question.

**Conclusion of the arc.** The pairing survives every attempted reduction — translation, reflection, conjugation, orbifold projection, orientation, and the scanned tilted family — with the last alternative excluded by finite-scan results (the structural exclusion originally claimed was refuted as stated; see Erratum v34.1). We conclude that within every construction scanned the generation pairing is **irreducible discrete data of the compactification**, of the same standing as discrete torsion or Wilson-line moduli; a general no-go for unpaired index-3 constructions remains open. "Why this pairing" is henceforth the same kind of question as "why this flux."

## 6. Limitations and outlook

- m_c/m_t remains the single factor-5 outlier of the MAP throughout (ratio 6–14 depending on the model) — the known weak spot of the M2 class.
- The scan window is {T¹, T², T³} × {Z₆, Z₁₂}, identical tori (N = 18, Q = 3), a single Gaussian Higgs profile, and CP phases unused. The full-9 evidence excludes T³×Z₁₂ (a 5×10⁹-term triple sum) explicitly.
- The flattening statement of §5c is a theorem *of the scanned family* ((f₊−f₋)² ≥ 1, window-free within the family), not of primality: the bare prime-Pfaffian inference is refuted by an exact counterexample (Erratum v34.1). Both the number theory of exact lattice degeneracy (which (flux, N) protect the index) and the metric-compatibility conditions under which a family-external equal-scale flux could realize two-tower depth are open and, we believe, interesting in their own right.
- The upstream discrete data — fluxes, Wilson lines, and now the pairing — await a common selection principle (moduli stabilization / vacuum selection), which lies beyond the present methods.

*Addendum (v17.1).* Both deferred questions — the unused CP phases and the selection principle — are taken up in the companion paper [paper/cp-complex-structure-full.md]: CP violation proves to *require complex structure* (the rectangular geometry's J = 0 is a structural zero, reversed by +306 nats once J is admitted as an observable, with the |V_td| holdout miss of this paper's model turning into a 5% hit), and the first surviving data-blind selection principle reduces, on dissection, to a measure correction on the configuration space.

## 7. Reproducibility

Everything is deterministic: no random Yukawa coefficients anywhere, fixed seeds where sampling occurs (Bayesian baselines), a dependency-free Rust implementation with a single shared numerical library, and built-in [PASS]/[FAIL] verification in every program (314 PASS / 0 FAIL at the time of writing). Three practices proved load-bearing and we commend them: (i) *claim ledgers* — every statement carries a grade (C0 premise … C5 interpretation), evidence files, and explicit limitations, machine-validated in CI; (ii) *independent reimplementation* — a figure-generation cross-check in a different language found the label knife edge of §2.4; (iii) *certified truncation* — evidence sums too large to enumerate are bounded rigorously, with the remainder interval reported (widths ~10⁻¹⁰ nats here). Errata are appended, never overwritten: two of our own intermediate claims (a Wilson-refinement equivalence conjecture and an early tolerance choice) were corrected in the record by later versions.

## Erratum (v34.1): the prime-Pfaffian argument is refuted as stated

The originally published Abstract and §5c asserted:

> "index 3, being prime, forces one magnetic eigenplane to flatten, capping the attainable hierarchy at the single-tower floor"

and concluded that the unpaired index-3 alternative "is excluded structurally". This inference is **false**. An independent cross-model clean-room replication (Quantum-Relational-Network-FollowUp, final report 2026-08-02, its claim YUK-005; verdict on the paper set: *Partially Replicated*) exhibited the exact integer counterexample

```text
F = ⎡  0   1   1   1 ⎤
    ⎢ −1   0   1  −1 ⎥      Pf(F) = 3,   FᵀF = 3·I₄
    ⎢ −1  −1   0   1 ⎥  ⇒   f₊ = f₋ = √3  (equal skew singular values)
    ⎣ −1   1  −1   0 ⎦
```

verified in this repository by exact i64 arithmetic (`v341_yukawa_erratum` [E1]; ledger entry QRN-YUK-034, status refuted_as_stated). Primality fixes the *product* f₊f₋ = |Pf| but not the *inequality* of the factors, because skew singular values of an integer alternating form need not be integers. Primality of the Pfaffian alone therefore does not imply unequal magnetic skew scales.

**What survives, sharpened.** Within the four-parameter family actually scanned (two coordinate-plane fluxes plus one tilt pair, Pf = Q₁Q₂ + ts, i.e. F₁₃ = F₂₄ = 0), the flattening *is* a theorem: the equal-scale condition Σ_{i<j}F_ij² = 2|Pf| = 6 has **no** integer solution in the family (window-free enumeration — Σ = 6 bounds every entry by 2), hence Pf = ±3 forces (f₊−f₋)² = Σ − 6 ≥ 1, with the minimal gap f₊ − f₋ = 1 attained at (2,1,1,1)-type points and the seven scanned fluxes at gap² ∈ {4, 13, 25, 49} (`v341_yukawa_erratum` [E3]). The counterexample has all six components non-zero and lies outside the family: the finite-scan observations (QRN-YUK-015, now split into a family-scoped claim) were correct; only their promotion to a universal primality argument was wrong.

**The general lesson.** Pf(F) is a basis-invariant algebraic/topological datum of the integer alternating form (invariant under unimodular change of lattice basis, Pf(SᵀFS) = det(S)·Pf(F)); the skew singular values are *metric-dependent geometric* data. The same form with Pf = 3 realizes both equal scales (F above, Σ = 6) and unequal scales (SᵀFS with a unit shear, Σ = 8, f₊−f₋ = √2) — `v341_yukawa_erratum` [E2]. Deriving skew-scale statements from the Pfaffian requires an explicit metric/lattice-compatibility bridge; promotions of topological invariants to metric anisotropies without such a bridge are henceforth flagged in this program's type ledger.

**Corrections applied (v34.1).** Abstract and §5c re-scoped to the scanned family; "excluded structurally" withdrawn in favor of "excluded within the scanned family, general exclusion open"; QRN-YUK-015 split (family-scoped finite result retained; universal inference registered as QRN-YUK-034, refuted_as_stated); QRN-META-013 corrected accordingly. The FollowUp counterexample is registered as an exact regression in `sim/src/bin/v341_yukawa_erratum.rs`.

## References

*(書誌は 2026-07-05 に Web 照合済み。[1] は正典として据え置き)*

[1] C. D. Froggatt and H. B. Nielsen, "Hierarchy of quark masses, Cabibbo angles and CP violation," Nucl. Phys. B 147 (1979) 277.
[2] C. Bachas, "A way to break supersymmetry," arXiv:hep-th/9503030.
[3] D. Cremades, L. E. Ibáñez and F. Marchesano, "Computing Yukawa couplings from magnetized extra dimensions," JHEP 05 (2004) 079, arXiv:hep-th/0404229.
[4] R. Blumenhagen, B. Körs, D. Lüst and S. Stieberger, "Four-dimensional string compactifications with D-branes, orientifolds and fluxes," Phys. Rept. 445 (2007) 1–193.
[5] L. E. Ibáñez and A. M. Uranga, *String Theory and Particle Physics: An Introduction to String Phenomenology*, Cambridge University Press (2012).
[6] R. Trotta, "Bayes in the sky: Bayesian inference and model selection in cosmology," Contemp. Phys. 49 (2008) 71–104.
[7] N. Haba and H. Murayama, "Anarchy and hierarchy: An approach to study models of fermion masses and mixings," Phys. Rev. D 63 (2001) 053010, arXiv:hep-ph/0009174.
[8] 本計画のリポジトリ: Quantum Relational Network (claims.yml, results/, figures/) — 全数値の一次ソース。
