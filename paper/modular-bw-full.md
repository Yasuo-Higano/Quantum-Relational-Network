# Exact Brillouin-Zone Moment Formulas for the Lattice Bisognano–Wichmann Normalization of Three-Dimensional Staggered Fermions

**Draft v1 (v25.2, frozen).** 対象誌: J. Stat. Mech. / PRB。全数値の一次ソースは `results/` (証明書 `results/v252_bz_certificate.json`)、主張の等級と限界は [claims.yml](../claims.yml) (QRN-C0-007, QRN-GRAV-030..034)。**本稿は v25.2 凍結版から生成される固定物であり、主張の拡張は行わない。**

---

## Abstract

The ground-state entanglement Hamiltonian (EH) of free staggered lattice fermions in 3+1 dimensions, for a half-space cut, is known numerically to take the Bisognano–Wichmann (BW) linear form with a direction-dependent normalization: the nearest-neighbour bond coefficients grow as πξ/λ with λ_x ≠ λ_⊥, a lattice (regulator) property rather than a physical temperature. We close the analytic problem of this normalization. First, an exact block reduction maps the half-space problem, transverse momentum by transverse momentum, onto semi-infinite one-dimensional staggered chains with mass μ(k_y,k_z) = √(cos²k_y + cos²k_z), verified by an operator-level sum rule to 4.5×10⁻⁷. Second, the EH of the semi-infinite staggered chain is exactly known (Eisler, J. Stat. Mech. (2025) 013101): H_EH = 4κK(κ′)·T with κ = 1/√(1+μ²). Translating to bond-gradient normalization gives a single prefactor function for both the hopping and the mass channel, g(μ) = g_m(μ) = (2/π)κK(κ′) = 1/AGM(1, √(1+μ²)). The direction-dependent normalizations are then exact Brillouin-zone moments, 1/λ_x = ⟨g(μ)⟩_BZ and 1/λ_⊥ = ⟨2cos²k_y · g(μ)⟩_BZ, giving λ_x = 1.185467287349258… and λ_⊥ = 1.229428764341310…, in agreement with direct 3D lattice measurements to 6×10⁻⁷ and 1×10⁻⁶. Third, we prove the sign of the anisotropy: since g is strictly decreasing and E[cos²k_y] = 1/2, the difference r_⊥ − r_x = 2 Cov(cos²k_y, g) is negative by the Chebyshev correlation inequality, so λ_⊥ > λ_x is forced by the geometry of the channel weights; interval arithmetic with outward rounding certifies λ_⊥ − λ_x ≥ 0.04391. We delimit the claims explicitly: the 1D closed form is prior work; the contributions here are the exact 3D reduction, the BZ moment formulas with certified enclosures, the anisotropy theorem, and a machine-validated boundary between the two. The normalization λ is not an entanglement temperature (in 1+1D the same estimator gives λ = 1 from the first lattice layer) and its area-entropy analogue is not a Newton constant.

## 1. Introduction

For a relativistic QFT and a half-space (Rindler wedge), the Bisognano–Wichmann theorem states that the modular (entanglement) Hamiltonian is the boost generator, H = 2π∫ x T₀₀ dx: a linear deformation of the physical Hamiltonian. On a lattice, the discretized version of this form is known to describe integrable chains in their gapped phase through the corner transfer matrix (CTM) correspondence, and, for free-fermion chains, exactly solvable commuting-operator structures exist for several geometries. For *gapped* chains the BW form holds up to a *nonuniversal, dispersion-dependent prefactor* [Eisler 2024/2025].

In a numerical program on 3+1D staggered fermions (quantum relational network repository, v22–v25), the half-space EH was computed directly from lattice correlation matrices with double-double precision and clamp-calibrated modular kernels. Three facts emerged: (i) the BW linear *form* survives to five digits near the cut; (ii) the normalization is direction dependent, λ_x = 1.185468 vs λ_⊥ = 1.229430 (transverse), stable across scaling windows; (iii) the same estimator applied to two 1+1D discretizations (XX and Wilson) gives λ = 1.0000–1.0002 from the first layer, identifying λ as a property of the 3+1D staggered regulator, not a continuum temperature. An earlier numerological reading λ = 32/27 was rejected by direct large-N computation (deviation 187× the convergence residual), and the associated "lattice entanglement temperature" and "bare Newton constant" interpretations were retired. This paper reports the analytic closure of what remains: the exact origin of λ_x, λ_⊥.

## 2. Exact reduction to one-dimensional staggered chains

The Hamiltonian is the standard staggered discretization on an L³ lattice (open in x, periodic in y,z; hopping ±1/2; η-phases η_y = (−1)^x, η_z = (−1)^{x+y}). Because the transverse part is x-independent, a transverse Fourier–taste rotation block-diagonalizes the problem: for each transverse momentum pair (k_y, k_z) in the halved Brillouin zone, the two taste bands reduce exactly to two decoupled 1D staggered chains with masses ±μ,

> μ(k_y, k_z) = √(cos²k_y + cos²k_z),  H_chain = Σ (1/2)(c†_x c_{x+1} + h.c.) + μ Σ (−1)^x c†_x c_x,

and charge conjugation c → (−1)^x c† maps K(−μ) = −ΣK(μ)Σ (Σ = diag(−1)^x). The half-space EH bond profile is therefore an exact BZ average of chain EH profiles:

> A_x(ξ) = ⟨K^{1D}_μ(ξ, ξ+1)⟩_BZ (uniform),  K_y(ξ) = ⟨(2cos²k_y/μ)·K^{1D}_μ(ξ,ξ)⟩_BZ (mass channel).

This reduction is verified as an operator-level **sum rule** at N = 32: the direct 3D half-space scan equals the 1D chain sum with max|Δ| = 4.5×10⁻⁷ / 1.1×10⁻⁷ / 1.2×10⁻⁷ (A_x/K_y/B_z) inside the certified window ξ ≤ 5 (the floor is the clamp-boundary resolution of the double-double kernel, not the algebra). [claims QRN-GRAV-030]

Defining the normalizations by K_NN = g(μ)·πξ + b and (−1)^x K_xx = g_m(μ)·2πξ_s·μ + b_m (ξ_s = ξ − 1/2 the half-site coordinate), the 3D normalizations become

> **1/λ_x = ⟨g(μ)⟩_BZ,  1/λ_⊥ = ⟨2cos²k_y · g_m(μ)⟩_BZ.**

## 3. The chain prefactor is exactly known

The EH of the semi-infinite staggered chain was solved exactly by Eisler [J. Stat. Mech. (2025) 013101, arXiv:2410.16433, eqs. (64), (67), (72)]: the correlation matrix commutes with the tridiagonal matrix T (t_m = −m/2, d_m = (−1)^m μ(m−1/2)) and

> H_EH = 4κ K(κ′) · T,  κ = 1/√(1+μ²),  κ′ = |μ|/√(1+μ²),

K being the complete elliptic integral of the first kind. Our chain differs by the local unitary c_n → (−1)ⁿc_n (hopping sign), which leaves the staggered term and the bond-gradient magnitudes invariant. Matching normalizations gives a single function for both channels:

> **g(μ) = g_m(μ) = (2/π) κ K(κ′) = (2/π)·K(κ′)/√(1+μ²) = 1/AGM(1, √(1+μ²)) = (2/π)∫₀^{π/2} dθ/√(1+μ²cos²θ).**

The equality g_m ≡ g — observed numerically at the 10⁻⁵ level before the identification — is thus a theorem: both channels carry the same nonuniversal prefactor. Notably g(1) = 1/AGM(1,√2) is Gauss's constant. Small-μ expansion: g = 1 − μ²/4 + 9μ⁴/64 − 25μ⁶/256 + … (an analytic series in μ²; earlier conjectured μ²ln μ terms do not exist). Large-μ: g = (2/π|μ|)·ln(4|μ|)·(1+O(1/μ²)). A pre-registered screening had rejected the spacing-type candidate π/(2K(κ′)) (drift 0.16% → 3.8% over μ ∈ [0.5, √2]); the exact solution explains why: the EH level spacing factorizes as [prefactor 4κK(κ′)]×[T-spacing π/(2κK(κ))], and the screening candidate was the spacing function, not the operator prefactor. Verification battery (21 PASS): three-representation agreement to 3.4×10⁻¹⁴, Gauss-constant anchor to 1 ulp, match to the lattice extraction max|ĝ/g − 1| = 1.18×10⁻⁵ over the pre-registered window μ ≥ 0.5 (mass channel: 4.4×10⁻⁷), and a mutation layer (five deliberate errors, all detected). [QRN-C0-007, QRN-GRAV-031]

## 4. Brillouin-zone moments: certified values

With g known, r_x = 1/λ_x and r_⊥ = 1/λ_⊥ are parameter-free numbers. Three implementations sharing no discretization or interpolation:

- **A** (no AGM): 3D torus midpoint rule of the single integrand [1+(cos²k_y+cos²k_z)cos²θ]^{−1/2};
- **B**: 2D BZ midpoint rule with f64 AGM;
- **R**: interval Riemann sum with outward rounding (n = 2¹⁵ cells on [0,π/2]²; interval π; cosine via remainder-padded Taylor series — no libm trust; interval AGM; hierarchical summation; exact power-of-two division).

A and B agree to 3.9×10⁻¹⁶ and both lie inside the rigorous enclosures of R:

| quantity | floating (A=B) | certified enclosure (R) |
|---|---|---|
| r_x | 0.843549215293854 | [0.843545503984956, 0.843552926638375] |
| r_⊥ | 0.813385882130202 | [0.813356957837832, 0.813414806840221] |
| **λ_x** | **1.185467287349258** | [1.185462071698428, 1.185472502995920] |
| **λ_⊥** | **1.229428764341310** | [1.229385046338885, 1.229472484821826] |

The direct 3D lattice measurements (λ_x = 1.185468, λ_⊥ = 1.229430, N ≤ 192 block method) agree to 6.0×10⁻⁷ and 1.0×10⁻⁶. The division of labour is explicit: *proof-grade* digits come from the interval enclosure (width O(h): 7.4×10⁻⁶ / 5.8×10⁻⁵), *reference* digits from two independent double-precision routes. The certificate JSON records formulas, domain, normalization, methods, enclosures, and SHA-256 hashes. [QRN-GRAV-032, QRN-GRAV-033]

## 5. The anisotropy is a theorem

**Theorem.** r_⊥ < r_x, i.e. λ_⊥ > λ_x.

*Proof.* Under the uniform BZ measure, Y = cos²k_y and Z = cos²k_z are independent with E[Y] = 1/2, and h(s) = g(√s) is strictly decreasing (g′(μ) = −μ⟨cos²θ(1+μ²cos²θ)^{−3/2}⟩ < 0). Then r_⊥ − r_x = ⟨2Y·h(Y+Z)⟩ − ⟨h(Y+Z)⟩ = 2Cov(Y, h(Y+Z)). Conditioning on Z (E[Y|Z] = E[Y] by independence), Cov(Y, h(Y+Z)) = E_Z[Cov(Y, h(Y+Z)|Z)], and for each z the function h(Y+z) is strictly decreasing in Y, so each conditional covariance is negative by the Chebyshev correlation inequality. ∎

Interval arithmetic seals the inequality numerically: r_x − r_⊥ ≥ 0.0301307 > 0, hence **λ_⊥ − λ_x ≥ 0.04391** (certified). Physically: the transverse channel weights 2cos²k_y are largest exactly where the effective mass μ is largest and the prefactor g smallest — a negative covariance between channel weight and mass-dependent prefactor. The direction of the anisotropy is thus geometric necessity, not a Lorentz-breaking surprise.

## 6. Finite-size behaviour and instrument floor

Redefining the finite-chain numerics as a *convergence test toward the exact g* (rather than a discovery tool) separates physics from instrument. Measured at μ ∈ {0.5, 1, √2}, N = 64…512: the residual ĝ_N/g − 1 is **N-independent** (fitted exponent p ≈ 0), of size ≤ 8.1×10⁻⁶ — an instrument floor set by the fit window (ξ* = 5–6, the κ-budget edge of the clamped double-double kernel), not a finite-size effect; the mass channel reads the exact value to 5.3×10⁻⁷. The gapped/critical crossover is organized by the single variable x = N·arsinh(μ) (ξ_corr = 1/arsinh μ): for x ≥ 25 the extraction is converged (|d| ≤ 6×10⁻⁶ observed), while for x ≤ 2 the chain is effectively massless and the extraction reads the massless value — quantifying why small-μ points must be excluded from prefactor fits at fixed N. Gates are fail-closed: an insufficient trust window returns "no measurement" rather than a silently extrapolated fit, and a deliberately degraded kernel (coarse clamp pair) is verified to trigger exclusion. [QRN-GRAV-034]

## 7. Scope, prior art, and negative statements

**Prior art.** The 1D closed form g(μ) — the exact EH of the semi-infinite staggered chain — is Eisler's result, not ours. We claim: (i) the exact 3D→1D reduction with its operator-level sum rule; (ii) the BZ moment formulas 1/λ_x = ⟨g⟩, 1/λ_⊥ = ⟨2cos²k_y·g⟩ and their certified evaluations; (iii) the anisotropy theorem; (iv) the identification g_m ≡ g as the two-channel consequence of the single prefactor.

**Negative statements (established in the same program, kept as constraints).** λ is a regulator normalization: direction-, discretization- and operator-normalization-dependent. It is not an entanglement temperature (the same estimator yields λ = 1 in two 1+1D discretizations from the first layer); λ = 32/27 is excluded (N = 96 direct reading, deviation 187× the residual); the associated area-coefficient combination G_entropy^bare(a; regulator, direction scheme) varies by 3.7% between direction schemes and is not a Newton constant. Nothing in this paper depends on, or supports, an induced-gravity interpretation.

**Verification taxonomy.** Of the 42 PASS gates behind this paper: mathematical identities (three-representation agreement, sum rule), known-value benchmarks (Gauss constant, Eisler prefactor), internal consistency (independent implementations, interval containment), negative controls (mutation layers, degraded-kernel exclusion, weight/measure mutants), and convergence/floor measurements. **None are frozen novel predictions** — this paper closes a normalization analytically; it does not predict new data.

## 8. Reproducibility

All computations run from a dependency-free Rust repository (fixed seeds, built-in PASS/FAIL against exact values, machine-validated claim ledger with explicit assumptions and falsifiers; interval arithmetic implemented from scratch with outward rounding and libm-independent π/cos). The frozen artifact set and its SHA-256 manifest are produced by `v252_manifest`; the certificate is `results/v252_bz_certificate.json`.

## References

1. J. J. Bisognano, E. H. Wichmann, J. Math. Phys. 16, 985 (1975); 17, 303 (1976).
2. V. Eisler, "On the Bisognano-Wichmann entanglement Hamiltonian of nonrelativistic fermions", J. Stat. Mech. (2025) 013101 [arXiv:2410.16433].
3. I. Peschel, J. Phys. A 36, L205 (2003) — reduced density matrices from correlation functions.
4. R. J. Baxter, *Exactly Solved Models in Statistical Mechanics* (Academic Press, 1982) — corner transfer matrices.
5. I. Peschel, M. Kaulke, Ö. Legeza, Ann. Phys. (Leipzig) 8, 153 (1999); V. Eisler, I. Peschel, J. Phys. A 50, 284003 (2017) — EH of free-fermion chains.
6. G. Wong, I. Klich, L. A. Pando Zayas, D. Vaman, JHEP 2013(12), 20 — entanglement temperature.
7. Repository: Quantum-Relational-Network (docs/uft-v22.2…v25.2; claims.yml; results/v252_bz_certificate.json).
