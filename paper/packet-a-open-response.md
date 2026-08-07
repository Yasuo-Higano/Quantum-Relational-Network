# Packet A — Open Signed-Response Quotient (adversary packet, v35.4)

- packet: A (Open Response)
- discipline: paper-closed — no source code, no recorded output values.
  Every claim below carries: status / domain / claim (exact quantifiers) /
  input-output types / gauge / falsifier / forbidden interpretation.
- adversary task: independent derivation, counterexample search, scope
  narrowing. `RefutedAsStated` and `Inconclusive` are first-class outcomes and
  will be accepted as-is.
- expected_outputs_in_packet: false (constructions are given symbolically; the
  adversary computes their own values).

## Domain shared by all claims

Finite-mode, number-conserving, quasi-free, Markovian open fermion systems:

```text
H = Σ_ab h_ab c†_a c_b   (h Hermitian),
loss jumps  L_k = Σ_a l_{k,a} c_a,     gain jumps  G_m = Σ_a g_{m,a} c†_a.
```

Normal covariance C_ab = ⟨c†_b c_a⟩. Node structure is a declared partition of
modes into disjoint index sets; P_i is the orthogonal projector onto node i.
Signed probes prepare C±(0) = C0 ± ε P_i with P_i P_j = 0 for i ≠ j.
The registered observable family is n_j(t) = Tr(P_j C(t)).

A first derivation task (OQ-0, calibration): show that C(t) obeys a closed
affine ODE  Ċ = X C + C X† + Y  and derive X and Y from (h, l, g) in the
covariance convention above, including which of the two dissipative matrices
enters transposed. This is deliberately left to the adversary; sign and
transpose conventions are a known failure mode.

## OQ-1 — Affine cancellation (theorem candidate)

- status: theorem candidate.
- claim: for every X, Y and every pair of initial conditions C⁺(0), C⁻(0),
  the difference D(t) = C⁺(t) − C⁻(t) obeys the homogeneous equation
  Ḋ = X D + D X†; with C±(0) = C0 ± εP_i, D(t) = 2ε e^{Xt} P_i e^{X†t}.
  The affine term Y cancels identically (not approximately) from the signed
  difference at all orders.
- input/output: (two evolved covariance trajectories) → (difference trajectory
  independent of Y and of C0).
- gauge: invariant under all gauges of OQ-7.
- falsifier: any (X, Y, C0, ε, P_i, t) with D(t) ≠ 2ε e^{Xt} P_i e^{X†t}.
- forbidden interpretation: single-sided (unsigned) responses do NOT cancel Y;
  the claim is about the signed difference only.

## OQ-2 — Effective-drift curvature law (theorem candidate)

- status: theorem candidate.
- claim: for i ≠ j (disjoint node projectors), the first time-derivative of
  Δn_j(t) = n_j⁺(t) − n_j⁻(t) vanishes at t = 0, and

  ```text
  ( d²/dt² Δn_j )(0) / (4ε)  =  ‖ P_j X P_i ‖_F²
  ```

  where X is the effective drift of OQ-0. In closed systems (X = −ih) this
  reduces to ‖P_j h P_i‖_F²; in open systems the curvature reads the
  effective drift, not the Hamiltonian.
- input/output type: OpenSignedCovarianceProbe → EffectiveDriftTopology.
- gauge: the statistic is invariant under OQ-7 gauges.
- falsifier: any model in the domain where the identity fails exactly
  (analytically) or beyond discretization error (numerically).
- forbidden interpretation: the right-hand side is NOT ‖P_j h P_i‖² in
  general; promoting it to a Hamiltonian reading without a dissipator
  certificate is the no-go OQ-3.

## OQ-3 — Hamiltonian-promotion no-go (theorem candidate)

- status: theorem candidate; the registered contract is named
  `CurvatureOnlyOpenResponse` = the family {w_ji} of OQ-2 statistics only.
- claim: there exist two models in the domain,
  - M_coh: two nodes, h with a single off-diagonal edge, no dissipation;
  - M_loss: h = 0, one collective loss vector supported evenly (with opposite
    signs) on the two nodes,
  such that every registered statistic of `CurvatureOnlyOpenResponse`
  coincides on M_coh and M_loss, while the Hamiltonian edge sets differ
  (M_coh has an h-edge, M_loss has none). Consequently no function from the
  contract output to Hamiltonian support is correct on both models.
- concrete pair (for refutation attempts): M_coh: h = ((0,1),(1,0)), l = 0.
  M_loss: h = 0, single loss vector l = √2 · (1, −1).
- input/output: EffectiveDriftTopology ↛ HamiltonianTopology (forbidden map).
- unlock (positive side): see OQ-4.
- falsifier: (i) a statistic inside `CurvatureOnlyOpenResponse` distinguishing
  the pair; (ii) an error in the pair (e.g. the two effective drifts do not
  produce equal curvature tables); (iii) a proof that the contract as stated
  is empty or ill-posed.
- forbidden interpretation: the no-go does NOT claim the two models are
  indistinguishable by their full time series (they are not); it is scoped to
  the curvature-only contract. It does NOT forbid promotion when a dissipator
  locality certificate is available (OQ-4).

## OQ-4 — Dissipator-locality unlock (theorem candidate)

- status: theorem candidate.
- claim: if the total dissipative drift matrix Γ (the non-Hamiltonian part of
  X, i.e. X = −ih − Γ/2) is node-block-diagonal, then P_j X P_i = −i P_j h P_i
  for i ≠ j exactly, hence the OQ-2 statistic equals ‖P_j h P_i‖²_F and the
  Hamiltonian reading is licensed. A `DissipatorLocalityCertificate` is a
  machine-checked bound on the cross-node blocks of Γ, not a declaration.
- falsifier: a block-diagonal-Γ model where the reduction fails; or a proof
  that the certificate as specified is not decidable from the model data.
- forbidden interpretation: certificate by declaration; extension to
  non-quasi-free dissipators.

## OQ-5 — Charge-attribution no-go (theorem candidate)

- status: theorem candidate; registered contract = the single scalar
  S = d/dt Tr C at a declared state (one state, one time).
- claim: charge non-conservation (S ≠ 0) does not identify Hamiltonian
  pairing. There exist (a) a pure-loss model (no pairing) and (b) a closed
  BdG pairing model evaluated at a pure BCS state (no dissipation) with the
  same value of the registered statistic. Positive side: for a closed
  number-conserving Nambu Hamiltonian (pairing block Δ = 0) the statistic
  vanishes identically at every state; hence under a verified zero-dissipation
  certificate, S ≠ 0 implies Δ ≠ 0.
- concrete pair (for refutation attempts): (a) two modes, h = 0, loss vector
  with squared amplitude 8 on mode 1 only, state C = diag(1, 0);
  (b) two modes, h = 0, Δ = ((0,1),(−1,0)), pure BCS state with C = I/2 and
  anomalous correlator of maximal modulus; scale conventions must be declared
  by the adversary and calibrated against a dense (2^N) computation.
- input/output: ChargeNonconservingResponse ↛ HamiltonianPairingWitness
  (forbidden map); unlock via DissipativeChargeConservationCertificate.
- falsifier: a proof that within the stated contract the two constructions
  cannot be made to coincide; or that the positive side fails (a Δ = 0 closed
  Nambu model with non-vanishing statistic).
- forbidden interpretation: the full charge trajectory (all t) may distinguish
  loss from pairing; the no-go is scoped to the single-point value. In the
  quasi-free lane with linear jumps, "charge-conserving dissipation" means
  zero dissipation (linear jumps change charge by ±1); quadratic dephasing
  lies outside this lane.

## OQ-6 — Finite-shot robust promotion (instrument-spec claim)

- status: specification claim (semantics), not a physics claim.
- claim: the finite-shot curvature reader must (i) gate every shot series
  through the correlation-granularity gates (ordered shots: split-half and
  transition-count; batch data cannot reach an iid certificate; aggregate
  counts are never assessed), (ii) estimate curvature by the two-scale
  Richardson combination [8Δ̂(δ) − Δ̂(2δ)]/(8εδ²), (iii) carry a registered
  bias bound of the form const · R⁴ δ² with R a declared drift-norm bound and
  the constant derived from the t⁴ remainder, and (iv) answer only when the
  whole propagated confidence interval lies on one side of the declared
  threshold (Straddled and InsufficientObservation are correct non-answers;
  correlated shots are OutOfDomain).
- falsifier: a registered configuration in which the reader as specified
  answers and is wrong with probability exceeding its declared level under
  the iid model; or an inconsistency in the bias-bound derivation.
- forbidden interpretation: point estimates or marginal intervals as
  certificates (forbidden transformations 22/23 of the finite-data core).

## OQ-7 — Gauge family of the full-time quotient (conjecture, open)

- status: conjecture / open classification. Adversary contribution requested.
- known invariances (theorem candidates): the full-time registered response
  table Φ_ji(t) = Tr(P_j e^{Xt} P_i e^{X†t}) is invariant under
  (i) local phase gauge X ↦ D X D† with D = diag(e^{iφ_a}),
  (ii) global frequency shift X ↦ X + iωI,
  (iii) complex conjugation X ↦ conj(X).
- open question: for generic X, is the group generated by (i)–(iii) the FULL
  observational equivalence of the full-time contract, or are there larger
  identifications? Deliverables accepted: a counterexample pair (X, X') with
  identical response tables not related by (i)–(iii), or a proof sketch of
  completeness for stated genericity assumptions, or `Inconclusive` with the
  obstruction named.
- forbidden interpretation: none of (i)–(iii) changes Hamiltonian support;
  OQ-7 therefore does NOT by itself provide a full-time-contract
  Hamiltonian-support no-go (see OQ-8).

## OQ-8 — Scope frontier: full-time Hamiltonian identifiability (open)

- status: open question (the packet's sharpest attack surface).
- question: under the FULL-time contract (all t, all registered node pairs),
  is Hamiltonian support identifiable within the domain, given that the
  curvature-only contract is not (OQ-3)? Either a full-time counterexample
  pair (equal response tables, different h-support) or a positive
  identifiability argument narrows the boundary of OQ-3's scope discipline.
- forbidden interpretation: OQ-3 must not be cited as if it covered the
  full-time contract; conversely a positive answer here would not weaken
  OQ-3 (different contracts).
