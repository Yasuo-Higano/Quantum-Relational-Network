# Operational Locality from Certified Quantum Interfaces: Identifiability, Contextual Descent, and Finite-Data Guarantees — Core Specification

**OCS version 1.0 (spec-frozen at v34.2).** 作業名:「証明付き量子インターフェースからの操作的局所性」。

**Purpose.** This document is a *paper-closed* specification: it defines every object,
procedure, theorem, normative test instance, and verdict vocabulary needed to
re-implement the operational core of this program **without reading any repository
source code and without access to any recorded output values**. It exists because an
independent cross-model clean-room replication (2026-08-02, verdict *Partially
Replicated*) correctly found that the previously supplied papers did not close the
operational core: many claims were reproducible only by repository replay. This
specification is the interface for all future clean-room implementations; the
repository's own implementation is one instance of it, not its definition.

**Discipline.**
- `expected_outputs_in_spec: false` — this document never states measured output
  values. It states *procedures*, *exact identities*, *decision bars*, and
  *qualitative verdicts required at normative instances*. If an implementation and
  this spec disagree, and the spec is found ambiguous, the spec version is bumped
  (OCS 1.0 → 1.1); results are never used to retro-fit the text.
- The companion machine-readable manifest `operational-core-closure.yml` records,
  per claim, whether it is closed by this paper (`paper_closed`) or only replayable
  from the repository (`repository_replay_only`). Its integrity, and this file's
  SHA-256 freeze, are machine-audited.

## 0. Reproducibility grades

Every claim in this program is tagged with the strongest grade it has actually met:

1. **repository replay** — re-running committed code reproduces committed outputs.
2. **paper/spec-closed** — an implementation written from this spec alone (no source,
   no outputs) can execute the procedure and reach the same *verdicts*.
3. **clean-room** — grade 2 achieved by an implementation whose authorship shares no
   code and no numerical kernel with the repository (a shared human operator is
   allowed but must be disclosed; cross-model clean-room is this grade, not grade 4).
4. **organizationally external** — grade 3 achieved by a different human author /
   organization, with protocol frozen before the run and results (including
   failures) public.

Grade 4 count for this program is currently zero and nothing in this spec changes
that; only the external-replication ledger may.

## 1. Kinematics: systems, roles, and declared operations

**S1 (system).** A system is a finite-dimensional complex Hilbert space H = C^d with
its full operator algebra M_d. All procedures below are exact finite-dimensional
linear algebra; the Hilbert–Schmidt (HS) inner product ⟨A,B⟩ = Tr(A†B) and Frobenius
norm ‖·‖_F are the default metric on operators.

**S2 (roles are types).** Operational primitives carry one of four roles, which are
*distinct types with no implicit conversion*:
- **Preparation**: a density matrix ρ (Hermitian, PSD, Tr ρ = 1; convex mixing
  allowed within a declared preparation family).
- **ControlGenerator**: a Hermitian generator G (the only role admissible as input
  to factorization recovery).
- **MeasurementEffect**: an operator E with 0 ≤ E ≤ I. Effects form an operator
  system, **not** an algebra: no product closure may be assumed or required.
- **DriftGenerator**: an uncontrolled Hermitian generator (present but not
  commandable).
Role qualification happens at construction; a role-mixed input to any recovery
procedure is a **construction-time rejection** (not an Abstain verdict): reasons
`role_mixed_recovery`, `no_declared_contexts`, `context_coverage_incomplete`.

**S3 (grading is a type).** A net is either *ordinary* (parity-even lane) or
*Z2-graded* (fermionic lane). An ordinary net rejects odd (fermionic) primitives at
construction (`grading` mismatch). In the graded lane the certified bracket is the
anticommutator on odd–odd pairs and the commutator otherwise. This gate exists
because independent fermionic modes have odd generators whose *ordinary* commutators
are large (Jordan–Wigner strings), which fabricates a complete graph out of an empty
one; the graded bracket of independent odd generators is exactly zero.

**S4 (declared vs accessible).** A *DeclaredOperation* is a labeled matrix with a
role intent. It has **no** operational standing. The only path to an
*AccessibleOperation* is through a certificate gate (§2). There is no default
coercion Declared → Accessible (Forbidden Transformation FT-14).

## 2. Certificates: the gate from declared to accessible

**C1 (certificate anatomy).** Every certificate carries: (i) its **target binding**
= a collision-resistant hash of the exact target matrix (SHA-256 of a canonical
serialization); (ii) its **scope** = the calibration data set / time window / noise
assumption under which it was produced; (iii) its **content** (below). Using a
certificate whose target hash does not match the operation it is attached to is a
construction-time rejection (`certificate_target_mismatch`) — certificates are
non-transferable (no reuse across targets, FT "origin binding").

**C2 (IndependentAddressabilityCertificate).** Input: command family
{u_1..u_k} with declared target generators {G_1..G_k} (one-to-one). Procedure:
HS-normalize the response directions actually generated per command; form the k×k
response matrix M_{jk} = ⟨Ĝ_j, D̂_k⟩ (HS inner products of unit-normalized targets
and generated directions); qualification requires:
- rank(M) = k with minimum singular value **σ_min ≥ 0.5** (frozen bar) —
  mathematical decomposability of a span does *not* qualify (FT-15: a device whose
  only knob drives G_1+G_2 has rank 1 and cannot certify {G_1, G_2});
- cross-talk = max_{j≠k} |M_{jk}| with **interval semantics**: qualified if the
  whole cross-talk interval ≤ 0.1 (frozen bar); rejected if entirely above;
  **Straddled** (abstain, no forced verdict) if the interval crosses the bar.
Failure vocabulary: `degenerate_targets`, `command_target_count_mismatch`,
`insufficient_command_rank`, `sigma_below_bar`, `crosstalk_excess`,
`crosstalk_margin_straddled`.

**C3 (SynthesisCertificate).** Input: a base of already-accessible generators and a
target G. Content: an explicit word (sequence of scaled sums and Lie brackets) whose
machine re-execution reproduces G with residual ≤ **1e-9** (frozen bar) relative to
‖G‖_F. The Lie closure is computed by iterated brackets with rank-revealing span
growth; if the target is not in the closure: `no_synthesis_path`; if the word's
residual exceeds the bar: `synthesis_residual_excess`. Accessibility is a relation
to an interface: the same G_1 can be unreachable from {G_1+G_2} and synthesizable
from {G_1+G_2, G_3} — reachability claims always name their base.

**C4 (TomographyCertificate).** Input: an informationally complete preparation
family (the ρ's span Hermitian d×d space) and measured frequencies for a declared
effect. Content: the linear-inversion reconstruction; qualification requires
reconstruction residual ≤ bar and the reconstructed E to satisfy 0 ≤ E ≤ I
(`not_informationally_complete`, `tomography_residual_excess`,
`effect_qualification_failed`).

**C5 (origin).** Every AccessibleOperation records its origin as exactly one of
{DirectlyCalibrated, Synthesized, TomographicallyInferred} with the corresponding
certificate embedded. Origin coverage of admitted operations is a scoring metric
(coverage 1.0 required in holdouts); raw promotions (admissions without origin) must
be zero.

## 3. Contexts: role-typed, witnessed compatibility

**X1 (four context types).**
- **ControlContext**: a declared subset of control primitives jointly executable.
- **MeasurementContext**: a set of effects with a **JointMeasurementCertificate** —
  an explicit joint POVM {R_ω} with R_ω ≥ 0, Σ_ω R_ω = I, whose marginals reproduce
  the member effects (`joint_candidate_not_positive`, `joint_sum_not_identity`,
  `marginal_mismatch`).
- **PreparationFamily**: a convex-reachability certificate — explicit convex weights
  reproducing the family's declared states (`weights_invalid`, `mixture_mismatch`).
- **DriftRegime**: a stability certificate — drift variation within the declared
  window ≤ bar (`drift_regime_unstable`).

**X2 (joint measurability is wider than commutativity).** Compatibility of effects
is certified by a joint POVM, never inferred from commutation. Normative instance
N3 (§10) has non-commuting unsharp qubit effects that *are* jointly measurable at
noise η = 0.6 and *not* at η = 0.8: the qualification criterion for the unbiased
orthogonal pair is η_x² + η_z² ≤ 1 (Busch), and the canonical joint candidate is
R_{ab} = (I + a η_x σ_x + b η_z σ_z)/4, a,b ∈ {±1}, whose positivity is exactly
this criterion. Conversely, algebraic commutation certificates do **not** create a
joint context: a net whose declared contexts are all singletons yields
Abstain(`operational_compatibility_unwitnessed`) even if all pairs certifiably
commute (FT-12).

## 4. Resource budgets: a five-component partial order

**R1 (budget).** A ResourceBudget is the 5-tuple (max_duration, max_amplitude,
max_bandwidth, max_depth, max_error), each a non-negative finite real. The order is
**componentwise only**; there is no total order and no scalarization. Collapsing the
budget to any weighted scalar is forbidden (FT-18): it can reverse factorization
readings by making incomparable budgets comparable.

**R2 (profile).** Given a finite declared grid B of budgets, the
*operational factorization profile* is the map b ↦ reading(b) where reading(b) is
the recovery verdict (§5) of the sub-net of operations whose budgets are ≤ b
(componentwise). Vocabulary per grid point: `no_accessible_operations` (honest
record of insufficiency), a construction rejection, or a reading.

**R3 (promotion rule, FT-17).** A profile point may be promoted to a *stable* local
factorization only if its (reading + gauge-orbit) equivalence class contains a
comparable chain of length ≥ 2 in the grid. Single-point (transient) classes are
never promoted — including the top of the grid, which is transient *relative to the
declared grid*. Transient promotions must be zero.

**R4 (invariance).** The profile is invariant under componentwise strictly monotone
reparameterizations of the budget axes.

## 5. Factorization recovery: marked nets, verdicts, orbits

**F1 (marking carries locality).** The input to recovery is a *marked operational
net*: the set of accessible control primitives, their certified pairwise
(graded) commutator **intervals**, and the declared contexts. The generating matrix
set used by recovery must be *identically* the net's admitted primitives (no side
channel input). A single closed global algebra carries no factorization: closures of
site-local and of Fourier-transformed generator families can be the same full M_d
(erasure no-go T1); there is no map GlobalClosure → factorization (FT-11).

**F2 (frozen decision procedure).** With commutation threshold **τ = 1e-3** (frozen):
1. Qualification: control role only; ≥ 1 declared context; every primitive belongs
   to some context (else the three construction rejections of S2).
2. Graph: vertices = primitives; edge iff the certified bracket-norm interval lies
   entirely above τ; non-edge iff entirely ≤ τ; if any interval straddles τ →
   Abstain(`commutator_margin_straddled`).
3. Components of the non-commutation graph.
4. **Witness gate**: every unordered pair of components must share at least one
   declared context (joint-context witness). Otherwise
   Abstain(`operational_compatibility_unwitnessed`) — commutation certificates alone
   never pass this gate (FT-12).
5. Closure: per component, the complex associative *-closure of its primitives (span
   growth by products, rank-revealing). Factor test: the closure must be a full
   matrix algebra of some dimension d_i on its support
   (else Abstain(`component_not_factor`) or, if generators are simply too few,
   Abstain(`insufficient_operational_generators`)).
6. Center: compute the center of the joint closure. Trivial center and Π d_i = d →
   **Exact** up to local unitary × component permutation, with local_dims {d_i}.
   Nontrivial center → decompose H ≅ ⊕_α (C^{m_α} ⊗ C^{n_α}) and return
   **SuperselectionSectors** {(m_α, n_α)} — forcing a tensor product across a
   nontrivial center is forbidden.
7. If distinct qualified constructions yield readings in different gauge orbits →
   **EquivalenceClassOnly** (no tie-breaking, ever).

**F3 (gauge orbit and matching).** The Exact reading is defined up to local unitary
× permutation. Concretely, a factorization is compared as the *set* of factor
subalgebras: for factors A, B (as subspaces of traceless operators under HS), the
overlap is ‖Π_A Π_B‖_F² / dim(A) with Π the HS-orthogonal projector; two readings
match iff some dimension-preserving permutation matches all components with overlap
≥ the frozen orbit bar **0.9** (equal subalgebras give 1 up to numerical dust;
genuinely different markings, e.g. site vs Fourier, fall far below). Non-matching
Exact readings coexist only as EquivalenceClassOnly.

**F4 (verdict vocabulary, frozen).** `exact_up_to_local_unitary_and_permutation
{local_dims}`, `superselection_sectors {(m_α,n_α)}`, `equivalence_class_only`,
`abstain(reason)` with reasons `insufficient_operational_generators`,
`commutator_margin_straddled`, `grading_mismatch`, `component_not_factor`,
`operational_compatibility_unwitnessed`.

## 6. Contextual descent: charts, glue, cocycles

**G1 (chart).** A chart is a declared subset of the net's primitives. Chart-local
recovery runs F2 on the sub-net, with commutator certificates *inherited* from the
global net (never recomputed locally) and the witness gate evaluated within the
chart. A chart speaks only about its own support.

**G2 (glue).** Charts A, B overlap on shared primitives. The glue condition is that
the factor subalgebras induced on the overlap **match** (orbit overlap = 1 within
the exact-certificate bar 1e-9). A consistent covering atlas (all pairwise
matchings hold, factor dimensions multiply to d, cross-chart factor pairs
witnessed) yields **GluedExact**.

**G3 (glue theorem, T5 — positive).** For a consistent covering atlas, the glued
reading equals the direct global recovery (F2 on the full net): same verdict, same
local_dims, same gauge orbit. Implementations must verify this equality on their
own instances; the spec fixes the *procedure*, not the numbers.

**G4 (cocycle failure).** If an entangler-twisted chart breaks matching on an
overlap, and no covering consistent sub-atlas exists, the verdict is
Abstain(`glue_inconsistent`) — even when *every* chart is locally Exact
(FT-19: ChartLocalFactorization ↛ GlobalFactorization). If several maximal
consistent covering atlases exist with non-matching global readings →
`equivalence_class_only {n_consistent_atlases}`. Remaining reasons:
`chart_failed(local failure)`, `coverage_incomplete`, `compatibility_unwitnessed`.

## 7. Graded lane: Majorana frames, charge witnesses, Dirac modes

**M1 (frame qualification).** A Majorana frame is 2N Hermitian operators γ_a with
{γ_a, γ_b} = 2δ_ab (CAR bar 1e-9). Rejections: `not_hermitian`, `car_violation`,
`odd_count`.

**M2 (Dirac pairing no-go, T6a).** Odd CAR data is invariant under the orthogonal
group O(2N) acting on frames; the pairing of Majoranas into complex modes is *not a
function of CAR data*. Without a witness the reading is **MajoranaFrameOnly**
(the O(2N) orbit); there is no map MajoranaFrame → ComplexModeFactorization
(FT-20). Normative check: two different pairings both satisfy mode-CAR exactly
while assigning different occupation to the same state — CAR cannot choose.

**M3 (charge witness recovery, T6b — positive).** Let Q be a declared even witness
whose adjoint action closes linearly on the frame:
J_{ba} = ⟨γ_b, i[Q, γ_a]⟩ / ‖γ‖² (HS, with ‖γ‖² = Tr γ² = d). Qualification: the
residual of i[Q,γ_a] off span_R{γ_b} ≤ bar (`witness_not_linear_on_frame`
otherwise), J real antisymmetric, and **J² = −I** within 1e-9
(`complex_structure_unresolved` otherwise — e.g. a partial charge witness).
Then the ±i eigenplanes of J define N complex modes a_j (unique up to U(N) gauge);
implementations must verify mode-CAR and Σ_j a_j†a_j = Q on their instances.

**M4 (no fabrication).** The ungraded recovery of §5 applied to a graded net must
abstain (`component_not_factor` after graded bracket evaluation), not fabricate
locality from Jordan–Wigner strings; the witness path M3 is the only gate to modes.

## 8. Response probes are types (the replication-driven separations)

**P1 (SignedInitialCovarianceProbe — one-particle lane).** Contract: quadratic
(one-particle) dynamics with single-particle Hamiltonian h on C^n; probe states are
covariances Γ± = Γ0 ± εP_i with Γ0 = I/2 and P_i, P_j orthogonal projectors
(source, target; P_iP_j = 0); readout n_j(t) = Tr(Γ(t)P_j),
Γ(t) = e^{iht} Γ e^{−iht}. Exact identity (T8, all ε ≠ 0, exact in finite dims):

  ( n̈_j[Γ+](0) − n̈_j[Γ−](0) ) / (4ε) = ‖ P_j h P_i ‖_F².

**P2 (HamiltonianQuench is a different experiment, T9).** Quenching the generator
(h ± εP_i) with a *fixed* state is not the same probe and does not satisfy P1's
identity. Normative instance N1: n = 2, h = σ_x (coupling weight ‖P_2 h P_1‖² = 1);
with Γ0 = I/2 the quench response vanishes identically (the state commutes with
every generator) while the signed-probe response equals 1. Required verdicts: the
probe-type field is part of the observation contract; applying a P1-certified law to
quench data is OutOfDomain (FT: SignedInitialCovarianceProbe ↛ HamiltonianQuench).
The phrase "local perturbation" without a probe type is not admissible in any
contract.

**P3 (NumberConservingResponse vs pairing, T10).** Contract: many-body (Fock)
dynamics; probe states are *product* density matrices with mode occupations
p_i = 1/2 ± ε (source mode) and 1/2 elsewhere; readout ⟨n_j(t)⟩; response =
(n̈_j⁺(0) − n̈_j⁻(0))/(4ε). For number-conserving quadratic H (hopping t between
modes 1,2) the response equals |t|²; adding a *diagonal* density–density term
preserves it. A pairing term breaks the law: normative instance N2, two modes,

  H = c₁†c₂ + c₂†c₁ + Δ (c₁†c₂† + c₂c₁),  source mode 1, target n₂:
  response = **1 − Δ²** exactly
  (hopping sector contributes +1 via [T,[T,n₂]] = −2(n₁−n₂)·(sign),
   pairing sector −Δ² via [P,[P,n₂]] = 2Δ²(n₁+n₂−1), cross terms vanish).

Certificates for number-conserving response laws must carry the charge witness
(‖[H, Σn]‖ = 0); applying them to a Hamiltonian that fails the witness is
OutOfDomain (FT: NumberConservingResponse ↛ BCS/PairingResponse).

## 9. Structured backends: scaling lanes with scope discipline

**B1 (Pauli GF(2) lane).** Operators are ± Pauli strings encoded as (x|z) ∈
GF(2)^{2n}; commutation is the symplectic form ω(v,w) = x_v·z_w + z_v·x_w (mod 2);
the closure of a generator span V has complex dimension 2^{dim V}; the **center is
the radical** of ω restricted to V (radical dimension r ⇒ 2^r superselection
labels). The F2 procedure runs verbatim on this encoding (edges = ω ≠ 0), and on
any instance small enough for both lanes the structured verdicts must **equal** the
dense verdicts, including sector structure. Non-Pauli sums (e.g. X₁+X₂) are not
encodable: constructing them in this lane is a scope error (`not_pauli_string`).

**B2 (Majorana quadratic lane).** Quadratic Hamiltonians are real antisymmetric
matrices A (H = (i/4) Σ A_{ab} γ_a γ_b); blocks = connected components of the
support graph of A; the Lie closure of a block of 2m Majoranas is so(2m);
correspondence principle to dense: the *-closure of a 2m-Majorana quadratic family
has complex dimension **2^{2m−1}** (even Clifford algebra — never full M_d: parity
superselection). Non-antisymmetric input: `not_antisymmetric`.

**B3 (ScopeExceeded is a correct verdict, FT-21).** The dense lane refuses
instances with dim > 4096 (frozen bar): the verdict `dimension_too_large_for_dense`
is a *correct answer*, not a failure. Structured certificates at large n must never
silently fall back to dense computation.

## 10. Normative instances (constructions only — no output values)

- **N1 (quench null).** n = 2, h = σ_x, P_i = |1⟩⟨1|, P_j = |2⟩⟨2|, Γ0 = I/2.
  Required: P1 response = coupling weight (= 1); quench response = 0; the two
  contracts are distinguished by type, not by fit.
- **N2 (pairing).** Two JW fermion modes, H as in P3 with Δ ∈ (0,1); product probe
  states p₁ = 1/2 ± ε, p₂ = 1/2. Required: response = 1 − Δ² exactly (any ε), and
  = 1 for Δ = 0 with or without a diagonal V n₁n₂ term.
- **N3 (Busch pair).** Qubit effects E±^x = (I ± η σ_x)/2, E±^z = (I ± η σ_z)/2.
  Required: the canonical joint candidate R_{ab} = (I + aησ_x + bησ_z)/4 qualifies
  (all R_{ab} ≥ 0, Σ = I, marginals exact) iff 2η² ≤ 1; in particular qualification
  at η = 0.6 and rejection (`joint_candidate_not_positive`) at η = 0.8. Commutation
  of the pair is not required and does not substitute for the certificate.
- **N4 (tied control).** Two-qubit device whose only command drives G = X₁+X₂.
  Required: addressability certification of targets {X₁, X₂} fails with rank 1 < 2
  (`insufficient_command_rank`); the honest net (admitting only X₁+X₂) abstains in
  recovery — a mathematical decomposition of the span is not device capability.
- **N5 (erasure).** d = 8; site-local complete generator family vs its global
  Fourier conjugate. Required: identical *-closures (full M₈); recovery readings in
  non-matching orbits (overlap far below the 0.9 bar) — marking, not closure,
  carries the factorization.

## 11. No-gos and positive theorems (with falsifiers)

| id | statement (scope: finite dims, procedures of this spec) | falsifier |
|---|---|---|
| T1 | Global-algebra erasure no-go: the closed global algebra's isomorphism class carries no factorization; distinct markings can share one closure. | a map from closures to factorizations reproducing F2 on N5-type instances |
| T2 | Marked recovery (positive): F2 recovers Exact/Sectors/EquivClass/Abstain as specified; witness gates are necessary (removing them changes verdicts on N5/FT-12 instances). | an F2-qualified instance where the procedure's verdict is wrong (e.g. fabricated Exact) |
| T3 | Controller-free decomposition no-go: no map from (H_drift, state, global algebra) alone to accessible operations/factorization — four interfaces over identical (H, drift, ρ) must yield non-equivalent readings ([d1,d2,d3]-orbit-α / same-dims-orbit-β / coarser / Abstain). | a controller-free construction reproducing interface-dependent readings |
| T4 | Tied-control no-go: rank-deficient command Jacobians never certify independent addressability (N4). | a rank-1 device instance passing C2 for two targets |
| T5 | Glue theorem (positive): consistent covering atlases reproduce direct global recovery (verdict, dims, orbit). | a consistent atlas whose glued reading differs from direct recovery |
| T6 | Graded boundary: (a) no witness ⇒ MajoranaFrameOnly (O(2N)-invariance of CAR); (b) qualified charge witness ⇒ complex modes up to U(N). | (a) a CAR-functional selecting a pairing; (b) a qualified witness instance where J fails J²=−I or modes fail CAR/Σa†a=Q |
| T7 | Resource discipline: profiles are componentwise-monotone-invariant; scalarization can flip readings; transient (chain < 2) classes are never stable. | a scalarization provably neutral on all profiles, or a stable single-point class |
| T8 | Signed-covariance identity: (n̈⁺−n̈⁻)/(4ε) = ‖P_j h P_i‖_F² exactly (P1 contract). | any P1-contract instance violating the identity beyond numerical dust |
| T9 | Probe separation: quench ≠ signed covariance (N1 null). | a derivation of P1's identity for quench data within P1's contract |
| T10 | Response-law scope: number-conserving laws fail under pairing (N2: 1−Δ²); certificates must carry the charge witness. | a charge-witness-free extension reproducing N2 |

**Forbidden interpretations (apply to every claim above).** No statement here is a
bridge law, a statement about nature, or an empirical prediction; physical scope is
finite-dimensional models and laboratory-style interfaces. The registered
bridge-law ledger remains empty; natural-observation hits remain zero.

## 12. Adjudication order (fail-closed, total precedence)

For any single evaluation, the first applicable verdict from the top wins:

1. **OutOfDomain / scope error** — input outside the registered contract (wrong
   grading lane, non-Pauli in the GF(2) lane, non-antisymmetric in the quadratic
   lane, dense dim > 4096, probe-type mismatch, rejected noise model).
2. **Construction-time rejection** — role mixing, empty/incomplete contexts,
   certificate target mismatch, rank-deficient addressability, degenerate targets.
3. **InsufficientObservation** — missing channels/generators/witnesses
   (`insufficient_operational_generators`, `no_synthesis_path`, missing charge
   witness, `coverage_incomplete`).
4. **Straddled** — any decision interval crossing its frozen bar (commutator τ,
   cross-talk, glue overlap): abstain, never force.
5. **EquivalenceClassOnly** — multiple qualified non-matching orbits/atlases.
6. **SuperselectionSectors** — positive reading across a nontrivial center.
7. **Exact / GluedExact / ComplexModeFactorization** — full qualification.

Forced answers (emitting 5–7 when 1–4 applies) are scoring errors of the *selective
risk* type; abstentions where 5–7 is achievable are *answerable recall* errors. Both
are first-class scoring targets; abstention on a non-identifiable instance is the
**correct** output, not a failure.

## 13. Frozen bars (complete list)

| bar | value | used in |
|---|---|---|
| commutator threshold τ | 1e-3 | F2 graph |
| addressability σ_min | ≥ 0.5 | C2 |
| cross-talk | ≤ 0.1 (interval semantics) | C2 |
| synthesis residual | ≤ 1e-9 | C3 |
| CAR / witness / J / glue-match residual | ≤ 1e-9 | M1/M3/G2 |
| gauge-orbit matching | ≥ 0.9 (match) | F3 |
| stability chain (promotion) | length ≥ 2 | R3 |
| dense scope | dim ≤ 4096 | B3 |

## 14. What this spec does not close (honest boundary)

Holdout generation/scoring harnesses (HOLD-5..9), the external-replication campaign
layer, homology/VR readouts, Gaussian-oracle lanes, and everything tagged
`repository_replay_only` in the closure manifest are *not* closed by this paper.
Finite-data (shot-noise) certificate semantics — confidence *sets* over interfaces
replacing exact certificates — are deliberately out of OCS 1.0 scope: they are the
subject of the finite-data extension (fourth no-go and robust promotion), which
will be specified only after its theorems are proven. This spec freezes the exact
(noiseless-certificate) semantics that any finite-data extension must reduce to.

## Appendix A. Terminology map (repository ↔ spec; informative, not normative)

Repository module names (operational_net / laboratory_interface / resource_profile /
contextual_factorization / graded_recovery / structured_backend) implement §5,
§1–§4, §4, §6, §7, §9 respectively; the repository's verdict enums match the
vocabulary strings quoted in this spec verbatim. An independent implementation
needs none of these names — conformance is defined by verdict agreement on the
procedures and normative instances of this document.
