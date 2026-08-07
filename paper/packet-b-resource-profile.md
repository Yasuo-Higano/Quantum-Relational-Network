# Packet B — Set-Valued Resource Profile (adversary packet, v35.4)

- packet: B (Resource Profile)
- discipline: paper-closed — no source code, no recorded output values.
  Claims carry: status / domain / claim (exact quantifiers) / falsifier /
  forbidden interpretation.
- adversary task: independent derivation, counterexample search, scope
  narrowing. `RefutedAsStated` and `Inconclusive` are first-class outcomes.
- expected_outputs_in_packet: false.

## Domain shared by all claims

A parameter space Θ; a poset of resource budgets. For each budget b:
an equivalence relation ~_b on Θ (identifiable granularity), with quotient map
q_b : Θ → Θ/~_b; for b ≤ b' the finer relation refines the coarser
(~_{b'} ⊆ ~_b), equivalently q_b factors as q_b = r_{b'b} ∘ q_{b'}.
From data D, each budget produces a confidence set C_b(D) ⊆ Θ with marginal
coverage ≥ 1 − α. The set-valued reading is Q_b(D) = q_b(C_b(D)).

## RP-1 — Non-functoriality no-go (theorem candidate)

- status: theorem candidate.
- claim: coverage alone does not imply samplewise nesting, and without
  nesting the profile b ↦ Q_b(D) fails naturality: there exist coverage-valid
  families with r_{b'b}(Q_{b'}(D)) ⊄ Q_b(D) for some realization D. A minimal
  counterexample exists with |Θ| = 2, two-point data, identity refinement,
  and both budget families exactly at their nominal coverage level.
- quantifier structure: ∃ (families, D) violating the containment — the no-go
  refutes "for all coverage-valid families, naturality holds".
- falsifier: a proof that marginal coverage implies samplewise nesting for
  arbitrary (or any natural nontrivial class of) confidence-set families.
- forbidden interpretation: this does not say refinement is impossible —
  it says refinement is a property of the CONSTRUCTION (see RP-2), not of
  coverage. It also does not claim nesting violations are rare or frequent;
  rates are an empirical matter outside this packet.

## RP-2 — Conditional lax refinement (theorem candidate)

- status: theorem candidate.
- claim: if q_b = r ∘ q_{b'} (refinement of equivalences) and the confidence
  sets are samplewise nested (∀θ: θ ∈ C_{b'}(D) ⇒ θ ∈ C_b(D)), then
  r(Q_{b'}(D)) ⊆ Q_b(D). Pure logic: no finiteness, measurability, or
  coverage assumptions are needed for this implication.
- construction note: samplewise nesting can be bought, e.g. by simultaneous
  error allocation plus running intersection
  C̃_{b'}(D) = C_{b'}(D) ∩ C̃_b(D); the coverage of the intersected family is
  controlled by a union bound over the allocated levels.
- adversary question (RP-2q): the packet claims that restriction to a common
  master dataset alone (cumulative data, per-budget intervals at full level,
  no intersection) does NOT guarantee nesting. Verify or refute with your own
  interval construction.
- falsifier: a counterexample to the implication itself (nesting + refinement
  yet containment fails).
- forbidden interpretation: "functor" language for arbitrary profiles; the
  result is a conditional (lax) statement.

## RP-3 — Margin stability (theorem candidate)

- status: theorem candidate.
- claim: for interval readings [lo, hi] (lo ≤ hi) against a decision boundary
  τ with the three-valued verdict (edge / no-edge / straddled), any endpoint
  perturbation strictly smaller than the verdict's margin leaves the verdict
  unchanged, where the margin is: lo − τ for edge; the distance by which hi
  clears τ for no-edge; the smaller of the two boundary distances for
  straddled.
- falsifier: a sub-margin perturbation that flips a verdict.
- forbidden interpretation: stability without a margin premise (see RP-4).

## RP-4 — No global Lipschitz stability (theorem candidate)

- status: theorem candidate.
- claim: there is no constant L such that the discrete verdict distance is
  bounded by L times the endpoint perturbation uniformly over all intervals:
  at the boundary, pairs with arbitrarily small perturbation and verdict
  distance one exist at every scale.
- falsifier: a global Lipschitz bound valid on all of interval space.
- forbidden interpretation: this does not contradict RP-3 (margin-conditional
  stability); it shows the margin premise cannot be dropped. Non-answers
  (straddled) near boundaries are semantics, not failure.

## RP-5 — Chain-of-length-two promotion as a corollary (theorem candidate)

- status: theorem candidate (derivation of a previously frozen operational
  rule from RP-1/RP-2).
- claim: (a) soundness — under refinement + nesting, if the finer budget's
  reading is a singleton {v'}, then r(v') belongs to the coarser budget's
  reading set: promoting a verdict that persists across two consecutive
  budgets can never contradict the coarser reading. (b) necessity of the
  two-chain premise — a singleton reading at one budget can evaporate at the
  next budget under coverage-valid independent sampling; a minimal
  two-hypothesis counterexample exists. Together these derive the rule
  "stable iff it persists over a chain of length ≥ 2; single-point readings
  are transient and must not be promoted".
- falsifier: (a) a nested-refinement instance where a persistent singleton
  contradicts the coarser set; or (b) a proof that singleton readings are
  automatically persistent under coverage alone.
- forbidden interpretation: the corollary licenses promotion only along
  chains where nesting was constructed; it does not license promotion across
  independently-sampled budgets.

## RP-6 — Adversary question: alternative nesting constructions (open)

- status: open design question (input to a future finite-data specification).
- question: besides running intersection with union-bound allocation, which
  constructions give samplewise nesting with better width/coverage trade-offs
  in this setting (e.g. anytime-valid confidence sequences, simultaneous
  allocation over the whole budget chain)? Deliverables: a construction with
  a stated coverage guarantee and a nesting proof, or `Inconclusive`.
- forbidden interpretation: none (design space question).
