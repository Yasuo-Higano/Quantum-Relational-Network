# One-page challenge: independently test a frozen topology-readout protocol (D2-R)

**TL;DR** — We publish a frozen, machine-checkable protocol that claims:
*given a declared node factorization*, the spatial adjacency of a quadratic
fermionic network can be read from signed density-response curvature, with a
fail-closed verdict semantics (the reader must answer "cannot be determined"
when the data does not qualify). We are looking for **one** independent team to
implement it from the paper-closed spec alone and report the result —
**success or failure, both are publishable outcomes for us**.

## What exactly is claimed (and what is not)

- Claimed: an exact response law (curvature of signed probes reads
  `|h_ij|^2`-type block weights), a set of no-go theorems (what canNOT be
  promoted from the declared observation contract), and a fail-closed
  finite-data certificate semantics — all under a **given factorization**.
- NOT claimed: discovery of network topology without declared structure,
  emergent geometry, gravity, or any natural-world observation. We also do not
  claim novelty for "inferring topology from local measurements" per se —
  the specific frozen response law + no-go + fail-closed semantics is the unit
  under test.

## What you get

- `paper/operational-core-spec.md` (OCS-1.0) — a paper-closed spec: no source
  code, no recorded output values needed. sha256-pinned.
- `reproducer/campaigns/d2r-v1/` — frozen preregistration schema, ambiguity
  ledger, valid/invalid report fixtures. sha256-pinned
  (`MANIFEST.sha256`, pinned by
  `af967b3eb9a34511bd93785c05b4b6ca1899cd9459f0f9e2680219767e22de28`).
- We do **not** hand over our numerical kernel, implementation skeleton, or
  translation code — a clean-room implementation is the point.

## What we ask

1. Preregister (JSON conforming to `PREREGISTRATION.schema.json`) **before**
   running.
2. Implement from the spec in any language; run the frozen protocol.
3. Report via GitHub issue/PR — including failures. One valid FULL report
   suffices (frozen promise: we will not move the goalposts to two).
4. Questions are answered only as public, non-normative clarifications in the
   ambiguity ledger (the frozen text does not change).

- **QUAL** (practice tier): reproduce the fixture calibration — good for a
  weekend evaluation, does not count as external replication.
- **FULL**: the end-to-end unit — counts as the first external replication if
  the six independence conditions hold (different author/org, independent
  repository, no shared numerical kernel, protocol frozen before run, commit
  hash recorded, results public including failures).

## Honest status

External replications to date: **0**. Recorded experimental runs: **0**.
This challenge is exactly our attempt to change that, one report at a time.

Repository: https://github.com/Yasuo-Higano/Quantum-Relational-Network
