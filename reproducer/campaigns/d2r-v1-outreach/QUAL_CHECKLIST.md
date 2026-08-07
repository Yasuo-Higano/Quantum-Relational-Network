# D2R-QUAL checklist (practice tier — not an external replication)

QUAL is a low-burden qualification pass: you verify that you can read the
frozen materials and reproduce the *calibration fixtures*, without committing
to the full end-to-end unit. Completing QUAL does not move any replication
counter; it exists to de-risk a later FULL attempt.

- [ ] Read `reproducer/SPEC.md` (general rules; failures are submitted in the
      same format as successes).
- [ ] Read `paper/operational-core-spec.md` (OCS-1.0) — confirm the sha256 pin
      matches the one recorded in the repository.
- [ ] Read `reproducer/campaigns/d2r-v1/REPLICATION_NOTE.md` (frozen scope:
      what one valid report does and does not unlock).
- [ ] Validate `INVALID_REPORT_FIXTURES/prereg_valid_minimal.json` as
      **conforming** against `PREREGISTRATION.schema.json` with your own
      tooling.
- [ ] Confirm the two `prereg_invalid_*.json` fixtures are **rejected** by
      your tooling (missing independence declaration / capability inflation).
- [ ] Check `AMBIGUITIES.yml` for known clarifications before asking new
      questions; new questions go to a public issue (answers are appended as
      non-normative clarifications only).
- [ ] Decide: proceed to FULL (preregister first), or stop here (please tell
      us why — "too much burden" is a valid, recordable outcome that helps us
      design D2R-v2).
