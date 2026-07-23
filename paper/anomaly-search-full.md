# A Machine-Certified Exhaustive Classification of Minimal Anomaly-Free Chiral Spectra in Bounded SU(3)×SU(2)×U(1) Domains

**Draft v2 (v14.2 起草, v25.2 で改題・主張域限定).** 対象誌: SciPost Physics / CPC / PRD。骨子は [anomaly-search.md](anomaly-search.md)、全数値の一次ソースは `results/` と `certificates/`、主張の等級と限界は [claims.yml](../claims.yml)。**改題の記録 (PROMPT/7)**: 旧題 "The Standard Model as the Unique Minimal Anomaly-Free Chiral Spectrum" は無条件の一意性を示唆する — 一般の U(1) anomaly 方程式には無限の解族が知られており (Costa–Dobrescu–Fox 型の一般解構成、および近年の分類研究)、本稿の定理は「明示された有界領域・電荷量子化・多重度上限の下」でのみ成立する。定理の資格をタイトル段階で限定し、仮定を §0 に前置した。

---

## Abstract

We present an exhaustive search for anomaly-free chiral fermion spectra of SU(3)×SU(2)×U(1) gauge theories within explicitly bounded domains (representations up to the color octet and weak triplet, hypercharges |Y| ≤ 3, up to 8 multiplets), subject to the five perturbative anomaly-cancellation conditions, the Witten SU(2) anomaly, chirality, and charge under all factors. Within every domain scanned, the minimal solution is the 15-component Standard Model generation; it is unique at 15 components and isolated above them. A negative-control map identifies which conditions carry which conclusion: minimality is enforced by chirality, all-factor charge, the Witten anomaly, and the cubic SU(3) anomaly, while uniqueness is enforced by the cubic U(1) anomaly — the three linear anomaly conditions are redundant within the domain, itself a finding. The search is implemented three times with disjoint algorithms (depth-first, linear-elimination grid, meet-in-the-middle — algorithmic diversity, not independent replication; see §4), certified by SHA-256 canonical-form certificates, and the two core domains are verified as theorems in Lean 4 via native_decide with a poisoned-fuel design that makes completeness of the enumeration part of the theorem statement. Extending the gauge group: the only rank-2 chiral extension in the window is B−L (requiring ν_R, unique as a charge plane), no rank-3 chiral extension exists (0 solutions against 355 in the vectorlike control), and among exceptional groups only E₆ admits chiral matter — its 27 reducing exactly to one SM generation plus vectorlike remainder. Consistency alone thus counts the forces and fixes the matter, within windows whose boundaries we state precisely. All results are reproducible from a dependency-free Rust repository with a machine-validated claim ledger.

---

## 0. Scope and standing assumptions (theorem preconditions)

Every "minimality", "uniqueness", and "isolation" statement in this paper is a theorem **under the following explicitly stated hypotheses**, and nothing beyond them is claimed:

1. **Bounded search domain.** Color representations in {1, 3, 3̄, 6, 6̄, 8}, weak representations in {1, 2, 3}, hypercharges |y| = |6Y| ≤ 18 (i.e. |Y| ≤ 3), at most 8 multiplets (domain-dependent; each domain's exact bounds are machine-readable in `certificates/v62_domains.json`). No claim is made outside these windows; the EFT reading (heavy vectorlike matter decouples) is an assumption, not a result of the scan.
2. **Charge quantization.** Hypercharges are integers in units of Y/6 after gcd rescaling; the classification is of *rational* charge assignments up to normalization. Irrational or unquantized charges are out of scope. Consequently, the known infinite families of anomaly-free U(1) charge assignments — the general solution of the U(1)³ + gravitational system [6] and subsequent classification studies — coexist with our results: they live outside the bounded, all-factor-charged, chiral window classified here. Our uniqueness claims are **classifications within bounds**, not global uniqueness theorems.
3. **Multiplicity bounds and multiset semantics.** Spectra are multisets of at most k multiplets; repeated multiplets are allowed up to the component cap. Uniqueness is stated up to the canonical equivalences (permutation, global conjugation, U(1) sign, gcd rescaling).
4. **Vectorlike pairs.** The chirality condition C excludes multiplets appearing together with their conjugates. All statements therefore concern *chiral content*: any spectrum may be extended by arbitrary vectorlike pairs without affecting anomaly cancellation, and the negative-control run with C removed (§3) quantifies exactly how minimality and uniqueness collapse when they are admitted.
5. **Sterile states.** States neutral under all three factors ("all-factor charge" condition F) are excluded from the baseline and reintroduced explicitly on the ν_R robustness axis; the gravitational anomaly A3 is imposed in both variants.
6. **Observational vs. mathematical exclusion.** The eight 16-component neighbours of the SM are excluded by *observation* (fractionally charged hadrons etc.), not by group theory; we mark every such exclusion as observational and keep it out of the theorem statements.
7. **Trust base.** Rust exhaustive enumeration (exact integer arithmetic) certified by SHA-256 canonical-form certificates; two core domains re-proved in Lean 4 (`native_decide`; kernel + native evaluator trusted). The three Rust implementations are *algorithmically diverse but not independently replicated* — see §4.

## 1. Introduction

Why this matter content? The Standard Model (SM) generation — fifteen Weyl components in five multiplets with hypercharges (1,1)₋₆ ⊕ (1,2)₃ ⊕ (3,1)₋₂ ⊕ (3,1)₄ ⊕ (3̄,2)₋₁ (in units of Y/6) — satisfies a tightly coupled system of anomaly-cancellation conditions. That the SM is anomaly-free is textbook material; that hypercharge is *forced* by anomaly cancellation under stated assumptions has a long literature. This paper asks and answers a computational version of the strongest form of the question:

> Within an explicitly bounded space of chiral gauge theories, is the SM generation the **minimal** consistent spectrum? Is it **unique** at its size? Is it **isolated** above it? And which consistency conditions carry each of these properties?

Our answers: yes, yes, yes — with the entire component-count spectrum {15: 1, 16: 8, 17: 1, 18: 18, 22: 2, 24: 459} enumerated in the widest domain, and with a *control map* attributing minimality and uniqueness to specific conditions. The method is exhaustive integer enumeration: no sampling, no heuristics, exact arithmetic throughout, with completeness itself made checkable (§4).

Three methodological features distinguish this from a table of known facts. First, *independent triple implementation*: the same domains are enumerated by three algorithms sharing no code path, and their canonical-form solution sets agree byte-for-byte (SHA-256 certificates in the repository). Second, *negative controls*: each consistency condition is switched off in turn, and the search is re-run to see what breaks — the instrument is demonstrated to be able to fail. Third, *machine verification*: the two core domains are re-proved inside Lean 4, with a fuel-poisoning construction that converts "the enumeration did not silently skip a branch" into part of the theorem.

## 2. Setup

A candidate spectrum is a multiset of at most k multiplets (R_c, R_w)_{y} with R_c ∈ {1, 3, 3̄, 6, 6̄, 8}, R_w ∈ {1, 2, 3}, and integer y = 6Y in a stated window. The consistency conditions are:

- A1–A3 (linear): [SU(3)]²U(1), [SU(2)]²U(1), grav²U(1);
- A4: [U(1)]³; A5: [SU(3)]³;
- W: Witten's SU(2) global anomaly (even number of weak doublets, counted with color dimension);
- C: chirality (no multiplet together with its conjugate; no completely real content);
- F: charged under all three factors (some color, some weak, some hypercharge non-trivially present).

Domains (machine-readable in `certificates/v62_domains.json`): the core domain D₁ (small representations, |Y| ≤ 3/2, ≤ 5 multiplets, ≤ 15 components), the extension domain D₂ (≤ 6 multiplets, ≤ 24 components), and robustness axes: large representations (6, 6̄, 8; weak triplets), |Y| ≤ 2 and ≤ 3, up to 8 multiplets, and an added right-handed neutrino. The physical reading of the window is an EFT assumption — heavy vectorlike matter decouples (verified separately at M^{−1.91} in this program) — so the theorems are of the form "within this window," stated without apology and without extrapolation.

## 3. Results: minimality, uniqueness, isolation

**Minimality and uniqueness.** In every domain scanned, the minimal chiral solution has 15 components and is exactly the SM generation, unique up to the canonical equivalences (permutation, global conjugation, U(1) sign, and charge rescaling). This holds with large representations admitted, with |Y| up to 3, with up to 8 multiplets, and with ν_R present (where the 16-component SM+ν_R appears as expected at the next size).

**The spectrum.** The widest domain yields component counts {15: 1, 16: 8, 17: 1, 18: 18, 22: 2, 24: 459}: one solution at 15; eight 16-component neighbors, each excluded empirically by fractionally charged hadrons or similar observational inputs (this exclusion is observational, not group-theoretic — we mark it as such); a single 17-component solution containing a weak triplet; and nothing between 18 and 22. The SM is not merely minimal; it sits in a *gap*.

**The control map.** Re-running with individual conditions removed:

| condition removed | minimal size | # at 15 | reading |
|---|---|---|---|
| chirality | 10 | 59 | minimality *and* uniqueness collapse (vectorlike pairs flood in) |
| all-factor charge | 5 | 1 | minimality collapses (pure-U(1) towers) |
| Witten SU(2) | 12 | 1 | minimality collapses |
| [SU(3)]³ | 12 | 6 | both degrade |
| [SU(3)]²U(1) | 15 | 1 | **redundant in-domain** |
| [SU(2)]²U(1) | 15 | 1 | **redundant in-domain** |
| grav²U(1) | 15 | 1 | **redundant in-domain** |
| [U(1)]³ | 15 | 9 | uniqueness collapses (nine 15-component solutions) |

The linear anomaly conditions, so prominent in textbook derivations, are *implied* by the remaining conditions within the domain; minimality is carried by chirality + all-factor charge + Witten + cubic SU(3), and uniqueness by the cubic U(1). The baseline solution set is contained in every control's solution set (monotonicity check), so the controls genuinely relax rather than shift the problem.

## 4. Methods: three implementations, certificates, machine verification

**Triple implementation.** (i) Depth-first with pruning over multiplet stacks; (ii) linear-elimination over a hypercharge grid; (iii) meet-in-the-middle over half-assignments with exact integer keys (the anomaly conditions split additively). All arithmetic is exact (i64/i128); canonical forms quotient permutations, conjugation, U(1) flips, and gcd rescaling. The three solution sets agree exactly in every domain; their canonicalized listings are certified by SHA-256 (`certificates/`), and CI recomputes and compares the hashes.

*Honesty note on independence.* These three implementations are **algorithmically diverse, not independently replicated**: they were produced within a single research program, by the same author/toolchain lineage, and a shared conceptual misunderstanding could in principle propagate into all three (and into the Lean encoding of the domain). What the triple agreement rules out is implementation-level error, not specification-level error. The domain specification and certificates exist precisely to make third-party replication cheap: reimplementing any domain from `certificates/v62_domains.json` alone and comparing SHA-256 canonical listings is an afternoon's work, and we solicit exactly that as the missing verification layer.

**Machine verification in Lean 4.** The core domains are re-proved as theorems of the form `scanResult = 4` (the number of canonical-orbit representatives, SM's orbit): the enumeration is a structurally recursive walk whose fuel is *poisoned* — running out of fuel contributes +10⁹ and an out-of-orbit hit contributes +10⁶ to the result — so the equality theorem simultaneously asserts correctness *and completeness* of the walk. `native_decide` evaluates the closed decision procedure (991 s for the small-representation domain; 9,734 s for the large-representation domain, 1.7×10¹⁰ recursive calls). The trust base is Lean's kernel plus the native evaluator, which we state explicitly; the Rust and Lean enumerations are independent codings of the same domain definition.

**Runtime prediction.** Lean native_decide wall time is proportional to the total number of recursive calls (~0.6 µs/call), measured by instrumenting the identical recursion in Rust first. We recommend this practice: enumeration-completeness proofs are affordable exactly when their cost can be predicted before writing them.

## 5. Climbing the gauge group: the U(1) staircase and the exceptional ceiling

**Second U(1).** With the SM matter fixed and a second hypercharge Y′ scanned over |6Y′| ≤ 12: without ν_R every solution is proportional to Y; with ν_R the solutions span exactly the {Y, B−L} plane. A structure-free rank-2 scan (both charges free, ≤ 5–6 multiplets) finds *no* chiral rank-2 solution below 16 components, and at 16 exactly one charge *plane* — {Y, B−L} on SM+ν_R — after Plücker-coordinate classification of planes (GL(2,ℤ) lattice classification overcounts bases of the same plane; we document this initially-wrong equivalence and its fix). **The minimal way to add a second force is B−L, and it forces the right-handed neutrino.**

**Third U(1).** A rank-3 scan (19 exact conditions: 9 linear, 10 cubic monomials) over the stated window finds **zero** chiral solutions; the vectorlike control finds 355, so the instrument is not blind (claim QRN-GAUGE-014). The staircase ends at two steps: Y, then B−L, then nothing.

**Exceptional groups.** G₂, F₄, E₇, E₈ admit only self-conjugate representations; we verify W = −W as exact integer weight-multiset identities for 7, 14, 26, 52, 56, 133, 248 (with Weyl-closure checks of each construction, and SU(3)'s complex 3 as the negative control), so their chiral cores are empty. Among exceptional groups only E₆ has complex representations, and the chiral core of its 27 — after stripping conjugate pairs and singlets exactly — is one SM generation (claims QRN-GAUGE-015, v8.2). The staircase and the exceptional ceiling point at the same spectrum.

## 6. Limitations

The windows are finite and stated: large-|Y| rank-3 territory, larger multiplet counts, and huge representations remain open; the 16-component neighbors are excluded by observational inputs, not consistency; the exceptional-group statement verifies representative low-dimensional representations and cites the standard self-conjugacy theorems for the universal part; the v4.3-domain Lean formalization at ~10¹¹ calls awaits an array-based enumeration design. None of these limitations, if lifted, can *remove* solutions — every extension so far has only confirmed the gap around the SM.

## 7. Reproducibility

Dependency-free Rust (std only), fixed seeds, exact integer arithmetic in all searches, [PASS]/[FAIL] verification embedded in every program (314 PASS / 0 FAIL across the repository at time of writing), SHA-256 certificates for all solution sets, machine-validated claim ledger (94 claims with grades, evidence files, and limitations), and CI that rebuilds, re-runs the audit suite, re-validates certificates, and re-checks both Lean theorems.

## References

*(書誌は 2026-07-05 に Web 照合済み。[1] は正典として据え置き)*

[1] E. Witten, "An SU(2) anomaly," Phys. Lett. B 117 (1982) 324.
[2] C. Q. Geng and R. E. Marshak, "Uniqueness of quark and lepton representations in the standard model from the anomalies viewpoint," Phys. Rev. D 39 (1989) 693; J. A. Minahan, P. Ramond and R. C. Warner, "A comment on anomaly cancellation in the standard model," Phys. Rev. D 41 (1990) 715.
[3] K. S. Babu and R. N. Mohapatra, "Is there a connection between quantization of electric charge and a Majorana neutrino?," Phys. Rev. Lett. 63 (1989) 938.
[4] R. Foot, G. C. Joshi, H. Lew and R. R. Volkas, "Charge quantization in the standard model and some of its extensions," Mod. Phys. Lett. A 5 (1990) 2721.
[5] B. C. Allanach, B. Gripaios and J. Tooby-Smith, "Anomaly cancellation with an extra gauge boson," Phys. Rev. Lett. 125 (2020) 161601.
[6] D. B. Costa, B. A. Dobrescu and P. J. Fox, "General solution to the U(1) anomaly equations," Phys. Rev. Lett. 123 (2019) 151601, arXiv:1905.13729.
[7] The Lean 4 theorem prover (the trust base of native_decide is stated explicitly in §4).
[8] 本計画のリポジトリ: Quantum Relational Network (claims.yml, results/, certificates/, proofs/) — 全数値・全証明書の一次ソース。
