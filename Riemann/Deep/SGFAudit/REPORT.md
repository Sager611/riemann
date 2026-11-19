AUDIT PASSED — every claim of small-gap-floor holds as written (Lemma 1, Theorems A and B, both Arb certificates re-verified independently at 256 bits, Lean statements faithful); one routine step is missing for CGdL's broader test class (added, reviewer CORRECT), plus four wording fixes.

# SGFAudit: line-by-line audit of discoveries/small-gap-floor.md (sgfaudit, 19 Nov 2026)

Scope: the page's Claim, Method, Certificate and Lean paragraphs, against Lamzouri 2611.14486, Chirre–Gonçalves–de Laat (CGdL)
1810.08843 and Durkan–Pearce-Crump (DPC) 2610.28527 (TeX read at the cited lines), Riemann/Deep/ExploreDPP/ (DERIVATION.md,
verify_dual.py, logs/, Lean/). Full chain with source lines: [FULLPROOF.md](FULLPROOF.md). Not rebuilt: Lean (brief: no Lean).

## Items (verdict, evidence)
| # | Claim item (page) | Verdict | Evidence |
|---|---|---|---|
| 1 | D(λ) definition | holds | 2611.14486 l.183, multiplicity l.177 |
| 2 | Lamzouri's criterion as quoted | holds | Thm Main2 l.392–397 verbatim (plus 0 < λ ≤ 1, harmless); FT e^{−2πiξu} l.307–309 |
| 3 | never proves D(λ) > 0 for λ ≤ 0.775, any θ < 4/7 | holds | Lemma 1 line by line (FULLPROOF §3); every freedom of Main2 covered (§2); CertB PASS twice |
| 4 | θ = 1 LP never reaches λ ≤ 0.595 | holds | CertB1 PASS twice; is CGdL's functional (l.527–562) |
| 4a | ... "the method of CGdL" | holds after a step | CGdL's class A_LP (l.420–430: continuous L¹, ĝ ≥ 0) is wider than C_c²/Schwartz; mollification closes it (§3.4) |
| 5 | 0.7851 within 0.0101 of the floor | holds | 0.7851 − 0.775; upper end SmallGap (float sanity §8) |
| 6 | float floor ≈ 0.784 | holds (float) | DERIVATION λ*(X) fit; 0.78 float-feasible, certification fails (verify_t47exact_078.log) |
| 7 | θ = 1: floor ≥ 0.595, CGdL best 0.6039 | holds | CGdL Thm 4 l.317–325 (RH + N* ~ N, pairs N(λ, T)) |
| 8 | DPC: D(λ) > 0 ∀λ > 1/2 under RH | holds | 2610.28527 Thm thm:positive-proportion l.61–69; log γ_n vs log T conversion (FULLPROOF I6) |
| 9 | RH case not capped by the functional | holds | follows from 4 and 8 |
| 10 | Theorem A (1/(1+θ); 7/11, 1/2) | holds | §3.1; Lean theoremA' hypothesis-free; series = triangle form to 10⁻⁶ |
| 11 | Lean sentence | holds | statements match (§5); recorded axioms propext, Classical.choice, Quot.sound; no sorry/axiom/native_decide |

## Certificate audit (FULLPROOF §4)
- verify_dual.py: exact rationals, Arb outward rounding, closed cells covering [0, 1/(2h)] from α = 0, θ a grid point with the
  inclusive left branch, valid Taylor minorant, M2 ≥ |P''|. Correct.
- Independent re-verification (indep_verify.py, own method: per-term sin_cos_pi of exact rationals, centred third-order adaptive
  cells, 256 bits) on the Lean lists certB/certB1 (identical to the JSON after rounding, cmp_cert.log): both PASS; true minima of G
  5.76·10⁻⁴ (θ = 4/7, α ≈ 1.2551) and 6.86·10⁻⁴ (θ = 1, α ≈ 1.9539). Their worst-cell bounds 0.000204450606588984… and
  0.000238928568865544… reproduced to 40 digits; Horner vs direct sums agree to 3.5·10⁻⁷⁴ at 56 random/edge points.

## Falsification (FULLPROOF §6): no counterexample
Poisson identity with a C² B-spline whose support contains m/h (10⁻¹⁴); Lamzouri's own 0.7892 kernel against CertB (identity to
10⁻¹⁶, its positive ∫f^dσ comes from the atom at exactly u = 0.775 where f^ = +0.0094); SmallGap's 0.7851 kernel has f^(0.775) > 0;
θ = 1/(2h) boundary (a = 1/2) fine; f ≥ 0 is load-bearing and is a hypothesis of both Lamzouri and CGdL.

## Corrections (none changes a number or a conclusion)
1. Claim, "the method of Chirre–Gonçalves–de Laat": CGdL's class A_LP (1810.08843 l.420–430) is only continuous L¹ with ĝ ≥ 0; the page,
   DERIVATION §1 and Lean (Admissible, AdmissibleDecay) cover C_c² and O(|x|⁻²)-decaying kernels. Add: "for all of A_LP by Gaussian
   mollification (SGFAudit FULLPROOF §3.4)". Reviewer: CORRECT.
2. Claim, "Under RH the θ = 1 functional stops at 0.595": read "cannot reach λ ≤ 0.595; its floor lies in [0.595, 0.6039]". Also
   CGdL's 0.6039 is for pairs N(λ, T) under RH and N* ~ N (l.317–325), as the Review paragraph already says.
3. Claim, "the unconditional record λ ≥ 0.7851": dated; the record is now λ > 7/11 (RECORDS, other methods). Read "smallgap's
   0.7851 (the method's best)".
4. verify_dual.py takes θ from the JSON float 0.5714285714285714 (< 4/7) unless argv[3] = 4/7 is passed; the logged run and the
   rerun command (ExploreDPP/REPORT.md l.57) pass it. DERIVATION §3 "K = 235 (η on (0, 26])": 235h = 26.018. Cosmetic.
Not audited: the Before and Novelty paragraphs beyond CGdL (Lagarias–Rodgers, CCCM, Inoue–Kobayashi–Toma), SmallGap's exact
certificate for the upper end 0.7851 (float sanity only, FULLPROOF §8), DPC's Thm B proof (dpcaudit).

## What survives
Everything: Theorem A (floor ≥ 7/11 for θ < 4/7, 1/2 at θ = 1, hypothesis-free in Lean), Theorem B (floors ≥ 0.775 and ≥ 0.595,
certified twice by different methods), the floor intervals [0.775, 0.7851] and [0.595, 0.6039].

## Dependents (rg -n 'small-gap-floor', plus 0.775 mentions)
DISCOVERIES.md l.14; RECORDS.md l.25 ("Method floor"); TREE.md l.12 (N8 Barriers); discoveries/small-gap-hardcore.md l.5;
Riemann/Tree/R442.lean (the round), R448.lean l.4 (escape routes); Riemann/Deep/DPCAudit/REPORT.md l.56; Riemann/Deep/DPCUncond/
REPORT.md l.25; NOGO.md l.615 (explore34), l.633 (explore40). None needs a change beyond corrections 1–3 where they repeat the wording.

## Files (Riemann/Deep/SGFAudit/)
FULLPROOF.md (chain with source lines); indep_verify.py + indep_verify.log (CertB, CertB1 at 256 bits); cmp_cert.py + cmp_cert.log
(Lean lists = JSON); falsify.py + falsify.log (F1–F6); smallgap_sanity.py + smallgap_sanity.log. Runs under .scratch/memcap.sh 290:
peaks 15 MB (indep_verify, footprint 14.8 MB, 61 s), 21 MB (falsify; first attempt killed at 382 MB by numpy's leggauss(4000)
eigenproblem, rewritten with a 120-node composite rule), 14 MB (smallgap_sanity). One process at a time, single-threaded.

RH relevance: N8
