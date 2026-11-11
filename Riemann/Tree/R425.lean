/-!
# Round 425 — node N8 (critical-line zeros), round 12: the Lean-kernel-checked arena bound κ ≈ 0.673638877 for distinct zeros on the line (Bei Yi, Riemann.fail, Zenodo 24097425) authenticated: frozen zeta23 statement with classical inputs only, standard axioms, Lean and nanoda acceptance; it beats 0.6735382 for zeros on the line, not for simple ones; no newer N8 item; quantum defect theory, well testing and power grids give no route

* Node: N8 (twelfth round; R421 before; gaplambda and smallgap agents active on gaps). Round type: [check] a record found by
  harvest6 (arena acceptance 9 Nov), added to the N8 line as [check] the same day.
* [check] Statement: github.com/josusanmartin/riemann, challenge/templates/Challenge.Candidate.lean sections 12-20 (and CONTRIBUTING.md):
  dyadic and cumulative (kappa - eps) Ncount <= N0star for the zeros of zeta, N0star = distinct zeros on the critical line,
  Ncount = all nontrivial zeros with multiplicity; the only hypotheses are zeta23's PaperInputs: EF (Weil's explicit formula),
  RvM (Riemann-von Mangoldt local count), cheb (Chebyshev-Mertens), MV (Montgomery-Vaughan Hilbert inequality), Gamma
  (Stirling facts) (our copy: Riemann/Deep/WeilEF/Z23/Hypotheses.lean l.149-160); no RH and no candidate hypotheses: HOLDS.
* [check] Submission: Zenodo 24097425 v3.0 (T022, d250; 5,447 points, 28,421 leaves; Solution_T022.lean 1,983,753 bytes):
  AXIOMS_T022.txt lists exactly propext, Classical.choice, Quot.sound; VERIFICATION_T022.json verifier exit 0; attestation
  'kernel-verified' (Lean and nanoda, 9 Nov 2026) with the exact numerator and denominator; accepted on the arena: HOLDS.
* Comparison: kappa = 0.67363887666504 exceeds 1669159/2478195 = 0.6735382 (2610.08965) by 1.0068e-4 for distinct zeros on
  the line (hence for zeros on the line), but not for simple zeros on the line, where 2610.08965 stands; the method is zeta23's
  Hilbert-space/pair-correlation bound with a ten-point (TEN) certificate refined by root coverage, not a mollifier, which is
  why it passes the basic method's 2 - C_MT = 0.672501. TREE N8 line updated (no longer [check]).
* zeta23 trust boundary: all five PaperInputs fields are classical published inputs (Weil 1952 / Iwaniec-Kowalski, Titchmarsh,
  Montgomery-Vaughan 2007, Montgomery-Vaughan 1974); none is a 2025-2026 preprint: HOLDS.
* [check] Newest N8 items: arXiv API (proportion on the line, positive proportion) nothing dated 1 Nov or later; Zenodo since
  6 Nov: 24104077, 24101490, 24097425, 24041188: no other record: HOLDS.
* Exploration: multichannel quantum defect theory (unitary K-matrix, Lu-Fano eigenphase interlacing: no channel structure
  tied to primes), pressure-transient well testing (Hankel moment positivity of a relaxation measure holds for every positive
  measure, so the prime-power measure's Hankel minors carry no RH information; Weil's form is not a moment positivity) and
  power grids (Kron reduction and swing-equation passivity are generic Laplacian positivity): no route.
-/

namespace Fayed.R425
end Fayed.R425
