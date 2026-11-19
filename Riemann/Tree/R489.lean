/-!
# Round 489 — node N5 (Nyman–Beurling / Báez-Duarte), round 16, Zenodo 19 Nov [check]: Munawer's paper M proves only the sufficiency direction (E_dir(N) → 0 ⇒ RH), the direction we machine-checked, while its text calls the matrix criterion 'strictly equivalent'; RiemannBlueprint v3.0.14 states Hardy's theorem in Lean without hypotheses; Skrynnik claims no proof; all eight scouts timed out under load

* Node: N5 (16th round; R476 before). Round type: [check] of the Zenodo records of 19 Nov (o44_zsweep.py from 2026-11-19, 9 records;
  texts pre-extracted in .scratch/r489/txt). Machine load ≈ 11 (five Orca agents): all eight scouts hit the 2-min limit; the checks
  below are the orchestrator's (lesson 228).
* [check] Munawer, Zenodo 24357107 (paper M, v3.0.0): MatrixFactorizationAndHalfObservation.lean proves
  riemannHypothesis_of_unitDefectMatrix_vanishing_unconditional (l.852: if inf{αᵀM_Nα : 1ᵀα = 1} → 0 then _root_.RiemannHypothesis)
  and the ground-state variant (l.2059), with M_N = unitDefectMatrix (l.639, entries unitDefectEntry, UnitPartitionHomogenization.lean
  l.50, built from Vasyunin's cotangent Gram kernel, paper l.139). Both are the sufficiency direction, as in our machine-checked
  RH_of_NBApprox (Riemann/Deep/NymanBeurling). The paper (l.382–384) calls RH 'strictly equivalent' to E_dir(N) = 1/(1ᵀM_N⁻¹1) → 0
  via its Theorems 2.1(3) and 5.2(3); the converse is not in the Lean file. No RH claim: KNOWN direction, nothing to refute.
* [check] Munawer, RiemannBlueprint v3.0.14 (Zenodo 24378371, 19 Nov 19:44): infinite_criticalLine_zeros
  (MATHLIB4_PR_STAGING_HardyCriticalLineZeros.lean l.221–223) states {t : ℝ | riemannZeta (criticalLine t) = 0}.Infinite with no
  hypotheses (via hardyZ sign changes); no sorry or axiom declaration in the Hardy files (orchestrator rg; not rebuilt): Hardy's 1914
  theorem as a hypothesis-free Lean statement (N8 side item; a formalization of a classical result). Paper L v3.0.25 (24378373):
  version of the R488 item.
* [check] Skrynnik (Zenodo 24356722, 24378143): the texts say the work 'can become a proof of the Riemann hypothesis; the missing step
  is …' (Fractal_Ladder_EN l.41): no claimed proof, nothing to refute. Koshelenko (24369317, mean square 1/π² of the error term of
  Jacobi-symbol sums): off the tree. Qin (24364735, claimed Goldbach proof) and Hackbarth (24359794): PDF downloads failed (file
  names with spaces/brackets) → next round if fetched.
* Exploration: Sinkhorn/Schrödinger bridges, Kolmogorov complexity, Markov-chain mixing: the three scouts timed out (untested).
* Scouts: 8, none answered (load).
-/

namespace Fayed.R489
end Fayed.R489
