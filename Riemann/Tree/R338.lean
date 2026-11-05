/-!
# Round 338 — node N13 (Hilbert-Polya / spectral), round 5: no new proved spectral realization; Aksman's GRH proof (Zenodo 23839754) breaks at (79); 2609.25564's open question answered (U_E off the unit circle for E >= 4)

* Node: N13 Hilbert-Polya and spectral approaches (fifth round; R325 before).
* [check] arXiv (Oct-Nov 2026): no proved spectrum = zeta zeros; 2610.08437 (Berry-Keating oscillator view) is formal. Zenodo:
  23873342 (2 Nov), 23691133, 23603376, 23635192 are new but claim no complete proof (conditional or finite/local positivity).
* Refutation (orchestrator check of the PDF): Aksman, Zenodo 23839754 (31 Oct 2026, 'unconditional resolution' of GRH via
  Bianchi orbifolds), Thm 8.1: (76)-(78) give only <(B* + B + I)P0 f, P0 f> = 0 on the zero mode and a fixed nonzero bound
  e^{-pi sqrt 163} = 3.81e-18 on the other modes; the exact identity A*G + GA = 0 (79), used in (80) and (84), does not follow
  (R223 had left this step unexamined). REFUTATIONS.
* Harvest: 2606.06604 (conclusion l.1695-1697) leaves open B^{phi=p} under descent to R+*/p^Z and its link to elliptic curves;
  2609.00299 states no open question. Half-plane substitution in 2402.13082 (scout timed out): not examined.
* Exploration and computation: arXiv 2609.25564 (finite Rodriguez-Villegas approximants Z_E -> xi) leaves open whether the
  numerators U_E have unit-circle zeros. They do not for every even E = 4..40 (max ||q| - 1| = 0.056 at E = 4, 0.496 at E = 40;
  E = 4 checked independently: B_4 P has zeros 2.856 +- 6.108i). Z_E itself has all zeros on Re s = 1/2 for E <= 24, but so do
  the approximants of H_lambda for lambda = -1 (non-real zeros), and the smallest zero (7.9-9.3) is far from 14.13: the Hurwitz
  criterion 'Z_E real-rooted for all E' is not probeable ([RVApprox](Riemann/Deep/RVApprox/README.md)). Golyshev's canonical
  strip (Rodriguez-Villegas mechanism): no idea beyond this. Zeta functions of groups / Igusa local zeta: only finite-field
  (Deligne) line statements, no idea. Spin glasses (Guerra interpolation): no transfer to Weil's form (1706.08462 is only a
  random Euler product analogy), no idea.
* N1 bookkeeping this round: chekaudit SOUND (Chekmachev 1.925 / 2.131 even = N1 record, caveat: c = 45, 47, 71 on the
  author's code); weilooc's independent certificate of 1.6 committed (reviewer: correct; not a record).
-/

namespace Fayed.R338
end Fayed.R338
