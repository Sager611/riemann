/-!
# Round 354 — node N9 (moments, Lindelof, subconvexity), round 6: records stand; conductor-aspect mu line from the half-plane and Petrow-Young (routine); Bluteau's local-Lindelof AFE is false

* Node: N9 moments, Lindelof, subconvexity (sixth round; R343 before).
* [check] state of the art: no improvement of 13/84, of the explicit 0.611 t^(1/6) log t (2602.05614) or of the exponent pairs and
  moment records; 2606.27323 is a sixth-moment LOWER bound (34.4 c_3 T log^9 T); Zenodo: nothing newer. HOLDS.
* Adversarial: 2602.05614 Thm 1 at small t: max |zeta(1/2+it)|/(t^(1/6) log t) = 0.5074 at t = 45.6 on [3, 300] (mpmath, step 0.01 to 30,
  0.1 beyond; r(3) = 0.4082), below 0.611: HOLDS.
* Computation (routine, COROLLARIES): the Lean-checked half-plane for every Dirichlet L-function gives Lindelof in the conductor
  aspect for Re s > B (Titchmarsh 14.2, uniform in q(|t|+1)); Rademacher's Phragmen-Lindelof with Petrow-Young's Weyl bound for all
  primitive chi (1908.10346 Thm 1.1, TeX l.193-199) gives L(sigma+it, chi) << (q(|t|+1))^((B-sigma)/(6(B-1/2))+eps) on [1/2, B], below
  Weyl + convexity (1-sigma)/3 and Burgess 3(1-sigma)/8 on all of (1/2, 1).
* Adversarial (Zenodo): Bluteau 17575894 (Nov 2025, local Lindelof-type bounds): the opening short AFE has the wrong second endpoint
  term (must be (-1)^N s/(4(1-2^(1-s))(N+1)^(s+1)), Boole summation); remainder decays like N^(-1.6), not N^(-2.6), at s = 0.6+20i
  (rechecked at 40 digits; corrected formula N^(-3.6)). REFUTATIONS.
* Harvest: no 2023-2026 result conditional on Lindelof off the critical line that the half-plane would make unconditional (beyond R314/R343).
* Exploration: Diophantine approximation / Littlewood conjecture (2607.11178 shrinking targets, no zeta bridge), Khovanov homology
  (2304.02200: superpolynomials only conjecturally tied to L-functions), minimal surfaces / Willmore / harmonic maps (no zeta link): no new route.
-/

namespace Fayed.R354
end Fayed.R354
