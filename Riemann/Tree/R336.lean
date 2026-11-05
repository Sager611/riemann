/-!
# Round 336 — node N10 (explicit psi, pi - li, Mertens), round 5: no explicit power-saving Mertens bound exists yet (candidate: mertens78); explicit 1/zeta bounds only in classical regions; theta(x) < x to 1.39e17 stands

* Node: N10 explicit prime and Mertens bounds (fifth round; R323 before).
* State of the art (Mertens): no explicit bound |M(x)| <= C x^theta with theta < 1, C and x_0 explicit, is published, even
  after the 2026 half-plane (2604.24505 says 'explicit error term' only; 2608.01498 weights prime sums; 2208.06141 gives
  sub-power savings). The explicit records stay Ramare's 0.013 x/log x (x >= 1,078,853) and Cohen-Dress-El Marraki's x/4345.
  Candidate open: an explicit power saving from the Lean-checked half-plane (Orca brief mertens78).
* Inputs for it: explicit 1/zeta bounds are published only in the classical, Littlewood and Korobov-Vinogradov regions
  (2208.06141 Thms 1-3: orders (log t)^{11/12}, (log t)^{11/12}(log log t)^{-3/4}, (log t)^{2/3}(log log t)^{1/4});
  2602.06199 is conditional on Re s = 1; no half-plane version: KNOWN.
* [check] theta(x) < x for 0 < x <= 1.39e17 (Platt-Trudgian 2016) not extended: HOLDS.
* Harvest: 2609.23222 again states no open computational question.
* Computation scout (Mertens numerics to 2e8) timed out.
* Exploration: Bayesian / statistical tests of zeros (2101.10336, 1509.08703) give heuristics, no falsifiable node constant
  (no idea); Landau sampling density and Gabor frames (1901.00356, 2609.03418) control average density only (no idea);
  arithmetic topology / Alexander polynomials (2511.03446: unit-circle roots from cyclotomic structure) fail the DH test.
-/

namespace Fayed.R336
end Fayed.R336
