/-!
# Round 361 — node N10 (explicit psi and pi - li bounds), round 7: records stand; our short-interval prime hand-overs and Hurst's Mertens oscillations rechecked; pi(x) < li(x) is proved only to 1e19 and RH to 3e12 could reach about 3.6e21 in principle

* Node: N10 explicit psi and pi - li bounds (seventh round; R352 before).
* [check] state of the art: the Oct-Nov listings and Zenodo give no new record; 2610.27830 Thm 2 (pi(n^2 + n) > pi(n^2) claimed
  equivalent to pi(x + sqrt x) > pi(x)) repeats the gap recorded in R121: HOLDS.
* Computation (sign ranges): theta(x) < x and pi(x) < li(x) are proved for 1 <= x <= 10^19 (Buthe 1511.02032 Thm 1, TeX l. 39, 42);
  no newer range found; the scaling sqrt(X) log X ~ T with RH to T = 3e12 suggests X ~ 3.6e21 in principle, if the explicit zero
  sums can be evaluated (candidate).
* Adversarial (our prime in (x, x + x^{7/8}] for every x >= 1, R323): the hand-overs overlap (1.12016e18 < 4e18 for the
  maximal-gap tables, 1476^{8/7} = 4186.1326 < 4186.14, r(e^{41.56}) = 0.99992 < 1) and 2609.23222 Thm 2 holds for x > x_0 >= e^3
  as used: HOLDS.
* [check] (Mertens oscillations): limsup M(x)/sqrt(x) > 1.826054 and liminf < -1.837625 (Hurst 1610.08551: lattice dimension 800,
  nu = 17,000, 14,400 zeros) stand; beating 2 would need weeks of reduction (each of Hurst's runs took about 35 days): HOLDS.
* Harvest: Kim-Nguyen 2502.21021 (l. 847) expect parameter optimization to lower their first-counterexample bound exp(1.96e19) for
  |M(x)| < sqrt(x) (candidate) and (l. 855-862) ask to enumerate q_N(x) candidates; 2607.07566 (l. 1087) asks for a practical
  benchmark of Hirsch-Kessler-Mendlovic against Helfgott-Thompson for computing M(x).
* Exploration: ML-guided program search (AlphaEvolve, FunSearch: only the sign-uncertainty constant 0.3279, 2511.02864 Problem 11,
  R323; no RH-node constant), cryptographic lattice reduction (BKZ 2.0, G6K: no application to Mertens or Skewes lattices found;
  fpylll does not build here), tensor networks and structured matrices (prolate/DPSS Toeplitz matrices have commuting tridiagonals,
  2006.00427, but Weil Galerkin matrices have no known Toeplitz or H-matrix structure): no new route.
-/

namespace Fayed.R361
end Fayed.R361
