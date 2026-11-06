/-!
# Round 390 — node N4 (Robin / Lagarias), round 9: records stand (arXiv since 28 Oct: only 2610.39120, seen, a numerical CA2 study; Zenodo, HAL nothing new); Sole-Planat's Dedekind Psi criterion holds on our Nicolas primorials (routine); the Robin ratio is not monotone along colossally abundant numbers (first 3000); 2610.39120's questions are asymptotic; positive mass, Brascamp-Lieb and Tutte polynomials give no route

* Node: N4 Robin / Lagarias (ninth round; R380 before).
* [check] state of the art: arXiv API (Robin inequality, colossally abundant, superabundant; Nicolas + primorial; sum of divisors
  + RH) since 28 Oct 2026: only 2610.39120 (seen R380; Thm th:wu, RH iff G(n) < e^gamma for every CA2 number n > 5040; a numerical
  study quoting Morrill-Platt's 10^(10^13) range, TeX l. 56); Zenodo ('"Robin inequality"', '"superabundant"', '"colossally
  abundant"', '"Nicolas inequality"') after 2026-11-02 and HAL after 2026-10-15: nothing: HOLDS.
* Routine corollary (COROLLARIES): Sole-Planat (1011.1825) prove RH iff Psi(N_k)/(N_k log log N_k) > e^gamma/zeta(2) for all
  primorials N_k, k >= 3 (Psi the Dedekind function); as Psi(N)/N = (N/phi(N)) prod_{p <= p_k}(1 - p^-2) and the product exceeds
  1/zeta(2), the criterion holds for every primorial with p_k <= 1.78e48 by our Nicolas range (the scout timed out; derived and
  the statement checked in the abstract here).
* Computation: the first 3000 colossally abundant numbers (Alaoglu-Erdos critical values eps(p, a) = log((p^(a+1) - 1)/(p^(a+1) -
  p))/log p; the first ten match A004490): the Robin ratio R = sigma(n)/(e^gamma n log log n) falls from 1.00556 at 5040 to 0.98325
  at 55440 and 0.97305 at 720720 and keeps falling before it rises (726 decreases in 3000 steps); max R after 5040 is 1 - 6.08e-4
  at index 2347, later values are lower, so R is not monotone along CA numbers, consistent with the oscillating residual of
  2610.39120; the monotonicity premise of the task was not a published claim: HOLDS (nothing refuted).
* Harvest: 2610.39120 asks whether Q(n) <= 3/2 for all large CA numbers under RH (Question q:bound, TeX l. 278-284) and conjectures
  that Ramanujan's constants are attained along CA numbers (conj:sharpness): asymptotic, not finite computations: NONE.
* Adversarial: no new Zenodo or HAL Robin-based RH claim to test: NOT TESTABLE.
* Exploration: the positive mass theorem and the Penrose inequality (Bochner and energy integrands are nonnegative, while Weil's
  prime and archimedean terms have opposite signs; 2610.27903), Brascamp-Lieb heat-flow monotonicity (forward heat flow on
  log-concave or real-rooted data bounds Lambda from below only; R44, R71, ExploreMonotone) and Tutte polynomials (deletion-
  contraction sign control has no counterpart for Jensen coefficients; 1812.01814): no route.
-/

namespace Fayed.R390
end Fayed.R390
