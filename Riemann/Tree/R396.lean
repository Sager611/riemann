/-!
# Round 396 — node N9 (moments, Lindelof, subconvexity), round 10: records stand (0.568 t^{1/6} log t, 54.6 t^{27/164}, mu(1/2) <= 13/84); the quoted inputs of both N9 discoveries authenticated in the TeX (Patel-Yang's Lemma 3.5 omission, the per-block hypotheses of 2602.05614's lemmas); scan of |zeta(1/2+it)| on [3, 600] stays below both constants; two harvest candidates (de Faveri-Pandey's savings, an exponent-database closure rerun); fair division, stochastic control and information geometry give no route

* Node: N9 moments, Lindelof, subconvexity (tenth round; R386 before), used to authenticate the inputs of our two N9 discoveries.
* [check] arXiv (axq.py: critical line, zeta, explicit / Weyl / subconvexity) and Zenodo since 25 Oct: nothing beats 0.568 t^{1/6}
  log t, 54.6 t^{27/164} or mu(1/2) <= 13/84; Zenodo 23769727 (Hardy Gram-lattice means) only quotes Bourgain's 13/84: HOLDS.
* Patel-Yang 2302.13444: Lemma 3.5's statement (TeX l. 863-892) lacks the factor (2 pi h)^{17/328} that its proof produces
  (l. 927-933); Theorem 1 states 66.7 t^{27/164} for t >= 3 (l. 60-64), as recorded in REFUTATIONS: HOLDS.
* 2602.05614: Lemmas 2.5, 2.7 and 2.9 allow per-block K, r, c and a per-interval L (TeX l. 400-417, 455-470, 570-589), as our
  0.568 uses (Explicit611 REPORT l. 22-25, 46-50); Theorem 1's 0.611 and 0.566 at l. 155-161: HOLDS.
* Computation (orchestrator; the scout timed out): mpmath scan of |zeta(1/2+it)| for 3 <= t <= 600, step 0.05: max of
  |zeta|/(t^{1/6} log t) is 0.5074 at t = 45.60 (bound 0.568), max of |zeta|/t^{27/164} is 2.985 at t = 480.40 (bound 54.6): HOLDS.
* Harvest (candidates, not started): optimise de Faveri-Pandey's savings delta_k and range nu for the sixth and twelfth moments
  (2610.01035, TeX l. 74; a finite parameter search over their proof); rerun the exponent-database polytope closure of Tao-Trudgian-
  Yang (2501.16779, l. 501-502) to test for an unreported exponent pair or beta improvement (exact arithmetic).
* Exploration: fair division (Sperner and KKM give envy-free partitions, 2412.13340, not positivity of a quadratic form for all
  test functions), stochastic control (KYP, Popov and Riccati positivity need a rational positive-real transfer function; Weil's
  kernel has digamma and infinitely many delay terms) and information geometry (Fisher convexity of log zeta holds only for real
  s > 1, already R352 and R372): no route.
-/

namespace Fayed.R396
end Fayed.R396
