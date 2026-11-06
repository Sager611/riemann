/-!
# Round 376 — node N9 (moments, Lindelof, subconvexity), round 8: our 0.568 record holds and its block lemma passes a random test; Patel-Yang's sub-Weyl constant 66.7 rests on removable uniform relaxations binding at log t in [540, 565] (candidate subweyl); mapping class groups, homogenization and optimal stopping give no route

* Node: N9 moments, Lindelof, subconvexity (eighth round; R366 before; record now ours: 0.568 t^{1/6} log t, explicit611).
* [check] arXiv API and Zenodo: nothing better than 0.4821 / 0.568 t^{1/6} log t, Patel-Yang's 66.7 t^{27/164} or the
  explicit strip bounds; 2609.00537 seen before: HOLDS.
* Adversarial: the per-block bound of our certificate (Lemmas 2.5-2.9 of 2602.05614) against true block sums at 60 random
  (t, r), t in [2e6, 1e9]: max |S|/B_r = 0.1749 (t = 4,427,285.48, K = 165, r = 4), independent mpmath check agrees: HOLDS.
* Harvest: Patel-Yang 2302.13444 (sub-Weyl 66.7 t^{27/164} for t >= 3): uniform relaxations at TeX l. 557-576, 869-883,
  917-941 (floor-to-geometric block bounds, a global h_3, R(t) <= R(t_1), triangle and square-root subadditivity, a trivial
  S_1 bound); 66.70 binds at log t in [540, 545] and [560, 565] (table rows l. 1189, 1197); per-interval recomputation can
  remove them (candidate subweyl, N9).
* Computation: Phragmen-Lindelof interpolation between 0.568 t^{1/6} log t on Re s = 1/2 and 0.75 log t on Re s = 1 gives
  C(sigma) t^{(1-sigma)/3} log t with C = 0.6005, 0.6348, 0.6711, 0.7094 at sigma = 0.6, 0.7, 0.8, 0.9 (with 0.611:
  0.6366, 0.6632, 0.6910, 0.7199); a routine corollary once an explicit PL lemma with log factors is checked (not recorded).
* Wild: a fully explicit Bourgain 13/84 bound needs explicit decoupling constants for every dyadic range: not feasible soon.
* Exploration: mapping class groups (Ruelle and Selberg continuations known; dilatation growth has no arithmetic bridge),
  homogenization (Bloch approximants from truncated Euler products do not converge to the zero measure in the strip),
  optimal stopping (no canonical filtration for psi, M or Robin's sum): no route.
-/

namespace Fayed.R376
end Fayed.R376
