/-!
# Round 343 — node N9 (moments, Lindelof, subconvexity), round 5: records stand; explicit Weyl-type record 0.611 t^{1/6} log t (2602.05614) tracked; half-plane line beats every exponent-pair mu(sigma) bound; Zenodo 23178474 vacuous

* Node: N9 moments, Lindelof, subconvexity (fifth round; R327 before).
* [check] Zenodo and arXiv since R327: mu(1/2) <= 13/84 and the 2610.01035 sixth and twelfth moment exponents stand (HOLDS).
* Explicit pointwise records, now tracked in TREE N9: |zeta(1/2+it)| <= 0.611 t^{1/6} log t for t >= 3 and 0.566 t^{1/6} log t
  for t >= 8.97e17 (arXiv 2602.05614 Thm 1, Feb 2026; TeX read; before 0.618, Hiary-Patel-Yang 2207.02366; 0.63, Hiary
  1507.01261); explicit sub-Weyl 66.7 t^{27/164} for t >= 3 (Patel-Yang 2302.13444 Thm 1, sharper for t >= e^105).
* Computation (completes R327): lower convex envelope of exponent-pair bounds mu(l - k) <= k (A/B orbit of (0, 1) to depth 8,
  Huxley, Bourgain, Tao-Trudgian-Yang (89/1282, 997/1282)) with mu(1/2) <= 13/84, against the same with mu(B) = 0,
  B = 0.87338196661320752: the half-plane line wins on all of (1/2, 1); gaps 0.11379 vs 0.11331 at 0.6, 0.07281 vs 0.07186
  at 0.7, 0.04183 vs 0.03042 at 0.8, maximum 0.02255 at sigma = B (routine, COROLLARIES).
* Adversarial: Zenodo 23178474 (Fox, Lindelof for X_0(143), Lean) proves theta_Lind < 0.143 by norm_num on the hard-coded
  theta_Lind := 0.055; no L-function occurs (REFUTATIONS).
* Harvest: Gonek's conjecture I_{-k}(alpha, T) asymp alpha^{-k^2} for 1/log T <= alpha <= 1 (2610.27591, TeX l.200; proved
  there as an upper bound for 0 < k < 1/2, l.231-254). The half-plane gives it for every fixed k only for alpha >= B - 1/2 + delta,
  the trivial bounded range (Titchmarsh 14.x): nothing new.
* Exploration: spin glasses (2411.15377: ultrametric Gaussian kernels only represent zeta models), Kac-Rice (heuristic counts
  for random models, no transfer to zeta or to fixed Jensen coefficients; R318), Painleve and Riemann-Hilbert problems
  (0805.0725 is a contour-integral analogy; no tau-function realization of xi): no idea.
-/

namespace Fayed.R343
end Fayed.R343
