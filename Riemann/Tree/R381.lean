/-!
# Round 381 — node N7 (zero-free regions and zero-density), round 9: B_2 + 2.2e-254 stands (tracker PRs #113-#119 weaker or unaccepted; Hexagon, Zenodo, HAL nothing new); Lee's k-th power rectangles (2602.14340) are superseded by the half-plane (k = 8, already 2610.20058); the hybrid twelfth moment stays the missing input; MDPs, general equilibrium and dimers give no route

* Node: N7 zero-free regions and zero-density constants (ninth round; R373 before).
* [check] tracker (github.com/GTimG/QuasiRiemannTracker, API): PR #119 merged 0.8749570194201 (weaker, same as #113), #114 open
  0.96875, #117 closed (0.9672, number fields, not accepted); no accepted value beats Bebel's 0.8733819666132075...: HOLDS.
* [check] off arXiv: hexagonmath.org newest items unrelated; Zenodo ('"zero-free region"', '"zero density"', '"quasi-Riemann"'):
  23935840 gives no global zero-free bound, 23866253 proves no new region, 23975506, 23971353, 23933662, 23932214, 23881837
  are not records; HAL: none: HOLDS.
* Adversarial (recompute B_2): the brief sent the scout to arXiv 2611.00082, an unrelated LDPC paper; Liu's B_2 is in
  Hexagon 2611.00082 (hexagonmath.org numbering): NOT TESTABLE this round (redo from the Hexagon PDF).
* Harvest: Lee 2602.14340 (primes between consecutive k-th powers for every n, k >= 86 unconditionally; zero-free rectangles
  such as 0.99168 <= sigma <= 1, H_R < |t| <= e^264.334 would give k = 70, TeX l. 182) is superseded by any half-plane below
  1 - 1/8: the 7/8 half-plane gives a prime in (n^8, n^8 + n^7] for all n, already shown in 2610.20058 (seen) and implied by
  our ShortPrimes corollary; k = 7 would need a half-plane below 6/7 < B_2: KNOWN. 2610.36069 (McShane-Rivin conjecture on
  Markov numbers, l. 607) is off-tree.
* Wild: a hybrid twelfth moment over q <= T^{1/2} feeding Heath-Brown's detector would give N(sigma, T) << T^{2(1-sigma)+0.01}
  for 0.86 <= sigma <= B_2 and could cross the contact; this is the old hybrid12s candidate (hard), no estimate in sight.
* Exploration: Markov decision processes (Bellman contraction needs a discount < 1; MDP LP duality certifies finite
  occupation-measure positivity only), general equilibrium (Hicksian/gross-substitute stability is a finite-dimensional local
  condition), dimers and matching polynomials (Heilmann-Lieb real-rootedness passes to limits only as in Polya's program; a
  matching-polynomial sequence tending to Xi would be a Jensen-type approximation, N6, with no construction): no route.
-/

namespace Fayed.R381
end Fayed.R381
