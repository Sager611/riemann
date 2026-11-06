/-!
# Round 379 — node N10 (explicit psi, pi - li and Mertens bounds), round 8: records stand (the Mertens counterexample bound exp(3.485e18) is already on the node line); Fiori-Jaskari's and Chirre-Helfgott's explicit psi bounds are stated for 7/8 or with a non-optimal zero-free region, so the Lean-checked edge gives routine sharpenings (short-interval primes candidate, low); 2608.01498's Mertens-type bound holds to 1e7 and is sharp at x = 113; sign uncertainty, Kuramoto locking and dilation networks give no route

* Node: N10 explicit psi, pi - li and Mertens bounds (eighth round; R361 before).
* [check] state of the art: the Oct-Nov math.NT listings, Zenodo ('"Mertens function"', '"prime number theorem"') and HAL: the
  scout's 'record' (first Mertens counterexample below exp(3.485e18), Bi-Hou HAL hal-05838582 and Gallego Zenodo 23700212) is
  already on the N10 line (MertensX); nothing else new: HOLDS.
* Harvest: Fiori-Jaskari 2609.23222 state their quasi-theorem (TeX l. 196-223) for Re s > 7/8 with constants C(sigma, x0), and
  Chirre-Helfgott 2512.15709 (final remarks, l. 5041-5042) say their use of zero-free regions is not optimal and that T should be
  updated with higher RH verification. With the Lean-checked edge 0.87338196661320752 our ShortPrimes corollary (a prime in
  (x, x + x^{7/8}] for all x >= 1) should become (x, x + x^{theta}] for theta just above 0.8734 (routine; candidate short0873, low).
* Computation: the Mertens78 certificate hardcodes the edge B = 7/8 (cert.py l. 8, opt.py l. 4, lp.py l. 6; a batch takes about
  12 min); the gap theta - B = 0.063 comes from Rademacher's convexity bound and the Harnack/LP step, so the new edge would lower
  theta = 0.9382 by only about 0.0016 (routine, not run).
* Adversarial: 2608.01498 Thm 1.1(ii), |sum_{p <= x} 1/p - log log x - M| <= 3.69/(log x)^3 for all x >= 2, holds on [2, 1e7]
  including the left limits at primes (664,579 primes; max ratio 0.9999850 at x = 113, left-limit max 0.747, both rerun here),
  so its constant is sharp at x = 113: HOLDS.
* Wild: the explicit Mobius mean |sum_{n <= x} mu(n)/n| <= 3.17e7 x^{-0.0618} proposed by the scout is already ours (COROLLARIES):
  KNOWN.
* Exploration: sign uncertainty principles (2003.10771: eventual sign conditions on f and f-hat, not Weil's support and positivity
  constraints), Kuramoto locking (stability is positive semidefiniteness of a Laplacian with nonnegative cos weights; Weil's form
  has the prime terms with the wrong sign, which is the whole problem), dilation networks (2510.18132's multiplicative near-
  orthogonality, seen, gives widths, i.e. lower bounds, not d_N -> 0): no route.
-/

namespace Fayed.R379
end Fayed.R379
