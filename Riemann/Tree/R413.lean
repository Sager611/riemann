/-!
# Round 413 — node N6 (Jensen-polynomial hyperbolicity), round 10: Kim-Lee's d^{1/2+eps} large-n wedge is ineffective only through Y.-O. Kim's two 1996 theorems on nonreal zeros of derivatives (Orca candidate: an explicit version from our strip, RH height and saddle code); no effective derivative-zero rate and no fully explicit Taylor-coefficient asymptotics of xi are known; HAL and OpenAlex give nothing new; the Dirichlet-L corollary extends from q <= 125,000 to q <= 400,000 with Platt's full heights; line spectral pairs, hydrodynamic stability and catastrophe optics give no route

* Node: N6 Jensen-polynomial hyperbolicity (tenth round; R393 before). Round type: effectivity harvest, [check], input check.
* Harvest (Kim-Lee 2105.05386): Theorem 1 (TeX l.132-140) N(f; d) = O(d^{c/2}) for every c > rho (order rho < 2, zeros in the
  strip), so N(Xi; d), N(Xi_0; d) = O(d^c) for c > 1/2; Theorem 2 (l.148-165) gives n_1, n_2 only existentially, from Y.-O. Kim,
  PAMS 124 (1996) 819-830, Thm 2, and J. Korean Math. Soc. 33 (1996) 657-667, Thm 1 (k = 1/2, g(z) = f_0(z - 1/4)); Theorem 3 turns
  the zero location into hyperbolicity. The paper does not say whether Kim's n_1 is intrinsically ineffective: UNCLEAR, Orca
  candidate kimeff (explicit n_2 from the Lean-checked strip, RH to height T and our SavvaK saddle code).
* Search (arXiv API, web): no effective version of Kim's theorems and no explicit rate for nonreal zeros of derivatives of Xi:
  NONE FOUND.
* Search: no fully explicit asymptotics of the Taylor coefficients of xi valid for all n >= n0 (GORTTW 1910.01227 Thm 1.1:
  n >= c e^{d/2}, c unspecified; 1902.06330: O(log log n / log n) only): NONE FOUND.
* [check] HAL ('Jensen polynomials', 'Turan inequalities', 'Laguerre inequalities': only Planat 2026), OpenAlex (no work cites
  Kim-Lee), 2610.38870 (d - C sqrt(d log d) real zeros, not all-n hyperbolicity): nothing beats d <= 1.6141e25 or
  n^2 log(n+2) >= 24 d^3: HOLDS.
* Input check (Platt 1305.3087, TeX l.95-103 and Theorem l.867-871): GRH for all primitive characters mod q <= 400,000 to height
  max(10^8/q, 7.5e7/q + 200) for even q and max(10^8/q, 3.75e7/q + 200) for odd q; R345 used q <= 125,000 only. Orchestrator
  (exact fractions, Riemann/Deep/Jensen9/dirichlet_q400k.py): J^{d,n}(Xi_chi) hyperbolic for all n and d <= floor((1 + r^2)^2/(4r^2)),
  r = a/t0(q), for all 202,644 real primitive moduli q <= 400,000; minimum 154,735 at q = 399,999 (t0 = 293.75); values for
  q <= 125,000 unchanged (q = 3: 1,992,465,552,200,840). Routine, COROLLARIES.
* Exploration: line spectral pairs (unit-circle interlacing and Schur parameters of LPC polynomials) certify only finite
  self-reciprocal truncations, with no passage to the Laguerre-Polya class; hydrodynamic stability (Howard's semicircle theorem
  needs a sign-definite profile; the explicit-formula kernel has signed prime terms) could at best reproduce a strip; catastrophe
  optics (fold catastrophes of heat-flow zero collisions) is local and qualitative, no bound on Lambda: no route.
-/

namespace Fayed.R413
end Fayed.R413
