# RVApprox: zeros of the finite Rodriguez-Villegas approximants of arXiv 2609.25564 (R338)

arXiv 2609.25564 builds U_E(q) = sum_j (-1)^j C(E,j) P(j/E) q^j from the theta-logistic profile P and finite transforms Z_E -> xi;
its concluding section leaves open whether the zeros of U_E lie on the unit circle.

- rv.py (mpmath, 50 digits; int Phi = xi(1/2) checked to 12 digits): U_E has zeros off the unit circle for every even E = 4..40;
  max | |q| - 1 | = 0.0562 (E = 4), 0.2055 (E = 10), 0.3965 (E = 26), 0.4963 (E = 40); E = 2 is on the circle (exactly: P(0) > P(1/2)).
  Independent check of E = 4: B_4 P(x) = R(p), p = x(1 - x), has roots p = 32.0103 +- 28.7832 i (not real >= 1/4), i.e. x = 2.856 +- 6.108 i.
- zEl.py <lambda> <Emax> (kernel Phi(u) e^{lambda u^2}): Z_E has all zeros on Re s = 1/2 for E <= 24 at lambda = 0, -0.3 and -1
  (H_lambda has non-real zeros for lambda < 0), so small-E real-rootedness is generic; the smallest zero (7.9-9.3 for E <= 6) is far
  from 14.1347 (Bernstein convergence O(1/E)), so 'Z_E real-rooted for all E' (sufficient for RH by Hurwitz) cannot be probed numerically.
  Above E = 24 at 50 digits the leading coefficients lose precision (spurious off-line roots); use dps >= 3E there.

- lp.py <riemann|gauss> <dps> <param> <Emax> (trapezoid quadrature, Z_E(1/2 + it) as a polynomial in u = t^2; R341): the map
  kernel -> Z_E does NOT preserve the Laguerre-Polya class: the Gaussian kernel (transform e^{-t^2/4}, no zeros) gives non-real
  zeros from E = 16 on (u = 11.68 +- 0.52i), the Gaussian-cos kernel b = 1 (transform cos(t) e^{-t^2/4}, real zeros) from E = 40;
  Riemann's kernel stays real-rooted to E = 90 for lambda = 0 and 0.5 (smallest zero 13.01 at E = 90, slowly -> 14.1347).
  So RH does not imply 'Z_E real-rooted for all E' by any LP-preservation argument, and the Hurwitz criterion is unproductive.
