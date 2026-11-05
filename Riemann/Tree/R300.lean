/-!
# Round 300 — node N8 (critical-line proportion), second round: 67.354% re-verified from the TeX constants; Levinson's 0.3428 reproduced; nothing newer

* Node: N8 critical-line proportion (second round).
* State of the art: 2610.08965 (>= 67.354% of zeros simple and on the line) is still v1 and still the record (Oct/Nov
  2026 listings; 2610.23467 counts zeros simple OR on the line, not both).
* Adversarial recheck of our record: the TeX coefficients of Theorem thm:line match the ancillary JSON; theta = 17/20,
  1 - sum |a'_j| = 49333/62500 > 0.789, ratio 1669159/2478195 = 0.67354 (exact arithmetic): HOLDS.
* Computation: Levinson's method with theta = 1/2, P(x) = x, Q(x) = 1 - x (Young 1002.4403, Thm 1) gives
  R* = 1.28691, c = 2.32986, kappa = 0.342762: reproduced.
* Wild: explicit Selberg-type density near the line (1910.08274 Thm 1, sigma <= 0.831) allows c <= (0.831 - 1/2) log T
  ~ 9.50 at T = 3*10^12, while the 7/8 half-plane matters only for c >= (3/8) log T ~ 10.76: no gain. No new
  100%-on-the-line claim to refute.
* Exploration: Elias-Williamson-type induction (Hodge-Riemann at one step gives hard Lefschetz at the next, by
  deformation; 2408.10155) would need a graded decomposition of finite Weil windows (N1, idea, no candidate);
  reflection positivity on R_+ forces Laplace mixtures of positive measures (1312.6161 Thm 2.2), a test class for
  Weil's distribution (N1, idea); Duursma code zeta functions: monotone coefficients plus Enestrom-Kakeya and Rouche
  put roots on the circle (0704.3903), and code-RH forces formal self-duality (1505.01993) (N12/N6, weak).
-/

namespace Fayed.R300
end Fayed.R300
