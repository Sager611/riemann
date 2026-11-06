/-!
# Round 387 — node N5 (Nyman-Beurling/Baez-Duarte), round 10: records stand (Zenodo, HAL nothing after 4 Nov; 2611.03924 still v1); Burnol's bound is a liminf only, so 2611.03924's finite value 0.0451 at N = 5120 is consistent; the B_2 zeta_K range is already recorded; our NBLp Theorem B already gives a power rate N^{-eta} for eta < 1/p - B with a non-explicit constant (explicit C needs effective 1/zeta bounds); quantum algorithms, spherical designs and curvature-dimension spaces give no route

* Node: N5 Nyman-Beurling/Baez-Duarte (tenth round; R377 before).
* [check] Zenodo ("Nyman-Beurling", "Baez-Duarte", "Nyman Beurling"): no item after 2026-11-04; HAL newest 2026-08-16
  (hal-05718412); arXiv 2611.03924 is still v1: HOLDS.
* Adversarial: Burnol's Thm 1.3 (math/0103058, TeX l. 651-660) is only liminf_{lambda -> 0} D(lambda)^2 log(1/lambda) >=
  sum m_rho^2/|rho|^2, with no finite-N error; 2611.03924 (TeX l. 192-200, 947-970) uses the same H^2(C_{1/2}) distance and
  reports d_5120^2 = 0.005283, d^2 log N = 0.0451216 < C = 0.0461914: a finite value below a liminf bound is consistent: HOLDS.
* Routine check: the zeta_K criterion for 1 < p <= 1/B_2 = 1.1449744, conditional on Liu's Hexagon 2611.00082, is already in
  COROLLARIES (NBLpK); the general Hecke L-function version is not established there: KNOWN.
* Harvest: 2611.03924 (TeX l. 923-927, 1002-1006) repeats questions already recorded (R340, R367); no newer unseen NB paper:
  NONE.
* Wild: an explicit L^p rate: our Riemann/Deep/NBLp/PROOF.md Theorem B already gives d_{N,p} <= C(p, B, eta) N^{-eta} for every
  eta < 1/p - B (p = 8/7: eta < 0.00162), stronger than any power of log N, but C is not explicit; making it explicit needs
  effective lower bounds for 1/zeta on compact strips and explicit Borel-Caratheodory constants: low.
* Exploration: quantum algorithms (phase estimation reads a spectrum only after a zeta-linked unitary is supplied: circular),
  spherical designs and Delsarte LP bounds (no map from Weil test functions to spherical harmonics), metric measure spaces
  with Ricci bounds (an RCD structure on a prime/Weil space would encode Weil positivity rather than derive it): no route.
-/

namespace Fayed.R387
end Fayed.R387
