/-!
# Round 321 — node N6 (Jensen polynomials), round: a plan for the first explicit large-n threshold; the Dirichlet Jensen corollary fails; interval hyperbolicity has no vertex test

* Node: N6 Jensen polynomials of xi.
* State of the art: no explicit large-n threshold N(d) found (two scouts timed out). Nearby 2026 results: 2608.19160 proves the
  Csordas-Dimitrov conjecture (log of the first Laguerre expression of the kernel is strictly concave, Csordas' Open Problem 4.14);
  2610.04972 proves strict infinite log-concavity of xi's Taylor coefficients (Zhu's conjecture); 2611.02944 proves PF_{k*} for
  k* = 1,810,540,931,075 from the RH height (already recorded).
* Plan for an explicit n_0(d) (Orca candidate jensen-n0): GORTTW 1910.01227 Thm 1 (hyperbolic for n >= c e^{d/2}, c
  inexplicit) plus explicit Plancherel-Rotach errors for Hermite polynomials (2109.01054 Thm 2; exact midpoint minima of H_d
  are positive for d = 10, 20, 50) or a quantitative finite free CLT (2310.15485 Thm 1.1: Levy error < C_d/sqrt n, then a
  Rouche argument once n > c_d^{2d} K(d)/eps^2).
* Wild: the Dirichlet-L analogue of our all-n range fails as stated: Chasse's sector lemma needs positive Taylor coefficients,
  and a real primitive chi gives only real (signed) ones (as R313 found).
* Exploration: hyperbolicity of interval polynomials is not vertex-determined. Cubic counterexample: all 8 vertex discriminants
  are >= 29.9 while the interior point (4.460887, 0.136860, 0.022797) has discriminant -7.4956. So there is no Kharitonov-type
  certificate (no idea). Two scouts timed out.
-/

namespace Fayed.R321
end Fayed.R321
