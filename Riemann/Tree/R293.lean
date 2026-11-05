/-!
# Round 293 — node N2 (de Bruijn-Newman), second round: a newer half-plane alone does not move 0.13963; barrier rows must be re-optimised

* Node: N2 de Bruijn-Newman constant (second round).
* Correction to R288: substituting q_0 = (2*0.874058504 - 1)^2 = 0.55968 (< 9/16) into 2610.12803's Table 2 leaves every
  later row unchanged: the first gate merely gains margin (+0.0244), and the algebraic bound stays Lambda <= 0.13963. The
  R288 values (0.13958 / 0.13941 / 0.13924) assumed the barrier is re-optimised for the new q_0; that needs regenerating the
  24 rows, the gates and the 4,800 Bernstein coefficients (the paper has no public code; scalar certificates take seconds to
  minutes, the Polymath boundary computation stays external). Sensitivity under that re-optimisation model:
  dLambda/dtheta = 0.405 at 7/8; Lambda(0.87) = 0.13759, Lambda(0.86) = 0.13350, Lambda(13/15) = 0.13623 (interpolated, not
  certified).
* Adversarial: a Bernstein gate coefficient (row 21, l = 4, a = 3, b = 4) recomputed exactly: 0.0013004990 > 13/10000.
* Wild: for an exact complex pair the backward-heat strip y(t) = sqrt(b^2 - 2t) does not depend on the height, so
  height-local strips give no faster shrinking.
* Exploration ideas: finite-free inverse heat root ODE r_i' = p''(r_i)/(2 p'(r_i)) as a Wasserstein/entropy flow on finite
  Xi-zero boxes (N2); parity splitting of centrosymmetric Legendre matrices to halve certified eigenproblems (N1/N2);
  Kharitonov-sector vertex tests for robustness of Jensen heat flows (N2/N6).
-/

namespace Fayed.R293
end Fayed.R293
