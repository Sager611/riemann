/-!
# Round 289 — node N5 (Nyman-Beurling / Baez-Duarte): Burnol's bound made precise; d_N computed to N = 200

* Node: N5 Nyman-Beurling / Baez-Duarte.
* State of the art: Burnol (math/0103058, Thm 1): liminf d_N^2 log N >= sum over zeros ON the critical line of m_rho^2/|rho|^2,
  unconditionally; under RH this equals C = 2 + gamma - log 4 pi = 0.046191417932. d_N^2 ~ C/log N under RH plus simple
  zeros (Baez-Duarte-Balazard-Landreau-Saias, 1211.5191). Newer: 2611.03924 (zero-sensitive d_n(alpha), exact zeta-free
  series to 10^6), 2606.22562 (power-law refinement only with simple separated zeros), 2405.06349 (Gram matrices,
  d_N >> (log N)^{-1/2}), 2610.03568 (invariant-subspace criterion). TREE.md N5 made precise.
* Computation: d_N^2 log N = 0.05252, 0.04818, 0.04558, 0.04620, 0.04670 for N = 10, 20, 50, 100, 200 (two independent
  scripts), near C; Vasyunin's Gram formula checked against piecewise quadrature for j, k <= 6 (error < 2e-15).
* Wild: the 7/8 half-plane gives only M(x) << x^{7/8+eps}, not enough for any power decay of d_N.
* Exploration ideas: Beurling-Malliavin density for the dilation system {1/(kx)} after a Mellin transform (N5); matrix-free
  CG/Lanczos with harmonic-sum Gram kernels to reach N = 10^4..10^5 (N5); singularity-subtracted observables for Mayer's
  transfer operator (N5/N13).
-/

namespace Fayed.R289
end Fayed.R289
