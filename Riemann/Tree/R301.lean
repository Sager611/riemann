/-!
# Round 301 — node N5 (Nyman-Beurling), second round: two routine corollaries - explicit unconditional liminf d_N^2 log N >= 0.046191417926, and density in L^{8/7}

* Node: N5 Nyman-Beurling / Baez-Duarte (second round).
* State of the art: Burnol (math/0103058 Thm 1.3): liminf D(lambda)^2 log(1/lambda) >= sum over zeros on the line of
  m_rho^2/|rho|^2 (multiplicities squared), hence the same for d_N; TREE.md N5 correct. No 2025-2026 d_N record
  (2607.12084: damped Baez-Duarte approximants, global boundedness open).
* Height-tied corollary: with RH to T = 3,000,175,332,800 (Platt-Trudgian 2004.09765 Thm 1; multiplicities only help
  since m^2 >= m) and C = sum_rho 1/(rho(1 - rho)) = 2 + gamma - log 4 pi unconditionally, the tail
  sum_{|gamma| > T} 1/gamma^2 <= 5.7*10^-12 (Trudgian's N(T) bound) gives liminf d_N^2 log N >= 0.046191417926
  unconditionally (COROLLARIES; routine).
* Wild corollary: Beurling (1955): the span of rho_a, 0 < a <= 1, is dense in L^p(0,1) iff zeta != 0 for Re s > 1/p;
  with the 2026 half-plane it is dense in L^{8/7}(0,1) (COROLLARIES; routine). The natural-number version is open for
  p != 2 (Baez-Duarte math/0505453 remark).
* Computation: independent Gram-matrix check gives d_N^2 log N = 0.05252, 0.04818, 0.04558 at N = 10, 20, 50
  (TREE.md values confirmed).
* Exploration: Ehrhart polynomials - reflexive reciprocity L(-z-1) = L(z) plus sign patterns of h* force roots on
  Re z = -1/2 or real (1503.05739, 2610.40496) (N6, weak); Muntz theory - near-best Muntz approximants need
  exponentially large coefficients (2208.04176 Thm 2), suggesting Gram ill-conditioning as the mechanism behind the
  logarithmic decay of d_N (N5, idea); Biane-Pitman-Yor-type laws with E[X^s] = |D|^{-s/2} xi_K(s) (2411.08863)
  turn the de Bruijn-Newman flow into adding Gaussian noise to log X (N2, idea; only lambda_1, lambda_2 > 0 proved).
-/

namespace Fayed.R301
end Fayed.R301
