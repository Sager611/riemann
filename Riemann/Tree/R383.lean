/-!
# Round 383 — node N3 (Li/Keiper), round 8: records stand (Zenodo, HAL: nothing after 4 Nov); Platt 2016 is still the best rigorous GRH height for all Dirichlet characters (2610.29018 adds q <= 1e9 for real even characters, height 32 only); Brown's converse at N = 106.25 T^2 returns only the input region; Cesaro-averaged Li positivity is an RH-equivalent with a larger range but no progress; Hele-Shaw flow, rough paths and nodal domains give no route

* Node: N3 Li/Keiper (eighth round; R368 before).
* [check] Zenodo ("Li coefficients", "Li criterion", "Keiper") and HAL: no item after 2026-11-04, no range beyond our
  106.25 T^2 or Fabbian's 77.03 T^2, no new RH claim via Li's criterion (the newest Zenodo Li items are already in NOGO): HOLDS.
* Computation ([check] of the DirLi input): Platt's Math. Comp. 2016 verification (all primitive chi mod q <= 400,000,
  t0(q) = 10^8/q for q <= 125,000) is still the best rigorous GRH height for all characters; arXiv 2610.29018 (seen, R133)
  verifies real even primitive chi for q <= 10^9 only to height 32, so the DirLi ranges stand: HOLDS.
* Adversarial: Brown's converse (J. Number Theory 111 (2005), Thm 2: lambda_n >= 0 for n <= N excludes zeros from an explicit
  region) evaluated at N = 106.25 T^2 reaches only |gamma| <= T and Re rho <= 0.87338196661320752, the region we fed in; no free
  zero-free region and no inconsistency: HOLDS.
* Harvest: 2301.05779 (TeX l. 1230) asks for the relation between lambda_n and its G_n unconditionally or when RH fails; a
  computation tests only finite cases: not settleable.
* Wild: Cesaro means A_N = (1/N) sum_{n<=N} lambda_n: an off-line zero gives finitely many terms z^n with |z| = r > 1 and
  arg z not 0 (no real zeros), so A_N, like lambda_n, is negative infinitely often; RH iff A_N > 0 for all N, and the geometric
  sum damps the worst term by about gamma/N, a range gain of about 2 ln T/(2B - 1) = 76.9 in the coefficient of T^2 (scout
  estimate); not in 1304.7895 (coefficientwise criteria only); but every such range is a corollary of RH to T plus the
  half-plane, so a self-made criterion with a larger range is no progress: DEAD.
* Exploration: Hele-Shaw flow and quadrature domains (Richardson moments are domain moments with no positivity; identifying
  them with zero power sums is circular), rough paths (expected-signature kernels are positive definite, but no path law
  encodes Weil test functions), nodal domains (Courant, Pleijel and Sturm oscillation need the self-adjoint operator first):
  no route.
-/

namespace Fayed.R383
end Fayed.R383
