/-!
# Round 311 — node N5 (Nyman-Beurling), third round: Baez-Duarte's natural-number criterion extended to L^p, 1 < p < 2, with the endpoint (agent nblp); reviewed CORRECT, novelty not found; unconditional for p <= 1.14497

* Node: N5 Nyman-Beurling / Baez-Duarte (third round; verification of the Orca agent nblp).
* Theorem A (nblp, Riemann/Deep/NBLp/PROOF.md): for 1 < p < 2, zeta(s) != 0 on Re s > 1/p iff 1_(0,1) is in the
  L^p(0,1)-closure of span{rho_{1/k}: k >= 2} (equivalently the span is dense in the step-function subspace S_p; literal
  density in L^p(0,1) is false). Theorem B: if zeta != 0 on Re s > theta and p < 1/theta, the Riesz-Mobius approximants
  1 + sum_{k<=N} mu(k)(1 - k/N) rho_{1/k} have L^p norm << N^{-eta}, eta < 1/p - theta. Unconditional for
  1 < p <= 8/7 (7/8 half-plane) and, with the Lean-checked Re s > 0.87338196661320752, for 1 < p <= 1.14497.
* Review: independent reviewer CORRECT (endpoint sufficiency: Lemma A4 residue signs give B(eps) = O_p(eps), A5 gives
  O_p(eps^{1-1/p}), A6 finite-span approximation); SIGNIFICANT if new. Adversarial scout: the endpoint proof uses fixed
  eps-strips (sigma' = 1/p - eps/2, sigma' + eps > 1/p), N -> infinity before eps -> 0, nothing on the line Re s = 1/p:
  HOLDS.
* Novelty: no sufficiency result for p != 2 found (arXiv phrasings; 2606.16097 Prop. 3 is necessity only);
  Dimitrov-Oliveira 1608.07887 Thm 2 gives L^p only for theta < 1/p strictly (open range, routine); Gallardo-Gutierrez-
  Seco 2206.02654 Thm A necessity only; Zenodo nothing. Logged in DISCOVERIES (N5).
* Own check: checks.py rerun (14 s): ||chi + H_eps||_{1.1} = 0.360, 0.205, 0.110, 0.057 for eps = 0.4, 0.2, 0.1, 0.05;
  agent's run: approximant norms 0.34595 (N = 10) and 0.05923 (N = 100) at p = 1.1. The independent norm scout timed out.
* Wild: a Dirichlet-L extension needs character-dependent gamma factors and twisted Mobius identities at the endpoint
  (idea). Exploration: type/cotype multipliers (2601.04803), H^2 kernels for Nyman-Beurling (2610.03568 Thm 6.1) and
  L^p dilation completeness (1612.06791, 2511.18188) suggest L^p tools but no criterion for p > 2 (N5, ideas).
-/

namespace Fayed.R311
end Fayed.R311
