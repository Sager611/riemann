/-!
# Round 294 — node N4 (Robin/Lagarias), second round: the range is computational (Morrill-Platt, Briggs' algorithm), 21-free is Axler; the 7/8 half-plane extends neither

* Node: N4 Robin / Lagarias (second round; lesson: look for height-tied ranges first).
* State of the art: the range 5041 <= n <= 10^{10^{13.11485}} is Morrill-Platt (arXiv 1809.10813, Theorem "Platt"):
  Briggs' colossally-abundant algorithm rerun in 100-bit interval arithmetic; the endpoint (largest prime
  p = 29,996,208,012,611) is computational, not tied to the RH height T = 3*10^12. T only fixes Buthe's range
  B = 2.169*10^25 for |theta(x) - x| < sqrt(x) log^2 x/(8 pi) (arXiv 1410.7015: 4.92 sqrt(x/log x) <= T), used for the
  t-free families. The 21-free family is Axler (arXiv 2110.13478; Morrill-Platt had 20-free). TREE.md N4 corrected.
* Computation: at p = 29,996,208,012,611 the Fiori-Jaskari bound from the 7/8 half-plane has relative error 0.0208
  (Buthe's bound: about 7*10^-6), so it extends neither the range nor the t-free family (Q_21 = 1.00067891,
  Q_22 = 1.00067915; t = 22 not reached). Open: whether Re rho <= 7/8 for |gamma| > T enlarges Buthe's B (a scout's
  heuristic B ~ T^{8/3} is unverified; next round, node N10).
* Wild: Nicolas' primorial criterion has a direct-arithmetic verified range (ours k ~ 10^7, R119), not height-tied.
* Adversarial: 2609.26787 (Mishra-Sarkar, finite-truncation Robin) strong-tail lemma rechecked (R_21 = 6.358*10^-15 <
  6.4*10^-15); no new Oct/Nov 2026 Robin or Lagarias claim to refute.
* Exploration: homogeneous dynamics - no unconditional closed-horocycle rate beyond L^{-1/2} (2604.08753's
  T^{-1/4+eps} concerns affine SL2 x (R^2)^k orbits): none. Information theory - monotone q-Fisher information
  (2601.12944) as a heat-flow diagnostic for de Bruijn-Newman interpolations (N2/N6; weak, the flow direction is wrong).
  Circuit complexity - Liouville's Boolean degree/sparsity bounds (2411.05471) follow from lambda(2n) = -lambda(n)
  alone, so only classes containing prefix-sum tests could carry RH information (N10 barrier).
-/

namespace Fayed.R294
end Fayed.R294
