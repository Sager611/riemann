/-!
# Round 177 — three RH-adjacent claims triaged (two on Siegel zeros, one zero-proportion record); a 4D packing proof refuted

* 2610.08965 (Knausgård, v1 6 Oct 2026) claims, unconditionally, ≥ 83.900% of the zeta zeros distinct and ≥ 67.354%
  simple and on the critical line (TeX l. 95–96, 172–178). Its percentages follow from its formulas (zero separation
  θ = 17/20, two Fourier windows; our mpmath: 0.8390046 and 0.6735383), but the method rests on an "unconditional"
  pair-correlation or energy input (l. 780–808) and on a computer search that its Lean file takes as axioms (l. 69–74).
  Deep audit running.
* 2610.30476 (Guo, v1 21 Oct 2026) claims an elementary proof that L(s, χ) has no real zero with (1 − β) log q ≤ 1/100
  for every quadratic χ: one wrong derivative (l. 1164–1168, harmless); every displayed numerical certificate
  recomputes. Deep audit running.
* 2610.31402 (Banks, v1 22 Oct 2026) claims an effective exclusion of Landau–Siegel zeros over every fixed base field,
  resting on an "OpenAI argument" and a Masser-type zero estimate (l. 183–214, 271–280). Deep audit running.
* Correct as tested: 2610.12816 (Robin's criterion via reciprocal LCM powers; κ = 10⁵, 5041 ≤ n ≤ 10⁶, max ratio
  0.98547 at n = 10080); Zenodo 23583993 (the N = 64 rule is primitive, period 2⁶⁴ − 1). 2609.18414 is a valid but
  weaker bounded-intervals result, not a record.
* Refuted (Orca audit): the claimed proof in 2609.25120 v2 that D₄ is an optimal packing in R⁴ (REFUTATIONS.md).
-/

namespace Fayed.R177
end Fayed.R177
