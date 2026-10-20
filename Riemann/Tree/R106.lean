/-!
# Round 106 — what a uniform 7/8 zero-free half-plane settles (scouts, corrected by the orchestrator)

| item | before → after | verdict |
|------|----------------|---------|
| ψ(x) − x | Vinogradov–Korobov x·exp(−c(log x)^{3/5}(loglog x)^{−1/5}) → O(x^{7/8+ε}) | NEW: first power saving (a scout said "weaker": wrong, x^{7/8} = x·e^{−(log x)/8}) |
| M(x) | same VK shape → O(x^{7/8+ε}) | NEW |
| primes in short intervals | BHP x^{0.525} already better than x^{7/8+ε} | weaker than known |
| zero density | DH known for σ ≥ 25/32 (Bourgain); now N(σ,T) = 0 for σ > 7/8 | partial |
| Lindelöf μ(σ) | Littlewood's three-circles argument should give μ(σ) = 0 for σ > θ (scouts disagree: σ > 15/16 or no control) | check Titchmarsh §14.2 |
| Landau–Siegel zeros | if the half-plane is uniform for every Dirichlet L-function (abstracts say so), no Siegel zeros: effective Siegel bound, effective class numbers | SETTLES OPEN PROBLEM (if uniform) |
| Linnik's constant | a half-plane σ > θ gives about p ≪ q^{1/(1−θ)+ε} = q^{8+ε}; Xylouris L = 5 is better | weaker than known |

Lean targets once the port lands: no_siegel_zero (direct corollary of the Dirichlet version),
ψ(x) = x + O(x^{7/8+ε}), M(x) = O(x^{7/8+ε}).
-/
