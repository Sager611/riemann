/-!
# Round 428 — node N10 (explicit ψ, π − li and Mertens bounds), round 12: the inputs of mertensch's forming Mertens record authenticated: Chirre–Helfgott 2511.14736v1 (the only version) Corollary mertensimple, |M(x)| ≤ 3x/(π·10^10) + 11.39√x and |m(x)| ≤ 3/(π·10^10) + 11.39/√x for x ≥ 1, Theorem 1.1's identity, the 'crude estimate' remark, T = 10^10 + 1 and 11.35051; Rademacher's bound as quoted by Kadiri–Lumley–Ng; the comparison column reproduced; no pre-emption; radiocarbon dating, exterior ballistics and computational neuroscience give no route

* Node: N10 explicit ψ, π − li and Mertens bounds (twelfth round; R410 before). Round type: authenticate the quoted inputs of a
  forming record (mertensch: |M(x)| ≤ 9.77·10^4 x^0.9382 for every x ≥ 1, Riemann/Deep/MertensCH, agent still working).
* [check] Chirre–Helfgott arXiv 2511.14736 (v1 only): Corollary mertensimple (Moebius.tex l.390–395): for x ≥ 1,
  |M(x)| ≤ 3x/(π·10^10) + 11.39√x and |m(x)| ≤ 3/(π·10^10) + 11.39/√x, m(x) = Σ_{n≤x} μ(n)/n: HOLDS (the m(x) bound added
  to TREE N10). Theorem 1.1 (l.267–286): the identity with δ = π/(2T), a_∞ = sup|a_n|, |ι| ≤ δ, w = coth(δ(s − σ)) −
  tanh(δ(1 − σ)) and |ε| ≤ (π/4)(a_∞/L + a_∞/L² + I)/T² + 2a_∞/x; the remark "1/T comes from the crude estimate
  |Σ_{N<n≤N+H} μ(n)| ≤ H applied to H ≫ N/T" (l.5149–5150); T = 10^10 + 1 and δΣ|coth(δρ)/ζ′(ρ)| = 11.35051… (l.382,
  398–399): HOLDS. Its general bounded-coefficient theorem (mainthmA: T ≥ 4π, x > e²T) is a framework, not a record.
* [check] Rademacher's convexity bound as we quote it (Math. Z. 72 (1959), Thm 4): Kadiri–Lumley–Ng arXiv 2101.12263 TeX
  l.396–403, |ζ(s)| ≤ 3|1+s|/|1−s|·(|1+s|/2π)^{(1−σ+η)/2} ζ(1+η) for −1/2 ≤ −η ≤ σ ≤ 1 + η ≤ 3/2: HOLDS (the scout found no
  source in time; checked by the orchestrator). Platt–Trudgian's height 3,000,175,332,800 was authenticated in R405.
* Computation: the "published bounds combined" column of MertensCH/REPORT.md reproduced (compare.py: 10^5.96 at θ = 0.9382,
  10^5.26 at 0.9409, gains 9.35x and 11.9x); a scout using only Mertens78's all-x constants got 10^6.04 / 10^5.33, because
  the column also uses Mertens78's stated constants for x ≥ x_a: no flaw. Float J_L0 = 2∫_0^1000 |1/ζ(0.85 + it)|/|s| dt =
  17.483 (Simpson, step 0.005; mpmath agrees at three points) is below the Arb bound 19.596: HOLDS.
* State of the art: arXiv (2611.00525, 2610.11487, 2610.09567: Möbius-weighted GCD sums and short-interval cancellation),
  Zenodo and HAL (hal-05838582 is the counterexample range): no explicit power saving for M(x) besides ours: no pre-emption.
* Exploration: radiocarbon dating (Bateman chains: positive exponential sums, Prony identifiability), exterior ballistics
  (dissipative drag ODEs, spin-stability eigenvalues) and computational neuroscience (neural-field Krein–Rutman principal
  eigenvalue only): no route.
-/

namespace Fayed.R428
end Fayed.R428
