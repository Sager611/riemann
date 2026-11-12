/-!
# Round 451 — node N3 (Li/Keiper), round 12, falsification of a fresh discovery of ours (lirate's Li-coefficient remainder: mean square N^{3/2}√log N under RH, Ω((n log n)^{1/4}) unconditionally, smoothed mean 7/4, sharp weighted Arias de Reyna criterion): Lemma A's constant −1/4 and the trend's 3/4 − 1/(24n) recomputed independently; Theorem 4 holds for four more windows; the mean 7/4 is Voros's heuristic correction 2Z(0) (math/0404213, math/0506326), now proved; no mean-square law or Ω-result in seven sources; annihilating pairs, hypersingular Riesz energy and the Cesàro operator transfer no tool

* Node: N3 (12th round; R436 before). Round type: falsification and adversarial novelty of a fresh discovery of ours
  (discoveries/li-remainder.md, agent lirate, 11 Nov 2026), aimed at its constants, its smoothed-mean theorem and its prior art.
* [check] Lemma A (Riemann/Deep/LiRate/proofs.md l.13–28): λ_n^(0) = ∫_0^∞ 2(1 − cos nθ(t)) A′(t) dt, θ(t) = 2 arctan(1/(2t)),
  A′ = θ_RS′/π, recomputed by the orchestrator's own mpmath quadrature in v = θ: λ_n^(0) − t(n) + 1/4 + 1/(24n) = 8.39·10⁻⁵
  (n = 10), 3.81·10⁻⁸ (n = 50); ∫_0^∞ (A′ − log(t/2π)/(2π)) dt = −0.125: HOLDS (the scout timed out).
* [check] Trend λ̄_n (Li coefficients of π^(−s/2)Γ(1 + s/2); closed form 1 − (n/2)(γ + log 4π) + Σ_{j≥2} (−1)^j C(n,j)(1 − 2^(−j))ζ(j)
  equals the derivative formula at n = 1, 2, 3): λ̄_n − t(n) − 3/4 + 1/(24n) = 3.81·10⁻⁸, 4.20·10⁻⁹, 5.21·10⁻¹⁰, 6.51·10⁻¹¹ at
  n = 50, 100, 200, 400, i.e. ≈ 1/(240n³) as in Voros's Eq. (TS): HOLDS.
* [falsify] Theorem 4 (every φ ∈ C_c^∞(0, ∞), φ ≥ 0) on λ_n, n ≤ 10⁵ (orchestrator rerun of fft_lambda.py, r = 0.999, 2¹⁷ points,
  360 bits, peak 98 MB, λ_100000 identical): bumps on (0.1, 1), (0.9, 1), (0.5, 1) and x⁵·bump(0.3, 1): worst deviation from 7/4
  at N = 10⁵ is 0.0182 (narrow bump; allowance 30 log N/Σφ = 0.156), shrinking with N; λ_n − λ̄_n has smoothed mean 1: HOLDS.
* [known] Prior art (TeX): Voros math/0404213 l.275–282 and math/0506326 l.444–450 (from his math/0104051 Eq. (41)): the
  correction δλ_n = 2Z(0) = +7/4 improves the agreement 'in the mean', 'although this term should not be asymptotically
  meaningful': the mean 7/4 was a heuristic and Theorem 4 proves it; the trend constant 3/4 is Voros's Eq. (TS). No mean-square
  law, Ω-result or proof of the mean in Voros (also 1403.4558), Lagarias math/0404394, Maślanka math/0402168, Flajolet–Vepstas
  math/0611332 (only Coffey's inequality A_n(1,2) ≥ (n/2) log n + (γ − 1)n/2 + 1/2) or Coffey math-ph/0505052 (PDF); Arias de
  Reyna's ℓ² equivalence is his Thm 1, Cor 2, Thm 6 (attribution correct). Page, DISCOVERIES line and TREE N3 now say
  'Voros's heuristic, proved'; the SIGNIFICANT verdict rested on the Ω-result and the weighted threshold, which stand.
* [aside, N8] Riemann.fail arena: Bei Yi's beiyi-ten-d400-t024 (κ ≈ 0.673638979, 12 Nov) and Knausgård's knausgard-ten-win10-lp
  (κ = 33366251109982218737921064894121261303861213265473/49531331836748571742426383080768585205078125000000 ≈ 0.673639288,
  12 Nov: cutting-plane LP over a 16-harmonic window and 10-point weights) authenticated (metadata, Lean + nanoda, standard
  axioms); TREE N8 updated.
* Exploration: Logvinenko–Sereda / Amrein–Berthier–Benedicks (N5): Mellin dilation closure is not a Fourier spectral inequality;
  Hardin–Saff hypersingular Riesz energy (N1): normalised zero ordinates are not E_2 minimisers (N = 3: 10.93 > 9), no bound;
  Brown–Halmos–Shields Cesàro spectrum (N13): |1/ρ − 1| ≤ 1 for all ρ is RH restated: no transfer.
-/

namespace Fayed.R451
end Fayed.R451
