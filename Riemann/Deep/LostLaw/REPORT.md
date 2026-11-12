PARTIAL: the asymptotic lost-sample law is not refuted; its exact finite-μ equality is rigorously falsified at μ = 10 — certified (Arb) two-sided brackets of log(λ/λ′) contain 4g₂(y) for all single lost samples tested at μ = 12 (7) and 16 (8) (y = 0.14–1.94) and for 4 of 8 at μ = 10, the other four (y = 0.48–0.97) lying 0.003–0.028 below; five lost samples cost Σ4g₂ + [1.02, 1.44] (μ = 10), + [0.62, 1.48] (μ = 16), and {γ₁, γ₁₄, γ₃₅} costs certifiably more than its three single costs (μ = 10); the merged-pair cost and the crossover across δ* are certified; converged Galerkin shows the multi-sample excess decaying (fixed five-point set +1.22 → +0.29, μ = 10 → 30) and the single error κ(μ) = −0.12 ± 0.08 on 8 ≤ μ ≤ 50 (band means −0.18, −0.11, −0.07 ≈ −0.5/√μ, equal for all y to ±0.02, oscillating as zeros cross t = 4πμ). Proved: an exact reduction to full-form data, log(w₀/(1 − θ)) ≤ log(λ/λ′) ≤ log(w₀/(1 − θ)) + log(1 + ρ/w₀) + log(1/(1 − r)); certified: MoriRate's sufficient target 1 − θ ≥ e^{−4g₂} fails at μ = 16, 30 (1 − θ ≤ e^{−178}, e^{−353}). Open: the e^{o(1)} bridge, μ = 20 certificates (order 1100 > memory cap), convergence of C(μ) = −log λ − 4πμ + 5 log μ ∈ [−15.56, −14.70] (5 ≤ μ ≤ 45).

# LostLaw — the lost-sample law for λ_min: falsification at larger μ, what is provable, the O(1) term (lostlaw, Orca, 12 Nov 2026)

Law (MoriRate): deleting one zero at height γ₀ = 2πμy multiplies λ_min of the windowed Weil form (even sector, support [−a, a],
a = ½ log μ, Q = Σ_ρ|û(γ)|² under RH) by e^{−4g₂(y)}, g₂(y) = arccosh(2/|y|) (|y| < 2), 0 (|y| ≥ 2); a pair closer than δ* ≈ e^{−3.2μ} is
one lost sample; costs add. λ = λ_min(Q), λ′ = λ_min of the model form (zeros deleted or moved), Σ4g₂ = the law's prediction.
Code (this folder, copied from .scratch/lostlaw; MoriRate untouched): gal.py (Galerkin: ChekAudit's closed-form matrix in Arb, Ritz on span{x, R⁻¹b, R⁻²b}, Arb
Rayleigh quotients), certmulti.py (Chekmachev's certificate for model forms), tailritz.py (upper bound for λ beyond the solvable N),
theta_rig.py; derivations DERIVATION.md. All runs single-threaded, sequential, under an RSS watchdog (kill at 290 MB); kept runs ≤ 280 MB
(/usr/bin/time -l). Over the cap: gal.py N = 800 at 416 bits (369 MB, self-abort, replaced by N ≤ 520) and certmulti.py at μ = 20 (killed at
291 MB in assembly: certificate order ≈ 1100).

## Short answer
1. (i) The asymptotic law is not falsified. Certified brackets (§1) contain 4g₂ for every single deletion and merged pair at μ = 12 and 16
   (width 0.51, 0.79–0.86) and for 4 of 8 singles at μ = 10 (width 0.35 with 3× Chekmachev's nodes); at μ = 10 the law is certifiably inexact (single costs 0.003–0.028 below
   4g₂ at y = 0.48, 0.76, 0.79, 0.97; three to five lost samples cost Σ4g₂ + 0.17 … 1.44; F({γ₁, γ₁₄, γ₃₅}) ≥ 19.078 > F(γ₁) + F(γ₁₄) + F(γ₃₅)
   ≤ 19.020, F = log(λ/λ′), so the cost is certifiably super-additive; a resolved pair raises λ), as an o(1) law may be. Converged Galerkin: the excess of
   fixed multi-point sets decays with μ; the single error κ(μ) shows no failure of o(1) and no clean decay up to μ = 50 (band means
   ≈ −0.5/√μ, the size of MoriRate Lemma B's tip term O(μ^{−1/2} log μ)). μ = 20 was not certifiable under the cap.
2. (ii) Proved (DERIVATION.md §1–4): Proposition 1 (exact identity) reduces the cost to the ground state's energy share w₀ at ±γ₀ and the
   extremal share 1 − θ = inf Q′/Q, both of the undeleted form, up to log(1 + ρ/w₀) + log(1/(1 − r)) ≈ 0.01 + r. MoriRate's scale-free
   sufficient target for the upper half fails at every computed window by 9–337 nats (certified), the margin ≈ 2μP(y₀) − 4g₂(y₀), P the
   Green potential of the deficit, whose rate −log(1 − θ) ≈ 2μP(y₀) is matched to 5 nats. Each half of the law needs a sampling-side and a
   ground-state-side estimate sharp to e^{o(1)}; neither is proved.
3. (iii) C(μ) = −log λ − 4πμ + 5 log μ: certified [−15.555, −15.436] (μ = 5), [−15.518, −15.367] (10), [−15.408, −15.165] (12),
   [−15.401, −15.009] (16), [−15.138, −14.697] (45); uncertified Galerkin −15.35 … −15.62 on 5 ≤ μ ≤ 40 with ±0.1 fluctuations (fit with 4π
   fixed: exponent 4.88 instead of 5). A limit near −15.3 ± 0.2 fits; a drift ≤ 0.5 is not excluded. The L-function comparison was not possible (§4).

## 1. Certified brackets
Lower end log(ε_full/X′), upper end log(X_full/ε′): ε = certmulti.py lower bound (rigorous Lemma 5.5 check at the reported ε), X = Arb
Rayleigh quotient of an explicit vector (gal.py N = 520; for X_full tailritz.py). Zeros: 179-digit values in certified balls (sign change
of Hardy's Z; the deleted columns use the balls). excess = bracket − Σ4g₂; logs of the certified endpoints in double precision (margins
≥ 0.003 ≫ rounding, zero-ball radii ≤ 1e−140).

1a. μ = 10, K2 = 603 (3 × Chekmachev's): λ ∈ [1.2546e−43, 1.4593e−43] (width 0.151). 20 sets 196 s, 228 MB; A-sets 54 s, 202 MB.
| deletion (y) | certified log(λ/λ′) | Σ4g₂ | Galerkin | excess |
|---|---|---|---|---|
| γ₁ (0.225) | [11.272, 11.617] | 11.500 | 11.453 | [−0.228, +0.117] |
| γ₂ (0.335) | [9.579, 9.924] | 9.896 | 9.761 | [−0.317, +0.028] |
| γ₄ (0.484) | [8.038, 8.383] | 8.386 | 8.219 | [−0.348, −0.003] |
| γ₉ (0.764) | [6.100, 6.445] | 6.467 | 6.282 | [−0.367, −0.022] |
| γ₁₄ (0.968) | [5.043, 5.388] | 5.416 | 5.224 | [−0.374, −0.028] |
| γ₃₅ (1.781) | [1.661, 2.015] | 1.966 | 1.843 | [−0.305, +0.049] |
| γ₄₀ (1.957) | [0.826, 1.184] | 0.839 | 1.008 | [−0.013, +0.345] |
| γ₁, γ₂ | [21.366, 21.753] | 21.396 | 21.548 | [−0.030, +0.357] |
| γ₁₄, γ₁₅ | [10.359, 10.741] | 10.519 | 10.540 | [−0.160, +0.223] |
| γ₁, γ₂, γ₃ | [30.754, 31.109] | 30.586 | 30.935 | [+0.168, +0.523] |
| γ₁, γ₁₄, γ₃₅ | [19.078, 19.455] | 18.882 | 19.259 | [+0.196, +0.573] |
| A4 = γ₄, γ₉, γ₁₄, γ₂₁ (y 0.48–1.26) | [24.827, 25.243] | 24.406 | 25.030 | [+0.421, +0.838] |
| A5 = A4 + γ₂₇ (1.51) | [28.580, 29.005] | 27.561 | 28.783 | [+1.019, +1.444] |
| γ₁₀ → γ₉ + δ, δ = 1e−40 … 1e−15 | [5.942, 6.287] | 6.310 (γ₁₀ lost) | 6.123 | [−0.368, −0.023] |
| γ₁₀ → γ₉ + δ, δ = 1e−10, 1e−5 | [−0.341, −0.008] | 0 (resolved) | −0.159 | the resolved pair raises λ |

1a′. μ = 12, K2 = 529 (2.2 ×): λ ∈ [3.1064e−54, 3.9631e−54] (width 0.244; tailritz 520 → 1000). 10 sets 139 s, 222 MB.
| deletion (y) | γ₁ (0.187) | γ₂ (0.279) | γ₆ (0.498) | γ₁₂ (0.749) | γ₁₉ (1.004) | γ₂₇ (1.255) | γ₄₆ (1.787) | γ₁, γ₁₉, γ₄₆ | γ₆, γ₁₂, γ₁₉ |
|---|---|---|---|---|---|---|---|---|---|
| certified | [11.917, 12.430] | [10.231, 10.744] | [7.822, 8.336] | [6.097, 6.613] | [4.781, 5.300] | [3.689, 4.212] | [1.475, 2.015] | [19.516, 20.078] | [19.774, 20.321] |
| Σ4g₂ | 12.233 | 10.634 | 8.266 | 6.555 | 5.249 | 4.166 | 1.933 | 19.415 | 20.070 |
All singles contain 4g₂ (ε/Galerkin 0.74–0.76 is too loose to see κ(12) = −0.17); the spread triple costs Σ4g₂ + [0.10, 0.66].

1b. μ = 16, K2 = 427 (Chekmachev's): λ ∈ [1.6557e−75, 2.4495e−75] (width 0.392). 20 sets 800 s, 280 MB; three sets 180 s, 220 MB.
| deletion (y) | certified | Σ4g₂ | Galerkin |
|---|---|---|---|
| γ₁ (0.141) | [12.896, 13.683] | 13.388 | 13.313 |
| γ₂ (0.209) | [11.217, 12.005] | 11.794 | 11.635 |
| γ₁₀ (0.495) | [7.669, 8.459] | 8.294 | 8.087 |
| γ₃₀ (1.008) | [4.586, 5.386] | 5.232 | 5.004 |
| γ₆₉ (1.790) | [1.309, 2.157] | 1.921 | 1.727 |
| γ₇₇ (1.942) | [0.563, 1.426] | 0.972 | 0.980 |
| γ₁, γ₂ / γ₃₀, γ₃₁ / γ₁, γ₃₀ | [24.896, 25.711] / [9.908, 10.712] / [18.286, 19.095] | 25.181 / 10.354 / 18.619 | 25.314 / 10.325 / 18.703 |
| γ₁, γ₂, γ₃ / γ₁, γ₃₀, γ₆₉ / γ₁₀, γ₃₀, γ₅₀ | [36.106, 36.936] / [20.509, 21.357] / [16.697, 17.540] | 36.276 / 20.540 / 17.014 | 36.524 / 20.927 / 17.115 |
| γ₁, γ₁₀, γ₃₀, γ₅₀, γ₆₉ | [32.943, 33.804] | 32.323 | 33.361 |
| γ₂₀ → γ₁₉ + δ: 1e−40, 1e−30 / 1e−20 / ≥ 1e−15 | [5.812, 6.606] / [0.505, 1.096] / [−0.500, +0.310] | 6.448 / – / 0 | 6.229 / 0.923 / −0.082 |
The crossover δ* lies in (1e−30, 1e−15), certified, with a certified intermediate cost at δ = 1e−20 (MoriRate's fit: δ* ≈ 1.2e−21).

1c. μ = 20: not certified. The certificate has order ≈ 1100 (|K₊| = 788, |G| = 8032; CCMFMM's particle count 1101), a dense Arb Gram
matrix of 1.25·10⁶ balls; killed at 291 MB in assembly. ExploreFMM's O(n)-memory structured LDL (64 MB at μ = 20) would need deletion
columns as extra borders (CCMFMM's border mechanism); not done.

Tightness: ε/Galerkin = 0.61–0.67 with Chekmachev's nodes, 0.79–0.83 with 3×, nearly the same for the full and every model form, so
brackets shrink only with more nodes; tailritz.py lowers the upper bound for λ by 0.018 (μ = 5), 0.030 (10), 0.026 (16) nats.

## 2. Converged Galerkin: the error term of the law
log-ratios converge fast in N (μ = 10: N = 200 vs 400 ≤ 0.004; μ = 16: 400 vs 520 ≤ 0.015; μ = 20: 300 vs 480 ≤ 0.006), λ itself slowly
(2–3 % per doubling). κ(μ) = mean of log(λ/λ′) − 4g₂(y) over single deletions with 0.3 ≤ y ≤ 1.5 (spread = max − min):
| μ | 8 | 10 | 12 | 14 | 16 | 18 | 20 | 22 | 24 | 25 | 26 | 28 | 30 | 35 | 40 | 45 | 50 |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| κ | −0.232 | −0.184 | −0.173 | −0.099 | −0.225 | −0.145 | −0.082 | −0.179 | −0.178 | +0.078 | −0.143 | −0.091 | −0.129 | −0.042 | −0.129 | −0.040 | −0.014 |
| spread | 0.04 | 0.06 | 0.07 | 0.05 | 0.05 | 0.03 | 0.02 | 0.03 | 0.02 | 0.01 | 0.02 | 0.04 | 0.03 | 0.03 | 0.01 | 0.01 | 0.00 |
(N = 200–520; μ = 35–50 with 250/310-digit zeros, 9–12 deletions each; μ = 40 at N = 400 vs 460 moves κ by 0.002.)
- y-independence: within one μ ≥ 18 all single offsets agree to ±0.02 (relative law λ(Z∖γ)/λ(Z∖γ′) = e^{4g₂(y′) − 4g₂(y)} to 2–4 %);
  the first zero is +0.12…+0.15 above κ.
- Mean κ over μ ∈ [8, 16], [18, 26], [28, 50]: −0.183, −0.108, −0.074 (≈ −0.63, −0.51, −0.46 times 1/√μ at the band centres): slow decay,
  roughly −0.5/√μ; one positive value (μ = 25).
- Edge oscillation: κ is smooth in μ and moves as zeros cross t = 4πμ: μ = 10.00 … 10.24: −0.187 (flat), then −0.146, −0.074, +0.026 after
  the edge passes γ₄₂ = 127.52; μ = 20.00 … 20.30: −0.078 → +0.026 → +0.002 around γ₁₁₀ = 253.07; μ = 30.00 … 30.10: −0.121 → −0.145 → −0.124
  (peak-to-peak 0.21, 0.10, 0.025 over 1.5, 1.7, 0.8 zero spacings). Pairs oscillate with peak-to-peak 0.08–0.15 (μ ≈ 10, 20), 0.02 (μ ≈ 30).
- Several deletions: excess of fixed configurations over Σ4g₂: A5 (y = 0.5, 0.75, 1, 1.25, 1.5) +1.22, +0.71, +0.67, +0.62, +0.29 at
  μ = 10, 16, 20, 25, 30; A4 (first four) +0.62, +0.41, +0.24, +0.36, +0.13, +0.14, +0.14 (μ = 10 … 40); A3 (0.5, 1, 1.5) +0.15, +0.12, +0.13,
  +0.19, −0.03; A2 (0.5, 1.5) +0.09, 0.00, +0.05, +0.08, −0.05. Relative to the measured single costs the A2 interaction is +0.44, +0.45, +0.20,
  −0.07, +0.21: κ is paid about once per configuration, not per deletion; what remains decays (consistent with MoriRate Lemma B's additivity).
- Pairs: δ ≪ δ* costs the deletion of the moved zero (6 digits), δ ≫ δ* costs −0.06 … −0.16 (helps) at every μ; at μ = 30, δ = 1e−40 is
  already partly resolved (5.98 against 6.37): δ* ≈ 2e−41, as MoriRate's fit predicts.

## 3. Proved statements (DERIVATION.md)
- Proposition 1 (exact; Z ∖ {±γ₀} a uniqueness set for the form domain, true for ζ and every model here): with θ = 2⟨e₀, A⁻¹e₀⟩,
  w₀ = 2|ψ̂₁(γ₀)|²/λ, r = λ′/λ, ρ = λΣ_{i≥2}2c_i²/(E_i(E_i − λ′)) ≤ λ/(E₂ − λ): 1 − θ = r[w₀/(1 − r) + ρ], hence
  log(w₀/(1 − θ)) ≤ log(λ/λ′) ≤ log(w₀/(1 − θ)) + log(1 + ρ/w₀) + log(1/(1 − r)). Numbers: ρ/w₀ ≈ 0.011; ln w₀ = −164.7,
  ln(1 − θ) = −178.0 (μ = 16, γ₁); −337.3, −353.2 (μ = 30).
- MoriRate's sufficient target fails where tested (theta_rig.py: θ ≥ 2(b·u)²/(uᵀMu), Arb, certified zero balls):
| μ | y | 1 − θ ≤ | e^{−4g₂(y)} | margin (nats) |
|---|---|---|---|---|
| 16 | 0.141 / 0.495 / 1.008 / 1.790 | e^{−177.9} / e^{−130.7} / e^{−72.2} / e^{−11.3} | e^{−13.4} / e^{−8.3} / e^{−5.2} / e^{−1.9} | 164 / 122 / 67 / 9 |
| 30 | 0.075 / 0.502 / 1.005 / 1.511 | e^{−353.1} / e^{−242.4} / e^{−133.7} / e^{−48.7} | e^{−15.9} / e^{−8.2} / e^{−5.2} / e^{−3.1} | 337 / 234 / 128 / 46 |
  The margin ≈ 2μP(y) − 4g₂(y) grows linearly in μ; an asymptotic refutation needs a proof of 1 − θ ≤ e^{−cμ} (reviewer).
- Rate of 1 − θ: −log(1 − θ(γ₀)) = 2μP(y₀) + d, P(y₀) = ∫_{−2}^{2} g_Ω(t, y₀)(−log|t|)dt (P(0) = 2π, P(0.5) = 4.013, P(1) = 2.152,
  P(1.5) = 0.750), d ∈ [−5.3, +1.8] at μ = 16, 30; for Mori's proxy P(y₀) is a rigorous lower bound (his certificate at base point y₀).
- What each half needs: lower half = a test vector with Q′(u)/Q(u) ≤ w₀e^{−4g₂+o(1)} plus a lower bound for w₀; upper half = a lower bound
  for 1 − θ (sampling/multiplier side) plus an upper bound for w₀ (ground-state side), each sharp to e^{o(1)}.

## 4. The O(1) term C(μ) = −log λ − 4πμ + 5 log μ
| μ | certified λ ≥ (C ≤) | Arb λ ≤ (C ≥) | C |
|---|---|---|---|
| 5 | 8.1509e−18 (Chekmachev's nodes) | 9.1780e−18 (tailritz 200 → 700) | [−15.555, −15.436] |
| 10 | 1.2546e−43 (3× nodes) | 1.4593e−43 (tailritz 520 → 1000) | [−15.518, −15.367] |
| 12 | 3.1064e−54 (2.2× nodes) | 3.9631e−54 (tailritz 520 → 1000) | [−15.408, −15.165] |
| 16 | 1.6557e−75 | 2.4495e−75 (tailritz 520 → 900) | [−15.401, −15.009] |
| 45 | 1.1506e−231 (Chekmachev) | 1.7885e−231 (ChekAudit, N = 1200) | [−15.138, −14.697] |
Uncertified Galerkin values (λ_N ≥ λ, so lower values of C): −15.62 (8), −15.55 (12, 14), −15.40 (18), −15.45 (20–24), −15.61 (25),
−15.41 (26, 28), −15.35 (30), −15.44 (35), −15.45 (40) (μ = 45, 50 at N ≈ 0.8 × 2.7μ log μ: −15.62, −15.89, too coarse for λ, fine for κ);
they move with μ by ±0.1, like κ. Fits to these 15 values (μ = 8–40):
log λ_N + 4πμ = 4.88 log μ + 15.83 (max residual 0.16); free slope 12.556 (4π = 12.566), exponent 4.68; exponent 5 fixed: C = −15.48 ± 0.08. Chekmachev's published ε_c give −14.95 (23),
−15.00 (25), −14.93 (29), −14.77 (47); his ε_c keep a margin (our certificate at his node count gives 9.4 % more at μ = 5), so MoriRate's
−14.9 ± 0.4 was biased upward. Verdict: C(μ) ∈ [−15.56, −14.70] on 5 ≤ μ ≤ 45, best estimate −15.3 ± 0.2 and flat on 5 ≤ μ ≤ 40; convergence and
the exponent 5 (rather than 4.8–5.0) are not decided. L-functions: not done — ExploreAntenna's conductor-5 forms (L(s, χ₅) for a complex odd
character, Davenport–Heilbronn) are Galerkin-only and unconverged at μ_eff = c/5 ≤ 9.6, and the χ₅ control has no e^{−4πc/5} depth at all
(λ_min = 5.5e−11 at c = 48), so no constant can be compared until that is explained.

## Review
reviewer (subagent, argument pasted inline): Proposition 1 GAP for general Z (a deletion can destroy uniqueness: ordinates ±(n + ½)π with
multiplicities (n + 1)², a = 1), CORRECT under the uniqueness hypothesis now stated; "MoriRate's sufficient target is false" GAP as an
asymptotic claim (finite-μ values do not refute an o(1) estimate; now "fails where tested", asymptotic version open); certificate changes
(Gram linear in ε plus an exact second-order ball, ε-free weights from fixed certified d⁰_j, extrapolated ε then a rigorous check) CORRECT
for 0 ≤ ε < d_min, ε ≤ ¼ (here ε ≤ 1e−17, d_min ≥ 1.7e−3). Significance ROUTINE; weakest point: promoting finite-μ discrepancies to
asymptotic statements (the first line keeps them finite-μ). Second reviewer (tailritz + finite-μ claims pasted inline): tailritz.py upper
bound CORRECT (any exact trial vector; orthonormal basis {x/|x|} ∪ tail); single-deletion falsification of the exact equality at μ = 10
CORRECT given zero balls in the deleted columns (used) and endpoint margins; "super-additive" needed certified single costs, not Σ4g₂ —
now stated only for {γ₁, γ₁₄, γ₃₅}, whose three singles are certified; wording "exact finite-μ equality falsified, asymptotic law not
refuted" adopted. Not reviewed: the a-posteriori Woodbury mode (same ε as the standard mode at μ = 5, 10; used for no reported number).

## Novelty
.scratch/lst.py 'Weil|quadratic form|smallest eigenvalue|Mori|lost sample|sampling.*zero|Christoffel|Paley-Wiener' (arXiv NT, CA, CV, SP,
PR, FA, CO, 2026-10/11) and Zenodo "Weil quadratic form": 33 + 156 titles; abstracts checked by a scout for 2610.32704 (Suzuki), 2610.19132
(Bondarenko–Ortega-Cerdà–Seip, point-evaluation extremals in PW), Zenodo 24144544, 24126058, 24120157: no deletion law, no error term of a
deletion cost, no Green-potential rate of 1 − θ, no certified model-form brackets (24126058 has a λ⁵ prefactor in a prolate-leakage bound).
New here: certified model-form brackets, κ(μ) data, certified finite-μ super-additivity and resolved-pair gain, Proposition 1. All ROUTINE
(reviewer); nothing for DISCOVERIES; RECORDS: certified ε at Chekmachev's node count beats his printed ε₅, ε₁₀, ε₁₆ by 9–10 % (his margin
convention), not a record.

RH relevance: N1 (rate of Weil positivity; the lost-sample law behind MoriRate's barrier survives certified tests at μ = 10, 12, 16, its
error term measured to μ = 50).

## Files
Riemann/Deep/LostLaw (copied from .scratch/lostlaw after the runs): gal.py, certmulti.py, tailritz.py, theta_rig.py, zhp.py,
galerkin_fast.py + galerkin_x.py (ChekAudit copies), certlow.py (MoriRate copy, imported for Sym/scan), zeros_dps179.txt (MoriRate's
zeros), sets_*.json (model forms), chain2–8.sh (the sequential runs as executed, with the RSS watchdog; paths point to .scratch/lostlaw),
logs/ (one JSON log + /usr/bin/time per run; chain.status; scan caches), REPORT.md, DERIVATION.md. Scratch state: .scratch/lostlaw/STATE.md.
Rerun from this directory: uv run --no-project --with mpmath --with python-flint python3 gal.py 16 520 512 2.3 sets_c16.json rig
(36 s, 230 MB); K2MUL=3 … certmulti.py 10 416 sets_c10.json logs/gal_c10_N520.log (196 s, 228 MB); tailritz.py 10 520 1000 416 (37 s,
222 MB); theta_rig.py 16 300 512 1,10,30,69 (5 s, 89 MB).
