# The Li-coefficient remainder: sharp mean-square order, an Ω-bound and a sharpened Li criterion

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** Let λ_n = Σ_ρ [1 − (1 − 1/ρ)^n] (Li: RH ⟺ λ_n > 0 for all n) and R(n) = λ_n − (n/2)(log n − 1 − log 2π + γ).
1. Under RH, Σ_{n≤N} R(n)² ≥ (c₁ − o(1)) N^{3/2}(log N)^{1/2} with c₁ = (3/5)^{1/2}(4/5)²/(2π√2) = 0.0558 (RH only: no pair correlation, no simplicity), and Σ_{N<n≤2N} R(n)² ≪ N^{3/2}(log N)^{3/2}; with Montgomery's pair-correlation theorem the upper bound is N^{3/2}(log N)^{1/2} log log N.
2. Unconditionally R(n) = Ω((n log n)^{1/4}). Under RH limsup |R(n)|/(n log n)^{1/4} ≥ 0.289; if RH fails, |R(n)| ≥ e^{cn} infinitely often.
3. RH ⟺ Σ_m m^δ y_m² < ∞ for some (equivalently every) δ ∈ [0, 1/2), where Keiper's λ_m/m = ½(log m + γ − log 2π − 1) + y_m. Under RH the weight m^{1/2} is never summable.
4. Under RH the smoothed mean of R(n) is exactly 7/4 (Voros's heuristic constant 2Z(0), proved): Σ φ(n/N) R(n) = (7/4) Σ φ(n/N) + O(log N) for φ ∈ C_c^∞(0, ∞).
5. Hence Merlini–Sala–Sala's heuristic conjecture |λ_tiny(n)| ≤ a·log n (arXiv 1903.06233, "better O(n^ε)") is false for every a, and O(n^ε) fails for every ε < 1/4. Numerically a = 2 first fails at n = 5083.

**Before.** Under RH only the upper bound R(n) = O(√n log n) was known: Lagarias (math/0404394, Thm 1.1, automorphic L; he remarks that the arithmetic part may be smaller) and Mazhouda–Voros (1506.01755, Thm 4.1). Known equivalences: RH ⟺ R(n) = o(n) (Voros, math/0506326) and RH ⟺ (y_m) ∈ ℓ² (Arias de Reyna, Funct. Approx. 45, 2011). No lower bound or mean-square law was known. The mean 7/4 was known only heuristically: Voros (math/0404213 l.275–282, math/0506326 l.444–450, from his math/0104051 Eq. (41)) noted that the correction δλ_n = 2Z(0) = +7/4 from the next pole of his Mellin integral improves the fit to the data in the mean, but expected it not to count asymptotically; the trend's t(n) + 3/4 − Σ B_{2k}n^{1−2k}/(4k) is his Eq. (TS), unconditional. Our R436 read the residual as a mean-zero oscillation.

**Method.**
- Lemma A (RH): an explicit zero-sum formula for λ_n with the exact smooth constant 7/4: an endpoint term 2 plus −1/4 from the −π/8 constant of θ (the Fejér kernel gives t(n) + O(1/n)).
- Thm 1: a two-sided weight φ = ψ*ψ with nonnegative Fourier series, λ_{−n} = λ_n, and Cauchy–Schwarz against W(n) = Σ_{γ≤T} cos nθ_γ. Every off-diagonal pair of zeros then counts with the favourable sign; T² log(T/2π) = (3/5)N∫φ is optimal.
- Thm 2: Montgomery's theorem below √(N/log N); Schur's test and Selberg's ∫S² above.
- Barrier: the pointwise exponent lies in [1/4, 1/2]. Lowering 1/2 needs cancellation in Σ_{γ≤√n} e^{inθ_γ}, whose phases are ≈ n/γ ≥ √n, and RH gives none (a von Koch analogue).

**Numerics.** λ_n for n ≤ 10⁵ by a ball-arithmetic DFT of log ξ(1/(1 − z)) on |z| = 0.999 (360 bits, error bound 5.5·10⁻⁵¹) and on 0.998. λ_100000 = 462580.782406902231409416038, equal to Johansson (1309.2877) to 27 digits.
- Σ_{n≤N} R²/(c₁N^{3/2}√log N) = 2.12, 1.96, 2.46 at N = 10³, 10⁴, 10⁵.
- The smoothed mean of R on (N/2, N) is 1.7443 (N = 10⁴), 1.7512 (N = 3·10⁴) and 1.7500 (N = 10⁵).
- max |R|/(n log n)^{1/4} = 1.558, at n = 79302.
- The prime-side line −π^{−1/2}Λ(m)m^{−1/2}(log m)^{−3/4} n^{1/4} cos(2√(n log m) − 3π/4) at m = 2 has amplitude 0.1820 (predicted 0.1820) and phase 0.787 (π/4 = 0.785).

**Checks.**
- The orchestrator reran fft_lambda.py (r = 0.999, 114 MB): all 1,991 sample values are identical, and λ_100000 is as above. From the table it recomputed R(n) at n = 10², …, 10⁵, the ratios, the smoothed mean and the maximum.
- Checked in Merlini–Sala–Sala 1903.06233 (PDF, §2 and conclusion): "we may strength the conjecture to |λ tiny(n)| = o(n) = a·log(n) with a~2"; "a more strong but more gratuitous conjecture that the tiny oscillations may be bounded in absolute value by a·log(n) (better O(n^ε))".
- R451 falsification: Lemma A's constant recomputed by independent quadrature (λ_n^{(0)} − t(n) + 1/4 + 1/(24n) = 8.4·10⁻⁵ at n = 10, 3.8·10⁻⁸ at n = 50; ∫r = −1/8); the trend equals t(n) + 3/4 − 1/(24n) + 1/(240n³) + … to n = 400 (closed form = derivative formula at n = 1, 2, 3); Theorem 4 holds for bumps on (0.1, 1) and (0.9, 1) and for x⁵·bump(0.3, 1) (worst deviation 0.018 at N = 10⁵, inside the O(log N) allowance), and λ_n − λ̄_n has smoothed mean 1.

**Review.**
- Agent's reviewers: Thm 1 CORRECT; a Lemma A gap was repaired (exact c₀ = −1/4); the Thm 3 constant was fixed; Thm 4 and Cor 4 CORRECT.
- Orchestrator's reviewer (statements and proof idea pasted inline; it read proofs.md): Thm 1 with Lemma A, Cor 1, Cor 3 and Thm 4 all CORRECT. Verdict SIGNIFICANT: "substantive Ω-result and sharp weighted-ℓ² threshold, not merely routine manipulation". The two mean-square bounds differ by log log N, so the exponent is sharp but the exact order is not fixed.

- Line-by-line audit (liaudit, 12 Nov 2026): AUDIT PASSED — complete proofs of Lemma 0, Lemma A, Theorems 1–4, Corollaries 1–4, the Ω-result and the weighted Li criterion ([FULLPROOF](../Riemann/Deep/LiAudit/FULLPROOF.md)); sketch steps filled in (the trend lemma by a distributional Fourier transform and Poisson, the order of limits in Theorem 1, the log singularity in Lemma A, Montgomery's weight in Theorem 2, the |G_N″| bound in Theorem 4); six documentation slips, fixed here and in proofs.md; numerics: λ_n for n ≤ 10⁵ recomputed from Odlyzko's zeros with a tail bound (agrees with the FFT to 4.9·10⁻⁶), Theorem 1's ratio 1.472–2.694 over every 100 ≤ N ≤ 10⁵ (orchestrator rerun of lemmaA.py identical, 38 MB) ([LiAudit](../Riemann/Deep/LiAudit/REPORT.md)).
**Novelty.**
- arXiv API (Li coefficients, Keiper–Li, Li criterion, with mean, square, Omega, oscillation, remainder, fluctuation; 40 abstracts) and the Oct–Nov 2026 monthly listings: nothing.
- OpenAlex: the 7 works citing Arias de Reyna 2011 do not treat the size of the remainder.
- HAL (Li coefficients, Keiper–Li): Mazhouda's HDR, Maślanka's analytic extension, Jonot and Voros; none has these results.
- Zenodo "Li coefficients" (69 records): nothing. GXL evidence API (agent): no support.
- Sources read by the agent: Lagarias, Mazhouda–Voros, Voros, Coffey (per-term Laguerre asymptotic, |S_γ(n) + n| = O(n^{1/4}), a conjecture |S_Λ| = O(n^{1/2+ε})), Maślanka, Johansson, 0912.2390 (conjectures O(n^{1/2+ε})), Sekatskii, Odžak–Smajlović and Merlini–Sala–Sala. None has an Ω-result, a mean-square bound or the mean 7/4.
- Our files: TREE N3 (R430, R436) listed only O(√n log n).
- R451 (adversarial, TeX): Voros math/0404213, math/0506326 and 1403.4558, Lagarias math/0404394, Maślanka math/0402168, Flajolet–Vepstas math/0611332 (Coffey's inequality only) and Coffey math-ph/0505052 (PDF) have no mean-square law, no Ω-result and no proof of the mean; Voros has the heuristic 7/4 (Before); Arias de Reyna's ℓ² equivalence is his Thm 1, Cor 2, Thm 6.

**Caveat.**
- Pointwise only (0.202 + o(1))·√n log n/√(log log n) is proved, under RH.
- The constant of the mean square (2√2/(3π) under GUE pair correlation) is a heuristic prediction.
- Positivity of λ_n is not advanced. The pointwise-exponent 'barrier' of proofs.md is a heuristic remark; proved is the range [1/4, 1/2] for the pointwise exponent under RH.
- Lean: Theorem 1 and Corollaries 1–4 in Lean 4 + Mathlib (LiRate.thm1_of_lemmas, li_cor1–4, li_rh_iff; inputs as named hypotheses: Li's formula, Riemann–von Mangoldt, Theorem 2's bound, ¬RH ⇒ radius < 1, the Fourier estimates of the weights; sorry-free, standard axioms; reviewer: faithful with caveats; [LEAN.md](../Riemann/Deep/LiRate/LEAN.md)).

**RH relevance.** N3, the rate of Li's criterion: the sharp mean-square exponent of the remainder, an unconditional Ω-bound, and the sharp weighted form of Arias de Reyna's equivalence.

**Files.** [LiRate report](../Riemann/Deep/LiRate/REPORT.md); in Riemann/Deep/LiRate/: proofs.md, derivation.md, fft_lambda.py, analyze.py, tiny.py, logs/.
