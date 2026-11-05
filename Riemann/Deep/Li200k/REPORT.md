COROLLARY: with the zero-free half-plane Re s > 7/8, lambda_n > 0 for all n <= 952,080,909,647,480,761,027,133,440 ≈ 9.52·10^26 = 105.77 T^2 (T = 3,000,175,332,800); routine, the a = 1/2 range is KNOWN (Fabbian, Zenodo 22839423: 77.03 T^2 certified, 78.72 T^2 continuation)

# li200k report (node N3, Li / Keiper–Li coefficients), 2026-11-05

## Short answer

- **Known part (checked ourselves).** RH verified up to T = 3,000,175,332,800 (Platt–Trudgian) already gives
  λ_n > 0 far beyond the direct-computation record 10^5 (Johansson). Fabbian (Zenodo 10.5281/zenodo.22839423,
  19 Sep 2026, unrefereed) proves it explicitly: Theorem 1, n ≤ H²(201/200)^871 = 693,318,936,098,417,738,514,687,162
  (77.03 T²); Remark B.7, a single-implementation continuation to 708,581,529,455,821,832,027,338,082 (78.72 T²).
  Oesterlé (2000, unpublished) and Brown (2005, Thm 2, proof has a gap) stated T²-scale ranges earlier.
- **Our independent proof and certificate** of the same statement (zeros above T only known to lie in 0 < β < 1)
  gives 708,505,587,700,793,981,366,435,840 = 78.714 T², which matches Fabbian's continuation to 0.011%. It is a second
  implementation of his number, not a new result.
- **Possibly new, routine: the 7/8 half-plane.** If zeros above T satisfy |β − 1/2| ≤ 3/8 (OpenAI's half-plane
  Re s > 7/8, unrefereed, Lean proof builds in Riemann/Deep/QRH), the same argument certifies
  **λ_n > 0 for n ≤ 952,080,909,647,480,761,027,133,440 = 105.774 T²**: 1.373 times Fabbian's Theorem 1 and 1.344 times
  his continuation. The range scales like 1/a for strip half-width a, so 4/3 was expected. Nobody has an a = 3/8 range
  (search below); Fabbian's §8.2 only says that such information "would extend the range".
- **RH relevance (N3): none beyond the inputs.** Finite ranges come from the verified height plus the strip width and
  carry no information on zeros above T (Fabbian, Prop. 8.1). The 2·10^5 direct computation of the brief was skipped.

## Theorem

Let a(t) ∈ (0, 1/2] with a(t)/t² and a(t)²/t³ nonincreasing, such that every zero ρ = β + iγ with γ ≥ T has
|β − 1/2| ≤ a(γ), and suppose β = 1/2 for 0 < γ ≤ T. Then λ_n > 0 for 1 ≤ n ≤ N(a):

| a | input on zeros above T | N(a) | N/T² | status |
|---|---|---|---|---|
| 1/2 | 0 < β < 1 | 708,505,587,700,793,981,366,435,840 | 78.714 | known (Fabbian: 77.026 thm, 78.722 continuation) |
| 1/2 − 1/(5.558691 log t) | Mossinghoff–Trudgian–Yang zero-free region (refereed) | 717,734,280,405,823,147,676,270,592 | 79.739 | routine, +1.3% over 78.72: insignificant |
| **3/8** | **half-plane Re s > 7/8** | **952,080,909,647,480,761,027,133,440** | **105.774** | **the corollary** |

Conditional half-planes (verify.py --conditional): Liu 0.874957… → 105.787 T², Chang–Wang 0.87446205769353 → 105.930 T².

## Proof

Notation: λ_n = Σ_ρ [1 − (1 − 1/ρ)^n] (symmetric summation; Li 1997, Bombieri–Lagarias 1999), g(x) = 2 − 2cos x,
θ(t) = 2 arctan(1/(2t)), N(t) = #{ρ: 0 < γ < t} with multiplicity, Ñ(t) = (t/2π) log(t/2πe) + 7/8,
|N − Ñ| ≤ E(t) = 0.111 log t + 0.275 log log t + 2.450 + 0.2/t for t ≥ e (Trudgian 2014, Thm 1 and (3.11)).

1. **Per zero.** An on-line zero 1/2 + iγ and its conjugate contribute g(nθ(γ)) ≥ 0. An off-line zero ρ = 1/2 + δ + iγ
   comes with ρ* = 1 − ρ̄ = 1/2 − δ + iγ; z_ρ = re^{iφ}, z_{ρ*} = r^{−1}e^{iφ}. The quartet ρ, ρ*, ρ̄, ρ̄* contributes
   4 − 2(r^n + r^{−n}) cos nφ = 2 × [2 − (r^n + r^{−n}) cos nφ] ≥ 2 × [g(nθ(γ)) − q_n(γ)], with
   q_n(t) = 4 sinh²(n a(t)/(2t²)) + 4n a(t)²/t³, using ℓ = −log r ≤ δ/γ² ≤ a/γ², |φ − θ(γ)| ≤ 2δ²/γ³, |g′| ≤ 2.
   **Both ρ and ρ* have height γ > 0, so each carries one g-term and one q-term** (counting the quartet once would
   lose half of the loss and give false counterexamples). Off-line zeros have γ > T, hence
   λ_n ≥ M_n − Q_n, M_n = Σ_{γ>0} g(nθ(γ)) (ordinates only, no RH), Q_n = Σ_{γ≥T} q_n(γ).
2. **Loss.** q_n decreases in t; Stieltjes against N gives Q_n ≤ ∫_T^∞ q_n (Ñ′ + E′) dt + 2q_n(T)E(T), E′ ≤ 0.13/t,
   bounded by an upper Riemann sum on [T, 100T] plus an analytic tail.
3. **Gain** (terms of M_n are ≥ 0, so low zeros may be dropped). For n < 10^4: keep γ ≥ t_n = ½cot(π/2n), where
   G(t) = g(nθ(t)) falls from 4 to 0, and M_n ≥ (n/π)[K₀(log(n/2π) − π²/4n²) − K₁] − 8E(t_n) − (0.111 + 0.275/log t_n)n²/(2t_n²),
   K₀ = Si(π) − 2/π, K₁ = ∫₀^π (1 − cos u) log u · u⁻² du. For n ≥ 10^4: keep γ ≥ √n, cut where nθ(t) ∈ πℤ (G monotone
   on each piece, each piece costs ≤ 8E); the integral is ≥ (n/π)[(log(n/2π) − 1/(4n))(π/2 − 2/U) − (π/2)(1 − γ_E)],
   U = nθ(√n), which is the RH asymptotic (n/2)(log n − 1 − log 2π + γ_E) up to lower-order terms.
4. **Small n.** Hardy's Z changes sign on [14.13, 14.14], so a zero γ* there is on the line; for n ≤ 88,
   nθ(γ*) ∈ (0.0706, 6.23) and λ_n ≥ g(nθ(γ*)) − Q_n ≥ 0.00335 − 2·10^−23.
5. **Chain.** On integer blocks [n₁, n₂] the gain is taken at its worst end and Q_n is increasing in n; verify.py
   certifies gain > loss on about 95 contiguous blocks from 89 to N(a), every comparison as a ball sign. ∎

Reviewer subagent (verdict only, two rounds): Steps 1–4 correct; the first round flagged that blocks below 10^4 needed
bound 3a (omitted from the description), the second round: "CORRECT as a certification scheme".

## Inputs

| input | used for | status |
|---|---|---|
| Platt–Trudgian, Bull. LMS 53 (2021) | RH for 0 < γ ≤ 3,000,175,332,800 | refereed |
| Trudgian, JNT 134 (2014), arXiv 1208.5846 | E(t) | refereed |
| OpenAI, "The Quasi-Riemann Hypothesis: a zero-free half-plane Re s > 7/8" (30 Sep 2026) | a = 3/8 | unrefereed; Lean proof builds here |
| Mossinghoff–Trudgian–Yang, arXiv 2212.06867 | a(t) of the side line | refereed |

## Certificate (Riemann/Deep/Li200k/verify.py; python-flint 0.9.0, FLINT 3.6; 39 s, peak 36 MB, one thread)

Part 0: the per-zero inequality of step 1 on 3000 random off-line zeros (min slack 0.0037). Part 1: the zero at
14.13–14.14 and n ≤ 88. Part 3: λ_n for n ≤ 600 computed directly (λ₁ = 0.0230957089661, as Johansson) lies above
our bound (ratio ≤ 0.84); at n = 10^5 our bound 457,145.2 is below Johansson's λ = 462,580.78. Part 2: the block
chains for the three rows of the table.

## Novelty check

Zenodo API ("Li coefficients", "Keiper-Li", "Li's criterion", with 7/8 / half-plane / zero-free; all Fabbian records);
Fabbian 22839423 and 23027727 full text (no 7/8, no 3/8, no strip computation); Polson–Sokolov, "Computing the Riemann
Hypothesis" Vol. I–II (Zenodo 23891944/23891946, 2 Nov 2026): chapter 23 lists the consequences of the 7/8 half-plane,
none for Li coefficients, and §I.23.6 treats λ_n only for n ≤ 100; arXiv math.NT and math.CV listings 2026-10 and
2026-11 (no Li-coefficient titles); GXL evidence API (2 claims: only the T²-scale statements of Oesterlé/Voros/Maślanka,
nothing with a strip). Not found elsewhere.

## Classification and suggested repo edits (orchestrator)

- Routine deduction from the 7/8 half-plane → COROLLARIES.md: "λ_n > 0 for n ≤ 9.52·10^26 = 105.77 T² (7/8 half-plane)".
- TREE N3 record: "λ_n > 0 for n ≤ 6.93·10^26 (Fabbian 2026, explicit proof from RH to T; 7.09·10^26 continuation,
  reproduced here); with 7/8: 9.52·10^26 (ours); direct computation 10^5 (Johansson)". Barrier: ranges ≈ (T²/a)·log T;
  nothing above T. R282 ("RH-to-height statements do not give a larger unconditional range") and NOGO R118 are
  superseded.

## Files

- Riemann/Deep/Li200k/REPORT.md — this report.
- Riemann/Deep/Li200k/verify.py — certificate (39 s).
