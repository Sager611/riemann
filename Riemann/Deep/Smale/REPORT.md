NO DISCOVERY: the disproof is sound and only existential. No false step was found, the (9/16, 3/8) certificate checks out exactly, and the earlier disproof (GPT-6 Astra's Lean proof, github.com/tadamcz/mean-value-problem) is also non-explicit. No explicit P with S(P) > 1 was found: random-start optimisation never beats (d−1)/d for d ≤ 48, and our potential-theory estimates put every construction from known pairs at degree ≳ 10⁴.

# Smale's mean value constant (arXiv 2610.21472): audit and counterexample hunt (2026-10-22)
Details: [NOTE.md](NOTE.md). Scripts: smale_exact.py (sympy, 4 s), smale_numerics.py (numpy/scipy, each mode < 4 s, < 120 MB). No Lean.

## (0) The earlier disproof
- **Source:** the paper's [Lean] reference is github.com/tadamcz/mean-value-problem (Sept 2026), written by GPT-6 Astra in an Epoch AI run and checked by Comparator. It is a Lean 4 proof that some P has |P(c)/c| > 1 at every critical point; the model uses the lemniscate {|z⁴+64i| ≤ √4097}, z₀ = 33/32+2i and a Blaschke factor.
- **Not explicit:** r, the base polynomial and n are all existential (`disk_exclusion`, `exists_base_polynomial`, `∀ᶠ n`); the README says "very large, unspecified degree". Its approximation rate is only e^{−3·10⁻⁵ N}, from the Green value (1/8)log(4097/4096) at the conformal map's branch points.
- **Nothing else** turned up on arXiv (API), openai/math (749 preprints) or Zenodo.

## (1) Audit of 𝓜 ≥ 𝓑 and of the pair (9/16, 3/8): no false or unjustified step
- **The chain 𝓑 ≤ 𝓐 ≤ 𝓜 holds:** regularise, take P̃ = wB(W⁻¹w)/B(0), approximate by Runge, then use P_n = (p + w qⁿ)/(1+q(0)ⁿ) with a Hilbert lemniscate; the critical points are trapped in the collar, where |qⁿ| = O(1/n). Only implicit details are missing: Lemma 2.1 also needs t > B(0), and (T5) needs a single-curve lemniscate (standard via Faber polynomials).
- **The pair checks out exactly:** H matches the paper and 1/B(0)² = 8081/6656. The Schur–Cohn constants δ₁, δ₂, δ₃ (degrees 6, 12, 24) and |lead T₁|² have no root in [−1,1], |lead T₂|² vanishes only at c = 0, δ₁·2³² = D₁, and the values at c = 0, ±1/√2 equal the paper's. The margin is thin: min |root of P_c| = 1.0022, and H' has a zero at |z| = 1.0039.
- **Admissible but not regular:** the model has a critical point on ∂Ω, at w* = H(1)/N, with |f/w| = 1/B(0); Lemma 2.1 covers this. 𝓜 > 1 is independently confirmed by the Lean proof.

## (2)–(3) Explicit P: none
- **Direct optimisation** (SLP on the critical points, exact quadrature): best S = (d−1)/d for d = 4–12; 0.9375, 0.9500, 0.9526 and 0.9655 for d = 16, 24, 32, 48.
- **Exact-certificate family** P = wR(w^k), P' = T(w^k): the critical points are known exactly and S(P) = min|R(t_l)|. The k → ∞ limit is Smale's problem for u·e^{L(u)}, impossible for deg L = 2; optima found are Smale's extremal for N ≤ 20 orbits.
- **Other searches, all negative:** relaxed pairs on a disk with polynomial f of degree ≤ 16 (max 1.000); Astra-type lemniscates |zⁿ−1| < ρ with n ≤ 6 (no admissible pair with A ≥ 1.003); polynomials of degree 4–6 with mountain-pass level R*(p) > 1 (none).
- **Side result (numerical, uncertified):** a new degree-one pair. W = e^{3.1z}, b = −0.1814+0.6221i, zero a = −0.0555−0.8659i gives 1/B(0) = 1.0316, above the paper's discarded 1.0226 and 1.0255.
- **Cost of the indirect route:** all critical points on one curve need n ≥ max(−σ/ω). That is 9.3·10⁶ for the paper's pair at 1.01 on Γ (its base point sits at the centre of a spiral, where the equilibrium density is 2·10⁻⁶) and 4.7·10⁴ for the e^{3.1z} pair at 1.003.
- **Answer to (3):** the indirect argument does give polynomials for all large degrees, and Lean confirms 𝓜 > 1, but neither gives an explicit one. A degree-independent certificate exists (Rouché on two curves plus a collar bound, NOTE §4); what is missing is a pair (p, q) of moderate degree with positive margins.

Suggested TREE line: "Smale mean-value (Orca): 2610.21472 sound (exact Schur–Cohn re-check); earlier disproof = Astra Lean (non-explicit); no explicit P (S ≤ (d−1)/d to d = 48; known pairs need degree ≳ 10⁴)."
