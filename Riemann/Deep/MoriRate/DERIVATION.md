# MoriRate — derivation of the lost-sample law (heuristic, first-order LP sensitivity of Mori's Theorem 4)

Setting (Mori 24116425 v4, Lemma 2.1; explore19 wmat.py has the same normalisation). Window [λ⁻¹, λ], a = log λ, μ = λ², s = log μ.
Under RH, for real u ∈ C_c^∞(−a, a): W(u) = Σ_ρ m_ρ |û(γ_ρ)|², û of exponential type a. So λ_min is the lower sampling
constant (times s) of the zero ordinates (both signs, with multiplicity) for PW_a, normalised by ‖u‖₂².

1. Deficit. In v = γ/2π the zero density is log|v| (Riemann–von Mangoldt), the Nyquist density of PW_a is s, so the deficit
   density is (s − log|v|)₊ = log₊(μ/|v|). With y = v/μ it is μ·(−log|y|)₊ dy: about 2μ samples missing (∫_{−1}^{1} −log|y| = 2).
2. Cost (Mori §1.5, Prop. 5.1, Thm 4). Virtual nodes fill the deficit; each is removed by a multiplier Φ_j whose zero density is
   ϑ + σ*(v/R_ϑ), R_ϑ = μ(1 + O(1/s)); cost 2 log sup|Φ_j| = 2 R_ϑ p(σ*) = 4πμ(1 + O(1/s)), with
   P(constraint) = inf{p(σ) = sup_x U_σ(x) − U_σ(0) : σ ≥ constraint, ∫dσ ≤ 0}, P(−log|y| dy) = 2π, extremal σ* (Joukowski).
3. Dual certificate (Mori Thm 4(i), Lemma 3.2). ω₂ = harmonic measure at 0 of Ω = ℂ ∖ ((−∞, −2] ∪ [2, ∞)), g₂ = U_{ω₂} − log|y|
   = Green function of Ω with pole at 0. With z = y/2 and the conformal map z ↦ z/(1 + √(1 − z²)) of the doubly slit plane onto the
   disc, g₂(y) = log((2 + √(4 − y²))/|y|) for |y| < 2 and g₂ = 0 for |y| ≥ 2. p(σ) ≥ ∫ g₂ dσ ≥ ∫ g₂ (−log|y|) dy = 2π
   (mr_numbers.py: the integral equals 2π to 15 digits; ∫_0^2 g₂ = π).
4. A lost sample. Deleting the zero 1/2 + iγ₀ (or merging it into a neighbour: a double zero, or a pair closer than δ*) removes the
   two samples ±γ₀, i.e. adds unit masses at ±y₀ = ±γ₀/(2πμ) to the deficit: constraint −log|y| dy + μ⁻¹(δ_{y₀} + δ_{−y₀}).
   First order (LP envelope, ω₂ stays dual-optimal): P increases by 2 g₂(y₀)/μ, so −log λ_min increases by 2μ · 2g₂(y₀)/μ:

       log(λ_min / λ_min^{(−ρ₀)}) = 4 g₂(γ₀/2πμ) + o(1),   i.e.  λ_min^{(−ρ₀)} = λ_min · (y₀/(2 + √(4 − y₀²)))⁴ · e^{o(1)}.

   Zeros above 4πμ (y₀ ≥ 2) cost nothing at first order (g₂ = 0 on supp ω₂). The factor 4 = 2 (ρ₀ and ρ̄₀) × 2 (energy |F|²).
5. Several lost samples. The value P(constraint) is convex in the constraint (inf of a convex functional over a feasible set that
   moves linearly), so by the certificate P ≥ 2π + μ⁻¹ Σ_i 2g₂(y_i): costs are at least additive in the multiplier problem.
   A band of lost samples above 4πμ that pushes the distinct-point density below s is not a small perturbation (the constraint
   becomes positive there) and costs a positive multiple of the new deficit (exterior collapse; tested at c = 8).
6. Close pair at distance δ. Let F₀ be the ground state of the model with γ₂ deleted (F₀(γ₁) ≈ 0, so it is also the merged
   ground state). Moving γ₂ to γ₁ + δ adds 2F₀(γ₁ + δ)² ≈ 2F₀′(γ₁)²δ²: λ(δ) = λ_merged + Kδ², K ≤ 2F₀′(γ₁)² (coupling to the
   next near-null modes halves it at μ = 8: measured 2.30·10⁻¹³ vs first-order 4.60·10⁻¹³), until Kδ² reaches λ_resolved − λ_merged,
   i.e. δ* = √((λ_resolved − λ_merged)/K) ≈ e^{−3.24μ+3.6} (μ = 6–14, y ≈ 0.75); beyond δ* the pair is resolved, costs nothing, helps.

Consequences (rate). If ζ (under RH) had k(μ) lost samples in |γ| < 4πμ with Σ 4g₂(y_i) ≥ εμ for infinitely many μ, then
−log λ_min ≥ (4π + ε + o(1))μ along them. A fixed proportion p of double zeros below 4πμ gives Σ 4g₂ ≈ 4p ∫_0^2 g₂(y) μ log(μ y) dy
≈ 4πpμ log μ: the rate becomes superlinear. Known multiplicity results (≥ 19/27 simple under RH; ≥ 67.25 % simple and ≥ 83.62 %
distinct unconditionally) allow p ≈ 0.15 (doubles), so they cannot give the limit 4π; even p = ε/log T, invisible to every asymptotic
proportion, raises the rate to ≈ 4π(1 + ε).

Status: steps 2–3 are Mori's theorems; the lower bound P ≥ 2π + μ⁻¹Σ2g₂(y_i) in step 5 is rigorous for the multiplier problem;
steps 4–6 for the spectral minimum are a conjecture (reviewer: no two-sided comparison between multiplier cost and λ_min), supported
by the numerics and the certified instances in REPORT.md.
