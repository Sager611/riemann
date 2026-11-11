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

## (ii') Two-sided first-order law: what is proved (follow-up, 11 Nov 2026)

Setting. Even sector (odd: sin in place of cos). Z ⊂ ℝ a symmetric locally finite multiset (zeros of ζ under RH, or a model),
Q(u) = Σ_{γ∈Z} m_γ|û(γ)|² on L²(−a, a), a = ½ ln μ, û(t) = ∫u e^{itx}; A = its Friedrichs operator, λ = λ_min(A) > 0.
For γ₀ ∈ Z, e₀ = cos(γ₀x)·1_{(−a,a)}, so û(γ₀) = ⟨u, e₀⟩ for real even u; deleted form Q′ = Q − 2⟨·, e₀⟩², λ′ = λ_min(Q′).
y₀ = γ₀/2πμ, g₂(y) = log((2 + √(4 − y²))/|y|) (|y| < 2), 0 (|y| ≥ 2). Conjectured law: log(λ/λ′) = 4g₂(y₀) + o(1).

Lemma A (exact, every Z, every μ).
(a) 0 ≤ λ′ ≤ λ, since Q′ = Σ_{Z∖{±γ₀}} m|û|² ≥ 0 and Q′ ≤ Q.
(b) Let θ = sup_u 2⟨u, e₀⟩²/Q(u) = 2⟨e₀, A⁻¹e₀⟩ ∈ (0, 1]. Then λ′ ≥ (1 − θ)λ. Proof: Q′(u) = Q(u) − 2⟨u, e₀⟩² ≥ (1 − θ)Q(u) ≥ (1 − θ)λ‖u‖².
(c) λ′ is the least x ∈ [0, λ) with 2⟨e₀, (A − x)⁻¹e₀⟩ = 1 if such x exists, and λ′ = λ otherwise (Weinstein–Aronszajn for a
    rank-one decrement; the resolvent exists only below λ).
(d) 1 − θ = inf_u Q′(u)/Q(u) = 1/(1 + 2C), C = sup_u ⟨u, e₀⟩²/Q′(u): the share of zero-side energy that must stay off ±γ₀.
So the upper half of the law, log(λ/λ′) ≤ 4g₂(y₀) + o(1), follows from the scale-free estimate 1 − θ ≥ e^{−4g₂(y₀) − o(1)}; the lower
half, log(λ/λ′) ≥ 4g₂(y₀) − o(1), needs test vectors u with Q′(u) ≤ e^{−4g₂(y₀) + o(1)} λ‖u‖².

Lemma B (the multiplier proxy, two-sided; rigorous). Mori's Theorem 4 problem: P(c) = inf{p(σ) : σ ≥ c, ∫dσ ≤ 0,
∫(1 + |log|y||)d|σ| < ∞}, p(σ) = sup_x U_σ(x) − U_σ(0), c₀ = max(−log|y|, −b)dy, b ∈ [log 2, ∞], P(c₀) = 2π (Mori). A lost sample
at ±y₀ is the constraint c₀ + Δ, Δ = μ⁻¹(δ_{y₀} + δ_{−y₀}). For 0 < y₀ < 2 fixed:
    2π + 2g₂(y₀)/μ ≤ P(c₀ + Δ) ≤ 2π + 2g₂(y₀)/μ + K(y₀) μ^{−3/2} log μ   (μ ≥ μ₀(y₀)).
Lower (all μ): Mori's certificate. ω₂ = harmonic measure of Ω = ℂ ∖ ((−∞, −2] ∪ [2, ∞)) at 0, g₂ = U_{ω₂} − log|y| (Mori Lemma 3.2).
For admissible σ: p(σ) ≥ ∫(U_σ − U_σ(0))dω₂ = ∫g₂ dσ ≥ ∫g₂ dc₀ + ∫g₂ dΔ = 2π + 2g₂(y₀)/μ (g₂ ≥ 0, σ − c₀ − Δ ≥ 0, and
∫g₂(−log|y|)dy = 2π; on |y| < 2, max(−log|y|, −b) = −log|y| because b ≥ log 2).
Upper: β_y = harmonic measure of Ω at y (a probability measure on K = ℝ∖(−2, 2) with finite log-moment, density ≍ (|t| − 2)^{−1/2} at
the tips); as in Mori's Lemma 3.2, g_Ω(x, y) = U_{β_y}(x) − log|x − y| (both sides harmonic on Ω∖{y}, same logarithmic pole, 0 on K,
→ 0 at ∞). Put σ_Δ = σ* + μ⁻¹Σ_± (δ_{±y₀} − β_{±y₀}). Then U_{σ_Δ}(x) − U_{σ_Δ}(0) = [U_{σ*}(x) − U_{σ*}(0)] − μ⁻¹Σ_±[g_Ω(x, ±y₀) −
g_Ω(0, ±y₀)] ≤ 2π + 2g₂(y₀)/μ, with equality on K (Mori Thm 4(ii): U_{σ*} − U_{σ*}(0) ≤ 2π, = 2π on K; g_Ω ≥ 0, = 0 on K;
g_Ω(0, ±y₀) = g₂(y₀) by symmetry). The balayage identity U_{β_y} = log|· − y| on K (no constant) was also checked numerically by
the disc map w = (z/2)/(1 + √(1 − z²/4)) and the Poisson kernel (balayage_check.py: equal to 15 digits at y = 0.5, 1.3). ∫dσ_Δ = ∫dσ* ≤ 0. The constraint holds on (−2, 2) (σ_Δ = −log|y| + the atoms) and on K except where
σ*'s slack arccosh(|y|/2) ≈ √(|y| − 2) is below μ⁻¹ dβ/dy ≍ μ⁻¹(|y| − 2)^{−1/2}, i.e. on |y| − 2 ≤ c(y₀)/μ, of mass O(μ^{−3/2});
add that mass back there and remove the same mass from 3 ≤ |y| ≤ 4 (slack ≥ 0.55 for every b ≥ log 2). The combined correction has
mass zero, so its potential is O(μ^{−3/2} log μ) uniformly on ℝ (each piece alone grows like m·log|x|), and p rises by at most
K μ^{−3/2} log μ. For finitely many fixed interior points y₁, …, y_m the construction is linear: P = 2π + (2/μ)Σ g₂(y_i) +
O_{y₁..y_m}(μ^{−3/2} log μ) — additive, not only super-additive (constants depend on the locations). Times Mori's scale
2R_ϑ = 2μe^{2/s}: the proxy cost of one lost sample is 4g₂(y₀)(1 + O(1/log μ)) + O(μ^{−1/2} log μ). This closes the reviewer's gap (1) for the proxy (envelope theorem replaced by an explicit
balayage measure), not for λ_min.

Proposition C (λ_min itself: certified two-sided brackets; first zero γ₁ = 14.1347…, enclosed by a certified sign change of Hardy's Z).
Lower bounds are Chekmachev's certificate (Thm 5.10; ChekAudit's re-implementation) with the deleted zero as one more negative rank-one
column √2·e₀ (Q′ − ε ≥ D + UMU* − BB* − 2e₀e₀* and the Woodbury test λ_max(B′*P⁻¹B′) ≤ 1, Lemma 5.5 Cholesky in Arb: certlow.py);
upper bounds are Arb Rayleigh quotients of explicit vectors (wmat.py Galerkin ground states: certdel.py, ubfull.py).
| μ | y₁ | λ (even) | λ′ = λ_min(Q − 2|û(γ₁)|²) | log(λ/λ′) | 4g₂(y₁) |
|---|---|---|---|---|---|
| 5 | 0.4499 | [7.4475e−18, 1.0128e−17] | [1.3883e−21, 1.7545e−21] | [8.354, 8.895] | 8.688 |
| 10 | 0.2250 | [9.1338e−44, 1.7279e−43] | [1.0130e−48, 1.8225e−48] | [10.822, 12.047] | 11.500 |
(λ lower bounds ε₅, ε₁₀ certified by Chekmachev and re-certified in ChekAudit, even sector.) Both brackets contain 4g₂(y₁): at these two
windows |log(λ/λ′) − 4g₂(y₁)| ≤ 0.34 (μ = 5) and ≤ 0.68 (μ = 10), rigorously. Finite μ only.

What remains open. (1) No O(1)-precise asymptotic of λ_min is known even without deletion (log λ_min = −4πμ + o(μ) only under Mori's
hypotheses), and the law is an O(1) statement. (2) Certificate half: Mori's chain λ_min ≥ s·A(Λ) ≥ exp(−2R_ϑP − E), E = O((A_μ + s³)log s).
With Lemma B, an extra virtual node at ±y₀ raises 2R_ϑP by 4g₂(y₀) + o(1), but a node's placement and elimination and the Voss–Hall
factor are controlled only up to O(1)–O(log s), and the bound is relative to Mori's lower bound for Z, not to λ_min(Z). Needed: a
sampling constant sharp to e^{o(1)}, i.e. log λ_min(Z) = −2R_ϑ P(c_Z) + κ(μ) + o(1) with κ insensitive to one point. (3) Test-vector
half: needs u with Q′(u) ≤ e^{−4g₂ + o(1)}λ_min(Z)‖u‖² for all large μ; the CCM/Mori prolate vector is zero-blind (its zero-side energy
sits at y ≈ 2.5–3, so deleting a central zero changes its Rayleigh quotient by a negligible amount), and no zero-adaptive test vector
with controlled Rayleigh quotient is known. Sufficient (not equivalent) target for the upper half via Lemma A(b), (d): the scale-free extremal value 1 − θ(γ₀) ≥ e^{−4g₂ − o(1)}.

Review of (ii′) (reviewer subagent, argument pasted inline): Lemma A GAP, fixed above ((c): λ′ = λ when the secular equation has no root
below λ). Lemma B: lower bound and no-constant balayage identity CORRECT; repairable gaps fixed above (zero-mass tip correction, slack
0.55 for b = log 2, location-dependent constants, O(g₂/s) from 2R_ϑ). Proposition C: CORRECT as a technique (implementation details:
Arb balls throughout, zero ball radius 10⁻⁵⁶ (μ = 5), 10⁻¹⁰⁴ (μ = 10), ChekAudit's Lemma 5.5 Cholesky check). Significance: ROUTINE as
progress toward the spectral law; the explicit proxy sensitivity is useful. 'What remains open': accurate. Weakest point: no
deletion-stable, e^{o(1)}-sharp bridge from the multiplier optimum to the spectral bottom.
