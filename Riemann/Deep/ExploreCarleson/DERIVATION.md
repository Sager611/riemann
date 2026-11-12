# explore42 derivation: the continuous Nyman–Beurling distance in a Sonine space, Riesz geometry of the zero evaluators, Carleson constant of shifted zeros

Notation as explore38 ([DERIVATION](../ExploreChristoffel/DERIVATION.md)): K = L²(0,∞), Mellin on Re s = ½ with dτ/2π, (f,g) linear in f; χ₁ = (1 + log t)1_(0,1];
A = inverse Mellin of Z(s) = (s−1)ζ(s)/s², supported in (0,1]; D_θf(t) = f(t/θ); C_λ = closed span{D_θA : λ ≤ θ ≤ 1}; D(λ) = dist(χ₁, C_λ) ≤ d_N for λ = 1/N.
M f(t) = t⁻¹∫₀ᵗf (multiplier 1/(1−s)), M′f(t) = ∫_t^∞ f(u)du/u (multiplier 1/s), (1−M)⁻¹ = (1−M)* = 1−M′ (all unitary on K).
J f(t) = f(1/t)/t (Mellin f̂(1−s)), C f(v) = 2∫₀^∞ cos(2πvt) f(t) dt (unitary, C² = 1). CJ has multiplier 2(2π)^{−s}Γ(s)cos(πs/2) = 1/χ(s).
V = explore38's unitary, multiplier (s/(1−s))³/χ(s), so V = −(1−M)³CJ and V⁻¹ = −JC(1−M′)³. Q_λ = projection onto L²[λ,∞). ℓ := log(λ/v).

## 1. C_λ sits in a Sonine space (Lemma 1)
VA = JA (V = Z̄/Z on the line), supported in [1,∞); V commutes with D_θ, so VD_θA = D_θJA is supported in [θ,∞) ⊂ [λ,∞). Also D_θA ⊂ L²(0,1). Hence
  C_λ ⊂ W_λ := {f ∈ L²(0,1) : Vf = 0 on (0,λ)},
a Sonine-type space (Burnol's L_a with a = N^{−1/2}, up to the finite-dimensional bookkeeping of the (1−M)³ factor). Orthogonally,
  D(λ)² = F(λ) + ‖P_{W_λ⊖C_λ} χ₁‖²,   F(λ) := dist²(χ₁, W_λ)   (F depends only on the gamma factor, not on ζ's zeros).

## 2. The ζ-free floor (Lemma 2)
W_λ^⊥ = closure(L²(1,∞) + V⁻¹L²(0,λ)) = L²(1,∞) ⊕ closure(R), R := P_(0,1)V⁻¹L²(0,λ). Since χ₁ ⊥ L²(1,∞),
  F = sup_{g∈L²(0,λ)} |(χ₁, V⁻¹g)|² / ‖P_(0,1)V⁻¹g‖².
Substitute h = (1−M′)³g. For f ∈ L²(0,λ): (1−M′)f stays in L²(0,λ) while (1−M)f = −t⁻¹∫₀^λ f on (λ,∞); hence h ↦ g = (1−M)³h is unitary from
H₃ := {h ∈ L²(0,λ) : h ⊥ 1, ℓ, ℓ²} onto L²(0,λ) ((1−M′)ℓ^k = ℓ^k − ℓ^{k+1}/(k+1) and (1−M′)t^p = (1 + 1/p)t^p − λ^p/p give the three conditions). Then
  V⁻¹g = −JCh,  ‖P_(0,1)V⁻¹g‖² = ‖P_(1,∞)Ch‖² = ‖h‖² − ‖P_(0,1)Ch‖² = ((I − T)h, h),  (χ₁, V⁻¹g) = ((1−M′)³Vχ₁, h) = −(q, h),
  q := CJχ₁,  q(v) = 2∫₁^∞ cos(2πvt)(1 − log t) dt/t = −2[(1 + log a)Ci(a) + J(a)],  a = 2πv,  J(a) = ∫_a^∞ cos u · log u · du/u,
  T(v,v′) = sin 2π(v−v′)/(π(v−v′)) + sin 2π(v+v′)/(π(v+v′))  on (0,λ)²  (T = P_(0,λ)CP_(0,1)CP_(0,λ): the even prolate concentration operator, band (0,1)).
So F(λ) = (P_{H₃}q, (P_{H₃}(I−T)P_{H₃})⁻¹P_{H₃}q) ≥ ε(λ) := ‖P_{H₃}q‖² = ‖(1−Q_λ)Vχ₁‖² = ∫₀^λ|Σ_j b_j(2πt)^{2j}|²dt (explore38 §2; the two computations agree to 7 digits).
F − ε is the gain from the time–frequency defect of T (top eigenvalue μ₀ = 0.038 at N = 2, 5.2·10⁻⁵ at N = 10). Nyström, graded nodes v = λx³, M = 40 vs 80: 7 digits.

## 3. The zero part (Lemma 3) and the joint bound
Burnol's vectors Y_ρ = V⁻¹Q_λVψ_ρ (ψ_w = t^{−w}1_(0,1], limit w → ρ) lie in C_λ^⊥ but not in W_λ. For f ∈ W_λ, (f, P_WY_ρ) = (Vf, Q_λVψ_ρ) = (f, ψ_ρ):
P_WY_ρ is the evaluator of W_λ at ρ, and P_WY_ρ ⊥ C_λ. By duality, W_λ ⊖ C_λ ⊇ closed span of the evaluators, with equality iff C_λ = {f ∈ W_λ : f vanishes at every zero}
(spectral synthesis). Burnol's Theorem A (math/0203120: for a < 1, L_a = \bar P_a ⊕ Y_a, P_a the co-Poisson subspace of sums of g ∈ L¹(a, 1/a), Y_a spanned by the zero
evaluators) is that statement in his extended Sonine space L_a. After J and dilation by a = λ^{1/2}, W_λ = {F = 0 on (0,a), CF ∈ span{1, log, log²} on (0,a)}, a 3-dimensional
extension of de Branges' K_a (L_a is a 2-dimensional one), so the transfer is finite-dimensional bookkeeping, not done here. Rigorously
  D(λ)² ≥ F(λ) + c_Wᴴ G_W⁻¹ c_W   (evaluators of W_λ at any finite set of zeros, with multiplicity), with equality over all zeros under synthesis.
Computable lower bound without Theorem A: S := W^⊥ + span{Y_ρ : ρ ∈ R} ⊂ C_λ^⊥ contains L²(1,∞), so ‖P_Sχ₁‖² = ‖P_{R + span y}χ₁‖², y_ρ = P_(0,1)Y_ρ.
Write Y_ρ = ψ_ρ + JCh_ρ with h_ρ = (1−M′)³[1_(0,λ)(V(ρ)t^{−ρ} + F_ρ(t))] (F_ρ = Σ a_j(ρ)(2πt)^{2j}, explore38 §2); h_ρ is t^{−ρ} times a constant plus
Σ_j t^{2j}-terms plus a polynomial of degree ≤ 2 in ℓ (exact algebra above), so Ch_ρ(u) = 2Σ_n (−1)^n(2πu)^{2n}/(2n)! ∫₀^λ t^{2n}h_ρ exactly
(∫₀^λ t^m ℓ^k = λ^{m+1}k!/(m+1)^{k+1}). Since Y_ρ ⊥ V⁻¹L²(0,λ):
  (y_j, y_i) = H_ij − ∫₀¹ Ch_j \overline{Ch_i} du,   (y_j, P_(0,1)V⁻¹g) = ∫₀^λ (CP_(0,1)Ch_j)(v) \overline{h(v)} dv,   (χ₁, y_i) = c_i,
with H, c explore38's exact finite-λ Gram and scalar products. The block system (R in the H₃ basis, y over K zero pairs) is solved by a Schur complement
(s15_joint.py). Its value is a lower bound for D² for every finite K; by Theorem A it tends to D² as K → ∞.

## 4. Carleson / Shapiro–Shields reading
Normalize: G̃ = Δ^{−1/2}HΔ^{−1/2}, Δ = diag Φ_N′(γ), c̃ = Δ^{−1/2}c, so |c̃|² = Σ_ρ |c_ρ|²/Φ_N′(γ) ≈ Σ_ρ 1/(|ρ|²(log N + log(γ/2π))) and
  |c̃|²/λ_max(G̃) ≤ cᴴH⁻¹c ≤ |c̃|²/λ_min(G̃).
Shapiro–Shields/Carleson: λ_min > 0 uniformly requires separation of the points in the space's metric, here σ = Φ_N′(γ)·Δγ/2π (gaps in Clark units of the phase Φ_N);
λ_max bounded requires a Carleson box condition (bounded counts in windows 2π/Φ_N′). Uniform separation (inf σ > 0, which GUE predicts false) plus the
box condition would make the evaluators a Riesz sequence and give D(λ)² ≤ F(λ) + |c̃|²/λ_min ≍ 1/log(1/λ); without it the upper side is void.

## 5. N8 side: Carleson's constant of the shifted zeros
Half-plane Re s > ½, b_z(s) = (s − z)/(s + z̄ − 1). For ε > B − ½ (B the zero-free abscissa; under RH any ε > 0) every z = ρ + ε has Re z > ½ and the reflections
1 − z̄ = ρ′ − ε (ρ′ = 1 − ρ̄ a zero) are the zeros of ξ(s + ε). Hence ∏_ρ b_{ρ+ε} = ξ(s−ε)/ξ(s+ε) up to a unimodular constant (both are unimodular on the
line, the ratio has no zeros or poles in the half-plane, and log|∏ b| ~ −ε log(σ/2π) matches ξ(σ−ε)/ξ(σ+ε) ~ (σ/2π)^{−ε}, so no singular factor), and
  δ_n(ε) := ∏_{m≠n} |b_{z_m}(z_n)| = (2β_n + 2ε − 1) |ξ′(ρ_n)| / |ξ(ρ_n + 2ε)|   (simple ρ_n).
For fixed ε, δ_n ≍ |ζ′(ρ_n)|/|ζ(ρ_n+2ε)| · (γ_n/2π)^{−ε} → 0: the shifted zeros are never an H^∞-interpolating sequence (Carleson). At ε = a/log(T/2π)
δ_n is a functional of the unfolded configuration (∏ (1 + a²/(π²x²))^{−1/2}, x unfolded distances): a local gap statistic.
