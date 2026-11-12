# ExplorePotential (explore35): derivations

Zeros ρ = β + iγ of ζ (or of the Davenport–Heilbronn function f), counted with multiplicity; δ = β − 1/2; 0 < γ ≤ T.
S(x, T) = Σ_{0<γ≤T} x^{iγ} (ordinates only), L(x, T) = Σ_{0<γ≤T} x^ρ, x > 1; e(u) = e^{2πiu}; x_k = e^{2πkα}, so Σ_γ e(kαγ) = S(x_k, T).

## 1. Tool: Erdős–Turán
For points u_1, …, u_N in ℝ/ℤ and every K ≥ 1 (Montgomery, Ten Lectures, ch. 1, with Vaaler's constants):
  D := sup_I |#{n: u_n ∈ I} − N|I|| ≤ N/(K + 1) + 3 Σ_{k=1}^{K} |Σ_n e(k u_n)|/k.
Applied to u_γ = αγ the exponential sums are S(x_k, T), the inputs the brief names (Landau–Gonek).

## 2. Lemma 1 (pairing identity; exact, unconditional)
  S(x, T) = x^{−1/2} L(x, T) − 2 E_off(x, T),   E_off(x, T) = Σ_{δ>0, 0<γ≤T} x^{iγ} (cosh(δ log x) − 1).
Proof. ρ ↦ 1 − ρ̄ maps zeros to zeros (functional equation and real coefficients; for f too) and keeps γ. A zero with
δ ≠ 0 pairs with 1/2 − δ + iγ: the pair gives √x x^{iγ}(x^δ + x^{−δ}) in L and 2x^{iγ} in S; zeros with δ = 0 give √x x^{iγ}
and x^{iγ}. Subtract. ∎
The same split with on-line ordinates only: S_on(x, T) = x^{−1/2} L(x, T) − Σ_{δ>0} 2 x^{iγ} cosh(δ log x)  (on-line share).

## 3. Input: Landau–Gonek (Gonek 1993, unconditional)
  L(x, T) = −(T/2π) Λ(x) + O(x log(2xT) log log 3x) + O(log x · min(T, x/⟨x⟩)) + O(log(2T) · min(T, 1/log x)),
⟨x⟩ = distance from x to the nearest prime power other than x, Λ(x) = 0 off the prime powers. For f the same contour
argument (right edge beyond the abscissa of absolute convergence of f′/f) gives L_f(x, T) = −(T/2π) c(x) + O_x(log T) with
−f′/f = Σ c(n) n^{−s}, i.e. Σ_{d|n} c(d) a(n/d) = a(n) log n, a = (1, κ, −κ, −1, 0) mod 5, κ = (√(10 − 2√5) − 2)/(√5 − 1).

## 4. Lemma 2 (Selberg's density bounds the off-line term)
For 2 ≤ x ≤ T^{1/8}:  E_off(x, T) ≪ T (log x)²/log T.
Proof. Selberg (1946): N(1/2 + u, T) ≪ T^{1 − u/4} log T for 0 ≤ u ≤ 1/2. By parts,
Σ_{δ>0}(cosh(δ log x) − 1) = ∫_0^{1/2} log x · sinh(u log x) N(1/2 + u, T) du ≪ T log T log x ∫_0^∞ sinh(u log x) e^{−u(log T)/4} du
= T log T (log x)²/((log T/4)² − (log x)²) ≪ T (log x)²/log T. ∎

Corollary 3 (first-order law, unconditional). S(x, T) = −(T/2π) Λ(x)/√x + O(T (log x)²/log T) + O(√x log(2xT) log log 3x) + O(x^{−1/2} log 2T/log x).
With N(T) ~ (T/2π) log(T/2π): the Fourier coefficients of the counting measure of {αγ} are
  N(T)^{−1} S(x_k, T) = −Λ(x_k)/(√x_k log(T/2π)) + O_{α,k}(1/log² T).
Under RH, E_off = 0 and only Gonek's error remains. Corollary 4 (RH): with K = ⌊(log T)/(2πα)⌋ in §1,
D/N(T) ≤ (2πα + 3 Σ_{k≤K} Λ(x_k)/(k √x_k) + o(1))/log T: Hlawka's 1/log T (1975), with this constant
(Gonek's error enters S as O(√x log(2xT) log log 3x), at most O(√T log T log log T) for x ≤ T).

## 5. What the statistic sees (the barrier)
(i) Moved pair (γ fixed, δ: 0 → δ₀): S unchanged; L changes by 2√x x^{iγ}(cosh(δ₀ log x) − 1) = O_x(1); invisible against T.
(ii) Near-line families (δ log x ≪ 1): cosh − 1 = O(δ² log² x) per zero (it grows like x^δ/2 once δ log x ≳ 1); any off-line set
    allowed by Lemma 2 shifts S by o(T).
(iii) A family at fixed distance with ≍ T zeros and coherent phases shifts S (and the on-line share) by ≍ T. Selberg excludes
    it for ζ (T^{1−δ/4} log T = o(T)). For f, N(σ, T) ≍ T at fixed σ > 1/2 (almost periodicity; f has zeros in σ > 1,
    Davenport–Heilbronn 1936), so Lemma 2 fails and such a family can exist.
So, given the primes, the first-order ordinate law constrains off-line zeros only through E_off. Every hypothetical ζ zero set that
obeys Selberg's density obeys it, so the law certifies nothing beyond zero density (N7). This is a statement about this law, not about
every ordinate statistic: L itself depends on β. S_on is not fixed by Landau plus Selberg, and whether another method could fix it
unconditionally is open.
(iv) Locality. ET at level K resolves {αγ} at 1/K, a statistic aggregated over all heights. Counts in [t, t + h] (N8: gaps,
multiplicities) need α ≲ 1/h, where x_k → 1 and S is N(T) and S(T) (Backlund), or x ≥ T, where Gonek's √x log(xT) exceeds
the main term T Λ(x)/√x as soon as x ≳ T/log T. So ET + Landau–Gonek never reach N8's local quantities.

## 6. DH on-line share (adversarial measurement)
By §2 for f: S_on + Σ_{δ>0} 2x^{iγ}cosh(δ log x) = −(T/2π) c(x)/√x + O_x(log T). The on-line share is S_on divided by the right side.
For ζ under RH (verified height) the share is 1 up to O(log T). For f (dh_ft.py, dh_offsearch.py): 2,406 on-line zeros and 109
off-line right members below 2500 (218 zeros, count deficit 217.4); with them the sum rule closes (residual ≤ 2.1 for x ≤ 40 at
T = 420, …, 2500). The on-line share at T = 2500: −1.86 (x = 2), −1.20 (3), 0.59 (4), 0.38 (6), 0.61 (9), 0.87 (11),
falling with T (x = 6: 0.66, 0.57, 0.52, 0.39, 0.38 at T = 420, 1000, 1500, 2000, 2500). The off-line zeros are phase-locked to
the leading Dirichlet terms (f = 1 + κ2^{−s} − κ3^{−s} − 4^{−s} + …): mean cos(γ log 2) = −0.73, mean cos(γ log 3) = +0.70,
i.e. 2^{−iγ} ≈ −1 and 3^{−iγ} ≈ +1, where those terms oppose the 1.

## 7. Phase-locked ζ family (model)
Take F = {γ < 2·10⁴: cos(γ log 2) < −0.9} (3,470 ordinates, 15.4%) as split double zeros 1/2 ± δ + iγ (zeta_family.py). Then S is
unchanged, and x^{−1/2}L changes by 2E_off: at x = 2, |2E_off| = 0.16, 16.1, 146 for δ = 0.01, 0.1, 0.3, against Landau's
error scale max |S − P| = 6.5. The on-line ordinates lose 2Σ_F x^{iγ} = −4.30 × the Landau term at every δ. So the primes see the family only through
(δ log x)² (visible here for δ ≳ 0.064, i.e. δ log T ≳ 0.6), while the on-line share sees it at first order.
