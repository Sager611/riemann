# DERIVATION (aclgh2): models and the onset identity

Notation of [../ACLGFals/DERIVATION.md](../ACLGFals/DERIVATION.md): x = log u, M g(iτ) = ∫ g e^{iτx} dx; Γ_t = (1 − P_{1/t}) U P_t ≅ Hankel
operator with symbol φ_t = t^{2iτ} u, u = u_∞ ∏_{p∈S} u_p; Γ_t e_n = κ_n g_n, f_n = U^* g_n, a_n = (e_n + f_n)/√(2 + 2κ_n);
m_n(R) = (1/2π) ∫ χ_R |M a_n|²; T(n₀, R) = Σ_{n>n₀} κ_n m_n(R) (ACLG eq. H2ii, TeX l. 1257); κ_0 ≥ κ_1 ≥ …

## 1. Cut model (code/semi_pair2.py; aclgfals's Nyström model)
Kernel 2t C_S(t²uw) on L²(0,1), C_S(z) = Σ_{n≤nmax} γ(n) cos(2πnz) (n = ∏ p^{k_p}, k_p ≥ −1). Changes to semi_pair.py: C_S tabulated
(32 points per top wavelength, 8-point Lagrange; |Δκ| ≤ 7·10⁻¹¹ against the direct sum), scale t, two-ended eigen-list cut to its
complete range, one vector matrix. Why it fails for ≥ 3 primes: the cut n ≤ nmax removes depths |x| > log nmax, and the near-self-dual
vectors live there: G(u) ~ Σ_{n S-smooth} h(nu) ~ (log 1/u)^r (r = |S| − 1), centroid depth ≈ 2r (band-limited model: 4.9, 6.1, 7.3,
8.3 for r = 2…5; 21.7 for r = 25). Operator error Σ_{n>nmax} |γ(n)| n^{−1/2} ~ nmax^{−1/2} (log nmax)^{r−1}.

## 2. Band-limited model (code/laurent.py, laurent_count.py, centroid.py)
- Grid x_j = (j + ½)Δx, Δx = π/τ_max; band-limited functions ↔ samples (×√Δx). U ↔ Laurent matrix (c_{i−j}),
  c_m = (1/2π) ∫_{−π}^{π} φ̃(θ/Δx) e^{−imθ} dθ (FFT, real); unitary on ℓ²(ℤ) since |φ̃| = 1. P_t ↔ j ≤ −1, so Γ_t ↔ H_{il} = c_{i+l+1},
  i, l < N = depth/Δx; H real symmetric, H v = μ v: κ = |μ|, e = v (x = −(l + ½)Δx), g = sign(μ) v (x = (i + ½)Δx).
  (P, U^*(1 − P)U) is an exact pair of projections; u_p enters exactly (no Laurent cut) on |τ| ≤ τ_1.
- Phase closure: φ̃ = e^{iψ̃}, ψ = −2ϑ(τ) + 2τ log t + Σ_p 2 atan2(a sin τℓ, 1 − a cos τℓ), a = p^{−1/2}, ℓ = log p (= arg u_p(iτ));
  ψ̃ = ψ on |τ| ≤ τ_1 = b₁τ_max (b₁ = 0.75), smooth blend on [τ_1, (1 + b₁)τ_max/2] into the line through (τ_1, ψ(τ_1)) and
  (τ_max, πK), K the nearest integer (slope ≈ ψ'(τ_1) < 0); odd and smooth on the circle. Negative slope = negative delay, so the
  closure adds no coupled directions (S = {∞, 2}: N(κ > ½) = 4 at τ_max = 400, 800, 1180).
- Mellin data: M e(iτ) = √Δx Σ_l e_l e^{−iτ(l+½)Δx}, M g(iτ) = √Δx Σ_i g_i e^{iτ(i+½)Δx}, M f = conj(φ) M g (|τ| ≤ 2R < τ_1).
- S = {∞}: κ_n − (1 − δ_n(1)) = O(Δx²) (largest relative error over κ_0…κ_4: 3.7·10⁻² at τ_max = 100, 2.3·10⁻³ at 400).

## 3. Onset = counting function at an angle threshold
T(n₀, R) = W(κ_{n₀}, R), W(κ, R) := Σ_{κ_n < κ} κ_n m_n(R) (ties aside). Put
κ♮(S, R, τ) := min{κ_n : Σ_{κ_{n'} ≤ κ_n} κ_{n'} m_{n'}(R) > τ}. Then the onset n*(S, R, τ) = min{n₀ : T(n₀, R) ≤ τ} is the index of
κ♮: n* = N_S(κ♮), N_S(κ) := #{n : κ_n > κ}. W sees only directions with Mellin mass at |τ| ≤ 2R (resolved and converged in §2);
N_S counts all directions, mostly at high Mellin frequency. As κ_n ≤ 1, n* = N_S(κ♮) ≥ Σ_{κ_n > κ♮} κ_n².

## 4. Semiclassical count of coupled directions
Delay D_t(τ) = ψ'(τ) = 2 log t − 2ϑ'(τ) + Σ_p δ_p(τ log p), δ_p(θ) = 2 log p Re(a e^{iθ}/(1 − a e^{iθ})) = 2 log p Σ_{k≥1} p^{−k/2} cos kθ,
range [−2 log p/(√p + 1), 2 log p/(√p − 1)]. U shifts a wave packet at frequency τ by D_t(τ) in x; Γ_t couples {−D_t(τ) < x < 0};
Weyl–Landau count N_sc = (1/2π) ∫ (D_t)_+ dτ. As −2ϑ'(τ) = −log(|τ|/2π) + O(τ⁻²), coupling at large τ needs Σ_p δ_p > log(|τ|/2πt²):
joint near-resonances τ log p ≡ 0 (mod 2π), up to the horizon 2πt² e^{M_S}, M_S = Σ_p 2 log p/(√p − 1) (~ 4λ for p ≤ λ²).
Identity: (τ log p mod 2π)_p equidistributes (Kronecker–Weyl), so with X = Σ_p δ_p(θ_p), θ_p iid uniform,
(1/2π) ∫ (log(2πt²/|τ|) + X)_+ dτ = (1/π) E ∫_0^{2πt²e^X} log(2πt² e^X/τ) dτ = 2t² E e^X = 2t² F_S,  F_S := ∏_{p∈S} E_θ e^{δ_p(θ)}.
E_θ e^{δ_p} = 1 + (log p)²/(p − 1) + O((log p)³ p^{−3/2}) (1.0085278 by quadrature, 1.0085019 by the second-order term, p = 9973),
hence log F_S = Σ_{p≤x} (log p)²/p + O(1) = ½ log² x + O(log x) (partial summation from Mertens); numerically 94.76 against
½ log² x = 95.43 at x = 10⁶. Direct integration of (1/2π)∫ D_+ (code/semicl.py): 4.89, 11.24, 24.65 for S = {∞,2}, {∞,2,3},
{∞,2,3,5}, against 2F_S = 4.67, 11.19, 24.95; S = {∞}: 2.25 (exact count of κ > ½: 3). Exact pocket scan (code/pockets.py):
for S = {∞,2,…,11}, Σ over the 970 pockets on 0 ≤ τ ≤ 5·10⁶ of A/2π = 48.66, twice that 97.3 against 2F_S = 96.9.
Per pocket (shifted band, code/window_count.py; e^{iτ_c x} commutes with P): one direction per sign of τ with κ² ≈ A/2π for thin pockets
(A < 2π), so Σκ² over coupled directions ≈ N_sc; fat pockets (only at τ ≈ 0) give ≈ A/2π directions with κ ≈ 1.
