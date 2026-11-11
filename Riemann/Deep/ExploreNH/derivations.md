# explore24 derivations (non-Hermitian physics → RH tree)

Notation: H_t(z) = ∫_0^∞ e^{tu²}Φ(u) cos(zu) du, H_0(z) = ξ(1/2 + iz/2)/8, ∂_t H = −∂_z² H, zeros of H_0 at z = 2γ.

## D1. Exceptional-point (EP) normal form of a de Bruijn–Newman collision
At a double real zero (x_c, t_c): H = H_z = 0, H_zz = 2a ≠ 0, so H_t = −H_zz = −2a and
H ≈ a[(z − x_c)² − 2(t − t_c)], z_± = x_c ± √(2(t − t_c)).
Gap² = 8(t − t_c): a square-root (Puiseux) EP whose coefficient 8 is universal (independent of a, of the height and of the other zeros).
With the other zeros as a bath (E(x) = 2Σ_k 1/(x − z_k), S = Σ_k (x̄ − z_k)^{−2}):
dδ/dt = 4/δ + E(z_1) − E(z_2) ≈ 4/δ − 2Sδ, so d(δ²)/dt = 8 − 4Sδ² and
- Sδ₀² < 2: real EP at t_c = −(1/4S) log(1/(1 − Sδ₀²/2)) = −δ₀²/8 − Sδ₀⁴/32 − …;
- Sδ₀² > 2: complex EPs t_c = (1/4S)[log(Sδ₀²/2 − 1) ± iπ(2k + 1)] (in this two-body model).
Csordas–Smith–Varga (Constr. Approx. 10 (1994)): λ = ((1 − 5Δ²g/4)^{4/5} − 1)/(8g) = −Δ²/8 − Δ⁴g/64 − …, g ≈ 2S: the same expansion to second order. So the EP time of a Lehmer pair is the CSV bound.

## D2. Re-entry EP of an off-line pair (forward flow)
For z = x + iy with partner z̄: ẏ = −1/y − 2yS′, S′ = Σ_k |z − z_k|^{−2}, so d(y²)/dt = −2 − 4S′y² and
t_e = (1/4S) log(1 + 2Sy₀²) ≤ y₀²/2 (de Bruijn). This is Planat 2609.37164's inequality and explore15's screening factor 2Sδ²/log(1 + 2Sδ²).
An off-line pair at t = 0 is a broken-PT state whose EP lies at t_e > 0, hence Λ ≥ t_e. Nothing forces t_e ≤ 0.

## D3. Krein / Pontryagin index of Weil's form
Weil's form W(f) = Σ_ρ f̂(γ_ρ) · conj f̂(γ̄_ρ) (γ = (ρ − 1/2)/i; W(f) ≥ 0 for all f ⟺ RH).
- Real γ contributes |f̂(γ)|² ≥ 0.
- A conjugate pair γ, γ̄ contributes 2 Re(u v̄), u = f̂(γ), v = f̂(γ̄): a 2×2 Hermitian block of signature (1, 1).
So on every finite-dimensional test space the number of negative squares is at most the number of off-line pairs it sees. It equals that number when the evaluations at those points are independent (Paley–Wiener interpolation).
The test-function space with W is a Pontryagin space Π_κ, κ = number of off-line pairs. "Unbroken PT" means a positive metric, κ = 0, which is Weil positivity (N1).
Same count for −ξ′/ξ: it is a generalized Nevanlinna function of class N_κ (Krein–Langer 1977), κ = number of zeros of H_0 in the upper half-plane. κ = 0 is Lagarias' Herglotz/positivity criterion.
BBM's metric η = Δ†Δ is positive for free and induces the dilation generator with spectrum ℝ (2607.19067, R308). Its positivity carries no RH content.

## D4. Graphical Krein signature and the collision rule
Write a real-on-ℝ function F = 2 Re E with E = |E|e^{iφ}. Its zeros are where cos φ = 0, and the signature of a zero is sign φ′.
- Two adjacent zeros can merge only if φ′ vanishes between them (or |E| → 0). Hence collisions need opposite signatures (Krein–MacKay rule).
- If E is Hermite–Biehler (|E(z)| > |E^#(z)| on Im z > 0), then arg E is strictly decreasing on ℝ (φ′ < 0 everywhere; de Branges' phase −arg E increases). Every zero then has the same signature and persists (strong stability).
- For E = H + iH′ every simple zero has φ′ = −1, the same for all zeros. A collision of H is a real zero of E, where the metric degenerates.
- E is HB exactly when H ∈ LP: Im(H′/H) < 0 on the upper half-plane gives |H + iH′| > |H − iH′|. The heat flow carries E_t = e^{−t∂²}E_0 = H_t + iH_t′ and preserves HB only forward in t.
- So "a signature that forbids collisions for t > t_0" exists ⟺ H_{t_0} ∈ LP ⟺ Λ ≤ t_0: a restatement.
- The scalar picture has no protection: consecutive real zeros of a real function always carry opposite graphical signatures (alternating slopes), so every adjacent pair may collide.

## D5. EPs are invisible to spectral invariants
- The Taylor coefficients a_k(t) = ((−1)^k/(2k)!) ∫_0^∞ u^{2k} e^{tu²} Φ(u) du are entire in t.
- The power sums s_k(t) = Σ_j z_j(t)^{−2k} are real-analytic for all real t (H_t(0) > 0).
- Log-potentials Σ_j log|1 − x/z_j(t)| are real-analytic away from x ∈ zero set.
- At an EP the two Puiseux branches are swapped by monodromy, so every symmetric function stays single-valued and bounded, hence analytic.
- Hence Λ is not a singularity of any holomorphic spectral invariant (no Lee–Yang/Fisher-type non-analyticity). It is a signature threshold only: the Hermite–Grommer Hankel forms of (s_k(t)) lose positivity, with negative index = N_+(t), non-increasing in t (Pólya; de Bruijn 1950). This is why analyticity arguments (R48) cannot see Λ.

## D6. Scattering / CPA / lasing reading (N7, N8)
- Modular surface (Pavlov–Faddeev, Lax–Phillips): S(s) = Λ(2s − 1)/Λ(2s), Λ(s) = π^{−s/2}Γ(s/2)ζ(s), unimodular on Re s = 1/2, poles (resonances) at s = ρ/2.
- A pole or zero on the unitary axis (lasing threshold or coherent perfect absorber) would need ζ(1 + 2ir) = 0. This is excluded by the prime number theorem, which Sarnak (2004) proves this way via Eisenstein series.
- RH ⟺ every resonance has the same decay rate 1/4 ("a single common lasing threshold"). That is a restatement. Davenport–Heilbronn has no automorphic scattering realization, and resonance positions are unconstrained by the formalism.

## D7. Skin effect / generalized Brillouin zone (N8, N13)
- An OBC Toeplitz spectrum is real when an imaginary gauge transformation (Hatano–Nelson, z → rz) makes the symbol real on a circle. Reason: T_n(a) is then similar to a Hermitian matrix (2303.15186 gives the precise criterion).
- ξ is real on Re s = 1/2 by the functional equation, but so is Davenport–Heilbronn. A real symbol forces real zeros only through a determinant/operator realization (Hilbert–Pólya).
- The equimodular (Schmidt–Spitzer, Beraha–Kahane–Weiss) principle "zeros where two dominant terms balance" is Polymath15's large-height mechanism for t > 0: the dominant AFE half is locally HB. It fails at t = 0, where no two terms dominate.
