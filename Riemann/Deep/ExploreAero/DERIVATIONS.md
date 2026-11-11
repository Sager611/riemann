# ExploreAero derivations (explore26)

## D1 Airfoil operator = function of the dilation generator (Koppelman–Pincus 1959, Mellin form)
T f(x) = (1/π) PV∫_{-1}^{1} f(t) dt/(t − x) on L²(−1,1). Put u = (1+x)/(1−x) ∈ (0,∞), F(v) = f(t(v))/(1+v). Then
t − x = 2(v − u)/((1+v)(1+u)), dt = 2dv/(1+v)², so (Tf)(x(u))/(1+u) = (1/π) PV∫_0^∞ F(v) dv/(v − u) = (H₊F)(u) and
‖f‖² = 2‖F‖²: T ≅ H₊ unitarily. From PV∫_0^∞ v^{s−1}/(v − u) dv = −π u^{s−1} cot πs (0 < Re s < 1), on Re s = 1/2:
**T = i tanh(πH)**, H = −i(u d/du + 1/2) (Berry–Keating xp, essentially self-adjoint on L²(0,∞)); generalized eigenfunctions
g_τ(x) = ((1+x)/(1−x))^{iτ}/√(1−x²), T g_τ = i tanh(πτ) g_τ (checked to 12 digits, code/t1_airfoil.py). In x, u d/du = ((1−x²)/2) d/dx.
- Kutta mode: g₀ = 1/√(1−x²) (∉ L²) is the generalized eigenvector at τ = 0, i.e. s = 1/2. Inversion: with G = γ/(1+v), H₊G = Φ gives
  G(v) = ∫_{Re s = c} −φ(s) tan(πs) v^{s−1} ds/2πi, 0 < c < 1; tan πs has its pole at s = 1/2. c < 1/2: no residue as v → ∞, γ bounded at
  the trailing edge x = 1 (Kutta); c > 1/2: γ ~ √((1+x)/(1−x)) there (bounded at the leading edge instead). Söhngen–Tricomi's weights
  √((1∓x)/(1±x)) = u^{∓1/2} implement the two sides; the two inverses differ by the residue at s = 1/2 = circulation × g₀; L² sits on
  the pole (Re s = 1/2), hence no bounded inverse.
- L^p: Mellin line Re s = 1/p, symbol −cot πs vanishes only at s = 1/2 → T is onto with kernel span{g₀} for 1 < p < 2, injective with
  codim-1 range for p > 2, not Fredholm at p = 2 (Widom 1960, Okada–Elliott 1991). This is Beurling's L^p/Nyman mechanism with one
  "zero" at s = 1/2 (N5 already has the arithmetic version: ζ ≠ 0 on Re s > 1/p ⟺ 1 ∈ L^p-closure).
- Skewness: T* = −T, so ⟨f, Tf⟩ = 0 for real f: no finite-Hilbert-transform operator contributes to a real quadratic form.
- Even partner: in y = log u the kernel of T is PV 1/(2π sinh((w−y)/2)). Its even (finite-part) twin 1/(2 sinh(|y|/2)) has symbol
  Re ψ(1/2+iτ) − ψ(1/2): ∫_0^∞(1 − cos τt)/(2 sinh(t/2)) dt = Re ψ(1/2+iτ) − ψ(1/2), ∫_0^∞ sin τt/(2 sinh(t/2)) dt = Im ψ(1/2+iτ) =
  (π/2) tanh πτ. So the causal kernel 1_{t>0}/(2 sinh(t/2)) has symbol −ψ(1/2 − iτ) + const: its odd part is the airfoil operator, its
  even part is the archimedean Weil term of the complex place Γ_C(s). For ζ (real place) the kernel in Q is
  e^{x/2}/sinh x = 1/(2 sinh(x/2)) + 1/(2 cosh(x/2)) (from ψ(z) = −γ + ∫_0^∞(e^{−t} − e^{−zt})/(1−e^{−t})dt at z = 1/4 + ir/2, t = 2x).
  The airfoil/archimedean pair is a Kramers–Kronig (Hilbert) pair — the R13.6 relation at the archimedean place.
- Multhopp/Legendre–Galerkin truncation: ⟨p_m, T p_n⟩ = −(2/π)√((2m+1)(2n+1))/2 · (1 − (−1)^{m+n})/((m−n)(m+n+1)) (Neumann
  T P_n = −(2/π)Q_n), real antisymmetric. Resolution ~K^{−2} at the edges = a log-window in u: box-regularized xp, so N(τ) is linear
  in τ with slope ∝ log K; Riemann–von Mangoldt needs (τ/2π) log(τ/2πe).

## D2 Lifting line, Glauert, Munk; the Weil wing
- Prandtl: induced drag D = ⟨Γ, |D|Γ⟩ (|D| = (−Δ)^{1/2} restricted to the span, Γ extended by 0), lift ∝ ∫Γ. Glauert x = a cos θ,
  Γ = Σ A_k sin kθ: H(Γ') maps sin kθ to k sin kθ/sin θ (Glauert's integral (1/π)∫_0^π cos kφ/(cos φ − cos θ)dφ = sin kθ/sin θ), so
  ⟨sin kθ, |D| sin lθ⟩ = (π/2)k δ_kl, lift ∝ A₁: D ∝ Σ k A_k² ≥ A₁² — Munk 1921: elliptic loading is optimal. FT of sin kθ:
  π i^{k−1} k J_k(ar)/r; pole functional ⟨sin kθ, cosh(x/2)⟩ = 2πk I_k(a/2) (≈ a single mode for small a).
- Weil (even f on [−a, a], c = e^{2a}): Q(f) = 2ℓ(f)² + R(f), ℓ(f) = ⟨f, cosh(x/2)⟩, R = A − P (archimedean minus primes). Hence
  Q ≥ 0 ⟺ R ≥ 0 on ker ℓ and D_min := min{R(f) : ℓ(f) = 1} ≥ −2 ("minimum drag at unit lift ≥ −2"). Odd f: Q = −2s(f)² + R,
  s(f) = ⟨f, sinh(x/2)⟩ ("negative lift", roll), so Q ≥ 0 ⟺ R > 0 and Munk's maximal efficiency 2 sᵀR^{−1}s ≤ 1, optimiser R^{−1}s.
  This is the rank-one (Schur-complement) treatment of the pole term (Yoshida 1992; Connes–Consani 2006.13771; CCM 2511.22755).
- Munk's proof needs (i) a drag operator diagonal in the basis that contains the lift functional as one mode. Weil's archimedean symbol
  θ′(r) = (1/2)(Re ψ(1/4 + ir/2) − log π) ~ (1/2) log(r/2π) is of order 0 (log), negative for |r| < 6.29, and not diagonal in any
  Chebyshev basis (only the log kernel, the inverse direction, is: ∫ −log|x−t| T_n(t)dt/√(1−t²) = (π/n)T_n(x)); plus the prime
  shifts. Order 0 also means the form domain contains jumps, so no √ tip law: the optimum's tip load vanishes only logarithmically.
- Arch matrix in x-space (code/weilleg.py): (−γ − log π)g(0) + ∫_0^{2a}[g(0)e^{−2x} − e^{−x/2}g(x)] 2dx/(1−e^{−2x}) − g(0) log(1 − e^{−4a}),
  g = f ⋆ f̃; primes −2Σ_{k<c} Λ(k)k^{−1/2} g(log k).
- Moved pair: replacing ±γ by the quartet ±γ ± iδ changes Q by −2 f̂(γ)² + 4 Re f̂(γ + iδ)² (even real f).

## D3 Flutter ↔ de Bruijn–Newman (no new numerics; explore24 has them)
- Coupled-mode flutter: two real frequencies meet and leave the axis (Krein collision). H_t: ż_k = 2Σ_{j≠k} 1/(z_k − z_j), so a close
  pair obeys d(gap²)/dt = 8 + O(gap² × external field): gap² = 8(t − t_c) is the Zimmerman–Weissenburger flutter margin (a discriminant,
  linear in the parameter); its extrapolated onset is the Csordas–Smith–Varga Lehmer-pair bound (explore24: −7.11203·10⁻⁴ ODE vs
  −7.11206·10⁻⁴ CSV, Lehmer's pair at 7005). Lower bounds on Λ only; Λ ≥ 0 (Rodgers–Tao).
- Energy method (Nissim's aerodynamic energy, Lyapunov): for an isolated conjugate pair d(Im z)²/dt = −2 → de Bruijn's strip theorem
  Λ ≤ y_max²/2 (with the Lean half-plane 0.2788, N2); V-g / p-k sweeps = argument principle (R404).
- Theodorsen C(k) = H₁⁽²⁾/(H₁⁽²⁾ + iH₀⁽²⁾) = K₁/(K₀ + K₁)(ik): 1/2 < Re C ≤ 1 ⟺ |H₀⁽²⁾| < |H₁⁽²⁾| (Nicholson), no zeros of
  K₀ + K₁ on the physical sheet: a passive (positive-real) transfer function; the ξ analogue is Lagarias' Re ξ′/ξ > 0 (explore22).
  Wagner/Küssner: branch cut on the negative real p-axis ⇒ 1 − φ completely monotone (R409, R419).
