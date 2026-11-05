# explore10 derivations (2026-11-05)

Notation. Zero ρ = 1/2 + α + i x enters every explicit formula as H(x − iα), with H(z) = ∫ G(u) e^{izu} du,
G even, real, supported in [−R, R]. Its mirror 1 − ρ̄ enters as H(x + iα) = conj H(x − iα), so a pair weighs
2 Re H(x − iα) = 2 ∫ G(u) cosh(αu) cos(xu) du.

## D1. Trace method on surfaces (what is used)

Selberg: Σ_j ĥ(r_j) = (g−1)∫ ĥ(r) r tanh(πr) dr + Σ_γ ℓ(γ₀) h(ℓ(γ)) / (2 sinh(ℓ(γ)/2)), r_j ∈ ℝ ∪ i[−1/2, 1/2]
(self-adjointness). With ĥ = φ̂² (entire extension; φ even, real, φ ≥ 0, supp ⊂ [−R/2, R/2]) every spectral term
is ≥ 0, ĥ(iβ) = (∫φ cosh βu)² increases in β, an exceptional r = iα weighs ≈ e^{αR}, so P(λ₁ ≤ 1/4 − α²) ≤ E[Σ_{j≥1} ĥ(r_j)] / ĥ(iα) (Markov). The right side
needs: (a) the identity term ≈ g/R² (Plancherel density r tanh πr vanishes at r = 0), so R ≥ (1/α) log g;
(b) E[geodesic side] at lengths up to A log g: Mirzakhani's integration formula for simple curves, the
local-type densities V_g^T(ℓ) with 1/g-expansions whose coefficients are Friedman–Ramanujan
(p(ℓ)e^ℓ + O(ℓ^c e^{ℓ/2})), the operator D = 1/4 − d²/dℓ² applied to the test function (D annihilates sinh(ℓ/2),
so it kills λ₀ = 0; since (1/4 − d²)(r(ℓ)e^{ℓ/2}) = −(r'' + r')e^{ℓ/2}, D^m kills, after the weight ℓ/(2 sinh(ℓ/2)),
the FR main terms p(ℓ)e^ℓ for m > deg p + 1, by integration by parts), and the
removal of tangles (Moebius inversion, 2401.01601). Petri's Bourbaki survey 2607.06331 §(trace method):
"using the fact that h_R is non-negative on R ∪ iR, all the contributions to the spectral side ... are
non-negative"; 3/16 − ε comes from R = 4 log g with simple geodesics only (Wu–Xue, Lipnowski–Wright).
Randomness enters only (b). Positivity is deterministic.

## D2. Family first moment for primitive characters mod a prime q

For primitive χ of parity a (IK Thm 5.12 form):
Σ_ρ H(γ_ρ) = (1/2π)∫ H(r)[log(q/π) + Re ψ((1/2 + a + ir)/2)] dr − Σ_n Λ(n) n^{−1/2} [χ(n) + χ̄(n)] G(log n).
Average over the q − 2 primitive χ mod q: Σ*_χ χ(n) = (q−1)·1_{n≡1 (q)} − 1_{q∤n}. Hence
  Z(q, R) := E_χ Σ_ρ H(γ_ρ) = arch(q, R) − (2/(q−2)) · D,   D = (q−1) S₁ − S₀,
  S₁ = Σ_{n≡1 (q)} Λ(n) n^{−1/2} G(log n),  S₀ = Σ_{q∤n} Λ(n) n^{−1/2} G(log n).
Test function: G(u) = (1 − |u|/R)₊, H(r) = R (sin(rR/2)/(rR/2))² ≥ 0 on ℝ and on iℝ (the AM class),
(1/2π)∫H = 1, R = A log q.
* The missing principal character contributes −S₀ ≈ −(4/R)e^{R/2}: the analogue of λ₀ = 0. It is cancelled by
  the AP primes (q−1)S₁: the analogue of "simple geodesics cancel λ₀".
* For R < log(q+1) (A < 1; in fact up to the least prime POWER ≡ 1 mod q, e.g. 8 for q = 7): S₁ = 0 exactly
  (orthogonality is exact: the family is "tangle-free" below length log q).
* Unconditional upper bound (sign of the prime side): S₁ ≥ 0, so Z ≤ Z_up := arch + 2S₀/(q−2), and
  2S₀/(q−2) ≈ 8 q^{A/2−1}/(A log q): harmless up to A = 2, then grows like q^{A/2−1}.
* Given ZERO-WISE positivity (every single zero contributes ≥ 0, which holds when all zeros lie on the critical
  line or on the real segment, as the spectrum of a surface does), Markov gives
  for the box B_α = [1/2 + α, 1] × [−π/(3R), π/(3R)]: a pair in B_α weighs ≥ cos(π/3)·2H(iα) = H(iα), so
  P_χ(zero in B_α) ≤ Z / H(iα),  H(iα) = R (sinh(αR/2)/(αR/2))².
  At A = 2 (provable range): P ≲ (log q)·α²R·e^{−αR} = 2α²(log q)² q^{−2α}: → 0 for every α > 0.
  So for this "almost all near the real axis" statement the averaged "geometric" side is not the bottleneck;
  zero-wise positivity is. (Whether (q−1)S₁ cancels S₀ beyond A = 2 is not proved; it is not needed here.)
* Zero-wise positivity for one fixed kernel is weaker than GRH, but for the whole AM class (all R) it fails at
  every non-real off-line zero: with w = α + ix, ∫₀^R (1 − u/R)e^{wu}du = (e^{wR} − wR − 1)/(w²R), so the triangle
  kernel's pair weight 2Re H(x − iα) is 2e^{αR} Re(e^{ixR}/w²)/R + O(1 + 1/|w|): for x ≠ 0 its sign rotates with R at amplitude
  → ∞ (T2 numerics), so it forces every
  non-real zero onto the critical line: GRH up to real zeros. The random model never supplies it (on surfaces it
  is free). This is the whole transfer gap.

## D3. Harnack cap for kernels that are positive on the whole strip

To drop unknown zeros pointwise, without GRH, a kernel needs Re H ≥ 0 wherever zeros may lie: on the closed strip
S = {|Im z| ≤ 1/2}.
Characterisation (G even, real, continuous, compactly supported). Re H(x + iy) = FT[G(u) cosh(yu)](x). On the
boundary line y = 1/2 this is ≥ 0 for all x iff k := G·cosh(u/2) is positive definite (Bochner; continuity matters:
changing G at one point changes nothing in H). So every strip-positive kernel is G = k/cosh(u/2) with k
positive definite (the Poitou–Odlyzko form), and conversely: for |y| < 1/2, cosh(yu)/cosh(u/2) has the positive
Fourier density 4π cos(πy) cosh(πx)/(cosh 2πx + cos 2πy); at |y| = 1/2 the ratio is 1 (transform 2πδ₀); so
Re H ≥ 0 on S (Schur product).
Resolution. On the line Re H(x) = (k̂ * sech(π·))(x) (FT of sech(u/2) is 2π sech(πx)): a positive mixture of
translates of sech(πx), whatever the support R. Measured half-width at half maximum 0.730, 0.540, 0.450, 0.430
for R = 5, 10, 40, 160 against 0.419 for sech(πx); H(0) → 2π k(0) (t2c_out.txt).
Cap. Re H is positive harmonic in S and S is invariant under real shifts. The map z ↦ tanh(πz/2) sends S onto the
unit disc and x − iα to a point at Möbius distance tan(πα/2) from the image of x. Harnack:
  Re H(x − iα) ≤ tan(π/4 + πα/2) · Re H(x)   for all x ∈ ℝ, 0 ≤ α < 1/2,
independent of R and of the support; the constant is sharp for positive harmonic functions (Poisson kernel at the
boundary point x − i/2). An off-line zero is worth at most tan(π/4 + πα/2) on-line zeros at the same height
(α = 1/4: 2.414; α = 0.3085: 3.224).
Consequence for a first moment over characters mod q (k regular, e.g. k = φ*φ̃ with φ ∈ C¹_c, so the explicit
formula applies). Unconditionally Z = E_χ Σ_ρ Re H ≥ k(0)(log(q/π) + ψ(1/4)) − O(k(0)R/q): Re ψ(σ + it) ≥ ψ(σ) and
H ≥ 0 on the line bound the archimedean term below uniformly in k; the prime terms are ≤ Σ Λ(n)·2k(0)/n per class
(|k| ≤ k(0)), and Brun–Titchmarsh makes their family average O(k(0)R/q). Meanwhile a bad pair weighs ≤ 2 tan(π/4 + πα/2)·sup_x Re H(x)
≤ 4π k(0) tan(π/4 + πα/2) (Re H = k̂ * sech(π·) ≤ ∫k̂ = 2π k(0)). So the Markov bound is ≥ (log q − 5.4 − O(R/q))/(4π tan(π/4 + πα/2)) for every support R ≤ q:
it grows with q (main term 0.46 at q = 10⁶, 2.28 at q = 10³⁰ for α = 1/4) and never isolates a member.
Near Re s = 1 (α = 1/2 − δ) the cap is ≈ 2/(πδ): the classical width-1/log q, pole-anchored mechanism (R30).

## D4. Strong convergence needs a non-amenable limit

Random covers: std∘ρ_n for uniform ρ_n ∈ Hom(Γ, S_n) converges strongly to the regular representation of Γ
(Bordenave–Collins; CGTV 2405.16026; MPvH 2504.08988); the limit norm is the Kesten/McKean value (λ₀(ℍ) = 1/4),
which is below the trivial value because Γ is non-amenable. Arithmetic analogue: the direct sum of all characters
mod q (the regular representation of (ℤ/q)^*), restricted to the primes p ≤ y (y fixed), converges strongly as
q → ∞ to the regular representation of ℤ^{π(y)} (the vectors (χ(p))_{p≤y} equidistribute on T^{π(y)}; a single
sampled character does not converge). The limit norm of Σ_p (z_p + z̄_p) is 2π(y), the trivial value, because
ℤ^{π(y)} is amenable (Brooks: amenable covers have no uniform gap). Strong convergence itself holds; a gap needs a
non-amenable limit. For growing y the observed decay of the family maximum follows, heuristically, the union
bound over q − 2 characters at Steinhaus size √(π(y) log q); the moments are Steinhaus-exact (up to the principal
character, off-diagonal average −1/(q−2)) only for products below q: the large-sieve range. No RH content.
