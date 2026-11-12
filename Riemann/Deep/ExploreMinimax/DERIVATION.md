# ExploreMinimax — derivations (explore36)

## D1. Prime data = zero symbol on |y| < L
Weil's explicit formula, g ∈ C_c^∞(−L, L), ĝ(z) = ∫ g(y) e^{izy} dy, ordinates γ_ρ = (ρ − ½)/i (complex off the line):
Σ_ρ ĝ(γ_ρ) = ĝ(i/2) + ĝ(−i/2) − Σ_{n < e^L} Λ(n) n^{−1/2} (g(log n) + g(−log n)) + Arch(g).
So {Λ(n) : n < X = e^L} ⟺ the distribution G_Z(y) = Σ_ρ e^{iγ_ρ y} on (−L, L). A quadruple ½ ± δ ± iγ̄ contributes
4 cos(γ̄y) cosh(δy); an on-line pair 2 cos(γ₁y) + 2 cos(γ₂y).

## D2. Weil windows and the Le Cam distance
Q_Z(f) = ∫∫ G_Z(x − x') f(x) f̄(x') dx dx', f ∈ C_c^∞(−a, a), 2a ≤ L (N1's window c = e^{2a} uses primes ≤ c, so X = c);
Q_Z itself is unbounded on L². For two configurations differing in finitely many zeros, ΔG = G_{Z₁} − G_{Z₀} is a smooth
function and Q_{Z₁} − Q_{Z₀} extends to the bounded convolution operator T_{ΔG} on L²(−a, a).
Le Cam distance at window c: η_c(Z₀, Z₁) = ‖T_{ΔG}‖ = sup_{‖f‖=1} |Q_{Z₁}(f) − Q_{Z₀}(f)|.
Two-point lemma (deterministic): if η_c(Z₀, Z₁) ≤ m, every argument whose input is the Weil form on window c and whose
conclusion is invariant under every perturbation ΔQ with |ΔQ(f)| ≤ m‖f‖² (a certificate with margin m) concludes the same
for both. This is a robustness implication; the exact prime data of Z₀ and Z₁ still differ (D3).
For ΔG concentrated at frequency ±γ̄ (|γ̄|a ≫ 1), ‖T_{ΔG}‖ = ‖T_h‖ + O(1/(γ̄a)) with the demodulated
h(y) = e^{−iγ̄y} ΔG₊(y) (two nearly orthogonal components); lecam1.py computes ‖T_h‖ by Gauss–Legendre Nyström.

## D3. Finite changes are never exactly invisible
Σ_j c_j e^{iz_j y} (distinct z_j, finitely many) is entire in y; vanishing on an interval forces c_j = 0. So exact
indistinguishability needs infinitely many changed zeros (Beurling/DH-type); every finite construction has η > 0.

## D4. Why L* = ½ log(t/2π)
Moving ordinate γ_k by ε adds ε·iy e^{iγ_k y} to ΔG. The M neighbours of γ̄ (density D = log(t/2π)/2π) span a band of
width W = M/D; functions on (−L, L) with spectrum in that band have ≈ WL/π degrees of freedom (Landau). The count is
overdetermined when M > ML/(πD), i.e. L < πD = ½ log(t/2π) = L*. This is a heuristic for approximate cancellation (finite M
never cancels exactly, D3); the sharp threshold is a conjecture supported by the numerics, not a theorem. Lattice case exactly (explore30): an off-line pair lattice of step Δ is
invisible on |y| < 2π/Δ; with the zero density 2/Δ = D this is L < πD = L*. X* = e^{L*} = √(t/2π) = the Riemann–Siegel
main-sum length (Z(t) = 2Σ_{n ≤ √(t/2π)} n^{−1/2} cos(θ(t) − t log n) + R).
Bounded moves add a second condition: the merged-pair signal 2cosh(δy) − e^{iθ₁y} − e^{iθ₂y} grows like e^{δL}; with
|ε_k| ≤ 1 spacing it is cancelled only for δL ≲ 1.5 (measured, REPORT S2).

## D5. Exact global pair (Euler products)
A(s) = 1 + b p^{−s} + p^{1−2s} = (1 − αp^{−s})(1 − βp^{−s}), α + β = −b, αβ = p, b real, b ≠ −(p + 1) (that b makes A vanish
at s = 0, 1 and cancels ζ's pole, which changes the pole terms of the explicit formula). p^s A(s) = p^s + b + p^{1−s} is invariant
under s ↦ 1 − s, so F_{p,b} = ζA satisfies ξ_F(s) := p^s π^{−s/2} Γ(s/2) ζ(s) A(s) = ξ_F(1 − s): ζ's gamma factor,
conductor p², root number 1, Euler product. Coefficients c_n = Σ_{d|n} a_A(d) = 1 and Λ_F(n) = Λ(n) for n < p
(Λ_A lives on powers of p: Λ_A(p^k) = −(α^k + β^k) log p). Zeros of A: p^{s−½} = w, √p(w + 1/w) = −b: on the line iff
|b| ≤ 2√p, else Re s = ½ ± log|w|/log p (p = b = 7: 0.902588; p = 11, b = 11.5: 0.978492). Poisson summation over the two
vertical lattices: Σ_{A-zeros} ĝ(γ) = log p Σ_k (w^k + w^{−k}) g(k log p) (k = 0 term 2g(0)); for supp g ⊂ (−log p, log p)
this is 2 log p · g(0) for every b. Hence F_{p,b} and F_{p,b'} have identical prime data below p and
Q_{F_{p,b}} = Q_ζ + 2 log p ‖f‖² on every window c < p. On windows just above p, two bumps at distance log p give
Q_A(f) = log p·g(0)(2 − b/√p): negative iff b > 2√p. Hamburger (1921): a Dirichlet series absolutely convergent in Re s > 1, with
ζ's gamma factor and functional equation (conductor 1) and meromorphic continuation of finite order with finitely many
poles, is cζ, so no such pair exists for ζ itself; F_{p,b} violates only the
Selberg-class Euler-product growth axiom (b(p^k) ≪ p^{kθ}, θ < ½) at p.

## D6. Assouad
Hypercube ω ∈ {0, 1}^K (candidate pair j off-line iff ω_j = 1), all 2^K configurations η-close to Z₀ when the
compensation works for the full block: an η-robust estimator of #{off-line zeros in the block} computed from primes ≤ X
outputs one value on the whole cube, so its worst-case error is ≥ K/2 pairs; an η-robust estimate of which pairs are
off-line has expected Hamming error ≥ K/2 under the uniform cube (Assouad). Fano needs a noise model;
with exact data the mutual information is the entropy of the data map, i.e. counting — not used.
