# explore11 derivations (D1–D5)

Notation: χ_12 = (12/·); ψ = ψ_4ψ_5, the even primitive character mod 20 with ψ_5(2) = i, so ψ(3) = i; ε = τ(ψ)/√20; φ = (1+√5)/2; T_k = Glaisher's numbers 1, 23, 1681, 257543, 67637281, ….

## D1. Strange identities: coefficients, Stokes constants, zero sensitivity
- For an f-periodic even χ with Σχ = 0, A(t) = Σ_{n≥1} nχ(n) e^{−n²t/D} has Mellin transform Γ(s) D^s L(2s−1, χ). That transform is entire apart from the poles of Γ, so A(t) ~ Σ_k (−1)^k L(−2k−1, χ) (t/D)^k / k!.
- Strange identity (Kontsevich–Zagier; Hikami for T(2,2m+1)): X(e^{−t}) = −½ e^{c²t/D} A(t) as formal power series.
  - Trefoil: X = Σ_k (q;q)_k, D = 24, c = 1, χ = χ_12.
  - T(2,5): X = Σ_{k₂≥k₁≥0} (q;q)_{k₂} q^{k₁²+k₁} [k₂ choose k₁]_q, D = 40, c = 3, χ = Im ψ.
  - Both checked exactly to t¹³ (t1). The partner with q^{k₁²} does not match −½·(±1, ±9 function) and is not used.
- FE (χ real, even, primitive, root number 1): L(−2k−1, χ) = (−1)^{k+1} (f/π)^{2k+3/2} k! Γ(k+3/2) π^{−1} L(2k+2, χ). Checked for k = 0, 1, 5, 20 (t2).
  - For χ_12 this gives L(2k+2, χ_12) = 2π T_k / ((12/π)^{2k+3/2} k! Γ(k+3/2)).
  - The k-th coefficient is therefore ∝ (2k+1)!(f/2π)^{2k+2} Σ_m χ(m) m^{−2k−2}. The Borel transform is a sum of simple-pole terms at points ∝ m², with weights χ(m). So the Stokes constants are the Dirichlet coefficients (the single-tower case of Fantini–Rella 2404.11550).
  - Peeling T_0…T_139 recovers χ_12(1…13) exactly (t2).
- Zero sensitivity. Move a zero ρ = ½+iγ and its conjugate to the FE-symmetric quadruple ½±δ±iγ. This multiplies Λ by R(s) = (1 − δ²/(s−ρ)²)(1 − δ²/(s−ρ̄)²).
  - At s = 2k+2: R − 1 ≈ −2δ² Re(2k+3/2−iγ)^{−2} = O(δ²/k²).
  - For the lowest zero of L(s, χ_12), γ_1 = 3.80462763305087: relative change of T_k ≈ +0.087δ² (k = 0), −0.016δ² (k = 2), −0.011δ² (k = 5), −0.0011δ² (k = 20) (t2).
  - So T_0 = 1 would become 1.00087 at δ = 0.1. Integrality pins the function (R37 rigidity), not the class of its zero locations (D3).

## D2. Character content of torus knots and Brieskorn spheres
- Torus knot T(s,t): χ_{2st} is +1 at ±(st−s−t), −1 at ±(st+s−t), 0 elsewhere. These are 4 classes, all units mod 2st.
  - A Dirichlet character is nonzero on all φ(2st) units. φ(2st) = 4 ⇔ 2st ∈ {8, 10, 12} ⇔ (s,t) = (2,3).
  - So every torus knot except the trefoil carries a combination of ≥ 2 characters (t5):
    - T(2,5): ψ, ψ̄ (conductor 20)
    - T(3,4): (12/·), (24/·)
    - T(2,7): 3 characters, conductor 28
    - T(3,5): 2 characters, conductor 15
    - T(4,5): 4 characters
- Brieskorn sphere Σ(p₁,p₂,p₃): Hikami's χ_{2P}, n ≡ P(1 + Σε_i/p_i) ↦ −ε₁ε₂ε₃, has 8 classes. Since P ≥ 30, φ(2P) ≥ 16 > 8: never a single character (Σ(2,3,5): two odd characters of conductor 60).
- Zeros in Re s > 1. Saias–Weingartner (Acta Arith. 2009): a periodic Dirichlet series with no zeros in Re s > 1 equals P(s)L(s,χ), with P a Dirichlet polynomial.
  - If L(s, χ_{2st}) = P(s)L(s,χ), then a(ℓ) = p(1)χ(ℓ) for primes ℓ > deg P.
  - χ_{2st} vanishes on some unit classes and not on others. By Dirichlet's theorem this forces p(1) = 0 and p(1) ≠ 0, a contradiction.
  - Hence for every torus knot ≠ trefoil (and every Brieskorn sphere), the L-series of the strange identity has zeros in Re s > 1, infinitely many by almost periodicity.

## D3. T(2,5) and the Lee–Yang S-matrix
- ε(ψ) = (2/√5)(sin 72° + i sin 36°), so tan arg ε = sin 36°/sin 72° = 1/φ.
- Let θ_χ(τ) = Σ_{n≥1} χ(n) q^{n²/40} and Θ = (θ_{Re ψ}, θ_{Im ψ}). Hecke's θ_ψ(−1/τ) = ε (τ/i)^{½} θ_ψ̄(τ) gives Θ(−1/τ) = (τ/i)^{½} S Θ(τ), with S = [[c, s], [s, −c]] and c + is = ε.
  - This is the M(2,5) (Lee–Yang) S-matrix. Fitted numerically, the residual is 4·10⁻¹⁶ (t3c).
- S is real, symmetric and orthogonal, with eigenvalues ±1. Its eigenvectors (cos θ/2, sin θ/2) and (−sin θ/2, cos θ/2), θ = arg ε, give
  - G₊(s) = Σ 2Re(e^{−iθ/2}ψ(n)) n^{−s} and H₋(s) = Σ 2Im(e^{−iθ/2}ψ(n)) n^{−s},
  - with Λ(1−s) = ±Λ(s), where Λ(s) = (20/π)^{s/2}Γ(s/2)·(series).
  - These are Davenport–Heilbronn-type combinations, of the same kind as Nakamura's f(s,χ) (2008.02570, infinitely many zeros off the line for non-real χ).
  - Saias–Weingartner applies directly: |a(ℓ)| takes the two values 2cos(θ/2) and 2sin(θ/2) on primes in different classes, so a(ℓ) = p(1)χ(ℓ) is impossible. Hence there are infinitely many zeros in Re s > 1.
- First zeros right of the line (t3b, mp.dirichlet check at 40 digits):
  - H₋: 1.195270648642 + 5.539314643116i (in Re s > 1)
  - G₊: 0.652400414407 + 122.846724764454i; 0.5897221094 + 374.87624634i
  - The knot's own series L(s, Im ψ): 0.838544451285 + 8.292146000549i
  - T(3,4): 0.986306813311 + 6.982295921607i
  - T(2,7): 0.820173484809 + 9.766689580103i
  - Argument-principle counts in [0.51, 2.6] × [0.5, 60]: 17 (T(2,5)), 10 (H₋), 16 (T(3,4)), 12 (T(2,7)); G₊: 1 in [0.51, 2.6] × [0.5, 123.8].
- Poincaré sphere: χ_60 = −Re ψ′, where ψ′ is odd and primitive mod 60 and ε(ψ′) = 0.525731 − 0.850651i (golden ratio again). Its H₋ vanishes at 0.604987666972 + 19.746196916245i (t3d, t3e).
- Principle: a unitary, real-symmetric S-matrix on a vector of theta series yields ±1 eigencombinations with scalar FEs. Unless the vector is one-dimensional with multiplicative coefficients, these are Davenport–Heilbronn functions. Unitarity builds DH; it does not exclude it.

## D4. Riesz criterion in Glaisher numbers (known)
- R(x) = Σ_{k≥0} (−x)^{k+1}/(k! L(2k+2, χ_12)) = −x Σ_n μ(n)χ_12(n) n^{−2} e^{−x/n²}, and L(2k+2, χ_12) is the rational multiple of π^{−2k−1/2} and T_k from D1.
- GRH(χ_12) ⇔ R(x) = O(x^{1/4+ε}). This is Riesz 1916 (NOGO R13.4 ≡ RH; Riesz sufficiency machine-checked for ζ in this repo).
- A zero at height γ enters with weight |Γ(1−ρ/2)| ≈ e^{−πγ/4}.

## D5. Kashaev invariant of 4₁
- ⟨4₁⟩_N = Σ_{k<N} Π_{j≤k} (2 sin(πj/N))² > 0.
- log⟨4₁⟩_N = NV/2π + (3/2)log N − ¼ log 3 + 0.5544/N + O(N⁻²) (t4, N ≤ 3200), with V = vol(4₁) = 2D(e^{iπ/3}) = (3√3/2) L(2, χ_{−3}) = 2.029883212819307.
- The 1/N term matches Garoufalidis–Zagier's 11/(72√−3) at h = 2πi/N: 11·2π/(72√3) = 0.5542. The only L-value is at s = 2, inside absolute convergence.
