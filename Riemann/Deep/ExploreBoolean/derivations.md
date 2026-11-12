# ExploreBoolean — derivations (explore33, 12 Nov 2026)

Notation. x ≥ 2, P_x = primes ≤ x, n = π(x), cube {0,1}^{P_x} ≅ subsets S, m_S = Π_{p∈S} p. Knapsack down-set D_x = {S : Σ_{p∈S} log p ≤ log x} ≅ squarefree m ≤ x; f_x = 1_{D_x}, α = |D_x|/2^n = Q(x)/2^n. Integer Walsh coefficients F(T) := 2^n \hat f(T) = Σ_{S∈D} (−1)^{|S∩T|}. Generalized primes q_1 < q_2 < … give D_q(x) = {S : Π_{i∈S} q_i ≤ x}, M_q(x) = Σ_{S∈D_q(x)} (−1)^{|S|}.

## D1 Spectrum = all sign-twisted squarefree sums
F(T) = Σ_{m≤x} μ²(m) g_T(m), g_T completely multiplicative, g_T(p) = −1 (p ∈ T), +1 (p ∉ T). So F(∅) = Q(x), F(P_x) = M(x) (μ = (−1)^ω on squarefree m), and Parseval Σ_T F(T)² = 2^n Q(x): the RMS coefficient is √Q(x) ≍ √x (the Rademacher random multiplicative function in L²). Hence

  RH ⟺ |F(P_x)| ≤ x^ε · RMS_T F(T) for every ε > 0, x ≥ x_0(ε):

the all-minus point is not an outlier among the 2^n Rademacher points. The class C_LTF⁺ (positive-weight threshold down-sets) is closed under coordinate permutations; f_x ∈ C_LTF⁺ is also shifted (Björner) and monotone.

## D2 The four tools at the top coefficient (exact statements: O'Donnell 2105.10386)
(a) Hypercontractivity. ρ^n |\hat f([n])| ≤ ‖T_ρ f‖_2 ≤ ‖f‖_{1+ρ²} = α^{1/(1+ρ²)}. B(ρ) = ρ^{−n} α^{1/(1+ρ²)}: (ρ/n) d log B/dρ = −1 + 2(ρ²/(1+ρ²)²)(ln 2 − ln N/n) ≤ −1 + (ln 2)/2 < 0, so the minimum is at ρ = 1: |M| ≤ √(2^n Q). Parseval; ≫ Q.
(b) Level-k (k ≤ 2 ln(1/α)): W^{≤k} ≤ ((2e/k) ln(1/α))^k α². At k = n (valid iff (2 ln 2 − 1) n ≥ 2 ln Q, first at n = 20): |M| ≤ Q ((2e/n) ln(1/α))^{n/2} → Q (2e ln 2)^{n/2} = Q·1.941^n. Any inequality valid for all Boolean f must exceed its value on Hamming balls, whose top coefficient is ±C(n−1, r) = |D|(1 − 2r/n + O((r/n)²)) (→ |D|(1 − 2c) for r = cn, c < 1/2; D6): no level inequality can beat |M| ≲ |D| at k = n.
(c) Peres. F_± = 1 − 2f is an LTF; NS_δ ≤ √(1/⌊1/δ⌋) (O'Donnell's proof); W^{≥1/δ} ≤ (2/(1 − e^{−2})) NS_δ ≤ 3 NS_δ. δ = 1/n: \hat F_±([n])² ≤ 3/√n, \hat F_±([n]) = −2\hat f([n]): |M| ≤ (√3/2) 2^n n^{−1/4}. Additive in μ: void when α ≪ n^{−1/4}.
(d) KKL / Friedgut. Lower bounds on influences and junta approximation; a junta h on J ⊊ [n] has \hat h([n]) = 0, so |\hat f([n])| ≤ ‖f − h‖_2 ≤ √α at best. For f_x: Inf_p = 2b_p/2^n with b_p = #{m ∈ D: p∤m, pm > x} = Q − 2c_p (Chow parameter = boundary count, monotone f); b_p = Q − 2 for every p > x/2, I[f]/α ≈ 2n: total influence ≈ maximal. With absolute error ε Friedgut only says F_± is ε-close to a constant (I[F_±] ≈ 2nα ≪ ε); with error εα relative to the measure the junta size exp(O(n/ε)) is vacuous.
(e) Topology (Björner 1101.5704). D_x shifted ⇒ wedge of spheres, Σβ_k = #{S ∈ D: 2 ∉ S, S∪{2} ∉ D} = #{odd squarefree m ∈ (x/2, x]} = (2/π²)x + O(x^θ). |M| ≤ Σβ: the identity M(x) = Σ_{x/2<m≤x, 2∤m} μ(m). Only bound below trivial, by the factor 1/3.

## D3 Biased measure makes f_x balanced, at an exponential price
Product measure π with Pr[p ∈ S] = q_p := 1/(p+1): π(S) = C_x/m_S, C_x = Π_{p≤x} p/(p+1), E_π f_x = C_x Σ_{m≤x} μ²(m)/m → e^{−γ} (Dickman law of log m_S/log x). Orthonormal basis φ_p(1) = √p, φ_p(0) = −1/√p. The sign weight (−1)^{|S|}/π(S) = Π_p g_p with E_π g_p = 0 and g_p = −((p+1)/√p) φ_p, so

  M(x) = Π_{p≤x} (−(p+1)/√p) · \tilde f_π(P_x),   Π (p+1)/√p = e^{θ(x)/2} Π(1 + 1/p).

Checked exactly for x = 10, 30, 60, 100 (biased_check.py; |\tilde f_π| = 1.3·10^{−19} at x = 100). Peres applies with the same bound (it holds for LTFs in product spaces, O'Donnell Ch. 8 exercise); KKL and Friedgut hold only with losses in the smallest bias λ = 1/(x+1) (general hypercontractivity). Either way RH needs |\tilde f_π(P_x)| ≤ x^{1/2+ε} e^{−θ(x)/2}: no Fourier-tail bound sees e^{−x/2}. For general biases the amplification is Π 1/√(q_p(1−q_p)) ≥ 2^n, minimal at q = 1/2: no product basis helps.

## D4 Symmetrization: the top coefficient sees only the ω-profile
F(P_x) is invariant under permuting coordinates, so F(P_x) = Σ_k (−1)^k A_k, A_k = #{S ∈ D: |S| = k} = π_k^*(x): M(x) = Σ_k (−1)^k π_k^*(x) = P_x(−1), P_x(z) = Σ_{m≤x} μ²(m) z^{ω(m)} (explore1's prime-factor gas). Level profile: for T random with Pr[p ∈ T] = α, the coefficient is a Selberg–Delange sum with mean sign z = 1 − 2α:

  F(T) ≈ G_T x (log x)^{−2α}/Γ(1 − 2α)   (main term; vanishes at α = 1/2 and α = 1).

Measured at x = 10^6 (median of 400 samples): α = 0.75: −4834 (pred. −5490 G_T, 1/Γ(−1/2) = −0.282); α = 0.9: −1468 (pred. −1544 G_T); α = 1/2: median 4, RMS 761 = √Q·0.976. The top coefficient sits at the second zero of 1/Γ(1 − 2α): Mertens' cancellation is the Selberg–Delange zero at z = −1, i.e. Perron's formula with the truncated Euler product Π_{p≤x}(1 − p^{−s}).

## D5 Barrier theorem: a planted zero is invisible to bounded-level Walsh data
Fix 1/2 < β < 1, γ > 0, ρ = β + iγ, u_0 ≥ 2^{1/(1−β)}. Let q_i be the primes ≤ u_0, then the points Π_c^{−1}(k − 1/2), k > π(u_0), where Π_c(t) = π(u_0) + ∫_{u_0}^t (1 − 2u^{β−1} cos(γ log u)) du/log u (density ≥ 0).
(a) P(s) := Σ q_i^{−s} = −log(s−1) + log(s−ρ) + log(s−ρ̄) + H(s), H analytic on σ > 0: the continuous part is E_1((s−1)log u_0) − E_1((s−ρ)log u_0) − E_1((s−ρ̄)log u_0) (E_1(z) = −γ − log z + entire); the rounding error |Π_d − Π_c| ≤ 1/2 gives s∫(Π_d − Π_c)u^{−s−1}du, analytic on σ > 0. With Σ_{k≥2} P(ks)/k analytic on σ > 1/2:

  Π_i (1 − q_i^{−s}) = (s − 1) e^{A(s)} / ((s − ρ)(s − ρ̄)),   A analytic on σ > 1/2.

e^{A} ≠ 0, so the Mellin transform s∫ M_q(x) x^{−s−1} dx has poles at ρ, ρ̄; if M_q(x) = O(x^{β−ε}) the transform would be analytic on σ > β − ε, so M_q(x) ≠ O(x^{β−ε}) for every ε > 0.
(b) For every finite T, Σ_{S} (−1)^{|S∩T|} m_S^{−s} = Π_{i∈T} (1 − q_i^{−s})/(1 + q_i^{−s}) · (s−ρ)(s−ρ̄) e^{B(s)}/(s−1), B analytic on σ > 1/2: every bounded-level Walsh coefficient has, on σ > 1/2, the single pole s = 1, exactly as for the primes (Π_{p∈T}(…)·ζ(s)/ζ(2s)). On σ > 1/2 the planted ρ creates no singularity of these transforms: it appears in the zero factor (s−ρ)(s−ρ̄) (no x^ρ term), inside the analytic zero-free factor e^{B(s)} (through H and the k ≥ 2 terms) and in the residue at s = 1 (the density constant); its first singularities come from P(ks), k ≥ 2, at s = ρ/2, ρ̄/2, terms of order x^{β/2} < x^{1/2}.
(c) Numbers (x ≤ 10^8, β = 0.8, γ = 10, u_0 = 40): fit M_q(x)/x^{0.8} = 0.0201 cos(10 log x − 1.456), R² = 0.988 on [10^5, 10^8], amplitude 0.0200 / 0.0201 on the two halves, best grid (β, γ) = (0.79, 10.00); |M_q|/√x decade maxima 0.74, 1.22, 2.43, 5.04 (×10^{0.3} per decade). Boolean statistics vs the primes at 10^8: N/x 0.5893 vs 0.6079, W¹/(2N² ln(1/α)) 0.7213 vs 0.7213 (→ 1/(2 ln 2) for every sparse knapsack set), I/(2nα) 1.0000 vs 1.0000, MaxInf/α 2 vs 2, Σβ/N 0.3333 vs 0.3333, mean ω 2.860 vs 2.888, var ω 1.191 vs 1.212; relative Hamming distance |D_x Δ D_q|/|D_x| = 0.054. Top: |M|/√N = 5.73 vs 0.25.

## D6 Equal weights: no cancellation
Weights 1 + ε_i (ε_i small, generic), threshold 2 + η: D = {∅} ∪ singletons ∪ {pairs with ε_i + ε_j ≤ η}, an LTF of every size N ∈ [n+1, n+1+C(n,2)]; Σ_{S∈D}(−1)^{|S|} = N − 2n. With n = π(x), N = Q(x): top = Q − 2π(x) ≈ Q(1 − 3.29/log x) (0.81 Q at 10^8). Hence any bound depending only on (n, |D|) and membership in C_LTF⁺ is ≥ |D| − 2n.

## D7 Perturbed weights
- q_p = p + 1: Π(1 − (p+1)^{−s}) = E(s)/ζ(s), E = Π (1 − (p+1)^{−s})/(1 − p^{−s}) analytic and zero-free on σ > 0 (Σ|s|p^{−σ−1} < ∞). Off-line ζ-zeros force M_{p+1} = Ω(x^{β−ε}); M_{p+1}(x) ≪ x^{1/2+ε} ⟹ RH. Yet |D_x Δ D_{p+1}| = 0.224 |D_x|.
- Shift q_p = p e^{−ε}: Σ_m μ(m)(m e^{−εω(m)})^{−s} = Π(1 − e^{εs}p^{−s}) = ζ(s)^{−e^{εs}} G(s), G = Π(1 − zp^{−s})(1 − p^{−s})^{−z}, z = e^{εs}, analytic zero-free on σ > 1/2 (ε < log 2). Selberg–Delange (exponent analytic in s): M ≈ (G(1)/Γ(−e^ε)) x (log x)^{−1−e^ε}. ε = 0.25: lead 585, 4113, 30316 at 10^6, 10^7, 10^8; observed 578, 5114, 36003 (second-order term ≈ 0.73 log log x/log x and √x-size ζ-oscillation).
- iid q_p = p e^{δ_p}, δ_p bounded (Selberg–Delange heuristic, not proved): E M_q(x) = Perron integral of Π(1 − K(s)p^{−s}), K(s) = E e^{−sδ} (Pakianathan–Winfree 1104.4324, Thm ExpTop); ≈ (G/Γ(−K(1))) x (log x)^{−1−K(1)}. δ ~ U[−0.5, 0.5]: K(1) = 1.0422, lead 10940 at 10^8; one realization 6635 (√x = 10^4: invisible below ~10^{10}). Valid as stated when K(1) ∉ ℕ and G(1) ≠ 0, e.g. every small symmetric nondegenerate δ (K(1) = E e^{−δ} ∈ (1, 2), all factors 1 − K/p ≠ 0). If K(1) = k ∈ ℕ, 1/Γ(−k) = 0; with X = e^{−δ}, K'(1) = E[X log X] ≥ k log k (Jensen), strictly > 0 for nondegenerate δ also when k = 1, so the first term is K'(1)(s−1)^{k+1} log(s−1) → x/(log x)^{k+2}, unless k is prime (then the factor 1 − K(1)/p vanishes at p = k, G(1) = 0, and further terms are needed). Exact prime density (K(1) = 1) therefore still gives E M_q ≍ x/(log x)³. In the generic case E M_q ≍ x/(log x)^{O(1)}: generic weights lose the power saving; log p is not generic.

## D8 Noise stability sees only t ≥ 0 (heuristic obstruction, not a theorem)
Stab_ρ[f] = Σ_k ρ^k W^{=k} = 2^{−n}((1+ρ)/2)^n Σ_{m,m'∈D} t^{d(m,m')}, t = (1−ρ)/(1+ρ) ∈ [0, ∞) for ρ ∈ (−1, 1], d = ω(mm'/(m,m')²); M² = Σ_{m,m'} (−1)^{d} is the same pair series at t = −1 (ρ → ∞: \hat f([n])² = lim Stab_ρ/ρ^n). The generic route from data on [−1, 1] to the leading coefficient is Chebyshev's inequality: Stab_ρ known to ±η certifies \hat f([n])² only to ±η 2^{n−1}, so certifying |M| ≤ √x this way needs η ≤ x 2^{1−3n} (relative precision 2^{1−2n}). This bounds what the data certify by that route; it is not a proof that constrained stability data cannot do better (the proved barriers are D5, D6).
