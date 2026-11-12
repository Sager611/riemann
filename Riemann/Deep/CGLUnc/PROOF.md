# CGLUnc — statements and proofs (cglunc, 12 Nov 2026)

Notation as in REPORT.md: γ_ρ = (ρ−1/2)/i; reflection ρ ↦ 1−ρ̄ sends γ_ρ ↦ conj γ_ρ (equal multiplicity);
x-units Δ = (γ_ρ − γ_ρ')L/2π; K(x) = ∫K̂(α)e^{2πiαx}dα (K̂ even, smooth, K̂(α)e^{c|α|} ∈ L¹ for all c);
S(K) = Σ_{ρ,ρ'} m_ρ m_ρ' K(Δ_{ρρ'}) over zeros in the window. "Bandwidth-one data" = tr G̃ and ‖G̃_ψ‖²_F (all windows ψ on
[−1/2,1/2], with cross terms), i.e. S(K) for supp K̂ ⊂ [−1,1].

## Proposition U (what survives; known)
(i) F(α) ≥ 0 for all real α and F(α) = T^{−2α}(log T + O(1)) + α + O((log T)^{−1/2}) on [0,1], zeros complex
(Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh 2306.04799, Thm 1, TeX l.140–156; positivity = reflection symmetry:
Σ m m' e^{iξ(γ_ρ − γ_ρ')} = |Σ m e^{iξγ_ρ}|² because conj Σ m e^{iξγ_ρ} = Σ m e^{−iξγ_ρ}).
(ii) Hence for K̂ ≤ 0 on |α| ≥ 1: S(K)·w ≤ N(K̂(0) + ∫_{−1}^{1}K̂|α| + o(1)) unconditionally (CGL's upper bound).

## Proposition O (the CGL lower bound cannot be made unconditional; configuration-wise frame)
Kernel classes (K̂ even, K̂ ≤ 0 on |α| ≥ 1, K̂ ≢ 0 on |α| > 1):
 (C1) band-limited: K̂ ∈ C^k_c([−A, A]), A > 1 the edge of supp K̂ (so ∫_{|α|>A−ε}|K̂| > 0 for every ε > 0), k ≥ 4;
 (C2) CGL's own class: K̂(α) = f(rα), f = p(x)e^{−πx²}, p a polynomial, f eventually ≤ 0 (A_LP with the SDP ansatz).
(Both have all exponential moments, which S(K) = N∫K̂F needs for complex zeros at every T.) Put
I₂(y) = ∫_{|α|>1}|K̂|cosh(4παy) and Γ_k(η) = Σ_{j≤k}(1+|η|)^{k−j}∫|∂^jK̂(α)|e^{2π|α||η|}dα, so that
|K(d + iη)| ≤ c_k Γ_k(η)(1+|d|)^{−k} (k integrations by parts, x-units). Explicitly:
 (C1) Γ_k(η) ≪ (1+η)^k e^{2πAη} and I₂(y) ≫_ε e^{4π(A−ε)y};
 (C2) Γ_k(η) ≪ (1+η)^{m} e^{πη²/r²} and I₂(y) ≫ (1+y)^{−m} e^{4πy²/r²} (Laplace at α = 2y/r², where f(rα) < 0).
In both, (1+y)^k Γ_k(y + y₀)/I₂(y) → 0 exponentially in y (C1: ≪ (1+y)^k e^{−2π(A−2ε)y}; C2: ≪ e^{−3πy²/r² + O(y)}),
Γ_k(2y)/I₂(y) = e^{o(y)} (C1; ≪ (1+y)^k e^{4πεy} for every ε) resp. O((1+y)^{m'}) (C2), and e^{4πy}/I₂(y) → 0.
Claim: let Z be any configuration (reflection-symmetric, locally finite, N(t,t+1) ≪ log t) whose zeros all have depth
≤ y₀ (y₀ fixed; e.g. on-line configurations, AF's ceiling law, the TEN duals of Prop D), with N zeros in the window and
‖G̃_ψ(Z)‖²_F = O(N) for the bandwidth-one windows ψ in use. For every M > 0, θ ∈ [0, 1/2] and T large there is Z' ⊃ Z,
obtained by adding n_p ≍ N^θ off-line pairs (θ = 0: a bounded number), with S_{Z'}(K) ≤ S_Z(K) − (M − o(1))N,
N(Z') − N(Z) = O(N^θ), N0*(Z') = N0*(Z), N0s(Z') = N0s(Z), and tr G̃_ψ, ‖G̃_ψ‖²_F, tr(G̃_ψG̃_ψ') within o(N) of Z's.
Consequence (scope): no inequality S(K) ≥ Σ_{on-line ρ}c(m_ρ) − C·(bandwidth-one data) − o(N) with c(1) > 0 holds for
all configurations; every shallow (in particular every on-line) adversary of a configuration-wise certificate built on
bandwidth-one data and counts can be modified at o(N) cost to satisfy S(K) ≤ N·U(K), so the constraint cannot lower any
bound such adversaries give (e.g. AF's 0.6818287). Nothing is claimed about arguments with further arithmetic input.
Proof. A pair {ρ, 1−ρ̄}, m = 1, at depth y = (β−1/2)L/2π has self-terms 2K(0) + 2K(2iy) in S, and
K(2iy) = ∫K̂cosh(4παy) ≤ ‖K̂‖₁cosh(4πy) − I₂(y). Choose y = y(T) with I₂(y) ≍ MN^{1−θ} (I₂ is continuous and
increasing). Then y → ∞: in C1, c_εe^{4π(A−ε)y} ≤ I₂(y) ≤ ‖K̂‖₁cosh(4πAy) gives y = (1−θ)log N/(4πA') with
A' ∈ [A−ε, A] + o(1), so y ≍ log N, y ≤ L/(4π) (β ≤ 1) and e^{4πy} ≪ N^{(1−θ)/(A−ε)}; in C2, y ≍ (log N)^{1/2} and
e^{4πy} = N^{o(1)}. Insert n_p = ⌈MN/I₂(y)⌉ ≍ N^θ pairs, mutually ≥ D = N/(2n_p) apart (x-units). An x-interval of
length 1 holds ≪ log T zeros of Z, so each pair's cross terms with Z total ≪ Γ_k(y + y₀) log T = o(I₂(y)) (y ≫ log log T),
and the pair–pair terms total ≪ n_p²Γ_k(2y)D^{−k} ≪ n_p^{k+1}N^{−k}(n_pI₂(y))e^{o(y)} ≍ M N^{1−k(1−θ)+θ+o(1)} = o(N) (k ≥ 4).
Hence S_{Z'} − S_Z ≤ −2n_pI₂(y)(1 − o(1)) + o(N) ≤ −(2M − o(1))N (n_pI₂(y) ≥ MN).
Bandwidth-one data: with E the padding's part of G̃_ψ, Poisson–Gabor gives tr E = 2n_p; AF's Paley–Wiener bound (PW)
for the C² window φ gives an HS kernel ≪ e^{2π(|η|+|η'|)}(1+|d|)^{−4} between points at depths η, η' and horizontal
distance d, so ‖E‖²_F ≪ n_p e^{4πy}(1 + n_pD^{−4}) ≪ N^{θ + (1−θ)/(A−ε)} (C1) resp. N^{θ+o(1)} (C2) = o(N) (A − ε > 1),
and |‖G̃+E‖² − ‖G̃‖²| ≤ 2‖G̃‖‖E‖ + ‖E‖² = o(N); likewise for tr(G̃_ψG̃_ψ'). N grows by 2n_p; on-line counts unchanged. ∎
Zero density (what would block the padding). It places ≍ N^θ zeros at β − 1/2 = 2πy/L.
 (C1) β − 1/2 → (1−θ)/(2A'): a bound N(σ,T) ≪ T^{1−a(σ−1/2)}log T allows this iff θ ≤ 1 − a(1−θ)/(2A'), i.e. for every
 θ ∈ [0, 1/2] when a < 2A; a slope a > 2A excludes it (a = 2A undecided). Known: a = 8/7 (Conrey; Baluyot's thesis, θ < 4/7).
 Even θ = 0 works: ONE pair at β ≈ 1/2 + 1/(2A) lowers S(K) by MN at bandwidth-one cost N^{1/(A−ε)} = o(N).
 (C2) |β − 1/2| ≍ (log T)^{−1/2}: compatible with every fixed slope a; since θ = 0 is allowed, only a box hypothesis with
 no exceptions (every zero of the window with |β − 1/2| = o((log T)^{−1/2})) blocks it — the thin-box hypotheses
 |β − 1/2| ≤ b/log T of Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh and Goldston–Suriajaya are of this kind.

## Lemma P (the rank–trace frame has only positive-definite second moments)
If M_γ (γ ∈ ℝ) is a translation-covariant family of real symmetric matrices (e.g. Σ_j a_j v^{(j)}_γ v^{(j)ᵀ}_γ for
modulated windows on the Gabor grid), κ(γ − γ') := tr(M_γ M_γ') satisfies Σ c_i c̄_l κ(γ_i − γ_l) = ‖Σ c_i M_γi‖²_F ≥ 0,
so κ is positive definite and κ̂ ≥ 0 (Bochner). Then ∫_{|α|>1}κ̂F ≥ 0: F ≥ 0 bounds the unknown part of ‖Σ m M_ρ‖² from
below, while the rank–trace inequality needs it from above. CGL's improvement requires κ̂ < 0 beyond 1. ∎

## Proposition D (TEN-class dual)
Class (scalar affine-defect certificates). Window W: ψ = cos(√2 s) + Σ_{j≤J}c_j cos(2πjs) on [−1/2,1/2] (real c; mixtures
below), the arena's normalisation mtParams T = atAM(paramsOf stdProfile 1) T (Solution.lean l.4001, λ = 1), x = xnorm =
L(T)(Im ρ − T)/2π (l.4003), retained = simple on-line zeros of I' with L² ≤ x ≤ d − L² (l.4006–4008), Γ = gram (l.4013),
D∘ = Dcirc = Σ_i(λ_i(Γ) − 1)² (l.4019). A certificate of the class is any proof, from any local certificate (any block
size n, pair weights, gap pressures, floors, reflection rows), of the arena's ActualDefectEstimate (l.16111–16115):
for every ZeroConfig Z and δ > 0, eventually in T, (c − 4(n−1)δ)·#retained − B·span(xret) − (n−1)c ≤ D∘, with
0 ≤ c < 1; actual_defect_estimate (l.16204–16224) proves it for every Z : ZeroConfig (no PaperInputs). ZeroConfig is
not defined in Solution.lean but imported from the pinned zeta23 (imports l.223–246); its hypotheses (zeta23
Zeta23/Defs.lean l.120–130): carrier ⊂ ℂ, mult ≥ 1 on carrier (l.125), 0 ≤ Re ρ ≤ 1 (l.126), reflect ρ = 1 − ρ̄ ∈ carrier
(l.127, reflect l.108) with equal multiplicity (l.128), finitely many points with T₁ < Im ρ ≤ T₂ (l.130). The O(1) term
(n−1)c is independent of T and Z. Combined with AF's chain N0s ≥ #retained ≥ H_W N + D∘ − o(N), κ = (H_W − B)/(1 − c).
Dual lemma. Let X ⊂ [0,S) be p points (density ν = p/S, distinct mod S) and e_W(X) = lim (1/len)Σ_{i≠j}k_W(x_i−x_j)² its
energy per unit length (Poisson below). If e_W(X) ≤ ν − H_W, every certificate of the class has κ ≤ ν.
Proof. Let U(γ) = (γ/2π)(log(γ/2π) − 1), so U′(γ) = L(γ)/2π, and let Z_X have carrier {1/2 + iγ : γ ≥ γ₀,
U(γ) ∈ SZ + X}, mult ≡ 1. Z_X is a ZeroConfig: Re = 1/2 (l.126), reflect(1/2 + iγ) = 1/2 + iγ (l.127, l.128), U is
increasing so windows are finite (l.130). In the window at height T, xnorm = (U(γ) − U(T))(1 + O(1/L)) uniformly, so the
retained set is M(T) ≍ d/S periods of X under a map with derivative 1 + O(1/L): #retained/span → ν. Entries: kernel_limit_
uniform (l.13885–13889) gives |Γ_ij − k_W(x_i − x_j)| ≤ δ eventually; AF's Poisson–Gabor identity and |(φ²)^(Δ)| ≤
‖(φ²)′‖₁/|Δ| give |Γ_ij| ≤ C/(1 + |x_i − x_j|) and |Γ_ii − 1| ≪ L⁻⁴ uniformly in T (retained points are ≥ L² from the grid
ends), a summable majorant at bounded density; dominated convergence gives D∘/span → e_W(X). Divide the estimate by
span ≍ d → ∞ (the O(1) vanishes), let T → ∞ (M(T) → ∞), then δ → 0: cν − B ≤ e_W(X) ≤ ν − H_W, so
H_W − B ≤ (1 − c)ν and κ ≤ ν (c < 1). ∎
Poisson. Σ_n k²(y + nS) = S⁻¹Σ_m K̂(m/S)e^{2πimy/S} (K̂ = ψ⋆ψ/(∫ψ)² continuous, supported in [−1,1], k² = O(x⁻²)), so
e_W(X) = S⁻²[K̂(0)p² + 2Σ_{1≤m<S}K̂(m/S)|Σ_{x∈X}e^{2πimx/S}|²] − p/S (checked against a real-space sum to 2e−7).
Corollary D1 (arena window, Arb). X = (.xx)^{12}(.xx.xx.xxx) relaxed, S = 46, p = 31, positions rationalised to 10⁻⁹:
ν − H_W − e_W(X) = 7.1796e−5 > 0 (certify_dual_TEN.log; H_W = 0.67206299948881203266318 = the Lean constant). κ ≤ 31/46.
Corollary D2 (the cosine family, Arb). Write c = (1, c′), e_W(X; c) = cᵀQ_Xc − ν (Q_X real symmetric, (J+1)×(J+1)),
H(c) = H_MT − Σ_j w_j c_j², w_j = (1/2 − 1/(4π²j²))/(2sin²(1/√2)) (|u−v| is diagonal on cos(2πjs)). With the augmented
normalised Gram matrix G = [[1, c′ᵀ],[c′, c′c′ᵀ]] (G ⪰ 0, G₀₀ = 1): ν − H(c) − e_W(X;c) = tr(M_X G),
M_X = [[2ν − H_MT − Q₀₀, −qᵀ],[−q, diag(w) − Q′]]. Twelve configurations X_i of density ν = 153/227 and rational λ_i ≥ 0,
Σλ_i = 1 (exact, verify_cover.py) make M(λ) = Σλ_iM_{X_i} positive definite (Arb LDLᵀ pivots > 0, min 3.85e−4) for
J = 16, 48, and for J = 96 (orchestrator rerun, min pivot 3.8468e−4, 265 MB); J = 96 contains every J ≤ 96.
One common configuration. Let X̄_B be the S-periodic concatenation (S = B·227·Σa_i) of a_i blocks of B periods of X_i,
λ_i = a_i/Σa_j. Since |κ_j(x)| ≤ C_J/(1 + |x|) for the basis kernels (j ≤ J), the pairs straddling a block boundary
contribute ≤ C_J log(227B)/(227B) per unit length to each entry, so ‖Q_{X̄_B} − Σλ_iQ_{X_i}‖_op ≤ (J+1)C_J log(227B)/(227B)
→ 0 and M_{X̄_B} ≻ 0 for B ≥ B₀(J). Then tr(M_{X̄_B}G) ≥ λ_min(M_{X̄_B})·tr G ≥ λ_min > 0 for every G ⪰ 0 with G₀₀ = 1,
so the single configuration X̄_B satisfies the dual lemma for every window of the family, for every mixture
(G = Σμ_k c^{(k)}c^{(k)ᵀ}, Σμ_k = 1, energy Σμ_k e_k) and for the best of several windows: κ ≤ 153/227 = 0.6740088.
(Order: B fixed while T → ∞ in the dual lemma; B₀ comes from the margin, no second limit is needed.)
Complex coefficients only in Hermitian form: G = cc*, energies |k_c|² = Σ c_j c̄_k κ_jκ_k (D∘ = ‖Γ − I‖²_F = Σ_{i≠j}|Γ_ij|² + Σ_i|Γ_ii − 1|², the diagonal correction o(span); H(c) = H_MT − Σ_j w_j|c_j|²); tr(MG) > 0
holds there too; complex-symmetric k_c² is not covered (and AF's Hermitian Gram kernel is the transform of |φ|², so
complex profiles do not arise).
Not covered: windows outside this basis (other frequencies or profiles), other supports or scales (λ ≠ 1, other
intervals, other tapers), estimates proved only for zeta-like configurations (extra zeta-specific restrictions, e.g.
the explicit-formula data of other windows), deductions outside the scalar affine-defect class (non-affine defect
bounds, Q-side defects, cross terms). Covered: every certificate shape, block size and weight set that yields the same
universal affine-defect estimate.
Repairs (orchestrator referee): (1) class stated as the arena's ActualDefectEstimate for every ZeroConfig, ZeroConfig's
fields quoted, Z_X built and checked, O(1) uniform, limit order T → ∞ then δ → 0 (B fixed); (2) tr(MG) > 0 only for
G ⪰ 0 with G₀₀ = 1 (augmented normalised Gram), certified rational λ_i ≥ 0, Σλ_i = 1; (3) one common configuration via
operator-norm convergence of Q_{X̄_B} and the PD margin; (4) complex coefficients only via Hermitian G and |k|²;
(5) family narrowed to J ≤ 96 (J = 96 by the orchestrator) and uncovered freedoms listed. CLOSED re-review (orchestrator, 12 Nov 2026): items (1)–(3), (5) CLOSED; (4) STILL OPEN on the display D∘ = Σ|Γ_ij|², corrected here with the referee's own wording (D∘ = ‖Γ − I‖²_F, H(c) with |c_j|²); SIGNIFICANT as a tight finite-dimensional TEN-method barrier, not an unrestricted TEN ceiling.
