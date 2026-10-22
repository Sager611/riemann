# RSFam: Rankin–Selberg twists by class group characters: large sieve, zero density, and Hecke-eigenvalue conditions on primes in ideal classes (proofs)

Date 2026-10-21. Brief: `.scratch/orca/rsfam-brief.md`. Scratch: `.scratch/rsfam/` (`params.py` → `params.out`: every parameter
inequality of §5 at sample A, and the minorants of Corollary D1; `thresh.py`, `thresh2.py`, `hyp.py` → `hyp.out`: thresholds; `axq*.py`
→ `axq1.out`…`axq4.out`: 35 arXiv API queries; `gxl.py` → `gxl.out`, `gxl_*.json`: 4 evidence-API claims; `src_*`: TeX of HT, BTZ, LOT, PT).
References: [C] = [Combos NOTE](../Combos/NOTE.md) (Lemmas 2–5), [R] = [RCApps NOTE](../RCApps/NOTE.md) (§0 notation, (I1)–(I6),
§2), [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md) (Proposition 1, Theorems 2′, 4), [D2] = [Density2 NOTE](../Density2/NOTE.md) (Lemma 5).
Nothing here uses the 2026 half-planes. Paper proofs, not refereed, no Lean.
Repaired 2026-10-21 after review (brief `.scratch/orca/rsfamfix-brief.md`, scratch `.scratch/rsfamfix/`): five items, listed in
§10 with the places changed; `params2.py` → `params2.out` re-checks §5 with the corrected margins and the minorant constants.
Second novelty audit: [NOVELTY.md](NOVELTY.md).

## Verdicts

| Item | Result | Status |
|---|---|---|
| (a) large sieve | Lemma A: Σ_χ\|Σ_𝔞 λ_m(𝔞)b_𝔞χ(𝔞)\|² ≤ (30/ε)log²(3\|D\|)(N + \|D\|^{1/2})Σ m(𝔞)^{1+2ε}d_{m+1}(𝔞)²\|b_𝔞\|², λ_m the base-change coefficients of (R1′) (not λ_{Sym^m f}(N𝔞)) | routine ([C] Lemma 2 with \|λ_m(𝔞)\| ≤ d_{m+1}(𝔞), (R1′)) |
| (b) second moment | Lemma B: Σ_χ\|L(½+iv, Sym^m f × θ_χ)M_V\|² ≪ \|D\|^{(m+1)/2+ε}V^{1+ε}(\|v\|+3)^{m+1+ε}, uniformly in V ≥ 1, v ∈ R; saves h ≍ \|D\|^{1/2} over convexity summed | routine; proof repaired (cutoff enlarged by (\|D\|V(\|v\|+3))^ε); used only with V = 1, in (c) and (c′) |
| (c) per-D density | Theorem C: Σ_{χ∈Ĝ}N(σ;t, Sym^m f × θ_χ) ≪ (\|D\|^{m+2}T^{2m+2})^{(1−σ)/(3−2σ)+ε}, 0.6 ≤ σ ≤ 1, every m ≥ 1, zeros in the unit window \|γ − t\| ≤ 1; all \|γ\| ≤ T costs a factor T (C′). For f (m = 1): \|D\|^{3(1−σ)/(3−2σ)} per window, nontrivial for σ > 3/4 | new statement; best published for this family: \|D\|^{(m+7)(1−σ)/(3−2σ)} (Pascadi–Thorner 2025, Thm 1.3), nontrivial only for σ > (2m+11)/(2m+12) |
| (c′) density over D | Proposition F: Σ_{\|D\|~X}Σ_{χ²≠1}N(σ,T; Sym^m f × θ_χ) ≪ (X^{m+10}T^{2m+10})^{(1−σ)/(3−2σ)+ε} | new statement; Pascadi–Thorner give X^{(m+11)(1−σ)/(3−2σ)} |
| (d) application | Theorem D: for every A > 3, all but O(X^{1−c}) fundamental D < 0, \|D\| ≤ X: every ideal class has degree-one primes 𝔭 with N𝔭 = p ∈ (x, x+x\|D\|^{−κ}], λ_f(p) > 0, and others with λ_f(p) < 0, for all x ≥ \|D\|^A. Positivity comes from the m = 0 lower bound ≥ ν_DH/8 together with m = 1, 2 errors ≤ η₀ν_DH (§5, Positivity). Any set S admitting a Chebyshev minorant of degree ≤ 2 with positive mean: A > 3; degree 1 (S ⊃ (α, 2], α < 0): A > 9/4. Same for D > 0 (Theorem D⁺, indefinite forms) | **new** (no statement found with an explicit power of \|D\|, on average or not; a per-D bound with an inexplicit power of order 10^{10} looks derivable from published inputs, [NOVELTY.md](NOVELTY.md)) |
| thresholds | A₀ = 3/2, 9/4, 3, (M+7)/3 for minorant degree M = 0, 1, 2, ≥ 2; the published densities would give 8/3, 3, 10/3, (M+8)/3 | M ≥ 3 now unconditional: Thorner (IMRN 2025, arXiv 2404.06482, Thm 1.1 with n = 0) excludes Landau–Siegel zeros of Sym^m f ⊗ χ for every m (§7) |
| below 3 for the sign | at m = 2, A = 3, σ = 5/6 the trivial per-D count, Theorem C, Proposition F and the Markov bar all equal X^{3/2} | NOGO with present tools (§7, §9) |

Effective constants: c, κ, X₀ are effective but not numerical, because the constants of (J4) (Brumley–Thorner–Zaman) and (J5)
(no Landau–Siegel zeros for the Dirichlet twists Sym^m f ⊗ χ_d) are not explicit. The thresholds A₀ are exact.

## 0. Setup

- f ∈ S_k(Γ₀(N)): a holomorphic newform, k ≥ 2, trivial nebentypus, without CM; λ_f(n) := a_f(n)n^{−(k−1)/2}, so
  λ_f(p) = 2cos θ_p ∈ [−2, 2] for p ∤ N (Deligne). For E/Q without CM take f = f_E: a_E(p) = √p λ_f(p).
- π_f: the automorphic representation of f; Π_m := Sym^m π_f (m ≥ 1), cuspidal on GL_{m+1}(A_Q) (Gelbart–Jacquet m = 2;
  Kim–Shahidi m = 3; Kim m = 4; Newton–Thorne, Publ. IHÉS 2021, all m; Theorems D, D⁺ use m ≤ 2, Theorem D_M of §7 uses m ≤ M). Π_{m,v} = Sym^m π_{f,v}
  under the local Langlands correspondence, so Π_m is tempered at every place (π_{f,p} is tempered at every p by Deligne and
  Carayol's local–global compatibility) and self-dual (trivial central character). Π_0 := 𝟙. For p ∤ N the trace of Π_{m,p} is
  U_m(cos θ_p) = sin((m+1)θ_p)/sin θ_p (Chebyshev); U_1(x/2) = x, U_2(x/2) = x² − 1.
- 𝔇, 𝔇(X) = {D ∈ 𝔇 : X/2 < |D| ≤ X}, K = Q(√D), G = Cl(K), h = |G|, Ĝ, Ĝ₂ (genus characters), Ĝ_c = Ĝ ∖ Ĝ₂, 𝔉(D), χ_d,
  ε_d(C) and genus theory (G): as in [R] §0 (D fundamental, so every χ ∈ Ĝ is unramified). d_k(𝔞), m(𝔞): as in [C] §0.
- θ_χ := Σ_𝔞 χ(𝔞)q^{N𝔞} (Hecke; weight 1, level |D|, nebentypus χ_D), π_χ its automorphic representation of GL₂(A_Q):
  cuspidal iff χ ∈ Ĝ_c; π_χ = χ_d ⊞ χ_{D/d} if χ = χ_{d,D/d} ∈ Ĝ₂. π_χ ≅ π_χ̄ ≅ π̃_χ, since χ̄(𝔞) = χ(𝔞̄).
- L_m(s,χ) := L(s, Π_m × π_χ) = L(s, Sym^m f ⊗ θ_χ); for χ = χ_{d,D/d} it is L(s, Π_m ⊗ χ_d)L(s, Π_m ⊗ χ_{D/d}).
  L_1(s,χ) is the Rankin–Selberg L-function of the brief, of degree 4 and conductor ≍ |D|²; there are h ≍ |D|^{1/2±ε} of them.
- N_m(σ;t,χ) := #{ρ = β+iγ : L_m(ρ,χ) = 0, β ≥ σ, |γ − t| ≤ 1}; N_m(σ,T;χ): the same with |γ| ≤ T;
  N_D^{(m)}(σ;t) := Σ_{χ∈Ĝ}N_m(σ;t,χ); N_D^{c,m}(σ,T) := Σ_{χ∈Ĝ_c}N_m(σ,T;χ). Zeros with multiplicity. T := |t| + 3.
- C(ε), c(ε), C, c: effective constants depending on ε, f, m only (and on A, P where stated).

Inputs besides [C], [R], [R2] (exact statements used):
- (J1) Rankin–Selberg theory (Jacquet–Piatetski-Shapiro–Shalika; Mœglin–Waldspurger; Gelbart–Shahidi, JAMS 2001): for cuspidal
  π, π′ on GL_n, GL_{n′} with π′ ≇ π̃|·|^{iu}, L(s, π × π′) is entire of order 1, bounded in vertical strips of bounded width,
  with a functional equation of the usual shape.
- (J2) Local factors: under the local Langlands correspondence for GL_n (Harris–Taylor, Henniart) the Rankin–Selberg factors
  are Artin factors of tensor products, which are inductive.
- (J3) Humphries–Thorner, Compos. Math. 160 (2024) (arXiv 2103.05634v2), Theorem 4.1, Corollary 4.3 (cor:MVT) and inequality (bound_int)
  of §5, with n = 2, F = Q, θ₂ = 0. Statement used: let 𝒮 be a finite set of cuspidal representations of GL₂(A_Q) satisfying
  GRC, Q ≥ max_{π∈𝒮} C(π), π′ ∈ 𝔉_{n′} with π′ ≇ π̃ for all π ∈ 𝒮, T ≥ 3, ε > 0, V := C(π′)^εQ^{2+ε}T^{3+ε}|𝒮|. For
  0 ≤ σ ≤ 1 and every Y ≥ 3 with log Y ≍ log V:
  Σ_{π∈𝒮}N_{π×π′}(σ,T) ≪_{ε,n′} V^ε[Y^{2(1−σ)} + Y^{1/2−σ}(V·I₂)^{1/2}],  I₂ := Σ_{π∈𝒮}∫_{−T}^{T}|L^{ur}(½+iv, π×π′)|²dv.
  (The mean value threshold of their Theorem 4.1 does not depend on π′; the κ-term of (bound_int) vanishes since π′ ≇ π̃.
  HT bound I₂ by convexity; we do not.) This is (I3) of [R] with π′ = Π_m.
- (J4) Brumley–Thorner–Zaman, J. Eur. Math. Soc. 24 (2022) (arXiv 1804.06402v3), Theorems 1.2 and 1.4 (labels thm:lfzde_RS,
  thm:Landau_Siegel; other papers cite the published numbering), with the family of
  primitive Dirichlet characters (d = 1; Hypothesis 1.1 holds by GRC) and π₀ = Π_m (GRC): with c_m := 10^7(m+1)^4 and an
  effective B_m, for Q, T ≥ 1 and ½ ≤ σ ≤ 1,
  Σ_{χ primitive, cond(χ) ≤ Q} N(σ,T; Π_m ⊗ χ) ≤ B_m ν^♯ (C(Π_m)QT)^{c_m(1−σ)},
  where ν^♯ = 1, or ν^♯ = min{1, (1−β_ψ)log(QT)} for any real primitive ψ of conductor ≤ 2Q with a real zero β_ψ (log-free;
  the second form is their Landau–Siegel repulsion).
- (J5) No Landau–Siegel zeros for these twists. Iwaniec–Kowalski Thm 5.10 gives the region below for all zeros except possibly
  one real simple zero when Π_m⊗χ is self-dual; that zero is excluded for f⊗χ (non-dihedral cuspidal GL₂) by Hoffstein–
  Ramakrishnan (IMRN 1995), and for Sym²f⊗χ by Goldfeld–Hoffstein–Lieman (Ann. of Math. 140 (1994), appendix; χ = 1) and Banks
  (Duke Math. J. 87 (1997)); see [BTZ] §2.2. So for m ∈ {1, 2} there is an effective c₁ > 0 such that for every primitive
  Dirichlet character χ, L(s, Π_m ⊗ χ) ≠ 0 for Re s ≥ 1 − c₁/log(C(Π_m⊗χ)(|Im s| + 3)). The same holds for every m ≥ 1
  by Thorner (IMRN 2025, arXiv 2404.06482, Theorem 1.1 and Corollary 1.2 with F = Q, n = 0; f is regular algebraic and
  non-dihedral): L(σ+it, Sym^m f ⊗ χ) ≠ 0 for σ ≥ 1 − c/(m²[m^{y}log(mC(f)) + log(C(χ)(3+|t|))]), c absolute, y ∈ {0, 2}.
  Pinned 2026-10-21 in the TeX of v2 ([ChebST FIX.md](../ChebST/FIX.md)): Theorem 1.1 (l. 97–118), Corollary 1.2 and the Remark after
  it (n = 0: C(π′) → 1, y_{π′} → 0; l. 123–132), and, in the proof, eq. (5.1) (l. 600–604): L(σ, Sym^m(π) ⊗ ξ) ≠ 0 for
  σ ≥ 1 − c/(m^{2+y}log(mC(π)) + m²log C(ξ)), every Hecke character ξ (ξ = χ|·|^{it}, C(ξ) ≤ (4/3)C(χ)(3+|t|)); it covers m = 1, 2
  as well. Thorner, Res. Math. Sci. 2021 (arXiv 2002.10450) is a different paper (its Cor. ZFR allows one real simple zero); not used.
  (J5) concerns only Dirichlet twists, i.e. the genus characters; the twists by χ ∈ Ĝ_c are treated in §5.6.
- (J6) Explicit formula ([R] (I6), Iwaniec–Kowalski Thm 5.11): for an entire L-function L of degree 2m+2 as in (R2) below,
  w ∈ C_c^∞([x, x+H]), x ≥ 2: Σ_n Λ_L(n)w(n) = −Σ_ρ W(ρ) + R, W(s) := ∫w(u)u^{s−1}du, |R| ≤ C(H/x + x^{−1/2})log(q x).

## 1. The twisted family

**(R1) Dirichlet series.** For m ≥ 0 and χ ∈ Ĝ,
  L_m(s,χ) = Σ_𝔞 λ_m(𝔞)χ(𝔞)N𝔞^{−s} = Π_𝔭 Π_{i=1}^{r_𝔭}(1 − α_{m,i}(𝔭)χ(𝔭)N𝔭^{−s})^{−1}  (Re s > 1),
with r_𝔭 ≤ m+1, |α_{m,i}(𝔭)| ≤ 1 and (α_{m,i}(𝔭))_i independent of χ; the Euler factors at 𝔭 and 𝔭̄ coincide and have real
coefficients. Hence λ_m(𝔞) is real, λ_m(𝔞̄) = λ_m(𝔞), |λ_m(𝔞)| ≤ d_{m+1}(𝔞), and
−(L_m′/L_m)(s,χ) = Σ_𝔞 Λ_K(𝔞)b_m(𝔞)χ(𝔞)N𝔞^{−s}, b_m(𝔭^k) := Σ_i α_{m,i}(𝔭)^k, |b_m| ≤ m+1. For a degree-one 𝔭 with
p = N𝔭 ∤ N: (α_{m,i}(𝔭))_i = (e^{i(m−2j)θ_p})_{j=0}^{m}, so b_m(𝔭) = U_m(λ_f(p)/2) (b_1(𝔭) = λ_f(p), b_2(𝔭) = λ_f(p)² − 1).
λ_0(𝔞) = 1, and L_0(s,χ) = L(s,χ).

*Proof.* Local, by (J2). Let τ be the Weil–Deligne representation of Π_{m,p}. Then π_{χ,p} ↔ ⊕_{w|p}Ind_{W_{K_w}}^{W_{Q_p}}χ_w
and L(s, τ ⊗ Ind χ_w) = L(s, τ|_{W_{K_w}} ⊗ χ_w). As χ_w is unramified, L(s, τ|_{W_{K_w}} ⊗ χ_w) = P_w(χ(𝔭_w)N𝔭_w^{−s})^{−1},
where P_w(X)^{−1} is the L-polynomial of τ|_{W_{K_w}}, independent of χ, of degree dim(τ|_{W_{K_w}})^{I_{K_w}} ≤ m+1, with roots
of modulus ≤ 1 (τ is pure of weight 0: Π_{m,p} is tempered). Conjugate places give conjugate restrictions, hence equal
factors; self-duality of τ gives real coefficients. For split p ∤ N, τ is unramified with Frobenius eigenvalues e^{i(m−2j)θ_p}.
|λ_m(𝔭^k)| = |h_k(α)| ≤ binom(k+m, m) = d_{m+1}(𝔭^k). □
(Concretely: at p | D, p ∤ N, π_{χ,p} = ν ⊞ νη_p with ν unramified, ν(p) = χ(𝔭), η_p the local component of χ_D; at inert p,
π_{χ,p} is unramified with Satake parameters ±1, since χ((p)) = 1.)

**(R1′) The coefficients explicitly.** λ_m(𝔞) are the coefficients of the base change of Sym^m f to K, an L-function over K;
they are multiplicative on ideals of O_K and are *not* λ_{Sym^m f}(N𝔞) in general. For p ∤ N put α_j := e^{i(m−2j)θ_p}
(0 ≤ j ≤ m), so λ_{Sym^m f}(p^k) = h_k(α_0, …, α_m) (complete symmetric polynomial) and λ_{Sym^m f}(p) = U_m(λ_f(p)/2). By (R1):
- p = 𝔭𝔭̄ split: λ_m(𝔭^a𝔭̄^b) = λ_{Sym^m f}(p^a)λ_{Sym^m f}(p^b);
- p = 𝔭² ramified (p | D): λ_m(𝔭^k) = λ_{Sym^m f}(p^k);
- p inert, 𝔭 = (p), N𝔭 = p²: λ_m(𝔭^k) = h_k(α_0², …, α_m²) = Σ_{i=0}^{2k}(−1)^iλ_{Sym^m f}(p^i)λ_{Sym^m f}(p^{2k−i})
  (the local factor is Π_j(1 − α_jp^{−s})^{−1}(1 + α_jp^{−s})^{−1} = Π_j(1 − α_j²N𝔭^{−s})^{−1});
- p | N: λ_m(𝔭^k) = h_k(α_{m,1}(𝔭), …, α_{m,r_𝔭}(𝔭)), r_𝔭 ≤ m+1, |α_{m,i}(𝔭)| ≤ 1 (proof of (R1)).
For m = 1: λ_1(𝔭𝔭̄) = λ_f(p)² but λ_f(N(𝔭𝔭̄)) = λ_f(p²) = λ_f(p)² − 1; at inert p, λ_1((p)) = λ_f(p)² − 2 but λ_f(N(p)) = λ_f(p)² − 1.
What holds is the identity of Dirichlet series over Q, Σ_{N𝔞=n}λ_m(𝔞)χ(𝔞) = a_{Π_m×π_χ}(n), and for m = 1 Rankin's identity
  L_1(s,χ) = L(2s, χ_D)·Σ_𝔞 λ_f(N𝔞)χ(𝔞)N𝔞^{−s}   (up to Euler factors at p | N),
checked prime by prime: at split p both sides equal (1 − X²)^{−1}Σ_{a,b}λ_f(p^{a+b})z^a z̄^b X^{a+b} with X = p^{−s}, z = χ(𝔭); at inert p,
(1 + X²)^{−1}Σ_aλ_f(p^{2a})X^{2a} = Π_j(1 − α_j²X²)^{−1}; at p | D, p ∤ N, both are Σ_kλ_f(p^k)χ(𝔭)^kX^k. So the brief's coefficients
λ_f(N𝔞)χ(𝔞) are those of L_1(s,χ)/L(2s,χ_D); for m ≥ 2 there is no such identity with λ_{Sym^m f}(N𝔞).
*Bounds* (Deligne; temperedness at p | N). In every case λ_m(𝔭^k) = h_k(β_1, …, β_r) with r ≤ m+1 and |β_i| ≤ 1, hence
|λ_m(𝔭^k)| ≤ binom(k+r−1, r−1) ≤ binom(k+m, m) = d_{m+1}(𝔭^k), and |λ_m(𝔞)| ≤ d_{m+1}(𝔞) for every ideal 𝔞 (d_{m+1} is the divisor
function of O_K; at inert p it is d_{m+1}((p)^k) = binom(k+m, m), not d_{m+1}(p^{2k})). Likewise |μ_m(𝔭^k)| = |e_k(β)| ≤
binom(m+1, k) ≤ d_{m+1}(𝔭^k) (§3), |b_m(𝔭^k)| ≤ m+1, and |λ_f(N𝔞)| ≤ d(N𝔞) ≤ d(𝔞)² ≤ d_4(𝔞). These bounds are the only properties of
the coefficients used in Lemmas A, B and Theorem C.

**(R2) Analytic properties.** Let m ≥ 1. For every χ ∈ Ĝ, L_m(s,χ) is entire of order 1, of polynomial growth in vertical
strips, and Λ_m(s,χ) := q_m^{s/2}γ_m(s)L_m(s,χ) = ε_m(χ)Λ_m(1−s,χ), with |ε_m(χ)| = 1, γ_m(s) = Π_{j=1}^{2m+2}Γ_R(s+μ_j),
Re μ_j ≥ 0, and q_m = q_m(D), γ_m independent of χ, C^{−1}|D|^{m+1} ≤ q_m(D) ≤ C|D|^{m+1}.

*Proof.* χ ∈ Ĝ_c: (J1) with π = π_χ cuspidal and Π_m ≇ π̃_χ|·|^{iu} (different GL for m ≥ 2; for m = 1, f is not dihedral and has
weight k ≠ 1). χ ∈ Ĝ₂: a product of two twists of the cuspidal Π_m. Self-duality: Π̃_m ≅ Π_m, π̃_χ ≅ π_χ. Conductor exponents,
using (R1)'s local description: p ∤ ND: 0; p | D, p ∤ N: (m+1)a(η_p) = (m+1)v_p(D); p | N, p ∤ D: 2a(Π_{m,p}); p | (N,D):
a(Π_{m,p}⊗ν) + a(Π_{m,p}⊗νη_p) = a(Π_{m,p}) + a(Π_{m,p}⊗η_p), which does not depend on the unramified ν; all bounded in terms of
f and m. At ∞, π_{χ,∞} ↔ Ind_{W_C}^{W_R}1 for every χ (θ_χ has weight 1 for χ ∈ Ĝ_c; sgn ⊞ 1 for χ ∈ Ĝ₂), so γ_m does not
depend on χ. □

**(R3) Approximate functional equation** (Iwaniec–Kowalski Thm 5.3, Prop 5.4). For m ≥ 1, s = ½ + iv, χ ∈ Ĝ:
L_m(s,χ) = S₁(χ) + ε_m(χ,s)S₂(χ), |ε_m(χ,s)| = 1, S₁ = Σ_𝔞 λ_m(𝔞)χ(𝔞)N𝔞^{−s}V_s(N𝔞/√q_m),
S₂ = Σ_𝔞 λ_m(𝔞)χ(𝔞)N𝔞^{s−1}V_{1−s}(N𝔞/√q_m), |V_s(y)| ≤ C_B(1 + y/(|v|+3)^{m+1})^{−B}.
(The dual coefficients are the conjugates Σ_{N𝔞=n}λ_m(𝔞)χ̄(𝔞) = Σ_{N𝔞=n}λ_m(𝔞̄)χ(𝔞̄); reindex 𝔞 ↦ 𝔞̄.)

**(R4) Zero counts** (Iwaniec–Kowalski Prop 5.7). N_m(0;t,χ) ≤ C log(|D|T); hence N_D^{(m)}(σ;t) ≤ C h log(|D|T)
≤ C|D|^{1/2}log³(3|D|T) for every σ ([C] (F3)).

## 2. Large sieve and second moments

**Lemma A (weighted large sieve, twisted).** For m ≥ 0, 0 < ε ≤ ½, N ≥ 1 and complex b_𝔞 supported on N𝔞 ≤ N:
  Σ_{χ∈Ĝ}|Σ_𝔞 λ_m(𝔞)b_𝔞χ(𝔞)|² ≤ (30/ε)log²(3|D|)(N + |D|^{1/2})Σ_𝔞 m(𝔞)^{1+2ε}d_{m+1}(𝔞)²|b_𝔞|².

*Proof.* [C] Lemma 2 holds for arbitrary complex c_𝔞 supported on N𝔞 ≤ N (here f = 1, so on all ideals); take
c_𝔞 = λ_m(𝔞)b_𝔞 and use |λ_m(𝔞)| ≤ d_{m+1}(𝔞) from (R1′). Nothing else about the coefficients enters. □

Remarks. (i) The weight m(𝔞)^{1+2ε} cannot be dropped ([C] Remark after Lemma 2). (ii) For coefficients bounded by divisor
functions, [C] Lemma 3 makes the weight cost N^{2ε}: the h ≍ |D|^{1/2} characters behave like characters of a modulus |D|^{1/2},
for the twisted coefficients as well. (iii) The brief's coefficients λ_f(N𝔞)χ(𝔞) differ from λ_1(𝔞)χ(𝔞) ((R1′)): up to Euler factors at p | N,
L_1(s,χ) = L(2s,χ_D)Σ_𝔞 λ_f(N𝔞)χ(𝔞)N𝔞^{−s} (Rankin's identity, (R1′)). Lemma A also holds with λ_f(N𝔞) in place of λ_m(𝔞) and
d(N𝔞) in place of d_{m+1}(𝔞) ([C] Lemma 2 again), and [C] Lemma 3 applies to it because d(N𝔞) ≤ d(𝔞)² ≤ d_{16}(𝔞) (check at
split, inert and ramified prime powers). (iv) General-family large
sieves (Humphries–Thorner, Pascadi–Thorner) have threshold Q²|𝒮| ≍ |D|^{5/2} for this family, i.e. |D|² times worse.

**Lemma B (mollified second moment).** Let m ≥ 1, ε > 0, V ≥ 1, v ∈ R, and |μ(𝔡)| ≤ d_{m+1}(𝔡); M_V(s,χ) :=
Σ_{N𝔡≤V}μ(𝔡)χ(𝔡)N𝔡^{−s}. Then, with C(ε) depending only on ε, f, m (uniform in V ≥ 1 and v ∈ R),
  Σ_{χ∈Ĝ}|L_m(½+iv,χ)M_V(½+iv,χ)|² ≤ C(ε)|D|^{(m+1)/2+ε}V^{1+ε}(|v|+3)^{m+1+ε},
and Σ_{χ∈Ĝ}|M_V(½+iv,χ)|² ≤ C(ε)(V + |D|^{1/2})(V|D|)^{ε}.

*Proof.* Put s = ½ + iv, L₀ := √q_m(|v|+3)^{m+1} and 𝔏 := 3|D|V(|v|+3) ≥ 9. By (R3), |L_mM_V|² ≤ 2|S₁M_V|² + 2|S₂M_V|², and
S_jM_V = Σ_𝔟 c_{j,𝔟}χ(𝔟) with c_{j,𝔟} independent of χ and, by (R1′) and the bound on V_s in (R3),
  |c_{j,𝔟}| ≤ C_B N𝔟^{−1/2} Σ_{𝔞𝔡=𝔟, N𝔡≤V} d_{m+1}(𝔞)d_{m+1}(𝔡)(1 + N𝔞/L₀)^{−B} ≤ C_B d_{2m+2}(𝔟)N𝔟^{−1/2}.
Cut at P := VL₀𝔏^ε (the old proof cut at VL₀|D|^ε and claimed a tail O(|D|^{−8}); that tail is not uniform in V and v, since it
carries a factor ≍ (VL₀)^{1/2}).
*Tail.* If N𝔟 > P, every term of c_{j,𝔟} has N𝔞 > P/V = L₀𝔏^ε. The number of ideals of norm n is ≤ d(n), d_{m+1}(𝔞) ≤
d_{m+1}(N𝔞)², and d(n)d_{m+1}(n)² ≤ C(m)n^{1/2}. So for B ≥ 2, with Y := L₀𝔏^ε,
  Σ_{N𝔟>P}|c_{j,𝔟}| ≤ C_B Σ_{N𝔡≤V}d_{m+1}(𝔡)N𝔡^{−1/2} · L₀^B Σ_{n>Y}C(m)n^{−B} ≤ C(m,B)·V·L₀(L₀/Y)^{B−1} ≤ C(m,B)𝔏^{m+3−ε(B−1)},
using V ≤ 𝔏 and L₀ ≤ C|D|^{(m+1)/2}(|v|+3)^{m+1} ≤ C𝔏^{m+1} (R2). With B := 1 + (m+13)/ε the tail is ≤ C(ε)𝔏^{−10} for every
χ, uniformly in V and v (a dyadic decomposition of N𝔞 > Y gives the same).
*Main part.* Lemma A (m = 0) and [C] Lemma 3 (k = 2m+2) on N𝔟 ≤ P:
  Σ_χ|Σ_{N𝔟≤P}c_{j,𝔟}χ(𝔟)|² ≤ (30/ε)log²(3|D|)(P + |D|^{1/2})C_B²Σ_{N𝔟≤P}m(𝔟)^{1+2ε}d_{2m+2}(𝔟)²N𝔟^{−1} ≤ C(ε)log²(3|D|)P^{1+2ε},
as P ≥ L₀ ≥ √q_m ≥ c|D|^{1/2}. With the tail (h ≤ |D|): Σ_χ|S_jM_V|² ≤ C(ε)P^{1+3ε} + C(ε)|D|𝔏^{−20} ≤ C(ε)P^{1+3ε}. Finally
P ≤ CV|D|^{(m+1)/2}(|v|+3)^{m+1}𝔏^ε and 𝔏^{O(ε)} = (|D|V(|v|+3))^{O(ε)}; rename ε. The second claim is Lemma A (m = 0) with
b_𝔡 = μ(𝔡)N𝔡^{−1/2−iv} (N𝔡 ≤ V) and [C] Lemma 3; there is no tail. □

*Uses (re-checked).* Only V = 1 (M_V ≡ 1) and the second claim are used: Theorem C class II integrates
Σ_χ|L_m(½+iv,χ)|² ≤ C(ε)|D|^{(m+1)/2+ε}(|v|+3)^{m+1+ε} against e^{−|v−t|} over all v ∈ R, which needs exactly this uniformity in v
(∫e^{−|v−t|}(|v|+3)^{(m+1+ε)/2}dv ≤ C(m)T^{(m+1+ε)/2}); Proposition F uses it on |v| ≤ T (I₂); Theorem D⁺ uses it for each
archimedean type. The statement did not change, so every bound downstream is unchanged.

Compare: convexity gives |L_m(½+iv,χ)|² ≪ |D|^{(m+1)/2+ε} for each χ, i.e. h·|D|^{(m+1)/2} for the sum; Lemma B saves h. Lindelöf
on average would be h|D|^ε (§9).

## 3. Per-D zero density

**Theorem C.** Let m ≥ 1, ε > 0, D ∈ 𝔇, t ∈ R, T = |t|+3 and 0.6 ≤ σ ≤ 1. Then, for the zeros in the unit window
|γ − t| ≤ 1,
  N_D^{(m)}(σ;t) ≤ C(ε)(|D|^{m+2}T^{2m+2})^{(1−σ)/(3−2σ)+ε}.
**(C′) All heights.** For T ≥ 3: Σ_{χ∈Ĝ}N_m(σ,T;χ) ≤ C(ε)·T·(|D|^{m+2}T^{2m+2})^{(1−σ)/(3−2σ)+ε}. Cover |γ| ≤ T by the ≤ 2T+1
windows |γ − k| ≤ 1 (k ∈ Z, |k| ≤ T), each with |k| + 3 ≤ 2T. The factor T is the number of windows: Theorem C is a unit-window
estimate and does not absorb it. Every use below says which form it uses.
For m = 1 (f ⊗ θ_χ, degree 4): ≤ C(ε)(|D|³T⁴)^{(1−σ)/(3−2σ)+ε} per window; this is below the trivial per-window bound |D|^{1/2}
exactly for σ > 3/4 (general m: σ > (2m+1)/(2m+2)); for |γ| ≤ T both bounds carry the factor T.

*Proof.* [R2] Proposition 1 with the following substitutions; we write it out. Assume |D| ≥ D₀(ε) (else use (R4)).
Put ε₀ := ε/10, 𝔔 := |D|^{m+2}T^{2m+2}, V := |D|^{1/2}, Y := 𝔔^{1/(2(3−2σ))}, Z := 2Y log²Y. Then Y ≥ |D|^{(m+2)/(2(3−2σ))}
≥ V, because m+2 ≥ 3 ≥ 3−2σ. Let μ_m(𝔡) be the coefficients of Π_𝔭Π_i(1 − α_{m,i}(𝔭)N𝔭^{−s}), so 1/L_m(s,χ) =
Σ_𝔡μ_m(𝔡)χ(𝔡)N𝔡^{−s} and |μ_m(𝔡)| ≤ d_{m+1}(𝔡) (|e_k(α)| ≤ binom(m+1,k) ≤ binom(m+k,k)). Let M_V be as in Lemma B with
μ = μ_m, and a_𝔞 := Σ_{𝔡|𝔞, N𝔡≤V}μ_m(𝔡)λ_m(𝔞𝔡^{−1}): a_{O_K} = 1, a_𝔞 = 0 for 1 < N𝔞 ≤ V, |a_𝔞| ≤ d_{2m+2}(𝔞), and
L_m(s,χ)M_V(s,χ) = Σ_𝔞a_𝔞χ(𝔞)N𝔞^{−s} (Re s > 1). For each χ with a zero in [σ,1] × [t−1,t+1] fix one such zero ρ_χ; 𝒳 is
the set of these χ, and N_D^{(m)}(σ;t) ≤ C log(|D|T)|𝒳| by (R4).
*Detector.* [C] Lemma 5 holds verbatim for L_m in place of L(·,χ): L_m(·,χ) is entire of polynomial growth in vertical strips
(R2), so the contour shift to Re w = ½ − β ∈ [−½, −1/10] crosses only w = 0, where the residue L_m(ρ)M_V(ρ) vanishes. So
𝒳 = 𝒳_I ∪ 𝒳_II with (I) |Σ_{N𝔞>V}a_𝔞χ(𝔞)N𝔞^{−ρ}e^{−N𝔞/Y}| ≥ ¼ or (II) ∫|L_mM_V(½+i(γ_χ+u),χ)|e^{−|u|}du ≥ c₁Y^{σ−½}.
*Class II.* Sum (II) over 𝒳_II, substitute v = γ_χ + u (e^{−|v−γ_χ|} ≤ e^{1−|v−t|}) and apply Cauchy–Schwarz at each v over
all of Ĝ:
c₁Y^{σ−½}|𝒳_II| ≤ e∫e^{−|v−t|}(Σ_χ|L_m(½+iv,χ)|²)^{1/2}(Σ_χ|M_V(½+iv,χ)|²)^{1/2}dv ≤ C(ε)(|D|^{(m+1)/2}T^{m+1}·|D|^{1/2})^{1/2}
|D|^{2ε₀}T^{ε₀}, by Lemma B (V = 1 for L_m; the second claim for M_V). So |𝒳_II| ≤ C(ε)Y^{½−σ}𝔔^{1/4}(|D|T)^{2ε₀}.
*Class I.* As in [C] Theorem 1 / [R2] Proposition 1: the tail N𝔞 > Z contributes ≤ 1/16 (|a_𝔞| ≤ d_{2m+2}(𝔞) ≤ d_{2m+2}(N𝔞)²);
dyadic blocks j with 2^{j+1} > V and 2^{j−1} < Z (J ≤ 2log Z of them); removal of β by Mellin inversion; Sobolev on unit
intervals; Lemma A (m = 0) with M := 2^{j+1} ≥ V = |D|^{1/2}, so the constant is ≤ (60/ε₀)log²(3|D|)M, and [C] Lemma 3 with
|b_𝔞| ≤ d_{2m+2}(𝔞)2^{−(j−1)σ}: Σ_χ sup_{t′∈I_τ}|P_j(χ,t′)|² ≤ C(ε)log⁴(3|D|Y)M^{2−2σ+2ε₀}. Hence |𝒳_I| ≤ C(ε)(|D|T)^{ε₀}
Z^{2−2σ+2ε₀}.
*Conclusion.* Y^{2−2σ} = 𝔔^{(1−σ)/(3−2σ)} and Y^{½−σ}𝔔^{1/4} = 𝔔^{(1−2σ)/(4(3−2σ)) + 1/4} = 𝔔^{(1−σ)/(3−2σ)}. Rename ε. □

*Comparison.* For 𝒮 = {π_χ : χ ∈ Ĝ_c} (n = 2, Q ≍ |D|, |𝒮| ≍ |D|^{1/2}) and π₀ = Π_m (n₀ = m+1), Pascadi–Thorner
(arXiv 2508.14888, Theorem 1.3, second line) give Σ_{π∈𝒮}N_{π×π₀}(σ,T) ≪ (Q^{2n+n₀}|𝒮|⁴T^{2(m+3)+4})^{(1−σ)/(3−2σ)+o(1)} =
(|D|^{m+7}T^{2m+10})^{(1−σ)/(3−2σ)+o(1)}; Humphries–Thorner Thm 1.1 give (|D|^{2+6.15·max(4,2m+2)})^{1−σ+ε}; the log-free bounds
(BTZ Thm 1.2, PT Thm 1.5) have exponents 10^7(2m+2)^4(1−σ) and 543·max{8,(m+1)³}(1−σ). PT count all |γ| ≤ T, so the comparison is with (C′):
T(|D|^{m+2}T^{2m+2})^{e} against (|D|^{m+7}T^{2m+10})^{e}, e = (1−σ)/(3−2σ). In the D-aspect (T fixed, or T ≤ |D|^{o(1)}) (C′) lowers
the |D|-exponent from (m+7)e to (m+2)e for every σ ∈ [0.6, 1), and moves the start of the nontrivial range from
σ > (2m+11)/(2m+12) (13/14 for f) to σ > (2m+1)/(2m+2) (3/4 for f). In T, (C′) has T^{1+(2m+2)e} against T^{(2m+10)e}: better
for σ < 5/6, worse for σ > 5/6. The gain has two sources: the per-D sieve threshold |D|^{1/2}
instead of Q²|𝒮| = |D|^{5/2}, and Lemma B instead of convexity. Only the D-aspect is claimed.
(Kowalski–Michel's 2002 family theorem needs automorphy of Π_m ⊠ π_χ and has exponent c(1−σ)/(2σ−1) with an unspecified c
larger than the family-size exponent; not computed here.)

## 4. Density over all discriminants

**Proposition F.** Let m ≥ 1, X, T ≥ 3, 0 ≤ σ ≤ 1, ε > 0. Then
  Σ_{D∈𝔇(X)} N_D^{c,m}(σ,T) ≤ C(ε)(X^{m+10}T^{2m+10})^{(1−σ)/(3−2σ)+ε}.

*Proof.* (J3) with 𝒮 := {π_χ : D ∈ 𝔇(X), χ ∈ Ĝ_c(D)} and π′ := Π_m. π_χ = π_{χ′} iff χ′ ∈ {χ, χ̄} (π_χ determines its
conductor |D| and its central character χ_D, hence K), so the left side is 2Σ_{π∈𝒮}N_{π×Π_m}(σ,T). GRC holds for π_χ ([R]
Prop 1); C(π_χ) ≤ 12|D|, so Q := 12X; |𝒮| ≤ Σ_{|D|≤X}h(D) ≤ CX^{3/2}log X; Π_m ≇ π̃ for π ∈ 𝒮 (R2). Hence
V ≤ C(ε)X^{7/2+2ε}T^{3+ε}. For I₂: |L^{ur}| ≤ 4^{(m+1)ω(DN)}|L| on Re s = ½ (the removed Euler factors have ≤ 2m+2 roots of
modulus ≤ 1), so Lemma B (V = 1, M ≡ 1; uniform in v, integrated over |v| ≤ T) gives I₂ ≤ C(ε)Σ_{D∈𝔇(X)}|D|^{(m+1)/2+ε}T^{m+2+ε} ≤ C(ε)X^{(m+3)/2+2ε}T^{m+2+ε}.
(HT's convexity bound is |𝒮|·X^{(m+1)/2}T^{m+2} ≍ X^{(m+4)/2}.) Take Y := (V·I₂)^{1/(3−2σ)}; then log Y ≍ log V and
Y^{½−σ}(VI₂)^{1/2} = Y^{2−2σ} ≤ C(ε)(X^{m+10}T^{2m+10})^{(1−σ)/(3−2σ)+O(ε)}. Rename ε. □

Comparison: Pascadi–Thorner Theorem 1.3 with Q ≍ X, |𝒮| ≍ X^{3/2}: (X^{4+(m+1)+6}T^{…})^{(1−σ)/(3−2σ)} = X^{(m+11)(1−σ)/(3−2σ)}.
The saving X^{(1−σ)/(3−2σ)} is exactly Lemma B against convexity. For m = 0 this is [R] Proposition 1 (X^{10}).

## 5. Primes in ideal classes with a Hecke-eigenvalue condition

**Theorem D.** Let P(x) = c₀ + c₁x + c₂(x² − 1) with c₀ > 0, and let S ⊂ [−2,2] satisfy P ≤ 1_S on [−2,2]. Put A₀ := 3/2 if
c₁ = c₂ = 0; A₀ := 9/4 if c₂ = 0 ≠ c₁; A₀ := 3 if c₂ ≠ 0. For every A > A₀ there are effective c, κ, δ > 0 and X₀
(depending on f, A, P) such that for X ≥ X₀ all D ∈ 𝔇 with |D| ≤ X outside a set of size ≤ X^{1−c} satisfy: for every
C ∈ G, every x ≥ |D|^A and every H ∈ [x|D|^{−κ}, x],
  Σ_{𝔭∈C, N𝔭 = p, p ∤ ND, x < p ≤ x+H, λ_f(p) ∈ S} log p ≥ c₀ ν_D H/(32 h(D)),        (5.1)
where ν_D is as in [R] Theorem 2 with this δ: ν_D = 1 unless the exceptional character χ₁ of [R] (I4) for Q₀ = |D|^δ exists
and equals χ_d for some d ∈ 𝔉(D), and always ν_D ≥ c′|D|^{−δ/2}/log|D|.

**Corollary D1.** Let E/Q have no CM and f = f_E. For every A > 3, all but O(X^{1−c}) D ∈ 𝔇 with |D| ≤ X satisfy: every
primitive positive definite binary quadratic form of discriminant D represents primes p, p′ ∤ ND with p, p′ ≤ 2|D|^A and
a_E(p) > 0 > a_E(p′); indeed ≫ ν_D H/(h log x) such primes of each sign in every (x, x+H] with x ≥ |D|^A, H ≥ x|D|^{−κ}.
The same holds (A > 3) for a_E(p) ≠ 0, |a_E(p)| < β√p (β > 1), |a_E(p)| > β√p (β < 1), a_E(p) > α√p (α < ½), and (A > 9/4)
for a_E(p) > −α√p and a_E(p) < α√p with α > 0. In particular, for all but O(N^{1−c}) squarefree n ≤ N with n ≡ 1, 2 (mod 4):
a prime p = a² + nb² ≤ 2(4n)^A with a_E(p) > 0 (and one with a_E(p) < 0).

*Proof of D1 from D.* Minorants, all checked in `params.py` and `params2.py` (P ≤ 1_S on [−2,2], Sato–Tate mean = c₀;
c₁, c₂, η₀ and the constant of (5.1) in the table after the proof of Theorem D):
S = (0,2]: P = x(x+2)/8 = 1/8 + x/4 + (x²−1)/8; S = [−2,0): x(x−2)/8; S = [−2,2]∖{0}: x²/4 = ¼ + (x²−1)/4;
S = (α,2], α < 0: (x−α)/(2−α), c₀ = |α|/(2+|α|), c₂ = 0; S = (α,2], 0 ≤ α < ½: (x+2)(x−α)/(4(2−α)), c₀ = (1−2α)/(4(2−α));
S = (−β,β), β > 1: 1 − x²/β², c₀ = 1 − 1/β²; S = {|x| > β}, β < 1: (x²−β²)/(4−β²), c₀ = (1−β²)/(4−β²); [−2,α) by x ↦ −x.
A degree-one 𝔭 ∈ C = [𝔞_Q]^{−1} with N𝔭 = p ∤ D gives a representation of p by Q ([C] (F1), (F2)); each p lies under ≤ 2
such 𝔭. For x² + ny² take D = −4n and the principal class ([R] Corollary 4). □

Remark (degree 1 cannot see the sign). If P = c₀ + c₁x ≤ 1_{(0,2]} then P ≤ 0 on [−2,0], so c₀ = P(0) ≤ 0. A minorant with
c₁ ≠ 0, c₂ = 0 exists iff S contains (α,2] or [−2,α′) with α < 0 < α′. So the sign needs Sym²f, and A > 3.

**Proposition T.** Let m ∈ {1,2}, A₁ := 9/4, A₂ := 3, A > A_m, η₀ ∈ (0,1]. There are effective δ_m = δ_m(f,A,η₀) > 0, c > 0
such that for every δ ∈ (0, δ_m], with κ := δ/4, and X ≥ X₀(f,A,η₀,δ), there is E^{(m)}(X) ⊂ 𝔇(X) with #E^{(m)}(X) ≤ X^{1−δ/8}
such that for D ∈ 𝔇(X) ∖ E^{(m)}(X), C ∈ G, x ≥ |D|^A and H ∈ [x|D|^{−κ}, x]:
  |h S_C^{(m)}| ≤ η₀ν_D H,  S_C^{(m)} := Σ_{𝔞∈C}Λ_K(𝔞)b_m(𝔞)w(N𝔞),
with w(u) := w₀((u−x)/H) as in [R] §2.3 (w₀ ∈ C_c^∞((0,1)), 0 ≤ w₀ ≤ 1, w₀ = 1 on [¼,¾]) and ν_D defined with δ.

*Proof of Theorem D.* m = 0: [R2] Theorem 2′ (A ≥ 1.55 numbers, then any A′ ≥ A) holds with every δ ≤ 1/2840 and κ = δ/4
(every condition in [R] §2.1–2.5 and [R2] Thm 2′ is an upper bound on δ, κ, ε₀, ε₁), and its proof gives h S_C^{(0)} ≥ ν_D H/8
outside its exceptional set ([R] §2.5). Take η₀ := c₀/(16(|c₁|+|c₂|)) (if c₁ = c₂ = 0 only m = 0 is needed), δ := min{1/2840,
δ_1, δ_2} over the m that occur, κ := δ/4, and remove the union over dyadic X′ ≤ X of the three exceptional sets
(≤ X^{1−c}, c := δ/9). Then P_𝔞 := c₀ + c₁b₁(𝔞) + c₂b₂(𝔞) satisfies
h Σ_{𝔞∈C}Λ_K(𝔞)P_𝔞w(N𝔞) = c₀hS_C^{(0)} + c₁hS_C^{(1)} + c₂hS_C^{(2)} ≥ c₀ν_DH/8 − (|c₁|+|c₂|)η₀ν_DH = c₀ν_DH/16.
For degree-one 𝔭 with p = N𝔭 ∤ ND, P_𝔭 = P(λ_f(p)) by (R1). All other 𝔞 with N𝔞 ∈ (x, x+H] (prime powers, inert and
ramified primes, primes above p | N) have |P_𝔞| ≤ |c₀| + 2|c₁| + 3|c₂| and number ≤ C(H/√x + ω(ND) + log x)log x; since
x ≥ |D|² (A > 9/4 when c₁ or c₂ ≠ 0; if c₁ = c₂ = 0 the claim is [R2] Thm 2′), h ≤ |D|^{1/2+o(1)} and ν_D ≥ |D|^{−δ}, their
contribution is ≤ c₀ν_DH/(32h) for |D| ≥ D₀. So
Σ_{p∤ND} P(λ_f(p))log p·w(p) ≥ c₀ν_DH/(32h) over degree-one 𝔭 ∈ C; as P(λ) ≤ 1 for λ ∈ S, P(λ) ≤ 0 otherwise, and
0 ≤ w ≤ 1_{(x,x+H]}, (5.1) follows. □

**Positivity: what c₀ > 0 does and does not give.** c₀ > 0 is the Sato–Tate mean of P (∫U_m dμ_ST = 0 for m ≥ 1), but
Theorem D never uses equidistribution. Positivity of the minorant sum needs three inputs for the same D, C, x, H:
(i) the m = 0 lower bound hS_C^{(0)} ≥ ν_DH/8 ([R] §2.5 with [R2] Theorem 2′); a Landau–Siegel zero of a genus factor χ_{d₀},
d₀ | D, shrinks this main term to ν_DH/8, with ν_D as small as c′|D|^{−δ/2}/log|D|;
(ii) the m = 1, 2 sums as errors relative to that shrunken main term: |hS_C^{(m)}| ≤ η₀ν_DH, η₀ := c₀/(16(|c₁| + |c₂|))
(Proposition T); this is why (e) of §5.4 needs the Landau–Siegel repulsion ν^♯ ≤ 2ν_D of (J4): without it (e) gives only
E_e ≤ η₀H/8, which is not ≤ η₀ν_DH/8 when ν_D is small;
(iii) the non-degree-one ideals (≤ c₀ν_DH/(32h) after division by h).
Then c₀ν_DH/8 − (|c₁| + |c₂|)η₀ν_DH = c₀ν_DH/16 and c₀ν_DH/(16h) − c₀ν_DH/(32h) = c₀ν_DH/(32h), which is (5.1). With c₀ > 0 but
without (i) or (ii) at this size nothing follows. The constants of Corollary D1 (`params2.py` re-checks P ≤ 1_S, the mean c₀, η₀
and the constant in (5.1); α, β as in Corollary D1):

| S | P | c₀ | c₁ | c₂ | η₀ | (5.1): Σ log p ≥ |
|---|---|---|---|---|---|---|
| λ > 0 | x(x+2)/8 | 1/8 | 1/4 | 1/8 | 1/48 | ν_DH/(256h) |
| λ < 0 | x(x−2)/8 | 1/8 | −1/4 | 1/8 | 1/48 | ν_DH/(256h) |
| λ ≠ 0 | x²/4 | 1/4 | 0 | 1/4 | 1/16 | ν_DH/(128h) |
| λ > α, α < 0 | (x−α)/(2−α) | \|α\|/(2+\|α\|) | 1/(2+\|α\|) | 0 | \|α\|/16 | \|α\|ν_DH/(32(2+\|α\|)h) |
| λ > α, 0 ≤ α < ½ | (x+2)(x−α)/(4(2−α)) | (1−2α)/(4(2−α)) | 1/4 | 1/(4(2−α)) | (1−2α)/(16(3−α)) | (1−2α)ν_DH/(128(2−α)h) |
| \|λ\| < β, 1 < β ≤ 2 | 1 − x²/β² | 1 − 1/β² | 0 | −1/β² | (β²−1)/16 | (1−β^{−2})ν_DH/(32h) |
| \|λ\| > β, β < 1 | (x²−β²)/(4−β²) | (1−β²)/(4−β²) | 0 | 1/(4−β²) | (1−β²)/16 | (1−β²)ν_DH/(32(4−β²)h) |

### 5.1 Parameters for Proposition T
The statement for A implies it for every A′ > A (same δ, κ, exceptional set), so assume A_m < A ≤ A_m + ¼.
- m = 2: σ₀ := 5/6, s₁ := (A−3)/6, ε′ := s₁/4. For σ ≤ 5/6: A(1−σ) − ½ ≥ s₁. φ₂(σ) := (12/(3−2σ) − A)(1−σ) has
  φ₂′ = A − 12/(3−2σ)² ≤ 3.25 − 27/4 < 0 on [5/6, 1] and φ₂(5/6) = (9−A)/6, so φ₂ ≤ 1 − 4ε′ on [σ₀, 1].
- m = 1: σ₁ := 3/2 − 3/(2A), τ := (A − 9/4)/42, σ₀ := σ₁ − τ (so 0.82 ≤ σ₀ ≤ 0.9), ε′ := (A − 9/4)/4,
  s₁ := (1−σ₀)(A − 3/(3−2σ₀)) ≥ (A − 9/4)/300. φ₁(σ) := (11/(3−2σ) − A)(1−σ): φ₁(σ₁) = 4 − 4A/3, |φ₁′| ≤ 14 on [0.6,1], so
  φ₁(σ₀) ≤ 13/4 − A = 1 − 4ε′, and φ₁′ = A − 11/(3−2σ)² < 0 on [σ₀,1] (3 − 2σ₀ ≤ 1.36); so φ₁ ≤ 1 − 4ε′ on [σ₀,1]. Per-D margins: for
  σ ∈ [0.6, ¾], A(1−σ) − ½ ≥ 1/16; for σ ∈ [¾, σ₀], (1−σ)(A − 3/(3−2σ)) ≥ s₁ (both factors positive, decreasing).
- Common: c_m := 10^7(m+1)^4 and B_m from (J4); c₁ from (J5), c₁′ := c₁/(m+3); σ_g := 1 − 1/(2c_m);
  δ_m := min{1/2840, s₁/10, ε′/10, A/(8c_m), δ_e}, δ_e := A/(1.3(c_m + c₁′^{−1}log(64B_m/η₀)));
  for δ ≤ δ_m: κ := δ/4, ε₀ := ε₁ := δ/100, σ_j := σ₀ + jε′/A (0 ≤ j ≤ J, J minimal with σ_J ≥ 1 − ε′/A),
  T_X := 2X^{κ+ε₁}, Q₀(X) := (X/2)^δ; for a given D: Q₀ := |D|^δ, T₁ := 2|D|^{κ+ε₁} ≤ T_X, 𝔉_s := {d ∈ 𝔉(D) : |d| ≤ Q₀},
  𝔉_l := 𝔉(D) ∖ 𝔉_s.
`params.py` checks all of these and the margins of §5.4 at A ∈ {3.0001, 3.01, 3.05, 3.1, 3.25} (m = 2) and A ∈ {2.2501, 2.26,
2.3, 2.4, 2.5} (m = 1), with δ the explicit part of the minimum (δ_e only makes δ smaller).

### 5.2 Exceptional sets (X a power of 2)
- E₁^{(m)}(X): D with N_D^{c,m}(σ_j, T_X) > X^{A(1−σ_j)−2ε′} for some j. By Proposition F (ε = ε₀) and Markov, using
  φ_m(σ_j) ≤ 1 − 4ε′ and (2m+10)(1−σ)/(3−2σ) ≤ 2 for σ ≥ 0.8:
  #E₁ ≤ Σ_j X^{2ε′−A(1−σ_j)}C(X^{m+10}T_X^{2m+10})^{(1−σ_j)/(3−2σ_j)+ε₀} ≤ CJX^{2ε′+φ_m(σ_j)+2(κ+ε₁)+14ε₀} ≤ CJX^{1−ε′}.
- E₂^{(m)}(X): D such that some d ∈ 𝔉(D) with |d| > Q₀(X) has a zero of L(s, Π_m⊗χ_d) with β > σ_g, |γ| ≤ T_X. By (J4)
  (ν^♯ = 1, Q = 2Z, T = T_X) at most B_m(2C(Π_m)ZT_X)^{1/2} primitive χ of conductor ≤ 2Z have such a zero; each d with
  |d| ∈ (Z, 2Z] divides ≤ 2X/Z elements of 𝔇(X). Summing over dyadic Z ∈ [Q₀(X), X]:
  #E₂ ≤ CXQ₀(X)^{−1/2}T_X^{1/2}log X ≤ X^{1−δ/4} (κ + ε₁ < δ/2).
- E^{(m)}(X) := E₁ ∪ E₂: size ≤ X^{1−δ/8} for X ≥ X₀ (δ ≤ ε′/10).

### 5.3 Decomposition
Fix D ∈ 𝔇(X) ∖ E^{(m)}(X) with |D| ≥ D₀, C, x, H. By orthogonality, (G) and (R1),
  h S_C^{(m)} = Σ_{χ∈Ĝ}χ̄(C)S_χ^{(m)},  S_χ^{(m)} := Σ_n Λ_{L_m(·,χ)}(n)w(n);
for χ = χ_{d,D/d} ∈ Ĝ₂, S_χ^{(m)} = ψ_w(Π_m⊗χ_d) + ψ_w(Π_m⊗χ_{D/d}), so the genus part is Σ_{d∈𝔉(D)}ε_d(C)ψ_w(Π_m⊗χ_d).
Apply (J6) to every term (m ≥ 1: no main terms). |W(ρ)| ≤ Hx^{β−1}, and |W(ρ)| ≤ C_kHx^{β−1}(2x/(H|ρ|))^k. Let E_a, …, E_f
bound the parts: (a) zeros of L_m(·,χ), χ ∈ Ĝ_c, β ≥ σ₀, |γ| ≤ T₁; (b) zeros of L_m(·,χ), all χ ∈ Ĝ, β < σ₀, |γ| ≤ T₁;
(c) all zeros with |γ| > T₁; (d) zeros of L(s, Π_m⊗χ_d), d ∈ 𝔉_l, |γ| ≤ T₁; (e) zeros of L(s, Π_m⊗χ_d), d ∈ 𝔉_s,
|γ| ≤ T₁; (f) remainders. Every zero lies in at least one part; overlaps only enlarge the bound.

### 5.4 Bounds (|D| ≥ D₀; x = |D|^A y, y ≥ 1)
- (a) (N_D^{c,m}(σ,T_X) and Proposition F count all |γ| ≤ T_X; no window count enters here.)
  D ∉ E₁ and x^{σ−1} ≤ |D|^{−A(1−σ)} ≤ 2^AX^{−A(1−σ)}: E_a ≤ HΣ_{j<J}N_D^{c,m}(σ_j,T_X)x^{σ_{j+1}−1} ≤ J2^AHX^{−ε′};
  N_D^{c,m}(σ_J,T_X) < 1 since A(1−σ_J) ≤ ε′.
- (b) Cover |γ| ≤ T₁ by the unit windows |γ − k| ≤ 1, k ∈ Z, |k| ≤ T₁: at most 2T₁ + 1 ≤ 3T₁ windows, each with |k| + 3 ≤ 2T₁.
  Per window, with N(σ) := N_D^{(m)}(σ;k): Σ_{β<σ₀}x^{β−1} ≤ x^{−½}N(0) + log x∫_{½}^{σ₀}x^{σ−1}N(σ)dσ, and
  x^{σ−1}log x ≤ |D|^{−A(1−σ)}(A log|D| + 10) for σ ≤ σ₀ ≤ 0.9. (R4): N(σ) ≤ |D|^{½+o(1)}(|k|+3)^{o(1)}.
  m = 2: each window ≤ |D|^{½−A/2+o(1)} + |D|^{½−A/6+o(1)} = |D|^{−s₁+o(1)}; times 3T₁ windows: E_b ≤ H|D|^{−s₁+κ+ε₁+o(1)}.
  m = 1: σ ∈ [½, 0.6]: |D|^{½−0.4A+o(1)} ≤ |D|^{−0.4}; [0.6, ¾]: |D|^{−1/16+o(1)}; [¾, σ₀]: Theorem C per window (ε = ε₀),
  |D|^{−(1−σ)(A−3/(3−2σ))+3ε₀}(|k|+3)^{2/3+4ε₀} ≤ |D|^{−s₁+3ε₀}(2T₁)^{2/3+4ε₀}. Summing over the ≤ 3T₁ windows (this factor T₁ is
  the window count of (C′), not part of Theorem C): E_b ≤ H|D|^{−s₁+3ε₀+(5/3+4ε₀)(κ+ε₁)+o(1)} + H|D|^{−1/16+κ+ε₁+o(1)}.
- (c) Each window holds ≤ C h log(|D|(|k|+3)) zeros (R4); with k ≥ 2/ε₁ in the decay bound and x/H ≤ |D|^κ: E_c ≤ H|D|^{−1}.
- (d) D ∉ E₂ (|d| > |D|^δ ≥ Q₀(X), T₁ ≤ T_X): these zeros have β ≤ σ_g, and number ≤ 2^{ω(D)}CT₁log(|D|T₁):
  E_d ≤ H|D|^{−A/(2c_m)+κ+ε₁+o(1)}.
- (e) These are zeros of L(s, Π_m⊗χ) for primitive χ of conductor ≤ Q₀. L := log(3C(Π_m)Q₀T₁) ≤ 1.3δ log|D|, ℓ := log x
  ≥ A log|D| ≥ 2c_mL. (J5): none has β > 1 − c₁′/L. (J4) with Q = 3Q₀, T = T₁ and ψ = χ₁ if [R]'s exceptional χ₁ for Q₀
  exists: N(σ) := Σ_{cond χ≤3Q₀}N(σ,T₁;Π_m⊗χ) ≤ B_mν^♯e^{c_mL(1−σ)}. Here ν^♯ ≤ min{1, 1.3(1−β₁)δ log|D|} ≤ 1.3ν_D if χ₁ = χ_{d₀}
  with d₀ ∈ 𝔉(D); ν^♯ ≤ 1 = ν_D if χ₁ exists but is not a genus factor of D; ν^♯ = 1 = ν_D if χ₁ does not exist. So
  ν^♯ ≤ 2ν_D, and (using log(3Q₀T₁) ≤ 1.3δ log|D|)
  Σ x^{β−1} ≤ x^{−½}N(½) + ℓ∫_{½}^{1−c₁′/L}x^{σ−1}N(σ)dσ ≤ B_mν^♯[|D|^{−A/2+0.7c_mδ} + (ℓ/(ℓ−c_mL))e^{−(c₁′/L)(ℓ−c_mL)}]
  ≤ 2B_mν_D[|D|^{−A/2+A/10} + 2e^{−c₁′(A/(1.3δ)−c_m)}] ≤ η₀ν_D/8 for |D| ≥ D₀, by δ ≤ δ_e. So E_e ≤ η₀ν_DH/8.
- (f) Σ|R| ≤ Ch(H/x + x^{−½})log(|D|x) ≤ H|D|^{−1} (x ≥ |D|², H ≥ x|D|^{−κ}).

### 5.5 Conclusion
By §5.1 (checked in `params.py`, `params2.py`): ε′ > δ/2; m = 2: s₁ − κ − ε₁ > δ/2; m = 1: s₁ − 3ε₀ − (5/3+4ε₀)(κ+ε₁) > δ/2 and
1/16 − κ − ε₁ > δ/2; A/(2c_m) − κ − ε₁ > δ/2. So E_a + E_b + E_c + E_d + E_f ≤ H|D|^{−δ/2−η₁} for some η₁ > 0, which is
≤ η₀ν_DH/8 for |D| ≥ D₀ because ν_D ≥ c′|D|^{−δ/2}/log|D| ([R] (I5)). With (e): |hS_C^{(m)}| ≤ η₀ν_DH/4 < η₀ν_DH. □

### 5.6 Exceptional zeros: which are excluded, and by what
- *Genus characters* χ = χ_{d,D/d}: L_m(s,χ) = L(s,Π_m⊗χ_d)L(s,Π_m⊗χ_{D/d}) are Dirichlet twists of Π_m. Only these are covered by
  (J5) (Hoffstein–Ramakrishnan for f⊗χ_d, Banks for Sym²f⊗χ_d, Thorner 2025 for every m) and by the repulsion ν^♯ of (J4). They
  enter only (d), (e) and E₂: |d| > Q₀ through E₂ (a count over D), |d| ≤ Q₀ through (J5) and (J4).
- *Non-genus characters* χ ∈ Ĝ_c: L_m(s,χ) = L(s, Π_m × π_χ) with π_χ cuspidal and dihedral; over K it is L(s, BC_K(Π_m) ⊗ χ). Its
  coefficients are real ((R1)), so self-duality does not rule out a real zero near 1, and (J5) says nothing about it. The proof
  never excludes such zeros by a zero-free region. Every zero of L_m(s,χ), χ ∈ Ĝ_c, with β ≥ σ₀ and |γ| ≤ T_X, a real exceptional
  zero included, is counted by N_D^{c,m}(σ₀, T_X); the D for which one of the counts N_D^{c,m}(σ_j, T_X) exceeds its Markov bar
  X^{A(1−σ_j)−2ε′} form E₁^{(m)}(X), the exceptional set of Proposition F (≤ CJX^{1−ε′} discriminants per dyadic X). For
  D ∉ E₁^{(m)}(X), N_D^{c,m}(σ_J, T_X) = 0, so there is no zero at all with β ≥ σ_J = 1 − O(ε′/A), |γ| ≤ T_X, and between σ₀ and σ_J
  the count is at most the bar, which (a) absorbs. So possible exceptional zeros of non-genus twists are removed by discarding the
  D in E₁^{(m)}(X). They do not enter ν_D, which concerns only the Dirichlet characters χ_d (m = 0, [R] (I4)). Zeros with β < σ₀,
  genus or not, are bounded per D in (b) by (R4) and Theorem C.
- *Remark (not used).* Such zeros do not exist anyway: by Ramakrishnan–Wang (m = 1: f is non-dihedral and f ≁ π_χ) and Luo
  (m = 2: π_χ is dihedral), as listed in Thorner–Zhao (arXiv 2601.04189, Proposition 1.1 (4), (6)), and for every m by their
  Theorem 1.7 (π_χ ⊗ χ_D = π_χ while Sym^m f ⊗ χ_D ≠ Sym^m f, so L(s, Sym^m f × π_χ) has no exceptional zero). Even with
  this, E₁ stays necessary: the zeros with σ₀ ≤ β ≤ 1 − c/log(|D|T) need a count, and the only per-D log-free count (BTZ,
  exponent ≈ 10^{10}(1−σ)) is far too weak for A near 3.

## 6. Positive discriminants

**Theorem D⁺.** Theorem D, Corollary D1 (with "indefinite" for "positive definite", h⁺(D) for h(D), narrow classes, and
a² − nb², n ≡ 2, 3 (mod 4), D = 4n) hold for fundamental D > 0, with the same A₀.

*Proof.* The dictionary of [R2] Theorem 4: G⁺, h⁺, Ĝ⁺_c, Ĝ⁺₂; π_χ (χ ∈ Ĝ⁺_c) is a Maass form of eigenvalue ¼ with
π_{χ,∞} = sgn^a ⊞ sgn^a, tempered, self-dual (𝔞𝔞^c = (N𝔞), N𝔞 > 0). Changes here: (R1) is unchanged (χ unramified at finite
places, χ((p)) = 1 at inert p), and so are the local formulas and bounds of (R1′); in (R2)–(R3) the gamma factor depends on
a ∈ {0,1}, so Lemma B (uniform in v) is applied to the two groups
of χ separately; Lemma A holds with [D2] Lemma 5 (𝔔 = 9D, M = 0) in place of [C] Lemma 2; (R4) with h⁺ ≤ 5D^{1/2}log²(3D).
Theorem C (needed for A₀ = 9/4) follows as in §3 (the detector of [R2] Prop 1 covers narrow class characters). Proposition F:
|𝒮⁺| ≤ CX^{3/2}, C(π_χ) ≤ 12D ([R2] Prop 2), I₂ by Lemma B. Genus theory for narrow classes and the factorization of
L_m(s, χ_{d,D/d}) are as in [R2] Theorem 4; (J4), (J5) concern the same twists L(s,Π_m⊗χ_d). The m = 0 input is [R2]
Theorem 4, and the representation step is [R2] Corollary 5. □

## 7. Thresholds and what fixes them

With per-D numerator a and family numerator b (exponents a(1−σ)/(3−2σ) per D and b(1−σ)/(3−2σ) over D), the per-D range is
where A(1−σ) > min{½, a(1−σ)/(3−2σ)}, and Markov needs φ(σ) = (b/(3−2σ) − A)(1−σ) < 1 beyond it (`thresh.py`, `thresh2.py`):

| minorant degree M (top twist m = M) | a, b (this note) | A₀ | with the published densities (a = m+7, b = m+11) |
|---|---|---|---|
| 0 (no condition) | 2, 10 ([R2], [R]) | 3/2 | 8/3 |
| 1 (λ_f(p) > α, α < 0) | 3, 11 | 9/4 | 3 |
| 2 (sign, nonvanishing, …) | 4, 12 | 3 | 10/3 |
| M ≥ 2 | M+2, M+10 | (M+7)/3 | (M+8)/3 |

- For m = 2 and A > 3 the per-D trivial count already covers σ < 1 − 1/(2A), and Proposition F covers the rest:
  φ₂(1 − 1/(2A)) = 6/(A+1) − ½ < 1 ⟺ A > 3. Theorem C is needed only for M = 1.
- At m = 2, A = 3, σ = 5/6 four exponents coincide: X times the trivial per-D count (X^{3/2}), X times Theorem C
  (X^{1+4·(1/6)/(4/3)} = X^{3/2}), Proposition F (X^{12·(1/6)/(4/3)} = X^{3/2}) and the Markov bar X^{1+A/6} = X^{3/2}.
- Payoffs (`hyp.py`, m = 2): A₀ = b/3 − 1 for 6 ≤ b ≤ 12: each factor X^{1/2} saved in V·I₂ of (J3) lowers A₀ by 1/3. A family large
  sieve with threshold |𝒮| = X^{3/2} (instead of X^{7/2}) gives b = 8, A₀ = 5/3; with also Lindelöf on average for the
  degree-6 family (I₂ ≍ X^{3/2}) b = 6, A₀ = 1 (GRH: |D|^{1+ε} for every D).
- All exponents here are per unit window (Theorem C) or for bounded T (Proposition F); in §5 the T₁ ≤ 2|D|^{κ+ε₁} windows cost
  |D|^{O(κ+ε₁)} = |D|^{O(δ)}, which goes into the margins and does not move A₀.
- M ≥ 3 (short intervals of the Sato–Tate parameter), now unconditional.
  **Theorem D_M.** Let M ≥ 2, P = Σ_{m=0}^{M}c_mU_m(x/2) ≤ 1_S on [−2,2], c₀ > 0, c_M ≠ 0. Then Theorem D (and Theorem D⁺)
  holds for every A > (M+7)/3, with η₀ := c₀/(16Σ_{m≥1}|c_m|).
  *Proof.* Proposition T for m = 1 (as proved) and for every 2 ≤ m ≤ M with A > A_m := (m+7)/3, exactly as for m = 2
  (A₂ = 3): σ₀ := 1 − 1/(2A_m), s₁ := (A − A_m)/(2A_m), ε′ := s₁/4. Then φ_m(σ₀) = 1 − (A − A_m)/(2A_m) = 1 − 4ε′, φ_m′ = A − (m+10)/(3−2σ)² < 0 on [σ₀, 1]
  for A ≤ A_m + ¼, the T_X-exponent (2m+10)(1−σ)/(3−2σ) is ≤ 3 there, and the E₁ margin is 3(κ+ε₁) + (m+11)ε₀ ≤ ε′
  (`params2.py`, m ≤ 20). In (b) the constant 10 becomes 2A_m, since 1 − σ ≥ 1/(2A_m). Automorphy of Sym^m f: Newton–Thorne;
  (J3), (J4) need only GRC; (J5) for m ≥ 3: Thorner (IMRN 2025, arXiv 2404.06482, Theorem 1.1 and Corollary 1.2 with F = Q,
  n = 0; eq. (5.1) of its proof), whose region c/(m²[m^{y}log(mC(f)) + log(C(χ)(3+|t|))]) is uniform in χ, so c₁′ := c/(2m²) for |D| ≥ D₀. Positivity
  as in §5 with Σ_{m≥1}|c_m|. □ (The first version of this note made M ≥ 3 conditional on (J5); Thorner's theorem removes that.)

## 8. Novelty audit (2026-10-21)
- arXiv API, 35 queries (`axq1.out`…`axq4.out`): Sato–Tate × quadratic form / ring class / class group / least prime / imaginary
  quadratic ideal class / split primes / almost all discriminants; Chebotarev–Sato–Tate (joint, effective); Rankin–Selberg ×
  class group (characters) × zero(s) / moment; zero density × Rankin–Selberg × family; zero density × class group; theta
  series × zero density; ring class characters × L-functions; Hecke eigenvalues × binary quadratic form; sign changes ×
  binary quadratic; least prime × Hecke eigenvalue × class; titles of HT, BTZ, LOT, Thorner, CPS, Nordentoft.
- Relevant hits, read at abstract level or in TeX (`src_*`): Pascadi–Thorner 2508.14888 (Theorems 1.3, 1.5: general families;
  weaker here, §3–§4), Humphries–Thorner 2103.05634, Brumley–Thorner–Zaman 1804.06402 (inputs), Lemke Oliver–Thorner 1505.03122
  and Thorner 2002.10450 (least prime with θ_p ∈ I in the level aspect; joint Sato–Tate of two forms), Chen–Park–Swaminathan
  1507.07122 (GRH, level aspect), Nordentoft 2105.09130 and Burungale–Tian 1712.01465 (nonvanishing of L(f⊗χ,½) over class
  group characters), Van Order 1903.06722 (class group twists of GL_n, central values), Vaishya et al. 2204.08310, 2401.18055,
  2601.16352 (sign changes of λ_f(n) at integers represented by one fixed form with h(D) = 1), Yang 2610.22728 (Motohashi formula
  for class group L-functions), Bucur–Kedlaya 1301.0139, Amri 2203.11498, Lee 2607.18832, Yang 2609.07080, Malik–Prabhu
  2012.12534 (Chebotarev–Sato–Tate for a fixed field), and the least-prime-in-a-class papers of [R2] §6 (Ditchen, Sardari,
  Gaudet, Zaman), none with a Hecke-eigenvalue condition.
- GXL evidence API, 4 claims (`gxl.out`): almost-all-D primes with a_E(p) > 0 represented by every form; zero density for
  L(f × θ_χ) over class group characters; least prime ideal in a class with λ_f(p) in an interval polynomial in |D|; effective
  Chebotarev–Sato–Tate in ring class fields. Support only for the unrestricted least-prime results (Ditchen, Sardari, Gaudet)
  and for fixed-field effective Sato–Tate (Bucur–Kedlaya, abelian varieties); nothing for the twisted density or for
  Hecke-eigenvalue conditions in growing class groups.
- Conclusion: Theorem C, Proposition F and Theorems D, D⁺ were not found. Per D, the joint (class, Sato–Tate) least-prime
  problem has, as far as found, no stated polynomial bound at all; fixed-field Chebotarev–Sato–Tate results are exponential in
  h(D) when applied to ring class fields.
- Second audit after the review (2026-10-21, [NOVELTY.md](NOVELTY.md), `.scratch/rsfamfix/`): 65 more arXiv API queries and 2 author lists, OpenAlex
  (25 searches, 15 of them replacing arXiv queries refused by the API rate limit, and the work lists of Thorner, Wong, Zaman,
  Lemke Oliver), Zenodo (15 queries, 2025–2026), Hexagon (all 323
  works), openai/math (749 preprints), 5 GXL claims, and the TeX of Thorner 2404.06482, Thorner–Zhao 2601.04189, Akbary–Wong
  2103.15713. Still nothing with an ideal-class (or form) condition and a sign or Sato–Tate condition, uniform in D. One caveat:
  a per-D bound |D|^{C} with C of order 10^{10} looks derivable from published results (BTZ Theorems 1.2 and 1.4 for the dihedral
  family; Ramakrishnan–Wang, Luo, Thorner–Zhao, Hoffstein–Ramakrishnan, Banks and Thorner 2025 for exceptional zeros;
  Thorner–Zaman for m = 0) by Linnik's method; we did not write it out and found it stated nowhere. Theorem D's content is the
  exponent 3 (9/4 one-sided) for almost all D, with short intervals.

## 9. Dead ends (NOGO candidates)
- Sign with a degree-one minorant: impossible (Remark in §5). Sign needs the degree-6 family Sym²f ⊗ θ_χ, so A₀ = 3.
- Below A = 3 for the sign: the four-way coincidence of §7. The per-D second moment |D|^{3/2} for Sym²f ⊗ θ_χ is the large
  sieve diagonal at AFE length |D|^{3/2} ≫ h; class II is sharp for it ([R2] §4 argument).
- Subconvexity (Michel–Venkatesh, L(½, f⊗θ_χ) ≪ |D|^{½−δ}) in class II: sup × mean square gives |𝒳_II| ≤ Y^{1−2σ}|D|^{3/2−2δ},
  the square of the Cauchy–Schwarz bound up to |D|^{−2δ}; worse whenever the density bound is ≥ 1. No use.
- Hölder with fourth moments (|D|^{m+1} from Lemma A): exponent (m+2)(1−σ)/(2−σ) > (m+2)(1−σ)/(3−2σ). Worse.
- Per-D sieve inside the family detector: reproduces X times Theorem C (the sum of per-D bounds), no gain over the endpoints.
- Lindelöf on average for Σ_χ|L(½, f⊗θ_χ)|² via Waldspurger/Gross periods: equals a fourth moment of toric periods over the
  class group, i.e. averaged mixing (Michel–Venkatesh conjecture; Khayutin, Blomer–Brumley–Khayutin under GRH). Out of reach.
- (Withdrawn 2026-10-21: "intervals needing M ≥ 3 are blocked by Landau–Siegel zeros of Sym^m f ⊗ χ_d". Thorner, IMRN 2025,
  excludes those zeros for every m; Theorem D_M of §7 is unconditional.)
- Lead (not pursued): a log-free form of Theorem C (no |D|^ε loss near σ = 1) together with the exceptional-zero results of §5.6
  would give a per-D theorem with a moderate exponent instead of the almost-all statement.

## 10. Repairs after review (2026-10-21)
The reviewer judged the first version FLAWED as written, repairable. Items and where they are fixed:
1. *Unit windows.* Theorem C counts zeros with |γ − t| ≤ 1; (C′) states the all-heights version with the extra factor T. Uses:
   §3 comparison (PT count |γ| ≤ T; T-exponents compared), §5.4(b) (≤ 3T₁ windows written out; the m = 1 margin is now
   s₁ − 3ε₀ − (5/3+4ε₀)(κ+ε₁) > δ/2), §7 (per window, resp. bounded T), REPORT.
2. *Lemma B tail.* The old tail bound O(|D|^{−8}) was not uniform in V and v. New cutoff P = VL₀(3|D|V(|v|+3))^ε with
   B = 1 + (m+13)/ε gives a tail ≤ C(ε)(|D|V(|v|+3))^{−10} for every χ; the statement is unchanged and now proved uniformly in
   V ≥ 1, v ∈ R. Uses re-checked: only V = 1 and the mollifier bound are used (Theorem C class II, Proposition F, Theorem D⁺).
3. *Coefficients.* (R1′) gives the base-change coefficients λ_m(𝔞) prime by prime (split, ramified, inert, p | N), shows where they
   differ from λ_{Sym^m f}(N𝔞), gives Rankin's identity for the brief's λ_f(N𝔞)χ(𝔞), and proves |λ_m(𝔞)| ≤ d_{m+1}(𝔞),
   |μ_m(𝔞)| ≤ d_{m+1}(𝔞), |λ_f(N𝔞)| ≤ d_4(𝔞); Lemma A uses nothing else.
4. *Exceptional zeros.* §5.6: (J5) and the repulsion of (J4) cover only the genus (Dirichlet) twists; possible exceptional zeros of
   the non-genus twists are removed by discarding the exceptional set E₁ of Proposition F (for D ∉ E₁ there is no non-genus zero
   with β ≥ σ_J at all). (Remark: such zeros do not exist anyway, by Ramakrishnan–Wang, Luo, Thorner–Zhao; not used.)
5. *Minorants.* §5, "Positivity": c₀ > 0 gives positivity only together with the m = 0 lower bound ν_DH/8 and m = 1, 2 errors
   ≤ η₀ν_DH (which is why (e) needs the repulsion ν^♯ ≤ 2ν_D); table of c₀, c₁, c₂, η₀ and the constant of (5.1) for every set of
   Corollary D1. `params2.py` → `params2.out`: 85 checks, all OK (minorants, margins, δ_e, the Lemma B exponent, Theorem D_M).
Found while repairing: Thorner (IMRN 2025) excludes Landau–Siegel zeros of Sym^m f ⊗ χ for every m, so the thresholds
(M+7)/3 for minorant degree M ≥ 3 are unconditional (Theorem D_M, §7); the §9 entry saying otherwise is withdrawn.
Citation check (2026-10-21, [ChebST FIX.md](../ChebST/FIX.md)): a reviewer did not find the (J5) region in Thorner arXiv 2002.10450,
which is not the cited paper. In the TeX of the cited one (arXiv 2404.06482v2 = IMRN 2025, rnaf307) it is eq. (5.1), for every m ≥ 1
and every Hecke character; (J5) is now pinned to it. No statement or constant changes.
