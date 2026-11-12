# CMAudit derivations (2026-11-12)

Re-derivations behind [REPORT.md](REPORT.md). [CMQ] = [CMQuartic NOTE](../CMQuartic/NOTE.md), [CM2] = [CM2 NOTE](../CM2/NOTE.md), [CM3] = [CM3 NOTE](../CM3/NOTE.md); "l." = NOTE line. Notation as there: F totally real, n = [F:Q]; K/F CM, c = conjugation; 𝔇 = 𝔇_{K/F}, D = N𝔡_{K/F}, d_K = d_F²D; 𝔤(𝔞) = O_F-content.

## D1. [CMQ] Lemma 1 (l.64–75)
- O_{K,𝔭} is free of rank 2 over the DVR O_{F,𝔭} and 1 ∉ πO_{K,𝔭}, so 1 extends to a basis {1, y} and O_{K,𝔭} = O_{F,𝔭}[y]; for a monogenic order 𝔇O_{K,𝔭} = (f′(y)) = (y − yᶜ).
- x = ββ′ᶜ ∈ 𝔠𝔠ᶜ = N_{K/F}(𝔠)O_K; locally x = g(a + by), x − xᶜ = gb(y − yᶜ). So (Δ) ⊆ N_{K/F}(𝔠)𝔇 and |NΔ| ≥ N(N_{K/F}(𝔠)O_K)·N𝔇 = (N𝔠)²D.
- σ_v∘c = σ̄_v: σ_vΔ = z − z̄, z = σ_vβ·σ̄_vβ′; |N_{K/Q}Δ| = Π_v|σ_vΔ|². Δ = 0 ⟺ β′/β ∈ F. Correct.
- Remark 2.3 (l.111): β′ = 1, 𝔠 = O_K: β − βᶜ ∈ 𝔇∖0, |σ_v(β − βᶜ)| ≤ 2|σ_vβ|, so Nβ ≥ D/4ⁿ = d_K/(4ⁿd_F²). Same bound as Miller–Peckner's key lemma (TeX l.653–661: p ≥ Δ/(4^N D_{K⁺}²), written there only for O_K = O_{K⁺}[√β]).

## D2. [CMQ] Lemma 2 (l.79–99) and §7.2 (l.299–308)
- (i): λO_F = 𝔪𝔫^{−1} coprime, 𝔫O_K𝔞′ = 𝔪O_K𝔞; 𝔭 | 𝔫 ⇒ v_𝔓(𝔞) ≥ e(𝔓/𝔭) for every 𝔓 | 𝔭 ⇒ 𝔭O_K | 𝔞. So distinct O_F-primitive ideals of one class have generator ratios ∉ F, for every choice of generators.
- Original form: generators of 𝔞𝔠 form a coset of O_K^× ⊇ O_F^×, so one with ℓ(β) in a fundamental domain 𝒫 of λ(O_F^×) exists; |σ_vβ|² ≤ (N·N𝔠)^{1/n}e^{r_F}; two ideals in one of kⁿ angle boxes contradict Lemma 1; kⁿ ≤ (1 + y)ⁿ ≤ 2^{n−1}(1 + yⁿ). Q(√5): 2(4π)²ε_F²d_F = 2·157.91·2.618·5 = 4134. Correct.
- Strengthened form: Π_v|σ_vβ| = (Nβ)^{1/2}, so two ideals in one box give N𝔠√D < (4π/k)ⁿN·N𝔠 and π_C(N) ≤ (1 + 4π(N/√D)^{1/n})ⁿ ≤ 2 + 2(4π)²d_FN/√d_K (= 2 + 1579N/√d_K for Q(√5)). No unit, no regulator. Correct.
- Constant lost: mean count per class ≈ (2π)ⁿR_KN/(w_Kζ_F(2)√d_K) (Remark 2.1). Q(√5), R_K = 2 log ε_F, w_K = 2: mean ≈ 16.42N/√d_K against bound ≈ (4π)²·5·N/√d_K = 789.6N/√d_K: ratio 0.0208. Measured max_Cπ_C/bound at N = 16√d_K: 0.0205–0.0220 on all 13 Q(√5) fields.

## D3. [CMQ] Lemmas 3, 3′ (l.118–135)
- 𝔞 = 𝔤O_K·𝔟 with 𝔟 primitive is unique: v_𝔭(𝔤(𝔞)) = min_{𝔓|𝔭}⌊v_𝔓(𝔞)/e(𝔓/𝔭)⌋ (split: min of the two valuations; inert: v; ramified: ⌊v/2⌋). χ ∈ 𝒵_K kills 𝔤O_K.
- Sharp constant: with w_𝔞 := N𝔤(𝔞)^{1+2ε} the dual Gram G_{χχ′} = Σ_𝔞 w_𝔞^{−1}χ(𝔞)χ̄′(𝔞) over Ĉl_K is a convolution on Cl_K with eigenvalues h_K·Σ_{𝔞∈C}w_𝔞^{−1}. So the best constant is exactly h_K·max_CΣ_{𝔞∈C,N𝔞≤N}w_𝔞^{−1} (for 𝒵_K: |𝒵_K|·max over cosets of j(Cl_F)). Direct diagonalisation agrees (EIG lines in checks/, h = 11 and 36: relative difference ≤ 10^{−35}).
- Lemma 3 proof: Cauchy–Schwarz in 𝔤 (ζ_F(1+2ε)), orthogonality on Cl_K/j(Cl_F), Cauchy–Schwarz inside a coset (≤ |j(Cl_F)|·max_Cπ_C(N) primitive ideals), |𝒵_K||j(Cl_F)| = h_K. Correct.
- Fact (c): h_K = w_K√d_Kκ_K/((2π)ⁿR_K) (r₁ = 0, r₂ = n); Louboutin κ_K ≤ (e log d_K/(2(2n−1)))^{2n−1}; Friedman R_K ≥ 0.2. Uniform in K. Correct.
- Large-N limit: Σ_{𝔞∈C}N𝔤^{−2} ≈ (κ_KN/h_K)·ζ_F(4)/ζ_F(2); Q(√5): 1.010/1.157 = 0.873. Measured K/((N + √d_K)κ_K) at N = 16√d_K: 0.82–0.88 (Q(√5); Q(√2) 0.70–0.74, Q(√10) 0.60–0.62).
- Lemma 3′: N𝔤^{1+2ε}/N𝔞 = N𝔤^{2ε−1}/N𝔟 ≤ M^εN𝔤^{−1}N𝔟^{−1}; d_k(𝔤O_K) ≤ d_k(𝔤)² prime by prime (ramified: binom(2a+k−1,k−1) ≤ binom(a+k−1,k−1)²). Correct.

## D4. [CMQ] Theorem A (l.139–168)
- Gamma factor Γ_C(s)ⁿ for every χ ∈ Ĉl_K (finite-order characters of C^× are trivial), conductor d_K, analytic conductor d_K(|t|+3)^{2n}, AFE length √d_KUⁿ.
- (2.1′) (l.153–157): the dual AFE term for χ = ψ₀ξ, ξ ∈ 𝒵_K, has χ̄(𝔞)χ(𝔨) = ψ̄₀(𝔞)ψ₀(𝔨)ξ(𝔞ᶜ𝔨), linear in ξ with ψ₀ fixed; the root number has modulus 1 and drops out after |LM| ≤ |D₁M| + |D₂M|. Lemma 3 at length X√d_KUⁿ with Lemma 3′ gives d_K^ε(X√d_KUⁿ)^{1+ε}. Correct.
- Class II: |𝒳_II| ≤ Y^{1/2−σ}(Σ_χ∫|L|²)^{1/2}(Σ_χ∫|M_V|²)^{1/2} ≈ Y^{1/2−σ}(√d_KTⁿ)^{1/2}(√d_K)^{1/2} = Y^{1/2−σ}d_K^{1/2}T^{n/2}. Class I: Lemma 3 at M ≥ V = √d_K has constant ≤ 2(B/ε)M, giving Z^{2−2σ+ε}. Balance Y^{3/2−σ} = d_K^{1/2}T^{n/2}: Y = (d_KTⁿ)^{1/(3−2σ)}, count (d_KTⁿ)^{2(1−σ)/(3−2σ)}. Correct.
- Mollifier content: squarefree 𝔨 may have content (𝔭O_K = 𝔓𝔓ᶜ), covered by Lemma 3′; in [CM2] Lemma 5(c) the pieces M_𝔱 run over squarefree relative norms, hence primitive ideals ([R3] Lemma 3(c), expansion via θ_p). Correct.
- σ-range: e(σ) ≥ ½ ⟺ σ ≤ ½, so for σ ≤ ½ the bound follows from (b) (c₀h_K log ≤ d_K^{1/2+ε}T^ε). On ½ < σ < 0.6 it is not proved ([C] Lemma 5 needs σ ≥ 0.6); routine with the detector of [CM2] Lemma 6 (σ ≥ ½ + 2δ). The discovery file states claim 2 with no σ-range.
- Comparison (l.166–168): Pascadi–Thorner Thm 1.3 (TeX v4 l.374–381) with n = 1 over K, Q ≍ d_K, |S| = h_K: (d_K^{−2}d_K³d_K²)^{(1−σ)/(3−2σ)} = d_K^{3(1−σ)/(3−2σ)}. Correct.

## D5. [CMQ] Proposition C (l.178–189)
- π_χ = AI(χ), χ ≠ χᶜ: cuspidal, archimedean type D₁ (μ = 0, 1; Γ_R(s)Γ_R(s+1) = Γ_C(s)), C(π_χ) = d_F²·D·Π_v(3+0)(3+1) = 12ⁿd_F²D, GRC, finite-order central character. Each π comes from ≤ 3 quadratic extensions with 2 characters each: multiplicity ≤ 6.
- HT (TeX: cor:MVT l.841–850, bound_int l.949–952): V = Q^{4m²θ_m+m+ε}T^{m²[F:Q]/2+1+ε}|𝒮|; θ_m = 0 under GRC gives Q²T^{2n+1}|𝒮| (m = 2), as quoted. Base field arbitrary.
- V ≍ X²·T^{2n+1}·X^{3/2}; I₂ ≪ X^{3/2+ε}T^{n+1+ε} ((2.1′) at X = 1; Euler factors at 𝔡 cost ≤ 2^{1/δ₀}). Y = (VX^{3/2}T^{n+1})^{1/(3−2σ)} gives Y^{1/2−σ}(VI₂)^{1/2} = Y^{2(1−σ)} = (X⁵T^{3n+2})^{2(1−σ)/(3−2σ)} = (XT^{(3n+2)/5})^{10(1−σ)/(3−2σ)}. Correct.

## D6. [CMQ] Lemma 4 (l.193–209) and (5.2) (l.240–247)
- Gal(M/F) is an extension of Z/2 by the cyclic Gal(M/K) with trivial action (χ∘c = χ), hence abelian; inertia has order ≤ 2; for odd 𝔭 the characters of order ≤ 2 of O_𝔭^× are 1 and the Legendre symbol; N(O_𝔓^×) = (O_𝔭^×)² for ramified odd 𝔭; kernel {1, η} by class field theory; Ind(ψ∘N) = ψ ⊕ ψη. Correct.
- ε_ψ(C) = ψ̄(N𝔞_C) does not depend on 𝔞_C ∈ C prime to 𝔡: a ratio λ has Nλ ≫ 0 and Nλ ∈ (O_𝔭^×)² at 𝔭 | 𝔡. η ∈ 𝒢(K)∖𝒢₀ (conductor 𝔡), so its zeros fall under E₂. Correct.
- Numerics: |Cl_K/(1−c)Cl_K| = |𝒢(K)|/2 in 26/26 fields over Q(√5), Q(√2), Q(√3) (h⁺ = 2), Q(√10) (h = h⁺ = 2), 𝒢(K) computed independently from bnrinit(F, [𝔡, [1,1]]) (checks/out_genus.txt).

## D7. [CMQ] Theorem B (l.211–263)
- A = 1.55, σ₀ = 0.84: φ(σ₀) = (10/1.32 − 1.55)·0.16 = 0.96412, ε′ = 0.00897, s₁ = 0.16·(1.55 − 1.51515) = 0.005576, as stated.
- E₁: X^{φ(σ_j)+2ε′+ε₀+(16/7)(κ+ε₁)} ≤ X^{1−2ε′+0.0028} ≤ X^{1−ε′} at κ = 10^{−3}. E_b exponent: −s₁ + ε₀ + (1 + 8/9)(κ+ε₁) = −0.0034 (κ = 10^{−3}), −0.0045 (κ = 5·10^{−4}). E_d: −0.005A + 1.1κ < 0. E₂: 1 − 0.965δ₀ + 0.01δ₀ + 0.011δ₀ ≤ 1 − 0.9δ₀.
- (e): Σx^{β−1} ≤ x^{σ₀−1}N_M(σ₀,T₁) + y∫_{c_M/L}^{1/2}C_Me^{−u(y−B_ML)}du ≤ o(1) + 2C_Me^{−c_My/(2L)} for y ≥ 2B_ML; κ₀ as stated gives ≤ 1/32. TZ Theorem 3.2 confirmed in TeX (l.1137–1146): {e^{O(n_K)}D_K²QT^{n_K+2}}^{81(1−σ)}, ½ ≤ σ < 1, T ≥ max{n_K^{5/6}D_K^{−4/(3n_K)}Q^{−4/(9n_K)}, 1}.
- (f): h_K(H/x + x^{−1/2}) ≤ HD^{−1} for x ≥ D^{3/2}, H ≥ xD^{−κ}. Prime powers: ≤ C(Hx^{−1/2} + log x)log^{2n+1}x ≤ H/(16h_K). Correct.

## D8. [CM2]
- Lemma 4 units (l.107–112): "mu ≫ 0" depends only on u mod U² (squares are totally positive), so R(m) is well defined; each U-orbit meets 𝒜′ in one pair ±α; ker(N on O_𝔣^×) = μ(O_𝔣) (CM), so [O_𝔣^× : μU] = [N(O_𝔣^×) : U²] = q and [O_𝔣^× : U] = q|μ(O_𝔣)|/2; |𝒜/O_𝔣^×| = |𝒜′|/ω. Correct.
- Classes (l.101–106): S^∨ = j^{−1}Σb̄T_𝔫 from [𝔞]^{−1} = [𝔞ᶜ]j_𝔣(N𝔞)^{−1}; j = j_𝔣(𝔯)^{−1} since 𝔫𝔯 = (ν_𝔫) is prime to 𝔣. Hecke relations (Lemma 3) re-derived in the split, inert, ramified and 𝔭 | 𝔣 cases. Correct.
- Lemma 7 (l.180–187): sub-box differences k ∈ 𝔞 have Π|k_v| ≤ 2ⁿΠZ_v/mⁿ ≤ N𝔞/2, so k = 0; only Π_vZ_v enters because every nonzero k ∈ 𝔞 has Π|k_v| ≥ N𝔞. Numerics: 800 random boxes in ideals of Q(√5), Q(√2), Q(√10), skew Z₁/Z₂ up to 10^{±6}, Z₁Z₂/N𝔞 ∈ {0.3, 1, 4, 16}: max count/bound 0.167.
- Lemma 9 (l.212–229): ∫f = 2ⁿ√M·N𝔶²v·K_n(N𝔶²v) (t_v = 2√μ_vτ_v, Y/(4ⁿM) = N𝔶²v); with √M = N𝔩√M₀ the main term is (2^{1−n}/√d_F)(μ(𝔠)/N𝔠)(X′/√M₀)ΣA_μ(𝔶)K_n(N𝔶²v). Local factor at 𝔭 | 𝔠: −N𝔭^{−1} (Jacobsthal, k = 0) and 0 (k ≥ 1: linear character sum over the lifts). Dual lattice 𝔮^{−1}𝔇_F^{−1}, covolume N𝔮√d_F, R ≍ X′/(N𝔠N𝔰d_F√M₀). Correct.
- Lemma 10 (l.231–252): (a) Hensel; (b) y² ≡ μ mod 𝔭^m forces y ≡ ±y₀ mod 𝔭^{m−e₂}, so A ≤ 2N𝔭^{2e₂}; the dependence on μ mod squares follows from C_k(x,μ) ⟺ C_k(λx,λ²μ) alone (the sentence on the generated group is vacuous). Correct.
- Lemma 11 (l.282–297): Σ_{u∈U/U²}Ξ(u) = 2ⁿ1[Ξ|_U = 1] ([U:U²] = 2ⁿ); 1 + 8O_𝔭 ⊆ (O_𝔭^×)² at dyadic 𝔭; Θ̃ is multiplicative because g_{𝔞𝔞′} = g_𝔞g_𝔞′Πγ_i^{c_i}·unit and Θ̃(𝔤_i)^{h_i} = Ξ(γ_i). Correct.
- §4.7 (l.299–323): Q² = E₀²/N(T₂T₃T₄T₅)² ≤ E₀²/(NκNκ′) ≤ L/(NκNκ′) when E₀² ≤ L; Σ_{T₂}N(T₂)^{−1/2}A_{T₂κ₀}^{1/2} ≤ (log)^{1/2}(Σd²|β|²)^{1/2}. Lemma 12 by duality and convexity at Re s = −ε. Correct.
- Lemma 13 (l.327–343): the residue at s′ = 1 is ∝ ∫W(v)v^{−σ−it′}dv = 0 by (Z) ((Z) re-derived by u = 4v). Lemma 14 (l.345–371): orthogonality in q through the quadratic characters χ_x, x reduced mod squares to x ∈ O_F with N(x) ≪ N𝔞N𝔞′. Correct.

## D9. [CM3]
- Lemma P_F (l.30–63): (a) 𝔭 | 𝔣 with v_𝔭(𝔫𝔫′) = 1 gives zero ([CM2] Lemma 2(c)); (b) {α ∈ O_{𝔭,f} : 𝔭² | Nα} = πO_{𝔭,f−1} (α = a + π^fβ, Nα ≡ a² mod π); (c) shift by λ ∈ F^× with 𝔩 = λ𝔩′, 𝔩′ prime to 𝔣: N(α/λ)O_F = 𝔩′²𝔐, and orbits ↔ integral 𝔟 ∈ I_K(𝔣_𝔢) of norm 𝔐 in the class j(𝔩′)^{−1} = j(A[𝔢]^{−1})j(𝔠)^{−1}; (d) Hecke relations in Z[G_{K,𝔣_𝔢}] and S^∨ = j′^{−1}S̃. Correct; every block is ω·Σ_C|S_C|² ≥ 0.
- Theorem 1 (l.83–102), 𝔴 = O_F: N𝔠√M₀ = √(N𝔫N𝔫′) ≤ 4N, so C₃·4N·X^δ ≤ X for N ≤ X^{1−4δ}; M₀ ≤ 16N² ≤ X²; E₀ = 2^{n+3}N/(N𝔠√X) ≥ 2v^{−1/2} (v^{−1/2} = 2ⁿ√M₀/√X ≤ 2^{n+2}N/(N𝔠√X)); E₀²/L = 4^{n+2}N/(N𝔠X) ≤ 1. Correct.
- Theorem 2 (l.108–117): (1−5δ)(½+δ−σ) ≤ ½ − σ + 3δ for σ ≤ 9/10; X^{3/2}V^{1−2σ} = 2^{2σ−1}X^{2−σ}. Correct.
- Theorem 3 (l.121–156): φ = 2 − σ − A(1−σ) = 1 − (A−1)(1−σ) ≤ 1 − (A−1)/10 = 1 − 8ε′ on [σ₀, 9/10]; F₂(9/10) = 5/6; #E₁ ≤ JX^{1−6ε′+ε′/10}; (b) needs (A−1)/2 − 6δ > 0. Correct.

## D10. Thresholds
- F₁ = F₂: 4σ² − 40σ + 31 = 0, σ_c = 5 − √69/2 = 0.846688; 1 + A₀(1−σ_c) = F₂(σ_c) = 1.17334 gives A₀ = 1.130662 = (3+√69)/10.
- 2 − σ = F₂: 2σ² + 3σ − 4 = 0, σ = (√41−3)/4 = 0.8508. (13−8σ)/(7−2σ) − (2−σ) = (2σ−1)(1−σ)/(7−2σ) ≥ 0.
