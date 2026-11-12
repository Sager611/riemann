# RC5Audit follow-up derivations (2026-11-12, orchestrator 31)

Scope: (a) [R2] = [RCApps2 NOTE](../RCApps2/NOTE.md) Propositions 1–2 and their use in RCApps5; (b) the Humphries–Thorner input (I3) of [R] = [RCApps NOTE](../RCApps/NOTE.md); (c) [R4] = [RCApps4 NOTE](../RCApps4/NOTE.md) as RCApps5 §6 and `discoveries/least-primes-forms.md` use it. [C] = Combos NOTE, [D2] = Density2 NOTE, [R3] = RCApps3 NOTE, [R5] = RCApps5 NOTE. Line numbers are those of each NOTE.

## F1. [R2] Proposition 1 (l.28–56): re-derived, correct

Hypotheses ([D2] Prop 6, [D2] l.131–143): Z-trivial family 𝒵 on I_K(f) with (LS) Σ_ξ|Σc_𝔞ξ(𝔞)|² ≤ (B/ε)(N+Λ)‖c‖²_ε, B ≤ 50log²(3𝔔), Λ ≤ 2𝔔^{1/2}; ℱ = ψ₀𝒵; 𝔔 = |D|(M₀+3)²; gamma shifts |μ_j| ≤ M₀+1.
- |ℱ|: (LS) with c = δ_{O_K} (N = 1, ‖c‖² = 1, ε = ½) gives |ℱ| ≤ 2B(1+Λ). Small 𝔔: |ℱ|·c₀log(16𝔔T²) ≤ C(ε)T^ε.
- Parameters: V = 𝔔^{1/2}, Y = (𝔔T)^{1/(3−2σ)} ≥ (𝔔T)^{1/2} ≥ V since 3 − 2σ ≤ 2; Z = 2Y log²Y.
- Detector: [C] Lemma 5 (refereed by RCAudit) with mollifier length X := V ≤ Y; valid for ℱ by [D2] Prop 6 step 3 (L(s,χ) = E_χ L(s,χ*) entire for χ ≠ 1). Zeros of L(s,χ) and L(s,χ*) with β ≥ 0.6 coincide (E_χ vanishes only on Re s = 0).
- Class II: e^{−|v−γ_χ|} ≤ e^{1−|v−t|}; Cauchy–Schwarz over all of ℱ∖{1}; Σ_χ|L(½+iv)|² ≤ C𝔔^{1/2+O(ε₀)}(|v|+3)^{1+ε₀} ([D2] (2.1) at mollifier length 1); Σ_χ|M_V|² ≤ (B/ε₀)(V+Λ)Σ_{N𝔡≤V}m(𝔡)^{1+2ε₀}N𝔡^{−1} ≤ C𝔔^{1/2+ε₀} ([C] Lemma 3). So |𝒳_II| ≤ C Y^{1/2−σ}𝔔^{1/2+ε₀}T^{1/2+ε₀}.
- Class I: a_𝔞 = 0 for 1 < N𝔞 ≤ V, so every block has M = 2^{j+1} > V = 𝔔^{1/2} ≥ Λ/2 and the (LS) constant is ≤ 3(B/ε₀)M: no k-th powers needed. Sobolev and |b_𝔞| ≤ d(𝔞)(M/4)^{−σ}: Σ_χ sup|P_j|² ≤ C log⁴·M·M^{−2σ}·M^{1+2ε₀}; |𝒳_I| ≤ C(𝔔T)^{ε₀}Z^{2−2σ+2ε₀}.
- Exponents: Y^{2−2σ} = (𝔔T)^{(2−2σ)/(3−2σ)} = (𝔔T)^{e(σ)}; Y^{1/2−σ}(𝔔T)^{1/2} = (𝔔T)^{(1/2−σ+3/2−σ)/(3−2σ)} = (𝔔T)^{e(σ)}.
- e(σ) < 2σ(1−σ) on (½, 1): 2σ(1−σ)/e(σ) = σ(3−2σ) and σ(3−2σ) − 1 = (2σ−1)(1−σ) > 0. In family-size terms (|ℱ| ≍ 𝔔^{1/2}) it is |ℱ|^{4(1−σ)/(3−2σ)}, Montgomery's argument with large-sieve constant M + |ℱ|.
- Theorem 1′(a) (ring class characters of imaginary quadratic orders, (LS) = [C] Lemma 2, (2.1) = [C] Lemma 4 at X = 1): all inputs are refereed (RCAudit), so it holds as stated. Theorem 1′(c) uses [D2] Lemma 5 (h⁺π⁺_C ≤ (h⁺log η + h⁺)(1 + 2.351N/D^{1/2}) by the hyperbolic count [D2] Lemma 3); argument checked, not tested numerically.
- Role in [R5]: none load-bearing. [R5] cites it only in the comparison at l.139; Theorems 3, 3⁺, 3^∠ of [R5] count zeros with β < ½ + 4δ trivially ([C] Thm 1(b); h⁺ per window for D > 0). It was load-bearing only in [R2] Theorems 2′, 4 (A > 3/2), now superseded.

## F2. [R2] Proposition 2 (l.112–122): correct; one cosmetic constant

- χ ∈ Ĝ⁺_c: 𝔞𝔞^c = (N𝔞) with N𝔞 > 0 totally positive ⇒ χ^c = χ̄; AI_{K/Q}(χ) cuspidal ⟺ χ ≠ χ^c ⟺ χ² ≠ 1; π_χ = π_{χ′} ⟺ χ′ ∈ {χ, χ̄}.
- ∞-type: χ((α)) = sgn(Nα)^a for α ≡ 1 mod^× f ⇒ π_∞ = π(sgn^a, sgn^a), Laplace eigenvalue ¼, tempered; finite Satake parameters χ(𝔭), χ(𝔭^c) or ±1: GRC.
- Conductor: in HT's normalisation C(π) = N𝔮_π Π_j(3 + |μ_j|) (TeX eqn:analytic_conductor_def), μ = (0,0) for a = 0 and (1,1) for a = 1, so C(π_χ) ≤ 16D, not 12D. Q := 16X changes only constants.
- |𝒮| ≤ Σ_{D≤X}h⁺ ≤ CX^{3/2}: h⁺ ≤ 2h, hR = ½√D L(1,χ_D), R ≥ log((1+√D)/2).
- |L^{ur}(½+iv,π_χ)| ≤ 2^{ω(D)}|L(½+iv,χ)| (ramified factors 1 − χ(𝔭_p)p^{−s}); I₂ ≤ C X^{3/2+2ε}T^{2+ε} from [D2] (2.1) for Ĝ⁺ (two archimedean types; Prop 6 allows χ-dependent shifts).
- Exponent: V·I₂ ≪ X^{5+O(ε)}T^{5+O(ε)} and Y = (VX^{3/2}T²)^{1/(3−2σ)} give (XT)^{10(1−σ)/(3−2σ)+ε}, as in [R] Prop 1.

## F3. How [R] Prop 1, [R2] Props 1–2 enter [R5]

| [R5] step | uses | where |
|---|---|---|
| Thm 3 step 2′, σ ≤ 9/10 | [R] Prop 1 at σ = 9/10 for the zeros with β > 9/10 (X^{5/6+ε} < X^{2−σ}) | l.159–160 |
| Thm 3 step 2′, σ > 9/10 | F* = F₂ = [R] Prop 1 | l.161 |
| Thm 3 step 3 | [C] Thm 1(b) (trivial count), not [R2] Prop 1 | l.164 |
| Thm 3⁺ | F₂ from [R2] Prop 2; trivial count h⁺ ≤ 5D^{1/2}log²(3D) per window | l.213–214 |
| Thm 3^∠ | [R] Prop 1′ (CM forms of weight \|ℓ\|+1, exponent F₂ + O(κ+ε₂)) | l.219–221 |

[R] Prop 1′ re-checked: V·I₂′ ≍ (X²M⁴)(T³)(X^{3/2}M)·X^{3/2}MT² = X⁵M⁶T⁵ ([R] l.167–169); M ≤ 6X^{κ+ε₂}.

## F4. Humphries–Thorner (I3) in the TeX (arXiv 2103.05634, e-print of 13 Mar 2024 = Compositio Math. 160 (2024) 1041–1072)

- thm:ZDE (published Theorem 1.2): Σ_{π∈𝒮(Q)}N_{π×π′}(σ,T) ≪ (|𝒮(Q)|⁴(C(π′)QT^{[F:Q]})^{6.15max{n²,n′n}})^{1−σ+ε}, for arbitrary 𝒮 ⊆ 𝔉_n. [R] does not use this; it uses the intermediate inequality below.
- thm:main_theorem_1 / prop:pre_large_sieve: Σ_{π∈𝒮(Q)}|Σ_{N𝔫∈(x,xe^{1/T}]}μ_{π×π′}(𝔫)b(𝔫)|² ≪ Q^ε(x/T + Q^{4n²θ_n+n}T^{n²[F:Q]/2+ε}|𝒮(Q)|)Σλ_{π′×π̃′}|b|², "θ_n … the best known exponent towards GRC for π ∈ 𝒮". θ_n enters only through (LRS_2) |α_{π×π′}(𝔭)| ≤ N𝔭^{θ_n+θ_{n′}}, Re μ ≥ −θ_n−θ_{n′}, for pairs in 𝒮 × 𝒮, i.e. in the ramified Euler factors on Re s = 1/(4log(ex)) in (cahrs) and in (preconvex). If every π ∈ 𝒮 satisfies GRC at all places, these hold with θ = 0: θ₂ = 0 is legitimate, as [R] says.
- cor:MVT (published Cor 4.4): if Y ≥ e, X ≥ Q^{4n²θ_n+n+ε}T^{n²[F:Q]/2+1+ε}|𝒮(Q)| and log Y ≍ log X, then Σ_π∫_{−T}^{T}|Σ_{N𝔫>X}μ(𝔫)N𝔫^{−1−1/log Y−iv}|²dv ≪ C(π′)^εQ^ε log X and Σ_π∫|Σ_{N𝔫≤X}μ(𝔫)N𝔫^{−1/2−iv}|²dv ≪ C(π′)^εQ^εX log X. At n = 2, F = Q, π′ = 1, θ = 0: X ≥ Q^{2+ε}T^{3+ε}|𝒮| = [R]'s V.
- §5, eqn:bound_int: Σ_π N(σ,T) ≪_ε X^ε(Y^{2(1−σ)}∫Σ_π|1 − LM_X(1+1/log Y+iv)|²dv + Y^{1/2−σ}∫Σ_π|LM_X(½+iv)|dv + Y^{1−σ}). Its derivation (zero_detect_1–3, boxes of height 2(log Y)², T ≥ 2(log Y)²) uses Y only through log Y ≍ log X; the κ term vanishes for π′ = 1 (cuspidal L(s,π) has no pole). HT then apply Cauchy–Schwarz to the second integral, bound I₂ := Σ_π∫|L^{ur}(½+iv,π)|²dv by preconvexity and choose Y; [R] keeps I₂ and chooses its own Y. With cor:MVT this is exactly [R]'s (I3): V^ε[Y^{2(1−σ)} + Y^{1/2−σ}(V·I₂)^{1/2} + Y^{1−σ}].
- Hypotheses against use. 𝒮: any finite set of distinct cuspidal π on GL₂(A_Q) ([R]: π_χ, multiplicity ≤ 2 handled by the factor 2); GRC: dihedral weight 1 (D < 0), Maass ¼ (D > 0), CM weight ≥ 3 (sectors), all with unimodular Satake parameters and tempered ∞-type; Q ≥ max C(π): 12|D|, 16D, 12X(M+4)²; log Y ≍ log V: 1 ≤ 3 − 2σ ≤ 3; T_X = 2X^{κ+ε₁} ≥ 2(log Y)²; implied constants depend on n, n′, [F:Q], ε only, so uniform in D, χ, T. Cosmetic: HT's boxes reach |v| ≤ T + (log Y)², written ∫_{−T}^{T}; use 2T.
- Verdict: (I3) is quoted correctly, and [R]'s use is within its hypotheses.

## F5. [R4] as used by [R5] §6 and the indefinite-forms corollary

- Lemma 1⁺: (γ) ↦ e(ku_γ/R⁺) is well defined on P⁺_{K,Z}(f) (good generators differ by η^Z, u shifts by R⁺Z), multiplicative, extends to I_K(f) (S¹ divisible); λ_k((m)) = 1, λ_k(𝔞̄) = λ̄_k(𝔞) since (N𝔞) has the generator N𝔞 > 0. Correct.
- Lemma 2⁺(i): (T_p^{(k)})² = T_{p²}^{(k)} + (Δ/p) uses λ_k((p)) = 1 and [(p)] = 1 in G⁺ (p ≫ 0); [1](T_M^{(k)}) = 1[(M,f)=1]Σ_{γ mod η^Z}e(ku_γ/R⁺); Poisson Σ_kŵ(k/R)e(ku/R) = RΣ_j w(u−jR) re-derived; hence Σ_kŵ(k/R⁺)Σ_C|S_C^{(k)}|² = R⁺Q_Δ(b), Q_Δ ≥ 0 (ŵ = ν̂² ≥ 0), and the k = 0 term gives Σ_χ|Σbλ_χ|² ≤ ŵ(0)^{−1}h⁺R⁺Q_Δ. Correct; [R4] §5 tests it to 8.9·10⁻¹⁶.
- Lemma 2⁺(ii): Dirichlet's formula h(O)R(O) = h(O_K)R(O_K)fΠ(1 − (d_K/p)/p), h⁺R⁺ = 2hR in both cases (unit of norm −1 or not). Correct.
- Lemma P⁺ ([R5] l.175–182): steps (a), (b) of Lemma P carry over (T² ≡ 4M mod p²; γ/e ≫ 0, u_{γ/e} = u_γ), each block is Lemma 2⁺(i) for the non-square order Δ_e. Correct.
- Theorem 1⁺ at w = 1 ([R5] l.184–191): Lemma 5⁺ (|u| ≤ 1 ⇒ T ≤ 2√M cosh 1, y² ≤ 4(T²−4M)/X′ ≤ 4·4M sinh²1/X′ = 22.1M/X′ ≤ 23M/X′); square Δ = a² give ≪ M^ε triples each, total ≪ X^δN‖b‖² ≤ X^{1+δ}‖b‖²; Lemma 6⁺ re-derived: t = 2√M(1+s)^{1/2}, arccosh√(1+s) = arcsinh√s, ∫f = √M·v·K⁺(v), factor 2 from ±y with T > 0 only; W smooth because w(x) = w₀(x²) and arccosh² is analytic at 1, and f vanishes near t = 2√M (Φ = 0 there); |W^{(i)}| ≪ (2√M)^{−i} ≤ (2λ₁)^i as y²X′ ≤ 23M < 32M; error ≪ √M(96c√M/(2π|k|X′))^i. Hypothesis X ≥ 128c√MX^δ holds as c√M ≤ 4N ≤ 4X^{1−4δ}. Correct.
- Theorem 2⁺: [R3] §2, §4 with Ĝ⁺; L_p(s,χ)^{−1} unchanged (χ((p)) = 1); Lemma 3(c) via [D2] Lemma 5; class II via [D2] (2.1) for Ĝ⁺ (stated on Re s = ½, used at ½+δ: same proof, routine, as for D < 0). Correct as used with Y = X^{1−5δ}.
- Theorem 3⁺: the dictionary of [R2] Theorem 4 (narrow genus theory: Ĝ⁺₂ = {χ_{d,D/d}}, L = L(s,χ_d)L(s,χ_{D/d}); (I6) for Hecke characters of K; E₂ with d | D; Σ|R| ≤ C(h⁺ + 2^{ω(D)})(H/x + x^{−1/2})log(Dx)); [R5]'s statement keeps H ∈ [xD^{−κ}, x], as the [R4] review asked (the [R4] statement says only H ≥ xD^{−κ}).
- Corollary 5 ([R2] l.147–156): N(ua + v(−b+√D)/2) = a(au² − buv + cv²) = aQ(u,−v) re-derived; α ≫ 0 ⇒ Q(u,−v) = N𝔭 = p > 0; ≤ 2 prime ideals per p. x² − ny²: D = 4n fundamental for squarefree n ≡ 2, 3 (4), principal narrow class. Correct.
- Sectors ([R4] §2, [R5] Lemma P^∠, Thms 1^∠–3^∠): ψ_ℓ well defined for |Δ| > 4 and ℓ even; the block value [1](S^∠S^∠*) depends only on ψ_ℓ on principal ideals (independent of the extension); Δ_e ∈ {−3, −4}: Σ_{u∈O^×}u^ℓ = w_e1[w_e | ℓ]; Theorem 1^∠: 64c(1+2ℓ²)√MX^δ ≤ X for |ℓ| ≤ X^δ; Theorem 3^∠ step 3: (2M+1)h zeros per window, 4Aδ + 2κ + ε₁ + ε₂ ≤ 5.32δ ≤ 6δ; Markov over ≤ 2M_X+1 values of ℓ costs X^{κ+ε₂} ≤ X^δ. Correct.
- Not re-refereed here: [R4] §3–§4 (a_E signs, sharpening; not used), and [D2] Lemmas 4 (M > 0) and 5 beyond the argument check above.

## N3. D > 0 numerical test (`fampos.gp`, PARI 2.19; narrow class groups bnrinit(bnf, [1, [1,1]]))

Family 𝔇⁺(X) (fundamental D ∈ (X/2, X]), classes from bnrisprincipal of one prime above each split or ramified p ≥ 19 ([𝔭̄] = −[𝔭] since 𝔭𝔭̄ = (p), p ≫ 0); LHS = Σ_D h⁺Σ_C|S_C|², same detector b as `fam.c` (σ = 5/6, V = √X/2, Y = X), admissible n. R_fam := LHS/(Σ_D h⁺·Σ|b|²) (actual family size), R_X := LHS/(X^{3/2}Σ|b|²) (Theorem 1⁺ normalisation); Σ_D h⁺R⁺ = 5420, 16109, 46207, 131782 for X = 1000, 2000, 4000, 8000 (0.17–0.18 X^{3/2}).

| X | Σ_D h⁺ (Σ_D h⁺R⁺) | N | R_fam killed | R_fam pole | R_X killed | R_X pole | χ = 1 share, killed |
|---|---|---|---|---|---|---|---|
| 1000 | 454 (5420) | 128 / 256 / 512 | 1.13 / 1.44 / 1.41 | 15.3 / 24.1 / 39.4 | 0.016 / 0.021 / 0.020 | 0.22 / 0.35 / 0.57 | 29–50 % |
| 2000 | 1033 (16109) | 128 / 512 / 1024 | 1.34 / 1.52 / 1.26 | 15.1 / 41.4 / 69.3 | 0.015 / 0.018 / 0.015 | 0.17 / 0.48 / 0.80 | 29–51 % |
| 4000 | 2268 (46207) | 256 / 1024 / 2048 | 1.15 / 1.77 / 1.33 | 22.6 / 60.5 / 110.9 | 0.010 / 0.016 / 0.012 | 0.20 / 0.54 / 0.99 | 23–56 % |
| 8000 | 5290 (131782) | 256 / 1024 / 4096 | 1.37 / 2.10 / 1.23 | 19.7 / 48.2 / 179.3 | 0.010 / 0.016 / 0.009 | 0.15 / 0.36 / 1.33 | 32–59 % |

- Theorem 1⁺'s normalisation: the pole-killed R_X stays at 0.009–0.021 for N from 2.9√X to X/2 (at N = X/8: 0.016, 0.017, 0.016, 0.016 for X = 1000 … 8000; at N = X/2 it falls, 0.020 → 0.009). The pole-present R_X grows almost linearly in N (0.15 → 1.33, ×9 for N ×16, at X = 8000) and with X at fixed N/X = ½ (0.57, 0.80, 0.99, 1.33, slope ≈ 0.41).
- Family-size normalisation: the pole-killed sum is 1.1–2.1 times Σh⁺·Σ|b|² for N ≥ 2.9√X, the family size of the narrow class characters, far below Theorem 1⁺'s X^{3/2} (Σ_D h⁺ = 0.45–0.66X here). Its slow rise at N = X/8 (1.13 → 2.10) is the trivial character, about half of the sum. With a mollifier of length V = √X/2 that term can reach X·V·Σ|b|² ≍ X^{3/2}Σ|b|², which Theorem 1⁺ allows. At N ≤ 2√X it is 3–8 (χ = 1 carries 75–90 %, as for D < 0).
- Memory 24 MB, 241 s for X = 8000; gp under a perl fork-and-kill timer.
