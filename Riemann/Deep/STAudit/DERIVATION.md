# STAudit — derivations (2026-11-12)

Line numbers: [S] = [STLemP NOTE](../STLemP/NOTE.md), [R3] = [RCApps3 NOTE](../RCApps3/NOTE.md), [CS] = [ChebST NOTE](../ChebST/NOTE.md), [RS] = [RSFam NOTE](../RSFam/NOTE.md). Exponents are powers of X; X^{O(δ)} dropped where marked. Scripts: `.scratch/staudit/` (copied to [checks/](checks/)).

## D1. Theorem 1′ ([S] 72–91) is [R5] Theorem 1 with Lemma 9 replaced by (B2)
Where [R3] §3.4–§3.8 (at w = 1, the RC5Audit-refereed form) touches b:
- Lemma 6 (Poisson in T, [R3] 187–205): needs (c, M) = 1, M ≤ X′², X′ ≥ 64c√M X′^δ. With n = cn₂ squarefree, c√M = √(nn′) ≤ 4N ≤ 4X′^{1−4δ}, so 256NX′^δ ≤ X′ once X′^{3δ} ≥ 256. Arbitrary admissible b.
- Lemma 7: M = n₂n₂′ odd, cube-free: squarefree support. Arbitrary b.
- §3.5: S(2^ate√v) = 0 unless 2^ate ≤ E₀ = 23N/(c√X′) (√v = √X′/(2√M) ≥ c√X′/(8N), S(u) = 0 for u ≥ √8). Uses only supp b ⊂ [N/16, 4N].
- y = 0 (Lemma 5(a)): ≤ 4X′Σ|β|². Arbitrary b.
- §3.6 ([R3] 255–270): Heath-Brown's sieve with arbitrary coefficients γ_m and A_κ = Σ_{κ|n}d(n)|β(n)|²; needs E₀ + L/κ ≤ 2L/κ ⟸ E₀² ≤ L ⟺ N ≤ 4cX′/529. Then ≪ X′^{8δ}X′Σd⁴|b|², and (B1) gives Σd⁴|b|² ≪ X′^{δ}N^{1−2σ} (N ≤ X′).
- §3.7 ([R3] 272–286): b enters only through B_{κ,q}(1) = ±Σ_{(m,κ)=1}b(cκm)ψ(κm)(m/q) and B′_{κ′,q}(1) (with b̄). (B2) with cκ for κ gives |B| ≤ C_bX′^{3δ}N^{1/2−σ}(cκ)^{−1/2}d(cκ)^{C_b}𝔙^{1/2}𝓛_b(q,ψ). b̄: ψ mod 8 and (·/q) are real, so |Σb̄(κm)ψ̃(m)| = |Σb(κm)ψ̃(m)|: the same bound. The weights (2^ate)^{−1}, |ϑ| ≤ 2^{ω(T₄)}T₄T₅, (κκ′)^{−1/2} = (T₂T₃)^{−1/2}(T₄T₅)^{−1} and Σ_{T₁T₅=q}T₁^{−1}T₅^{−2} = q^{−1}σ_{−1}(q) are as in [R3] l.278–281; d(·)^{C_b} costs (log X′)^{C}. The q-sum: Σ_{q≤E₀}q^{−1}σ_{−1}(q)𝓛_b² ≤ Σ_{dyadic Q≤E₀}Q^{−1}·C X′^{3δ}Q^{1+θ}·C log log ≪ X′^{4δ}E^θ. Hence |𝔐_c| ≪ X′^{13δ}N^{1−2σ}𝔙E^θ/c, main part of Q_c ≪ c^{−2}√X′𝔙E^θN^{1−2σ}.
- Sum over c ≤ 4N: Σ|Q_c| ≪ X′^{15δ}(X′ + √X′𝔙E^θ)N^{1−2σ}. ✔ (reviewer run A: CORRECT, ROUTINE.)
- For [R3]'s detectors 𝔙 = V, θ = 0, 𝓛 = Jutila: [R5] Theorem 1. For arbitrary b ≥ 0, (B2) at κ = q = 1 fails (pole), as the NOTE says (l.67–70); RC5Audit measured the resulting X·N growth.

## D2. Lemma L ([S] 95–107)
G_Δ = I_K(p₀)/P_{K,Z}(p₀) → G_D = I_K(p₀)/P_K(p₀) surjective; χ ↦ χ∘π injective Ĝ(D) → Ĝ_Δ. For (n, p₀) = 1, invertible O_Δ-ideals of norm n ↔ O_K-ideals of norm n, so λ_{χ∘π}(n) = λ_χ(n). Adding the other ψ ∈ Ĝ_Δ: Σ_{χ∈Ĝ(D)}|Σbλ_χ|² ≤ h(Δ)Σ_{C∈G_Δ}|S_C^{(Δ)}(b)|². Lemma P (refereed) for f(Δ) = p₀: Q′_Δ(b) = Σ_C|S_C^{(Δ)}(b)|² + (w_D/2)Σ_C|S_C^{(D)}(b_{p₀})|², b_{p₀} ≡ 0 as p₀ ∤ n on supp b. D ↦ p₀²D injective, p₀²D ≡ D (8), Φ(|p₀²D|/X′) = Φ(|D|/X) = 1; all other terms of 𝒬′_{X′} are ≥ 0 (Lemma P). h(p₀²D) = h(D)(p₀ − (D/p₀)) (D < −4; also p₀ | D) ≤ 2CX′^{1/2}log X′. ✔ Corollary L: X′ = p₀²X ∈ [N^{1+5δ}, 4N^{1+5δ}], (1+5δ)(1−4δ) ≥ 1. ✔

## D3. Lemma H, Lemma W ([S] 134–184)
- Lemma H: χ := ψ(·/q) differs from its primitive χ* only at 2 ((·/q) primitive, q odd squarefree); π⊗χ* cuspidal (Sym²f cuspidal, f non-CM); AFE with |ε(s)| ≍ (𝔮(|u|+3)^d)^{−δ}; Mellin separation of V uniform in q; odd n = n₁n₂²: χ*(n) = ψ*(n₁)(n₁/q)1[(n₂,q)=1]; Cauchy–Schwarz in n₂; Heath-Brown per n₂. Gives (Q(|u|+3))^{d/2+δ}: θ = 0 (m = 1), θ = ½ (m = 2). Standard. ✔ (S₂ there has the same X^{O(δ)} bookkeeping as in D4; absorbed.)
- Lemma W(i): a_V(κm) splits over k₁ | κ, k₂ | m (κm squarefree); λ_m(κm) = λ_m(κ)λ_m(m); Σ_{m′}λ_m(m′)ψ̃(m′)m′^{−z} = L(z, Sym^m f⊗ψ̃)·(finite and convergent Euler products); h_p = (1+λ_m(p)ψ̃(p)x)Π_i(1−α_iψ̃(p)x) = 1 + O(x²) since λ_m(p) = Σα_i; |h_p−1| ≤ 2·4^{m+1}p^{−1−2δ} ≤ ½ for p > P_m = 64(m+1)². L(z, Sym^m f⊗ψ̃) entire for m ≤ 2 and every ψ̃, the principal one included (f non-CM). Shift to Re(s+u) = ½+δ: no pole; |κ^{−s}(N/κ)^u| = N^{½+δ−σ}κ^{−½−δ}. Σ_q𝓛² by Cauchy–Schwarz in r and Lemma H. ✔
- Lemma W(ii): s₀ < 2P: N = 2^{⌈log₂(P/s₀)⌉} ∈ [P/s₀, 2P/s₀) (= 1 if P < s₀ < 2P), supp ⊂ (P/(4s₀), 4P/s₀) ⊂ [N/8, 4N]; line Re u = ½+δ−σ ∈ [0, 2δ]. s₀ ≥ 2P: N = 1, b = φ₂(s₀/P)1_{k=1}, κ = 1 only. ✔ (the 2026-10-22 repair is correct).

## D4. Theorem MS ([S] 188–220)
- p₀ ∈ [X^{(m−1)/4+6δ}, 2X^{(m−1)/4+6δ}], X′ = p₀²X ∈ [X^{(m+1)/2+12δ}, 4X^{(m+1)/2+12δ}]. AFE length √q_m(|v|+3)^{m+1} ≤ CX^{(m+1)/2+3δ} (q_m ≍ |D|^{m+1}, independent of χ; genus χ give L(s,Π_m⊗χ_{d₁})L(s,Π_m⊗χ_{d₂}), entire); cut at P₀ = X^{(m+1)/2+4δ}.
- n = s·k (k := Π_{p∥n, p∉𝔅𝔓}p): unique, coprime; a_χ multiplicative; a_χ(k) = λ_m(k)λ_χ(k) for (m,𝒫)-admissible k (split, inert, ramified primes checked). 𝔖 = Σ_{s∈𝒮,s≤4P}d_{2m+2}(s)s^{−½} ≤ C(f,m)(log X)^{C} (squarefull part: Π_{p≤4P}(1+C/p)).
- Blocks: s ≤ 4P ≤ 8P₀ ≤ X′, |v+τ| ≤ 2X^δ ≤ X′^δ, N ≤ max{1, 2P/s} ≤ 4X^{(m+1)/2+4δ} ≤ X′^{1−4δ} (m = 1: X^{1+8δ−48δ²}; m = 2: X^{3/2+6δ−48δ²}). Lemma L + Theorem 1′ with 𝔙 = 1, θ = (m−1)/2, E ≤ 23X′^{1/2}: each block's family sum ≤ CX′^{1/2+21δ}·X′. |τ| > X^δ: per-D sieve × (1+|τ|)^{−A+2}, A = 20/δ: ≪ X^{−5}.
- **Finding (bookkeeping, ROUTINE).** "S₂ is treated identically" (l.217–220) omits that the outer weight becomes d_{2m+2}(s)s^{−½+δ}. Repair: Cauchy–Schwarz with w₂(s) := d_{2m+2}(s)s^{−½+δ}, 𝔖₂ := Σw₂ ≤ (4P)^δ𝔖; inner family sums ≤ CX′^{3/2+21δ}N^{2δ} with N ≤ 4P/s (N = 1 blocks included); w₂(s)N^{2δ} ≤ d_{2m+2}(s)s^{−½}(4P)^{2δ}. So Σ_{D,χ}|B^{(2)}|² ≤ C(4P)^{3δ}(log X)^{2C}X′^{3/2+21δ}, and (8P₀)^{3δ}|ε_m|² ≤ CX^{(3/2)(m+1)δ+12δ²−(2m+2)δ} ≤ C for δ ≤ (m+1)/24. The S₂ bound equals the S₁ bound; Theorem MS's X^{3(m+1)/4+60δ} is unchanged. (The NOTE's own count, "N^{4δ}", already reaches X^{(2m+2)δ}: it over-counts the inner factor by what it omits outside.) Reviewer run B: "GAP — SIGNIFICANT … the claimed exponent survives"; reviewer run C on this repair: "GAP — ROUTINE … MS's exponent is unchanged".
- Exponent: X′^{3/2+21δ} ≤ 16X^{3(m+1)/4+50δ}. m = 1: X^{3/2+O(δ)} = family size ✔; m = 2: X^{9/4} ✔. Critical line (σ = ½ in both Lemma W(ii) sums, |ε_m| = 1): ≪_ε X^{3(m+1)/4+ε} ✔.

## D5. Theorem 2₁ ([S] 238–246), Proposition F₂ ([S] 252–268)
- Class II: Σ_D|𝒳_II| ≤ CY^{½+δ−σ}∫e^{−|v−t|}(ΣΣ|L_m|²)^{½}(ΣΣ|M̃|²)^{½}dv; ΣΣ|M̃|² ≤ CX^{3/2}(log X)^C ([CS] F1(c), V ≤ X^{1/2}/2, per-D); Theorem MS for |v| ≤ X^δ, per-D F1(d) beyond. m = 1: X^{3/2+31δ}Y^{½+δ−σ}, Y = X^{1−5δ}: ≤ X^{2−σ+35δ}. m = 2: X^{15/8+31δ}Y^{½+δ−σ}.
- Class I, m = 1: V = X^{1/2}/2, 2Z ≤ X^{1−4δ}, Lemma W(i) (𝔙 = V, θ = 0), p₀ = 1: X^{3/2+21δ}N^{1−2σ}; Σ_{N>V/2}: X^{2−σ+23δ}. ✔ Same structure as the refereed [R5] Theorem 2.
- Class I, m = 2: N ≤ X^{1−4δ}: ½ + max{1, ¼ + v + ν/2} + ν(2η−1), max over [v, 1] is R(v) = 3/2 − v(1−2η) (second branch ≤ R ⟺ v ≤ ½); X^{1−4δ} < N < 2Z: lift, X′ ≍ N^{1+5δ}, E ≤ 23N^{(1−5δ)/2}: (N^{3/2} + VN^{5/4})N^{2η−1}, exponent max{ν/2, v+ν/4} + 2νη ≤ y/2 + 2yη (v = y/4). Balance c₂ = 15/8 − y(½−η) = y(½+2η): y = 15/(8(1+η)), F₂ = 15(1+4η)/(16(1+η)); R(y/4) ≤ F₂ ⟺ η ≥ 1/14. 𝒫: ≤ 2log Z primes ≥ X^{δ/4}, Σp^{−1/2} ≤ 1. ✔ (`thr_check.py`, exact.)
- Independent optimisation (`thr_opt2.py`, free y, v, same sieves): E*(η) = F₂(η) to 4 digits on [1/14, 0.24]; threshold max_η (E*−1)/η = 2.0365 at η = 0.1752.

## D6. Thresholds ([S] 276–348)
- T¹ (A > 1): E₁ by Theorem 2₁ for σ_j ≤ 9/10 (exponent 1 + 2ε′ + 41δ − (A−1)(1−σ_j) ≤ 1 − ε′) and [RS] Prop F (X^{11e(σ)}, 11e ≤ 11/12) above; (b) below σ₀ = ½+4δ: H|D|^{½−A(½−4δ)+κ+ε₁+o(1)} ≤ H|D|^{−(A−1)/2+6δ}; (f), non-degree-one ideals need x ≥ |D|^{1+3δ}. ✔
- T² (A > A₂*): f(η) = (F₂−1)/η = (44η−1)/(16η(1+η)), max at η* = (1+3√5)/44 ∈ [1/14, 1/(2A₂*)], f(η*) = (23−3√5)/8 = 2.0364745; σ₀ = 1 − 1/(2A₂*) = 0.75448 ≥ ¾; (13/14, 1): 12e ≤ ¾. ✔ (exact)
- T₃ (A > 3, PT): class II X^{9/4}Y^{½−σ} with Y = X^{3/(1+2η)} is X^{¾+6e(σ)}; class I max{1+3η, ¾ + 6η/(1+2η)}; at η = 1/(2A): ¼ + 3/(A+1) < 1 ⟺ A > 3. ✔ 𝒜(4..8) = 11/3: class II needs y > 2A/(A−1) > 11/4 for A < 11/3, and every sieve at N = X^{11/4} costs ≥ ¾ + 5.5η. ✔ 𝒜(M ≥ 8) = (M+3)/3 ([CS] route H′, (m+6)/(A+1) ≤ 3). ✔
- Theorem D′/Cor D′1: degree-2 minorant of 1_{x>0} on [−2,2]: x(x+2)/8 (mean ⅛); degree-1 impossible (P(0) ≤ 0 forces c₀ ≤ 0); 1_{x>−α}: (α + x)/(2+α). ✔
- Discrepancy: n = ⌈(3A−3)/2⌉ − 1 gives 5/(n+1) ≤ 10/(3A−3). ✔ A(J) ≤ max{11/3, 25.2/L + 1.7} from M*(J) ≤ 75.4/L + 2 ✔; numerically (`minor.py`, 4812 intervals) the least n with μ(J) > λ_n(α) + λ_n(β) has nL ≤ 10.0 (≈ 3π at an end), so the stated bound holds with room (≈ 3.8×).

## D7. Imported statements, read in the TeX
- PT 2508.14888v2: thm:large_sieve2 (l.288–296) C(N,𝒮) ≪ (NQ)^{o(1)}(N + Q^n|𝒮|); cor:large_sieve2 (l.334–344) C^∞(N,𝒮) ≪ (NQ)^{o(1)}(N|𝒮|^{(k−1)/k} + Q^{n/k}|𝒮|), C^∞ normalised by N^{1/2}max|a|. n = 2, k = 2, |𝒮| ≍ X^{3/2}, Q ≤ 12X: (LS3) NX^{3/4} + X^{5/2} ✔, (LS2) N + X^{7/2} ✔; no GRC, no coprimality. Preprint.
- HT 2103.05634v2: thm:main_theorem_1 (l.717–730): (x/T + Q^{4n²θ_n+n}T^{n²/2+ε}|𝒮(Q)|), gcd(𝔫, 𝔮_π𝔮_π′) = 1; θ₂ = 0 for θ_χ (GRC); eqn:bound_int (l.949–954) with cor:MVT: Y^{2−2σ}·(…) + Y^{½−σ}(I₂·X_mol)^{½}: [RS] Prop F's (V·I₂)^{2e(σ)} ✔.
- BTZ 1804.06402v3: thm:lfzde_RS (l.334–339) (C(π₀)QT)^{10⁷(m₀m)⁴(1−σ)} under hyp (GRC on average) ✔ c_m = 10⁷(m+1)⁴; thm:Landau_Siegel (l.373–379) factor min{1, (1−β_χ)log QT} ✔ (ν^♯).
- Thorner 2404.06482v2 eq:SymL (l.600–604): L(σ, Sym^m(π)⊗ξ) ≠ 0 for σ ≥ 1 − c/(m^{2+y_π}log(mC(π)) + m² log C(ξ)), every ξ (twist |·|^{it} gives all heights) ✔ (J5).
- [R5] Lemma P, Corollary P, Theorem 1 (w = 1), Lemma P⁺, Theorem 1⁺; [R4] Lemma 2⁺: refereed AUDIT PASSED (RC5Audit REPORT, REPORT2).

## D8. D > 0 ([S] 120–124, 350–356)
Lemma L⁺: G⁺_{p₀²D} → G⁺_D surjective; h⁺R⁺(p₀²D) = h⁺R⁺(D)(p₀ − (D/p₀)); [R4] Lemma 2⁺ and [R5] Lemma P⁺ (refereed) make every non-square order nonnegative; Theorem 1⁺'s main term is [R3] §3.7, so (B2) applies verbatim; square Δ ≪ X′^δN‖b‖². Theorem MS for D > 0: two archimedean types treated separately (the gamma factor of Π_m × π_χ depends on it). Not re-derived beyond this reduction. Sectors not claimed.
