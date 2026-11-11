# AUDIT — Lamzouri, arXiv 2610.26841v2, Theorem 1.1 (isolated zeros)

Source: IsolatedZeros2.tex of v2 (26 Oct 2026; v2 is the latest version on arxiv.org/abs/2610.26841 on 6 Nov 2026). Line numbers below are TeX lines.
Verdict: NO FLAW. The chain is correct and unconditional. One cited input (Lem:Conrey, the constant 0.0806) has no published proof; it follows from a published theorem of Conrey by a short argument (DERIVATION §3), which also improves it. The O(τ²N) term of Pro:IsolatedSMALL is correct but not explicit; DERIVATION §2 makes it explicit.

## Chain

| step | statement (TeX lines) | inputs | status | needs RH? |
|---|---|---|---|---|
| Thm 1.1 (1.5) | ∃ c > 0: liminf B_0(T,c)/N(T) ≥ 0.6552 (l.211–216) | Pro:IsolatedSMALL, Lem:ConstructionKernel, Lem:Conrey, α_0 | correct; c not explicit | no |
| Thm 1.1 (1.6) | liminf B_0(T,1/6)/N(T) ≥ 0.50875 (l.218–220, proof l.704–736) | Pro:IsolatedZeros with K_ε(2τ_0) ≥ 5/6, Lem:ConstructionKernel | correct (τ_0 = 0.169641840… > 1/6; K_ε(2τ_0) ≥ 5/6 argument l.716–725 checked) | no |
| Pro:IsolatedZeros | #isolated ≥ 5/2 N − 3/2 ΣK² if κ_0 ≥ 5/6 (l.275–468) | Hilbert space: Bessel, Parseval, Gram–Schmidt with real coefficients | correct; every displayed step re-derived (inner products real l.371–374; y² ≥ 4y−4, y² ≥ 2y−1; dimension count l.366; pair/triple norms l.428–436; κ_1 ≥ 2κ_0² − 1 l.454–457; final linear algebra l.461–466) | no |
| Pro:IsolatedSMALL | #isolated ≥ 2N − ΣK² − Σ_ℓ max(1 − 4K(iy_ℓ)²/(1+K(2iy_ℓ)), 0) + O(τ²N) (l.511–573) | same as Pro:IsolatedZeros plus S_0 | correct; O(τ²N) made explicit: cluster loss λ_C ≈ 12.2τ² per clustered point and θ̄ = π²τ²/2 in the witness bound (DERIVATION §2.2–2.4); its proof drops 2b ≥ 0 at l.524–528 (harmless) | no |
| Lem:ConstructionKernel (= La 2609.02882, Lemma 3.2) | ΣK_ε(i(ρ−ρ')L)² = (C_{η_ε} + O_ε(1/√log T)) (T/2π) log T, |C_{η_ε} − C_MT| < ε (l.670–700) | BGST, Acta Arith. 214 (2024) Lemma 5: unconditional, all zeros with multiplicity | La is a 2026 preprint; BGST refereed; our ZetaZeros audit re-derived C(f) = ∫f² + ∬|u−v|ff from BGST. C(f_0) = C_MT checked here in closed form and by Arb quadrature | no |
| Lem:Conrey | Σ_{β>1/2, γ≤T}(β − ½) ≤ (0.0806 + o(1))T (l.757–767) | Conrey, Bull. AMS 20 (1989) 79–81, Theorem 3 — an announcement; the paper says "Conrey did not include a proof" (l.765) and sketches it via BCHB | GAP in the literature, not a flaw: 0.0806 = (1/4π) log(1 + 1/θ) at θ = 4/7, rounded up (0.0805006…); this follows from Conrey's published Theorem 2 (Crelle 399 (1989)) via the functional equation and Littlewood's lemma (DERIVATION §3, R → 0+). BCHB's published 0.0845 = (1/4π) log(26/9) is the same argument with θ = 9/17 | no |
| (Eq:BoundSumOffLineZeros) | Σ_{β>1/2} max(H_ε, 0) ≤ α_ε L Σ(β − ½) (l.797–801) | max(H,0) ≤ α y for y > 0, sum over distinct pairs ≤ sum with multiplicity | correct | no |
| α_0 (l.802–812) | α_0 = max_{y>0} H_0(y)/y < 0.21462, "using Mathematica" | H_0 < 0 on (0, 1/10], H_0(y)/y < 1/5 < H_0(2)/2 for y ≥ 5 | certified: α_0 ∈ [0.214612853119, 0.214612853122] (Arb, Taylor-model branch and bound), H_0(2) = 0.40925570… > 0.4092 | no |
| l.816 | 2 − C_MT − 0.0806 α_0 > 0.655202 | C_MT = 1.3274992963…, α_0 | certified ≥ 0.655202908 (margin 9.1e-7) | no |

Which inputs need RH or Montgomery's conjecture: none. The pair-correlation input is BGST's unconditional theorem over all zeros (on and off the line, with multiplicity); the off-line zeros enter only through the conjugate-pair structure of Z_T and through Lem:Conrey.

## Points checked in detail

1. τ–c relation. Z_T elements are γL (L = log T/2π) for zeros on the line; |x − Re z| > τ ⇔ |γ − γ'| > 2πτ/log T. Pro:IsolatedSMALL is applied with τ = c (l.776). The counted set I_1 (isolated simple real points with no other point within τ in real part) gives c-isolated zeros except possibly the top zero below T and n = 1: B_0(T,c) ≥ |I_1| − 2.
2. Hidden c-dependence. None besides the O(τ²N): Lem:ConstructionKernel's constant and Lem:Conrey are c-free; the kernel K_ε is chosen independently of c (l.776). The implied constant in O(τ²) depends on ∫u²η² ≤ ¼, uniformly in ε (so the order "c small, then ε small" at l.776 is harmless).
3. liminf versus finite T. All statements are liminf statements; the o(1) of Lem:ConstructionKernel is O_ε(1/√log T), not explicit. No finite-T claim is made here or by the paper.
4. Multiplicity. Zeros are counted with multiplicity throughout (l.180, l.577); isolated zeros are simple. In Pro:IsolatedSMALL the sum is over distinct conjugate pairs; (SMallGapsProof) replaces it by a sum over zeros with β > ½ with multiplicity, which is larger (correct direction).
5. log T versus log γ_n. Thm 1.1 defines B_0 with 2πc/log T; the abstract with 2πc/log γ_n. Harmless for liminf statements (DERIVATION §7.1), but the two definitions differ at fixed c.
6. The ε → 0 limits (|C_{η_ε} − C_MT| < ε, α_ε = α_0 + O(ε)) are correct: f_ε → f_0 in L¹, L², and K_ε → K_0 locally uniformly (DERIVATION §1.3, §4).
7. Dirichlet analogue (l.966–1032; Lem:GaPa, Lem:ConstructionKernelDirichlet, LemmaWu, Lem:AnalogueConreyL): read, not re-derived; Wu's 1 + 1/θ gives (1/4π)log(11/4) = 0.0805 as for ζ. Pro:BombieriHejhal (l.868–964) not audited. Wu's Theorem 1.1 (arXiv 1802.09704, (1.10)) prints (h,k)^{α+β} where its own (1.12) needs (h,k)^{1+α+β} (typo in Wu, harmless here).

## Remarks on the paper's Remark 1.2 (l.222–227)

Both suggested improvements exist, and the second is large. (a) Kernel: within cos(wx) kernels the gain is only +1e-4 at small c with the paper's inputs (DERIVATION §5). (b) Off-line term: the bottleneck is how Lem:Conrey is used, not the mollifier length. Littlewood's lemma at σ = ½ + 2πρ/log T (Conrey's published mean value moved there by the functional equation) bounds Σ_{β>1/2}((β − ½)log T/2π − ρ)^+ by B(ρ)N with B(1.2) ≤ 1.441e-5; since the deficiency 1 − 4K(iy)²/(1+K(2iy)) is negative for y < 1.22, the whole off-line term drops below 1e-5 and the proportion rises from 0.6552 to 0.672487 at c = 0.001, within 1.4e-5 of the method's simple-and-on-line constant 2 − C_MT (DERIVATION §3.7).

## Pro:BombieriHejhal and Theorem 1.4 (follow-up audit, l.824–964 and l.966–1027)

Verdict: correct as a reduction; it rests on Bombieri–Hejhal's Theorem B and Lemmas 7/9 (Duke 1995), cited, not re-derived; they need only Selberg's zero-density estimate N(σ,T,χ) ≪ T^{1−c(σ−½)} log T (l.878–880), classical for fixed χ.
| step (lines) | check |
|---|---|
| F, θ_j, Λ_F(s) = conj Λ_F(1 − s̄) (l.840–855) | e^{iθ_j}Λ(s,χ_j) = e^{−iθ_j}Λ(1−s,χ̄_j) as |ε_χ| = 1; b_j real; so G(t) = e^{iϑ(t)}F(½+it) is real (l.940–945) ✓ |
| H_j(T): gap ≤ 2π/((log T)(log log T)^{1/4}) (l.862–865) | o(average gap): for any fixed c, H_j ⊂ non-c-isolated zeros ∪ O(1) for large T ✓ |
| Lem:SelbergCLT, Lem:Separation (l.873–896) | joint CLT, independent limits for χ_1 ≠ χ_2; exceptional set E of measure ≪ δT ✓ (CLT taken from BH Thm B) |
| Lem:BHLem9 (l.898–910) | taken from BH Lemma 9 (via their Lemma 7) ✓ as cited |
| good ℓ, count of bad ℓ (l.913–927) | meas((y_ℓ, y_ℓ+√δM/log T) ∩ A_ℓ) ≥ (√δ − δ)M/log T ≥ √δM/(2 log T) for δ ≤ 1/4; bad ℓ ≪ √δ T log T/M ✓ |
| B_ℓ by Markov (l.929–935) | exceptional u-measure ≪ M√(log M)/(δ²(log log T)^{1/3}) ✓ |
| dominance (l.937–949) | δ√(log log T) − 2·O((log log T)^{1/3}) ≥ (δ/2)√(log log T) ✓ |
| zeros in I_ℓ, k_ℓ (l.950–953) | all zeros of L_1 (on and off the line, multiplicity) ✓ |
| τ-isolated zeros, τ = (log log T)^{−1/4} (l.953–956) | interior zeros only; simple, on the line; each neighbouring gap has u-length > 2π(log log T)^{−1/4}, larger than the exceptional measure, so it contains a point of B_ℓ ✓ (one point per gap keeps the intervals disjoint) |
| sign change (l.957–961) | G ≈ b_1 e^{i(ϑ+θ_1)}L_1 on B_ℓ; the Hardy-type function of L_1 changes sign at a simple on-line zero and is of constant sign between consecutive ordinates ✓ |
| count (l.961–964) | ≥ k_ℓ − 2 − #(H_{j_ℓ} ∩ I_ℓ); summing, δ = M^{−1/3}, M large ✓ |
| Thm 1.4 (l.1018–1027) | dyadic summation; |H_j| ≤ (1 − δ_j + o(1))(T/2π) log T; result δ_1 + δ_2 − 1 ✓ |
Inputs for δ_j: Lem:GaPa (Glasnik Mat., refereed), Lem:ConstructionKernelDirichlet (from La's preprint), and for the off-line term our hinge lemma for L(s,χ) from Wu's Proposition 3.1 (DERIVATION §8.2) instead of Lem:AnalogueConreyL.

F1 (added by the orchestrator from LamzAudit2, 11 Nov 2026): 2610.26841 l.622–701 quotes Lamzouri's Lemma 3.2 as |C_{η_ε} − C_MT| < ε for any admissible cutoff ψ_ε; this is false for steep cutoffs, but only an O(ε) bound is used there and here, so nothing breaks ([LamzAudit2](../LamzAudit2/REPORT.md) §3, cutoff.py).
