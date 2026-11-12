# DPC 2610.28527 — statements and inputs

Source: arXiv e-print v2 `Gaps_for_Dirichlet.tex` (1092 lines; line numbers below). v1 `main.tex` (1090 lines, 20 Oct 2026) differs only in the acknowledgements (v2 L972–974: ChatGPT used for ideas, proofreading, "to run the optimisations found in the large gaps calculations", writing; paper to be split into a small-gap and a large-gap paper). E-print = TeX only: no code, no coefficient data.

Notation: χ primitive mod Q (Q=1: ζ); γ_n positive ordinates **counted with multiplicity** (L47, L68); μ, λ = liminf/limsup of (γ_{n+1}−γ_n)log γ_n/2π (L47–50); h = 2πφ/log T; V_φ(t) = π⁻¹(log ξ_χ(½+i(t+h/2)) − log ξ_χ(½+i(t−h/2))), boundary values from Re s>½ (L163); L = log T, N = ⌊T/L²⁰⌋, Y = √T/L¹⁰ (L237); E_T = average on [T,2T] with weight |R_T|²exp(−(t−3T/2)²L²/2T²) (L242); K_{kq,q} = −2iΛ(k)k^{−1/2}(log k)^{−1}sin(πφ log k/L) (L254); A_φf(x) = iφf(x) + (2i/π)∫₀ˣ sin(πφu)/u f(x−u)du (L398).

## Main results

| # | TeX | Statement | Hypothesis |
|---|-----|-----------|------------|
| Thm 1 (A) | 53–55 | μ ≤ ½ for every primitive χ | none; L57: off-line zeros pair (ρ, 1−ρ̄) at one ordinate ⇒ zero gaps ⇒ μ = 0, else all zeros above some Y₀ on the line |
| Thm 2 (B) | 61–67 | ∀ε ∃c_ε: #{n: T≤γ_n≤γ_{n+1}≤2T, gap·log γ_n/2π < ½+ε} ≥ c_ε T log T (zero gaps counted, L69) | RH for ζ |
| Thm 3 (C) | 75–77 | λ > 3.63846723 | RH for L(s,χ) |

## Lemmas and propositions

| Item | TeX | Claim | Hypotheses | Inputs |
|------|-----|-------|------------|--------|
| Lemma 1 logarithmic-boundary | 158–196 | boundary values a.e.; Im V_φ = #zeros in (t−h/2,t+h/2); sup_ε∫_u^{u+1}\|log ξ(½+ε+it)\|^j ≪ T^{j+1} | zeros above Y₀ ≥ 2 on line | Hadamard (order 1) [IK 5.6, Titchmarsh 2.12]; O(log v) zeros per unit interval [Titchmarsh 9.2–9.4; Bennett et al.] |
| Lemma 2 arithmetic-sums | 198–229 | Σμ²(n)/n f f̄(+v) and Mertens limits, uniform in v | — | squarefree sums, Mertens |
| Prop 1 finite-moments | 262–388 | E_T[V^j] = ⟨(iφ+K/π)^j c,c⟩/‖c‖² + o(1), uniform in c ≠ 0; index sets: integers ≤ N, or reduced m/n with m,n ≤ Y | Lemma 1 hyp. | Stirling; Euler series at Re s = 3/2; Gaussian Fourier transform; \|log(kq/q′)\| ≥ 1/2B, B = T/L²⁰; Schur test ‖\|K\|‖ = O(1) |
| Prop 2 integer-limit | 404–480 | c_n = μ(n)n^{−½}f(log n/log N): lim E_T[V^j] = ⟨A_φ^j f,f⟩/‖f‖²; spec A_φ = {iφ} | Lemma 1 hyp. | Lemma 2 |
| Lemma 3 two-lines | 494–546 | B bounded, spec B ⊂ {0<Im w<1}, ‖g‖=1: real signed ν on Im w ∈ {0,1} with moments ⟨B^jg,g⟩, densities p₀(x) = e^{πx}Im⟨(e^{πx}−e^{πB})⁻¹g,g⟩, p₁ with −e^{πx}; unique among positive measures with all moments finite | — | Riesz–Dunford calculus |
| Lemma 4 two-line-limit | 548–577 | positive measure on the two lines with moments ⟨A_φ^jf,f⟩/‖f‖² | every interval of length 2πφ/log T in [T−1,2T+1] has ≤ 1 zero, all large T | Prop 2 |
| Proof Thm 1 | 580–601 | e^{πA}f = e^{ia}(f + 2ia∫₀ˣe^{ia(x−y)}f) (L585); g = e^{−iax}−e^{ia}e^{iax}; p₀(0) = sin a sin 2a/(a‖g‖²) < 0, a = πφ ∈ (π/2,π) (L596) | — | Laplace transform of sin(au)/u |
| Lemma 5 proportion-absolute-moments | 609–671 | ∫_T^{2T}\|V_φ\|^p ≪ T (p>1); prime-sum moments; ∫\|R_T\|^p ≪ T for length T^θ, primes > T^{η/2} | RH(ζ) | Rodgers CMH 92 (2017) Thm 1.1 (L629); Hilbert transform L^p (L641); explicit formula [IK 108–109, Rodgers Thm 2.1] (L643); Newton identities |
| Lemma 6 proportion-rough-limit | 686–741 | lim_{η↓0}lim_T of G_T-weighted mean of P(V_φ) = ⟨P(A_φ)f,f⟩/‖f‖² on [0,θ]; denominator ~ cT | RH(ζ) | Prop 1 averaged over centres; Dickman-type limit of the rough-number measure |
| Proof Thm 2 | 743–788 | θ ∈ (1/φ−1,1): p₀(0) = sin(aθ)sin(a(1+θ))/(a∫₀^θ\|g\|²) < 0; P = q²+ib₀ with Re P ≥ 0 on both lines, Re⟨P(A)g,g⟩ < 0 (L760); Hölder ⇒ meas E_T ≫ T ⇒ ≫ T log T pairs | RH(ζ) | Lemmas 3, 5, 6 |
| Coefficients | 795–806 | c_{m/n} = ℓ^{Ω(mn)}(mn)^{−½}H(m-data)H(n-data), m,n ≤ Y squarefree coprime; data = (log a/log Y, Σ_{p\|a}(log p/log a)^e, e = 5/4, 3/2, 2, 3) | — | — |
| Lemma 7 large-marked-law | 811–855 | limit law ℓ²v^{ℓ²−1}dv dP, P with Poisson–Dirichlet(ℓ²) partition moments (L813–820); pair weights factor | — | Selberg–Delange [Tenenbaum II.5] |
| Prop 3 large-marked-average | 882–903 | I > 0 ⇒ lim E_T[Im V_φ] = φ − (4ℓ/π)J(φ)/I | RH(L(s,χ)) | Props 1, Lemma 7 |
| Numerics | 905–966 | space v^r t^j, r ≤ 12, \|j\| ≤ 5 (1638 coeffs); J_K closed form (L922–931); remainder (L934); ℓ = 39327/20000, φ = 3.63846723 (L941); Arb: Σ\|b\| < 7.46135, 10³⁴I = 1.000000000000000020968965813, 10³⁴J₃₂ = 1.053854765588860867393555243 (±1e−27) (L943–947); margin > 2.15088826449e−9 (L953); \|J−J₃₂\| < 4.273e−44 (L957); final > 1.0811e−9 (L961–966) | — | Arb; coefficients b not published |
| Proof Thm 3 | 967–969 | mean count < 1 ⇒ empty interval; continuity in φ | RH(L(s,χ)) | Prop 3 |

## Inputs not used

Montgomery pair correlation (RH form or the all-zeros form of Baluyot–Goldston–Suriajaya–Turnage-Butterbaugh 2306.04799), zero density, Gonek's lemma, mollifier theorems beyond the diagonal evaluation in Prop 1.
