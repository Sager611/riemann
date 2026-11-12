AUDIT PASSED: all three theorems of Durkan–Pearce-Crump 2610.28527v2 check out: every decisive step was re-derived, the small-gap method provably stops exactly at the Alternative-Hypothesis value ½, and the λ > 3.63846723 margin (+2.1509e−9) was reproduced and Arb-certified with our own coefficients.

# Audit of arXiv 2610.28527 (Durkan, Pearce-Crump, 20/23 Oct 2026, unrefereed)

Files: `STATEMENTS.md` (every statement and input, with TeX lines), `DERIVATION.md` (re-derivations D1–D10), scripts `opcheck.py`, `largegap.py`, `cert.py`, `models.py`, `zfinite.py`. The TeX is v2 `Gaps_for_Dirichlet.tex`. v1 differs only in the acknowledgements, which state that ChatGPT ran the large-gap optimisation and that the paper will be split in two.

## Verdicts

| Claim | TeX | Inputs | Verdict |
|---|---|---|---|
| A: μ ≤ ½ unconditionally, every primitive χ | 53–55; proof 57, 158–601 | Hadamard product, O(log T) zeros per unit interval, Stirling, Euler series at Re s = 3/2, Mertens and squarefree sums. No pair correlation, zero density or RH-moment bounds | **Correct.** "Unconditional" only through the multiplicity convention (L47, L68): an off-line pair ρ, 1−ρ̄ is a zero gap. The real content is μ ≤ ½ whenever all high zeros lie on the line. Under RH this is the result of Bondarenko–Heap 2610.26038, by a related route |
| B: under RH, at least c_ε T log T normalised gaps < ½+ε in [T,2T] | 61–67; proof 603–788 | A's inputs, plus RH, Rodgers 2017 Thm 1.1 (tail bounds for zero counts), the L^p bound for the Hilbert transform, and the explicit formula | **Correct.** One cosmetic slip at L665: the bound should read ⌊2θ/η⌋ prime factors above T^{η/2}, which is harmless. This gives D(λ) > 0 for every λ > ½, up from Chirre–Gonçalves–de Laat's 0.6039 |
| C: under RH for L(s,χ), λ > 3.63846723 | 75–77; proof 790–969 | Prop 1 on reduced fractions, Selberg–Delange, Poisson–Dirichlet(ℓ²) moments, and a 1638-dimensional quadratic-form optimisation | **Correct.** The reduction to the finite problem was re-derived. The finite problem was then solved and certified independently of the authors' coefficients, which are not published and turn out not to be needed |

Reviewer: one reviewer subagent was given the decisive chain of A with exact TeX quotes. Verdict: "GAP — ROUTINE". Its only objection is that uniformity in c requires coefficients coprime to Q, and the TeX already requires this at L241 ("integers q≤N with (q,Q)=1, or reduced fractions m/n … (mn,Q)=1"), so there is no gap. The reviewer explicitly confirmed the Prop 1 estimates, Lemma 4 (two-line-limit) and Lemma 3 (two-lines).

## Decisive steps (details in DERIVATION.md)

- **Prop 1** (D3). Moment asymptotics hold for polynomials of length T/L²⁰, uniformly in the coefficients. The method: shift to Im t = −1, expand the Euler series, then move each oscillating exponential, together with the gamma-factor difference G but never V, to Im t = ∓H with H = T/L. Since |log(kq/q′)| ≥ L²⁰/2T, the off-diagonal terms are e^{−L¹⁹/4}-small. Every inequality was re-derived.
- **Lemma 4** (D5). If each window holds at most one zero, then N ∈ {0,1}. In that case the holomorphic moments E[V^j] determine every E[X^{2j}] by triangular induction (the top coefficient is 2j+1). This gives tightness and uniform integrability for free, which is why A needs no absolute-moment input.
- **Lemma 3 and the contradiction** (D6). Checked numerically: expm(πA) matches the closed form at L585, the moments of ν match ⟨A^jg,g⟩, and at φ = 0.6, p₀(0) = −0.105301 agrees with sin a sin 2a/(a‖g‖²).
- **Where ½ enters** (D7, our addition). Im(r−C)⁻¹ = (r−C)⁻¹(Im C)(r−C*)⁻¹, so every limit measure is positive exactly when Im e^{πA_φ} ⪰ 0. We computed ⟨Im C f,f⟩ = sin a‖F‖² + a cos a|∫F|² + a sin a⟨SF,F⟩, where S has spectrum {2/((2k+1)π)}. Hence Im C ⪰ 0 if and only if φ ≤ ½. The smallest eigenvalue of Im C is +0.462, 0, −0.476 and −1.414 at φ = 0.3, 0.5, 0.6 and 0.75. The method therefore yields μ ≤ ½ and can never yield μ < ½; it needs φ > ½ precisely through sin 2πφ < 0 in A and sin πφ(1+θ) < 0 with θ < 1 in B.

## Model tests

| Model | Shares | Outcome |
|---|---|---|
| Alternative Hypothesis (Lagarias–Rodgers), gaps in ½ℤ | Montgomery's F(α) = \|α\| on (0,1) | F is 2-periodic, so the pair weight at distance ½ is w₁ = ½ − 2/π² ≈ 0.297. About 30% of gaps equal ½, hence μ = ½ and D(½+ε) > 0, consistent with A and B. By D7 the argument cannot prove μ < ½, so it proves nothing false about AH. For C, AH would need gaps ≥ 4, which it allows |
| Picket fence Λ = γ(s)+γ(1−s) | functional equation, real on the line, all high zeros on the line, μ = λ = 1 | With no Euler product the prime operator K is 0. The moments are (iφ)^j for every weight and the law is positive. models.py gives E[V^j] = (iφ)^j for j ≤ 5 and sign Z = (−1)^N. No false claim results |
| Davenport–Heilbronn (`ExploreHelio/results/dh_zeros.json`: 9 off-line zeros up to 420) | functional equation, no Euler product | Off-line pairs share ordinates, so the reduction gives μ = 0, which is true. The analytic core is never reached |
| Actual zeta zeros, T = 2000 (zfinite.py) | everything | Prop 1 tested as a finite-T identity: weighted E_T[V^j] against Σ_r C(j,r)γ_{j−r}⟨(K/π)^r c,c⟩/‖c‖², for c = δ₁ and Möbius polynomials of length 20. Agreement is 1e−5 at j = 1, 1e−4 at j = 2 and 1e−3 at j = 3, 4, within the quadrature noise from the log singularities. A and B are asymptotic and have no finite check |

## Reproduction of 3.63846723 (D9)

Setup: ℓ = 39327/20000 and φ = 3.63846723. I = A_v ⊗ A_t, where A_t holds the PD(ℓ²) moments computed with the exponential formula (sanity check: with every exponent 1 all moments equal 1). J_K is evaluated as on TeX L922–931 with K = 40. Both forms are whitened by Cholesky at 1000 bits.

| Space | Dim | max J/I | 1−φ+(4ℓ/π)J/I |
|---|---|---|---|
| r ≤ 12, \|j\| ≤ 3 | 455 | 1.053854347236 | −1.045e−6 |
| r ≤ 12, \|j\| ≤ 4 | 910 | 1.053854726398 | −9.597e−8 |
| r ≤ 10, \|j\| ≤ 5 | 1386 | 1.053854765283 | +1.386e−9 |
| r ≤ 12, \|j\| ≤ 5 (paper) | 1638 | 1.053854765588851 | **+2.150864e−9** (paper: 2.15088826449e−9) |

- **Arb certificate** (cert.py). We took our top eigenvector, rounded it to exact dyadic coefficients, and evaluated I and J₄₀ in ball arithmetic, including the L934 remainder (remainder/I = 5.1e−28). Result: rigorous margin **[2.15089187038150e−9 ± 6e−25] > 0**. This matches the paper's coefficient scale: Σ|b|/√I = 7.46e17, i.e. Σ|b| ≈ 7.46 at I = 1e−34.
- **What the paper leaves out.** The coefficients b, the 10⁸⁰-denominator rationals and the Arb script are not in the e-print. None of them are needed: the optimum of the stated finite problem can be recomputed from the stated formulas.
- **The degree-5 space matters.** Prime-sum degree 4 fails. C rests on the five-variable H.

## Cross-checks

- Bondarenko–Heap 2610.26038 (19 Oct 2026, 6 pages, RH): μ ≤ ½ via "resonance to a quotient of shifted zeta functions … Möbius resonator". This is the same object as DPC's Z = e^{πV} = ξ(+)/ξ(−), which is real with sign (−1)^N. Two independent proofs of the ½ bound.
- Inoue 2604.05733: μ < 0.50895 under RH. This is superseded and consistent.
- L83 attributes to [IKT] a large-gap restriction of 3.022 for the Montgomery–Odlyzko criterion with polynomials of length ≤ T. C does not contradict it: DPC's weight is R², with R = Σ_d μ(d)|A_d|² over reduced fractions, which is not |A|² for an integer polynomial A of length ≤ T. Not checked: the bibliography entry for [IKT] is "Explicit extreme values of the argument of the Riemann zeta-function", so the 3.022 attribution may be a mis-citation.

## For the orchestrator

- The TREE.md N8 line and `discoveries/small-gap-floor.md` ("DPC prove D(λ) > 0 for every λ > ½ under RH") are supported. R123's open concern, the unpublished 1638-coefficient computation, is resolved: it was reproduced and certified.
- NOGO candidate: the Möbius-coefficient two-line moment method gives positive limit measures for every f exactly when Im e^{πA_φ} ⪰ 0, i.e. φ ≤ ½. So it cannot go below ½, in line with the Alternative Hypothesis (D7).
- C is a constant record (λ > 3.63846723 under RH), so it belongs in RECORDS.md.
- Resources: one single-threaded job at a time, longest run 225 s. **Memory-cap breach:** the first, unchunked zfinite.py run reached a max RSS of 503 MB (9 s), above the 300 MB cap. The kept script is chunked and aborts itself above 280 MB; its runs peaked at 156–281 MB. All other runs stayed under 300 MB (largegap.py 267 MB, cert.py 174 MB, models.py 213 MB). No worktrees, no commits, no Lean.

RH relevance: N8.
