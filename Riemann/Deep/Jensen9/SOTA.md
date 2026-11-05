# State of the art: Jensen polynomials of Riemann's xi (node N6), checked 2026-11-05

**Short answer.** Our TREE record "d ≤ 8 for all n (GORZ 2019)" is wrong by 24 orders of magnitude. Since Chasse
(thesis 2011, paper 2013) the range of degrees that are hyperbolic for *every* shift n is tied to the height T up
to which RH is verified: d ≤ about T². With Platt–Trudgian's T = 3,000,175,332,800 the proved range before this
round was **d ≤ T² = 9,001,052,027,541,590,755,840,000 ≈ 9.0·10^24** (Kim–Lee 2021; GORTTW Thm 1.2; stated by
Dumas 2026). No explicit N(d) is known for any degree beyond that range; the best large-n wedge is Kim–Lee's
ineffective N(d) = O_ε(d^{1/2+ε}).

Notation: 8ξ(1/2+z) = (−1+4z²)Λ(1/2+z) = Σ γ(n) z^{2n}/n!, J^{d,n}(X) = Σ_j C(d,j) γ(n+j) X^j. These are the
Jensen polynomials of the n-th derivative of Θ(w) = ξ(1/2+√w) (Kim–Lee: Ξ_0(−w)); rescaling by 8 is harmless.

## 1. Degree ranges proved for all n ≥ 0

| year | source | proved range | method | status |
|---|---|---|---|---|
| 1986 | Csordas–Norfolk–Varga, Trans. AMS 296 | d = 2 | Turán inequalities | published |
| 2011 | Dimitrov–Lucas, Proc. AMS 139 | d = 3 | higher Turán inequalities | published |
| 2011/2013 | Chasse, PhD thesis (U. Hawaii) Thm 177 / Complex Var. Elliptic Equ. 58 (2013) Thm 1.8 | d ≤ 2·10^17 | Obreschkoff sector lemma (thesis Thm 171) + RH to height 545,439,823; all n because the sector class is closed under differentiation (stated after Thm 177) | published |
| 2019 | Griffin–Ono–Rolen–Zagier (GORZ), PNAS 116, arXiv 1902.07321, Thm 2 | d ≤ 8 | Hermite model; numerics for n ≤ 10^6 and an effective bound "0 < C(n,j) < 14.25 for n ≥ 7" that is not shown | **proof sketch only**; already weaker than Chasse (its footnote cites Chasse only for n = 0, but Chasse's remark after Thm 177 covers all derivatives) |
| 2021 | Kim–Lee, arXiv 2105.05386 (not in a journal per OpenAlex) | d ≤ T²(1+1/(4T²))², i.e. 9·10^24 for T ≥ 3·10^12 | Obreschkoff + Gauss–Lucas (remark after their Thm 4) | arXiv |
| 2022 | Griffin–Ono–Rolen–Thorner–Tripp–Wagner (GORTTW), Adv. Math. 397, 108186, arXiv 1910.01227, Thm 1.2 / Cor 1.3 | d ≤ ⌊T⌋²; stated 9.36·10^20 (Platt 2017, T = 3.06·10^10) | Chasse's lemma | published |
| 2026 | Dumas, Zenodo 10.5281/zenodo.22812117 (v1.1) | d ≤ 9·10^24 (as d+1 for proper position) | GORTTW Thm 1.2 + Platt–Trudgian | Zenodo |
| exact | GORTTW Thm 1.2 with T = 3,000,175,332,800 | **d ≤ 9,001,052,027,541,590,755,840,000** | | the bar for this round |

Inputs: Platt–Trudgian, Bull. LMS 53 (2021), arXiv 2004.09765, Thm 1: RH is true up to height 3,000,175,332,800,
the lowest 12,363,153,437,138 zeros lie on the line. FLINT 3.6 recomputes N(T) = 12,363,153,437,138 by Turing's
method in 10 s (verify.py part C). No later published extension of the height (R108 agrees).

Second family (same criterion for Ξ itself, J(Ξ^{(n)}; d)): Kim–Lee Thm 4, d ≤ 1 + 4T² ≈ 3.6·10^25.

## 2. Large-n results (n ≥ N(d)); none is explicit beyond the all-n range

| source | statement | constant | note |
|---|---|---|---|
| GORZ 2019 Thm 1 | every d: n ≥ N(d) | ineffective | Hermite limit |
| GORZ 2019 Thm 2 sketch | d ≤ 8: M_ε < 10^6 | explicit, details unpublished | superseded |
| GORTTW 2022 Thm 1.1 | n ≥ c·e^{d/2} | c not computed (O_C constants, M_C) | |
| Kim–Lee 2021 Thm 1 | N(d) = O(d^{c}) for every c > 1/2 | ineffective (Kim 1996 on non-real zeros of derivatives) | best shape |
| O'Sullivan, Res. Math. Sci. 8 (2021), doi 10.1007/s40687-020-00240-5, arXiv 2007.13582 | modified P^{d,n} (Hermite basis): hyperbolic if n/log²n ≥ d^{3/4}/2, d large | not explicit | not J^{d,n}; asymptotics of γ(n) |
| Holland, arXiv 2608.08682 (9 Aug 2026) | n³log²(n+2) ≥ K d^5 ⇒ d distinct zeros; semicircle law | K absolute, not explicit | wedge weaker than Kim–Lee (does not cite Kim–Lee or Chasse) |
| Savva, Zenodo 10.5281/zenodo.22293641 (Sep 2026) | n² log(n+2) ≥ K d³ | not explicit | AI-written, no human review; also weaker than Kim–Lee |

So the best explicit N(d): N(d) = 0 for d ≤ D (all-n range); for larger d nothing explicit is published.

## 3. Other 2023–2026 work on N6

- Liu–Mao, arXiv 2610.04972: the ξ Taylor coefficients are strictly infinitely log-concave (Zhu's conjecture); interval arithmetic + closure argument. Consistent with RH, does not give hyperbolicity.
- Ono, arXiv 2610.38870: unshifted Jensen polynomials of Dedekind zeta functions have ≥ d − C_K√(d log d) real zeros for large d.
- Schieber, Zenodo 10.5281/zenodo.22936341: Marchenko–Pastur/semicircle laws for the zeros, conditional on hyperbolicity.
- Farmer, arXiv 2008.07206: Jensen polynomials are not a plausible route to RH (the asymptotics are universal).
- Zenodo has many 2026 deposits (degree 4, "degree 14", "n ≥ 1 for all d", RH claims); all are below Chasse's 2011 range or unrefereed RH claims.

## 4. Barriers

- The all-n range comes from the angle of the worst possible zero of Θ: d ≤ 1/sin²δ (Obreschkoff). The bound is sharp for the class: one zero pair at angle δ already fails at the next degree (verify.py part D). Beyond it one needs a larger verified height or a narrower zero-free strip; D grows like T²/(4a²) for a strip |Re ρ − 1/2| ≤ a.
- All large-n results have ineffective constants, and n ≫ d^{1/2+ε} is the best shape; uniformity in d with n = 0 is RH.

## 5. Search log

GORZ, GORTTW, Kim–Lee, O'Sullivan, Holland, Liu–Mao, Platt–Trudgian: TeX sources from arXiv e-print. Chasse thesis:
ScholarSpace PDF (DSpace REST endpoint). Dumas: Zenodo TeX. OpenAlex lookups; arXiv listings math.NT, math.CV,
math.CA for 2026-10 and 2026-11 (titles); Zenodo search "Jensen polynomials" (102 records); GXL evidence API (2 claims:
returned GORTTW and Kim–Lee only).
