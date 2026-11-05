# Counterexamples to thirteen published conjectures of Z.-W. Sun

[Back to DISCOVERIES](../DISCOVERIES.md)

Elementary and small; each certified by exhaustive search or exact arithmetic. "1309.1679" is Sun's "Some new problems in additive combinatorics" (Nanjing Univ. J. Math. Biquarterly 36 (2019); arXiv 1309.1679), with printed conjecture numbers.

| Conjecture | Counterexample |
|---|---|
| 1309.1679, Conj. 3.2 | Z/6Z, (a₁, …, a₄) = (0, 1, 3, 5): after the prescribed first element, a coset of order 3 forces two equal differences |
| 1309.1679, Conj. 3.4(i) | every 5-subset of Z/7Z: no circular order makes all a_i + 2a_{i+1} distinct |
| 1309.1679, Conj. 3.10(ii) | {2, 3, 4, 6}, since 2·6 = 3·4 |
| arXiv 1811.10503, Conj. 4.9 (also MathOverflow 316142) | n = 3 and 4: F₁₁ with A = {0, 1, 2, 5} (\|S(A)\| = 9 < 11), F₁₃ with A = {0, 1, 3, 4} (10 < 11); the n = 2 case (3-term progressions) was posted first elsewhere (see below) |
| MathOverflow 316143 (the bound) | λ = (1, 1, 2), A = {0, 1, 2, 3}: the sums are {3, 4, 5, 7, 8, 9}, so 6 < 7 over Q |
| OEIS A228766 (general form, in a comment) | Z/6Z, {0, 1, 3, 4} |
| OEIS A228728 (comment) | {0, 1, 2, 3, 4, 5, 8} |
| Wang–Sun, arXiv 2210.12044 = Electron. J. Combin. 33(2) (2026) #P2.64, Conj. 1.1, inequality (C) | Z/7, n = 4, A = ({0, 1}, {0, 1, 2, 4, 5}, {2, 4, 5}, {0, 1, 4, 5}): bound 7, but C misses 2, so \|C\| = 6 |
| OEIS A229827, Conj. (i), general part | Frobenius: F₈ with P(x) = x² or x⁴ (all 5,040 circular orders give ≤ 6 distinct sums of 8), F₉ with P(x) = x³ (≤ 7 of 9); q = 8 also follows from Bors–Wang, arXiv 2207.09642, Thm 1.1 and Table 1, which do not mention Sun |
| arXiv 1208.2683, Conj. 2.1 | n = 9 at the threshold, every 0 < α < 0.0023364: first order in α, θ(23)/9 > log 23 − 1 |
| 1208.2683, Conj. 3.1 (Lucas part) | (A, B) = (3, −26), n = 5: u₆⁷⁰ < u₅⁴² u₇³⁰, so the ratio sequence is not decreasing |
| arXiv 1211.1588, Conj. 4.30(iii) | m = 1, n = 2: 11² · 3³ = 3267 > 5⁵ = 3125 (the only failure for m = 1, n < 34833) |
| arXiv 1402.6641, Conj. 2.25 | n = 5 (also 6, 7, 13, 17, 21): prime-ideal counts 3, 4, 6, 8, 8; Remark 2.25 counts each split prime once |

Likely misprint, listed for completeness: arXiv 1809.07766, Conj. 6.8 fails at n = 7 (det = −2 ≠ 0); with "prime > 3" for "prime > 7" it holds for odd n < 160.

**Excluded as not new.**
- Conjecture 3.8 of 1309.1679 fails exactly at (q, a₀) = (11, 2); first posted by an anonymous Zenodo preprint (10.5281/zenodo.23307915, 2026-10-12); we confirmed it is the only exception for prime powers q ≤ 10⁵.
- 3.3(ii) is corrected in the printed version; A229827(ii) fails on {0, 1, 2}, as Sun's own Remark 1.3 notes.
- Wang–Sun inequality (L) (A₁ = A₃ = {0, 1, 2}, A₂ = {0, 2}): reported for sizes (3, 2, 3) by H. M. Daşdemir (github.com/darkness-38/counter-hunter, 2026-10-17).
- 1811.10503 Conj. 4.9 for n = 2: the-omega-institute issue #16084 (13–14 Oct 2026, Lean refutation).
- arXiv 2607.07638 (Domb numbers) at p = 3: K. Calderon, arXiv 2610.04772v1, Remark rem:p3.

**Review.** Every certificate correct (reviewer agents); all minor except the two Wang–Sun and A229827 items, judged significant; one entry for all.

**Novelty.** Adversarial checks (OEIS histories, MathOverflow, arXiv versions and citing papers, Zenodo, GitHub, Sun's homepage): [SunPrim](../Riemann/Deep/SunPrim/NOVELTY.md), [SunCx](../Riemann/Deep/SunCx/NOVELTY.md), [SunAudit](../Riemann/Deep/SunAudit/NOVELTY.md); novelty is provisional.

**RH relevance.** None (side result).

**Files.** [SunPrim](../Riemann/Deep/SunPrim/), [SunCx](../Riemann/Deep/SunCx/), [SunAudit](../Riemann/Deep/SunAudit/) (REPORT.md, verify.py)

