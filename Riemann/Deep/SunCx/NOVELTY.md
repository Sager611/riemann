(1) KNOWN: H. M. Daşdemir, github.com/darkness-38/counter-hunter, paper/counter_hunter_verification.tex (dated and committed 2026-10-17): "the companion linear inequality (1) was disproved for (|A_1|,|A_2|,|A_3|) = (3,2,3)" (no sets, no source given); part (C) also fails (Z/7, n = 4): NEW ; (2) NEW: no earlier report found; for q = 4, 8 it is a two-line corollary of Bors–Wang, arXiv 2207.09642 (2022), Theorem 1.1 and Table 1 (all complete mappings of F_4 and F_8 are even), which does not mention Sun

# Adversarial novelty check and classification (sunnov3, 2026-11-05)

## Verdicts

**(1) Wang–Sun, restricted sumsets.** EJC 33(2) (2026) #P2.64, Conjecture 1 = arXiv 2210.12044 (only v1, 2022-10-21), Conjecture 1.1 = OEIS A357130, Conjecture 2 (2022-09-16). Inequality (L) fails for A = ({0,1,2}, {0,2}, {0,1,2}): L = {2,3,4}, 3 < 4. A public GitHub preprint (above) already says (L) "was disproved" for the size pattern (3,2,3). It names neither the sets nor a source. In Z every (3,2,3) counterexample with elements in [0,10] is an affine image of ours (table below), so this is the same counterexample. Novelty is lost; record it as KNOWN.

New here, not in that preprint (which searched (C) for n = 4, 5, 6 and found nothing):
- **(C) is false too:** in Z/7 with n = 4, take A = ({0,1}, {0,1,2,4,5}, {2,4,5}, {0,1,4,5}). Then S = 14 and the bound is min{7, 14 − 8 + 1} = 7, but C misses 2 (the integer sums miss 2, 9 and 16), so |C| = 6. Among all 4-tuples in Z/7 this is the only example up to x ↦ ux + v, rotation and reversal (336 tuples).
- **(L) fails for every odd n ≥ 3:** put n − 3 copies of {0,1} in front of the example. The prefix must alternate, so L = (n−3)/2 + {2,3,4}, 3 < 4.
- **(L) fails in Z/p for every prime p ≥ 5 in a second way:** take ({0,1}, {0,…,h}, {2} ∪ {4,…,h}) with h = (p+3)/2. Then S = p + 4 and the bound is p, but 4 ∉ L: the integer sums lie in [2, p+4], a sum of 4 would need a_1 = a_2 or a_2 = a_3, and the sum p + 4 = 1 + h + h needs a_2 = a_3. Checked for p ≤ 101; |L| = p − 1.
- **Classification:** see below. No counterexample has equal sizes, even n, or p(G) ≤ 3.

**(2) Sun, OEIS A229827, Conjecture (i), general part (2013).** For every finite field F_q and every P of degree < q − 1 not of the form c − x, some circular order of F_q has the sums P(a_i) + a_{i+1} pairwise distinct. NEW, and wider than the candidate in the brief:
- q = 4: fails for P = ωx + c and ω²x + c (ω² = ω + 1). This is the smallest counterexample.
- q = 8: fails for exactly the even permutation polynomials, 20152 of them besides the excluded x + c. They include every affine bijection ax^4 + bx^2 + cx + d other than x + d, e.g. x², x⁴, λx with λ ∉ F_2, and nonlinear ones such as x⁵ + x³ + x and x⁶ + x³ + x².
- q = 9: fails for exactly 648 odd permutation polynomials: c x³ + d with c ∈ {±1, ±i}, λx + d with λ primitive, and two nonlinear orbits (x⁶ + x⁴ ± x³ + x²).
- No other failure for q ≤ 9, and none found for q = 16 … 81.

Sun's monomial form fails exactly at (q, k) = (8, 2), (8, 4), (9, 3) for q ≤ 81. That form is P = x^k with 0 ≤ k ≤ q − 2, except 2 | q and k = 1 (OEIS revisions #16–#36 of 2013-09-30/10-01, and Number Theory List message 2013-10-01). The already-refuted part (ii) also fails with F = A = F_9 and P = (1+i)x, since deg P = 1 < p − 1.

## Searched (all 2026-11-05)
- **Sources.** arXiv 2210.12044 TeX (v1 only). EJC published PDF (submitted 2025-01-31, published 2026-06-19): Conjecture 1 unchanged, with no remark or erratum; the article page has no errata either. Note the numbering: the EJC Theorems 2, 4, 5 are arXiv Theorems 1.1, 1.2, 1.3.
- **Sun's own view.** His Graz talk of 2026-05-07 (imsc.uni-graz.at/AlgNTh/slides/Slides_070526.pdf, slide 49) still states the L/C conjecture as open.
- **OEIS.** Full revision histories: A357130 (52 revisions, 6 pages) and A229827 (49 revisions, 5 pages). Neither has a correction note; post-2013 edits of A229827 are terms and links (Dobbelaere 2021: P = x² on Z/n, a(n) > 0 iff gcd(n, 6) = 1 for n < 256). Entries citing A357130: none. Entries citing A229827: only A229878 (different problem). Cross-references A227456, A229082, A229141 were read.
- **Number Theory List (NMBRTHRY).** Sun's message "A combinatorial conjecture on finite fields" (2013-10-01) states the monomial version. The list indexes for 2013-09, 10 and 11 show no reply.
- **Citations.** OpenAlex W4307206708 (merged arXiv and EJC records): cited_by 0. Crossref 10.37236/13831: 0. DataCite 10.48550/arXiv.2210.12044: 0. Semantic Scholar (DOI): 0; the arXiv-id query was rate-limited. Google Scholar is not reachable.
- **arXiv.** Author query for Z.-W. Sun (237 papers): 2022–2026 combinatorics papers screened. Abstracts with "restricted sumset(s)" or "Erdős–Heilbronn" since 2022-10: 23 papers. TeX rg'd for 2210.12044, A357130, A229827, "two new kinds" and Wang/Sun:
  - 2610.25073, 2609.23343, 2609.10148, 2605.19542, 2511.23223, 2503.09121, 2411.02365, 2402.01471, 2305.11574, 2606.30506, 2606.18234, 2410.00143, 2602.10402, 2610.06289, 2610.19965, 2309.14979, 2105.00140, 2207.09642.
  - No citation of either source. Of 1309.1679, all versions v1–v9 lack the A229827 conjecture.
- **Sun's homepage** (last modified 2026-08-04) and talks page. Slides read: "Progress on Questions in the book …" (2024), "Problems and Results on Permutations" (2019), "Determinants, Permutations and Additive Combinatorics" (2018), "Introduction to CNT (VI)" (2020). None discusses A229827 (i). The 2021 book is not legally online and was not checked.
- **MathOverflow** (Stack Exchange API). All 252 of Sun's questions: none on these problems after 2022; 363834 is a different problem. Searches for A357130, A229827, "restricted sumsets", the paper title and "circular permutation finite field" found nothing relevant.
- **GitHub code search** (about 25 queries).
  - darkness-38/counter-hunter: the source of (1). Repo created 2026-10-17; the claim and its (C) scan are in commit 25df6d38 (2026-10-17T18:42Z), files src/wang_sun_circular.rs and the paper's §"Wang & Sun (EJC 2026) Circular Restricted Sumset Breaker". The (C) scan covered fixed size patterns (3,2,3,2), (4,2,4,2) and so on, which miss (2,5,3,4).
  - graph-theory-AI/mathpocalypse-project: lists the paper with status "collected" and no verdict.
  - Other hits were unrelated, e.g. OEIS program repositories for A357130.
- **Zenodo API** (phrase queries "restricted sumset(s)", "circular permutation" "finite field", "Zhi-Wei Sun" counterexample, "two new kinds", "Han Wang" sumset, A357130, A229827): nothing on either conjecture.
- **GXL evidence API** (7 claims) and **Exa web search** (4 queries): nothing on either conjecture. They surfaced Bors–Wang 2207.09642: complete mappings of F_4 and F_8 are all even, and F_9 has both parities. It does not mention Sun or circular orders with P(a_i) + a_{i+1}.

## (1) Classification (classify_sumset.c, exhaustive; verify_sumset.py re-checks the starred rows in Python)
Counts are of violating tuples (A_1, …, A_n) with all |A_i| ≥ 2. In the S column, "8: 21" means 21 violations with S = Σ|A_i| = 8.

| G | n | (L) violations (by S) | (C) violations |
|---|---|---|---|
| Z/5 | 3* | 20 (8: 10, 9: 10) | 0 |
| Z/5 | 4*, 6 | 0 | 0 |
| Z/5 | 5* | 470 (12: 270, 13: 200) | 0 |
| Z/7 | 3* | 203 (8: 21, 11: 168, 12: 14) | 0 |
| Z/7 | 4* | 0 | **336 (14)** |
| Z/7 | 5 | 9471 (12: 1113, 14: 546, 15: 7224, 16: 588) | 0 |
| Z/11 | 3 | 495 (8: 55, 15: 440) | 0 |
| Z/11 | 4 (A_1 = {0,1} only) | 0 | 0 |
| Z/13 | 3 | 858 (8: 78, 17: 780) | 0 |
| Z, elements in [0,11) | 3 | 25, all affine images of ({0,1,2},{0,2},{0,1,2}) | 0 |
| Z, [0,7) | 4 | 0 | 0 |
| Z, [0,6) | 5 | 492 (S = 12 … 20) | 0 |
| Z/4, Z/2², Z/6, Z/8, Z/2×Z/4, Z/9, Z/3² (p(G) ≤ 3) | 3, 4 (and 5, 6 for the smallest) | 0 | 0 |

Classes for n = 3, up to x ↦ ux + v and reversal:
- S = 8, for every p ≥ 5 and for Z: ({0,1,2}, {0,2}, {0,1,2}). This is the smallest S, so it is the minimal counterexample.
- Z/5, S = 9: ({0,1}, Z/5, {2,4}).
- Z/7, S = 11: ({0,1}, {0..5}, {2,4,5}) and ({0,1,3}, {0..4}, {2,4,5}); S = 12: ({0,1,3}, {0..5}, {2,4,5}).
- Z/11, S = 15: two classes. Z/13, S = 17: three classes.
- Every class with S > 8 is capped: |L| = p − 1 < p, and the family above covers all p ≥ 5.

The examples sit outside the proved cases. Theorem 1.1 (EJC Theorem 2) needs |A_2| − |A_1| ∈ {0,1} in fields, or sorted sizes in torsion-free groups; our sizes are (3,2,3) and (2,(p+5)/2,(p−1)/2).

**Equal sizes.** Theorems 1.2–1.3 (EJC Theorems 4–5) prove (L) and, for even n, (C) in fields with p(F) large. In the small-p range they leave open, no violation was found:
- Z/5: k = 2, 3, 4 with n = 3, 4, 5, and k = 3 with n = 6, 7, 8.
- Z/7: k = 2…6 with n = 3; k = 3, 4, 5 with n = 4, 5; k = 4, 5 with n = 6; k = 3 with n = 7.
- Z/11: k = 5, 6, 7 with n = 3; k = 7, 8 with n = 4.

Torsion-free groups with equal sizes are covered by those theorems over Q.

**Even n.** No counterexample was found for (L) with n = 4 (Z/5, Z/7, Z/11 with A_1 = {0,1}, and Z on [0,7)) or n = 6 (Z/5).

## (2) Classification

**Symmetries (checked).** The problem depends only on the additive group (F_q, +) ≅ F_p^k.
- Translation and linear maps: if (a_i) works for P, then (g(a_i)) works for gPg⁻¹ for every g(x) = Mx + v with M ∈ GL_k(F_p), because the sums become M·s_i + 2v. This contains Sun's u·P(u⁻¹(x − v)) + v and conjugation by field automorphisms.
- Constants: P + c shifts every sum by c.
- Invariance of the hypotheses: deg P < q − 1 ⟺ Σ_x P(x) = 0, and the family c − x is mapped to itself.
- Inverse: for bijective P also P ↦ P⁻¹, by reversing (P(a_i)).
- So the conjecture concerns all functions with Σ P = 0, not only low-degree polynomials.

**Exhaustive table, q ≤ 9** (orbits under conjugation by AGL(k,p) and adding constants; every orbit decided by complete search):

| q | functions with Σ P = 0 | orbits | failing orbits besides c − x | failing functions | which |
|---|---|---|---|---|---|
| 2 | 2 | 1 | 0 | 0 | |
| 3 | 9 | 3 | 0 | 0 | |
| 4 | 64 | 6 | 1 | 8 | ωx + c, ω²x + c (the even permutations other than x + c) |
| 5 | 625 | 16 | 0 | 0 | |
| 7 | 117 649 | 477 | 0 | 0 | |
| 8 | 2 097 152 | 312 | 10 (sizes 896, 168, 1344, 10752, 336, 5376, 448, 448, 192, 192) | 20 152 | every even permutation other than x + c |
| 9 | 43 046 721 | 11 520 | 5 (sizes 108, 216, 216, 54, 54) | 648 | odd permutations: c x³ + d with c ∈ {±1, ±i} (108); λx + d with λ ∈ {1 ± i} (54) or {−1 ± i} (54); x⁶ + x⁴ + x³ + x² and x⁶ + x⁴ − x³ + x² orbits (216 each) |

Every failing P is a permutation; all non-bijective P work. The best circular order reaches 2 of 4 distinct sums (q = 4), 6 of 8 (q = 8) and 7 of 9 (q = 9).

**Why F_4, F_8, F_9 fail (proof; the finite facts are re-checked by verify_a229827.py).**
1. *Reduction.* Let P be a bijection and (a_i) a circular order, with π: a_i ↦ a_{i+1} the corresponding q-cycle. Put u = π∘P⁻¹. Then the sums P(a_i) + a_{i+1} = (u + id)(P(a_i)) are distinct iff u + id is a bijection, that is, iff u is a complete mapping of (F_q, +). So P works iff u∘P is a q-cycle for some complete mapping u. A q-cycle has sign (−1)^{q−1}.
2. *q = 4, 8.* Every complete mapping is even: for F_4 there are 8, all 3-cycles ωx + c and ω²x + c; for F_8 there are 384, all affine with linear part of order 7, and AGL(3,2) ≤ A_8 because GL(3,2) is simple and translations are even (Bors–Wang, Table 1). A q-cycle with q even is odd. So u∘P is never a q-cycle when P is even, and every even P fails. The computation shows that every odd P works, so the failures are exactly the even permutations; x + c are among them, which is why Sun had to exclude c − x. Examples: x² on F_8 has cycle type (3,3,1,1); λx is a 7-cycle on F_8^×; every element of AGL(3,2) is even. On F_4, x² is a transposition (odd), so it works.
3. *q = 9.* A 9-cycle is even, so u must have the sign of P. Of the 2241 complete mappings of F_9, exactly 108 are odd: x ↦ Mx + c with M ∈ GL(2,3) of order 8. Hence for an odd affine P, every candidate u∘P is affine. But AGL(2,3) has no 9-cycle: a 9-cycle has order 9, so its linear part T has 3-power order; T = I gives a translation (order ≤ 3), and T ≠ I gives f³(x) = x + (T − I)²c = x. So every odd affine P fails, for instance x³ (three transpositions) and λx with λ primitive (an 8-cycle). The 432 odd non-affine failures are finite facts: for each of them none of the 108 maps u∘P is a 9-cycle.
4. *Larger q.* For q = 16 (Bors–Wang give an odd complete mapping) and in every larger case tried, both parities occur and no obstruction appears.

**Larger q (randomized DFS; for P with small image an Euler-circuit construction):**
- With witnesses in a229827_witnesses.txt:
  - Frobenius x^{p^j}, every j, for q = 16, 25, 27, 32, 49, 64, 81: all work.
  - Every monomial x^k (k ≤ q − 2) for q ∈ {16, 25, 27, 32, 49, 64, 81}: all work.
  - λx with λ primitive, one per Galois orbit, for q = 16 … 81: all work.
- Computed, witnesses not stored:
  - all c·x^{p^j} (all c, j) for q = 16, 25, 27, 32, 49 (432 maps);
  - all F_p-linear maps, singular ones included, for q = 16 (65 535), 25 (624), 27 (19 682) and 49 (2 400);
  - 4000 random permutations and functions with Σ P = 0 for q = 16, 2000 for q = 25 and for q = 27, 600 for q = 32, 400 for q = 49.
  - All work.
- Small images need care: x^{(q−1)/d}-type P such as x¹² on F_25, or x^{16} on F_81, defeat plain DFS. The reduction then is to choose δ: F_q → Im P with y ↦ y + δ(y) bijective and with the right counts, then take an Euler circuit of the multigraph δ(y) → P(y). This always succeeded.

## Files
- verify_sumset.py: (1) certificates and counts (25 s).
- verify_a229827.py: (2) lemmas, q ≤ 8 full classification and witnesses (30 s, 85 MB).
- classify_sumset.c, classify_a229827.c: exhaustive runs, with expected outputs in the headers (q = 9: 10.5 s, 50 MB).
- a229827_witnesses.txt: 335 circular orders.
