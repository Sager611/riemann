# A counterexample to Sun's circular-permutation conjecture for primitive elements (q = 11, a₀ = 2)

Date 2026-10-22. Found by a counterexample-hunting scout (Riemann/Tree/R145.lean); proof and script below.

Numbering: conjectures of arXiv 1309.1679 are cited by their printed numbers (arXiv v9 = Nanjing Univ. J. Math. Biquarterly 36 (2019) 134–155). The TeX labels differ: printed 3.k is the k-th conjecture of Section 3, so label Conj3.3 → 3.2, Conj3.4 → 3.3, Conj3.5 → 3.4, Conj3.7 → 3.5 (Sun–Hou), Conj3.8 → 3.6, Conj3.9 → 3.7, Conj3.10 → 3.8, Conj3.11 → 3.9, Conj3.12 → 3.10, Conj3.13–Conj3.18 → 3.11–3.16. Our logs R146–R150 used the labels.

## Statement refuted
arXiv 2610.06289 (Conjecture 1.1, TeX l. 74–80; attributed to Z.-W. Sun, 2013, Conjecture 3.8): *Let F_q be a finite field with q > 7. For any given a₀ ∈ F_q there exists a circular permutation (a₁, …, a_{q−1}) of F_q^* such that a₀ + a₁a₂, a₀ + a₂a₃, …, a₀ + a_{q−1}a₁ are primitive elements of F_q.* Equivalently (the paper's own reformulation, l. 85–88): the graph G_{q,a₀} on F_q^* with xy an edge iff x ≠ y and a₀ + xy is primitive is Hamiltonian. The paper proves this for q > 18 888 871 (and for a₀ = 0, q > 4).

## Counterexample: q = 11, a₀ = 2
The primitive roots mod 11 are 2, 6, 7, 8. In G_{11,2}:
- N(2): 2 + 2y ∈ {2, 6, 7, 8} ⇔ 2y ∈ {0, 4, 5, 6} ⇔ y ∈ {0, 2, 8, 3}; excluding y = 0 (not in F_11^*) and y = 2 (loop): N(2) = {3, 8}.
- N(9): 9y ∈ {0, 4, 5, 6} ⇔ y ∈ {0, 9, 3, 8} (9⁻¹ = 5); excluding 0 and 9: N(9) = {3, 8}.

A Hamilton cycle must use both edges at every vertex of degree 2, so it contains 2–3, 2–8, 9–3, 9–8, i.e. the closed 4-cycle 2–3–9–8–2 on only 4 of the 10 vertices. Contradiction: G_{11,2} is not Hamiltonian, and no circular permutation exists for (q, a₀) = (11, 2).

Full neighbourhoods: 1: {4,5,6}, 2: {3,8}, 3: {2,5,9}, 4: {1,7}, 5: {1,3,10}, 6: {1,8,10}, 7: {4,10}, 8: {2,6,9}, 9: {3,8}, 10: {5,6,7}.

## Verification
`python3 verify.py` (< 1 s): recomputes the primitive roots and the neighbourhoods and confirms by exhaustive search that no Hamilton cycle exists. A scout also searched every prime q ≤ 61 and every a₀: (11, 2) is the only failure with q > 7 among primes up to 61.

## Status (2026-10-22)
- Original: Z.-W. Sun, arXiv 1309.1679, Conjecture 3.8 (TeX label Conj3.10, l. 629–633): every q > 7, every a₀, no exceptions; arXiv 2610.06289 cites it correctly as Conjecture 3.8. Neither source reports small-q checks.
- Reviewer: correct; significant if new.
- Wider search: no other failing (q, a₀) for primes q ≤ 61 or q = 8, 9, 16, 25, 27 (all a₀).
- Novelty so far: arXiv 1309.1679 v9 (2020) unchanged, with a q = 11 example for a₀ = −1 only; OEIS A229232 states the cases a₀ = ±1; nothing found on (11, 2). Pending: the citing papers ("Hamiltonicity in graphs defined by primes and primitive elements", Sun's later notes), adversarial novelty check.
- Novelty check (2026-10-22, see NOVELTY.md): KNOWN. Anonymous, "Sun's circular primitive-permutation conjecture: one exception and a complete classification", Zenodo 10.5281/zenodo.23307915 (2026-10-12), gives the same proof that G_{11,2} is not Hamiltonian and claims (11, 2) is the only exception. `census.py` confirms that (11, 2) is the only failure among all 2388 pairs (q, a₀) with 7 < q ≤ 128.

## More counterexamples in the same paper (arXiv 1309.1679 v9, Section 3)
Found by a scout (Riemann/Tree/R148.lean), checked exhaustively by `python3 verify2.py` (< 1 s). Statements quoted from the TeX (l. 512–516, 525–534, 541–555).

- **Conjecture 3.2** (label Conj3.3; 2013-09-02): *a₁, …, a_n distinct in a finite abelian group G; n ∤ |G|, or n even with cyclic Sylow 2-subgroup. Then some permutation (b₁, …, b_n) with b₁ = a₁ has the n − 1 differences b_i − b_{i+1} pairwise distinct.* Fails for G = Z/6Z, (a₁, …, a₄) = (0, 1, 3, 5) (4 ∤ 6, and the Sylow 2-subgroup is cyclic): after b₁ = 0 the last three entries are the coset 1 + {0, 2, 4}, whose two consecutive differences d, d′ lie in {2, 4} with d + d′ = b₂ − b₄ ≠ 0, forcing d = d′. Same in Z/9Z with (0, 1, 4, 7).
- **Conjecture 3.3(ii)** (label Conj3.4; 2013-09-03): *n > 2, A an n-subset of G with |G| odd and n ∤ |G|. Then some numbering has the n circular differences a₁ − a₂, …, a_n − a₁ pairwise distinct.* Fails for every 3-term progression {x, x + d, x + 2d} when 3 ∤ |G|: up to rotation the differences are (−d, −d, 2d) or (−2d, d, d). Smallest: G = Z/5Z, A = {0, 1, 2}.
- **Conjecture 3.4(i)** (label Conj3.5; 2013-09-20): *|A| = n > 3, G finite with 3 ∤ |G|. Then some numbering has the n circular sums a_i + 2a_{i+1} pairwise distinct.* Fails for G = Z/7Z and every 5-subset (21 of 21; e.g. A = {0, 1, 2, 3, 4}, all 120 numberings fail), and in Z/8Z for 8, 8, 4 of the subsets of sizes 4, 5, 6.
- Survive the small cases: 3.1, 3.3(i), 3.4(ii); Section 4 is not finitely testable.
- Status: pending adversarial novelty check (these are elementary; someone may have noticed).

## Two more false generalizations in Sun's OEIS comments (checked 2026-10-22; no correction note in OEIS)
Checked by `python3 verify3.py` (< 2 s).
- **A228766** (comment, general form): *if a₁, …, a_n are n > 2 distinct elements of a finite additive abelian group G with n odd or |G| not divisible by n, then some circular permutation has the n adjacent sums pairwise distinct.* Fails for G = Z/6Z, {0, 1, 3, 4} and G = Z/10Z, {0, 1, 5, 6} (n = 4); these are the only failing (m, n) with m = |G| ≤ 12 in Z/m (3 and 10 failing 4-subsets; 2 and 4 of them contain 0).
- **A229827 (ii)**: *F any field with |F| > 7, A ⊂ F finite with |A| = n > 2, P a polynomial (degree < p − 1 in characteristic p) not of the form c − x: some circular permutation has the n sums P(a₁) + a₂, …, P(a_n) + a₁ pairwise distinct.* Fails for P(x) = 2x and every 3-term progression {x, x + d, x + 2d}: the order (x, x + d, x + 2d) gives 2(x + d) + (x + 2d) = 2(x + 2d) + x, and the other order gives 2x + (x + 2d) = 2(x + d) + x. E.g. A = {0, 1, 2} in Q or in F₁₁.
- Related: A227399's "strong version" (any abelian group, n > 4) fails in Z/7Z like Conjecture 3.4(i); A185645(ii) is Conjecture 3.2. Already known: the finite-field form in A229910 (fails at q = 25, 121; Cohen–Oliveira e Silva–Sutherland–Trudgian, OEIS note by E. Rodríguez, Sep 2026).

## Two more (checked 2026-10-22; `python3 verify4.py`, < 1 s)
- **Conjecture 3.10(ii)** of arXiv 1309.1679 (label Conj3.12, TeX l. 653–673): *for any set A of n > 3 distinct nonzero reals there is a circular permutation whose n adjacent differences are pairwise distinct and whose n adjacent products are pairwise distinct, unless |A| = 4 and A = {±s, ±t}.* Fails for A = {2, 3, 4, 6} (2·6 = 3·4; every circular order repeats a difference or a product), also {−8, −2, 1, 4}, {−4, −1, 2, 8}, {−6, −4, −3, −2} (the only failures with n = 4 in [−8, 8]), and {±1/3, ±1, ±3} (n = 6). No failure with n = 5 in [−6, 6]; 3.10(i) survives all 25,776 integer sets tested.
- **A228728** (Sun, comment: "somewhat curious but we are unable to find a counterexample"): *if distinct reals a₁ < … < a_n admit a cyclic arrangement with pairwise distinct cyclic distances, then some arrangement with c₁ = a₁ and c_n = a_n has pairwise distinct cyclic distances |c₁ − c₂|, …, |c_n − c₁|.* Fails for {0, 1, 2, 3, 4, 5, 8}: the arrangement (0, 3, 4, 2, 8, 1, 5) has distances 3, 1, 2, 6, 7, 4, 5, but none of the 120 arrangements with ends 0 and 8 works. Among subsets of {0, …, 12} containing 0 with n ≤ 7 there are 4 such sets, all with n = 7.

## Sun, arXiv 1811.10503, Conjecture 4.9 (checked 2026-10-22; `python3 verify5.py`, < 1 s)
- *n ≥ 1, F a field with p(F) > n + 1 (the characteristic, or +∞), A ⊂ F finite with |A| ≥ n + δ_{n,3}: S(A) = {Σ_{k≤n} k a_k : a₁, …, a_n distinct in A} has |S(A)| ≥ min{p(F), (|A| − n)n(n+1)/2 + n(n² − 1)/6 + 1}.* Fails for n = 2 and every 3-term progression, e.g. A = {0, 1, 2}: S(A) = {a + 2b : a ≠ b} = {1, 2, 4, 5} has 4 elements, but the bound is min{p(F), 5} = 5 for every field with p(F) ≥ 5, including Q. (Sun added δ_{n,3} for n = 3 but not this n = 2 case.) The same statement for n = 2 (λ = (1, 2)) is MathOverflow question 316142 (Sun, 2018, no answers); it fails the same way.

## Sun, MathOverflow question 316143 (2018-11-25; no answers) (checked 2026-10-22; `python3 verify6.py`, < 1 s)
- *λ₁, …, λ_n (n ≥ 3) positive integers with gcd 1 and λ_{k+1} − λ_k ∈ {0, 1}; F a field with p(F) > n + 1; A ⊂ F finite with |A| ≥ n + δ_{n,3}: |{Σ λ_k a_k : a_k distinct in A}| ≥ min{p(F), (λ₁ + … + λ_n)(|A| − n) + Σ_{k ≤ n/2} (n + 1 − 2k)(λ_{n+1−k} − λ_k) + 1}.* The bound is the length of the integer interval [min, max] for A = {1, …, m}, and Sun's motivation assumes the sums fill that interval. They need not: for λ = (1, 1, 2) and A = {0, 1, 2, 3} the sums are {3, 4, 5, 7, 8, 9} (6 is missing), so 6 < 7 over Q, and the bound also fails in F_p for p = 5, 7, 11, 13. Other failures with A = {0, …, m − 1}: λ = (2, 3, 3) (m = 4, …, 7: two short), (3, 3, 4), (3, 4, 4).

## Novelty check and census (2026-10-22; details and sources in NOVELTY.md)
- Verdicts:
  - Known: 3.8 (Zenodo 10.5281/zenodo.23307915).
  - Known in substance: 3.3(ii). Sun excluded n = 3 in v4–v8, and his homepage copy of the journal version (196a.pdf, 2020-05-07) restricts (ii) to n = |G| − 1. That version is the known R-sequenceability of odd-order abelian groups.
  - Known in substance: A229827(ii), refuted by Sun's own Remark 1.3 of v9.
  - New: 3.2, 3.4(i), A228766, 3.10(ii), A228728 (a conjecture silently dropped from v9), 1811.10503 Conj. 4.9 (= MO 316142, unanswered), MO 316143.
- Version caveat: the homepage journal copy differs from arXiv v9 in exactly one place, 3.3(ii) ("If n ∤ |G|" in v9 becomes "If n = |G| − 1").
- `python3 census.py` (45 s, 21 MB; `--ext` adds Z/13–Z/16; `--primes 17` checks Z/17 for 3.3(ii) and 3.4(i)) checks all abelian groups of order ≤ 12 exhaustively, up to translations and automorphisms.
  - 3.2: a permutation with b₁ = a₁ exists iff neither A nor A∖{a₁} is a coset of a subgroup H with |H| ≥ 3 and Σ_{h∈H} h = 0. This is asserted for every subset and every a₁ (also Z/13–Z/16).
  - 3.3(ii): failures of the v9 version; the n = |G| − 1 version never fails.
  - 3.4(i) and the A228766 general form: failures tabulated.
  - 3.8: the Hamilton search over q ≤ 128 re-checks every cycle with separate polynomial arithmetic.
