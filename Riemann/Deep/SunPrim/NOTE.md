# A counterexample to Sun's circular-permutation conjecture for primitive elements (q = 11, a₀ = 2)

Date 2026-10-22. Found by a counterexample-hunting scout (Riemann/Tree/R145.lean); proof and script below.

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
- Original: Z.-W. Sun, arXiv 1309.1679, Conjecture 3.10 (TeX l. 629–633): every q > 7, every a₀, no exceptions. arXiv 2610.06289 cites it as Conjecture 3.8 (a numbering slip). Neither source reports small-q checks.
- Reviewer: correct; significant if new.
- Wider search: no other failing (q, a₀) for primes q ≤ 61 or q = 8, 9, 16, 25, 27 (all a₀).
- Novelty so far: arXiv 1309.1679 v9 (2020) unchanged, with a q = 11 example for a₀ = −1 only; OEIS A229232 states the cases a₀ = ±1; nothing found on (11, 2). Pending: the citing papers ("Hamiltonicity in graphs defined by primes and primitive elements", Sun's later notes), adversarial novelty check.

## More counterexamples in the same paper (arXiv 1309.1679 v9, Section 3)
Found by a scout (Riemann/Tree/R148.lean), checked exhaustively by `python3 verify2.py` (< 1 s). Statements quoted from the TeX (l. 512–516, 525–534, 541–555).

- **Conjecture 3.3** (2013-09-02): *a₁, …, a_n distinct in a finite abelian group G; n ∤ |G|, or n even with cyclic Sylow 2-subgroup. Then some permutation (b₁, …, b_n) with b₁ = a₁ has the n − 1 differences b_i − b_{i+1} pairwise distinct.* Fails for G = Z/6Z, (a₁, …, a₄) = (0, 1, 3, 5) (4 ∤ 6, and the Sylow 2-subgroup is cyclic): after b₁ = 0 the last three entries are the coset 1 + {0, 2, 4}, whose two consecutive differences d, d′ lie in {2, 4} with d + d′ = b₂ − b₄ ≠ 0, forcing d = d′. Same in Z/9Z with (0, 1, 4, 7).
- **Conjecture 3.4(ii)** (2013-09-03): *n > 2, A an n-subset of G with |G| odd and n ∤ |G|. Then some numbering has the n circular differences a₁ − a₂, …, a_n − a₁ pairwise distinct.* Fails for every 3-term progression {x, x + d, x + 2d} when 3 ∤ |G|: up to rotation the differences are (−d, −d, 2d) or (−2d, d, d). Smallest: G = Z/5Z, A = {0, 1, 2}.
- **Conjecture 3.5(i)** (2013-09-20): *|A| = n > 3, G finite with 3 ∤ |G|. Then some numbering has the n circular sums a_i + 2a_{i+1} pairwise distinct.* Fails for G = Z/7Z and every 5-subset (21 of 21; e.g. A = {0, 1, 2, 3, 4}, all 120 numberings fail), and in Z/8Z for 8, 8, 4 of the subsets of sizes 4, 5, 6.
- Survive the small cases: 3.1, 3.4(i), 3.5(ii); Section 4 is not finitely testable.
- Status: pending adversarial novelty check (these are elementary; someone may have noticed).
