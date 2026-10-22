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

## Status
To check: Sun's original wording (exceptions?), prime powers q = 8, 9, 16, 25, 27, …, novelty, reviewer.
