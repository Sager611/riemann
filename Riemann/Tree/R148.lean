/-!
# Round 148 — three more counterexamples in Sun's 2013 paper; novelty of (11, 2); wider search

* arXiv 1309.1679 v9 (Z.-W. Sun, 2013–2020), Section 3: three more conjectures fail on tiny cases (exhaustive checks,
  Riemann/Deep/SunPrim/verify2.py; proofs in its NOTE):
  - Conj. 3.3 (n ∤ |G|, or n even with cyclic Sylow 2-subgroup; b₁ prescribed): G = Z/6Z, (a₁,…,a₄) = (0,1,3,5): after
    b₁ = 0 the coset 1 + {0,2,4} fills the rest, and its two differences d, d′ ∈ {2,4} need d + d′ ≠ 0, so d = d′.
  - Conj. 3.4(ii) (|G| odd, n ∤ |G|, circular differences distinct): fails for every 3-term progression when 3 ∤ |G|
    (differences (−d,−d,2d) or (−2d,d,d)); smallest G = Z/5Z, A = {0,1,2}.
  - Conj. 3.5(i) (3 ∤ |G|, n > 3, circular sums a_i + 2a_(i+1) distinct): every 5-subset of Z/7Z fails; also some
    4-, 5-, 6-subsets of Z/8Z.
  3.1, 3.4(i), 3.5(ii) survive small cases; Section 4 is not finitely testable.
* Conj. 3.10 ((11, 2)): no further failure for primes 67 ≤ q ≤ 113 (1092 instances). Novelty: She, arXiv 2609.19114
  (Hamiltonicity of the additive graphs, q > 18,888,871) does not treat 3.10; Sun's 1405.0290 and 1811.10503 are
  silent; nothing found on (11, 2).
* Recent Sun papers: the conjecture of 2608.18197 (sums of two primes and Fibonacci numbers) holds for n ≤ 2000.
* Finders: 2609.08803 (no C₄-free graph in its case (4) for v ≥ 17; the authors checked v ≤ 16), 2609.16195 (polytope
  slicing questions, not finite). Two scouts timed out.
-/

namespace Fayed.R148
end Fayed.R148
