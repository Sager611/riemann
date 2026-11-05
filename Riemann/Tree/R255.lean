/-!
# Round 255 — Wang-Sun restricted-sumset conjecture fails for A_1 = A_3 = {0,1,2}, A_2 = {0,2} (candidate); other hunts clean

* Counterexample candidate (scout, orchestrator re-check by exhaustive search): Conjecture 1.1, inequality (L), of H. Wang and
  Z.-W. Sun, arXiv 2210.12044 (= Sun's Conjecture 2(ii) in OEIS A357130, 2022): |L(A_1,...,A_n)| >= min{p(G),
  |A_1|+...+|A_n| - 2n + 1 + (n mod 2)}, where L = sums a_1+...+a_n with a_i in A_i and a_i != a_{i+1}. In Z (or Z/p, p >= 5)
  take n = 3, A_1 = A_3 = {0,1,2}, A_2 = {0,2}: a_2 = 0 forces a_1, a_3 in {1,2}, a_2 = 2 forces a_1, a_3 in {0,1}, so
  L = {2,3,4} has 3 < 4 elements. Exhaustive: 20 violating triples in Z/5, 203 in Z/7; none for n = 4 in Z/5; the cyclic
  version (C) has no violation in these ranges. The paper's Theorem 1.1 needs |A_2| - |A_1| in {0,1} or sorted sizes, so
  only the conjecture fails. Reviewer and novelty checks running.
* Hold in all tested cases: Sun 1304.5988 Conj. 1.2 (n <= 30) and 1.4 (n <= 19); Sun 1807.01181 Conj. 1.1 (p <= 13);
  2508.01163 inertia conjectures (connected graphs n <= 8); OEIS A401001 and A401145 conjectures (extended ranges).
* Verified (others' results, correct): Lean repos of 2611.04172 and 2610.22410 (no sorry/axiom), 2610.33653 (Prop., random
  tests), 2611.01340 (Thm 1 labelling), 2606.17021 (Lemma 1, R(J_3, K_3) = 5).
* Not testable: Zenodo 23728847 (conditional on a hypothesis it states as unproved).
-/

namespace Fayed.R255
end Fayed.R255
