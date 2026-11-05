/-!
# Round 258 — finite-field Frobenius hunt: only known exceptions; a false step in a Collatz-permutation preprint

* False step (scout; orchestrator re-read the PDF): Zenodo 23909592 ("The original Collatz permutation of 1932", 3 Nov
  2026), Theorem 8.17: "this holds in particular if L <= m" rests on "L <= m implies 1/(2L^2) >= 1/(2m)", which needs
  L^2 <= m (L = m = 2: 1/8 < 1/4). Logged in REFUTATIONS.md (minor; the rest of the theorem stands).
* Known: Sun OEIS A229910 finite-field clause fails at q = 25 (and 121): Cohen, Oliveira e Silva, Sutherland, Trudgian;
  1711.06392 Cor. 2(i); our R149. A229827(i) at F_8, F_9 is our R257 candidate (scouts re-derived it).
* Rejected scout claim: 2602.17148's period conjecture "omits bijectivity": the paper restricts to bijective f (its search
  and observation are for bijective maps), so f = 0 is outside it.
* Verified (others' results, correct): 2611.05962 (values of f_0), 2611.05961 (22 classes, orders 5-7), 2611.04688
  (digraph example). Not testable: 2611.05433 Thm 6.
-/

namespace Fayed.R258
end Fayed.R258
