/-!
# Round 257 — Sun's OEIS conjectures: A229827(i) fails for Frobenius in F_8 and F_9 (candidate); A229827(ii) withdrawn (known)

* Counterexample candidate (scout found F_9; orchestrator exhaustive check found F_8 too): Z.-W. Sun, OEIS A229827
  Conjecture (i), general part (2013): for every finite field F_q and every P(x) of degree < q - 1 not of the form c - x,
  some circular permutation a_1..a_q of F_q has {P(a_1)+a_2, ..., P(a_q)+a_1} = F_q. It fails for the Frobenius map:
  q = 8 with P(x) = x^2 or x^4 (all 5,040 circular orders: at most 6 distinct sums), q = 9 with P(x) = x^3 (all 40,320:
  at most 7). Our R152 had tested only q <= 7. The Frobenius cases q = 4, 16, 25, 27, 49 do have such permutations (DFS).
  Script: Riemann/Deep/SunCx/verify_a229827.py.
* Withdrawn (R256): the failure of A229827(ii) at A = {0,1,2}, P(x) = 2x is already noted by Sun (arXiv 1309.1679,
  Remark after Theorem 1.3) and was logged in R149/R155 (discoveries/sun-counterexamples.md, "excluded").
* Rediscovered, already on that page: A228766 (Z/6, {0,1,3,4}), A227399 strong form (5-subsets of Z/7 = printed Conj. 3.4(i)).
* Hold in tested cases: A356714 Conj. 1-2 (n = 4, 5 in F_11, F_13, Z_3^2), A229878, A229232 (q = 8, 9, 17; n <= 8),
  A273287, A273488 (n <= 45), A276270 (n <= 100), A232891, A232894, A232898 (stated ranges).
* Scout error: "q = 8 passes for degree <= 3" (a search bug; x^2 fails).
-/

namespace Fayed.R257
end Fayed.R257
