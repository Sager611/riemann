/-!
# Round 197 — corrections to R196 and two open leads

* Correction (R196): the Padé "candidate counterexample" is not ours. 2610.39051 (Scardino) itself refutes Fasi–Higham
  Conjecture 3.2 (SIAM J. Matrix Anal. Appl. 40 (2019)) with c₃₈(66, 5) < 0, c₃₇(44, 7) < 0, c₄₂(45, 9) < 0; we only
  reproduced them (exact rationals and mpmath). Map (n ≤ 300, k ≤ 150): first failure for m = 5, 7, 9, 11 at
  k = 66, 44, 42, 43; none for m = 3; the least-modulus root of q_{k,m} is real throughout. Still open there: the
  diagonal k = m.
* Correction (R196): 2610.34300's conjecture N_A(4m² + 3) = 4 is certified only for m ≤ 6 and N_A(3^(2j+1)) = 4 only for
  j ≤ 1 (state repeats after 4 steps; pre-periods 20, 52, 207, 3217, 6218, 93520, 177479). For m = 7, …, 15 and j = 2, 3
  the sequence does not enter a 4-cycle within 3·10⁶ steps: open lead.
* No counterexample: 2610.31906 (Finsler systole on RP², round metrics only), 2610.38757 (eternally totally positive
  Toeplitz matrices of order 5, sampling), 2610.34912 (multivariate dynamical Chowla, Liouville sums to 1800²).
  The math.AC scout timed out.
-/

namespace Fayed.R197
end Fayed.R197
