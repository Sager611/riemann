/-!
# Round 205 — fresh open questions; an RH-related lead (Robin's G on colossally abundant numbers)

* Lead: 2610.39120 (Leshchiner–Musin–Thompson, CA2 numbers and RH) asks (Question 1): if N_i is CA3 and the next
  colossally abundant number is N_{i+1} = 2N_i, must G(N_{i+1}) > G(N_i), G(n) = σ(n)/(n log log n)? (Wu proved it
  for odd primes.) The scout's "No" (2520 → 5040) is invalid: 2520 and 720720 are CA2 (p < log n < q), not CA3.
  Heuristic: with ε = ε(2, a+1), the step raises G iff ε·log 2 > log(log(L + log 2)/log L), L = log N_i; CA3 gives L ≥ q
  with ε > ε(q, 1), which falls short only when q ≤ L < q + ≈0.15: possible for small a. Exact check of all 2-steps next.
* Holds (our check): 2610.33418 (supplementary difference sets (99; 43, 42; 36) and (85; 39, 34; 31): a D-optimal
  design of order 198).
* No counterexample / not answered: 2610.39599 (bracelet counts p ≤ 12), 2610.09474 (5 colours from two parameters),
  2610.37990 (period-24 spirals; L ≤ 3 exhaustive), 2610.41027 (Erdős–Straus grid, n < 1.01·10⁸), 2610.37300
  (qualitative), 2610.39705 (BCC optimal for E_p: Monte Carlo p = 6 … 50).
-/

namespace Fayed.R205
end Fayed.R205
