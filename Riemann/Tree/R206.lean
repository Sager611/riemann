/-!
# Round 206 — 2610.39120 Question 1 holds for every CA3 doubling step with log N ≤ 6·10⁸

* Question 1 (Robin's G at CA3 doubling steps N → 2N): exact check of every doubling step with exponent a + 1 ≤ 33
  (P ≤ 5.9·10⁸): the CA3 steps are a = 6, 7, 9–32 and all raise G (smallest margin 1.2·10⁻¹⁰ at a = 21; closest
  approach to the failure window L − Q = 0.427 at a = 6). The doubling steps that lower G start from CA1 (a = 1, 2) or
  CA2 (a = 3, 4, 5). A scout's "failure at a = 33" is an artifact: its sieve stopped at 1.05·10⁹ < P₁ = 1.14·10⁹
  (L − Q = −4.6·10⁷ contradicts CA3). Why p = 2 is delicate: the bound ε > ε(Q, 1) falls short at second order by
  (log 2/(2Q² log Q))(1 − log 2 − log 2/log Q), so a failure needs Q ≤ L < Q + ≈0.15 and ε near ε(Q, 1): heuristically
  rare and decreasing in a. The paper had not examined doubling steps numerically.
* No counterexample: 2610.40321 (dice, n = 11 attempt), 2610.38017 (Melas sequences, k = 16), 2610.26940 (touching
  cylinders, restricted search), 2610.39666 (torsional rigidity, ball equality only).
-/

namespace Fayed.R206
end Fayed.R206
