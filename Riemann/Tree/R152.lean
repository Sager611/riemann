/-!
# Round 152 — Sun–Zhao, MathOverflow, 1811.10503 (rest), OEIS sweeps

* Sun's MathOverflow question 316143 (2018, unanswered): the lower bound for |{Σ λ_k a_k : a_k distinct in A}|
  (λ positive, gcd 1, steps in {0, 1}, n ≥ 3) fails already for A = {0,1,2,3}, λ = (1,1,2): the sums are
  {3,4,5,7,8,9}, 6 < 7, over Q and F_5, F_7, F_11, F_13; also λ = (2,3,3), (3,3,4), (3,4,4) (Riemann/Deep/SunPrim/verify6.py).
* MathOverflow 316699 (sums Σ a_k a_π(k) over the nonzero squares cover F_q): fails for q = 5, 7, but the answers treat
  large q; not pursued. MO 316921, 301832, 303632: no failure in the tested ranges.
* Sun–Zhao linear extension of the Erdős–Heilbronn conjecture (arXiv 0810.0467, Conj. 1.1): 9,100,021 cases (p ≤ 13,
  n = 2, 3, |A| ≤ 6, all nonzero coefficients): no failure.
* 1811.10503: Sections 1–3 have no conjectures; Conj. 4.2–4.8 hold for n ≤ 9 or 10.
* OEIS: A229827(i) (q ≤ 7), A227399(i) (Z/m, m ≤ 12, n = 5, 6; n = 7 for m ≤ 9; integer sets): no failure. Two sweeps timed out.
-/

namespace Fayed.R152
end Fayed.R152
