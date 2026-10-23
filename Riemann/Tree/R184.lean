/-!
# Round 184 — the prefix-reversal conjecture of 2610.30892 is false (two independent checks); other hunts clean

* Counterexample, confirmed twice (own closure and sympy): arXiv 2610.30892, Conjecture 4 (TeX l. 268): "if m + k < n then
  ⟨r_n, r_m, r_k⟩ < S_n", with r_i reversing the first i positions (l. 87–99), 2 ≤ k < m < n. But ⟨r_7, r_4, r_2⟩ = S_7
  (order 5040) and ⟨r_11, r_6, r_4⟩ = S_11; both lie in the case n odd, m + k = n − 1, which the paper's theorem leaves
  open (l. 221–223). These are the only such triples with n ≤ 11. An Orca agent is characterizing the open case.
* Misread, no counterexample: OEIS A401236 (its terms lie in A373846 by definition; the reported "5" was an
  intermediate value).
* Already known: 2307.15996 (locked polyomino tilings: the connectivity conjecture fails at n = 6, 50 locked tilings
  among 80,092; the paper's own v5 states this).
* No counterexample: 2610.31392 (fpd versus FPD formulas, d, n ≤ 10⁴), 2610.30845 (bi-connected Yang–Baxter solutions,
  |X| ≤ 3: 1, 3, 12), 2601.21958 (linear systems on blow-ups of F_0, 174 cases), 2610.10665 (finite proxies only).
-/

namespace Fayed.R184
end Fayed.R184
