/-!
# Round 149 — Sun's OEIS generalizations; data checks; other Sun papers; thorough tests of 3.1, 3.4(i), 3.5(ii)

* Sun's OEIS comments (abelian groups, finite fields): two more false generalizations, verified by us
  (Riemann/Deep/SunPrim/verify3.py; no correction note in OEIS):
  - A228766 (general form: n odd or n ∤ |G|; circular adjacent sums distinct): fails for Z/6Z,
    {0,1,3,4} and Z/10Z, {0,1,5,6} (n = 4), the only failures with |G| ≤ 12.
  - A229827(ii) (any field with |F| > 7, |A| > 2, P not of the form c − x; P(a_i) + a_(i+1) distinct): fails for P(x) = 2x and
    every 3-term progression, e.g. A = {0,1,2} ⊂ Q.
  - A227399's strong version (n > 4, any abelian group) fails in Z/7Z like Conj. 3.5(i); A185645(ii) is Conj. 3.3.
  - Already known: A229910's finite-field form (q = 25, 121; Cohen–Oliveira e Silva–Sutherland–Trudgian; OEIS, Sep 2026).
  - Rejected scout claim: "A229232 has no cycle for q = 11": Sun's own example (1,9,2,4,5,8,10,3,6,7) works for a₀ = −1.
* OEIS data against "a(n) > 0 for n > N" conjectures: 100 of Sun's sequences, no contradiction.
* Sun's other papers: 1807.01181 Conj. 1.1 (F_p, p ≤ 11) and 1811.10503 Conj. 4.1 (Z/m, m ≤ 12): no failure.
* 1309.1679 Conj. 3.4(i) (odd m ≤ 15 and Z/3 × Z/3), 3.5(ii) (m ≤ 9, n = 4, 5), 3.1 (subsets of {1..12}, n ≤ 8): no failure.
* 2610.05776 (least N with a(N) = v is L_(4v−5)): holds for v ≤ 10. Large-q scout: arXiv API rate limit.
-/

namespace Fayed.R149
end Fayed.R149
