/-!
# Round 256 — Sun's OEIS conjecture A229827 (ii) fails for A = {0,1,2}, P(x) = 2x (candidate); other hunts clean

* Counterexample candidate (two scouts independently; orchestrator re-check): Z.-W. Sun, OEIS A229827 Conjecture (ii) (2013):
  for any field F with |F| > 7, finite A in F with |A| = n > 2 and P(x) not of the form c - x (deg P < p - 1 in characteristic
  p), some circular ordering of A has all n sums P(a_i) + a_{i+1} distinct. Over Q (or F_p, p >= 11) take A = {0,1,2},
  P(x) = 2x: the order (0,1,2) gives sums 1, 4, 4 and (0,2,1) gives 2, 5, 2. Also F_11, A = {0,3,6}, P(x) = x^2.
  Reviewer and novelty checks running. Wang-Sun 2210.12044 (R255): reviewer CORRECT and SIGNIFICANT; no prior refutation
  found; journal version EJC 33(2) (2026) #P2.64 keeps the conjecture. Scripts: Riemann/Deep/SunCx/.
* Hold in all tested cases: Sun OEIS A228772 conjecture (Z/m, m <= 15, and Z_3^2, Z_3 x Z_5, Z_5^2); Wang-Sun 0810.0467 Conj. 1.1
  (p <= 7, n <= 3, |A| <= 5); Sun 1901.04837 Conj. 5.1 (p = 3); 2610.36179 Conj. 1.3 (height-1 quivers); 2610.38864 Conj. 2.3.
* Already known or proved in the source: Zenodo 23960915 (in REFUTATIONS), 2610.19193 (lim inf conjecture, not finite),
  2610.08424 Conj. 53-54 (proved there for completely multiplicative functions).
* Cosmetic: 2610.40725 Fig. 1(c): all 45 "unit" edges have length 1.3 (uniform scaling, no mathematical error).
-/

namespace Fayed.R256
end Fayed.R256
