/-!
# Round 198 — 2610.34300 holds where tested; the diagonal Fasi–Higham case holds for small m; new record leads

* 2610.34300 (Finch): our C cycle finder (Brent, rolling hash; validated on k = 103, 147) gives N_A(199) = 4 and
  N_A(259) = 4; the paper had already checked k = 4m² + 3 for m ≤ 11 and 3^(2j+1) for j ≤ 2, so R197's lead is closed
  (and R196's "m ≤ 30" was a scout error).
* Fasi–Higham (2610.39051): no negative coefficient on the diagonal k = m for odd m ≤ 61, n ≤ 12m (mpmath, dps 2m + 40);
  none for m = 3, k ≤ 20000, n ≤ 3000; none for k ≤ m, m ≤ 151 (float scan, reliable only for m ≲ 50). The dominant
  root is real (mpmath; numpy.roots wrongly reported complex roots from m = 25). Diagonal pair amplitude 2|w₁|/w_r rises
  from 1.25 to 1.90 (m = 5 … 49): an Orca deep dive is open.
* Record leads (November): 2610.36648 (R(3,3,3,3) ≤ 61), 2610.39078 (b₄ ≤ 7), 2610.36675 (average kissing numbers).
* Hygiene: killed an orphaned scout python (catastrophic regex, 17.5 h at 100% CPU).
-/

namespace Fayed.R198
end Fayed.R198
