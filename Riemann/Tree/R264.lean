/-!
# Round 264 — Sun's 2026 Domb-number conjecture fails at p = 3 (minor); other Sun-related hunts hold or known

* Minor counterexample (scout; orchestrator re-check): Z.-W. Sun, arXiv 2607.07638 (v3, 21 Jul 2026), first conjecture of
  Section 4, part (i): for an odd prime p with (p/3) = 1 or (p/13) = 1, sum_{k<p} D(k) u_k == 0 (mod p^2) (D = Domb numbers,
  u_{n+1} = 11 u_n - u_{n-1}). At p = 3, (3/13) = 1 but 1*0 + 4*1 + 28*11 = 312 == 6 (mod 9). No other failure for p < 400
  (nor of the mod p^3 part or part (ii)), so only the smallest prime is wrong. Script: Riemann/Deep/SunCx/verify_2607.py.
* Already known: 1309.1679 printed Conj. 3.2 also fails in Z_3^2 (n = 4, a point plus a coset of order 3), the general form
  found by Sun38 and on the Sun page.
* Hold in tested cases: 2511.23223 (Sun's conjecture on restricted sums of squares, n <= 5000), 2603.20961 Conj. 1 in Z_23
  with |A| = 21 (one ordering found), Zenodo 23858164 Table 6 (no proof claimed), 2610.41971 (gr_2^3(K_3, K_4, K_4) = 5).
* Not testable: 2405.03626 (parameters a_p, b_p not given explicitly); two scouts hit arXiv rate limits.
-/

namespace Fayed.R264
end Fayed.R264
