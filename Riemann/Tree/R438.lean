/-!
# Round 438 — node N9 (moments, Lindelöf, subconvexity), round 13, falsification of the inputs of our sub-Weyl record 54.6 t^{27/164}: Patel–Yang 2302.13444's explicit k-th derivative test (k = 3, 4, 5, both readings of its recursion), improved Poisson formula and stationary-phase lemma hold at random, corner and adversarial parameters with wide margins; journal version J. Number Theory 262 (2024), no error reported by 12 citing works; Yang's source 2301.03165 still v1; Backus–Gilbert inversion, CLEAN deconvolution and Schur–Cohn tests transfer no tool

* Node: N9 moments, Lindelöf, subconvexity (13th round; R433 before). Round type: falsify the inputs of our record (the sub-Weyl one; Orca
  agents revaudit and pyaudit referee 2602.05614 and 2302.13444 line by line in parallel). Exact sums in numpy, ratios |LHS|/RHS, η over
  {0.01, …, 10} (the lemma holds for every η, so the largest ratio over η is tested).
* [falsify] Lemma kth_deriv_test (2302.13444v1 TeX l. 129–150, from Yang 2301.03165 Lem. 2.4–2.5), k = 3: 300 log phases (t/2π)log x and
  300 adversarial cubics c x³ + (p/q)x² + βx (c ≥ 4·10⁻¹¹, q ≤ 12): max ratio 0.2438 (log, N = 2) and 0.4109 (cubic): HOLDS (scout).
  k = 4, 5 (A_k, B_k by the recursion, both readings of its K): 120 log phases each, max 0.27; exact dyadic small-c polynomials c x^k + (p/q)x²
  (c from 10⁻¹² to 10⁻⁸): max 0.019; the scout's flagged case (c = 8.16·10⁻¹⁰, p/q = 7/2, N = 90,254) has |S| = 279.6 at 30 digits, ratio
  0.003: HOLDS (orchestrator; the kth45 scout timed out).
* [falsify] Lemma imp_poisson_summation_formula (l. 467–475; the B-process of our certificate): 300 random log-phase blocks on both sides of
  the stationary range, integer shifts, corners b − a = 1, 2: max ratio 0.0887 (mpmath confirms), smallest margin RHS − LHS = 9.08: HOLDS.
* [falsify] Lemma titchmarsh_lem_46 (l. 340–347): 200 stationary-phase integrals, corners with |f′(a)| or |f′(b)| small: max ratio 0.154
  (the scout's first maximum 42.2 was quadrature aliasing, removed by a 30-digit recheck): HOLDS.
* [check] 2302.13444 is still v1 on arXiv, published in J. Number Theory 262 (2024) (doi 10.1016/j.jnt.2024.03.011); OpenAlex lists 12
  citing works, none reports an error or a smaller explicit 27/164 constant for all t ≥ 3; Yang's derivative-test source is arXiv
  2301.03165, v1, no erratum found: HOLDS. Our Lemma 3.5 slip (R396) stands as the only error found.
* Exploration (tool transfers, tested): Backus–Gilbert minimum-variance kernels compute d_N for N5 (d_2² = 0.0906278455 on span{ρ_{1/2},
  ρ_{1/3}}) but d_N → 0 stays RH; CLEAN coordinate descent decreases the Weil energy only once Q ⪰ 0 is assumed; the Schur–Cohn test
  certifies finite polynomials while the mapped zeros lie on the boundary: no transfer.
-/

namespace Fayed.R438
end Fayed.R438
