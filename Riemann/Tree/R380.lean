/-!
# Round 380 — node N4 (Robin / Lagarias), round 8: records stand; 2601.22413's progressive-realization conjecture holds for all 97 superabundant N in (5040, 1e19) (saturated by level 4); a Robin counterexample needs more than 1.3e45 distinct prime factors (routine from our range); Zenodo 23240511 leaves its bridge open; mixing-time cutoff, Zipf laws and juggling give no route

* Node: N4 Robin / Lagarias (eighth round; R372 before).
* [check] state of the art (Zenodo- and HAL-weighted between listings): math.NT listings, Zenodo ('"Robin inequality"',
  '"superabundant"') and HAL: 2610.12816 and 2610.39120 claim no record; nothing new dated 20 Oct - 6 Nov 2026: HOLDS.
* Computation: 2601.22413 conjectures (TeX l. 868) that for every superabundant N > 5040, d^(r)(N) < N implies
  d^(r+1)(N) > d^(r)(N), where d^(r)(N) is the last divisor whose reciprocal prefix sum stays <= S_r log log N (l. 195-210;
  S_1..S_5 = 1, 1.5, 1.6962554944, 1.7590919512, 1.7760747282). Exact divisor lists for all 97 superabundant N in
  (5040, 1e19) (OEIS A004394 b-file): no failure; every N is saturated by level 4, so the test is weak (scout timed out;
  rerun here): HOLDS.
* Wild: Robin for n with omega(n) = 7 (proposed) is trivial: sigma(n)/n < prod_{p|n} p/(p-1), and with Rosser-Schoenfeld's
  prod_{p <= y} p/(p-1) < e^gamma log y (1 + 1/(2 log^2 y)) and our verified range X = 10^(10^46.79) (log log X = 108.572),
  Robin holds for every n > 5040 with omega(n) <= pi(1.41e47) >= 1.30e45 (routine; not recorded as a corollary).
* Adversarial: the newest Zenodo Robin framework, 23240511 (8 Oct), states its CA-to-admissible bridge as an open target
  (TeX sec. 'Coverage of potential counterexamples'); no claim to refute: NOT TESTABLE.
* Harvest: 2610.12816 (seen) asks to determine C_kappa (l. 659) and to control canonical endpoints (l. 1303); neither is a
  finite computation: NONE.
* Exploration: mixing times and cutoff (1904.03533: Ramanujan-graph cutoff sees only the extreme eigenvalues; the graph RH
  is Ihara's analogue), quantitative linguistics (Zipf and Simon laws fix tail exponents of zeta(s) for real s > 1 only), the
  mathematics of juggling (1310.2725: q-rook stationary laws, no Dirichlet series or Mobius link): no route.
-/

namespace Fayed.R380
end Fayed.R380
