/-!
# Round 367 — node N5 (Nyman-Beurling / Baez-Duarte), round 8: records stand; d_N numerics unchanged (N = 5120); 2611.03924 asks whether its critical-line obstruction is intrinsic to the subspaces; a Zenodo 'proof' is a finite truncation

* Node: N5 Nyman-Beurling / Baez-Duarte (eighth round; R350 before).
* [check] state of the art: the Oct-Nov listings and Zenodo give no new criterion, d_N bound or d_N computation (2610.00217,
  2610.10093, 2610.25002, 2610.41154, 2611.03437 are unrelated dilation or general papers): HOLDS.
* [check] numerics: 2611.03924 is v1 (2 Nov 2026) and 1211.5191 v1; the optimal d_N record stays N = 5120
  (d_N^2 log N = 0.0470 ... 0.0451; the 10^6 range is the Mobius-coefficient error only): HOLDS.
* Harvest: 2611.03924 (TeX l. 552-556) leaves open whether the critical-line obstruction of its Cor. sensor is intrinsic to the
  subspaces span{G_2, ..., G_n} or only to the Mobius coefficients mu(k) (testable by optimal-coefficient Gram computations);
  1211.5191 (l. 61-75) conjectures d_N^2 ~ C/log N under RH (C from the zeros; compare our numerics with the zero-derived C).
* Adversarial: Zenodo 23703001 (Sec. 5, 'FINALER NACHWEIS') computes a finite truncation (N = 4001, finite Ymax, L = 0.033 < 0.045)
  and itself says 'kein formaler Beweis': not a proof.
* Wild (Hardy space of Dirichlet series, dilation completeness): the cyclicity statements found (e.g. 2111.05398) are RH
  equivalences or generalizations; the strongest unconditional one remains our L^p criterion for 1 < p <= 1.14497 (R311): NO IDEA.
* Exploration: positive geometry and the amplituhedron (canonical forms have signed numerators; obstruction already in R347),
  black holes and quasinormal modes (open-system resonances; modular microstate counts do not locate zeros): no new route;
  Donaldson-Thomas wall-crossing: not tested (timeout, arXiv API throttled).
-/

namespace Fayed.R367
end Fayed.R367
