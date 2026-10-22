/-!
# Round 128 — verification of our own candidates: primary sources and record novelty (11 scouts, 2 reviewers)

* Weil positivity, half-width 1.55 (weilrec4): the reviewer says RIGOROUS (cert4.py, run4/ logs: margins
  1.215e−111 even, 4.004e−107 odd; Gram pieces each covered exactly once). Novelty, 8 sources checked
  (GitHub: Mori log 15 = half 1.354, Elston 1.495, jose-paul-b 0.969, sectio-aurea-q 0.9, davidangularme log 2;
  Zenodo: Watts, Piza, Ziad log 7, Mehennaoui): nothing larger. The scout's one "YES" (darreal44, µ = 11) is a
  47×47 truncated form, not a window certificate. → DISCOVERIES entry raised from 3/2 to 1.55.
  The log 5 run was stopped: its certificate step would need about 10 GB.
* Sato–Tate in every ideal class (chebst): reviewer UNSURE, no internal gap. Primary sources:
  - Pascadi–Thorner Cor 1.2: matches our use;
  - Humphries–Thorner Thm 4.1: matches our use;
  - "Thorner's zero-free region (J5)" is NOT in arXiv 2002.10450. That paper excludes all but one real simple
    zero for self-dual π × π′.
  → Orca chebstfix: repair, and check the earlier Frobenius-sign entry, which uses the same input.
* Suzuki screw-function certificate (Mittermeier, Zenodo 23377174): Ψ(t) > 0 for 0 < t ≤ log 10^10, with
  directed rounding. This is a pointwise RH-equivalent check, not a Weil-window (quadratic-form) certificate,
  so it is not comparable with our window record.
-/

namespace Fayed.R128
end Fayed.R128
