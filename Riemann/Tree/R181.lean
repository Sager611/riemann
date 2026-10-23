/-!
# Round 181 — Erdős-problem solutions checked; the multiple unicast counterexample is correct; Banks audited

* Correct (others' results, our own checks): Erdős #1147 (answer no for α = √2, built from (3 + 2√2)^(2j+1); checked at
  70 digits); Erdős #193 (an infinite S-walk in Z³ with no three collinear points; 16 displacement vectors; first 300
  points checked); 2610.10108 (Braverman–He: the 182-vertex PG(2, 9) network carries all 157 sessions at coding rate 1
  while fractional flow gives at most 147/157, so the undirected multiple unicast conjecture is false; our own exact F₉
  check, Orca); 2610.31402 (Banks: effective exclusion of Landau–Siegel zeros over a base field; every step checked
  against Masser 1981 and Stark 1974, Orca).
* Numerically supported: Erdős #1129 (Lebesgue-constant minimizers, n ≤ 6). Not checkable in 50 s: #318 (an
  8,801-job Lean build), #190 (statements match; the numbers check).
* {2,4}-distance codes, largest sizes: 4, 9, 16 (q = 3, n = 3, 4, 5); 4, 16 (q = 4, n = 3, 4).
* Wang–Zhi's Connes-embedding certificate is available only from the authors, so it cannot be checked.
* Targets for our own work (open Erdős problems with computational content): #170 (difference bases), #634, #451,
  #663, #854, #962, #1072, #1095.
* Sofa3D (Orca): certified V*(K₄) ≥ 1.8755015 (was 1.8558404); reviewer: significant, but floating-point rounding is not
  yet bounded; fix and novelty check running.
-/

namespace Fayed.R181
end Fayed.R181
