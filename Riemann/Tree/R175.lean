/-!
# Round 175 — October's claimed counterexamples triaged; an OEIS conjecture false as written; RamRepair under review

* Look correct (our partial checks): 2610.13848 (Banks–Ghomi–Qiu–Raffaelli, counterexample to Nirenberg's conjecture on
  closed asymptotic curves; closure integrals ≈ 6·10⁻¹⁵, sampled certificate min E = 0.1266 > 0); 2610.00153 (Cheng, Yau's
  Problem 52; |Rm(p_j)| = 12a(j+3) exactly); 2610.13979 (Guo, three conjectures of Luo; the two-tetrahedron L(3,1)
  certificate checks exactly).
* Doubtful, not checkable from the source: 2610.01536 (Wang–Zhi, an explicit degree-72 polynomial counterexample to
  Connes' embedding conjecture; its 18,326-term certificate is not included); 2610.31589 (Simon, thrackles have at most
  ⌊4(n−1)/3⌋ edges; rests on an unpublished computation for n ≤ 12); 2610.02809 (Dong–Cao, 0.6421 ≤ β_TSP ≤ 0.8810;
  the enumerations are only asserted; its final fraction is 0.88097907 < 0.8810); 2610.08581 (Gelander, continuous
  Aldous–Lyons; plausible, a compressed step at TeX l. 823–826).
* OEIS A401145: the comment "for n ≥ 7, a(n) > max(a(n−1), a(n+1))" is false as written at n = 8 (a(7), a(8), a(9) =
  3, 1, 3) and cannot hold for two consecutive n, so it is probably a typo; "a(n) > 0 for n ≠ 2" holds for n ≤ 112.
* Orca result under review: RamRepair claims Ramanujan's inequality π(x)² < (ex/log x) π(x/e) for every x ≥
  38,358,837,683, unconditionally. Four reviewers: the logic and coverage are correct; the certification has gaps
  (outward rounding of float conversions, the E2 bound for small arguments, smoothing constants for c = 75, 95, 105,
  upstream scripts not included). A fix pass and a clean rerun are running, with a novelty check.
-/

namespace Fayed.R175
end Fayed.R175
