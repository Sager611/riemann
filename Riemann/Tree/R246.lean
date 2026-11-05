/-!
# Round 246 — Krzyż n = 6, 7 computer-assisted proof audited: no flaw; hunts clean

* Audit (Orca "krzyz", Riemann/Deep/Krzyz/REPORT.md): Wikström's computer-assisted proof of Krzyż's bound for n = 6, 7
  (Zenodo 23964169 v2.0, arXiv 2610.15206): reduction re-derived by hand, all provided checks reproduce, 5,000+
  critical-stratum leaves re-excluded by an independent Arb/LP/Krawczyk prover, full coverage of the n = 6 critical tree
  verified. Not checkable: the Lei–Zhang entropy theorem, the n = 7 tree (54 GB logs withheld), unlogged smaller
  strata. No flaw found. Not our result.
* Minor (orchestrator check): 2611.04165's ancillary torus_ac_threshold.txt prints the "certified" x* interval as
  [0.74612462415580527164, same], a double-precision value 4.5·10⁻¹⁷ below the true root 0.74612462415580531679…;
  the a_c* = √x* bracket it actually uses is correct. Cosmetic, not a flaw in the theorem.
* Verified (others' results, correct): 2611.04523 (M₃(17) = 2¹⁶, even-weight words), 2611.03695 (its (1,3;2;4)₈
  construction has no rainbow triangle). Not testable: 2611.03778, 2611.01980, 2611.03835 (timed out),
  2208.11082 (no Collatz claim in Part II), Zenodo 23728890 (open hypothesis, no proof claimed).
-/

namespace Fayed.R246
end Fayed.R246
