/-!
# Round 234 — hunts clean; a scout's flaw claim was a pdftotext artifact

* Rejected (orchestrator check of the TeX source): 2611.04599 (Kaplansky's conjectures fail for torsion-free hyperbolic
  groups in characteristic p): Lemma 7.4 defines z([i:j:1]) = ⟨i − h[j=0]⟩, the centred representative (main.tex l. 781);
  pdftotext dropped the angle brackets, which produced the scout's "z([2:0:1]) = −3 at p = 3". The lemma and the
  paper's own p = 3 remark are consistent.
* Verified (others' results, correct): 2611.03527 (failure of the critical orthonormal Strichartz estimate, N = 100
  witness, two quadratures), 2610.28820 (all GF(2) generators n ≤ 4).
* No counterexample: 2611.01703 (von Neumann inequality sample), 2610.34053 (n = 2 case). Not testable: 2611.03607,
  2611.04699 (unspecified constants). Unverified scout report: Zenodo 23680304 (Hilbert–Pólya operator never built).
-/

namespace Fayed.R234
end Fayed.R234
