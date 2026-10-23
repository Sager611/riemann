/-!
# Round 171 — claimed counterexamples verified; one rests on a disconnected graph; the Catalan gap is already known

* Verified correct (others' counterexamples, our own code): 2610.06033 (Oxley's Conjecture 15.9.5: the binary matroid
  of the extended Golay code has circuit–cocircuit intersections of sizes {0, 2, 4, 6, 8, 12}, so 12 but not 10);
  2610.11642 (the Erdős–Sós bipartite-link bound n³/24: a 497-vertex 3-graph with bipartite links and 5,115,152 >
  5,115,144.71 edges); 2610.13987 (Anstee–Sali: forb(n, 1_{r−2} × C₄) = Θ(n^(r−1/2)) while X(F) = r).
* 2610.00086 ("A counterexample to Wormald's conjecture"): its 16-vertex cubic graph K_{3,3} ⊔ H has no partition into
  two isomorphic spanning linear forests (our exhaustive check of all 2²⁴ edge 2-colourings agrees), but the graph is
  disconnected, while the conjecture as quoted (Wormald, Problem 13, Ars Combin. 23A, 1987) is about connected cubic
  graphs. The original wording is checked next round.
* Known: the gap in Sun's Catalan-irrationality proof (2609.04176) was found by D. Wachs (arXiv 2609.22339,
  16 Sep 2026): the ledger omits the 2-adic term (2 log 2) B².
* Not verifiable in 50 s: 2610.03550 (strong Roberson; objects rebuilt, the exclusion theorem not rechecked); Zenodo
  23566520 (Lehmer's totient problem, ω(n) ≥ 16; arithmetic and history check out, the exhaustive search not rerun).
* 2609.10633 (Liu's triangle inequality): true on 10⁶ random and extreme triangles; the first proof step checks.
-/

namespace Fayed.R171
end Fayed.R171
