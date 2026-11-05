/-!
# Round 268 — Fasi-Higham off-diagonal thresholds, AI-found counterexamples, Lean repos, prime inequalities: known or correct

* Fasi-Higham off the diagonal (numerics, two precisions; the conjecture is already refuted by Scardino, 2610.39051): the
  smallest k with a non-positive coefficient of 1/q_{k,m} (k <= 400, n <= 3000) is k = 66 for m = 5 (first at n = 38),
  44 for m = 7 (n = 37), 42 for m = 9 (n = 41), 43 for m = 11 (n = 46); every larger k up to 400 fails too. Data for a later
  attempt at the positivity region (all k below these thresholds).
* Verified (others' results, correct): PatternBoost 2411.00566 (81-edge spanning subgraph of Q_6 of diameter 6), 2610.18636
  (its own refutation of Conj. 3.9-3.10, n = 4), 2511.13691 Conj. 1 (Bonse-type thresholds) and 2002.08022 Conj. 3.1-3.2
  (Robin's inequality at primorials), Zenodo 23421625 (elementary node criterion, no famous claim).
* Minor: GitHub Ghuffy1/erdos710-lean hard-codes the conjectured constant of Erdős #710 in CanonicalBridge.lean (not a proof).
* Not testable: no other Lean-linked papers from Sep-Nov 2026 beyond those already scanned.
-/

namespace Fayed.R268
end Fayed.R268
