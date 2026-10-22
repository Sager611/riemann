/-!
# Round 146 — Sun's conjecture: original wording and verdict; evolution algebras; four more hunts; STLemP review

* Sun's circular-permutation conjecture (R145 counterexample q = 11, a₀ = 2). The original is Conjecture 3.10 of
  Z.-W. Sun, arXiv 1309.1679 (TeX l. 629–633): every q > 7 and every a₀, no exceptions. 2610.06289 cites it as
  "Conjecture 3.8" (a numbering slip: 3.8 there is about quadratic residues). Neither source checks small q.
  Reviewer: the counterexample is CORRECT, significant if new. Prime powers and q ≤ 199: search timed out (rerun).
* Evolution algebras (2610.22248, l. 1351–1354): the problem names no field. The known non-singular case
  (Cadavid–Rodiño Montoya–Rodriguez 2020/21) is over R and uses positive real cube roots (their l. 994–1002); they
  conjectured the singular case from examples (their l. 961–965). The Q failure (P3) is noted nowhere we found.
  Over R the open part is the singular graphs. Independent Gröbner check: timed out (rerun).
* No counterexample: 2610.24758 necklace Hamiltonicity (n = 5, 7, 9, every ℓ); 2610.06527 induced-saturation
  dichotomy (all posets on ≤ 4 elements, n ≤ 4); 2610.07442 s₃₁₂(r,k,d) closed form (all r, k ≤ 4); 2610.01747
  t-intersecting partitions (k ≤ 4, n ≤ 16; the only failure, k = 3, n = 8, lies below n ≥ t·k²).
* STLemP (our candidate, Riemann/Deep/STLemP): reviewer INCOMPLETE with a repairable gap. Lemma W(ii) (NOTE l. 152–153,
  201–203) misses the k = 1 AFE block when s₀ ≥ 2P (N ≤ 1/2 makes log(2N) ≤ 0); fix by bounding k = 1 trivially or
  N = max(1, 2^⌈log₂(P/s₀)⌉). Lemma L, the class-𝔅 checks, the Theorem MS bookkeeping and the thresholds 1,
  (23 − 3√5)/8 and 3 check. Also NOTE l. 365 (an extra h) and the REPORT's hypotheses (weight ≥ 2, trivial
  nebentypus, minorant with positive Sato–Tate mean). → Orca repair.
-/

namespace Fayed.R146
end Fayed.R146
