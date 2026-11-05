/-!
# Round 253 — today's mailing and automated-conjecture hunts: nothing new (TxGraffiti harmonic counterexample already published)

* Rediscovered, already published: TxGraffiti's conjecture mu*(G) <= H(G) (arXiv 2507.17780, "open since 2023"; mu* = smallest
  maximal matching, H = harmonic index) fails for the subdivided star S(2,2,2,2,2) (11 vertices, mu* = 5, H = 100/21;
  orchestrator re-check). Already refuted by Biyikoglu and in 2606.15761 (same smallest tree counterexample).
* Conjectures the papers themselves already refute: 2209.00787 (exponential reduced Sombor; the scout's K_{1,4} case is an
  equality), 2509.12294 (extremal graphs for the diminished Sombor index differ from Movahedi et al.'s conjecture).
* Scout error: Zenodo 23934882 Eq. (11): the printed phase speed is ln 2 / (3 - 2 sqrt 2 cos(b ln 2)) (pdftotext jumbled it);
  it equals 0.43710 at b = pi/(3 ln 2), matching mpmath. No flaw.
* Verified (others' results, correct): 2611.02299 Thm 1, 2611.02295 Thm 1 (all solutions in the scanned ranges), 2611.04674
  (n <= 6), 2608.14159 Thm 2.3 (n < 500), 2610.08015 Lemma finite (p = 17, q = 13, m = 6), Sun 1710.08763 1-3-5 conjecture
  (n <= 10^4), Sun 2405.08552 Conj. 1.1 (p <= 43). OEIS A401496: closed form has a typo (108^(k-1) for 10*8^(k-1), 38^k for
  3*8^k); the recurrence fits the data.
* Not testable: 2611.04339 (theorem holds beyond an unspecified threshold).
-/

namespace Fayed.R253
end Fayed.R253
