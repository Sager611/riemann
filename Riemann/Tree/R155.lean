/-!
# Round 155 — checks of our candidates; weekly finders; RH triage; Zenodo; Sun determinants

* EvoAlg: the quoted BCS22 theorem matches arXiv 2006.14460 (Theorem 2.4(i),(ii), TeX l. 400–411, label nat-vec; any
  field); the NOTE's use (natural vector ⇒ rank 1) is justified, and its Lemma 1 proves it independently.
* Gale record: an independent scout recomputation from points.txt alone (own C++ DSATUR, exact integers) gives
  β₄²/D² = 7932058815501876/10000000211897483, i.e. 0.890621055636; reviewer CORRECT (fixing a K₄ is a valid reduction;
  β₄ = d_c exactly; the DRAT proof passes RUP on the regenerated CNF).
* Weekly finders: 2610.29197 (local dyadic conjecture for ideal clutters) survives 9,567 instances on 2–4 elements;
  2610.25661 (factorial separation) n ≤ 100; 2610.25271 disproves its own Properties C/D. Large-q finder: arXiv API empty.
* RH triage (math.GM, math.NT, past week): nothing to audit (2610.23985, RH for unitary groups over function fields, rests on
  asserted base change). Zenodo: 10.5281/zenodo.23440770 claims Ramanujan's inequality π(x)² < (ex/log x)π(x/e) for all
  x ≥ 38,358,837,683 UNCONDITIONALLY (Dudek–Platt needed RH) → Orca audit. Also two Lean "RH proofs" (23357512,
  23334982), screw-function certificates (23377174, 23544568), a Pólya–Vinogradov formalization (23563651).
* Sun determinants: arXiv 2606.03970's conjecture is for primes p ≡ 3 (mod 4); its natural F_q analogue fails at q = 27
  (determinant 0), which the authors do not claim. Not a counterexample.
* Novelty (Orca sunnov): NEW 1309.1679 printed 3.2, 3.4(i), 3.10(ii), OEIS A228766 (general form), A228728, 1811.10503
  Conj. 4.9, MO 316143; KNOWN 3.8 (anonymous Zenodo 10.5281/zenodo.23307915, 2026-10-12, with a claimed classification),
  3.3(ii) (corrected by Sun in the homepage journal version), A229827(ii) (Sun's own Remark 1.3).
* Novelty (Orca stlnov): STLemP's mean square, lift and thresholds NEW → DISCOVERIES updated.
-/

namespace Fayed.R155
end Fayed.R155
