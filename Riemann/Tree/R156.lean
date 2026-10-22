/-!
# Round 156 — checks of SunComplete; Ramanujan literature; Zenodo Lean "proofs" of RH; record finder; finders

* SunComplete novelty: the anonymous Zenodo 10.5281/zenodo.23307915 treats only Sun's Conj. 3.8 (it announces Sun–Hou
  as a "separate release"); no other 2026 record; She (arXiv 2609.19114) proves Sun–Hou for q > 18,888,871 and primes
  p > 61 only, so the prime powers 8 ≤ q ≤ 18,888,871 were open. Our completion (Theorem A) is new if correct.
  Certificate spot-check: timed out (rerun scoped).
* Ramanujan's inequality: Dudek–Platt (under RH) for x ≥ 38,358,837,682; unconditionally only x ≥ e^3915
  (Platt–Trudgian 2021), no later improvement found; RH-strength error bounds would be needed down to 3.8·10¹⁰.
* Zenodo Lean "proofs" of RH 10.5281/zenodo.23357512 and 23334982: refuted (no axioms or sorry, but the theorem is about
  a projected zero set Rec ZetaZeros; the missing identity ∀ q, Rec ZetaZeros q ↔ ZetaZeros q is equivalent to RH).
* Record finder: 2610.27175 (Hadamard matrices, 95 open orders ≤ 10,000), 2610.24025 (constant-weight codes, data
  public), 2610.29584, 2610.27133; Gale (2610.28828) is ours already.
* Finders: 2610.22663, 2610.27523 (double-cap conjecture proved there), 2610.22351, 2610.29197 (clutters, survived),
  2609.18522 (2-sumfree sequences eventually periodic: 100 starts, no failure), 2609.35252 (y² = x³ + k with small k:
  three finite exceptions up to 10⁶, no refutation).
* Sun batch reviewer: all seven certificates correct, minor; source snapshots and the printed-version comparison added
  (Riemann/Deep/SunPrim/SOURCES.md): one DISCOVERIES entry.
-/

namespace Fayed.R156
end Fayed.R156
