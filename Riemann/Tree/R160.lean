/-!
# Round 160 — checks around the Ramanujan audit; Bellotti; perfect cuboid; today's mailing

* Ramanujan audit (Padhi, Zenodo 23440770): reviewer NOT YET CONFIRMED: the ordered-vs-unordered counting objection is valid in
  principle, but our NOTE omits the full original summation; the repair inequality |X^w − Y^w| ≤ (1 + (Y/X)^{2ν})X^{Re w} is
  correct; one audit script has a bug (scanmod.py). Passed to Orca ramrepair.
* Bellotti, arXiv 2405.12545 (explicit zero density for ζ), read independently: (a) REAL FLAW at l. 855–872: the rectangles
  (l. 648–658) are one-sided, so the far off-diagonal sum needs a factor 2 (c₄ → 2c₄); (b) REAL FLAW at l. 933:
  |G(½ + it)| ≤ Σ θ_j θ_k [j,k]^{−1/2} drops signs. Both feed c₅, C₁, C₂, J and the table constants (l. 1060–1115), which are
  unsupported until re-derived. → after ramrepair, a dedicated check of Bellotti's stated constants.
* Literature: unconditionally Ramanujan's inequality is known for x ≤ e^103 down to 38,358,837,683 (Johnston 2022) and for
  x ≥ e^3158.442 (Fiori–Kadiri–Swidinsky, JNT 2021; reported by Axler 2022); the gap between is genuine.
* arXiv 2610.10643 (perfect cuboids force rank ≥ 2): no gap (Ingram's Lemma 7 gives good primitive divisors for every m ≥ 3);
  it does not claim non-existence.
* Sun–Hou, independent search: Hamilton cycles for 50 of the 58 primes 11 ≤ q ≤ 300 (8 timed out), none fail; q = 3, 5, 7 fail.
  EvoAlg 8-vertex numerics: test design flawed (determinant normalisation), inconclusive. Stabbing squares: the TeX has no
  instance data (only the Zenodo artifact).
* Today's mailing, audit targets: 2610.28874 (sporadic simplicial arrangement of 51 lines vs the Cuntz–Mücksch bound),
  2610.03229 (counterexamples to a published Hegselmann–Krause theorem), 2610.29043, 2610.29048, 2610.29115, 2608.24930
  (Erdős–Graham constant c = 8), 2608.30798, 2610.28774.
-/

namespace Fayed.R160
end Fayed.R160
