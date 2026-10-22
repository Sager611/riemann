/-!
# Round 159 — independent checks of our candidates; Zenodo triage; Kéri tables; finders; the Ramanujan audit

* Orca ramaudit: REFUTATION of Zenodo 10.5281/zenodo.23440770 v2 (Padhi, 2026-10-17, "Ramanujan's inequality holds for all
  x ≥ 38,358,837,683", unconditionally) as published: in the far off-diagonal term of its new zero-density bound (Lemma 6.4 with
  Lemma 6.2(v)) zeros are counted on one side only, losing a factor 2 (copied from Bellotti, arXiv 2405.12545, (4.13));
  corrected, the paper's own certified pipeline gives q = 1.03–1.07 > 1 near log x = 875 (published 0.99137). A repair
  (|X^w − Y^w| ≤ (1 + (Y/X)^{2ν})X^{Re w}) gives q ≤ 0.8957 on [57, 3159): plausibly true, other inputs unchecked → Orca
  ramrepair. Side finding: Bellotti (4.13), (4.15) share the two flaws. Best proved threshold stays e^3158.442 (Axler 2022).
  Pending review.
* SunComplete: an independent exact recheck of Lemma 6 (uniform coverage for q ≥ 30,140, all ω(q − 1)) finds no discrepancy.
* Zenodo 10.5281/zenodo.23307915 (Sun's Conj. 3.8 classified): character-sum sieve + Chvátal–Erdős down to 510 fields
  (largest 92,821), then heuristic cycles checked by a deterministic verifier; not independently certified; no Sun–Hou.
* Kéri's tables (2009): K₄(6,2) ∈ [32, 52], K₄(8,2) ∈ [240, 352]; 2610.03760 gives ≤ 48 and [251, 256].
* Finders: 2609.02041, 2609.10234 (log-concavity, m ≤ 100), 2609.17561: no counterexample. math.GM: 2610.10643 (perfect
  cuboid "rank obstruction") imports a stronger specialization of Ingram's Lemma 7 without proof (medium) → check.
* Timed out: Sun–Hou independent search (q ≤ 600), stabbing-squares experiment, EvoAlg on 8 vertices.
-/

namespace Fayed.R159
end Fayed.R159
