# Records: stronger counterexamples to Gale's conjecture in R⁴ and R⁵

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** Explicit point sets with β₅(V)/D(V) = 0.9164673596… in R⁴ (126 points) and β₆(V)/D(V) = 0.9283700398… in R⁵ (149 points), exact values (square roots of rationals): every partition of V into d + 1 parts has a part of diameter at least that fraction of D(V).

**Before.** 0.913148536990 (61 points, R⁴) and 0.925763634681 (87 points, R⁵): Voronov, arXiv 2610.28828 (20 Oct 2026) and its repository. Gale's 1953 conjecture says the ball is extremal; the regular-simplex bounds for the ball are σ₄ = 0.9128709 and σ₅ = 0.9238795, so both old and new sets are counterexamples, and ours are stronger ones.

**Method.** Our R³ method (interior refinement, core-SLP polishing of the critical non-colourable core, pruning) run as single-threaded chains; coordinates rounded to an integer grid.

**Certificate.** verify.py (13 s, standard library): exact integer distances; exhaustive searches that G(≥ dc²) is not (d+1)-colourable and an explicit (d+1)-colouring of G(> dc²); DRAT proofs are included as optional extras.

**Review.** Correct; significant as stronger counterexamples.

**Novelty.** Nothing above the 2610.28828 values found, checked 23 Oct 2026 ([novelty](../Riemann/Deep/Gale45/NOVELTY.md)).

**Files.** [note](../Riemann/Deep/Gale45/NOTE.md), [report](../Riemann/Deep/Gale45/REPORT.md), [R⁴ points](../Riemann/Deep/Gale45/points4.txt), [R⁵ points](../Riemann/Deep/Gale45/points5.txt), [verify.py](../Riemann/Deep/Gale45/verify.py)
