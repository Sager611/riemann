# Record: a larger body that moves through the 3D snaky corridor K₄

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** V*(K₄) ≥ 1.8755015: an explicit connected body of that volume moves continuously through the three-dimensional "snaky corridor" K₄ (the 3D moving-sofa problem of arXiv 2610.14051; Problem 6.64 of arXiv 2511.02864).

**Before.** V*(K₄) ≥ 1.8558404 (arXiv 2610.14051, Oct 2026); earlier 1.81. Larger numbers in the literature (about 1.88–1.89) are uncertified sampling estimates.

**Method.** The paper's own public motion (run w5_night2), refined at its tolerances, with translation and rotation run together between poses. The body is 195,220 sheared prisms over a quadtree of the cross-section, each checked against the whole continuous motion.

**Certificate.** verify.py with verify_core.c (98 s, one thread): every prism stays inside the corridor along the motion (chord lemma with an exact rotation-sag bound, exact slice test), the body is connected, and the volume is summed exactly. Rounding is bounded: exact or outward-rounded arithmetic, plus a worst-case bound of 1.9·10⁻¹⁴ for the floating-point core against a 5·10⁻¹⁰ margin.

**Review.** Correct and significant; an earlier gap (rounding not bounded) is closed.

**Novelty.** No certified bound above 1.8558404 found anywhere ([novelty](../Riemann/Deep/Sofa3D/NOVELTY.md)).

**Files.** [report](../Riemann/Deep/Sofa3D/REPORT.md), [note](../Riemann/Deep/Sofa3D/NOTE.md), [verify.py](../Riemann/Deep/Sofa3D/verify.py), [verify_core.c](../Riemann/Deep/Sofa3D/verify_core.c), [path](../Riemann/Deep/Sofa3D/w5_night2_path.json)
