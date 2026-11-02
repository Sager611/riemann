# Record: a smaller convex cover for Moser's worm problem

[Back to DISCOVERIES](../DISCOVERIES.md)

**Claim.** There is a convex quadrilateral of area 1785486893883814994120173447337/7277367935880000000000000000000 = 0.2453478936… that contains a congruent copy (rotations, translations and reflections allowed) of every planar arc of length 1. So the least area α of such a convex set satisfies α ≤ 0.2453478936.

**Before.** α ≤ 0.246332237201 (Keller, arXiv 2609.21968, Sept 2026; a quadrilateral certified by an exact case tree). Earlier refereed: π/12 = 0.2618. Lower bound 0.239.

**Method.** Keller's quadrilateral shrunk by the factor 499/500, re-certified with a new, larger finite case tree (29,137 leaves, 1,493 placements) found by an adaptive search; the proof method is Keller's support/length argument.

**Certificate.** verify_upper.py (exact rationals, standard library, about 8 minutes, one thread) accepts the compressed certificate in cert/; it is our own re-implementation of every finite check of Keller's proof, also accepts his published certificate, and rejects tampered certificates. NOTE.md maps each case of Keller's §4 to its check.

**Review.** Correct after two fixes (the case correspondence with Keller's §4 written out; an exact winding check replacing a floating-point one); significant.

**Novelty.** Nothing below 0.246332237201 found, rechecked 2 Nov 2026 ([novelty](../Riemann/Deep/Moser/NOVELTY.md)).

**Files.** [report](../Riemann/Deep/Moser/REPORT.md), [note](../Riemann/Deep/Moser/NOTE.md), [polygon](../Riemann/Deep/Moser/polygon.json), [verify_upper.py](../Riemann/Deep/Moser/verify_upper.py), [certificate](../Riemann/Deep/Moser/cert/)
