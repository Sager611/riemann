RECORD: α ≤ 1785486893883814994120173447337/7277367935880000000000000000000 = 0.2453478935812… (was 0.246332237201…)

Setting (Keller, arXiv 2609.21968): α = least area of a convex set containing a congruent copy
(rotations, translations and reflections allowed) of every planar arc of length 1.
Cover: Keller's quadrilateral shrunk by 499/500; side directions t = (-631/58, -1, 7/16, 37/56),
offsets h = (0, 0, 65260717/125000000, 222653301/500000000); area = Keller's area * (499/500)^2,
0.40 % smaller. Exact vertices: polygon.json.

Proof: Keller's support/length method (his Section 4) with a new finite case tree: 29 137 leaves,
1 493 orientations (Keller's 1 463 + 30 found by an adaptive search), exact integer leaf multipliers.
Certificate: cert/ (13.6 MB, compact format). NOTE.md, "Correspondence with Keller §4", maps each
step of his proof (TeX line numbers) to the check that covers it and re-derives why the case tree
and the four prefix families cover every arc.

Checks (exact rational arithmetic only, Python standard library):
- `python3 verify_upper.py cert` -> ACCEPT with the area above (2026-10-23 10:15 UTC, ~10 min on one
  low-power core, 70 MB). Its only floating-point step, a redundant atan2 winding test, was replaced
  that morning by an exact count; floats remain only in the printed decimal. Earlier versions, which
  differ only in that test, accepted Keller's own certificate and rejected 7 tampered ones.
- `python3 search/crosscheck_exact.py <keller-repo>/scripts cert` -> all 29 137 leaves also pass
  Keller's own exact leaf test (scripts/upper/exact.py, Geometry.valid, offsets patched).
- Not done: Lean kernel check (brief: no Lean builds).

Novelty (2026-10-23: arXiv API, math.MG/math.CO October listings, Zenodo): no newer upper bound.
Not a bound for Tao et al.'s constant 13a (no reflections), whose best upper bound is still π/12.

Room left: the hardest arcs found for this shape need scale >= 0.99615 (area ~0.2444). Scale 0.997
was stopped (35 900 open cells, >= 30 MB certificate); sparse leaf multipliers (~4x smaller
certificates) are the next step. Other quadrilaterals and pentagons gave no better estimate.
