RECORD: V*(K₄) ≥ 1.8755015 (was 1.8558404)

- What: an explicit connected body of volume ≥ 1.8755015 moves continuously through the 3D snaky corridor K₄
  (arXiv 2610.14051 / 2511.02864 Problem 6.64). Gain +0.0197 over the October 2026 record.
- Motion: the paper's own public run-II path (w5_night2), refined at the paper's tolerances (0.005, 0.05°);
  between poses translation and rotation run together (not the pipeline's x-y-z-then-rotate staircase).
  Orientations are stored as quaternions, so the rotations are exact rational matrices (no trigonometry).
- Body: 195,220 slanted prisms over a quadtree of the cross-section (cells down to 1/3200), whose end faces
  follow the tilted surfaces of the carved body.
- Rigor: each prism is checked against the whole continuous motion (chord lemma, exact slice test). No
  points-to-cubes counting rule (the old record's caveat). Connectivity checked; only the largest piece is claimed.
- Rounding (reviewer's gap, now closed): exact fractions for rotations, angles and the volume; outward-rounded
  intervals for geometry, start/end and connectivity; a written worst-case bound for the C core: total error
  ≤ 1.96e-14 against the 5e-10 margin (25,000× room), computed and enforced by verify.py itself.
- Verified: all 195,220 prisms pass, 90 s, 168 MB, one thread. Controls: 1,497 of 1,500 pushed-out prisms
  rejected; a duplicated prism rejected; a too-small margin rejected; Monte Carlo of the true motion: 0 escapes.
- Novelty: arXiv search (Oct 2026) finds nothing newer than 2610.14051 for this problem.
- Side findings: the published paths are tuned to the pipeline's sampling lattice; their true volumes are
  w5_night2 1.8811, w6_night3 1.8815, w5_night1 1.8782. The record path w6_night3 has slivers on lattice planes.
- Headroom: ~0.006 left for this path; optimizing a path against an unbiased objective should give more.
- Assumes only IEEE-754 double arithmetic (numpy, C with -ffp-contract=off) and Python's correct float rounding.
- Files (6.6 MB): NOTE.md, w5_night2_path.json, certificate_w5_night2.npz, verify.py, verify_core.c, build/.
- Check: python3 verify.py w5_night2_path.json certificate_w5_night2.npz
