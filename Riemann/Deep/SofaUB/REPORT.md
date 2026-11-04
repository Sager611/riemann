RECORD: V*(K_M) ≤ 3.1 for M > 2.2 (was 3.4083)

- What: every connected body that passes through the 3D snaky corridor K_M has volume ≤ 3.0999999989, for every
  middle length M > 2.09 (so also for M > 2.2 and the problem's M = 4). Previous upper bound 3.4083 (arXiv
  2610.14050, Oct 2026), computed in plain double precision. Lower bound stays 1.8755.
- Same reduction as the paper (a "first crossing" moment pins the body inside a tilted tube that meets two unit
  windows), but a new, exact way to bound tube volumes: slice the tube horizontally; each slice is a parallelogram,
  so the volume is an exactly integrable piecewise-quadratic function of height.
- Finding: the paper's tube family can never give less than 1+u (a vertical tube fits the whole stub). At their
  u = 2.2 the floor is 3.2. Only one other configuration competes (a shallow tube through the entry arm); the best
  any certificate of this kind can reach is about 3.073. Ours (3.1, at u = 2.09) is within 0.03 of that limit.
- Rigor: interval arithmetic in double precision with outward widening after every operation (rigorous enclosures),
  rotations parametrized rationally (no trigonometry), exact breakpoint handling. Not just a written error bound.
- Certificate: branch-and-bound tree with 16.1 million leaf boxes (1.5 MB compressed). verify.py replays it
  (6 min, 89 MB, one core) and cross-checks sampled leaves with exact fractions, Monte Carlo and an exact
  coverage proof. All pass.
- Reproduced the paper's own certificate from its archived files: 3.408292 exactly (needs a numpy-2 shim).
- Novelty: no newer upper bound in arXiv Oct-Nov 2026 listings or Zenodo.
- Next: going below ~3.07 needs coupling several moments of the motion (the vertical "rod" cannot turn the corner).
- Files: Riemann/Deep/SofaUB/ (NOTE.md, verify.py, sofaub.c, exact_check.py, mc_check.py, meta.json, codes.bin.gz).
- Check: cd Riemann/Deep/SofaUB && python3 verify.py
