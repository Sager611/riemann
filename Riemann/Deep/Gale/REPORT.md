RECORD: β₄/D ≥ 0.8906210556 for an explicit 170-point set in R³ (certified by exact integer arithmetic: two exhaustive backtracking searches, plus a CaDiCaL DRAT proof checked by drat-trim)

# Gale's conjecture in R³: stronger counterexample (2026-10-22)
Previous record: the 101-point set of arXiv 2610.28828, with exact β₄/D = √(39541040058709/50000010510015) = 0.8892809651 (R139). The ball value is σ₃ = √((3+√3)/6) = 0.8880738. Details are in [NOTE.md](NOTE.md).

## Result
- **The set.** [points.txt](points.txt) holds 170 points with coordinates in 10⁻⁸Z³.
- **The value.** β₄(V′)/D(V′) = √(7932058815501876/10000000211897483) = 0.890621055636…, exactly. That is +1.340·10⁻³ over the previous record and +2.547·10⁻³ over the ball.
- **The check.** [verify.py](verify.py) uses only the Python standard library and exact integers, and runs in about 1 s. It does four things:
  - finds the critical squared distance dc² by binary search;
  - proves G(≥dc²) (5538 edges) has no proper 4-colouring with two exhaustive searches: forward checking with a K4 fixed (29 nodes), and plain DSATUR with a different K4 fixed (946 nodes);
  - checks an explicit 4-colouring of G(>dc²), so β₄ = dc exactly;
  - compares the value exactly with σ₃ and with the old record.
- **Independent check.** CaDiCaL gives UNSAT/SAT on the two graphs. A DRAT proof [cert.drat](cert.drat) for the plain CNF (no symmetry breaking) is "s VERIFIED" by drat-trim. To redo it: `python3 verify.py --cnf g.cnf; drat-trim g.cnf cert.drat`.

## How it was found
- **Paper's method reproduced.** SAT search for the critical non-4-colourable edge prefix, then local max–min optimisation with that edge set fixed. Our local step is a sequential LP with conservative linearisation on a minimal UNSAT core, re-extracting the core after each step.
  - On the paper's set this reaches only 0.8892869.
  - The paper's other moves (noisy pruning, twin splitting, add/delete/jitter; about 2 CPU-h) gave nothing more. Random restarts reach at most 0.873.
- **New move: interior refinement, then prune.** Add 60–200 convex combinations of nearby points. This leaves D and the current critical value unchanged. Re-optimising on the enlarged set finds a better core. Then cut back to that core's vertices (the ratio cannot drop) and polish.
- **Progress.** Iterating the move (14 chains, about 10 CPU-h) gave 0.88933 → 0.89023 (118 pts) → 0.89042 (131) → 0.89059 (167) → 0.89062 (170).
- **Rounding.** Rounding to the 10⁻⁸ grid lost 1.9·10⁻⁸.

## Caveats
- The search had not converged: at N ≈ 170 a 30-min chain still gained about 3·10⁻⁵. More points or CPU time would likely give more.
- This is a lower bound β₃,₄ ≥ 0.8906210556 from one explicit set. There is no structural explanation, and nothing is in Lean.
- The set is far from spherical: distances from the centroid range over 0.484–0.589 D. It is highly balanced: 244 core edges sit at the critical length and 253 pairs at the diameter, each within 10⁻⁷.
