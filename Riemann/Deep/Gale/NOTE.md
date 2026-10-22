# Gale's conjecture in R³: a 170-point set with β₄/D = 0.8906210556…

Orca record hunt, 2026-10-22. The previous record is arXiv 2610.28828 (Voronov), a 101-point set with β₄/D = 0.8892809651. This note gives a 170-point set V′ ⊂ R³ with exact rational coordinates and

  **β₄(V′)/D(V′) = √(7932058815501876/10000000211897483) = 0.890621055636…**

This is 1.340·10⁻³ above the previous record and 2.547·10⁻³ above the ball value σ₃ = √((3+√3)/6) = 0.8880738. The check is exhaustive and uses only exact integer arithmetic. Files: points.txt (the set), verify.py (stdlib only, about 1 s), cert.drat (an optional DRAT proof).

## 1. Definitions and the two monotonicity facts we use
For a finite V with diameter D, β₄(V) is the minimum, over partitions V = A₁⊔A₂⊔A₃⊔A₄, of max diam Aᵢ. For a squared distance s, G(≥s) is the graph on V whose edges are the pairs with |x−y|² ≥ s. Let dc² be the largest squared distance s for which G(≥s) is not 4-colourable. Then β₄(V) = dc:
- every 4-partition is a 4-colouring of G(≥dc²), and none is proper, so some part has diameter ≥ dc;
- a proper 4-colouring of G(>dc²) is a partition whose parts all have diameter ≤ dc.

(This is the brief's definition: dc is the smallest pairwise distance d for which G(>d) is 4-colourable.)

Two facts make the search moves safe. Write t(V) = β₄(V)/D(V).
- **(M1) Adding points inside conv V does not lower t.** D is unchanged, because a set and its convex hull have the same diameter. β₄ cannot drop, because a 4-partition of the larger set restricts to one of V.
- **(M2) Removing points outside a non-4-colourable subgraph H ⊂ G(≥s) keeps t ≥ √s / D′, with D′ ≤ D.** H survives, and the diameter can only shrink.

## 2. The paper (source: arXiv e-print, main_en.tex)
The basic method:
1. Take random points.
2. Binary search, with a SAT solver at each step, for the shortest prefix E of the pairs sorted by decreasing length such that ([s],E) has χ ≥ 5.
3. Fix E and maximise min_E |vᵢ−vⱼ| / diam V locally (L-BFGS-B or SLSQP).

Refinements in the paper: k-core reduction, a noisy re-sorting of near-critical edges (parameters δ, η), vertex splitting into twins (single and paired), and "geometric refinement" that resamples a conjectured piecewise-smooth boundary. Reported values: 0.865 (basic), 0.872 (evolved code), 0.88928 (101 points). The best 50–250-point sets all share one structure: three long equatorial arcs, two triples of short arcs near the poles, and a few points in the polar trihedral angles. Imposing symmetry, or projecting radially onto a sphere, makes things worse. The author expected that further gains "require much greater computational resources or new ideas".

## 3. Our search (scripts lived in .scratch/gale, now deleted)
- **Core.** For the critical prefix we get an UNSAT core from CaDiCaL via PySAT, with one selector literal per edge. We shrink it to a minimal non-4-colourable edge set by deletion, removing the shortest edges first.
- **Local step (SLP).** Variables: displacements δ ∈ R^{3N} and s. Each LP maximises s subject to:
  - |p_ij|² + 2p_ij·(δᵢ−δⱼ) ≥ s on every core edge;
  - |p_ij|² + 2p_ij·(δᵢ−δⱼ) + 12r² ≤ 1 on every pair with |p_ij| + 2√3 r ≥ 1;
  - ‖δ‖∞ ≤ r.

  Here p_ij = vᵢ−vⱼ and D = 1. Both linearisations are conservative, since |p+Δ|² = |p|² + 2p·Δ + |Δ|² and |Δ|² ≤ 12r². So an accepted step is a true improvement. The trust radius doubles after a successful step that reaches the box boundary and shrinks by a factor 0.3 after a failure. After each step we renormalise, recompute the critical value from scratch with SAT, and re-extract the core ("core exchange").
- **Paper's set.** Exact value 0.8892809651 (same as R139). The SLP on its core gives 0.8892868, so the 7-decimal table is not locally optimal. Core exchange and perturbations give 0.8892869. Basin hopping with the paper's moves (noisy pruning, twin splits, add/delete/jitter; about 2 CPU-h) gave no further gain. Random restarts with N = 50–100 reach at most 0.873; samples of the Reuleaux tetrahedron reach about 0.74.
- **New move: interior refinement, then prune.** Add 60–200 points, each either a midpoint (or a random 0.3–0.7 convex combination) of a pair at distance < r, or a random convex combination of a point and two neighbours within r, with r ∈ {0.03, 0.05, 0.08, 0.12}. By (M1) this changes nothing at first. But the threshold graph gains many new vertices, and core-SLP on the enlarged set finds a different, better core. We then prune to the vertices of a minimal core, which is safe by (M2), and polish again. A move is accepted only if it improves. The first application gave 0.8893313 (274 points), which pruned to 0.8893336 (100 points).
- **Chains of this move:**

  | generation | chains × time | best result |
  |---|---|---|
  | 1 | 4 × 40 min | 0.8904185 (131 pts) |
  | 2 | 5 × 60 min, N ≤ 195 | 0.8905864 (167 pts) |
  | 3 | 5 × 30 min, N ≤ 199 | 0.8906211 (170 pts) |

  Total about 10 CPU-h on a heavily loaded machine. An intermediate 118-point set had exact value 0.8902281369 (not kept).
- **Rationalisation.** We centred the set, scaled it to D ≈ 1 and rounded to the grid 10⁻⁸Z³. The exact value after rounding is 0.8906210556, a loss of 1.9·10⁻⁸ from the float optimum.

## 4. Certificate (all exact)
- **Distances.** verify.py computes all 14365 squared distances as integers (common denominator 10⁸):
  - D² = 10000000211897483·10⁻¹⁶;
  - dc² = 7932058815501876·10⁻¹⁶.
- **G(≥dc²) is not 4-colourable** (5538 edges). Two independent exhaustive searches show this:
  - forward-checking backtracking with the K4 {85,86,131,164} coloured 0,1,2,3 (allowed by colour symmetry): 29 nodes;
  - plain DSATUR backtracking without propagation, with a different K4 {74,85,86,120} fixed: 946 nodes.
- **G(>dc²) is 4-colourable** (5537 edges): an explicit proper 4-colouring is checked edge by edge, so β₄(V′) = dc exactly.
- **Comparisons**, exact: r = dc²/D² > 39541040058709/50000010510015 (the previous record), and r > (3+√3)/6, tested as 6r−3 > 0 and (6r−3)² > 3.
- **Independent checks** (tools in scratch, now deleted):
  - CaDiCaL 1.5.3 (PySAT) gives UNSAT for G(≥dc²) and SAT for G(>dc²) on the exact graphs.
  - CaDiCaL 1.9.5 produced a DRAT proof (cert.drat, 1651 lines) for the plain CNF: 680 variables, 22322 clauses, no symmetry breaking. drat-trim, built from github.com/marijnheule/drat-trim, reports "s VERIFIED" in 0.1 s.
  - To reproduce: `python3 verify.py --cnf g.cnf && drat-trim g.cnf cert.drat`.
- A float recomputation with numpy and SAT gives the same value, 0.890621055636.

## 5. What the set looks like
- 170 points. Distances from the centroid range over 0.484–0.589 D, and the circumradius is about 0.5775 D (Jung's bound is 0.6124 D), so the set is far from spherical. Some points sit in tight clusters (minimum spacing 7.9·10⁻⁴ D).
- It is a very balanced max–min configuration:
  - 244 of the 549 edges of a minimal non-4-colourable core are within 10⁻⁷ of the critical length;
  - 253 pairs are within 10⁻⁷ of the diameter;
  - the core has five hubs, of degrees 70, 69, 57, 36 and 33.
- The search had not converged. At the end, a 30-min chain still gained about 3·10⁻⁵ at N ≈ 170. Allowing more points, or more CPU time, very likely raises the value further. The method gives no evidence about the upper bound β₃,₄ ≤ 0.9654.

## 6. Status
β₃,₄ ≥ 0.8906210556 (certified). This is a lower bound from one explicit finite set, not a structural result. Nothing here is formalised in Lean.
