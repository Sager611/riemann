# Gale's conjecture in R⁶: a 99-point set with β₇/D = 0.9373384468…

Orca record hunt, 2026-10-23 to 2026-11-02. The previous record is the 64-point set `best_known/6d_64p_093696` of arXiv 2610.28828 (Voronov, v1, Table 2) and its repository github.com/vsvor/Gale-conjecture. This note gives a 99-point set V₆ ⊂ R⁶ with exact rational coordinates and

  **β₇(V₆)/D(V₆) = √(8786033541902487/9999999889285618) = 0.937338446836…**  (old 0.936967247964, ball σ₆ = √(7/8) = 0.935414347)

That is +3.712·10⁻⁴ over the old record and +1.924·10⁻³ over the ball value. Files: points6.txt (the set), verify.py (Python standard library only, exact integers, about 2 s), cert6.drat (DRAT proof, 185 KB), points56.txt (a symmetric 56-point counterexample, §4). Method and earlier records: ../Gale/NOTE.md (R³), ../Gale45/NOTE.md (R⁴, R⁵).

## 1. Definitions
For a finite V ⊂ R⁶ with diameter D, β₇(V) is the least possible largest part diameter over partitions of V into 7 parts. G(≥s) joins the pairs with squared distance ≥ s. If dc² is the largest squared distance with G(≥dc²) not 7-colourable, then β₇(V) = dc: every 7-partition is a non-proper colouring of G(≥dc²), and a proper 7-colouring of G(>dc²) is a partition whose parts all have diameter ≤ dc. Gale (1953) conjectured that the ball is extremal; the regular-simplex partition gives β₇(B⁶) ≤ σ₆·D with σ₆² = 7/8, so any V with β₇/D > σ₆ is a counterexample.

## 2. The old record, recomputed exactly (task 1)
- points.csv (17-digit decimals, read as exact rationals): D² = 10000000000000002601437547173996163·10⁻³⁴, dc² = 8779076237580130897884591797000439·10⁻³⁴, so β₇/D = 0.936967247964416770… G(≥dc²) has 1139 edges and is not 7-colourable (CaDiCaL), G(>dc²) is 7-colourable. Both integers and the edge count agree with the repository's certificate.json (which notes that its UNSAT results are not proof-logged).
- Unchanged on 2026-11-02: arXiv abs page (v1 only), repository (last commit bd6dc959, 2026-10-22).

## 3. Search (task 2)
- **Structure of the old set.** It is a perturbed "simplex product", the family of the repository's simplex_product_seed.py: R⁶ = A ⊕ B with two orthogonal 3-spaces, poles at the vertices of a regular tetrahedron in A and one in B (8 poles, pole–pole distance ≈ D inside a tetrahedron), 4 points around each pole (offsets ≈ 0.1 D in the other 3-space), 16 single points near the midpoints of the A-pole/B-pole pairs, and 8 extra points.
- **Core-SLP polish** (as in R³–R⁵: SAT critical value, deletion-minimal non-7-colourable core, LP steps that lengthen the shortest core edges while keeping all pairs ≤ D, core exchange). The old set alone rises to 0.9372106; it was not locally optimal.
- **Interior refinement chains** (add convex combinations, outward pushes and twins, polish, prune to the core): 0.9372223 (66 points), then every LP stalled, also with cores biased toward the new points. Random kicks of size 10⁻³ D: +2·10⁻⁸.
- **What worked: SAT-screened vertex splitting** (the paper's twin-splitting idea, combined with our core-SLP). For a core vertex v and two of its critical core edges e = va, f = vb, replace v by twins at the same position: v keeps all its core edges except e, the new twin v′ all except f. If the split core is still not 7-colourable (one SAT call), one LP step can lengthen e and f independently. Positive-gain splits are accepted; then polish and prune. In later rounds up to 6 positive splits on distinct vertices are combined (re-checked by SAT) before one LP step. Every accepted step is safe: the split core is a non-7-colourable subgraph of G(≥s) for the new set, and the conservative LP keeps all distances ≤ D.

  | step | points | β₇/D |
  |---|---|---|
  | paper, best_known/6d_64p_093696 | 64 | 0.9369672480 |
  | core-SLP polish | 64 | 0.9372106014 |
  | refinement chain, kick | 66 | 0.9372225445 |
  | single splits (9 rounds) | 75 | 0.9372395847 |
  | multi-splits (4 rounds) | 95 | 0.9373313334 |
  | multi-splits (1 round) | 99 | 0.9373384441 |
  | rounded to the grid 10⁻⁸ Z⁶, exact | 99 | 0.9373384468 |

  About 3 CPU-hours in total, one single-threaded process at a time, on a loaded machine.
- **Dead ends** (all float values):
  - exact T_d × T_d symmetry (every point x·uᵢ ⊕ y·vⱼ for a few profile points (x, y), uᵢ and vⱼ the tetrahedron directions): best 0.9368472 with 56 points, and this optimum is rigid (free polish does not move it);
  - random symmetric structures with other splits and orbits (about 1100 short trials): ≤ 0.93685 for A ⊕ B = R³ ⊕ R³, 0.93234 for R² ⊕ R⁴, 0.91137 for R¹ ⊕ R⁵; adding orbits along the 2-fold axes did not help;
  - random restarts from the noisy symmetric set plus twins: 0.919–0.933.

## 4. Certificate (task 3, all exact)
`python3 verify.py` (stdlib only, about 2 s) checks points6.txt:
- all 4851 squared distances as integers (common denominator 10⁸); D² = 9999999889285618·10⁻¹⁶, dc² = 8786033541902487·10⁻¹⁶;
- G(≥dc²) (2439 edges) is not 7-colourable: forward-checking backtracking with the K5 {1,25,27,33,39} fixed (1419 nodes) and, independently, plain DSATUR with the K5 {20,25,29,43,63} fixed (16351 nodes);
- G(>dc²) (2438 edges) has an explicit proper 7-colouring, checked edge by edge, so β₇(V₆) = dc exactly;
- exact comparisons: dc²/D² > 7/8 and > the old record's fraction (§2).
- **Independent check.** `python3 verify.py points6.txt --cnf g6.cnf --sym && drat-trim g6.cnf cert6.drat`: the CNF "G(≥dc²) is 7-colourable" (693 variables, 17177 clauses, unit clauses fix the K5's colours, allowed by colour symmetry) is UNSAT by CaDiCaL 1.5.3 (PySAT); its DRAT proof cert6.drat (12514 lines) is "s VERIFIED" by drat-trim.
- **Side result: a symmetric 56-point counterexample.** points56.txt is the T_d × T_d-symmetric optimum of §3 rounded to 10⁻⁸ Z⁶: poles ρuᵢ and ρvⱼ with ρ = √(3/8) (pole tetrahedra of edge D), the 16 points αuᵢ − βvⱼ, their 16 mirror images −βuᵢ + αvⱼ, and the 16 points γ(uᵢ + vⱼ), with α = 0.59986445, β = 0.10724233 (α + β = 1/√2, these pairs are diameters), γ = 0.42771320. `python3 verify.py points56.txt` certifies β₇/D = √(2925609174002721/3333333325758482) = 0.9368472417 > σ₆ (2 s; no DRAT file). It is below the old record but uses fewer points (56 against 64) and is described by three numbers.

## 5. What the set looks like
99 points, all in the final minimal core (874 edges, five hubs of degrees 57, 56, 50, 49, 48); distances from the centroid 0.530–0.637 D; minimum spacing 8.6·10⁻⁴ D (several split twins stay close). The simplex-product skeleton of §3 is still visible; the split twins thicken the pole patches.

## 6. Status and caveats
β₆,₇ ≥ 0.9373384468 (certified), a lower bound from one explicit set; no structural result. The search had not converged: the last multi-split rounds still gained 0.7–2.6·10⁻⁵ each. The same splitting step may also help the R³–R⁵ records (not tried). Nothing is formalised in Lean.

