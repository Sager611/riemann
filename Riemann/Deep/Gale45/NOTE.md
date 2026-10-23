# Gale's conjecture in R⁴ and R⁵: stronger counterexamples (Orca record hunt, 2026-10-23)

Previous records: the 61-point set in R⁴ and the 87-point set in R⁵ of arXiv 2610.28828 (Voronov, v1, 20 Oct 2026, Table 2) and its repository github.com/vsvor/Gale-conjecture (`best_known/`). This note gives a 126-point set V₄ ⊂ R⁴ and a 149-point set V₅ ⊂ R⁵ with exact rational coordinates and

  **β₅(V₄)/D(V₄) = √(8399124363373578/10000000179265409) = 0.916467359637…**  (old 0.913148536990, ball σ₄ = 0.912870929)

  **β₆(V₅)/D(V₅) = √(2154677368122705/2500000047209683) = 0.928370039894…**  (old 0.925763634680, ball σ₅ = 0.923879533)

That is +3.319·10⁻³ and +2.606·10⁻³ over the old records, and +3.596·10⁻³ and +4.491·10⁻³ over the ball values. Files: points4.txt, points5.txt (the sets), verify.py (stdlib only, exact integers, about 15 s for both), cert4.drat, cert5.drat (optional DRAT proofs). The method and the R³ record (170 points, 0.8906210556) are in ../Gale/NOTE.md.

## 1. Definitions
For a finite V ⊂ R^d with diameter D and k = d+1, β_k(V) is the minimum over partitions V = A₁⊔…⊔A_k of max diam Aᵢ. G(≥s) is the graph on V joining pairs with |x−y|² ≥ s. If dc² is the largest squared distance with G(≥dc²) not k-colourable, then β_k(V) = dc: every k-partition is a non-proper colouring of G(≥dc²), and a proper k-colouring of G(>dc²) is a partition with all parts of diameter ≤ dc. Adding points inside conv V keeps D and cannot lower β_k; deleting points outside a non-k-colourable subgraph of G(≥s) keeps β_k ≥ √s while D can only shrink. Both facts make the search moves below safe.

## 2. What the repository numbers are (task 1)
- **The quantity.** `validation/example_verification.json` and `docs/history/base-search.md` report R_k(X) = β_k(X)/diam X with k = d+1 parts: a lower bound from an UNSAT witness graph and an upper bound from an explicit partition, both checked in integer arithmetic, bracketed to 12 decimals. The repository itself says "Graph UNSAT results are not independently proof-logged".
- **Ball value.** Gale (1953) conjectured that the ball is extremal: β_{n+1}(V) ≤ β_{n+1}(B^n)·D(V) for every V ⊂ Rⁿ. The paper's §1 gives the regular-simplex partition bound β_{n+1}(B^n) ≤ σ_n, with σ_n² = (n+1)/(n+2) for even n and σ_n² = 1/2 + ½√((n−1)/(n+3)) for odd n. So σ₄ = √(5/6) = 0.9128709292 and σ₅ = √(1/2 + √2/4) = cos(π/8) = 0.9238795325. The true ball values are known only for n = 2, 3; any V with β_{n+1}(V)/D(V) > σ_n is a counterexample.
- **examples/d4_m61** (61 points, 0.913148536990) has the same exact witness integers as `best_known/4d_61p_091314` (the `examples/` directory itself is not shipped in the repository zip; `docs/assembly.md` says d4_m61 "lacked a standard candidate directory"). It beats σ₄ by 2.8·10⁻⁴, so it is a counterexample in R⁴.
- **examples/d5_m32** (32 points, 0.892304779206) and examples/d4_m32 (0.887138937436) are "small smoke tests" of the basic search (docs/history/base-search.md). They are below σ₅ = 0.92388 and σ₄ = 0.91287, so they are **not** counterexamples. The R⁵ counterexample is `best_known/5d_87p_092576` (87 points, 0.925763634681, 1.9·10⁻³ above σ₅).
- **Where posted.** Only in the paper's TeX (Table 2, "Current records": n = 4: 0.91314…, 61 points; n = 5: 0.92576…, 87 points; n = 6: 0.93696…, 64 points; "most likely far from optimal") and in the repository. Checked 2026-10-23: arXiv abs page (v1 only), arXiv math.MG and math.CO listings for 2026-10 (no other Gale paper), Zenodo (nothing), the GXL evidence index (no Gale-type result in R⁴ or R⁵), repository commits (unchanged since bd6dc959, 2026-10-22 12:55 UTC). A wider adversarial check is in NOVELTY.md.

## 3. Search (task 2)
Same machinery as in R³ (../Gale/NOTE.md §3), with k = d+1:
- **Core-SLP polish.** SAT (CaDiCaL via PySAT) finds the critical prefix of pairs sorted by length; a deletion-minimal non-k-colourable core is extracted; LP steps (HiGHS) maximise the shortest core distance with all near-diameter pairs kept ≤ D (linearised, trust region); after a few steps the critical value is recomputed exactly and the core exchanged.
- **Interior refinement, then prune.** Add 30–140 convex combinations of nearby points (midpoints, random 0.3–0.7 combinations, random triples within radius 0.05–0.3). This changes neither D nor β at first, but core-SLP on the enlarged set finds a better core. Then cut back to the vertices of the final minimal core and polish again. Accept only on improvement.
- **Outward variant.** Half of the new points are pushed out along the ray from the centroid, by up to 3 %, as far as every distance stays ≤ D (then D is unchanged). With radii up to 0.3 this variant broke a long stall in R⁴: 0.91471 → 0.91647 in three chains of 40–50 minutes, in bursts between idle generations.
- **First observation.** Polishing the published sets alone: on the R⁴ 61-point set core-SLP finds no improving step (0.9131485); the R⁵ 87-point set rises to 0.92630 (305 accepted steps), so the published R⁵ set was not even locally optimal.
- **Chains** (single-threaded, under heavy machine load; about 8.6 chain-hours in R⁴ and 6.6 in R⁵). Best values along the way:

  | R⁴ (points) | R⁵ (points) |
  |---|---|
  | 0.9131485 (61, paper) | 0.9257636 (87, paper) |
  | 0.9134060 (71) | 0.9267090 (107) |
  | 0.9145067 (79) | 0.9272544 (113) |
  | 0.9147087 (93) | 0.9276714 (122) |
  | 0.9153830 (97) | 0.9279318 (135) |
  | 0.9158143 (100) | 0.9280423 (138) |
  | 0.9162111 (120) | 0.9282348 (144) |
  | 0.9164674 (126) | 0.9283701 (149) |

- **Rationalisation.** Centre, scale to D ≈ 1, round to the grid 10⁻⁸ Z^d; the exact value after rounding is within 2·10⁻⁸ of the float optimum.

## 4. Certificate (task 3, all exact)
`python3 verify.py` (stdlib only) checks both files:
- all squared distances as integers (common denominator 10⁸); binary search for dc² over the distinct squared distances;
- **R⁴, 126 points:** D² = 10000000179265409·10⁻¹⁶, dc² = 8399124363373578·10⁻¹⁶. G(≥dc²) (3197 edges) is not 5-colourable: forward-checking backtracking with the K4 {33,43,47,94} fixed (181 nodes) and, independently, plain DSATUR with the K4 {7,30,43,47} fixed (3245 nodes). G(>dc²) (3196 edges) has an explicit proper 5-colouring, checked edge by edge.
- **R⁵, 149 points:** D² = 10000000188838732·10⁻¹⁶, dc² = 8618709472490820·10⁻¹⁶. G(≥dc²) (4863 edges) is not 6-colourable: forward checking with the K5 {9,36,51,55,71} fixed (11877 nodes), plain DSATUR with the K4 {1,7,66,114} fixed (83620 nodes). G(>dc²) (4862 edges) has an explicit proper 6-colouring.
- exact comparisons with σ_d (even d: r > (d+1)/(d+2); odd d: 2r−1 > 0 and (2r−1)² > (d−1)/(d+3)) and with the previous records (recomputed exactly from the repository's points.csv).
- **Independent check.** `python3 verify.py points4.txt --cnf g4.cnf --sym` writes the CNF "G(≥dc²) is 5-colourable" (630 variables, 16115 clauses) with unit clauses fixing the K4's colours (allowed by colour symmetry). CaDiCaL 1.5.3 (PySAT) proves it UNSAT; its DRAT proof cert4.drat is "s VERIFIED" by drat-trim. The same for R⁵ (894 variables, 29332 clauses, cert5.drat).

## 5. Caveats
- Lower bounds β₅/D ≥ 0.9164673596 (R⁴) and β₆/D ≥ 0.9283700398 (R⁵) from explicit sets. The searches had not converged: both landscapes alternate long stalls with sudden jumps. No structure theory, nothing in Lean.
- The comparison is with the records in arXiv 2610.28828 and its repository as of 2026-10-23; R⁶ (0.93696, 64 points) was not attempted.
