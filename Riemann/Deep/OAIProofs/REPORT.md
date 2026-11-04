NO FLAW FOUND: OpenAI's proofs of Seymour's second-neighbourhood conjecture, Saxl's conjecture and the tensor-square conjecture (plus the Hill and Zarankiewicz crossing-number proofs): every step read; lemmas tested on small cases; all hold.

# Audit: OpenAI's claimed proofs of famous combinatorial conjectures (agent "oaiproofs", Nov 2026)
Source: github.com/openai/math, preprints/<name>/build (TeX, single-file downloads). Method: every lemma that has a finite instance was retested with our own code, often in a form stronger than the paper needs (the pruning lemma for every binary relation, cyclic-module supports for the exact generator, the crossing bound for every normalized drawing, not only minimum ones). Exact integer or rational arithmetic, one single-threaded job at a time. Scripts are in this folder; their output lines are in OUTPUT.txt.

| preprint | lemmas tested | ranges | verdict |
|---|---|---|---|
| Seymour 2nd neighbourhood (23 Sep) | Pruning lemma: matching bound nu <= abs(Z)+abs(H), Konig witness, exact optimum | all 512 relations on 3 points (incl. loops) x all R, C: 134,217,728 instances; exact optimum for 200k; random on 4-6 points (3.3M) | holds; the lemma is equivalent to the matching bound, which is tight in 32% of cases |
| | Prop. 3.1 (no graph with abs(U)+abs(F^2 U) < 2 abs(F U)); Seymour's statement; subset-deficit property | every oriented graph on n <= 7 vertices (2,142,288 at n = 7, via nauty) | none exists; the closest graphs miss by exactly 1 |
| | Section 3 steps C1-C6 (proper columns, diagonal kept and conflict-free, H outside F1P u F2Q, strict gain, maximiser inequality) | all compatible pairs, n <= 4; 6.3M random pairs, n = 5, 6 | all hold |
| | Lemma 1.2 step inequalities; blow-up D[T_m] bounds | 277k steps (n <= 9), 98k blow-ups | all hold |
| Saxl (24 Sep) | Theorem: g(rho_m, rho_m, lam) > 0 | m <= 8 (N = 36, 17,977 targets) | holds |
| | Cyclic Thm 3.1: the single cyclic module W_m holds every irreducible | m <= 4 (S_10 enumerated) | holds |
| | Band quotient: unique high set, factorization, signs; path contraction; segment table; Example 5.4; Q-map | m <= 9 / 5 / 7 / all / m <= 4 | all exact |
| | Two-band Prop.: every eta with <= 4 rows, by brute-force pairing; independent isotypic check | m <= 8; m <= 5 | holds (exactly the <= 4-row shapes) |
| Universal tensor squares (24 Sep) | Finite range: universal self-conjugate square exists exactly when n is not 2, 4 or 9 | n <= 24, independent characters | matches |
| | Odd-path support; one extra copy (cyclic supports) | k <= 9; abs(E) <= 9 | hold (one-extra support is larger than claimed) |
| | Candidate L; capacity rectangle; large-M bounds; 2,832 intermediate assertions; balance packing and Young inequalities; odd strip; two-path Pieri step; two-extra LR splits | n <= 30000; n <= 26; M <= 2000; exact; M <= 400 / 16; c <= 80; q <= 14; r <= 11 | all hold; failing (M,r) are exactly the 19 listed exceptions |
| Hill: crossing number of K_n (23 Sep) | interpolation and order lemmas; cycle pairing; injectivity; isotropy; bound abs(T)/2 >= Z(n) for every normalized drawing | s <= 5; straight and bent drawings of K_5, K_7, K_9; hill-climbs | all hold; climbs reach Z(n) = 1, 9, 36, never below |
| Zarankiewicz: crossing number of K_{m,n} (23 Sep) | subspace bound (sum of dim(U_i cap V_j) <= floor(m^2/4) N); rank detector; additive relation; rank sum | exhaustive Q^2 (m <= 5), Q^3 (m <= 3), climbs Q^3, Q^4; K_{3,5}, K_{4,5}, K_{5,5}, K_{3,7} | all hold; bound attained, never exceeded |

## Hand check
Every step of the first three papers was read against the TeX: Seymour (subset deficit, blow-up, pruning via generic ranks NL = AB and the kernel-support lemma, extremal contradiction, weighted and directed-cycle corollaries); Saxl (Young's rule, dominated case, sector induction, band factorization, path/segment algebra, strip reduction, induction); and for the tensor-square paper the candidate, band, balance and capacity sections. Its 64-degree modular character computation was not re-run; n <= 24 was recomputed independently. No false step, gap or quantifier slip found. The crossing-number papers prove a stronger statement than claimed: in every normalized drawing, the number of independent edge pairs with nonzero signed crossing count is at least Z(n) (resp. d_m d_n). Our searches never beat it.

## Lean (statements read; not built)
- SeymourSecondNeighborhood, Saxl, UniversalTensorSquares, CompleteCrossing and BipartiteCrossing Comparator statements are faithful. Oriented means loopless and asymmetric, with the exact N2. Kronecker is dim Hom into the tensor square of Specht modules realized as polytabloid orbit spans. The exceptions 2, 4, 9 are excluded. Drawings use injective continuous arcs, proper (chart) crossings, no triple points, and count crossing points.
- Solutions exist: SecondNeighborhood (43 files), Saxl (38), UniversalSquare (538), CompleteCrossing (46), Crossing (339, for the bipartite paper). All five were scanned with rg: no sorry, admit, native_decide, axiom, extern or unsafe, and no unusual options.

## Prior reports
openai/math has issues and discussions disabled; no public refutation of these papers found.

## Caveat
One run (saxl_cyclic S2-S4) peaked at 1.35 GB because the numpy batch was too large. cyc.py now caps the batch at about 50 MB; every other run stayed under 330 MB.
