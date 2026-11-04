# Gale's conjecture fails in R⁷ and R⁸, and in every dimension n ≥ 5: a simplex-product family with a short proof

Orca record hunt, 2026-11-02 to 2026-11-04. Gale (1953) conjectured that the ball is the hardest set to cut: β_{n+1}(V) ≤ β_{n+1}(Bⁿ)·D(V) for every bounded V ⊂ Rⁿ. Voronov (arXiv 2610.28828, v1, 20 Oct 2026; github.com/vsvor/Gale-conjecture) found computer counterexamples for n = 3–6 (our records for n = 4, 5, 6: ../Gale45, ../Gale6); the paper has no n = 7 example (its 7-row is commented out) and asks (Discussion, Problem): *prove that a counterexample exists for every n ≥ 3, describe it geometrically, and prove it without exhaustive enumeration of graph colourings.*

**Results.**
1. **R⁷:** an 84-point set with rational coordinates (points7.txt), **β₈/D = √(8900856812318928/10000000213244525) = 0.9434435130… > σ₇ = 0.9419651451** (+1.478·10⁻³). The clean 69-point member of the family below gives 0.9433983553 with exact parameters (points7_sym69.txt, on the 10⁻⁸ grid 0.9433983441); the 84-point set is that one plus twins and a free polish.
2. **R⁸:** the 85-point member of the family (points8.txt), **β₉/D = √(9023994457279939/10000000157507010) = 0.9499470677… > σ₈ = 0.9486832981** (+1.264·10⁻³).
3. **Every n ≥ 5:** an explicit five-orbit set Vₙ in Rⁿ = R^⌊n/2⌋ ⊕ R^⌈n/2⌉ (two regular simplices) with β_{n+1}(Vₙ)/D(Vₙ) > σₙ. The colouring part is a short lemma (§3), proved by hand; the metric part is an exact check of 36 rational numbers per dimension, done for every n from 5 to 1000 (family.py, 1.6 s). For larger n Borsuk's conjecture itself fails (counterexamples are known for every n ≥ 64), so with the computer-certified sets for n = 3, 4 (../Gale, ../Gale45, both with DRAT proofs), **Gale's conjecture fails in every dimension n ≥ 3.** This answers the paper's problem for n ≥ 5; n = 3, 4 still need a computer (the family gives only 0.8752 < σ₃ and 0.9106 < σ₄).

## 1. Definitions
For finite V ⊂ Rⁿ with diameter D, β_k(V) is the least possible largest part diameter over partitions of V into k parts. G(≥s) joins the pairs at squared distance ≥ s. If G(≥s) has no proper (n+1)-colouring, every partition into n+1 parts has a part with two points at distance ≥ √s, so β_{n+1}(V) ≥ √s. If dc² is the largest squared distance with G(≥dc²) not (n+1)-colourable, a proper (n+1)-colouring of G(>dc²) is a partition with parts of diameter ≤ dc, so β_{n+1}(V) = dc. The regular-simplex partition of the ball gives β_{n+1}(Bⁿ) ≤ σₙ with σₙ² = (n+1)/(n+2) (n even) and σₙ² = 1/2 + ½√((n−1)/(n+3)) (n odd): σ₇ = 0.9419651451, σ₈ = 0.9486832981. Any V with β_{n+1}(V)/D(V) > σₙ is a counterexample.

## 2. The family
Rⁿ = R^{m₁} ⊕ R^{m₂}, m₁ = ⌊n/2⌋, m₂ = ⌈n/2⌉; u₀…u_{m₁} and v₀…v_{m₂} are unit vectors of regular simplices (uᵢ·u_k = −1/m₁, v_j·v_l = −1/m₂ for different indices). With 8 parameters p = (r_A, r_B, a, b, c, d, e, f):

  Pᵢ = r_A uᵢ (poles of the first simplex), Q_j = r_B v_j (poles of the second), A_ij = a uᵢ − b v_j (near Pᵢ, pushed away from Q_j), B_ij = −c uᵢ + d v_j (near Q_j, pushed away from Pᵢ), C_ij = e uᵢ + f v_j (between Pᵢ and Q_j),

i = 0…m₁, j = 0…m₂: (m₁+1) + (m₂+1) + 3(m₁+1)(m₂+1) points (n = 7: 69, n = 8: 85). The five points Pᵢ, A_ij, C_ij, Q_j, B_ij of one cell lie in the plane of uᵢ and v_j, roughly on one circular arc. A squared distance depends only on the two orbits and on which indices agree, so all distances are 36 explicit quadratic forms in p ("classes"); for rational p they are exact rationals. At the optimum the poles are at distance D within each simplex (r_A, r_B = circumradii of unit-edge simplices), the binding short edges are P_i B_ij, Q_j A_ij, A_ij C_kl and B_ij C_kl (i ≠ k, j ≠ l), and the pairs at distance D are P–P, Q–Q, A_ij B_ij, A_ij C_kj (i ≠ k) and, depending on n, A_ij A_kl (i ≠ k, j ≠ l) or B_ij C_il (j ≠ l).

  | n | split | points | β_{n+1}/D (this family, exact check) | σₙ | excess |
  |---|---|---|---|---|---|
  | 5 | 2+3 | 43 | 0.9246685234 | 0.9238795325 | +7.89e-04 |
  | 6 | 3+3 | 56 | 0.9368472505 | 0.9354143466 | +1.43e-03 |
  | 7 | 3+4 | 69 | 0.9433983553 | 0.9419651451 | +1.43e-03 |
  | 8 | 4+4 | 85 | 0.9499470859 | 0.9486832980 | +1.26e-03 |
  | 9 | 4+5 | 101 | 0.9542066051 | 0.9530206138 | +1.19e-03 |
  | 10 | 5+5 | 120 | 0.9585019831 | 0.9574271077 | +1.07e-03 |
  | 11 | 5+6 | 139 | 0.9614803843 | 0.9605087856 | +9.72e-04 |
  | 12 | 6+6 | 161 | 0.9645018124 | 0.9636241116 | +8.78e-04 |
  | 16 | 8+8 | 261 | 0.9724289336 | 0.9718253158 | +6.04e-04 |
  | 20 | 10+10 | 385 | 0.9774452783 | 0.9770084209 | +4.37e-04 |
  | 32 | 16+16 | 901 | 0.9853861827 | 0.9851843661 | +2.02e-04 |
  | 63 | 31+32 | 3233 | 0.9923378425 | 0.9922760123 | +6.18e-05 |
  | 64 | 32+32 | 3333 | 0.9924551619 | 0.9923953268 | +5.98e-05 |
  | 100 | 50+50 | 7905 | 0.9951120711 | 0.9950859653 | +2.61e-05 |
  | 1000 | 500+500 | 754005 | 0.9995011632 | 0.9995008734 | +2.90e-07 |

The excess decays like ≈ 0.29/n² and stays positive at every n checked (min 2.9·10⁻⁷ at n = 1000). Parameters: family_params.txt (12 decimals, read as exact rationals).

## 3. The colouring lemma (no computer needed)
Let p, q ≥ 1, rows i ∈ I (|I| = p), columns j ∈ J (|J| = q). H(p, q) is the graph on Pᵢ, Q_j, A_ij, B_ij, C_ij with edges
- PᵢP_k, A_ijA_kl, P_kA_ij, P_kC_ij, A_ijC_kl whenever i ≠ k;
- Q_jQ_l, B_ijB_kl, Q_lB_ij, Q_lC_ij, B_ijC_kl whenever j ≠ l;
- C_ijC_kl whenever i ≠ k and j ≠ l;
- PᵢB_ij and Q_jA_ij for all i, j.

**Lemma.** H(p, q) has no proper colouring with p + q − 1 colours.

*Proof.* Suppose it has one. Let Xᵢ = {Pᵢ} ∪ {A_ij : j} (row blocks), Y_j = {Q_j} ∪ {B_ij : i} (column blocks).
(a) Two different row blocks are completely joined, and so are two different column blocks; so a colour class meets at most one row block and at most one column block.
(b) C_ij is adjacent to all of X_k (k ≠ i) and all of Y_l (l ≠ j); so the C-vertices of a class meeting Xᵢ lie in row i, those of a class meeting Y_j lie in column j.
(c) The C-vertices of a class lie in one row or one column (two C's in different rows and columns are adjacent; if a class held C_ij, C_il with j ≠ l and some C_km with k ≠ i, then m = j and m = l).
Call a class X(i) if it meets Xᵢ but no column block, Y(j) if it meets Y_j but no row block, XY(i,j) if it meets Xᵢ and Y_j (its only possible C-vertex is then C_ij), L if it meets no block (its C-vertices lie in a line). Let R₀ be the rows without an X(i) class, C₀ the columns without a Y(j) class, r₀ = |R₀|, c₀ = |C₀|, and z, ℓ the numbers of XY and L classes. There are at least p − r₀ X-classes and q − c₀ Y-classes, so
  z + ℓ ≤ r₀ + c₀ − 1.  (1)
For i ∈ R₀ every vertex of Xᵢ lies in an XY(i,·) class, so each row of R₀ and likewise each column of C₀ contains an XY class; hence z ≥ max(r₀, c₀), and (1) forces r₀, c₀ ≥ 1. For i ∈ R₀, j ∈ C₀ the vertex C_ij lies in the class XY(i,j) or in an L class, by (b). Let ℓ₁ L-classes lie in rows and ℓ₂ in columns. If ℓ₁ ≥ r₀ then z + ℓ ≥ c₀ + r₀, against (1); so ℓ₁ < r₀ and likewise ℓ₂ < c₀. At least (r₀ − ℓ₁)(c₀ − ℓ₂) cells of R₀ × C₀ meet no L class, each needs its own XY class, and
  z + ℓ ≥ (r₀ − ℓ₁)(c₀ − ℓ₂) + ℓ₁ + ℓ₂ ≥ r₀ + c₀ − 1  (xy ≥ x + y − 1 for x, y ≥ 1).
So everything is tight: all XY classes sit on the cells of R₀ × C₀ that no L class meets, one per cell. A row of R₀ met by an L class would have no XY class, so ℓ₁ = 0, likewise ℓ₂ = 0; then z = r₀c₀ = r₀ + c₀ − 1, i.e. (r₀ − 1)(c₀ − 1) = 0, with exactly one XY class S_ij for each cell of R₀ × C₀ and no other.
If r₀ = 1, R₀ = {i*}: for j ∈ C₀ the only class meeting Y_j is S_{i*j}, so B_{i*j} ∈ S_{i*j}; but P_{i*} lies in some S_{i*j} as well, and P_{i*}B_{i*j} is an edge. If c₀ = 1, C₀ = {j*}: symmetrically Q_{j*} and A_{ij*} share a class S_{ij*}, and Q_{j*}A_{ij*} is an edge. Contradiction. ∎

**Theorem.** Take p = m₁ + 1 rows and q = m₂ + 1 columns, so p + q − 1 = n + 1. If all pairs of Vₙ in the 17 classes of H (AA, BB, AC, BC split by whether the other index agrees) have squared distance ≥ t, then H ⊆ G(≥t) has no proper (n+1)-colouring and β_{n+1}(Vₙ)/D(Vₙ) ≥ √(t/D²). family.py checks t/D² > σₙ² exactly for every n = 5…1000. The class A_ijB_ij is not used (it sits at distance D anyway).
Checks of the lemma itself: family.py --lemma (exhaustive search, standard library) for (p, q) up to (4, 5); CaDiCaL for the 17-class graph up to (5, 7) (n = 10), and for the full 18-class graph up to (6, 7) (n = 11).

## 4. Certificates
`python3 verify.py` (standard library, exact integers; about 1 min) checks points7.txt and points8.txt; `python3 verify.py points7_sym69.txt` the 69-point set. For each file: all squared distances as integers (common denominator 10⁸), D², binary search for dc² over the distinct squared distances, G(≥dc²) not (n+1)-colourable by two independent exhaustive searches (forward checking, plain DSATUR, different cliques fixed), an explicit proper colouring of G(>dc²) checked edge by edge (so β = dc exactly), and the exact comparison with σₙ.

  | file | points | D²·10¹⁶ | dc²·10¹⁶ | edges of G(≥dc²) | search nodes (FC / DSATUR) | β/D |
  |---|---|---|---|---|---|---|
  | points7.txt | 84 | 10000000213244525 | 8900856812318928 | 1765 | 1327 / 23097 | 0.94344351302 |
  | points7_sym69.txt | 69 | 10000000157507010 | 8900004496494898 | 1406 | 1792 / 28190 | 0.94339834409 |
  | points8.txt | 85 | 10000000157507010 | 9023994457279939 | 2195 | 169807 / 801496 | 0.94994706774 |

Independent DRAT check: `python3 verify.py points7.txt --cnf g7.cnf --sym && drat-trim g7.cnf cert7.drat` (CNF "G(≥dc²) is 8-colourable": 672 variables, 14210 clauses, unit clauses fix a K6's colours; CaDiCaL 1.5.3 via PySAT proves UNSAT; cert7.drat, 188 KB, 11721 lines: "s VERIFIED"). The same for points8.txt (cert8.drat, 1.2 MB: "s VERIFIED"). For every n, the family needs no SAT at all: family.py checks the 36 class values exactly (all 996 rows in 1.6 s) and the lemma does the rest.

## 5. How it was found
- The R⁶ structure (../Gale6/NOTE.md §3: two orthogonal tetrahedra, 4 points per pole patch, 16 pair points) was made exactly symmetric in R³ ⊕ R⁴ and scaled; a parametric core-SLP (8 parameters; SAT critical value, minimal non-8-colourable core, LP step on the parameters) went 0.9226 → 0.9433984 in 7 steps, 1 minute. The same seed in R⁴ ⊕ R⁴ gave 0.9499253 (R⁸), in R⁴ ⊕ R⁵ 0.9542066 (R⁹).
- Reading off which distance classes are edges at the optimum gave the graph H; it is the same for n = 7, 8, 9, which suggested the lemma. With the lemma the optimisation is a 9-variable nonlinear program over 36 class values (SLSQP), so any n is cheap. The asymmetric optimum in R⁴ ⊕ R⁴ (0.9499471) beats the A↔B-symmetric one (0.9499253).
- Other splits are worse: R² ⊕ R⁵ 0.93683, R¹ ⊕ R⁶ 0.91560 (n = 7); unbalanced splits lose for every n ≤ 12.
- Dead ends for pushing R⁷ above the family value: extra symmetric orbits (never enter the core), free core-SLP polish and kicks (the 69-point optimum is rigid), symmetric orbit splits (separating two binding classes makes the graph colourable; the others give zero LP gain), adding single points to kill the 8-colourings of G(>dc²) (best new point reaches only 0.862 against 0.890 needed). What worked a little: random twins (copies of 4–16 points shifted by 0.003–0.02·D) plus a kick and a free polish (1 of 6 trials left the symmetric basin: 0.9434120, after a longer polish 0.9434306), then SAT-screened vertex splitting as in ../Gale6 (6 rounds: 0.9434306 → 0.9434435, 84 points; the last two rounds gained nothing). About 1.5 CPU-hours.

## 6. Novelty check (2026-11-04)
arXiv 2610.28828 has v1 only; its repository is unchanged since bd6dc959 (2026-10-22), with a simplex-product seed generator for dimensions 5–8 (extrapolated from a fitted R⁵ shape) but no 7- or 8-dimensional set and no proof. arXiv math.MG, math.CO, math.GM, math.OC listings for October and November 2026: the only Gale paper is Voronov's. The paper (§1) says counterexamples to Gale were known only through Borsuk counterexamples, none below dimension 9; the 9-dimensional one is OpenAI's compact set of rank-one projectors (github.com/openai/math, item 156, "Borsuk's conjecture fails in dimension nine"); openai/math lists nothing on Gale's conjecture. Borsuk counterexamples are known for every n ≥ 64 (Bondarenko 2013, Jenrich–Brouwer 2014; arXiv 2305.06283 abstract). So the sets here are the first in R⁷, R⁸ and R¹⁰…R⁶³, and the first uniform construction with a conceptual proof.

## 7. Caveats
- σₙ is the simplex-partition value; the true ball value β_{n+1}(Bⁿ) is known only for n = 2, 3, but it is ≤ σₙ, so exceeding σₙ suffices.
- The excesses are small (1.4·10⁻³ at n = 7, ≈ 0.29/n² for large n). The exact check covers 5 ≤ n ≤ 1000; beyond, an asymptotic expansion is not done (Borsuk covers n ≥ 64 anyway).
- Lean: nothing formalised. The lemma's proof is short and finite; the metric part is rational arithmetic, so a formalisation looks feasible.


## Files
- points7.txt (84 points in R⁷), cert7.drat; points7_sym69.txt (69 points, the family member for n = 7); points8.txt (85 points in R⁸), cert8.drat.
- verify.py: exact certificate check for finite sets (default: points7.txt and points8.txt), writes the CNF for drat-trim.
- family.py, family_params.txt: the family for every n = 5…1000, exact check; `--lemma` re-checks the colouring lemma by exhaustive search; `--n N` prints all 36 class values for one N.
- family_opt.py (numpy and scipy): re-derives one row of family_params.txt (SLSQP over the 36 class values, multistart; reproduces the rows for n = 7, 8 digit for digit).
