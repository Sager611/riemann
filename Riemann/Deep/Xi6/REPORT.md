RECORD: ξ(K₆) ≥ 33 (previous proved bound 25)

# ξ(K₆) ≥ 33: a 3-colouring of the edges of K₃₃ in which every K₆ sees all three colours

Agent xi6, 2026-11-05. Status: final for this run.

## Short answer

- New lower bound **ξ(K₆) ≥ 33**, by an explicit colouring of K₃₃ with a cyclic symmetry of order 11
  (file K33_Z11x3.txt), checked over all C(33,6) = 1,107,568 six-sets by an independent brute-force checker.
- The best previously proved lower bound we found is **25** (the (k−1)² grid colouring). The claimed 31 of
  arXiv 2611.03878 is false as proved (REFUTATIONS.md). The value 31 is nevertheless true (any 31 vertices
  of our K₃₃; also several other colourings of K₃₁ found here), but **no circulant colouring of K₃₁ works**:
  the exhaustive search over all colourings of the 15 difference classes of Z₃₁ has no solution.
- Upper bound: **ξ(K₆) ≤ 103** (from R(K₆, K₆−e) ≤ 104). So 33 ≤ ξ(K₆) ≤ 103.
- Side result: **ξ(K₇) ≥ 55** (circulant on Z₅₅, file K55_circulant_K7.txt); the published bound was 31.

## 1. The problem

ξ(G) (Almási–Axenovich, arXiv 2512.08053; also 2611.03878) is the largest n such that some red-blue-green
colouring of the edges of K_n gives every copy of G all three colours. For G = K₆ this means: split the edges
of K_n into three graphs G₀, G₁, G₂ such that no union Gᵢ ∪ Gⱼ contains K₆ (equivalently, every Gᵢ has
independence number at most 5). In older terms ξ(K_k) = R^3_2(K_k) − 1, the weakened (2-chromatic) Ramsey
number of Chung–Liu (1978) and Harborth–Möller (1999), and ξ(K_k) + 1 = R(k; 3, 2), the set-colouring
Ramsey number.

## 2. The colouring of K₃₃

Vertices (x, i) with x ∈ Z₁₁, i ∈ {0, 1, 2}; in the file, vertex 11·i + x. The shift x ↦ x + 1 (on all three
blocks at once) preserves the colouring. With d = y − x (mod 11):

- inside block i, the colour of {(x,i),(y,i)} is Wᵢ[min(d, 11 − d)], listed for 1..5:
  W₀ = 0 0 1 2 0,  W₁ = 2 2 1 2 0,  W₂ = 1 0 0 2 0;
- across blocks i < j, the colour of {(x,i),(y,j)} is X_ij[d], listed for d = 0..10:
  X₀₁ = 1 0 0 2 0 1 0 0 1 2 2,  X₀₂ = 1 2 1 1 1 2 0 2 2 2 0,  X₁₂ = 2 0 2 1 0 1 1 1 0 2 0.

Colour classes: 198, 154 and 176 edges. Colour degrees per vertex: {9,10,13}, {9,11,12} or {10,10,12}.
Every K₇ ⊂ K₃₃ contains a K₆, so the same colouring gives ξ(K₇) ≥ 33 as well (superseded by §7).

## 3. Verification

verify.py (Python standard library only; shares no code with the search) reads the matrix, checks symmetry
and that entries are 0/1/2, then runs explicit nested loops over every 6-subset (no pruning) and ORs the
15 edge colours:

    python3 verify.py K33_Z11x3.txt
    n=33 k=6 subsets checked=1107568 (=C(33,6)) two-coloured=0 first_bad=None
    VALID: every K_6 uses all three colours

Controls: the K₃₁ cubic-residue colouring of 2611.03878 gives two-coloured=186 (first bad set
{0,1,2,4,8,16}, as in REFUTATIONS.md); the K₁₉ cubic-residue colouring passes for k = 5; the K₃₁ colouring
passes for k = 7 (consistent with Budden–Stender–Zhang); a K₅₅ colouring with K₇ on {0..6} forced
monochromatic fails (4549 bad 7-sets). For k ≠ 6 and n > 30 verify.py uses a complete depth-first
enumeration that counts the completions of an already three-coloured prefix at once (still every
k-subset is accounted for; the total is asserted to equal C(n,k)).

Also verified (k = 6): K26_circulant.txt (largest circulant), K28_Cayley_D14.txt (largest Cayley colouring found).

## 4. Previous bounds and literature checked

Lower bounds known before this run:
- ξ(K_k) ≥ (k−1)²: grid colouring (same row red, same column blue, rest green). For k = 6: **25**
  (stated in 2611.03878, Thm 3.1 and Cor. 3.1). Exact: ξ(K₃) = 4, ξ(K₄) = 9 (Chung–Liu 1978);
  19 ≤ ξ(K₅) ≤ 28 (Ortlieb, KIT bachelor thesis 2017, quoted in 2512.08053).
- Budden–Stender–Zhang, "Weakened Ramsey numbers and their hypergraph analogues", Integers 17 (2017) #A23
  (Zenodo 10.5281/zenodo.10636186): cubic-residue colourings give R^3_2(K₅) > 19 and R^3_2(K₇) > 31,
  nothing for K₆ (the cubic-residue graph mod 31 has the bigger clique number 6); they write that little
  is known for n ≥ 5.
- Almási–Axenovich, arXiv 2512.08053 (Dec 2025): values for K₃, K₄, K₅ only, plus ξ(K_t) ≥ 1.224^t.
- arXiv 2611.03878 (Nov 2026): ξ(K₆) ≥ 31, false as proved (REFUTATIONS.md).
- Radziszowski, Small Ramsey Numbers, DS1.18 (April 2026): explicitly excludes weakened and
  set-colouring Ramsey numbers. Set-colouring papers (2206.11371, 2212.06802, 2305.14132, 2505.20652,
  2609.35079) are asymptotic or about stars/paths/matchings.
- Search: evidence API (2 claims), arXiv API (8 queries), OpenAlex (6 queries + author list of M. Budden),
  Zenodo. Not checked (not open): M. Budden, *Weakened Ramsey Theory* (CRC Press 2026), Ch. 3; the KIT
  theses of Ortlieb (2017), Lee (2022), Almási (2023).

## 5. Upper bounds

1. ξ(k) < R(k,k): merge red and blue; a red/blue K_k or a green K_k is two-coloured.
   R(6,6) ≤ 160 (DS1.18, Table Ib) gives ξ(K₆) ≤ 159.
2. ξ(k) < R(K_k, K_k − e) (Almási, quoted in 2611.03878): merge red and blue into yellow. A yellow K_k is
   two-coloured; a green K_k − e plus its missing edge is a K_k with green and one other colour.
   R(K₆, K₆ − e) ≤ 104 (DS1.18, Table IIIa) gives **ξ(K₆) ≤ 103**.
3. Neighbourhood recursion: in a colour-c neighbourhood the two unions containing c must be K₅-free, so
   F(a,b,c) ≤ 1 + F(a,b−1,c−1) + F(a−1,b,c−1) + F(a−1,b−1,c) (F = largest n with Gᵢ ∪ Gⱼ free of the
   given cliques). With Chung–Liu's exact small values this gives only ξ(K₆) ≤ 208; not competitive.

## 6. What was searched (k = 6 unless stated)

Solvers: own exhaustive C backtracking (circulants), CaDiCaL 1.5.3/3.0 and Kissat 4.0.4 via python-sat,
own WalkSAT-style local search. Each SAT variable is an orbit of edges under the imposed symmetry group;
each constraint is the set of edge orbits met by a six-set (only inclusion-minimal sets kept). UNSAT results
are solver verdicts without proof certificates.

| structure | result |
|---|---|
| circulant on Z_n, exhaustive | solutions for every n = 20..26; none for n = 27..46 (incl. 31) |
| Cayley colourings, all groups of orders 27, 28, 29, 30, 31, 33, 34, 35, 37, 38, 39; 7 abelian + ~15 non-abelian of order 32; ~12 of order 36; ~13 of order 40; Z₇², Z₂⁶, Z₄³, Z₈², Z₃⁴ | among groups of order ≥ 27 only the dihedral group D₁₄ (**n = 28**) works; all others UNSAT (Z₅², Z₂₅, D₁₃ also work, n = 25, 26) |
| Z_m acting with r orbits (+ f fixed points) | SAT: 9×3+2 (29), 10×3 (30), 10×3+1 (31), **11×3 (33)** |
| same, UNSAT | 14×2+1, 15×2, 14×2+2, 15×2+1, 16×2 (780 s), 10×3+2, 15×2+2, 16×2+1, **11×3+1 (34)**, 17×2 (34, 1479 s), 16×2+2 (34), 10×3+4 (34), 11×3+2 (35), 17×2+1 (35), **12×3 (36, 1046 s)**, 12×3+1 (37) |
| same, undecided (time limit) | 8×4+2 (34, 1100 s), 7×5 (35, 600 s), 9×4 (36, 300 s) |
| Z_m × 3 with an extra colour-rotating symmetry (orbit cycling, a ↦ a + 1 on colours) | SAT 27 (9×3), 31 (10×3+1); UNSAT 33, 36, 39 |
| Z_m × r with the reflection x ↦ −x | UNSAT for n = 33..39 tried |
| one-vertex extensions of 33 distinct 33-vertex solutions | all UNSAT (3.7–5.5 thousand clauses each) |
| local search on 34–36-vertex symmetric instances | best: 2 violated constraints (Z₁₇ × 2, later proved UNSAT) |

The 34-vertex case Z₁₁ × 3 + 1 was closed by case analysis: the fixed vertex sends one colour to each
block; patterns (0,0,0) and (0,0,1) are UNSAT directly; for (0,1,2) each block's inner colouring must be
one of 30 (6 up to multipliers), and all 180 block-0/block-1 pairs are UNSAT. The 36-vertex Z₁₂ × 3 run
used two valid symmetry reductions (colours numbered in order of first appearance in W₀; blocks 1 and 2
shifted so that X₀₁[0], X₀₂[0] carry the smallest colour present).

So within all cyclic block structures we could decide, 33 is the maximum; n = 34–36 stays open.

## 7. Side results

- **K₇**: circulant colourings in which every K₇ sees three colours exist for every n = 31..51 and for
  n = 55 (Z₅₅, colour of difference ±d given by s[d] for d = 1..27, s = 012201221020121102111000210);
  none for n = 52, 53, 54, 56, 57, 58, 59, 60, 63; n = 61, 62, 64 undecided. verify.py confirms K₅₅
  (all C(55,7) = 202,927,725 seven-sets). So **ξ(K₇) ≥ 55**, against 31 published (Budden–Stender–Zhang
  2017). Larger k (e.g. K₈, published 43) were not searched; circulants should improve them as well.
- **K₅**: nothing beyond 19 found: circulants n = 20..30 (exhaustive, none), 23 Cayley groups of
  orders 20..27, cyclic block structures n = 20..28 (all UNSAT or undecided). Full SAT without symmetry
  reproduces ξ(K₄) = 9 (n = 9 SAT, n = 10 UNSAT in 0.2 s) but cannot even refind the K₁₉ colouring in
  300 s, so ξ(K₅) = 19 versus 20 is out of reach this way.

## 8. Files (Riemann/Deep/Xi6/)

- REPORT.md (this file); verify.py (independent checker).
- K33_Z11x3.txt (record), K28_Cayley_D14.txt, K26_circulant.txt (k = 6); K55_circulant_K7.txt (k = 7).
- search/: the search code (C: cay.c exhaustive circulants, hyp.c minimal constraint sets, ls2.c local
  search; Python: sym.py, grp.py, groups.py, batch.py, fam.py, ext.py, ext_many.py, cases.py, cube34.py,
  sb3.py, sbr.py, fullsat.py; run.sh, runf.sh). Needs clang and python-sat (pip install python-sat);
  e.g. (in search/, after clang -O2 -o hyp hyp.c) python3 fam.py P:11:3:0 sat 600 refinds a
  33-vertex colouring in about a minute (CaDiCaL); results go to search/sol/.

Colouring format: comment lines start with '#', then n, then n rows of n characters (colour 0/1/2 of
edge ij, '.' on the diagonal).

