# Moser's worm problem: smaller convex universal cover (record hunt, 2026-10-22/23)

## Result
alpha <= 1785486893883814994120173447337/7277367935880000000000000000000 = 0.2453478935812...
(Keller, arXiv 2609.21968: 0.246332237201...). Cover = Keller's quadrilateral K0 scaled by 499/500
(t = (-631/58, -1, 7/16, 37/56), h = (0, 0, 65260717/125000000, 222653301/500000000)); polygon.json.
Alpha here = Keller's (congruent copies, reflections allowed); not Tao et al.'s constant 13a.

## Keller's method (paper Sec. 4), as re-implemented here
* K = {x : n_i . x <= h_i}, n_i = n(t_i) rational unit normals, h_0 = h_1 = 0. Translation criterion:
  Q fits iff Psi_j(Q) = a_j0 h_Q(n_0) + a_j1 h_Q(n_1) + h_Q(n_j) <= h_j for all j >= 2.
* Uncovered arc -> convex polygon Q with shortest visiting path <= 1; normal form with bottom edge
  [P-, P+], top vertex P0. Support contacts in kR + kL = 60 rational directions; for each listed
  orientation one escape inequality holds ("vrot" orientations fill gaps).
* Leaf: multipliers on support, bottom-edge and inherited escape rows give f with sum f = 0,
  |sum_I f| <= B on the prefix family F(R) and B = sum mu h_j > 0 => visiting length > 1.
  Tree nodes branch on escape alternatives or split the contact range R.

## Correspondence with Keller §4
Sources: Keller, arXiv 2609.21968v1, TeX file `sections/upper.tex` (= Section 4; "u.N" = its line N) and
`verify_upper.py` ("V.N" = its line N, version of 2026-10-23 10:05 UTC). Paper numbering: Theorem 4.1
(upper bound), Lemma 4.2 (finite witness), Lemma 4.3 (normal form), Lemma 4.4 (translation criterion);
equations (4.1) K, (4.2) escape, (4.3) boundary order, (4.4) support, (4.5) length, (4.6) contact pair,
(4.7) leaf checks. Keller's proof of Theorem 4.1 (u.442-457) is a chain; every finite input of it is checked:

| Keller's step (quoted or summarised, location) | Check in verify_upper.py |
|---|---|
| (4.1) side data (u.17); "uniquely determined positive rational coefficients a_j0, a_j1" with a_j0 n_0 + a_j1 n_1 + n_j = 0 (u.187-191) | check_shape (V.53-107): rational unit normals (V.60), h_0 = h_1 = 0 < h_j (V.62-64), each turn left (V.68), total turning exactly 2π (V.69-75, exact), a_j0, a_j1 > 0 and the identity (V.78-85), every side active (V.89-100), exact area (V.101-103) |
| Lemma 4.2 (u.67-75): "If K does not cover a unit arc, it does not cover some convex polygon Q with m ≥ 3 vertices ... ℓ_Q ≤ 1. Every prefix {Q_1,...,Q_j} is a boundary interval." Needs a unit segment in K to exclude collinear Q (u.51-52, u.95-96) | K has diameter ≥ 1 (V.104-106); the rest of the lemma uses no data |
| Lemma 4.3 (u.121-127): "After a rigid motion, Q has a bottom edge [P−,P+] on the x-axis and a topmost vertex P0 such that x(P−) < x(P+), τ(P−) < τ(P0) < τ(P+)"; the motion may be a reflection (u.152-153) | no data; the reflection is why the bound is for congruent copies with reflections |
| (4.3) (u.231-241): rational unit directions u^R_i in the open right half-plane "strictly ordered from downward to upward", u^L_j in the open left half-plane "strictly ordered from upward to downward"; weak cyclic order P−, P+, R_1..R_kR, P0, L_1..L_kL | directions M n_i of the non-"vrot" orientations (V.152-159), split by sign of x, sorted, strictness checked (V.160-165), kR, kL match the header (V.166-168); labels 0 = P−, 1 = P+ (direction (0,−1)), 2 = P0 ((0,1)), then R_1..R_kR, then L_1..L_kL (V.169-171): the order (4.3) |
| (4.4) (u.251-255): u_λ·(P_λ − P_η) ≥ 0 (λ, η ∈ Λ) and x(P+) − x(P−) ≥ 0 | support row (λ, η), λ ≠ η: f_λ += μ u_λ, f_η −= μ u_λ (V.286-292); bottom row: f_{P+} += μ e_x, f_{P−} −= μ e_x (V.293-295); right side 0 |
| Lemma 4.4 (u.198-208): for every M ∈ O(2) an uncovered Q satisfies a_j0 h_Q(Mn_0) + a_j1 h_Q(Mn_1) + h_Q(Mn_j) > h_j for j = 2 or j = 3 | every orientation exactly orthogonal (V.110-131: "flush", "rot", "vrot", optional mirror); escape row (M, j) with weights (a_j0, a_j1, 1) and right side h_j (V.189-209) |
| Linear escape inequalities (u.257-274): u = Σ κ_λ u_λ with κ ≥ 0; "A direction already chosen uses a one-term expansion, with P+ used for the downward direction" | one-term rows must hit a contact direction (V.198-201), (0,−1) ↦ P+ (V.171); "vrot" directions are expanded into their two angular neighbours with exact nonnegative coefficients (V.175-184); each row sums to 0 (V.206-208) |
| Covering all cases (u.398-406): root without assumptions on [0,kR]×[0,kL]; a node branches into both escape alternatives of one M or covers its rectangle by child rectangles; "Each leaf uses only the inherited escape inequalities" | tree walk (V.217-248): one child per alternative, child r = side 2 + r (V.226-232); child cells inside the parent and covering all its integer pairs (V.233-244); root (0, kR, 0, kL) (V.248); a leaf may use row (M, r) only if assumed on its branch (V.296-300) |
| (4.7), certifying a leaf (u.408-433): μ ≥ 0, B = h_2 Σ_{E_2} μ + h_3 Σ_{E_3} μ, B ≥ 1, Σ f = 0, \|Σ_I f\|² ≤ 1 for I ∈ F(R) | positive integer multipliers (V.282-284), B (V.304), B > 0 (V.307-308), Σ f = 0 (V.310-311), \|Σ_I f\|² ≤ B² (V.313-336): (4.7) after dividing every μ by B |
| (4.5) (u.283-305): \|w_r\| ≤ 1 for all r gives Σ f·P ≤ ℓ_Q | no data |
| F(R), "Which prefixes can occur" (u.315-390) | the four loops below |

| Stage (u) | Keller's set and range | verify_upper.py |
|---|---|---|
| r < τ(P−) (u.372-373) | Λ_r = {L_a..L_b}, q_lo+1 ≤ a ≤ b ≤ kL | V.323-325: a = q_lo+1..kL, b = a..kL |
| τ(P−) ≤ r < τ(P0) (u.374-375) | Λ_r = {P−, L_j..L_kL}, 1 ≤ j ≤ q_hi+1 | V.326-328: j = 1..q_hi+1 |
| τ(P0) ≤ r < τ(P+) (u.376-377) | Λ_r^c = {P+, R_1..R_s}, p_lo ≤ s ≤ kR | V.329-330: s = p_lo..kR |
| τ(P+) ≤ r (u.378-379) | Λ_r^c = {R_{a+1}..R_b}, 0 ≤ a ≤ b ≤ p_hi | V.331-333: a = 0..p_hi, b = a..p_hi |

Empty sets (allowed, u.384) have sum 0 and pass because B > 0, so V skips empty stage-1 sets. Checking
Λ_r^c instead of Λ_r is equivalent because Σ f = 0 (u.310-313). k = 4 here; the k-gon path is unused.

Why these checks exhaust all arcs (re-derivation of u.315-390 and u.442-457):
1. If a unit arc is not covered, Lemma 4.2 gives an uncovered convex polygon Q, ℓ_Q ≤ 1, whose visited
   prefixes are boundary intervals; Lemma 4.3 puts it in normal form. Take its contacts, (p, q) by (4.6);
   (4.4) holds for Q.
2. Tree path: at an orientation node Lemma 4.4 gives a true alternative; at a split node (p, q) lies in a
   cell. So Q reaches a leaf whose assumed escape rows hold and whose rectangle R contains (p, q).
3. The stages are exhaustive: by Lemma 4.3, τ(P−) < τ(P0) < τ(P+), so each r = 0..m lies in exactly one
   of r < τ(P−), [τ(P−), τ(P0)), [τ(P0), τ(P+)), r ≥ τ(P+).
4. Each prefix Λ_r and its complement are arcs of the boundary cycle (4.3), and [P−, P+] is an edge.
   (a) τ(L_j) ≤ τ(P0) ≤ τ(R_i): if τ(R_i) < τ(P0), the prefix at max(τ(R_i), τ(P−)) contains P− and R_i
   but not P+ or P0, which no arc does; the prefix at τ(P0) contains P−, P0 but not P+, so it contains
   P0, L_1..L_kL, P−. (b) The prefix at τ(P−) avoids P+ and P0, so its L-labels are a final block
   L_{q+1}..L_kL; the prefix at τ(P+) − 1 contains P0, P− but not P+, so its R-labels are a final block
   R_{p+1}..R_kR. This is (4.6). Hence: stage 1, the arc avoids P−, P+, P0 and holds no R-label, so
   Λ_r = {L_a..L_b} with a ≥ q+1 ≥ q_lo+1; stage 2, the arc ends at P− and stops before P0, so
   Λ_r = {P−, L_j..L_kL} ⊇ {P−, L_{q+1}..L_kL}, j ≤ q+1 ≤ q_hi+1; stage 3, the complement starts at P+,
   stops before P0, holds no L, so it is {P+, R_1..R_s} ⊇ {P+, R_1..R_p}, s ≥ p ≥ p_lo; stage 4, the
   complement avoids P−, P+, P0 and every L and lies in R_1..R_p, so it is {R_{a+1}..R_b}, b ≤ p ≤ p_hi.
   Empty and full prefixes give w_r = 0.
5. So \|w_r\| = \|Σ_{Λ_r} f\| = \|Σ_{Λ_r^c} f\| ≤ B for every r; with every μ divided by B, (4.5) gives
   Σ f·P ≤ ℓ_Q ≤ 1, while the leaf's rows give Σ f·P > 1 (B > 0 forces a positive multiplier on a strict
   escape row). Contradiction, so K covers every unit arc.
The only case splits in Keller's §4 are the escape alternatives (Lemma 4.4), the (p, q) rectangle and the
four stages; V checks the first two for completeness, and the third is exhaustive by step 3.

## Files
* verify_upper.py: independent exact checker (stdlib, Fractions, no floating point) written from the paper's
  statements: shape, orientations, expansions, escape rows, tree coverage, leaf conditions in
  homogeneous form. Accepts Keller's certificate (12 663 leaves) and cert/ (29 137 leaves).
* cert/: header.json (t, h, 1 493 orientation specs, kR = kL = 30), tree.json.gz,
  multipliers-NNN.jsonl.xz (lines [leaf, den, key_deltas, integer weights]), manifest.json (sha256).
* search/: searchk.py (tree search driver; --jobs 1 runs in-process; --init-tree reuses an old tree;
  --resume from checkpoint), modelk.py (Clarabel proposals, adaptive orientations, exact rounding,
  sparse re-solves, row dropping), exactk.py, crosscheck_exact.py (leaves vs Keller's own exact.py),
  tamper_tests.py, compact_cert.py, certstats.py, polyinfo.py, adversary.py/adv2.py/adv3.py (hardest
  arcs), shapeopt*.py/localshape.py (shape search), K0_s998.json (config of the record run),
  requirements.txt (numpy 2.0.2, scipy 1.13.1, clarabel 0.11.1; the checkers need only Python).

## History
1. Keller's certificate re-verified (ACCEPT, 208 s). Its leaf programs have little slack
   (min optimum 1.0000079): the tree must be regrown for a smaller K.
2. Scale 0.9995: 14 986 leaves (from Keller's tree, adaptive). Scale 0.998: 29 137 leaves, 27 more
   orientations, 43 300 expansions in 23 min with 3 workers (from the 0.9995 tree). Record.
3. Re-check 2026-10-23: verify_upper.py ACCEPT (~12 min, one low-power core, 56 MB); final re-run
   09:14 UTC on the shipped cert/: ACCEPT, same area; polygon.json checked exactly against cert/header.json;
   crosscheck_exact.py OK (29 137 leaves pass Keller's exact.Geometry.valid; ~6 min).
   Review fixes (10:05 UTC): the floating-point atan2 winding test in check_shape replaced by an
   exact count of angle descents (rejects star-pentagon normals, accepts Keller's K and ours);
   correspondence section added; final verify_upper.py run on cert/: see REPORT.md.
4. Hardest arcs for K0: splayed 3-segment "U" (legs ~0.275, base ~0.45), rho = 0.99615: this shape
   cannot certify below scale ~0.9962 (area ~0.2444). Other quadrilaterals and pentagons:
   estimated limits 0.2445-0.2452, no gain.
5. Scale 0.997 (area 0.24486): grown from the 0.998 tree; after ~66 000 expansions 29 700 leaves done,
   35 900 cells open and flat (single-core rate ~300 expansions/min under the CPU limits). It needs
   >= 65 000 leaves = >= 30 MB in this format (budget 20 MB) and est. 4-11 more CPU-hours: stopped.
6. Certificate size: 461 bytes per leaf after xz; ~572 rows per leaf with ~2.4-digit weights (dense
   interior-point solutions). A leaf only needs about 2 rows per contact point (Caratheodory), so
   sparse basic re-solves should shrink certificates ~4x and make scales 0.997-0.9965 affordable.

## Novelty check (2026-10-23)
arXiv API ("worm problem"; worm + convex), math.MG and math.CO October 2026 listings, Zenodo: no
upper bound newer than 2609.21968. New on Zenodo today: lower bound 0.2414 for the rotation-only
version (Tao et al. constant 13a, whose best listed upper bound is still pi/12).

## Lower bound
Not attempted: it needs a new certified minimisation of hull areas over arc placements (the
rotation-only proof above used 434 million leaves, 9.9 GB). The upper bound was the easier target.
