NO RECORD: every Leech-based component in 25–48D (and the P48p class for 49–60D) was tried; all are at proven ceilings, already lattice-coset shells, or locally saturated; the one new lever (r = 4 "tesseract" caps) did not pack.

# Kissing2 — extending the Leech-sublattice equator trick (kiss2, 5 Nov 2026)

## Short answer
- No new record. The 38D equator swap (τ₃₈ ≥ 624 020, Riemann/Deep/Kissing) is the only place where a large
  search-built set of Leech-shell lines could be replaced by a sublattice shell.
- In every other Leech-based record (25–31D), the equator is the Leech minimal shell itself. That shell cannot grow.
  The other layers either sit at proven ceilings or are built from cosets of L or 2L already.
- New structure found (not a record): the current 26–31D records (PackingStar plus Kravatskiy's refills) are built
  from **squares** of Leech minimal vectors. Each deleted equator vector lies in 4 caps, and **8 is the proven
  maximum**. A layout with 8 everywhere (tesseracts) would earn 1.5× the gain at the same per-line density,
  but in every test tesseracts packed far worse than squares.

Sources pinned: KissingNumbers main = 1f9e23d (4 Nov 21:04 UTC); Cohn's table, read 4 Nov (28D 221 988,
31D 251 188, 43D 2 919 782); PackingStar 477df34 (28 Oct release). The brief said dim28-joint-leech-layers holds the
28D record of 221 988. That package has 205 232; the 221 988 is dim27-31-joint-refills.

## 1. Survey of Leech-based packages (components with a shell or angle condition)
| package (dims) | record? | Leech-shell component and constraint | replacement tried | outcome / bound |
|---|---|---|---|---|
| dim38-joint-latitudes (38) | yes | equator: norm-48 lines, pairwise abs(q·q') ≤ 24 | min-48 sublattice M of index 3⁶ | done earlier: 624 020 (Kissing/) |
| dim38-latitude-extension (38) | superseded | the same 13 710-point equator | the same swap | 605 346 + 12 498 = 617 844 < 624 020 |
| dim38-leech-large-codimension (38) | superseded | 644 classes (lines with abs(cos) ≤ 1/4) partition all 98 280 lines; axis = τ₁₄ | — | at its own ceiling 3τ₂₄ + τ₁₄ = 591 612 |
| dim28-joint-leech-layers (28) | superseded (205 232) | equator = L_min minus owners; 3 caps per owner (class condition); 256 auxiliary | class via sublattice cosets (§5) | cosets cannot beat 248 lines; family ≪ 221 988 |
| norm8/frame/sqrt-layer packages (28–31) | superseded | frame-half heads (norm-8 frames) | — | superseded by the PackingStar family (§2) |
| dim25-collective-deformation (25) | yes | Leech shell + 971 deformed circles + 102 moved points | none applies (k = 1: no orthogonal tail lines, so no squares) | — |
| PackingStar family: dim26-30-asymmetric-exchanges, dim27-31-joint-refills, dim27-31-algebraic-completion (26–31) | yes | equator L_min∖D; caps on Leech edge midpoints; mid layer on frames; poles | free additions, tesseracts, completions (§2–4) | no gain; ceilings in §3 |
Non-Leech record packages (32–35 ERS/codes, 36–37 Coxeter–Todd, 39 codes, 40–46 P48n sections, 49–60 P48p caps): §6.

## 2. The PackingStar family decoded (26–31D; scripts/decode_ps.py, exact)
Scale: squared norm 32, product ≤ 16, body = first 24 coordinates, tail = the last k = d−24.
Layers in the current records:

| dim | equator (32,0) | caps (24,8) | mid (16,16) | poles (0,32) | deleted D | caps/abs(D) | squares | square groups with 32 squares |
|---|---|---|---|---|---|---|---|---|
| 26 (201 136) | 192 154 | 8 786 | 192 | 4 | 4 406 | 1.994 | 1 088 | 34 of 34 |
| 27 (209 825) | 183 543 | 25 982 | 288 | 12 | 13 017 | 1.996 | 3 233 | 98 of 120 |
| 28 (221 988) | 172 444 | 48 368 | 1 152 | 24 | 24 116 | 2.006 | 5 958 | 179 of 205 |
| 29 (231 697) | 163 504 | 66 228 | 1 925 | 40 | 33 056 | 2.004 | 8 214 | 240 of 264 |
| 30 (247 972) | 150 240 | 92 788 | 4 884 | 60 | 46 320 | 2.003 | 11 480 | 348 of 362 |
| 31 (251 196) | 148 688 | 95 832 | 6 592 | 84 | 47 872 | 2.002 | 11 916 | 366 of 376 |

Facts, all checked exactly in every record:
1. The equator is a subset of L_min, with the Golay code read off the configuration itself.
2. Every cap body b is the midpoint (a+a')/2 of a Leech **edge** (a·a' = 16). The coset b+L has minimum 8, reached
   only at ±(a−a')/2. Hence u·b ≤ 16 for every minimal u except the two endpoints, and **exactly the endpoints are
   deleted**. Endpoints are found fast from the class of 2b mod 2L.
3. Cap tails are 2·D_k roots (norm 8). Constraints: b·b' ≤ 8 on the same tail line, ≤ 12 for lines at 60°, ≤ 16 for
   orthogonal lines. The last needs a shared endpoint to fail.
4. Gain identity: cap gain = Σ over a in D of (r_a/2 − 1), where r_a is the number of caps containing a.
   Every record has r_a = 4 almost everywhere. The edge graph is a disjoint union of **squares**
   a, a−u, a−u−v, a−v (u ⊥ v). Each square has 4 deleted vectors and 8 caps, two orthogonal tail lines, gain +4.
5. Squares come in **groups** sharing a direction pair (u,v). For orthogonal minimal u, v there are exactly
   44 = 2·22 admissible centres ±f_i/2, an orthogonal frame of (u,v)⊥ (checked on 40 random pairs). The records use
   32 of the 44 in almost every group (34 of 34 in 26D, 366 of 376 in 31D).
6. Mid layer = frame halves. On one norm-16 tail, b·b' ≤ 0, so at most 48 bodies (an obtuse set in R²⁴).
   28D has 24 tails × 48 = 1 152: **saturated**. Poles ≤ τ_k: saturated in 27, 28, 29D (12, 24, 40); 26D has 4 of τ₂ = 6. In 30D and 31D,
   D₆/D₇ poles (60, 84) sit below E₆/E₇ (72, 126). E₆/E₇ have no integral coordinates at norm 32 and would have to
   clear the mid-layer tails at 45°. Not pursued.

## 3. Ceiling inside the family (new, routine)
**Lemma.** A deleted vector a lies in at most 8 caps.
*Proof.* Every such cap p = (b,t) has p·(a,0) = 24. Hence p = (3/4)(a,0) + q with q ⊥ (a,0) and abs(q)² = 14.
From p·p' ≤ 16 we get q·q' ≤ −2. Vectors of norm 14 with pairwise products ≤ −2 number at most 1 + 14/2 = 8.
Equality is a regular simplex. ∎

So caps ≤ 4·abs(D), and cap gain ≤ 3·abs(D). The records are at caps ≈ 2·abs(D) (gain ≈ abs(D)). Equality holds for a
**tesseract**: the 16 vectors (±u₁±u₂±u₃±u₄)/2 for an orthogonal minimal frame whose half-sum lies in L (a D₄ minus
one of its frames). It carries all 32 edges, one tail frame bijectively, 64 caps, 16 deletions: +48.
It is internally consistent. It needs 4 orthogonal tail lines, so k ≥ 4 (28–31D).
Counted: 7 172 550 tesseracts through each minimal vector.

Why it is worth testing: at the records' density (~2 000 midpoints per tail line), tesseracts would give 1.5× the
cap gain. That would be 28D ≈ 234 000 and 31D ≈ 275 000.

Why it fails:
- **Shared directions kill it.** Tesseracts sharing a direction u on one line put ±u/2 into every midpoint on their
  other lines. That forces the face centres on those lines to be an obtuse set in u⊥ (at most 46), and allows at
  most 7 tesseracts per shared direction (orthogonal triples). Sharing two directions allows at most one completion
  per group.
- **Automatically compatible families are tiny.** They are: orthogonal 4-spaces (≤ 6); tesseracts inside one
  orthogonal 24-frame (≤ 46 by pair packing; the coordinate frame and the three E₈-block frames each give only 6); D₄'s whose
  cross products are all 0 or ±8. The last condition is a dense conflict graph: by a heuristic Hurwitz count,
  each quaternionic line conflicts with about 4 600 of 8 190.
- **Measured** (scripts/tess_greedy.py and its square analogue, 28D, 300–600 s randomized greedy):

| greedy packing | units | gain |
|---|---|---|
| squares | 202 | +808 |
| tesseracts, one tail frame | 9 | +432 |
| tesseracts, three tail frames | 12–13 | +576 to +624 |

  Random packing is far below PackingStar for both.
- **Completing PackingStar's groups** (scripts/tess_completion.py, D28_221864). Any two squares of a group, with
  their antipodes, form a tesseract whose 16 missing edges cost no deletion. But all 42 480 candidate completions
  conflict with ≥ 378 existing points (median 1 120), spread over ≥ 4 other full groups.

## 4. Free additions on the current records (exact, none found)
- Free squares (scripts/freesq.py): all 44 squares of every direction pair meeting at a vertex, with 4 equator
  vertices, under every orthogonal tail-line pair (one line per direction). Candidates checked: 27D 1 846, 28D 2 478, 29D 3 978, 30D 5 540,
  31D 5 140. **0 free.**
  An earlier "2 free squares in 28D" was invalid: both directions on one tail line give b·b' + t·t = 24 at the
  shared vertex.
- Free mid points, 31D: 2 064 frame-mate bodies × 574 norm-16 tails: **0 free**.
- Signed squares (scripts/freesq_signed.py): every cap is an individual point, with any tails and any of the
  4 edges empty; the constraints inside the square are kept. Best net gain 0 in 28D (2 476 of 2 478 candidates
  take no cap at all) and −4 in 26D, 27D and 31D (all 1 054, 1 846 and 5 140 candidates fully blocked).

## 5. The class problem (old cap family, 25–31D and 38D large codimension)
A class is a set of minimal lines with abs(cos) ≤ 1/4. Known best: 248 lines. LP bound: 425 lines.
- **F₂-linear classes have at most 24 lines.** Take a coset of M with 2L ⊂ M ⊂ L and M of minimum ≥ 48.
  The character sum of the 98 280 minimal classes is two-valued: 4072 on minimal classes, −24 elsewhere
  (Parseval checked). This forces abs(T₂ ∩ V^⊥) = 3·2^(15−d) − 24, and every coset then holds at most 24 lines,
  i.e. one frame (scripts/f2_class_bound.py).
- **F₃ (Eisenstein) cosets** always contain the zero-sum triples {u, ωu, ω²u}. Taken as line systems they are never
  classes. Taken as vector sets they give about 270 (KNOWLEDGE.md §155.1), far below 496.
- So sublattice cosets cannot attack 248.

## 6. Other dimensions
- **36–37D (Coxeter–Todd):** points are triples of K₁₂ blocks labelled in K₁₂/θK₁₂ = F₃⁶. This is already a mod-3
  code construction; its "equator" is a set of cells, not a Leech shell. Nothing to swap.
- **40D (affine shells):** the core is the norm-54 shell of one class mod 3Λ, with differences of norm ≥ 54. This is
  the same trick, already used. Measured slack: C–C products 12 against a bound of 15, B–C 15 against 16.4,
  A–C 18 against 20.1. A further layer of other norms might fit. Not attempted: it needs enumeration of the K₄₀
  shells.
- **41–46D (P48n heads):** body families are P48n shells by construction. No search-built Leech component.
- **49–60D (P48p classes, base 7 783 lines):** an F₂-linear class has at most max(−f̂) vectors.
  Exact design moments of the 11-design give f̂ on the classes of minimal norm 6, 8 and 10:
  +1 833 472, +380 416, +64 000. Their distributions are fully determined; norm 6 gives 23 766 960 / 12 608 784 /
  1 678 887 / 36 848 / 1.
  Norms 12–28 leave one free parameter. The range for norm 12 is [−4 608, 396 800]; for norm 24 it is
  [−51.25M, 38.1M] (scripts/p48_fourier.py, scripts/p48_fhat_range.py).
  A coset of an index-2^k minimum-8 sublattice (2P ⊂ M) holds 52 416 000/2^k vectors on average.
  Beating 7 783 lines on average needs k ≤ 11: normalized minimum ≥ 5.82 against P48p's 6, which is implausible.
  At k = 13 (normalized minimum 5.50) one exceptional coset would need 2.4× the average.
  F₃ cosets (θP ⊂ M) contain ω-triples, so they keep one vector per orbit. Beating 7 783 then needs index ≤ 3⁷.
  In Leech, avoiding subspaces reach only half dimension; if P48p behaves the same, the index is ≥ 3¹².
  Not computed (the 26.2M-line shell exceeds the 500 MB budget). Lead left open.

## 7. Near misses
None above a record. Closest structural miss: tesseract completion. Each completion adds +32 caps at no deletion,
but is blocked by ≥ 378 conflicts. A design that leaves room for completions would need to give up about
1 120 / 32 ≈ 35× more than it gains.

## Files (Riemann/Deep/Kissing2/scripts; NumPy, exact integers; inputs are downloaded, not stored)
- decode_ps.py: layers, edge-midpoint check, r_a histogram, squares and groups (10–30 s).
- freesq.py: free squares (20 s–6 min).
- freesq_signed.py: free squares with signed caps (13–64 s).
- tess_completion.py: completion conflicts (16 min).
- tess_greedy.py: square and tesseract greedy.
- f2_class_bound.py: the 24-line bound.
- p48_fourier.py, p48_fhat_range.py: P48p character sums.

Inputs:
- github.com/alexlegeartis/KissingNumbers @1f9e23d: verifications/improved/dim27-31-joint-refills/data/coordinates{28,29,30}_*.npz
  and dim27-31-algebraic-completion/data/coordinates{27,31}_*.npz.
- github.com/CDM1619/PackingStar @477df34: 25D-31D/new_bounds_20261028/D28_221864_integer_coordinates.npz.

