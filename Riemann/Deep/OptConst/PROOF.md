# C_41b <= 19/10: an upper bound for the ambidextrous moving sofa constant

Constant 41b of [teorth/optimizationproblems](https://teorth.github.io/optimizationproblems/constants/41b.html):
the largest area of a sofa that can pass a Z-shaped hallway of width 1 (one right-angled turn to the right and,
far away, one to the left).  Known: lower bound 1.64495521 (Romik 2016); upper bounds 2√2 (Hammersley-type),
2.37 (Kallus–Romik 2018, the published bound for the one-sided constant 41a, which is ≥ 41b) and 2.2195 (via Baek's
announced, unverified proof that Gerver's sofa is optimal for 41a).  Kallus and Romik (2018, end of §5) explicitly
proposed adapting their method to the ambidextrous problem; we found no such adaptation in the literature.

**Theorem.** Every ambidextrous sofa has area at most 19/10.  This holds for the Z-hallway formulation of the
41b page and also for the formulation "the same placed shape can be moved around a right turn and around a left
turn" (Romik 2016; Kallus–Romik 2018).  The proof uses no result about rotation paths (no Gerver reduction).

## Notation
* H = R × [0,1], V = [0,1] × R, L0 = (−∞,1]×[0,1] ∪ [0,1]×(−∞,1] (the hallway of Kallus–Romik: entry along H
  from the left, exit downward), Rot_ψ = rotation by ψ about the origin, C(ψ,u) = Rot_ψ(L0 + u).
* In rotated coordinates p = x cos ψ + y sin ψ, q = −x sin ψ + y cos ψ:
  C(ψ,u) = P ∪ Q with P = {p ≤ u1+1, u2 ≤ q ≤ u2+1}, Q = {u1 ≤ p ≤ u1+1, q ≤ u2}.
* T = 19/10 and ε ∈ (0, π/4) with sin ε = 1/T.  Windows Q_j = [jπ/2 + ε, (j+1)π/2 − ε] (mod 2π), j = 0..3.
* A sofa is a connected compact set S of positive area; motions are continuous families of orientation-preserving
  rigid motions g_τ(z) = Rot_{θ(τ)} z + c(τ), with θ the continuous real lift, θ(0) = 0.
* Z-hallway: horizontal segment (heading east), clockwise corner, long vertical segment, counterclockwise corner,
  horizontal segment; the corners are at distance ≥ 10·diam(S) ("sufficiently far apart").  The mirror-image Z
  (left turn first) is equivalent by reflecting the whole picture.

## Step 1 (sofa frame)
Place the sofa frame so that S ⊂ H at time 0 (first horizontal segment).  If at time τ the moved sofa lies in a
world set X, then S ⊂ Rot_{−θ(τ)}(X) + t.  The first corner (clockwise) is L0, so near it S ⊂ C(−θ, u).  The second
corner (counterclockwise, down→right) is, as a set, Rot_π(L0) + c, so near it S ⊂ C(π − θ, u).

## Step 2 (strip bounds, after Gerver)
Two width-1 strips at angle γ meet in area 1/|sin γ|.  Assume area(S) > T.  (a) When the sofa lies in the vertical
segment, S ⊂ H ∩ Rot_{−θ}(V + t), so 1/|cos θ| > T: θ ∈ (−π/2 + kπ − ε, −π/2 + kπ + ε) for a unique integer k
(unique because ε < π/4).  (b) When it lies in the last horizontal segment, S ⊂ H ∩ Rot_{−θ}(H + t), so
θ ∈ (mπ − ε, mπ + ε).

## Step 3 (two full windows are always available)
Say window j is available if for every ψ ∈ Q_j there is u with S ⊂ C(ψ,u).
Let r(τ) be a point of the moved sofa and M the (closed, nonempty) set of times at which r is on the midline of
the vertical segment; on M, (a) defines k(τ), locally constant.  Let τ' = min M, τ_end = first time the sofa lies in
the last segment, τ'' = max M ∩ [0, τ_end].  While r is above the midline the sofa is inside L0 (world), while below
it is inside the second corner's hallway.  By the intermediate value theorem θ attains every value between its
values at two times.
* k constant on M (= k).  k ≤ −1: on [0,τ'] θ covers [−3π/2+ε, 0], so −θ covers Q_0 and Q_2.  k ≥ 2: −θ covers
  Q_3 and Q_1.  k = 0: −θ covers Q_0; on [τ'',τ_end] θ goes from (−π/2 ± ε) to (mπ ± ε), so π − θ covers Q_2 (m ≥ 0)
  or Q_3 (m ≤ −1).  k = 1: −θ covers Q_3; π − θ covers Q_1 (m ≤ 0) or Q_0 (m ≥ 1).
* k not constant: there is an excursion (τc, τd), disjoint from M, with τc, τd ∈ M and k(τc) ≠ k(τd)
  (take τd = inf of the times with a different k than at τ', τc = last earlier time of M).  During it the sofa stays
  in one corner region and θ sweeps an interval of length ≥ π − 2ε between (−π/2 + kπ ± ε) values: −θ covers
  Q_0,Q_3 or Q_1,Q_2; π − θ covers Q_1,Q_2 or Q_3,Q_0.
So one of the window pairs {0,2}, {0,3}, {1,2}, {1,3} is available.  (Right/left-turn formulation: the left-turn
hallway M(L0) (M: y ↦ 1−y) equals Rot_{−π/2}(L0) + t, its hallway angles are −π/2 − θ; the same bookkeeping gives
{0,2}, {0,3}, {2,3}, {1,3}, or the single window Q_3 (both turns rotating "the wrong way").)

## Step 4 (compact domain)
Take an available angle ψ.  S ⊂ H ∩ C(ψ,u) = (H∩P) ∪ (H∩Q); H meets the full bands containing P and Q in sets of
x-extent (1+|cos ψ|)/|sin ψ| and (1+|sin ψ|)/|cos ψ|; S is connected, so after a horizontal shift (H is invariant)
S ⊂ K = [0,W] × [0,1] with W = 155/32 ≥ 29/6 (the value for cos = 21/29, sin = 20/29).
Clamping (proved by inspection of P and Q): with a = min_K p, b = max_K p, c = min_K q, d = max_K q, replacing u1 by
min(u1, b) and u2 by min(u2, d) changes K ∩ C(ψ,u) only by a null set, and u1 < a − 1 or u2 < c − 1 makes it empty.
Hence area(S) ≤ max over u_j ∈ [a_j − 1, b_j] × [c_j − 1, d_j] of area(K ∩ ∩_j C(ψ_j, u_j)) for any finite set of
angles ψ_j taken from available windows.

## Step 5 (certified branch and bound)
On a box B = [a1,b1] × [a2,b2] of translations, C(ψ,u) ⊂ Pt ∪ Qt with Pt = {p ≤ b1+1, a2 ≤ q ≤ b2+1},
Qt = {a1 ≤ p ≤ b1+1, q ≤ b2}.  The area of K ∩ ∩_j (Pt_j ∪ Qt_j) equals the inclusion–exclusion sum over
{Pt, Qt, Pt∩Qt}^k of signed areas of convex polygons.  We use the base angles φ = arctan(20/21) (cos 21/29,
sin 20/29) and π/2 − φ, both inside [ε, π/2 − ε] since 20/29, 21/29 ≥ 10/19.  For each case (the six pairs of
distinct windows, each with one corridor at jπ/2 + φ, and windows 1 and 3 alone with two corridors at jπ/2 + φ,
jπ/2 + π/2 − φ) a branch-and-bound tree (split the widest coordinate at its midpoint) covers the clamped domain;
verify.py recomputes every leaf bound exactly (integer homogeneous coordinates, Fractions) and checks ≤ 19/10:

| case | corridors (hallway angles) | leaves | max leaf bound |
|---|---|---|---|
| {0,2} natural | φ, π+φ | 513503 | 1.8999998225 |
| {0,1} {0,3} {1,2} {1,3} {2,3} | jπ/2+φ, j'π/2+φ | 2100–3006 | ≤ 1.8999824195 |
| Q_1 alone, Q_3 alone | jπ/2+φ, jπ/2+π/2−φ | 831, 833 | ≤ 1.8978989000 |

Contradiction with area(S) > 19/10.  ∎

## Remarks
* The relaxation is nearly tight for these angles: for the natural pair {0,2} there are translations with area
  1.8298 (for φ = π/4 the sup is 2√2 − 1 numerically), so 19/10 is close to what one angle per turn can give.
* With two angles per turn (36.87°, 53.13°) the floating-point search closes at 1.80 (28M leaves; not exactly
  verified, too large for the Python checker) and the relaxation value is ≥ 1.729; Kallus–Romik's convergence
  theorem suggests the method converges to the true constant (conjecturally Romik's 1.64495521).
* Rotating "the wrong way" through a corner (windows 1, 3) is very restrictive: two nearby angles already give ≤ 1.9.
