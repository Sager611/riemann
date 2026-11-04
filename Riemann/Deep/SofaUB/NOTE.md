# 3D snaky-corridor sofa: certified V*(K_M) <= 3.1 for every M > 2.09 (was 3.4083 for M > 2.2)

**Statement.** Let u = fl(2.09) = 2.08999999999999985789... (a double; c = 1 + u is exact). Every compact connected
body that weakly traverses K_M with M > u has volume <= 3.0999999989 < 3.1. Previous best upper bound: 3.4083 for
M > 2.2 (arXiv 2610.14050, Theorem C, double precision without rounding control). Best lower bound: 1.8755015
(Riemann/Deep/Sofa3D). Covers M = 4 (Problem 6.64 of arXiv 2511.02864).

**Reduction (2610.14050, Prop. 6.1; re-derived, correct).** For u in (sqrt2, M), c = 1+u: V* <= max{u, V_ii(u)},
V_ii(u) = sup vol(T ^ (S_X u Stub(c))) over unit-square tubes T with |w.yhat| >= mu0 = (u^2-2)/(u^2+2) meeting
W_1 = [0,1]x{1}x[0,1] and W_c = [0,1]x{c}x[0,1] (first-crossing event, deep point, entry tube, sharp tilt lemma).

**Finding: the limit of this family.** The vertical tube [0,1]xRx[0,1] is admissible with volume c, so
V_ii(u) >= 1+u: at the paper's u = 2.2 no bound below 3.2 is possible. Numerical maximization (exact volumes,
Nelder-Mead, 300 starts per u) finds one other competitor, a shallow tube raking the entry arm (theta ~ 64 deg from
yhat, w_z ~ 0, roll ~ 45 deg): S(u) = 3.231 (u=2.0), 3.068 (2.08), 3.031 (2.10), 2.998 (2.12), 2.895 (2.2).
Hence min_u max(1+u, S(u)) ~ 3.073 at u ~ 2.073 is the best any single-event certificate can reach. We certify
u = 2.09, T = 3.1. (The paper's top leaves, theta ~ 0.4 rad, are evaluator slack: true values there are ~2.8.)

**Slice formulation.** k = (w_x, w_z)/w_y. Horizontal slices of a tube: Sigma(y) = X + (y-ya)k + A_k R(t) Q0,
Q0 = [-1/2,1/2]^2, A_k = I + k k^T/(1+s), s = sqrt(1+|k|^2) (the oblique projection stretches by s along k),
R(t) = rotation with cos = (1-t^2)/(1+t^2), sin = 2t/(1+t^2), t in [0,1] (all rolls, rational: no trigonometry).
V = int_0^c area(Sigma(y) ^ R_y) dy with R_y = (-inf,1]x[0,1] (y <= 1), [0,1]^2 (1 < y <= c). Windows: Sigma(1),
Sigma(c) meet [0,1]^2. Tilt: |k|^2 <= 8u^2/(u^2-2)^2. z-mirror: k_z >= 0. Five parameters (kx, kz, t, X1, X2),
X = center of Sigma(ya), ya = 0.75. Near-vertical bound: V <= u + sqrt(1+|k|^2).

**Box bound (sofaub.c).** For a parameter box, every slice lies in the dilated center slice
D(y) = {x : |cross(x-X(y), g2)| <= det/2 + PA0 + PA1|y-ya|, |cross(g1, x-X(y))| <= det/2 + PB0 + PB1|y-ya|},
g1, g2 = edge vectors of the center slice, det = s. Push-outs: offset (support function of the X box), slope
(|y-ya| times that of the k box), roll (det(d + d^2/2)/2, d <= 2 r_t/(1+t_min^2)), stretch (first-order term
sum_a r_ka |cross(dA/dk_a v_j, g)| plus remainders D2|dk|^2|g|/2 and L_A|dk| d |g| with explicit bounds L_A, D2 on
||DA||, ||D^2 A|| over the box; Taylor along segments inside the box). So V <= int_0^c area(D(y) ^ R_y) dy.
A box is infeasible if D(1) or D(c) is separated from [0,1]^2 (separating axis) or min|k|^2 > K^2.

**Exact integration.** On each sub-range of {0, ya, 1, c} the vertices of D(y) move affinely, so
area(D(y) ^ R_y) is piecewise quadratic; it changes formula only when a vertex of D crosses an edge segment of R or
a corner of R crosses an edge segment of D. These events are enclosed by intervals; between them the 3-node
interpolatory rule is exact; each event interval (width ~1e-15) adds width x max area. Areas by Green's theorem
(clip each edge to R; right edge of R inside D), so no polygon combinatorics. Tilt, window and coverage tests and
the near-vertical bound are interval checks too.

## Rounding
IEEE-754 double, round to nearest, -ffp-contract=off, no fast-math; only + - * / sqrt (correctly rounded).
Every interval operation is evaluated in round-to-nearest and widened outward by |x| 2^-50 + 2^-1074, at least
4 ulp, more than the rounding error (<= 2^-53|x|(1+2^-52) for normal results, <= 2^-1075 for subnormal), so every
interval encloses the exact value. A division by an interval containing 0 sets an error flag: the leaf then uses
the near-vertical bound (never happened: 0 degenerate leaves). The split heuristic does not affect soundness.

## Certificate and verification
Root grid 16 x 8 x 4 x 4 x 4 = 8,192 boxes: kx in [-K,K], kz in [0,K], t in [0,1], X(ya) in the range forced by
both windows (slice half-extent <= s/sqrt2). Depth-first tree, one byte per node (0 leaf, d+1 split dim d at the
midpoint): 32,247,596 nodes, 16,127,894 leaves = 7,434,183 infeasible + 579 near-vertical + 8,693,132 slice
integral; max leaf bound 3.099999998842. codes.bin.gz 1.5 MB, sha256(raw) 3f888cf6...10017fd835.
python3 verify.py (stdlib only + cc): [1] exact constants (fractions); [2] compile, replay every leaf: PASS,
372 s, 89 MB, one background core; [3a] 120 sampled leaves re-integrated exactly in rational arithmetic by an
independent algorithm (all breakpoints, Simpson): C bound - exact >= 1.4e-12; [3b] Monte Carlo: random admissible
tubes in sampled boxes lie below the leaf bound, 536 samples in infeasible boxes are all non-admissible;
[3c] root coverage proved in exact rational arithmetic. Development tests: 2e5 area enclosures vs float clipping,
12,000 random tubes inside random boxes below the box bound, vertical tube volume exactly c.
Search: ./sofaub search 2.09 3.10 0.75 16 8 4 4 codes.bin (929 s). Also certified: u = 2.12, T = 3.13 (3.3M leaves).

## Reproduction of 2610.14050
Their archived frontier (data/wu3/bnb_u2.2_g0.02_final.npz: 49,513 leaves, 15,067 infeasible) tiles the parameter
box (volume ratio 1.0000000000); frontier max 3.408568; 50 leaves (top 10, 20 random feasible, 20 infeasible)
re-evaluated with their evaluator never exceed the stored bound; their polish (g = 0.01, top 27 leaves) gives
3.408292027 = Theorem C. Needs a shim: numpy >= 2 dropped np.cross for 2-vectors (their geometry.py uses it).

## Headroom and next steps
Single-event optimum ~3.073 (u ~ 2.073); 3.09 would need ~10x more leaves (first-order box bounds in 5D). The
vertical tube is a 1 x 1 x c rod that cannot turn the first corner; only coupling several events (e.g. axial
profiles at two first-crossing heights) can see that, which is the route to bounds well below 3.

**Novelty.** arXiv listings Oct-Nov 2026 (math.MG, cs.CG, math.CO, math.OC, math.NA, math.GM; "sofa" hits:
2610.13955, 2610.14050, 2610.14051 only) and Zenodo: no newer upper bound.

**Files.** NOTE.md, REPORT.md, meta.json, codes.bin.gz (tree), sofaub.c (search / verify / dump / roots),
verify.py, exact_check.py (exact rational integrator), mc_check.py (admissibility, tube volume, root coverage).
