# 3D snaky-corridor moving sofa: certified V*(K_4) >= 1.8755015

**Statement.** There is an explicit compact connected body B of volume >= 1.8755015382 and an explicit continuous
rigid motion carrying B from the entry arm to the exit arm of the snaky corridor K_4 while staying in K_4.
Hence V*(K_4) >= 1.8755015. Previous record: 1.8558404 (arXiv:2610.14051, Oct 2026), itself resting on the
pipeline's unproved points-to-cubes counting rule; before that 1.81 (arXiv:2511.02864, Problem 6.64).
The new bound does not use that rule: the body is an explicit union of polytopes, each checked against the whole
continuous motion, and every rounding error is bounded (section "Rounding").

**Corridor.** K = (-inf,0]x[0,1]x[0,1] u [0,1]x[0,5]x[0,1] u [0,1]x[5,6]x[0,inf) = P1 u P2 u P3 with the
boxes P1 = (-inf,1]x[0,1]^2, P2 = [0,1]x[0,6]x[0,1], P3 = [0,1]x[5,6]x[0,inf).

**Motion.** Built from the published run-II path of 2610.14051 (w5_night2_path.json, 512 poses incl. the fixed
ends (-3,0,0,0,0,0) and (0,6,8,90,90,90)), refined exactly as the pipeline does at its tolerances (0.005, 0.05 deg)
to 7352 poses. The certificate stores the pose table: translations T_k (bit-identical to the pipeline's
refinement, re-checked by verify.py) and float quaternions Q_k (from the Euler angles, Rx(roll)Ry(pitch)Rz(yaw)).
Orientation k is R(Q_k), the exact rotation matrix of Q_k/|Q_k| (a rational function of Q_k, so no trigonometry
enters the proof). Between consecutive poses: linear translation and uniform rotation about the fixed axis of
R(Q_k+1)R(Q_k)^T (the pipeline instead uses an x,y,z,rotate staircase; the simultaneous motion keeps ~0.003 more).

**Body.** 195,220 sheared prisms over a quadtree of the cross-section (y,z): base grid 100x100, up to 5
refinements (finest cell 1/3200), lattice shifted by 0.382/3200. Prism = {(x,y,z) : (y,z) in cell,
A(y,z) <= x <= B(y,z)} with A, B affine (a convex polytope; its end faces follow the tilted surfaces of the carved
body; an axis-aligned prism grid at 1/400 gives only 1.8638).

**Why it is rigorous.**
1. Chord lemma: during one step a body point q stays within |q| th^2/8 of the chord between its end positions
   (th = rotation angle; translation is linear; linear-interpolation error of the arc t -> e^{i t th}). If both end
   positions lie in the same box P_i shrunk by that amount, the whole trajectory lies in P_i, a subset of K.
2. Slices: for fixed s the slice {x = A + s(B-A)} is a parallelogram; its image lies in a box iff its 4 corner
   images do, and each corner image is affine in s. So the set of s passing a check in box i is an interval; a
   prism survives iff, for every step and the start pose, the union over boxes covers [0,1].
3. Geometry: B > A at all cell corners; cells are distinct dyadic squares (none inside another); prisms in one cell
   are disjoint; at the start pose every vertex has x <= 0 (entry arm), at the end pose z >= 1 (exit arm).
4. Connectedness: union-find over face-adjacent prisms that share a point (x-ranges overlap at an exact rational
   point of the shared face): 7 components; only the largest (volume >= 1.8755015382) is claimed.

## Rounding
Approach: (a) exact or outward-rounded arithmetic everywhere except the floating-point core; (b) a written
worst-case bound for the core (an interval version of the core would cost 3-4x its 78 s and break the 2-minute
budget). Model: IEEE-754 double, round to nearest, u = 2^-53, |fl(a op b) - a op b| <= u|a op b| + tau
(tau = 2^-1074, underflow); C compiled with -ffp-contract=off (no fused operations).

(a) In verify.py, exact (fractions.Fraction): R(Q_k) and the max error of its float roundings
(delta = 5.6e-17); the step angle via cos(th/2) = C/(|Q1||Q2|), C = Q1.Q2 > 0, and Lagrange's identity
tan^2(th/2) = S/C^2, S = sum_{i<j}(Q1_i Q2_j - Q1_j Q2_i)^2, so with th <= 2 tan(th/2) the sag bound
eps_k = QR S/(2C^2) is evaluated exactly and rounded up (max 5.5e-7); the final volume (each prism's B0-A0
rounded down times the exact area, exact rational sum, floor to 7 decimals). Outward-rounded intervals (np.nextafter
after every operation, valid since a rounded result is within half an ulp): prism corners and vertices, the bounds
QR = 3.5230003 >= |q| and |x| <= 4, |y|,|z| <= 1.0002, the errors of the corner values handed to the core
(<= 1.8e-15 for x, 3.7e-15 for B-A, 3.4e-16 for y,z), start/end containment, B > A, disjointness, connectivity.

(b) The core (verify_core.c) uses these bounded inputs (asserted at run time: also |T| <= 8, block margins <= 4).
1. Per corner, coordinate and map it computes a~ = ((R0 x + R1 y) + R2 z) + T and b~ = R0 (B-A). Against the exact
   a, b (exact R(Q), exact corners): |a - a~| <= E_a = delta(X+Y+Z+e) + (1+delta)(e_x+e_y+e_z) + g4(|R|(X+Y+Z+e)+|T|) + 3tau,
   |b - b~| <= E_b = delta(D+e_d) + (1+delta)e_d + u(D+e_d) + tau, g4 = 4u/(1-4u).
2. Endpoints: for b~ > 0 it keeps s <= s~ = fl(fl(H' - a~)/b~), H' = fl(H - f~). With fl(H'-a~) = (H'-a~)(1+d1),
   s~ = (..)(1+d2) + eta: a~ + b~ s~ = H' + (H'-a~)((1+d1)(1+d2)-1) + b~ eta, so for every s in [0,s~] n [0,1]:
   a + b s <= H' + E_div + E_a + E_b, E_div = (|H'|+|a~|)(2u+u^2) + |b~| tau. Same for b~ < 0 (s >= s~), b~ = 0
   (needs a~ <= H') and lower bounds; s~ = +-inf (overflow) gives the right one-sided answer; set operations on
   s-intervals are comparisons only.
3. Margin: f~ = fl(eps + ETA) and H' = fl(H - f~) give a + b s <= H - eps - (ETA - E_core),
   E_core = u(eps+ETA) + u(6+F) + E_div + E_a + E_b = 1.90e-14, ETA = 5e-10. So every check the core passes holds
   exactly, with the true sag bound.
4. Block skips use one representative map with margin F = (dev + max eps)(1+1e-12) + ETA. The exact distance of any
   map of the block from it is <= dev(1+2 g5) + 2 delta(X+Y+Z) + 5 tau, so F - (that) - eps >= ETA(1-u) - E_block,
   E_block = 1.96e-14 (computed like E_core).
verify.py evaluates E_core and E_block exactly from the measured input bounds and stops unless both are below
ETA/1000; they are ETA/25,491. Assumptions left: IEEE-754 double arithmetic in numpy and C, correctly rounded
float(Fraction) and sqrt (CPython, IEEE).

**Verification.** python3 verify.py w5_night2_path.json certificate_w5_night2.npz -> 195,220 / 195,220 prisms
pass every check, 7 components, largest >= 1.8755015382 (all prisms 1.8755016314), CERTIFIED V*(K_4) >= 1.8755015;
90 s, 168 MB, one background-throttled core. Negative controls: 1,500 prisms pushed out by 0.01 in x -> 1,497
rejected; a duplicated prism -> rejected (not disjoint); margin ETA = 1e-15 -> rejected by the error-bound check.
Monte Carlo of the true continuous motion (2,500 random points, half on prism faces, 17 time samples per step,
7,351 steps): 0 points leave K.

**Findings along the way.**
- The published paths are tuned to the pipeline's lattice. Unbiased (jittered random line) estimates of the true
  carved volume: w5_night2 1.8811, w6_night3 (the 1.8558 record path) 1.8815, w5_night1 1.8782; estimates on
  lattice-aligned lines differ by up to 0.01. w6_night3 carries comb-like slivers on the lattice plane y = 0.36
  (its gradient polish exploited the h = 0.001 lattice); certification of it is harder (1.8606 at a coarse setting).
- The pipeline's staircase costs about 0.003 against the simultaneous motion; finer tolerances add only ~0.0002.
- Remaining gap for this path: about 0.006 (true volume ~1.881). Headroom: polish a path against an unbiased
  objective (random line offsets) and certify it the same way.

**Files.** w5_night2_path.json (path), certificate_w5_night2.npz (6.7 MB, sha256 da63fdad...1212b: planes = 6
coefficients per prism; cells = level, j, k; meta = base grid, levels, offset, tolerances; T, Q = pose table),
verify.py + verify_core.c (verifier), build/ (construction: carve2.c, sofa2.py, adapt2.py, adapt4.py,
run_adapt4.py, mkcert.py, mctest.py).

**Rebuild** (about 9 min, < 250 MB): cd build && cc -O3 -o carve2 carve2.c -lm && mkdir -p tmp &&
python3 run_adapt4.py --path ../w5_night2_path.json --nc0 100 --L 5 --tau_rel 0.01 --chunks 3 --minvol 1e-8
--save tmp/out.npz && python3 mkcert.py tmp/out.npz ../w5_night2_path.json ../certificate_new.npz
