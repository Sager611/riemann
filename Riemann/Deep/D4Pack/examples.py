#!/usr/bin/env python3
"""examples.py -- explicit packings that break the reduction steps of arXiv 2609.25120v2.

Frame: 'native' = sqrt2*D4 scaled down by sqrt2, i.e. D4 roots alpha (norm^2 2) are the
contacts, packing distance sqrt2, volumes x4 give the unit-ball frame (24-cell = 8).
Unit-frame distance d <-> native d/sqrt2.  Voronoi halfspace of neighbour y: x.y <= |y|^2/2.
'Corner' = every neighbour pulled radially to packing distance (Lemma radial-new, TeX 1032-1040).
"""
import itertools, sys
from fractions import Fraction as Fr
import numpy as np
from voro4 import exact_volume, float_cell, d4_roots

R = [tuple(r) for r in d4_roots()]
SQ2 = np.sqrt(2.0)

def check_packing(Y, exact=True):
    """min squared distance among {0} u Y, native frame (must be >= 2)."""
    pts = [tuple(Fr(0) for _ in range(4))] + [tuple(Fr(v) for v in y) for y in Y] if exact else \
          [np.zeros(4)] + [np.asarray(y, float) for y in Y]
    best = None
    for p, q in itertools.combinations(pts, 2):
        d2 = sum((p[k] - q[k]) ** 2 for k in range(4))
        best = d2 if best is None or d2 < best else best
    return best

def true_cell_exact(Y):
    A = [[Fr(v) for v in y] for y in Y]
    b = [sum(Fr(v) ** 2 for v in y) / 2 for y in Y]
    return 4 * exact_volume(A, b)

def corner_float(Y):
    Y = np.asarray([[float(v) for v in y] for y in Y])
    U = Y / np.linalg.norm(Y, axis=1)[:, None]
    vol, _ = float_cell(U, np.full(len(U), SQ2 / 2))   # plane at native distance sqrt2/2
    return 4 * vol, U

def max_inner(U):
    G = U @ U.T
    np.fill_diagonal(G, -2)
    return G.max()

def report(name, Y, true_exact=True):
    d2 = check_packing(Y, exact=true_exact)
    dist = [float(np.sqrt(sum(float(v) ** 2 for v in y))) * SQ2 for y in Y]   # unit frame
    if true_exact:
        tv = true_cell_exact(Y)
        tvs = f"{float(tv):.6f} (= {tv})" if tv.denominator < 10**8 else f"{float(tv):.6f} (exact rational)"
    else:
        A = np.asarray(Y, float)
        tv = 4 * float_cell(A, 0.5 * (A * A).sum(1))[0]
        tvs = f"{tv:.6f} (float)"
    cv, U = corner_float(Y)
    print(f"\n[{name}]  neighbours {len(Y)}; unit-frame distances {min(dist):.4f}..{max(dist):.4f};"
          f" packing min dist^2 (native, need >= 2) = {float(d2):.6f}{' exact' if true_exact else ''}")
    print(f"   true cell volume {tvs}")
    print(f"   corner: max <u,u'> = {max_inner(U):.4f} (contact config needs <= 1/2), "
          f"{len(U)} directions, volume {cv:.6f}")
    return tv, cv

print("E0  sanity: D4 contacts ->", true_cell_exact([list(r) for r in R]))

# E1: shell localisation (Lemma shell-new, TeX 985-1002) is false:
# deletion of root alpha0 plus a centre at unit distance 2.9 > 2 sqrt2 on the ray of alpha0.
a0 = R[0]
Y1 = [list(r) for r in R[1:]] + [[Fr(29, 20) * v for v in a0]]      # native 1.45*alpha0 = unit 2.9
tv1 = true_cell_exact(Y1)
print(f"\n[E1 far centre cuts the cell]  deletion alone: {true_cell_exact([list(r) for r in R[1:]])};"
      f" + centre at unit distance 2.9 (> 2*sqrt2 = 2.8284): {tv1} = {float(tv1):.6f};"
      f" packing min dist^2 = {check_packing(Y1)}")

# E2: the authors' 49-ball example: 48 centres at unit distance 1/sin(22.5 deg) in the F4 directions.
Rr = 1 / np.sin(np.pi / 8)
dirs = [np.array(r, float) / SQ2 for r in R]
dirs += [np.eye(4)[i] * s for i in range(4) for s in (1, -1)]
dirs += [np.array(s) / 2 for s in itertools.product((1, -1), repeat=4)]
Y2 = [d * Rr / SQ2 for d in dirs]   # native
report("E2 49 balls, F4 directions, unit distance %.4f" % Rr, Y2, true_exact=False)

# E3: rational 49-ball variant: roots at 5/4*alpha (unit 2.5), short at 39/20*sigma (unit 2.758).
sig = [tuple(Fr(int(v)) for v in np.eye(4)[i] * s) for i in range(4) for s in (1, -1)]
sig += [tuple(Fr(v, 2) for v in s) for s in itertools.product((1, -1), repeat=4)]
Y3 = [[Fr(5, 4) * v for v in r] for r in R] + [[Fr(39, 20) * v for v in s] for s in sig]
report("E3 49 balls, rational radii", Y3)

# E4: a 23-direction pattern that is not a contact configuration and whose corner is below 8:
# delete adjacent roots (1,1,0,0),(1,0,1,0); add a centre at native (2/3)(2,1,1,0) (unit 2.3094).
rem = {(1, 1, 0, 0), (1, 0, 1, 0)}
Y4 = [list(r) for r in R if r not in rem] + [[Fr(4, 3), Fr(2, 3), Fr(2, 3), Fr(0)]]
report("E4 two roots replaced by one centre between them", Y4)

# E5: single-deviation theorem (TeX 323-331) via the corner: 24 root neighbours at unit 2.3,
# one deep-hole neighbour at unit 2.8 (native 1.98 e1): hypothesis holds, corner < 8.
Y5 = [[Fr(23, 20) * v for v in r] for r in R] + [[Fr(99, 50), 0, 0, 0]]
report("E5 24 roots at 2.3 + one non-root at 2.8", Y5)

