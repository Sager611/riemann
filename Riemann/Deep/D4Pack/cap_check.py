#!/usr/bin/env python3
"""cap_check.py -- thm:cap-inequality (TeX 6546): for B = {|x|_1 <= 2} in R^4 and unit u,
vol(B n {x.u >= 1}) <= 1/3, equality at vertex directions.  cap = 32/3 - vol(B n {x.u <= 1})."""
import itertools, sys
import numpy as np
from voro4 import float_cell
S = np.array(list(itertools.product((1., -1.), repeat=4)))
def cap(u):
    u = u / np.linalg.norm(u)
    v, _ = float_cell(np.vstack([S, u]), np.concatenate([np.full(16, 2.0), [1.0]]))
    return 32 / 3 - v
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
best = (-1, None)
U = list(rng.normal(size=(1500, 4)))
U += [np.eye(4)[0] + rng.normal(scale=s, size=4) for s in (1e-3, 1e-2, 5e-2, 0.1, 0.2) for _ in range(60)]
U += [np.array([1, 1, 0, 0.]), np.array([1, 1, 1, 0.]), np.array([1, 1, 1, 1.]), np.array([1, .5, 0, 0.])]
for u in U:
    c = cap(np.asarray(u, float))
    if c > best[0]: best = (c, u / np.linalg.norm(u))
print(f"vertex direction cap = {cap(np.eye(4)[0]):.12f} (1/3 = {1/3:.12f})")
print(f"{len(U)} directions: max cap = {best[0]:.12f} at {np.round(best[1], 4)}; exceeds 1/3: {best[0] > 1/3 + 1e-10}")

