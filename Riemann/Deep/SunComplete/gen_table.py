#!/usr/bin/env python3
"""gen_table.py -- regenerates the table of Lemma 6 (NOTE.md sec. 4): for each case (A/B) and
r = omega(q-1) <= 60, the core size t and the least q0 at which the clean criterion (CS) holds
with worst-case parameters (bisection; (CS) is monotone in q). Exact rational arithmetic (verify.CS_holds).
Runtime ~2-3 min."""
import os, sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from verify import CS_holds, P_r

def least_q0(case, r, t, hi=10**40):
    if not CS_holds(case, r, t, hi): return None
    lo = 2
    while hi - lo > 1:
        mid = (lo + hi)//2
        if CS_holds(case, r, t, mid): hi = mid
        else: lo = mid
    return hi

for case in 'AB':
    rows = []
    for r in range(1, 61):
        best = None
        for t in range(0, min(r, 8) + 1):
            q0 = least_q0(case, r, t)
            if q0 is not None and (best is None or q0 < best[1]): best = (t, q0)
        rows.append((r, best[0], best[1], best[1] <= P_r(r) + 1))
        print(case, "r=%2d t=%d q0=%d  q0<=P_r+1: %s" % rows[-1], flush=True)
    print(case, "uniform threshold (max q0 over r with q0 > P_r+1):", max(q0 for r, t, q0, cov in rows if not cov))
