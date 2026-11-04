#!/usr/bin/env python3
"""Optimise the 8 parameters of the simplex-product family V_n (see family.py) for one dimension n.
Needs numpy and scipy (SLSQP).  The output row can be appended to family_params.txt and checked exactly with
family.py.  Usage:  python3 family_opt.py N [starts]
Model: maximise t subject to (17 edge classes of the lemma) >= t and (all 36 classes) <= 1, where every class value
is the exact quadratic form of family.py evaluated in floating point."""
import math, sys
import numpy as np
from scipy.optimize import minimize
from fractions import Fraction
import family as fam

def fval(m1, m2, x, cl):
    return float(fam.class_d2(m1, m2, tuple(Fraction(v) for v in x), cl))

def cdist(m1, m2, x, cl):
    X, Y, s1, s2 = cl
    cf = {'P': (x[0], 0.0), 'Q': (0.0, x[1]), 'A': (x[2], -x[3]), 'B': (-x[4], x[5]), 'C': (x[6], x[7])}
    (a1, b1), (a2, b2) = cf[X], cf[Y]
    def part(u, v, rel, m):
        if rel == '-': return u*u + v*v
        g = 1.0 if rel == 'e' else -1.0/m
        return u*u + v*v - 2*u*v*g
    return part(a1, a2, s1, m1) + part(b1, b2, s2, m2)

def solve(n, starts=16, seed=0):
    m1, m2 = n//2, n - n//2
    cls = fam.classes()
    rA, rB = math.sqrt(m1/(2*(m1+1))), math.sqrt(m2/(2*(m2+1)))
    base = np.array([rA, rB, 0.99*rA, 0.14*rB, 0.15*rA, 0.98*rB, 0.62*rA, 0.78*rB])
    rng = np.random.default_rng(seed)
    best = None
    for s in range(starts):
        x0 = base*np.exp(rng.normal(0, 0.05 if s else 0.0, 8))
        z0 = np.r_[x0, 0.8]
        cons = [{'type': 'ineq', 'fun': (lambda z, c=c: cdist(m1, m2, z[:8], c) - z[8])} for c in fam.EDGES]
        cons += [{'type': 'ineq', 'fun': (lambda z, c=c: 1.0 - cdist(m1, m2, z[:8], c))} for c in cls]
        res = minimize(lambda z: -z[8], z0, constraints=cons, method='SLSQP', options={'maxiter': 1000, 'ftol': 1e-16})
        x = res.x[:8]
        r = min(cdist(m1, m2, x, c) for c in fam.EDGES)/max(cdist(m1, m2, x, c) for c in cls)
        if best is None or r > best[0]: best = (r, x)
    return m1, m2, best

if __name__ == '__main__':
    n = int(sys.argv[1]); starts = int(sys.argv[2]) if len(sys.argv) > 2 else 16
    m1, m2, (r, x) = solve(n, starts)
    row = '%d %d %d %s' % (n, m1, m2, ' '.join('%.12f' % v for v in x))
    p = tuple(Fraction(w) for w in row.split()[3:])
    rex, ok = fam.check(n, m1, m2, p)
    print('# float ratio %.12f; exact after rounding: ratio^2 = %s, above sigma_%d: %s' % (math.sqrt(r), rex, n, ok))
    print(row)

