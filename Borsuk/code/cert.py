"""Rigorous diameter bound for Omega ∩ cell (Borsuk P8).  Omega = balls B(c, sqrt(r2)) (rational c, r2) ∩ halfspaces;
cell = halfspaces.  Polytope enclosure K: cell and linear cuts + tangent halfspaces <a, x> <= <a, c> + U(r2 |a|^2) of
each ball (U = rational upper bound of the square root), refined by cutting planes at vertices outside a ball.
Final check in exact rationals: vertices by graph traversal (pivoting) from one vertex, after an outward perturbation of
the right-hand sides (K only grows, and becomes simple); diam(K)^2 = max over vertex pairs.  Returns that bound."""
import numpy as np
from fractions import Fraction as F
from math import isqrt
from scipy.spatial import HalfspaceIntersection
from scipy.optimize import linprog

def sqrt_up(q):                      # rational q >= 0 -> rational >= sqrt(q), denominator 2^30
    s = isqrt(int(q * (1 << 60)) + 1) + 1
    return F(s, 1 << 30)

def tangent(c, r2, a):                # a rational vector
    na2 = sum(x * x for x in a)
    return (a, sum(x * y for x, y in zip(a, c)) + sqrt_up(r2 * na2))

def rat(v, den=1 << 24): return [F(int(round(x * den)), den) for x in v]

def float_vertices(H):                # H: list of (a, b) meaning a.x <= b
    A = np.array([[float(x) for x in a] for a, b in H]); b = np.array([float(b) for a, b in H])
    res = linprog(np.zeros(A.shape[1]), A_ub=A, b_ub=b - 1e-9, bounds=[(None, None)] * A.shape[1], method='highs')
    # interior point: Chebyshev centre
    nrm = np.linalg.norm(A, axis=1)
    rc = linprog(np.r_[np.zeros(A.shape[1]), -1], A_ub=np.c_[A, nrm], b_ub=b, bounds=[(None, None)] * A.shape[1] + [(0, None)], method='highs')
    x0 = rc.x[:-1]
    hs = HalfspaceIntersection(np.c_[A, -b], x0)
    return hs.intersections

def enclose(balls, lin, rounds=12, tol=1e-12):
    H = list(lin)
    for c, r2 in balls:
        R = float(r2) ** 0.5
        for k in range(4):
            for s in (1, -1):
                a = [F(0)] * 4; a[k] = F(s); H.append(tangent(c, r2, a))
    for _ in range(rounds):
        V = float_vertices(H); added = 0
        for c, r2 in balls:
            cf = np.array([float(x) for x in c]); d = ((V - cf) ** 2).sum(1)
            out = np.flatnonzero(d > float(r2) * (1 + 1e-9))
            for i in out[np.argsort(-d[out])][:40]:
                H.append(tangent(c, r2, rat((V[i] - cf) / np.sqrt(d[i])))); added += 1
        if not added: break
    return H

def exact_diam2(H, eps_pert=F(1, 1 << 36)):
    m = len(H); A = [list(a) for a, b in H]; b = [bb + eps_pert * (k + 1) / m for k, (a, bb) in enumerate(H)]   # outward
    V0 = float_vertices(list(zip(A, b)))
    x = V0[0]; Af = np.array([[float(t) for t in a] for a in A]); bf = np.array([float(t) for t in b])
    act = list(np.argsort(np.abs(Af @ x - bf))[:4])
    def solve(S):
        M = [A[k][:] + [b[k]] for k in S]
        for col in range(4):
            p = next(r for r in range(col, 4) if M[r][col] != 0); M[col], M[p] = M[p], M[col]
            iv = 1 / M[col][col]; M[col] = [t * iv for t in M[col]]
            for r in range(4):
                if r != col and M[r][col] != 0:
                    f = M[r][col]; M[r] = [s - f * t for s, t in zip(M[r], M[col])]
        return [M[r][4] for r in range(4)]
    def inv_cols(S):
        M = [A[k][:] + [F(int(i == j)) for j in range(4)] for i, k in enumerate(S)]
        for col in range(4):
            p = next(r for r in range(col, 4) if M[r][col] != 0); M[col], M[p] = M[p], M[col]
            iv = 1 / M[col][col]; M[col] = [t * iv for t in M[col]]
            for r in range(4):
                if r != col and M[r][col] != 0:
                    f = M[r][col]; M[r] = [s - f * t for s, t in zip(M[r], M[col])]
        return [[M[r][4 + j] for r in range(4)] for j in range(4)]   # column j of A_S^{-1}
    S0 = tuple(sorted(int(k) for k in act)); v0 = solve(S0)
    assert all(sum(p * q for p, q in zip(A[k], v0)) <= b[k] for k in range(m)), 'start vertex infeasible'
    seen = {S0: v0}; todo = [S0]
    while todo:
        S = todo.pop(); v = seen[S]; cols = inv_cols(list(S))
        for jj, j in enumerate(S):
            d = [-t for t in cols[jj]]
            best = None
            for k in range(m):
                if k in S: continue
                ad = sum(p * q for p, q in zip(A[k], d))
                if ad > 0:
                    t = (b[k] - sum(p * q for p, q in zip(A[k], v))) / ad
                    if best is None or t < best[0]: best = (t, k, 1)
                    elif t == best[0]: best = (t, best[1], best[2] + 1)
            assert best is not None, 'unbounded'
            assert best[2] == 1, 'degenerate pivot'
            S2 = tuple(sorted([k for k in S if k != j] + [best[1]]))
            if S2 not in seen:
                seen[S2] = [p + best[0] * q for p, q in zip(v, d)]; todo.append(S2)
    Vs = list(seen.values())
    return max(sum((p - q) ** 2 for p, q in zip(v, w)) for i, v in enumerate(Vs) for w in Vs[i + 1:]), len(Vs)
