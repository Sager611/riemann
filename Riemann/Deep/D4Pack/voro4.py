#!/usr/bin/env python3
"""voro4.py -- volume of a 4-polytope P = {x : a_i.x <= b_i} with the origin inside.

Own code, no scipy.  Two modes:
  exact_volume(A, b): A, b rational (fractions.Fraction or int); vertices found
      with a float prefilter, then recomputed and checked EXACTLY; volume by a
      pulling triangulation (cone from a vertex over facets, recursively), exact.
  float_cell(A, b): float version returning volume, facet 3-volumes and facet
      centroids (for gradients); same combinatorics with a tolerance.
Voronoi cell of centre 0 with neighbours y_j: halfspaces x.y_j <= |y_j|^2/2.
"""
import itertools
from fractions import Fraction as Fr
import numpy as np


# ---------------------------------------------------------------- exact helpers
def _det_exact(M):
    M = [list(r) for r in M]
    n = len(M)
    det = Fr(1)
    for c in range(n):
        p = next((r for r in range(c, n) if M[r][c] != 0), None)
        if p is None:
            return Fr(0)
        if p != c:
            M[c], M[p] = M[p], M[c]
            det = -det
        piv = M[c][c]
        det *= piv
        for r in range(c + 1, n):
            if M[r][c] != 0:
                f = M[r][c] / piv
                for k in range(c, n):
                    M[r][k] -= f * M[c][k]
    return det


def _rank_exact(rows):
    M = [list(r) for r in rows]
    if not M:
        return 0
    n, m = len(M), len(M[0])
    rk = 0
    for c in range(m):
        p = next((r for r in range(rk, n) if M[r][c] != 0), None)
        if p is None:
            continue
        M[rk], M[p] = M[p], M[rk]
        piv = M[rk][c]
        for r in range(n):
            if r != rk and M[r][c] != 0:
                f = M[r][c] / piv
                for k in range(c, m):
                    M[r][k] -= f * M[rk][k]
        rk += 1
        if rk == n:
            break
    return rk


def _solve_exact(A4, b4):
    M = [list(A4[i]) + [b4[i]] for i in range(4)]
    for c in range(4):
        p = next((r for r in range(c, 4) if M[r][c] != 0), None)
        if p is None:
            return None
        M[c], M[p] = M[p], M[c]
        piv = M[c][c]
        for r in range(4):
            if r != c and M[r][c] != 0:
                f = M[r][c] / piv
                for k in range(c, 5):
                    M[r][k] -= f * M[c][k]
    return tuple(M[i][4] / M[i][i] for i in range(4))


def _aff_rank_exact(pts):
    if len(pts) <= 1:
        return 0
    p0 = pts[0]
    return _rank_exact([[p[k] - p0[k] for k in range(4)] for p in pts[1:]])


# ---------------------------------------------------------------- float vertex prefilter
def vertices_float(A, b, tol=1e-9, chunk=40000):
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    m = len(b)
    scale = max(1.0, float(np.abs(b).max()))
    out = []
    it = itertools.combinations(range(m), 4)
    while True:
        cmb = np.array(list(itertools.islice(it, chunk)), dtype=np.int64)
        if len(cmb) == 0:
            break
        As = A[cmb]
        bs = b[cmb]
        det = np.linalg.det(As)
        ok = np.abs(det) > 1e-10
        if not ok.any():
            continue
        x = np.linalg.solve(As[ok], bs[ok][..., None])[..., 0]
        feas = (x @ A.T <= b[None, :] + tol * scale).all(axis=1)
        out.append(x[feas])
    if not out:
        return np.zeros((0, 4))
    X = np.concatenate(out)
    if len(X) == 0:
        return X
    key = np.round(X / (1e-7 * scale)).astype(np.int64)
    _, idx = np.unique(key, axis=0, return_index=True)
    return X[np.sort(idx)]


# ---------------------------------------------------------------- exact volume
def exact_vertices(A, b):
    Af = [[Fr(v) for v in row] for row in A]
    bf = [Fr(v) for v in b]
    m = len(bf)
    Xf = vertices_float([[float(v) for v in r] for r in Af], [float(v) for v in bf])
    verts = {}
    for x in Xf:
        vals = np.array([sum(float(Af[i][k]) * x[k] for k in range(4)) - float(bf[i]) for i in range(m)])
        T = [i for i in range(m) if abs(vals[i]) <= 1e-6 * (1 + abs(float(bf[i])))]
        # choose 4 independent tight rows exactly
        chosen = []
        for i in T:
            if _rank_exact([Af[j] for j in chosen + [i]]) == len(chosen) + 1:
                chosen.append(i)
            if len(chosen) == 4:
                break
        if len(chosen) < 4:
            raise RuntimeError("degenerate float vertex: fewer than 4 independent tight constraints")
        xe = _solve_exact([Af[i] for i in chosen], [bf[i] for i in chosen])
        lhs = [sum(Af[i][k] * xe[k] for k in range(4)) for i in range(m)]
        if any(lhs[i] > bf[i] for i in range(m)):
            raise RuntimeError("float vertex not exactly feasible")
        Te = frozenset(i for i in range(m) if lhs[i] == bf[i])
        verts[xe] = Te
    return verts


def exact_volume(A, b, return_info=False):
    """Exact volume of {x: A x <= b} (bounded, origin interior); A, b rational."""
    verts = exact_vertices(A, b)
    V = list(verts.keys())
    T = [verts[v] for v in V]
    m = len(b)
    nv = len(V)
    on = [frozenset(k for k in range(nv) if i in T[k]) for i in range(m)]
    # facets: constraints whose vertex set has affine rank 3; dedupe by vertex set
    facets = {}
    for i in range(m):
        S = on[i]
        if len(S) >= 4 and S not in facets:
            if _aff_rank_exact([V[k] for k in sorted(S)]) == 3:
                facets[S] = i
    vol = Fr(0)
    v0 = 0
    for S in facets:
        if v0 in S:
            continue
        vF = min(S)
        # 2-faces of facet S
        two = set()
        for j in range(m):
            G = S & on[j]
            if len(G) >= 3 and G != S and G not in two:
                if _aff_rank_exact([V[k] for k in sorted(G)]) == 2:
                    two.add(G)
        for G in two:
            if vF in G:
                continue
            vG = min(G)
            edges = set()
            for k in range(m):
                E = G & on[k]
                if len(E) == 2:
                    edges.add(E)
            for E in edges:
                if vG in E:
                    continue
                a, c = sorted(E)
                P0 = V[v0]
                rows = [[V[q][t] - P0[t] for t in range(4)] for q in (vF, vG, a, c)]
                vol += abs(_det_exact(rows)) / 24
    if return_info:
        return vol, {"nverts": nv, "nfacets": len(facets)}
    return vol


# ---------------------------------------------------------------- float cell with gradient data
def _aff_rank_float(P, tol=1e-9):
    if len(P) <= 1:
        return 0
    D = P[1:] - P[0]
    s = np.linalg.svd(D, compute_uv=False)
    return int((s > tol * max(1.0, s[0] if len(s) else 1.0)).sum())


def float_cell(A, b, tol=1e-9):
    """Return (volume, dict facet_constraint_index -> (3-volume, centroid)). Origin must be interior."""
    A = np.asarray(A, dtype=float)
    b = np.asarray(b, dtype=float)
    m = len(b)
    X = vertices_float(A, b, tol=tol)
    nv = len(X)
    if nv < 5:
        return None, None
    scale = max(1.0, float(np.abs(b).max()))
    R = X @ A.T - b[None, :]
    tight = np.abs(R) <= 1e-7 * scale
    on = [frozenset(np.nonzero(tight[:, i])[0].tolist()) for i in range(m)]
    seen = {}
    facet_data = {}
    vol = 0.0
    for i in range(m):
        S = on[i]
        if len(S) < 4 or S in seen:
            continue
        idx = sorted(S)
        if _aff_rank_float(X[idx]) != 3:
            continue
        seen[S] = i
        vF = idx[0]
        two = set()
        for j in range(m):
            G = S & on[j]
            if len(G) >= 3 and G != S and G not in two:
                if _aff_rank_float(X[sorted(G)]) == 2:
                    two.add(G)
        area = 0.0
        cen = np.zeros(4)
        for G in two:
            if vF in G:
                continue
            gl = sorted(G)
            vG = gl[0]
            edges = set()
            for k in range(m):
                E = G & on[k]
                if len(E) == 2:
                    edges.add(E)
            for E in edges:
                if vG in E:
                    continue
                a, c = sorted(E)
                D = np.array([X[vG] - X[vF], X[a] - X[vF], X[c] - X[vF]])
                g = D @ D.T
                v3 = np.sqrt(max(np.linalg.det(g), 0.0)) / 6.0
                area += v3
                cen += v3 * (X[vF] + X[vG] + X[a] + X[c]) / 4.0
        if area > 0:
            cen /= area
        h = b[i] / np.linalg.norm(A[i])
        vol += h * area / 4.0
        facet_data[i] = (area, cen)
    return vol, facet_data


def voronoi_halfspaces(Y):
    Y = np.asarray(Y, dtype=float)
    return Y, 0.5 * (Y * Y).sum(axis=1)


def voronoi_volume_and_grad(Y):
    """Volume of the Voronoi cell of 0 among neighbours Y (rows) and its gradient in Y."""
    A, b = voronoi_halfspaces(Y)
    vol, fd = float_cell(A, b)
    if vol is None:
        return None, None
    G = np.zeros_like(A)
    for i, (area, cen) in fd.items():
        y = A[i]
        G[i] = area * (y - cen) / np.linalg.norm(y)
    return vol, G


# ---------------------------------------------------------------- standard configurations
def d4_roots():
    R = []
    for i in range(4):
        for j in range(i + 1, 4):
            for si in (1, -1):
                for sj in (1, -1):
                    v = [0, 0, 0, 0]
                    v[i] = si
                    v[j] = sj
                    R.append(tuple(v))
    return R  # 24 integer roots, norm^2 = 2


if __name__ == "__main__":
    # self-test: D4 cell.  Native frame: neighbours alpha (norm^2 2), halfspace x.alpha <= 1,
    # volume 2; unit-ball frame multiplies volumes by (sqrt 2)^4 = 4.
    R = d4_roots()
    v = exact_volume([list(r) for r in R], [1] * 24)
    print("24-cell (native) exact volume:", v, " unit frame:", 4 * v)
    # deletion of one root: 25/3 in the unit frame
    v2 = exact_volume([list(r) for r in R[1:]], [1] * 23)
    print("deletion (unit frame):", 4 * v2)
    # cube [-1,1]^4 sanity
    cube = []
    for k in range(4):
        for s in (1, -1):
            e = [0] * 4
            e[k] = s
            cube.append(e)
    print("cube [-1,1]^4:", exact_volume(cube, [1] * 8), "(expect 16)")
    # cross-polytope |x|_1 <= 2: 2^4 * 2^4/24 = 32/3
    cp = [list(s) for s in itertools.product((1, -1), repeat=4)]
    print("cross-polytope |x|_1<=2:", exact_volume(cp, [2] * 16), "(expect 32/3)")
    # float check
    fv, _ = float_cell(np.array(R, float), np.ones(24))
    print("24-cell float (native):", fv)
    Y = 2 * np.array(R, float) / np.sqrt(2)  # unit frame contacts at distance 2
    vv, G = voronoi_volume_and_grad(Y)
    print("D4 Voronoi (unit frame) float:", vv, " max|grad tangential|:",
          max(abs(G[i] @ (np.eye(4) - np.outer(Y[i], Y[i]) / 4)).max() for i in range(24)))

