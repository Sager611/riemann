"""Cutting-plane shape optimization for covers of the paper's type (estimates only).

shape params (k sides): theta0 (deg, angle of n0), thetas of n_2..n_{k-1} (deg),
h ratios h_j/h_2 (j >= 3). n1 = (0,-1), h0 = h1 = 0.
objective = area(K) * max_family rho_K(Z)^2   (scale-invariant)
"""
import json
import sys
import time

import numpy as np
from scipy.optimize import minimize  # pyright: ignore[reportMissingImports]

from adversary import Shape, polyline


def unit_angle(deg):
    a = np.radians(deg)
    return np.array([np.cos(a), np.sin(a)])


class FShape(Shape):
    """Shape from float normals/offsets."""

    def __init__(self, normals, H, ngrid=1440):
        self.N = np.array(normals, dtype=float)
        self.H = np.array(H, dtype=float)
        k = len(H)
        A = np.column_stack([self.N[0], self.N[1]])
        self.a = np.array([np.linalg.solve(A, -self.N[j]) for j in range(2, k)])
        self.k = k
        self.set_grid(ngrid)


def make_shape(p, k, ngrid=720):
    th0 = p[0]
    ths = list(p[1:k - 1])
    ratios = [1.0] + list(p[k - 1:])
    normals = [unit_angle(th0), np.array([0.0, -1.0])] + [unit_angle(t) for t in ths]
    H = [0.0, 0.0] + ratios
    return FShape(normals, H, ngrid)


def feasible(p, k):
    th0 = p[0]
    ths = list(p[1:k - 1])
    if not (90 < th0 < 270):
        return False
    lo, hi = th0 - 180, 90
    if any(not (lo + 0.05 < t < hi - 0.05) for t in ths):
        return False
    if any(b <= a + 0.05 for a, b in zip(ths, ths[1:])):
        return False
    if any(r <= 0 for r in p[k - 1:]):
        return False
    return True


def poly_area(shape):
    N, H = shape.N, shape.H
    k = len(H)
    V = []
    for i in range(k):
        a, b = N[i - 1], N[i]
        M = np.array([a, b])
        V.append(np.linalg.solve(M, [H[i - 1], H[i]]))
    V = np.array(V)
    # check all sides active: each vertex satisfies other constraints
    for i in range(k):
        s = N @ V[i] - H
        s[[i - 1, i]] = 0
        if (s > 1e-12).any():
            return None
    x, y = V[:, 0], V[:, 1]
    return 0.5 * float(np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y))


def objective(p, k, family, refine=False):
    if not feasible(p, k):
        return 1e9
    S = make_shape(p, k)
    A = poly_area(S)
    if A is None or A <= 0:
        return 1e9
    r = max(S.rho(Z, refine=True) for Z in family)
    return A * r * r


def to_params(Z):
    d = np.diff(Z, axis=0)
    L = np.linalg.norm(d, axis=1)
    keep = L > 1e-9
    d, L = d[keep], L[keep]
    phi = np.arctan2(d[:, 1], d[:, 0])
    phi = phi - phi[0]
    m = len(L)
    return np.concatenate([phi[1:], np.log(L / L.sum())]), m


def reopt_arc(S, Z, iters=800):
    x0, m = to_params(Z)
    if m < 1:
        return S.rho(Z), Z
    if m == 1:
        return S.rho(Z), Z
    f = lambda x: -S.rho(polyline(x, m), refine=True)
    res = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=iters, xatol=1e-8, fatol=1e-10, adaptive=True))
    Zn = polyline(res.x, m)
    a, b = S.rho(Zn), S.rho(Z)
    return (a, Zn) if a > b else (b, Z)


def adversary(S, m, starts, rng, iters=2500):
    best = (-1, None)
    for _ in range(starts):
        x0 = np.concatenate([rng.uniform(-np.pi, np.pi, m - 1), rng.normal(0, 0.5, m)])
        f = lambda x: -S.rho(polyline(x, m), refine=False)
        res = minimize(f, x0, method="Nelder-Mead", options=dict(maxiter=iters, xatol=1e-8, fatol=1e-10, adaptive=True))
        Z = polyline(res.x, m)
        v = S.rho(Z)
        if v > best[0]:
            best = (v, Z)
    return best


def canonical():
    fam = [np.array([[0, 0], [1.0, 0]])]
    for a in (0.15, 0.2, 0.25, 0.3):
        s = 1 - 2 * a
        fam.append(np.array([[0, 0], [a, 0], [a, s], [0, s]]))       # U
        fam.append(np.array([[0, 0], [a, 0], [a, s], [2 * a, s]]))   # Z
    for ang in (60, 90, 120, 150):
        t = np.radians(ang)
        fam.append(np.array([[0.5 * np.cos(t), 0.5 * np.sin(t)], [0, 0], [0.5, 0]]))  # V
    for beta in (90, 180, 270, 360):
        b = np.radians(beta)
        R = 1 / b
        tt = np.linspace(0, b, 13)
        fam.append(np.column_stack([R * np.cos(tt), R * np.sin(tt)]))
    return fam


if __name__ == "__main__":
    k = int(sys.argv[1])
    p0 = np.array(json.loads(sys.argv[2]), dtype=float)
    seed_files = sys.argv[3].split(",") if len(sys.argv) > 3 and sys.argv[3] else []
    rounds = int(sys.argv[4]) if len(sys.argv) > 4 else 6
    out = sys.argv[5] if len(sys.argv) > 5 else "shapeopt_out.json"
    rng = np.random.default_rng(k)
    family = canonical()
    for f in seed_files:
        try:
            family.append(np.array(json.load(open(f))["Z"]))
        except Exception:
            pass
    p = p0
    t0 = time.time()
    hist = []
    for rnd in range(rounds):
        S = make_shape(p, k)
        vals = [S.rho(Z) for Z in family]
        mx = max(vals)
        active = [Z for Z, v in zip(family, vals) if v >= 0.97 * mx]
        res = minimize(objective, p, args=(k, active), method="Nelder-Mead",
                       options=dict(maxiter=500, xatol=1e-6, fatol=1e-10, adaptive=True))
        p = res.x
        S = make_shape(p, k)
        A = poly_area(S)
        vals = [S.rho(Z) for Z in family]
        mx = max(vals)
        newfam = []
        for Z, v in zip(family, vals):
            if v >= 0.97 * mx:
                v2, Zn = reopt_arc(S, Z, iters=400)
                newfam.append(Zn)
            else:
                newfam.append(Z)
        fresh = []
        for m in (3, 4, 6):
            v, Z = adversary(S, m, 3, rng, iters=1000)
            fresh.append((v, Z))
        fam_r = max(S.rho(Z) for Z in newfam)
        fr = max(v for v, _ in fresh)
        family = newfam + [Z for v, Z in fresh if v > 0.98 * max(fam_r, fr)]
        rtrue = max(fam_r, fr)
        est = A * rtrue ** 2
        hist.append(dict(round=rnd, params=p.tolist(), area_unscaled=A, rho=rtrue, est_area=est, obj=res.fun, n_family=len(family), n_active=len(active), secs=round(time.time() - t0, 1)))
        print(json.dumps(hist[-1]), flush=True)
        json.dump(dict(k=k, hist=hist, family=[Z.tolist() for Z in family]), open(out, "w"))
