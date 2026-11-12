"""Enumerate S-periodic configurations seeded by subsets of Z_S (necklaces), relax positions, record (nu, e).
Lower convex envelope of (nu, e) = mixture energies (block interfaces cost O(log B), negligible per length)."""
import sys, itertools, json
import numpy as np
from scipy.optimize import minimize
from tenwin import TEN_C, MT_C, H_of
from periodic import setup, energy


def necklaces(S):
    seen = set()
    for mask in range(1, 1 << S):
        bits = tuple((mask >> i) & 1 for i in range(S))
        rots = [bits[i:] + bits[:i] for i in range(S)]
        r = min(rots)
        if r in seen:
            continue
        seen.add(r)
        yield r


def relax(c, S, x0, cache={}):
    key = (id(c), S)
    if key not in cache:
        cache[key] = setup(c, S)
    ms, kh, k0 = cache[key]
    r = minimize(lambda x: energy(x, S, ms, kh, k0), x0, jac=True, method="L-BFGS-B",
                 options={"maxiter": 2000, "gtol": 1e-14, "ftol": 1e-16})
    e0 = energy(x0, S, ms, kh, k0)[0]
    return min(r.fun, e0), (np.sort(np.mod(r.x, S)) if r.fun < e0 else x0)


def envelope(pts):
    pts = sorted(pts)
    hull = []
    for p in pts:
        while len(hull) >= 2 and (hull[-1][0] - hull[-2][0]) * (p[1] - hull[-2][1]) - (hull[-1][1] - hull[-2][1]) * (p[0] - hull[-2][0]) <= 0:
            hull.pop()
        hull.append(p)
    return hull


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "TEN"
    Smax = int(sys.argv[2]) if len(sys.argv) > 2 else 12
    c = TEN_C if which == "TEN" else (MT_C if which == "MT" else np.array(json.load(open(which))["c"]))
    H = H_of(c)
    pts = []
    for S in range(1, Smax + 1):
        for nk in necklaces(S):
            p = sum(nk)
            nu = p / S
            if nu < 0.55 or nu > 1.0:
                continue
            x0 = np.array([i for i in range(S) if nk[i]], dtype=float)
            e, x = relax(c, S, x0)
            pts.append((nu, e, S, "".join("x" if b else "." for b in nk), [round(v, 6) for v in x]))
    pts.sort()
    hull = envelope([(a, b) for a, b, *_ in pts])
    print("H_W =", H)
    print("hull:")
    for (nu, e) in hull:
        lab = [q for q in pts if abs(q[0] - nu) < 1e-12 and abs(q[1] - e) < 1e-15][0]
        print(f"  nu={nu:.6f} e={e:.6e} S={lab[2]} pattern={lab[3]} x={lab[4]}")
    # root of nu = H + ebar(nu) along the hull
    for (a, b), (a2, b2) in zip(hull, hull[1:]):
        fa, fb = a - H - b, a2 - H - b2
        if fa < 0 <= fb:
            t = -fa / (fb - fa)
            print("root nu* =", a + t * (a2 - a), "between", a, a2)
    json.dump(pts, open(f"necklace_{which.replace('.json', '')}_{Smax}.json", "w"))
