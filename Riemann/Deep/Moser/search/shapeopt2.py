"""Trust-region cutting-plane shape optimization (estimates of true optimal area)."""
import json
import sys
import time

import numpy as np
from scipy.optimize import minimize  # pyright: ignore[reportMissingImports]

from shapeopt import make_shape, poly_area, feasible, reopt_arc, adversary, canonical
from adversary import polyline


def obj_family(p, k, fam):
    if not feasible(p, k):
        return 1e9
    S = make_shape(p, k)
    A = poly_area(S)
    if A is None or A <= 0:
        return 1e9
    return A * max(S.rho(Z) for Z in fam) ** 2


def true_estimate(p, k, fam, rng, fresh_starts=3):
    S = make_shape(p, k)
    A = poly_area(S)
    vals = [S.rho(Z) for Z in fam]
    mx = max(vals)
    newfam = []
    for Z, v in zip(fam, vals):
        if v >= 0.985 * mx:
            _, Zn = reopt_arc(S, Z, iters=500)
            newfam.append(Zn)
    for m in (3, 4, 5, 6, 8):
        v, Z = adversary(S, m, fresh_starts, rng, iters=1200)
        newfam.append(Z)
    rr = [S.rho(Z) for Z in newfam]
    return A * max(rr) ** 2, newfam, max(rr)


if __name__ == "__main__":
    k = int(sys.argv[1])
    p = np.array(json.loads(sys.argv[2]), dtype=float)
    seeds = [f for f in sys.argv[3].split(",") if f]
    rounds = int(sys.argv[4])
    out = sys.argv[5]
    delta = float(sys.argv[6]) if len(sys.argv) > 6 else 1.0
    rng = np.random.default_rng(1000 + k)
    fam = canonical()
    for f in seeds:
        try:
            d = json.load(open(f))
            if "Z" in d:
                fam.append(np.array(d["Z"]))
            for Z in d.get("family", []):
                fam.append(np.array(Z))
        except Exception:
            pass
    t0 = time.time()
    est, newarcs, r = true_estimate(p, k, fam, rng)
    fam += newarcs
    best = (est, p.copy())
    hist = [dict(round=-1, params=p.tolist(), est=est, rho=r, delta=delta)]
    print(json.dumps(hist[-1]), flush=True)
    scale = np.array([1.0] * (k - 1) + [0.02] * (k - 3))  # degrees for angles, ratio units
    for rnd in range(rounds):
        pb = best[1]
        S = make_shape(pb, k)
        vals = [S.rho(Z) for Z in fam]
        mx = max(vals)
        act = [Z for Z, v in zip(fam, vals) if v >= 0.95 * mx]

        def f(x):
            q = pb + x * scale
            if np.max(np.abs(x)) > delta:
                return 1e9
            return obj_family(q, k, act)

        sim = [np.zeros(len(pb))]
        for i in range(len(pb)):
            e = np.zeros(len(pb))
            e[i] = delta / 2
            sim.append(e)
        res = minimize(f, np.zeros(len(pb)), method="Nelder-Mead",
                       options=dict(maxiter=400, initial_simplex=np.array(sim), xatol=1e-5, fatol=1e-10, adaptive=True))
        q = pb + res.x * scale
        est, newarcs, r = true_estimate(q, k, act, rng)
        fam += newarcs
        if est < best[0]:
            best = (est, q.copy())
            delta = min(delta * 1.5, 8.0)
            acc = True
        else:
            delta = delta * 0.5
            acc = False
        # prune family: keep arcs relevant at best shape
        S = make_shape(best[1], k)
        vals = [S.rho(Z) for Z in fam]
        mx = max(vals)
        fam = [Z for Z, v in zip(fam, vals) if v >= 0.93 * mx]
        hist.append(dict(round=rnd, params=q.tolist(), est=est, model=res.fun, rho=r, accepted=acc,
                         best_est=best[0], best=best[1].tolist(), delta=delta, nfam=len(fam), secs=round(time.time() - t0)))
        print(json.dumps(hist[-1]), flush=True)
        json.dump(dict(k=k, hist=hist, best=best[1].tolist(), best_est=best[0], family=[Z.tolist() for Z in fam]), open(out, "w"))
        if delta < 0.01:
            break

