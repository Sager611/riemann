"""Low-energy periodic configurations for the TEN-class dual.
Energy per unit length of a S-periodic configuration x_1..x_p (density nu = p/S), kernel k_W^2 with Fourier
transform Khat supported in [-1,1] (Poisson summation is exact and finite):
  e(x) = S^-2 * sum_{|m|<S} Khat(m/S) |sum_j exp(2 pi i m x_j / S)|^2 - p/S.
Dual: if e(x) <= nu - H_W then no TEN-type certificate for window W gives kappa > nu."""
import sys, json
import numpy as np
from scipy.optimize import minimize
from tenwin import TEN_C, MT_C, H_of, Khat

rng = np.random.default_rng(12345)


def setup(c, S):
    M = int(np.ceil(S)) - 1 if float(S).is_integer() else int(np.floor(S))
    ms = np.arange(1, M + 1)
    ms = ms[ms / S < 1]
    kh = Khat(ms / S, c)
    k0 = Khat(0.0, c)[0]
    return ms, kh, k0


def energy(x, S, ms, kh, k0):
    p = len(x)
    ph = np.exp(2j * np.pi * np.outer(ms, x) / S)
    mu = ph.sum(axis=1)
    e = (k0 * p * p + 2 * np.sum(kh * np.abs(mu) ** 2)) / S ** 2 - p / S
    g = (2 * 2 * np.real(np.conj(mu)[:, None] * (2j * np.pi * ms[:, None] / S) * ph) * kh[:, None]).sum(axis=0) / S ** 2
    return e, g


def best(c, p, S, tries=40, seeds=None):
    ms, kh, k0 = setup(c, S)
    bestv, bestx = 1e9, None
    starts = []
    # structured starts: subsets of the integer lattice (S integer) and random
    if float(S).is_integer():
        for _ in range(tries // 2):
            starts.append(np.sort(rng.choice(int(S), size=p, replace=False)).astype(float) + rng.normal(0, 0.05, p))
    while len(starts) < tries:
        starts.append(np.sort(rng.uniform(0, S, p)))
    for x0 in starts:
        r = minimize(lambda x: energy(x, S, ms, kh, k0), x0, jac=True, method="L-BFGS-B",
                     options={"maxiter": 3000, "gtol": 1e-13, "ftol": 1e-15})
        if r.fun < bestv:
            bestv, bestx = r.fun, np.sort(np.mod(r.x, S))
    return bestv, bestx


if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "TEN"
    c = TEN_C if which == "TEN" else MT_C
    H = H_of(c)
    out = []
    for (p, S) in [(2, 3), (4, 6), (6, 9), (8, 12), (10, 15), (12, 18), (14, 21), (16, 24), (20, 30), (24, 36),
                   (27, 40), (29, 43), (31, 46), (33, 49), (35, 52), (37, 55), (39, 58)]:
        nu = p / S
        v, x = best(c, p, S)
        out.append((p, S, nu, v, nu - H, v <= nu - H))
        print(f"p={p:3d} S={S:3d} nu={nu:.6f} e_min={v:.6e} nu-H={nu-H:.6e} dual_ok={v <= nu - H}", flush=True)
