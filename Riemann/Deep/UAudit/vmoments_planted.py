#!/usr/bin/env python3
"""uaudit N3b: Step 5c on the interval J = [T, T+400]: actual zeros, and with planted off-line pairs.
A fraction p of the zeros in [T/2,5T/2] is replaced by off-line pairs (1/2 +- d + i gamma).
Re V = pi^-1 sum_on log|(t+h/2-g)/(t-h/2-g)| + pi^-1 sum_pairs log(((t+h/2-g)^2+d^2)/((t-h/2-g)^2+d^2))
(exact Hadamard terms of the zeros in [T/2,5T/2]; zeros farther than 60 from J enter through a smooth
interpolated far part; the O(1) remainder of zeros outside [T/2,5T/2] and Gamma is dropped, see vmoments.py (1)).
Im V: sharp count, pair counted twice. Im Vhat = F~: pair Poisson-smeared (PROOF 6b).
"""
import sys, json, time, math
import numpy as np
from scipy.integrate import quad

ZF = sys.argv[2] if len(sys.argv) > 2 else "../ExploreHelio/results/zeta_zeros.npy"
Z = np.sort(np.load(ZF))
T = 2000.0
L = math.log(T)
H = 2 * math.pi * 0.7 / L
ZW = Z[(Z >= T / 2) & (Z <= 2.5 * T)]
A, B = T, T + 400.0
W = 60.0
rng = np.random.default_rng(11)
TINY = 1e-300

def mom(fun, sing, m):
    pts = np.unique(np.concatenate([[A, B], sing[(sing > A) & (sing < B)]]))
    tot = 0.0
    for i in range(len(pts) - 1):
        if pts[i + 1] - pts[i] > 1e-12:
            tot += quad(lambda t: abs(fun(t)) ** m, pts[i], pts[i + 1], limit=50, epsrel=1e-5, epsabs=1e-10)[0]
    return tot

def run(p, d, ms=(2, 4, 6)):
    mask = rng.random(len(ZW)) < p
    g_on, g_off = ZW[~mask], ZW[mask]
    near = lambda g: g[(g > A - W) & (g < B + W)]
    far = lambda g: g[(g <= A - W) | (g >= B + W)]
    n_on, n_off, f_on, f_off = near(g_on), near(g_off), far(g_on), far(g_off)
    grid = np.linspace(A, B, 801)
    farv = np.array([np.sum(np.log(np.abs((t + H / 2 - f_on) / (t - H / 2 - f_on)))) +
                     np.sum(np.log(((t + H / 2 - f_off) ** 2 + d * d) / ((t - H / 2 - f_off) ** 2 + d * d))) for t in grid])
    def ReV(t):
        v = np.sum(np.log(np.abs(t + H / 2 - n_on) + TINY) - np.log(np.abs(t - H / 2 - n_on) + TINY))
        v += np.sum(np.log(((t + H / 2 - n_off) ** 2 + d * d) / ((t - H / 2 - n_off) ** 2 + d * d)))
        return (float(v) + float(np.interp(t, grid, farv))) / math.pi
    def ImVhat(t):
        v = float(np.sum(np.abs(t - n_on) < H / 2))
        v += (2 / math.pi) * float(np.sum(np.arctan((t + H / 2 - g_off) / d) - np.arctan((t - H / 2 - g_off) / d)))
        return v
    def ImV(t):
        return float(np.sum(np.abs(t - n_on) < H / 2)) + 2.0 * float(np.sum(np.abs(t - n_off) < H / 2))
    sing = np.concatenate([n_on - H / 2, n_on + H / 2])
    sing2 = np.concatenate([sing, n_off - H / 2, n_off + H / 2])
    res = {}
    for m in ms:
        res[m] = dict(ReV=mom(ReV, sing, m) / (B - A), ImVhat=mom(ImVhat, sing, m) / (B - A), ImV=mom(ImV, sing2, m) / (B - A))
    return int(mask.sum()), res

if __name__ == "__main__":
    t0 = time.time()
    out = dict(T=T, L=L, h=H, J=[A, B], rows=[])
    for (p, d) in ((0.0, 1.0), (0.05, 1 / L), (0.05, 3 / L), (0.05, 0.45), (0.3, 1 / L), (0.3, 0.45)):
        n, res = run(p, d)
        out["rows"].append(dict(p=p, d=d, nplanted=n, moments=res))
        print("p=%.2f d=%.3f n=%d " % (p, d, n) + "  ".join("m=%d: |ReV|^m=%.4f ImVhat^m=%.4f ImV^m=%.4f" % (
            m, r["ReV"], r["ImVhat"], r["ImV"]) for m, r in res.items()) + "  (%.0fs)" % (time.time() - t0), flush=True)
    json.dump(out, open(sys.argv[1] if len(sys.argv) > 1 else "vmoments_planted.json", "w"), indent=1)
    print("time %.1f s" % (time.time() - t0))
