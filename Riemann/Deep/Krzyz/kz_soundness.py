#!/usr/bin/env python3
"""Soundness tests for kz_check.py: find genuine solutions of the reduced moment system
(KKT points of n-1 atom configurations) numerically, check |a_n| = e^{-t} C (Prop. central-coeff)
against a direct Taylor-coefficient computation, and check that the checker never excludes
a box containing a solution when the coefficient threshold is lowered below its |a_n|."""
import sys, math, json, random
import numpy as np
from scipy.optimize import least_squares
import kz_check as K
from flint import arb

def F(x, n):
    N = n - 1
    th, q, lam = x[:N], x[N], x[N + 1:]
    out = []
    for k in range(1, n + 1):
        Tk = math.cosh(k * math.acosh(max(q, 1.0)))
        z = sum((lam[j] - 1.0 / k) * complex(math.cos(k * th[j]), math.sin(k * th[j])) for j in range(N)) - Tk / k
        out += [z.real, z.imag]
    return np.array(out)

def an_direct(th, lam, n):
    # Taylor coefficient of exp(-sum lam_j (1+a_j z)/(1-a_j z)) via log-series recursion
    N = len(th)
    h = [0j] * (n + 1)
    for m in range(1, n + 1):
        h[m] = -2 * sum(lam[j] * complex(math.cos(m * th[j]), math.sin(m * th[j])) for j in range(N))
    f = [0j] * (n + 1); f[0] = 1
    for k in range(1, n + 1):
        f[k] = sum(j * h[j] * f[k - j] for j in range(1, k + 1)) / k
    return abs(f[n]) * math.exp(-sum(lam))

def C_float(th, q):
    p = np.array([1.0 + 0j])
    for t in th:
        p = np.convolve(p, np.array([-complex(math.cos(t), math.sin(t)), 1.0]))
    s0 = float(np.sum(np.abs(p) ** 2)); s1 = sum((p[k + 1] * np.conj(p[k])).real for k in range(len(p) - 1))
    return 2 * q * s0 - 2 * s1

def find(n, tries, seed):
    rng = random.Random(seed); N = n - 1; sols = []
    for _ in range(tries):
        x0 = [rng.uniform(0, 2 * math.pi) for _ in range(N)] + [rng.uniform(1, 1.25)] + [rng.uniform(0.02, 1.5) for _ in range(N)]
        lb = [-10] * N + [1.0] + [0.0] * N; ub = [20] * N + [1.4] + [K.THI[n]] * N
        r = least_squares(F, x0, args=(n,), bounds=(lb, ub), xtol=1e-15, ftol=1e-15, gtol=1e-15, max_nfev=4000)
        if np.max(np.abs(r.fun)) < 1e-11 and min(r.x[N + 1:]) > 1e-6:
            th = [t % (2 * math.pi) for t in r.x[:N]]; lam = list(r.x[N + 1:])
            order = sorted(range(N), key=lambda j: th[j]); th = [th[j] for j in order]; lam = [lam[j] for j in order]
            if th[0] + th[-1] > 2 * math.pi:
                th2 = [(2 * math.pi - t) % (2 * math.pi) for t in th]; order = sorted(range(N), key=lambda j: th2[j])
                th = [th2[j] for j in order]; lam = [lam[j] for j in order]
            q = float(r.x[N]); t = sum(lam)
            a1 = an_direct(th, lam, n); a2 = math.exp(-t) * C_float(th, q)
            sols.append(dict(theta=th, q=q, lam=lam, t=t, an_direct=a1, an_formula=a2))
    # dedupe
    out = []
    for s in sorted(sols, key=lambda s: -s['an_direct']):
        if all(abs(s['an_direct'] - o['an_direct']) > 1e-9 or abs(s['t'] - o['t']) > 1e-7 for o in out):
            out.append(s)
    return out

def box_around(s, rad):
    lo = [t - rad for t in s['theta']] + [max(1.0, s['q'] - rad)]
    hi = [t + rad for t in s['theta']] + [s['q'] + rad]
    return lo, hi

if __name__ == '__main__':
    n = int(sys.argv[1]); tries = int(sys.argv[2]) if len(sys.argv) > 2 else 60
    sols = find(n, tries, 12345 + n)
    report = []
    for s in sols[:4]:
        lo, hi = box_around(s, 1e-3)
        ok_true, nodes1, reasons1, _ = K.check_cell(n, lo, hi, limit=400)
        saved = K.TWO_OVER_E
        K.TWO_OVER_E = arb(s['an_direct'] * 0.9)          # threshold below the true |a_n|: must NOT exclude
        ok_low, nodes2, reasons2, _ = K.check_cell(n, lo, hi, limit=150)
        K.TWO_OVER_E = saved
        report.append(dict(an=s['an_direct'], an_formula_relerr=abs(s['an_formula'] - s['an_direct']) / s['an_direct'],
                           t=s['t'], q=s['q'], excluded_with_true_threshold=ok_true, reasons_true=reasons1,
                           excluded_with_lowered_threshold=ok_low, reasons_low=reasons2, nodes_low=nodes2))
    print(json.dumps(dict(n=n, solutions_found=len(sols), top=[round(s['an_direct'], 9) for s in sols[:8]], tests=report), indent=1))
