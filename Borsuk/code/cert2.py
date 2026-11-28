"""Certifier v2 (Borsuk P8): diam(K) = max over unit u of the width h_K(u) + h_K(-u) for convex compact K.
K = ∩ balls B(c_i, sqrt(r2_i)) ∩ halfspaces {a_j.x <= b_j}.  Upper bound of h_K(u) by Lagrange duality: for any
lam >= 0 (sum > 0), mu >= 0:  h_K(u) <= |u + 2 Σ lam_i c_i - Σ mu_j a_j|^2 / (4 Σ lam_i) - Σ lam_i (|c_i|^2 - r2_i) + Σ mu_j b_j.
Multipliers from a float optimisation, bound evaluated in exact rationals (u rational, not necessarily unit: then the
bound is for h_K(u) with |u| = |u|, handled by dividing by a rational lower bound of |u|).  Covering of S^3: the
boundary of the cube [-1,1]^4 (8 facets) subdivided adaptively; a sub-cube with centre v0 (|v0| >= 1) and half side s
covers directions within chord delta = s*sqrt(3) of u0 = v0/|v0| (|v| >= 1 on the cube surface), and
|h(u) - h(u0)| <= Rmax * |u - u0| with Rmax = max |x| over K (<= circumradius).  Also max distance between two convex
sets K1, K2: max_u h_K1(u) + h_K2(-u), same covering.  All final comparisons exact."""
import numpy as np
from fractions import Fraction as F
from math import isqrt
from scipy.optimize import minimize

def sqrt_up(q):
    s = isqrt(int(q * (1 << 62)) + 1) + 1
    return F(s, 1 << 31)
def sqrt_down(q):
    s = isqrt(int(q * (1 << 62)))
    return F(max(s - 1, 0), 1 << 31)

class Body:
    def __init__(s, balls, halfs):
        s.balls = [(tuple(F(x) for x in c), F(r2)) for c, r2 in balls]
        s.halfs = [(tuple(F(x) for x in a), F(b)) for a, b in halfs]
        s.C = np.array([[float(x) for x in c] for c, r2 in s.balls]); s.R2 = np.array([float(r2) for c, r2 in s.balls])
        s.A = np.array([[float(x) for x in a] for a, b in s.halfs]).reshape(-1, 4); s.B = np.array([float(b) for a, b in s.halfs])
        s.warm = {}
    def dual_float(s, u, x0=None):
        nb, nh = len(s.balls), len(s.halfs)
        def g(z):
            lam = np.exp(z[:nb]); mu = np.exp(z[nb:]) if nh else np.zeros(0)
            w = u + 2 * lam @ s.C - (mu @ s.A if nh else 0)
            return w @ w / (4 * lam.sum()) - lam @ ((s.C ** 2).sum(1) - s.R2) + (mu @ s.B if nh else 0)
        z0 = x0 if x0 is not None else np.r_[np.zeros(nb), -3 * np.ones(nh)]
        r = minimize(g, z0, method='L-BFGS-B', bounds=[(-30, 10)] * (nb + nh))
        return r.x, r.fun
    def h_ub(s, ur, key=None):
        """rigorous upper bound of h_K(u) for rational vector ur (any length)"""
        u = np.array([float(x) for x in ur])
        z, _ = s.dual_float(u, s.warm.get(key))
        if key is not None: s.warm[key] = z
        nb = len(s.balls)
        lam = [F(np.exp(t)).limit_denominator(1 << 40) for t in z[:nb]]; mu = [F(np.exp(t)).limit_denominator(1 << 40) for t in z[nb:]]
        lam = [max(l, F(1, 1 << 40)) for l in lam]
        w = [ur[k] + 2 * sum(l * c[k] for l, (c, r2) in zip(lam, s.balls)) - sum(m * a[k] for m, (a, b) in zip(mu, s.halfs)) for k in range(4)]
        L = sum(lam)
        return sum(t * t for t in w) / (4 * L) - sum(l * (sum(x * x for x in c) - r2) for l, (c, r2) in zip(lam, s.balls)) + sum(m * b for m, (a, b) in zip(mu, s.halfs))

def facets_init(n=2):
    out = []
    for ax in range(4):
        for sg in (1, -1):
            for idx in np.ndindex(*([n] * 3)):
                c = [F(-1) + F(2 * i + 1, n) for i in idx]; v = c[:ax] + [F(sg)] + c[ax:]
                out.append((tuple(v), F(1, n), ax))
    return out

def max_pair_width(K1, K2, Rmax, target, max_boxes=200000):
    """certify max_{x in K1, y in K2} |x - y| < target (target rational); returns (ok, worst bound, boxes)"""
    todo = facets_init(2); worst = F(0); n = 0
    while todo:
        v, s, ax = todo.pop(); n += 1
        if n > max_boxes: return False, None, n
        nv2 = sum(x * x for x in v); nlo = sqrt_down(nv2)          # |v| >= nlo >= 1
        hu = K1.h_ub(v, key=('a', v)) + K2.h_ub(tuple(-x for x in v), key=('b', v))   # bound for |v| * width(u0)
        wid = hu / nlo                                              # >= width in direction u0 = v/|v|
        delta = s * sqrt_up(F(3))                                   # chord radius of the covered directions
        bound = wid + 2 * Rmax * delta
        if bound < target: worst = max(worst, bound); continue
        if s < F(1, 1 << 12): return False, bound, n
        h = s / 2
        for d in np.ndindex(2, 2, 2):
            off = [h if t else -h for t in d]; c = list(v)
            j = 0
            for k in range(4):
                if k != ax: c[k] = v[k] + off[j]; j += 1
            todo.append((tuple(c), h, ax))
    return True, worst, n
