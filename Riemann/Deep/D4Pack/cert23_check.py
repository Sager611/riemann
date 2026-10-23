#!/usr/bin/env python3
"""cert23_check.py -- independent sample check of the m=23 certificate (thm:certificate,
TeX 5081-5097) from the supplementary file certificate_d8.txt (GitHub creelie/d4-voronoi-cells,
multi_cap/continuation_out/).  Own implementation of omega, G_k, Y_k, S_k from the TeX.
Checks: (1) F_k positive definite, f_k >= 0 (exact LDL over Q); (2) B exact vs 8 - A_*;
(3) condition (C) 1000(w(u)+w(v)+w(t)) - P(u,v,t) >= 0 at many admissible triples (float).
Not a proof of (C); a search for a violating triple."""
import sys, itertools
from fractions import Fraction as Fr
import numpy as np
from numpy.polynomial import legendre as L

path = sys.argv[1]
lines = open(path).read().split("\n")
f = None; Fm = []; cur = None
for ln in lines:
    if ln.startswith("f = "):
        f = [float(x) for x in ln[4:].split()]
    elif ln.startswith("F") and "(" in ln:
        cur = []; Fm.append(cur)
    elif cur is not None and ln.strip() and ln.startswith("  "):
        cur.append([float(x) for x in ln.split()])
d = len(f) - 1
assert d == 8 and len(Fm) == 9 and all(len(Fm[k]) == d - k + 1 for k in range(9))

def ldl_pd(M):
    A = [[Fr(x) for x in r] for r in M]; n = len(A)
    for c in range(n):
        if A[c][c] <= 0: return False
        for r in range(c + 1, n):
            q = A[r][c] / A[c][c]
            for k in range(c, n): A[r][k] -= q * A[c][k]
    return True
sym = all(Fm[k][i][j] == Fm[k][j][i] for k in range(9) for i in range(len(Fm[k])) for j in range(len(Fm[k])))
pd = all(ldl_pd(F) for F in Fm)
print(f"(1) F_k symmetric: {sym}; positive definite (exact LDL over Q): {pd}; f_k>=0 for k>=1: {all(x >= 0 for x in f[1:])}")

N = 23
fq = [Fr(x) for x in f]
f1 = sum(fq)                       # G_k(1) = 1
F111 = sum(Fr(x) for r in Fm[0] for x in r)
B = (Fr(N, 2) * (N * fq[0] - f1) - Fr(N) * F111 / (6 * (N - 2))) / 1000
import mpmath as mp
mp.mp.dps = 30
rs = mp.atan(1 / mp.sqrt(2))
Astar = 9 * mp.pi**2 / 8 - 207 * mp.pi / 8 * rs + 253 * mp.pi / (12 * mp.sqrt(2))
print(f"(2) B = {mp.nstr(mp.mpf(B.numerator)/B.denominator, 13)} (paper 0.0929000002); 8 - A_* = {mp.nstr(8-Astar, 13)}; B > 8-A_*: {mp.mpf(B.numerator)/B.denominator > 8 - Astar}")

# omega(u) on [1/3, 1/2] (zero below), omega(gamma) = int_{gamma/2}^{r_*} Lam(r,gamma) tan r sec^4 r dr
r_star = float(rs)
def Lam(r, g):
    return 4*np.pi*((r - g/2)/2 - (np.sin(2*r) - np.sin(g))/4 - np.tan(g/2)/2*(np.sin(r)**2 - np.sin(g/2)**2))
xg, wg = np.polynomial.legendre.leggauss(40)
ug = np.linspace(1/3, 0.5, 4001)
gg = np.arccos(ug)
wtab = np.zeros_like(ug)
for i, g in enumerate(gg):
    a, b = g / 2, r_star
    if b <= a: continue
    r = (b - a)/2*xg + (a + b)/2
    wtab[i] = ((b - a)/2*wg*Lam(r, g)*np.tan(r)/np.cos(r)**4).sum()
def omega(u):
    u = np.asarray(u, float)
    return np.where(u > 1/3, np.interp(u, ug, wtab), 0.0)
print(f"    omega(1/2) = {omega(0.5):.10f} (paper 0.00144540...)")

def cheb(x, n):          # T_0..T_n at x  -> (len(x), n+1)
    T = np.zeros((len(x), n + 1)); T[:, 0] = 1
    if n >= 1: T[:, 1] = x
    for i in range(2, n + 1): T[:, i] = 2*x*T[:, i-1] - T[:, i-2]
    return T
def gk(x, n):            # G_k = U_k/(k+1)
    U = np.zeros((len(x), n + 1)); U[:, 0] = 1
    if n >= 1: U[:, 1] = 2*x
    for i in range(2, n + 1): U[:, i] = 2*x*U[:, i-1] - U[:, i-2]
    return U / np.arange(1, n + 2)
legc = [L.leg2poly([0]*k + [1]) for k in range(d + 1)]   # power-basis coefficients of P_k
def Wk(k, a, b, t):      # ((1-a^2)(1-b^2))^{k/2} P_k((t-ab)/sqrt(..)) as a polynomial
    c = legc[k]; s = t - a*b; AB = (1 - a*a)*(1 - b*b); out = np.zeros_like(a)
    for p in range(k + 1):
        if c[p] != 0:     # p = k - 2m
            m = (k - p)//2; out += c[p] * s**p * AB**m
    return out
Fa = [np.array(F) for F in Fm]
def Fpoly(u, v, t):
    Tc = {}
    tot = np.zeros_like(u)
    for (a, b, c) in itertools.permutations((0, 1, 2)):
        X = (u, v, t)
        A, Bv, Cv = X[a], X[b], X[c]
        TA = cheb(A, d); TB = cheb(Bv, d)
        for k in range(d + 1):
            n = d - k + 1
            q = ((TA[:, :n] @ Fa[k]) * TB[:, :n]).sum(1)
            tot += q * Wk(k, A, Bv, Cv) / 6
    return tot
fv = np.array(f)
def fpoly(x):
    return gk(x, d) @ fv
def Q(u, v, t):
    one = np.ones_like(u)
    P = fpoly(u) + fpoly(v) + fpoly(t) + Fpoly(u, v, t) + (Fpoly(one, u, u) + Fpoly(one, v, v) + Fpoly(one, t, t)) / (N - 2)
    return 1000*(omega(u) + omega(v) + omega(t)) - P

rng = np.random.default_rng(7)
def admissible(u, v, t):
    return (u <= .5) & (v <= .5) & (t <= .5) & (u >= -1) & (v >= -1) & (t >= -1) & (1 + 2*u*v*t - u*u - v*v - t*t >= -1e-15)
worst = (np.inf, None); tested = 0
for batch in range(12):
    n = 40000
    u = rng.uniform(-1, .5, n); v = rng.uniform(-1, .5, n)
    s = np.sqrt(np.clip((1 - u*u)*(1 - v*v), 0, None))
    if batch % 3 == 0:     # interior
        t = u*v + s*rng.uniform(-1, 1, n)
    elif batch % 3 == 1:   # coplanar boundary
        t = u*v + s*rng.choice([-1., 1.], n)
    else:                  # near 60 degrees in one or more pairs
        u = .5 - rng.exponential(0.02, n); v = np.where(rng.uniform(size=n) < .5, .5 - rng.exponential(.02, n), v)
        s = np.sqrt(np.clip((1 - u*u)*(1 - v*v), 0, None)); t = u*v + s*rng.uniform(-1, 1, n)
    ok = admissible(u, v, t); u, v, t = u[ok], v[ok], t[ok]
    q = Q(u, v, t); tested += len(q)
    i = int(np.argmin(q))
    if q[i] < worst[0]: worst = (q[i], (u[i], v[i], t[i]))
vals = [-1, -.5, 0, .5]
dl = [(a, b, c) for a in vals for b in vals for c in vals if admissible(np.array([a]), np.array([b]), np.array([c]))[0]]
qd = Q(*[np.array(z, float) for z in zip(*dl)])
q16 = Q(np.array([-1/16]), np.array([-1/16]), np.array([-1/16]))[0]
print(f"(3) {tested} admissible triples: min Q = {worst[0]:.3e} at {tuple(round(x, 4) for x in worst[1])};"
      f" {len(dl)} root-system triples: min Q = {qd.min():.3e}; Q(-1/16,-1/16,-1/16) = {q16:.3e}")
print("    no violation found" if min(worst[0], qd.min(), q16) >= 0 else "    VIOLATION FOUND")

