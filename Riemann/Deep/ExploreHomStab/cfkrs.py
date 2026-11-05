# explore9: CFKRS recipe polynomial for the family chi_{8d}, d odd squarefree (Soundararajan's family).
#   sum*_{d<=X} L(1/2, chi_8d)^k  ~  sum*_{d<=X} Q_k(log(8d/pi))
# Q_k(L) = ((-1)^{k(k-1)/2}/k!) [prod z_j^{2k-1}] { exp((L/2) sum z) F(z) V(z) },
#   F = A(z) prod_{i<=j} zt(z_i+z_j) prod_j g(z_j)^{-1/2},  zt(s) = s zeta(1+s),
#   g(z) = Gamma((1/2-z)/2)/Gamma((1/2+z)/2),  V = prod_{i<j} (z_i-z_j)^2 (z_i+z_j),
#   A = prod_{i<=j}(1-2^{-1-z_i-z_j}) prod_{p odd} [1 + p/(p+1)(H_e(x)-1)] prod_{i<=j}(1-x_i x_j), x_i = p^{-1/2-z_i}.
# Taylor coefficients of F by FFT on a polycircle (radii r_j); primes p > P0 via the p^{-2} term of log A_p and
# prime sums tau_m = sum_{p>P0} (log p)^m p^{-2} (to P1, then PNT integral).
# Usage: uv run --no-project --with numpy --with scipy --with mpmath python3 cfkrs.py k [N] [r] [P0]
import sys, json
from math import factorial
from itertools import product
import numpy as np, mpmath as mp
from scipy.special import gammaincc, exp1, gammaln

k = int(sys.argv[1]); N = int(sys.argv[2]) if len(sys.argv) > 2 else 32
r = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
P0 = int(sys.argv[4]) if len(sys.argv) > 4 else 20000
P1 = 10**7
mp.mp.dps = 30

def sieve(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0]
pr = sieve(P1)
small = pr[(pr > 2) & (pr <= P0)].astype(float)
big = pr[pr > P0].astype(float)

# ---- polynomials as dicts {exponent tuple: coef}
def padd(a, b, c=1.0):
    out = dict(a)
    for e, v in b.items(): out[e] = out.get(e, 0) + c * v
    return out
def pmul(a, b):
    out = {}
    for e1, v1 in a.items():
        for e2, v2 in b.items():
            e = tuple(x + y for x, y in zip(e1, e2)); out[e] = out.get(e, 0) + v1 * v2
    return out
def unit(i, k):
    return tuple(1 if j == i else 0 for j in range(k))
def hpoly(m, k):
    return {e: 1.0 for e in product(range(m + 1), repeat=k) if sum(e) == m}
def Spoly(m, k):
    out = {}
    for i in range(k):
        for j in range(i, k):
            e = [0] * k; e[i] += m; e[j] += m; e = tuple(e); out[e] = out.get(e, 0) + 1.0
    return out

# ---- tail: sum_{p>P0} p^{-2} c2(p^{-z}),  c2 = h4 - h2 - h2^2/2 - S2/2
h2, h4 = hpoly(2, k), hpoly(4, k)
c2 = padd(padd(padd(h4, h2, -1.0), pmul(h2, h2), -0.5), Spoly(2, k), -0.5)
M = 90
lb = np.log(big)
tau = np.zeros(M)   # tau[m] = sum_{p>P0} p^{-2} (log p)^m / m!
for m in range(M):
    tau[m] = np.sum(np.exp(m * np.log(lb) - gammaln(m + 1) - 2 * lb))
    # PNT tail beyond P1: int_{log P1}^inf v^{m-1} e^{-v} dv / m!  (s = 2)
    x0 = np.log(P1)
    tail = exp1(x0) if m == 0 else gammaincc(m, x0) / m
    tau[m] += tail
def T2(w):
    out = np.zeros_like(w); pw = np.ones_like(w)
    for m in range(M):
        out += tau[m] * pw; pw = pw * (-w)
    return out

# ---- grid
radii = [r * (1.0 - 0.07 * j) for j in range(k)]
om = np.exp(2j * np.pi * np.arange(N) / N)
Z = np.meshgrid(*[rad * om for rad in radii], indexing='ij')
logF = np.zeros(Z[0].shape, dtype=complex)
# p = 2 factor of A
for i in range(k):
    for j in range(i, k):
        logF += np.log(1 - 2.0 ** (-1 - Z[i] - Z[j]))
# odd primes <= P0, chunked
for p in small:
    lp = np.log(p)
    x = [np.exp(-(0.5 + Zi) * lp) for Zi in Z]
    Pm = np.ones_like(Z[0]); Pp = np.ones_like(Z[0])
    for xi in x:
        Pm = Pm / (1 - xi); Pp = Pp / (1 + xi)
    He = 0.5 * (Pm + Pp)
    R = 1 + (p / (p + 1)) * (He - 1)
    C = np.ones_like(Z[0])
    for i in range(k):
        for j in range(i, k):
            C = C * (1 - x[i] * x[j])
    logF += np.log(R) + np.log(C)
# tail primes
for e, v in c2.items():
    w = sum(e[j] * Z[j] for j in range(k))
    logF += v * T2(w)
# zt(s) = s zeta(1+s) = 1 + sum_{n>=0} (-1)^n gamma_n s^{n+1}/n!
gam = [float(mp.stieltjes(n)) for n in range(40)]
def zt(s):
    out = np.ones_like(s); pw = s.copy()
    for n in range(40):
        out += ((-1) ** n) * gam[n] / factorial(n) * pw; pw = pw * s
    return out
for i in range(k):
    for j in range(i, k):
        logF += np.log(zt(Z[i] + Z[j]))
# g(z)^{-1/2}
def logg(z):
    zz = z.ravel()
    v = np.array([complex(mp.loggamma((0.5 - complex(t)) / 2) - mp.loggamma((0.5 + complex(t)) / 2)) for t in zz])
    return v.reshape(z.shape)
for i in range(k):
    ax = radii[i] * om
    lg = logg(ax)
    shape = [1] * k; shape[i] = N
    logF += -0.5 * lg.reshape(shape)
F = np.exp(logF)
f = np.fft.fftn(F) / N ** k
# divide by prod r_j^{b_j}
for i in range(k):
    shape = [1] * k; shape[i] = N
    f = f / (radii[i] ** np.arange(N)).reshape(shape)

# ---- assemble Q_k
V = {tuple([0] * k): 1.0}
for i in range(k):
    for j in range(i + 1, k):
        dif = padd({unit(i, k): 1.0}, {unit(j, k): 1.0}, -1.0)
        sm = padd({unit(i, k): 1.0}, {unit(j, k): 1.0})
        V = pmul(V, pmul(pmul(dif, dif), sm))
S1 = {unit(i, k): 1.0 for i in range(k)}
deg = k * (k + 1) // 2
target = tuple([2 * k - 1] * k)
pref = (-1) ** (k * (k - 1) // 2) / factorial(k)
coef = []
G = dict(V)
for m in range(deg + 1):
    tot = 0j
    for e, v in G.items():
        b = tuple(t - x for t, x in zip(target, e))
        if min(b) < 0 or max(b) >= N: continue
        tot += v * f[b]
    coef.append(pref * tot / (2 ** m * factorial(m)))
    G = pmul(G, S1)
c = [complex(x) for x in coef]
print(json.dumps({"k": k, "N": N, "r": r, "P0": P0, "coef_re": [x.real for x in c], "coef_im_max": max(abs(x.imag) for x in c)}))
