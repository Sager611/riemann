#!/usr/bin/env python3
"""selbaudit check 5: exact finite sums vs the main-term profile used in §8 (q = 1, d = delta = 1):
S(v) = sum_{n<xi} M(n) tau(n) n^{-1-iv} exactly; main: s = -Delta/((1-theta) 2 sqrt(pi) sqrt(v) log xi), Delta = Delta_2(X) - Delta_2(theta X), X = v log xi.
d_exact(X) := -S(v)(1-theta)2 sqrt(pi) sqrt(v)/X  vs  d(X;0) = [Delta_2(X) - Delta_2(theta X)]/X; and J over x in [1,B] with both (B <= X range of Lemma 12).
uv run --no-project --with numpy --with mpmath python3 exactsums.py"""
import math, numpy as np, mpmath as mp
def tau_pk(a):
    r = 1.0
    for i in range(1, a + 1): r *= (i - 1.5)/i
    return r
def tau_arr(n):
    t = np.ones(n + 1); t[0] = 0; spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, n + 1):
        if spf[i] == 0: spf[i::i][spf[i::i] == 0] = i
    for m in range(2, n + 1):
        p = spf[m]; r = m; a = 0
        while r % p == 0: r //= p; a += 1
        t[m] = t[r]*tau_pk(a)
    return t
psi = lambda Y: complex(mp.sqrt(mp.pi)*mp.exp(-1j*mp.pi/4)*mp.erf(mp.exp(1j*mp.pi/4)*mp.sqrt(Y)))
D2 = lambda Y: -2*math.sqrt(Y)*complex(math.cos(Y), -math.sin(Y)) - (1 + 2j*Y)*psi(Y)
TH = 0.15
for xi in [10**4, 10**6]:
    L = math.log(xi); t = tau_arr(xi - 1); n = np.arange(1, xi); M = (np.log(xi/n) - np.maximum(np.log(xi**TH/n), 0))/((1 - TH)*L)
    w = t[1:]*M/n; ln = np.log(n)
    Xmax = math.sqrt(L)   # Lemma 12 range v <= 1/sqrt(log xi)  <=>  X <= sqrt(log xi)
    rows = []; xs = np.linspace(0.05, Xmax, 40); dE = []; dM = []
    for X in xs:
        v = X/L; S = np.sum(w*np.exp(-1j*v*ln))
        de = -S*(1 - TH)*2*math.sqrt(math.pi)*math.sqrt(v)/X*math.sqrt(L)/math.sqrt(L)
        de = -S*(1 - TH)*2*math.sqrt(math.pi)*math.sqrt(L)*math.sqrt(X)/X/ math.sqrt(L) * math.sqrt(L)/math.sqrt(L)
        # s = -Delta/((1-th)2 sqrt(pi) sqrt(v) L)  =>  Delta/X = -s (1-th) 2 sqrt(pi) sqrt(v) L / X
        de = -S*(1 - TH)*2*math.sqrt(math.pi)*math.sqrt(v)*L/X
        dm = (D2(X) - D2(TH*X))/X
        dE.append(de); dM.append(dm)
    dE = np.array(dE); dM = np.array(dM); rel = np.abs(dE - dM)/np.abs(dM)
    print("xi=1e%d: X in [0.05, %.2f]: |d_exact - d_main|/|d_main| max %.3f, median %.3f; at X=%.2f: exact %.4f%+.4fi main %.4f%+.4fi" % (round(math.log10(xi)), Xmax, rel.max(), np.median(rel), xs[-1], dE[-1].real, dE[-1].imag, dM[-1].real, dM[-1].imag))
    # J on [1, B] with rho = 2.0 (worst case), both profiles, trapezoid on the 40-point grid restricted to X >= 1 is too coarse: use 400 points
    xg = np.linspace(1.0, Xmax, 400); h = xg[1] - xg[0]
    de = np.array([-np.sum(w*np.exp(-1j*(X/L)*ln))*(1 - TH)*2*math.sqrt(math.pi)*math.sqrt(X/L)*L/X for X in xg])
    dm = np.array([(D2(X) - D2(TH*X))/X for X in xg])
    for rho in [2.0, 6.0]:
        B = xg[-1]; f = lambda d: np.sum(((B - xg)*np.exp(1j*rho*xg)*d*d)[:-1] + ((B - xg)*np.exp(1j*rho*xg)*d*d)[1:])*h/2
        print("   rho=%.0f, B=%.2f: J_exact=%.4f%+.4fi  J_main=%.4f%+.4fi" % (rho, B, f(de).real, f(de).imag, f(dm).real, f(dm).imag), flush=True)
