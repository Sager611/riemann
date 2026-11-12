#!/usr/bin/env python3
"""selbaudit check 2 (DERIVATION §6 + §10): exact Sigma_1 = sum_{h,k} chi0(hk) gamma(h)gamma(k)(h,k)/(hk) at finite xi,
against the bound 16 P * int_0^2 Phi_BB^2 /(pi^2 (1-theta)^2) / log xi with P = P1eff (~1, §10) and P = P1 (11.543, absolute).
Phi_BB(t) = min(1, (1 - t/2)/(1 - theta)) exactly, so int Phi_BB^2 = 2 theta + 2(1 - theta)/3.
uv run --no-project --with numpy python3 sigma1.py"""
import math, sys, time
import numpy as np
T0 = time.time()
def tau_pk(a):
    r = 1.0
    for i in range(1, a + 1): r *= (i - 1.5)/i
    return r
def tau_arr(n):
    t = np.ones(n + 1); t[0] = 0
    spf = np.zeros(n + 1, dtype=np.int64)
    for i in range(2, n + 1):
        if spf[i] == 0: spf[i::i][spf[i::i] == 0] = i
    for m in range(2, n + 1):
        p = spf[m]; r = m; a = 0
        while r % p == 0: r //= p; a += 1
        t[m] = t[r]*tau_pk(a)
    return t
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0]
def sigma1(xi, th, q=1):
    t = tau_arr(xi); L = math.log(xi)
    a = np.arange(1, xi); M = (np.log(xi/a) - np.maximum(np.log(xi**th/a), 0))/((1 - th)*L)
    f = np.zeros(xi); f[1:] = t[1:xi]*M
    if q > 1: f[np.arange(xi) % q == 0] = 0.0          # chi0(n): the convolution of chi0-supported sequences
    N = (xi - 1)**2
    g = np.zeros(N + 1)
    for x in range(1, xi):
        if f[x] != 0.0: g[x*np.arange(1, xi)] += f[x]*f[1:]
    g[1:] /= np.arange(1, N + 1)
    A = g                                              # in place: A[d] = sum_{d|n} g[n]
    for p in primes_upto(N):
        k = 1; hi = N//p
        while hi >= 1:
            lo = N//(p**(k + 1)) if p**(k + 1) <= N else 0
            A[lo + 1:hi + 1] += A[(lo + 1)*p:hi*p + 1:p]
            k += 1; hi = lo
    phi = np.arange(N + 1, dtype=np.int64)
    for p in primes_upto(N): phi[p::p] -= phi[p::p]//p
    S = float(np.sum(phi[1:]*A[1:]**2))
    return S, A, phi
# brute-force validation at tiny xi
xi = 12; th = 0.3
S, A, phi = sigma1(xi, th)
t = tau_arr(xi); a = np.arange(1, xi); M = (np.log(xi/a) - np.maximum(np.log(xi**th/a), 0))/((1 - th)*math.log(xi))
f = np.zeros(xi); f[1:] = t[1:xi]*M; N = (xi - 1)**2; gam = np.zeros(N + 1)
for x in range(1, xi):
    for y in range(1, xi): gam[x*y] += f[x]*f[y]
bf = sum(gam[h]*gam[k]*math.gcd(h, k)/(h*k) for h in range(1, N + 1) if gam[h] for k in range(1, N + 1) if gam[k])
print("validation xi=12: transform %.12f  brute force %.12f" % (S, bf)); assert abs(S - bf) < 1e-12*max(1, abs(bf))
P1 = 11.542935; P1eff = 1.000003001
for th in [0.15, 0.05]:
    i0 = 2*th + 2*(1 - th)/3; c = 16*i0/(math.pi**2*(1 - th)**2)
    print("theta=%.2f: bound Sigma1*log(xi) <= %.4f (P1eff), %.3f (P1)" % (th, c*P1eff, c*P1))
    for q in [1, 13]:
        for xi in [100, 300, 1000, 2000]:
            S, A, phi = sigma1(xi, th, q)
            fac = q/(q - 1) if q > 1 else 1.0   # bound carries q/phi(q)
            print("  q=%-2d xi=%-5d Sigma1*log(xi) = %.5f   ratio to P1eff-bound %.3f, to P1-bound %.4f   (%.0fs)" % (q, xi, S*math.log(xi), S*math.log(xi)/(fac*c*P1eff), S*math.log(xi)/(fac*c*P1), time.time() - T0), flush=True)
