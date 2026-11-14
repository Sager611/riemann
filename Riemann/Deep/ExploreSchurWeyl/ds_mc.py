#!/usr/bin/env python3
"""Monte Carlo sanity check of the Diaconis-Shahshahani moments used as constraints (explore50).
Haar U(N) by QR of complex Ginibre (phase-corrected); E p_lam conj(p_mu) vs delta_{lam,mu} z_lam for |lam|=|mu| <= N+1."""
import sys, math, numpy as np
def parts(n, mx=None):
    mx = n if mx is None else mx
    if n == 0: yield (); return
    for k in range(min(n, mx), 0, -1):
        for r in parts(n-k, k): yield (k,) + r
def zl(lam):
    z = 1
    for j in set(lam):
        a = lam.count(j); z *= j**a*math.factorial(a)
    return z
rng = np.random.default_rng(7)
for N in (3, 4):
    S = 200000; P = np.zeros((S, N+2), dtype=complex)
    for c in range(0, S, 20000):
        G = (rng.standard_normal((20000, N, N)) + 1j*rng.standard_normal((20000, N, N)))/np.sqrt(2)
        Q, R = np.linalg.qr(G); d = np.diagonal(R, axis1=1, axis2=2); Q = Q*(d/np.abs(d))[:, None, :]
        ev = np.linalg.eigvals(Q)
        for k in range(N+2): P[c:c+20000, k] = (ev**k).sum(axis=1)
    worst_in, rows_out = 0.0, []
    for n in range(1, N+2):
        L = list(parts(n))
        for i, lam in enumerate(L):
            for mu in L[i:]:
                v = np.prod([P[:, k] for k in lam], axis=0)*np.conj(np.prod([P[:, k] for k in mu], axis=0))
                m = v.mean(); se = v.std()/np.sqrt(S); t = zl(lam) if lam == mu else 0
                dev = abs(m - t)/max(se, 1e-12)
                if n <= N: worst_in = max(worst_in, dev)
                else: rows_out.append((lam, mu, m.real, t, dev))
    print(f"N={N}: |lam|<=N: max |MC - DS|/se = {worst_in:.2f} over all pairs")
    for lam, mu, m, t, dev in sorted(rows_out, key=lambda r: -r[4])[:3]:
        print(f"   |lam|=N+1 {lam}|{mu}: MC {m:.3f} vs Wick {t}  ({dev:.1f} se)")
