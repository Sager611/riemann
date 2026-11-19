#!/usr/bin/env python3
"""tail_check.py (ppraudit) -- falsification / recomputation of the tail (DERIVATION sec 5) inputs.
(i) TRUE G_q(xi_q) (exact enumeration of squarefree d < xi with P+(d) < q) vs tail.py's certified G_lo, every prime
    q in [50, 1e5] at u = u1 and u1 + 0.1; true Part A vs tail.py's partA bound.
(ii) q > Q2 = 1e5 (first 300 primes) at u = u1, u1 + 1: true Phi_q - 1 vs psi(a_q).
(iii) constants: c_sigma, h_Q'(a_cap) > 0, u1, L1, m_delta, K_B, R_B(c), (1 + eps_V) I.
Usage: python3 tail_check.py"""
import sys, math, time
import numpy as np
import mpmath as M
sys.path.insert(0, "../PrimRootNum"); sys.dont_write_bytecode = True
from fractions import Fraction as Fr
import tail as TL
M.mp.prec = 128
t0 = time.time()
s = 2.3; Q2 = 100000
T = TL.Tail(Fr(23, 10), w0=50, Q2=Q2)
NMAX = 2 * 10 ** 6
# sieve: phi, squarefree, largest prime factor for d < NMAX
phi = np.arange(NMAX, dtype=np.float64); sq = np.ones(NMAX, dtype=bool); lpf = np.ones(NMAX, dtype=np.int64)
isp = np.ones(NMAX, dtype=bool); isp[:2] = False
for p in range(2, NMAX):
    if isp[p]:
        if p * p < NMAX: isp[p * p::p] = False
        phi[p::p] *= (1 - 1.0 / p); lpf[p::p] = p
        if p * p < NMAX: sq[p * p::p * p] = False
w = np.where(sq, 1.0 / np.maximum(phi, 1), 0.0); w[0] = 0.0
primes = np.nonzero(isp)[0]
def Gtrue(q, xi):
    n = int(math.ceil(xi))                     # d < xi  <=> d <= ceil(xi) - 1
    if n > NMAX: return None
    sel = w[1:n] if q > n else w[1:n][lpf[1:n] < q]
    return float(sel.sum()) + 0.0
u1 = math.log(104743)
print("[i] G_lo vs true G at u = u1, u1 + 0.1 for primes q in [50, 1e5]")
for u in (u1, u1 + 0.1):
    worst = 9.0; worstq = None; trueA = 0.0; nq = 0; skipped = 0
    for q in T.PA:
        lx = (s * u - math.log(q)) / 2; xi = math.exp(lx)
        G = Gtrue(q, xi)
        if G is None: skipped += 1; continue
        glo = float(T.G_lo(q, TL.A(Fr(lx).limit_denominator(10 ** 15))).upper())
        r = glo / G
        if r > worst or worstq is None: worst, worstq = r, (q, round(xi), glo, G)
        trueA += (1.0 / G - float(T.V[q].mid())) / q; nq += 1
    print("   u = %.4f: %d primes (%d skipped, xi > NMAX); max G_lo/G_true = %.6f at (q, xi, G_lo, G) = %s (must be <= 1); true Part A = %.6g" % (u, nq, skipped, worst, worstq, trueA))
uu = TL.A(Fr(u1).limit_denominator(10 ** 15))
pa = float(T.partA(uu, uu + TL.A(Fr(1, 10))).upper())
print("   tail.py partA(u1, u1 + 0.1) upper bound = %.6g (true value above must be smaller)" % pa)
print("[ii] true Phi_q - 1 vs psi(a_q) for the first 300 primes q > 1e5")
for u in (u1, u1 + 1.0):
    worst = 0.0; wq = None
    for q in [int(x) for x in primes[primes > Q2][:300]]:
        lx = (s * u - math.log(q)) / 2; xi = math.exp(lx); a = lx / math.log(q)
        G = Gtrue(q, xi)
        V = float(T.V[max(k for k in T.V)].mid())  # placeholder, replaced below
        V = math.exp(float(np.sum(np.log1p(-1.0 / primes[primes < q]))))
        true = 1.0 / (V * G) - 1
        ps = float(T.psi(Fr(a).limit_denominator(10 ** 12)).upper())
        if true / ps > worst: worst, wq = true / ps, (q, round(xi), true, ps)
    print("   u = %.3f: max (Phi_q - 1)/psi(a_q) = %.4f at (q, xi, true, psi) = %s (must be <= 1)" % (u, worst, wq))
print("[iii] constants")
cs = M.mpf(1)
for p in primes: cs *= 1 + M.mpf(2) / (int(p) * (int(p) - 1))
print("   c_sigma: prod_{p < 2e6} = %s, + tail exp(2/(2e6-1)) = %s (cert <= 3.2795975)" % (M.nstr(cs, 12), M.nstr(cs * M.e ** (M.mpf(2) / (NMAX - 1)), 12)))
LQ = M.log(Q2); etaB = 1 / LQ ** 2 + M.mpf(1) / Q2
hp = lambda a: 1 - M.log(a) - etaB - 1 / (a * LQ)
print("   h_Q'(2.609) = %s (> 0 needed), argmax h_Q = %s" % (M.nstr(hp(M.mpf(2609) / 1000), 8), M.nstr(M.findroot(hp, 2.6), 10)))
P = [int(x) for x in primes[:10001]]
L1 = M.fsum(M.log(q) for q in P); u1m = M.log(P[-1])
print("   p_10001 = %d, u1 = %s, L1 = log(P_10001 + 1) ~ theta(p_10001) = %s (cert 104403.76)" % (P[-1], M.nstr(u1m, 10), M.nstr(L1, 10)))
kz = 1 - 1 / u1m; ell1 = M.log(L1 / kz)
print("   delta*ell(L1) = %s (>= 1), m_delta = %s (cert 0.272694)" % (M.nstr(ell1 / 10, 6), M.nstr(L1 ** (M.mpf(1) / 10) / ell1, 8)))
KB = (1 / M.mpf(0.65) ** 2) * (1 / (2 * LQ ** 2)) * (1 + 2 / LQ ** 2)
print("   K_B = %s (cert 0.00906308)" % M.nstr(KB, 8))
print("   acap used = %s; Psi(1) = psi(0.65) = %s (cert 1.75049)" % (T.acap, M.nstr(M.mpf(str(float(T.psi(Fr(13, 20)).upper()))), 8)))
print("done %.0f s" % (time.time() - t0))
