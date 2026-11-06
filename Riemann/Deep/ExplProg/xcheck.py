#!/usr/bin/env python3
"""xcheck.py -- independent float cross-check of the closed-form zero sums in cert.py:
upper Stieltjes sums  N+(5/7) (G(0) - G(5/7)) + sum_i Nbmor(t_{i+1}) (G(t_i) - G(t_{i+1}))  with the UNdecomposed BMOR bound
max(0, (t/pi) log(qt/2 pi e) + 1/4 + 0.22737 l + 2 log(1+l) - 1/2), l = log(q(t+2)/2pi), and the Rosser count for zeta.
The closed forms (decomposed, hence larger) must be >= these sums (safe direction) and close to them for L >= 9
(at L = log 9 the decomposition 2 log(1+l) <= 2 log(1+L) + 2 c+/(1+L) and the dropped max(0, .) cost more).
Run: uv run --no-project --with python-flint --with mpmath python3 xcheck.py"""
import math, sys
sys.path.insert(0, '.')
import cert
from cert import arb, A, Fr
def nb(t, L):
    l = L + math.log((t + 2) / (2 * math.pi))
    v = (t / math.pi) * (L + math.log(t / (2 * math.pi * math.e))) + 0.25 + 0.22737 * l + 2 * math.log1p(l) - 0.5
    return max(0.0, v)
def nz(t):
    return (t / math.pi) * math.log(t / (2 * math.pi * math.e)) + 4.926 + 0.274 * math.log(t) + 0.886 * math.log(math.log(t))
def G(t, n0, k1, k2):
    return n0 if t == 0 else min(n0, k1 / t, k2 / t / t)
def stieltjes(n0, k1, k2, L=None, zeta=False, r=1.0005, T=1e9):
    t = 14.0 if zeta else 5 / 7
    g = G(t, n0, k1, k2)
    s = 0.0 if zeta else nb(t, L) * (n0 - g)          # zeros below 5/7 at height 0; zeta has none below 14
    while t < T:
        t2 = t * r; g2 = G(t2, n0, k1, k2)
        s += (nz(t2) if zeta else nb(t2, L)) * (g - g2)
        t, g = t2, g2
    s += (nz(t) if zeta else nb(t, L)) * g          # crude tail
    return s
K = cert.KER
(n0, k1, k2), (m0, k1L, k2L) = K.ks(arb(80))
f = lambda x: float(x.mid())
tests = []
for L in (math.log(9), 9.21, 50.0):
    zc = f(cert.zsum_dir(arb(L), n0, k1, k2)); zs = stieltjes(f(n0), f(k1), f(k2), L)
    tests.append(("Linnik kernel, L=%.2f" % L, zc, zs))
for d in (Fr(1, 2), Fr(1, 10), Fr(1, 1000)):
    (WB, A1, A2), _ = cert.trap_consts(A(d), arb(1) / 8)
    zc = f(cert.zsum_dir(arb(12.0), WB, A1, A2)); zs = stieltjes(f(WB), f(A1), f(A2), 12.0)
    tests.append(("trapezoid delta=%s, L=12" % d, zc, zs))
    zc = f(cert.zsum_zeta(WB, A1, A2)); zs = stieltjes(f(WB), f(A1), f(A2), zeta=True)
    tests.append(("trapezoid delta=%s, zeta" % d, zc, zs))
Gk, _ = cert.lp_kernels(Fr(7, 10))
zc = f(cert.zsum_dir(arb(100.0), *Gk)); zs = stieltjes(*[f(v) for v in Gk], 100.0)
tests.append(("L'/L kernel h=0.7, L=100", zc, zs))
ok = True
for (name, zc, zs) in tests:
    rel = (zc - zs) / zs
    flag = "ok" if -1e-6 < rel < (0.30 if "L=2.20" in name else 0.02) else "CHECK"
    ok &= flag == "ok"
    print("%-32s closed form %12.5f   Stieltjes %12.5f   (closed - sum)/sum = %+.5f  %s" % (name, zc, zs, rel, flag))
print("ALL OK" if ok else "SOME CHECKS FLAGGED")
