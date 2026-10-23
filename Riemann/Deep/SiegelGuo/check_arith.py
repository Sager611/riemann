#!/usr/bin/env python3
"""Data tests of Proposition arith:main (TeX l.417-838) of arXiv 2610.30476.
(a) S_-(x,chi) <= 11/20 log x + 29/60 log q + 6 for all fundamental discriminants |D| <= DMAX, x <= XMAX
    (sup over x is attained at inert primes; report the largest excess over 11/20 log x + 29/60 log q).
(b) For p0 >= max(10000, 4q) inert: exact H = sum_P H_p log p versus the counting bound
    n log n + sum log w_{r,s}, and the paper's displayed lower/upper forms (arith:H-lower), (arith:H-upper)."""
import numpy as np, math, sys
DMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 10000
XMAX = 10**6
s = np.ones(XMAX + 1, dtype=np.bool_); s[:2] = False
for i in range(2, int(XMAX**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.nonzero(s)[0].astype(np.int64); del s
lp = np.log(P.astype(np.float64))
def sqfree(m):
    m = abs(m)
    for p in range(2, int(m**0.5) + 1):
        if m % (p * p) == 0: return False
    return True
def fundamentals(dmax):
    for D in range(-dmax, dmax + 1):
        if D in (0, 1): continue
        if D % 4 == 1 and sqfree(D): yield D, D
        elif D % 4 == 0 and (D // 4) % 4 in (2, 3) and sqfree(D // 4): yield D, D // 4
def chi_primes(D):
    """Kronecker (D/p) for all p in P (vectorised Euler criterion for odd p)."""
    out = np.zeros(len(P), dtype=np.int64)
    if D % 2 != 0:
        out[0] = 1 if D % 8 in (1, 7) else -1
    Po = P[1:]
    base = np.mod(D, Po)
    e = (Po - 1) // 2
    res = np.ones_like(Po); b = base.copy()
    while np.any(e > 0):
        odd = (e & 1) == 1
        res = np.where(odd, (res * b) % Po, res)
        b = (b * b) % Po
        e >>= 1
    leg = np.where(base == 0, 0, np.where(res == 1, 1, -1))
    out[1:] = leg
    return out
worst = (-1e9, None); count = 0; viol = 0
top = []
for D, d in fundamentals(DMAX):
    q = abs(D); chi = chi_primes(D)
    inert = chi == -1
    Sm = np.cumsum(np.where(inert, lp / P, 0.0))
    exc = (Sm - 0.55 * lp)[inert]
    if exc.size == 0: continue
    j = int(np.argmax(exc)); e = float(exc[j]) - 29 / 60 * math.log(q)
    count += 1
    if e > 6: viol += 1
    top.append((e, D, int(P[inert][j])))
    if e > worst[0]: worst = (e, D, int(P[inert][j]))
top.sort(reverse=True)
print(f"(a) {count} fundamental discriminants |D|<={DMAX}, x<=1e6: violations of arith:target = {viol}")
print("    largest excess S_-(x) - 11/20 log x - 29/60 log q (bound 6):", [(round(e, 3), D, x) for e, D, x in top[:6]])
# ---------- (b) exact H at real parameters
def Fk(k, s0):
    return k * (k - 1) // 2 + (s0 - k) * k if k < s0 else s0 * (s0 - 1) // 2
def Hp(p, p0, s0):
    a, b = divmod(p0, p)
    return sum(p * Fk(k, s0) for k in range(a)) + b * Fk(a, s0)
print("(b) exact exponent mass H versus counting bound (all logs natural)")
for D, d in [(-3, -3), (-4, -1), (5, 5), (-163, -163), (-1003, -1003), (8, 2), (-8, -2), (-4 * 2003, -2003)]:
    q = abs(D); chi = chi_primes(D)
    cand = [int(p) for p, c in zip(P, chi) if c == -1 and p >= max(10000, 4 * q) and p % 2]
    p0 = cand[0]
    s0 = p0 // 30; n = p0 * s0; m = math.isqrt(p0 - 1) + 1; b0 = 2 * math.sqrt(p0 * q)
    Pset = [int(p) for p, c in zip(P, chi) if c == -1 and p % 2 and p <= p0]
    H = sum(Hp(p, p0, s0) * math.log(p) for p in Pset)
    U = sum(math.log(p) / p for p in Pset)
    th = sum(math.log(p) for p in P if p <= p0)
    upper = n * math.log(n) + n * (s0 - 1) * math.log(p0 / b0) + p0 * s0 * (p0 + s0 - 2) / 2 * math.log(b0)
    nrm = p0 * s0 * (p0 - 1) / 2
    up_disp = 11 / 20 * math.log(p0) + 29 / 60 * math.log(q) + 29 / 30 * math.log(2) + 4 * math.log(p0) / (p0 - 1)
    lo_disp = U - U / (2 * s0) - th / (p0 - 1) - 7 * (2 * p0 - 1) / (120 * s0)
    ok = (H <= upper) and (H / nrm >= lo_disp) and (upper / nrm <= up_disp + 1e-12)
    print(f"  D={D:6d} p0={p0} s0={s0} m={m}: H/norm={H/nrm:.4f}  counting/norm={upper/nrm:.4f}  "
          f"(arith:H-lower)={lo_disp:.4f} (arith:H-upper)={up_disp:.4f}  U={U:.4f}  consistent={ok}")

