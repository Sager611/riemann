#!/usr/bin/env python3
"""Siegel-zero model sanity check for the arithmetic step (TeX l.640-772).
Model: every odd prime <= p0 inert (the pattern a Siegel zero forces). Compare the exact exponent mass
H = sum_p H_p log p with the counting bound n log n + sum log w (q = 3). If the model ever exceeds the
bound, the counting inequality genuinely excludes it -- so the arithmetic step is load-bearing, and it
holds for real characters only through the Frobenius divisibility (verified in check_inject.py)."""
import numpy as np, math
N = 2 * 10**7
s = np.ones(N + 1, dtype=np.bool_); s[:2] = False
for i in range(2, int(N**0.5) + 1):
    if s[i]: s[i*i::i] = False
P = np.nonzero(s)[0]; del s
def Hp(p, p0, s0):
    a, b = divmod(p0, p)
    k = np.arange(a, dtype=np.float64)
    F = np.where(k < s0, k * (k - 1) / 2 + (s0 - k) * k, s0 * (s0 - 1) / 2)
    Fa = a * (a - 1) / 2 + (s0 - a) * a if a < s0 else s0 * (s0 - 1) / 2
    return p * F.sum() + b * Fa
q = 3
inert3 = lambda p: p % 3 == 2           # (-3/p) = -1 iff p = 2 mod 3
for p0 in [10007, 100043, 1000037, 10000079, 19999999]:
    while not (p0 % 3 == 2 and all(p0 % k for k in range(2, math.isqrt(p0) + 1))): p0 += 1
    s0 = p0 // 30; n = p0 * s0; b0 = 2 * math.sqrt(p0 * q); nrm = p0 * s0 * (p0 - 1) / 2
    cnt = (n * math.log(n) + n * (s0 - 1) * math.log(p0 / b0) + p0 * s0 * (p0 + s0 - 2) / 2 * math.log(b0)) / nrm
    odd = P[(P >= 3) & (P <= p0)]
    Hm = sum(Hp(int(p), p0, s0) * math.log(p) for p in odd) / nrm
    Ht = sum(Hp(int(p), p0, s0) * math.log(p) for p in odd if inert3(int(p))) / nrm
    print(f"p0={p0:>9}: counting bound/norm={cnt:7.3f}  H/norm all-inert model={Hm:7.3f} ({'EXCEEDS' if Hm > cnt else 'below'})  true chi_-3={Ht:6.3f}")

