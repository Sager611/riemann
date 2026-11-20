#!/usr/bin/env python3
"""Test A (remark item 2, classical input of the core): at a CM point, derivatives D^a of the
Weierstrass function are divisible by p^(a/p + a/p^3 + ... - O(log a)) for p INERT in K
(supersingular reduction, height-2 formal group), but not for p split.
Curves: y^2 = 4x^3 - 4x (g2=4, g3=0, CM by Z[i]) and y^2 = 4x^3 - 4 (g2=0, g3=4, CM by Z[w]).
wp(z) = z^-2 + sum_k c_k z^(2k); D^(2k) of the regular part at 0 is (2k)! c_k."""
import sys, math
from flint import fmpq
K = int(sys.argv[1]) if len(sys.argv) > 1 else 600
def coeffs(g2, g3, K):
    c = [None, fmpq(g2, 20), fmpq(g3, 28)]
    for k in range(3, K + 1):
        s = fmpq(0)
        for m in range(1, k - 1):
            s += c[m] * c[k - 1 - m]
        c.append(fmpq(3, (2 * k + 3) * (k - 2)) * s)
    return c
def vp(n, p):
    if n == 0: return 10**9
    v = 0
    while n % p == 0: n //= p; v += 1
    return v
def vq(x, p):
    return vp(int(x.p), p) - vp(int(x.q), p)
def vfact(n, p):
    s, q = 0, p
    while q <= n: s += n // q; q *= p
    return s
for name, g2, g3, inert, split in [("Z[i]: y^2=4x^3-4x", 4, 0, [3, 7, 11, 19, 23], [5, 13, 17, 29]),
                                   ("Z[w]: y^2=4x^3-4", 0, 4, [5, 11, 17, 23], [7, 13, 19, 31])]:
    c = coeffs(g2, g3, K)
    print("==", name, "k <=", K)
    for p in inert + split:
        kind = "inert" if p in inert else "split"
        worst_odd, worst_all, worst_zero = 10**9, 10**9, 10**9
        rows = []
        for k in range(1, K + 1):
            if c[k] == 0: continue
            a = 2 * k
            d = vq(c[k], p) + vfact(a, p)          # v_p(D^a wp(0))
            odd = a * p / (p * p - 1)                # a/p + a/p^3 + ...
            allp = a / (p - 1)                       # a/p + a/p^2 + ... ~ v_p(a!)
            L = math.log(a) / math.log(p)
            worst_odd = min(worst_odd, (d - odd) / max(L, 1))
            worst_all = min(worst_all, (d - allp) / max(L, 1))
            worst_zero = min(worst_zero, d / max(L, 1))
            if k in (K // 4, K // 2, K): rows.append((a, d, round(odd, 1), round(allp, 1)))
        print(f"  p={p:3d} {kind}: min_k (v - [a/p+a/p^3+..])/log_p a = {worst_odd:8.2f} | min (v - a/(p-1))/log_p a = {worst_all:8.2f} | min v/log_p a = {worst_zero:6.2f} | (a, v, odd, all) {rows}")
