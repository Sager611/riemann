#!/usr/bin/env python3
"""For sampled generic conics B_delta(A,q) = -4n (delta = [[al,be],[g,e]] in SL2(Z)), count ALL integer points (A,q) with
A = 4a, q odd, Q <= q < 2Q, 1/2 <= a/q <= 2 (the box of the second spacing problem), exactly (python ints).
Compares with the Bombieri-Pila bound B^{1/2} for a conic in a box of side B = 8Q. usage: conic_points.py dumpfile"""
import sys, math
L = open(sys.argv[1]).read().split('\n'); Q = int(L[0].split()[1])
rows = [list(map(int, l.split())) for l in L[1:] if l.strip()]
cnt = []
for al, be, g, e, n in rows:
    a2, a1, a0 = al*g, 2*be*g, be*e      # B(A,q) = a2 A^2 + a1 A q + a0 q^2 ; solve B = -4n for A, each q
    c = 0
    for q in range(Q | 1, 2*Q, 2):
        # a2 A^2 + (a1 q) A + (a0 q^2 + 4n) = 0
        D = (a1*q)**2 - 4*a2*(a0*q*q + 4*n)
        if D < 0: continue
        s = math.isqrt(D)
        if s*s != D: continue
        for num in (-a1*q + s, -a1*q - s):
            if a2 != 0 and num % (2*a2) == 0:
                A = num // (2*a2)
                if A % 4 == 0 and 2*A >= 4*q and A <= 8*q and math.gcd(A//4, q) == 1: c += 1
            if s == 0: break
    cnt.append(c)
cnt.sort()
print(f"Q={Q} sampled generic conics={len(cnt)} points in box per conic: min={cnt[0]} median={cnt[len(cnt)//2]} max={cnt[-1]} mean={sum(cnt)/len(cnt):.3f}; Bombieri-Pila sqrt(8Q)={math.sqrt(8*Q):.1f}")
from collections import Counter; print("histogram", sorted(Counter(cnt).items()))
