#!/usr/bin/env python3
"""Check of Proposition 'alternating products' (Section 4) on all rational-type solutions.

x in Z[i][C_N] whose character values depend only on the character order d: x(chi_d) = alpha_d.
Such x has coefficients in Z[i] iff  k | sum_{j|k} mu(k/j) G(j)  for every k | N,  G(k) = alpha_{N/k}
(x = sum_h a_h J_h, J_h = sum of the subgroup of order h, a_k k = sum_{j|k} mu(k/j) G(j)).
x x^* = u^2  <=>  |alpha_d|^2 = u^2 for all d.
Delta(x) = prod_{S subset I} alpha_{d(S)}^{(-1)^{|S|}},  d(S) = prod_{p in S} p.
Paper: if N = u^2 then Delta(x) is a root of unity of p-power order for every p | u (so = 1 here, values in Q(i)).
Control: with N = rad(u) (scalar exponent 2 > group exponent 1) Delta = -1 must occur.
"""
import sys
from fractions import Fraction
from itertools import product

def divisors(n): return [d for d in range(1, n+1) if n % d == 0]
def mu(n):
    r, m, p = 1, n, 2
    while p*p <= m:
        if m % p == 0:
            m //= p
            if m % p == 0: return 0
            r = -r
        p += 1
    return -r if m > 1 else r
def primes(n): return [p for p in range(2, n+1) if n % p == 0 and all(p % q for q in range(2, p))]
def gauss_of_norm(n):
    out = []
    m = int(n**0.5) + 1
    for a in range(-m, m+1):
        for b in range(-m, m+1):
            if a*a + b*b == n: out.append((a, b))
    return out
def gmul(a, b): return (a[0]*b[0]-a[1]*b[1], a[0]*b[1]+a[1]*b[0])
def gdiv_exact_is_one(num, den):
    return num == den

def run(u, N):
    I = primes(u)
    divs = divisors(N)
    vals = gauss_of_norm(u*u)
    order = sorted(divs)
    sub = {k: [j for j in divs if k % j == 0] for k in divs}
    sols = []
    G = {}
    def rec(idx):
        if idx == len(order):
            sols.append(dict(G)); return
        k = order[idx]
        for v in vals:
            G[k] = v
            sr = si = 0
            for j in sub[k]:
                m = mu(k // j)
                if m: sr += m*G[j][0]; si += m*G[j][1]
            if sr % k == 0 and si % k == 0:
                rec(idx+1)
        del G[k]
    rec(0)
    bad = []
    from math import prod
    for G in sols:
        alpha = {d: G[N // d] for d in divs}
        num = (1, 0); den = (1, 0)
        for r in range(len(I)+1):
            from itertools import combinations
            for S in combinations(I, r):
                d = prod(S) if S else 1
                if r % 2 == 0: num = gmul(num, alpha[d])
                else: den = gmul(den, alpha[d])
        if num != den: bad.append((alpha, num, den))
    return I, len(sols), bad

if __name__ == '__main__':
    ok = True
    for (u, N) in [(3, 9), (5, 25), (9, 81), (15, 225), (21, 441), (3*5*7, 0)]:
        if N == 0: continue
        I, n, bad = run(u, N)
        print(f"u={u:3d} N=u^2={N:4d} primes {I}: {n} integral rational x with xx*=u^2; Delta != 1 for {len(bad)}")
        ok &= (len(bad) == 0)
    # control: group exponent below scalar exponent
    for (u, N) in [(3, 3), (5, 5), (15, 15), (15, 45)]:
        I, n, bad = run(u, N)
        ex = bad[0] if bad else None
        msg = (f"e.g. alpha={ {d: a for d, a in ex[0].items()} }, Delta = {ex[1]}/{ex[2]}") if ex else ""
        print(f"control u={u} N={N}: {n} integral x with xx*=u^2 in Z[i][C_{N}]; Delta != 1 for {len(bad)}  {msg}")
    print("SUMMARY alternating products:", "consistent" if ok else "INCONSISTENT")

