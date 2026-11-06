"""Auxiliary primes of OpenAI's factoring paper (4 Oct 2026, eq. (auxiliary)):
   l prime, l not in {2,3,p,q}, l = 1 mod 12q, p^((l-1)/q) != 1 mod l.
Checks the bound of DERIVATION Thm B:  l_q(p) <= (0.49 (log 3 + D(q,p)))^32,
   D(q,p) = q phi(12q) log(12q) + phi(12q) (q log q + (q-1) log p)   (upper bound for log Delta_L),
and the simplified form l_q(p) <= (1.96 q^2 log(12 q^2 p))^32.
Part 1: every prime 5 <= p <= PMAX with p not dividing 6q, q in Q.  Part 2: adversarial p = 1 mod (first k such l).
Run: uv run --no-project --with sympy python3 aux.py
"""
import math, sys
from sympy import isprime, primerange, totient
PMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 200000
Q = [2, 3, 5, 7, 11, 13]
def D(q, p):
    f = int(totient(12 * q))
    return q * f * math.log(12 * q) + f * (q * math.log(q) + (q - 1) * math.log(p))
def log10_bound(q, p):
    return 32 * math.log10(0.49 * (math.log(3) + D(q, p)))
def log10_simple(q, p):
    return 32 * math.log10(1.96 * q * q * math.log(12 * q * q * p))
def least_aux(p, q, cap=10**9):
    l = 12 * q + 1
    while l < cap:
        if l != p and isprime(l) and pow(p, (l - 1) // q, l) != 1:
            return l
        l += 12 * q
    return None
print("Part 1: all primes 5 <= p <= %d, p not dividing 6q" % PMAX)
for q in Q:
    mx, mp, worst_slack = 0, 0, 1e9
    cnt = 0
    for p in primerange(5, PMAX + 1):
        if (6 * q) % p == 0:
            continue
        cnt += 1
        l = least_aux(p, q)
        assert l is not None
        if l > mx:
            mx, mp = l, p
        slack = log10_bound(q, p) - math.log10(l)
        assert slack > 0 and log10_simple(q, p) >= log10_bound(q, p) - 1e-9
        worst_slack = min(worst_slack, slack)
    print("q=%2d: %6d primes; max l_q(p) = %d at p = %d; ratio l/(q^2 log p) = %.2f; min log10(bound/l) = %.1f"
          % (q, cnt, mx, mp, mx / (q * q * math.log(mp)), worst_slack))
print("Part 2: adversarial p = 1 + m*prod(first k primes l = 1 mod 12q) (p prime), so p is a q-th power mod each")
for q in [2, 3, 5, 7]:
    ls = [l for l in range(12 * q + 1, 10**6, 12 * q) if isprime(l)]
    for k in [2, 4, 8, 12]:
        P = 1
        for l in ls[:k]:
            P *= l
        m = 1
        while True:
            p = 1 + m * P
            if (6 * q) % p != 0 and isprime(p):
                break
            m += 1
        l = least_aux(p, q)
        print("q=%d k=%2d: p ~ 10^%.1f, l_q(p) = %d (> l_k = %d), q^2 log p = %.0f, log10 bound = %.1f, log10 simple = %.1f"
              % (q, k, math.log10(p), l, ls[k - 1], q * q * math.log(p), log10_bound(q, p), log10_simple(q, p)))
        assert l > ls[k - 1] and math.log10(l) < log10_bound(q, p)
print("ALL CHECKS PASSED")
