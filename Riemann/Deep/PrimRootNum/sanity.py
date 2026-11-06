#!/usr/bin/env python3
"""
sanity.py -- independent checks of the two sieve lemmas (not part of the proof; DERIVATION.md sec 3).
 (S1) For the DFS weight sets D(y) of sieve.py: sum_{d | n, d in D} mu(d) <= [n = 1] for EVERY squarefree n | P_w
      (brute force over all 2^w subsets, w <= 13, many y), and T, C agree with sieve.py.
 (S2) Shifted-sieve main term: for random sets of w primes l_1 < ... < l_w (l_j >= p_j),
      sum_{d | P_w} lambda_d / d'  >=  V_Ups * T_w / V_w,   d' = image of d under p_j -> l_j.
Run: python3 sanity.py
"""
import random, itertools
from fractions import Fraction as Fr
from sieve import primes_first, enumerate_weights, cumulative, ONE

def dset(y, P):
    """explicit D(y) as bitmasks over indices of P (same rule as sieve.py)"""
    out = [0]                                  # d = 1
    def rec(mask, prod, li, r):
        for j in range(li):
            q = P[j]
            pos = r + 1                        # position of the new prime
            if pos % 2 == 0 and not (prod * q ** 3 < y):
                continue
            m2 = mask | (1 << j); out.append(m2)
            rec(m2, prod * q, j, r + 1)
    rec(0, 1, len(P), 0)
    return out

def mu_of(mask):
    return -1 if bin(mask).count("1") % 2 else 1

def check_S1(w, y):
    P = primes_first(w)
    D = dset(y, P)
    bad = 0
    for n in range(1 << w):
        s = sum(mu_of(d) for d in D if d & ~n == 0)
        if s > (1 if n == 0 else 0):
            bad += 1
    # compare T, C with sieve.py
    T = sum(Fr(mu_of(d), eval_prod(d, P)) for d in D); C = len(D)
    Tq, Cq = enumerate_weights(y, P); Tc, Cc = cumulative(Tq, Cq)
    assert Cc[-1] == C, (w, y, Cc[-1], C)
    assert Tc[-1] <= T * ONE <= Tc[-1] + 4 * C, (w, y)      # exact lower bound within rounding
    return bad, C, float(T)

def eval_prod(mask, P):
    p = 1
    for j in range(len(P)):
        if mask >> j & 1:
            p *= P[j]
    return p

def check_S2(w, y, trials=300, seed=1):
    P = primes_first(w)
    D = dset(y, P)
    Vw = Fr(1)
    for q in P:
        Vw *= Fr(q - 1, q)
    T = sum(Fr(mu_of(d), eval_prod(d, P)) for d in D)
    rnd = random.Random(seed)
    pool = primes_first(200)
    worst = None
    for _ in range(trials):
        # random increasing l_j >= p_j
        while True:
            L = sorted(rnd.sample(pool, w))
            if all(L[j] >= P[j] for j in range(w)):
                break
        Vu = Fr(1)
        for l in L:
            Vu *= Fr(l - 1, l)
        lhs = sum(Fr(mu_of(d), eval_prod(d, L)) for d in D)
        rhs = Vu * T / Vw
        r = lhs - rhs
        if worst is None or r < worst:
            worst = r
    return worst

def check_S3(w, y):
    """plain DFS listing every element of D(y) (no bulk leaves), exact rational T and count C, against sieve.py"""
    P = primes_first(w)
    T = Fr(0); C = 0
    stack = [(1, len(P), 0)]
    while stack:
        prod, li, r = stack.pop()
        T += Fr(1 if r % 2 == 0 else -1, prod); C += 1
        for j in range(li):
            q = P[j]
            if (r + 1) % 2 == 0 and not (prod * q ** 3 < y):
                continue
            stack.append((prod * q, j, r + 1))
    Tq, Cq = enumerate_weights(y, P); Tc, Cc = cumulative(Tq, Cq)
    ok = (Cc[-1] == C) and (Tc[-1] <= T * ONE <= Tc[-1] + 4 * C)
    return ok, C, float(T), Tc[-1] / ONE

if __name__ == "__main__":
    for (w, y) in [(20, 10 ** 4), (20, 10 ** 6), (50, 10 ** 4), (50, 10 ** 5), (100, 10 ** 5), (100, 3 * 10 ** 5), (200, 10 ** 5)]:
        ok, C, T, Ts = check_S3(w, y)
        print("S3 w=%3d y=%7d: plain DFS C=%d T=%.10f, sieve.py T_lo=%.10f  %s" % (w, y, C, T, Ts, "OK" if ok else "MISMATCH"), flush=True)
        assert ok
    tot_bad = 0
    for w in range(1, 14):
        for y in (3, 10, 30, 100, 300, 1000, 3000, 10 ** 4, 10 ** 5, 10 ** 6, 10 ** 9):
            bad, C, T = check_S1(w, y)
            tot_bad += bad
        print("S1 w=%2d: violations %d (last y: C=%d, T=%.5f)" % (w, tot_bad, C, T), flush=True)
    print("S1 total violations:", tot_bad)
    for w in (1, 2, 3, 4, 5, 6, 8, 10):
        for y in (10, 100, 1000, 10 ** 4, 10 ** 6):
            wst = check_S2(w, y, trials=200, seed=w * 1000 + y % 997)
            print("S2 w=%2d y=%7d: min(lhs - rhs) = %.3e  (>= 0 required)" % (w, y, float(wst)), flush=True)
            assert wst >= 0
    print("sanity OK")
