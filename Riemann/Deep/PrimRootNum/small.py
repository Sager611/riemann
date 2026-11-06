#!/usr/bin/env python3
"""
small.py -- exact least prime primitive roots g*(p) for all primes 3 <= p <= N (default 10^7), and adversarial large p.
Checks g*(p) <= (c log p)^E for the constants given on the command line, records the record values of g*(p).
Run: uv run --no-project --with numpy python3 small.py N c E      (one thread, ~0.2 GB for N = 10^7)
"""
import sys, math, time
import numpy as np

def spf_table(N):
    spf = np.zeros(N + 1, dtype=np.int32)
    spf[1] = 1
    for i in range(2, int(N ** 0.5) + 1):
        if spf[i] == 0:
            spf[i] = i
            blk = spf[i * i::i]
            mask = blk == 0
            blk[mask] = i
            spf[i * i::i] = blk
    rest = np.nonzero(spf == 0)[0]
    spf[rest] = rest
    spf[0] = 0
    return spf

def distinct_primes(n, spf):
    out = []
    while n > 1:
        l = int(spf[n]); out.append(l)
        while n % l == 0:
            n //= l
    return out

SMALLP = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97]

def gstar(p, ls, primes_iter=None):
    """least prime r (r != p) that is a primitive root mod p; ls = distinct primes of p - 1"""
    exps = [(p - 1) // l for l in ls]
    cand = primes_iter if primes_iter is not None else SMALLP
    for r in cand:
        if r == p:
            continue
        if all(pow(r, e, p) != 1 for e in exps):
            return r
    return None

def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0]

if __name__ == "__main__":
    N = int(float(sys.argv[1])) if len(sys.argv) > 1 else 10 ** 7
    c = float(sys.argv[2]) if len(sys.argv) > 2 else None
    E = float(sys.argv[3]) if len(sys.argv) > 3 else None
    t0 = time.time()
    spf = spf_table(N)
    pr = primes_upto(N)
    small_list = [int(x) for x in pr[:2000]]
    rec = 0; records = []; worst_ratio = 0.0; worst = None; fails = 0; cnt = 0
    maxw = 0
    for p in pr:
        p = int(p)
        if p < 3:
            continue
        ls = distinct_primes(p - 1, spf)
        maxw = max(maxw, len(ls))
        g = gstar(p, ls, small_list)
        if g is None:
            raise SystemExit("no prime primitive root found below %d for p = %d" % (small_list[-1], p))
        cnt += 1
        if g > rec:
            rec = g; records.append((p, g))
        if c is not None:
            # g <= (c log p)^E  <=>  log g <= E log(c log p); compare with a safety margin 1e-9
            lhs = math.log(g); rhs = E * math.log(c * math.log(p))
            if lhs > rhs - 1e-9:
                fails += 1
                if fails <= 10:
                    print("FAIL p=%d g*=%d bound=%.4g" % (p, g, math.exp(rhs)))
            ratio = lhs / rhs if rhs > 0 else float("inf")
            if ratio > worst_ratio:
                worst_ratio = ratio; worst = (p, g)
    print("N=%d primes checked=%d max omega(p-1)=%d time %.1fs" % (N, cnt, maxw, time.time() - t0))
    print("records of g*(p): " + ", ".join("%d:%d" % r for r in records))
    if c is not None:
        print("bound (%.4f log p)^%.4f: fails=%d, worst log g*/log bound = %.4f at p=%d (g*=%d)" % (c, E, fails, worst_ratio, worst[0], worst[1]))
