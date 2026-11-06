#!/usr/bin/env python3
"""mc.py -- Monte Carlo sanity check (not part of the proof).
Draw iid primes from (lo, hi] (log-weighted, as the law nu), and record the least N for which the first N distinct draws
cover every 2 <= q <= Q (subset products of the drawn primes not dividing q).  Compares the typical N with log2(Q)
and with the certified bound N(Q) of cert.py.  Run: python3 mc.py  (pure Python; < 2 min)."""
import random, math, sys, bisect

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]

def covered_upto(seq, Q):
    """least prefix length m of seq covering all 2 <= q <= Q, or None"""
    worst = 0
    for q in range(2, Q + 1):
        units = sum(1 for x in range(1, q) if math.gcd(x, q) == 1)
        if units == 1: continue
        S = bytearray(q); S[1] = 1; mem = [1]; m = None
        for i, p in enumerate(seq):
            if q % p == 0: continue
            pm = p % q; new = []
            for x in mem:
                zz = (x * pm) % q
                if not S[zz]: S[zz] = 1; new.append(zz)
            mem += new
            if len(mem) == units: m = i + 1; break
        if m is None: return None
        worst = max(worst, m)
    return worst

if __name__ == "__main__":
    random.seed(20261106)
    P = primes_upto(2 * 10 ** 6)
    pool = [p for p in P if p > 10 ** 6]
    w = [math.log(p) for p in pool]; cw = []
    t = 0
    for x in w: t += x; cw.append(t)
    for Q, trials in ((100, 40), (300, 30), (1000, 12)):
        res = []
        for _ in range(trials):
            seq = []
            while len(seq) < 80:
                p = pool[bisect.bisect_left(cw, random.random() * t)]
                if p not in seq: seq.append(p)
            res.append(covered_upto(seq, Q))
        res.sort()
        print("Q=%d (log2 Q = %.2f): least N over %d random draws from (1e6, 2e6]: min %d, median %d, max %d"
              % (Q, math.log2(Q), trials, res[0], res[len(res) // 2], res[-1]), flush=True)
