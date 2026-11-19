#!/usr/bin/env python3
"""mc_indep.py (unitaudit) -- independent rerun of UnitCover mc.py's Monte Carlo (not part of the proof): least N such
that N uniform random primes from (1e6, 2e6] (distinct) cover every 2 <= q <= Q; own numpy closure, own seed."""
import numpy as np, math, random, statistics
s = np.ones(2 * 10**6 + 1, bool); s[:2] = False
for i in range(2, int((2 * 10**6) ** .5) + 1):
    if s[i]: s[i * i::i] = False
PR = [int(p) for p in np.nonzero(s)[0] if p > 10**6]
rng = random.Random(4242)
for Q, trials in ((100, 40), (300, 30), (1000, 12)):
    units = {q: np.array([math.gcd(x, q) == 1 for x in range(q)]) for q in range(3, Q + 1)}
    ar = {q: np.arange(q, dtype=np.int64) for q in range(3, Q + 1)}
    Ns = []
    for _ in range(trials):
        seq = rng.sample(PR, 60)
        S = {q: np.zeros(q, bool) for q in range(3, Q + 1)}
        for q in S: S[q][1] = True
        open_q = set(S)
        for n, p in enumerate(seq, 1):
            for q in list(open_q):
                S[q] = S[q] | S[q][(ar[q] * pow(p, -1, q)) % q]
                if S[q].sum() == units[q].sum(): open_q.discard(q)
            if not open_q: Ns.append(n); break
    print("Q=%d (log2 Q = %.2f): least N over %d draws: min %d, median %g, max %d" % (Q, math.log2(Q), trials, min(Ns), statistics.median(Ns), max(Ns)), flush=True)
