#!/usr/bin/env python3
"""grid.py -- choose, for every w <= Wmax, the sieve level y (from a fixed grid) minimising (C_w - 1)/T_w.
Writes logs/sieve_best.json: {w: [y, T_numerator(str), C]}.  Exact integers from sieve.py (unit 2^-64)."""
import sys, json, time
from sieve import primes_first, enumerate_weights, cumulative, ONE
Wmax = int(sys.argv[1]); kmax = int(sys.argv[2]) if len(sys.argv) > 2 else 82
P = primes_first(Wmax)
best = [None] * Wmax
t00 = time.time()
for k in range(4, kmax + 1):
    y = int(round(10 ** (k / 8)))
    t0 = time.time()
    Tq, Cq = enumerate_weights(y, P)
    T, C = cumulative(Tq, Cq)
    for i in range(Wmax):
        if T[i] <= 0:
            continue
        r = (C[i] - 1) / T[i]
        if best[i] is None or r < best[i][3]:
            best[i] = (y, T[i], C[i], r)
    print("k=%d y=%d %.1fs (total %.0fs)" % (k, y, time.time() - t0, time.time() - t00), file=sys.stderr, flush=True)
out = {"Wmax": Wmax, "pmax": P[-1], "K": 64, "best": [[b[0], str(b[1]), b[2]] for b in best]}
json.dump(out, open(sys.argv[3] if len(sys.argv) > 3 else "logs/sieve_best.json", "w"))
ys = sorted(set(b[0] for b in best))
print("distinct y: %d" % len(ys), ys, file=sys.stderr)
for w in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000):
    if w <= Wmax:
        b = best[w - 1]
        print("w=%d p_w=%d y=%d T=%.6g C=%d (C-1)/T=%.5g  /p_w^2=%.4f" % (w, P[w - 1], b[0], b[1] / ONE, b[2], b[3] * ONE, b[3] * ONE / P[w - 1] ** 2), file=sys.stderr)
