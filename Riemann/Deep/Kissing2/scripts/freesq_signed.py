#!/usr/bin/env python3
"""freesq_signed.py -- like freesq.py, but every cap of a candidate square is an individual signed point:
each of the 4 edges independently takes 0, 1 (t) or 2 (+-t) caps, any tails; constraints inside the square:
a body carries only t or {t,-t}; a u-edge cap and a v-edge cap (they always share a vertex, b.b' = 16) need
t.t' <= 0; parallel edges (b.b' = 8) are free.  Reports the best net gain caps - 4 over all candidate squares.
Usage: python3 freesq_signed.py config.npz
"""
import sys, itertools, collections, random, time
import numpy as np
from decode_ps import golay_from_config, leech_min, hnf, canon, cl, load_points

def main(path):
    t0 = time.time()
    P = load_points(path)
    if (P * P).sum(1)[0] == 128: P //= 2
    n, d = P.shape; k = d - 24
    B, T = P[:, :24], P[:, 24:]; bn = (B * B).sum(1)
    U = leech_min(golay_from_config(B[bn == 32]))
    key = {tuple(x): i for i, x in enumerate(U.tolist())}
    inE = np.zeros(len(U), dtype=bool)
    for x in B[bn == 32].tolist(): inE[key[tuple(x)]] = True
    random.seed(1); H2 = 2 * hnf(U[random.sample(range(len(U)), 400)])
    cls = {tuple(r): i for i, r in enumerate(canon(U, H2).tolist())}
    ci = np.flatnonzero(bn == 24); Bc = B[ci]
    w = [cls[tuple(r)] for r in canon(2 * Bc, H2).tolist()]
    ends = np.array([(key[tuple(((2 * b + U[j]) // 2).tolist())], key[tuple(((2 * b - U[j]) // 2).tolist())])
                     for b, j in zip(Bc, w)])
    adj = collections.defaultdict(set)
    for x, y in ends: adj[x].add(y); adj[y].add(x)
    pairs = set()
    for a, nb in adj.items():
        for b1, b2 in itertools.combinations(sorted(nb), 2):
            d1, d2 = U[a] - U[b1], U[a] - U[b2]
            if d1 @ d2 == 0: pairs.add(tuple(sorted([cl(d1), cl(d2)])))
    tails = np.array([t for t in itertools.product([-2, 0, 2], repeat=k) if sum(1 for c in t if c) == 2])
    TT = tails @ tails.T
    neg = {i: int(np.flatnonzero((tails == -tails[i]).all(1))[0]) for i in range(len(tails))}
    X = np.flatnonzero(bn < 32); Xb, Xt = B[X], T[X]; XT = Xt @ tails.T
    best = (-99, None); hist = collections.Counter(); nchk = 0
    for d1, d2 in pairs:
        u, v = np.array(d1), np.array(d2)
        Zi = np.flatnonzero(((U @ u) == 16) & ((U @ v) == 16))
        Mu, Mv = Xb @ u, Xb @ v
        for f in 2 * U[Zi] - u - v:
            vs = [key[tuple(((f + s1 * u + s2 * v) // 2).tolist())] for s1 in (1, -1) for s2 in (1, -1)]
            if not inE[vs].all(): continue
            nchk += 1
            Mf = Xb @ f
            ok = []
            for m2 in (Mf + Mv, Mf - Mv, Mf + Mu, Mf - Mu):          # 2b.x for u-edges (2) then v-edges (2)
                rows = np.flatnonzero(m2 > 16)                    # conflict needs 2b.x + 2t.tx > 32 with t.tx <= 8
                ok.append(set(np.flatnonzero(~((m2[rows, None] + 2 * XT[rows]) > 32).any(0)).tolist()))
            def options(s):
                o = [()] + [(t,) for t in s] + [(t, neg[t]) for t in s if neg[t] in s and t < neg[t]]
                return o
            ou = [options(ok[0]), options(ok[1])]; ov = [options(ok[2]), options(ok[3])]
            top = 0
            for o1 in ou[0]:
                for o2 in ou[1]:
                    A = o1 + o2
                    tot = len(A)
                    for s in (ok[2], ok[3]):
                        good = [t for t in s if all(TT[t, a] <= 0 for a in A)]
                        g = set(good)
                        tot += 2 if any(neg[t] in g for t in good) else (1 if good else 0)
                    top = max(top, tot)
            hist[top] += 1
            if top - 4 > best[0]: best = (top - 4, (d1, d2, tuple(f.tolist())))
    print(path, 'squares checked', nchk, 'max caps per square histogram', sorted(hist.items()),
          'best net gain', best[0], '(%.0fs)' % (time.time() - t0))

if __name__ == '__main__':
    main(sys.argv[1])
