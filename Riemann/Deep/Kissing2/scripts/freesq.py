#!/usr/bin/env python3
"""freesq.py -- free SQUARES for a PackingStar-family configuration (k = 3..7).

A square = 4 Leech minimal vectors a, a-u, a-u-v, a-v (u, v orthogonal minimal vectors, a.u = a.v = 16) carrying
8 caps: the 4 edge midpoints, each with both signs of a tail line (2*D_k root), the two directions on two
ORTHOGONAL tail lines.  Adding it deletes the 4 equator vectors and adds 8 caps (+4).  Candidates: all squares
(44 per pair) for every direction pair (u, v) that already meets at a vertex of the configuration.  A candidate
is free if its 4 vertices are on the equator and its 8 caps have inner product <= 16 (norm-32 scale) with every
non-equator point.  Result on the current 27D-30D records: 0 free squares (REPORT.md).
Usage: python3 freesq.py config.npz
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
    lines = sorted(set(cl(t) for i in range(k) for j in range(i + 1, k) for si in (2, -2) for sj in (2, -2)
                       for t in [tuple(si if c == i else sj if c == j else 0 for c in range(k))]))
    L = np.array(lines); LL = L @ L.T
    orthl = [(a, b) for a in range(len(L)) for b in range(len(L)) if a != b and LL[a, b] == 0]
    X = np.flatnonzero(bn < 32); Xb, Xt = B[X], T[X]; XT = np.abs(Xt @ L.T)
    nfree = 0; nchk = 0
    for d1, d2 in pairs:
        u, v = np.array(d1), np.array(d2)
        Zi = np.flatnonzero(((U @ u) == 16) & ((U @ v) == 16))
        Mu, Mv = Xb @ u, Xb @ v
        for f in 2 * U[Zi] - u - v:
            vs = [key[tuple(((f + s1 * u + s2 * v) // 2).tolist())] for s1 in (1, -1) for s2 in (1, -1)]
            if not inE[vs].all(): continue
            nchk += 1
            Mf = Xb @ f
            mu = np.maximum(Mf + Mv, Mf - Mv); mv = np.maximum(Mf + Mu, Mf - Mu)   # 2 b.x for the 4 midpoints
            badu = ((mu[:, None] + 2 * XT) > 32).any(0); badv = ((mv[:, None] + 2 * XT) > 32).any(0)
            nfree += sum(1 for a, b in orthl if not badu[a] and not badv[b])
    print(path, 'direction pairs', len(pairs), 'squares with 4 equator vertices', nchk,
          'free (square, line-pair) options', nfree, '(%.0fs)' % (time.time() - t0))

if __name__ == '__main__':
    main(sys.argv[1])
