#!/usr/bin/env python3
"""tess_completion.py -- can PackingStar's square groups be completed to tesseracts?  (28D, about 16 min)

A full group (u,v; 16 orthogonal frame vectors f_i) holds 32 squares at centres +-f_i/2.  The squares at
+-f_i/2, +-f_j/2 are the 16 vertices of the tesseract {u, v, (f_i+f_j)/2, (f_i-f_j)/2}; its 16 missing edges
(directions (f_i +- f_j)/2, on the two remaining lines of the tail frame) add 32 caps and delete nothing.
For every group, pair (i,j) and line assignment this counts the existing caps/mid points that conflict.
Measured on PackingStar's D28_221864 (commit 477df34): 42480 candidates, every one conflicts with >= 378 points
spread over >= 5 other groups (median 1120).  Also: within one group two completions always collide
((u+v)/2 components give b.b' = 16 > 8), so at most one per group.
Input: https://raw.githubusercontent.com/CDM1619/PackingStar/477df34cf2a09bb9476b735a685b538a82da8c55/25D-31D/new_bounds_20261028/D28_221864_integer_coordinates.npz
Usage: python3 tess_completion.py D28_221864_integer_coordinates.npz
"""
import sys, collections, itertools, random, time
import numpy as np
from decode_ps import golay_from_config, leech_min, hnf, canon, cl, load_points

def main(path):
    t0 = time.time()
    P = load_points(path)
    if (P * P).sum(1)[0] == 128: P //= 2
    B, T = P[:, :24], P[:, 24:]; bn = (B * B).sum(1)
    U = leech_min(golay_from_config(B[bn == 32]))
    key = {tuple(x): i for i, x in enumerate(U.tolist())}
    random.seed(1); H2 = 2 * hnf(U[random.sample(range(len(U)), 400)])
    cls = {tuple(r): i for i, r in enumerate(canon(U, H2).tolist())}
    ci = np.flatnonzero(bn == 24); Bc = B[ci]
    w = [cls[tuple(r)] for r in canon(2 * Bc, H2).tolist()]
    ends = np.array([(key[tuple(((2 * b + U[j]) // 2).tolist())], key[tuple(((2 * b - U[j]) // 2).tolist())])
                     for b, j in zip(Bc, w)])
    edge = {}
    for r, (x, y) in enumerate(ends): edge.setdefault((min(x, y), max(x, y)), set()).add(cl(T[ci[r]]))
    adj = collections.defaultdict(list)
    for (x, y) in edge: adj[x].append(y); adj[y].append(x)
    seen = set(); groups = collections.defaultdict(list); vgroup = {}
    for s in list(adj):
        if s in seen or len(adj[s]) != 2: continue
        st = [s]; seen.add(s); vs = []
        while st:
            x = st.pop(); vs.append(x)
            for y in adj[x]:
                if y not in seen: seen.add(y); st.append(y)
        if len(vs) != 4 or any(len(adj[x]) != 2 for x in vs): continue
        ds = {}
        for a in vs:
            for b in adj[a]:
                if a < b: ds[cl(U[a] - U[b])] = next(iter(edge[(a, b)]))
        if len(ds) != 2: continue
        (d1, l1), (d2, l2) = sorted(ds.items())
        if np.dot(l1, l2) != 0: continue
        groups[(d1, l1, d2, l2)].append((vs, U[vs].sum(0)))
    full = [g for g, v in groups.items() if len(v) == 32]
    gid = {g: i for i, g in enumerate(full)}
    for g in full:
        for vs, _ in groups[g]:
            for x in vs: vgroup[x] = gid[g]
    owner = np.full(len(P), -1)
    for r, (x, y) in enumerate(ends): owner[ci[r]] = vgroup.get(x, -1)
    nq = np.flatnonzero(bn < 32); Xb, Xt, Xo = B[nq], T[nq], owner[nq]
    k = P.shape[1] - 24
    lines = sorted(set(cl(t) for t in itertools.product([-2, 0, 2], repeat=k) if sum(1 for c in t if c) == 2))
    vals = []; ngr = []
    for g in full:
        d1, l1, d2, l2 = g; u, v = np.array(d1), np.array(d2)
        rem = [np.array(l) for l in lines if np.dot(l, l1) == 0 and np.dot(l, l2) == 0]
        fk = np.array(list({cl(zs // 2): zs // 2 for _, zs in groups[g]}.values()))
        Mu, Mv, MF = Xb @ u, Xb @ v, Xb @ fk.T
        Tr = [np.abs(Xt @ r) for r in rem]
        for i in range(len(fk)):
            for j in range(i + 1, len(fk)):
                for asg in (0, 1):
                    bad = np.zeros(len(Xb), dtype=bool)
                    for wdir, tr in ((0, Tr[asg]), (1, Tr[1 - asg])):
                        sg = -1 if wdir == 0 else 1
                        for s1, s2, s3 in itertools.product((1, -1), repeat=3):
                            bad |= (2 * (s1 * Mu + s2 * Mv) + s3 * (MF[:, i] + sg * MF[:, j]) + 4 * tr) > 64
                    bad &= (Xo != gid[g])
                    vals.append(int(bad.sum())); ngr.append(len(set(Xo[bad].tolist()) - {-1}))
    print('full groups', len(full), 'candidates', len(vals), 'conflicts min', min(vals), 'median', sorted(vals)[len(vals) // 2],
          'groups hit min', min(ngr), '(%.0fs)' % (time.time() - t0))

if __name__ == '__main__':
    main(sys.argv[1])
