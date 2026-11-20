#!/usr/bin/env python3
# alpdep: analysis of the lazy BFS result: shortest name paths from the roots to each target, and per target the
# "path cone" (reached decls from which the target is reachable inside the reached subgraph), with the decls
# closest to the target listed (distance-to-target <= K).
import sys, numpy as np
OUT = sys.argv[1]; targets = [t for t in sys.argv[2].split(',') if t]; K = int(sys.argv[3]); conet = set(sys.argv[4].split(','))
dist = np.load(OUT + '/dist.npy'); par = np.load(OUT + '/par.npy'); eoff = np.load(OUT + '/eoff.npy'); elen = np.load(OUT + '/elen.npy')
store = np.fromfile(OUT + '/estore.bin', dtype=np.int32); N = len(dist)
name_id = {}
with open(OUT + '/names.tsv') as fh:
    for ln in fh:
        p = ln.split('\t', 2)
        if p[1] in targets: name_id.setdefault(p[1], []).append(int(p[0]))
R = np.nonzero(dist >= 0)[0]
print('reached', len(R), 'of', N, '; max dist', int(dist.max()))
lens = elen[R].astype(np.int64); st = eoff[R]
tot = int(lens.sum()); rel = np.arange(tot) - np.repeat(np.cumsum(lens) - lens, lens)
dst = store[np.repeat(st, lens) + rel]; src = np.repeat(R, lens).astype(np.int32)
print('edges among reached', tot)
need = set(); out = []
for t in targets:
    ids = name_id.get(t, [])
    best = None
    for i in ids:
        if dist[i] >= 0 and (best is None or dist[i] < dist[best]): best = i
    if best is None:
        out.append(('NOT', t, ids, None, None)); continue
    p = [int(best)]
    while par[p[-1]] >= 0: p.append(int(par[p[-1]]))
    p = p[::-1]; need.update(p)
    cone = None
    if t in conet:
        bd = np.full(N, -1, dtype=np.int32); bd[best] = 0; F = np.array([best]); d = 0
        inF = np.zeros(N, dtype=bool)
        while len(F) and d < 200:
            inF[:] = False; inF[F] = True
            s = src[inF[dst]]; s = np.unique(s[bd[s] < 0])
            d += 1; bd[s] = d; F = s
        cone = bd
        near = np.nonzero((bd >= 0) & (bd <= K))[0]; need.update(near.tolist())
    out.append(('YES', t, ids, p, cone))
info = {}
with open(OUT + '/names.tsv') as fh:
    for ln in fh:
        q = ln.rstrip('\n').split('\t')
        if int(q[0]) in need: info[int(q[0])] = q
nf = set(int(info[i][3]) for i in info); files = {}
with open('/Users/adriansager/playground/riemann/.scratch/alpcore/files.txt') as fh:
    for k, ln in enumerate(fh):
        if k in nf: files[k] = ln.strip()[2:]
def show(i):
    q = info[i]; return '%s [%s] %s:%s' % (q[1], q[2], files[int(q[3])], q[4])
for (st_, t, ids, p, cone) in out:
    if st_ == 'NOT':
        print('TARGET', t, 'ids', ids[:4], 'NOT REACHED'); continue
    print('TARGET', t, 'dist', len(p) - 1, 'path:')
    for i in p: print('   ', show(i))
    if cone is not None:
        c = np.nonzero(cone >= 0)[0]
        print('  path-cone size', len(c), '(reached decls that reach it); roots in cone:', [int(x) for x in c if dist[x] == 0])
        near = sorted(np.nonzero((cone >= 1) & (cone <= K))[0].tolist(), key=lambda x: (cone[x], x))
        for i in near: print('    d2t=%d dist=%d %s' % (cone[i], dist[i], show(i)))
