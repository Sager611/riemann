#!/usr/bin/env python3
# alpdep: forward BFS over the name-level graph from the archive roots; shortest paths to targets;
# per-target "path cone" (= forward-reachable nodes that can reach the target).
import sys, re, numpy as np
OUT = sys.argv[1]; roots = sys.argv[2].split(','); targets = sys.argv[3].split(',')
cnt = np.fromfile(OUT + '/ecount.bin', dtype=np.int32)
N = len(cnt)
off = np.zeros(N + 1, dtype=np.int64); np.cumsum(cnt, out=off[1:])
tg = np.memmap(OUT + '/etarget.bin', dtype=np.int32, mode='r')
assert len(tg) == off[-1], (len(tg), off[-1])
# name lookup: stream names.tsv for given full names
want = set(roots) | set(targets)
name_id = {}
with open(OUT + '/names.tsv') as fh:
    for ln in fh:
        p = ln.split('\t', 2)
        if p[1] in want: name_id.setdefault(p[1], []).append(int(p[0]))
missing = [n for n in want if n not in name_id]
print('N', N, 'edges', int(off[-1]), 'missing names', missing)
par = np.full(N, -2, dtype=np.int32); dist = np.full(N, -1, dtype=np.int32)
F = np.array(sorted(i for r in roots for i in name_id.get(r, [])), dtype=np.int64)
par[F] = -1; dist[F] = 0; d = 0
while len(F):
    lens = cnt[F].astype(np.int64); tot = int(lens.sum())
    if tot == 0: break
    st = np.repeat(off[F], lens); rel = np.arange(tot) - np.repeat(np.cumsum(lens) - lens, lens)
    nb = np.asarray(tg[st + rel]); src = np.repeat(F, lens)
    m = par[nb] == -2; nb = nb[m]; src = src[m]
    u, ix = np.unique(nb, return_index=True)
    par[u] = src[ix]; d += 1; dist[u] = d; F = u.astype(np.int64)
reach = dist >= 0
print('reachable', int(reach.sum()), 'max depth', d)
np.save(OUT + '/dist.npy', dist); np.save(OUT + '/par.npy', par)
# collect names for paths
need = set()
paths = {}
for t in targets:
    best = None
    for i in name_id.get(t, []):
        if dist[i] >= 0 and (best is None or dist[i] < dist[best]): best = i
    if best is None: paths[t] = None; continue
    p = [best]
    while par[p[-1]] >= 0: p.append(int(par[p[-1]]))
    paths[t] = p[::-1]; need.update(p)
info = {}
with open(OUT + '/names.tsv') as fh:
    for ln in fh:
        p = ln.rstrip('\n').split('\t')
        if int(p[0]) in need: info[int(p[0])] = p
files = {}
needf = set(int(info[i][3]) for i in info)
with open('/Users/adriansager/playground/riemann/.scratch/alpcore/files.txt') as fh:
    for k, ln in enumerate(fh):
        if k in needf: files[k] = ln.strip()[2:]
for t in targets:
    p = paths[t]
    ids = name_id.get(t, [])
    if p is None:
        print('TARGET', t, 'ids', ids[:5], 'NOT REACHABLE'); continue
    print('TARGET', t, 'dist', len(p) - 1)
    for i in p:
        q = info[i]; print('   %s [%s] %s:%s' % (q[1], q[2], files[int(q[3])], q[4]))
