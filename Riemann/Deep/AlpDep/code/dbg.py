#!/usr/bin/env python3
# alpdep debug: for given full names (or suffixes), print id, dist, file:line and the names of its out-edges.
import sys, numpy as np
OUT = sys.argv[1]; qs = sys.argv[2].split(','); mx = int(sys.argv[3]) if len(sys.argv) > 3 else 60
dist = np.load(OUT + '/dist.npy'); eoff = np.load(OUT + '/eoff.npy'); elen = np.load(OUT + '/elen.npy')
store = np.fromfile(OUT + '/estore.bin', dtype=np.int32)
hit = {}
with open(OUT + '/names.tsv') as fh:
    for ln in fh:
        p = ln.rstrip('\n').split('\t')
        for q in qs:
            if p[1] == q or p[1].endswith('.' + q): hit.setdefault(q, []).append(p)
need = set()
for q, ps in hit.items():
    for p in ps:
        i = int(p[0])
        if eoff[i] >= 0: need.update(store[eoff[i]:eoff[i] + elen[i]].tolist())
nm = {}
with open(OUT + '/names.tsv') as fh:
    for ln in fh:
        p = ln.split('\t', 3)
        if int(p[0]) in need: nm[int(p[0])] = p[1]
for q in qs:
    for p in hit.get(q, []):
        i = int(p[0])
        print('%s id=%d dist=%d kind=%s f%s:%s edges=%s' % (p[1], i, dist[i], p[2], p[3], p[4], 'unprocessed' if eoff[i] < 0 else elen[i]))
        if eoff[i] >= 0:
            print('   ->', ' '.join(nm.get(j, str(j)) + '(%d)' % dist[j] for j in store[eoff[i]:eoff[i] + elen[i]].tolist()[:mx]))
