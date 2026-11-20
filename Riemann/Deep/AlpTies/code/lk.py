#!/usr/bin/env python3
# usage: python3 lk.py max regex1 regex2 ...  -> per regex: dist kind file:line name (on-cone first)
import sys, re, numpy as np
mx = int(sys.argv[1]); rxs = [re.compile(a) for a in sys.argv[2:]]
D = np.load('/Users/adriansager/playground/riemann/.scratch/alpdep/dist.npy', mmap_mode='r')
hits = [[] for _ in rxs]
with open('/Users/adriansager/playground/riemann/.scratch/alpdep/names.tsv', encoding='utf-8') as f:
    for line in f:
        p = line.rstrip('\n').split('\t')
        if len(p) < 6: continue
        for k, rx in enumerate(rxs):
            if rx.search(p[1]):
                i = int(p[0]); hits[k].append((int(D[i]) if i < len(D) else -9, p[2], int(p[3]), p[4], p[1]))
need = set(h[2] for hs in hits for h in hs)
fl = {}
with open('/Users/adriansager/playground/riemann/.scratch/alpcore/files.txt', encoding='utf-8') as f:
    for k, line in enumerate(f):
        if k in need: fl[k] = line.strip()
for k, hs in enumerate(hits):
    hs.sort(key=lambda h: (h[0] < 0, h[0]))
    print('##', sys.argv[2+k], len(hs), 'hits;', sum(1 for h in hs if h[0] >= 0), 'on cone')
    for h in hs[:mx]:
        print(h[0], h[1], fl.get(h[2], h[2]) + ':' + h[3], h[4])
