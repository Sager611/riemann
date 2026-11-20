#!/usr/bin/env python3
# alpaut lookup: lk.py REGEX [REGEX...] -> full name, kind, file:start-end, alpdep dist (-1 = off the cone)
# regex is matched (re.search) against the full declaration name; anchor with $ for suffixes.
import sys, re, numpy as np
S = '/Users/adriansager/playground/riemann/.scratch/'
pats = [re.compile(p) for p in sys.argv[1:]]
dist = np.load(S + 'alpdep/dist.npy', mmap_mode='r')
hits = []
with open(S + 'alpdep/names.tsv') as fh:
    for ln in fh:
        p = ln.rstrip('\n').split('\t')
        if any(q.search(p[1]) for q in pats): hits.append(p)
need = {int(h[3]) for h in hits}
files = {}
with open(S + 'alpcore/files.txt') as fh:
    for i, ln in enumerate(fh):
        if i in need: files[i] = ln.strip()[2:]
for h in hits[:400]:
    print(h[1], h[2], '%s:%s-%s' % (files.get(int(h[3]), '?'), h[4], h[5]), 'dist=%d' % int(dist[int(h[0])]))
print('#hits', len(hits))
