#!/usr/bin/env python3
# alprows: look up declarations by exact name or dotted suffix in alpdep's names.tsv; print kind, file:line range, cone dist.
# usage: lookup.py name1 name2 ...   (or -f file with one name per line)
import sys, numpy as np
S = '/Users/adriansager/playground/riemann/.scratch/'
args = sys.argv[1:]
if args and args[0] == '-f':
    args = [l.strip() for l in open(args[1]) if l.strip()]
want = set(args)
dist = np.load(S + 'alpdep/dist.npy', mmap_mode='r')
hits = []
with open(S + 'alpdep/names.tsv') as fh:
    for ln in fh:
        q = ln.rstrip('\n').split('\t')
        n = q[1]
        last = n.rsplit('.', 1)[-1]
        if n in want or last in want or any(n.endswith('.' + w) for w in want if '.' in w):
            hits.append(q)
nf = set(int(q[3]) for q in hits); files = {}
with open(S + 'alpcore/files.txt') as fh:
    for k, l in enumerate(fh):
        if k in nf: files[k] = l.strip()[2:]
for q in sorted(hits, key=lambda q: q[1]):
    i = int(q[0])
    print('%s [%s] dist=%d %s:%s-%s' % (q[1], q[2], int(dist[i]), files[int(q[3])], q[4], q[5]))

