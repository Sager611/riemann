# q.py: query alpdep's name table. usage: q.py REGEX [maxhits]
# prints: dist  kind  fullname  file:start-end   (dist = alpdep BFS distance from the roots, -1 = off the cone)
import sys, re, numpy as np
R = '/Users/adriansager/playground/riemann/.scratch/'
rx = re.compile(sys.argv[1]); mx = int(sys.argv[2]) if len(sys.argv) > 2 else 200
dist = np.load(R + 'alpdep/dist.npy', mmap_mode='r')
hits = []
with open(R + 'alpdep/names.tsv') as fh:
    for ln in fh:
        p = ln.rstrip('\n').split('\t')
        if rx.search(p[1]):
            hits.append(p)
            if len(hits) >= mx: break
need = {int(p[3]) for p in hits}
fn = {}
with open(R + 'alpcore/files.txt') as fh:
    for i, ln in enumerate(fh):
        if i in need: fn[i] = ln.strip()
for p in hits:
    print(int(dist[int(p[0])]), p[2], p[1], '%s:%s-%s' % (fn.get(int(p[3]), '?'), p[4], p[5]))
