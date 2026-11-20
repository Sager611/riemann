# path.py: shortest BFS path (alpdep's par.npy) from the roots to each target full name. usage: path.py NAME1,NAME2,...
import sys, numpy as np
R = '/Users/adriansager/playground/riemann/.scratch/'
targets = sys.argv[1].split(',')
dist = np.load(R + 'alpdep/dist.npy', mmap_mode='r'); par = np.load(R + 'alpdep/par.npy', mmap_mode='r')
tid = {}
with open(R + 'alpdep/names.tsv') as fh:
    for ln in fh:
        p = ln.split('\t', 2)
        if p[1] in targets: tid[p[1]] = int(p[0])
chains = {}
need = set()
for t, i in tid.items():
    ch = [i]
    while dist[ch[-1]] > 0 and len(ch) < 80:
        ch.append(int(par[ch[-1]]))
    chains[t] = ch; need.update(ch)
nm = {}
with open(R + 'alpdep/names.tsv') as fh:
    for ln in fh:
        p = ln.split('\t', 2)
        if int(p[0]) in need: nm[int(p[0])] = p[1]
for t in targets:
    if t not in tid: print(t, 'NOT FOUND'); continue
    ch = chains[t]
    print(t, 'dist', int(dist[tid[t]]))
    print('  ' + ' < '.join(nm[i].split('.')[-1] + '[' + str(int(dist[i])) + ']' for i in ch))
