#!/usr/bin/env python3
# Textual import cone of src/Main.lean from an rg dump of '^import' lines (path:import M1 M2 ...).
import sys, collections, os
DUMP, SRC, OUT = sys.argv[1], sys.argv[2], sys.argv[3]
ids = {}; names = []
def mid(n):
    i = ids.get(n)
    if i is None:
        i = len(names); ids[n] = i; names.append(n)
    return i
edges = collections.defaultdict(list)
proj = set()
with open(DUMP, encoding='utf-8') as f:
    for line in f:
        path, _, rest = line.partition(':import ')
        m = path[4:-5].replace('/', '.')  # strip 'src/' and '.lean'
        a = mid(m); proj.add(a)
        for d in rest.split('--')[0].split():
            edges[a].append(mid(d))
# project files with no imports at all still count as project modules
start = mid('Main')
seen = set([start]); q = collections.deque([start])
while q:
    a = q.popleft()
    for b in edges.get(a, ()):
        if b not in seen:
            seen.add(b); q.append(b)
cone_proj = sorted(names[i] for i in seen if os.path.exists(os.path.join(SRC, *names[i].split('.')) + '.lean'))
cps = set(cone_proj)
cone_ext = sorted(names[i] for i in seen if names[i] not in cps)
with open(os.path.join(OUT, 'cone_files.txt'), 'w') as f: f.write('\n'.join(cone_proj) + '\n')
with open(os.path.join(OUT, 'cone_external_imports.txt'), 'w') as f: f.write('\n'.join(cone_ext) + '\n')
tops = collections.Counter(x.split('.')[0] for x in cone_proj)
nb = sum(os.path.getsize(os.path.join(SRC, *x.split('.')) + '.lean') for x in cone_proj)
ne = sum(len(edges.get(ids[x], ())) for x in cone_proj)
print('project modules with imports:', len(proj), 'all module names seen:', len(names))
print('cone project modules:', len(cone_proj), 'bytes:', nb, 'import edges inside cone files:', ne)
print('cone by top dir:', dict(tops))
print('external modules in cone:', len(cone_ext), collections.Counter(x.split('.')[0] for x in cone_ext))
