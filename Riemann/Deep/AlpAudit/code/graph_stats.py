#!/usr/bin/env python3
# Cross-check of shosonoda's import-graph figures from our rg dump of '^import' lines:
# internal (project→project) import count, acyclicity, longest import chain from Main.
import sys, collections
DUMP = sys.argv[1]
ids = {}; names = []
def mid(n):
    i = ids.get(n)
    if i is None:
        i = len(names); ids[n] = i; names.append(n)
    return i
raw = []
proj = set()
with open(DUMP, encoding='utf-8') as f:
    for line in f:
        path, _, rest = line.partition(':import ')
        a = mid(path[4:-5].replace('/', '.')); proj.add(a)
        for d in rest.split('--')[0].split():
            raw.append((a, mid(d)))
internal = [(a, b) for a, b in raw if b in proj]
dedup = set(internal)
print('edges total', len(raw), 'internal', len(internal), 'internal distinct', len(dedup), 'external', len(raw) - len(internal))
succ = collections.defaultdict(list); indeg = collections.Counter()
for a, b in dedup:
    succ[a].append(b); indeg[b] += 1
# Kahn topological order on the project graph (edges a -> b mean a imports b)
q = collections.deque(n for n in proj if indeg[n] == 0)
order = []
indeg2 = dict(indeg)
while q:
    n = q.popleft(); order.append(n)
    for b in succ.get(n, ()):
        indeg2[b] -= 1
        if indeg2[b] == 0: q.append(b)
print('acyclic', len(order) == len(proj), 'roots (imported by nobody)', sum(1 for n in proj if indeg[n] == 0), [names[n] for n in proj if indeg[n] == 0][:5])
depth = {}
for n in reversed(order):
    depth[n] = 1 + max((depth[b] for b in succ.get(n, ())), default=0)
m = ids['Main']
print('longest import chain from Main (modules, project only):', depth[m], '=> edges', depth[m] - 1)
