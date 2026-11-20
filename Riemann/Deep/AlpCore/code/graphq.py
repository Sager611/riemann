#!/usr/bin/env python3
# Query the import graph built by importgraph.py.
# usage: graphq.py reach | path <mod> [from] | rdeps <mod> | deps <mod> | closure <mod> | grepmods <regex> [reachable-only]
import sys, array, re
from collections import deque
OUT = "/Users/adriansager/playground/riemann/.scratch/alpcore"
mods = [l.rstrip("\n") for l in open(OUT + "/mods.txt")]
idx = {m: i for i, m in enumerate(mods)}
raw = array.array("I")
with open(OUT + "/edges.bin", "rb") as fh:
    raw.frombytes(fh.read())
E = len(raw) // 2
S, D = raw[:E], raw[E:]
nl = array.array("I")
with open(OUT + "/nlines.bin", "rb") as fh: nl.frombytes(fh.read())
N = len(mods)
def adj(rev=False):
    a = [[] for _ in range(N)]
    for s, d in zip(S, D):
        if rev: a[d].append(s)
        else: a[s].append(d)
    return a
def bfs(a, src):
    par = {src: -1}; q = deque([src])
    while q:
        u = q.popleft()
        for v in a[u]:
            if v not in par: par[v] = u; q.append(v)
    return par
cmd = sys.argv[1]
if cmd == "reach":
    a = adj(); par = bfs(a, idx["Main"])
    from collections import Counter
    c = Counter(); L = Counter()
    for v in par:
        k = ".".join(mods[v].split(".")[:2]); c[k] += 1; L[k] += nl[v]
    print("reachable from Main:", len(par), "of", N, "lines", sum(nl[v] for v in par))
    for k, n in c.most_common(40): print(f"{n:7d} {L[k]:9d} {k}")
    with open(OUT + "/reach.txt", "w") as fo:
        for v in sorted(par): fo.write(mods[v] + "\n")
elif cmd == "path":
    a = adj(); src = idx[sys.argv[3]] if len(sys.argv) > 3 else idx["Main"]
    par = bfs(a, src); t = idx.get(sys.argv[2])
    if t is None or t not in par: print("unreachable", sys.argv[2]); sys.exit()
    p = []
    while t != -1: p.append(mods[t]); t = par[t]
    for m in reversed(p): print(m)
elif cmd in ("rdeps", "deps"):
    a = adj(rev=(cmd == "rdeps")); t = idx[sys.argv[2]]
    for v in a[t]: print(mods[v])
elif cmd == "closure":
    a = adj(); par = bfs(a, idx[sys.argv[2]])
    print(len(par), "modules", sum(nl[v] for v in par), "lines")
elif cmd == "grepmods":
    rx = re.compile(sys.argv[2], re.I)
    reach = None
    if len(sys.argv) > 3:
        reach = set(l.strip() for l in open(OUT + "/reach.txt"))
    for m in mods:
        if rx.search(m) and (reach is None or m in reach): print(m)
