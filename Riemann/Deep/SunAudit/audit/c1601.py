import sys, itertools
sys.path.insert(0, '.')
from lib import *
from math import gcd
def lcm(a, b): return a * b // gcd(a, b)

def check_group(mods, node_limit=2*10**7):
    els = list(itertools.product(*[range(m) for m in mods]))
    idx = {e: i for i, e in enumerate(els)}
    n = 1
    for m in mods: n = lcm(n, m)
    def mul(k, x): return tuple((k * xi) % m for xi, m in zip(x, mods))
    def add(x, y): return tuple((a + b) % m for a, b, m in zip(x, y, mods))
    zero = els[0]
    allowed = {s: [x for x in els if mul(s, x) != zero] for s in range(1, n)}
    addtab = [[idx[add(x, y)] for y in els] for x in els]
    nodes = [0]
    found = []
    # DFS over s=1..n-1 ; Sigma as frozenset of indices
    def dfs(s, Sig, seq):
        nodes[0] += 1
        if nodes[0] > node_limit: raise TimeoutError
        if s == n:
            if len(Sig) < n:
                found.append(list(seq)); return True
            return False
        for x in allowed[s]:
            xi = idx[x]
            S2 = Sig | {addtab[y][xi] for y in Sig}
            if len(S2) >= n: continue
            seq.append(x)
            if dfs(s + 1, S2, seq): return True
            seq.pop()
        return False
    try:
        dfs(1, frozenset([0]), [])
        status = 'exhaustive'
    except TimeoutError:
        status = 'node-limit'
    return n, status, nodes[0], found

for mods in [[2,2],[3,3],[2,4],[2,2,4],[4,4],[3,3,3],[2,6],[5,5],[3,6],[2,2,6],[2,8],[3,9],[6,6],[2,10],[4,8],[2,12],[7,7],[2,14],[5,10],[3,12],[4,12],[2,2,2,4],[2,16],[3,15],[2,18]]:
    n, status, nodes, found = check_group(mods)
    print(mods, 'exp', n, status, nodes, 'COUNTEREXAMPLE' if found else 'none', found[:1], flush=True)
