import sys
from math import gcd
sys.path.insert(0, ".scratch/prefixlean")
from words import rev, perm_of, cycles, bfs_pair, lcm

def search(n, m, k, maxlen=15):
    gens = [rev(n, n), rev(n, m), rev(n, k)]
    best = None
    # DFS over reduced words (no immediate repeats)
    stack = [((), tuple(range(n)))]
    while stack:
        w, p = stack.pop()
        if w:
            lens = [len(c) for c in cycles(p)]
            if lens.count(2) == 1 and all(l % 2 == 1 for l in lens if l != 2):
                e = 1
                for l in lens:
                    if l != 2: e = lcm(e, l)
                cost = len(w) * e
                if best is None or cost < best[0]:
                    best = (cost, list(w), e)
        if len(w) < maxlen:
            for s in range(3):
                if w and w[-1] == s: continue
                q = tuple(gens[s][p[x]] for x in range(n))  # apply s after p
                stack.append((w + (s,), q))
    return best

for (n, m, k) in [(23,18,4), (23,16,6), (23,12,10)]:
    cost, hseq, e = search(n, m, k)
    gens = [rev(n, n), rev(n, m), rev(n, k)]
    tau = perm_of(n, gens, hseq * e)
    moved = [x for x in range(n) if tau[x] != x]
    t1, t2 = moved
    a = k + 1
    c = bfs_pair(n, gens, (t1, t2), (0, a))
    seq = c[::-1] + hseq * e + c
    p = perm_of(n, gens, seq)
    target = tuple(a if x == 0 else 0 if x == a else x for x in range(n))
    assert p == target
    chunks = []
    for i in range(0, len(seq), 20):
        ch = seq[i:i+20]
        chunks.append((len(ch), sum(d * 3**j for j, d in enumerate(ch))))
    open(".scratch/prefixlean/chunks_%d_%d_%d.txt" % (n, m, k), "w").write(", ".join("(%d, %d)" % t for t in chunks))
    print(n, m, k, "cost", cost, "h", len(hseq), "e", e, "c", len(c), "L", len(seq))
