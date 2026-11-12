import random, sys
from collections import deque
from math import gcd

def rev(n, i):
    return tuple((i - 1 - x) if x < i else x for x in range(n))

def apply_seq(n, gens, seq, x):
    for s in seq:
        x = gens[s][x]
    return x

def perm_of(n, gens, seq):
    return tuple(apply_seq(n, gens, seq, x) for x in range(n))

def cycles(p):
    n = len(p); seen = [False]*n; out = []
    for i in range(n):
        if not seen[i]:
            c = []; j = i
            while not seen[j]:
                seen[j] = True; c.append(j); j = p[j]
            out.append(c)
    return out

def lcm(a, b): return a*b//gcd(a, b)

def find_transposition(n, gens, rng, tries=100000):
    best = None
    for _ in range(tries):
        L = rng.randint(5, 60)
        seq = [rng.randrange(3) for _ in range(L)]
        p = perm_of(n, gens, seq)
        cs = cycles(p)
        lens = [len(c) for c in cs]
        if lens.count(2) == 1 and all(l % 2 == 1 for l in lens if l != 2):
            e = 1
            for l in lens:
                if l != 2: e = lcm(e, l)
            cost = L * e
            if best is None or cost < best[0]:
                best = (cost, seq, e)
    return best

def bfs_pair(n, gens, src, dst):
    # find seq c with c(src[0]) = dst[0], c(src[1]) = dst[1]
    start = src; prev = {start: None}; q = deque([start])
    while q:
        u = q.popleft()
        if u == dst:
            seq = []
            while prev[u] is not None:
                u0, s = prev[u]; seq.append(s); u = u0
            return seq[::-1]
        for s in range(3):
            v = (gens[s][u[0]], gens[s][u[1]])
            if v not in prev:
                prev[v] = (u, s); q.append(v)
    return None

def main(n, m, k, seed=1):
    rng = random.Random(seed)
    gens = [rev(n, n), rev(n, m), rev(n, k)]
    a = k + 1
    cost, hseq, e = find_transposition(n, gens, rng)
    tau = perm_of(n, gens, hseq * e)
    moved = [x for x in range(n) if tau[x] != x]
    assert len(moved) == 2
    t1, t2 = moved
    c = bfs_pair(n, gens, (t1, t2), (0, a))
    seq = c[::-1] + hseq * e + c
    p = perm_of(n, gens, seq)
    target = tuple(a if x == 0 else 0 if x == a else x for x in range(n))
    assert p == target, (p, target)
    code = 0
    for j, s in enumerate(seq):
        code += s * 3**j
    return len(seq), code, len(hseq), e, len(c)

if __name__ == "__main__":
    for (n, m, k) in [(23,18,4), (23,16,6), (23,12,10)]:
        best = None
        for seed in range(1, 3):
            r = main(n, m, k, seed)
            if best is None or r[0] < best[0]: best = r
        L, code, lh, e, lc = best
        print(n, m, k, "L", L, "h", lh, "e", e, "c", lc)
        with open(".scratch/prefixlean/word_%d_%d_%d.txt" % (n, m, k), "w") as f:
            f.write("%d %d\n" % (L, code))
