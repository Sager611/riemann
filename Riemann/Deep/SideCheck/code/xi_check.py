#!/usr/bin/env python3
"""sidecheck: independent check of a 3-edge-colouring of K_n (every K_k sees all 3 colours).
Method (different from Riemann/Deep/Xi6/verify.py, which brute-forces all k-sets):
for each pair of colours {a,b}, build the union graph G_a u G_b as bitmasks and compute its clique
number exactly by a pivoting Bron-Kerbosch search; also count its k-cliques by ordered extension.
Every K_k sees all 3 colours  <=>  every pair union is K_k-free  <=>  each clique number <= k-1."""
import sys, itertools
from math import comb

def parse(fn):
    lines = [l.rstrip('\n') for l in open(fn)]
    body = [l.strip() for l in lines if l.strip() and not l.lstrip().startswith('#')]
    n = int(body[0]); rows = body[1:]
    assert len(rows) == n, (len(rows), n)
    col = [[None]*n for _ in range(n)]
    for i, r in enumerate(rows):
        assert len(r) == n
        for j, ch in enumerate(r):
            if i == j: assert ch == '.'; continue
            assert ch in '012'; col[i][j] = int(ch)
    for i in range(n):
        for j in range(n):
            if i != j: assert col[i][j] == col[j][i]
    return n, col

def clique_number(adj, n):
    best = 0
    def bk(R, P, X):
        nonlocal best
        if P == 0 and X == 0:
            best = max(best, R); return
        if R + bin(P).count('1') <= best: return
        PX = P | X
        # pivot: vertex of P|X with most neighbours in P
        u = max((v for v in range(n) if PX >> v & 1), key=lambda v: bin(P & adj[v]).count('1'))
        cand = P & ~adj[u]
        while cand:
            v = (cand & -cand).bit_length() - 1
            cand &= cand - 1
            bk(R + 1, P & adj[v], X & adj[v])
            P &= ~(1 << v); X |= 1 << v
    bk(0, (1 << n) - 1, 0)
    return best

def count_cliques(adj, n, k):
    cnt = 0
    def ext(size, P):
        nonlocal cnt
        if size == k: cnt += 1; return
        while P:
            v = (P & -P).bit_length() - 1
            P &= P - 1
            ext(size + 1, P & adj[v])
    # only higher-numbered neighbours to count each clique once
    hi = [adj[v] & ~((1 << (v + 1)) - 1) for v in range(n)]
    def ext2(size, P):
        nonlocal cnt
        if size == k: cnt += 1; return
        while P:
            v = (P & -P).bit_length() - 1
            P &= P - 1
            ext2(size + 1, P & hi[v])
    for v in range(n): ext2(1, hi[v])
    return cnt

def main():
    fn = sys.argv[1]; k = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    n, col = parse(fn)
    sizes = [sum(1 for i in range(n) for j in range(i+1, n) if col[i][j] == c) for c in range(3)]
    print(f"file={fn} n={n} k={k} class sizes={sizes} total={sum(sizes)} (=C(n,2)={comb(n,2)})")
    ok = True
    for a, b in [(0, 1), (0, 2), (1, 2)]:
        adj = [sum(1 << j for j in range(n) if j != i and col[i][j] in (a, b)) for i in range(n)]
        w = clique_number(adj, n)
        ck = count_cliques(adj, n, k); ck1 = count_cliques(adj, n, k - 1)
        print(f"  union of colours {a},{b}: clique number {w}; #K_{k} = {ck}; #K_{k-1} = {ck1}")
        ok &= (w <= k - 1 and ck == 0)
    # also: each colour class must be non-empty and every vertex degree profile
    degs = sorted(set(tuple(sorted(sum(1 for j in range(n) if j != i and col[i][j] == c) for c in range(3))) for i in range(n)))
    print("  colour-degree profiles:", degs)
    print("RESULT:", "VALID (every K_%d sees all three colours)" % k if ok else "INVALID")
    return 0 if ok else 1

if __name__ == '__main__':
    sys.exit(main())
