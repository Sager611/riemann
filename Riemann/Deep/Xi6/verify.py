#!/usr/bin/env python3
"""verify.py - independent exhaustive checker (stdlib only).

Usage: python3 verify.py <colouring-file> [k]   (default k = 6)

Colouring file: lines starting with '#' are comments; first other line is n;
then n rows of n characters: row i, column j = colour (0, 1 or 2) of edge {i, j},
'.' on the diagonal.

Checks that the matrix is symmetric with entries in {0,1,2}, then enumerates ALL
C(n,k) k-subsets of the vertex set (k = 6, or n <= 30: explicit loop over every subset, no pruning;
otherwise a complete depth-first enumeration that counts the completions of already 3-coloured prefixes at once) and counts the ones whose
C(k,2) edges use fewer than 3 colours.  Exit code 0 iff that count is 0.
"""
import sys, itertools, math

def read(fn):
    rows = [l.strip() for l in open(fn) if l.strip() and not l.startswith('#')]
    n = int(rows[0]); rows = rows[1:1 + n]
    assert len(rows) == n, "need n rows"
    M = [[0] * n for _ in range(n)]
    for i, r in enumerate(rows):
        assert len(r) == n, f"row {i} has length {len(r)}"
        for j, ch in enumerate(r):
            if i == j:
                assert ch == '.', f"diagonal ({i},{i}) must be '.'"
            else:
                assert ch in '012', f"bad colour {ch!r} at ({i},{j})"
                M[i][j] = 1 << int(ch)          # colour c -> bit c
    for i in range(n):
        for j in range(n):
            assert M[i][j] == M[j][i], f"not symmetric at ({i},{j})"
    return n, M

def main():
    fn = sys.argv[1]; k = int(sys.argv[2]) if len(sys.argv) > 2 else 6
    n, M = read(fn)
    total = bad = 0; first = None
    if k == 6:
        # explicit nested loops over a<b<c<d<e<f, OR of the 15 edge-colour bits
        for a in range(n):
            Ma = M[a]
            for b in range(a + 1, n):
                o2 = Ma[b]; Mb = M[b]
                for c in range(b + 1, n):
                    Mc = M[c]; o3 = o2 | Ma[c] | Mb[c]
                    for d in range(c + 1, n):
                        Md = M[d]; o4 = o3 | Ma[d] | Mb[d] | Mc[d]
                        for e in range(d + 1, n):
                            Me = M[e]; o5 = o4 | Ma[e] | Mb[e] | Mc[e] | Md[e]
                            for f in range(e + 1, n):
                                total += 1
                                if (o5 | Ma[f] | Mb[f] | Mc[f] | Md[f] | Me[f]) != 7:
                                    bad += 1
                                    if first is None: first = (a, b, c, d, e, f)
    elif n <= 30:
        for S in itertools.combinations(range(n), k):
            total += 1
            o = 0
            for x, y in itertools.combinations(S, 2): o |= M[x][y]
            if o != 7:
                bad += 1
                if first is None: first = S
    else:
        # depth-first over increasing vertex sequences, OR-ing edge colours; once a prefix already
        # shows all 3 colours every completion does too, so its C(.,.) completions are counted at once.
        # Still exhaustive: every k-subset is accounted for exactly once.
        sys.setrecursionlimit(10000)
        def dfs(pref, o, start):
            nonlocal total, bad, first
            r = k - len(pref)
            if r == 0:
                total += 1
                if o != 7:
                    bad += 1
                    if first is None: first = tuple(pref)
                return
            if o == 7:
                total += math.comb(n - start, r); return
            for x in range(start, n - r + 1):
                Mx = M[x]; o2 = o
                for y in pref: o2 |= Mx[y]
                pref.append(x); dfs(pref, o2, x + 1); pref.pop()
        dfs([], 0, 0)
    assert total == math.comb(n, k)
    print(f"n={n} k={k} subsets checked={total} (=C({n},{k})) two-coloured={bad} first_bad={first}")
    print("VALID: every K_%d uses all three colours" % k if bad == 0 else "INVALID")
    sys.exit(0 if bad == 0 else 1)

if __name__ == '__main__':
    main()

