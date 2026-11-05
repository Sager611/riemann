import sys, itertools, random
sys.path.insert(0, '.')
from lib import *
part = sys.argv[1]

def S_of_A(A, n, add, mulk, zero):
    # set of sums sum_{k=1}^n k a_k over injective sequences from A
    A = list(A)
    m = len(A)
    # DP over masks: dict mask -> set of partial sums
    from collections import defaultdict
    cur = {0: {zero}}
    for k in range(1, n + 1):
        nxt = defaultdict(set)
        for mask, sums in cur.items():
            for i in range(m):
                if mask >> i & 1: continue
                t = mulk(k, A[i])
                s2 = nxt[mask | (1 << i)]
                for s in sums: s2.add(add(s, t))
        cur = nxt
    res = set()
    for sums in cur.values(): res |= sums
    return res

if part == '49p':
    # F = Z/p, p > n+1, |A| >= n + delta(n,3); normalize 0,1 in A (affine invariance)
    out = []
    for p in [5, 7, 11, 13]:
        add = lambda a, b, p=p: (a + b) % p
        mulk = lambda k, a, p=p: (k * a) % p
        for n in range(1, p - 1):
            if not p > n + 1: continue
            bads = []
            for size in range(max(2, n + (1 if n == 3 else 0)), p + 1):
                bound = min(p, (size - n) * n * (n + 1) // 2 + n * (n * n - 1) // 6 + 1)
                others = list(range(2, p))
                cnt = 0
                for rest in itertools.combinations(others, size - 2):
                    A = (0, 1) + rest
                    if size > 9 and p == 13 and n > 6: break
                    cnt += 1
                    s = S_of_A(A, n, add, mulk, 0)
                    if len(s) < bound:
                        bads.append((p, n, A, len(s), bound)); break
            out.append((p, n, bads[:2]))
            print(p, n, bads[:2], flush=True)
if part == '49z':
    # F = Q: A subset of {0..L} containing 0 and with gcd 1 (affine normalization)
    from math import gcd
    L = int(sys.argv[2]); maxsize = int(sys.argv[3])
    add = lambda a, b: a + b
    mulk = lambda k, a: k * a
    for n in range(1, maxsize + 1):
        bads = []
        for size in range(n + (1 if n == 3 else 0), maxsize + 1):
            bound = (size - n) * n * (n + 1) // 2 + n * (n * n - 1) // 6 + 1
            for rest in itertools.combinations(range(1, L + 1), size - 1):
                g = 0
                for x in rest: g = gcd(g, x)
                if g != 1: continue
                A = (0,) + rest
                s = S_of_A(A, n, add, mulk, 0)
                if len(s) < bound:
                    bads.append((n, A, len(s), bound)); break
        print('Q n=%d' % n, bads[:2], flush=True)
if part == '410':
    def ok_numbering(A, add3):
        A = list(A); n = len(A)
        # cyclic arrangement starting with A[0]; DFS with partial distinctness
        seq = [A[0]]; used = [False] * n; used[0] = True
        sums = set()
        def dfs():
            if len(seq) == n:
                # closing sums: (a_{n-1}, a_n, a_1), (a_n, a_1, a_2)
                s1 = add3(seq[-2], seq[-1], seq[0]); s2 = add3(seq[-1], seq[0], seq[1])
                return s1 not in sums and s2 not in sums and s1 != s2
            for i in range(1, n):
                if used[i]: continue
                x = A[i]
                if len(seq) >= 2:
                    s = add3(seq[-2], seq[-1], x)
                    if s in sums: continue
                    sums.add(s)
                used[i] = True; seq.append(x)
                if dfs(): return True
                seq.pop(); used[i] = False
                if len(seq) >= 2: sums.discard(add3(seq[-2], seq[-1], x))
            return False
        return dfs()
    groups = []
    for m in range(4, int(sys.argv[2]) + 1):
        groups.append(('Z/%d' % m, [m]))
    groups += [('Z/3xZ/3', [3, 3]), ('Z/3xZ/9', [3, 9]), ('Z/5xZ/5', [5, 5]), ('Z/3^3', [3, 3, 3]), ('Z/2xZ/2', [2, 2]), ('Z/2xZ/4', [2, 4]), ('Z/2^3', [2, 2, 2]), ('Z/2xZ/6', [2, 6]), ('Z/4xZ/4', [4, 4])]
    maxsz = int(sys.argv[3])
    for name, mods in groups:
        els = list(itertools.product(*[range(m) for m in mods]))
        zero = els[0]
        add3 = lambda a, b, c, mods=mods: tuple((x + y + z) % m for x, y, z, m in zip(a, b, c, mods))
        bad = []
        others = els[1:]
        G = len(els)
        sizes = list(range(4, min(G, maxsz) + 1))
        if G > maxsz: sizes += [G]  # full group too
        for size in sizes:
            for rest in itertools.combinations(others, size - 1):
                A = (zero,) + rest
                if not ok_numbering(A, add3):
                    bad.append(A); break
        print(name, 'sizes', sizes[0], '..', sizes[-1], 'bad:', bad[:3], flush=True)
