import sys, itertools, random
sys.path.insert(0, '.')
from lib import *
part = sys.argv[1]

# ---------- groups ----------
class AbGroup:
    def __init__(self, mods):
        self.mods = mods
        self.els = list(itertools.product(*[range(m) for m in mods]))
    def power(self, x, k):
        return tuple((k * xi) % m for xi, m in zip(x, self.mods))
    def add(self, x, y):
        return tuple((a + b) % m for a, b, m in zip(x, y, self.mods))
    def order_set(self):
        s = set()
        for x in self.els:
            k = 1
            while self.power(x, k) != self.power(x, 0): k += 1
            s.add(k)
        return s

class SemiDirect:
    """Z/p x| Z/q with generator acting by multiplication by r (r^q = 1 mod p)"""
    def __init__(self, p, q, r):
        assert pow(r, q, p) == 1 and r % p != 1
        self.p, self.q, self.r = p, q, r
        self.els = [(a, b) for a in range(p) for b in range(q)]
        self.e = (0, 0)
    def mul(self, x, y):
        return ((x[0] + pow(self.r, x[1], self.p) * y[0]) % self.p, (x[1] + y[1]) % self.q)
    def power(self, x, k):
        res = self.e
        for _ in range(k): res = self.mul(res, x)
        return res

class Heis:
    def __init__(self, p):
        self.p = p
        self.els = list(itertools.product(range(p), repeat=3))
        self.e = (0, 0, 0)
    def mul(self, x, y):
        p = self.p
        return ((x[0] + y[0]) % p, (x[1] + y[1]) % p, (x[2] + y[2] + x[0] * y[1]) % p)
    def power(self, x, k):
        res = self.e
        for _ in range(k): res = self.mul(res, x)
        return res

def orders(G):
    s = set()
    e = G.power(G.els[0], 0)
    for x in G.els:
        k = 1
        while G.power(x, k) != e: k += 1
        s.add(k)
    return s

def has_good_order(G, A, n, pw):
    # pw[x][k] = x^k ; find bijection positions 1..n -> A with pw distinct
    A = list(A)
    used = [False] * n
    vals = set()
    def dfs(k):
        if k > n: return True
        for i in range(n):
            if used[i]: continue
            v = pw[A[i]][k]
            if v in vals: continue
            used[i] = True; vals.add(v)
            if dfs(k + 1): return True
            used[i] = False; vals.discard(v)
        return False
    return dfs(1)

if part == '41':
    tests = []
    for m in [5, 7, 11, 13]:
        tests.append(('Z/%d' % m, AbGroup([m])))
    tests += [('Z/25', AbGroup([25])), ('Z/5xZ/5', AbGroup([5, 5])), ('Z/35', AbGroup([35])), ('Z/49', AbGroup([49])), ('Z/7xZ/7', AbGroup([7, 7])),
              ('Z/11:Z/5', SemiDirect(11, 5, 3)), ('Heis(5)', Heis(5)), ('Z/5^3', AbGroup([5, 5, 5])), ('Z/31:Z/5', SemiDirect(31, 5, 2))]
    for name, G in tests:
        os_ = orders(G)
        nmax = 1
        while not any(2 <= o <= nmax + 2 for o in os_) and nmax + 1 < len(G.els): nmax += 1
        # nmax = largest n with no element order in 2..n+1
        nmax = max(n for n in range(1, len(G.els) + 1) if not any(2 <= o <= n + 1 for o in os_))
        nmax = min(nmax, len(G.els))
        pw = {x: [G.power(x, k) for k in range(nmax + 2)] for x in G.els}
        lim = {'Z/49': 5, 'Z/7xZ/7': 5, 'Heis(5)': 3, 'Z/5^3': 3, 'Z/31:Z/5': 3}.get(name, nmax)
        res = []
        for n in range(2, min(nmax, lim) + 1):
            cnt = 0; bad = []
            if name in ('Z/49', 'Z/7xZ/7') and n == 5:
                it = (tuple(random.sample(G.els, n)) for _ in range(300000)); mode = 'sample300k'
            else:
                it = itertools.combinations(G.els, n); mode = 'all'
            for A in it:
                cnt += 1
                if not has_good_order(G, A, n, pw):
                    bad.append(A)
                    if len(bad) > 3: break
            res.append((n, mode, cnt, bad[:3]))
        print(name, 'orders', sorted(os_), 'n<=', nmax, res, flush=True)

if part == '41s':
    random.seed(7)
    tests = [('Z/29:Z/7', SemiDirect(29, 7, 7), [4, 5]), ('Z/43:Z/7', SemiDirect(43, 7, 4), [4, 5]), ('Heis(7)', Heis(7), [4, 5]), ('Z/121', AbGroup([121]), [6, 7, 8, 9]), ('Z/11xZ/11', AbGroup([11, 11]), [6, 7, 8, 9])]
    for name, G, ns in tests:
        os_ = orders(G)
        pw = {x: [G.power(x, k) for k in range(max(ns) + 2)] for x in G.els}
        out = []
        for n in ns:
            assert not any(2 <= o <= n + 1 for o in os_)
            bad = []
            for _ in range(int(sys.argv[2])):
                A = tuple(random.sample(G.els, n))
                if not has_good_order(G, A, n, pw): bad.append(A); break
            out.append((n, bad))
        print(name, sorted(os_), out, flush=True)
