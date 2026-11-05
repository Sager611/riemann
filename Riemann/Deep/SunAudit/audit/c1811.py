import sys, itertools, random
from fractions import Fraction
sys.path.insert(0, '.')
from lib import *
part = sys.argv[1]

def bip_match(n, ok):
    """perfect matching rows 1..n -> cols 1..n with ok(r,c); returns match or None"""
    adj = {r: [c for c in range(1, n + 1) if ok(r, c)] for r in range(1, n + 1)}
    matchc = {}
    def aug(r, seen):
        for c in adj[r]:
            if c in seen: continue
            seen.add(c)
            if c not in matchc or aug(matchc[c], seen):
                matchc[c] = r
                return True
        return False
    sys.setrecursionlimit(10000)
    for r in range(1, n + 1):
        if not aug(r, set()):
            return None
    return {r: c for c, r in matchc.items()}

if part == '43':
    PR = set(primes_upto(200000))
    bad1 = [n for n in range(1, 301) if bip_match(n, lambda k, s: k * s + 1 in PR) is None]
    bad2 = [n for n in range(1, 301) if bip_match(n, lambda k, s: k * s - 1 in PR) is None]
    print('4.3(i) no perm n<=300:', bad1)
    print('4.3(ii) no perm n<=300:', bad2, '(conj allows n<=2)')
if part == '44':
    PR = set(primes_upto(300000))
    bad1 = [n for n in range(1, 301) if bip_match(n, lambda k, s: k*k + k*s + s*s in PR) is None]
    bad2 = [n for n in range(1, 301) if bip_match(n, lambda k, s: k*k + s*s in PR) is None]
    print('4.4(i) no perm n<=300:', bad1)
    print('4.4(ii) no perm n<=300:', bad2, '(conj: only n=7)')
if part == '45':
    P = primes_upto(30000); PR = set(primes_upto(70000))
    pp = [None] + P
    bad = [n for n in range(1, 301) if bip_match(n, lambda k, s: pp[k] + pp[s] + 1 in PR) is None]
    print('4.5 no perm n<=300:', bad)
if part == '46':
    P = primes_upto(2000); pp = [None] + P
    def isp(x): return is_prime(x)
    res = []
    for n in range(1, 11):
        ok = [[(k > 0 and s > 0 and isp(pp[k] * pp[s] - 2)) for s in range(n + 1)] for k in range(n + 1)]
        even = odd = False
        # DFS over permutations with parity tracking
        used = [False] * (n + 1)
        perm = []
        def dfs(k):
            global even, odd
            if even and odd: return
            if k > n:
                inv = sum(1 for i in range(n) for j in range(i + 1, n) if perm[i] > perm[j])
                if inv % 2 == 0: even = True
                else: odd = True
                return
            for s in range(1, n + 1):
                if not used[s] and ok[k][s]:
                    used[s] = True; perm.append(s)
                    dfs(k + 1)
                    perm.pop(); used[s] = False
        dfs(1)
        res.append((n, even, odd))
    print('4.6 (n, even exists, odd exists):', res)
if part == '42':
    def per_mod(M, mod):
        n = len(M)
        # Ryser with Gray code
        total = 0
        rowsums = [0] * n
        prev = 0
        for g in range(1, 1 << n):
            gray = g ^ (g >> 1)
            diff = gray ^ prev
            j = diff.bit_length() - 1
            if gray & diff:
                for i in range(n): rowsums[i] += M[i][j]
            else:
                for i in range(n): rowsums[i] -= M[i][j]
            prev = gray
            prod = 1
            for i in range(n):
                prod = prod * rowsums[i] % mod
                if prod == 0: break
            bits = bin(gray).count('1')
            if (n - bits) % 2: total -= prod
            else: total += prod
        return total % mod
    out = []
    for n in range(2, int(sys.argv[2]) + 1):
        if n % 4 == 2: continue
        m = n - 1
        mod = n * n
        M = [[pow(i, j - 1, mod) for j in range(1, m + 1)] for i in range(1, m + 1)]
        v = per_mod(M, mod)
        out.append((n, v % n == 0, (v % (n*n) == 0) if n % 3 == 1 else None))
        print(n, v % n == 0, (v % (n*n) == 0) if n % 3 == 1 else '-', flush=True)
if part == '47':
    # (4.3) n>6 path sum 1/(a+b) = 1; (4.4) n>7 cyclic; (4.5) n>7 sum 1/(a^2-b^2) = 0
    def search(n, cyclic, f, target, limit_nodes=5*10**6):
        # DFS; f(a,b) Fraction; returns example or None (exhaustive unless node limit)
        best = None
        used = [False] * (n + 1)
        seq = []
        nodes = [0]
        def dfs(s):
            nodes[0] += 1
            if nodes[0] > limit_nodes: raise TimeoutError
            if len(seq) == n:
                tot = s + (f(seq[-1], seq[0]) if cyclic else 0)
                return list(seq) if tot == target else None
            for x in range(1, n + 1):
                if used[x]: continue
                if cyclic and not seq and x != 1: continue
                ns = s + (f(seq[-1], x) if seq else 0)
                used[x] = True; seq.append(x)
                r = dfs(ns)
                seq.pop(); used[x] = False
                if r: return r
            return None
        return dfs(Fraction(0)), nodes[0]
    f1 = lambda a, b: Fraction(1, a + b)
    f2 = lambda a, b: Fraction(1, a*a - b*b)
    for n in range(2, int(sys.argv[2]) + 1):
        try:
            r1 = search(n, False, f1, 1)
        except TimeoutError: r1 = ('timeout',)
        try:
            r2 = search(n, True, f1, 1)
        except TimeoutError: r2 = ('timeout',)
        try:
            r3 = search(n, False, f2, 0)
        except TimeoutError: r3 = ('timeout',)
        print(n, '4.3:', r1[0], '| 4.4:', r2[0], '| 4.5:', r3[0], flush=True)
