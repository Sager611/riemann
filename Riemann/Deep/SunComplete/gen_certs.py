#!/usr/bin/env python3
"""gen_certs.py -- regenerates certA.json and certB.json.gz (NOTE.md sec. 5); deterministic, ~1-2 min.
Part 1: finite fields F_p[x]/(f) (f = first monic polynomial in itertools.product order with x of order q-1).
Part 2: certificates for the exceptional sets EA, EB of verify.py (She / alternating / geometric / explicit cycles)."""

# ---- Part 1. Honest finite-field arithmetic for F_q, q = p^k, elements encoded as integers 0..q-1
# (base-p digits = coefficients of 1, x, ..., x^(k-1) in F_p[x]/(f)).
from math import gcd

def factor_distinct(n):
    ps = []; d = 2
    while d*d <= n:
        if n % d == 0:
            ps.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: ps.append(n)
    return ps

def is_prime(n):
    if n < 2: return False
    d = 2
    while d*d <= n:
        if n % d == 0: return False
        d += 1
    return True

def prime_power(q):
    for p in range(2, q+1):
        if q % p == 0:
            k = 0; x = q
            while x % p == 0: x //= p; k += 1
            return (p, k) if x == 1 else None
    return None

class GF:
    def __init__(self, q):
        pk = prime_power(q)
        assert pk is not None, q
        self.p, self.k = pk; self.q = q
        p, k = self.p, self.k
        self.pw = [p**i for i in range(k+1)]
        if k == 1:
            self.f = None
            g = self._prim_root_mod_p(p)
            exp = [1]*(q-1)
            for i in range(1, q-1): exp[i] = exp[i-1]*g % p
            self.g = g
        else:
            self.f = self._find_prim_poly()
            exp = [0]*(q-1); exp[0] = 1
            for i in range(1, q-1): exp[i] = self._mulx(exp[i-1])
            assert self._mulx(exp[q-2]) == 1
            self.g = p  # the element x
        assert len(set(exp)) == q-1 and 0 not in exp  # field + x (resp. g) primitive
        self.exp = exp
        self.log = [None]*q
        for i, e in enumerate(exp): self.log[e] = i
        n = q-1
        self.isprim = bytearray(q)
        for i in range(n):
            if gcd(i, n) == 1: self.isprim[exp[i]] = 1
        # quadratic character via log parity (q odd)
    def _prim_root_mod_p(self, p):
        if p == 2: return 1
        ls = factor_distinct(p-1)
        for g in range(2, p):
            if all(pow(g, (p-1)//l, p) != 1 for l in ls): return g
    def digits(self, a):
        d = []
        for _ in range(self.k): d.append(a % self.p); a //= self.p
        return d
    def undigits(self, d):
        a = 0
        for c in reversed(d): a = a*self.p + c
        return a
    def add(self, a, b):
        if self.k == 1: return (a+b) % self.p
        p = self.p; r = 0; m = 1
        for _ in range(self.k):
            r += ((a % p + b % p) % p)*m; a //= p; b //= p; m *= p
        return r
    def neg(self, a):
        if self.k == 1: return (-a) % self.p
        p = self.p; r = 0; m = 1
        for _ in range(self.k):
            r += ((-(a % p)) % p)*m; a //= p; m *= p
        return r
    def sub(self, a, b): return self.add(a, self.neg(b))
    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        return self.exp[(self.log[a] + self.log[b]) % (self.q-1)]
    def inv(self, a): return self.exp[(-self.log[a]) % (self.q-1)]
    def _mulx(self, a):
        # multiply by x modulo monic f (f given as list of k+1 coeffs, f[k]=1)
        p, k, f = self.p, self.k, self.f
        d = self.digits(a)
        top = d[k-1]
        nd = [0] + d[:k-1]
        for i in range(k): nd[i] = (nd[i] - top*f[i]) % p
        return self.undigits(nd)
    def _find_prim_poly(self):
        p, k, q = self.p, self.k, self.q
        from itertools import product
        ls = factor_distinct(q-1)
        for cs in product(range(p), repeat=k):
            f = list(cs) + [1]
            if f[0] == 0: continue
            self.f = f
            # order of x must be exactly q-1
            a = 1; ok = True
            seen_one = False
            # compute x^((q-1)/l) via repeated multiplication is O(q); fine for our sizes
            powx = [1]
            for i in range(1, q):
                a = self._mulx(a)
                if a == 1 and i < q-1: ok = False; break
            if ok and a == 1: return f
        raise ValueError



# ---- Part 2. certificates

# Computational verification for the exceptional sets E_A, E_B (q > 7).
import sys, json, random, time, gzip, os
from math import gcd
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import verify   # exceptional sets EA, EB (re-derived by verify.py, Prop. 7)

T0 = time.time()
EA = [q for q in verify.EA if q > 7]
EB = [q for q in verify.EB if q > 7]

def ham_dfs(V, adj, limit_nodes=5*10**6):
    n = len(V)
    nb = [[j for j in range(n) if j != i and adj(V[i], V[j])] for i in range(n)]
    s = min(range(n), key=lambda i: len(nb[i]))
    path = [s]; used = [False]*n; used[s] = True
    cnt = [0]
    sys.setrecursionlimit(10000)
    def rec():
        cnt[0] += 1
        if cnt[0] > limit_nodes: raise TimeoutError
        if len(path) == n: return s in nb[path[-1]]
        end = path[-1]
        cand = [j for j in nb[end] if not used[j]]
        cand.sort(key=lambda j: sum(1 for k in nb[j] if not used[k]))
        for j in cand:
            used[j] = True; path.append(j)
            if rec(): return True
            path.pop(); used[j] = False
        return False
    ok = rec()
    return ([V[i] for i in path] if ok else None), cnt[0]

def ham_posa(V, adj, seed, restarts=40):
    n = len(V); rng = random.Random(seed)
    for att in range(restarts):
        start = V[rng.randrange(n)]
        path = [start]
        unv = [v for v in V if v != start]; rng.shuffle(unv)
        steps = 0; maxsteps = 50*n + 2000
        while steps < maxsteps:
            steps += 1
            end = path[-1]
            if unv:
                hit = -1
                for idx in range(len(unv)):
                    if adj(end, unv[idx]): hit = idx; break
                if hit >= 0:
                    w = unv[hit]; unv[hit] = unv[-1]; unv.pop(); path.append(w); continue
            elif adj(end, path[0]):
                return path
            L = len(path)
            cands = [i for i in range(L-2) if adj(end, path[i])]
            if not cands or rng.random() < 0.1:
                path.reverse(); continue
            i = cands[rng.randrange(len(cands))]
            path[i+1:] = path[i+1:][::-1]
    return None

def check_cycle(F, cyc, kind, a0=None):
    q = F.q
    if kind == 'A':
        assert sorted(cyc) == list(range(q))
        for i in range(q):
            assert F.isprim[F.add(cyc[i], cyc[(i+1) % q])], ('A', q, i)
    else:
        assert sorted(cyc) == list(range(1, q))
        N = q-1
        for i in range(N):
            assert F.isprim[F.add(a0, F.mul(cyc[i], cyc[(i+1) % N]))], ('B', q, a0, i)
    return True

# ---------------- (A) ----------------
certA = {}
for q in EA:
    F = GF(q)
    rec = None
    if F.k == 1:
        p = q
        xs = [x for x in range(p) if F.isprim[x] and F.isprim[(x+1) % p] and F.isprim[(x+2) % p]]
        if xs:
            x = xs[0]; m = (p-1)//2; h = x*pow(2, -1, p) % p
            seq = {}
            for kk in range(0, m+1): seq[2*kk+1] = (h + 2*kk + 1) % p
            for kk in range(1, m+1): seq[2*kk] = (h - 2*kk + 1) % p
            cyc = [seq[i] for i in range(1, p+1)]
            check_cycle(F, cyc, 'A')
            rec = {'type': 'she', 'x': x}
    if rec is None:
        V = list(range(q))
        adj = lambda a, b: F.isprim[F.add(a, b)] == 1
        cyc = None
        if q <= 40:
            cyc, _ = ham_dfs(V, adj)
        if cyc is None:
            cyc = ham_posa(V, adj, seed=q)
        assert cyc is not None, ('A search failed', q)
        check_cycle(F, cyc, 'A')
        rec = {'type': 'cycle', 'f': F.f, 'cycle': cyc}
    certA[str(q)] = rec
print("(A) done:", len(certA), "fields;", sum(1 for v in certA.values() if v['type']=='she'), "by She construction;",
      [int(k) for k, v in certA.items() if v['type']=='cycle'], "by explicit cycle; %.1fs" % (time.time()-T0), flush=True)
json.dump(certA, open(os.path.join(HERE, 'certA.json'), 'w'))

# ---------------- (B) ----------------
def addv(F, a0, arr):
    if F.k == 1: return (arr + a0) % F.p
    p = F.p; res = np.zeros_like(arr); a = a0; x = arr.copy()
    for i in range(F.k):
        d = x % p; e = a % p
        res += ((d + e) % p) * F.pw[i]
        x //= p; a //= p
    return res

certB = {}
failures = []
stats = {'alt': 0, 'geo': 0, 'search': 0}
for q in EB:
    F = GF(q); n = q-1
    expA = np.array(F.exp, dtype=np.int64)
    logA = np.array([(-1 if v is None else v) for v in F.log], dtype=np.int64)
    isprimA = np.frombuffer(bytes(F.isprim), dtype=np.uint8)
    prims = [z for z in range(q) if F.isprim[z]]
    recs = []
    for a0 in range(1, q):
        cert = None
        if q % 2 == 1:
            m = n//2
            for u in range(1, q):
                if F.log[u] % 2 == 0 or not F.isprim[F.add(a0, u)]: continue
                lu = F.log[u]
                for w in prims:
                    v = F.sub(w, a0)
                    if v == 0 or F.log[v] % 2 == 0: continue
                    j = ((F.log[v] - lu) % n)//2
                    if gcd(j, m) == 1:
                        cert = ('alt', u, F.exp[(F.log[v]-lu) % n]); break
                if cert: break
            if cert:
                _, u, r = cert; lu, lr = F.log[u], F.log[r]
                # cycle: a_{2j+1} = r^j (j=0..m-1), a_{2j} = u r^{1-j} (j=1..m)
                cyc = np.zeros(n, dtype=np.int64)
                jj = np.arange(m)
                cyc[0::2] = expA[(jj*lr) % n]
                cyc[1::2] = expA[(lu + (1 - (jj+1))*lr) % n]
        else:
            m = (n-1)//2
            sq = lambda z: F.exp[(F.log[z]*(q//2)) % n]   # square root in char 2
            gens = [z for z in range(1, q) if gcd(F.log[z], n) == 1]
            for u in range(1, q):
                if not F.isprim[F.add(a0, u)]: continue
                for s in gens:
                    if F.isprim[F.add(a0, F.mul(u, s))] and F.isprim[F.add(a0, F.mul(u, F.mul(s, s)))]:
                        cert = ('geo', u, s); break
                if cert: break
            if cert:
                _, u, s = cert; r = F.mul(s, s); t = F.mul(s, sq(u))
                lt, lr, lu = F.log[t], F.log[r], F.log[u]
                cyc = np.zeros(n, dtype=np.int64)
                cyc[0::2] = expA[(lt + np.arange(m+1)*lr) % n]
                jj = np.arange(1, m+1)
                cyc[1::2] = expA[(lu - lt + (1-jj)*lr) % n]
        if cert:
            assert sorted(cyc.tolist()) == list(range(1, q)), (q, a0, cert)
            prods = expA[(logA[cyc] + logA[np.roll(cyc, -1)]) % n]
            assert isprimA[addv(F, a0, prods)].all(), (q, a0, cert)
            stats[cert[0]] += 1
            recs.append([a0, cert[0], cert[1], cert[2]])
            continue
        # fallback: explicit search
        V = list(range(1, q))
        adj = lambda x, y: F.isprim[F.add(a0, F.mul(x, y))] == 1
        cyc = None; proven_none = False
        if q <= 40:
            cyc, nodes = ham_dfs(V, adj)
            if cyc is None: proven_none = True
        if cyc is None and not proven_none:
            cyc = ham_posa(V, adj, seed=1000003*q + a0)
        if cyc is None:
            failures.append((q, a0, 'none (exhaustive)' if proven_none else 'not found'))
            recs.append([a0, 'none'])
            continue
        check_cycle(F, cyc, 'B', a0)
        stats['search'] += 1
        recs.append([a0, 'cycle', cyc])
    certB[str(q)] = {'f': F.f, 'g': F.g, 'certs': recs}
    if q > 1000: print("  (B) q=%d done %.1fs" % (q, time.time()-T0), flush=True)
print("(B) done:", len(certB), "fields; stats", stats, "failures", failures, "%.1fs" % (time.time()-T0), flush=True)
json.dump(certB, gzip.open(os.path.join(HERE, 'certB.json.gz'), 'wt'))

