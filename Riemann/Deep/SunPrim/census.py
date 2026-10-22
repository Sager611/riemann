#!/usr/bin/env python3
"""Census for four conjectures of Z.-W. Sun, arXiv 1309.1679 (printed numbering of v9 / Nanjing Univ. J. Math. Biquarterly
36 (2019)), and the OEIS A228766 variant. Run: python3 census.py [--ext]   (about 1 min; --ext adds Z/13..Z/16, ~5 min)

Part 1, all abelian groups of order <= 12 (and Z/13..Z/16 with --ext); exhaustive, up to translations and automorphisms:
  3.2     path b_1 = a_1, ..., b_n with the n-1 differences b_i - b_{i+1} distinct; hyp: n !| |G|, or n even and Syl_2 cyclic
  3.3(ii) circular differences distinct; v9: |G| odd, n > 2, n !| |G|;  v4-v8: (n odd or n !| |G|), 3 < n < |G|;
          Sun's homepage copy of the journal version (196a.pdf, 2020-05-07): |G| odd, n = |G| - 1 (= R-sequenceability)
  3.4(i)  circular sums a_i + 2a_{i+1} distinct; hyp: |G| not divisible by 3, n > 3
  A228766 (= 3.3(i) of v1-v8) circular sums a_i + a_{i+1} distinct; hyp: n > 2, n odd or n !| |G|; v9 3.3(i): |G| odd, n > 2
  Asserts: every 3.2 failure has {a_1..a_n} or {a_2..a_n} a coset of a subgroup H with |H| >= 3 and sum(H) = 0 (and all
  such configurations fail); the homepage 3.3(ii) and v9 3.3(i) have no failures.
Part 2, Conjecture 3.8: for every prime power q <= 128 and a_0 in F_q, Hamilton cycle in G_{q,a0} (x ~ y iff x != y and
  a_0 + xy primitive). Cycles re-checked with separate polynomial arithmetic; non-Hamiltonian = complete search.
  Asserts that (11, 2) is the only failure with q > 7.
Option --primes 17,19: for Z/p, all subsets containing {0, 1} (enough by 2-transitivity of Aff(Z/p)), list the n with
  failures of the circular-difference (3.3(ii)) and circular x + 2y (3.4(i)) properties (Z/17: about 3 min)."""
import itertools, sys, resource, random

EXT = "--ext" in sys.argv
PRIMES = [int(t) for t in sys.argv[sys.argv.index("--primes") + 1].split(",")] if "--primes" in sys.argv else []

def memcheck():
    if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss > 250 * 1024 * 1024:
        sys.exit("memory limit")

# ---------------------------------------------------------------- Part 1: groups
class Group:
    def __init__(self, ns):
        self.ns, self.el = ns, list(itertools.product(*[range(n) for n in ns]))
        self.ix = {e: i for i, e in enumerate(self.el)}
        N = self.N = len(self.el)
        self.add = [[self.ix[tuple((a + b) % n for a, b, n in zip(x, y, ns))] for y in self.el] for x in self.el]
        self.neg = [self.ix[tuple((-a) % n for a, n in zip(x, ns))] for x in self.el]
        self.zero = 0
        self.sub = [[self.add[x][self.neg[y]] for y in range(N)] for x in range(N)]
        self.s2 = [[self.add[x][self.add[y][y]] for y in range(N)] for x in range(N)]
        self.syl2_cyclic = sum(1 for n in ns if n % 2 == 0) <= 1
        self.name = "x".join("Z%d" % n for n in ns)
        self.total = self.gsum(range(N))
        self.auts = self._auts()

    def gsum(self, S):
        r = 0
        for x in S:
            r = self.add[r][x]
        return r

    def mul(self, k, x):
        r = 0
        for _ in range(k):
            r = self.add[r][x]
        return r

    def _auts(self):
        gens = [self.ix[tuple(int(i == j) for j in range(len(self.ns)))] for i in range(len(self.ns))]
        cands = [[g for g in range(self.N) if self.mul(n, g) == 0] for n in self.ns]
        out = []
        for imgs in itertools.product(*cands):
            phi = [self.gsum(self.mul(c, g) for c, g in zip(x, imgs)) for x in self.el]
            if len(set(phi)) == self.N:
                out.append(phi)
        return out

    def canon(self, A, mark=None):
        best = None
        for phi in self.auts:
            B = [phi[a] for a in A]
            for t in ([self.neg[phi[mark]]] if mark is not None else [self.neg[b] for b in B]):
                C = tuple(sorted(self.add[b][t] for b in B))
                if best is None or C < best:
                    best = C
        return best

    def zero_sum_cosets(self):
        subs = set()
        for S in itertools.combinations(range(self.N), 3):
            H, grow = {0, *S}, True
            while grow:
                new = {self.add[x][y] for x in H for y in H} - H
                grow = bool(new)
                H |= new
            subs.add(frozenset(H))
        subs |= {frozenset(range(self.N))}
        cos = set()
        for H in subs:
            if len(H) >= 3 and self.gsum(H) == 0:
                cos |= {frozenset(self.add[x][h] for h in H) for x in range(self.N)}
        return cos

def circ_ok(G, A, f, kind):
    n, N = len(A), G.N
    if kind == "diff":                      # n nonzero distinct differences with sum 0
        if n > N - 1 or (n == N - 1 and G.total != 0) or (n == N - 2 and G.total == 0):
            return False
    if kind == "sum" and n == N and G.total != 0:
        return False                        # values = G, but sum of values = 2 sum(G) = 0
    s, rem = A[0], set(A[1:])
    def dfs(used, cur):
        if not rem:
            return not (used >> f[cur][s]) & 1
        for y in list(rem):
            v = f[cur][y]
            if (used >> v) & 1:
                continue
            rem.discard(y)
            ok = dfs(used | (1 << v), y)
            rem.add(y)
            if ok:
                return True
        return False
    return dfs(0, s)

def path_ok(G, A):                          # A canonical with marked point 0 = b_1
    if len(A) == G.N and G.total == 0:
        return False
    rem = set(A) - {0}
    def dfs(used, cur):
        if not rem:
            return True
        for y in list(rem):
            v = G.sub[cur][y]
            if (used >> v) & 1:
                continue
            rem.discard(y)
            ok = dfs(used | (1 << v), y)
            rem.add(y)
            if ok:
                return True
        return False
    return dfs(0, 0)

def part1():
    groups = [(m,) for m in range(2, 13)] + [(2, 2), (2, 4), (2, 2, 2), (3, 3), (2, 6)]
    if EXT:
        groups += [(m,) for m in range(13, 17)]
    hyps = {
        "3.2": lambda G, n: G.N % n != 0 or (n % 2 == 0 and G.syl2_cyclic),
        "3.3(ii) v9": lambda G, n: G.N % 2 == 1 and n > 2 and G.N % n != 0,
        "3.3(ii) v4-v8": lambda G, n: (n % 2 == 1 or G.N % n != 0) and 3 < n < G.N,
        "3.3(ii) homepage": lambda G, n: G.N % 2 == 1 and n > 2 and n == G.N - 1,
        "3.4(i)": lambda G, n: G.N % 3 != 0 and n > 3,
        "A228766": lambda G, n: n > 2 and (n % 2 == 1 or G.N % n != 0),
        "3.3(i) v9": lambda G, n: G.N % 2 == 1 and n > 2,
    }
    prob_of = {"3.2": "p33", "3.3(ii) v9": "diff", "3.3(ii) v4-v8": "diff", "3.3(ii) homepage": "diff",
               "3.4(i)": "s2", "A228766": "sum", "3.3(i) v9": "sum"}
    report = {h: [] for h in hyps}
    for ns in groups:
        G = Group(ns)
        tab = {"diff": G.sub, "sum": G.add, "s2": G.s2}
        cos = G.zero_sum_cosets()
        cache, fails = {}, {}
        for n in range(2, G.N + 1):
            for A in itertools.combinations(range(G.N), n):
                cA = G.canon(A)
                for kind in ("diff", "sum", "s2"):
                    if (kind, cA) not in cache:
                        cache[(kind, cA)] = circ_ok(G, list(cA), tab[kind], kind)
                    if not cache[(kind, cA)]:
                        fails.setdefault((kind, n), set()).add(cA)
                for a in A:
                    key = ("p33", G.canon(A, mark=a))
                    if key not in cache:
                        cache[key] = path_ok(G, key[1])
                    rep = key[1]
                    obstructed = frozenset(rep) in cos or (frozenset(rep) - {0}) in cos
                    assert cache[key] == (not obstructed), (G.name, rep)   # 3.2: fails iff coset obstruction
                    if not cache[key]:
                        fails.setdefault(("p33", n), set()).add(rep)
            memcheck()
        for h, test in hyps.items():
            for (kind, n), orbs in sorted(fails.items()):
                if kind == prob_of[h] and test(G, n):
                    show = lambda A: tuple(G.el[i] if len(ns) > 1 else i for i in A)
                    report[h].append("%s n=%d: %d orbit(s), e.g. %s" % (G.name, n, len(orbs), show(min(orbs))))
        print("done", G.name, flush=True)
    for h, lines in report.items():
        print("\n[%s] failures within the hypotheses:" % h, "none" if not lines else "")
        for l in lines:
            print("   ", l)
    assert not report["3.3(ii) homepage"] and not report["3.3(i) v9"]
    print("\n3.2: every failure has A or A minus a_1 a coset of a zero-sum subgroup of order >= 3 (asserted), and conversely.")

# ---------------------------------------------------------------- Part 2: Conjecture 3.8
def prime_powers(n):
    out = []
    for q in range(2, n + 1):
        p = next(d for d in range(2, q + 1) if q % d == 0)
        k, m = 0, q
        while m % p == 0:
            m, k = m // p, k + 1
        if m == 1:
            out.append((q, p, k))
    return out

def build_field(q, p, k):
    """Elements 0..q-1 as base-p digit vectors (digit i = coefficient of x^i); f monic of degree k with x of order q-1."""
    def mulx(a, f):
        c = a[-1]
        return [((a[i - 1] if i else 0) - c * f[i]) % p for i in range(k)]
    enc = lambda a: sum(c * p ** i for i, c in enumerate(a))
    if k == 1:
        g = next(g for g in range(1, q) if len({pow(g, i, q) for i in range(q - 1)}) == q - 1)
        exp, f = [pow(g, i, q) for i in range(q - 1)], None
    else:
        for f in itertools.product(range(p), repeat=k):
            f = list(f)
            if f[0] == 0:
                continue
            a, exp = [1] + [0] * (k - 1), []
            for _ in range(q - 1):
                exp.append(enc(a))
                a = mulx(a, f)
            if len(set(exp)) == q - 1 and enc(a) == 1:
                break
    log = {e: i for i, e in enumerate(exp)}
    def add(a, b):
        r, m = 0, 1
        while a or b:
            r, a, b, m = r + ((a + b) % p) * m, a // p, b // p, m * p
        return r
    addt = [[add(a, b) for b in range(q)] for a in range(q)]
    mul = lambda a, b: 0 if a == 0 or b == 0 else exp[(log[a] + log[b]) % (q - 1)]
    prim = {exp[i] for i in range(q - 1) if gcd(i, q - 1) == 1}
    return addt, mul, prim, f

def gcd(a, b):
    while b:
        a, b = b, a % b
    return a

def check_cycle_independently(q, p, k, f, a0, cyc):
    """Separate arithmetic: polynomials mod f (or integers mod p), primitivity via z^((q-1)/r) != 1."""
    if sorted(cyc) != list(range(1, q)):
        return False
    dec = lambda e: [(e // p ** i) % p for i in range(k)]
    def pmul(a, b):
        if k == 1:
            return [(a[0] * b[0]) % p]
        r = [0] * (2 * k - 1)
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
        for d in range(2 * k - 2, k - 1, -1):          # x^k = -sum f_i x^i
            c, r[d] = r[d], 0
            for i in range(k):
                r[d - k + i] = (r[d - k + i] - c * f[i]) % p
        return r[:k]
    def ppow(a, e):
        r, b = [1] + [0] * (k - 1), a
        while e:
            if e & 1:
                r = pmul(r, b)
            b, e = pmul(b, b), e >> 1
        return r
    one = [1] + [0] * (k - 1)
    rs = [r for r in range(2, q) if (q - 1) % r == 0 and all(r % s for s in range(2, r))]
    if k > 1:                                         # x has order q-1, so the quotient ring is a field
        xx = [0, 1] + [0] * (k - 2)
        if ppow(xx, q - 1) != one or any(ppow(xx, (q - 1) // r) == one for r in rs):
            return False
    for i in range(len(cyc)):
        z = [(u + v) % p for u, v in zip(dec(a0), pmul(dec(cyc[i]), dec(cyc[(i + 1) % len(cyc)])))]
        if not any(z) or any(ppow(z, (q - 1) // r) == one for r in rs):
            return False
    return True

def ham(V, adj, budget, rng):
    n = len(V)
    ix = {v: i for i, v in enumerate(V)}
    A = [sum(1 << ix[w] for w in adj[v]) for v in V]
    if n >= 3 and min(a.bit_count() for a in A) < 2:
        return None, True
    s = min(range(n), key=lambda i: A[i].bit_count())
    path, nodes = [s], [0]
    def feasible(U, end):
        m = U
        while m:
            u = (m & -m).bit_length() - 1
            m &= m - 1
            if (A[u] & U).bit_count() + ((A[u] >> end) & 1) + ((A[u] >> s) & 1) < 2:
                return False
        return True
    def dfs(U, end):
        nodes[0] += 1
        if nodes[0] > budget:
            raise TimeoutError
        if U == 0:
            return (A[end] >> s) & 1 == 1
        if not A[s] & U:
            return False
        cand, m = [], A[end] & U
        while m:
            cand.append((m & -m).bit_length() - 1)
            m &= m - 1
        cand.sort(key=lambda w: ((A[w] & U).bit_count(), rng.random()))
        for w in cand:
            U2 = U & ~(1 << w)
            if U2 and not feasible(U2, w):
                continue
            path.append(w)
            if dfs(U2, w):
                return True
            path.pop()
        return False
    U0 = ((1 << n) - 1) & ~(1 << s)
    try:
        ok = dfs(U0, s)
    except TimeoutError:
        return None, False
    return ([V[i] for i in path] if ok else None), True

def part2(qmax=128):
    sys.setrecursionlimit(10000)
    bad, total = [], 0
    for q, p, k in prime_powers(qmax):
        addt, mul, prim, f = build_field(q, p, k)
        for a0 in range(q):
            V = list(range(1, q))
            adj = {x: [y for y in V if y != x and addt[a0][mul(x, y)] in prim] for x in V}
            rng = random.Random(1000 * q + a0)
            for budget in (20000,) * 6 + (10 ** 7,):
                cyc, complete = ham(V, adj, budget, rng)
                if cyc or complete:
                    break
            if cyc:
                assert check_cycle_independently(q, p, k, f, a0, cyc), (q, a0)
            if q > 7:
                total += 1
                if not cyc:
                    bad.append((q, a0, "non-Hamiltonian (complete search)" if complete else "UNDECIDED"))
            memcheck()
    print("\n3.8: %d instances (q, a0) with 7 < q <= %d; failures: %s" % (total, qmax, bad))
    assert bad == [(11, 2, "non-Hamiltonian (complete search)")]

def prime_check(p):
    G = Group((p,))
    fails = {"diff": set(), "s2": set()}
    for n in range(3, p + 1):
        for rest in itertools.combinations(range(2, p), n - 2):
            A = [0, 1, *rest]
            for kind, tab in (("diff", G.sub), ("s2", G.s2)):
                if not circ_ok(G, A, tab, kind):
                    fails[kind].add(n)
        memcheck()
    print("Z/%d: n with failures: circular differences %s; circular x+2y %s" % (p, sorted(fails["diff"]), sorted(fails["s2"])))

if __name__ == "__main__":
    if PRIMES:
        for p in PRIMES:
            prime_check(p)
        sys.exit(0)
    part1()
    part2()
    print("all checks passed")
