#!/usr/bin/env python3
"""verify.py (stdlib only) for Riemann/Deep/Sun38.

Part A. Sun, arXiv 1309.1679, TeX label Conj3.10 = printed Conjecture 3.8 (2013-09-17):
  for q > 7 and a0 in F_q the graph G_{q,a0} on F_q^* (x ~ y iff x != y and a0 + x*y primitive)
  is Hamiltonian.  Checks every stored certificate in data/certs.txt.gz (all prime powers
  7 < q <= Q, every a0), materializes and checks the full circular permutation for q <= MAT,
  and re-proves the exception (11, 2).  Field arithmetic is rebuilt here from scratch
  (polynomial multiplication mod the stored modulus; the stored generator must have order q-1,
  which also proves that the quotient ring is a field).
  Options: --mat Q (materialize up to Q, default 1000), --recompute LO HI (independent search
  of certificates for all prime powers in [LO, HI], nothing stored), --sample K (recompute K
  random primes from the extension range, default 1).
Part B. Exhaustive re-check of the counterexamples to TeX Conj3.3, Conj3.4(ii), Conj3.5(i)
  (printed 3.2, 3.3(ii), 3.4(i)) and of their minimality (all abelian groups of order <= 8).
"""
import gzip, json, math, os, random, sys, itertools, time

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data')

def is_prime(m):
    if m < 2: return False
    d = 2
    while d * d <= m:
        if m % d == 0: return False
        d += 1
    return True

def prime_factors(m):
    out, d = [], 2
    while d * d <= m:
        if m % d == 0:
            out.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: out.append(m)
    return out

def prime_powers(lo, hi):
    res = []
    for q in range(max(lo, 2), hi + 1):
        p = prime_factors(q)[0]
        k, r = 0, q
        while r % p == 0: r //= p; k += 1
        if r == 1: res.append((q, p, k))
    return res

class GF:
    """F_q = F_p[x]/(f), f = x^k + sum_{i<k} f_i x^i; element <-> integer sum c_i p^i."""
    def __init__(self, p, k, fint):
        assert is_prime(p) and k >= 1
        self.p, self.k, self.q, self.n = p, k, p ** k, p ** k - 1
        self.f = [(fint // p ** i) % p for i in range(k)]
    def dig(self, e):
        d = []
        for _ in range(self.k): d.append(e % self.p); e //= self.p
        return d
    def und(self, d):
        e, m = 0, 1
        for c in d: e += c * m; m *= self.p
        return e
    def add(self, a, b):
        if self.k == 1: return (a + b) % self.p
        return self.und([(x + y) % self.p for x, y in zip(self.dig(a), self.dig(b))])
    def mul(self, a, b):
        p, k = self.p, self.k
        if k == 1: return a * b % p
        A, B = self.dig(a), self.dig(b)
        c = [0] * (2 * k - 1)
        for i, x in enumerate(A):
            if x:
                for j, y in enumerate(B): c[i + j] += x * y
        for i in range(2 * k - 2, k - 1, -1):
            t = c[i] % p
            if t:
                for j in range(k): c[i - k + j] -= t * self.f[j]
            c[i] = 0
        return self.und([v % p for v in c[:k]])
    def power(self, a, m):
        r = 1
        while m:
            if m & 1: r = self.mul(r, a)
            a = self.mul(a, a); m >>= 1
        return r
    def has_order_n(self, g):
        return self.power(g, self.n) == 1 and all(self.power(g, self.n // r) != 1 for r in prime_factors(self.n))
    def set_generator(self, g):
        E, e = [], 1
        for i in range(self.n):
            if i and e == 1: break
            E.append(e); e = self.mul(e, g)
        if e != 1 or len(E) != self.n or len(set(E)) != self.n or 0 in E:
            raise AssertionError('stored generator does not have order q-1 (or modulus not irreducible)')
        self.E, self.L = E, {v: i for i, v in enumerate(E)}
    def prim(self, z):
        return z != 0 and math.gcd(self.L[z], self.n) == 1

def cycle_from_token(n, tok):
    kind, rest = tok[0], tok[1:]
    if kind == 'x': return [int(v) for v in rest.split(',')]
    if kind == 'z':                       # Luo-She: (n-1, 0, n-2, 1, h, h-1, h+1, h-2, h+2, ..., 2, n-3)
        h = (n - 1) // 2
        cyc = [n - 1, 0, n - 2, 1, h]
        for j in range(1, h - 1): cyc += [h - j, h + j]
        return cyc
    a, b = (int(v) for v in rest.split(','))
    if kind == 'p':
        t, d = a, b
        cyc = []
        for m in range(n // 2): cyc += [(m * d) % n, (t - m * d) % n]
        return cyc
    if kind == 'u':
        t, d = a, b
        tp = (t + d) % n
        v0 = tp // 2                      # 2*v0 = tp (tp even), v0 in [0, n/2)
        cyc = []
        for m in range(n // 2): cyc += [(v0 + m * d) % n, (t - v0 - m * d) % n]
        return cyc
    if kind == 'a':
        t, e = a, b
        d = 2 * e % n
        inv2 = (n + 1) // 2
        v0 = (t + d) * inv2 % n
        cyc = []
        for m in range((n + 1) // 2):
            cyc.append((v0 + m * d) % n)
            if m < (n - 1) // 2: cyc.append((t - v0 - m * d) % n)
        return cyc
    raise ValueError(tok)

def check_conditions(F, a0, tok):
    """Lemma conditions (enough by the construction proved in REPORT.md)."""
    n, kind = F.n, tok[0]
    if kind == 'z':
        h = (n - 1) // 2
        return a0 == 0 and n % 2 == 1 and n >= 9 and all(math.gcd(s, n) == 1 for s in (n - 1, n - 2, n - 4, h + 1))
    a, b = (int(v) for v in tok[1:].split(','))
    z = lambda s: F.prim(F.add(a0, F.E[s % n]))
    if kind == 'p':
        return n % 2 == 0 and a % 2 == 1 and 0 < b < n and math.gcd(b, n) == 2 and z(a) and z(a + b)
    if kind == 'u':
        return (n % 2 == 0 and a % 2 == 1 and 0 < b < n and math.gcd(b, n) == 1
                and z(a) and z(a + b) and z(a + b + n // 2))
    if kind == 'a':
        return n % 2 == 1 and 0 < b < n and math.gcd(b, n) == 1 and z(a) and z(a + b) and z(a + 2 * b)
    return False

def check_cycle(F, a0, cyc, direct):
    n = F.n
    if sorted(cyc) != list(range(n)): return False
    E = F.E
    for i in range(n):
        x, y = cyc[i], cyc[(i + 1) % n]
        if x == y: return False
        prod = F.mul(E[x], E[y]) if direct else E[(x + y) % n]
        if not F.prim(F.add(a0, prod)): return False
    return True

def hamiltonian_bruteforce(F, a0):
    """Exhaustive DFS on field elements (small q only). Returns a cycle or None."""
    V = list(range(1, F.q))
    adj = {x: [y for y in V if y != x and F.prim(F.add(a0, F.mul(x, y)))] for x in V}
    start = V[0]
    path, seen = [start], {start}
    def dfs():
        if len(path) == len(V): return start in adj[path[-1]]
        for y in adj[path[-1]]:
            if y not in seen:
                seen.add(y); path.append(y)
                if dfs(): return True
                path.pop(); seen.discard(y)
        return False
    return (path[:] if dfs() else None), adj

def part_a(mat, recompute, sample):
    t0 = time.time()
    meta = json.load(open(os.path.join(DATA, 'meta.json')))
    Q = meta['Q_stored']
    exc_claim = json.load(open(os.path.join(DATA, 'exceptions.json')))
    claimed = sorted((e['q'], e['a0']) for e in exc_claim)
    seen_q, n_cert, n_mat, found_exc = [], 0, 0, []
    with gzip.open(os.path.join(DATA, 'certs.txt.gz'), 'rt') as fi:
        while True:
            hdr = fi.readline()
            if not hdr: break
            _, q, p, k, fint, g = hdr.split(); q, p, k, fint, g = map(int, (q, p, k, fint, g))
            toks = fi.readline().split()
            assert q == p ** k and is_prime(p) and len(toks) == q, ('bad field record', q)
            F = GF(p, k, fint); F.set_generator(g)
            for a0, tok in enumerate(toks):
                if tok[0] == 'n':
                    found_exc.append((q, a0)); continue
                if tok[0] == 'x':
                    assert check_cycle(F, a0, cycle_from_token(F.n, tok), True), ('bad cycle', q, a0)
                else:
                    assert check_conditions(F, a0, tok), ('bad certificate', q, a0, tok)
                    if q <= mat:
                        assert check_cycle(F, a0, cycle_from_token(F.n, tok), k == 1), ('bad cycle', q, a0, tok)
                        n_mat += 1
                n_cert += 1
            seen_q.append(q)
    expect = [q for q, _, _ in prime_powers(8, Q)]
    assert seen_q == expect, 'field coverage mismatch'
    assert sorted(found_exc) == claimed, ('exception list mismatch', found_exc, claimed)
    # re-prove every exception
    for q, a0 in claimed:
        p = prime_factors(q)[0]; k = round(math.log(q, p))
        assert k == 1, 'only prime-field exceptions expected'
        F = GF(p, 1, 0)
        g = next(h for h in range(2, q) if all(pow(h, F.n // r, q) != 1 for r in prime_factors(F.n)))
        F.set_generator(g)
        cyc, adj = hamiltonian_bruteforce(F, a0)
        assert cyc is None, ('exception has a Hamilton cycle!', q, a0)
        deg2 = {x: tuple(sorted(v)) for x, v in adj.items() if len(v) == 2}
        twins = [(x, y) for x in deg2 for y in deg2 if x < y and deg2[x] == deg2[y]]
        print('  exception (q, a0) = (%d, %d): exhaustive search finds no Hamilton cycle; degree-2 twins %s with N = %s'
              % (q, a0, twins, [deg2[x] for x, _ in twins]))
    print('Part A: %d fields 7 < q <= %d, %d Hamiltonian (q, a0) certified (%d materialized and edge-checked), '
          'exceptions %s  [%.1f s]' % (len(seen_q), Q, n_cert, n_mat, claimed, time.time() - t0))
    # extension range (not stored): summary + independent recomputation of a sample / a range
    ranges = []
    if recompute: ranges.append(tuple(recompute))
    ext = meta.get('extension')
    if ext:
        lo, hi = ext['lo'], ext['hi']
        rows = {}
        with gzip.open(os.path.join(DATA, 'ext_primes.txt.gz'), 'rt') as fi:
            for line in fi:
                v = list(map(int, line.split()))
                rows[v[0]] = v[1:]
        want = [q for q in range(lo, hi + 1) if is_prime(q)]
        assert all(q in rows and rows[q][0] == q and rows[q][1] == 0 for q in want), 'ext_primes coverage'
        pp = json.load(open(os.path.join(DATA, 'ext_pp_stats.json')))
        wantpp = [q for q, p, k in prime_powers(lo, hi) if k > 1]
        assert sorted(map(int, pp)) == wantpp and all(sum(pp[str(q)].values()) == q and set(pp[str(q)]) <= set('pauz')
                                                       for q in wantpp), 'ext_pp coverage'
        print('  extension %d..%d (summaries, certificates not stored): %d primes via ext.c, %d prime powers via gen.py, '
              'every a0 certified, no exception' % (lo, hi, len(want), len(wantpp)))
    if '--recompute-pp' in sys.argv:
        t1 = time.time(); cnt = 0
        for q, p, k in prime_powers(ext['lo'], ext['hi']):
            if k > 1: cnt += recompute_field(q, p, k)
        print('  recomputed all prime powers (k > 1) in [%d, %d]: %d graphs OK [%.1f s]' % (ext['lo'], ext['hi'], cnt, time.time() - t1))
    if ext and sample:
        rng = random.Random()
        cands = [q for q in range(ext['lo'], ext['hi'] + 1) if is_prime(q)]
        for q in rng.sample(cands, sample): ranges.append((q, q))
    for lo, hi in ranges:
        t1 = time.time(); cnt = 0
        for q, p, k in prime_powers(max(lo, 8), hi):
            cnt += recompute_field(q, p, k)
        print('  recomputed certificates for all a0, prime powers in [%d, %d]: %d graphs OK [%.1f s]'
              % (lo, hi, cnt, time.time() - t1))

def recompute_field(q, p, k):
    """Independent search (own primitive polynomial / root) for a certificate for every a0."""
    n = q - 1
    if k == 1:
        g = next(h for h in range(2, q) if all(pow(h, n // r, q) != 1 for r in prime_factors(n)))
        F = GF(p, 1, 0); F.set_generator(g)
    else:
        for fint in range(1, q):
            if fint % p == 0: continue
            F = GF(p, k, fint)
            if F.has_order_n(p):              # x has order q-1, so f is primitive (hence irreducible)
                F.set_generator(p); break
    unitn = [math.gcd(s, n) for s in range(n)]
    for a0 in range(q):
        inT = lambda s: F.prim(F.add(a0, F.E[s % n]))
        tok = None
        if n % 2 == 0:
            odds = []
            for s in range(1, n, 2):
                if inT(s):
                    hit = next((t for t in odds if unitn[s - t] == 2), None)
                    if hit is not None: tok = 'p%d,%d' % (hit, s - hit); break
                    odds.append(s)
                    if len(odds) > 200: break
        else:
            lim = min(n, 400)
            tok = next(('a%d,%d' % (t, e) for t in range(lim) if inT(t) for e in range(1, lim)
                        if unitn[e] == 1 and inT(t + e) and inT(t + 2 * e)), None)
            if tok is None and a0 == 0 and n >= 9: tok = 'z'
        if tok is None:
            raise AssertionError('no P2/AP3 certificate for (q, a0) = (%d, %d); use gen.py search' % (q, a0))
        assert check_conditions(F, a0, tok)
    return q

# ------------------------------------------------------------------ Part B
def abelian_groups(m):
    """All abelian groups of order m as lists of cyclic moduli (prime-power factors)."""
    def parts(e, mx=None):
        if e == 0: yield []; return
        mx = e if mx is None else mx
        for x in range(min(e, mx), 0, -1):
            for rest in parts(e - x, x): yield [x] + rest
    fac = []
    r = m
    for p in prime_factors(m):
        e = 0
        while r % p == 0: r //= p; e += 1
        fac.append((p, e))
    out = []
    for combo in itertools.product(*[list(parts(e)) for _, e in fac]):
        mods, two = [], 0
        for (p, _), part in zip(fac, combo):
            mods += [p ** x for x in part]
            if p == 2: two = len(part)
        out.append((mods, two <= 1))
    return out or [([], True)]

def group_elems(mods):
    return list(itertools.product(*[range(md) for md in mods])) if mods else [()]

def name(mods):
    return 'x'.join('Z/%d' % md for md in mods) if mods else 'trivial'

def search_orders(A, mods, mode, first=None):
    """mode 'path': permutation with b1 = first, the n-1 differences b_i - b_{i+1} distinct.
       mode 'cdiff': circular order with the n differences a_i - a_{i+1} distinct.
       mode 'csum2': circular order with the n sums a_i + 2 a_{i+1} distinct."""
    sub = lambda x, y: tuple((a - b) % md for a, b, md in zip(x, y, mods))
    s2 = lambda x, y: tuple((a + 2 * b) % md for a, b, md in zip(x, y, mods))
    op = sub if mode in ('path', 'cdiff') else s2
    A = list(A)
    start = first if mode == 'path' else A[0]
    rest = [a for a in A if a != start]
    used, path = set(), [start]
    def dfs(remaining):
        if not remaining:
            if mode == 'path': return True
            v = op(path[-1], start)
            return v not in used
        for i, y in enumerate(remaining):
            v = op(path[-1], y)
            if v in used: continue
            used.add(v); path.append(y)
            if dfs(remaining[:i] + remaining[i + 1:]): return True
            path.pop(); used.discard(v)
        return False
    return dfs(rest)

def census(maxorder, which):
    """Failures within the stated hypotheses, per (group, n)."""
    fails = []
    for m in range(1, maxorder + 1):
        for mods, syl2cyc in abelian_groups(m):
            G = group_elems(mods)
            for n in range(1, m + 1):
                if which == '3.3':
                    if n < 2 or not (m % n != 0 or (n % 2 == 0 and syl2cyc)): continue
                elif which == '3.4ii':
                    if not (m % 2 == 1 and n > 2 and m % n != 0): continue
                elif which == '3.5i':
                    if not (n > 3 and m % 3 != 0): continue
                bad = []
                for A in itertools.combinations(G, n):
                    if which == '3.3':
                        for a1 in A:
                            if not search_orders(A, mods, 'path', a1): bad.append((A, a1))
                    elif not search_orders(A, mods, 'cdiff' if which == '3.4ii' else 'csum2'):
                        bad.append(A)
                if bad: fails.append((m, name(mods), n, len(bad), bad[0]))
    return fails

def part_b():
    t0 = time.time()
    # Conj3.3 (printed 3.2): G = Z/6, (a1..a4) = (0, 1, 3, 5)
    A = [(0,), (1,), (3,), (5,)]
    assert 6 % 4 != 0 and not search_orders(A, [6], 'path', (0,))
    f33 = census(8, '3.3')
    first = min((m, n) for m, _, n, _, _ in f33)
    assert first == (6, 4), f33
    print('Part B, TeX Conj3.3 (printed 3.2): Z/6, (0,1,3,5), a1 = 0 has no valid permutation; '
          'smallest failing (|G|, n) = %s; failures |G| <= 8: %s' % (first, [(g, n, c) for _, g, n, c, _ in f33]))
    # Conj3.4(ii) (printed 3.3(ii) of arXiv v9): Z/5, A = {0, 1, 2}
    assert not search_orders([(0,), (1,), (2,)], [5], 'cdiff')
    f34 = census(9, '3.4ii')
    first = min((m, n) for m, _, n, _, _ in f34)
    assert first == (5, 3), f34
    print('Part B, TeX Conj3.4(ii) (printed 3.3(ii), v9 wording): Z/5, {0,1,2} fails; smallest failing (|G|, n) = %s; '
          'failures |G| <= 9 (odd): %s' % (first, [(g, n, c) for _, g, n, c, _ in f34]))
    # Conj3.5(i) (printed 3.4(i)): Z/7, n = 5 (all 21 subsets)
    G7 = [(i,) for i in range(7)]
    bad7 = [A for A in itertools.combinations(G7, 5) if not search_orders(A, [7], 'csum2')]
    assert len(bad7) == 21
    f35 = census(8, '3.5i')
    first = min((m, n) for m, _, n, _, _ in f35)
    assert first == (7, 5), f35
    print('Part B, TeX Conj3.5(i) (printed 3.4(i)): all 21 five-subsets of Z/7 fail; smallest failing (|G|, n) = %s; '
          'failures |G| <= 8: %s  [%.1f s]' % (first, [(g, n, c, b) for _, g, n, c, b in f35], time.time() - t0))

def main():
    args = sys.argv[1:]
    mat = int(args[args.index('--mat') + 1]) if '--mat' in args else 1000
    rec = [int(args[args.index('--recompute') + 1]), int(args[args.index('--recompute') + 2])] if '--recompute' in args else None
    sample = int(args[args.index('--sample') + 1]) if '--sample' in args else 1
    if '--only-b' not in args: part_a(mat, rec, sample)
    if '--only-a' not in args: part_b()
    print('ALL CHECKS PASSED')

if __name__ == '__main__':
    main()
