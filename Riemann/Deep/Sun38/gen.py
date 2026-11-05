#!/usr/bin/env python3
"""sun38 generator: Hamiltonicity of G_{q,a0} (Sun 1309.1679, TeX Conj3.10 = printed Conj 3.8).

G_{q,a0}: vertices F_q^*, x ~ y iff x != y and a0 + x*y is primitive.
In log coordinates (x = g^i, g a fixed primitive element, n = q-1) this is the Cayley sum graph
on Z/n: i ~ j iff i != j and i + j in T, where T = T(a0) = {s : a0 + g^s primitive}.

Certificates (cycle written as logs v_0..v_{n-1}; circular permutation a_k = g^{v_{k-1}}):
  p t,d  (n even) t odd, t in T, t+d in T, gcd(d,n)=2:  v_{2m}=m*d, v_{2m+1}=t-m*d.
  u t,d  (n even) t odd, t in T, t'=t+d even in T, gcd(d,n)=1, t'+n/2 in T:
         v_0 = t'/2 (the solution in [0,n/2)), v_{2m}=v_0+m*d, v_{2m+1}=t-v_0-m*d.
  a t,e  (n odd)  t, t+e, t+2e in T, gcd(e,n)=1:  d=2e, v_0=(t+d)/2 mod n,
         v_{2m}=v_0+m*d (m<=(n-1)/2), v_{2m+1}=t-v_0-m*d (m<=(n-3)/2).
  z      (a0 = 0, n = 2h+1 odd, n >= 9) Luo-She cycle (n-1, 0, n-2, 1, h, h-1, h+1, ..., 2, n-3);
         adjacent sums n-1, n-2, n-4, h+1, all coprime to n.
  x v0,v1,...  explicit cycle (logs) found by search.
  n<tag> non-Hamiltonian (details in exceptions.json).
Usage: gen.py QLO QHI OUTFILE   (all prime powers QLO <= q <= QHI)
"""
import sys, math, random, time, gzip, json, os

def is_prime(m):
    if m < 2: return False
    i = 2
    while i * i <= m:
        if m % i == 0: return False
        i += 1
    return True

def prime_factors(m):
    f, d = [], 2
    while d * d <= m:
        if m % d == 0:
            f.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: f.append(m)
    return f

def prime_powers(lo, hi):
    out = []
    for q in range(max(lo, 2), hi + 1):
        for p in prime_factors(q)[:1]:
            k, r = 0, q
            while r % p == 0: r //= p; k += 1
            if r == 1: out.append((q, p, k))
    return out

class Field:
    """F_q, q=p^k; element encoding e = sum c_i p^i  <->  sum c_i x^i (mod f)."""
    def __init__(self, p, k):
        self.p, self.k, self.q = p, k, p ** k
        q, n = self.q, self.q - 1
        self.n = n
        if k == 1:
            fs = prime_factors(n)
            g = 2 if q > 2 else 1
            while any(pow(g, n // r, q) == 1 for r in fs): g += 1
            self.g, self.fint = g, 0
            exp = [1] * n
            for i in range(1, n): exp[i] = exp[i - 1] * g % q
        else:
            # smallest monic f = x^k + sum_{i<k} c_i x^i (fint = sum c_i p^i) with x of order n
            pk1 = p ** (k - 1)
            for fint in range(1, q):
                if fint % p == 0: continue          # c_0 != 0
                negf = [(-((fint // p ** i) % p)) % p for i in range(k)]
                exp = [1] * n
                e, ok = 1, True
                for i in range(1, n):
                    top, rest = e // pk1, e % pk1
                    e = rest * p
                    if top:
                        # add top * (-f_lower) digitwise
                        r, m, ee = 0, 1, e
                        for j in range(k):
                            r += ((ee % p + top * negf[j]) % p) * m
                            ee //= p; m *= p
                        e = r
                    if e == 1: ok = False; break
                    exp[i] = e
                if ok:
                    # x^n must be 1
                    top, rest = e // pk1, e % pk1
                    e2 = rest * p
                    if top:
                        r, m, ee = 0, 1, e2
                        for j in range(k):
                            r += ((ee % p + top * negf[j]) % p) * m
                            ee //= p; m *= p
                        e2 = r
                    if e2 == 1:
                        self.fint, self.g = fint, p     # generator x (encoding p)
                        break
            else:
                raise RuntimeError('no primitive polynomial')
        self.exp = exp
        log = [-1] * q
        for i, e in enumerate(exp): log[e] = i
        assert all(l >= 0 for l in log[1:]), 'exp table not a bijection'
        self.log = log

    def add(self, a, b):
        p = self.p
        if self.k == 1: return (a + b) % p
        if p == 2: return a ^ b
        r, m = 0, 1
        while a or b:
            r += ((a % p + b % p) % p) * m
            a //= p; b //= p; m *= p
        return r

def find_cert(F, a0, unit, gcdv):
    """Return certificate token or None (T accessed lazily)."""
    n, exp, log = F.n, F.exp, F.log
    memo = {}
    def inT(s):
        s %= n
        v = memo.get(s)
        if v is None:
            z = F.add(a0, exp[s])
            v = z != 0 and unit[log[z]]
            memo[s] = v
        return v
    if n % 2 == 0:
        odds = []
        for s in range(1, n, 2):
            if inT(s):
                for t in odds:
                    if gcdv[s - t] == 2: return 'p%d,%d' % (t, s - t)
                odds.append(s)
                if len(odds) > 60: break
        # U2 fallback
        h = n // 2
        Tl = [s for s in range(n) if inT(s)]
        Te = [s for s in Tl if s % 2 == 0]
        To = [s for s in Tl if s % 2 == 1]
        for tp in Te:
            if not inT(tp + h): continue
            for t in To:
                d = (tp - t) % n
                if gcdv[d] == 1: return 'u%d,%d' % (t, d)
        return None
    else:
        big = n > 2000
        for t in range(n):
            if big and t > 400: break
            if not inT(t): continue
            for e in range(1, n):
                if big and e > 400: break
                if gcdv[e] == 1 and inT(t + e) and inT(t + 2 * e):
                    return 'a%d,%d' % (t, e)
        if a0 == 0 and n >= 9:
            return 'z'                  # Luo-She explicit cycle for a0 = 0, n odd (sums n-1, n-2, n-4, (n+1)/2)
        return None

# ---------------------------------------------------------------- search for leftovers
def graph_adj(F, a0, unit):
    n = F.n
    T = [s for s in range(n) if (lambda z: z != 0 and unit[F.log[z]])(F.add(a0, F.exp[s]))]
    adj = [sorted({(s - v) % n for s in T} - {v}) for v in range(n)]
    return T, adj

def posa(n, adj, rng, max_steps):
    adjset = [set(a) for a in adj]
    start = rng.randrange(n)
    path, pos = [start], [-1] * n
    pos[start] = 0
    for _ in range(max_steps):
        v = path[-1]
        cand = [w for w in adj[v] if pos[w] < 0]
        L = len(path)
        if cand:
            if n > 3000: w = rng.choice(cand)
            else: w = min(cand, key=lambda w: (sum(1 for x in adj[w] if pos[x] < 0), rng.random()))
            pos[w] = L; path.append(w); continue
        if L == n and path[0] in adjset[v]: return path
        ch = [w for w in adj[v] if 0 <= pos[w] < L - 2]
        if not ch: return None
        i = pos[rng.choice(ch)]
        seg = path[i + 1:]; seg.reverse(); path[i + 1:] = seg
        for j in range(i + 1, L): pos[path[j]] = j
    return None

def exact(n, adj, node_limit):
    """Exhaustive Hamilton-cycle search with pruning. Returns (status, cycle, nodes)."""
    am = [0] * n
    for v in range(n):
        for w in adj[v]: am[v] |= 1 << w
    if min(len(a) for a in adj) < 2: return ('non', None, 0)
    s = min(range(n), key=lambda v: (len(adj[v]), v))
    nodes = [0]
    path = [s]
    def feasible(U, end):
        base = U | (1 << end) | (1 << s)
        fe = fs = 0
        x = U
        while x:
            b = x & -x; u = b.bit_length() - 1; x ^= b
            a = am[u] & base
            c = a.bit_count()
            if c < 2: return False
            if c == 2:
                if (a >> end) & 1: fe += 1
                if (a >> s) & 1: fs += 1
        lim = 2 if end == s else 1
        if fe > lim or fs > lim: return False
        if U and not (am[end] & U): return False
        if U and not (am[s] & U): return False
        seen = fr = 1 << end
        while fr:
            nb = 0
            y = fr
            while y:
                b = y & -y; y ^= b
                nb |= am[b.bit_length() - 1]
            nb &= base & ~seen
            seen |= nb; fr = nb
        return seen == base
    def dfs(U, end):
        nodes[0] += 1
        if nodes[0] > node_limit: raise TimeoutError
        if not U: return (am[end] >> s) & 1 == 1
        base = U | (1 << end) | (1 << s)
        cands = []
        x = am[end] & U
        forced = None
        while x:
            b = x & -x; u = b.bit_length() - 1; x ^= b
            c = (am[u] & base).bit_count()
            if c == 2 and U != b: forced = u
            cands.append((c, u))
        if forced is not None: cands = [(0, forced)]
        cands.sort()
        for _, u in cands:
            U2 = U & ~(1 << u)
            if feasible(U2, u):
                path.append(u)
                if dfs(U2, u): return True
                path.pop()
        return False
    full = (1 << n) - 1
    try:
        U0 = full & ~(1 << s)
        ok = feasible(U0, s) and dfs(U0, s)
    except TimeoutError:
        return ('undecided', None, nodes[0])
    return ('ham', path[:], nodes[0]) if ok else ('non', None, nodes[0])

def check_cycle_logs(n, T, cyc):
    Ts = set(T)
    return sorted(cyc) == list(range(n)) and all(((cyc[i] + cyc[(i + 1) % n]) % n) in Ts and cyc[i] != cyc[(i + 1) % n] for i in range(n))

def nonham_reason(n, adj):
    if min(len(a) for a in adj) < 2: return 'deg<2'
    forced = {}
    for v in range(n):
        if len(adj[v]) == 2:
            for w in adj[v]:
                forced.setdefault(v, set()).add(w); forced.setdefault(w, set()).add(v)
    if any(len(s) > 2 for s in forced.values()): return 'forced3'
    # forced short cycle?
    seen = set()
    for v in forced:
        if v in seen: continue
        comp, stack = [], [v]
        seen.add(v)
        while stack:
            x = stack.pop(); comp.append(x)
            for y in forced[x]:
                if y not in seen: seen.add(y); stack.append(y)
        if all(len(forced[x]) == 2 for x in comp) and len(comp) < n: return 'forcedcycle%d' % len(comp)
    return 'search'

def main():
    lo, hi, out = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    nostore = '--nostore' in sys.argv
    onlypp = '--onlypp' in sys.argv          # skip prime fields (done by ext.c)
    rng = random.Random(20261105)
    t0 = time.time()
    exc = []
    stats = {}
    with gzip.open(out, 'wt') as fo:
        for q, p, k in prime_powers(max(lo, 8), hi):
            if onlypp and k == 1: continue
            F = Field(p, k)
            n = F.n
            unit = bytearray(1 if math.gcd(s, n) == 1 else 0 for s in range(n))
            gcdv = [math.gcd(d, n) for d in range(n)]
            toks = []
            cnt = {}
            for a0 in range(q):
                tok = find_cert(F, a0, unit, gcdv)
                if tok is None:
                    T, adj = graph_adj(F, a0, unit)
                    cyc = None
                    if min(len(a) for a in adj) >= 2:
                        for _ in range(20):
                            cyc = posa(n, adj, rng, 50 * n)
                            if cyc: break
                    if cyc is None:
                        st, cyc, nodes = exact(n, adj, 10 ** 7)
                        if st != 'ham':
                            reason = nonham_reason(n, adj)
                            exc.append({'q': q, 'a0': a0, 'status': st, 'nodes': nodes, 'reason': reason,
                                        'adj_logs': adj})
                            tok = 'n' + st
                    if cyc is not None:
                        assert check_cycle_logs(n, T, cyc)
                        tok = 'x' + ','.join(map(str, cyc))
                cnt[tok[0]] = cnt.get(tok[0], 0) + 1
                toks.append(tok)
            if not nostore:
                fo.write('F %d %d %d %d %d\n' % (q, p, k, F.fint, F.g))
                fo.write(' '.join(toks) + '\n')
            stats[q] = cnt
            if q < 200 or q % 100 < 3 or k > 1:
                print(q, p, k, cnt, round(time.time() - t0, 1), flush=True)
    base = out.replace('.txt.gz', '')
    with open(base + '.exc.json', 'w') as fe:
        json.dump(exc, fe)
    with open(base + '.stats.json', 'w') as fs:
        json.dump(stats, fs)
    print('done', lo, hi, 'exceptions:', [(e['q'], e['a0'], e['status'], e['reason'], e['nodes']) for e in exc],
          'time', round(time.time() - t0, 1))

if __name__ == '__main__':
    main()
