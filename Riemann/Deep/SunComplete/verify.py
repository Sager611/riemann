#!/usr/bin/env python3
"""verify.py -- re-checks everything in NOTE.md (Riemann/Deep/SunComplete), < 60 s.

  1. Lemma 6 (uniform criterion (CS)) at the tabulated points (exact rational arithmetic),
     and the tail estimate used for omega(q-1) > 60.
  2. Proposition 4 (exact per-q criterion) for every prime power 3 <= q < 30140; the
     failing q must be exactly the exceptional sets EA, EB below.
  3. Certificates for EA (Theorem A) and EB (Theorem B, every a0 != 0), from certA.json and
     certB.json.gz: all certificate conditions are checked with independent arithmetic, and
     full Hamilton cycles are rebuilt and checked for every (A)-instance and for a random sample
     of (B)-instances; (q,a0) = (11,2) is shown non-Hamiltonian by exhaustive search.
  4. a0 = 0 (Luo-She construction): exponent sequences for all prime powers 4 < q < 3000.
  5. (optional, --cross) recompiles crossA.c / crossB.c and compares with crossA.out / crossB.out.
All inequalities use fractions.Fraction; sqrt(q) is replaced by a rational upper bound.
"""
import sys, json, gzip, math, random, itertools, time, os, subprocess
from fractions import Fraction as Fr
from math import gcd, isqrt

T_START = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
def say(*a): print(*a, flush=True)

# ---------------------------------------------------------------- basic arithmetic
def primes_upto(N):
    s = bytearray([1])*(N+1); s[0] = s[1] = 0
    for i in range(2, isqrt(N)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(N+1) if s[i]]
PR = primes_upto(5000)
def factor_distinct(n):
    ps = []; d = 2
    while d*d <= n:
        if n % d == 0:
            ps.append(d)
            while n % d == 0: n //= d
        d += 1
    if n > 1: ps.append(n)
    return ps
def phi(n):
    r = n
    for p in factor_distinct(n): r = r//p*(p-1)
    return r
def prime_power(q):
    ps = factor_distinct(q)
    if len(ps) != 1: return None
    p = ps[0]; k = 0
    while q % p == 0: q //= p; k += 1
    return p, k
def sqrt_ub(q):
    r = isqrt(q)
    if r*r == q: return Fr(r)
    return Fr(isqrt(q*10**12)+1, 10**6)
def theta_of(ps):
    t = Fr(1)
    for p in ps: t *= Fr(p-1, p)
    return t

# ---------------------------------------------------------------- 1. Lemma 6 (CS)
def sieve_params(core, siev):
    s = len(siev); Wk = 2**len(core); thk = theta_of(core)
    if s == 0: return Wk, thk, Fr(1), Fr(1)
    delta = 1 - sum(Fr(1, p) for p in siev)
    if delta <= 0: return None
    return Wk, thk, delta, Fr(s-1)/delta + 2

def alpha_r(r): return theta_of(PR[:r])
def P_r(r):
    x = 1
    for p in PR[:r]: x *= p
    return x

def CS_holds(case, r, t, q):
    sp = sieve_params(PR[:t], PR[t:r])
    if sp is None: return False
    Wk, thk, dl, Dl = sp
    tau = thk*dl
    sq = sqrt_ub(q)
    E = (Wk*Dl - 1)*sq + 1
    if case == 'A': n, beta = q, 0
    else:
        n, beta = q-1, 1
        if not Fr(isqrt(q)) >= 1/tau: return False
    Nt = Fr(q-1) - Fr(beta)/tau
    Ap = 1/tau + E
    D = alpha_r(r)*(q-1) - 1
    if not (D > Ap and 2*Ap < n and Nt > E): return False
    return (D - Ap)*Nt*Nt*(1 - 2*Ap/n) > E*E*D

# (r, t, q0): (CS) holds at q0 with core = t smallest primes; found by bisection (monotone in q).
CS_TABLE = {
 'A': [(1,0,23),(2,1,146),(3,2,908),(4,2,3904),(5,2,11825),(6,2,30138),(7,2,65775),(8,2,133389),(9,2,249816),(10,2,432449),(11,2,730603),(12,2,1177203)] + [(r,3,None) for r in range(13,35)] + [(r,4,None) for r in range(35,61)],
 'B': [(1,0,25),(2,1,147),(3,2,909),(4,2,3905),(5,2,11827),(6,2,30139),(7,2,65776),(8,2,133390),(9,2,249817),(10,2,432451),(11,2,730605),(12,2,1177204)] + [(r,3,None) for r in range(13,35)] + [(r,4,None) for r in range(35,61)],
}
QCUT = 30140   # every prime power q >= QCUT satisfies (CS); every q < QCUT is checked individually
def check_lemma6():
    for case in 'AB':
        for (r, t, q0) in CS_TABLE[case]:
            if q0 is None: q0 = P_r(r) + 1
            assert CS_holds(case, r, t, q0), ('CS fails', case, r, t, q0)
            if r <= 6: assert q0 <= QCUT, (case, r, q0)
            else: assert q0 <= P_r(r) + 1, (case, r, q0)
    # tail r > 60 (NOTE, Lemma 6(c)): with q >= P_r + 1 >= (r+1)!, the quantity
    #   rho_r = 2^r (r+1)^2 / sqrt((r+1)!)  bounds A'/D etc.; need rho_61 tiny and rho decreasing.
    r = 61
    fact = math.factorial(r+1)
    rho = Fr(2**r*(r+1)**2) / Fr(isqrt(fact))
    assert rho < Fr(1, 10**20), float(rho)
    assert all(p >= j+2 for j, p in enumerate(PR[:200])), 'p_j >= j+1'
    say("[1] Lemma 6 (CS): verified at %d tabulated points (r <= 60) + tail bound; uniform range q >= %d" % (2*60, QCUT))

# ---------------------------------------------------------------- 2. Proposition 4 per q
def case_consts(case, q, Phi):
    if case == 'A':
        eps = 0 if q % 2 == 0 else 1
        return q, eps, 0, Phi - eps
    return q-1, 1, 1, Phi - 2
def crit_exact(case, q, Phi, core, siev):
    sp = sieve_params(core, siev)
    if sp is None: return False
    Wk, thk, dl, Dl = sp
    tau = thk*dl
    E = (Wk*Dl - 1)*sqrt_ub(q) + 1
    n, eps, beta, dmin = case_consts(case, q, Phi)
    Nt = Fr(q-1) - Fr(beta)/tau
    if not Nt > E: return False
    astar = n*(Fr(eps)/tau + E)/(Nt + E)
    C = math.floor(astar) - 1
    a1 = dmin - C + 1
    if a1 < 1: return False
    F = lambda a: Nt*Nt*a*(n - a - C) - E*E*(a + C)*(n - a)
    if Fr(a1) <= Fr(n, 2):
        if not (F(Fr(a1)) > 0 and F(Fr(n, 2)) > 0): return False
    return True
def crit_best(case, q):
    ps = factor_distinct(q-1); Phi = phi(q-1)
    for m in range(len(ps)+1):
        for S in itertools.combinations(ps, m):
            if crit_exact(case, q, Phi, [p for p in ps if p not in S], list(S)): return True
    return False

EA = [3,4,5,7,9,11,13,16,19,23,25,27,29,31,37,41,43,47,49,61,67,71,73,79,103,121,127,131,139,151,157,169,181,191,199,211,223,229,241,271,277,283,307,313,331,337,361,379,421,463,541,547,571,631,661,691,841,859,967,991,1021,1051,1093,1171,1231,1291,1321,1381,1429,1471,1531,1681,1951,2311,2731,3571,4621,4831,6091]
EB = sorted(EA + [8, 53, 1123, 1741])
def check_prop4():
    pp = sorted(x for p in primes_upto(QCUT) for x in [p**i for i in range(1, 20)] if x < QCUT)
    fa = [q for q in pp if q >= 3 and not crit_best('A', q)]
    fb = [q for q in pp if q >= 3 and not crit_best('B', q)]
    assert fa == EA, (fa, EA)
    assert fb == EB, (fb, EB)
    say("[2] Prop. 4 checked for all %d prime powers 3 <= q < %d: exceptions |EA| = %d, |EB| = %d (max %d)" % (len(pp)-1, QCUT, len(fa), len(fb), max(fb)))

# ---------------------------------------------------------------- finite fields (honest arithmetic)
class Field:
    """F_q = F_p[x]/(f), elements coded by base-p digits (coefficient of x^i = digit i).
    For k >= 2 the power table of x is built by honest multiplication by x and must have
    exactly q-1 distinct nonzero entries with x^(q-1) = 1: this proves f irreducible and x primitive."""
    def __init__(self, q, f=None, convention='C'):
        self.q = q; self.p, self.k = prime_power(q); self.n = q-1
        self.nps = factor_distinct(q-1)
        p, k = self.p, self.k
        self.pw = [p**i for i in range(k+1)]
        if k == 1:
            self.f = None; return
        if f is None: f = self.find_poly(convention)
        self.f = f
        exp = self.power_table(f)
        assert exp is not None, ('not a primitive polynomial', q, f)
        self.exp = exp
        self.log = [None]*q
        for i, e in enumerate(exp): self.log[e] = i
    def mulx(self, a, f):
        p, k = self.p, self.k
        d = [(a // self.pw[i]) % p for i in range(k)]
        top = d[k-1]; nd = [0] + d[:k-1]
        return sum(((nd[i] - top*f[i]) % p)*self.pw[i] for i in range(k))
    def power_table(self, f):
        exp = [1]; seen = {1}
        a = 1
        for i in range(1, self.q-1):
            a = self.mulx(a, f)
            if a == 0 or a in seen: return None
            seen.add(a); exp.append(a)
        return exp if self.mulx(a, f) == 1 else None
    def find_poly(self, convention):
        p, k = self.p, self.k
        for c in range(p**k):
            f = [(c // p**i) % p for i in range(k)] + [1]
            if f[0] == 0: continue
            if self.power_table(f) is not None: return f
    def add(self, a, b):
        if self.k == 1: return (a + b) % self.p
        p = self.p
        return sum((((a // w) + (b // w)) % p)*w for w in self.pw[:self.k])
    def neg(self, a):
        if self.k == 1: return (-a) % self.p
        p = self.p
        return sum(((-(a // w)) % p)*w for w in self.pw[:self.k])
    def mul(self, a, b):
        if a == 0 or b == 0: return 0
        if self.k == 1: return a*b % self.p
        return self.exp[(self.log[a] + self.log[b]) % self.n]
    def pow(self, a, e):
        if self.k == 1: return pow(a, e, self.p)
        if a == 0: return 0 if e else 1
        return self.exp[(self.log[a]*e) % self.n]
    def inv(self, a): return self.pow(a, self.n - 1)
    def is_prim(self, z):
        if z == 0: return False
        return all(self.pow(z, self.n // l) != 1 for l in self.nps)
    def is_square(self, z):          # q odd, z != 0
        return self.pow(z, self.n // 2) == 1
    def order_is(self, z, m):        # multiplicative order of z equals m (m | q-1)
        return self.pow(z, m) == 1 and all(self.pow(z, m // l) != 1 for l in factor_distinct(m))

# ---------------------------------------------------------------- 3. certificates
def check_cycle_A(F, cyc):
    q = F.q
    assert sorted(cyc) == list(range(q))
    for i in range(q): assert F.is_prim(F.add(cyc[i], cyc[(i+1) % q])), (q, i)
def check_cycle_B(F, a0, cyc):
    q = F.q; N = q-1
    assert sorted(cyc) == list(range(1, q))
    for i in range(N): assert F.is_prim(F.add(a0, F.mul(cyc[i], cyc[(i+1) % N]))), (q, a0, i)

def she_cycle(p, x):     # She's construction: x, x+1, x+2 primitive mod p
    m = (p-1)//2; h = x*pow(2, -1, p) % p; seq = {}
    for k in range(m+1): seq[2*k+1] = (h + 2*k + 1) % p
    for k in range(1, m+1): seq[2*k] = (h - 2*k + 1) % p
    return [seq[i] for i in range(1, p+1)]

def alt_cycle(F, u, r):  # q odd: a_{2j+1} = r^j, a_{2j} = u r^{1-j}
    m = F.n // 2; cyc = []
    for j in range(m):
        cyc.append(F.pow(r, j)); cyc.append(F.mul(u, F.pow(r, (1 - (j+1)) % F.n)))
    return cyc
def geo_cycle(F, u, s):  # q even: r = s^2, t = s*sqrt(u), a_{2j+1} = t r^j, a_{2j} = (u/t) r^{1-j}
    n = F.n; m = (n-1)//2
    r = F.mul(s, s); t = F.mul(s, F.pow(u, F.q // 2)); ut = F.mul(u, F.inv(t))
    cyc = [None]*n
    for j in range(m+1): cyc[2*j] = F.mul(t, F.pow(r, j))
    for j in range(1, m+1): cyc[2*j-1] = F.mul(ut, F.pow(r, (1-j) % n))
    return cyc

def ham_exists(V, adj):
    n = len(V); nb = [[j for j in range(n) if j != i and adj(V[i], V[j])] for i in range(n)]
    path = [0]; used = [False]*n; used[0] = True
    def rec():
        if len(path) == n: return 0 in nb[path[-1]]
        for j in nb[path[-1]]:
            if not used[j]:
                used[j] = True; path.append(j)
                if rec(): return True
                path.pop(); used[j] = False
        return False
    return rec()

def check_certificates(sample_B=400):
    rng = random.Random(20261022)
    cA = json.load(open(os.path.join(HERE, 'certA.json')))
    assert sorted(int(q) for q in cA) == [q for q in EA if q > 7]
    nshe = 0
    for qs, rec in cA.items():
        q = int(qs)
        if rec['type'] == 'she':
            assert prime_power(q)[1] == 1
            F = Field(q); x = rec['x']
            assert F.is_prim(x) and F.is_prim((x+1) % q) and F.is_prim((x+2) % q)
            check_cycle_A(F, she_cycle(q, x)); nshe += 1
        else:
            F = Field(q, rec['f']); check_cycle_A(F, rec['cycle'])
    say("[3a] (A): all %d fields of EA with q > 7 certified (%d by She's construction, %d explicit cycles)" % (len(cA), nshe, len(cA)-nshe))
    cB = json.load(gzip.open(os.path.join(HERE, 'certB.json.gz'), 'rt'))
    assert sorted(int(q) for q in cB) == [q for q in EB if q > 7]
    counts = {}; fullchecked = 0; nones = []
    for qs, rec in cB.items():
        q = int(qs); F = Field(q, rec['f'])
        a0s = [c[0] for c in rec['certs']]
        assert sorted(a0s) == list(range(1, q)), q
        for c in rec['certs']:
            a0, typ = c[0], c[1]
            counts[typ] = counts.get(typ, 0) + 1
            if typ == 'alt':
                u, r = c[2], c[3]
                assert q % 2 == 1 and u != 0 and not F.is_square(u)
                assert F.is_square(r) and F.order_is(r, F.n // 2)
                assert F.is_prim(F.add(a0, u)) and F.is_prim(F.add(a0, F.mul(u, r)))
                if rng.random() < sample_B / 59000 or q < 30:
                    check_cycle_B(F, a0, alt_cycle(F, u, r)); fullchecked += 1
            elif typ == 'geo':
                u, s = c[2], c[3]
                assert q % 2 == 0 and u != 0 and F.is_prim(s)
                assert all(F.is_prim(F.add(a0, F.mul(u, F.pow(s, e)))) for e in (0, 1, 2))
                check_cycle_B(F, a0, geo_cycle(F, u, s)); fullchecked += 1
            elif typ == 'cycle':
                check_cycle_B(F, a0, c[2]); fullchecked += 1
            else:
                nones.append((q, a0))
    assert nones == [(11, 2)], nones
    F = Field(11)
    assert not ham_exists(list(range(1, 11)), lambda x, y: F.is_prim((2 + x*y) % 11))
    say("[3b] (B): all %d instances (q in EB, q > 7, a0 != 0) certified: %s; %d full cycles rebuilt and checked;"
        " (11,2) has no Hamilton cycle (exhaustive)" % (sum(counts.values()) - 1, counts, fullchecked))

# ---------------------------------------------------------------- 4. a0 = 0 (Luo-She)
def luo_she_seq(n):
    if n % 2 == 0:
        seq = []
        for i in range(n//2): seq += [2*i, n-1-2*i]
        return seq
    h = (n-1)//2
    if h == 3: return [6, 0, 5, 1, 3, 2, 4]
    seq = [n-1, 0, n-2, 1, h]
    for i in range(1, h-1): seq += [h-i, h+i]
    return seq
def check_a0_zero(LIM=3000):
    cnt = 0
    for p in primes_upto(LIM):
        x = p
        while x < LIM:
            if x > 4:
                n = x - 1; s = luo_she_seq(n)
                assert sorted(s) == list(range(n)), n
                assert all(gcd(s[i] + s[(i+1) % n], n) == 1 for i in range(n)), n
                cnt += 1
            x *= p
    for q in [8, 9, 16, 25, 27, 49, 121]:
        F = Field(q, convention='C')
        g = next(z for z in range(1, q) if F.is_prim(z))
        cyc = [F.pow(g, e) for e in luo_she_seq(q-1)]
        check_cycle_B(F, 0, cyc)
    say("[4] a0 = 0: Luo-She exponent sequences valid for all %d prime powers 4 < q < %d (proof for all q > 4: NOTE sec. 6;"
        " realised in F_q for 7 fields)" % (cnt, LIM))

# ---------------------------------------------------------------- 5. cross-check outputs
def fnv64(h, v):
    for i in range(8):
        h ^= (v >> (8*i)) & 0xff; h = (h * 1099511628211) & 0xFFFFFFFFFFFFFFFF
    return h
def crossB_hash(q):
    """Python re-implementation of the canonical certificate search of crossB.c (all a0 in F_q)."""
    F = Field(q, convention='C')
    n = q - 1
    if F.k == 1:
        g = next(z for z in range(2, q) if F.is_prim(z)) if q > 2 else 1
        log = [None]*q; e = 1
        for i in range(n): log[e] = i; e = e*g % q
        exp = None
    else:
        log = F.log
    prims = [z for z in range(q) if z and gcd(log[z], n) == 1]
    isprim = bytearray(q)
    for z in prims: isprim[z] = 1
    h = 1469598103934665603
    for a0 in range(q):
        cu = cr = -1
        if q % 2:
            m = n // 2; na0 = F.neg(a0)
            for u in range(1, q):
                if log[u] % 2 == 0 or not isprim[F.add(a0, u)]: continue
                for w in prims:
                    v = F.add(w, na0)
                    if v == 0 or log[v] % 2 == 0: continue
                    d = (log[v] - log[u]) % n
                    if gcd(d // 2, m) == 1:
                        cu = u; cr = F.mul(v, F.inv(u)); break
                if cu >= 0: break
        else:
            for u in range(1, q):
                if not isprim[F.add(a0, u)]: continue
                for s in prims:
                    us = F.mul(u, s); us2 = F.mul(us, s)
                    if isprim[F.add(a0, us)] and isprim[F.add(a0, us2)]: cu = u; cr = s; break
                if cu >= 0: break
        h = fnv64(h, (a0 << 40) if cu < 0 else ((a0 << 40) | (cu << 20) | cr))
    return '%016x' % h
def check_cross(full=False):
    lines = {}
    for ln in open(os.path.join(HERE, 'crossB.out')):
        d = dict(kv.split('=') for kv in ln.split() if '=' in kv)
        lines[int(d['q'])] = d
    tot = len(lines)
    bad = {q: d.get('fail_a0') for q, d in lines.items() if d['failsB'] != '0'}
    assert set(bad) == {9, 11, 13, 16, 25, 31, 43, 64, 256, 1024, 4096, 16384}, bad
    assert all(d['fail_a0'] == '0,' for q, d in lines.items() if q in (16, 64, 256, 1024, 4096, 16384))
    she_missing = [q for q, d in lines.items() if d['k'] == '1' and q % 2 and d['sheX'] == '-1']
    assert she_missing == [13, 29, 61], she_missing
    for q in [8, 27, 97, 125, 256, 1009]:
        assert crossB_hash(q) == lines[q]['hash'], q
    linesA = [ln for ln in open(os.path.join(HERE, 'crossA.out'))]
    assert len(linesA) == 74 and all('cycle_ok=1' in ln for ln in linesA)
    say("[5] crossB.out: %d prime powers 8 <= q < %d, all a0 (46,731,497 instances): construction fails only at"
        " a0 != 0 for q in {9,11,13,25,31,43} (all in EB) and a0 = 0 for q = 4^j; hashes re-derived for 6 sample q;"
        " crossA.out: 74 non-prime q, all cycles ok" % (tot, QCUT))
    if full:
        for src in ['crossA', 'crossB']:
            subprocess.run(['cc', '-O2', '-o', '/tmp/' + src, os.path.join(HERE, src + '.c')], check=True)
            out = subprocess.run(['/tmp/' + src, str(QCUT)], capture_output=True, text=True, timeout=600).stdout
            assert out == open(os.path.join(HERE, src + '.out')).read(), src
        say("[5'] --cross: crossA/crossB recompiled and rerun; outputs identical")

if __name__ == '__main__':
    check_lemma6()
    check_prop4()
    check_certificates()
    check_a0_zero()
    check_cross(full='--cross' in sys.argv)
    say("ALL CHECKS PASSED in %.1f s" % (time.time() - T_START))

