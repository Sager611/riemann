"""liuaudit2 T12: Remark 10.9 item 2 / Prop 5.2 / Lemma 4.2 of Liu 2610.26835v2 over F = Q(zeta12):
the theta-type coefficient a_0 = conj(Lambda) gamma_2 (angular type +-1 at both complex places) has NO
Patterson-type main term, while the angular-type-0 cubic Gauss sums gamma_2 do (Kubota-Patterson: X^{5/6}).
Sums run over squarefree exterior ideals n of norm <= X, in a fixed class b in B_S^(3) (81 classes) of the
3-allowed generator a_n (B_S^(3) an isotropic complement of the unit image in F_S^x/F_S^x3, Lemma 3.4, 3-part).
gamma_2(n) = q_n^{-1/2} sum_{x mod n} chi_n(x)^2 e(x/n) on a_n; CRT with cubic reciprocity (exact for 3-allowed).
Prediction (Liu, l.5290-5296, 2892): rms_b |sum a_0| ~ X^{1/2} (no main term); rms_b |sum gamma_2| ~ X^{5/6}.
Run: uv run --no-project --with numpy python3 theta_bias_q12.py [X=300000] [out.npz]
"""
import sys, math, cmath, itertools, random, time
import numpy as np
from fractions import Fraction as Fr

XMAX = int(float(sys.argv[1])) if len(sys.argv) > 1 else 300000
t0 = time.time()
# ---------------- field Q(zeta12), x^4 = x^2 - 1 ----------------
def mul(a, b, mod=None):
    c = [0]*7
    for i in range(4):
        ai = a[i]
        if ai == 0: continue
        for k in range(4): c[i+k] += ai*b[k]
    for k in (6, 5, 4):
        ck = c[k]
        if ck: c[k-2] += ck; c[k-4] -= ck; c[k] = 0
    r = tuple(c[:4])
    if mod is not None: r = tuple(t % mod for t in r)
    return r
def add(a, b): return tuple(a[i]+b[i] for i in range(4))
def scal(s, a): return tuple(s*t for t in a)
ONE = (1, 0, 0, 0); X_ = (0, 1, 0, 0)
SQM3 = (-1, 0, 2, 0); PI2 = (1, 0, 0, 1); EPSF = (1, 1, 0, 0)
TX = {5: (0, -1, 0, 1), 7: (0, -1, 0, 0), 11: (0, 1, 0, -1)}
def gal(a, k):
    if k == 1: return a
    t = TX[k]; r = (a[0], 0, 0, 0); p = ONE
    for i in range(1, 4):
        p = mul(p, t); r = add(r, scal(a[i], p))
    return r
def conjprod(a): return mul(gal(a, 5), mul(gal(a, 7), gal(a, 11)))
def norm(a):
    c = mul(a, conjprod(a)); assert c[1] == c[2] == c[3] == 0; return c[0]
def powe(a, n, mod=None):
    r = ONE; b = a
    while n:
        if n & 1: r = mul(r, b, mod)
        b = mul(b, b, mod); n >>= 1
    return r

# ---------------- local cube classes at S (mod 8 = p2^6, mod 9 = p3^4) ----------------
def coset_table(units, H, m, ell):
    reps = [((1, 0, 0, 0), ())]; cur = set(H)
    for a in units:
        if a in cur: continue
        newreps = []
        for (r, v) in reps:
            for k in range(ell): newreps.append((mul(r, powe(a, k, m), m), v + (k,)))
        reps = newreps
        cur = set(mul(r, h, m) for (r, v) in reps for h in H)
        if len(cur) == len(units): break
    L = max(len(v) for r, v in reps)
    table = {}
    for (r, v) in reps:
        v = v + (0,)*(L - len(v))
        for h in H: table[mul(r, h, m)] = v
    return table
units8 = [a for a in itertools.product(range(8), repeat=4) if norm(a) % 2]
units9 = [a for a in itertools.product(range(9), repeat=4) if norm(a) % 3]
CU8 = coset_table(units8, set(powe(a, 3, 8) for a in units8), 8, 3)
CU9 = coset_table(units9, set(powe(a, 3, 9) for a in units9), 9, 3)
def classes3(a):
    """cube class of an integral S-integer generator a in F_S^x/F_S^x3: vector in F_3^8."""
    n = norm(a); v2 = 0; v3 = 0
    while n % 4 == 0 and n > 0: n //= 4; v2 += 1
    while n % 9 == 0 and n > 0: n //= 9; v3 += 1
    def divpow(a, pi, v):
        for _ in range(v):
            c = mul(a, conjprod(pi)); N = norm(pi)
            assert all(t % N == 0 for t in c), (a, pi)
            a = tuple(t//N for t in c)
        return a
    u2 = divpow(a, PI2, v2); u3 = divpow(a, SQM3, v3)
    return (v2 % 3,) + CU8[tuple(t % 8 for t in u2)] + (v3 % 3,) + CU9[tuple(t % 9 for t in u3)]

# ---------------- primes ----------------
def primes_upto(n):
    s = bytearray([1])*(n+1); s[0] = s[1] = 0
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]
def prim_root(l):
    q = l - 1; fac = []
    p = 2
    while p*p <= q:
        if q % p == 0:
            fac.append(p)
            while q % p == 0: q //= p
        p += 1
    if q > 1: fac.append(q)
    for g in range(2, l):
        if all(pow(g, (l-1)//p, l) != 1 for p in fac): return g

GRAM = [[4, 0, 2, 0], [0, 4, 0, 2], [2, 0, 4, 0], [0, 2, 0, 4]]   # Tr(x^(i-j))
def lll(B, delta=0.75):
    """float LLL for the trace form (4-dim, entries < 1e6: double precision is ample)."""
    B = [list(map(int, b)) for b in B]; n = len(B)
    def ip(u, v): return float(sum(u[i]*GRAM[i][j]*v[j] for i in range(4) for j in range(4)))
    def gso(B):
        Bs = []; mu = [[0.0]*n for _ in range(n)]; c = []
        for i in range(n):
            v = [float(x) for x in B[i]]
            for j in range(i):
                mu[i][j] = sum(B[i][a]*GRAM[a][b]*Bs[j][b] for a in range(4) for b in range(4)) / c[j]
                v = [v[k] - mu[i][j]*Bs[j][k] for k in range(4)]
            Bs.append(v); c.append(sum(v[a]*GRAM[a][b]*v[b] for a in range(4) for b in range(4)))
        return Bs, mu, c
    Bs, mu, c = gso(B); k = 1; it = 0
    while k < n and it < 10000:
        it += 1
        for j in range(k-1, -1, -1):
            q = int(round(mu[k][j]))
            if q:
                B[k] = [B[k][t] - q*B[j][t] for t in range(4)]
                Bs, mu, c = gso(B)
        if c[k] >= (delta - mu[k][k-1]**2)*c[k-1]:
            k += 1
        else:
            B[k], B[k-1] = B[k-1], B[k]
            Bs, mu, c = gso(B); k = max(k-1, 1)
    return [tuple(b) for b in B]
def find_generator(basis, target):
    R = lll(basis)
    for rng in ((-1, 0, 1), (-2, -1, 0, 1, 2)):
        for coeffs in itertools.product(rng, repeat=4):
            if not any(coeffs): continue
            a = tuple(sum(coeffs[i]*R[i][t] for i in range(4)) for t in range(4))
            if abs(norm(a)) == target: return a
    return None

primes = []
for l in primes_upto(XMAX):
    if l in (2, 3): continue
    if l % 12 == 1:
        g = prim_root(l); z = pow(g, (l-1)//12, l)
        for k in (1, 5, 7, 11):
            r = pow(z, k, l)
            assert (r**4 - r**2 + 1) % l == 0
            basis = [(l, 0, 0, 0), ((-r) % l, 1, 0, 0), ((-r*r) % l, 0, 1, 0), ((-r**3) % l, 0, 0, 1)]
            gen = find_generator(basis, l)
            assert gen is not None, (l, r)
            primes.append({'l': l, 'deg': 1, 'r': r, 'q': l, 'gen': gen})
    elif l*l <= XMAX:
        quads = []
        for b in range(l):
            for c in range(l):
                p = [1, 0, (-1) % l, 0, 1]
                for i in range(3):
                    q_ = p[i]
                    if q_: p[i+1] = (p[i+1] - q_*b) % l; p[i+2] = (p[i+2] - q_*c) % l; p[i] = 0
                if p[3] == 0 and p[4] == 0: quads.append((b, c))
        assert len(quads) == 2
        for (b, c) in quads:
            basis = [(l, 0, 0, 0), (0, l, 0, 0), (c, b, 1, 0), (0, c, b, 1)]
            gen = find_generator(basis, l*l)
            assert gen is not None
            primes.append({'l': l, 'deg': 2, 'bc': (b, c), 'q': l*l, 'gen': gen})
primes.sort(key=lambda P: (P['q'], P['l'], P.get('r', 0)))
print("prime ideals coprime to 6 with norm <= %d: %d (%.1fs)" % (XMAX, len(primes), time.time() - t0))

# ---------------- residue maps ----------------
def red1(P, a):
    l, r = P['l'], P['r']
    return (a[0] + a[1]*r + a[2]*r*r + a[3]*r*r*r) % l
class RF2:
    def __init__(self, P):
        self.l = P['l']; self.b, self.c = P['bc']; self.q = self.l**2
        l = self.l; self.one = (1, 0)
        q1 = self.q - 1
        fac = [p for p in primes_upto(q1) if q1 % p == 0]
        for g in [(u, v) for u in range(l) for v in range(l) if (u, v) != (0, 0)]:
            if all(self.pw(g, q1//p) != self.one for p in fac): break
        self.g = g; self.log = {}; y = self.one
        for m in range(q1): self.log[y] = m; y = self.mul(y, g)
        z6 = self.red((0, 0, 1, 0)); h = self.pw(g, q1//6)
        y = self.one; k0 = None
        for k in range(6):
            if y == h: k0 = k
            y = self.mul(y, z6)
        self.k0 = k0
    def mul(self, s, t):
        l, b, c = self.l, self.b, self.c
        u0 = s[0]*t[0]; u1 = s[0]*t[1] + s[1]*t[0]; u2 = s[1]*t[1]
        return ((u0 - u2*c) % l, (u1 - u2*b) % l)
    def pw(self, s, n):
        r = self.one
        while n:
            if n & 1: r = self.mul(r, s)
            s = self.mul(s, s); n >>= 1
        return r
    def red(self, a):
        l = self.l
        x2 = ((-self.c) % l, (-self.b) % l); x3 = self.mul((0, 1), x2)
        return ((a[0] + a[2]*x2[0] + a[3]*x3[0]) % l, (a[1] + a[2]*x2[1] + a[3]*x3[1]) % l)
    def tr(self, s): return (2*s[0] - self.b*s[1]) % self.l
    def chi_exp(self, s): return (self.log[s]*self.k0) % 6

LOG1 = {}
def chi_exp(Pi, a):
    """k in Z/6 with chi_P(a) = zeta6^k (a integral, P-unit)"""
    P = primes[Pi]
    if P['deg'] == 1:
        if Pi not in LOG1:
            l, r = P['l'], P['r']; g = prim_root(l)
            lg = np.zeros(l, dtype=np.int64); y = 1
            for m in range(l-1): lg[y] = m; y = y*g % l
            h = pow(g, (l-1)//6, l); k0 = next(k for k in range(6) if pow(r*r, k, l) == h)
            LOG1[Pi] = (lg, k0)
        lg, k0 = LOG1[Pi]; v = red1(P, a); assert v != 0
        return int(lg[v]*k0 % 6)
    if 'rf' not in P: P['rf'] = RF2(P)
    R = P['rf']; s = R.red(a); assert s != (0, 0)
    return R.chi_exp(s)

# ---------------- 3-part of Lemma 3.4 ----------------
UNITS = [X_, EPSF, PI2, SQM3]
ucl = [classes3(u) for u in UNITS]
def Hexp(ka, a, kb, b):
    if ka == 'u' and kb == 'u': return 0
    if ka == 'u' and kb == 'p': return (-chi_exp(b, a)) % 6
    if ka == 'p' and kb == 'u': return chi_exp(a, b) % 6
    if a == b: return 0
    return (chi_exp(a, primes[b]['gen']) - chi_exp(b, primes[a]['gen'])) % 6
def rank_mod(vecs, ell):
    M = [list(v) for v in vecs]; r = 0; nc = len(M[0]) if M else 0
    for c in range(nc):
        piv = next((i for i in range(r, len(M)) if M[i][c] % ell), None)
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]; iv = pow(M[r][c], -1, ell); M[r] = [t*iv % ell for t in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % ell:
                f = M[i][c]; M[i] = [(M[i][k] - f*M[r][k]) % ell for k in range(nc)]
        r += 1
    return r
def solve_mod(cols, rhs, ell):
    n = len(rhs); m = len(cols)
    A = [[cols[k][i] % ell for k in range(m)] + [rhs[i] % ell] for i in range(n)]
    r = 0; pivc = []
    for c in range(m):
        piv = next((i for i in range(r, n) if A[i][c] % ell), None)
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]; iv = pow(A[r][c], -1, ell); A[r] = [t*iv % ell for t in A[r]]
        for i in range(n):
            if i != r and A[i][c] % ell:
                f = A[i][c]; A[i] = [(A[i][k] - f*A[r][k]) % ell for k in range(m+1)]
        pivc.append(c); r += 1
    assert r == m
    y = [0]*m
    for i, c in enumerate(pivc): y[c] = A[i][m]
    return y
pcl = {}
def pclass(i):
    if i not in pcl: pcl[i] = classes3(primes[i]['gen'])
    return pcl[i]
vecs = list(ucl); chosen = []
assert rank_mod(vecs, 3) == 4
for i in range(len(primes)):
    if primes[i]['q'] > 3000: break
    v = pclass(i)
    if rank_mod(vecs + [v], 3) > len(vecs): vecs.append(v); chosen.append(i)
    if len(vecs) == 8: break
basis = [('u', u) for u in UNITS] + [('p', i) for i in chosen]
G3 = [[Hexp(basis[i][0], basis[i][1], basis[k][0], basis[k][1]) % 3 for k in range(8)] for i in range(8)]
C3 = vecs
assert rank_mod(G3, 3) == 8
A = [[G3[i][4+k] for k in range(4)] for i in range(4)]
M = [A[i][:] + [1 if i == k else 0 for k in range(4)] for i in range(4)]
for c in range(4):
    piv = next(i for i in range(c, 4) if M[i][c] % 3); M[c], M[piv] = M[piv], M[c]
    iv = pow(M[c][c], -1, 3); M[c] = [t*iv % 3 for t in M[c]]
    for i in range(4):
        if i != c and M[i][c] % 3:
            f = M[i][c]; M[i] = [(M[i][k] - f*M[c][k]) % 3 for k in range(8)]
Ainv = [row[4:] for row in M]
def Gv(u, v): return sum(u[i]*G3[i][k]*v[k] for i in range(8) for k in range(8)) % 3
fp = []
for jj in range(4):
    v = [0]*8
    for k in range(4): v[4+k] = Ainv[k][jj] % 3
    fp.append(v)
Fb = []
for jj in range(4):
    v = fp[jj][:]
    for i in range(jj): v[i] = (v[i] + Gv(fp[i], fp[jj])) % 3
    Fb.append(v)
assert all(Gv(Fb[i], Fb[k]) == 0 for i in range(4) for k in range(4))
NEWC = [[1 if t == i else 0 for t in range(8)] for i in range(4)] + Fb
def decomp(vec):
    z = solve_mod(NEWC, solve_mod(C3, vec, 3), 3)
    return z[:4], tuple(z[4:])
for i in range(len(primes)):
    P = primes[i]; a = P['gen']
    e, b = decomp(pclass(i))
    for k in range(4):
        for _ in range((-e[k]) % 3): a = mul(a, UNITS[k])
    e2, b2 = decomp(classes3(a))
    assert all(t == 0 for t in e2) and b2 == b
    P['agen'] = a; P['b3'] = b
print("3-allowed generators built (%.1fs)" % (time.time() - t0))

# ---------------- gamma_2(p), Lambda(p) ----------------
RESCACHE = {}
def cubic_classes_deg1(l):
    """x^((l-1)/3) mod l for x = 1..l-1 (cached for the 4 primes above l)"""
    if l in RESCACHE: return RESCACHE[l]
    RESCACHE.clear()
    X = np.arange(1, l, dtype=np.int64)
    e = (l - 1)//3; res = np.ones_like(X); base = X.copy()
    while e:
        if e & 1: res = res*base % l
        base = base*base % l; e >>= 1
    RESCACHE[l] = (X, res)
    return RESCACHE[l]
def gamma2(i):
    P = primes[i]; a = P['agen']; l = P['l']
    Na = norm(a); cp = conjprod(a)
    if P['deg'] == 1:
        m = Na // l; assert Na % l == 0 and m % l
        w = red1(P, cp) * pow(m % l, -1, l) % l
        assert w != 0
        X, res = cubic_classes_deg1(l)
        r = P['r']; r4 = pow(r, 4, l); r8 = r4*r4 % l
        mexp = np.where(res == 1, 0, np.where(res == r4, 1, np.where(res == r8, 2, -1)))
        assert (mexp >= 0).all()
        ph = np.exp(2j*np.pi*(mexp/3.0 - ((X*w) % l)/l))
        return complex(ph.sum())/math.sqrt(l)
    # degree 2: F_q = F_l[t]/(t^2 + b t + c), vectorized pair arithmetic (no log tables)
    b, c = P['bc']; m = Na // (l*l); assert Na % (l*l) == 0 and m % l
    assert all(x % l == 0 for x in cp)          # (cp) = l * P' : ell/a = (cp/l)/m
    cpl = tuple(x//l for x in cp)
    def redq(a_):
        x2 = ((-c) % l, (-b) % l)
        x3 = ((-c*x2[1]) % l, (x2[0] - b*x2[1]) % l)   # t * x2 = x2[0] t + x2[1] t^2
        return ((a_[0] + a_[2]*x2[0] + a_[3]*x3[0]) % l, (a_[1] + a_[2]*x2[1] + a_[3]*x3[1]) % l)
    wv = redq(cpl); im = pow(m % l, -1, l); wv = ((wv[0]*im) % l, (wv[1]*im) % l)
    assert wv != (0, 0)
    U, V = np.meshgrid(np.arange(l, dtype=np.int64), np.arange(l, dtype=np.int64), indexing='ij')
    U = U.ravel()[1:]; V = V.ravel()[1:]                 # all nonzero (u, v)
    def pmul(u1, v1, u2, v2):
        s0 = u1*u2; s1 = u1*v2 + v1*u2; s2 = v1*v2
        return (s0 - s2*c) % l, (s1 - s2*b) % l
    e = (l*l - 1)//3; ru = np.ones_like(U); rv = np.zeros_like(V); bu, bv = U.copy(), V.copy()
    while e:
        if e & 1: ru, rv = pmul(ru, rv, bu, bv)
        bu, bv = pmul(bu, bv, bu, bv); e >>= 1
    om = redq((-1, 0, 1, 0))                             # omega = x^4 = x^2 - 1
    om2 = pmul(np.array([om[0]]), np.array([om[1]]), np.array([om[0]]), np.array([om[1]]))
    o2 = (int(om2[0][0]), int(om2[1][0]))
    mexp = np.where((ru == 1) & (rv == 0), 0, np.where((ru == om[0]) & (rv == om[1]), 1, np.where((ru == o2[0]) & (rv == o2[1]), 2, -1)))
    assert (mexp >= 0).all()
    pu, pv = pmul(U, V, np.full_like(U, wv[0]), np.full_like(V, wv[1]))
    tr = (2*pu - b*pv) % l
    ph = np.exp(2j*np.pi*(mexp/3.0 - tr/l))
    return complex(ph.sum())/l
def Lambda(i):
    a = primes[i]['gen']
    m = mul(a, gal(a, 7)); assert m[1] == 0 and m[3] == 0
    a_, b_ = m[0] + m[2], m[2]
    cands = []
    for k in range(3):
        for s in (1, -1): cands.append((s*a_, s*b_))
        a_, b_ = -b_, a_ - b_
    good = [c for c in cands if c[0] % 3 == 1 and c[1] % 3 == 0]
    assert len(good) == 1
    z = good[0][0] + good[0][1]*cmath.exp(2j*math.pi/3)
    return z/abs(z)
G2 = np.zeros(len(primes), dtype=complex); LAM = np.zeros(len(primes), dtype=complex)
for i in range(len(primes)):
    G2[i] = gamma2(i); LAM[i] = Lambda(i)
print("gamma_2, Lambda for all primes (%.1fs); max | |gamma2|-1 | = %.1e; max |gamma2^3 + Lambda| = %.1e" % (
    time.time() - t0, float(np.max(np.abs(np.abs(G2) - 1))), float(np.max(np.abs(G2**3 + LAM)))))

def cub(i, j):   # exponent e in Z/3 with chi_{P_i}(a_j)^2 = zeta3^e  (chi_P(a) = zeta6^k  =>  chi_P(a)^2 = zeta3^k)
    return chi_exp(i, primes[j]['agen']) % 3
# (fixed 12 Nov follow-up: the first version used 2k mod 3, the exponent of chi^4, which conjugated every CRT pair phase)
SMALL = int(XMAX**0.5) + 1
sm = [i for i in range(len(primes)) if primes[i]['q'] <= SMALL]
nrec = bad = 0
for i, j in itertools.combinations(sm[:60], 2):
    nrec += 1
    if cub(i, j) != cub(j, i): bad += 1
print("cubic reciprocity on 3-allowed generators: %d pairs, %d failures" % (nrec, bad))
assert bad == 0

CHK = sorted(set([c for c in (1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6) if c <= XMAX] + [float(XMAX)]))
nb = 81
def bidx(b): return b[0] + 3*b[1] + 9*b[2] + 27*b[3]
S = {t: np.zeros((len(CHK), nb), dtype=complex) for t in ('g2', 'Lg2', 'Lbg2')}
CNT = np.zeros((len(CHK), nb))
qs = [P['q'] for P in primes]
Z3 = [1, cmath.exp(2j*math.pi/3), cmath.exp(4j*math.pi/3)]
CUBC = {}
def cubc(i, j):
    key = (i, j)
    if key not in CUBC: CUBC[key] = cub(i, j)
    return CUBC[key]
count = 0
REC_N = []; REC_B = []; REC_G = []; REC_L = []
def record(N, g, lam, b):
    global count
    count += 1
    REC_N.append(N); REC_B.append(bidx(b)); REC_G.append(g); REC_L.append(lam)
    k = next(t for t, c in enumerate(CHK) if N <= c)
    bi = bidx(b)
    S['g2'][k, bi] += g; S['Lbg2'][k, bi] += lam.conjugate()*g; S['Lg2'][k, bi] += lam*g; CNT[k, bi] += 1
sys.setrecursionlimit(10000)
def dfs(start, N, g, lam, b, members):
    for j in range(start, len(primes)):
        q = qs[j]
        if N*q > XMAX: break
        f = 1+0j
        for i in members: f *= Z3[(2*cubc(i, j)) % 3]   # (chi_i(a_j) chi_j(a_i))^2 = zeta3^(cub(i,j)+cub(j,i)) = zeta3^(2 cub)
        g2 = g*G2[j]*f; lam2 = lam*LAM[j]
        bj = primes[j]['b3']; b2 = ((b[0]+bj[0]) % 3, (b[1]+bj[1]) % 3, (b[2]+bj[2]) % 3, (b[3]+bj[3]) % 3)
        record(N*q, g2, lam2, b2)
        if j + 1 < len(primes) and N*q*qs[j+1] <= XMAX:
            dfs(j+1, N*q, g2, lam2, b2, members + [j])
record(1, 1+0j, 1+0j, (0, 0, 0, 0))
dfs(0, 1, 1+0j, 1+0j, (0, 0, 0, 0), [])
print("squarefree exterior ideals with norm <= %d: %d (%.1fs)" % (XMAX, count, time.time() - t0))
cum = {t: np.cumsum(S[t], axis=0) for t in S}; ccum = np.cumsum(CNT, axis=0)
print("%9s %8s | %12s %12s %12s | %10s %10s %10s" % ("X", "n/class", "rms g2", "rms Lbar g2", "rms L g2", "b=0: g2", "Lbar g2", "L g2"))
for k, c in enumerate(CHK):
    rms = {t: float(np.sqrt(np.mean(np.abs(cum[t][k])**2))) for t in cum}
    print("%9d %8.1f | %12.2f %12.2f %12.2f | %10.2f %10.2f %10.2f" % (c, ccum[k].mean(), rms['g2'], rms['Lbg2'], rms['Lg2'],
          abs(cum['g2'][k][0]), abs(cum['Lbg2'][k][0]), abs(cum['Lg2'][k][0])))
if len(CHK) >= 3:
    k1, k2 = len(CHK) - 3, len(CHK) - 1
    for t in ('g2', 'Lbg2', 'Lg2'):
        r1 = np.sqrt(np.mean(np.abs(cum[t][k1])**2)); r2 = np.sqrt(np.mean(np.abs(cum[t][k2])**2))
        print("growth exponent of rms_b |sum %s| between X=%d and %d: %.3f" % (t, CHK[k1], CHK[k2], math.log(r2/r1)/math.log(CHK[k2]/CHK[k1])))
    r1 = np.sqrt(ccum[k1].mean()); r2 = np.sqrt(ccum[k2].mean())
    print("random-walk reference sqrt(n/class): exponent %.3f" % (math.log(r2/r1)/math.log(CHK[k2]/CHK[k1])))
OUT = sys.argv[2] if len(sys.argv) > 2 else None
if OUT:
    np.savez_compressed(OUT, N=np.array(REC_N, dtype=np.int64), b=np.array(REC_B, dtype=np.int16), g=np.array(REC_G), lam=np.array(REC_L),
        pq=np.array(qs, dtype=np.int64), pb=np.array([bidx(P['b3']) for P in primes], dtype=np.int16), pg=G2, plam=LAM)
    print('saved', OUT)
from collections import defaultdict as _dd
_acc = _dd(complex)
for _t in range(len(REC_N)):
    if REC_B[_t] == 0: _acc[REC_N[_t]] += np.conj(REC_L[_t])*REC_G[_t]
_im = [abs(z.imag) for z in _acc.values()]
print("Galois check: class-0 per-norm sums of conj(Lambda) gamma_2 real for %d of %d norms (max |Im| %.1e)" % (sum(1 for x in _im if x < 1e-8), len(_im), max(_im)))
print("elapsed %.1fs" % (time.time() - t0))
