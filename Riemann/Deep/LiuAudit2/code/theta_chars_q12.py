"""liuaudit2: Lemmas 3.4, 3.5, 4.2 of Liu arXiv 2610.26835v2 over F = Q(zeta12).

F = Q(x), x^4 = x^2 - 1 (x = zeta12), O_F = Z[x] (class number 1), r_F = 2, [F:Q] = 4.
S = {two complex places, p2 = (1+i), p3 = (sqrt(-3))}, R = Z[x][1/6], U_S = <x, 1+x, 1+i, sqrt(-3)>.
Embedding j: x -> exp(2 pi i/12).  chi_p(u) = j(lift of u^((q-1)/6) mod p) (l.1880).

Tests (TeX line numbers of number_fields_v4.tex):
 T1  gamma_2(p)^3 = -Lambda(p), Lambda = Lambda_{E0} o N_{F/E0}             (Lemma 4.2, l.2378-2392)
 T2  gamma_1 gamma_2 = conj(chi_p(4)) gamma_3 gamma_2^3                     (l.2482-2486)
 T3  H_S nondegenerate on V_S = F_S^x/F_S^x6; unit image E_S injective, isotropic (Lemma 3.4, l.1908)
 T4  allowed generators (class in an isotropic complement B_S): chi_a(b) = chi_b(a)  (eq. l.1926)
 T5  chi_a(eps) = H_S(eps,a)^-1 for all generators a (eq. l.1928); xi_eps factors through B_S
 T6  gamma_3(a), a allowed, = lambda([a]) with lambda a quadratic character of B_S;
     G_0(p) = conj(chi_p(4)) gamma_3(p) a character of B_S of order | 6     (l.2440-2481)
 T7  a_0(n) gamma_1(n) = mu(n) G_0(n), mu(n) gamma_{-1}(n) = a_0(n) conj G_0(n),
     n squarefree allowed composite (Gauss sums by CRT, verified directly)  (Lemma 4.2, l.2356; §10 l.3675-3700)
 T8  gamma_j(eps h)-type unit substitution: G(U, eps h) = conj(xi_eps(U)) G(U,h) (Lemma 3.5, l.1990)
Run: uv run --no-project python3 theta_chars_q12.py [normbound]
"""
import sys, cmath, math, itertools, random, time
from fractions import Fraction as Fr

NB = int(sys.argv[1]) if len(sys.argv) > 1 else 3000
t0 = time.time()
random.seed(1)
FAIL = []
def check(name, ok, info=""):
    if not ok:
        FAIL.append((name, info))
    return ok

# ---------- field arithmetic ----------
def mul(a, b, mod=None):
    c = [0]*7
    for i in range(4):
        ai = a[i]
        if ai == 0: continue
        for k in range(4):
            c[i+k] += ai*b[k]
    for k in (6, 5, 4):
        ck = c[k]
        if ck:
            c[k-2] += ck; c[k-4] -= ck; c[k] = 0
    r = c[:4]
    if mod is not None:
        r = [t % mod for t in r]
    return tuple(r)

def add(a, b): return tuple(a[i]+b[i] for i in range(4))
def neg(a): return tuple(-t for t in a)
def scal(s, a): return tuple(s*t for t in a)
ONE = (1, 0, 0, 0); X = (0, 1, 0, 0)
I_ = (0, 0, 0, 1)                 # i = x^3
OMEGA = (-1, 0, 1, 0)             # omega = x^4 = x^2 - 1
SQM3 = (-1, 0, 2, 0)              # sqrt(-3) = 2 omega + 1 = 2x^2 - 1
PI2 = (1, 0, 0, 1)                # 1 + i
EPSF = (1, 1, 0, 0)               # 1 + x  (norm 1)

def powe(a, n, mod=None):
    r = ONE; b = a
    if n < 0: raise ValueError
    while n:
        if n & 1: r = mul(r, b, mod)
        b = mul(b, b, mod); n >>= 1
    return r

TX = {5: (0, -1, 0, 1), 7: (0, -1, 0, 0), 11: (0, 1, 0, -1)}   # images of x under x -> x^k
def gal(a, k):
    if k == 1: return a
    t = TX[k]; r = (a[0], 0, 0, 0); p = ONE
    for i in range(1, 4):
        p = mul(p, t)
        r = add(r, scal(a[i], p))
    return r
def norm(a):
    c = mul(mul(a, gal(a, 5)), mul(gal(a, 7), gal(a, 11)))
    assert c[1] == 0 and c[2] == 0 and c[3] == 0, c
    return c[0]
def inv(a):
    c = mul(gal(a, 5), mul(gal(a, 7), gal(a, 11)))
    n = norm(a)
    return tuple(Fr(t, n) for t in c)
def mulF(a, b):  # Fractions
    c = [Fr(0)]*7
    for i in range(4):
        for k in range(4):
            c[i+k] += Fr(a[i])*Fr(b[k])
    for k in (6, 5, 4):
        ck = c[k]
        if ck:
            c[k-2] += ck; c[k-4] -= ck; c[k] = Fr(0)
    return tuple(c[:4])
def normF(a):
    # norm of Fraction element: clear denominators
    d = 1
    for t in a: d = d*Fr(t).denominator//math.gcd(d, Fr(t).denominator)
    ai = tuple(int(Fr(t)*d) for t in a)
    return Fr(norm(ai), d**4)
def powF(a, n):
    if n < 0:
        a = invF(a); n = -n
    r = tuple(Fr(t) for t in ONE)
    for _ in range(n): r = mulF(r, a)
    return r
def invF(a):
    d = 1
    for t in a: d = d*Fr(t).denominator//math.gcd(d, Fr(t).denominator)
    ai = tuple(int(Fr(t)*d) for t in a)
    ia = inv(ai)
    return tuple(t*d for t in ia)

def emb(a, k=1):
    z = cmath.exp(2j*math.pi*k/12)
    return sum(complex(a[i])*z**i for i in range(4))

def vp(n, p):
    n = Fr(n)
    if n == 0: return 10**9
    v = 0; a, b = n.numerator, n.denominator
    while a % p == 0: a //= p; v += 1
    while b % p == 0: b //= p; v -= 1
    return v

def redmod(a, m):
    """reduce Fraction element (denominators prime to m) mod m"""
    out = []
    for t in a:
        t = Fr(t)
        out.append(t.numerator*pow(t.denominator, -1, m) % m)
    return tuple(out)

# ---------- local classes at S ----------
def build_tables():
    # mod 8 (p2^6) and mod 9 (p3^4)
    units8 = [a for a in itertools.product(range(8), repeat=4) if norm(a) % 2]
    units9 = [a for a in itertools.product(range(9), repeat=4) if norm(a) % 3]
    def subgroup(units, e, m):
        return set(powe(a, e, m) for a in units)
    def coset_table(units, H, m, ell):
        H = set(H); basis = []; span = {ONE[:]: ()}
        # span: element-class representatives -> coordinate vector
        span = {tuple(t % m for t in ONE): tuple()}
        reps = [((1, 0, 0, 0), ())]
        cur = set(H)
        for a in units:
            if a in cur: continue
            basis.append(a)
            newreps = []
            for (r, v) in reps:
                for k in range(ell):
                    newreps.append((mul(r, powe(a, k, m), m), v + (k,)))
            reps = newreps
            cur = set(mul(r, h, m) for (r, v) in reps for h in H)
            if len(cur) == len(units): break
        # pad coordinate vectors to same length
        L = len(basis)
        table = {}
        for (r, v) in reps:
            v = v + (0,)*(L - len(v))
            for h in H:
                table[mul(r, h, m)] = v
        assert len(table) == len(units), (len(table), len(units))
        return table, L
    SQ8 = subgroup(units8, 2, 8); CU8 = subgroup(units8, 3, 8)
    SQ9 = subgroup(units9, 2, 9); CU9 = subgroup(units9, 3, 9)
    T = {}
    T['sq8'], L1 = coset_table(units8, SQ8, 8, 2)
    T['cu8'], L2 = coset_table(units8, CU8, 8, 3)
    T['sq9'], L3 = coset_table(units9, SQ9, 9, 2)
    T['cu9'], L4 = coset_table(units9, CU9, 9, 3)
    return T, (len(units8), len(SQ8), len(CU8), L1, L2, len(units9), len(SQ9), len(CU9), L3, L4)

T, dims = build_tables()
print("local tables: |(O/8)^x|=%d |sq|=%d |cu|=%d dim2=%d dim3=%d; |(O/9)^x|=%d |sq|=%d |cu|=%d dim2=%d dim3=%d" % dims)
check("T3-local-dims", dims[3] == 5 and dims[4] == 1 and dims[8] == 1 and dims[9] == 5, str(dims))
INVPI2 = invF(PI2); INVPI3 = invF(SQM3)

def classes(a):
    """a: element (ints or Fractions), nonzero. returns (vec2 in F2^8, vec3 in F3^8)"""
    a = tuple(Fr(t) for t in a)
    N = normF(a)
    v2 = vp(N, 2); v3 = vp(N, 3)
    assert v2 % 2 == 0 and v3 % 2 == 0
    v2 //= 2; v3 //= 2
    u2 = mulF(a, powF(PI2, -v2)) if v2 else a
    u3 = mulF(a, powF(SQM3, -v3)) if v3 else a
    r8 = redmod(u2, 8); r9 = redmod(u3, 9)
    s8 = T['sq8'][r8]; c8 = T['cu8'][r8]; s9 = T['sq9'][r9]; c9 = T['cu9'][r9]
    vec2 = (v2 % 2,) + s8 + (v3 % 2,) + s9
    vec3 = (v2 % 3,) + c8 + (v3 % 3,) + c9
    return vec2, vec3

# ---------- primes ----------
def primes_upto(n):
    s = bytearray([1])*(n+1); s[0] = s[1] = 0
    for i in range(2, int(n**0.5)+1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n+1) if s[i]]

# residue field data for prime ideals: ('lin', ell, r) or ('quad', ell, (b, c)) meaning x^2 + b x + c
def factor_mod(ell):
    f = lambda t: (t**4 - t**2 + 1) % ell
    roots = [r for r in range(ell) if f(r) == 0]
    if roots:
        return [('lin', ell, r) for r in roots]
    quads = []
    for b in range(ell):
        for c in range(ell):
            # divide x^4 - x^2 + 1 by x^2 + b x + c mod ell; remainder zero?
            p = [1, 0, -1, 0, 1]  # coefficients high->low
            p = [t % ell for t in p]
            for i in range(3):
                q = p[i]
                if q:
                    p[i+1] = (p[i+1] - q*b) % ell
                    p[i+2] = (p[i+2] - q*c) % ell
                    p[i] = 0
            if p[3] == 0 and p[4] == 0:
                quads.append(('quad', ell, (b, c)))
    assert len(quads) == 2, (ell, quads)
    return quads

class RF:
    """residue field O/P; elements as ints (lin) or pairs (quad)."""
    def __init__(self, P):
        self.P = P; kind, ell, d = P; self.ell = ell
        self.q = ell if kind == 'lin' else ell*ell
        if kind == 'lin':
            self.r = d
        else:
            self.b, self.c = d
        # element enumeration and log tables
        self.zero = self.red((0, 0, 0, 0))
        self.one = self.red(ONE)
        els = self.elements()
        # find generator
        q1 = self.q - 1
        fac = [p for p in primes_upto(q1) if q1 % p == 0]
        for g in els:
            if g == self.zero: continue
            if all(self.pw(g, q1//p) != self.one for p in fac):
                break
        self.g = g
        self.log = {}; y = self.one
        for m in range(q1):
            self.log[y] = m; y = self.mul(y, g)
        t = self.red(X); z6 = self.mul(t, t)  # zeta6 = x^2
        h = self.pw(g, q1//6)
        k0 = None; y = self.one
        for k in range(6):
            if y == h: k0 = k
            y = self.mul(y, z6)
        assert k0 is not None
        self.k0 = k0   # chi(g) = zeta6^k0 (exponent)
    def elements(self):
        if self.P[0] == 'lin': return list(range(self.ell))
        return [(u, v) for u in range(self.ell) for v in range(self.ell)]
    def red(self, a):
        """reduce Fraction/int element of O_(P)"""
        ell = self.ell
        a = [Fr(t) for t in a]
        a = [t.numerator*pow(t.denominator, -1, ell) % ell for t in a]
        if self.P[0] == 'lin':
            r = self.r
            return (a[0] + a[1]*r + a[2]*r*r + a[3]*r*r*r) % ell
        b, c = self.b, self.c
        # powers of x mod x^2 + b x + c: x^2 = -b x - c
        x2 = ((-c) % ell, (-b) % ell)
        x3 = self.mul((0, 1), x2)
        u = (a[0] + a[2]*x2[0] + a[3]*x3[0]) % ell
        v = (a[1] + a[2]*x2[1] + a[3]*x3[1]) % ell
        return (u, v)
    def mul(self, s, t):
        ell = self.ell
        if self.P[0] == 'lin': return s*t % ell
        b, c = self.b, self.c
        u0 = s[0]*t[0]; u1 = s[0]*t[1] + s[1]*t[0]; u2 = s[1]*t[1]
        return ((u0 - u2*c) % ell, (u1 - u2*b) % ell)
    def pw(self, s, n):
        r = self.one
        while n:
            if n & 1: r = self.mul(r, s)
            s = self.mul(s, s); n >>= 1
        return r
    def tr(self, s):
        if self.P[0] == 'lin': return s % self.ell
        return (2*s[0] - self.b*s[1]) % self.ell
    def chi_exp(self, s):
        """exponent e in Z/6 with chi_P(s) = zeta6^e (None if s = 0)"""
        if s == self.zero: return None
        return (self.log[s]*self.k0) % 6

def resexp(RFp, a):
    """exponent of (a/P)_6 for element a (Fractions ok, P-integral)."""
    return RFp.chi_exp(RFp.red(a))

Z6 = [cmath.exp(2j*math.pi*k/6) for k in range(6)]

def gauss(RFp, gen, j):
    """gamma_j = q^(-1/2) sum_x chi_P(x)^j e(-Tr(x w)/ell), w = ell/gen mod P."""
    ell = RFp.ell
    w = RFp.red(tuple(ell*t for t in invF(gen)))
    assert w != RFp.zero
    s = 0j
    g = RFp.g; y = RFp.one; k0 = RFp.k0
    for m in range(RFp.q - 1):
        e = (m*k0*j) % 6
        s += Z6[e]*cmath.exp(-2j*math.pi*RFp.tr(RFp.mul(y, w))/ell)
        y = RFp.mul(y, g)
    return s/math.sqrt(RFp.q)

# ---------- find prime elements ----------
B = 6
plist = {}  # P -> generator (ints)
ells = [l for l in primes_upto(int(NB**1)) if l > 3 and (l % 12 == 1 or l*l <= NB)]
FACT = {l: factor_mod(l) for l in ells}
Ps = {}
for l in ells:
    for P in FACT[l]:
        q = l if P[0] == 'lin' else l*l
        if q <= NB: Ps[P] = None
want = len(Ps)
for a in itertools.product(range(-B, B+1), repeat=4):
    if not any(a): continue
    n = norm(a)
    if n > NB or n < 5: continue
    # prime norm or square of prime with l != 1 mod 12
    l = None
    if n in set(): pass
    rt = int(round(n**0.5))
    if rt*rt == n and rt % 12 != 1 and all(rt % p for p in range(2, int(rt**0.5)+1)) and rt > 3:
        l = rt
    elif n % 12 == 1 and all(n % p for p in range(2, int(n**0.5)+1)):
        l = n
    if l is None or l not in FACT: continue
    for P in FACT[l]:
        if P in Ps and Ps[P] is None:
            R = RF.__new__(RF); R.P = P; R.ell = l
            if P[0] == 'lin': R.r = P[2]
            else: R.b, R.c = P[2]
            if RF.red(R, a) == (0 if P[0] == 'lin' else (0, 0)):
                Ps[P] = a
missing = [P for P in Ps if Ps[P] is None]
print("prime ideals with norm <= %d: %d, generators found for %d" % (NB, want, want - len(missing)))
Ps = {P: g for P, g in Ps.items() if g is not None}
plist = sorted(Ps.items(), key=lambda t: (t[0][1] if t[0][0] == 'lin' else t[0][1]**2, t[0]))
RFs = {P: RF(P) for P, _ in plist}
print("residue fields built %.1fs" % (time.time() - t0))

# ---------- T1, T2 ----------
def Lambda(gen):
    """Lambda_{E0}(N_{F/E0} gen) with pi_+ = 1 mod 3 (gen in O_F, generating P exactly)."""
    m = mul(gen, gal(gen, 7))
    assert m[1] == 0 and m[3] == 0
    A, Bc = m[0] + m[2], m[2]   # A + B omega
    # associates: multiply by units +-1, +-omega, +-omega^2 ; omega^2 = -1 - omega
    def mulw(A, Bc):  # (A + B w) * w = A w + B w^2 = -B + (A - B) w
        return (-Bc, A - Bc)
    cands = []
    a_, b_ = A, Bc
    for k in range(3):
        for s in (1, -1):
            cands.append((s*a_, s*b_))
        a_, b_ = mulw(a_, b_)
    good = [c for c in cands if c[0] % 3 == 1 and c[1] % 3 == 0]
    assert len(good) == 1, (A, Bc, cands)
    A1, B1 = good[0]
    z = A1 + B1*cmath.exp(2j*math.pi/3)
    return z/abs(z)

nT1 = nT2 = 0; maxerr1 = maxerr2 = 0.0
gam = {}
for P, gen in plist:
    R = RFs[P]
    g1 = gauss(R, gen, 1); g2 = gauss(R, gen, 2); g3 = gauss(R, gen, 3)
    gam[P] = (g1, g2, g3)
    L = Lambda(gen)
    e1 = abs(g2**3 + L); maxerr1 = max(maxerr1, e1)
    check("T1", e1 < 1e-8, "%s %s %g" % (P, gen, e1)); nT1 += 1
    c4 = Z6[(-resexp(R, (4, 0, 0, 0))) % 6]
    e2 = abs(g1*g2 - c4*g3*g2**3); maxerr2 = max(maxerr2, e2)
    check("T2", e2 < 1e-8, "%s %g" % (P, e2)); nT2 += 1
    check("T2-g3real", abs(abs(g3) - 1) < 1e-8 and abs(g3.imag) < 1e-8, str(g3))
print("T1 gamma2^3 = -Lambda: %d primes, max err %.2e; T2 duplication: %d primes, max err %.2e  (%.1fs)" % (nT1, maxerr1, nT2, maxerr2, time.time() - t0))

# ---------- T3: Hilbert pairing on V_S via tame reciprocity ----------
UNITS = [X, EPSF, PI2, SQM3]
ucl = [classes(u) for u in UNITS]
pcl = {P: classes(g) for P, g in plist}

def rank_mod(vecs, ell):
    M = [list(v) for v in vecs]; r = 0; ncol = len(M[0]) if M else 0
    for c in range(ncol):
        piv = None
        for i in range(r, len(M)):
            if M[i][c] % ell: piv = i; break
        if piv is None: continue
        M[r], M[piv] = M[piv], M[r]
        iv = pow(M[r][c], -1, ell)
        M[r] = [t*iv % ell for t in M[r]]
        for i in range(len(M)):
            if i != r and M[i][c] % ell:
                f = M[i][c]
                M[i] = [(M[i][k] - f*M[r][k]) % ell for k in range(ncol)]
        r += 1
    return r

def solve_mod(Mcols, rhs, ell):
    """solve sum_k y_k Mcols[k] = rhs (vectors length n) mod ell; Mcols basis (n vectors)."""
    n = len(rhs)
    A = [[Mcols[k][i] % ell for k in range(len(Mcols))] + [rhs[i] % ell] for i in range(n)]
    m = len(Mcols); r = 0; pivc = []
    for c in range(m):
        piv = None
        for i in range(r, n):
            if A[i][c] % ell: piv = i; break
        if piv is None: continue
        A[r], A[piv] = A[piv], A[r]
        iv = pow(A[r][c], -1, ell)
        A[r] = [t*iv % ell for t in A[r]]
        for i in range(n):
            if i != r and A[i][c] % ell:
                f = A[i][c]
                A[i] = [(A[i][k] - f*A[r][k]) % ell for k in range(m+1)]
        pivc.append(c); r += 1
    assert r == m
    y = [0]*m
    for i, c in enumerate(pivc): y[c] = A[i][m]
    return y

rk2 = rank_mod([u[0] for u in ucl], 2); rk3 = rank_mod([u[1] for u in ucl], 3)
print("T3 unit image ranks: V2 %d, V3 %d (Lemma 3.4: 4 = |S_f| + [F:Q]/2)" % (rk2, rk3))
check("T3-unit-injective", rk2 == 4 and rk3 == 4)

# choose primes extending unit images to bases of V2 and V3
def extend(ell, idx):
    chosen = []; vecs = [u[idx] for u in ucl]
    for P, g in plist:
        v = pcl[P][idx]
        if rank_mod(vecs + [v], ell) > len(vecs):
            vecs.append(v); chosen.append(P)
        if len(vecs) == 8: break
    return chosen
B2P = extend(2, 0); B3P = extend(3, 1)
print("basis primes V2:", [(P[0], P[1]) for P in B2P], " V3:", [(P[0], P[1]) for P in B3P])

def Hexp(a_kind, a, b_kind, b):
    """exponent (Z/6) of H_S(a,b); kinds: ('u', unit element) or ('p', P)."""
    if a_kind == 'u' and b_kind == 'u': return 0
    if a_kind == 'u' and b_kind == 'p':
        return (-resexp(RFs[b], a)) % 6
    if a_kind == 'p' and b_kind == 'u':
        return resexp(RFs[a], b) % 6
    if a == b: return 0
    # H_S(p,p') = (p'/P) (p/P')^{-1}
    return (resexp(RFs[a], Ps[b]) - resexp(RFs[b], Ps[a])) % 6

def gram(ell, chosen, idx):
    basis = [('u', u) for u in UNITS] + [('p', P) for P in chosen]
    n = len(basis)
    G = [[0]*n for _ in range(n)]
    for i in range(n):
        for k in range(n):
            e = Hexp(basis[i][0], basis[i][1], basis[k][0], basis[k][1])
            G[i][k] = (e % 2) if ell == 2 else (e % 3)
    cols = [ucl[i][idx] for i in range(4)] + [pcl[P][idx] for P in chosen]
    return G, cols

G2, C2 = gram(2, B2P, 0); G3, C3 = gram(3, B3P, 1)
def det_mod(G, ell): return rank_mod(G, ell)
r2 = det_mod(G2, 2); r3 = det_mod(G3, 3)
anti2 = all((G2[i][k] + G2[k][i]) % 2 == 0 for i in range(8) for k in range(8))
anti3 = all((G3[i][k] + G3[k][i]) % 3 == 0 for i in range(8) for k in range(8))
print("T3 Gram ranks: h2 %d/8, h3 %d/8; antisymmetric %s %s" % (r2, r3, anti2, anti3))
check("T3-nondegenerate", r2 == 8 and r3 == 8 and anti2 and anti3)

# Gram in standard coordinates: G_std(c, d) = y(c)^T G y(d), y = basis coordinates
def coords(cols, vec, ell): return solve_mod(cols, vec, ell)
def pair(G, cols, ell, c, d):
    y = coords(cols, c, ell); z = coords(cols, d, ell)
    return sum(y[i]*G[i][k]*z[k] for i in range(8) for k in range(8)) % ell

# consistency of H_S with classes: test H_S(p, p') from residue symbols against Gram (non-basis primes)
nH = 0; badH = 0
sample = random.sample(plist, min(60, len(plist)))
for (P, g), (Q, h) in itertools.combinations(sample, 2):
    if P == Q: continue
    e = Hexp('p', P, 'p', Q)
    ok2 = (e % 2) == pair(G2, C2, 2, pcl[P][0], pcl[Q][0])
    ok3 = (e % 3) == pair(G3, C3, 3, pcl[P][1], pcl[Q][1])
    nH += 1
    if not (ok2 and ok3): badH += 1
for (P, g) in sample:
    for ui, u in enumerate(UNITS + [(4, 0, 0, 0)]):
        e = Hexp('u', u, 'p', P)
        cu = classes(u)
        ok2 = (e % 2) == pair(G2, C2, 2, cu[0], pcl[P][0])
        ok3 = (e % 3) == pair(G3, C3, 3, cu[1], pcl[P][1])
        nH += 1
        if not (ok2 and ok3): badH += 1
print("T3/T5 H_S(class,class) = tame-reciprocity value: %d pairs, %d mismatches" % (nH, badH))
check("T3-H-classfunction", badH == 0, str(badH))

# ---------- T4: Lagrangian complement B_S and allowed generators ----------
def lagrangian_complement(G, ell):
    # basis e_1..e_4 (units, isotropic), f_1..f_4 (primes). A_ik = G(e_i, f_k)
    A = [[G[i][4+k] % ell for k in range(4)] for i in range(4)]
    # inverse of A mod ell
    M = [A[i][:] + [1 if i == k else 0 for k in range(4)] for i in range(4)]
    for c in range(4):
        piv = next(i for i in range(c, 4) if M[i][c] % ell)
        M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], -1, ell); M[c] = [t*iv % ell for t in M[c]]
        for i in range(4):
            if i != c and M[i][c] % ell:
                f = M[i][c]; M[i] = [(M[i][k] - f*M[c][k]) % ell for k in range(8)]
    Ainv = [row[4:] for row in M]
    # f'_j = sum_k Ainv[k][j] f_k  (vectors in basis coordinates, length 8)
    fp = []
    for j in range(4):
        v = [0]*8
        for k in range(4): v[4+k] = Ainv[k][j] % ell
        fp.append(v)
    def Gv(u, v): return sum(u[i]*G[i][k]*v[k] for i in range(8) for k in range(8)) % ell
    for i in range(4):
        for j in range(4):
            assert Gv([1 if t == i else 0 for t in range(8)], fp[j]) == (1 if i == j else 0)
    fpp = []
    for i in range(4):
        v = fp[i][:]
        for j in range(i+1, 4):
            pass
        fpp.append(v)
    # x_{ij} = G(f'_i, f'_j) for i<j ; f''_j = f'_j + sum_i x_{ij} e_i
    for j in range(4):
        v = fp[j][:]
        for i in range(j):
            x = Gv(fp[i], fp[j])
            v[i] = (v[i] + x) % ell
        fpp[j] = v
    for i in range(4):
        for j in range(4):
            assert Gv(fpp[i], fpp[j]) == 0
    return fpp
F2 = lagrangian_complement(G2, 2); F3 = lagrangian_complement(G3, 3)
SHEAR = int(sys.argv[2]) if len(sys.argv) > 2 else 0
def shear(Fb, ell, rng):
    # other Lagrangian complements: f''_j + sum_i s_ij e_i with s symmetric (keeps isotropy, Lemma 3.4 'one may choose')
    S = [[0]*4 for _ in range(4)]
    for i in range(4):
        for j in range(i, 4):
            S[i][j] = S[j][i] = rng.randrange(ell)
    out = []
    for j in range(4):
        v = Fb[j][:]
        for i in range(4): v[i] = (v[i] + S[i][j]) % ell
        out.append(v)
    return out
if SHEAR:
    rng = random.Random(SHEAR)
    F2 = shear(F2, 2, rng); F3 = shear(F3, 3, rng)
    def Gv_(G, u, v, ell): return sum(u[i]*G[i][k]*v[k] for i in range(8) for k in range(8)) % ell
    assert all(Gv_(G2, F2[i], F2[j], 2) == 0 for i in range(4) for j in range(4))
    assert all(Gv_(G3, F3[i], F3[j], 3) == 0 for i in range(4) for j in range(4))
    print("B_S replaced by a sheared isotropic complement (seed %d)" % SHEAR)

def allowed_decomp(vec, cols, Fb, ell):
    """vec (standard coords) = sum_i a_i e_i + b ; return (a_1..a_4, b-coords in f'' basis)."""
    y = coords(cols, vec, ell)   # basis coordinates
    # express y in basis {e_1..e_4, f''_1..f''_4} (in basis coordinates)
    newcols = [[1 if t == i else 0 for t in range(8)] for i in range(4)] + Fb
    z = solve_mod(newcols, y, ell)
    return z[:4], z[4:]

def allowed(gen):
    c2, c3 = classes(gen)
    a2, b2 = allowed_decomp(c2, C2, F2, 2); a3, b3 = allowed_decomp(c3, C3, F3, 3)
    A = [((3*a2[k]) + (4*a3[k])) % 6 for k in range(4)]   # CRT: =a2 mod 2, =a3 mod 3
    for k in range(4): assert A[k] % 2 == a2[k] and A[k] % 3 == a3[k]
    eta = tuple(Fr(t) for t in ONE)
    for k in range(4):
        if A[k]: eta = mulF(eta, powF(tuple(Fr(t) for t in UNITS[k]), A[k]))
    agen = mulF(tuple(Fr(t) for t in gen), invF(eta))
    c2b, c3b = classes(agen)
    a2b, bb2 = allowed_decomp(c2b, C2, F2, 2); a3b, bb3 = allowed_decomp(c3b, C3, F3, 3)
    assert all(t == 0 for t in a2b) and all(t == 0 for t in a3b)
    return agen, (tuple(bb2), tuple(bb3))

AL = {}
for P, g in plist:
    AL[P] = allowed(g)
print("allowed generators computed for %d primes (%.1fs)" % (len(AL), time.time() - t0))

# T3b Lemma 3.5 (eq. v-unit-projection): pairing E x B_S -> mu6 perfect: N_S^-1 sum_eps xi_eps(b) conj xi_eps(c) = 1_{b=c}
# xi_eps(b) = H_S(eps,b)^-1; eps runs over all 6^4 = 1296 classes of U_S/U_S^6, b over B_S classes realized by allowed primes
realized = {}
for P, g in plist:
    realized.setdefault(AL[P][1], P)
eps_cls = []
for A in itertools.product(range(6), repeat=4):
    v2 = [sum(A[k]*ucl[k][0][i] for k in range(4)) % 2 for i in range(8)]
    v3 = [sum(A[k]*ucl[k][1][i] for k in range(4)) % 3 for i in range(8)]
    eps_cls.append((v2, v3))
def xi_exp(ec, P):
    a = AL[P][0]; c2, c3 = classes(a)
    h2 = pair(G2, C2, 2, ec[0], c2); h3 = pair(G3, C3, 3, ec[1], c3)
    e = (3*h2 + 4*h3) % 6          # exponent of H_S(eps, a) in Z/6 (CRT)
    return (-e) % 6
keys = list(realized.keys())[:24]
XI = {kk: [xi_exp(ec, realized[kk]) for ec in eps_cls] for kk in keys}
badorth = 0; northo = 0
for b in keys:
    for c in keys:
        s_ = sum(Z6[(XI[b][i] - XI[c][i]) % 6] for i in range(len(eps_cls)))/len(eps_cls)
        northo += 1
        if abs(s_ - (1 if b == c else 0)) > 1e-9: badorth += 1
print("T3b unit projection orthogonality over 1296 unit classes: %d (b,c) pairs of realized B_S classes (%d distinct classes realized by %d primes), %d failures" % (northo, len(realized), len(plist), badorth))
check("T3b", badorth == 0)

# T4 reciprocity for allowed generators
nrec = badrec = 0
sampP = random.sample(plist, min(150, len(plist)))
for (P, _), (Q, _) in itertools.combinations(sampP, 2):
    a, b = AL[P][0], AL[Q][0]
    e1 = resexp(RFs[Q], a); e2 = resexp(RFs[P], b)  # chi_b(a) = (a/Q), chi_a(b) = (b/P)
    nrec += 1
    if e1 != e2: badrec += 1
# control: arbitrary (non-allowed) generators
nctl = badctl = 0
for (P, g), (Q, h) in itertools.combinations(sampP[:60], 2):
    nctl += 1
    if resexp(RFs[Q], g) != resexp(RFs[P], h): badctl += 1
print("T4 chi_a(b) = chi_b(a), allowed: %d pairs, %d failures; control with O_F generators: %d of %d pairs fail" % (nrec, badrec, badctl, nctl))
check("T4-allowed-reciprocity", badrec == 0, str(badrec))

# T5: chi_a(eps) depends only on class of a in B_S; xi_eps(a) = H_S(eps,a)^-1
nT5 = badT5 = 0
for P, g in plist:
    a, (b2, b3) = AL[P]
    for u in UNITS + [(4, 0, 0, 0)]:
        e = resexp(RFs[P], u)
        cu = classes(u)
        h2 = pair(G2, C2, 2, cu[0], classes(a)[0]); h3 = pair(G3, C3, 3, cu[1], classes(a)[1])
        nT5 += 1
        if (-e) % 2 != h2 or (-e) % 3 != h3: badT5 += 1
print("T5 chi_a(eps) = H_S(eps,a)^-1 (eps in U_S gens and 4): %d checks, %d failures" % (nT5, badT5))
check("T5", badT5 == 0, str(badT5))

# ---------- T6: gamma_3 on allowed generators is a character of B_S (order 2); G_0 order | 6 ----------
g3a = {}; G0a = {}
for P, g in plist:
    a, (b2, b3) = AL[P]
    R = RFs[P]
    g3 = gauss(R, a, 3)
    g3a[P] = g3
    c4 = Z6[(-resexp(R, (4, 0, 0, 0))) % 6]
    G0a[P] = c4*g3
# fit lambda: gamma_3 = (-1)^{L.b2}
def fit_and_test(vals, keyf, ell, name):
    # vals: P -> complex root of unity of order | ell ; key: coordinate vector over F_ell (length 4)
    pts = [(keyf(P), vals[P]) for P, _ in plist]
    # exponent: value = exp(2 pi i e / ell)
    def ex(z):
        ang = cmath.phase(z)/(2*math.pi)*ell
        e = int(round(ang)) % ell
        assert abs(z - cmath.exp(2j*math.pi*e/ell)) < 1e-7, (z, ell)
        return e
    pts = [(k, ex(z)) for k, z in pts]
    # choose 4 independent keys
    basis = []; bvals = []
    for k, e in pts:
        if rank_mod(basis + [list(k)], ell) > len(basis):
            basis.append(list(k)); bvals.append(e)
        if len(basis) == 4: break
    if len(basis) < 4:
        return None
    # solve L . k = e
    Lv = solve_mod([[basis[i][c] for i in range(4)] for c in range(4)], bvals, ell) if False else None
    # direct: matrix rows = basis keys; solve K L = e
    M = [basis[i] + [bvals[i]] for i in range(4)]
    for c in range(4):
        piv = next(i for i in range(c, 4) if M[i][c] % ell)
        M[c], M[piv] = M[piv], M[c]
        iv = pow(M[c][c], -1, ell); M[c] = [t*iv % ell for t in M[c]]
        for i in range(4):
            if i != c and M[i][c] % ell:
                f = M[i][c]; M[i] = [(M[i][k] - f*M[c][k]) % ell for k in range(5)]
    Lv = [M[i][4] for i in range(4)]
    bad = sum(1 for k, e in pts if sum(Lv[i]*k[i] for i in range(4)) % ell != e)
    return Lv, bad, len(pts)

res6 = fit_and_test(g3a, lambda P: AL[P][1][0], 2, "lambda")
print("T6 gamma_3(a) = lambda([a]) quadratic character of B_S: L = %s, %d failures of %d primes" % (res6[0], res6[1], res6[2]))
check("T6-lambda", res6[1] == 0)
# G_0: order | 6 -> test 2-part and 3-part separately: G0^3 is +-1 char of b2, G0^2 is mu3 char of b3
res6b = fit_and_test({P: G0a[P]**3 for P in G0a}, lambda P: AL[P][1][0], 2, "G0^3")
res6c = fit_and_test({P: G0a[P]**2 for P in G0a}, lambda P: AL[P][1][1], 3, "G0^2")
print("T6 G_0 = conj(chi(4)) gamma_3 character of B_S: 2-part L=%s fails %d; 3-part L=%s fails %d (of %d)" % (res6b[0], res6b[1], res6c[0], res6c[1], res6b[2]))
check("T6-G0", res6b[1] == 0 and res6c[1] == 0)
# G_0 values are sixth roots of unity
check("T6-G0-order6", all(abs(G0a[P]**6 - 1) < 1e-8 for P in G0a))

# ---------- T7: composite identity a_0 gamma_1 = mu G_0, mu gamma_{-1} = a_0 conj(G_0) ----------
def gauss_composite(Plist, gen, j):
    """full Gauss sum mod n = gen (squarefree, primes Plist), computed by CRT decomposition of x mod n:
    gamma_j(n) = q_n^{-1/2} sum_{x mod n} chi_n(x)^j e(x/n), e(x/n) = prod_i psi_{P_i}(-x/n)."""
    tot = 1+0j
    for P in Plist:
        R = RFs[P]; ell = R.ell
        w = R.red(tuple(ell*t for t in invF(gen)))   # ell/n mod P
        # chi_n(x) = prod_i chi_{P_i}(x); for x = (x_P) via CRT the other factors chi_{P_i}(x) depend on x mod P_i.
        s = 0j; y = R.one
        for m in range(R.q - 1):
            e = (m*R.k0*j) % 6
            s += Z6[e]*cmath.exp(-2j*math.pi*R.tr(R.mul(y, w))/ell)
            y = R.mul(y, R.g)
        tot *= s/math.sqrt(R.q)
    return tot
def Lam_ideal(Plist):
    z = 1+0j
    for P in Plist: z *= Lambda(Ps[P])
    return z
nT7 = badT7 = 0; maxe7 = 0.0
small = [P for P, _ in plist[:40]]
for r in (2, 3):
    for combo in itertools.combinations(small, r):
        if r == 3 and random.random() > 0.15: continue
        gen = tuple(Fr(t) for t in ONE)
        for P in combo: gen = mulF(gen, AL[P][0])
        g1 = gauss_composite(combo, gen, 1); g2 = gauss_composite(combo, gen, 2)
        gm1 = gauss_composite(combo, gen, 5)  # gamma_{-1} = chi^{-1} = chi^5
        a0 = Lam_ideal(combo).conjugate()*g2
        mu = (-1)**r
        G0 = 1+0j
        for P in combo: G0 *= G0a[P]
        e1 = abs(a0*g1 - mu*G0); e2 = abs(mu*gm1 - a0*G0.conjugate())
        # CRT formula of Liu: gamma_j(c1 c2) = gamma_j(c1) gamma_j(c2) chi_c1(c2)^j chi_c2(c1)^j
        maxe7 = max(maxe7, e1, e2); nT7 += 1
        if e1 > 1e-7 or e2 > 1e-7: badT7 += 1
print("T7 a0*gamma1 = mu*G0 and mu*gamma_-1 = a0*conj(G0) on squarefree allowed n (2 and 3 primes): %d n, %d failures, max err %.1e" % (nT7, badT7, maxe7))
check("T7", badT7 == 0)

# T7b: direct check of the CRT factor for two primes against brute force over R/nR (small norms)
nB = badB = 0
for (P, Q) in itertools.combinations([P for P, _ in plist[:8]], 2):
    a, b = AL[P][0], AL[Q][0]
    n = mulF(a, b)
    for j in (1, 2, 3):
        lhs = gauss_composite([P, Q], n, j)
        rhs = gauss(RFs[P], a, j)*gauss(RFs[Q], b, j)*Z6[(j*(resexp(RFs[P], b) + resexp(RFs[Q], a))) % 6]
        nB += 1
        if abs(lhs - rhs) > 1e-8: badB += 1
print("T7b CRT factor chi_c1(c2)^j chi_c2(c1)^j: %d checks, %d failures" % (nB, badB))
check("T7b", badB == 0)

# ---------- T8: Lemma 3.5 unit substitution for prime U ----------
nT8 = badT8 = 0
for P, _ in plist[:30]:
    a = AL[P][0]; R = RFs[P]; ell = R.ell
    for u in UNITS:
        for h in [(1, 0, 0, 0), (2, 1, 0, 0), (0, 1, 1, 3)]:
            if R.red(h) == R.zero: continue
            def G(hh):
                w = R.red(tuple(ell*t for t in mulF(tuple(Fr(t) for t in hh), invF(a))))
                s = 0j; y = R.one
                for m in range(R.q - 1):
                    s += Z6[(m*R.k0) % 6]*cmath.exp(-2j*math.pi*R.tr(R.mul(y, w))/ell)
                    y = R.mul(y, R.g)
                return s
            lhs = G(mul(u, h)); xi = Z6[resexp(R, u) % 6]
            rhs = xi.conjugate()*G(h)
            nT8 += 1
            if abs(lhs - rhs) > 1e-7*math.sqrt(R.q): badT8 += 1
print("T8 G(U, eps h) = conj(xi_eps(U)) G(U, h): %d checks, %d failures" % (nT8, badT8))
check("T8", badT8 == 0)

# ---------- controls (the tests have power) ----------
# C1: conjugate convention gamma_2^3 = -conj(Lambda) fails
c1 = sum(1 for P, gen in plist if abs(gam[P][1]**3 + Lambda(gen).conjugate()) > 1e-8)
# C2: gamma_3 of O_F generators is NOT a character of the full square class (V_S is not isotropic)
pts = [(pcl[P][0], gam[P][2]) for P, _ in plist]
def ex2(z): return 0 if z.real > 0 else 1
rk_keys = rank_mod([list(k) for k, z in pts], 2)
rk_aug = rank_mod([list(k) + [ex2(z)] for k, z in pts], 2)
c2 = rk_aug - rk_keys   # 0 iff gamma_3 is a quadratic character of the square class of the O_F generator
print("controls: C1 gamma2^3 = -conj(Lambda) fails at %d of %d primes; C2 gamma_3(O_F generator) vs its V_2 class: rank(keys)=%d, rank(keys|value)=%d (not a character iff larger)" % (c1, len(plist), rk_keys, rk_aug))
check("C1-power", c1 > 0); check("C2-power", c2 > 0)
print("elapsed %.1fs" % (time.time() - t0))
print("FAILURES:", len(FAIL))
for f in FAIL[:20]: print("  ", f)
print("RESULT:", "PASS" if not FAIL else "FAIL")
