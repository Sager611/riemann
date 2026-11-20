# alpext P0: the field of the cubic class-group character psi of F = Q(sqrt-23) lies in no CM field (exact).
# (1) x^3-x-1 has one real root -> K3 = Q(alpha) has signature (1,1): neither totally real nor CM; subfields of CM fields
#     are totally real or CM (lemma in PROBES.md), so no CM field contains K3 or the Hilbert class field H = F(alpha).
# (2) complex conjugation acts on the roots as a transposition; Z(S3) = 1, so it is not central in Gal(H/Q) (H not CM).
# (3) first imaginary quadratic discriminant with a class-group character psi, psi^2 != 1 (psi^c = psi^-1 != psi): D = -23.
# (4) Frobenius types of x^3-x-1 mod p (p < 2e5) agree with (-23/p) and with the forms x^2+xy+6y^2, 2x^2+-xy+3y^2 (p < 3000).
import sympy as sp
from math import gcd, isqrt
x = sp.symbols('x'); f = x**3 - x - 1
print('disc(x^3-x-1) =', sp.discriminant(f, x), '; real roots:', len(sp.real_roots(f)), '; all roots:', len(sp.Poly(f).all_roots()))
# (2) complex conjugation on the three roots: real root fixed, complex pair swapped
rts = [complex(r.evalf(30)) for r in sp.Poly(f).all_roots()]
perm = []
for r in rts:
    c = r.conjugate(); k = min(range(3), key=lambda j: abs(rts[j]-c)); perm.append(k)
print('complex conjugation on roots (index map):', perm, '-> fixed points:', sum(1 for i,k in enumerate(perm) if i==k))
import itertools
S3 = list(itertools.permutations(range(3)))
comp = lambda a,b: tuple(a[b[i]] for i in range(3))
c = tuple(perm)
print('c central in S3:', all(comp(c,g)==comp(g,c) for g in S3), '; centre of S3:', [g for g in S3 if all(comp(g,h)==comp(h,g) for h in S3)])
# (3) class groups of imaginary quadratic discriminants -D, D <= 400, via reduced forms; exponent > 2 <=> psi with psi^2 != 1
def forms(D):
    out = []
    for a in range(1, isqrt(D//3)+1):
        for b in range(-a+1, a+1):
            if (b*b + D) % (4*a): continue
            cc = (b*b + D)//(4*a)
            if cc < a or (cc == a and b < 0): continue
            if gcd(gcd(a, abs(b)), cc) != 1: continue
            out.append((a, b, cc))
    return out
def compose(f1, f2, D):  # Gauss composition (Shanks/Cohen alg 5.4.7, primitive forms)
    a1,b1,c1 = f1; a2,b2,c2 = f2
    if a1 > a2: a1,b1,c1,a2,b2,c2 = a2,b2,c2,a1,b1,c1
    s = (b1+b2)//2; n = b2 - s
    def egcd(u,v):
        if v == 0: return (u,1,0)
        g,x1,y1 = egcd(v, u % v); return (g, y1, x1 - (u//v)*y1)
    if a2 % a1 == 0: y1 = 0; d = a1
    else:
        d, u, v = egcd(a2, a1); y1 = u
    if s % d == 0: y2 = -1; x2 = 0; d1 = d
    else:
        d1, x2, y2 = egcd(s, d); y2 = -y2
    v1 = a1//d1; v2 = a2//d1; r = (y1*y2*n - x2*c2) % v1
    b3 = b2 + 2*v2*r; a3 = v1*v2; c3 = (c2*d1 + r*(b2 + v2*r))//v1
    return reduce_form((a3, b3, c3), D)
def reduce_form(F, D):
    a,b,c = F
    while True:
        if c < a: a,b,c = c,-b,a
        if b > a or b <= -a:
            k = (a - b)//(2*a); b2 = b + 2*k*a; c = (b2*b2 + D)//(4*a); b = b2; continue
        if a == c and b < 0: b = -b
        return (a,b,c)
def is_fund(D):
    if D % 4 == 3:
        return all(D % (p*p) for p in range(3, isqrt(D)+1, 2))
    if D % 4 == 0:
        m = D//4
        return m % 4 in (1, 2) and all(m % (p*p) for p in range(3, isqrt(m)+1, 2))
    return False
first = None
for D in range(3, 401):
    if not is_fund(D): continue
    F = forms(D); h = len(F)
    I = reduce_form((1, D % 2, (D % 2 + D)//4), D)
    exp2 = all(reduce_form(compose(g, g, D), D) == I for g in F)
    if not exp2 and first is None: first = (D, h, F)
    if D <= 23 or (not exp2 and D <= 60): print(f'D=-{D}: h={h}, class group exponent <= 2: {exp2}')
print('first fundamental -D with a class-group character of order >= 3:', -first[0], 'h =', first[1], 'forms', first[2])
# (4) Frobenius types
def leg(a, p): r = pow(a % p, (p-1)//2, p); return -1 if r == p-1 else r
def pmulmod(A, B, p):  # polys as lists low->high, mod x^3 - x - 1
    r = [0]*5
    for i,ai in enumerate(A):
        for j,bj in enumerate(B): r[i+j] += ai*bj
    for k in (4, 3):  # x^3 = x + 1
        if r[k]: r[k-2] += r[k]; r[k-3] += r[k]; r[k] = 0
    return [v % p for v in r[:3]]
def xpow(p):
    R = [1,0,0]; B = [0,1,0]; e = p
    while e:
        if e & 1: R = pmulmod(R, B, p)
        B = pmulmod(B, B, p); e >>= 1
    return R
def nroots(p): return sum(1 for t in range(p) if (t*t*t - t - 1) % p == 0)
N = 200000
sieve = bytearray([1])*(N+1); sieve[0]=sieve[1]=0
for i in range(2, isqrt(N)+1):
    if sieve[i]: sieve[i*i::i] = bytearray(len(range(i*i, N+1, i)))
bad = 0; cnt = {'id':0,'tr':0,'3c':0}
for p in range(5, N+1):
    if not sieve[p] or p == 23: continue
    L = leg(-23, p); triv = xpow(p) == [0,1,0]
    typ = 'tr' if L == -1 else ('id' if triv else '3c')
    if L == -1 and triv: bad += 1
    cnt[typ] += 1
    if p < 3000:
        nr = nroots(p)
        if {'id':3,'tr':1,'3c':0}[typ] != nr: bad += 1
        rep0 = any((4*p - 23*y*y) >= 0 and isqrt(4*p-23*y*y)**2 == 4*p-23*y*y and (isqrt(4*p-23*y*y) + y) % 2 == 0 for y in range(0, isqrt(4*p//23)+1))
        # x^2+xy+6y^2 = p  <=>  (2x+y)^2 + 23 y^2 = 4p
        if L == 1 and rep0 != (typ == 'id'): bad += 1
print(f'Frobenius types for 5 <= p <= {N}, p != 23: {cnt}; mismatches = {bad}')
tot = sum(cnt.values()); print('shares id/tr/3c = %.4f %.4f %.4f (Chebotarev 1/6, 1/2, 1/3)' % (cnt['id']/tot, cnt['tr']/tot, cnt['3c']/tot))
