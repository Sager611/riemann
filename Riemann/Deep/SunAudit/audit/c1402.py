import sys, math
from math import gcd, isqrt
sys.path.insert(0, '.')
from lib import *
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
which = sys.argv[2] if len(sys.argv) > 2 else 'all'
PL = primes_upto(6 * 10**7)
PS = set(PL)
MAXP = PL[-1]
def isp(x):
    if x < 2: return False
    if x <= MAXP: return x in PS
    return mr(x)
def mr(n):
    if n < 2: return False
    for p in [2,3,5,7,11,13,17,19,23,29,31,37]:
        if n % p == 0: return n == p
    d = n - 1; s = 0
    while d % 2 == 0: d //= 2; s += 1
    for a in [2,3,5,7,11,13,17,19,23,29,31,37]:
        x = pow(a, d, n)
        if x == 1 or x == n - 1: continue
        for _ in range(s - 1):
            x = x * x % n
            if x == n - 1: break
        else: return False
    return True
PHI = list(range(4 * N + 10))
for i in range(2, len(PHI)):
    if PHI[i] == i:
        for j in range(i, len(PHI), i): PHI[j] -= PHI[j] // i
SIG = [0] * (4 * N + 10)
for i in range(1, len(SIG)):
    for j in range(i, len(SIG), i): SIG[j] += i
pk = [None] + PL  # p_k
def pi_(x):
    import bisect
    return bisect.bisect_right(PL, x)
def practical(n):
    if n == 1: return True
    if n % 2: return False
    m = n; p = 2; prod = 1
    sig = 1
    while m > 1:
        if m % p == 0:
            if p > 1 + sig: return False
            e = 1
            while m % p == 0: m //= p; e *= p
            sig *= (e * p - 1) // (p - 1)
        p += 1
        if p * p > m and m > 1:
            if m > 1 + sig: return False
            break
    return True
PRAC = [False] + [practical(n) for n in range(1, 4 * N + 10)]
def issq(x): return x >= 0 and isqrt(x) ** 2 == x
def tri(x): return x >= 0 and issq(8 * x + 1)
def is_pr_root(g, p):
    if p == 2: return g % 2 == 1
    g %= p
    if g == 0: return False
    m = p - 1; d = 2; fs = []
    while d * d <= m:
        if m % d == 0:
            fs.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: fs.append(m)
    return all(pow(g, (p - 1) // f, p) != 1 for f in fs)
R = {}
def conj(name, start, exc=(), lim=None):
    def deco(f):
        R[name] = (f, start, set(exc), lim)
        return f
    return deco
P_upto = lambda n: [p for p in PL if p < n] if n < 10**6 else None
import bisect
def primes_below(n): return PL[:bisect.bisect_left(PL, n)]


import array
FACT = [1]
for _i in range(1, 3000): FACT.append(FACT[-1] * _i)
PIARR = None
LIMPI = 6 * 10**7
def build_pi():
    global PIARR
    arr = array.array('i', [0]) * (LIMPI + 1)
    c = 0; j = 0
    nxt = PL[0]
    idx = 0
    # cumulative count
    flags = bytearray(LIMPI + 1)
    for p in PL: flags[p] = 1
    s = 0
    for x in range(LIMPI + 1):
        s += flags[x]; arr[x] = s
    PIARR = arr
build_pi()
def pi(x): return PIARR[x] if x <= LIMPI else None
TW = array.array('i', [0]) * (2 * 10**6 + 1)
s = 0
for x in range(2 * 10**6 + 1):
    if x >= 5 and isp(x) and isp(x - 2): s += 1
    TW[x] = s
def pi2(x): return TW[x]
pp = [1] + [0] * 3000
for k in range(1, 3001):
    for t in range(k, 3001): pp[t] += pp[t - k]
qq = [1] + [0] * 3000
for k in range(1, 3001):
    for t in range(3000, k - 1, -1): qq[t] += qq[t - k]
qb = [pp[i] - qq[i] for i in range(3001)]

@conj('2.1ia', 2)
def _(n): return any(isp(pi(k*n)) for k in range(1, n + 1))
@conj('2.1ib', 1)
def _(n): return any(isp(pi(k*n)) for k in range(1, int(3*math.sqrt(n) + 3) + 1) if k < 3*math.sqrt(n) + 3)
for d_, n0 in [(0, 5), (1, 3), (-1, 6)]:
    def mkd(d_):
        def f(n): return any(isp(k*k + k - 1) and isp(pi(k*n) + d_) for k in range(1, n))
        return f
    R['2.1ii d=%d' % d_] = (mkd(d_), n0 + 1, set(), None)
@conj('2.2i', 1)
def _(n): return n == 1 or any(pi(k*n) % n == 0 for k in range(1, pk[n]))
@conj('2.2ii', 1, lim=300)
def _(n): return len(set(pi(k*n) % n for k in range(1, 2*pk[n] + 1))) == n
@conj('2.3', 2, lim=400)
def _(n): return any(pi(k*n) % pi(j*n) == 0 for j in range(1, n) for k in range(j + 1, n + 1) if k % j == 1 % j and pi(j*n) > 0)
@conj('2.4i', 1)
def _(n): return any(issq(pi(k*n)) for k in range(1, pk[n]))
@conj('2.4ii', 1)
def _(n): return any(issq(pi2(k*n)) for k in range(1, n + 1))
@conj('2.5i', 6)
def _(n): return any(isp(k*n + pi(k*n)) for k in range(1, n))
@conj('2.5ii', 1, lim=600)
def _(n): return any(isp(pk[k*n] - pi(k*n)) for k in range(1, n + 1))
@conj('2.6i', 3)
def _(n): return any(isp(pi(pi((p - 1)*n))) for p in primes_below(n + 1))
@conj('2.6iia', 1)
def _(n): return any(issq(pi(pi(k*n))) for k in range(1, n + 1))
@conj('2.6iib', 1)
def _(n): return any(tri(pi(pi(k*n))) for k in range(1, (n + 1)//2 + 1))
@conj('2.7i', 5, lim=1500)
def _(n):
    # pi(kn)^{1/k} > pi((k+1)n)^{1/(k+1)} for all k=1..n
    return all(pi(k*n) > 0 and math.log(pi(k*n)) / k > math.log(pi((k+1)*n)) / (k + 1) for k in range(1, n + 1))
@conj('2.7ii', 1)
def _(n): return any(issq(pi((k+1)*n) - pi(k*n)) for k in range(0, n))
@conj('2.8ia', 4)
def _(n): return any(isp(pi(p*n) - pi((p-1)*n)) for p in primes_below(n))
@conj('2.8ib', 4)
def _(n): return any(isp(pi((p+1)//2*n) - pi((p-1)//2*n)) for p in primes_below(n + 1) if p > 2)
@conj('2.8ii', 4)
def _(n): return any(isp(pi(k*n) - pi((k-1)*n)) and isp(pi((k+1)*n) - pi(k*n)) for k in range(1, n))
@conj('2.9i', 2)
def _(n): return any(pi((k+2)*n) - pi((k+1)*n) == pi((k+1)*n) - pi(k*n) for k in range(1, n))
@conj('2.9ii', 5, lim=800)
def _(n): return any(pi((k+1)*n) - pi(k*n) == pi((k+2)*n) - pi((k+1)*n) == pi((k+3)*n) - pi((k+2)*n) for k in range(1, pk[n]))
def c210(n): return len(set(pi((k+1)*n) - pi(k*n) for k in range(0, n)))
@conj('2.10 (>= sqrt(n-1), eq only 2,26)', 1)
def _(n):
    c = c210(n)
    if c * c < n - 1: return False
    if c * c == n - 1 and n not in (2, 26): return False
    return True
@conj('2.11', 2)
def _(n): return any(pi(p + 2*n) - pi(p + n) == pi(p + n) - pi(p) > 0 for p in primes_below(pk[n] + 1))
@conj('2.12ia', 11)
def _(n): return any(isp(pi(k*(n-k))) for k in range(1, n))
@conj('2.12ib', 11)
def _(n): return any(isp(pi(k*k*(n-k))) for k in range(1, n))
@conj('2.12iia', 10)
def _(n): return any(isp(pi(2*k) - pi(k)) and isp(pi(2*(n-k)) - pi(n-k)) for k in range(1, n))
@conj('2.12iib', 5)
def _(n): return any(isp(pi(2*k) - pi(k)) and issq(pi(2*(n-k)) - pi(n-k)) for k in range(1, n))
@conj('2.13i', 6)
def _(n): return any(isp(pi2(k*(n-k))) for k in range(1, n))
@conj('2.13ii', 9)
def _(n): return any(isp(pi2(k*(n-k)) - 1) and isp(pi2(k*(n-k)) + 1) for k in range(1, n))
@conj('2.14i', 1, lim=200)
def _(n): return len(set(pi(k*k) % n for k in range(1, 2*pk[n+1] - 2))) == n
@conj('2.14ii', 3)
def _(n): return math.log(pi(n*n)) / n > math.log(pi((n+1)**2)) / (n + 1)
@conj('2.14iii', 1, [25, 35, 44, 46, 105])
def _(n): return any(isp(x) for x in range(pi(n*n), pi((n+1)**2) + 1))
@conj('2.15i', 9)
def _(n): return any(isp(pi(k)) and isp(pi(k*k)) for k in range(n + 1, 2*n))
@conj('2.16ia', 2)
def _(n): return any(isp(pi(n + k*k)) for k in range(1, n))
@conj('2.16ii', 5)
def _(n): return any(isp(n + pi(k*k)) for k in range(1, n + 1))
@conj('2.16iiia', 1, [1, 2, 3, 4, 6, 12])
def _(n): return any(isp(n*n + pi(k*k)) for k in range(2, n))
@conj('2.16iiib', 5)
def _(n): return any(isp(pi(n*n) + pi(k*k)) for k in range(2, n))
@conj('2.16iiic', 2)
def _(n): return any(isp(pi((k+1)**2) - pi(k*k)) and isp(pi(n*n) - pi(k*k)) for k in range(1, n))
@conj('2.17ia', 5)
def _(n): return any(isp(p*n + pi(p)) for p in primes_below(n))
@conj('2.17ib', 1)
def _(n): return any(isp(p*n + pi(p)) for p in primes_below(math.sqrt(2*n) * math.log(5*n)))
@conj('2.17ii', 3)
def _(n): return any(isp(2*pi(p) - (-1)**n) and isp(p*n + ((-1)**n - 3)//2) for p in primes_below(n + 1))
@conj('2.18i', 4)
def _(n):
    ops = [q for q in primes_below(n + 1) if q > 2]
    return any(isp(p) and p > 2 and p <= n for q in ops for p in [n - q + pi(q)])
@conj('2.18ii', 8)
def _(n): return any(isp(n + p - pi(p)) for p in primes_below(n))
@conj('2.19i', 3)
def _(n): return any(isp(pi(p)) and isp(2*n - p) for p in primes_below(2*n))
@conj('2.20ia', 5)
def _(n): return any(isp(pi(n-p)) and isp(2*pi(n-p) + 1) for p in primes_below(n))
@conj('2.20ib', 9)
def _(n): return any(isp(pi(n-p) - 1) and isp(pi(n-p) + 1) for p in primes_below(n))
@conj('2.20iia', 5)
def _(n): return any(isp(3*pi(n-p) - 1) and isp(3*pi(n-p) + 1) and isp(3*pi(n-p) + 5) for p in primes_below(n))
@conj('2.20iib', 9)
def _(n): return any(isp(3*pi(n-p) - 1) and isp(3*pi(n-p) + 1) and isp(3*pi(n-p) - 5) for p in primes_below(n))
SGc = array.array('i', [0]) * (10**6 + 1)
s = 0
for x in range(10**6 + 1):
    if isp(x) and isp(2*x + 1): s += 1
    SGc[x] = s
@conj('2.21i', 5)
def _(n): return any(isp(SGc[n-p]) and isp(2*SGc[n-p] + 1) for p in primes_below(n))
@conj('2.21ii', 13)
def _(n): return any(isp(pi2(n - p + 2) if False else sum(1 for q in []) ) for p in []) or any(isp(r) and isp(r + 2) for p in primes_below(n) for r in [TWc(n - p)])
def TWc(x):
    # number of q <= x with q, q+2 prime
    return pi2(x + 2)
@conj('2.22ia', 3)
def _(n): return any(issq(pi(n-p)) for p in primes_below(n))
@conj('2.22ib', 3)
def _(n): return any(tri(pi(n-p)) for p in primes_below(n))
@conj('2.22ii', 3)
def _(n): return any(issq(pi(n+p)) for p in primes_below(pk[n] + 1))
def iscube(x):
    r = round(x ** (1/3)) if x > 0 else 0
    return any(c**3 == x for c in (r-1, r, r+1) if c >= 0)
@conj('2.23i', 12)
def _(n): return any(issq(SGc[n-p]) for p in primes_below(n))
@conj('2.23ii', 54)
def _(n): return any(iscube(SGc[n-p]) for p in primes_below(n))
SF = array.array('i', [0]) * (4 * 10**6 + 1)
s = 0
mu2 = bytearray([1]) * (4 * 10**6 + 1)
for d in range(2, 2001):
    for t in range(d*d, 4 * 10**6 + 1, d*d): mu2[t] = 0
for x in range(1, 4 * 10**6 + 1):
    s += mu2[x]; SF[x] = s
@conj('2.24i', 2)
def _(n): return any(isp(SF[(p-1)//2*n]) for p in primes_below(2*n) if p > 2)
@conj('2.24ii', 4)
def _(n): return any(isp(SF[n-p]) for p in primes_below(n))
def gauss_ideals(x):
    # primes p=1 mod 4 contribute 2 ideals (norm p), p=2 one, p=3 mod 4 ideal (p) norm p^2
    c = 0
    for p in PL:
        if p > x: break
        if p == 2: c += 1
        elif p % 4 == 1: c += 2
        elif p * p <= x: c += 1
    return c
@conj('2.25', 5, lim=300)
def _(n): return any(isp(g) and g % 4 == 1 for k in range(1, n + 1) for g in [gauss_ideals(k*n)])
@conj('3.1', 3)
def _(n): return any(isp(2*n - q) and isp(pk[q + 2] + 2) for q in primes_below(2*n))
@conj('3.2', 3)
def _(n): return any(isp(pk[k] + 2) and isp(pk[pk[n-k]] + 2) for k in range(1, n))
@conj('3.3', 3)
def _(n): return any(isp(6*k-1) and isp(6*k+1) and isp(pk[n-k] + 2) for k in range(1, n))
@conj('3.4', 2)
def _(n): return any(isp(pk[k]**2 - 2) and isp(pk[n-k]**2 - 2) and isp(pk[pk[n-k]]**2 - 2) for k in range(1, n))
@conj('3.5ii', 1, lim=200)
def _(n): return any(pk[(k+2)*n] - pk[(k+1)*n] == pk[(k+1)*n] - pk[k*n] for k in range(1, 3*pk[n] + 9) if (k+2)*n < len(pk))
@conj('3.6ia', 7)
def _(n): return any(isp(pk[k*n] + 2) for k in range(1, n + 1))
@conj('3.6ib', 7, lim=500)
def _(n): return any(isp(pk[k*k*n] + 2) for k in range(1, n + 1) if k*k*n < len(pk))
@conj('3.6ic', 1)
def _(n): return any(isp(pk[k*n] + 2) for k in range(1, int(3*math.sqrt(n) + 6) + 1) if k < 3*math.sqrt(n) + 6)
@conj('3.6ii', 1)
def _(n): return any(isp(2*k + 1) and isp(pk[k*n]**2 - 2) for k in range(1, n + 1))
@conj('3.8ia', 2)
def _(n): return any(isp(k*pk[n-k] + 1) for k in range(1, n))
@conj('3.8ib', 3)
def _(n): return any(isp(k*pk[n-k] - 1) for k in range(1, n))
@conj('3.8ii', 6)
def _(n): return any(isp(pk[k]*pk[n-k] - 6) for k in range(1, n))
@conj('3.9a', 7)
def _(n): return any(isp(pk[k] + pk[n-k] - 1) for k in range(1, n))
@conj('3.9b', 7)
def _(n): return any(isp(pk[k]**2 + pk[n-k]**2 - 1) for k in range(1, n))
@conj('3.10i', 4)
def _(n): return any(isp(pk[k]**2 + 4*pk[n-k]**2) for k in range(1, n))
@conj('3.10iia', 11)
def _(n): return any(isp(pk[k]**3 + 2*pk[n-k]**3) for k in range(1, n))
@conj('3.10iib', 11)
def _(n): return any(isp(pk[k]**3 + 2*pk[n-k]**2) for k in range(1, n))
@conj('3.11ia', 1, [1, 2, 3, 6])
def _(n): return any(isp(pk[q]**2 + (pk[n] - 1)**2) for q in primes_below(n))
@conj('3.11ib', 1, [1, 2, 9])
def _(n): return any(isp((pk[q] - 1)**2 + pk[n]**2) for q in primes_below(n))
@conj('3.11ii', 2)
def _(n): return any(isp(pk[n]**3 + 2*pk[k]**3) for k in range(1, n))
@conj('3.12i', 8, lim=600)
def _(n): return any(isp(2**k + pk[n-k]) for k in range(1, n))
@conj('3.12ii', 4, lim=400)
def _(n): return any(isp(FACT[k] + pk[n-k]) for k in range(1, n))
from math import comb
@conj('3.13i', 3, lim=500)
def _(n): return any(isp(comb(2*k, k) + pk[n-k]) for k in range(1, n) if n - k > k)
@conj('3.13ii', 5, lim=400)
def _(n): return any(isp(pk[n] + comb(pk[k] - 1, (pk[k] - 1)//2)) for k in range(2, n) if k < math.sqrt(n) * math.log(n))
@conj('3.14', 20, lim=500)
def _(n): return any((8*PHI[k] + PHI[n-k]) % 8 == 0 and isp(comb(2*m, m) + pk[m]) for k in range(1, n) for m in [(8*PHI[k] + PHI[n-k])//8] if m >= 1)
@conj('3.15i', 10)
def _(n): return any(isp(k + pk[n-k]) and isp(pk[k + pk[n-k]] - (k + pk[n-k]) + 1) for k in range(1, n))
@conj('3.15ii', 2)
def _(n): return any(isp(pk[pk[k]] - pk[k] + 1) and isp(pk[pk[2*k+1]] - pk[2*k+1] + 1) and isp(pk[pk[n-k]] - pk[n-k] + 1) for k in range(1, n))
@conj('3.15iii', 1, lim=400)
def _(n): return any(isp(pk[pk[k]] - pk[k] + 1) and pk[k*n] < len(pk) and isp(pk[pk[k*n]] - pk[k*n] + 1) for k in range(1, n + 1))
QSET = set(q for q in PL[:200000] if isp(pk[q] - q + 1))
@conj('3.17-odd', 7)
def _(n):
    if n % 2 == 0: return True
    Q = sorted(x for x in QSET if x < n)
    Qs = set(Q)
    return any((n - a - b) in Qs for a in Q for b in Q if a + b < n)
@conj('3.18', 32)
def _(n): return any(PHI[n-k] % 2 == 0 and isp(q) and isp(pk[q] - q - 1) and isp(pk[q] - q + 1) for k in range(1, n - 2) for q in [PHI[k] + PHI[n-k]//2 + 1])
@conj('3.19', 38)
def _(n):
    for k in range(1, n):
        if PHI[n-k] % 3: continue
        q = PHI[k] + PHI[n-k]//3 + 1
        if not isp(q): continue
        r = pk[q] - q + 1
        if isp(r) and r < len(pk) and isp(pk[r] - r + 1): return True
    return False
@conj('3.21i', 6)
def _(n): return any(isp(2*k - 1) and isp(pk[k*n] + k*n) for k in range(1, n))
@conj('3.21ii', 9)
def _(n): return any(isp(pk[k*(k+1)//2] + PHI[n - k*(k+1)//2]) for k in range(1, n) if k*(k+1)//2 < n)
@conj('3.22i', 101)
def _(n): return any(isp(k*PHI[k] + pk[n - k*k]) for k in range(1, isqrt(n) + 1) if n - k*k >= 1)
@conj('3.22ii', 7, [18])
def _(n): return any(isp(SIG[k*k] + pk[n - k*k] - 1) for k in range(1, isqrt(n) + 1) if n - k*k >= 1)
@conj('3.23i', 14)
def _(n): return any(isp(q + 2) and isp(pk[n-q] + q + 1) for q in primes_below(n))
@conj('3.23ii', 1, [1, 2, 3, 4, 6, 12])
def _(n): return any(isp(3*(pk[n-q] + q) - 1) and isp(3*(pk[n-q] + q) + 1) for q in primes_below(n))
INVSET = set(k for k in range(1, 5000) if isp(pow(k, -1, pk[k])) )
PRSET = set(k for k in range(1, 5000) if is_pr_root(k, pk[k]))
@conj('3.24i', 4)
def _(n): return any(a in INVSET and (n - a) in INVSET for a in range(1, n))
@conj('3.24ii', 2)
def _(n): return any(a in PRSET and (n - a) in PRSET for a in range(1, n))
@conj('3.25a', 7)
def _(n): return any(is_pr_root(p*n, pk[n]) for p in primes_below(n))
@conj('3.25b', 7)
def _(n): return any(is_pr_root(q*(n-q), pk[n]) for q in primes_below(n))
@conj('4.1ia', 1, lim=400)
def _(n): return any(isp(pp[n] + k) for k in range(1, n + 1))
@conj('4.1ib', 1, lim=400)
def _(n): return any(isp(qq[n] + k) for k in range(1, n + 1))
@conj('4.1iia', 2, lim=400)
def _(n): return any(isp(pp[n] + pp[k] - 1) for k in range(1, n))
@conj('4.1iib', 2, lim=400)
def _(n): return any(isp(pp[k] + qq[n]) for k in range(1, n))
@conj('4.1iic', 8, lim=400)
def _(n): return any(isp(n + pp[k]) for k in range(1, n))
@conj('4.2ia', 3, lim=500)
def _(n): return any(isp(2*pp[n-q] + 1) for q in primes_below(n))
@conj('4.2ib', 4, lim=500)
def _(n): return any(isp(2*pp[n-q] - 1) for q in primes_below(n))
@conj('4.2iia', 3, lim=500)
def _(n): return any(isp(qq[n-p] + 1) for p in primes_below(n))
@conj('4.2iib', 7, lim=500)
def _(n): return any(isp(qq[n-p] - 1) for p in primes_below(n))
@conj('4.4ia', 4, lim=300)
def _(n): return any(isp(pp[n+k] + 1) for k in range(1, n + 1))
@conj('4.4ib', 16, lim=300)
def _(n): return any(isp(pp[n+k] - 1) for k in range(1, n + 1))
@conj('4.5a', 8, lim=300)
def _(n): return any(isp(p) and p <= 3000 and isp(qq[p] - 1) for k in range(1, n) for p in [pk[k] + PHI[n-k]])
@conj('4.5b', 8, [15], lim=300)
def _(n): return any(isp(p) and p <= 3000 and isp(qq[p] + 1) for k in range(1, n) for p in [pk[k] + PHI[n-k]])
@conj('4.9ii', 3, lim=500)
def _(n): return any(isp(qq[k] + qb[n-k]) for k in range(1, n))
@conj('4.9iii', 2, lim=500)
def _(n): return any(isp(2**k - 1 + qq[n-k]) for k in range(1, n))
PSET_part = sorted(set(pp[1:200]))
QSET_part = sorted(set(qq[1:300]))
@conj('4.10i (p-th prime)', 1)
def _(n):
    p = pk[n]
    return any(is_pr_root(g, p) for g in PSET_part if g < p)
@conj('4.10ii (p-th prime, p>3)', 3)
def _(n):
    p = pk[n]
    return any(is_pr_root(g, p) for g in QSET_part if g < p)


@conj('3.16', 1, lim=1)
def _(n): return True
@conj('4.3i', 128, lim=1200)
def _(n): return any(PHI[n-k] % 2 == 0 and isp(pp[k + PHI[n-k]//2]) for k in range(1, n - 2))
@conj('4.3ii', 728, lim=1200)
def _(n): return any(PHI[n-k] % 2 == 0 and isp(PHI[k] + PHI[n-k]//2 + 1) and isp(pp[PHI[k] + PHI[n-k]//2]) for k in range(1, n - 2))
@conj('4.4ii', 6, lim=500)
def _(n): return any(isp(qq[PHI[k]*PHI[n-k]//4] + 1) for k in range(3, n - 2) if (PHI[k]*PHI[n-k]) % 4 == 0 and PHI[k]*PHI[n-k]//4 <= 3000)
@conj('4.6', 60, lim=1200)
def _(n): return any(PHI[n-k] % 4 == 0 and isp(m - 1) and isp(m + 1) and isp(qq[m] + 1) for k in range(1, n) for m in [PHI[k] + PHI[n-k]//4] if m <= 3000)
@conj('4.7i', 128, lim=1000)
def _(n): return any(PHI[n-k] % 6 == 0 and isp(r) and isp(pp[r] + qq[r]) for k in range(1, n) for r in [PHI[k] + PHI[n-k]//6 + 1] if r <= 3000)
@conj('4.7ii', 18, lim=1000)
def _(n): return any((4*PHI[k] + PHI[n-k]) % 8 == 0 and isp(pp[m]**2 + qq[m]**2) for k in range(1, n) for m in [(4*PHI[k] + PHI[n-k])//8] if 1 <= m <= 3000)
@conj('4.8i', 99, lim=1000)
def _(n): return any((6*PHI[k] + PHI[n-k]) % 12 == 0 and isp(p) and isp(qb[p]) for k in range(1, n) for p in [(6*PHI[k] + PHI[n-k])//12 + 1] if p <= 3000)
@conj('4.8ii', 4, lim=600)
def _(n): return any(PHI[n-k] % 2 == 0 and isp(qq[m]**2 + qb[m]**2) for k in range(1, n - 2) for m in [k + PHI[n-k]//2] if m <= 3000)
@conj('4.9ia', 2, lim=300)
def _(n): return any(isp(k*pp[n]*qq[n]*qb[n] - 1) for k in range(1, n + 1))
@conj('4.9ib', 2, lim=300)
def _(n): return any(isp(2*pp[k]*pp[n]*qq[n]*qb[n] + 1) for k in range(1, n))
@conj('4.1iiia', 2, lim=300)
def _(n): return any(isp(k*pp[n]*(pp[n]-1) + 1) for k in range(1, n))
@conj('4.1iiib', 2, lim=300)
def _(n): return any(isp(pp[k]*pp[n]*(pp[n]-1) + 1) for k in range(1, n))
@conj('4.1iiic', 2, lim=300)
def _(n): return any(isp(pp[k]*pp[n]*(pp[n]+1) - 1) for k in range(1, n))

names = list(R) if which == 'all' else [w for w in R if w.startswith(which)]
for name in names:
    f, start, exc, lim = R[name]
    top = min(N, lim) if lim else N
    def safe(n):
        try:
            return bool(f(n))
        except Exception as e:
            if n >= start: print('  ERR', name, n, repr(e)[:80])
            return False
    fails = [n for n in range(1, top + 1) if not safe(n)]
    bad = [n for n in fails if n >= start and n not in exc]
    print(name, 'n<=%d' % top, 'COUNTEREXAMPLES:' if bad else 'ok', bad[:10], '| failures below threshold/exc:', [n for n in fails if n < start or n in exc][-6:], flush=True)
