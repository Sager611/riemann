import sys, math
from math import gcd, isqrt
sys.path.insert(0, '.')
from lib import *
N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
which = sys.argv[2] if len(sys.argv) > 2 else 'all'
PL = primes_upto(30 * N * N if N <= 300 else 4 * 10**7)
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

@conj('2.2', 7)
def _(n):
    n1 = 1 if n % 2 else 2
    return any(isp(n - (p*n1 - 1)) and isp(n + (p*n1 - 1)) for p in primes_below((n + n1 - 1) // n1) if p * n1 < n)
@conj('2.3', 6)
def _(n): return any(isp(6*n - p) and isp(6*n + p) for p in primes_below(n))
@conj('2.4', 12)
def _(n): return any(isp(p + 6) and isp(6*(n-p) - 1) and isp(6*(n-p) + 1) for p in primes_below(n))
@conj('2.5', 2)
def _(n): return any(isp(6*k - 1) and isp(12*k - 1) and isp(6*(n-k) - 1) and isp(6*(n-k) + 1) for k in range(1, n))
@conj('2.7i', 4)
def _(n): return any(isp(2*n - p) and isp(3*p - 10) and isp(3*(2*n - p) + 10) for p in primes_below(2*n))
@conj('2.7ii', 5, [76])
def _(n): return any(isp(2*n - p) and isp(3*p - 10) and isp(3*(2*n - p) - 10) for p in primes_below(2*n))
@conj('2.8', 9)
def _(n): return any(isp(q) and isp(2*n - 1 - 2*q) and isp((2*n - 1 - 2*q)**2 + 60*q*q) for q in range(2, n))
@conj('2.9', 4)
def _(n): return any(isp(2*p*p - 1) and isp(2*(n-p)**2 - 1) for p in primes_below(n))
@conj('2.10i', 4)
def _(n): return any(isp(k*n + 1) and isp(k*(n-k) - 1) for k in range(1, n))
@conj('2.10ii-odd', 3)
def _(n): return n % 2 == 0 or any(isp(k*n + 1) and isp(k*k + (n-k)**2) for k in range(1, n))
@conj('2.11a', 3)
def _(n): return any(isp(k*k*(n-k) - 1) for k in range(1, n))
@conj('2.11b', 5)
def _(n): return any(isp(k*k*(n-k) + 1) for k in range(1, n))
@conj('2.11c', 2)
def _(n): return any(isp((k*(n-k))**2 + k*(n-k) + 1) for k in range(1, n))
@conj('2.11d', 3)
def _(n): return any(isp((k*(n-k))**2 + k*(n-k) - 1) for k in range(1, n))
@conj('2.11e', 3)
def _(n): return any(isp(2*(k*(n-k))**2 - 1) for k in range(1, n))
@conj('2.12a', 2)
def _(n): return any(isp((p - 1)*n + 1) for p in primes_below(n + 1))
@conj('2.12b', 5)
def _(n): return any(isp(3*p + 8) and isp((p - 1)*n + 1) for p in primes_below(n))
@conj('2.12c', 6)
def _(n): return any(isp(3*p - 10) and isp((p - 1)*(n - p) - 1) for p in primes_below(n))
@conj('2.13', 2)
def _(n): return any(isp(2*p + 1) and isp((p - 1)**2 + (2*n - p)**2) for p in primes_below(2*n))
@conj('2.14i-odd', 3)
def _(n): return n % 2 == 0 or any(isp(x**4 + (n-x)**4) for x in range(1, n))
@conj('2.14ii', 11)
def _(n): return any(isp(p + 6) and isp(n*n + p*(n-p)) for p in primes_below(n))
@conj('2.14iii-odd', 3)
def _(n): return n % 2 == 0 or any(isp(k*k + (n-k)**2) and isp(k**3 + (n-k)**2) for k in range(1, n))
@conj('2.15', 2)
def _(n): return any(isp(2*n - p) and isp((p - 1)*(2*n - p + 1) - 1) for p in primes_below(2*n))
@conj('2.16', 8)
def _(n): return any(isp(2*p*(n-p) + 1) for p in primes_below(n))
@conj('2.17', 4)
def _(n): return any(isp((p + 1)*(n - p)//2 + 1) and ((p + 1)*(n - p)) % 2 == 0 for p in primes_below(n))
@conj('2.18', 1)
def _(n): return any(isp(n + k) and isp(n + k*k) for k in range(0, n + 1))
@conj('2.19', 1)
def _(n): return any(isp(n*n + k*k) for k in range(1, int(4 * math.sqrt(n + 1)) + 1))
@conj('2.20i', 6)
def _(n): return any(isp(p + 6) and isp(n + (n-p)**2) for p in primes_below(n))
@conj('2.20ii', 4)
def _(n): return any(isp(3*p - 4) and isp(n*n + (n-p)**2) for p in primes_below(n))
@conj('2.21i', 2)
def _(n): return any(isp(x + n*(n-x)) and isp(x*x + n*(n-x)**2) for x in range(1, n))
@conj('2.21iia', 3)
def _(n): return any(isp(n*x + (n-x)) and isp(n*x - (n-x)) for x in range(1, n))
@conj('2.21iib', 3)
def _(n): return any(isp(x*x + (n-2)*(n-x)**2) for x in range(1, n))
@conj('2.21iiia', 3)
def _(n): return any(isp(p**3 + n*(n-p)**2) for p in primes_below(n))
@conj('2.21iiib', 3)
def _(n): return any(isp(p + n*(n-p)) for p in primes_below(n))
@conj('2.22a', 2)
def _(n): return any(isp(2*x + 1) and isp(x*x + x + 1) and isp((n-x)**2 + (n-x) + 1) for x in range(1, n))
@conj('2.22b', 2)
def _(n): return any(isp(x*x + 1) and isp(4*(n-x)**2 + 1) for x in range(1, n))
@conj('2.22c', 2)
def _(n): return any(isp(4*x*x + 1) and isp(4*(n-x)**2 + 1) for x in range(1, n))
@conj('2.23', 1)
def _(n): return any(isp(x**3 + 2*(n-x)**3) for x in range(0, n + 1))
@conj('2.24', 5)
def _(n): return any(isp(2*n + k) and isp(2*n**3 + k**3) for k in range(1, n))
for m_ in [1, 2, 3, 4, 5, 6, 18]:
    def mk(m_):
        def f(n): return n % 2 == 0 or any(isp(x**m_ + 3*(n-x)**m_) for x in range(0, n + 1))
        return f
    R['2.25 m=%d odd n' % m_] = (mk(m_), 1, set(), 400 if m_ >= 5 else None)
@conj('3.1a', 2)
def _(n): return any(isp(x*(n-x)//gcd(x, n-x) + 1) for x in range(1, n))
@conj('3.1b', 4)
def _(n): return any(isp(x*(n-x)//gcd(x, n-x) - 1) for x in range(1, n))
@conj('3.2ia', 2)
def _(n): return any(isp(x*(x+1)//2 + (n-x)**2) for x in range(1, n))
@conj('3.2ib', 7)
def _(n): return any(isp(x*(x+1)//2 + (n-x)**4) for x in range(1, n))
@conj('3.2ii', 2)
def _(n): return any(isp(2*x + 1) and isp(x*(x+1)//2 + n - x) for x in range(1, n))
@conj('3.3', 1)
def _(n): return any(isp((n - x*(x+1)//2) * (n - x*(x+1)//2 + 1)//2 + 1) for x in range(0, n + 1) if x*(x+1)//2 <= n)
@conj('3.4', 3)
def _(n): return any(isp(2*k*(n - k*k) - 1) for k in range(1, isqrt(n) + 1) if n - k*k >= 1)
@conj('3.5', 2)
def _(n): return any(isp(2*(n - k*k)**2 - 1) for k in range(0, isqrt(n) + 1))
@conj('3.6i', 2, lim=600)
def _(n): return any(isp(n + 2**k - k) for k in range(1, n))
@conj('3.6ii', 4, lim=600)
def _(n): return any(isp(n + 2**k + k) for k in range(1, n))
@conj('3.7', 4)
def _(n):
    vals = [2**k - k for k in range(1, 40) if 2**k - k < n]
    return any(isp(n - a - b) for a in vals for b in vals)
for r_ in (1, 2):
    def mk2(r_, s):
        def f(n): return any(isp(2**k * (n-k)**r_ + s) for k in range(1, n))
        return f
    R['3.8 r=%d +1' % r_] = (mk2(r_, 1), 2, set(), 600)
    R['3.8 r=%d -1' % r_] = (mk2(r_, -1), 3, set(), 600)
FIB = [0, 1]
while len(FIB) < 3 * N + 10: FIB.append(FIB[-1] + FIB[-2])
LUC = [2, 1]
while len(LUC) < 3 * N + 10: LUC.append(LUC[-1] + LUC[-2])
@conj('3.9a', 2, lim=800)
def _(n): return any(isp(FIB[k] + n - k) for k in range(1, n))
@conj('3.9b', 2, lim=800)
def _(n): return any(isp(FIB[k] + 2*(n - k)) for k in range(1, n))
@conj('3.9c', 2, lim=800)
def _(n): return any(isp(FIB[k] + (n-k)*(n-k+1)) for k in range(1, n))
@conj('3.9ii', 2, lim=800)
def _(n): return any(isp(LUC[k] + n - k) for k in range(1, n))
FACT = [1]
for i in range(1, 3 * N): FACT.append(FACT[-1] * i)
@conj('3.10a', 2, lim=400)
def _(n): return any(isp(FACT[k]*(n-k) + 1) for k in range(1, n))
@conj('3.10b', 4, lim=400)
def _(n): return any(isp(FACT[k]*(n-k) - 1) for k in range(1, n))
@conj('3.10ii+', 2, lim=300)
def _(n): return any(isp(FACT[k]*n + 1) for k in range(1, n + 1))
@conj('3.10ii-', 2, lim=300)
def _(n): return any(isp(FACT[k]*n - 1) for k in range(1, n + 1))
@conj('3.17i', 7)
def _(n): return any(tri(pk[k] + pi_(n-k)) for k in range(1, n))
@conj('3.17ii', 11)
def _(n): return any(isp(pk[k] + pi_(n-k)) and isp(pk[k] + pi_(n-k) + 2) for k in range(1, n))
@conj('3.18', 3)
def _(n): return any(isp(pk[q] - q + 1) and isp(pk[pk[n-q]] - pk[n-q] + 1) for q in primes_below(n))
@conj('3.19', 7)
def _(n): return any(2*q < n and isp(pk[q] - q + 1) and isp(n - (1 + n % 2)*q) for q in primes_below(n))
@conj('3.20ia', 3)
def _(n): return any(isp(2*n - p) and isp(PHI[p+2] - 1) and isp(PHI[p+2] + 1) for p in primes_below(2*n))
@conj('3.20ib', 13)
def _(n): return any(isp(2*n - 1 - 2*p) and isp(PHI[p+1] - 1) and isp(PHI[p+1] + 1) for p in primes_below(n))
@conj('3.20ii', 24)
def _(n):
    c = 1 + n % 2
    return any(isp(n - c*p) and isp(PHI[p+1] - 1) and isp(PHI[n - c*p - 1] + 1) for p in primes_below(n // c + 1) if n - c*p >= 2)
@conj('3.21i', 12)
def _(n): return any(k != n-k and isp(PHI[k]-1) and isp(PHI[k]+1) and isp(PHI[n-k]-1) and isp(PHI[n-k]+1) for k in range(1, n))
@conj('3.21ii', 7)
def _(n): return any(isp(pk[k] + 2) and isp(PHI[n-k]-1) and isp(PHI[n-k]+1) for k in range(1, n))
@conj('3.21iiia', 6)
def _(n): return any(isp(pk[pk[k]] + 2) and isp(PHI[n-k]-1) and isp(PHI[n-k]+1) for k in range(1, n))
@conj('3.21iiib', 2)
def _(n): return any(isp(pk[pk[k]] + 2) and isp(6*(n-k)-1) and isp(6*(n-k)+1) for k in range(1, n))
@conj('3.21iv', 8, lim=1500)
def _(n): return any(pk[k] < len(pk) and pk[pk[k]] < len(pk) and isp(pk[pk[pk[k]]] - 2) and isp(PHI[n-k]-1) and isp(PHI[n-k]+1) for k in range(1, n))
@conj('3.21v', 9)
def _(n): return any(isp(3*k-1) and isp(3*k+1) and isp(PHI[n-k]-1) and isp(PHI[n-k]+1) for k in range(1, n))
@conj('3.21vi', 12)
def _(n): return any(isp(p+6) and isp(PHI[n-p]-1) and isp(PHI[n-p]+1) for p in primes_below(n))
@conj('3.22i', 2, [8])
def _(n): return any(isp(k + PHI[n-k]) and isp(2*n - k - PHI[n-k]) for k in range(1, n))
@conj('3.22ii', 3)
def _(n): return any(isp(k + PHI[n-k]) and isp(2*n + 1 - 2*(k + PHI[n-k])) for k in range(1, n))
@conj('3.23i', 2)
def _(n): return any(isp(6*k-1) and isp(6*k+1) and isp(k + PHI[n-k]) for k in range(1, n))
@conj('3.23ii', 4, [12])
def _(n): return any(isp(6*k-1) and isp(6*k+1) and PHI[n-k] % 2 == 0 and isp(k + PHI[n-k]//2) for k in range(1, n))
@conj('3.23iii', 6)
def _(n): return any(PHI[n-k] % 2 == 0 and issq(k + PHI[n-k]//2) for k in range(1, n))
KSET = set(k for k in range(1, 4*N) if isp(k*(k+1) - pk[k]))
@conj('3.24i', 4)
def _(n): return any(a in KSET and (n-a) in KSET for a in range(1, n))
@conj('3.24ii', 3)
def _(n): return any(a in KSET and tri(n-a) and n - a > 0 for a in range(1, n))
@conj('3.24iii', 4)
def _(n): return any(a in KSET and isp(n-a) and isp(pk[n-a] - (n-a) + 1) for a in range(1, n))
@conj('3.24iv', 8)
def _(n): return any(isp(pk[k] + PHI[n-k]) and pk[k] + PHI[n-k] < len(pk) and isp((pk[k] + PHI[n-k])*(pk[k] + PHI[n-k] + 1) - pk[pk[k] + PHI[n-k]]) for k in range(1, n))
@conj('3.25i', 8)
def _(n): return any(PHI[n-k] % 2 == 0 and isp(PHI[k] + PHI[n-k]//2 - 1) and PHI[k] + PHI[n-k]//2 - 1 > 2 and is_pr_root(2, PHI[k] + PHI[n-k]//2 - 1) for k in range(1, n))
@conj('3.25ii', 38)
def _(n): return any(isp(pk[k] + PHI[n-k]) and isp(2*(pk[k] + PHI[n-k]) + 1) and is_pr_root(2, pk[k] + PHI[n-k]) for k in range(1, n))
@conj('3.26i', 6, lim=500)
def _(n): return any(isp(2**PHI[k] + 2**PHI[n-k] - 1) for k in range(3, n - 2))
for a_, st in [(2, 6), (3, 6), (4, 6), (5, 6), (6, 6), (8, 6), (7, 18)]:
    def mk3(a_):
        def f(n): return any(PHI[n-k] % 2 == 0 and isp(a_**PHI[k] + a_**(PHI[n-k]//2) - 1) for k in range(3, n - 2))
        return f
    R['3.26ii a=%d' % a_] = (mk3(a_), st, set(), 250)
def iexp(k, m, dk, dm):
    num = PHI[k] * dm + PHI[m] * dk
    if num % (dk * dm): return None
    return num // (dk * dm)
def mk4(base, mult, dk, dm, add):
    def f(n):
        for k in range(1, n):
            e = iexp(k, n-k, dk, dm)
            if e is not None and isp(mult * base**e + add): return True
        return False
    return f
R['3.27ia'] = (mk4(2, 1, 2, 6, 3), 10, set(), 600)
R['3.27ib'] = (mk4(2, 1, 2, 6, -3), 14, set(), 600)
R['3.27iia'] = (mk4(2, 3, 2, 8, 1), 26, set(), 600)
R['3.27iib'] = (mk4(2, 3, 2, 12, -1), 15, set(), 600)
R['3.27iiia'] = (mk4(3, 2, 2, 12, 1), 27, set(), 400)
R['3.27iiib'] = (mk4(3, 2, 2, 12, -1), 38, set(), 400)
R['3.27iv'] = (mk4(2, 1, 1, 4, -5), 11, set(), 400)
@conj('3.29i', 9)
def _(n): return any(issq(PHI[i]*PHI[n-i]) for i in range(1, (n+1)//2) if i < n - i)
def phi_big(x):
    r = x; m = x; d = 2
    while d * d <= m:
        if m % d == 0:
            while m % d == 0: m //= d
            r -= r // d
        d += 1
    if m > 1: r -= r // m
    return r
@conj('3.29ii', 8)
def _(n): return any(issq(phi_big(i*(n-i)) + 1) for i in range(1, (n+1)//2) if i < n-i)
@conj('3.30i', 6)
def _(n): return any((PHI[k] + PHI[n-k]) % 2 == 0 and isp((PHI[k] + PHI[n-k])//2) for k in range(1, n))
@conj('3.30iia', 1, [1, 2, 3, 6])
def _(n): return any(issq(k*PHI[n-k] + 1) for k in range(1, n))
@conj('3.30iib', 5)
def _(n): return any(k < n-k and isp(k*PHI[n-k] - 1) and isp(k*PHI[n-k] + 1) for k in range(1, n))
@conj('3.30iii', 6)
def _(n): return any(isp(PHI[k]*PHI[n-k] - 1) and isp(PHI[k]*PHI[n-k] + 1) for k in range(1, n))
@conj('3.30iv', 5)
def _(n): return any(isp(k*PHI[k]*PHI[n-k] - 1) and isp(2*(k*PHI[k]*PHI[n-k] - 1) + 1) for k in range(1, n))
@conj('3.31i', 2)
def _(n): return any(isp(SIG[k*k] + PHI[n - k*k]) for k in range(1, isqrt(n) + 1) if k*k <= n - k*k)
@conj('3.31iia', 2)
def _(n): return any(isp(SIG[k]**2 + PHI[n-k]) for k in range(1, n))
@conj('3.31iib', 2)
def _(n): return any(isp(SIG[k] + PHI[n-k]**2) for k in range(1, n))
@conj('3.32', 3)
def _(n):
    if n % 2 == 0: return any(isp(n - SIG[k]) and n - SIG[k] > 2 for k in range(1, n))
    return any(isp(n - k*PHI[k]) for k in range(1, isqrt(n - 1) + 1) if k*k < n)
@conj('3.33ia', 9, [14])
def _(n): return any(jacobi(n, p) == 1 for p in PL[bisect.bisect_right(PL, n):bisect.bisect_left(PL, 2*n)] if p > 2) if True else False
@conj('3.33ib-nonsq', 1)
def _(n): return issq(n) or any(leg(n, p) == -1 for p in PL[bisect.bisect_right(PL, n):bisect.bisect_left(PL, 2*n)] if p > 2)
@conj('3.33iia', 6)
def _(n): return any(leg(2*n, p) == 1 for p in PL[bisect.bisect_right(PL, n):bisect.bisect_left(PL, 2*n)] if p > 2)
@conj('3.33iib', 7)
def _(n): return any(leg(-n, p) == -1 for p in PL[bisect.bisect_right(PL, n):bisect.bisect_left(PL, 2*n)] if p > 2)
@conj('3.34a', 1)
def _(n): return any(leg(n, p) == 1 for p in PL[bisect.bisect_right(PL, n*n):bisect.bisect_left(PL, (n+1)**2)] if p > 2)
@conj('3.34b', 2)
def _(n): return any(leg(n*(n+1), p) == 1 for p in PL[bisect.bisect_right(PL, n*n):bisect.bisect_left(PL, (n+1)**2)] if p > 2)
@conj('3.35', 400, [757, 1069, 1238])
def _(n):
    c = 1 + n % 2
    return any(q > 2 and isp(n - c*q) and n - c*q > 2 and leg(n - c*q, q) == 1 and leg(q, n - c*q) == 1 for q in primes_below(n // c + 1))
@conj('3.37', 6)
def _(n):
    c = 1 + n % 2
    for q in primes_below(n // 2 + 1):
        p = n - c*q
        if p > 2 and isp(p):
            if n % 2 and jacobi(q, n) == 1: return True
            if n % 2 == 0 and q > 2 and jacobi((q + 1)//2, n + 1) == 1: return True
            if n % 2 == 0 and q == 2 and jacobi((q + 1)//2 if (q+1) % 2 == 0 else 0, n + 1) == 1: return True
    return False
@conj('3.38i', 3)
def _(n): return any(isp(2*n - p) and PRAC[p+1] and PRAC[2*n - p - 1] for p in primes_below(2*n))
@conj('3.38ii', 9)
def _(n): return any(isp(2*n - 1 - 2*p) and PRAC[2*n - 1 - p] for p in primes_below(n))
def sand1(p): return isp(p) and PRAC[p-1] and PRAC[p+1]
def sand2(q): return PRAC[q] and isp(q-1) and isp(q+1)
@conj('3.39i', 9)
def _(n): return any(sand1(p) for p in range(n + 1, 2*n))
@conj('3.39ii', 7)
def _(n): return any(sand2(q) for q in range(n + 1, 2*n))
@conj('3.39iii', 232)
def _(n): return any(isp(p) and isp(p+2) and PRAC[p-1] and PRAC[p+1] for p in range(n + 1, 2*n - 1))
@conj('3.40i', 4)
def _(n): return any(sand1(p) and (isp(n-p) or PRAC[n-p]) for p in range(2, n))
@conj('3.43i-odd', 3)
def _(n): return n % 2 == 0 or any(isp(2*p + 1) and PRAC[n - p] for p in primes_below(n))
@conj('3.43ii', 3)
def _(n): return any(PRAC[q] and ((isp(n-q) and isp(n+q)) or (PRAC[n-q] and PRAC[n+q])) for q in range(1, n))
@conj('3.44i', 7)
def _(n): return any(n - p - 1 >= 1 and ((isp(n-p-1) and isp(n+p+1)) or (PRAC[n-p-1] and PRAC[n+p+1])) for p in primes_below(n))
@conj('3.44ii', 3)
def _(n): return any(n - p + 1 >= 1 and ((isp(n-p+1) and isp(n+p-1)) or (PRAC[n-p+1] and PRAC[n+p-1])) for p in primes_below(n))


if which == 'extra':
    out = []
    # 2.6
    L = 6000
    A = [x for x in range(1, L) if isp(6*x-1) and isp(6*x+1)]
    B = [x for x in range(1, L) if isp(6*x+1) and isp(6*x+5)]
    C = [x for x in range(1, L) if isp(2*x-3) and isp(2*x+3)]
    def sums(X, Y, top):
        s = bytearray(top + 1)
        for x in X:
            if x > top: break
            for y in Y:
                if x + y > top: break
                s[x + y] = 1
        return s
    top = 5000
    sAB = sums(A, B, top); sBC = sums(B, C, top); sAC = sums(A, C, top)
    s2A = sums(A, A, top); s2B = sums(B, B, top); s2C = sums(C, C, top)
    print('2.6 A+B missing in [2,%d]:' % top, [n for n in range(2, top + 1) if not sAB[n]])
    print('2.6 B+C missing in [5,%d]:' % top, [n for n in range(5, top + 1) if not sBC[n]])
    print('2.6 A+C missing in [5,%d]:' % top, [n for n in range(5, top + 1) if not sAC[n]], '(claim only 161)')
    print('2.6 2A last missing <=%d:' % top, [n for n in range(2, top + 1) if not s2A[n]][-3:], '(claim >=702 all)')
    print('2.6 2B last missing:', [n for n in range(2, top + 1) if not s2B[n]][-3:], '(claim >=492 all)')
    print('2.6 2C last missing:', [n for n in range(2, top + 1) if not s2C[n]][-3:], '(claim >=4006 all)')
    # 3.12(i): 1<m<n, m does not divide n: floor(kn/m) prime for some k in 1..n-1
    bad = [(m, n) for n in range(3, 400) for m in range(2, n) if n % m and not any(isp(k*n//m) for k in range(1, n))]
    print('3.12(i) failures n<400:', bad[:10])
    # 3.13(i): S={n: 2n+1, 2n^3+1 prime}; n>2 sum of three elements of S
    S = [x for x in range(1, 3000) if isp(2*x+1) and isp(2*x**3+1)]
    s2 = bytearray(6001)
    for a in S:
        for b in S:
            if a + b <= 6000: s2[a+b] = 1
    bad = [n for n in range(3, 2001) if not any(s2[n - c] for c in S if c < n)]
    print('3.13(i) failures n in [3,2000]:', bad)
    # 3.14: P={p: p, p+6, 3p+8 prime}; n>6: 2n+1=p+q+r with p,q,r in P and p+q+9 prime
    P = [p for p in PL[:5000] if isp(p+6) and isp(3*p+8)]
    Pset = set(P)
    bad = []
    for n in range(7, 1500):
        t = 2*n + 1; ok = False
        for p in P:
            if p >= t: break
            for q in P:
                if p + q >= t: break
                if (t - p - q) in Pset and isp(p + q + 9): ok = True; break
            if ok: break
        if not ok: bad.append(n)
    print('3.14 failures n in [7,1500):', bad)
    # 3.15: P'={p: p, 3p-4, 3p-10, 3p-14 prime}; n>17: 2n = p+q+r+s
    P2 = [p for p in PL[:5000] if isp(3*p-4) and isp(3*p-10) and isp(3*p-14)]
    s2 = bytearray(8001)
    for a in P2:
        for b in P2:
            if a + b <= 8000: s2[a+b] = 1
    pairs = [x for x in range(8001) if s2[x]]
    bad = [n for n in range(18, 2001) if not any(s2[2*n - x] for x in pairs if x < 2*n)]
    print('3.15 failures n in [18,2000]:', bad)
    # 3.16(ii): each n = x^2+y^2+p(p+e)/2
    vals = sorted(set(p*(p+e)//2 for p in PL[:2000] for e in (1, -1)))
    sq2 = bytearray(20001)
    for x in range(0, 142):
        for y in range(x, 142):
            if x*x + y*y <= 20000: sq2[x*x + y*y] = 1
    bad = [n for n in range(1, 20001) if not any(sq2[n - v] for v in vals if v <= n)]
    print('3.16(ii) failures n<=20000:', bad)
    # 3.16(i): n>1: x^2+y^2+phi(z^2), x<=y, y or z prime
    phis = {}
    for z in range(1, 200):
        phis.setdefault(z * PHI[z], []).append(z)
    bad = []
    for n in range(2, 4001):
        ok = False
        for v, zs in phis.items():
            if v > n: continue
            rest = n - v
            for x in range(0, isqrt(rest // 2) + 1):
                y2 = rest - x*x
                y = isqrt(y2)
                if y*y == y2 and x <= y and (isp(y) or any(isp(z) for z in zs)): ok = True; break
            if ok: break
        if not ok: bad.append(n)
    print('3.16(i) failures n in [2,4000]:', bad)
    # 4.17(i) [#87]: p>5: k^2<p with inverse of k^2 mod p prime
    bad = [p for p in PL if 5 < p < 200000 and not any(isp(pow(k*k, -1, p)) for k in range(1, isqrt(p - 1) + 1))]
    print('4.17(i) failures p<2e5:', bad)
    # 4.18(i) [#88]: n>2, n!=8: n=k+m, k!=m, p_k prim root mod p_m and vice versa
    bad = [n for n in range(3, 1500) if n != 8 and not any(k != n-k and is_pr_root(pk[k], pk[n-k]) and is_pr_root(pk[n-k], pk[k]) for k in range(1, n))]
    print('4.18(i) failures n in [3,1500):', bad)

names = list(R) if which == 'all' else [w for w in R if w.startswith(which)]
for name in names:
    f, start, exc, lim = R[name]
    top = min(N, lim) if lim else N
    fails = [n for n in range(1, top + 1) if not f(n)]
    bad = [n for n in fails if n >= start and n not in exc]
    print(name, 'n<=%d' % top, 'COUNTEREXAMPLES:' if bad else 'ok', bad[:10], '| failures below threshold/exc:', [n for n in fails if n < start or n in exc][-6:], flush=True)
