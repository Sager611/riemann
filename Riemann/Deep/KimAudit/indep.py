#!/usr/bin/env python3
"""kimaudit T1-T3: independent mpmath recomputation from DERIVATION formulas (not from certify.py).
Run: bash .scratch/memcap.sh 290 900 uv run --no-project --with mpmath python3 indep.py"""
import mpmath as mp
mp.mp.dps = 60
pi, mpf = mp.pi, mp.mpf
A = {'a': mpf('0.37338196661320752'), 'half': mpf(1)/2}
KAP = {'a': mpf('6.45'), 'half': mpf('3.59')}
D_all = mpf(16140857488599731839205676)
def Z(u): return u*(pi*mp.exp(2*u) - mpf(9)/4)
def uinv(n):
    lo, hi = mpf('1e-6'), mpf(1000)
    for _ in range(260):
        mid = (lo+hi)/2
        if Z(mid) < n: lo = mid
        else: hi = mid
    return (lo+hi)/2
def Rn(n): return mpf(9)/10*(n-mpf(1)/2)**2/(uinv(n)+1)**2
def Dn(n, a): R = Rn(n); return R**2/(4*a**2*(R+a**2))
def bis(f, lo, hi, it=200):   # f(lo) false, f(hi) true, f monotone
    for _ in range(it):
        mid = (lo+hi)/2
        if f(mid): hi = mid
        else: lo = mid
    return lo, hi
print('T1  Z(12) =', mp.nstr(Z(mpf(12)), 30), ' n_* = ceil =', int(mp.ceil(Z(mpf(12)))))
ns = int(mp.ceil(Z(mpf(12)))); assert ns == 998616377771
print('    u(n_*) =', mp.nstr(uinv(mpf(ns)), 20), '(>= 12 needed)')
a = A['a']
for n in [mpf(10)**13, mpf('4.6974172e13'), mpf('4.6974173e13')]:
    print('T2  n = %s: R_n = %s  D_n = %s  D_n - D_all = %s' % (mp.nstr(n, 9), mp.nstr(Rn(n), 15), mp.nstr(Dn(n, a), 20), mp.nstr(Dn(n, a)-D_all, 6)))
lo, hi = bis(lambda n: Dn(n, a) > D_all, mpf(10)**13, mpf(10)**14, 120)
print('    real crossover D_n = D_all at n =', mp.nstr(hi, 15), ' integer threshold', int(mp.ceil(hi)))
lo, hi = bis(lambda n: Dn(n, a) >= D_all+1, mpf(10)**13, mpf(10)**14, 120)
print('    smallest n with D_n >= D_all+1:', int(mp.ceil(hi)))
lo3, hi3 = bis(lambda n: n**2*mp.log(n+2) >= 24*(D_all+1)**3, mpf(10)**30, mpf(10)**45, 300)
print('    two-thirds wedge at d = D_all+1 needs n >=', mp.nstr(hi3, 8))
n30 = mpf(10)**30
print('    n = 1e30: D_n =', mp.nstr(Dn(n30, a), 8), ' two-thirds d <=', mp.nstr((n30**2*mp.log(n30+2)/24)**(mpf(1)/3), 8))
print('    d0 = (n_*^2 log(n_*+2)/24)^(1/3) =', mp.nstr((mpf(ns)**2*mp.log(ns+2)/24)**(mpf(1)/3), 10))
for key in ['a', 'half']:
    a = A[key]; kap = KAP[key]
    print('T3  a =', mp.nstr(a, 20), ' 3.6/(4a^2) =', mp.nstr(mpf('3.6')/(4*a**2), 12), ' stated', kap)
    worst = mpf(10)**9; worstu = mpf(-1)
    for k in range(0, 2001):
        n = mpf(ns)*mp.power(10, mpf(k)/7)
        r = Dn(n, a)/(kap*n**2/mp.log(n)**2)
        worst = min(worst, r)
        worstu = max(worstu, uinv(n)+1-mp.log(n)/2)
    print('    min_{n in [n_*, 1e297]} D_n/(kap n^2/log^2 n) =', mp.nstr(worst, 10), ' max (u+1 - log(n)/2) =', mp.nstr(worstu, 8))
    for c in ['0.51', '0.55', '0.6', '1']:
        c = mpf(c)
        f = lambda L: mp.exp((1-2*c)*L)*(L**2/kap)**c
        Ls = mp.log(ns); Lm = 2*c/(2*c-1)
        best = max(f(Ls + mpf(j)/100) for j in range(0, 40000))
        print('    c = %s: sup f (grid) = %s, closed form f(max(L*, 2c/(2c-1))) = %s' % (mp.nstr(c, 3), mp.nstr(best, 8), mp.nstr(f(max(Ls, Lm)), 8)))
    Cpage = {'a': ['7.6890', '0.87160', '0.069829', '1.1853e-10'], 'half': ['10.367', '1.2031', '0.099247', '2.1294e-10']}[key]
    bad = 0; mn = mpf(10)**9
    for c, C in zip(['0.51', '0.55', '0.6', '1'], Cpage):
        c, C = mpf(c), mpf(C)
        for j in range(0, 601):
            d = mp.power(10, mpf(j)/10)
            n = mp.ceil(max(mpf(ns), C*d**c))
            q = Dn(n, a)/d
            mn = min(mn, q)
            if q < 1: bad += 1
    print('    item-3 sampling (2404 (c,d) points, 1 <= d <= 1e60): failures', bad, ' min D_n/d =', mp.nstr(mn, 8))
