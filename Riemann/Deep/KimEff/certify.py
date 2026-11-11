#!/usr/bin/env python3
"""kimeff certificate: explicit Kim-Lee-shape hyperbolicity region for the Jensen polynomials of Riemann's Xi.
Run: uv run --no-project --with python-flint python3 certify.py [a]   (a = 0.37338196661320752 default, or 0.5)
Every inequality of DERIVATION.md Sec. 8-10 is evaluated in Arb ball arithmetic; any failure aborts (exit 1).
Strict inequalities are checked as 'x > 0' on balls (the whole ball is positive)."""
import sys
# python-flint is supplied at run time by uv (--with python-flint); it is not installed system-wide (no pip).
from flint import arb, ctx  # pyright: ignore[reportMissingImports]
ctx.prec = 256
A_STR = sys.argv[1] if len(sys.argv) > 1 else '0.37338196661320752'
a = arb(A_STR)
pi = arb.pi()
def ex(x): return arb(x).exp()
def lg(x): return arb(x).log()
def sq(x): return arb(x).sqrt()
FAIL = []
def check(name, cond_ball, show):
    ok = bool(cond_ball > 0)
    print(('PASS ' if ok else 'FAIL ') + name + ': ' + show)
    if not ok: FAIL.append(name)
from decimal import Decimal, Context, ROUND_CEILING, ROUND_FLOOR
def up(x, d=8):
    """decimal string with d significant digits that is >= every point of the ball x (directed rounding up)"""
    v = x.upper(); c = Context(prec=d, rounding=ROUND_CEILING)
    r = c.plus(Decimal(v.str(d + 12, radius=False)))
    while not bool(arb(str(r)) >= v): r = c.next_plus(r)
    return str(r)
def lo(x, d=8):
    """decimal string with d significant digits that is <= every point of the ball x (directed rounding down)"""
    v = x.lower(); c = Context(prec=d, rounding=ROUND_FLOOR)
    r = c.plus(Decimal(v.str(d + 12, radius=False)))
    while not bool(arb(str(r)) <= v): r = c.next_minus(r)
    return str(r)

# ---- fixed parameters (DERIVATION Sec. 7) ----
theta = arb(9)/10; thetaU = arb(49)/50; dR = arb(1); Delta = arb(17)/25; K = arb(4); ee = arb(1)/1000
u_star = arb(12)
beta = 3/(2*pi)
rho = sq(thetaU/theta)
Chat = theta*(2-theta)/(4*(1-theta)**2)
print('a =', A_STR, '| theta = 9/10, thetaU = 49/50, delta_R = 1, Delta = 17/25, K = 4, eps_eta = 1/1000, u_* = 12')

# ---- eps0 = sup_t phi_{>=2}/phi_1 = value at t = 0 (Sec. 5) ----
eps0 = arb(0)
for k in range(2, 12):
    k = arb(k)
    eps0 += k**2*(2*pi*k**2-3)/(2*pi-3)*ex(-pi*(k**2-1))
eps0 += ex(-300)   # tail k >= 12: sum k^2(2 pi k^2-3)/(2pi-3) e^{-pi(k^2-1)} < e^{-300}
print('eps0 <=', up(eps0))

# ---- Q2 (n-independent) ----
Q2 = (1-ex(-2*Delta)) - (theta/4 + theta**2/(16*(1-theta)))
check('Q2 left bracket (n-independent): 1-e^{-2Delta} - theta/4 - theta^2/(16(1-theta))', Q2, '>= ' + lo(Q2))
check('rho - 1 >= 1/(2 u_*) (F increasing on u >= u_*)', (rho-1) - 1/(2*u_star), 'rho >= ' + lo(rho, 10))

def Zf(u): return u*(pi*ex(2*u) - arb(9)/4)

def majorants(u):
    """Monotone majorants of DERIVATION Sec. 9 at a given u (n = Z(u))."""
    n = Zf(u)
    out = {'n': n}
    out['M1'] = 4 - 2*beta*u*ex(-2*u)/(1-beta) - u*(u+1)*a**2/(2*n)
    klow = 4*(n/u)*ex(-2*(Delta+2*ee)) - a**2/(2*n)
    eta_up = sq(K/klow)
    out['eta_up'] = eta_up
    out['M3a'] = ee - eta_up
    out['M3b'] = dR - 1/n - 2*eta_up
    out['Phi_up'] = 2*a*(1+Chat)*eta_up*(sq(theta) + a*(u+1)/n)
    out['E_up'] = a**2*rho**4*(theta*(u+1)**2/n + a**2*(u+1)**4/n**3)/(8*(1-thetaU)**2)
    y = (rho-1)*u + rho*dR                       # U - u
    F = (ex(2*y) - 1 - 2*y)/u
    G1 = thetaU/4 + a*rho**2*(u+1)*(sq(theta) + a*(u+1)/n)/(2*n)
    eps = 1/n + 2*eta_up
    uL = u - Delta
    dh1 = Delta**2/(uL*u) + (2/u)*(1+9*u/(4*n))*(Delta-(1-ex(-2*Delta))/2)
    dh2 = eps**2/u**2 + (1/u)*(1+9*u/(4*n))*(ex(2*eps)-1-2*eps)
    dh = dh1 + dh2   # >= max(dh1, dh2), both >= 0
    mLG = theta/4 + theta**2/(32*(1-theta))
    out['F'] = F; out['G'] = G1 + dh + mLG
    c6 = F - (G1 + dh + mLG)
    out['c6'] = c6
    out['Hcp_low'] = (n/u)*(2*(ex(2*rho)-1) - rho*theta/2) - a*rho*(sq(theta)+a)
    logratio = lg(1+eps0) - lg(1-beta*ex(-2*uL)) + lg(5*n/K)/2 - c6*n
    out['logratio'] = logratio
    out['final'] = 2*(out['Phi_up']).cos() - (1+eps0)*ex(out['E_up'])/(2*K) - ex(logratio)
    return out

m = majorants(u_star)
zf = m['n']
nstar_int = int(zf.floor().unique_fmpz()) + 1
print('n_* = ceil(Z(12)) =', nstar_int, ' (Z(12) in', zf.str(20), ')')
check('M1 right bracket Q1*u at u_*', m['M1'], '>= ' + lo(m['M1']))
check('M3a eta_up <= eps_eta', m['M3a'], 'eta_up <= ' + up(m['eta_up']))
check('M3b 1/n + 2 eta_up <= delta_R', m['M3b'], 'ok')
check('M4 phase Phi_up < pi/2', pi/2 - m['Phi_up'], 'Phi_up <= ' + up(m['Phi_up']))
check('M5 E_up < 1', 1 - m['E_up'], 'E_up <= ' + up(m['E_up']))
check('M6 c6 = F(u_*) - G_up(u_*) > 0', m['c6'], 'F >= ' + lo(m['F']) + ', G_up <= ' + up(m['G']) + ', c6 >= ' + lo(m['c6']))
check('M6 Hcp_low >= 1', m['Hcp_low'] - 1, 'Hcp_low >= ' + lo(m['Hcp_low'], 6))
check('M6 c6 > 1/(2 n_*) (logratio decreasing in n)', m['c6'] - 1/(2*m['n']), 'ok')
print('     logratio(u_*) <=', up(m['logratio'], 6))
# Hc concavity on [U, oo): Hc'' = -2n/t^2 - 4 pi e^{2t} + Wmax/(2b) < 0 since 4 pi e^{2t} >= 4n/u and
# Wmax/(2b) <= (R + vmax)/(2n) <= theta n/(2(u+1)^2) + a(sqrt(theta)/(u+1) + a/n);  ratio to 4n/u decreasing in u
n_ = m['n']; u_ = u_star
hcc = 4*n_/u_ - (theta*n_/(2*(u_+1)**2) + a*(sq(theta)/(u_+1) + a/n_))
check('Hc concave on [U,oo): 4n/u - Wmax/(2b) bound > 0 at u_* (then increasing)', hcc, '>= ' + lo(hcc, 6))
check('Phi_up <= pi/2 (cos decreasing on [0, pi/2], so cos(Phi_up) increases with u)', pi/2 - m['Phi_up'], 'ok')
check('FINAL 2cos(Phi_up) - (1+eps0)e^{E_up}/(2K) - e^{logratio} > 0', m['final'], '>= ' + lo(m['final']))

print('cross-check grid (majorants are monotone by Sec. 9; values shown, not needed for the proof):')
for uu in [12, 13, 14, 16, 20, 30, 50, 100, 300]:
    mm = majorants(arb(uu))
    print('   u=%d  n=%s  M1=%s eta_up=%s Phi_up=%s E_up=%s c6=%s final=%s' % (uu, mm['n'].str(4, radius=False), mm['M1'].str(5,radius=False), mm['eta_up'].str(3,radius=False), mm['Phi_up'].str(3,radius=False), mm['E_up'].str(3,radius=False), mm['c6'].str(4,radius=False), mm['final'].str(5,radius=False)))
    if not (mm['M1'] > 0 and mm['M3a'] > 0 and mm['c6'] > 0 and mm['final'] > 0): FAIL.append('grid u=%d' % uu)

# ---- Sec. 10: hyperbolic region ----
four_a2 = 4*a**2
kappa = (arb(36)/10)/four_a2           # D_n >= kappa (n-1/2)^2/log^2 n - 1/4
kap_dec = arb('6.45') if A_STR.startswith('0.3733') else arb('3.59')
print('4a^2 in', four_a2.str(20), ' 3.6/(4a^2) >=', lo(kappa, 8), ' stated constant', kap_dec.str(4, radius=False))
check('H3 constant: 3.6/(4a^2) > stated', kappa - kap_dec, 'ok')
N0 = arb(nstar_int)
H2 = (lg(N0)/2 - 1)*(pi/ex(2) - arb(9)/(4*N0)) - 1
check('H2 (log n/2 - 1)(pi/e^2 - 9/(4n)) >= 1 at n_* (increasing in n)', H2, '>= ' + lo(H2))
H3 = (kappa - kap_dec)*N0**2 - kappa*N0 - lg(N0)**2/4
check('H3 slack (kappa-c)n^2 - kappa n - log^2(n)/4 > 0 at n_* (increasing)', H3, 'ok')

D_all = arb('16140857488599731839205676')
def u_upper(nn):
    lo_, hi_ = arb(1), arb(400)
    for _ in range(300):
        mid = (lo_+hi_)/2
        if bool(Zf(mid) < nn): lo_ = mid
        elif bool(Zf(mid) > nn): hi_ = mid
        else: break
    return hi_      # Z(hi_) > nn  =>  u(nn) < hi_
def D_of_n(nn):
    uu = u_upper(nn)
    R = theta*(nn-arb(1)/2)**2/(uu+1)**2
    return R**2/(four_a2*(R + a**2)), R
for e in [12, 13, 14, 15, 20, 30, 40]:
    nn = arb(10)**e
    if bool(nn < N0): continue
    D, R = D_of_n(nn)
    print('   n = 1e%d: R_n >= %s, D_n >= %s, c n^2/log^2 n >= %s, two-thirds wedge d <= %s' % (e, lo(R, 5), lo(D, 5), lo(kap_dec*nn**2/lg(nn)**2, 5), up((nn**2*lg(nn+2)/24)**(arb(1)/3), 5)))
lo_n, hi_n = arb(10)**12, arb(10)**16
for _ in range(120):
    mid = ((lo_n+hi_n)/2).floor()
    if bool(D_of_n(mid)[0] > D_all): hi_n = mid
    else: lo_n = mid
print('crossover: D_n > D_all for all n >=', up(hi_n, 8), '(D_n increasing in n)')
d1 = D_all + 1
lo_n, hi_n = arb(10)**30, arb(10)**45
for _ in range(300):
    mid = (lo_n+hi_n)/2
    if bool(mid**2*lg(mid+2) >= 24*d1**3): hi_n = mid
    else: lo_n = mid
print('two-thirds wedge at d = D_all+1 needs n >=', lo(lo_n, 6))
lo2, hi2 = arb(10)**12, arb(10)**16
for _ in range(120):
    mid = ((lo2+hi2)/2).floor()
    if bool(D_of_n(mid)[0] >= d1): hi2 = mid
    else: lo2 = mid
print('ours at d = D_all+1: n >=', up(hi2, 6))
C23 = 24*kap_dec**3
print('containment of the two-thirds region for d > D_all: needs log(n+2) log^6(n) <= 24*kap^3 n^4, 24*kap^3 in', C23.str(12))
# for n >= 8: d/dlog n of log(n^4/(log(n+2) log^6 n)) = 4 - n/((n+2)log(n+2)) - 6/log n >= 4 - 1/log 10 - 6/log 8 > 0
check('ratio n^4/(log(n+2)log^6 n) increasing for n >= 8', 4 - 1/lg(arb(10)) - 6/lg(arb(8)), 'ok')
check('containment at n = 8 (hence all n >= 8)', C23*arb(8)**4 - lg(arb(10))*lg(arb(8))**6, 'ok')
# self-contained comparison (no all-n range): if d > d0 := (n_*^2 log(n_*+2)/24)^{1/3}, every (d,n) with n^2 log(n+2) >= 24 d^3
# has n > n_* (n^2 log(n+2) increasing), hence d <= 6.45 n^2/log^2 n by the containment inequality: our region contains it.
d0 = (N0**2*lg(N0+2)/24)**(arb(1)/3)
print('self-contained: for every d >', up(d0, 6), 'the two-thirds region at d lies inside ours (no all-n range used)')
# power forms: n >= C d^c and n >= n_*  =>  d <= (n/C)^{1/c} <= kap_dec n^2/log^2 n  iff  C >= f(n) := n^{1-2c} (log^2 n/kap_dec)^c;
# log f is concave in L = log n with maximum at L = 2c/(2c-1); so sup_{n >= n_*} f = f(n_*) if log n_* >= 2c/(2c-1), else f(e^{2c/(2c-1)}).
for c in ['0.51', '0.55', '0.6', '1']:
    cc = arb(c)
    Lm = 2*cc/(2*cc-1)
    L = lg(N0) if bool(lg(N0) >= Lm) else Lm
    Cc = ex((1-2*cc)*L)*(L**2/kap_dec)**cc
    print('   c = %s: J^{d,n} hyperbolic for every d >= 1 and every n >= max(n_*, C d^c) with C = %s (sup of f at log n = %s)' % (c, up(Cc, 5), L.str(5, radius=False)))
print()
print('ALL CHECKS PASS' if not FAIL else 'FAILURES: ' + ', '.join(FAIL))
sys.exit(1 if FAIL else 0)
