# indep.py (wsgsaudit): independent mpmath checks of WeilSimple's matrix entries, mu(k) and inner-product closed forms.
# Usage: indep.py c   (integer c = lambda^2). Imports the audited code from a /tmp copy only to compare against it.
import sys, time, json
sys.path.insert(0, '/tmp/wsgs_copy')
import mpmath as mp
from flint import arb, ctx, fmpq
ctx.prec = 200
from weil_matrix import WeilMat
import cert_simple as CS
mp.mp.dps = 40
cn = int(sys.argv[1]); c = mp.mpf(cn); L = mp.log(c); H = L/2
def vM(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n
            while m % p == 0: m //= p
            return mp.log(p) if m == 1 else mp.mpf(0)
atoms = [(mp.log(n), vM(n)/mp.sqrt(n)) for n in range(2, cn) if vM(n) != 0]   # n < c (n = c sits at C(L) = 0)
def bas(j, sec):
    w = 2*mp.pi*j/L
    if sec == 'even':
        N = 1/mp.sqrt(L) if j == 0 else mp.sqrt(2/L)
        return lambda x: N*mp.cos(w*x)
    return lambda x: mp.sqrt(2/L)*mp.sin(w*x)
def ip(f, g): return mp.quad(lambda x: f(x)*g(x), mp.linspace(-H, H, 9))
import numpy as np
GX, GW = np.polynomial.legendre.leggauss(80)
def gl(fun, a, b):   # 80-point Gauss-Legendre on [a, b], double precision (vectorized fun)
    x = (b - a)/2*GX + (a + b)/2; return (b - a)/2*np.dot(GW, fun(x))
Lf = float(L); Hf = Lf/2
def nb(j, sec):
    w = 2*np.pi*j/Lf
    if sec == 'even': N = 1/np.sqrt(Lf) if j == 0 else np.sqrt(2/Lf); return lambda x: N*np.cos(w*x)
    return lambda x: np.sqrt(2/Lf)*np.sin(w*x)
def ncorr(f, g, x):   # int f(y) g(y - x) dy, 0 <= x <= L, by quadrature (no closed form)
    return gl(lambda y: f(y)*g(y - x), x - Hf, Hf) if x < Lf else 0.0
def Qn(f, g, sec):   # CCM (3.10) in x = log u: pole + archimedean (Gauss integral of psi) + primes
    v = (lambda x: np.cosh(x/2)) if sec == 'even' else (lambda x: np.sinh(x/2))
    ipn = lambda a, b: sum(gl(lambda x: a(x)*b(x), lo, lo + Lf/4) for lo in (-Hf, -Hf/2, 0, Hf/2))
    pole = (2 if sec == 'even' else -2)*ipn(f, v)*ipn(g, v); fg = ipn(f, g)
    kern = lambda xs: np.array([(2*np.exp(-x)*fg - 2*np.exp(x/2)*ncorr(f, g, x))/(np.exp(x) - np.exp(-x)) for x in xs])
    arch = -(float(mp.log(mp.pi)) + float(mp.euler))*fg + sum(gl(kern, lo, lo + Lf/4) for lo in (0, Lf/4, Lf/2, 3*Lf/4))
    arch += 2*fg*float(-mp.log(1 - c**-2)/2)          # int_L^inf 2 e^-x fg/(e^x - e^-x) dx
    return pole + arch - 2*sum(float(w)*ncorr(f, g, float(ln)) for ln, w in atoms)
out = {}
t0 = time.time()
W = WeilMat(fmpq(cn))
# (A) entries by x-space quadrature (double precision) vs weil_matrix closed forms (Arb)
mx = 0
for sec, idx in (('even', [0, 1, 2, 5, 9]), ('odd', [1, 2, 5, 9])):
    T, I = W.block(max(idx), sec)
    for a in range(len(idx)):
        for b in range(a, len(idx)):
            i, j = idx[a], idx[b]
            q = Qn(nb(i, sec), nb(j, sec), sec)
            t = float(T[I.index(i), I.index(j)].mid().str(20, radius=False))
            mx = max(mx, abs(q - t))
            print(f"A {sec} ({i},{j}) quad {q:.15e} closed {t:.15e} diff {abs(q - t):.2e}", flush=True)
out['A_maxdiff'] = mp.nstr(mx, 3)
print('A done', round(time.time() - t0), 's', flush=True)
# (C) lattice identity with my own mu(k) vs weil_matrix on f = cos^2(pi x/L) (even) and sin(a x) + sin(2 a x)/2 (odd)
def mu(k):
    r = mp.pi*k/L
    v = -mp.log(mp.pi) + mp.re(mp.digamma(mp.mpf(1)/4 + 0.5j*r))
    v += 2*(-1)**k*sum(c**(-(2*n + mp.mpf(1)/2))*(2*n + mp.mpf(1)/2)/((2*n + mp.mpf(1)/2)**2 + r*r) for n in range(int(60/float(mp.log10(c))) + 3))
    return v - 2*sum(w*mp.cos(r*ln) for ln, w in atoms)
S = CS.SymZeta(fmpq(cn))
md = max(abs(mu(k) - mp.mpf(S.mu(k).mid().str(45, radius=False))) for k in (0, 1, 2, 3, 7, 40, 401, 12345))
out['mu_maxdiff'] = mp.nstr(md, 3); print('mu vs SymZeta.mu maxdiff', out['mu_maxdiff'], flush=True)
mp.mp.dps = 25
a1 = 2*mp.pi/L; K = 4001
lat_e = mu(0)*(L/2)**2 + 2*mu(2)*(L/4)**2
lat_o = 2*mu(2)*(L/2)**2 + 2*mu(4)*(L/4)**2
for k in range(1, K + 1, 2):
    r = mp.pi*k/L; s = (-1)**((k - 1)//2); m = mu(k)
    lat_e += 2*m*(s*a1*a1/(r*(r*r - a1*a1)))**2
    lat_o += 2*m*(6*s*a1**3/((r*r - a1*a1)*(r*r - 4*a1*a1)))**2
fe = lambda x: mp.cos(mp.pi*x/L)**2; fo = lambda x: mp.sin(a1*x) + mp.sin(2*a1*x)/2
lat_e = 2*ip(fe, lambda x: mp.cosh(x/2))**2 + lat_e/(2*L)
lat_o = -2*ip(fo, lambda x: mp.sinh(x/2))**2 + lat_o/(2*L)
Te, _ = W.block(1, 'even'); To, _ = W.block(2, 'odd')
ge = [mp.sqrt(L)/2, mp.sqrt(L/2)/2]; go = [mp.sqrt(L/2), mp.sqrt(L/2)/2]
m_ = lambda T, i, j: mp.mpf(T[i, j].mid().str(30, radius=False))
qe = sum(ge[i]*ge[j]*m_(Te, i, j) for i in range(2) for j in range(2))
qo = sum(go[i]*go[j]*m_(To, i, j) for i in range(2) for j in range(2))
out['LAT_even'] = [mp.nstr(lat_e, 16), mp.nstr(qe, 16)]; out['LAT_odd'] = [mp.nstr(lat_o, 16), mp.nstr(qo, 16)]
print('LAT even (lattice, matrix)', out['LAT_even'], 'odd', out['LAT_odd'], flush=True)
# (D) cert_simple inner-product closed forms vs quadrature (formulas copied from cert_simple.py l. 301-311)
mp.mp.dps = 30
dmax = 0
for sec in ('even', 'odd'):
    Nj = lambda j: 1/mp.sqrt(L) if (j == 0 and sec == 'even') else mp.sqrt(2/L)
    sgn = lambda j: 1 if j % 2 == 0 else -1; sgk = lambda k: 1 if ((k - 1)//2) % 2 == 0 else -1
    rk = lambda k: mp.pi*k/L; yv = -(L*L)/(4*mp.pi**2)
    if sec == 'even':
        A = lambda j: Nj(j)*(L/mp.pi)*sgn(j); bk = lambda k: 2*k*sgk(k); bv = -(L/mp.pi)*mp.sinh(L/4)
        inv = lambda k: 2*rk(k)*sgk(k)*mp.cosh(L/4)/(rk(k)**2 + 0.25); vv = mp.sinh(L/2) + L/2
        ph = lambda j: (lambda x: Nj(j)*mp.cos(2*mp.pi*j*x/L)); ek = lambda k: (lambda x: mp.cos(mp.pi*k*x/L)); v = lambda x: mp.cosh(x/2)
    else:
        A = lambda j: Nj(j)*sgn(j)*4*j*(L/mp.pi); bk = lambda k: sgk(k); bv = mp.sinh(L/4)
        inv = lambda k: sgk(k)*mp.cosh(L/4)/(rk(k)**2 + 0.25); vv = mp.sinh(L/2) - L/2
        ph = lambda j: (lambda x: Nj(j)*mp.sin(2*mp.pi*j*x/L)); ek = lambda k: (lambda x: mp.sin(mp.pi*k*x/L)); v = lambda x: mp.sinh(x/2)
    for j in ([0, 1, 3] if sec == 'even' else [1, 3]):
        for k in (1, 3, 7):
            dmax = max(dmax, abs(ip(ph(j), ek(k)) - A(j)*bk(k)/(k*k - 4*j*j)))
        dmax = max(dmax, abs(ip(ph(j), v) - A(j)*bv/(yv - 4*j*j)))
    for k in (1, 3, 7): dmax = max(dmax, abs(ip(ek(k), v) - inv(k)))
    dmax = max(dmax, abs(ip(v, v) - vv))
out['D_maxdiff'] = mp.nstr(dmax, 3); out['sec'] = round(time.time() - t0)
print(json.dumps(dict(c=cn, **out)), flush=True)
