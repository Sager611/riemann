# indep.py: independent mpmath recomputation (isoaudit). Formulas re-derived from DERIVATION.md, not from cert.py.
# Run: OMP_NUM_THREADS=1 uv run --no-project --with mpmath python3 indep.py
import mpmath as mp, sys, time
mp.mp.dps = 30
t0 = time.time()
PI = mp.pi
TH = mp.mpf(4)/7
def kern(w):
    Z = 2*mp.sin(w/2)/w
    f = lambda x: mp.cos(w*x)/Z
    def Kr(u):   # K(u) = 2 int_0^{1/2} f cos(2 pi u x), closed form (re-derived)
        a, b = (w+2*PI*u)/2, (w-2*PI*u)/2
        s = lambda z: mp.sin(z)/z if z != 0 else mp.mpf(1)
        return w*(s(a)+s(b))/(4*mp.sin(w/2))
    def Ki(y):   # K(iy) = 2 int_0^{1/2} f cosh(2 pi y x), closed form (re-derived)
        return (w/mp.sin(w/2))*(w*mp.sin(w/2)*mp.cosh(PI*y) + 2*PI*y*mp.cos(w/2)*mp.sinh(PI*y))/(w**2+4*PI**2*y**2)
    return f, Kr, Ki
def Cq(f):   # C(f) = int f^2 + iint |u-v| f f, by quadrature (definition, no closed form)
    I1 = mp.quad(lambda x: f(x)**2, [-0.5, 0.5])
    G = lambda u: mp.quad(lambda v: abs(u-v)*f(v), [-0.5, u, 0.5])
    I2 = mp.quad(lambda u: f(u)*G(u), [-0.5, 0, 0.5])
    return I1 + I2
def cstar(R, th=TH): return 1 + (1-mp.exp(-2*R))/(mp.exp(2*th*R)-1)
def Bh(rho): return mp.log(cstar(2*PI*rho))/(4*PI)
def lamC(Kr, c):
    k0, k1 = Kr(2*c), Kr(4*c)
    return max(2*(1-k0), (6-16*k0/3-8*k1/3)/3, 0)
def Hfun(Ki, c):
    thb = PI**2*c**2/2
    return lambda y: 1 - (1-thb)**2*4*Ki(y)**2/(1+Ki(2*y))
def first_root(F, a=mp.mpf('0.05'), b=mp.mpf(6), n=600):
    ys = [a + (b-a)*k/n for k in range(n+1)]
    for y0, y1 in zip(ys, ys[1:]):
        if F(y0) < 0 <= F(y1): return mp.findroot(F, (y0, y1), solver='illinois')
    return None
def supratio(F, rho, ystart, yend=40, step=mp.mpf('0.005')):
    best, yb = -mp.inf, None; y = ystart
    while y < yend:
        v = F(y)/(y-rho)
        if v > best: best, yb = v, y
        y += step
    g = lambda t: -F(t)/(t-rho)
    a, b = yb-step, yb+step; gr = (mp.sqrt(5)-1)/2
    for _ in range(80):
        c1, c2 = b-gr*(b-a), a+gr*(b-a)
        if g(c1) < g(c2): b = c2
        else: a = c1
    ym = (a+b)/2
    return F(ym)/(ym-rho), ym
def p_hinge(c, w=mp.sqrt(2), back=mp.mpf('0.02'), C=None, verbose=False):
    f, Kr, Ki = kern(w)
    if C is None: C = Cq(f)
    lam = lamC(Kr, c); H = Hfun(Ki, c); F = lambda y: H(y) - 3*lam
    rmax = first_root(F)
    if rmax is None: rmax = mp.mpf(6)
    rho = rmax - back
    a, ym = supratio(F, rho, rmax)
    B = Bh(rho)
    p = (2 - lam - C - a*B)/(1-lam)
    if verbose: print('  c=%s lam=%s rho_max=%s rho=%s alpha=%s (y*=%s) B=%s p=%s' % tuple(mp.nstr(v, 12) for v in (c, lam, rmax, rho, a, ym, B, p)))
    return p, dict(lam=lam, rmax=rmax, rho=rho, alpha=a, B=B, C=C)
