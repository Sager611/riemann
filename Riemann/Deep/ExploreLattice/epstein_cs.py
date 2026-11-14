# epstein_cs.py (explore48): Epstein zeta of the unimodular 2D lattice L_z, z = x + iy (y >= sqrt(3)/2), via Chowla-Selberg.
# E*(z,s) = pi^-s Gamma(s) (1/2) sum' y^s/|mz+n|^{2s}
#         = zeta*(2s) y^s + zeta*(2s-1) y^{1-s} + 4 sqrt(y) sum_{k>=1} k^{s-1/2} sigma_{1-2s}(k) K_{s-1/2}(2 pi k y) cos(2 pi k x),
# zeta*(w) = pi^{-w/2} Gamma(w/2) zeta(w).  For Q = a m^2 + b mn + c n^2, disc D: Lambda_Q(s) = (sqrt|D|/2pi)^s Gamma(s) Z_Q(s) = 2 E*(z,s),
# z = (-b + i sqrt|D|)/(2a).  The K-Bessel terms are the dual (Poisson) frequencies of the line sums along the shortest vector:
# the lattice-crypto 'smoothing' tail.  On the line s = 1/2 + it:
#   F(t) := E*(z,1/2+it) / (2 sqrt(y) |zeta*(1+2it)|) = cos(theta_a(t)) + beta(t),
#   theta_a = arg(zeta*(1+2it) y^{it}),  beta = 2 sqrt(cosh(pi t)) sum_k C_k(t) K_{it}(2 pi k y) cos(2 pi k x) / |zeta(1+2it)|,
#   C_k(t) = sum_{d|k} cos(t log(k/d^2))  (= k^{it} sigma_{-2it}(k), real).
import math, cmath
from flint import acb, arb, ctx

def divisors(k):
    ds = []
    i = 1
    while i * i <= k:
        if k % i == 0:
            ds.append(i)
            if i * i != k: ds.append(k // i)
        i += 1
    return ds

_DIV = {}
def C_k(k, t):
    ds = _DIV.get(k)
    if ds is None:
        ds = divisors(k); _DIV[k] = ds
    return sum(math.cos(t * math.log(k / (d * d))) for d in ds)

def kmax_needed(t, y, eps_exp=40.0):
    # K_{it}(2 pi k y) negligible (relative e^{-eps_exp}) once 2 pi k y * g(t/(2 pi k y)) > eps_exp; g(a)=sqrt(1-a^2)-a*acos(a)
    k = 1
    while True:
        u = 2 * math.pi * k * y
        if u > t:
            a = t / u
            g = math.sqrt(1 - a * a) - a * math.acos(a)
            if u * g > eps_exp: return k - 1
        k += 1
        if k > 100000: return k

def Kscaled(t, u, tol=1e-13):
    """e^{pi t/2} K_{it}(u) as float (t >= 0, u > 0); precision doubled until the arb radius < tol (absolute)."""
    old = ctx.prec
    p = int(2.3 * abs(t)) + 80
    for _ in range(10):
        ctx.prec = p
        v = acb(u).bessel_k(acb(0, t)) * (arb.pi() * abs(t) / 2).exp()
        if float(v.real.rad()) < tol: break
        p *= 2
    else:
        ctx.prec = old
        raise ArithmeticError('Kscaled: no convergence t=%g u=%g' % (t, u))
    r = float(v.real.mid())
    ctx.prec = old
    return r

def theta_a(t, y):
    # arg(zeta*(1+2it) y^{it}) continuous in t (using loggamma), + arg zeta(1+2it) principal branch (|arg zeta(1+2it)| < pi/2 roughly; fine for sign tests)
    s = acb(0.5, t)
    lg = s.lgamma()
    z1 = acb(1, 2 * t).zeta()
    th = float(lg.imag.mid()) - t * math.log(math.pi) + t * math.log(y) + float(z1.arg().mid())
    return th, abs(complex(float(z1.real.mid()), float(z1.imag.mid())))

def F_line(t, x, y, eps_exp=40.0, return_parts=False):
    th, az = theta_a(t, y)
    km = kmax_needed(t, y, eps_exp)
    S = 0.0
    for k in range(1, km + 1):
        ck = C_k(k, t)
        cx = math.cos(2 * math.pi * k * x)
        if ck == 0.0 or abs(cx) < 1e-15: continue
        S += ck * Kscaled(t, 2 * math.pi * k * y) * cx
    # 2 sqrt(cosh(pi t)) e^{-pi t/2} = sqrt(2) sqrt(1 + e^{-2 pi t})
    beta = 2 * math.sqrt((1 + math.exp(-2 * math.pi * t)) / 2) * S / az
    if return_parts: return math.cos(th) + beta, math.cos(th), beta, th
    return math.cos(th) + beta

def Estar(s, x, y, kextra=8, prec=None):
    """E*(z,s) as acb, complex s; normalized output G = e^{-i pi s/2} E* (magnitude O(poly)) returned as complex."""
    s = acb(s)
    t = abs(float(s.imag.mid()))
    old = ctx.prec
    ctx.prec = prec or (int(2.3 * t) + 100)
    s = acb(s.real.mid().str(30), s.imag.mid().str(30)) if False else s
    yy = arb(repr(y)) if isinstance(y, float) else arb(y)
    xx = arb(repr(x)) if isinstance(x, float) else arb(x)
    pi = arb.pi()
    def zstar(w):
        return (pi ** (-w / 2)) * (w / 2).gamma() * w.zeta()
    val = zstar(2 * s) * yy ** s + zstar(2 * s - 1) * yy ** (1 - s)
    km = kmax_needed(t, float(y), 60.0) + kextra
    nu = s - acb(0.5)
    tot = acb(0)
    for k in range(1, km + 1):
        sig = acb(0)
        for d in divisors(k):
            sig += acb(d) ** (1 - 2 * s)
        cx = (2 * pi * k * xx).cos()
        tot += acb(k) ** nu * sig * acb(2 * pi * k * yy).bessel_k(nu) * cx
    val += 4 * yy.sqrt() * tot
    G = (acb(0, -1) * pi * s / 2).exp() * val
    ctx.prec = old
    rad = max(float(G.real.rad()), float(G.imag.rad()))
    mag = abs(complex(float(G.real.mid()), float(G.imag.mid())))
    if rad > 1e-12 * max(1.0, mag):
        if (prec or 0) > 40000: raise ArithmeticError('Estar: no convergence at s=%s' % s)
        return Estar(s, x, y, kextra, (prec or (int(2.3 * t) + 100)) * 2)
    return complex(float(G.real.mid()), float(G.imag.mid()))
