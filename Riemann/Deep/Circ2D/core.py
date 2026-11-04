# core.py -- independent mpmath reconstruction of the planar certificate
# (OpenAI preprint "A sharp Fourier certificate for planar circle packing").
# Written from the TeX formulas only (not from the authors' verification code).
import mpmath as mp
import re, os

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.environ.get('CIRC2D_SRC', os.path.join(HERE, '..', '..', '..', '.scratch', 'circ2d', 'src'))
A = [0, 1, 3, 4, 7, 9]

def setup(dps=40):
    mp.mp.dps = dps

def consts():
    b = mp.sqrt(3) / 2
    h = mp.mpf(2) / 5
    B = mp.mpf(4) / 3
    return b, h, B

def laurent_P():
    """Coefficients P_l, l=-6..6, of P(s)=prod (2 sin(pi(s-a)/12))^2 in w^l, w=e^{i pi s/6}."""
    rho = mp.expjpi(mp.mpf(1) / 6)
    poly = {0: mp.mpc(1)}
    for a in A:
        f = {0: mp.mpc(2), 1: -rho ** (-a), -1: -rho ** a}
        new = {}
        for k1, v1 in poly.items():
            for k2, v2 in f.items():
                new[k1 + k2] = new.get(k1 + k2, 0) + v1 * v2
        poly = new
    return {l: poly.get(l, mp.mpc(0)) for l in range(-6, 7)}

def nodes(nmax):
    return [n for n in range(1, int(nmax) + 1) if n % 12 in A]

def Pval(s):
    r = mp.mpf(1)
    for a in A:
        r *= (2 * mp.sin(mp.pi * (s - a) / 12)) ** 2
    return r

def P_deriv(s, k, PL):
    """k-th derivative of P from the Laurent form (real part)."""
    tot = mp.mpc(0)
    for l, Pl in PL.items():
        t = mp.mpf(l) / 6
        tot += Pl * (1j * mp.pi * t) ** k * mp.expjpi(t * s)
    return tot

def node_QD(n, PL):
    P2 = P_deriv(n, 2, PL).real
    P3 = P_deriv(n, 3, PL).real
    Q = P2 / 2
    D = P3 / (3 * P2)
    return Q, D

def read_table(path=None):
    """Return dict n -> (c1,d1,c2,d2) as mpf, from data/coefficients.tex (scaled 1e-10)."""
    path = path or os.path.join(SRC, 'build', 'data', 'coefficients.tex')
    tab = {}
    for line in open(path):
        m = re.findall(r'-?\d+', line)
        if len(m) == 5:
            n = int(m[0])
            tab[n] = tuple(mp.mpf(int(x)) / mp.mpf(10) ** 10 for x in m[1:])
    return tab

def lists_from_table(tab):
    """Two coefficient lists (dict n->(c,d)) incl. node 0, and constants C_i."""
    L1 = {0: (mp.mpf(1), mp.mpf(44) / 100)}
    L2 = {0: (mp.mpf(0), mp.mpf(-368) / 1000)}
    for n, (c1, d1, c2, d2) in tab.items():
        L1[n] = (c1, d1)
        L2[n] = (c2, d2)
    C1 = mp.mpf(-13) / 1000
    C2 = mp.mpf(17) / 1000
    return (L1, C1), (L2, C2)

def R_val(L, C, s):
    r = mp.mpf(C)
    for n, (c, d) in L.items():
        u = s - n
        r += c / u ** 2 + d / u
    return r

def p_direct(L, C, s):
    """p(s)=P(s)R(s) evaluated directly (s must not be a node with nonzero coefficient)."""
    # accurate near nodes: P(s)/(s-n)^2 via sine product with the vanishing factor isolated
    return Pval(s) * R_val(L, C, s)

class Measure:
    """Folded spectral measure of p = P R on [0,1]: Gauss-Legendre on the six pieces + atoms."""
    def __init__(self, L, C, PL, gl_n=96, split=1):
        self.L, self.C, self.PL = L, C, PL
        xs, ws = gl_nodes(gl_n)
        self.pts = []  # (t, w*f(t))
        for j in range(1, 7):
            al_be = {}
            for n in L:
                al = mp.mpc(0); be = mp.mpc(0)
                for l in range(j, 7):
                    tl = mp.mpf(l) / 6
                    e = PL[l] * mp.expjpi(tl * n)
                    al += e; be += tl * e
                al_be[n] = (al, be)
            for q in range(split):
                lo = mp.mpf(j - 1) / 6 + mp.mpf(q) / (6 * split)
                hi = lo + mp.mpf(1) / (6 * split)
                mid = (lo + hi) / 2; half = (hi - lo) / 2
                for x, w in zip(xs, ws):
                    t = mid + half * x
                    f = mp.mpc(0)
                    for n, (c, d) in L.items():
                        al, be = al_be[n]
                        f += mp.expjpi(-t * n) * (-2 * mp.pi ** 2 * c * (be - t * al) + 2j * mp.pi * d * al)
                    self.pts.append((t, half * w * f))
        self.atoms = [(mp.mpf(0), C * PL[0])] + [(mp.mpf(j) / 6, 2 * C * PL[j]) for j in range(1, 7)]

def gl_nodes(n):
    xs, ws = [], []
    # Legendre roots by Newton from Chebyshev guesses (mpmath precision)
    for k in range(1, n + 1):
        x = mp.cos(mp.pi * (k - mp.mpf(1) / 4) / (n + mp.mpf(1) / 2))
        for _ in range(100):
            p0, p1 = mp.mpf(1), x
            for m in range(2, n + 1):
                p0, p1 = p1, ((2 * m - 1) * x * p1 - (m - 1) * p0) / m
            dp = n * (x * p1 - p0) / (x ** 2 - 1)
            dx = p1 / dp
            x -= dx
            if abs(dx) < mp.mpf(10) ** (-mp.mp.dps - 5):
                break
        p0, p1 = mp.mpf(1), x
        for m in range(2, n + 1):
            p0, p1 = p1, ((2 * m - 1) * x * p1 - (m - 1) * p0) / m
        dp = n * (x * p1 - p0) / (x ** 2 - 1)
        xs.append(x); ws.append(2 / ((1 - x ** 2) * dp ** 2))
    return xs, ws

def lam_z(t):
    b, h, B = consts()
    lam = 1j / (b * (t + 1j * h))
    z = -B / (t + 1j * h) - 1j * h
    return lam, z

def K_derivs(M, s, kmax=0):
    """[K^{(k)}(s)] for k<=kmax via the folded measure."""
    out = [mp.mpc(0)] * (kmax + 1)
    for t, g in list(M.pts) + list(M.atoms):
        lam, z = lam_z(t)
        base = g * lam * mp.exp(1j * mp.pi * z * s)
        iz = 1j * mp.pi * z
        for k in range(kmax + 1):
            out[k] += base
            base *= iz
    return [o.real for o in out]

def p_derivs(M, s, kmax=0):
    out = [mp.mpc(0)] * (kmax + 1)
    for t, g in list(M.pts) + list(M.atoms):
        base = g * mp.expjpi(t * s)
        it = 1j * mp.pi * t
        for k in range(kmax + 1):
            out[k] += base
            base *= it
    return [o.real for o in out]

