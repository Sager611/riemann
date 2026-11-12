# xform.py (explore39): Weil-type form on I = [-a, a] in x-space, parity-pure orthonormal Legendre basis.
#   Q(f) = Arch(f) + pole(f) - 2 sum_n w_n C_f(log n),  C_f(h) = int f(x) f(x-h) dx,
#   Arch(f) = (D0 + 2K(2a))|f|^2 + int_0^{2a} kappa(z) (2|f|^2 - 2 C_f(z)) dz   (Levy-Khintchine, explore32 S0),
#   zeta: kappa(z) = e^{z/2}/(e^z - e^{-z}), D0 = psi(1/4) - ln pi, pole(f) = 2 <f,e^{-x/2}><f,e^{x/2}>;
#   dh  : kappa(z) = e^{-z/2}/(e^z - e^{-z}), D0 = ln(5/pi) + psi(3/4), no pole, atoms b(n)/sqrt n (wmat.dh_coeffs).
# Exploration code, float64 (not a certificate).
import numpy as np, math
from numpy.polynomial import legendre as Lg
from scipy.special import digamma

def ppow(n):
    for p in range(2, n + 1):
        if n % p == 0:
            m = n; e = 0
            while m % p == 0: m //= p; e += 1
            return (p, e) if m == 1 else None
    return None

def dh_b(nmax):
    k = (math.sqrt(10 - 2*math.sqrt(5)) - 2)/(math.sqrt(5) - 1)
    per = [0.0, 1.0, k, -k, -1.0]
    a = [0.0] + [per[n % 5] for n in range(1, nmax + 1)]
    b = [0.0]*(nmax + 1)
    for n in range(2, nmax + 1):
        s = a[n]*math.log(n)
        for d in range(2, n):
            if n % d == 0: s -= a[d]*b[n//d]
        b[n] = s
    return b, k

class Form:
    def __init__(s, c, kind='zeta'):
        s.c = float(c); s.a = math.log(s.c)/2; s.kind = kind
        if kind == 'zeta':
            s.D0 = float(digamma(0.25)) - math.log(math.pi); s.pole = True; s.sg = +1
            s.atoms = [(math.log(n), math.log(ppow(n)[0])/math.sqrt(n)) for n in range(2, int(math.ceil(s.c))) if n < s.c and ppow(n)]
        elif kind == 'dh':
            s.D0 = math.log(5/math.pi) + float(digamma(0.75)); s.pole = False; s.sg = -1
            b, _ = dh_b(int(s.c) + 2)
            s.atoms = [(math.log(n), b[n]/math.sqrt(n)) for n in range(2, int(math.ceil(s.c))) if n < s.c and abs(b[n]) > 0]
        else: raise ValueError(kind)
    # kernel pieces
    def kap(s, z):
        z = np.asarray(z, float); return np.exp(s.sg*z/2)/(np.exp(z) - np.exp(-z))
    def zkap(s, z):   # z*kappa(z), analytic at 0
        z = np.asarray(z, float); zz = np.where(z == 0, 1.0, z)
        return np.where(z == 0, 0.5, 0.5*np.exp(s.sg*zz/2)*zz/np.sinh(zz))
    def K(s, r):      # int_r^inf kappa
        r = np.asarray(r, float)
        return 0.5*np.log(1/np.tanh(r/4)) + s.sg*np.arctan(np.exp(-r/2))
    # basis
    def degs(s, n, sector): return [2*k + (0 if sector == 'even' else 1) for k in range(n)]
    def B(s, x, n, sector):
        d = s.degs(n, sector); V = Lg.legvander(np.asarray(x, float)/s.a, d[-1])
        return V[:, d]*np.sqrt((2*np.array(d) + 1)/(2*s.a))
    def gauss(s, lo, hi, m):
        t, w = Lg.leggauss(m); return (hi + lo)/2 + (hi - lo)/2*t, (hi - lo)/2*w
    def Cmat(s, z, n, sector, m):   # C_ij(z) = int p_i(x) p_j(x - z) dx, 0 <= z <= 2a
        xs, ws = s.gauss(-s.a + z, s.a, m)
        return (s.B(xs, n, sector)*ws[:, None]).T @ s.B(xs - z, n, sector)
    def parts(s, n, sector, nz=None, nx=None):
        nz = nz or 2*n + 60; nx = nx or 2*n + 30
        I = np.eye(n); a = s.a
        zs, wz = s.gauss(0.0, 2*a, nz)
        A = (s.D0 + 2*s.K(2*a))*I
        for z, w in zip(zs, wz):
            Cz = s.Cmat(z, n, sector, nx); Cz = (Cz + Cz.T)/2
            A += w*s.zkap(z)*(2*I - 2*Cz)/z
        xs, ws = s.gauss(-a, a, nx + 20); Bx = s.B(xs, n, sector)
        em = (Bx*ws[:, None]).T @ np.exp(-xs/2); ep = (Bx*ws[:, None]).T @ np.exp(xs/2)
        P = (np.outer(em, ep) + np.outer(ep, em)) if s.pole else np.zeros((n, n))
        Pr = np.zeros((n, n))
        for h, w in s.atoms:
            Ch = s.Cmat(h, n, sector, nx); Pr += w*(Ch + Ch.T)
        return {'arch': A, 'pole': P, 'prime': Pr, 'em': em, 'ep': ep}
    def Q(s, n, sector, **kw):
        p = s.parts(n, sector, **kw); return p['arch'] + p['pole'] - p['prime'], p
    # reduced atom pieces on [0, a]: Dir(h) = int_h^a g(x) g(x-h), Ref(h) = int g(x) g(h - x) over x, h-x in [0,a]
    def atom_pieces(s, h, n, sector, m):
        a = s.a; Dm = np.zeros((n, n)); Rm = np.zeros((n, n))
        if h < a:
            xs, ws = s.gauss(h, a, m); Bx = s.B(xs, n, sector); Bh = s.B(xs - h, n, sector)
            Dm = (Bx*ws[:, None]).T @ Bh; Dm = (Dm + Dm.T)/2
        lo, hi = max(0.0, h - a), min(a, h)
        if hi > lo:
            xs, ws = s.gauss(lo, hi, m); Rm = (s.B(xs, n, sector)*ws[:, None]).T @ s.B(h - xs, n, sector); Rm = (Rm + Rm.T)/2
        return Dm, Rm
    def jc(s, z):   # continuous off-diagonal jump kernel j^c(z) = kappa(|z|) - 2 cosh(z/2) [pole]
        z = np.abs(np.asarray(z, float)); v = s.kap(np.where(z == 0, 1e-300, z))
        return v - (2*np.cosh(z/2) if s.pole else 0.0)
    def jev(s, x, y): return s.jc(x - y) + s.jc(x + y)
    def jodd(s, x, y): return s.jc(x - y) - s.jc(x + y)
    def flip(s, n, sector, m2=400, m1=None):
        # matrix F of the anti-Markov parts' doubled forms: Q_flip = Q - F (term-wise Beurling-Deny minorant)
        m1 = m1 or 2*n + 40; a = s.a; F = np.zeros((n, n)); info = {}
        if sector == 'even':
            xs, ws = s.gauss(0.0, a, m2); X, Y = np.meshgrid(xs, xs, indexing='ij')
            jm = np.maximum(-s.jev(X, Y), 0.0); info['jev_min'] = float(-jm.max()) if jm.max() > 0 else float(np.min(s.jev(X, Y)))
            Bx = s.B(xs, n, sector)*ws[:, None]
            F += 4*Bx.T @ jm @ Bx
        for h, w in s.atoms:
            Dm, Rm = s.atom_pieces(h, n, sector, m1)
            if sector == 'even':
                if w < 0: F += 2*(4*(-w)*Dm + 2*(-w)*Rm)
            else:
                if w < 0: F += 2*(4*(-w)*Dm)
                if w > 0: F += 2*(2*w*Rm)
        return F, info
    # pointwise operator for Barta: (H phi)(x), phi = B @ v, full interval, plus flip operator
    def Hphi(s, v, sector, xs, mq=80, flip=False, m2=600):
        a = s.a; n = len(v)
        phi = lambda x: s.B(np.atleast_1d(x), n, sector) @ v
        out = np.zeros(len(xs))
        xg, wg = s.gauss(-a, a, 200); pg = phi(xg)
        A = np.sum(wg*pg*np.exp(-xg/2)); Bp = np.sum(wg*pg*np.exp(xg/2))
        t, w = Lg.leggauss(mq)
        for k, x in enumerate(xs):
            px = phi(x)[0]
            val = (s.D0 + s.K(a - x) + s.K(a + x))*px
            for lo, hi in ((-a, x), (x, a)):
                if hi - lo < 1e-14: continue
                ys = (hi + lo)/2 + (hi - lo)/2*t; wy = (hi - lo)/2*w
                d = np.abs(x - ys); val += np.sum(wy*s.zkap(d)/d*(px - phi(ys)))
            if s.pole: val += np.exp(-x/2)*Bp + np.exp(x/2)*A
            for h, wt in s.atoms:
                for y in (x - h, x + h):
                    if abs(y) <= a: val -= wt*phi(y)[0]
            out[k] = val
        if flip:
            # subtract 2 * (reduced anti-Markov operator) at x in [0, a]
            ys, wy = s.gauss(0.0, a, m2); py = phi(ys)
            for k, x in enumerate(xs):
                corr = 0.0
                if sector == 'even':
                    jm = np.maximum(-s.jev(x, ys), 0.0); corr += 2*np.sum(wy*jm*py)
                for h, wt in s.atoms:
                    if wt < 0:   # direct atoms anti-Markov in both sectors
                        for y in (x - h, x + h):
                            if 0 <= y <= a: corr += 2*(-wt)*phi(y)[0]
                    if 0 <= h - x <= a and x < h:   # reflected atom: anti-Markov if w<0 (even) or w>0 (odd)
                        if sector == 'even' and wt < 0: corr += 2*(-wt)*phi(h - x)[0]
                        if sector == 'odd' and wt > 0: corr += 2*wt*phi(h - x)[0]
                out[k] -= corr
        return out, phi

def Hmat(F, n, sector, xs, mq=80, flip=False, m2=600):
    """rows: x in xs; cols: basis p_j; entries (H p_j)(x) [minus flip operator if flip]."""
    import numpy as np
    from numpy.polynomial import legendre as Lg
    a = F.a; Bf = lambda x: F.B(np.atleast_1d(np.asarray(x, float)), n, sector)
    xg, wg = F.gauss(-a, a, 200); Bg = Bf(xg)
    Avec = (Bg*wg[:, None]).T @ np.exp(-xg/2); Bvec = (Bg*wg[:, None]).T @ np.exp(xg/2)
    t, w = Lg.leggauss(mq); out = np.zeros((len(xs), n))
    for k, x in enumerate(xs):
        bx = Bf(x)[0]; row = (F.D0 + F.K(a - x) + F.K(a + x))*bx
        for lo, hi in ((-a, x), (x, a)):
            if hi - lo < 1e-14: continue
            ys = (hi + lo)/2 + (hi - lo)/2*t; wy = (hi - lo)/2*w; d = np.abs(x - ys)
            row = row + ((wy*F.zkap(d)/d)[:, None]*(bx[None, :] - Bf(ys))).sum(0)
        if F.pole: row = row + np.exp(-x/2)*Bvec + np.exp(x/2)*Avec
        for h, wt in F.atoms:
            for y in (x - h, x + h):
                if abs(y) <= a: row = row - wt*Bf(y)[0]
        out[k] = row
    if flip:
        ys, wy = F.gauss(0.0, a, m2); By = Bf(ys)
        for k, x in enumerate(xs):
            corr = np.zeros(n)
            if sector == 'even':
                jm = np.maximum(-F.jev(x, ys), 0.0); corr += 2*(wy*jm) @ By
            for h, wt in F.atoms:
                if wt < 0:
                    for y in (x - h, x + h):
                        if 0 <= y <= a: corr += 2*(-wt)*Bf(y)[0]
                if 0 <= h - x <= a and x < h:
                    if sector == 'even' and wt < 0: corr += 2*(-wt)*Bf(h - x)[0]
                    if sector == 'odd' and wt > 0: corr += 2*wt*Bf(h - x)[0]
            out[k] -= corr
    return out

def barta_lp(F, n, sector, flip=False, npts=400, tol=1e-6, v0=None):
    """best Barta bound max_{phi in span, phi>0 on grid} min_grid (H phi)/phi, by LP bisection; returns (bound, v, xs)."""
    import numpy as np
    from scipy.optimize import linprog
    a = F.a; u = (np.arange(npts) + 0.5)/npts; xs = a*np.sin(np.pi*u/2)
    if sector == 'odd': xs = xs[xs > 2e-3*a]
    HB = Hmat(F, n, sector, xs, flip=flip); BB = F.B(xs, n, sector)
    if sector == 'odd':   # odd: divide by x to keep the ratio well-scaled near 0 (both vanish linearly)
        HB = HB/xs[:, None]; BB = BB/xs[:, None]
    wq = np.gradient(xs)
    def feas(t):
        A_ub = -(HB - t*BB); b_ub = np.zeros(len(xs))
        A_ub = np.vstack([A_ub, -BB]); b_ub = np.concatenate([b_ub, -1e-9*np.ones(len(xs))])
        A_eq = (wq @ BB)[None, :]; b_eq = [1.0]
        r = linprog(np.zeros(n), A_ub=A_ub, b_ub=b_ub, A_eq=A_eq, b_eq=b_eq, bounds=[(None, None)]*n, method='highs')
        return r.status == 0, (r.x if r.status == 0 else None)
    lo, hi = -50.0, 50.0
    if v0 is not None:
        ph = BB @ v0; 
        if np.median(ph) < 0: v0 = -v0; ph = -ph
        if np.all(ph > 0): lo = float(np.min((HB @ v0)/ph))
    ok, vbest = feas(lo)
    if not ok: return None, None, xs
    for _ in range(60):
        mid = (lo + hi)/2; ok, v = feas(mid)
        if ok: lo, vbest = mid, v
        else: hi = mid
        if hi - lo < tol*max(1.0, abs(lo)): break
    return lo, vbest, xs

