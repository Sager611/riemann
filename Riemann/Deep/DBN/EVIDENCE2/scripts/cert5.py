# cert5.py -- paired-probe certificate (pairwise relaxation, split grid), Euler B&B, floor assembly.
import numpy as np, math, json, time
from fractions import Fraction as Fr
from flint import arb, fmpq, ctx  # pyright: ignore[reportMissingImports]
import lib5 as L
from lib5 import up, dn, fu, fd, aq, sum_up, U53, INF

# ---------------- interval arithmetic on numpy arrays (outward rounding) ----------------
class IA:
    __slots__ = ('lo', 'hi')
    def __init__(s, lo, hi=None):
        s.lo = np.asarray(lo, dtype=np.float64); s.hi = s.lo if hi is None else np.asarray(hi, dtype=np.float64)
    def __add__(s, o):
        o = _ia(o); return IA(dn(s.lo + o.lo), up(s.hi + o.hi))
    __radd__ = __add__
    def __sub__(s, o):
        o = _ia(o); return IA(dn(s.lo - o.hi), up(s.hi - o.lo))
    def __rsub__(s, o): return _ia(o) - s
    def __mul__(s, o):
        o = _ia(o); p = [s.lo*o.lo, s.lo*o.hi, s.hi*o.lo, s.hi*o.hi]
        return IA(dn(np.minimum(np.minimum(p[0], p[1]), np.minimum(p[2], p[3]))), up(np.maximum(np.maximum(p[0], p[1]), np.maximum(p[2], p[3]))))
    __rmul__ = __mul__
    def __neg__(s): return IA(-s.hi, -s.lo)
    def abs(s):
        lo = np.where(s.lo >= 0, s.lo, np.where(s.hi <= 0, -s.hi, 0.0))
        return IA(lo, np.maximum(np.abs(s.lo), np.abs(s.hi)))
    def sq(s):
        a = s.abs(); return IA(dn(a.lo*a.lo), up(a.hi*a.hi))
    def sqrt(s):
        return IA(dn(np.sqrt(np.maximum(s.lo, 0.0))), up(np.sqrt(np.maximum(s.hi, 0.0))))
def _ia(o): return o if isinstance(o, IA) else IA(o)
def isum_hi(x):
    """rigorous upper bound of the exact sum of an array of upper bounds (any sign)."""
    x = np.asarray(x, dtype=np.float64).ravel()
    n = x.size
    if n == 0: return 0.0
    e = float(np.sum(np.abs(x)))*(4.0*(n + 2)*U53)*(1 + 1e-12)
    return float(up(up(float(np.sum(x)) + up(e))))

# ---------------- interval trigonometry with arb ----------------
PI = arb.pi()
def cos_range(t1, t2):
    """[min, max] of cos on [t1, t2] (arb endpoints, t1 <= t2), returned as floats (outward)."""
    if (t2 - t1 > 2*PI): return -1.0, 1.0
    c1 = t1.cos(); c2 = t2.cos()
    lo = min(fd(c1), fd(c2)); hi = max(fu(c1), fu(c2))
    n1 = (t1/PI).floor(); n2 = (t2/PI).floor()
    a = n1.unique_fmpz(); b = n2.unique_fmpz()
    if a is None or b is None:
        lo2 = -1.0; hi2 = 1.0                           # undecided: be safe
        return lo2, hi2
    for n in range(int(a) + 1, int(b) + 1):              # multiples n*pi inside (t1, t2]
        if n % 2 == 0: hi = 1.0
        else: lo = -1.0
    return max(lo, -1.0), min(hi, 1.0)
def abscs_max(t1, t2):
    """(max |cos|, max |sin|) on [t1, t2]."""
    if (t2 - t1 > PI): return 1.0, 1.0
    s1_, c1_ = t1.sin_cos(); s2_, c2_ = t2.sin_cos()
    cm_ = max(fu(abs(c1_)), fu(abs(c2_))); sm_ = max(fu(abs(s1_)), fu(abs(s2_)))
    q1 = (2*t1/PI).floor().unique_fmpz(); q2 = (2*t2/PI).floor().unique_fmpz()
    if q1 is None or q2 is None: return 1.0, 1.0
    for j in range(int(q1) + 1, int(q2) + 1):
        if j % 2 == 0: cm_ = 1.0
        else: sm_ = 1.0
    return min(cm_, 1.0), min(sm_, 1.0)
def abscos_max(t1, t2):
    lo, hi = cos_range(t1, t2); return max(abs(lo), abs(hi))
def abssin_max(t1, t2):
    h = PI/2
    lo, hi = cos_range(t1 - h, t2 - h)                  # sin t = cos(t - pi/2)
    return max(abs(lo), abs(hi))

_FCACHE = {}
def clear_cache(): _FCACHE.clear()
# ---------------- bins and pair groups for one rectangle ----------------
class Geometry:
    def __init__(self, R, Psi):
        self.R = R
        m = R.masses(Psi)
        self.m = m
        lin = R.linked_refl
        fF = m['MF'] > 0
        kind = [np.zeros(fF.sum(), dtype=np.int64)]; kk = [m['kF'][fF]]; M = [m['MF'][fF]]; a = [m['aF'][fF]]
        llo = [m['lamF_lo'][fF]]; lhi = [m['lamF_hi'][fF]]
        fR = m['MR'] > 0
        AR = sum_up(m['MR']); BR = sum_up(up(m['MR']*m['lamR_hi']))
        if lin:
            kind.append(np.ones(fR.sum(), dtype=np.int64)); kk.append(m['kR'][fR]); M.append(m['MR'][fR]); a.append(m['aR'][fR])
            llo.append(m['lamR_lo'][fR]); lhi.append(m['lamR_hi'][fR])
        self.kind = np.concatenate(kind); self.k = np.concatenate(kk); self.M = np.concatenate(M); self.a = np.concatenate(a)
        self.lam_lo = np.concatenate(llo); self.lam_hi = np.concatenate(lhi)
        nb = self.M.size; self.nb = nb
        self.A_lin = sum_up(self.M); self.B_lin = sum_up(up(self.M*self.lam_hi))
        self.Dg = sum_up(up(up(self.a*self.M)*self.lam_hi))
        A_un = R.Aact; B_un = R.Bact
        if not lin:
            A_un = A_un + arb(AR); B_un = B_un + arb(BR)
        self.A_un = fu(A_un); self.B_un = fu(B_un)
        Un = float(up(up(self.B_un*up(up(1 + self.A_lin) + self.A_un)) + up(self.A_un*self.B_lin)))
        # boundary labels n* (one cutoff index can differ between the probes): unlinked correction
        lam_max = float(np.max(self.lam_hi)) if self.nb else 0.0
        bd = m['bd_amp']
        self.Un = float(up(Un + up(bd*up(up(lam_max*up(1 + self.A_lin + self.A_un)) + up(self.B_lin + self.B_un)))))
        self.Aplus = float(up(self.A_lin + self.A_un)); self.Bplus = float(up(self.B_lin + self.B_un))
        # pairs I <= J
        I, J = np.triu_indices(nb)
        self.I = I; self.J = J
        W = up(self.M[I]*self.M[J]); W = np.where(I == J, up(W*0.5), W)
        self.W = W
        tI = self.kind[I]; tJ = self.kind[J]
        typ = tI + tJ                                  # 0 FF, 1 FR, 2 RR
        kI = self.k[I]; kJ = self.k[J]
        kF = np.where(tI == 0, kI, kJ); kRr = np.where(tI == 0, kJ, kI)
        self.typ = typ
        self.delta = np.where(typ == 1, kF - kRr, np.abs(kI - kJ))
        key = np.where(typ == 1, kF + kRr, kI + kJ)
        # group ids (merged S-keys: FF keys >= 200 in pairs, FR/RR in triples)
        mfac = np.where(typ == 0, np.where(key >= 200, 2, 1), 3)
        k0 = np.where(mfac == 1, key, (key//mfac)*mfac)
        gkey = typ*10**8 + (k0 + 500000)*10 + mfac
        ug, gid = np.unique(gkey, return_inverse=True)
        self.gid = gid; self.ng = ug.size
        gt = ug//10**8; rest = ug - gt*10**8; gm = rest % 10; gk = rest//10 - 500000
        self.gtyp = gt; self.gk = gk
        w = 1.0/20
        L1lo = m['L1lo']; L1hi = m['L1hi']; om = m['omega']
        top = gk + gm + 1
        Slo = np.where(gt == 0, dn(gk*w), np.where(gt == 1, dn(L1lo + dn(gk*w)), dn(dn(2*L1lo) + dn(gk*w))))
        Shi = np.where(gt == 0, up(top*w), np.where(gt == 1, up(up(L1hi + om) + up(top*w)), up(up(2*L1hi + 2*om) + up(top*w))))
        self.Slo = dn(Slo); self.Shi = up(Shi)
        # unique (type, delta) table index per pair
        tdk = typ*10**6 + self.delta + 500000
        self.td_u, self.td_i = np.unique(tdk, return_inverse=True)
        self.td_t = self.td_u//10**6; self.td_d = self.td_u % 10**6 - 500000
        fac = 1.0 + 4.0*(W.size + 2)*U53
        self.W0 = up(np.bincount(gid, W, minlength=self.ng)*fac)
        self.fac = fac

    def shift_data(self, a_frac):
        """a = r*h (Fraction). Per-bin cos(a*lambda) ranges and per-group Wc, Ws (all rigorous)."""
        a = aq(Fr(a_frac)); R = self.R
        wid = arb(fmpq(2, 10**10))
        w = arb(fmpq(1, 20)); om = (R.L2 - R.L1 if R.L2 is not None else arb(0)) + arb(fmpq(L.SUB, L.NG))
        # per bin cos(a*lambda) range
        Clo = np.empty(self.nb); Chi = np.empty(self.nb)
        fc = _FCACHE.setdefault(('bin', Fr(a_frac)), {})
        for b in range(self.nb):
            k = int(self.k[b])
            if self.kind[b] == 0:
                if k not in fc: fc[k] = cos_range(a*k*w, a*(k + 1)*w)
                Clo[b], Chi[b] = fc[k]; continue
            else:
                t1 = a*(R.L1 + k*w) - wid; t2 = a*(R.L1 + (k + 1)*w + om) + wid
            Clo[b], Chi[b] = cos_range(t1, t2)
        # pair tables over unique (type, delta)
        half = a/2
        nt = self.td_u.size; cmx = np.empty(nt); smx = np.empty(nt)
        ffc = _FCACHE.setdefault(('ff', Fr(a_frac)), {})
        for i in range(nt):
            t_ = int(self.td_t[i]); dlt = int(self.td_d[i])
            if t_ == 0:
                if dlt not in ffc: ffc[dlt] = abscs_max(half*(dlt - 1)*w, half*(dlt + 1)*w)
                cmx[i], smx[i] = ffc[dlt]; continue
            elif t_ == 1:
                th1 = half*((dlt - 1)*w - R.L1 - om) - wid; th2 = half*((dlt + 1)*w - R.L1) + wid
            else:
                th1 = half*((dlt - 1)*w - om) - wid; th2 = half*((dlt + 1)*w + om) + wid
            cmx[i], smx[i] = abscs_max(th1, th2)
        cmax = cmx[self.td_i]; smax = smx[self.td_i]
        Wc = up(np.bincount(self.gid, up(self.W*cmax), minlength=self.ng)*self.fac)
        Ws = up(np.bincount(self.gid, up(self.W*smax), minlength=self.ng)*self.fac)
        return dict(Clo=Clo, Chi=Chi, Wc=Wc, Ws=Ws)

# ---------------- Psi: float (search) and interval (proof) versions ----------------
def psi_float(G, S, cm, cp, mu):
    """cm, cp scalars; mu array (n,) -> Psi values (n,)"""
    mu = np.atleast_1d(mu)[:, None]
    best = None
    for lam in (G.lam_lo, G.lam_hi):
        al = 2*cm - lam; be = mu*(2*cp - lam)
        for C in (S['Clo'], S['Chi']):
            v = al*al + be*be + 2*al*be*C
            best = v if best is None else np.maximum(best, v)
    s0 = np.sqrt(np.maximum(best, 0.0)) @ G.M
    P = np.maximum(np.abs((1 + mu)*G.Slo - 2*(cm + mu*cp)), np.abs((1 + mu)*G.Shi - 2*(cm + mu*cp)))
    Rr = np.maximum(np.abs(2*(cm - mu*cp) - (1 - mu)*G.Slo), np.abs(2*(cm - mu*cp) - (1 - mu)*G.Shi))
    Q = np.maximum(np.abs(2*cm - G.Slo) + mu*np.abs(2*cp - G.Slo), np.abs(2*cm - G.Shi) + mu*np.abs(2*cp - G.Shi))
    b = np.minimum(P*S['Wc'] + Rr*S['Ws'], Q*G.W0).sum(axis=1)
    mu = mu[:, 0]
    return (cm + mu*cp) - s0 - b - (1 + mu)*G.Dg - (1 + mu)*G.Un

def psi_ia(G, S, cm, cp, mu):
    """rigorous lower bound of Psi for exact float parameters cm, cp, mu."""
    cmI = IA(cm); cpI = IA(cp); muI = IA(mu)
    best = None
    for lam in (G.lam_lo, G.lam_hi):
        al = 2*cmI - IA(lam); be = muI*(2*cpI - IA(lam))
        for C in (S['Clo'], S['Chi']):
            v = al.sq() + be.sq() + 2*(al*be)*IA(C)
            best = v.hi if best is None else np.maximum(best, v.hi)
    g0 = up(np.sqrt(np.maximum(best, 0.0)))
    s0 = sum_up(up(g0*G.M))
    one = IA(1.0)
    def endmax(f):
        return np.maximum(f(IA(G.Slo)).abs().hi, f(IA(G.Shi)).abs().hi)
    P = endmax(lambda S_: (one + muI)*S_ - 2*(cmI + muI*cpI))
    Rr = endmax(lambda S_: 2*(cmI - muI*cpI) - (one - muI)*S_)
    Q = np.maximum((2*cmI - IA(G.Slo)).abs().hi + up(mu*(2*cpI - IA(G.Slo)).abs().hi),
                   (2*cmI - IA(G.Shi)).abs().hi + up(mu*(2*cpI - IA(G.Shi)).abs().hi))
    Q = up(Q)
    b1 = up(up(P*S['Wc']) + up(Rr*S['Ws'])); b2 = up(Q*G.W0)
    bsum = sum_up(np.minimum(b1, b2))
    pos = (cmI + muI*cpI).lo
    neg = float(up(up(bsum + s0) + up(up((1 + mu)*G.Dg)*(1 + 2*U53)) + up(up((1 + mu)*G.Un)*(1 + 2*U53))))
    neg = float(up(up(bsum + s0) + up(((one + muI)*IA(G.Dg)).hi) + up(((one + muI)*IA(G.Un)).hi)))
    return float(dn(float(pos) - neg))

def psi_ia_multi(G, S, cm, cp, mu):
    """rigorous lower bounds of Psi for arrays (K,) of exact float parameters."""
    cm = np.asarray(cm, dtype=np.float64)[:, None]; cp = np.asarray(cp, dtype=np.float64)[:, None]; mu = np.asarray(mu, dtype=np.float64)[:, None]
    cmI = IA(cm); cpI = IA(cp); muI = IA(mu)
    best = None
    for lam in (G.lam_lo, G.lam_hi):
        lamI = IA(lam[None, :])
        al = 2*cmI - lamI; be = muI*(2*cpI - lamI)
        for C_ in (S['Clo'], S['Chi']):
            v = al.sq() + be.sq() + 2*(al*be)*IA(C_[None, :])
            best = v.hi if best is None else np.maximum(best, v.hi)
    g0 = up(np.sqrt(np.maximum(best, 0.0)))
    prod = up(g0*G.M[None, :])
    n = prod.shape[1]
    s0 = up(up(np.sum(prod, axis=1)*(1.0 + 4.0*(n + 2)*U53)))
    one = IA(1.0)
    SloI = IA(G.Slo[None, :]); ShiI = IA(G.Shi[None, :])
    def endmax(f): return np.maximum(f(SloI).abs().hi, f(ShiI).abs().hi)
    P = endmax(lambda S_: (one + muI)*S_ - 2*(cmI + muI*cpI))
    Rr = endmax(lambda S_: 2*(cmI - muI*cpI) - (one - muI)*S_)
    Q = up(np.maximum(up((2*cmI - SloI).abs().hi + (muI*(2*cpI - SloI).abs()).hi),
                      up((2*cmI - ShiI).abs().hi + (muI*(2*cpI - ShiI).abs()).hi)))
    b1 = up(up(P*S['Wc'][None, :]) + up(Rr*S['Ws'][None, :])); b2 = up(Q*G.W0[None, :])
    bb = np.minimum(b1, b2); n2 = bb.shape[1]
    bsum = up(up(np.sum(bb, axis=1)*(1.0 + 4.0*(n2 + 2)*U53)))
    pos = (cmI + muI*cpI).lo[:, 0]
    tail = up(((one + muI)*IA(G.Dg)).hi[:, 0] + ((one + muI)*IA(G.Un)).hi[:, 0])
    neg = up(up(bsum + s0) + tail)
    return dn(pos - neg)

def best_mu(G, S, cm, cp, mu_hi=6.0):
    grid = np.linspace(0, mu_hi, 61)
    v = psi_float(G, S, cm, cp, grid)
    i = int(np.argmax(v))
    lo = grid[max(i - 1, 0)]; hi = grid[min(i + 1, len(grid) - 1)]
    g2 = np.linspace(lo, hi, 41); v2 = psi_float(G, S, cm, cp, g2); j = int(np.argmax(v2))
    return float(g2[j]), float(v2[j])

def one_point(G, S):
    f = lambda c: psi_float(G, S, c, c, np.array([0.0]))[0]
    lo, hi = 0.0, 6.0
    if f(hi) < 0: return None
    for _ in range(50):
        c = (lo + hi)/2
        if f(c) >= 0: hi = c
        else: lo = c
    c = float(up(hi*(1 + 1e-9)))
    while psi_ia(G, S, c, c, 0.0) < 0:
        c = float(up(c*(1 + 1e-6)))
    return c

def paired(G, S, Ctop, K=16, tol=1e-12):
    def feas(c, keep=False):
        grid = np.linspace(c, Ctop, K + 1); grid[-1] = Ctop
        mus = []
        for k in range(K):
            cm = 2*c - grid[k + 1]; cp = grid[k]
            if cm <= 0: return None
            mu, v = best_mu(G, S, cm, cp)
            if v < tol: return None
            mus.append(mu)
        return grid, mus
    lo, hi = 0.0, Ctop
    for _ in range(26):
        c = (lo + hi)/2
        if feas(c) is not None: hi = c
        else: lo = c
    # rigorous verification (raise c slightly until every grid point verifies)
    c = float(up(hi))
    for attempt in range(30):
        r = feas(c)
        if r is not None:
            grid, mus = r
            ok = all(psi_ia(G, S, float(2*c - grid[k + 1]), float(grid[k]), float(mus[k])) >= 0 for k in range(K))
            if ok: return c, [float(x) for x in grid], mus
        c = float(up(c + 1e-6*Ctop))
    return None

# ---------------- Euler part: paired supremum (grid + second-derivative bound) ----------------
EULER_N = 1024
_EULER_CACHE = {}
def euler_pair(rmin, rmax, a_frac):
    """Rigorous upper bound of sup_theta (1/2) sum_{+-} Re(E_h/E)(theta +- (a/2) log p), summed over p = 2,3,5.
    f_p(t) = (log p/4)[g(t+d)+g(t-d)], g(t) = (rho cos t - rho^2)/(1 - 2 rho cos t + rho^2) = sum_k rho^k cos kt.
    On each grid cell f <= max(endpoint values) + M2 h^2/8 with |f''| <= (log p/2) rho(1+rho)/(1-rho)^3."""
    a = aq(Fr(a_frac)); total = arb(0)
    hstep = PI/EULER_N
    for p in (2, 3, 5):
        lp = L.LOGP[p]; dlt = a/2*lp
        rho = rmax[p]; rh = rmax[p]      # sup_theta f is nondecreasing in rho (Poisson-kernel averaging)
        best = -1e300
        th = [hstep*i for i in range(EULER_N + 1)]
        for t_ in th:
            cp_ = (t_ + dlt).cos(); cm_ = (t_ - dlt).cos()
            v = lp/4*((rho*cp_ - rho*rho)/(1 - 2*rho*cp_ + rho*rho) + (rho*cm_ - rho*rho)/(1 - 2*rho*cm_ + rho*rho))
            best = max(best, fu(v))
        M2 = lp/2*rh*(1 + rh)/(1 - rh)**3
        total += arb(best) + M2*hstep*hstep/8
    return total

def euler_one(rmax):
    t = arb(0)
    for p in (2, 3, 5):
        t += L.LOGP[p]/2*rmax[p]/(1 - rmax[p])
    return t

def phys_err(R, Aplus, Bplus):
    En = R.En; Ey = R.Ey
    R0 = arb(fmpq(2, 10**6))*En; R1 = arb(fmpq(2, 10**4))*En + arb(fmpq(2, 10**6))*Ey
    g = 1 - arb(Aplus)
    assert Aplus < 0.8 and (g - R0) > 0
    return (R1 + arb(Bplus)*R0/(2*g))/(g - R0)

def floor_value(R, Gm, c_half_term, Eterm):
    """L1/4 - 1e-6 - physerr - Eterm - c/2  (arb)"""
    pe = phys_err(R, Gm.Aplus, Gm.Bplus)
    return R.L1/4 - arb(fmpq(1, 10**6)) - pe - Eterm - arb(c_half_term)/2

# ---------------- fast vectorised search ----------------
def psi_float_multi(G, S, cm, cp, mu):
    """cm, cp: (K,), mu: (K, n) -> Psi (K, n) in float (search only)."""
    cm = np.asarray(cm)[:, None, None]; cp = np.asarray(cp)[:, None, None]; mu3 = np.asarray(mu)[:, :, None]
    best = None
    for lam in (G.lam_lo, G.lam_hi):
        al = 2*cm - lam; be = mu3*(2*cp - lam)
        for C in (S['Clo'], S['Chi']):
            v = al*al + be*be + 2*al*be*C
            best = v if best is None else np.maximum(best, v)
    s0 = np.sqrt(np.maximum(best, 0.0)) @ G.M
    Slo = G.Slo; Shi = G.Shi
    t2 = 2*(cm + mu3*cp); m1 = 1 + mu3
    P = np.maximum(np.abs(m1*Slo - t2), np.abs(m1*Shi - t2))
    r0 = 2*(cm - mu3*cp); r1 = 1 - mu3
    Rr = np.maximum(np.abs(r0 - r1*Slo), np.abs(r0 - r1*Shi))
    Q = np.maximum(np.abs(2*cm - Slo) + mu3*np.abs(2*cp - Slo), np.abs(2*cm - Shi) + mu3*np.abs(2*cp - Shi))
    b = np.minimum(P*S['Wc'] + Rr*S['Ws'], Q*G.W0).sum(axis=2)
    mu2 = np.asarray(mu)
    return (cm[:, :, 0] + mu2*cp[:, :, 0]) - s0 - b - (1 + mu2)*(G.Dg + G.Un)

GR = (math.sqrt(5) - 1)/2
def eval_grid(G, S, c, Ctop, K, mu_init=None):
    grid = np.linspace(c, Ctop, K + 1); grid[-1] = Ctop
    cm = 2*c - grid[1:]; cp = grid[:-1]
    if np.any(cm <= 0): return -1e9, grid, None
    ar = np.arange(K)
    if mu_init is None:
        sc = np.linspace(0, 4, 9)
        v = psi_float_multi(G, S, cm, cp, np.repeat(sc[None, :], K, 0)); i = np.argmax(v, 1)
        a_ = sc[np.maximum(i - 1, 0)]; b_ = sc[np.minimum(i + 1, 8)]
        it = 11
    else:
        a_ = np.maximum(np.asarray(mu_init) - 0.25, 0.0); b_ = np.asarray(mu_init) + 0.25
        it = 9
    x1 = b_ - GR*(b_ - a_); x2 = a_ + GR*(b_ - a_)
    f = psi_float_multi(G, S, cm, cp, np.stack([x1, x2], 1)); f1 = f[:, 0]; f2 = f[:, 1]
    for _ in range(it):
        m = f1 < f2
        an = np.where(m, x1, a_); bn = np.where(m, b_, x2)
        x1n = np.where(m, x2, bn - GR*(bn - an)); x2n = np.where(m, an + GR*(bn - an), x1)
        fe = psi_float_multi(G, S, cm, cp, np.where(m, x2n, x1n)[:, None])[:, 0]
        f1n = np.where(m, f2, fe); f2n = np.where(m, fe, f1)
        a_, b_, x1, x2, f1, f2 = an, bn, x1n, x2n, f1n, f2n
    mub = np.where(f1 >= f2, x1, x2); best = np.maximum(f1, f2)
    return float(np.min(best)), grid, mub

def root_increasing(f, lo, hi, flo=None, fhi=None, tol=2e-6, maxit=40):
    """smallest x with f(x) >= 0 (f increasing), Illinois method; returns hi-side point with f >= 0."""
    if flo is None: flo = f(lo)
    if fhi is None: fhi = f(hi)
    if flo >= 0: return lo
    if fhi < 0: return None
    side = 0
    for _ in range(maxit):
        x = (lo*fhi - hi*flo)/(fhi - flo)
        x = min(max(x, lo + 0.05*(hi - lo)*0.0 + 1e-15), hi)
        fx = f(x)
        if fx >= 0:
            hi, fhi = x, fx
            if side == 1: flo /= 2
            side = 1
        else:
            lo, flo = x, fx
            if side == -1: fhi /= 2
            side = -1
        if hi - lo < tol: break
    return hi

def one_point_fast(G, S):
    f = lambda c: psi_float(G, S, c, c, np.array([0.0]))[0]
    hi = 6.0
    while f(hi) < 0 and hi < 200: hi *= 2
    c = root_increasing(f, 0.0, hi, tol=1e-9)
    if c is None: return None
    c = float(up(c*(1 + 1e-9)))
    for _ in range(200):
        if psi_ia(G, S, c, c, 0.0) >= 0: return c
        c = float(up(c*(1 + 1e-7)))
    return None

def bisect_root(f, lo, hi, flo, fhi, tol=2e-6, maxit=40):
    """f increasing, flo < 0 <= fhi: regula falsi (Illinois) with bisection fallback; returns hi."""
    side = 0
    for it in range(maxit):
        if hi - lo < tol: break
        if it % 3 == 2: x = (lo + hi)/2
        else:
            x = (lo*fhi - hi*flo)/(fhi - flo)
            x = min(max(x, lo + 0.02*(hi - lo)), hi - 0.02*(hi - lo))
        fx = f(x)
        if fx >= 0:
            hi, fhi = x, fx
            if side == 1: flo /= 2
            side = 1
        else:
            lo, flo = x, fx
            if side == -1: fhi /= 2
            side = -1
    return hi

def certify_grid(G, S, c, Ctop, K):
    v, grid, mub = eval_grid(G, S, c, Ctop, K)
    if v < 0 or mub is None: return None
    cm = 2*c - grid[1:]; cp = grid[:-1]
    lb = psi_ia_multi(G, S, cm, cp, mub)
    if np.all(lb >= 0): return [float(x) for x in grid], [float(x) for x in mub]
    return None

def paired_fast(G, S, Ctop, K=32, c_hint=None):
    st: dict = {'mu': None}
    def fK(c):
        v, grid, mub = eval_grid(G, S, c, Ctop, K, mu_init=st['mu'])
        if mub is not None: st['mu'] = mub
        return v
    if c_hint is not None:
        lo, hi = max(0.3*Ctop, c_hint - 0.006), min(Ctop, c_hint + 0.006)
    else:
        lo, hi = 0.55*Ctop, 0.9*Ctop
    flo = fK(lo)
    while flo >= 0 and lo > 0.3*Ctop:
        hi, lo = lo, max(0.3*Ctop, lo - 0.05*Ctop); flo = fK(lo)
    fhi = fK(hi)
    while fhi < 0:
        if hi >= Ctop: return None
        lo, flo = hi, fhi
        hi = min(Ctop, hi + 0.05*Ctop); fhi = fK(hi)
    c = bisect_root(fK, lo, hi, flo, fhi, tol=8e-6) if flo < 0 else lo
    for attempt in range(40):
        r = certify_grid(G, S, c, Ctop, K)
        if r is not None: return float(c), r[0], r[1]
        c = float(up(c + 2e-6))
    return None

def verify_cert(G, S, c, Ctop, K):
    """re-optimise mu for a prescribed (c, Ctop, K) and verify rigorously; returns (grid, mus) or None."""
    return certify_grid(G, S, c, Ctop, K)
