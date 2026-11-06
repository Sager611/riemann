#!/usr/bin/env python3
"""isozero certificate: explicit constants for Lamzouri's isolated-zero theorem (arXiv 2610.26841v2).

Run:  uv run --no-project --with python-flint --with mpmath python3 cert.py
All bounds are Arb balls; every printed lower bound is a lower endpoint, every upper bound an upper endpoint.
See DERIVATION.md for the formulas (section numbers are quoted in comments).
"""
import sys, math
from flint import arb, ctx

ctx.prec = 160
PI = arb.pi()
S2 = arb(2).sqrt()
A = 1 / S2                        # 1/sqrt(2)
SA, CA = A.sin(), A.cos()
COT = CA / SA
THETA = arb(4) / 7                # Conrey's mollifier length (limit theta -> 4/7)

def lo(x):  return x.lower()
def hi(x):  return x.upper()
def flo(x):   # float <= lower endpoint (directed rounding by one ulp)
    return math.nextafter(float(x.lower().mid()), -math.inf)
def fhi(x):   # float >= upper endpoint
    return math.nextafter(float(x.upper().mid()), math.inf)

# ---------- section 1: the Montgomery-Taylor kernel f0 and its transform K0 ----------
CMT = arb(1) / 2 + A * COT        # C_MT = 1/2 + cot(1/sqrt2)/sqrt2 = C(f0), DERIVATION 1.2

# Kernel family f_w(x) = cos(w x) / Z_w on [-1/2, 1/2], Z_w = 2 sin(w/2)/w (w = sqrt2: Montgomery-Taylor f0).
# Closed forms (DERIVATION 1.3): K_w(u), K_w(iy), C(f_w).
def K0_real(u, w=None):
    """K_w(u) = 2 int_0^{1/2} f_w(x) cos(2 pi u x) dx (u real), closed form via sinc."""
    if w is None:
        wp = S2 + 2 * PI * u
        wm = S2 - 2 * PI * u
        return ((wp / 2).sinc() + (wm / 2).sinc()) / (2 * S2 * SA)
    wp = w + 2 * PI * u
    wm = w - 2 * PI * u
    return w * ((wp / 2).sinc() + (wm / 2).sinc()) / (4 * (w / 2).sin())

def K0_imag(y, w=None):
    """K_w(iy) = 2 int_0^{1/2} f_w(x) cosh(2 pi x y) dx; for w = sqrt2 this is paper l.785."""
    py = PI * y
    if w is None:
        return (py.cosh() + S2 * py * COT * py.sinh()) / (1 + 2 * PI**2 * y**2)
    cw = (w / 2).cos() / (w / 2).sin()
    return w**2 * (py.cosh() + (2 * py / w) * cw * py.sinh()) / (w**2 + 4 * PI**2 * y**2)

def C_of(w):
    """C(f_w) = int f_w^2 + iint |u-v| f_w(u) f_w(v) du dv (closed form, DERIVATION 1.3)."""
    Z = 2 * (w / 2).sin() / w
    I = arb(1) / 2 + w.sin() / (2 * w)          # int_{-1/2}^{1/2} cos^2(w x) dx
    return I * (1 - 2 / w**2) / Z**2 + ((w / 2).sin() / w + 2 * (w / 2).cos() / w**2) / Z

def H(y, thb, w=None):
    """H_tau(y) = 1 - (1-thb)^2 * 4 K(iy)^2 / (1 + K(2iy)), DERIVATION 3.4 (thb = pi^2 tau^2 / 2)."""
    k1 = K0_imag(y, w)
    k2 = K0_imag(2 * y, w)
    return 1 - (1 - thb)**2 * 4 * k1**2 / (1 + k2)

def check_C_numeric():
    """Independent check of C(f0) = int f0^2 + iint |u-v| f0 f0 = C_MT (DERIVATION 1.2), by Arb quadrature."""
    from flint import acb
    c = 1 / (S2 * SA)
    def f0(x):  return c * (S2 * x).cos()
    def G(u):   # int |u-v| f0(v) dv, closed form
        return c * (SA / S2 + CA - (S2 * u).cos())
    I1 = acb.integral(lambda t, d: acb(c**2) * (acb(S2) * t).cos()**2, -0.5, 0.5).real
    I2 = acb.integral(lambda t, d: acb(c) * (acb(S2) * t).cos() * acb(c) * (acb(SA / S2 + CA) - (acb(S2) * t).cos()), -0.5, 0.5).real
    return I1, I2, I1 + I2

# ---------- section 2: certified sup of (H_tau(y) - s)^+ / y over y > 0 ----------
from flint import arb_series

def _g_series(y0, thb, s, n=3, w=None, rho=None):
    """Taylor series in t of g(y0 + t) = (H_tau(y0+t) - s)/(y0+t + rho (1 - H_tau(y0+t))); rho = 0 gives
    (H - s)/y. y0 may be a ball (then the coefficients enclose the derivatives over the ball)."""
    t = arb_series([y0, 1], prec=n)
    def cosh_sinh(z):
        e = z.exp(); ei = (-z).exp()
        return (e + ei) * arb('0.5'), (e - ei) * arb('0.5')
    if w is None:
        cw, ww = COT, S2
    else:
        cw, ww = (w / 2).cos() / (w / 2).sin(), w
    def K(z):
        ch, sh = cosh_sinh(PI * z)
        return ww**2 * (ch + (2 * PI / ww) * cw * z * sh) * (ww**2 + 4 * PI**2 * z * z).inv()
    k1 = K(t)
    k2 = K(2 * t)
    Hs = 1 - (1 - thb)**2 * 4 * k1 * k1 * (1 + k2).inv()
    if rho is None:
        return (Hs - s) * t.inv()
    return (Hs - s) * (t + rho * (1 - Hs)).inv()

def alpha_tilde(thb, s, ylo=arb('0.1'), yhi=None, tol=1e-11, maxit=400000, w=None, rho=None):
    """Upper bound for sup_{y>0} max(H_tau(y) - s, 0)/y (DERIVATION 3.6).
    (0, ylo]: H_tau < 0 <= s by K0(iy) >= 1, K0(2iy) <= cosh(2 pi y) (paper l.802-807).
    [yhi, inf): (H - s)^+/y < (1 - s)^+/yhi since H < 1.
    [ylo, yhi]: branch and bound with second-order Taylor models (Lagrange remainder enclosed by Arb)."""
    tail0 = 1 - (1 - thb)**2 * 4 / (1 + (2 * PI * ylo).cosh())
    assert hi(tail0) < 0, "small-y region not certified"
    if yhi is None:
        yhi = arb(60)
    tailinf = max(fhi(1 - s), 0.0) / float(yhi.mid()) * (1 + 1e-12)
    best_lo = 0.0
    work = [(float(ylo.mid()), float(yhi.mid()))]
    ubs = [tailinf]
    it = 0
    while work:
        it += 1
        if it > maxit:
            raise RuntimeError("B&B did not converge")
        a, b = work.pop()
        m = 0.5 * (a + b); r = max(m - a, b - m) * (1 + 1e-15)
        gp = _g_series(arb(m), thb, s, 2, w, rho)
        g0, g1 = gp[0], gp[1]
        if flo(g0) > best_lo:
            best_lo = flo(g0)
        gb = _g_series(arb.union(arb(a), arb(b)), thb, s, 3, w, rho)
        c2 = gb[2]
        ub = fhi(g0) + fhi(abs(g1)) * r + max(fhi(c2), 0.0) * r * r
        ub *= (1 + 1e-14)
        if ub <= 0:
            continue
        if ub <= best_lo + tol:
            ubs.append(ub); continue
        if r < 1e-9:
            ubs.append(ub); continue
        work.append((a, m)); work.append((m, b))
    return max(max(ubs), 0.0), max(best_lo, 0.0), it

# ---------- section 3: the off-line moment constant D (Littlewood + Conrey Thm 2) ----------
def cstar(R, theta=THETA):
    """inf over P of e^{-2R} c_1(P,1,R,theta) = e^{-2R} + (1-e^{-2R})/2 * (1 + coth(theta R)), DERIVATION 4.3."""
    e = (-2 * R).exp()
    return e + (1 - e) / 2 * (1 + (theta * R).coth())

def Fmom(R, kappa, theta=THETA):
    return R * (1 - kappa) + cstar(R, theta).log()

def D_shift(kappa, theta=THETA):
    """(1/4pi) min_R [R(1-kappa) + log c*(R)]: float search, Arb evaluation at the chosen R (any R is valid)."""
    k = float(kappa.lower().mid()) if isinstance(kappa, arb) else float(kappa)
    f = lambda r: float(Fmom(arb(r), arb(k), theta).mid())
    a, b = 1e-3, 4.0
    g = (math.sqrt(5) - 1) / 2
    for _ in range(200):
        c1, c2 = b - g * (b - a), a + g * (b - a)
        if f(c1) < f(c2): b = c2
        else: a = c1
    R = arb(round(0.5 * (a + b), 6))
    kap = kappa if isinstance(kappa, arb) else arb(kappa)
    val = Fmom(R, kap.lower(), theta) / (4 * PI)        # kappa rounded down -> larger (safe) D
    return val.upper(), R

D_PAPER = arb('0.0806')
D_R0 = ((1 + 1 / THETA).log() / (4 * PI)).upper()       # (1/4pi) log(11/4): the paper's Lemma, re-derived

# ---------- section 4: the proportion bound ----------
def tau_data(tau, w=None):
    tau = arb(tau) if not isinstance(tau, arb) else tau
    k00 = K0_real(tau, w)       # min of K on [0, tau]   (K decreasing on [0,1], DERIVATION 3.1)
    k0 = K0_real(2 * tau, w)    # kappa_0 = min on [0, 2tau]
    k1 = K0_real(4 * tau, w)    # kappa_1 = min on [0, 4tau]
    eps0 = 1 - k0
    betaT = (6 - 16 * k0 / 3 - 8 * k1 / 3) / 3
    lamC = max(fhi(2 * eps0), fhi(betaT), 0.0)
    thb = PI**2 * tau**2 / 2
    assert fhi(4 * tau) <= 1.0                       # K0 decreasing on [0, 4 tau] (DERIVATION 3.1)
    assert flo(2 * k00**2 - 1) >= 0                  # R_2-witness term nonnegative (DERIVATION 3.3)
    assert flo(1 - thb) > 0
    return dict(k00=k00, k0=k0, k1=k1, eps0=eps0, betaT=betaT, lamC=lamC, thb=thb)

_AT_CACHE = {}
def at_cached(tau_str, lam):
    key = (tau_str, round(lam, 15))
    if key not in _AT_CACHE:
        td = tau_data(arb(tau_str))
        _AT_CACHE[key] = alpha_tilde(td['thb'], 3 * arb(lam))
    return _AT_CACHE[key]

def p_from(lam, at_hi, D, C=CMT):
    """p = (2 - lam - C - D * alpha_tilde) / (1 - lam) (DERIVATION 3.7); lower endpoint."""
    lamA = arb(lam)
    return flo((2 - lamA - C - D * arb(at_hi)) / (1 - lamA))

def best_bound(tau_str, Ds, C=CMT, full=True):
    """For each D in Ds: max over lam in a candidate set (all >= lamC) of p. Returns list of (p, lam, at_hi)."""
    td = tau_data(arb(tau_str))
    lamC = td['lamC']
    cands = [lamC * (1 + 1e-12) + 1e-300]
    if full:
        cands += [lamC + (1 - lamC) * t for t in (0.002, 0.005, 0.01, 0.02, 0.04, 0.08)]
    res = []
    for D in Ds:
        best = None
        for lam in cands:
            if lam >= 0.999: continue
            ah, al, it = at_cached(tau_str, lam)
            p = p_from(lam, ah, D, C)
            if best is None or p > best[0]:
                best = (p, lam, ah)
        res.append(best)
    return res, td

# ---------- section 4b: joint use of (E) for the off-line fraction (DERIVATION 3.6) ----------
def D_of_R(R, C, theta=THETA):
    """(1/4pi)[R (C - 1) + log c*(R)]: the moment budget with off-line fraction <= C - 1 (kappa = 2 - C); upper endpoint."""
    return ((R * (C - 1) + cstar(R, theta).log()) / (4 * PI)).upper()

def joint_p(cstr, R, w=None):
    """p = (2 - lam - C - alpha_hat(R) D(R)) / (1 - lam), alpha_hat = sup (H - 3 lam)^+ / (y + (R/4pi)(1 - H))."""
    td = tau_data(arb(cstr), w)
    C = CMT if w is None else C_of(w)
    lam = td['lamC'] * (1 + 1e-12) + 1e-300
    assert lam <= 1.0                                # needed for the off-line budget M_S <= sum K^2 - N - ... (DERIVATION 3.5)
    Ra = arb(R)
    rho = Ra / (4 * PI)
    ah, al, it = alpha_tilde(td['thb'], 3 * arb(lam), w=w, rho=rho)
    D = D_of_R(Ra, C)
    lamA = arb(lam)
    p = (2 - lamA - C - D * arb(ah)) / (1 - lamA)
    return flo(p), ah, fhi(D), lam

def joint_best(cstr, w=None, Rlo=0.6, Rhi=1.8, iters=14):
    g = (math.sqrt(5) - 1) / 2
    a, b = Rlo, Rhi
    cache = {}
    def f(R):
        R = round(R, 6)
        if R not in cache: cache[R] = joint_p(cstr, R, w)
        return cache[R][0]
    for _ in range(iters):
        c1, c2 = b - g * (b - a), a + g * (b - a)
        if f(c1) > f(c2): b = c2
        else: a = c1
    R = round(0.5 * (a + b), 4)
    return (R,) + joint_p(cstr, R, w)

# ---------- section 4c: hinge bound (DERIVATION 3.7): Littlewood at every height rho ----------
def B_hinge(rho, theta=THETA):
    """B(rho) = (1/4pi) log c*(2 pi rho): limsup (1/N) sum_{beta>1/2} (y - rho)^+ <= B(rho), y = (beta-1/2) log T/(2pi)."""
    return (cstar(2 * PI * rho, theta).log() / (4 * PI)).upper()

def H_upper_on(thb, s, a, b, w=None, maxit=200000):
    """Certified upper bound for max_{y in [a,b]} (H_tau(y) - s) (b <= 60), by bisection with Taylor models of H."""
    best_hi = -math.inf
    work = [(a, b)]
    it = 0
    while work:
        it += 1
        if it > maxit: raise RuntimeError("H_upper_on did not converge")
        lo_, hi_ = work.pop()
        m = 0.5 * (lo_ + hi_); r = max(m - lo_, hi_ - m) * (1 + 1e-15)
        hm = H(arb(m), thb, w) - s
        Hb = H(arb.union(arb(lo_), arb(hi_)), thb, w) - s       # interval enclosure of H - s over the box
        ub = fhi(Hb)
        if ub < 0:
            best_hi = max(best_hi, ub); continue
        if fhi(hm) >= 0:
            return fhi(hm)                                    # a point with H - s >= 0: not negative here
        if r < 1e-7:
            return ub
        work.append((lo_, m)); work.append((m, hi_))
    return best_hi

def rho_max(thb, s, w=None):
    """Largest rho (to 1e-4) with H_tau(y) - s < 0 certified on (0, rho]: small-y bound on (0, 0.1], boxes on [0.1, rho]."""
    tail0 = 1 - (1 - thb)**2 * 4 / (1 + (2 * PI * arb('0.1')).cosh())
    assert hi(tail0) < 0
    lo_, hi_ = 0.1, 3.0
    for _ in range(15):
        m = 0.5 * (lo_ + hi_)
        if H_upper_on(thb, s, 0.1, m, w) < 0: lo_ = m
        else: hi_ = m
    rho = math.floor(lo_ * 1e4) / 1e4
    assert H_upper_on(thb, s, 0.1, rho, w) < 0
    return rho

def alpha_rho(thb, s, rho, w=None, yhi=60.0, tol=1e-9, maxit=400000):
    """Upper bound for sup_{y > rho} max(H_tau(y) - s, 0)/(y - rho), given H_tau - s < 0 on (0, rho]."""
    best_lo = 0.0
    ubs = [max(fhi(1 - s), 0.0) / (yhi - rho) * (1 + 1e-12)]   # tail [yhi, inf): H < 1
    work = [(rho, yhi)]
    it = 0
    while work:
        it += 1
        if it > maxit: raise RuntimeError("alpha_rho did not converge")
        a, b = work.pop()
        num = H(arb.union(arb(a), arb(b)), thb, w) - s
        if fhi(num) <= 0:
            continue
        if a <= rho:
            if b - a < 1e-12: raise RuntimeError("H - s >= 0 at rho")
            m = 0.5 * (a + b); work.append((a, m)); work.append((m, b)); continue
        m = 0.5 * (a + b); r = max(m - a, b - m) * (1 + 1e-15)
        # g(y) = (H - s)/(y - rho), Taylor model built from H's series
        def gser(y0, n):
            tt = arb_series([y0, 1], prec=n)
            Hs = _H_series(y0, thb, n, w)
            return (Hs - s) * (tt - rho).inv()
        g0s = gser(arb(m), 2)
        g0, g1 = g0s[0], g0s[1]
        if flo(g0) > best_lo: best_lo = flo(g0)
        c2 = gser(arb.union(arb(a), arb(b)), 3)[2]
        ub = (fhi(g0) + fhi(abs(g1)) * r + max(fhi(c2), 0.0) * r * r) * (1 + 1e-14)
        if ub <= 0: continue
        if ub <= best_lo + tol or r < 1e-9:
            ubs.append(ub); continue
        work.append((a, m)); work.append((m, b))
    return max(max(ubs), 0.0), best_lo, it

def _H_series(y0, thb, n=3, w=None):
    t = arb_series([y0, 1], prec=n)
    if w is None:
        cw, ww = COT, S2
    else:
        cw, ww = (w / 2).cos() / (w / 2).sin(), w
    def cosh_sinh(z):
        e = z.exp(); ei = (-z).exp()
        return (e + ei) * arb('0.5'), (e - ei) * arb('0.5')
    def K(z):
        ch, sh = cosh_sinh(PI * z)
        return ww**2 * (ch + (2 * PI / ww) * cw * z * sh) * (ww**2 + 4 * PI**2 * z * z).inv()
    k1 = K(t); k2 = K(2 * t)
    return 1 - (1 - thb)**2 * 4 * k1 * k1 * (1 + k2).inv()

def hinge_p(cstr, w=None, back=(0.02, 0.05, 0.1)):
    """p = (2 - lam - C - alpha_rho B(rho))/(1 - lam), rho = rho_max - back (best of a few)."""
    td = tau_data(arb(cstr), w)
    C = CMT if w is None else C_of(w)
    lam = td['lamC'] * (1 + 1e-12) + 1e-300
    assert lam < 1.0
    s = 3 * arb(lam)
    r0 = rho_max(td['thb'], s, w)
    best = None
    for bk in back:
        rho = round(r0 - bk, 4)
        if rho <= 0.1: continue
        ar, al, it = alpha_rho(td['thb'], s, rho, w)
        Bv = B_hinge(arb(rho))
        lamA = arb(lam)
        p = flo((2 - lamA - C - Bv * arb(ar)) / (1 - lamA))
        if best is None or p > best[0]:
            best = (p, rho, ar, fhi(Bv), lam, r0)
    return best

def prop21(C=CMT):
    """Lamzouri Prop 2.1 route (Thm 1.1(1.6)): 5/2 - 3/2 C, valid for tau <= tau0."""
    return 5 / arb(2) - 3 * C / 2

def tau0():
    """first positive root of K0(2 tau) = 5/6 (paper l.709-711), by bisection with Arb sign checks."""
    a, b = 0.1, 0.25
    for _ in range(80):
        m = 0.5 * (a + b)
        v = K0_real(arb(2 * m)) - arb(5) / 6
        if flo(v) > 0: a = m
        elif fhi(v) < 0: b = m
        else: break
    return a, b

def main():
    import os, time
    t0 = time.time()
    os.makedirs("logs", exist_ok=True)
    logf = open("logs/cert.log", "w")
    def P(*s):
        line = " ".join(str(x) for x in s)
        print(line); logf.write(line + "\n"); logf.flush(); sys.stdout.flush()
    P("== isozero certificate (python-flint", __import__('flint').__version__, ", prec", ctx.prec, ")")
    P("C_MT =", CMT)
    I1, I2, Cn = check_C_numeric()
    P("check C(f0) by quadrature: int f0^2 =", I1, "; iint|u-v|f0f0 =", I2, "; sum =", Cn, "; overlaps C_MT:", Cn.overlaps(CMT))
    P("5/2 - 3/2 C_MT =", prop21(), " (paper 0.508751)")
    a, b = tau0()
    P("tau0 in [%.9f, %.9f] (paper 0.169641...)" % (a, b))
    P("H0(2) =", H(arb(2), arb(0)), " (paper: > 0.4092)")
    ah, al, it = alpha_tilde(arb(0), arb(0))
    P("alpha_0 in [%.12f, %.12f]  (B&B boxes %d; paper: < 0.21462)" % (al, ah, it))
    P("paper l.816: 2 - C_MT - 0.0806 alpha_0 >= %.9f (paper: > 0.655202)" % flo(2 - CMT - D_PAPER * arb(ah)))
    P("D_R0 = log(11/4)/(4 pi) <= %.12f (paper's Lem:Conrey constant 0.0806, rounded up)" % fhi(D_R0))
    # ----- off-line moments -----
    P("== first moment D(kappa) = (1/4pi) min_R [R(1-kappa) + log c*(R)] (DERIVATION 3.4-3.5), theta -> 4/7")
    for name, kk in (("kappa = 0", arb(0)), ("kappa = 5/12 (PRZZ 2020)", arb(5) / 12), ("kappa = 2 - C_MT (from (E))", 2 - CMT)):
        Dk, Rk = D_shift(kk)
        [(pk, lk, ak)], _ = best_bound('0.001', [Dk], full=False)
        P("  %-30s R = %.6f  D <= %.9f   p(c=0.001) >= %.9f" % (name, float(Rk.mid()), fhi(Dk), pk))
    Rj, pj, ahj, Dj, lj = joint_best('0.001')
    P("  joint (DERIVATION 3.6), c = 0.001: R = %.4f, D(R) <= %.9f, alpha_hat <= %.9f, p >= %.9f" % (Rj, Dj, ahj, pj))
    P("== hinge moments B(rho) = (1/4pi) log c*(2 pi rho) >= limsup (1/N) sum_{beta>1/2} (y - rho)^+ (DERIVATION 3.7)")
    for rho in ('0.25', '0.5', '1', '1.2', '1.5', '2'):
        P("  rho = %-4s  B(rho) <= %.6e" % (rho, fhi(B_hinge(arb(rho)))))
    # ----- table -----
    cs = ['0.001', '0.002', '0.003', '0.005', '0.0075', '0.01', '0.015', '0.02', '0.03', '0.05', '0.075', '0.1', '0.125', '0.15', '0.16', '0.1666666666666666666666666667']
    P("== table: liminf B_0(T,c)/N(T) >= ;  kernel f_0 (Montgomery-Taylor)")
    P("   A: paper inputs (D = 0.0806); B: D = log(11/4)/4pi; C: joint (R = 1.19); H: hinge (new, headline)")
    P("   c        eps0       lamC      | A        | B        | C        | H        rho     alpha_rho  B(rho)    | Pro:IsolatedZeros")
    for c in cs:
        res, td = best_bound(c, [D_PAPER, D_R0], full=False)
        (pA, lA, aA), (pB, lB, aB) = res
        pC = joint_p(c, 1.19)[0]
        pH, rho, ar, Bv, lam, r0 = hinge_p(c)
        P("   %-8s %.3e  %.3e | %.6f | %.6f | %.6f | %.6f %.4f %.6f %.3e | %.6f" % (c[:8], fhi(td['eps0']), td['lamC'], pA, pB, pC, pH, rho, ar, Bv, flo(prop21())))
    # ----- explicit c0 -----
    def cmax(pfun, target, lo_=0.0005, hi_=0.12, steps=20):
        for _ in range(steps):
            m = 0.5 * (lo_ + hi_)
            if pfun('%.10f' % m) >= target: lo_ = m
            else: hi_ = m
        c0 = math.floor(lo_ * 1e5) / 1e5
        return c0, pfun('%.5f' % c0)
    def pA_(cstr):
        [(p, l, a_)], _ = best_bound(cstr, [D_PAPER], full=False); return p
    def pH_(cstr):
        return hinge_p(cstr, back=(0.02,))[0]
    P("== explicit c0: largest c (to 1e-5) with certified p >= target")
    for name, pf, targets, steps in (("A (paper inputs)", pA_, (0.655, 0.6552), 20), ("H (hinge, new)", pH_, (0.6552, 0.66, 0.665, 0.67, 0.672, 0.6724), 14)):
        for target in targets:
            c0, p0 = cmax(pf, target, steps=steps)
            P("  %-17s target %.4f: c0 = %.5f, p(c0) >= %.6f" % (name, target, c0, p0))
    # ----- kernel optimisation inside the family f_w = cos(w x)/Z_w (DERIVATION 5) -----
    P("== kernel f_w = cos(w x)/Z_w: per-c optimisation of w^2 (golden section), certified at the chosen w")
    def pA_w(cstr, w2):
        w = arb(w2).sqrt()
        td = tau_data(arb(cstr), w)
        lam = td['lamC'] * (1 + 1e-12) + 1e-300
        ah_, al_, it_ = alpha_tilde(td['thb'], 3 * arb(lam), w=w)
        lamA = arb(lam)
        return flo((2 - lamA - C_of(w) - D_PAPER * arb(ah_)) / (1 - lamA))
    def pH_w(cstr, w2):
        return hinge_p(cstr, arb(w2).sqrt(), back=(0.02,))[0]
    for cstr in ('0.001', '0.01', '0.05', '0.1', '0.1666666666666666666666666667'):
        for name, f in (("A (paper inputs)", lambda w2: pA_w(cstr, w2)), ("H (hinge, new)", lambda w2: pH_w(cstr, w2))):
            a_, b_ = 1.3, 6.0
            g = (math.sqrt(5) - 1) / 2
            for _ in range(12 if name.startswith("H") else 18):
                c1, c2 = b_ - g * (b_ - a_), a_ + g * (b_ - a_)
                if f(c1) > f(c2): b_ = c2
                else: a_ = c1
            w2 = round(0.5 * (a_ + b_), 4)
            pw, pm = f(w2), f(2.0)
            P("  %-17s c=%-8s w^2=%.4f: C(f_w)<=%.9f p>=%.6f   (f_0: p>=%.6f, gain %.2e)" % (name, cstr[:8], w2, fhi(C_of(arb(w2).sqrt())), pw, pm, pw - pm))
    P("elapsed %.1f s" % (time.time() - t0))
    logf.close()
    return 0

if __name__ == "__main__":
    sys.exit(main())
