#!/usr/bin/env python3
"""selbergn4: certified explicit Selberg-method proportion kappa_F of critical-line zeros of
F = sum_{j<=N} c_j eps_j q^{s/2-1/4} L(s,chi_j)  (distinct even primitive characters, one modulus).

Part 1 reproduces Dousselin (arXiv 2311.10285, Theorem 'application') in Arb.
Part 2 evaluates the improved chain of DERIVATION.md (all constants rounded upward, Arb balls).
Run:  uv run --no-project --with python-flint --with mpmath --with numpy python3 selbergn4.py [--search]
"""
import sys, time, math, json
from flint import arb, acb, ctx
import numpy as np
ctx.prec = 128
T0 = time.time()
PI = arb.pi(); IU = acb(0, 1); HALF = arb(1)/2
def log(*a):
    print(*a, flush=True)
def up(x):
    return float(x.upper()) if isinstance(x, arb) else float(abs(x).upper())
def lo(x):
    return float(x.lower())
def ball(a, b):
    return arb.union(arb(a), arb(b))
def sqrt_iv(a, b):
    a = max(a, 0.0); return arb.union(arb(a).sqrt(), arb(b).sqrt())
def lin(a, b, w):
    n = max(1, int(math.ceil((b - a)/w))); return [a + (b - a)*i/n for i in range(n + 1)]

# ---------------------------------------------------------------- Part 0: Euler products
def primes_upto(n):
    s = np.ones(n + 1, dtype=bool); s[:2] = False
    for i in range(2, int(n**0.5) + 1):
        if s[i]: s[i*i::i] = False
    return [int(p) for p in np.nonzero(s)[0]]
YP = 10**6
PR = primes_upto(YP)
def eulerprod(f, ctail):
    P = arb(1)
    for p in PR: P *= 1 + f(p)
    return arb.union(P, P*(arb(ctail)/(YP - 1)).exp())   # tail: f(p) <= ctail/(p-1)^2 for p > 10^6
P1 = eulerprod(lambda p: arb(3*p*p - 3*p + 1)/arb((p - 1)**3*p), "3.00001")
P2 = eulerprod(lambda p: arb(5*p**5 - 6*p**4 + 5*p*p - 4*p + 1)/arb((p - 1)**5*p*(p + 1)), "5.00001")
P1u = arb(P1.upper()); P2u = arb(P2.upper())
# Small primes p <= 10^6 treated with signed local factors (DERIVATION §10): their Euler factors become exactly 1,
# only prod_{p > 10^6} (1 + f1(p)) <= exp(3.00001/(10^6-1)) remains.
P1eff = arb((arb("3.00001")/(YP - 1)).exp().upper())
log("P1 =", P1.str(10), " P2 =", P2.str(10), " P1 upper %.10f" % up(P1))

# ---------------------------------------------------------------- Part 1: Dousselin reproduction
import mpmath as mp
GAM = arb(0.25).gamma()/arb(0.75).gamma()
def rho_f(th):
    return float(mp.findroot(lambda x: -1 + 2*th*x + mp.e**(x*(1 - th))*(-1 + 2*x), 0.7))
def rho2_f(th, v, kap):
    a = mp.sqrt(mp.pi*kap*v); b = mp.gamma(0.25)/mp.gamma(0.75)
    f = lambda x: 2*th*x*(a + b*mp.sqrt(x)) - 2*a - b*mp.sqrt(x) + mp.e**(x*(1 - th))*(2*x*(a + b*mp.sqrt(x)) - 2*a - b*mp.sqrt(x))
    return float(mp.findroot(f, 0.7))
def C6_D(th, kap, x, c):   # Lemma A bound, valid for every c > 0
    th = arb(th); c = arb(c); x = arb(x)
    return 1/(1 - th)*(c.exp() + (c*th).exp())/(2*(PI*kap*c).sqrt())*((x/c).sqrt()*(PI*kap).sqrt() + GAM)
def C4_D(t):
    th = arb(t)
    b1v = max(1, abs(1 - 4*t)) + 3*t*t**0.5; b2v = 1 - t + 3*t*(1 - t)**0.5
    b = (arb(max(1, abs(1 - 4*t))) + 3*th*th.sqrt()) if b1v <= b2v else (1 - th + 3*th*(1 - th).sqrt())
    a = 3*(1 - th).sqrt()
    m = a if float(a.mid()) >= float(b.mid()) else b
    return 1/(3*(1 - th)*PI.sqrt())*(arb(9)/5*(1 + th**arb(2.5)) + 2*m)
def dousselin(th, kap=arb(1)/8, nrec=100, s2fix=False):
    tha = arb(th); r = arb(rho_f(th))
    C5 = 1/(1 - tha)*(r.exp() + (r*tha).exp())/(2*(PI*kap*r).sqrt())*GAM
    c3fac = (1/(8*kap) + 2) if s2fix else (1/(8*kap) + arb(3)/2)
    C3 = c3fac*((r.exp() + (r*tha).exp())/((1 - tha)*(PI*r).sqrt())*GAM)**4*P1u
    C2 = 6*(C3 + 1 + 2*C3.sqrt())/(tha*kap)**2
    C4 = C4_D(th); h = 1/kap/nrec; I0 = arb(0); I1 = arb(0)
    for k in range(nrec):
        xk = h*k; xk1 = h*(k + 1)
        c6 = C6_D(th, kap, xk1, rho2_f(th, float(xk1.mid()), float(kap.mid())))
        c7 = (HALF + 2*kap)*c6**2 + 2*C4*c6*kap.sqrt()
        I0 += c7*h; I1 += c7*(xk1**2 - xk**2)/2
    E = arb(32)/3
    K1 = P1*E/(PI.sqrt()*(1 - tha))
    K2 = kap*P2*I0 + P1*(kap/(HALF - 2*kap) + E**2/(4*(1 - tha)**2*PI) + E/(PI.sqrt()*(1 - tha))*(kap.log() - kap))
    K3 = -E**2*P1/(4*(1 - tha)**2*PI*kap)
    K4 = -kap*P2*I1 + P1*((4*kap - HALF)/(HALF - 2*kap)**2 - E**2/(4*(1 - tha)**2*PI*kap)*(1 + kap.log()) + E/(PI.sqrt()*(1 - tha)*kap))
    return dict(C5=C5, C2=C2, K1=K1, K2=K2, K3=K3, K4=K4)
def kappa_dousselin(N, A, d):
    A = arb(A)
    C1 = 8*d['C5']**2*(d['K1']*A*A.log() + d['K2']*A + d['K3']*A.log() + d['K4'])
    if N == 1: return 2*PI*(1/(2*A) - (C1.sqrt() + d['C2'].sqrt())**2/A**3)
    return 2*PI*(1/(2*A) - 4*N*(C1 + d['C2'])/A**3)
TABLE = [(1, "29056699.107509706", 0.011, 5.45e-8), (2, "212583177.09901848", 0.0016, 7.38e-9), (3, "319102776.4709714", 0.0014, 4.91e-9),
         (4, "425715589.6389222", 0.0013, 3.68e-9), (5, "532459869.61320543", 0.0012, 2.94e-9), (10, "1067086846.4520979", 0.001, 1.46e-9)]

# ---------------------------------------------------------------- Part 2: special functions
JS = 40
def Eser(Y, ymax):
    """E(Y) = Delta_2(Y)/sqrt(Y) = sum_j (-iY)^j/(j!(j^2-1/4)) for a real ball 0 <= Y <= ymax <= 3, with tail bound."""
    Y = acb(Y); s = acb(0); term = acb(1)
    for j in range(JS):
        s += term/(arb(j*j) - arb(1)/4); term = term*(-IU*Y)/(j + 1)
    t = arb(ymax)**JS/arb(JS + 1).gamma()*2/(arb(JS*JS) - arb(1)/4)
    return s + acb(arb(0, t), arb(0, t))
def IEser(Y, ymax):
    """IE(Y) = Im E(Y)/Y = sum_k (-1)^{k+1} Y^{2k}/((2k+1)!((2k+1)^2-1/4))."""
    s = arb(0); term = arb(1)
    for k in range(JS//2):
        j = 2*k + 1; s += (-1)**(k + 1)*term/(arb(j*j) - arb(1)/4); term = term*Y*Y/((j + 1)*(j + 2))
    t = arb(ymax)**JS/arb(JS + 1).gamma()*2/(arb(JS*JS) - arb(1)/4)
    return s + arb(0, t)
W8 = acb.exp_pi_i(acb(arb(1)/4)); PSI_INF = PI.sqrt()*acb.exp_pi_i(acb(arb(-1)/4))
def psi_cf(Y):
    """psi(Y) = int_0^Y e^{-ix} x^{-1/2} dx = sqrt(pi) e^{-i pi/4} erf(e^{i pi/4} sqrt(Y))."""
    return PSI_INF*(W8*acb(arb(Y).sqrt())).erf()
def D2_cf(Y):
    """Delta_2(Y) = -2 sqrt(Y) e^{-iY} - (1+2iY) psi(Y)."""
    Y = arb(Y); return -2*Y.sqrt()*acb.exp(-IU*Y) - (1 + 2*IU*Y)*psi_cf(Y)
C1C = 2*PI.sqrt()*acb.exp_pi_i(acb(arb(-3)/4))    # c1 = 2 sqrt(pi) e^{-3 i pi/4} = lim Delta_2(X)/X
def selfcheck():
    for y in ["0.01", "0.5", "1", "2.5"]:
        a = Eser(arb(y), 2.5)*arb(y).sqrt(); b = D2_cf(arb(y))
        assert up(abs(a - b)) < 1e-30, (y, a, b)
    # Delta_2 against the defining double integral (mpmath quadrature, float check only)
    for y in [0.7, 3.0]:
        ref = mp.quad(lambda x: -2/mp.sqrt(x) + mp.quad(lambda t: (mp.e**(-1j*t) - 1)/t**1.5, [0, x]), [0, y])
        assert abs(complex(ref) - complex(float(D2_cf(arb(y)).real.mid()), float(D2_cf(arb(y)).imag.mid()))) < 1e-8
    # F4 representation against Im(e^{-iXu} G)/X
    X = arb("0.5"); u = arb("0.3"); l1 = 1 - u; G = l1.sqrt()*Eser(X*l1, 1.0)
    d = (acb.exp(-IU*X*u)*G).imag/X; F = (X*u).cos()*(l1*l1.sqrt()*IEser(X*l1, 1.0)) - u*(X*u).sinc()*G.real
    assert up(abs(d - F)) < 1e-30
    log("selfcheck: series = closed form, Delta_2 = definition (quadrature), F4 representation: OK")

# ---------------------------------------------------------------- Part 3: improved constants
def g_funcs(th, b=1.0, NX=128, NU1=40, NU2=600):
    """On X-intervals of [0,b]: g6 = sup_u |sqrt(l1)E(Xl1) - sqrt(l2)E(Xl2)|, g4 = sup_u |F4(X,u)|;
    l1 = 1-u, l2 = (th-u)^+, F4 = Im(e^{-iXu}[...])/X = cos(Xu)(l1^{3/2}IE(Xl1)-l2^{3/2}IE(Xl2)) - u sinc(Xu) Re[...]."""
    ub = [(th*i/NU1, th*(i + 1)/NU1, True) for i in range(NU1)] + [(th + (1 - th)*i/NU2, th + (1 - th)*(i + 1)/NU2, False) for i in range(NU2)]
    xs = [b*k/NX for k in range(NX + 1)]; g6 = []; g4 = []
    for k in range(NX):
        X = ball(xs[k], xs[k + 1]); m6 = 0.0; m4 = 0.0
        for (ua, ubb, low) in ub:
            u = ball(ua, ubb); l1 = ball(max(1 - ubb, 0), 1 - ua); s1 = sqrt_iv(1 - ubb, 1 - ua)
            Gc = s1*Eser(X*l1, xs[k + 1]); Gi = l1*s1*IEser(X*l1, xs[k + 1])
            if low:
                l2 = ball(max(th - ubb, 0), th - ua); s2 = sqrt_iv(th - ubb, th - ua)
                Gc = Gc - s2*Eser(X*l2, xs[k + 1]); Gi = Gi - l2*s2*IEser(X*l2, xs[k + 1])
            m6 = max(m6, up(abs(Gc)))
            Xu = X*u; m4 = max(m4, up(abs(Xu.cos()*Gi - u*Xu.sinc()*Gc.real)))
        g6.append(m6); g4.append(m4)
    return xs, g6, g4
def Gsin(y):
    """int_0^y |sin t| dt (y >= 0 arb), continuous and increasing."""
    n = math.floor(float((y/PI).mid())); return 2*n + 1 - (y - n*PI).cos()
def mavg(xs, kap):
    """upper bound on (1/(2X^2)) int_{X/(2kap)}^{(1/(2kap)+2)X} |sin y| dy on each X-interval; <= 1/(2kap)+1 always."""
    c = 1/(2*kap) + 2; cap = up(1/(2*kap) + 1); out = []
    for k in range(len(xs) - 1):
        if xs[k] == 0: out.append(cap); continue
        v = (Gsin(c*arb(xs[k + 1])) - Gsin(arb(xs[k])/(2*kap)))/(2*arb(xs[k])**2)
        out.append(min(cap, up(v)))
    return out
class LV:  # theta-independent large-v ingredients
    pass
def lv_common():
    L = LV()
    r0 = 0.0
    for pts, f in [(lin(0.0, 1.0, 0.0025), lambda Y: acb(Y.sqrt())*Eser(Y, 1.0) - C1C*acb(Y)), (lin(1.0, 50.0, 0.005), lambda Y: D2_cf(Y) - C1C*acb(Y))]:
        for i in range(len(pts) - 1): r0 = max(r0, up(abs(f(ball(pts[i], pts[i + 1])))))
    L.R0max = max(r0, up(PI.sqrt() + 4/arb(50).sqrt()))          # |R0(Y)| <= sqrt(pi) + 4/sqrt(Y) for Y >= 50
    da = 0.0
    for pts, f in [(lin(0.0, 1.0, 0.0025), lambda Y: acb(Y.sqrt())*Eser(Y, 1.0)), (lin(1.0, 20.0, 0.002), lambda Y: D2_cf(Y)/acb(Y))]:
        for i in range(len(pts) - 1): da = max(da, up(abs(f(ball(pts[i], pts[i + 1])))))
    L.DA = max(da, up(abs(C1C) + arb(L.R0max)/20))                # sup_Y |Delta_2(Y)|/max(1,Y)
    em = 0.0; pts = lin(0.0, 1.0, 0.0025)
    for i in range(len(pts) - 1): em = max(em, up(abs(Eser(ball(pts[i], pts[i + 1]), 1.0))))
    L.Emax = em; L.DAp = max(em, L.DA)                               # Dm(u) <= sqrt(1-u) * DAp  (u >= theta)
    ps = 0.0
    pts = lin(0.0, 50.0, 0.005)
    for i in range(len(pts) - 1): ps = max(ps, up(abs(psi_cf(ball(pts[i], pts[i + 1])))))
    L.Psimax = max(ps, up(PI.sqrt() + 2/arb(50).sqrt()))           # |psi(Y)| <= sqrt(pi) + 2/sqrt(Y)
    # W(s) = int_s^inf |psi(Y)| Y^-2 dY on a grid; |psi(Y)| <= 2 sqrt(Y)
    g = [1e-8*(1e8)**(i/800) for i in range(800)] + lin(1.0, 400.0, 0.02)
    ubp = [min(2*math.sqrt(g[i + 1]), up(abs(psi_cf(ball(g[i], g[i + 1]))))) for i in range(len(g) - 1)]
    W = [0.0]*len(g); W[-1] = up((PI.sqrt() + 2/arb(400).sqrt())/400)
    for i in range(len(g) - 2, -1, -1): W[i] = W[i + 1] + ubp[i]*(1/g[i] - 1/g[i + 1])*(1 + 1e-12)
    L.grid = g; L.W = W
    return L
def Ipsi(L, a, b):
    """sup over l in [a,b] of l * int_l^inf |psi(Y)| Y^-2 dY = sup_l int_1^inf |psi(xl)| x^-2 dx."""
    g = L.grid; W = L.W; m = 0.0
    if a < g[0]: m = max(m, 4*math.sqrt(g[0]) + g[0]*W[0])
    for i in range(len(g) - 1):
        if g[i + 1] < a or g[i] > b: continue
        m = max(m, min(g[i + 1], b)*W[i])
    return m
def lv_theta(L, th, xb=1.0):
    """D1max = sup_u |d(1,u)|, Dmax = sup_{x>=1,u} |d(x,u)|, Dpint = sup_u int_1^inf |d'(x,u)| dx, d(x,u) = [Delta_2(x l1)-Delta_2(x l2)]/x."""
    assert xb == 1.0
    d1 = 0.0; pts = lin(0.0, 1.0, 1/1000)
    for j in range(len(pts) - 1):
        ua, ubb = pts[j], pts[j + 1]
        for (a, b) in ([(ua, ubb)] if (ubb <= th or ua >= th) else [(ua, th), (th, ubb)]):
            v = acb(sqrt_iv(1 - b, 1 - a))*Eser(ball(max(1 - b, 0), 1 - a), 1.0)
            if b <= th: v = v - acb(sqrt_iv(th - b, th - a))*Eser(ball(max(th - b, 0), th - a), 1.0)
            d1 = max(d1, up(abs(v)))
    dB = 0.0; xs = lin(1.0, 10.0, 0.005); NUB = 8
    for i in range(len(xs) - 1):
        x = ball(xs[i], xs[i + 1])
        for j in range(NUB):
            ua, ubb = th*j/NUB, th*(j + 1)/NUB
            dd1 = D2_cf(x*ball(1 - ubb, 1 - ua)); y2hi = (th - ua)*xs[i + 1]
            l2 = ball(max(th - ubb, 0), th - ua)
            dd2 = acb(sqrt_iv(max(th - ubb, 0)*xs[i], y2hi))*Eser(x*l2, 1.0) if y2hi <= 1.0 else D2_cf(x*l2)
            dB = max(dB, up(abs((dd1 - dd2)/acb(x))))
    dB = max(dB, up(abs(C1C)*(1 - arb(th)) + 2*arb(L.R0max)/10))
    Dmax = max(L.DA, dB, d1)
    # Dpint, case u in [0,th]: direct, x in [1,200], tail 2 Psimax/200
    pB = 0.0; xs = lin(1.0, 10.0, 0.01) + lin(10.0, 200.0, 0.05)[1:]
    for j in range(NUB):
        ua, ubb = th*j/NUB, th*(j + 1)/NUB; s = 0.0
        l1 = ball(1 - ubb, 1 - ua); l2 = ball(max(th - ubb, 0), th - ua)
        for i in range(len(xs) - 1):
            x = ball(xs[i], xs[i + 1])
            s += up(abs(psi_cf(x*l1) - psi_cf(x*l2)))*(1/xs[i] - 1/xs[i + 1])*(1 + 1e-12)
        pB = max(pB, s + 2*L.Psimax/200)
    Dpint = max(Ipsi(L, 0.0, 1 - th), pB)
    return d1, Dmax, Dpint
class Consts:
    pass
def improved_consts(th, kap, L, NX=128, NU1=40, NU2=600):
    c = Consts(); c.th = th; c.kap = kap
    tha = arb(th); ka = arb(kap)
    c.C5p = 2/(PI*ka*(1 - tha)).sqrt()
    xs, g6, g4 = g_funcs(th, 1.0, NX, NU1, NU2)
    m = mavg(xs, ka)
    Ssum = arb(0); Sx = arb(0)
    for k in range(NX):
        w = arb(m[k])*arb(g6[k])**2 + 2*arb(g4[k])*arb(g6[k])
        Ssum += w*(arb(xs[k + 1]) - arb(xs[k])); Sx += w*(arb(xs[k + 1])**2 - arb(xs[k])**2)/(2*ka)
    c.Sa = Ssum; c.Sb = -Sx                       # S_small(A) <= Sa*A + Sb
    d1, Dmax, Dpint = lv_theta(L, th)
    c.D1 = arb(d1); c.Dmax = arb(Dmax); c.Dpint = arb(Dpint)
    r = 1/(HALF - 2*ka)
    # L(A) = r[(A kap - 1) D1^2 + 2 A kap Dmax Dpint] + kap r^2 2 Dmax (Dmax + Dpint)
    c.La = r*ka*(c.D1**2 + 2*c.Dmax*c.Dpint); c.Lb = -r*c.D1**2 + ka*r*r*2*c.Dmax*(c.Dmax + c.Dpint)
    c.pref = 8*P1u/(PI**2*ka*(1 - tha)**3)
    c.Ka = c.pref*(c.La + c.Sa); c.Kb = c.pref*(c.Lb + c.Sb)   # C1(A) <= Ka*A + Kb
    c.C3 = (2*ka + 20*ka**2)*c.C5p**4*P1u
    c.C2 = 4*(c.C3.sqrt() + 1)**2/(tha*ka)**2
    c.g6 = (min(g6), max(g6)); c.g4 = (min(g4), max(g4)); c.mint = sum(m[k]*(xs[k+1]-xs[k]) for k in range(NX))
    return c
def kappa_new(N, A, c):
    A = arb(A); C1 = c.Ka*A + c.Kb
    if not (lo(C1) > 0 and lo(A) > 1/c.kap): return arb(0)
    return 2*PI*(1/A - N*(C1.sqrt() + c.C2.sqrt())**2/A**3)
def best_A(N, c):
    Ka, Kb, C2 = float(c.Ka.mid()), float(c.Kb.mid()), float(c.C2.mid())
    f = lambda A: (1/A - N*(math.sqrt(max(Ka*A + Kb, 0)) + math.sqrt(C2))**2/A**3)
    a, b = 1/c.kap + 1, 1e9
    for _ in range(200):   # golden section on log A
        m1 = a*(b/a)**0.382; m2 = a*(b/a)**0.618
        if f(m1) < f(m2): a = m1
        else: b = m2
    return float("%.6g" % math.sqrt(a*b))


# ---------------------------------------------------------------- Part 3b: Phi-refined constants (u + u' >= t weights)
def ubox_list(th, NU1, NU2):
    return [(th*i/NU1, th*(i + 1)/NU1, True) for i in range(NU1)] + [(th + (1 - th)*i/NU2, th + (1 - th)*(i + 1)/NU2, False) for i in range(NU2)]
def Bhat_up(th, ua, ub_, low):
    # Bhat(u) = (sqrt(1-u) - sqrt((th-u)^+))/sqrt(1-th) <= 1; upper bound on the box
    v = arb(1 - ua).sqrt() - (arb(max(th - ub_, 0)).sqrt() if low else 0)
    return min(1.0, up(v/(1 - arb(th)).sqrt()))
def make_SB(ubs, Bh):
    order = np.argsort(ubs); us = np.asarray(ubs)[order]; bs = np.asarray(Bh)[order]
    suff = np.maximum.accumulate(bs[::-1])[::-1]
    suff = np.append(suff, 0.0)
    def SB(r):
        idx = np.searchsorted(us, r, side='left'); return suff[idx]
    return SB
def phi_vec(f, ubs, SB, tg):
    """Phi_f(t_m) = max_i f_i * SB(t_m - u_i^+) on the t-grid (upper bound, decreasing in t)."""
    f = np.asarray(f); ubs = np.asarray(ubs)
    R = tg[:, None] - ubs[None, :]
    return np.max(f[None, :]*SB(R), axis=1)*(1 + 1e-12)
def improved_consts2(th, kap, L, NX=128, NU1=40, NU2=600, NT=200, P1x=None):
    if P1x is None: P1x = P1u
    c = Consts(); c.th = th; c.kap = kap; tha = arb(th); ka = arb(kap); c.P1x = P1x
    c.C5p = 2/(PI*ka*(1 - tha)).sqrt()
    U = ubox_list(th, NU1, NU2); ubs = [b for (_, b, _) in U]
    Bh = [Bhat_up(th, a, b, lw) for (a, b, lw) in U]
    SB = make_SB(ubs, Bh)
    tg = np.array([2*m/NT for m in range(NT + 1)]); dt = 2/NT
    xs = [k/NX for k in range(NX + 1)]
    # small v
    Sa = arb(0); Sb = arb(0); cc = 1/(2*ka)
    for k in range(NX):
        X = ball(xs[k], xs[k + 1]); S6 = []; S4 = []
        for (ua, ubb, low) in U:
            u = ball(ua, ubb); l1 = ball(max(1 - ubb, 0), 1 - ua); s1 = sqrt_iv(1 - ubb, 1 - ua)
            Gc = s1*Eser(X*l1, xs[k + 1]); Gi = l1*s1*IEser(X*l1, xs[k + 1])
            if low:
                l2 = ball(max(th - ubb, 0), th - ua); s2 = sqrt_iv(th - ubb, th - ua)
                Gc = Gc - s2*Eser(X*l2, xs[k + 1]); Gi = Gi - l2*s2*IEser(X*l2, xs[k + 1])
            S6.append(up(abs(Gc))); Xu = X*u; S4.append(up(abs(Xu.cos()*Gi - u*Xu.sinc()*Gc.real)))
        P6 = phi_vec(S6, ubs, SB, tg); P4 = phi_vec(S4, ubs, SB, tg)
        # first term: (1/(2X)) int_0^2 Phi6^2 |sin(X/(2kap) + tX)| dt ; second: int_0^2 Phi6 Phi4 dt
        w1 = 0.0
        for m in range(NT):
            capv = up(cc + arb(tg[m + 1]))            # |sin y|/X <= y/X
            if xs[k] > 0:
                ylo = arb(xs[k])*(cc + arb(tg[m])); yhi = arb(xs[k + 1])*(cc + arb(tg[m + 1]))
                n1 = math.floor(float((ylo/PI - HALF).mid())); n2 = math.floor(float((yhi/PI - HALF).mid()))
                ms = 1.0 if n2 > n1 else max(up(abs(ylo.sin())), up(abs(yhi.sin())))
                capv = min(capv, ms/xs[k])
            w1 += P6[m]**2*capv*dt/2
        w2 = float(np.sum(P6[:-1]*P4[:-1]))*dt
        w = arb(w1*(1 + 1e-10)) + arb(w2*(1 + 1e-10))
        Sa += w*(arb(xs[k + 1]) - arb(xs[k])); Sb -= w*(arb(xs[k + 1])**2 - arb(xs[k])**2)/(2*ka)
    c.Sa = Sa; c.Sb = Sb
    # large v: S1(u) at X = 1, Dm(u), Dp(u)
    S1 = []; Dm = []; Dp = []
    DB = 0.0; xsB = lin(1.0, 10.0, 0.005); NUB = 8
    for i in range(len(xsB) - 1):
        x = ball(xsB[i], xsB[i + 1])
        for j in range(NUB):
            ua, ubb = th*j/NUB, th*(j + 1)/NUB
            dd1 = D2_cf(x*ball(1 - ubb, 1 - ua)); y2hi = (th - ua)*xsB[i + 1]; l2 = ball(max(th - ubb, 0), th - ua)
            dd2 = acb(sqrt_iv(max(th - ubb, 0)*xsB[i], y2hi))*Eser(x*l2, 1.0) if y2hi <= 1.0 else D2_cf(x*l2)
            DB = max(DB, up(abs((dd1 - dd2)/acb(x))))
    DB = max(DB, up(abs(C1C)*(1 - tha) + 2*arb(L.R0max)/10))
    pB = 0.0; xsP = lin(1.0, 10.0, 0.01) + lin(10.0, 200.0, 0.05)[1:]
    for j in range(NUB):
        ua, ubb = th*j/NUB, th*(j + 1)/NUB; s = 0.0; l1 = ball(1 - ubb, 1 - ua); l2 = ball(max(th - ubb, 0), th - ua)
        for i in range(len(xsP) - 1):
            x = ball(xsP[i], xsP[i + 1]); s += up(abs(psi_cf(x*l1) - psi_cf(x*l2)))*(1/xsP[i] - 1/xsP[i + 1])*(1 + 1e-12)
        pB = max(pB, s + 2*L.Psimax/200)
    for (ua, ubb, low) in U:
        v = acb(sqrt_iv(1 - ubb, 1 - ua))*Eser(ball(max(1 - ubb, 0), 1 - ua), 1.0)
        if low: v = v - acb(sqrt_iv(th - ubb, th - ua))*Eser(ball(max(th - ubb, 0), th - ua), 1.0)
        S1.append(up(abs(v)))
        if low: Dm.append(max(DB, S1[-1])); Dp.append(pB)
        else: Dm.append(max(math.sqrt(1 - ua)*L.DAp, S1[-1])); Dp.append(Ipsi(L, max(1 - ubb, 0), 1 - ua))
    PS1 = phi_vec(S1, ubs, SB, tg); PDm = phi_vec(Dm, ubs, SB, tg); PDp = phi_vec(Dp, ubs, SB, tg)
    avg = lambda arr: arb(float(np.sum(arr[:-1]))*dt/2*(1 + 1e-10))
    # log(xi)/Lambda <= kap*r(t), r(t) = 1/(1/2 + (t-2)kap)  (Lambda >= (1/2)log T + (t_d - 2) log xi since delta_1, delta_2 <= xi);
    # r is decreasing, so r(t_m) bounds it on the cell [t_m, t_{m+1}]
    rt = 1/(0.5 + (tg - 2)*kap)*(1 + 1e-12)
    r = 1/(HALF - 2*ka)
    assert abs(rt[0] - float(r.mid())) < 1e-6*float(r.mid()) + 1e-9
    aS1 = avg(rt*PS1**2); aX = avg(rt*2*PDm*PDp); a2 = avg(rt**2*(PS1**2 + PDm**2 + 2*PDm*PDp))
    c.La = ka*(aS1 + aX); c.Lb = -aS1 + ka*a2
    # Note -r*aS1 is fine: the true term is (A kap - 1)*(S1 S1 average) with A kap > 1.
    c.pref = 8*P1x/(PI**2*ka*(1 - tha)**3)
    c.Ka = c.pref*(c.La + c.Sa); c.Kb = c.pref*(c.Lb + c.Sb)
    PBB = phi_vec(Bh, ubs, SB, tg)
    i0 = arb(float(np.sum(PBB[:-1]**2))*dt*(1 + 1e-10))
    i1 = arb(float(np.sum(PBB[:-1]**2*(6*(tg[:-1] + dt/2) + 4)))*dt*(1 + 1e-10))   # (6t+4) on cell: use cell value at right end below
    i1 = arb(float(np.sum(PBB[:-1]**2*(6*tg[1:] + 4)))*dt*(1 + 1e-10))
    c.C3 = c.C5p**4*P1x*(ka*i0 + ka**2*i1)
    # ||L eta^2 - 1||^2 = ||L eta^2||^2 - 2 Re int L eta^2 + T <= (C3 - 1 + o(1)) T, since int_T^{2T} L eta^2 = T + o(T) (DERIVATION §3)
    c.C2 = 4*(c.C3 - 1)/(tha*ka)**2
    c.diag = dict(aS1r=float(aS1.mid()), aXr=float(aX.mid()), DB=DB, pB=pB, i0=float(i0.mid()), i1=float(i1.mid()))
    return c

if __name__ == "__main__":
    KAP = 0.1249999          # kappa < 1/8 (mean value theorem for |I|^2 needs 4 kappa < 1/2)
    log("== Part 1: Dousselin, Theorem 'application', reproduced in Arb (P1, P2 with rigorous tails) ==")
    DT = {}
    for (N, A, th, printed) in TABLE:
        d = dousselin(th); DT[N] = (A, th, d)
        k = kappa_dousselin(N, A, d); k2 = kappa_dousselin(N, A, dousselin(th, s2fix=True))
        log("N=%-3d theta=%-7g A=%-20s kappa >= %.6e (printed %.3g); with Sigma_2' bound doubled: %.6e" % (N, th, A, lo(k), printed, lo(k2)))
        assert abs(lo(k)/printed - 1) < 0.006
    selfcheck()
    L = lv_common()
    log("R0max %.6f  DA %.6f  Emax %.6f  DAp %.6f  Psimax %.6f  Ipsi[0,1] %.6f  (%.0fs)" % (L.R0max, L.DA, L.Emax, L.DAp, L.Psimax, Ipsi(L, 0, 1), time.time() - T0))
    if "--search" in sys.argv:
        for (th, kk) in [(t, k) for t in [0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.5] for k in [0.1249999, 0.11]]:
            c = improved_consts2(th, kk, L, NX=32, NU1=10, NU2=150, NT=100, P1x=P1eff)
            row = [(N, lo(kappa_new(N, best_A(N, c), c))) for N in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 100, 1000]]
            log("th=%.3f kap=%.4f Ka=%.1f Kb=%.1f C2=%.4g C3=%.4g La=%.3f Sa=%.3f | " % (th, kk, float(c.Ka.mid()), float(c.Kb.mid()), float(c.C2.mid()), float(c.C3.mid()), float(c.La.mid()), float(c.Sa.mid())) + " ".join("N%d:%.4g" % r for r in row))
        sys.exit(0)
    log("== Part 2a: Dousselin's constants C1, C2 with the new counting (1/A) and per-j combination ((sqrtC1+sqrtC2)^2) ==")
    for N in [2, 3, 4, 5, 10]:
        A0, th, d = DT[N]; best = (0, 0)
        for Af in [float(A0)*f for f in [0.05, 0.08, 0.1, 0.12, 0.15, 0.2, 0.3]]:
            A = arb(Af); C1 = 8*d['C5']**2*(d['K1']*A*A.log() + d['K2']*A + d['K3']*A.log() + d['K4'])
            k = lo(2*PI*(1/A - N*(C1.sqrt() + d['C2'].sqrt())**2/A**3))
            if k > best[0]: best = (k, Af)
        log("N=%-3d theta=%g: kappa >= %.4e at A=%.4g  (x%.1f over Dousselin)" % (N, th, best[0], best[1], best[0]/lo(kappa_dousselin(N, A0, d))))
    log("== Part 2b: + exact main terms (C5', C6', C4' via Delta_2), local factors cancelled (P1), |sin|-average, integration by parts for v > 1/log xi; sup-norm weights ==")
    cB = improved_consts(0.1, KAP, L, NX=64, NU1=20, NU2=300)
    log("theta=0.1: Ka=%.2f Kb=%.2f La=%.4f Sa=%.4f C3=%.3f C2=%.4e  D1=%.5f Dmax=%.5f Dpint=%.5f" % (up(cB.Ka), up(cB.Kb), up(cB.La), up(cB.Sa), up(cB.C3), up(cB.C2), up(cB.D1), up(cB.Dmax), up(cB.Dpint)))
    for N in [4, 5, 10]:
        A = best_A(N, cB); log("   N=%d: kappa >= %.4e (A=%g)" % (N, lo(kappa_new(N, A, cB)), A))
    log("== Part 2c: + weights restricted to u+u' >= log d/log xi (Lemma Phi), absolute local factors at every prime (Euler product P1) ==")
    c = improved_consts2(0.1, KAP, L, NX=128, NU1=40, NU2=600, NT=200, P1x=P1u)
    log("theta=0.100 kappa=%.7f: Ka=%.3f Kb=%.3f La=%.5f Sa=%.5f C3=%.4f C2=%.6e pref=%.4f (%.0fs)" % (KAP, up(c.Ka), up(c.Kb), up(c.La), up(c.Sa), up(c.C3), up(c.C2), up(c.pref), time.time() - T0))
    for N in [4, 5, 10]:
        A = best_A(N, c); log("   N=%d: kappa >= %.4e (A=%g)" % (N, lo(kappa_new(N, A, c)), A))
    log("== Part 2d: certified chain (DERIVATION.md): + signed local factors at the primes p <= 10^6 (Euler product -> P1eff = %.9f) ==" % up(P1eff))
    CS = {}
    for th in [0.03, 0.05, 0.075, 0.1, 0.125, 0.15, 0.175, 0.2, 0.25]:
        c = improved_consts2(th, KAP, L, NX=128, NU1=40, NU2=600, NT=200, P1x=P1eff); CS[th] = c
        log("theta=%.3f kappa=%.7f: C5'=%.6f Ka=%.4f Kb=%.4f La=%.5f Sa=%.5f C3=%.5f C2=%.6e pref=%.5f %s (%.0fs)" % (th, KAP, up(c.C5p), up(c.Ka), up(c.Kb), up(c.La), up(c.Sa), up(c.C3), up(c.C2), up(c.pref), json.dumps({k: round(v, 5) for k, v in c.diag.items()}), time.time() - T0))
    log("== RESULT: kappa_F >= (N, theta, A, bound, Dousselin's, ratio) ==")
    printed = {1: 5.45e-8, 2: 7.38e-9, 3: 4.91e-9, 4: 3.68e-9, 5: 2.94e-9, 10: 1.46e-9}
    RES = {}
    for N in range(1, 11):
        best = (0, None, None)
        for th, c in CS.items():
            A = best_A(N, c); k = lo(kappa_new(N, A, c))
            if k > best[0]: best = (k, th, A)
        RES[N] = best
        pd = printed.get(N)
        log("N=%-3d theta=%.3f A=%-10g kappa_F >= %.4e   %s" % (N, best[1], best[2], best[0], ("Dousselin %.3g, ratio %.0f" % (pd, best[0]/pd)) if pd else ""))
    for N in [100, 1000]:
        best = max(((lo(kappa_new(N, best_A(N, c), c)), th) for th, c in CS.items()))
        log("N=%-5d kappa_F >= %.4e (theta=%.3f); Dousselin's table: %s" % (N, best[0], best[1], {100: "1.45e-10", 1000: "1.43e-11"}[N]))
    log("== Large N: N*kappa_F >= 2 pi [1/alpha - (sqrt(Ka alpha) + sqrt(C2/N1))^2/alpha^3] for all N >= N1 (increasing in N; Kb <= 0) ==")
    for th, c in CS.items():
        assert up(c.Kb) <= 0
        out = []
        for N1 in [11, 100, 1000, 10**4, 10**6]:
            bestc = 0
            for al in [up(c.Ka)*f for f in [1.6, 1.8, 2.0, 2.2, 2.5, 3.0, 4.0]]:
                a = arb(al); v = lo(2*PI*(1/a - (( c.Ka*a).sqrt() + (c.C2/N1).sqrt())**2/a**3))
                bestc = max(bestc, v)
            out.append("N1=%g: %.4e" % (N1, bestc))
        log("theta=%.3f: c(N1) = " % th + ", ".join(out) + ";  limit pi/(2 Ka) = %.4e" % lo(PI/(2*c.Ka)))
    for N in [4, 10]:
        assert RES[N][0] > 10**5*printed[N]
    json.dump({str(N): {"kappa_lower": RES[N][0], "theta": RES[N][1], "A": RES[N][2]} for N in RES}, open("logs/result.json", "w"), indent=1)
    log("all assertions passed; total %.0fs" % (time.time() - T0))
