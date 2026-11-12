#!/usr/bin/env python3
"""gl2ceil certificate: rigorous upper bound for the one-piece Levinson-Conrey ceiling
   kappa*(xi) = sup_{R>0, P, Q} 1 - log c(P,Q,R,xi)/R      (DERIVATION.md of Degree2Mol, s2-3; REPORT.md here)
Reduction: c = (1+w(1)^2)/2 + alpha*W0 + gamma*W1, w = e^{Ry}Q, alpha = C_P/xi, gamma = xi*B_P; inf over P = inf over the sinh family
P_lam (lam >= 0); for fixed (R, lam) the inner minimum M(R,lam) over Q (Q(0)=1, Q(y)+Q(1-y)=beta) is bounded below by the Fenchel dual
   D(sigma, mu; R, alpha, gamma) = 1/2 + M - 2 sigma(0) - (2 sigma(1) + M e^{-R})^2/2 - (1/alpha) int (sigma' + mu e^{-Ry})^2 - (1/gamma) int sigma^2
for ANY continuous piecewise-polynomial sigma on [0,1] and ANY symmetric (mu(y) = mu(1-y)) piecewise-polynomial mu, M = int mu.
Cells (R, lam) in [Rc-d, Rc+d] x [lam_a, lam_b]: fixed (sigma, mu) built from the float optimum at the cell centre; monotone bound in lam
(C_lam increasing, B_lam decreasing), second-order Taylor bound in R with explicit remainder; everything in Arb (python-flint).
Tails: R <= R_tail (closed form at R=0 + slope bound), R >= 28 and lam >= LAM (derivative-free bound 'crude').
Run: uv run --no-project --with numpy --with scipy --with python-flint python3 ceil.py <xi_num> <xi_den> <kappa_up> [Rmax] [LAM]
"""
import numpy as np, sys, time, resource, json, os
from scipy.integrate import solve_ivp
from numpy.polynomial import chebyshev as Ch
from numpy.polynomial import legendre as Lg
from flint import arb, arb_poly, ctx, fmpq
ctx.prec = 256
T0 = time.time()
def memMB(): return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/2**20   # macOS: bytes
# ---------------- parameters ----------------
XN, XD = int(sys.argv[1]), int(sys.argv[2]); XIq = fmpq(XN, XD); XI = arb(XIq); xi = XN/XD
KUP = fmpq(*map(int, sys.argv[3].split('/'))) if '/' in sys.argv[3] else None
if KUP is None:
    from fractions import Fraction; f = Fraction(sys.argv[3]); KUP = fmpq(f.numerator, f.denominator)
RMAX = float(sys.argv[4]) if len(sys.argv) > 4 else 28.0
LAM = float(sys.argv[5]) if len(sys.argv) > 5 else 3.0
ONEmK = arb(1) - arb(KUP)
HL = 2.0**-int(os.environ.get('HLBITS', '16'))   # lam grid step (dyadic)
KP, DEG = 32, 12                  # pieces, degree
LOG = []
def log(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); LOG.append(s)
# ---------------- sinh family ----------------
def PAB_f(lam):
    if lam < 1e-3: return 1 + lam**4/45, 1/3 - 2*lam**2/45
    s = np.sinh(lam); return lam**2*(0.5 + np.sinh(2*lam)/(4*lam))/s**2, (np.sinh(2*lam)/(4*lam) - 0.5)/s**2
CBc = {}
def CB_arb(i):   # exact-ish balls of C_lam, B_lam at lam = i*HL
    if i in CBc: return CBc[i]
    if i == 0: v = (arb(1), arb(1)/3)
    else:
        l = arb(i)*arb(HL); s = l.sinh(); s2 = (2*l).sinh()/(4*l)
        v = (l*l*(arb(1)/2 + s2)/(s*s), (s2 - arb(1)/2)/(s*s))
    CBc[i] = v; return v
# ---------------- float primal optimum (y-form EL ODE) and dual functions ----------------
NG = 300
xg_, wg_ = Lg.leggauss(NG); YQ = 0.5 + (xg_+1)/4; WQ = wg_/4
XC = (np.cos(np.pi*(np.arange(DEG+1)+0.5)/(DEG+1)) + 1)/2
def build(R, lam):
    Cl, Bl = PAB_f(lam); al = Cl/xi; ga = xi*Bl; w2 = Cl/(Bl*xi*xi); k = R*R - w2
    ys_s = np.concatenate([(kk + XC)/KP for kk in range(KP)])
    ys_m = np.concatenate([0.5 + (j + XC)/KP for j in range(KP//2)])
    yref = np.where(ys_s >= 0.5, ys_s, 1 - ys_s)
    ys = np.unique(np.concatenate([[0.5, 1.0], YQ, yref, ys_m]))
    def f(y, Y):
        th = np.tanh(R*(2*y-1)); g = 1/(1+np.exp(2*R*(2*y-1)))
        return [Y[1], k*g - 2*R*th*Y[1] - k*Y[0], Y[3], -2*R*th*Y[3] - k*Y[2]]
    sol = solve_ivp(f, (0.5, 1.0), [0.5, 0.0, 0.0, 1.0], method='DOP853', rtol=1e-13, atol=1e-14, t_eval=ys)
    up, dup, ph, dph = sol.y
    a = (1 - up[-1])/ph[-1]; u = up + a*ph; du = dup + a*dph; v = -ph/ph[-1]; dv = -dph/ph[-1]
    i = np.searchsorted(ys, YQ); uq, duq, vq, dvq = u[i], du[i], v[i], dv[i]
    e1 = np.exp(2*R*YQ); e2 = np.exp(2*R*(1-YQ))
    def J(beta):
        Q = beta*uq + vq; dQ = beta*duq + dvq
        W0 = np.sum(WQ*(e1*Q**2 + e2*(beta-Q)**2)); W1 = np.sum(WQ*(e1*(R*Q+dQ)**2 + e2*(R*(beta-Q)+dQ)**2))
        return 0.5*(1 + np.exp(2*R)*(beta-1)**2) + al*W0 + ga*W1
    j0, j1, j2 = J(0.0), J(1.0), J(2.0); A2 = (j2 - 2*j1 + j0)/2; A1 = (j1 - j0) - A2; beta = -A1/(2*A2); Jb = J(beta)
    S = beta*u + v; dS = beta*du + dv
    g = 1/(1+np.exp(2*R*(2*ys-1))); d2S = k*(beta*g - S) - 2*R*np.tanh(R*(2*ys-1))*dS
    def Qd(y):
        up_ = y >= 0.5; j = np.searchsorted(ys, np.where(up_, y, 1-y))
        return np.where(up_, S[j], beta - S[j]), dS[j], np.where(up_, d2S[j], -d2S[j])
    Q0, Q1, Q2 = Qd(ys_s); sig = ga*np.exp(R*ys_s)*(R*Q0 + Q1)
    Q0, Q1, Q2 = Qd(ys_m); mm = np.exp(2*R*ys_m)*(al*Q0 - ga*(R*R*Q0 + 2*R*Q1 + Q2))
    Psub = np.polynomial.Polynomial([-1.0, 2.0])
    def fit(vals, n):
        out = []
        for ii in range(n):
            cc = Ch.chebfit(2*XC-1, vals[ii*(DEG+1):(ii+1)*(DEG+1)], DEG)
            out.append(np.polynomial.Polynomial(Ch.cheb2poly(cc))(Psub).coef)
        return out
    return fit(sig, KP), fit(mm, KP//2), beta, Jb
# ---------------- rigorous dual quantities ----------------
def moments(a, n):      # I_j(a) = int_0^1 x^j e^{a x} dx = 1F1(j+1; j+2; a)/(j+1)
    return [a.hypgeom_1f1(j+1, j+2)/(j+1) for j in range(n)]
I0 = [arb(1)/(j+1) for j in range(4*DEG+8)]
def dot(p, I):
    c = p.coeffs(); s = arb(0)
    for j, cj in enumerate(c): s += cj*I[j]
    return s
XPOLY = arb_poly([0, 1]); ONEMX = arb_poly([1, -1])
def rig(sc, mc, Rc, Rlo):
    """Rc, Rlo: fmpq (dyadic). Returns dict of arb quantities for the dual built from float coefficient lists sc (KP pieces), mc (KP/2)."""
    s = []; prev = None
    for kk in range(KP):
        c = [arb(float(x)) for x in sc[kk]]
        if kk > 0: c[0] = prev
        p = arb_poly(c); s.append(p); prev = p(arb(1))
        assert prev.rad() == 0, "continuity value not exact"
    m = [arb_poly([arb(float(x)) for x in mc[j]]) for j in range(KP//2)]
    mt = [m[kk-KP//2] if kk >= KP//2 else m[KP//2-1-kk](ONEMX) for kk in range(KP)]
    R = arb(Rc); RL = arb(Rlo); n = 2*DEG+4
    I1 = moments(-R/KP, n); I2 = moments(-2*R/KP, n); I3 = moments(-2*RL/KP, n)
    K2 = S1 = M = X = Xp = Xpp = Y = Yp = Ypp = Ylo = arb(0)
    for kk in range(KP):
        sk = s[kk]; dsk = sk.derivative(); mk = mt[kk]
        yk = (arb(kk) + XPOLY)*(arb(1)/KP); yk2 = yk*yk
        e1 = (-R*kk/KP).exp(); e2 = (-2*R*kk/KP).exp(); e3 = (-2*RL*kk/KP).exp()
        K2 += dot(sk*sk, I0)/KP; S1 += dot(dsk*dsk, I0)*KP; M += dot(mk, I0)/KP
        sm = dsk*mk; m2 = mk*mk
        X += e1*dot(sm, I1); Xp -= e1*dot(yk*sm, I1); Xpp += e1*dot(yk2*sm, I1)
        Y += e2*dot(m2, I2)/KP; Yp -= 2*e2*dot(yk*m2, I2)/KP; Ypp += 4*e2*dot(yk2*m2, I2)/KP; Ylo += e3*dot(m2, I3)/KP
    s0 = s[0](arb(0)); s1 = s[-1](arb(1))
    eR = (-R).exp(); eL = (-RL).exp()
    t = 2*s1 + M*eR
    K0 = arb(1)/2 + M - 2*s0 - t*t/2; K0p = t*M*eR
    F0 = M*M*eL*eL + (abs(2*s1) + abs(M)*eL)*abs(M)*eL
    K1 = S1 + 2*X + Y; K1p = 2*Xp + Yp
    F1 = 2*Ylo + 2*(S1.sqrt() + Ylo.sqrt())*Ylo.sqrt()          # sup |K1''| on [Rlo, oo)
    K1pp = 2*Xpp + Ypp                                              # K1''(Rc)
    G1 = 2*S1.sqrt()*Ylo.sqrt() + 8*Ylo                            # sup |K1'''| on [Rlo, oo)
    return dict(K0=K0, K0p=K0p, F0=F0, K1=K1, K1p=K1p, F1=F1, K2=K2, K1pp=K1pp, G1=G1)
# ---------------- cell checks ----------------
def target(Rq):   # exp(R (1 - kappa_up))
    return (arb(Rq)*ONEmK).exp()
def check_lam(q, Rc, d, ia, ib):
    """True iff for every lam-subcell of [ia*HL, ib*HL] (recursively bisected down to HL) the Taylor lower bound beats the target
    at both R endpoints. Returns (ok, worst_margin_float)."""
    Ra, Rb = Rc - d, Rc + d; da = arb(d)
    tA, tB = target(Ra), target(Rb)
    worst = [1e300]
    def L_end(x, z, sgn):
        # third-order Taylor in R: K0 to 2nd order (|K0''| <= F0), K1 to 3rd order (K1''(Rc) exact, |K1'''| <= G1);
        # the quadratic coefficient a2 = -(x/2) K1''(Rc) is replaced by min(lower(a2), 0) so that lower(R) - target(R) is concave
        a2 = -x*q['K1pp']/2; a2l = a2.lower(); a2c = a2l if a2l < 0 else arb(0)
        # one concave lower bound lower(eps) (same formula at both endpoints): lower - target is concave, so endpoints suffice
        return (q['K0'] + sgn*q['K0p']*da - x*(q['K1'] + sgn*q['K1p']*da) - z*q['K2'] + a2c*da*da - da*da/2*q['F0'] - x*q['G1']*da*da*da/6)
    def rec(i, j):
        Ci, _ = CB_arb(i); _, Bj = CB_arb(j)
        x = XI/Ci; z = 1/(XI*Bj)
        la = L_end(x, z, -1); lb = L_end(x, z, +1)
        mA = la - tA; mB = lb - tB
        if mA > 0 and mB > 0:
            worst[0] = min(worst[0], float((mA/tA).lower()), float((mB/tB).lower())); return True
        if j - i <= 1: return False
        h = (i + j)//2
        return rec(i, h) and rec(h, j)
    return rec(ia, ib), worst[0]
def crude(Rq, Cq):   # 1/2 + e^{2R}K/(e^{2R}+2K), K = (C/xi) e^R gd(R)/(4R); increasing in R and C
    R = arb(Rq); gd = 2*((R/2).tanh()).atan(); K = Cq/XI*R.exp()*gd/(4*R); e = (2*R).exp()
    return arb(1)/2 + e*K/(e + 2*K)
def crude_cell(Ra, Rb, ia):
    C, _ = CB_arb(ia); return bool(crude(Ra, C) > target(Rb))

# ---------------- driver ----------------
def frac_str(q): return f"{q.p}/{q.q}"
def tails():
    """Analytic tails. (a) R <= R_tail: m(0) >= 2 - 4xi/(3xi+1) for every lam (closed form at R=0 with C>=1, CB>=1/4, coth>=1),
    and d log c/dR >= -xi*sqrt(B/C) >= -xi/sqrt(3); so kappa <= 1 + xi/sqrt3 - log(m0)/R <= kappa_up for R <= R_tail.
    (b) R >= RMAX: kappa <= [log(4R xi/gd(RMAX)) + 2e^{-RMAX}]/R, decreasing for R >= e gd(RMAX)/(4 xi)."""
    m0 = arb(2) - 4*XI/(3*XI + 1); s = XI/arb(3).sqrt()
    Rt = m0.log()/(1 + s - arb(KUP))
    # strip lam >= LAM: m(0) >= 2 - 1/(3/4 + C_LAM/(4 xi)), slope <= xi sqrt(B_LAM/C_LAM)
    CL, BL = CB_arb(round(LAM/HL)); m0L = arb(2) - 1/(arb(3)/4 + CL/(4*XI)); sL = XI*(BL/CL).sqrt()
    RtL = m0L.log()/(1 + sL - arb(KUP))
    RM = arb(RMAX); gd = 2*((RM/2).tanh()).atan()
    hR = ((4*RM*XI/gd).log() + 2*(-RM).exp())/RM
    mono = bool(RM > arb(1).exp()*gd/(4*XI))
    return dict(m0=m0.str(12), R_tail=Rt.str(12), R_tail_lower=float(Rt.lower()), m0_strip=m0L.str(12), R_tail_strip_lower=float(RtL.lower()),
                h_RMAX=hR.str(12), h_ok=bool(hR < arb(KUP)), mono_ok=mono)
def run(tmax):
    t0 = time.time()
    tl = tails(); log("tails:", json.dumps(tl))
    R0 = fmpq(int(float(os.environ.get('R0', '0.5'))*64), 64)
    assert tl['R_tail_lower'] >= float(R0.p)/float(R0.q) and tl['h_ok'] and tl['mono_ok'], "tails fail"
    # strip lam >= LAM: crude with C_LAM on [Rs, RMAX], Rs = dyadic <= max(R_tail, R_tail_strip)
    iL = round(LAM/HL); CL, _ = CB_arb(iL)
    Rs = fmpq(int(max(tl['R_tail_lower'], tl['R_tail_strip_lower'])*64), 64)
    st = [(Rs, fmpq(int(RMAX*64)+1, 64))]; nst = 0; strip_ok = True
    while st:
        Ra, Rb = st.pop()
        if bool(crude(Ra, CL) > target(Rb)): nst += 1; continue
        if Rb - Ra < fmpq(1, 2**20): strip_ok = False; log("STRIP FAIL", float(Ra), float(Rb)); break
        Rm = (Ra + Rb)/2; st += [(Ra, Rm), (Rm, Rb)]
    log(f"strip lam >= {LAM}: {'OK' if strip_ok else 'FAIL'} with {nst} R-intervals ({time.time()-t0:.1f}s)")
    # core: R in [1/2, RMAX] x lam in [0, LAM]
    Q = []
    nR = int(np.ceil((RMAX - float(R0.p)/float(R0.q))*2)); nl = int(LAM*4)
    c0, c1 = float(os.environ.get('CORE_R0', '-1')), float(os.environ.get('CORE_R1', '1e9'))
    for a in range(nR):
        Rc = R0 + fmpq(2*a + 1, 4)
        if not (c0 <= float(Rc.p)/float(Rc.q) < c1): continue
        for b in range(nl):
            Q.append((Rc, fmpq(1, 4), b*(iL//nl), (b+1)*(iL//nl)))
    stats = dict(crude=0, dual=0, split=0, builds=0); worst = (1e300, None); fails = []; cert = []
    while Q:
        if time.time() - t0 > tmax: log("TIME LIMIT; remaining cells", len(Q)); break
        if memMB() > 280: log("MEMORY LIMIT", memMB()); break
        Rc, d, ia, ib = Q.pop()
        Ra, Rb = Rc - d, Rc + d
        if crude_cell(Ra, Rb, ia): stats['crude'] += 1; cert.append((frac_str(Rc), frac_str(d), ia, ib, 'crude')); continue
        ok = False
        if d <= fmpq(1, 4):
            lc = (ia + ib)/2*HL; Rcf = float(Rc.p)/float(Rc.q)
            try:
                sc, mc, beta, J = build(Rcf, lc); stats['builds'] += 1
                q = rig(sc, mc, Rc, Ra)
                ok, w = check_lam(q, Rc, d, ia, ib)
            except Exception as e:
                log("build/rig error", Rcf, lc, repr(e)); ok = False
            if ok:
                stats['dual'] += 1; cert.append((frac_str(Rc), frac_str(d), ia, ib, 'dual'))
                if w < worst[0]: worst = (w, (float(Rc.p)/float(Rc.q), float(d.p)/float(d.q), ia*HL, ib*HL))
                continue
        stats['split'] += 1
        lw = (ib - ia)*HL; df = float(d.p)/float(d.q)
        if d < fmpq(1, 2**14) and ib - ia <= 1: fails.append((float(Rc.p)/float(Rc.q), df, ia, ib)); log("CELL FAIL", float(Rc.p)/float(Rc.q), df, ia*HL, ib*HL); break
        if (2*df*13 > lw and d >= fmpq(1, 2**14)) or ib - ia <= 1:
            Q += [(Rc - d/2, d/2, ia, ib), (Rc + d/2, d/2, ia, ib)]
        else:
            im = (ia + ib)//2; Q += [(Rc, d, ia, im), (Rc, d, im, ib)]
        if (stats['split'] % 200) == 0:
            log(f"  progress: {stats} queue={len(Q)} worst={worst[0]:.3e} t={time.time()-t0:.0f}s mem={memMB():.0f}MB")
    done = not Q and not fails and strip_ok
    log(f"RESULT xi={XN}/{XD} kappa_up={frac_str(KUP)}: {'CERTIFIED' if done else 'NOT CERTIFIED'}; {stats}; worst rel margin {worst}; time {time.time()-t0:.0f}s; peak mem {memMB():.0f}MB")
    return done, cert
if __name__ == '__main__':
    tmax = float(os.environ.get('TMAX', '1500'))
    done, cert = run(tmax)
    tag = f"{XN}_{XD}_{sys.argv[3].replace('/','_')}"
    with open(os.path.join(os.path.dirname(os.path.abspath(__file__)), 'logs', f'ceil_{tag}.json'), 'w') as fh:
        json.dump(dict(xi=f"{XN}/{XD}", kappa_up=sys.argv[3], certified=done, RMAX=RMAX, LAM=LAM, KP=KP, DEG=DEG, HL=HL, cells=cert, log=LOG), fh)
