
# Untrusted optimiser: joint (P,Q,R) optimum of Conrey's functional at zeta-normalised length xi.
# Output: params.json with exact decimal strings (R, lam, h_n), used by cert.py.
import numpy as np, sys, json, time
from numpy.polynomial import chebyshev as C
from numpy.polynomial import legendre as Lg
NQ = 800
xg, wg = Lg.leggauss(NQ); yg = (xg + 1)/2; wg = wg/2
def basis(R, M):
    z = 1 - 2*yg; E = np.exp(R*yg)
    Phi = [E]; dPhi = [R*E]
    for n in range(1, M+1):
        c = np.zeros(2*n); c[2*n-1] = 1.0
        T = C.chebval(z, c); dT = C.chebval(z, C.chebder(c))*(-2.0)
        Phi.append(E*(T - 1)); dPhi.append(R*E*(T-1) + E*dT)
    Phi = np.array(Phi); dPhi = np.array(dPhi)
    end = np.array([np.exp(R)] + [-2*np.exp(R)]*M)
    return (Phi*wg) @ Phi.T, (dPhi*wg) @ dPhi.T, end
def solveQ(G0, G1, end, a, b, beta1):
    H = a*G0 + b*G1 + 0.5*np.outer(end, end)
    m = H.shape[0]-1
    if not beta1:
        h = np.linalg.solve(H[1:,1:], -H[1:,0])
    else:  # constraint end.v = 0  (Q(1)=0)
        K = np.zeros((m+1, m+1)); K[:m,:m] = H[1:,1:]; K[:m,m] = end[1:]; K[m,:m] = end[1:]
        rhs = np.concatenate([-H[1:,0], [-end[0]]])
        h = np.linalg.solve(K, rhs)[:m]
    v = np.concatenate([[1.0], h])
    return h, 0.5 + v @ H @ v, v
def PAB(lam):
    s = np.sinh(lam)
    return lam**2*(0.5 + np.sinh(2*lam)/(4*lam))/s**2, (np.sinh(2*lam)/(4*lam) - 0.5)/s**2
def best(R, xi, M, Popt=True, beta1=False):
    G0, G1, end = basis(R, M)
    if not Popt:
        h, c, v = solveQ(G0, G1, end, 1.0/xi, xi/3.0, beta1); return c, h, None
    lam = 0.75; cp = None
    for it in range(500):
        A, B = PAB(lam)
        h, c, v = solveQ(G0, G1, end, A/xi, xi*B, beta1)
        W0 = v@G0@v; W1 = v@G1@v
        lam = xi*np.sqrt(W1/W0)
        A, B = PAB(lam)
        c2 = 0.5*(1 + (end@v)**2) + (A/xi)*W0 + xi*B*W1
        if cp is not None and abs(cp - c2) < 1e-15*c2: break
        cp = c2
    return c2, h, lam
def kap(R, xi, M, **kw):
    c, h, lam = best(R, xi, M, **kw); return 1 - np.log(c)/R, h, lam
def golden(f, a, b, tol=1e-9):
    g = (np.sqrt(5)-1)/2
    x1 = b - g*(b-a); x2 = a + g*(b-a); f1 = f(x1); f2 = f(x2)
    while b - a > tol:
        if f1 > f2: b, x2, f2 = x2, x1, f1; x1 = b - g*(b-a); f1 = f(x1)
        else: a, x1, f1 = x1, x2, f2; x2 = a + g*(b-a); f2 = f(x2)
    return (a+b)/2
from fractions import Fraction as Fr
settings = {'AT21_uncond': Fr(25,256) - Fr(1,10**7), 'RSC': Fr(1,8) - Fr(1,10**7), 'Bernard_uncond': Fr(5,54) - Fr(1,10**7),
            # arXiv 2607.00282 Cor cor:fullresultpgl3per uses kappa(theta) = sup over (P,Q,R) at theta -> 1/6- (PGL3 twists), 1/2- (PGL2 twists)
            'family_PGL3': Fr(1,6) - Fr(1,10**7), 'family_PGL2': Fr(1,2) - Fr(1,10**7), 'family_PGL2_simple': Fr(1,2) - Fr(1,10**7)}
M0 = int(sys.argv[1]) if len(sys.argv) > 1 else 30
out = {}
for name, xiF in settings.items():
    xi = float(xiF)
    M = 1 if name.endswith('_simple') else M0   # simple zeros need linear Q (Heath-Brown; 2607.00282 RatioGL3.tex l.185)
    rows = {}
    for lab, kw, lo, hi in [('Px_beta1', dict(Popt=False, beta1=True), 0.3, 20), ('Px_free', dict(Popt=False, beta1=False), 0.3, 20),
                             ('Popt_beta1', dict(Popt=True, beta1=True), 0.3, 20), ('joint', dict(Popt=True, beta1=False), 0.3, 25)]:
        Rb = golden(lambda R: kap(R, xi, M, **kw)[0], lo, hi)  # golden() maximises
        k, h, lam = kap(Rb, xi, M, **kw)
        rows[lab] = (Rb, k, lam)
        print(f"{name:15s} xi={xi:.7f} {lab:11s} R*={Rb:.6f} kappa={k:.9f} lam={lam}", flush=True)
    # paper's c0=c1 construction: P=x, beta=1, R = sqrt(3/5)/xi
    Rp = np.sqrt(0.6)/xi; kp = kap(Rp, xi, M, Popt=False, beta1=True)[0]
    if name in ('RSC', 'Bernard_uncond'):
        xiP = {'RSC': Fr(1,8), 'Bernard_uncond': Fr(5,54)}[name]; RpF = Fr(round(np.sqrt(0.6)/float(xiP)*10**10), 10**10)
        cP, hP, _ = best(float(RpF), float(xiP), M, Popt=False, beta1=True)
        # enforce Q(1)=0 exactly: adjust last coefficient so that sum h = 1/2 (done exactly in cert.py)
        out['paper_'+name] = dict(xi=f"{xiP.numerator}/{xiP.denominator}", R=str(RpF), lam=None, h=[repr(float(x)) for x in hP], kappa_float=float(1-np.log(cP)/float(RpF)), beta1=True)
    print(f"{name:15s} paper c0=c1 (P=x, beta=1, R=sqrt(3/5)/xi={Rp:.5f}): kappa={kp:.6f}", flush=True)
    Rb, k, lam = rows['joint']
    Rr = Fr(round(Rb*10**8), 10**8); lr = Fr(round(lam*10**8), 10**8)
    kk, h, lam2 = kap(float(Rr), xi, M)
    G0, G1, end = basis(float(Rr), M)
    A, B = PAB(float(lr)); hh, c, v = solveQ(G0, G1, end, A/xi, xi*B, False)
    out[name] = dict(xi=f"{xiF.numerator}/{xiF.denominator}", R=str(Rr), lam=str(lr), h=[repr(float(x)) for x in hh], kappa_float=float(1-np.log(c)/float(Rr)), Q1=float(1-2*sum(hh)))
    print(f"   rounded R={Rr} lam={lr}: kappa_float={1-np.log(c)/float(Rr):.10f} Q(1)={1-2*sum(hh):.3e}", flush=True)
json.dump(out, open('params.json','w'), indent=1)
