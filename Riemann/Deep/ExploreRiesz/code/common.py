# explore41 common: unfolded zeros of zeta, Davenport-Heilbronn (DH), CUE, Poisson; normalized Gram of exponentials on L^2(-pi, pi).
import json, numpy as np, mpmath as mp
PI = np.pi
ROOT = '/Users/adriansager/playground/riemann/Riemann/Deep/ExploreHelio/results/'
def theta(t):
    t = np.asarray(t, float)
    return t/2*np.log(t/(2*PI)) - t/2 - PI/8 + 1/(48*t) + 7/(5760*t**3) + 31/(80640*t**5)
def dtheta(t):  # theta'(t)
    t = np.asarray(t, float)
    return 0.5*np.log(t/(2*PI)) - 1/(48*t**2) - 7/(1920*t**4)
def zeta_unfolded():
    g = np.load(ROOT + 'zeta_zeros.npy')
    lam = theta(g)/PI + 1.5          # lam_n = n - S(gamma_n) (midpoint convention), reference lattice Z
    return g, lam
def vartheta_dh(t):
    return np.array([float(mp.im(mp.loggamma(mp.mpf(0.75) + 0.5j*mp.mpf(x)))) + 0.5*x*np.log(5/PI) for x in np.atleast_1d(t)])
def dvartheta_dh(t):
    return np.array([0.5*np.log(5/PI) + 0.5*float(mp.re(mp.digamma(mp.mpf(0.75) + 0.5j*mp.mpf(x)))) for x in np.atleast_1d(t)])
def dh_unfolded(offline=True):
    d = json.load(open(ROOT + 'dh_zeros.json'))
    on = np.array(d['online'])
    pts = list(vartheta_dh(on)/PI + 0.5 + 0j)   # lam = n + dev, dev ~ 0 at the head
    if offline:
        for b, g in d['offline']:
            x = vartheta_dh([g])[0]/PI + 0.5; L = dvartheta_dh([g])[0]/PI
            dd = (b - 0.5)*L                     # unfolded height of the pair
            pts += [x + 1j*dd, x - 1j*dd]
    pts = np.array(pts); o = np.argsort(pts.real, kind='stable')
    return pts[o], d
def cue(N, rng):
    Z = (rng.standard_normal((N, N)) + 1j*rng.standard_normal((N, N)))/np.sqrt(2)
    Q, R = np.linalg.qr(Z); dg = np.diag(R); Q = Q*(dg/np.abs(dg))
    ang = np.sort(np.angle(np.linalg.eigvals(Q)))
    return (ang + PI)*N/(2*PI)
def sincc(w):
    w = np.asarray(w, complex); out = np.ones_like(w)
    m = np.abs(w) > 1e-12
    out[m] = np.sin(PI*w[m])/(PI*w[m])
    return out
def gram(y):
    # normalized Gram of e^{i y_j u} in L^2(-pi, pi) with (1/2pi) du; y complex allowed
    y = np.asarray(y, complex)
    G = sincc(y[:, None] - np.conj(y)[None, :])
    d = np.sqrt(np.real(np.diag(G)))
    return G/(d[:, None]*d[None, :])
def riesz(y):
    ev = np.linalg.eigvalsh(gram(y))
    return ev[0], ev[-1]
