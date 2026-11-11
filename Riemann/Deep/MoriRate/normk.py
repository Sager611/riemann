# normk.py (morirate audit): ||k||^2/H^2 for the unsmoothed CCM/Mori vector k = E(h)|[-a,a], h = I4 h0 - I0 h4, at large mu
# (input [5, Thm D]: ||k||^2 >= 0.109 H^2 for mu >= 50, used in paper III Prop. 3.3(c) and Thm G). float64 (norms only).
import sys, math, json
import numpy as np
from numpy.polynomial import legendre as Lg
from scipy.linalg import eigh_tridiagonal

def run(mu):
    c = 2*math.pi*mu; lam = math.sqrt(mu); a = math.log(lam)
    kmax = int(1.6*c) + 200
    if kmax % 2: kmax += 1
    ks = np.arange(0, kmax+1, 2, dtype=float)
    D = ks*(ks+1) + c*c*(2*ks*(ks+1) - 1)/((2*ks+3)*(2*ks-1))
    E = c*c*(ks[:-1]+1)*(ks[:-1]+2)/((2*ks[:-1]+3)*np.sqrt((2*ks[:-1]+1)*(2*ks[:-1]+5)))
    w, V = eigh_tridiagonal(D, E, select='i', select_range=(0, 2))
    coefs = {}
    for n, idx in ((0, 0), (4, 2)):
        v = V[:, idx].copy()
        if v[0] < 0: v = -v
        full = np.zeros(kmax+1); full[0::2] = v*np.sqrt(ks + 0.5)
        coefs[n] = (full, v[0]*math.sqrt(2))
    I0 = math.sqrt(lam)*coefs[0][1]; I4 = math.sqrt(lam)*coefs[4][1]
    H2 = I0*I0 + I4*I4
    x24, w24 = np.polynomial.legendre.leggauss(24)
    bps = sorted(set([-a, a] + [a - math.log(m) for m in range(2, int(mu) + 1) if a - math.log(m) > -a + 1e-12]))
    ys, ws = [], []
    for p, q in zip(bps[:-1], bps[1:]):
        npc = max(1, int(math.ceil((q - p)/0.01)))
        for j in range(npc):
            lo = p + (q-p)*j/npc; hi = p + (q-p)*(j+1)/npc
            ys.append((lo+hi)/2 + (hi-lo)/2*x24); ws.append((hi-lo)/2*w24)
    ys = np.concatenate(ys); ws = np.concatenate(ws)
    ey = np.exp(ys)
    kv = np.zeros_like(ys)
    m = 1
    while True:
        u = m*ey/lam
        mask = u < 1
        if not mask.any(): break
        uu = u[mask]
        h = (I4*Lg.legval(uu, coefs[0][0]) - I0*Lg.legval(uu, coefs[4][0]))/math.sqrt(lam)
        kv[mask] += h
        m += 1
    kv *= np.exp(ys/2)
    nk2 = float(np.sum(ws*kv*kv))
    print(json.dumps({'mu': mu, 'c': c, 'nodes': int(len(ys)), 'mmax': m-1, 'H2': H2, 'normk2': nk2, 'ratio': nk2/H2, 'chi_over_c2': [float(w[0]/c**2), float(w[2]/c**2)]}), flush=True)

if __name__ == '__main__':
    for mu in [float(x) for x in sys.argv[1:]]:
        run(mu)
