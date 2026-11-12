"""One-piece Ritz optimiser, copied from Riemann/Deep/Degree2Mol/opt2.py (basis, solveQ, PAB, best; untrusted, float)."""
import numpy as np
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
def solveQ(G0, G1, end, a, b):
    H = a*G0 + b*G1 + 0.5*np.outer(end, end)
    h = np.linalg.solve(H[1:,1:], -H[1:,0])
    v = np.concatenate([[1.0], h])
    return h, 0.5 + v @ H @ v, v
def PAB(lam):
    s = np.sinh(lam)
    return lam**2*(0.5 + np.sinh(2*lam)/(4*lam))/s**2, (np.sinh(2*lam)/(4*lam) - 0.5)/s**2
def best(R, xi, M=30):
    G0, G1, end = basis(R, M)
    lam = 0.75; cp = None
    for it in range(500):
        A, B = PAB(lam)
        h, c, v = solveQ(G0, G1, end, A/xi, xi*B)
        W0 = v@G0@v; W1 = v@G1@v
        lam = xi*np.sqrt(W1/W0)
        A, B = PAB(lam)
        c2 = 0.5*(1 + (end@v)**2) + (A/xi)*W0 + xi*B*W1
        if cp is not None and abs(cp - c2) < 1e-15*c2: break
        cp = c2
    return c2, h, lam
