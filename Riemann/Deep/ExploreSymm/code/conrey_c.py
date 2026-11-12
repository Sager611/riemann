# Conrey's c(P,Q,R) for P(x)=x, Q(x)=1-b x, theta given; kappa >= 1 - log(c)/R
import numpy as np
from scipy.optimize import minimize
from numpy.polynomial.legendre import leggauss
xg, wg = leggauss(60); u = (xg+1)/2; wu = wg/2
U, Vv = np.meshgrid(u, u, indexing='ij'); W = np.outer(wu, wu)
def c(b, R, th):
    Q = 1 - b*Vv; dQ = -b
    integ = (Q*(1 + R*th*U) + U*th*dQ)**2      # d/dx[e^{R th x}(x+u)Q(v+th x)] at x=0
    return 1 + (1/th)*np.sum(W*np.exp(2*R*Vv)*integ)
for th in (0.5, 4/7):
    f = lambda p: -(1 - np.log(c(p[0], p[1], th))/p[1])
    r = minimize(f, [1.0, 1.3], method='Nelder-Mead', options={'xatol':1e-8,'fatol':1e-10})
    b, R = r.x
    print(f"theta={th:.4f} b={b:.5f} R={R:.5f} kappa={-r.fun:.6f} c={c(b,R,th):.6f} h=0.5logc={0.5*np.log(c(b,R,th)):.6f}")
    for RR in (0.8, 1.0, 1.3, 1.6, 2.0, 2.5, 3.0):
        print(f"   R={RR}: c={c(b,RR,th):.5f}  h/R={0.5*np.log(c(b,RR,th))/RR:.5f}")
