
import numpy as np, sys, time
from numpy.polynomial import chebyshev as C
from numpy.polynomial import legendre as Lg
np.set_printoptions(precision=6, linewidth=150)

# Gauss-Legendre nodes on [0,1]
NQ = 600
xg, wg = Lg.leggauss(NQ); yg = (xg + 1)/2; wg = wg/2

def basis(R, M):
    """phi_0 = e^{Ry}, phi_n = e^{Ry}(T_{2n-1}(1-2y)-1), n=1..M ; values and derivatives at nodes, and at y=1"""
    z = 1 - 2*yg
    E = np.exp(R*yg)
    Phi = [E]; dPhi = [R*E]
    for n in range(1, M+1):
        c = np.zeros(2*n); c[2*n-1] = 1.0
        T = C.chebval(z, c); dT = C.chebval(z, C.chebder(c))*(-2.0)
        Phi.append(E*(T - 1)); dPhi.append(R*E*(T-1) + E*dT)
    Phi = np.array(Phi); dPhi = np.array(dPhi)
    # w(1) coefficients: e^R * (1, T(-1)-1 = -2, ...)
    end = np.array([np.exp(R)] + [-2*np.exp(R)]*M)
    G0 = (Phi*wg) @ Phi.T; G1 = (dPhi*wg) @ dPhi.T
    return G0, G1, end

def solveQ(G0, G1, end, a, b):
    # minimize 1/2 + 1/2 (end.(1,h))^2 + a (1,h)G0(1,h) + b (1,h)G1(1,h)
    H = a*G0 + b*G1 + 0.5*np.outer(end, end)
    Hh = H[1:,1:]; g = H[1:,0]
    h = np.linalg.solve(Hh, -g)
    v = np.concatenate([[1.0], h])
    c = 0.5 + v @ H @ v
    return h, c, v

def PAB(lam):
    # P = sinh(lam x)/sinh(lam): A=int P'^2, B=int P^2
    s = np.sinh(lam)
    A = lam**2*(0.5 + np.sinh(2*lam)/(4*lam))/s**2
    B = (np.sinh(2*lam)/(4*lam) - 0.5)/s**2
    return A, B

def best_c(R, th, M, mode='joint', iters=200):
    G0, G1, end = basis(R, M)
    if mode == 'Px':
        A, B = 1.0, 1.0/3
        h, c, v = solveQ(G0, G1, end, A/th, th*B)
        return c, h, None
    lam = 1.0
    cprev = None
    for it in range(iters):
        A, B = PAB(lam)
        h, c, v = solveQ(G0, G1, end, A/th, th*B)
        W0 = v@G0@v; W1 = v@G1@v
        lam = th*np.sqrt(W1/W0)
        # optimal-P closed form for this Q
        copt = 0.5*(1 + (end@v)**2) + np.sqrt(W0*W1)/np.tanh(lam)
        if cprev is not None and abs(cprev - copt) < 1e-14*copt: break
        cprev = copt
    return copt, h, lam

def kappa(R, th, M, mode='joint'):
    c, h, lam = best_c(R, th, M, mode)
    return 1 - np.log(c)/R, h, lam

def golden(f, a, b, tol=1e-7):
    g = (np.sqrt(5)-1)/2
    x1 = b - g*(b-a); x2 = a + g*(b-a); f1 = f(x1); f2 = f(x2)
    while b - a > tol:
        if f1 > f2: b, x2, f2 = x2, x1, f1; x1 = b - g*(b-a); f1 = f(x1)
        else: a, x1, f1 = x1, x2, f2; x2 = a + g*(b-a); f2 = f(x2)
    return (a+b)/2

mode = sys.argv[1]; th = eval(sys.argv[2]); M = int(sys.argv[3]); Rlo, Rhi = float(sys.argv[4]), float(sys.argv[5])
t0 = time.time()
Rb = golden(lambda R: kappa(R, th, M, mode)[0], Rlo, Rhi, 1e-6)
k, h, lam = kappa(Rb, th, M, mode)
print(f"mode={mode} th={th:.6f} M={M} R*={Rb:.6f} kappa={k:.8f} lam={lam} time={time.time()-t0:.1f}s")
