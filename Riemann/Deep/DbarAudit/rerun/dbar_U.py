# explore59 dbar_U.py : band-limited zero-interpolation cost U_Z(L) for the continuous Nyman-Beurling problem
# V_L = {E(t)=int_[0,L] e^{-iut} dmu(u)} in L^2(w dt), w=1/(2pi(1/4+t^2)); <e^{-iut},e^{-ivt}> = e^{-|u-v|/2}.
# K(z,w) = (1/2a)[(a^2+wb^2)(1-e^{-iL(z-wb)})/(i(z-wb)) + (a-i wb) + (a+i wb) e^{-iL(z-wb)}], a=1/2, wb=conj(w).
import sys, numpy as np, time
a = 0.5
def Kmat(Z, W, L):
    Z = np.asarray(Z, complex)[:, None]; Wb = np.conj(np.asarray(W, complex))[None, :]
    D = Z - Wb
    small = np.abs(D) < 1e-12
    Ds = np.where(small, 1.0, D)
    ph = np.exp(-1j * L * D)
    I = np.where(small, L, (1 - ph) / (1j * Ds))
    return ((a * a + Wb * Wb) * I + (a - 1j * Wb) + (a + 1j * Wb) * ph) / (2 * a)
def Kblock(Z, W, L, B=200):
    Z = np.asarray(Z, complex); out = np.empty((len(Z), len(W)), complex)
    for i in range(0, len(Z), B):
        out[i:i + B] = Kmat(Z[i:i + B], W, L)
    return out
def U(Zs, L, ret=False):
    G = Kblock(Zs, Zs, L)
    one = np.ones(len(Zs), complex)
    c = np.linalg.solve(G, one)
    u = float(np.real(np.vdot(one, c)))
    if ret: return u, c, np.linalg.cond(G)
    return u
zz = np.loadtxt('/Users/adriansager/playground/riemann/Riemann/Deep/WangAudit/zeros1250.txt')[:, 1]
mode = (sys.argv[1] if len(sys.argv) > 1 else 'all') if __name__ == '__main__' else 'none'
t0 = time.time()
if mode in ('check', 'all'):
    # reproducing-property check: <k_w1,k_w2>_{L2(w)} by quadrature vs K(w2,w1)
    L = 3.0
    w1, w2 = 14.134725, 21.02204 - 0.3j
    t = np.linspace(-4000, 4000, 1600001)
    k1 = Kmat(t, [w1], L)[:, 0]; k2 = Kmat(t, [w2], L)[:, 0]
    wt = 1 / (2 * np.pi * (0.25 + t * t))
    ip = np.trapezoid(k1 * np.conj(k2) * wt, t)
    print('CHECK <k_w1,k_w2> quad =', ip, ' K(w2,w1) =', Kmat([w2], [w1], L)[0, 0], ' (tail ~1/4000)')
    ip2 = np.trapezoid(np.abs(k2) ** 2 * wt, t)
    print('CHECK ||k_w2||^2 quad =', ip2, ' K(w2,w2) =', Kmat([w2], [w2], L)[0, 0])
if mode in ('U', 'all'):
    Csum = lambda g: float(np.sum(2.0 / (0.25 + g ** 2)))
    print('C = 2+gamma-log(4pi) =', 2 + 0.5772156649015329 - np.log(4 * np.pi))
    for K in (100, 400, 1000):
        g = zz[:K]; Zs = np.concatenate([g, -g])
        CK = Csum(g)
        for L in (1.0, 2.0, 3.912, 5.0, 7.155, 10.0, 15.0, 20.0, 40.0):
            u, c, cond = U(Zs, L, True)
            diag = float(np.sum(1.0 / ((0.25 + Zs ** 2) * L + 1)))
            burn = float(np.sum(1.0 / ((L + np.log(np.abs(Zs) / (2 * np.pi))) * (0.25 + Zs ** 2))))
            print(f'K={K:5d} L={L:6.3f}  L*U={L*u:.7f}  L*diag={L*diag:.7f}  L*Bdiag={L*burn:.7f}  C_K={CK:.7f}  cond={cond:.2e}')
print('time', round(time.time() - t0, 1), 's')
