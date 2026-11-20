# e1phase spec.py — Weyl-law model of the Kuznetsov/spectral-large-sieve route with the Bessel phase in the test function.
# In the regime Kuznetsov argument x << |t| (here x ~ T^0.1, |t| ~ Z ~ T^0.4), J_{2it}(x) ~ (x/2)^{2it}/Gamma(1+2it), so
# |h(t)|/cosh(pi t) ~ |phi~(2it)| |2t|^{-1/2}/sqrt(2 pi) with phi~ the Mellin transform of the test function phi(u) = w(u) exp(i Z sqrt(u)).
# The spectral large sieve majorizes sum_j |h(t_j)| |sum c_m rho_j(m)|^2/cosh by (Weyl density t/6) => mass M(Z) = int |phi~(2it)| |t|^{1/2} dt / 6.
import numpy as np
u = np.linspace(1, 2, 6001); du = u[1]-u[0]; w = np.sin(np.pi*(u-1))**2
for Z in [0, 10, 30, 100, 300, 1000]:
    tmax = Z + 60; t = np.arange(-tmax, tmax, 0.05)
    f = w*np.exp(1j*Z*np.sqrt(u))/u
    ph = np.array([np.sum(f*np.exp(2j*tt*np.log(u)))*du for tt in t])
    a = np.abs(ph); M = np.sum(a*np.sqrt(np.abs(t)))*0.05/6; L1 = np.sum(a)*0.05
    c = np.sum(a*np.abs(t))/np.sum(a)
    print('Z=%5d  int|phi~(2it)|dt=%.4f (x sqrt(1+Z): %.4f)  sup|phi~|=%.4f  mean|t|=%.1f  Weyl mass M=%.4f  M/(1+Z)=%.4f' % (
        Z, L1, L1/np.sqrt(1+Z), a.max(), c, M, M/(1+Z)), flush=True)
