#!/usr/bin/env python3
# qcheck.py (robinchaudit) -- (a) direct quadrature cross-check of fals.py's MGF/Taylor evaluation of F_phi(rho);
# (b) Lemma 13 against the exact D = u/L - log(1 + log(1+u)/L); (c) MAD_k formula against quadrature of |s - k/2| f_S(s).
#   uv run --no-project --with mpmath python3 qcheck.py
import math, time, itertools
import mpmath as mp
import importlib.util
spec = importlib.util.spec_from_file_location('fl', __file__.replace('qcheck.py', 'fals_lib.py')); fl = importlib.util.module_from_spec(spec); spec.loader.exec_module(fl)
mp.mp.dps = 60; t0 = time.time()

def FS(t, k):   # Irwin-Hall CDF
    return sum((-1)**j*math.comb(k, j)*(t - j)**k for j in range(int(mp.floor(t)) + 1))/mp.factorial(k)
def fS(t, k):   # Irwin-Hall density
    return sum((-1)**j*math.comb(k, j)*(t - j)**(k - 1) for j in range(int(mp.floor(t)) + 1))/mp.factorial(k - 1)

def Fphi_quad(s, Lp, lam, k):
    a = mp.mpf(lam)/k; tot = mp.mpc(0)
    for j in range(k):
        f = lambda u: mp.exp(s*u)*(1/u + 1/u**2)*FS((u - Lp)/a, k) if k > 1 else mp.exp(s*u)*(1/u + 1/u**2)*(u - Lp)/a
        tot += mp.quad(f, mp.linspace(Lp + j*a, Lp + (j + 1)*a, 5), method='gauss-legendre')
    return tot + fl.Fy(s, Lp + lam)

print('(a) F_phi: MGF/Taylor (fals.py) vs direct quadrature')
for k, gam, L, beta, fac, cen in ((2, '3000000000001', '111.1', '0.999', math.pi, True), (40, '3000000000001', '16.11809565', '0.5', 2.0, True),
                                  (1, '14.134725', '16.11809565', '0.87338196661320752', 0, False), (5, '1e16', '111.1', '0.87338196661320752', 4.0, True)):
    gam = mp.mpf(gam); L = mp.mpf(L); beta = mp.mpf(beta)
    lam = fac*k/gam if fac else mp.mpf('1e-3'); c = lam/2 if cen else 0; s = beta - 1 + 1j*gam
    A, _ = fl.Fphi(s, L - c, lam, k); Bq = Fphi_quad(s, L - c, lam, k)
    print(f'  k={k} gamma={mp.nstr(gam, 6)} L={mp.nstr(L, 6)} beta={mp.nstr(beta, 5)}: |MGF - quad|/|quad| = {mp.nstr(abs(A - Bq)/abs(Bq), 3)}; {time.time()-t0:.0f}s', flush=True)

print('(b) Lemma 13: max D_true/bound')
mp.mp.dps = 120; worst = (0, None)
for eta, f, L in itertools.product(('1e-12', '1e-8', '1e-3', '0.1', '0.3', '0.5'), ('1', '-1', '0.5', '-0.5', '0.99', '-0.99'), ('20', '23.03', '60', '111.1')):
    eta = mp.mpf(eta); u = eta*mp.mpf(f); L = mp.mpf(L)
    D = u/L - mp.log(1 + mp.log(1 + u)/L)
    b = eta**2/(2*L)*(1 + 2*eta/(3*(1 - eta)) + (1 + eta)**2/L + 2*eta*(1 + eta)**3/(3*L**2*(1 - eta*(1 + eta)/L)))
    if D/b > worst[0]: worst = (D/b, (float(eta), float(f), float(L)))
print(f'  worst ratio {mp.nstr(worst[0], 8)} at (eta, u/eta, L) = {worst[1]}')

print('(c) MAD_k: formula vs quadrature of |s - k/2| f_S(s)')
mp.mp.dps = 80
for k in (3, 12, 40):
    q = sum(mp.quad(lambda t: abs(t - mp.mpf(k)/2)*fS(t, k), [j, mp.mpf(j) + mp.mpf(1)/2, j + 1], method='gauss-legendre') for j in range(k))
    print(f'  k={k}: formula {mp.nstr(fl.MADk(k), 20)}  quad {mp.nstr(q, 20)}; {time.time()-t0:.0f}s', flush=True)
