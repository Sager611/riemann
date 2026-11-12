# Independent lambda_n, n <= NMAX, from Odlyzko's first 100000 zeros (zeros1, |err| <= 3e-9) under RH:
# lambda_n = sum_{gamma<=T} g_n(gamma) + int_T^inf g_n A' dt - g_n(T) S(T) + E_n,
# |E_n| <= max_{t>=T}|S_1(t)| (4n/T^2 + 2n^2/(3T^3))  (E_n = S_1(T) g_n'(T) + int_T^inf S_1 g_n'').
# usage: python3 zeros_lambda.py zeros1 NMAX out.npy
import sys, time, numpy as np, mpmath as mp
t0 = time.time()
g = np.loadtxt(sys.argv[1]); NMAX = int(sys.argv[2]); out = sys.argv[3]
assert g.size == 100000 and np.all(np.diff(g) > 0)
T = 0.5*(g[-1] + g[-2]); N0 = g.size
th = 2*np.arctan(1/(2*g))
C = np.zeros(NMAX+1)
z = np.exp(1j*th); w = np.ones_like(z)
for n in range(1, NMAX+1):
    if n % 1024 == 0: w = np.exp(1j*n*th)
    else: w = w*z
    C[n] = w.real.sum()
print('zero sums done %.1fs' % (time.time()-t0), flush=True)
mp.mp.dps = 30
A_T = mp.siegeltheta(T)/mp.pi + 1
S_T = float(N0 - A_T)
thT = 2*np.arctan(1/(2*T))
# tail int_T^inf 2(1-cos n theta(t)) A'(t) dt, t = T/w^6, Gauss-Legendre in w
xw, ww = np.polynomial.legendre.leggauss(400)
wq = 0.5*(xw+1); wwq = 0.5*ww
tq = T/wq**6; jac = 6*T/wq**7
Ap = (0.5*np.log(tq/(2*np.pi)) - 1/(48*tq**2))/np.pi
thq = 2*np.arctan(1/(2*tq))
base = Ap*jac*wwq
tail = np.zeros(NMAX+1)
for a in range(1, NMAX+1, 5000):
    ns = np.arange(a, min(a+5000, NMAX+1))
    tail[ns] = (4*np.sin(0.5*np.outer(ns, thq))**2*base).sum(axis=1)  # 2(1-cos x) = 4 sin^2(x/2), no cancellation
# mpmath check of the tail at a few n
for n in sorted(set(m for m in (10, 1000, 5083, NMAX) if m <= NMAX)):
    f = lambda t: 2*(1-mp.cos(n*2*mp.atan(1/(2*t))))*(0.5*mp.log(t/(2*mp.pi)) - 1/(48*t**2))/mp.pi
    q = mp.quad(f, [T, 2*T, 10*T, 100*T, mp.inf])
    print('tail check n=%d quad=%.12e GL=%.12e diff=%.2e' % (n, float(q), tail[n], float(q)-tail[n]))
ns = np.arange(NMAX+1)
gT = 4*np.sin(0.5*ns*thT)**2
lam = 2*N0 - 2*C + tail - gT*S_T
lam[0] = 0.0
S1max = 0.059*np.log(T*10) + 2.067
err = S1max*(4*ns/T**2 + 2*ns**2/(3*T**3))
np.save(out, np.vstack([lam, err]))
print('T=%.6f S(T)=%.6f; lam_1=%.15f lam_100=%.12f lam_NMAX=%.9f err_NMAX=%.2e; %.1fs' % (T, S_T, lam[1], lam[100], lam[NMAX], err[NMAX], time.time()-t0))
