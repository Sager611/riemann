# window_count.py: the band-limited model on a shifted band tau_c + [-tmax, tmax] (modulation by e^{i tau_c x} commutes with P):
# complex Hankel matrix c_{i+l+1}, c_m = (1/2pi) int phi~(tau_c + theta/dx) e^{-i m theta} dtheta, phi~ = exp(i psi~),
# psi~ = psi(tau_c + .) - psi(tau_c) on |tau'| <= tau1, closed by two lines of one negative slope s (V_R - V_L = 2 pi K) on the circle.
# Prints the singular values above 1e-3 and the one-sided semiclassical count (1/2pi) int_{|tau'|<=tau1} D_+ in the window.
# Usage: window_count.py primes tau_c tmax depth [b1]
import os
for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS'): os.environ[k] = '1'
import sys, json, math, numpy as np, scipy.linalg as sla
from scipy.special import loggamma, digamma
primes = [int(q) for q in sys.argv[1].split(',')]; tc = float(sys.argv[2]); tmax = float(sys.argv[3]); depth = float(sys.argv[4])
b1 = float(sys.argv[5]) if len(sys.argv) > 5 else 0.75
dx = math.pi/tmax; N = int(round(depth/dx)); tau1 = b1*tmax; tau2 = 0.5*(1 + b1)*tmax
def psi(tau):
    out = -2*(np.imag(loggamma(0.25 + 0.5j*tau)) - 0.5*tau*math.log(math.pi))
    for p in primes:
        a = p**-0.5; L = math.log(p); out = out + 2*np.arctan2(a*np.sin(tau*L), 1 - a*np.cos(tau*L))
    return out
def D(tau):
    out = -(np.real(digamma(0.25 + 0.5j*tau)) - math.log(math.pi))
    for p in primes:
        a = p**-0.5; L = math.log(p); c = np.cos(tau*L); out = out + 2*L*(a*c - a*a)/(1 - 2*a*c + a*a)
    return out
p0 = float(psi(np.array([tc]))[0])
pc = lambda tp: psi(tc + tp) - p0
pR = float(pc(np.array([tau1]))[0]); pL = float(pc(np.array([-tau1]))[0])
s0 = min(float(D(np.array([tc + tau1]))[0]), float(D(np.array([tc - tau1]))[0]), -0.5)
K = round((pR - pL + 2*s0*(tmax - tau1))/(2*math.pi))
s = (2*math.pi*K - (pR - pL))/(2*(tmax - tau1))
def step(x):
    x = np.clip(x, 0, 1); a = np.exp(-1/np.maximum(x, 1e-300)); b = np.exp(-1/np.maximum(1 - x, 1e-300))
    return np.where(x <= 0, 0.0, np.where(x >= 1, 1.0, a/(a + b)))
def psit(tp):
    a = np.abs(tp); bl = step((a - tau1)/(tau2 - tau1))
    lin = np.where(tp > 0, pR + s*(tp - tau1), pL + s*(tp + tau1))
    return (1 - bl)*pc(tp) + bl*lin
Nth = 1 << int(math.ceil(math.log2(4*N + 4096)))
k = np.arange(Nth); th = 2*math.pi*k/Nth; th = np.where(th >= math.pi, th - 2*math.pi, th)
c = np.fft.fft(np.exp(1j*psit(th/dx)))/Nth
H = np.empty((N, N), dtype=np.complex64); ar = np.arange(N)
for i0 in range(0, N, 128):
    i1 = min(N, i0 + 128); H[i0:i1] = c[1 + np.arange(i0, i1)[:, None] + ar[None, :]]
sv = sla.svdvals(H, overwrite_a=True, check_finite=False); del H
tg = np.arange(-tau1, tau1, 0.005) + 0.0025
nsc = float(np.clip(D(tc + tg), 0, None).sum()*0.005/(2*math.pi))
print(json.dumps({'primes': primes, 'tau_c': tc, 'tmax': tmax, 'tau1': tau1, 'N': N, 'slope': s, 'N_sc_window': round(nsc, 3),
                  'counts': {str(x): int((sv > x).sum()) for x in (0.9, 0.5, 0.3, 0.1, 0.05, 0.01)}, 'sum_k2': round(float((sv.astype(float)**2).sum()), 4), 'top': [round(float(x), 4) for x in sv[:12]]}))
