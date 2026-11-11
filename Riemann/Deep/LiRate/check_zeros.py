# Independent check of the FFT lambda_n via the zero sum (Odlyzko's first 10^5 zeros, |err| <= 3e-9):
# lambda_n = sum_{gamma<=H} g_n(gamma) + int_H^inf g_n A' dt - g_n(H) S(H) - int_H^inf S g_n'  (last two O(n^2/H^2))
import math, numpy as np, mpmath as mp
z = np.loadtxt('../../../.scratch/lirate/zeros1')
H = (z[-1] + z[-2]) / 2
def theta(t): return 2 * np.arctan(1 / (2 * t))
lam = {}
for l in open('logs/lam_r999.txt'):
    if l.startswith('#'): continue
    a, b = l.split()[:2]; lam[int(a)] = float(b)
mp.mp.dps = 20
def Ap(t): return (mp.re(mp.digamma(mp.mpf(1)/4 + 1j*t/2))/2 - mp.log(mp.pi)/2)/mp.pi
thRS = lambda t: mp.im(mp.loggamma(mp.mpf(1)/4 + 1j*t/2)) - t/2*mp.log(mp.pi)
SH = len(z) - (float(thRS(H))/math.pi + 1)
print(f'H = {H:.3f}, S(H) = {SH:.4f}')
for n in (100, 500, 1000, 2000):
    s = np.sum(2 * (1 - np.cos(n * theta(z))))
    tail = float(mp.quad(lambda t: 2*(1-mp.cos(n*2*mp.atan(1/(2*t))))*Ap(t), [H, 10*H, 100*H, mp.inf]))
    est = s + tail
    print(f'n={n}: zero-sum {est:.6f}  FFT {lam[n]:.6f}  diff {est-lam[n]:+.2e}  (boundary term ~ {2*(1-math.cos(n*theta(H)))*abs(SH):.1e})')
