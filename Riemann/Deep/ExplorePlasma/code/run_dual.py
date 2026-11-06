import sys, time; sys.path.insert(0, __file__.rsplit('/', 1)[0])
# T3: dual plasma f0' = Xi: C(x) = 2 pi i Xi+(x) = -2 pi S(x) + i pi Xi(x); neutral modes at -gamma_n, k_n^2 = 2 pi S(gamma_n);
# Penrose winding (number of unstable modes) of C around k^2 for Xi and Xi_off (delta = 0.25) on x in [-60, 60]
import numpy as np, mpmath as mp
from onesided import make_onesided, Phi_np
from offline_np import make as make_off
t0 = time.time()
F = make_onesided(Phi_np)
gam = [float(mp.zetazero(n).imag) for n in range(1, 11)]
print("k_n^2 = 2 pi S(gamma_n), n = 1..10:", ' '.join('%.5f' % (2*np.pi*s) for s in F(np.array(gam)).imag))
u, F0, F1, F2, info = make_off(0.25)
h = u[1] - u[0]; wtr = np.full(u.size, h); wtr[0] = wtr[-1] = h/2
def curve(K, xs, chunk=200):
    out = np.empty(xs.size, dtype=complex)
    for i in range(0, xs.size, chunk):
        E = np.exp(1j*np.multiply.outer(xs[i:i+chunk], u))
        out[i:i+chunk] = 2j*np.pi*2*(E*(K*wtr)).sum(-1)
    return out
xs = np.linspace(-60, 60, 120001)
Cz = curve(Phi_np(u), xs); Co = curve(F0, xs)
def winding(C, w):
    v = np.concatenate([C, C[:1]]) - w
    return int(np.round(np.sum(np.angle(v[1:]/v[:-1]))/(2*np.pi)))
ks = [0.30, 0.25, 0.20, 0.17, 0.15, 0.13, 0.12, 0.11]
print('k^2          :', ' '.join('%5.2f' % k for k in ks))
print('winding Xi   :', ' '.join('%5d' % winding(Cz, k) for k in ks))
print('winding Xioff:', ' '.join('%5d' % winding(Co, k) for k in ks))
def crossings(C):
    s = np.sign(C.imag); idx = np.where((s[1:]*s[:-1] < 0) & (xs[1:] < 0))[0]
    return [(round(float(xs[i]), 3), round(float(C[i].real), 5)) for i in idx]
print('positive-axis crossings (x, Re C), Xi   :', crossings(Cz)[-7:])
print('positive-axis crossings (x, Re C), Xioff:', crossings(Co)[-7:])
print('time %.1f' % (time.time() - t0))
