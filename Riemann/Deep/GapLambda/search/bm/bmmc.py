"""Monte Carlo of the raw 9-dim integrand of BM Thm smoothed (l.307-311), t3,t4 sampled (no analytic step),
versus the quadrature of bmfloat (t3,t4 integrated analytically). P = 1 - 5.8x + 6.4x^2, theta=1/4, kappa=3.18, v=1.26."""
import sys, time
sys.path.insert(0, '.')
import numpy as np
from numpy.polynomial import Polynomial, Legendre
import bmfloat as bf
th, k, v, r = 0.25, 3.18, 1.26, 1
P = Polynomial([1, -5.8, 6.4])
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 5)
N, reps = 400000, 40
acc0, acc1 = [], []
t0 = time.time()
for _ in range(reps):
    x = rng.random(N)
    R = 1 - x
    def tri():
        p, q = rng.random(N), rng.random(N); f = p + q > 1; p[f], q[f] = 1 - p[f], 1 - q[f]
        return R * p, R * q
    x1, x2 = tri(); x3, x4 = tri()
    vol = R ** 2 / 2 * R ** 2 / 2          # density correction: uniform x, uniform triangles of area R^2/2
    t1, t2, t3, t4 = rng.random((4, N))
    a = 1 - th * (x1 + x3); b = 1 - th * (x2 + x4)
    Phi = th * k * np.pi * (x2 - x4 - (x3 - x4) * t3 + (x1 - x2) * t4) - k * np.pi * (a * t1 - b * t2) * (t3 - t4)
    g = np.exp(1j * Phi) * a * b * (th * (x1 - x2) + a * t1 - b * t2) * (th * (x3 - x4) + a * t1 - b * t2) \
        * P(1 - x - x1 - x2) * P(1 - x - x3 - x4) * vol
    Q = v - th * (x + x1 + x2 + x3 + x4) - a * t1 - b * t2
    acc0.append(np.mean(g).real); acc1.append(np.mean(g * Q * Q).real)
m0, s0 = np.mean(acc0), np.std(acc0) / np.sqrt(reps); m1, s1 = np.mean(acc1), np.std(acc1) / np.sqrt(reps)
M, _ = bf.mats(th, r, k, 3, 10)
Pl = P.convert(kind=Legendre, domain=[0, 1]).coef
q0 = Pl @ M[0] @ Pl; q1 = Pl @ (v * v * M[0] - 2 * v * M[1] + M[2]) @ Pl
print('MC c0 %.6e +- %.1e  quad %.6e  z=%.2f' % (m0, s0, q0, (m0 - q0) / s0))
print('MC c1 %.6e +- %.1e  quad %.6e  z=%.2f' % (m1, s1, q1, (m1 - q1) / s1))
print('MC h1 = %.4f, quad h1 = %.6f  (%.0fs, %d samples)' % (m0 / (k * k * m1), q0 / (k * k * q1), time.time() - t0, N * reps))
