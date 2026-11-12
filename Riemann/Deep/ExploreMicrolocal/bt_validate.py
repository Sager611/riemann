# explore47 validation (double precision, scipy): explicit-formula identity for the Beurling–Selberg majorant/minorant
# R^± of [-t,t] (type Δ, primes n < e^Δ) on Odlyzko zeros1 (first 10^5 zeros).
import sys, os, math, time
import numpy as np
from scipy.special import polygamma, psi
import mpmath as mp
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from btlib import prime_powers, load, Rpm
t0 = time.time()
def H(x, sg):
    ax = np.abs(x); s2 = (np.sin(np.pi * x) / np.pi) ** 2
    core = np.sign(x) * (1 - s2 * (polygamma(1, ax) + polygamma(1, 1 + ax) - 2 / ax))
    return core + sg * s2 / np.where(ax == 0, 1, x * x)
def R(u, t, Dt, sg):
    d = Dt / (2 * np.pi)
    return 0.5 * (H(d * (t + u), sg) + H(d * (t - u), sg))
xs = np.array([-3.7, -1.2, -0.3, 0.2, 0.9, 1.5, 4.25, 10.1, 1e3 + 0.37])
print('H sanity', bool(np.all(H(xs, 1) >= np.sign(xs) - 1e-12) and np.all(H(xs, -1) <= np.sign(xs) + 1e-12)), flush=True)
_, _, vals = load('zeros1'); g = np.array([float(v) for v in vals])
n, lam, base = prime_powers(3000); logn = np.log(n.astype(float))
w = lambda u: 0.5 * (np.real(psi(0.25 + 0.5j * u)) - np.log(np.pi))
for (t, Dt) in [(100.0, 3.0), (100.0, 5.0), (1000.0, 4.0), (1000.0, 7.0)]:
    d = Dt / (2 * np.pi)
    for sg in (1, -1):
        lhs = 2 * np.sum(R(g, t, Dt, sg))
        # ∫(R−χ)w over u>0 (×2), split at t; dense grid near t, geometric far grid
        u1 = np.linspace(0, t, 400001)[:-1]; u2 = t + np.concatenate([np.linspace(0, 400 / d, 400001)[1:], 400 / d * np.geomspace(1, 1e6, 200001)[1:]])
        f1 = (R(u1, t, Dt, sg) - 1) * w(u1); f2 = R(u2, t, Dt, sg) * w(u2)
        tail = 0.0
        integ_E = 2 * (np.trapezoid(f1, u1) + np.trapezoid(f2, u2)) / np.pi
        th = float(mp.siegeltheta(t)); smooth = 2 * th / np.pi + integ_E
        pole = float(mp.re(Rpm(0.5j, t, Dt, sg) + Rpm(-0.5j, t, Dt, sg)))
        m = logn < Dt; e = logn[m] / Dt
        Rh = (2 * np.pi / Dt) * ((1 / np.pi + (1 - e) / np.tan(np.pi * e)) * np.sin(t * logn[m]) + sg * (1 - e) * np.cos(t * logn[m]))
        prime = float(np.sum(lam[m] / np.sqrt(n[m]) * Rh)) / np.pi
        rhs = pole + smooth - prime
        S = np.sum(g <= t) - th / np.pi - 1
        tau = integ_E / 2
        print(f't={t:.0f} Δ={Dt} {"+" if sg>0 else "-"}: LHS={lhs:.6f} RHS={rhs:.6f} diff={lhs-rhs:.1e} pole-2={pole-2:.1e} '
              f'tau={tau:.5f} vs ±log(t/2π)/(2Δ)={sg*math.log(t/2/math.pi)/(2*Dt):.5f} O={-prime/2:.4f} S={S:.4f} '
              f'bound={tau - prime/2 + (pole-2)/2:.4f}', flush=True)
print('time', round(time.time() - t0, 1), flush=True)
