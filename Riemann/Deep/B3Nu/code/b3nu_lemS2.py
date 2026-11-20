# b3nu Lemma S falsification WITH complex shifts: exact T_+-(n) from B3Cert's code/b3c_fals.py (imported read-only: complex beta in the
# profile, Hankel a_+- for complex d by Laguerre quadrature) vs the leading term
#   lead_s = s*(i/2)*H/sqrt((1+v)(2+v)) * e^{-i s pi/4} * e^{-y*^2+(1-2beta)y*/H} * xi*^{-beta-z/2} * X*^{-1/2} * a_s(X*, d) * e^{i Phi_s(x)}.
import sys, math, numpy as np
sys.path.insert(0, 'Riemann/Deep/B3Cert/code'); sys.path.insert(0, 'Riemann/Deep/B3Nu/code')
from b3c_fals import setup, Tpm, apm
from b3nu_T import Phis
def lead(T, H, x, s, al, be):
    z, d = al - be, al + be
    b = 2 * np.pi * x / T; v = (b + s * np.sqrt(b * b + 4 * b)) / 2; ys = 0.5 * H * np.log1p(v)
    xis = T / (2 * np.pi * (1 + v)); Xs = 4 * np.pi * np.sqrt(x * xis)
    ap, am = apm(np.atleast_1d(Xs), d); a = ap if s == 1 else am
    return s * 0.5j * H / np.sqrt((1 + v) * (2 + v)) * np.exp(-1j * s * np.pi / 4) * np.exp(-ys ** 2 + (1 - 2 * be) * ys / H) * xis ** (-be - z / 2) * Xs ** -0.5 * a * np.exp(1j * Phis(T, x, s))
for T, th in [(1e6, 0.45), (1e7, 0.45), (1e6, 0.49)]:
    H = T ** th; tau = T / H ** 2; K = T / H
    for (al, be) in [(0.05 + 0.1j, -0.03 + 0.07j), (0.2 - 0.1j, 0.03 + 0.2j)]:
        w, Phi, N = setup(T, H, be)
        worst = 0.0; mx = 0.0
        for r in np.exp(np.linspace(math.log(1e-3), math.log(3), 24)):
            x = r * tau; tp, tm = Tpm(T, H, x, al, be, w, Phi); sc = H * (8 * np.pi * x * T) ** -0.25
            lp = lead(T, H, x, 1, al, be)[0]; lm = lead(T, H, x, -1, al, be)[0]
            worst = max(worst, abs(tp - lp) / sc, abs(tm - lm) / sc); mx = max(mx, abs(tp) / sc, abs(tm) / sc)
        print(f"T={T:.0e} th={th} K={K:.0f} shifts=({al},{be}): max|T_s - lead_s|/(H X0^-1/2) = {worst:.2e}; K*err = {K*worst:.3f}; max|T_s|/(H X0^-1/2) = {mx:.4f}", flush=True)
