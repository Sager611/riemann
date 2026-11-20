# ceilaudit B4: bandwidth freedom lambda < 1 in zeta23/TEN (X = e^L, L = lambda log(T/2pi), lambda <= 1 admitted by the lemmas).
# R_lambda(psi) = A/lambda + lambda*B, A = int psi^2/(int psi)^2, B = iint |u-v| psi psi/(int psi)^2 (psi_0 = 1 gives 1/lambda + lambda/3).
import numpy as np
from numpy.polynomial.legendre import leggauss
x, w = leggauss(400); s = x / 2; ws = w / 2
for name, psi in [("psi_0 = 1", np.ones_like(s)), ("MT cos(sqrt2 s)", np.cos(np.sqrt(2) * s))]:
    I = ws @ psi; A = (ws @ psi ** 2) / I ** 2
    B = (ws * psi) @ np.abs(s[:, None] - s[None, :]) @ (ws * psi) / I ** 2
    print(f"{name}: A = {A:.7f}, B = {B:.7f}, R_1 = {A + B:.7f} (H = 2 - R = {2 - A - B:.7f}), lambda* = sqrt(A/B) = {np.sqrt(A / B):.4f} (> 1: lambda = 1 optimal)")
    for lam in [0.999, 0.99, 0.95, 0.9]:
        print(f"   lambda = {lam}: H loss {A / lam + lam * B - A - B:.2e} (first order (A - B)(1 - lambda) = {(A - B) * (1 - lam):.2e})")
print("TEN's whole defect gain over H on any window <= 0.6740088 - 0.6725007 = 1.5e-3, so lambda < 1 can help only for 1 - lambda < 2.2e-3 and")
print("only if the defect gain grows faster than (A - B) = 0.685 per unit decrease of lambda; the rank cap d = lambda N (AF 7(e)) also caps kappa <= lambda.")
