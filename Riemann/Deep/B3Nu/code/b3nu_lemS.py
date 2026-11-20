# b3nu Lemma S falsification: exact T_sigma(x) (x = n/pq, zero shifts) vs the 2D stationary-phase leading term
#   L_sigma(x) = sigma*(i/2) * H / sqrt((1+v)(2+v)) * e^{-i sigma pi/4} e^{-y*^2 + y*/H} X*^{-1/2} a_sigma(X*) e^{i Phi_sigma(x)},
#   v^2 = (2 pi x/T)(1+v), sign v = sigma, y* = (H/2)log(1+v), X* = 4 pi sqrt(x T/(2 pi (1+v)))   (chi(xi*/T) = 1 here).
# Reports max |T - L| / (H X0^{-1/2}) over r = x/tau in [1e-3, 3]; Lemma S predicts O(L^C/K).
import math, numpy as np
from b3nu_T import setup, Tbr, Phis
def lead(T, H, x, s):
    b = 2 * np.pi * x / T; v = (b + s * np.sqrt(b * b + 4 * b)) / 2; ys = 0.5 * H * np.log1p(v)
    Xs = 4 * np.pi * np.sqrt(x * T / (2 * np.pi * (1 + v)))
    a = 1 - 1j * s / (8 * Xs) - 9 / (128 * Xs ** 2)
    return s * 0.5j * H / np.sqrt((1 + v) * (2 + v)) * np.exp(-1j * s * np.pi / 4) * np.exp(-ys ** 2 + ys / H) * Xs ** -0.5 * a * np.exp(1j * Phis(T, x, s))
for T, th in [(1e5, 0.45), (1e6, 0.45), (1e7, 0.45), (1e6, 0.40), (1e6, 0.49), (1e7, 0.49)]:
    H = T ** th; tau = T / H ** 2; K = T / H
    w, Phi, ws, Ps = setup(T, H)
    r = np.exp(np.linspace(math.log(1e-3), math.log(3), 60)); x = r * tau
    sc = H * (8 * np.pi * x * T) ** -0.25   # H X0^{-1/2}, X0 = 2 sqrt(2 pi) sqrt(x T)
    errs = []; mx = []
    for s in (1, -1):
        ex = Tbr(T, H, x, ws, Ps, s); le = lead(T, H, x, s)
        errs.append(np.max(np.abs(ex - le) / sc)); mx.append(np.max(np.abs(ex) / sc))
    print(f"T={T:.0e} th={th} K={K:.0f} tau={tau:.2f}: max|T_s - lead_s|/(H X0^-1/2) = {errs[0]:.2e} (s=+1), {errs[1]:.2e} (s=-1); K*err = {K*max(errs):.2f}; max|T_s|/(H X0^-1/2) = {max(mx):.4f}", flush=True)
