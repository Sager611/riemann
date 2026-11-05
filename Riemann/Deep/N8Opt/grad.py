# grad.py -- first-order sensitivity of the 8-point LP value s(f) to the window coefficients a_j
# (envelope theorem on the converged cutting-plane LP of lp8.py; floating point; proposes only)
import json, sys
import numpy as np
from scipy.optimize import linprog
cuts = np.load(sys.argv[1]); G = cuts['G']; mu = cuts['duals']; z = cuts['z']; s = float(cuts['s'])
d = json.load(open('src/anc/certificate_line/inputs.json'))
a = np.array(d['window_coefficients'], float) / 1e6
J = 14; aa = np.zeros(J); aa[:len(a)] = a
NP = 8; NG = 7; pairs = [(i, j) for i in range(NP) for j in range(i + 1, NP)]
PI = np.array([p[0] for p in pairs]); PJ = np.array([p[1] for p in pairs])
iA, iK, iW, iB = 0, 1, 1 + NP, 1 + NP + len(pairs)
x = z[iK:iK + NP].sum(); W = z[iW:iW + len(pairs)]
def K(t, coef):
    s_ = np.sinc(t)
    for j, aj in enumerate(coef, 1): s_ = s_ + aj / 2 * (np.sinc(t - j) + np.sinc(t + j))
    return s_
y = np.concatenate([np.zeros((G.shape[0], 1)), np.cumsum(G, axis=1)], axis=1)
D = y[:, PJ] - y[:, PI]
Kd = K(D, aa)
act = mu > 1e-12
print('cuts %d, active duals %d, sum of duals %.4f, s_LP %.7f, x %.6f' % (G.shape[0], act.sum(), mu.sum(), s, x))
grad = []
for j in range(1, J + 1):
    dK = (np.sinc(D - j) + np.sinc(D + j)) / 2
    dPhi = float(mu @ ((2 * Kd * dK) @ W))
    dC = aa[j - 1] + ((-1) ** j / j**2 - aa[j - 1] / (2 * j**2)) / np.pi**2
    g = (dPhi - dC) / (1 - x)
    grad.append(g)
    print('a_%-2d = %+.6f  dPhi/da %+.5f  dC/da %+.5f  ds/da %+.5f' % (j, aa[j - 1], dPhi, dC, g))
grad = np.array(grad)
print('|grad| = %.5f ; a step of length 0.002 along it gains ~%.2e in s (first order)' % (np.linalg.norm(grad), 0.002 * np.linalg.norm(grad)))

# admissible separation range for the line certificate of this window
t = np.linspace(0, 1.2, 120001); K2 = K(t, aa) ** 2
bad = np.nonzero(K2 < x)[0]
th_max = t[bad[0]] if len(bad) else 1.2
fhalf = 1 + sum(aa[j - 1] * (-1) ** j for j in range(1, J + 1))
th = np.linspace(0.3, 1.0, 70001); lam_need = 2 * fhalf / (th * (1 + np.sinc(th)))
th_min = th[np.nonzero(lam_need < 2)[0][0]]
print('theta: K^2 >= x on [0,theta] needs theta < %.4f; majorant with Lambda < 2 at u=1/2 needs theta > %.4f (f(1/2) = %.6f)' % (th_max, th_min, fhalf))
