# nb2fals f6: data-free projection of the GP posterior sd of the constant for windows u in [log 100, u_hi],
# trend = 1 (C3/L) or 2 (C3/L + C4/L^2); hyperparameters from the N <= 4000 REML fits (f3 (a)); K = 50 lines; grid step 0.02.
import math, numpy as np
from scipy.linalg import cho_factor, cho_solve, toeplitz
Z = np.load('/Users/adriansager/playground/riemann/.scratch/nb2fals/zeros100.npz'); gam = Z['gam']; zp = np.abs(Z['zp']); K = 50; d = 0.02
i, j = np.triu_indices(K, 1); om = gam[j] - gam[i]; la = np.log(gam[0] * gam[1] / (gam[i] * gam[j])); lb = np.log((gam[1] - gam[0]) / om)
HP = {'B (alpha=1, beta=1)': (0.706, 7.84, 1.62e-4, 1.0, 1.0), 'fitted (alpha=1.05, beta=0.80)': (0.518, 7.75, 1.62e-4, 1.05, 0.80)}
def sdc(hp, uh, trend):
    tp, ts, s0, al, be = hp; ug = np.arange(math.log(100), uh + 1e-12, d); n = len(ug); lag = np.arange(n) * d
    s2 = (tp ** 2) * np.exp(2 * (al * la + be * lb)) / (gam[0] * gam[1]) ** 2
    Km = toeplitz(s2 @ np.cos(np.outer(om, lag)))
    e = np.exp(-ug / 2); Km += toeplitz(((ts / (gam[:K] ** 2 * zp[:K])) ** 2) @ np.cos(np.outer(gam[:K], lag))) * np.outer(e, e)
    Km[np.diag_indices(n)] += s0 ** 2
    X = np.column_stack([np.ones(n)] + [ug ** (-p) for p in range(1, trend + 1)]); cf = cho_factor(Km, lower=True)
    return math.sqrt(np.linalg.inv(X.T @ cho_solve(cf, X))[0, 0])
for name, hp in HP.items():
    for trend in (1, 2):
        row = [(uh, sdc(hp, uh, trend)) for uh in (math.log(4000), math.log(5500), 10.0, 12.0, 15.0, 20.0, 30.0, 45.0)]
        print(f'shape {name}, trend {"C3/L" if trend == 1 else "C3/L+C4/L^2"}: ' + ', '.join(f'N=e^{uh:.2f}: {s:.4f}' for uh, s in row), flush=True)
