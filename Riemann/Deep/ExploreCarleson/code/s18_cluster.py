# explore42 s18: is the zero form controlled by separation (Carleson/Shapiro-Shields lower Riesz bound) or by cluster multiplicity (box constant)?
# Main-term Burnol Gram (explore38), c = (rho-1)/rho^2. B/|c~|^2 vs lmin for (a) two coalescing zeros, (b) first 400 zero pairs with inserted near-duplicates.
import math
import numpy as np
zc = np.load('/Users/adriansager/playground/riemann/.scratch/explore42/zeros_cache.npz'); g = zc['g']
rng = np.random.default_rng(3)
def theta(t): return t / 2 * np.log(t / (2 * np.pi)) - t / 2 - np.pi / 8 + 1 / (48 * t) + 7 / (5760 * t ** 3)
def dtheta(t): return 0.5 * np.log(t / (2 * np.pi)) - 1 / (48 * t ** 2) - 7 / (1920 * t ** 4)
def Phi(t, L): return t * L + 2 * theta(np.abs(t)) * np.sign(t) + 6 * np.arctan(2 * t)
def dPhi(t, L): return L + 2 * dtheta(np.abs(t)) + 3 / (0.25 + t ** 2)
def form(xpos, L):
    x = np.concatenate([xpos, -xpos])
    P = Phi(x, L); dP = dPhi(x, L)
    d = x[None, :] - x[:, None]
    with np.errstate(divide='ignore', invalid='ignore'):
        G = (1 - np.exp(1j * (P[:, None] - P[None, :]))) / (1j * d)
    G[np.arange(len(x)), np.arange(len(x))] = dP
    H = G.T
    rho = 0.5 + 1j * x; c = (rho - 1) / rho ** 2
    B = float(np.real(np.vdot(c, np.linalg.solve(H, c))))
    ct = float(np.sum(np.abs(c) ** 2 / dP))
    s = 1 / np.sqrt(dP); ev = np.linalg.eigvalsh(s[:, None] * H * s[None, :])
    return B / ct, ev[0], ev[-1]
for N in [10, 1280, 10 ** 10]:
    L = math.log(N)
    g1 = g[0]; dp = dPhi(np.array([g1]), L)[0]
    out = []
    for sig in [4, 2, 1, 0.5, 0.25, 0.1, 0.03, 0.01, 0.001]:
        r, lmin, lmax = form(np.array([g1, g1 + sig * 2 * np.pi / dp]), L)
        out.append(f'{sig}:{r:.3f}/{lmin:.1e}')
    print(f'(a) N={N:.0e} two zeros at normalized distance sigma -> B/|c~|^2 / lmin: ' + '  '.join(out), flush=True)
for N in [10, 1280, 10 ** 10]:
    L = math.log(N)
    base = g[:200]
    r0 = form(base, L)
    rows = [f'zeta: {r0[0]:.4f} (lmin {r0[1]:.3f}, lmax {r0[2]:.2f})']
    for m, frac, sig in [(2, 0.3, 0.05), (2, 1.0, 0.05), (3, 0.3, 0.05), (2, 1.0, 0.005)]:
        idx = rng.choice(200, int(frac * 200), replace=False)
        extra = []
        for j in idx:
            dp = dPhi(np.array([base[j]]), L)[0]
            for k in range(1, m):
                extra.append(base[j] + k * sig * 2 * np.pi / dp)
        x = np.sort(np.concatenate([base, np.array(extra)]))
        r = form(x, L)
        # weighted cluster count: sum m^2|c~|^2 / sum m|c~|^2 for the clusters created
        rows.append(f'+{m-1} dup at {int(frac*100)}% (sigma {sig}): {r[0]:.4f} (lmin {r[1]:.1e}, lmax {r[2]:.2f})')
    print(f'(b) N={N:.0e} ' + ' | '.join(rows), flush=True)
