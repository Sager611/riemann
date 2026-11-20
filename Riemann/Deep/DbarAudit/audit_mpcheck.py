# audit_mpcheck.py (dbaraudit): mpmath (dps 30) vs double for the cost ratio U_all / sum_g U_g on a reduced chain configuration
# (first 60 zeros, every 8th zero n >= 20 -> m = 6 zeros at spacing 0.99 L^{-3/4}), L = 40, 160: is a ratio near 1 (or < 1) real?
import numpy as np, mpmath as mp, importlib.util, time
spec = importlib.util.spec_from_file_location('dU', '/Users/adriansager/playground/riemann/Riemann/Deep/DbarAudit/rerun/dbar_U.py')
dU = importlib.util.module_from_spec(spec); spec.loader.exec_module(dU)
mp.mp.dps = 30; t0 = time.time()
def Kmp(z, w, L):  # explore59's kernel K(z,w), a = 1/2, real nodes
    a = mp.mpf(1) / 2; d = z - w
    I = mp.mpf(L) if d == 0 else (1 - mp.expj(-L * d)) / (1j * d)
    return ((a * a + w * w) * I + (a - 1j * w) + (a + 1j * w) * mp.expj(-L * d)) / (2 * a)
def cost_mp(Z, L):
    n = len(Z); G = mp.matrix(n, n)
    for i in range(n):
        for j in range(n): G[i, j] = Kmp(Z[i], Z[j], L)
    one = mp.matrix([1] * n); c = mp.lu_solve(G, one)
    return mp.re(sum(c[i] for i in range(n)))
def cost_d(Z, L):
    G = dU.Kblock(np.array([float(z) for z in Z]), np.array([float(z) for z in Z]), L); G = 0.5 * (G + G.conj().T)
    one = np.ones(len(Z), complex); return float(np.real(np.vdot(one, np.linalg.solve(G, one))))
for L in (40, 160):
    step = mp.mpf(0.99) * mp.mpf(L) ** mp.mpf(-0.75)
    groups = []
    for n, g in enumerate(dU.zz[:60]):
        g = mp.mpf(g)
        for s in (1, -1):
            if n >= 20 and n % 8 == 0: groups.append([s * (g + k * step) for k in range(6)])
            else: groups.append([s * g])
    Z = [z for gr in groups for z in gr]
    ua, us = cost_mp(Z, L), sum(cost_mp(gr, L) for gr in groups)
    da, ds = cost_d(Z, L), sum(cost_d(gr, L) for gr in groups)
    print(f'L={L} n={len(Z)}  mp: L*U_all={mp.nstr(L*ua,12)} ratio={mp.nstr(ua/us,12)}   double: L*U_all={L*da:.12f} ratio={da/ds:.12f}', flush=True)
print('time', round(time.time() - t0, 1), 's')
