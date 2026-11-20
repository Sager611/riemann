# audit_chain.py (dbaraudit): full OU-norm cost L*U for RH-allowed worst cases not in dbarnb_falsify.py:
# (a) chains: every 8th zero (n >= 20) of the first 300 replaced by m zeros at spacing 0.99*s0, s0 = L^{-3/4} (one proof-group of
#     width (m-1)*0.99*s0, resolved by the band since s0*L = L^{1/4} > 1);  (b) dense clusters: every 2nd zero of the first 150
#     replaced by an m-cluster at spacing 0.3/L (unresolved).  Prints L*U_all, L*sum_g U_g and the ratio.
import sys, numpy as np, time, importlib.util
spec = importlib.util.spec_from_file_location('dU', '/Users/adriansager/playground/riemann/Riemann/Deep/DbarAudit/rerun/dbar_U.py')
dU = importlib.util.module_from_spec(spec); spec.loader.exec_module(dU)
Kblock, zz = dU.Kblock, dU.zz
def cost(G):
    one = np.ones(G.shape[0], complex); return float(np.real(np.vdot(one, np.linalg.solve(G, one))))
t0 = time.time()
for case, m, Kz, every in (('chain', 6, 300, 8), ('chain', 8, 300, 8), ('dense', 4, 150, 2)):
    for L in (10.0, 20.0, 40.0, 80.0, 160.0):
        step = 0.99 * L ** -0.75 if case == 'chain' else 0.3 / L
        pts, grp, gid = [], [], 0
        for n, g in enumerate(zz[:Kz]):
            if n >= 20 and n % every == 0:
                for s in (1, -1):
                    for k in range(m): pts.append(s * (g + k * step)); grp.append(gid)
                    gid += 1
            else:
                for s in (1, -1): pts.append(s * g); grp.append(gid); gid += 1
        Zs = np.array(pts); grp = np.array(grp)
        G = Kblock(Zs, Zs, L); G = 0.5 * (G + G.conj().T)
        uall = cost(G); usum = sum(cost(G[np.ix_(np.where(grp == q)[0], np.where(grp == q)[0])]) for q in range(gid))
        print(f'{case} m={m} L={L:6.1f} n={len(Zs)} step*L={step*L:.3f}  L*U_all={L*uall:.6f}  L*sum_g U_g={L*usum:.6f}  ratio={uall/usum:.5f}  cond={np.linalg.cond(G):.1e}', flush=True)
print('time', round(time.time() - t0, 1), 's')
