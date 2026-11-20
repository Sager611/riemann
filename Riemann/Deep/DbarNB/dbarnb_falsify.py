# dbarnb_falsify.py : falsification tests for Lemma R'' / R' (dbarnb, 20 Nov 2026)
# Uses explore59's kernel (dbar_U.Kmat). Modes: dual | clusters | orth
import sys, numpy as np, time
sys.argv = sys.argv[:1] + sys.argv[1:]
mode = sys.argv[1] if len(sys.argv) > 1 else 'dual'
import importlib.util
spec = importlib.util.spec_from_file_location('dU', '/Users/adriansager/playground/riemann/Riemann/Deep/DbarNB/dbar_U.py')
dU = importlib.util.module_from_spec(spec); spec.loader.exec_module(dU)
Kblock, zz = dU.Kblock, dU.zz
t0 = time.time()
def NGram(Zs, ell):  # <e_k, e_j>_N, e_j = e^{i gamma_j u}: int_0^ell (g_j g_k + 1/4) e^{i(g_k-g_j)u} du + 1/2 + 1/2 e^{i(g_k-g_j) ell}
    gj = Zs[:, None]; gk = Zs[None, :]; D = gk - gj
    small = np.abs(D) < 1e-12; Ds = np.where(small, 1.0, D)
    I = np.where(small, ell, (np.exp(1j * D * ell) - 1) / (1j * Ds))
    return (gj * gk + 0.25) * I + 0.5 + 0.5 * np.exp(1j * D * ell)
def cost(G):
    one = np.ones(G.shape[0], complex); return float(np.real(np.vdot(one, np.linalg.solve(G, one))))
if mode == 'dual':
    Zs = np.concatenate([zz[:40], -zz[:40], [100.0, 100.0 + 0.3 / 20, 100.0 + 0.7 / 20], [-100.0, -100.0 - 0.3 / 20, -100.0 - 0.7 / 20]])
    for L in (5.0, 20.0, 80.0):
        GK = Kblock(Zs, Zs, L); GN = NGram(Zs, L)
        print(f'DUAL L={L:5.1f} max|K-N|/max|K| = {np.max(np.abs(GK - GN)) / np.max(np.abs(GK)):.2e}  L*U_K={L*cost(GK):.10f}  L*U_N={L*cost(GN):.10f}')
if mode in ('clusters', 'orth'):
    m = int(sys.argv[2]) if len(sys.argv) > 2 else 4
    hstep = float(sys.argv[3]) if len(sys.argv) > 3 else 0.3   # cluster spacing h/L (unresolved at every L)
    Kz = 300; every = 8
    base = zz[:Kz]
    for L in (10.0, 20.0, 40.0, 80.0, 160.0):
        pts, grp = [], []
        gid = 0
        for n, g in enumerate(base):
            if n >= 20 and n % every == 0:
                for s in (1, -1):
                    for k in range(m): pts.append(s * (g + k * hstep / L)); grp.append(gid)
                    gid += 1
            else:
                for s in (1, -1): pts.append(s * g); grp.append(gid); gid += 1
        Zs = np.array(pts); grp = np.array(grp)
        G = Kblock(Zs, Zs, L)
        G = 0.5 * (G + G.conj().T)
        uall = cost(G)
        usum = 0.0
        B = np.zeros_like(G)
        for q in range(gid):
            idx = np.where(grp == q)[0]
            sub = G[np.ix_(idx, idx)]; usum += cost(sub); B[np.ix_(idx, idx)] = sub
        line = f'm={m} h={hstep} L={L:6.1f} n={len(Zs)}  L*U_all={L*uall:.6f}  L*sum_g U_g={L*usum:.6f}  ratio={uall/usum:.5f}'
        if mode == 'orth':
            import scipy.linalg as sl
            ev = sl.eigh(B, G, eigvals_only=True)
            line += f'  sup_p sum_g N(p_g)/N(p)={ev.max():.5f}  inf={ev.min():.5f}'
        print(line, flush=True)
print('time', round(time.time() - t0, 1), 's')
