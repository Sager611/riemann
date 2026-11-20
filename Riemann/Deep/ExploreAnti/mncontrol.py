# explore55 control: Milinovich-Ng J_{-1} >= |sum M(1-rho)/zeta'(rho)|^2 / sum|M(rho)|^2 (features M(rho)); does adding conj features help?
import numpy as np, sys
from moll import feats
from rs import zprime
g_all = np.loadtxt(sys.argv[1])
for (a, b) in [(10000, 20000), (50000, 100000)]:
    g = g_all[a:b]; zp = zprime(g); T = g[-1]
    for th in [0.5, 0.9]:
        D, Mp, Mm = feats(g, zp, T**th, 3); U = 1/D; J = np.sum(abs(U)**2)
        def proj(V):
            bb = np.conj(V).T @ U; G = np.conj(V).T @ V; return float(np.real(np.conj(bb) @ np.linalg.solve(G, bb)))/J
        print(f"zeros[{a}:{b}] theta={th}: J_-1/((3/pi^3)(T2-T1))={J/((3/np.pi**3)*(g[-1]-g[0])):.3f}  MN(Mp)={proj(Mp):.4f}  conj-only(Mm)={proj(Mm):.4f}  both={proj(np.concatenate([Mp, Mm], 1)):.4f}  MN-pred theta/(1+theta)={th/(1+th):.4f}", flush=True)
