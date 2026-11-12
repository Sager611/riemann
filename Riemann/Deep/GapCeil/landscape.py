"""sup_v lambda_max(kappa)/kappa^2 on L2(phi) for a range of kappa (r=1, th=1/4)."""
import sys, time, resource
import numpy as np
import gcore as G
th = 0.25
t0 = time.time()
ks = [float(s) for s in sys.argv[1].split(',')]
N = int(sys.argv[2]) if len(sys.argv) > 2 else 16
for kap in ks:
    c = kap * np.pi
    D = int(min(90, max(36, 2.2 * c * th + 30)))
    K = 128 if D > 50 else 64
    R = 0.3 if kap > 5 else 0.4
    out = []
    for Nn in (N // 2, N):
        M0, MS, MSS, B = G.mats(c, th, Nn, D, R=R, K=K)
        lam, v = G.best_v(M0, MS, MSS, lo=0.5, hi=2.5)
        out.append((lam, v))
    # multiplication-part (essential) ratio: rho0(0)/min_X rho1(X,0) at best v
    g0, g1, g2 = [g[0, 0].real for g in G.taylor2(c, 2, R=R, K=K)]
    X = np.linspace(0, th, 201)
    v = out[-1][1]
    ess = max(g0 / ((v - X) ** 2 * g0 - 2 * (v - X) * g1 + g2))
    vs = np.linspace(0.5, 2.5, 2001)
    ess_inf = min(max(g0 / ((vv - X) ** 2 * g0 - 2 * (vv - X) * g1 + g2)) for vv in vs[:1]) if False else g0 / (g2 - g1 ** 2 / g0)
    print('kappa=%.4f D=%d: N=%d lam/k^2=%.6f (v=%.4f) | N=%d lam/k^2=%.6f (v=%.4f) | lam=%.4f ; ess ratio(v)=%.4f, 1/Var_w=%.4f [%.0fs]' % (
        kap, D, N // 2, out[0][0] / kap ** 2, out[0][1], N, out[1][0] / kap ** 2, out[1][1], out[1][0], ess, ess_inf, time.time() - t0), flush=True)
print('peak %.0f MB' % (resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 2 ** 20))
