"""Lattice-node patch LP (feasibility, float64) for the c = 47 forced nodes beyond K2.
Lattice P = m L (m = 3, 4): nodes t_k = 2 pi k / P; guard band L <= |x| <= P/2 (d = m/2 - 1).
Thresholds: sublattice k = 0 mod m keeps exact weights (needs mu >= 0, Chekmachev's even rule); the other m-1 sublattices
are thresholded at tau = m deff/(m-1) with deff = 1/2 (so the default positive part is ||f||^2/2 as in Chekmachev).
Patch near each forced node t_j: h(t) = sum_i [a_i cos((t-t_j)x_i) + b_i sin((t-t_j)x_i)] D sinc^2((t-t_j)D/2),
x_i = L + i (P/2 - L)/(M+1); this is the FT of the real even kappa(x) = sum_i hat_i(|x|) [a_i cos(t_j x) + b_i cos(t_j x + pi/2 sgn x)],
supported in the guard band.  LP: maximize z = min over window nodes of (mu_k + h_k - tau_k), z <= zcap, |h_k| <= eta on the
window's outer zone.  Usage: patchlat.py m start stop [M] [eta] [W] [deff]  -> json lines"""
import math, sys, json, time, resource
import numpy as np
from scipy.optimize import linprog
sys.path.insert(0, '.')
from wsym import mu

c = 47; L = math.log(c)
m = int(sys.argv[1]); i0 = int(sys.argv[2]); i1 = int(sys.argv[3])
M = int(sys.argv[4]) if len(sys.argv) > 4 else 12
eta = float(sys.argv[5]) if len(sys.argv) > 5 else 0.5
W = float(sys.argv[6]) if len(sys.argv) > 6 else 25.0
P = m * L; d = m / 2 - 1
DEFF = float(sys.argv[7]) if len(sys.argv) > 7 else 0.5
tau_sub = m * DEFF / (m - 1)   # deff = (1/m) sum_r tau_r
xs = L + np.arange(1, M + 1) * (P / 2 - L) / (M + 1)
D = (P / 2 - L) / (M + 1)
G, NMAX, ZCAP = 25.0, 8, 1.0

def env(u):
    v = u * D / 2
    out = np.ones_like(v)
    nz = np.abs(v) > 1e-12
    out[nz] = (np.sin(v[nz]) / v[nz]) ** 2
    return D * out

F = np.load('/tmp/e27_forced47.npy')
t_all = np.sort(F[0])
cl = []; cur = [t_all[0]]
for t in t_all[1:]:
    if t - cur[-1] < G:
        cur.append(t)
    else:
        cl.append(cur); cur = [t]
cl.append(cur)
t0 = time.time()
for ci in range(i0, min(i1, len(cl))):
    C = np.array(cl[ci])
    ka = int(math.ceil((C[0] - W) * P / (2 * math.pi))); kb = int(math.floor((C[-1] + W) * P / (2 * math.pi)))
    k = np.arange(ka, kb + 1)
    t = 2 * math.pi * k / P
    mk = mu(t, c)
    tau = np.where(k % m == 0, 0.0, tau_sub)
    nforced = int((mk < tau).sum())
    rec = {'ci': ci, 't0': float(C[0]), 'n2L': int(len(C)), 'nP': int(len(k)), 'forcedP': nforced, 'minmu': float(mk.min())}
    if nforced == 0:
        rec['z'] = float((mk - tau).min()); print(json.dumps(rec)); continue
    if len(C) > NMAX:
        rec['z'] = None; rec['skip'] = 'cluster'; print(json.dumps(rec)); continue
    cols = []
    for tj in C:
        u = t - tj
        e = env(u)
        cols.append(e[:, None] * np.cos(np.outer(u, xs)))
        cols.append(e[:, None] * np.sin(np.outer(u, xs)))
    Bm = np.hstack(cols); n = Bm.shape[1]
    outer = np.min(np.abs(t[:, None] - C[None, :]), axis=1) > W - 6
    A1 = np.hstack([-Bm, np.ones((len(t), 1))])
    A2 = np.hstack([Bm[outer], np.zeros((outer.sum(), 1))]); A3 = np.hstack([-Bm[outer], np.zeros((outer.sum(), 1))])
    A_ub = np.vstack([A1, A2, A3]); b_ub = np.concatenate([mk - tau, np.full(outer.sum(), eta), np.full(outer.sum(), eta)])
    cobj = np.zeros(n + 1); cobj[-1] = -1.0
    r = linprog(cobj, A_ub=A_ub, b_ub=b_ub, bounds=[(-1e3, 1e3)] * n + [(None, ZCAP)], method='highs')
    rec['z'] = float(r.x[-1]) if r.status == 0 else None
    if r.status == 0:
        rec['hmax'] = float(np.max(np.abs(Bm @ r.x[:n])))
    print(json.dumps(rec))
sys.stderr.write(f'done {i0}-{i1} in {time.time()-t0:.1f}s maxrss {resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1e6:.0f} MB\n')
