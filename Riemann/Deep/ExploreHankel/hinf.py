# explore13 T1 test: the sup-norm (AAK/Nehari-type) Baez-Duarte problem
#   e_N = min_{c in R^N} sup_tau |1 + F(s) sum_k c_k k^{-s}| / |s|,  s = 1/2 + i tau,
# for F = zeta and F = Davenport-Heilbronn f.  At every zero rho on the line the error is 1/|rho|, so e_N >= 1/|rho_1|.
# LP: modulus by a 32-gon (true value <= LP value / cos(pi/32)); cutting planes on a fine grid tau in [0, TMAX].
# Usage: python3 hinf.py zeta|dh N1,N2,... TMAX STEP
import sys, math
import numpy as np
import mpmath as mp
from scipy.optimize import linprog
mp.mp.dps = 15
which = sys.argv[1]; Ns = [int(x) for x in sys.argv[2].split(',')]
TMAX = float(sys.argv[3]); STEP = float(sys.argv[4])
kap = (mp.sqrt(10 - 2 * mp.sqrt(5)) - 2) / (mp.sqrt(5) - 1)
A = {1: mp.mpf(1), 2: kap, 3: -kap, 4: mp.mpf(-1)}
def F(s):
    if which == 'zeta':
        return mp.zeta(s)
    return mp.power(5, -s) * sum(A[r] * mp.zeta(s, mp.mpf(r) / 5) for r in A)
taus = np.arange(0.0, TMAX + 1e-9, STEP)
cache = f"/Users/adriansager/playground/riemann/.scratch/explore13/F_{which}_{TMAX}_{STEP}.npy"
try:
    Fv = np.load(cache)
except FileNotFoundError:
    Fv = np.array([complex(F(mp.mpc(0.5, float(t)))) for t in taus]); np.save(cache, Fv)
s = 0.5 + 1j * taus; w = 1.0 / np.abs(s)
M = 32; phis = 2 * np.pi * np.arange(M) / M
def err(c, N):
    k = np.arange(1, N + 1)
    D = np.exp(-np.outer(s, np.log(k))) @ c
    return np.abs(1 + Fv * D) * w
for N in Ns:
    k = np.arange(1, N + 1)
    Z = Fv[:, None] * np.exp(-np.outer(s, np.log(k)))  # (J, N)
    act = list(range(0, len(taus), max(1, int(0.25 / STEP))))
    for it in range(30):
        rows = []; rhs = []
        Za = Z[act]; wa = w[act]
        for ph in phis:
            e = np.exp(1j * ph)
            rows.append(np.hstack([np.real(e * Za), -1.0 / wa[:, None]])); rhs.append(-np.cos(ph) * np.ones(len(act)))
        Aub = np.vstack(rows); bub = np.concatenate(rhs)
        cost = np.zeros(N + 1); cost[-1] = 1.0
        res = linprog(cost, A_ub=Aub, b_ub=bub, bounds=[(None, None)] * N + [(0, None)], method='highs')
        c = res.x[:N]; tlp = res.x[-1]
        e = err(c, N); jmax = int(np.argmax(e))
        if e[jmax] <= tlp / math.cos(math.pi / M) * (1 + 1e-9):
            break
        # add the worst violators
        viol = np.where(e > tlp / math.cos(math.pi / M))[0]
        new = viol[np.argsort(-e[viol])][:200]
        act = sorted(set(act) | set(new.tolist()))
    print(f"{which} N={N:3d}  LP={tlp:.5f}  true sup on grid={e.max():.5f} at tau={taus[jmax]:.3f}  iters={it+1}  |c|_1={np.abs(c).sum():.2f}", flush=True)
