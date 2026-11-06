#!/usr/bin/env python3
"""lincomb3 numerics step 1: Hardy Z-functions Z_j(t) of three odd primitive characters mod 7 on a grid.
Z_j(t) = e^{i theta_Arb(t,chi_j)} L(1/2+it, chi_j) (python-flint/Arb), real. Grid t = T0 + k*H, k < K.
Run: OMP_NUM_THREADS=1 perl -e 'alarm 1800; exec @ARGV' uv run --no-project --with python-flint --with numpy python3 zgrid.py <n> [Tmax]
"""
import sys, time, numpy as np
from flint import dirichlet_char, ctx
ctx.prec = 45
q = 7
n = int(sys.argv[1]); Tmax = float(sys.argv[2]) if len(sys.argv) > 2 else 2400.0
T0, H = 5.0, 0.03
K = int((Tmax - T0) / H) + 1
chi = dirichlet_char(q, n)
assert chi.is_primitive() and chi.parity() == 1
Z = np.empty(K)
t0 = time.time()
for k in range(K):
    Z[k] = float(chi.hardy_z(T0 + k * H).real.mid())
    if k % 20000 == 0:
        print(n, k, K, round(time.time() - t0, 1), flush=True)
np.save('Z_q7_n%d.npy' % n, Z)
print('done', n, K, round(time.time() - t0, 1))
