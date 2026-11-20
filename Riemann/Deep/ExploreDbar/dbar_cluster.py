# explore59 dbar_cluster.py : clusters of size m=3 (at gamma_1) and m=2 (at gamma_5), internal spacing h -> cost -> m^2/|rho|^2
import numpy as np
from dbar_U import U, zz
g = zz[:400].copy()
w = lambda x: 1.0 / (0.25 + x * x)
base = float(np.sum(2 * w(g)))
pred = base + 2 * (9 - 1) * w(g[0]) + 2 * (4 - 1) * w(g[4])
print(f'C_K={base:.6f}  predicted with m=3 at g1, m=2 at g5: {pred:.6f}')
for h in (1e-2, 1e-3):
    for L in (20.0, 40.0, 80.0):
        gg = np.concatenate([[g[0] - h, g[0], g[0] + h, g[4] - h / 2, g[4] + h / 2], np.delete(g, [0, 4])])
        Zs = np.concatenate([gg, -gg])
        u, c, cond = U(Zs, L, True)
        print(f'h={h:g} L={L:5.1f}  L*U={L*u:.6f}  ratio to prediction {L*u/pred:.5f}  cond {cond:.1e}')
