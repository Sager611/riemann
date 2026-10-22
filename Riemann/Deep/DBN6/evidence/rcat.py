# Kernel-catalogue lever on Raymond's H(7/8) chain (non-rigorous): dense rho grid, all 10 shifts with floors.
import numpy as np, importlib.util
spec = importlib.util.spec_from_file_location('rm', 'rmodel_lib.py'); rm = importlib.util.module_from_spec(spec); spec.loader.exec_module(rm)
from fractions import Fraction as Fr
shifts = sorted({float(Fr(k[2:].replace('_','/'))) for k in rm.cells[0][0][1].keys()})
rhos = np.round(np.arange(0.05, 0.4601, 0.005), 4)
cat = []
for r in shifts:
    for rho in rhos:
        lam = rm.lam_table(r, rho, [0])[0]
        w = np.linspace(0, rho*rho, 2000); y = np.sqrt(w)
        M = float(np.max(rm.K(1, r, y)))
        cat.append((r, float(rho), 0.9995*min(lam, 2.0), M*1.000001))
print('catalogue size', len(cat), 'shifts', shifts)
orig = rm.kern
rm.kern = cat
v1 = rm.run(False)
rm.kern = orig + cat
v2 = rm.run(False)
print('dense catalogue only: %.12f ; union with Raymond 57: %.12f ; Raymond: 0.116397726408291' % (v1, v2))
