# pyright: reportMissingImports=false
# python-flint lives in a local venv (see README.md), not the system interpreter.

import sys, time, math
from flint import ctx
import weil
ctx.prec = 1000
bands = [
 (80,  ['0.5','0.6','0.683147','log2','0.703147','0.8','0.9','1.0','1.088612','log3','1.108612','1.2','1.3','1.376294','log4','1.396294']),
 (140, ['1.5','log5','1.7','1.8','1.9','log7','2.0']),
 (200, ['log8','log9','2.25','log11','2.5']),
 (240, ['log13','2.75','3.0']),
]
out = open('table.txt', 'a')
for N, Ls in bands:
    for s in Ls:
        t0 = time.time(); L = weil.parseL(s)
        M = weil.weil_matrix(L, N)
        le = weil.lam_min(weil.sub(M, list(range(0, N, 2)))); lo = weil.lam_min(weil.sub(M, list(range(1, N, 2))))
        line = f"{s}\t{float(L):.6f}\tN={N}\teven={le.str(6, radius=False)}\todd={lo.str(6, radius=False)}\t{time.time()-t0:.0f}s"
        print(line, flush=True); out.write(line + "\n"); out.flush()
